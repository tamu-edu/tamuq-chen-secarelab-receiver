# ============================================================================
# Gate-free v20 nuisance stress screen at one structural cell.
#
# Deliberately extends below the physical identification brackets. The
# selected values are diagnostic compensators, never validated parameters.
# ============================================================================

using Statistics

include("profile_2D_v20_ua_plant.jl")

const STRESS_FELT_K_2D_v20 = (0.40, 0.55, 0.70)
const STRESS_FELT_CP_2D_v20 = (0.40, 0.55, 0.75)
const STRESS_POWER_GROUPS_2D_v20 = (
    (name="456", id="E69", index=1, levels=(0.85, 0.95, 1.05)),
    (name="304", id="E74", index=2, levels=(0.90, 1.05, 1.23)),
    (name="256", id="E79", index=3, levels=(0.55, 0.70, 0.84)),
)

function _stress_parameters_2D_v20(
    exponent, ratio;
    felt_k=0.70, felt_cp=0.75,
    powers=(1.05, 1.23, 0.84),
)
    nu50 = ratio * UA_PROFILE_NU50_REFERENCE_2D_v20
    return parameters_2D_v20(
        mesh=:nominal,
        source_model=:near_deep,
        deep_fraction=0.90,
        deep_length_m=0.12,
        nu_prefactor=nu50 / 50.0^exponent,
        reynolds_exponent=exponent,
        distributed_tube_flange_h_W_m2_K=0.0,
        probe_capacity_areal_J_m2_K=1200.0,
        probe_stem_conductance_areal_W_m2_K=0.0,
        felt_conductivity_scale=felt_k,
        felt_heat_capacity_scale=felt_cp,
        felt_contact_scale=0.30,
        power_scales=powers,
        t3_location=:receiver_136,
    )
end

function _stress_cooling_score_2D_v20(case)
    rows = 2:size(case.model, 1)
    threshold = first(case.inputs.times) +
        0.9 * (last(case.inputs.times) - first(case.inputs.times))
    start = something(
        findfirst(time -> time >= threshold, case.inputs.times),
        size(case.model, 1),
    )
    final = max(2, start):size(case.model, 1)
    side_error = case.model[rows, 1:3] .-
                 case.observed[rows, 1:3]
    T2_error = case.model[rows, 7] .-
               case.observed[rows, 7]
    side_curve = mean(
        _profile_rho_2D_v20(value / 35.0)
        for value in side_error
    )
    T2_curve = mean(
        _profile_rho_2D_v20(value / 15.0)
        for value in T2_error
    )
    side_level = mean(
        _profile_rho_2D_v20((
            mean(case.model[final, sensor]) -
            mean(case.observed[final, sensor])
        ) / 30.0) for sensor in 1:3
    )
    T2_level = _profile_rho_2D_v20((
        mean(case.model[final, 7]) -
        mean(case.observed[final, 7])
    ) / 10.0)
    return (
        objective=0.70 * (
            0.60side_curve + 0.40T2_curve
        ) + 0.30 * (
            0.60side_level + 0.40T2_level
        ),
        side_rmse_K=sqrt(mean(abs2, side_error)),
        side_final_bias_K=mean(
            case.model[final, 1:3] .-
            case.observed[final, 1:3],
        ),
        T2_rmse_K=sqrt(mean(abs2, T2_error)),
        T2_final_bias_K=mean(
            case.model[final, 7] .-
            case.observed[final, 7],
        ),
        T3_diagnostic_rmse_K=sqrt(mean(abs2,
            case.model[rows, 6] .-
            case.observed[rows, 6],
        )),
    )
end

function stress_screen_2D_v20_nuisance(
    exponent, ratio; max_points=31, suffix="",
)
    mkpath(OUTPUT_DIR_2D_v20)
    cooling_inputs = [
        case_inputs_2D_v20(
            id; cooling=true, max_points=max_points,
        ) for id in ("C69", "C80")
    ]
    felt_rows = NamedTuple[]
    felt_candidates = NamedTuple[]
    for felt_k in STRESS_FELT_K_2D_v20,
        felt_cp in STRESS_FELT_CP_2D_v20
        println("v20 gate-free felt k=$felt_k cp=$felt_cp")
        flush(stdout)
        p = _stress_parameters_2D_v20(
            exponent, ratio; felt_k=felt_k, felt_cp=felt_cp,
        )
        cases = [
            simulate_case_2D_v20(
                input, p; initialization=:side_T2_only,
                reltol=5e-4, abstol=1e-4, dtmax=120.0,
            ) for input in cooling_inputs
        ]
        score = _aggregate_plant_score_2D_v20(
            _stress_cooling_score_2D_v20.(cases),
        )
        flow = all(
            all(case.result.flow_solver_converged)
            for case in cases
        )
        row = merge((
            reynolds_exponent=exponent,
            Nu50_ratio=ratio,
            felt_conductivity_scale=felt_k,
            felt_heat_capacity_scale=felt_cp,
            flow_converged=flow,
        ), score)
        push!(felt_rows, row)
        push!(felt_candidates, (
            felt_k=felt_k, felt_cp=felt_cp,
            score=score, flow=flow, cases=cases,
        ))
    end
    sort!(felt_rows; by=row -> row.objective)
    felt = first(sort(
        felt_candidates; by=item -> (
            item.flow ? 0 : 1, item.score.objective,
        ),
    ))
    _write_namedtuples_csv_2D_v20(joinpath(
        OUTPUT_DIR_2D_v20,
        "stress_nuisance_felt$(suffix)_2D_v20.csv",
    ), felt_rows)

    selected_powers = [1.05, 1.23, 0.84]
    heating_cases = Any[]
    power_rows = NamedTuple[]
    for group in STRESS_POWER_GROUPS_2D_v20
        input = case_inputs_2D_v20(
            group.id; max_points=max_points,
        )
        candidates = NamedTuple[]
        for scale in group.levels
            powers = copy(selected_powers)
            powers[group.index] = scale
            println(
                "v20 gate-free power $(group.name)=$scale",
            )
            flush(stdout)
            p = _stress_parameters_2D_v20(
                exponent, ratio;
                felt_k=felt.felt_k,
                felt_cp=felt.felt_cp,
                powers=Tuple(powers),
            )
            case = simulate_case_2D_v20(
                input, p; reltol=5e-4, abstol=1e-4,
                dtmax=120.0,
            )
            score = _plant_case_score_2D_v20(case)
            flow = all(case.result.flow_solver_converged)
            row = merge((
                reynolds_exponent=exponent,
                Nu50_ratio=ratio,
                power_group=group.name,
                power_scale=scale,
                flow_converged=flow,
            ), score)
            push!(power_rows, row)
            push!(candidates, (
                scale=scale, score=score,
                flow=flow, case=case,
            ))
        end
        winner = first(sort(
            candidates; by=item -> (
                item.flow ? 0 : 1, item.score.objective,
            ),
        ))
        selected_powers[group.index] = winner.scale
        push!(heating_cases, winner.case)
    end
    _write_namedtuples_csv_2D_v20(joinpath(
        OUTPUT_DIR_2D_v20,
        "stress_nuisance_power$(suffix)_2D_v20.csv",
    ), power_rows)
    heating = _aggregate_plant_score_2D_v20(
        _plant_case_score_2D_v20.(heating_cases),
    )
    cooling = felt.score
    summary = (
        reynolds_exponent=exponent,
        Nu50_ratio=ratio,
        felt_conductivity_scale=felt.felt_k,
        felt_heat_capacity_scale=felt.felt_cp,
        power_scale_456=selected_powers[1],
        power_scale_304=selected_powers[2],
        power_scale_256=selected_powers[3],
        heating_objective=heating.objective,
        heating_side_rmse_K=heating.side_rmse_K,
        heating_T2_rmse_K=heating.T2_rmse_K,
        heating_T3_diagnostic_rmse_K=
            heating.T3_diagnostic_rmse_K,
        cooling_objective=cooling.objective,
        cooling_side_rmse_K=cooling.side_rmse_K,
        cooling_T2_rmse_K=cooling.T2_rmse_K,
        cooling_T3_diagnostic_rmse_K=
            cooling.T3_diagnostic_rmse_K,
        flow_converged=(
            felt.flow && all(
                all(case.result.flow_solver_converged)
                for case in heating_cases
            )
        ),
        original_heating_gate=(
            heating.side_rmse_K < 60.0 &&
            heating.T2_rmse_K < 15.0
        ),
        original_cooling_gate=(
            cooling.side_rmse_K < 35.0 &&
            cooling.T2_rmse_K < 12.0
        ),
    )
    _write_namedtuples_csv_2D_v20(joinpath(
        OUTPUT_DIR_2D_v20,
        "stress_nuisance_selected$(suffix)_2D_v20.csv",
    ), [summary])
    println("v20 gate-free nuisance best: ", summary)
    return (
        summary=summary, felt_rows=felt_rows,
        power_rows=power_rows,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    exponent = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 0.50
    ratio = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 2.00
    points = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 31
    suffix = length(ARGS) >= 4 ? ARGS[4] : ""
    stress_screen_2D_v20_nuisance(
        exponent, ratio; max_points=points, suffix=suffix,
    )
end
