# ============================================================================
# Bounded T3-free nuisance screen at one structural (source, Nu50, n) cell.
#
# Cooling selects felt k/Cp from side and T2 only. Each lamp group then
# selects its power scale from one representative heating case. T3, T9, and
# T10 are diagnostic only.
# ============================================================================

using Statistics

include("profile_2D_v20_ua_plant.jl")

const FELT_K_LEVELS_2D_v20 = (0.70, 1.00, 1.30)
const FELT_CP_LEVELS_2D_v20 = (0.75, 1.00, 1.50)
const POWER_LEVELS_2D_v20 = (
    (name="456", id="E69", index=1,
     levels=(1.05, 1.195, 1.34)),
    (name="304", id="E74", index=2,
     levels=(1.23, 1.405, 1.58)),
    (name="256", id="E79", index=3,
     levels=(0.84, 0.975, 1.11)),
)

function _cooling_case_score_2D_v20(case)
    rows = 2:size(case.model, 1)
    final_start = findfirst(
        time -> time >= first(case.inputs.times) +
                0.9 * (
                    last(case.inputs.times) -
                    first(case.inputs.times)
                ),
        case.inputs.times,
    )
    final = max(2, something(
        final_start, size(case.model, 1),
    )):size(case.model, 1)
    side_error = case.model[rows, 1:3] .-
                 case.observed[rows, 1:3]
    T2_error = case.model[rows, 7] .-
               case.observed[rows, 7]
    side_curve = mean(
        _profile_rho_2D_v20(value / 35.0)
        for value in side_error
    )
    side_level = mean(
        _profile_rho_2D_v20((
            mean(case.model[final, column]) -
            mean(case.observed[final, column])
        ) / 30.0) for column in 1:3
    )
    T2_curve = mean(
        _profile_rho_2D_v20(value / 15.0)
        for value in T2_error
    )
    T2_level = _profile_rho_2D_v20((
        mean(case.model[final, 7]) -
        mean(case.observed[final, 7])
    ) / 10.0)
    curve = 0.60side_curve + 0.40T2_curve
    level = 0.60side_level + 0.40T2_level
    return (
        objective=0.70curve + 0.30level,
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

function _nuisance_parameters_2D_v20(
    exponent, ratio;
    felt_k=0.70, felt_cp=0.75,
    powers=(1.05, 1.23, 0.84),
)
    nu50 = ratio * UA_PROFILE_NU50_REFERENCE_2D_v20
    prefactor = nu50 / 50.0^exponent
    return parameters_2D_v20(
        mesh=:nominal,
        source_model=:near_deep,
        deep_fraction=0.90,
        deep_length_m=0.12,
        nu_prefactor=prefactor,
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

function screen_2D_v20_plant_nuisance(
    exponent, ratio; max_points=31, suffix="",
)
    mkpath(OUTPUT_DIR_2D_v20)
    felt_rows = NamedTuple[]
    felt_results = NamedTuple[]
    cooling_inputs = [
        case_inputs_2D_v20(
            id; cooling=true, max_points=max_points,
        ) for id in ("C69", "C80")
    ]
    for felt_k in FELT_K_LEVELS_2D_v20,
        felt_cp in FELT_CP_LEVELS_2D_v20
        println("v20 nuisance felt k=$felt_k cp=$felt_cp")
        flush(stdout)
        p = _nuisance_parameters_2D_v20(
            exponent, ratio; felt_k=felt_k, felt_cp=felt_cp,
        )
        cases = [
            simulate_case_2D_v20(
                input, p; initialization=:side_T2_only,
                reltol=5e-4, abstol=1e-4,
                dtmax=120.0,
            ) for input in cooling_inputs
        ]
        score = _aggregate_plant_score_2D_v20(
            _cooling_case_score_2D_v20.(cases),
        )
        row = merge((
            reynolds_exponent=exponent,
            Nu50_ratio=ratio,
            felt_conductivity_scale=felt_k,
            felt_heat_capacity_scale=felt_cp,
        ), score)
        push!(felt_rows, row)
        flow_converged = all(
            all(case.result.flow_solver_converged)
            for case in cases
        )
        push!(felt_results, (
            felt_k=felt_k, felt_cp=felt_cp,
            score=score, cases=cases,
            flow_converged=flow_converged,
        ))
    end
    sort!(felt_rows; by=row -> row.objective)
    felt_winner = first(sort(
        felt_results; by=item -> (
            item.flow_converged ? 0 : 1,
            item.score.objective,
        ),
    ))
    _write_namedtuples_csv_2D_v20(joinpath(
        OUTPUT_DIR_2D_v20,
        "plant_nuisance_felt$(suffix)_2D_v20.csv",
    ), felt_rows)

    selected_powers = [1.05, 1.23, 0.84]
    power_rows = NamedTuple[]
    selected_heating_cases = Any[]
    power_boundary = false
    for group in POWER_LEVELS_2D_v20
        input = case_inputs_2D_v20(
            group.id; max_points=max_points,
        )
        candidates = NamedTuple[]
        for scale in group.levels
            powers = copy(selected_powers)
            powers[group.index] = scale
            println(
                "v20 nuisance power $(group.name) scale=$scale",
            )
            flush(stdout)
            p = _nuisance_parameters_2D_v20(
                exponent, ratio;
                felt_k=felt_winner.felt_k,
                felt_cp=felt_winner.felt_cp,
                powers=Tuple(powers),
            )
            case = simulate_case_2D_v20(
                input, p; reltol=5e-4, abstol=1e-4,
                dtmax=120.0,
            )
            score = _plant_case_score_2D_v20(case)
            flow_converged = all(
                case.result.flow_solver_converged,
            )
            row = merge((
                reynolds_exponent=exponent,
                Nu50_ratio=ratio,
                power_group=group.name,
                power_scale=scale,
                flow_converged=flow_converged,
            ), score)
            push!(power_rows, row)
            push!(candidates, (
                scale=scale, score=score, case=case,
                flow_converged=flow_converged,
            ))
        end
        winner = first(sort(
            candidates; by=item -> (
                item.flow_converged ? 0 : 1,
                item.score.objective,
            ),
        ))
        selected_powers[group.index] = winner.scale
        push!(selected_heating_cases, winner.case)
        power_boundary |= (
            winner.scale == first(group.levels) ||
            winner.scale == last(group.levels)
        )
    end
    sort!(power_rows; by=row -> (
        row.power_group, row.objective,
    ))
    _write_namedtuples_csv_2D_v20(joinpath(
        OUTPUT_DIR_2D_v20,
        "plant_nuisance_power$(suffix)_2D_v20.csv",
    ), power_rows)

    heating_score = _aggregate_plant_score_2D_v20(
        _plant_case_score_2D_v20.(selected_heating_cases),
    )
    cooling_score = felt_winner.score
    felt_boundary = (
        felt_winner.felt_k == first(FELT_K_LEVELS_2D_v20) ||
        felt_winner.felt_k == last(FELT_K_LEVELS_2D_v20) ||
        felt_winner.felt_cp == first(FELT_CP_LEVELS_2D_v20) ||
        felt_winner.felt_cp == last(FELT_CP_LEVELS_2D_v20)
    )
    heating_gate = heating_score.side_rmse_K < 60.0 &&
                   heating_score.T2_rmse_K < 15.0
    cooling_gate = cooling_score.side_rmse_K < 35.0 &&
                   cooling_score.T2_rmse_K < 12.0
    flow_converged = felt_winner.flow_converged &&
        all(
            all(case.result.flow_solver_converged)
            for case in selected_heating_cases
        )
    structural_boundary = exponent in (
        first(UA_PROFILE_EXPONENTS_2D_v20),
        last(UA_PROFILE_EXPONENTS_2D_v20),
    ) || ratio in (
        first(UA_PROFILE_RATIOS_2D_v20),
        last(UA_PROFILE_RATIOS_2D_v20),
    )
    summary = (
        reynolds_exponent=exponent,
        Nu50_ratio=ratio,
        felt_conductivity_scale=felt_winner.felt_k,
        felt_heat_capacity_scale=felt_winner.felt_cp,
        power_scale_456=selected_powers[1],
        power_scale_304=selected_powers[2],
        power_scale_256=selected_powers[3],
        heating_objective=heating_score.objective,
        heating_side_rmse_K=heating_score.side_rmse_K,
        heating_T2_rmse_K=heating_score.T2_rmse_K,
        heating_T3_diagnostic_rmse_K=
            heating_score.T3_diagnostic_rmse_K,
        cooling_objective=cooling_score.objective,
        cooling_side_rmse_K=cooling_score.side_rmse_K,
        cooling_T2_rmse_K=cooling_score.T2_rmse_K,
        cooling_T3_diagnostic_rmse_K=
            cooling_score.T3_diagnostic_rmse_K,
        felt_boundary=felt_boundary,
        power_boundary=power_boundary,
        structural_boundary=structural_boundary,
        flow_converged=flow_converged,
        heating_gate=heating_gate,
        cooling_gate=cooling_gate,
        plant_admissible=(
            heating_gate && cooling_gate &&
            !felt_boundary && !power_boundary &&
            !structural_boundary && flow_converged
        ),
    )
    _write_namedtuples_csv_2D_v20(joinpath(
        OUTPUT_DIR_2D_v20,
        "plant_nuisance_selected$(suffix)_2D_v20.csv",
    ), [summary])
    println("v20 nuisance selected: ", summary)
    return (
        summary=summary, felt_rows=felt_rows,
        power_rows=power_rows,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    exponent = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 1.25
    ratio = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1.35
    points = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 31
    suffix = length(ARGS) >= 4 ? ARGS[4] : ""
    screen_2D_v20_plant_nuisance(
        exponent, ratio; max_points=points, suffix=suffix,
    )
end
