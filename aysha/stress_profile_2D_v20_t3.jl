# ============================================================================
# Gate-free T3 observer stress profile for a fixed stress-test plant.
# ============================================================================

using Statistics

include("profile_2D_v20_t3_observer.jl")

const STRESS_T3_CAPACITIES_2D_v20 = (
    1200.0, 3000.0, 4500.0, 6000.0, 10000.0, 30000.0,
)
const STRESS_T3_STEMS_2D_v20 = (0.0, 20.0)
const STRESS_T3_NU50_REFERENCE_2D_v20 =
    3.07786e-4 * 50.0^1.44346

function stress_profile_2D_v20_t3(
    exponent, ratio, felt_k, felt_cp,
    power456, power304, power256;
    max_points=61, suffix="",
)
    mkpath(OUTPUT_DIR_2D_v20)
    nu50 = ratio * STRESS_T3_NU50_REFERENCE_2D_v20
    p = parameters_2D_v20(
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
        power_scales=(power456, power304, power256),
        t3_location=:receiver_136,
    )
    cases = Dict{String,Any}()
    for id in ("C69", "C80", "C81")
        println("v20 gate-free T3 plant solve: $id")
        flush(stdout)
        input = case_inputs_2D_v20(
            id; cooling=true, max_points=max_points,
        )
        cases[id] = simulate_case_2D_v20(
            input, p; initialization=:side_T2_only,
            reltol=5e-4, abstol=1e-4, dtmax=120.0,
        )
    end
    rows = NamedTuple[]
    for location in T3_LOCATIONS_2D_v20,
        capacity in STRESS_T3_CAPACITIES_2D_v20,
        stem in STRESS_T3_STEMS_2D_v20
        training_scores = [
            _observer_case_score_2D_v20(
                cases[id],
                _observer_prediction_2D_v20(
                    cases[id], location, capacity, stem,
                ).T3,
            ) for id in ("C69", "C80")
        ]
        train = _aggregate_observer_score_2D_v20(
            training_scores,
        )
        validation_prediction = _observer_prediction_2D_v20(
            cases["C81"], location, capacity, stem,
        )
        validation = _observer_case_score_2D_v20(
            cases["C81"], validation_prediction.T3,
        )
        push!(rows, (
            reynolds_exponent=exponent,
            Nu50_ratio=ratio,
            felt_conductivity_scale=felt_k,
            felt_heat_capacity_scale=felt_cp,
            power_scale_456=power456,
            power_scale_304=power304,
            power_scale_256=power256,
            location=String(location),
            global_z_mm=1000.0 *
                validation_prediction.T3_global_z_m,
            capacity_areal_J_m2_K=capacity,
            stem_conductance_areal_W_m2_K=stem,
            training_objective=train.objective,
            training_T3_rmse_K=train.T3_rmse_K,
            training_max_T3_rmse_K=train.max_T3_rmse_K,
            training_final_bias_K=train.final_bias_K,
            training_t90_mae_s=train.t90_mae_s,
            validation_T3_rmse_K=validation.T3_rmse_K,
            validation_final_bias_K=validation.final_bias_K,
            validation_t90_error_s=validation.t90_error_s,
        ))
    end
    sort!(rows; by=row -> row.training_objective)
    output = joinpath(
        OUTPUT_DIR_2D_v20,
        "stress_t3_profile$(suffix)_2D_v20.csv",
    )
    _write_namedtuples_csv_2D_v20(output, rows)
    winner = first(rows)
    _write_namedtuples_csv_2D_v20(joinpath(
        OUTPUT_DIR_2D_v20,
        "stress_t3_selected$(suffix)_2D_v20.csv",
    ), [winner])
    println("v20 gate-free T3 best: ", winner)
    return (winner=winner, rows=rows, cases=cases)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 7 || error(
        "usage: exponent ratio felt_k felt_cp p456 p304 p256 [points] [suffix]",
    )
    exponent = parse(Float64, ARGS[1])
    ratio = parse(Float64, ARGS[2])
    felt_k = parse(Float64, ARGS[3])
    felt_cp = parse(Float64, ARGS[4])
    p456 = parse(Float64, ARGS[5])
    p304 = parse(Float64, ARGS[6])
    p256 = parse(Float64, ARGS[7])
    points = length(ARGS) >= 8 ? parse(Int, ARGS[8]) : 61
    suffix = length(ARGS) >= 9 ? ARGS[9] : ""
    stress_profile_2D_v20_t3(
        exponent, ratio, felt_k, felt_cp,
        p456, p304, p256;
        max_points=points, suffix=suffix,
    )
end
