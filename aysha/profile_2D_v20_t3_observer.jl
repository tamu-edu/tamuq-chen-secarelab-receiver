# ============================================================================
# V20 discrete T3-location and probe-response profile.
#
# The plant trajectory is solved once per cooling experiment without using
# T3. Location, probe capacity, and stem conductance are post-processed.
# ============================================================================

using Statistics

include("run_2D_v20.jl")

const T3_LOCATIONS_2D_v20 = (
    :receiver_136, :receiver_exit, :rear_003,
)
const T3_CAPACITIES_2D_v20 = (
    200.0, 600.0, 1200.0, 2400.0, 3000.0,
)
const T3_STEMS_2D_v20 = (
    0.0, 20.0, 60.0, 120.0, 200.0,
)

function _observer_parameters_2D_v20(
    old::V20.ModelParameters2D, location, capacity, stem,
)
    old19 = old.base
    old_probe = old19.outlet_observation
    probe = V20.OutletObservationParameters2D(
        capacity_areal_J_m2_K=capacity,
        stem_conductance_areal_W_m2_K=stem,
        emissivity=old_probe.emissivity,
        probe_diameter_m=old_probe.probe_diameter_m,
        axial_position_m=old_probe.axial_position_m,
    )
    base = V20.V19.ModelParameters2D(
        base=old19.base,
        integrated_exchange=old19.integrated_exchange,
        outlet_observation=probe,
        rear_tube_flange=old19.rear_tube_flange,
    )
    return V20.ModelParameters2D(
        base=base,
        t3_location=V20.T3LocationParameters2D(
            location=location,
        ),
    )
end

function _observer_prediction_2D_v20(
    case, location, capacity, stem,
)
    parameters = _observer_parameters_2D_v20(
        case.result.parameters, location, capacity, stem,
    )
    result = V20.SimulationResult2D(
        case.result.base_result, parameters,
        case.result.diagnostics, case.result.ode_solution,
    )
    return V20.sensor_predictions2D(
        result; initial_T3=case.observed[1, 6],
    )
end

_rho_2D_v20(value) =
    2.0 * (sqrt(1.0 + Float64(value)^2) - 1.0)

function _observer_case_score_2D_v20(
    case, prediction,
)
    observed = case.observed[:, 6]
    rows = 2:length(observed)
    final = max(2, floor(Int, 0.9length(observed))):length(observed)
    error = prediction .- observed
    curve = mean(
        _rho_2D_v20(error[index] / 25.0)
        for index in rows
    )
    final_bias = mean(error[final])
    level = _rho_2D_v20(final_bias / 20.0)
    return (
        objective=0.70curve + 0.30level,
        T3_rmse_K=sqrt(mean(abs2, error[rows])),
        final_bias_K=final_bias,
        t90_error_s=V20.get_t90_2D(
            case.inputs.times, prediction,
        ) - V20.get_t90_2D(
            case.inputs.times, observed,
        ),
    )
end

function _aggregate_observer_score_2D_v20(scores)
    return (
        objective=mean(score.objective for score in scores),
        T3_rmse_K=mean(score.T3_rmse_K for score in scores),
        max_T3_rmse_K=maximum(
            score.T3_rmse_K for score in scores
        ),
        final_bias_K=mean(
            score.final_bias_K for score in scores
        ),
        max_abs_final_bias_K=maximum(
            abs(score.final_bias_K) for score in scores
        ),
        t90_mae_s=mean(
            abs(score.t90_error_s) for score in scores
        ),
    )
end

function profile_2D_v20_t3_observer(
    ; mesh=:nominal, max_points=61,
)
    mkpath(OUTPUT_DIR_2D_v20)
    p = parameters_2D_v20(
        mesh=mesh,
        source_model=:near_deep,
        deep_fraction=0.90,
        deep_length_m=0.12,
        nu_prefactor=3.80e-4,
        reynolds_exponent=1.54,
        distributed_tube_flange_h_W_m2_K=0.0,
        probe_capacity_areal_J_m2_K=1000.0,
        probe_stem_conductance_areal_W_m2_K=20.0,
        felt_conductivity_scale=0.70,
        felt_heat_capacity_scale=0.75,
        felt_contact_scale=0.30,
        power_scales=(1.05, 1.23, 0.84),
    )
    cases = Dict{String,Any}()
    for id in ("C69", "C80", "C81")
        println("v20 T3 plant solve: $id")
        flush(stdout)
        inputs = case_inputs_2D_v20(
            id; cooling=true, max_points=max_points,
        )
        cases[id] = simulate_case_2D_v20(
            inputs, p; initialization=:side_T2_only,
            reltol=5e-4, abstol=1e-4, dtmax=120.0,
        )
    end

    rows = NamedTuple[]
    candidates = NamedTuple[]
    for location in T3_LOCATIONS_2D_v20,
        capacity in T3_CAPACITIES_2D_v20,
        stem in T3_STEMS_2D_v20
        train_scores = [
            _observer_case_score_2D_v20(
                cases[id],
                _observer_prediction_2D_v20(
                    cases[id], location, capacity, stem,
                ).T3,
            ) for id in ("C69", "C80")
        ]
        validation_prediction = _observer_prediction_2D_v20(
            cases["C81"], location, capacity, stem,
        )
        validation = _observer_case_score_2D_v20(
            cases["C81"], validation_prediction.T3,
        )
        train = _aggregate_observer_score_2D_v20(train_scores)
        capacity_interior =
            capacity != first(T3_CAPACITIES_2D_v20) &&
            capacity != last(T3_CAPACITIES_2D_v20)
        stem_interior =
            stem != first(T3_STEMS_2D_v20) &&
            stem != last(T3_STEMS_2D_v20)
        gate = train.T3_rmse_K <= 30.0 &&
               validation.T3_rmse_K <= 30.0 &&
               abs(validation.final_bias_K) <= 20.0 &&
               abs(validation.t90_error_s) <= 500.0 &&
               capacity_interior && stem_interior
        row = (
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
            capacity_interior=capacity_interior,
            stem_interior=stem_interior,
            observer_gate=gate,
        )
        push!(rows, row)
        push!(candidates, (
            location=location, capacity=capacity, stem=stem,
            score=train.objective, row=row,
        ))
    end

    null_rows = NamedTuple[]
    for location in T3_LOCATIONS_2D_v20
        train_scores = NamedTuple[]
        validation_score = nothing
        for id in ("C69", "C80", "C81")
            prediction = _observer_prediction_2D_v20(
                cases[id], location, 1000.0, 20.0,
            )
            score = _observer_case_score_2D_v20(
                cases[id], prediction.T3_gas_raw,
            )
            if id == "C81"
                validation_score = score
            else
                push!(train_scores, score)
            end
        end
        train = _aggregate_observer_score_2D_v20(train_scores)
        push!(null_rows, (
            location=String(location),
            training_T3_rmse_K=train.T3_rmse_K,
            training_final_bias_K=train.final_bias_K,
            validation_T3_rmse_K=
                validation_score.T3_rmse_K,
            validation_final_bias_K=
                validation_score.final_bias_K,
            validation_t90_error_s=
                validation_score.t90_error_s,
        ))
    end
    sort!(rows; by=row -> row.training_objective)
    sort!(null_rows; by=row -> row.training_T3_rmse_K)
    _write_namedtuples_csv_2D_v20(joinpath(
        OUTPUT_DIR_2D_v20,
        "t3_observer_profile_2D_v20.csv",
    ), rows)
    _write_namedtuples_csv_2D_v20(joinpath(
        OUTPUT_DIR_2D_v20,
        "t3_location_raw_gas_null_2D_v20.csv",
    ), null_rows)
    winner = first(sort(
        candidates; by=item -> item.score,
    ))
    open(joinpath(
        OUTPUT_DIR_2D_v20,
        "t3_observer_selected_2D_v20.txt",
    ), "w") do io
        println(io, "location=", winner.location)
        println(io, "capacity_areal_J_m2_K=", winner.capacity)
        println(io, "stem_conductance_areal_W_m2_K=", winner.stem)
        for name in propertynames(winner.row)
            println(io, name, "=", getproperty(winner.row, name))
        end
        println(io, "any_observer_gate=",
                any(row.observer_gate for row in rows))
    end
    println("v20 T3 observer winner: ", winner.row)
    println("v20 any observer gate: ",
            any(row.observer_gate for row in rows))
    return (
        winner=winner, rows=rows, null_rows=null_rows,
        cases=cases,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    mesh = length(ARGS) >= 1 ? Symbol(ARGS[1]) : :nominal
    points = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 61
    profile_2D_v20_t3_observer(
        ; mesh=mesh, max_points=points,
    )
end
