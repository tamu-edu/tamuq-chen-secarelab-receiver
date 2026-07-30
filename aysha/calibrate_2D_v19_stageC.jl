# ============================================================================
# calibrate_2D_v19_stageC.jl
# Stage C: rear tube/flange topology, then physical T3 probe post-processing.
# ============================================================================

using Statistics

include("calibrate_2D_v19_stageB.jl")

function selected_ua_2D_v19()
    values = _read_keyvalues_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19, "stageB_selected_2D_v19.txt",
    ))
    return (
        prefactor=parse(Float64, values["nu_prefactor"]),
        exponent=parse(Float64, values["reynolds_exponent"]),
    )
end

function _stageC_parameters_2D_v19(
    source, ua, tube_h;
    capacity=1000.0, stem=20.0, mesh=:screen,
)
    return parameters_2D_v19(
        mesh=mesh,
        source_model=source.model,
        deep_fraction=source.fraction,
        deep_length_m=source.length_m,
        nu_prefactor=ua.prefactor,
        reynolds_exponent=ua.exponent,
        probe_capacity_areal_J_m2_K=capacity,
        probe_stem_conductance_areal_W_m2_K=stem,
        rear_tube_flange_contact_h_W_m2_K=tube_h,
        power_scales=(1.195, 1.405, 0.975),
        felt_conductivity_scale=1.0,
        felt_heat_capacity_scale=1.0,
    )
end

function cooling_plant_score_2D_v19(cases)
    side_losses = Float64[]
    felt_losses = Float64[]
    side_rmse = Float64[]
    felt_rmse = Float64[]
    for case in cases
        rows = 2:size(case.model, 1)
        side_error = case.model[rows, 1:3] .-
                     case.observed[rows, 1:3]
        felt_error = case.model[rows, 7] .-
                     case.observed[rows, 7]
        push!(side_losses, mean(
            rho_2D_v19(value / 35.0) for value in side_error
        ))
        push!(felt_losses, mean(
            rho_2D_v19(value / 15.0) for value in felt_error
        ))
        push!(side_rmse, sqrt(mean(abs2, side_error)))
        push!(felt_rmse, sqrt(mean(abs2, felt_error)))
    end
    return (
        objective=0.60 * mean(side_losses) +
                  0.40 * mean(felt_losses),
        side_rmse_K=mean(side_rmse),
        felt_rmse_K=mean(felt_rmse),
    )
end

function _replace_probe_2D_v19(result, capacity, stem)
    old = result.parameters
    probe = V19.OutletObservationParameters2D(
        capacity_areal_J_m2_K=capacity,
        stem_conductance_areal_W_m2_K=stem,
        emissivity=old.outlet_observation.emissivity,
        probe_diameter_m=
            old.outlet_observation.probe_diameter_m,
        axial_position_m=
            old.outlet_observation.axial_position_m,
    )
    parameters = V19.ModelParameters2D(
        base=old.base,
        integrated_exchange=old.integrated_exchange,
        outlet_observation=probe,
        rear_tube_flange=old.rear_tube_flange,
    )
    return V19.SimulationResult2D(
        result.base_result, parameters,
        result.diagnostics, result.ode_solution,
    )
end

function probe_score_2D_v19(cases, capacity, stem)
    losses = Float64[]
    rmses = Float64[]
    final_biases = Float64[]
    t90_errors = Float64[]
    for case in cases
        result = _replace_probe_2D_v19(
            case.result, capacity, stem,
        )
        prediction = V19.sensor_predictions2D(
            result; initial_T3=case.observed[1, 6],
        ).T3
        observed = case.observed[:, 6]
        rows = 2:length(prediction)
        final = max(2, floor(Int, 0.9 * length(prediction))):length(prediction)
        curve = mean(
            rho_2D_v19(
                (prediction[index] - observed[index]) / 25.0
            ) for index in rows
        )
        level = rho_2D_v19((
            mean(prediction[final]) - mean(observed[final])
        ) / 20.0)
        push!(losses, 0.70 * curve + 0.30 * level)
        push!(rmses, sqrt(mean(abs2, prediction .- observed)))
        push!(final_biases,
              mean(prediction[final]) - mean(observed[final]))
        push!(t90_errors,
            V19.get_t90_2D(case.inputs.times, prediction) -
            V19.get_t90_2D(case.inputs.times, observed))
    end
    return (
        objective=mean(losses),
        T3_rmse_K=mean(rmses),
        final_bias_K=mean(final_biases),
        t90_mae_s=mean(abs, t90_errors),
        max_abs_final_bias_K=maximum(abs, final_biases),
    )
end

function calibrate_2D_v19_stageC()
    source = selected_source_2D_v19()
    ua = selected_ua_2D_v19()
    training_inputs = [
        case_inputs_2D_v19(
            id; cooling=true, max_points=61,
        ) for id in ("C69", "C80")
    ]
    plant_rows = NamedTuple[]
    plant_results = NamedTuple[]
    for tube_h in (0.0, 25.0, 100.0, 400.0)
        p = _stageC_parameters_2D_v19(
            source, ua, tube_h,
        )
        cases = _simulate_set_2D_v19(training_inputs, p)
        score = cooling_plant_score_2D_v19(cases)
        push!(plant_results, (
            tube_h=tube_h, cases=cases, score=score,
        ))
        push!(plant_rows, merge((
            tube_flange_contact_h_W_m2_K=tube_h,
        ), score))
        println("v19-C plant h=$tube_h score=", score)
        flush(stdout)
    end
    plant = first(sort(
        plant_results; by=x -> x.score.objective,
    ))
    _write_namedtuples_csv_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19,
        "stageC_tube_flange_ablation_2D_v19.csv",
    ), plant_rows)

    probe_rows = NamedTuple[]
    probe_results = NamedTuple[]
    capacities = (200.0, 600.0, 1200.0, 2400.0, 3000.0)
    stems = (0.0, 20.0, 60.0, 120.0, 200.0)
    for capacity in capacities, stem in stems
        score = probe_score_2D_v19(
            plant.cases, capacity, stem,
        )
        push!(probe_results, (
            capacity=capacity, stem=stem, score=score,
        ))
        push!(probe_rows, merge((
            capacity_areal_J_m2_K=capacity,
            stem_conductance_areal_W_m2_K=stem,
        ), score))
    end
    probe = first(sort(
        probe_results; by=x -> x.score.objective,
    ))
    _write_namedtuples_csv_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19,
        "stageC_probe_profile_2D_v19.csv",
    ), probe_rows)

    holdout_input = case_inputs_2D_v19(
        "C81"; cooling=true, max_points=61,
    )
    holdout_case = simulate_case_2D_v19(
        holdout_input,
        _stageC_parameters_2D_v19(
            source, ua, plant.tube_h;
            capacity=probe.capacity, stem=probe.stem,
        );
        reltol=5e-4, abstol=1e-4, dtmax=120.0,
    )
    holdout = probe_score_2D_v19(
        [holdout_case], probe.capacity, probe.stem,
    )
    plant_holdout = cooling_plant_score_2D_v19([holdout_case])
    capacity_interior = probe.capacity != first(capacities) &&
                        probe.capacity != last(capacities)
    stem_interior = probe.stem != first(stems) &&
                    probe.stem != last(stems)
    pass = probe.score.T3_rmse_K <= 30.0 &&
        holdout.T3_rmse_K <= 30.0 &&
        holdout.max_abs_final_bias_K <= 20.0 &&
        holdout.t90_mae_s <= 500.0 &&
        capacity_interior && stem_interior
    open(joinpath(
        OUTPUT_DIR_2D_v19, "stageC_selected_2D_v19.txt",
    ), "w") do io
        println(io, "rear_tube_flange_contact_h_W_m2_K=",
                plant.tube_h)
        println(io, "probe_capacity_areal_J_m2_K=",
                probe.capacity)
        println(io, "probe_stem_conductance_areal_W_m2_K=",
                probe.stem)
        println(io, "training_T3_rmse_K=",
                probe.score.T3_rmse_K)
        println(io, "holdout_T3_rmse_K=", holdout.T3_rmse_K)
        println(io, "holdout_final_bias_K=",
                holdout.final_bias_K)
        println(io, "holdout_t90_mae_s=", holdout.t90_mae_s)
        println(io, "holdout_side_rmse_K=",
                plant_holdout.side_rmse_K)
        println(io, "holdout_felt_rmse_K=",
                plant_holdout.felt_rmse_K)
        println(io, "stageC_pass=", pass)
    end
    println("v19-C plant=", plant.tube_h,
            " probe=", probe, " holdout=", holdout,
            " pass=", pass)
    return (
        plant=plant, probe=probe, holdout=holdout,
        plant_holdout=plant_holdout, pass=pass,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    calibrate_2D_v19_stageC()
end
