# ============================================================================
# run_2D_v9_full.jl - Leakage-free thermal/hydraulic calibration and validation
# ============================================================================

using LinearAlgebra
using Statistics

include("run_2D_v9.jl")

const FULL_OUTPUT_DIR_2D_v9 = joinpath(@__DIR__, "summaries", "2D_v9")
const TRAIN_HEATING_2D_v9 = (
    "E67", "E69", "E71",
    "E72", "E74", "E76",
    "E77", "E79", "E81",
)
const VALIDATION_HEATING_2D_v9 = (
    "E68", "E70", "E73", "E75", "E78", "E80",
)
const VALIDATION_COOLING_2D_v9 = ("C69", "C80", "C81")
const SENSOR_NAMES_FULL_2D_v9 = (
    :T8, :T12_perim, :T11_perim, :T9_core, :T10_core, :T3, :T2,
)
const PARAMETER_NAMES_FULL_2D_v9 = (
    "A_Nu", "B_Re", "k_scale_r", "k_scale_z", "sigma_beam_mm",
    "spillage_fraction", "scale_456", "scale_304", "scale_256",
    "beta_rad", "receiver_felt_contact_resistance_m2K_W", "beta_opt",
)

function with_mesh_full_2D_v9(
    p; nr_rec, nr_felt, nr_case, nz, nz_rear,
)
    g0 = p.geometry
    geometry = Geometry2D(
        length=g0.length,
        receiver_width=g0.receiver_width,
        receiver_radius=g0.receiver_radius,
        insulation_outer_radius=g0.insulation_outer_radius,
        casing_outer_radius=g0.casing_outer_radius,
        channel_width=g0.channel_width,
        wall_thickness=g0.wall_thickness,
        channel_count=g0.channel_count,
        nodes_r_rec=nr_rec,
        nodes_r_felt=nr_felt,
        nodes_r_case=nr_case,
        nodes_z=nz,
        rear_tube_length=g0.rear_tube_length,
        rear_tube_nodes=nz_rear,
        rear_tube_inner_radius=g0.rear_tube_inner_radius,
        rear_tube_outer_radius=g0.rear_tube_outer_radius,
        t3_distance_from_receiver=g0.t3_distance_from_receiver,
    )
    return ModelParameters2D(
        geometry, p.solid, p.heat_transfer, p.losses, p.optics, p.hydraulics,
    )
end

calibration_mesh_full_2D_v9(p) = with_mesh_full_2D_v9(
    p; nr_rec=10, nr_felt=5, nr_case=2, nz=30, nz_rear=20,
)

function run_full_case_2D_v9(
    simulation_id, p=default_parameters2D(); is_cooling=false,
)
    data = is_cooling ? measurements_cooling : measurements
    hydraulic = run_hydraulic_case_2D_v9(
        simulation_id, p; cooling=is_cooling,
    )
    result = hydraulic.result
    predictions = sensor_predictions2D(result)
    model = hcat(
        predictions.T8, predictions.T12, predictions.T11,
        predictions.T9, predictions.T10, predictions.T3, predictions.T2,
    )
    experiment = hcat(
        Float64.(observation(data, simulation_id, "_T8")),
        Float64.(observation(data, simulation_id, "_T12")),
        Float64.(observation(data, simulation_id, "_T11")),
        Float64.(observation(data, simulation_id, "_T9")),
        Float64.(observation(data, simulation_id, "_T10")),
        Float64.(observation(data, simulation_id, "_T3")),
        Float64.(observation(data, simulation_id, "_T2")),
    )
    return (
        result=result,
        model=model,
        experiment=experiment,
        times=hydraulic.times,
        observed_dp1=hydraulic.observed_dp1,
    )
end

function run_calibration_case_full_2D_v9(
    simulation_id, p=default_parameters2D(); is_cooling=false,
)
    return run_full_case_2D_v9(
        simulation_id, calibration_mesh_full_2D_v9(p);
        is_cooling=is_cooling,
    )
end

function pilot_seed_full_2D_v9(p_base)
    # The v8 apparent-Nu pilot fitted B_Re=1.11359. Rescale A_Nu to the fixed
    # independently reported B_Re=1.44 at representative Re≈75.
    theta = [
        5.0e-4,
        1.44,
        0.05368421283813515,
        0.9086130005080822,
        14.0,
        0.10,
        0.8580655609964761,
        0.776185848837943,
        0.51672161518935,
        54.081387226214,
        0.001,
        96.10955829209641,
    ]
    return unpack_parameters2D(theta, p_base)
end

function with_minor_loss_full_2D_v9(p, coefficient)
    h0 = p.hydraulics
    hydraulics = HydraulicParameters2D(
        standard_pressure=h0.standard_pressure,
        standard_temperature=h0.standard_temperature,
        atmospheric_pressure=h0.atmospheric_pressure,
        mass_flow_scale=h0.mass_flow_scale,
        dp1_zero_offset_mbar=h0.dp1_zero_offset_mbar,
        hydraulic_resistance_scale=h0.hydraulic_resistance_scale,
        minor_loss_coefficient=Float64(coefficient),
    )
    return ModelParameters2D(
        p.geometry, p.solid, p.heat_transfer, p.losses, p.optics, hydraulics,
    )
end

function hot_excess_dp_basis_full_2D_v9(result, p)
    grid = Receiver2D_v9.build_grid2D(p)
    R_rec = grid.r_faces[grid.nr_rec + 1]
    psi = [
        1.0 - p.heat_transfer.c_radial_flow *
        (grid.r_centers[i] / R_rec)^2 for i in 1:grid.nr_rec
    ]
    psi_sum = sum(
        psi[i] * grid.area_flow[i] for i in 1:grid.nr_rec
    )
    fractions = [
        psi[i] * grid.area_flow[i] / psi_sum for i in 1:grid.nr_rec
    ]
    rho_ref = air_density(
        p.hydraulics.standard_temperature,
        p.hydraulics.atmospheric_pressure,
    )
    basis = zeros(length(result.time))
    for k in eachindex(result.time)
        for i in 1:grid.nr_rec
            mdot_ring = result.mass_flow_kg_s[k] * fractions[i]
            u_ref = mdot_ring / (rho_ref * grid.area_flow[i])
            q_excess = 0.5 * (
                result.gas_density[i, end, k] *
                result.gas_velocity[i, end, k]^2 -
                rho_ref * u_ref^2
            )
            basis[k] += fractions[i] * q_excess / 100.0
        end
    end
    return basis
end

function fit_hot_excess_loss_full_2D_v9(training_cases, p)
    numerator = 0.0
    denominator = 0.0
    rows = NamedTuple[]
    for id in training_cases
        case = run_full_case_2D_v9(id, p; is_cooling=false)
        basis = hot_excess_dp_basis_full_2D_v9(case.result, p)
        residual = case.observed_dp1 .- case.result.dp1_prediction_mbar
        weight = 1.0 / length(basis)
        numerator += weight * dot(basis, residual)
        denominator += weight * dot(basis, basis)
        push!(rows, (simulation_id=id, case=case, basis=basis))
    end
    raw_coefficient = numerator / max(denominator, eps(Float64))
    return (
        coefficient=max(0.0, raw_coefficient),
        raw_coefficient=raw_coefficient,
        rows=rows,
    )
end

function parameter_bound_diagnostics_full_2D_v9(p)
    theta = pack_parameters2D(p)
    rows = NamedTuple[]
    for index in FIT_INDICES_2D
        span = UB_2D[index] - LB_2D[index]
        bound_distance = min(
            (theta[index] - LB_2D[index]) / span,
            (UB_2D[index] - theta[index]) / span,
        )
        push!(rows, (
            index=index,
            name=PARAMETER_NAMES_FULL_2D_v9[index],
            value=theta[index],
            lower=LB_2D[index],
            upper=UB_2D[index],
            normalized_bound_distance=bound_distance,
            near_bound=bound_distance < 0.01,
        ))
    end
    return rows
end

function evaluate_full_2D_v9(p, fitted_hot_loss)
    thermal_rows = NamedTuple[]
    dp_rows = NamedTuple[]
    ordering_rows = NamedTuple[]
    transient_rows = NamedTuple[]
    velocity_min = Inf
    velocity_max = -Inf

    groups = (
        (TRAIN_HEATING_2D_v9, :heating_train, false),
        (VALIDATION_HEATING_2D_v9, :heating_validation, false),
        (VALIDATION_COOLING_2D_v9, :cooling_validation, true),
    )

    for (ids, phase, cooling) in groups
        for id in ids
            case = run_full_case_2D_v9(id, p; is_cooling=cooling)
            basis = hot_excess_dp_basis_full_2D_v9(case.result, p)
            dp_base = case.result.dp1_prediction_mbar
            dp_augmented = dp_base .+ fitted_hot_loss .* basis
            dp_error_base = dp_base .- case.observed_dp1
            dp_error_augmented = dp_augmented .- case.observed_dp1

            velocity_min = min(velocity_min, minimum(case.result.gas_velocity))
            velocity_max = max(velocity_max, maximum(case.result.gas_velocity))

            for (j, sensor) in enumerate(SENSOR_NAMES_FULL_2D_v9)
                error = case.model[:, j] .- case.experiment[:, j]
                push!(thermal_rows, (
                    simulation_id=id,
                    phase=phase,
                    sensor=sensor,
                    rmse_K=sqrt(mean(abs2, error)),
                    bias_K=mean(error),
                    final_error_K=error[end],
                    t90_error_s=get_t90_2D(case.times, case.model[:, j]) -
                                  get_t90_2D(case.times, case.experiment[:, j]),
                    shape_loss=normalized_slope_mse_2D(
                        case.model[:, j], case.experiment[:, j],
                    ),
                ))
            end

            push!(dp_rows, (
                simulation_id=id,
                phase=phase,
                base_rmse_mbar=sqrt(mean(abs2, dp_error_base)),
                augmented_rmse_mbar=sqrt(mean(abs2, dp_error_augmented)),
                base_bias_mbar=mean(dp_error_base),
                augmented_bias_mbar=mean(dp_error_augmented),
                base_final_error_mbar=dp_error_base[end],
                augmented_final_error_mbar=dp_error_augmented[end],
            ))

            if !cooling
                mid_model = case.model[end, 2] - case.model[end, 4]
                mid_experiment = case.experiment[end, 2] -
                                 case.experiment[end, 4]
                deep_model = case.model[end, 3] - case.model[end, 5]
                deep_experiment = case.experiment[end, 3] -
                                  case.experiment[end, 5]
                push!(ordering_rows, (
                    simulation_id=id,
                    phase=phase,
                    mid_model_K=mid_model,
                    mid_experiment_K=mid_experiment,
                    mid_sign_correct=mid_model * mid_experiment > 0.0,
                    deep_model_K=deep_model,
                    deep_experiment_K=deep_experiment,
                    deep_sign_correct=deep_model * deep_experiment > 0.0,
                ))
            end

            for k in eachindex(case.times)
                push!(transient_rows, (
                    simulation_id=id,
                    phase=phase,
                    time_s=case.times[k],
                    flow_lpm=standard_mass_flow2D(1.0, p.hydraulics) > 0 ?
                        case.result.mass_flow_kg_s[k] /
                        standard_mass_flow2D(1.0, p.hydraulics) : 0.0,
                    dp1_observed_mbar=case.observed_dp1[k],
                    dp1_base_mbar=dp_base[k],
                    dp1_augmented_mbar=dp_augmented[k],
                    T8_model_K=case.model[k, 1],
                    T8_experiment_K=case.experiment[k, 1],
                    T12_model_K=case.model[k, 2],
                    T12_experiment_K=case.experiment[k, 2],
                    T11_model_K=case.model[k, 3],
                    T11_experiment_K=case.experiment[k, 3],
                    T9_model_K=case.model[k, 4],
                    T9_experiment_K=case.experiment[k, 4],
                    T10_model_K=case.model[k, 5],
                    T10_experiment_K=case.experiment[k, 5],
                    T3_model_K=case.model[k, 6],
                    T3_experiment_K=case.experiment[k, 6],
                    T2_model_K=case.model[k, 7],
                    T2_experiment_K=case.experiment[k, 7],
                ))
            end
        end
    end
    return (
        thermal_rows=thermal_rows,
        dp_rows=dp_rows,
        ordering_rows=ordering_rows,
        transient_rows=transient_rows,
        velocity_min=velocity_min,
        velocity_max=velocity_max,
    )
end

function write_named_tuples_full_2D_v9(path, rows)
    isempty(rows) && error("Cannot write empty rows to $path")
    names = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(names, ','))
        for row in rows
            println(io, join((getproperty(row, name) for name in names), ','))
        end
    end
    return path
end

function phase_summary_full_2D_v9(thermal_rows, phase)
    rows = filter(row -> row.phase == phase, thermal_rows)
    return (
        mean_sensor_rmse_K=mean(row.rmse_K for row in rows),
        median_sensor_rmse_K=median(row.rmse_K for row in rows),
        steady_mae_K=mean(abs(row.final_error_K) for row in rows),
        t90_mae_s=mean(abs(row.t90_error_s) for row in rows),
    )
end

function dp_phase_summary_full_2D_v9(dp_rows, phase)
    rows = filter(row -> row.phase == phase, dp_rows)
    return (
        mean_base_rmse_mbar=mean(row.base_rmse_mbar for row in rows),
        mean_augmented_rmse_mbar=mean(
            row.augmented_rmse_mbar for row in rows
        ),
        mean_base_bias_mbar=mean(row.base_bias_mbar for row in rows),
        mean_augmented_bias_mbar=mean(
            row.augmented_bias_mbar for row in rows
        ),
    )
end

function main_full_2D_v9()
    mkpath(FULL_OUTPUT_DIR_2D_v9)
    calibration_dp1 = cold_t0_dp1_calibration_2D_v9()
    p_hydraulic = with_t0_hydraulics_2D_v9(
        default_parameters2D(), calibration_dp1,
    )
    p_seed = pilot_seed_full_2D_v9(p_hydraulic)
    max_iters = parse(Int, get(ENV, "V9_MAX_ITERS", "120"))
    max_time = parse(Float64, get(ENV, "V9_MAX_TIME_S", "3600"))

    seed_objective = Receiver2D_v9.loss_function_2D(
        pack_parameters2D(p_seed),
        (TRAIN_HEATING_2D_v9, run_calibration_case_full_2D_v9),
    )
    println("v9_seed_objective=$seed_objective")
    println("v9_training_cases=$(collect(TRAIN_HEATING_2D_v9))")
    println("v9_validation_heating=$(collect(VALIDATION_HEATING_2D_v9))")
    flush(stdout)

    calibration = calibrate2D(
        run_calibration_case_full_2D_v9;
        heating_cases=TRAIN_HEATING_2D_v9,
        max_iters=max_iters,
        max_time=max_time,
        p_init=p_seed,
    )
    p_thermal = calibration.parameters
    hot_fit = fit_hot_excess_loss_full_2D_v9(
        TRAIN_HEATING_2D_v9, p_thermal,
    )
    p_fitted = with_minor_loss_full_2D_v9(
        p_thermal, hot_fit.coefficient,
    )
    evaluation = evaluate_full_2D_v9(p_thermal, hot_fit.coefficient)
    bound_rows = parameter_bound_diagnostics_full_2D_v9(p_fitted)

    write_named_tuples_full_2D_v9(
        joinpath(FULL_OUTPUT_DIR_2D_v9, "thermal_metrics_2D_v9.csv"),
        evaluation.thermal_rows,
    )
    write_named_tuples_full_2D_v9(
        joinpath(FULL_OUTPUT_DIR_2D_v9, "dp1_transient_metrics_2D_v9.csv"),
        evaluation.dp_rows,
    )
    write_named_tuples_full_2D_v9(
        joinpath(FULL_OUTPUT_DIR_2D_v9, "steady_ordering_2D_v9.csv"),
        evaluation.ordering_rows,
    )
    write_named_tuples_full_2D_v9(
        joinpath(FULL_OUTPUT_DIR_2D_v9, "transient_predictions_2D_v9.csv"),
        evaluation.transient_rows,
    )
    write_named_tuples_full_2D_v9(
        joinpath(FULL_OUTPUT_DIR_2D_v9, "parameter_bounds_2D_v9.csv"),
        bound_rows,
    )

    theta = pack_parameters2D(p_fitted)
    parameter_rows = [
        (
            index=i,
            name=PARAMETER_NAMES_FULL_2D_v9[i],
            value=theta[i],
            fitted=i in FIT_INDICES_2D,
        ) for i in eachindex(theta)
    ]
    write_named_tuples_full_2D_v9(
        joinpath(FULL_OUTPUT_DIR_2D_v9, "parameters_fitted_2D_v9.csv"),
        parameter_rows,
    )

    summaries = Dict(
        phase => phase_summary_full_2D_v9(
            evaluation.thermal_rows, phase,
        ) for phase in (
            :heating_train, :heating_validation, :cooling_validation,
        )
    )
    dp_summaries = Dict(
        phase => dp_phase_summary_full_2D_v9(
            evaluation.dp_rows, phase,
        ) for phase in (
            :heating_train, :heating_validation, :cooling_validation,
        )
    )
    mid_correct = count(row -> row.mid_sign_correct, evaluation.ordering_rows)
    deep_correct = count(
        row -> row.deep_sign_correct, evaluation.ordering_rows
    )
    near_bounds = [row.name for row in bound_rows if row.near_bound]

    summary_path = joinpath(
        FULL_OUTPUT_DIR_2D_v9, "full_test_summary_2D_v9.txt",
    )
    open(summary_path, "w") do io
        println(io, "training_cases=$(collect(TRAIN_HEATING_2D_v9))")
        println(
            io,
            "validation_heating_cases=$(collect(VALIDATION_HEATING_2D_v9))",
        )
        println(
            io,
            "validation_cooling_cases=$(collect(VALIDATION_COOLING_2D_v9))",
        )
        println(io, "seed_objective=$seed_objective")
        println(io, "fitted_objective=$(calibration.objective)")
        println(io, "return_code=$(calibration.retcode)")
        println(io, "evaluations=$(calibration.evaluations)")
        println(io, "parameters=$theta")
        println(
            io,
            "hot_excess_loss_coefficient=$(hot_fit.coefficient)",
        )
        println(
            io,
            "hot_excess_loss_raw_coefficient=$(hot_fit.raw_coefficient)",
        )
        for phase in (
            :heating_train, :heating_validation, :cooling_validation,
        )
            println(io, "$(phase)_thermal=$(summaries[phase])")
            println(io, "$(phase)_dp1=$(dp_summaries[phase])")
        end
        println(
            io,
            "mid_ordering_correct=$mid_correct/$(length(evaluation.ordering_rows))",
        )
        println(
            io,
            "deep_ordering_correct=$deep_correct/$(length(evaluation.ordering_rows))",
        )
        println(io, "near_bound_parameters=$near_bounds")
        println(
            io,
            "local_velocity_range_m_s=$(evaluation.velocity_min),$(evaluation.velocity_max)",
        )
    end

    println("v9_fitted_objective=$(calibration.objective)")
    println("v9_return_code=$(calibration.retcode)")
    println("v9_evaluations=$(calibration.evaluations)")
    println("v9_hot_excess_loss_coefficient=$(hot_fit.coefficient)")
    println("v9_thermal_summaries=$summaries")
    println("v9_dp1_summaries=$dp_summaries")
    println(
        "v9_spatial_ordering=mid $mid_correct/$(length(evaluation.ordering_rows)), ",
        "deep $deep_correct/$(length(evaluation.ordering_rows))",
    )
    println("v9_near_bound_parameters=$near_bounds")
    println("v9_summary=$summary_path")
    return (
        calibration=calibration,
        parameters=p_fitted,
        evaluation=evaluation,
        hot_fit=hot_fit,
        summary_path=summary_path,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_full_2D_v9()
end
