# ============================================================================
# Staged v11 calibration
#
# Stage 1: fit profile-sensitive beta_opt, S_h, and k_scale_z using only the
#          nine established heating training cases and temperature
#          differences/flow ordering.
# Stage 2: fit the three irradiance-group power scales using only absolute
#          temperature levels in those same training cases.
# Validation: six untouched heating cases, three cooling cases, and a
#             representative nominal-mesh comparison.
# ============================================================================

using LinearAlgebra
using Statistics
using Optim

include("run_2D_v11_inheritance_sensitivity.jl")

const TRAIN_HEATING_2D_v11 = (
    "E67", "E69", "E71",
    "E72", "E74", "E76",
    "E77", "E79", "E81",
)
const VALID_HEATING_2D_v11 = (
    "E68", "E70", "E73", "E75", "E78", "E80",
)
const VALID_COOLING_2D_v11 = ("C69", "C80", "C81")
const SENSOR_NAMES_STAGED_2D_v11 = (
    :T8, :T12, :T11, :T9, :T10, :T3, :T2,
)
const PROFILE_LOWER_2D_v11 = [10.0, 0.05, 0.20]
const PROFILE_UPPER_2D_v11 = [300.0, 2.00, 2.00]
const POWER_LOWER_2D_v11 = [0.30, 0.30, 0.30]
const POWER_UPPER_2D_v11 = [2.20, 2.20, 1.50]

function with_exact_optics_2D_v11(
    p;
    beta_opt=p.optics.extinction_coefficient,
    scales=(
        p.optics.scale_456,
        p.optics.scale_304,
        p.optics.scale_256,
    ),
)
    o0 = p.optics
    optics = OpticalParameters2D(
        absorbed_fraction=o0.absorbed_fraction,
        extinction_coefficient=Float64(beta_opt),
        beam_radius_sigma=o0.beam_radius_sigma,
        spillage_fraction=o0.spillage_fraction,
        front_deposition_fraction=o0.front_deposition_fraction,
        scale_456=Float64(scales[1]),
        scale_304=Float64(scales[2]),
        scale_256=Float64(scales[3]),
    )
    return ModelParameters2D(
        p.geometry, p.solid, p.heat_transfer,
        p.losses, optics, p.hydraulics,
    )
end

function staged_base_2D_v11(; nominal_mesh=false)
    fitted = unpack_parameters2D(
        fitted_v9_theta_2D_v11(), default_parameters2D(),
    )
    base = nominal_mesh ? fitted : sensitivity_mesh_2D_v11(fitted)
    graetz = with_heat_model_2D_v11(
        base;
        model=:graetz,
        nu_fd=3.61,
        temperature_exponent=0.0,
        scale=1.0,
    )
    return with_hydraulics_2D_v11(
        graetz;
        offset=-0.5428447496201336,
        resistance_scale=0.9701974617588522,
        old_hot_K=0.0,
        distribution=:equal_pressure,
        groove_diameter=13.0e-3,
        groove_K=184.15800344228472,
        max_iterations=16,
        tolerance=1.0e-5,
    )
end

function staged_parameters_2D_v11(
    base,
    profile_theta,
    power_theta,
)
    beta_opt, heat_scale, axial_scale = profile_theta
    thermal = with_heat_model_2D_v11(
        base;
        model=:graetz,
        nu_fd=3.61,
        temperature_exponent=0.0,
        scale=heat_scale,
    )
    thermal = with_solid_sensitivity_2D_v11(
        thermal; axial_scale=axial_scale,
    )
    return with_exact_optics_2D_v11(
        thermal;
        beta_opt=beta_opt,
        scales=Tuple(power_theta),
    )
end

function observed_final_vector_2D_v11(id; cooling=false)
    data = cooling ? measurements_cooling : measurements
    return [
        Float64(observation(data, id, "_$(sensor)")[end])
        for sensor in SENSOR_NAMES_STAGED_2D_v11
    ]
end

function compact_heating_prediction_2D_v11(id, p)
    case = operating_case_2D_v11(id, p)
    result = simulate2D(
        p,
        case.op,
        [case.times[1], case.times[end]];
        initial_temperature=case.ambient[1],
    )
    pred = sensor_predictions2D(result)
    values = [
        getproperty(pred, sensor)[end]
        for sensor in SENSOR_NAMES_STAGED_2D_v11
    ]
    return (
        simulation_id=id,
        flow_lpm=mean(case.flow),
        irradiance_kW_m2=case.nominal / 1000.0,
        model=values,
        observed=observed_final_vector_2D_v11(id),
        dp1_model_mbar=result.dp1_prediction_mbar[end],
        dp1_observed_mbar=observed_final_2D_v11(id, "_DP1"),
    )
end

function compact_training_predictions_2D_v11(p)
    return compact_predictions_2D_v11(TRAIN_HEATING_2D_v11, p)
end

function compact_predictions_2D_v11(ids, p)
    outputs = Vector{Any}(undef, length(ids))
    Threads.@threads for index in eachindex(ids)
        outputs[index] = compact_heating_prediction_2D_v11(
            ids[index], p,
        )
    end
    return outputs
end

function profile_residuals_2D_v11(outputs)
    residuals = Float64[]
    axial_model = Dict{String, Float64}()
    axial_observed = Dict{String, Float64}()
    for output in outputs
        m = output.model
        e = output.observed
        # Side-wall axial profile: 5 -> 58 -> 107 mm.
        push!(residuals, ((m[2] - m[1]) - (e[2] - e[1])) / 50.0)
        push!(residuals, ((m[3] - m[2]) - (e[3] - e[2])) / 40.0)
        # Core axial profile: 58 -> 107 mm.
        push!(residuals, ((m[5] - m[4]) - (e[5] - e[4])) / 40.0)
        # Rear-tube relation to the mean rear receiver temperature.
        push!(
            residuals,
            (
                (m[6] - 0.5 * (m[3] + m[5])) -
                (e[6] - 0.5 * (e[3] + e[5]))
            ) / 60.0,
        )
        axial_model[output.simulation_id] = m[2] - m[1]
        axial_observed[output.simulation_id] = e[2] - e[1]
    end
    # Explicitly constrain within-irradiance flow ordering independently of
    # each group's absolute profile level.
    available_ids = Set(keys(axial_model))
    for (_, ids) in HEATING_GROUPS_2D_v11
        train_ids = [
            id for id in ids
            if id in TRAIN_HEATING_2D_v11 && id in available_ids
        ]
        length(train_ids) >= 2 || continue
        model_mean = mean(axial_model[id] for id in train_ids)
        observed_mean = mean(axial_observed[id] for id in train_ids)
        for id in train_ids
            push!(
                residuals,
                (
                    (axial_model[id] - model_mean) -
                    (axial_observed[id] - observed_mean)
                ) / 25.0,
            )
        end
    end
    return residuals
end

function level_residuals_2D_v11(outputs)
    residuals = Float64[]
    for output in outputs
        m = output.model
        e = output.observed
        # Radially averaged receiver temperatures prevent the known side/core
        # observation mismatch from dominating the power calibration.
        model_levels = (
            m[1],
            0.5 * (m[2] + m[4]),
            0.5 * (m[3] + m[5]),
            m[6],
            m[7],
        )
        observed_levels = (
            e[1],
            0.5 * (e[2] + e[4]),
            0.5 * (e[3] + e[5]),
            e[6],
            e[7],
        )
        append!(
            residuals,
            (model_levels .- observed_levels) ./ 100.0,
        )
    end
    return residuals
end

function normalized_parameter_penalty_2D_v11(
    theta,
    reference,
    scales;
    weight,
)
    return weight * mean(abs2, (theta .- reference) ./ scales)
end

function radical_inverse_2D_v11(index, base)
    result = 0.0
    factor = 1.0 / base
    value = index
    while value > 0
        result += factor * (value % base)
        value = div(value, base)
        factor /= base
    end
    return result
end

function quadratic_features_2D_v11(x)
    return [
        1.0,
        x[1], x[2], x[3],
        x[1]^2, x[2]^2, x[3]^2,
        x[1] * x[2],
        x[1] * x[3],
        x[2] * x[3],
    ]
end

function fit_profile_stage_2D_v11(base, power_theta)
    trace = NamedTuple[]
    evaluations = Ref(0)
    reference = [110.0, 1.0, 1.0]
    ranges = PROFILE_UPPER_2D_v11 .- PROFILE_LOWER_2D_v11
    evaluated_theta = Vector{Vector{Float64}}()
    evaluated_loss = Float64[]
    function evaluate(theta)
        candidate = clamp.(
            Float64.(theta),
            PROFILE_LOWER_2D_v11,
            PROFILE_UPPER_2D_v11,
        )
        evaluations[] += 1
        p = staged_parameters_2D_v11(base, candidate, power_theta)
        outputs = compact_training_predictions_2D_v11(p)
        residuals = profile_residuals_2D_v11(outputs)
        loss = mean(abs2, residuals) +
               normalized_parameter_penalty_2D_v11(
                   candidate, reference, ranges; weight=0.01,
               )
        push!(evaluated_theta, candidate)
        push!(evaluated_loss, loss)
        push!(trace, (
            stage="profile",
            evaluation=evaluations[],
            objective=loss,
            beta_opt=candidate[1],
            heat_transfer_scale=candidate[2],
            axial_conductivity_scale=candidate[3],
        ))
        println(
            "profile eval $(evaluations[]): loss=$(round(loss,digits=6)) " *
            "theta=$(round.(candidate,digits=5))",
        )
        flush(stdout)
        return loss
    end

    design = Vector{Vector{Float64}}()
    push!(design, [110.0, 1.0, 1.0])
    push!(design, [110.0, 0.50, 1.0])
    push!(design, [160.0, 1.0, 1.4])
    push!(design, [80.0, 1.5, 0.6])
    for index in 1:18
        normalized = [
            radical_inverse_2D_v11(index, 2),
            radical_inverse_2D_v11(index, 3),
            radical_inverse_2D_v11(index, 5),
        ]
        push!(
            design,
            PROFILE_LOWER_2D_v11 .+ normalized .* ranges,
        )
    end
    foreach(evaluate, design)

    normalized_design = [
        (theta .- PROFILE_LOWER_2D_v11) ./ ranges
        for theta in evaluated_theta
    ]
    feature_matrix = reduce(
        vcat,
        permutedims(quadratic_features_2D_v11(x))
        for x in normalized_design
    )
    coefficients = feature_matrix \ evaluated_loss
    grid_values = range(0.0, 1.0; length=26)
    best_x = [0.0, 0.0, 0.0]
    best_surrogate = Inf
    for x1 in grid_values, x2 in grid_values, x3 in grid_values
        x = [x1, x2, x3]
        value = dot(quadratic_features_2D_v11(x), coefficients)
        if value < best_surrogate
            best_surrogate = value
            best_x .= x
        end
    end
    surrogate_theta =
        PROFILE_LOWER_2D_v11 .+ best_x .* ranges
    evaluate(surrogate_theta)
    for index in eachindex(surrogate_theta), direction in (-1.0, 1.0)
        local_theta = copy(surrogate_theta)
        local_theta[index] += direction * 0.08 * ranges[index]
        evaluate(local_theta)
    end
    best_index = argmin(evaluated_loss)
    return (
        theta=evaluated_theta[best_index],
        objective=evaluated_loss[best_index],
        trace=trace,
        evaluations=evaluations[],
        converged=true,
        surrogate_coefficients=coefficients,
        surrogate_theta=surrogate_theta,
        surrogate_prediction=best_surrogate,
    )
end

function fit_power_stage_2D_v11(base, profile_theta, initial_power)
    trace = NamedTuple[]
    evaluations = Ref(0)
    ranges = POWER_UPPER_2D_v11 .- POWER_LOWER_2D_v11
    power = collect(initial_power)
    group_losses = Float64[]
    for (group_index, (_, all_ids)) in enumerate(HEATING_GROUPS_2D_v11)
        ids = [id for id in all_ids if id in TRAIN_HEATING_2D_v11]
        candidates = unique(sort(vcat(
            collect(range(
                POWER_LOWER_2D_v11[group_index],
                POWER_UPPER_2D_v11[group_index];
                length=7,
            )),
            initial_power[group_index],
        )))
        tested_scales = Float64[]
        tested_losses = Float64[]
        function evaluate_scale(scale)
            trial_power = copy(power)
            trial_power[group_index] = clamp(
                scale,
                POWER_LOWER_2D_v11[group_index],
                POWER_UPPER_2D_v11[group_index],
            )
            evaluations[] += 1
            p = staged_parameters_2D_v11(
                base, profile_theta, trial_power,
            )
            outputs = compact_predictions_2D_v11(ids, p)
            level = level_residuals_2D_v11(outputs)
            profile = profile_residuals_2D_v11(outputs)
            loss = mean(abs2, level) +
                   0.20 * mean(abs2, profile) +
                   0.01 * (
                       (
                           trial_power[group_index] -
                           initial_power[group_index]
                       ) / ranges[group_index]
                   )^2
            push!(tested_scales, trial_power[group_index])
            push!(tested_losses, loss)
            push!(trace, (
                stage="power",
                group=group_index,
                evaluation=evaluations[],
                objective=loss,
                tested_scale=trial_power[group_index],
            ))
            println(
                "power group $group_index eval $(evaluations[]): " *
                "loss=$(round(loss,digits=6)) " *
                "scale=$(round(trial_power[group_index],digits=5))",
            )
            flush(stdout)
            return loss
        end
        foreach(evaluate_scale, candidates)
        design = hcat(
            ones(length(tested_scales)),
            tested_scales,
            tested_scales .^ 2,
        )
        coefficients = design \ tested_losses
        if coefficients[3] > 0.0
            vertex = clamp(
                -coefficients[2] / (2.0 * coefficients[3]),
                POWER_LOWER_2D_v11[group_index],
                POWER_UPPER_2D_v11[group_index],
            )
            evaluate_scale(vertex)
        end
        best_index = argmin(tested_losses)
        power[group_index] = tested_scales[best_index]
        push!(group_losses, tested_losses[best_index])
    end
    return (
        theta=power,
        objective=mean(group_losses),
        trace=trace,
        evaluations=evaluations[],
        converged=true,
    )
end

function staged_response_vector_2D_v11(
    base,
    theta,
)
    profile_theta = theta[1:3]
    power_theta = theta[4:6]
    outputs = compact_training_predictions_2D_v11(
        staged_parameters_2D_v11(base, profile_theta, power_theta),
    )
    return vcat(
        profile_residuals_2D_v11(outputs),
        level_residuals_2D_v11(outputs),
    )
end

function identifiability_staged_2D_v11(base, theta)
    lower = vcat(PROFILE_LOWER_2D_v11, POWER_LOWER_2D_v11)
    upper = vcat(PROFILE_UPPER_2D_v11, POWER_UPPER_2D_v11)
    ranges = upper .- lower
    response0 = staged_response_vector_2D_v11(base, theta)
    jacobian = zeros(length(response0), length(theta))
    for index in eachindex(theta)
        step = 0.01 * ranges[index]
        low = copy(theta)
        high = copy(theta)
        low[index] = max(lower[index], theta[index] - step)
        high[index] = min(upper[index], theta[index] + step)
        response_low = staged_response_vector_2D_v11(base, low)
        response_high = staged_response_vector_2D_v11(base, high)
        delta = high[index] - low[index]
        jacobian[:, index] .= (
            response_high .- response_low
        ) ./ delta .* ranges[index]
    end
    singular_values = svdvals(jacobian)
    cutoff = 1.0e-3
    rank = count(
        value -> value > singular_values[1] * cutoff,
        singular_values,
    )
    norms = [norm(view(jacobian, :, i)) for i in axes(jacobian, 2)]
    correlation = zeros(length(theta), length(theta))
    for i in axes(correlation, 1), j in axes(correlation, 2)
        correlation[i, j] = dot(
            view(jacobian, :, i),
            view(jacobian, :, j),
        ) / max(norms[i] * norms[j], eps(Float64))
    end
    return (
        jacobian=jacobian,
        singular_values=singular_values,
        rank=rank,
        cutoff=cutoff,
        condition_number=singular_values[1] /
            max(singular_values[end], eps(Float64)),
        column_norms=norms,
        correlation=correlation,
        maximum_correlation=maximum(
            abs(correlation[i, j])
            for i in axes(correlation, 1), j in axes(correlation, 2)
            if i != j
        ),
    )
end

function cooling_t0_dict_2D_v11(id)
    return Dict(
        sensor => Float64(
            observation(measurements_cooling, id, "_$(sensor)")[1],
        )
        for sensor in SENSOR_NAMES_STAGED_2D_v11
    )
end

function full_case_staged_2D_v11(id, p; cooling=false)
    data = cooling ? measurements_cooling : measurements
    conditions = cooling ? simulation_conditions_cooling :
                 simulation_conditions
    times = Float64.(observation_time(data, id))
    flow = Float64.(observation(data, id, "_flow"))
    inlet = Float64.(observation(data, id, "_Tin"))
    ambient = Float64.(observation(data, id, "_Tamb"))
    nominal = Float64(conditions[id][Io])
    scale = nominal >= 400000.0 ? p.optics.scale_456 :
            (nominal >= 280000.0 ? p.optics.scale_304 :
             p.optics.scale_256)
    irradiance = cooling ? zeros(length(times)) :
                 fill(scale * nominal, length(times))
    op = OperatingCondition2D(
        irradiance=linear_history(times, irradiance),
        flow_lpm=linear_history(times, flow),
        inlet_temperature=linear_history(times, inlet),
        ambient_temperature=linear_history(times, ambient),
    )
    if cooling
        grid = Receiver2D_v11.build_grid2D(p)
        initial = build_initial_state_2D(
            grid, p, cooling_t0_dict_2D_v11(id), ambient[1],
        )
    else
        initial = ambient[1]
    end
    result = simulate2D(
        p, op, times; initial_temperature=initial,
    )
    pred = sensor_predictions2D(result)
    model = hcat(
        (getproperty(pred, sensor) for
         sensor in SENSOR_NAMES_STAGED_2D_v11)...,
    )
    observed = hcat(
        (
            Float64.(observation(data, id, "_$(sensor)"))
            for sensor in SENSOR_NAMES_STAGED_2D_v11
        )...,
    )
    dp1_observed = Float64.(observation(data, id, "_DP1"))
    return (
        simulation_id=id,
        cooling=cooling,
        times=times,
        flow=flow,
        model=model,
        observed=observed,
        dp1_model=result.dp1_prediction_mbar,
        dp1_observed=dp1_observed,
        result=result,
    )
end

function validation_metrics_staged_2D_v11(cases)
    sensor_rows = NamedTuple[]
    case_rows = NamedTuple[]
    transient_rows = NamedTuple[]
    for case in cases
        phase = case.cooling ? "cooling_validation" :
                (case.simulation_id in TRAIN_HEATING_2D_v11 ?
                 "heating_training" : "heating_validation")
        for (index, sensor) in enumerate(SENSOR_NAMES_STAGED_2D_v11)
            model = case.model[:, index]
            observed = case.observed[:, index]
            push!(sensor_rows, (
                phase=phase,
                simulation_id=case.simulation_id,
                sensor=String(sensor),
                rmse_K=sqrt(mean(abs2, model .- observed)),
                steady_error_K=model[end] - observed[end],
                t90_error_s=get_t90_2D(case.times, model) -
                    get_t90_2D(case.times, observed),
                shape_loss=normalized_slope_mse_2D(model, observed),
            ))
        end
        m = case.model[end, :]
        e = case.observed[end, :]
        push!(case_rows, (
            phase=phase,
            simulation_id=case.simulation_id,
            mean_flow_lpm=mean(case.flow),
            model_T8_K=m[1],
            observed_T8_K=e[1],
            model_T12_K=m[2],
            observed_T12_K=e[2],
            model_T11_K=m[3],
            observed_T11_K=e[3],
            model_T9_K=m[4],
            observed_T9_K=e[4],
            model_T10_K=m[5],
            observed_T10_K=e[5],
            model_T3_K=m[6],
            observed_T3_K=e[6],
            model_T2_K=m[7],
            observed_T2_K=e[7],
            model_T12_minus_T8_K=m[2] - m[1],
            observed_T12_minus_T8_K=e[2] - e[1],
            model_T12_minus_T9_K=m[2] - m[4],
            observed_T12_minus_T9_K=e[2] - e[4],
            model_T11_minus_T10_K=m[3] - m[5],
            observed_T11_minus_T10_K=e[3] - e[5],
            dp1_model_mbar=case.dp1_model[end],
            dp1_observed_mbar=case.dp1_observed[end],
        ))
        if case.simulation_id in ("E67", "E68", "C81")
            for index in eachindex(case.times)
                push!(transient_rows, (
                    phase=phase,
                    simulation_id=case.simulation_id,
                    time_s=case.times[index],
                    model_T8_K=case.model[index, 1],
                    observed_T8_K=case.observed[index, 1],
                    model_T12_K=case.model[index, 2],
                    observed_T12_K=case.observed[index, 2],
                    model_T11_K=case.model[index, 3],
                    observed_T11_K=case.observed[index, 3],
                    model_T9_K=case.model[index, 4],
                    observed_T9_K=case.observed[index, 4],
                    model_T10_K=case.model[index, 5],
                    observed_T10_K=case.observed[index, 5],
                    model_T3_K=case.model[index, 6],
                    observed_T3_K=case.observed[index, 6],
                    model_T2_K=case.model[index, 7],
                    observed_T2_K=case.observed[index, 7],
                ))
            end
        end
    end
    return (
        sensor_rows=sensor_rows,
        case_rows=case_rows,
        transient_rows=transient_rows,
    )
end

function aggregate_validation_staged_2D_v11(sensor_rows, case_rows)
    rows = NamedTuple[]
    for phase in unique(row.phase for row in sensor_rows)
        sensors = filter(row -> row.phase == phase, sensor_rows)
        cases = filter(row -> row.phase == phase, case_rows)
        push!(rows, (
            phase=phase,
            mean_sensor_rmse_K=mean(row.rmse_K for row in sensors),
            median_sensor_rmse_K=median(row.rmse_K for row in sensors),
            steady_mae_K=mean(abs(row.steady_error_K) for row in sensors),
            t90_mae_s=mean(abs(row.t90_error_s) for row in sensors),
            axial_rmse_K=sqrt(mean(
                (
                    row.model_T12_minus_T8_K -
                    row.observed_T12_minus_T8_K
                )^2 for row in cases
            )),
            mid_radial_rmse_K=sqrt(mean(
                (
                    row.model_T12_minus_T9_K -
                    row.observed_T12_minus_T9_K
                )^2 for row in cases
            )),
            deep_radial_rmse_K=sqrt(mean(
                (
                    row.model_T11_minus_T10_K -
                    row.observed_T11_minus_T10_K
                )^2 for row in cases
            )),
            dp1_rmse_mbar=sqrt(mean(
                (row.dp1_model_mbar - row.dp1_observed_mbar)^2
                for row in cases
            )),
            dp1_bias_mbar=mean(
                row.dp1_model_mbar - row.dp1_observed_mbar
                for row in cases
            ),
        ))
    end
    return rows
end

function mesh_confirmation_staged_2D_v11(
    profile_theta,
    power_theta,
)
    rows = NamedTuple[]
    for (mesh, base) in (
        ("sensitivity", staged_base_2D_v11(nominal_mesh=false)),
        ("nominal", staged_base_2D_v11(nominal_mesh=true)),
    )
        p = staged_parameters_2D_v11(base, profile_theta, power_theta)
        for id in ("E67", "E72", "E77")
            case = full_case_staged_2D_v11(id, p; cooling=false)
            m = case.model[end, :]
            push!(rows, (
                mesh=mesh,
                simulation_id=id,
                T8_K=m[1],
                T12_K=m[2],
                T11_K=m[3],
                T9_K=m[4],
                T10_K=m[5],
                T3_K=m[6],
                T2_K=m[7],
                dp1_mbar=case.dp1_model[end],
            ))
        end
    end
    return rows
end

function write_staged_parameters_2D_v11(
    path,
    profile_fit,
    power_fit,
)
    names = (
        "beta_opt_1_m",
        "heat_transfer_scale",
        "axial_conductivity_scale",
        "scale_456",
        "scale_304",
        "scale_256",
    )
    values = vcat(profile_fit.theta, power_fit.theta)
    rows = [
        (index=i, parameter=names[i], value=values[i])
        for i in eachindex(values)
    ]
    return write_namedtuple_csv_2D_v11(path, rows)
end

function main_staged_calibration_2D_v11()
    BLAS.set_num_threads(1)
    mkpath(OUTPUT_DIR_2D_v11)
    base = staged_base_2D_v11()
    initial_power = [
        base.optics.scale_456,
        base.optics.scale_304,
        base.optics.scale_256,
    ]

    println(
        "Staged v11 calibration with $(Threads.nthreads()) Julia threads",
    )
    println("Warm-up training solve...")
    flush(stdout)
    warm = compact_training_predictions_2D_v11(
        staged_parameters_2D_v11(
            base, [110.0, 1.0, 1.0], initial_power,
        ),
    )
    println(
        "Warm-up complete; profile objective=",
        mean(abs2, profile_residuals_2D_v11(warm)),
    )
    flush(stdout)

    println("Stage 1: profile-sensitive transport fit")
    profile_fit = fit_profile_stage_2D_v11(base, initial_power)
    println("Stage 2: irradiance-group power fit")
    power_fit = fit_power_stage_2D_v11(
        base, profile_fit.theta, initial_power,
    )

    profile_trace_path = joinpath(
        OUTPUT_DIR_2D_v11, "staged_profile_trace_2D_v11.csv",
    )
    power_trace_path = joinpath(
        OUTPUT_DIR_2D_v11, "staged_power_trace_2D_v11.csv",
    )
    write_namedtuple_csv_2D_v11(profile_trace_path, profile_fit.trace)
    write_namedtuple_csv_2D_v11(power_trace_path, power_fit.trace)
    write_staged_parameters_2D_v11(
        joinpath(
            OUTPUT_DIR_2D_v11,
            "staged_parameters_2D_v11.csv",
        ),
        profile_fit,
        power_fit,
    )

    fitted_theta = vcat(profile_fit.theta, power_fit.theta)
    println("Identifiability audit at staged fit")
    flush(stdout)
    ident = identifiability_staged_2D_v11(base, fitted_theta)
    parameter_names = (
        "beta_opt", "S_h", "k_scale_z",
        "scale_456", "scale_304", "scale_256",
    )
    ident_rows = [
        (
            parameter=parameter_names[i],
            column_norm=ident.column_norms[i],
        )
        for i in eachindex(parameter_names)
    ]
    write_namedtuple_csv_2D_v11(
        joinpath(
            OUTPUT_DIR_2D_v11,
            "staged_identifiability_columns_2D_v11.csv",
        ),
        ident_rows,
    )
    correlation_rows = NamedTuple[]
    for i in eachindex(parameter_names), j in eachindex(parameter_names)
        push!(correlation_rows, (
            parameter_i=parameter_names[i],
            parameter_j=parameter_names[j],
            correlation=ident.correlation[i, j],
        ))
    end
    write_namedtuple_csv_2D_v11(
        joinpath(
            OUTPUT_DIR_2D_v11,
            "staged_identifiability_correlation_2D_v11.csv",
        ),
        correlation_rows,
    )

    p_fitted = staged_parameters_2D_v11(
        base, profile_fit.theta, power_fit.theta,
    )
    println("Running full training and held-out validation cases")
    flush(stdout)
    all_cases = Any[]
    for id in vcat(
        collect(TRAIN_HEATING_2D_v11),
        collect(VALID_HEATING_2D_v11),
    )
        println("  heating $id")
        flush(stdout)
        push!(
            all_cases,
            full_case_staged_2D_v11(id, p_fitted; cooling=false),
        )
    end
    for id in VALID_COOLING_2D_v11
        println("  cooling $id")
        flush(stdout)
        push!(
            all_cases,
            full_case_staged_2D_v11(id, p_fitted; cooling=true),
        )
    end
    metrics = validation_metrics_staged_2D_v11(all_cases)
    aggregate = aggregate_validation_staged_2D_v11(
        metrics.sensor_rows, metrics.case_rows,
    )
    write_namedtuple_csv_2D_v11(
        joinpath(
            OUTPUT_DIR_2D_v11,
            "staged_sensor_metrics_2D_v11.csv",
        ),
        metrics.sensor_rows,
    )
    write_namedtuple_csv_2D_v11(
        joinpath(
            OUTPUT_DIR_2D_v11,
            "staged_case_metrics_2D_v11.csv",
        ),
        metrics.case_rows,
    )
    write_namedtuple_csv_2D_v11(
        joinpath(
            OUTPUT_DIR_2D_v11,
            "staged_transients_2D_v11.csv",
        ),
        metrics.transient_rows,
    )
    write_namedtuple_csv_2D_v11(
        joinpath(
            OUTPUT_DIR_2D_v11,
            "staged_validation_summary_2D_v11.csv",
        ),
        aggregate,
    )

    println("Running nominal-mesh confirmation")
    flush(stdout)
    mesh_rows = mesh_confirmation_staged_2D_v11(
        profile_fit.theta, power_fit.theta,
    )
    write_namedtuple_csv_2D_v11(
        joinpath(
            OUTPUT_DIR_2D_v11,
            "staged_mesh_cases_2D_v11.csv",
        ),
        mesh_rows,
    )

    summary_path = joinpath(
        OUTPUT_DIR_2D_v11,
        "staged_calibration_summary_2D_v11.txt",
    )
    open(summary_path, "w") do io
        println(io, "training_cases=$(collect(TRAIN_HEATING_2D_v11))")
        println(io, "validation_heating=$(collect(VALID_HEATING_2D_v11))")
        println(io, "validation_cooling=$(collect(VALID_COOLING_2D_v11))")
        println(io, "profile_theta=$(profile_fit.theta)")
        println(io, "profile_objective=$(profile_fit.objective)")
        println(io, "profile_evaluations=$(profile_fit.evaluations)")
        println(io, "profile_converged=$(profile_fit.converged)")
        println(io, "power_theta=$(power_fit.theta)")
        println(io, "power_objective=$(power_fit.objective)")
        println(io, "power_evaluations=$(power_fit.evaluations)")
        println(io, "power_converged=$(power_fit.converged)")
        println(io, "identifiability_singular_values=$(ident.singular_values)")
        println(io, "identifiability_rank=$(ident.rank)/6")
        println(io, "identifiability_condition=$(ident.condition_number)")
        println(io, "identifiability_max_correlation=$(ident.maximum_correlation)")
        foreach(row -> println(io, row), aggregate)
    end

    println("Staged calibration complete")
    println("profile_theta=$(profile_fit.theta)")
    println("power_theta=$(power_fit.theta)")
    println(
        "identifiability rank=$(ident.rank)/6, " *
        "condition=$(ident.condition_number), " *
        "max correlation=$(ident.maximum_correlation)",
    )
    println("validation summaries:")
    foreach(println, aggregate)
    return (
        profile_fit=profile_fit,
        power_fit=power_fit,
        identifiability=ident,
        validation=aggregate,
        mesh=mesh_rows,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_staged_calibration_2D_v11()
end
