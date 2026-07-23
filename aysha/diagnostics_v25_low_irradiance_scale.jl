include("1D_v25.jl")

const FITTED_PARAMS_v25_LOW_IRRADIANCE = [
    4.917603371707668,
    0.5127592508804074,
    1.0 / 3.0,
    1.0,
    0.0,
    1.0,
    2.0034190016068183,
    1.9980406077869783,
    0.9523485964817118,
    16.044477788186573,
    214.29323236142696,
    7.591183315111345,
    184.67243519237965,
    0.5251021463644217,
    3.629367998027823,
    0.9993097154917425,
    0.10138177146664679,
    3.336402359922935,
    0.0,
    0.999767609832555,
    0.0,
]

const LOW_IRRADIANCE_KEYS_v25 = [key for key in sim_key_heat if simulation_conditions[key][Io] == 256000.0]

function heating_loss_subset_v25(p, keys; nodes=11)
    return loss_cases_v25(p, keys; is_cooling=false, nodes=nodes)
end

function steady_biases_v25(p, keys; nodes=15)
    errors = Dict(sensor => Float64[] for sensor in (:T8, :T12, :T11, :T9, :T10, :T3, :T2))
    for simulation_id in keys
        outputs, _, experiment = solve_case_v25(p, simulation_id; nodes=nodes)
        push!(errors[:T8], outputs[end, 1] - experiment[end, 1])
        push!(errors[:T12], outputs[end, 2] - experiment[end, 2])
        push!(errors[:T11], outputs[end, 3] - experiment[end, 3])
        push!(errors[:T9], outputs[end, 4] - experiment[end, 4])
        push!(errors[:T10], outputs[end, 5] - experiment[end, 5])
        push!(errors[:T3], outputs[end, 6] - experiment[end, 6])
        push!(errors[:T2], outputs[end, 7] - experiment[end, 7])
    end
    return errors
end

function mean_error_v25(errors, sensor)
    return sum(errors[sensor]) / length(errors[sensor])
end

function mean_abs_error_v25(errors, sensor)
    return sum(abs, errors[sensor]) / length(errors[sensor])
end

function print_row_v25(label, p)
    low_errors = steady_biases_v25(p, LOW_IRRADIANCE_KEYS_v25)
    all_errors = steady_biases_v25(p, sim_key_heat)
    println(join((
        label,
        p[9],
        loss_heating_v25(p; nodes=11) + loss_cooling_v25(p; nodes=11),
        loss_heating_v25(p; nodes=11),
        heating_loss_subset_v25(p, LOW_IRRADIANCE_KEYS_v25; nodes=11),
        mean_error_v25(low_errors, :T8),
        mean_error_v25(low_errors, :T12),
        mean_error_v25(low_errors, :T11),
        mean_error_v25(low_errors, :T9),
        mean_error_v25(low_errors, :T10),
        mean_error_v25(low_errors, :T3),
        mean_abs_error_v25(low_errors, :T8),
        mean_abs_error_v25(low_errors, :T12),
        mean_abs_error_v25(low_errors, :T11),
        mean_abs_error_v25(low_errors, :T9),
        mean_abs_error_v25(low_errors, :T10),
        mean_abs_error_v25(low_errors, :T3),
        mean_error_v25(all_errors, :T8),
        mean_error_v25(all_errors, :T12),
        mean_error_v25(all_errors, :T9),
        mean_error_v25(all_errors, :T3),
    ), ','))
end

println("label,scale_256,objective_total,objective_heating,objective_256_heating,mean_256_T8_err,mean_256_T12_err,mean_256_T11_err,mean_256_T9_err,mean_256_T10_err,mean_256_T3_err,mae_256_T8,mae_256_T12,mae_256_T11,mae_256_T9,mae_256_T10,mae_256_T3,mean_all_T8_err,mean_all_T12_err,mean_all_T9_err,mean_all_T3_err")

for factor in (0.70, 0.85, 1.00, 1.10, 1.20, 1.35, 1.50, 1.75, 2.00)
    p = copy(FITTED_PARAMS_v25_LOW_IRRADIANCE)
    p[9] = FITTED_PARAMS_v25_LOW_IRRADIANCE[9] * factor
    print_row_v25("scale256_factor_$(factor)", p)
end
