include("1D_v25.jl")

const FITTED_PARAMS_v25_LOW_IRRADIANCE_BREAKDOWN = [
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

const LOW_IRRADIANCE_KEYS_v25_BREAKDOWN =
    [key for key in sim_key_heat if simulation_conditions[key][Io] == 256000.0]

const SENSOR_LABELS_v25_BREAKDOWN = (:T8, :T12, :T11, :T9, :T10, :T3, :T2)

function print_breakdown_row_v25(factor, p; nodes=11)
    sensor_losses = zeros(length(SENSOR_LABELS_v25_BREAKDOWN))
    ordering = 0.0
    steady_errors = zeros(length(SENSOR_LABELS_v25_BREAKDOWN))
    for simulation_id in LOW_IRRADIANCE_KEYS_v25_BREAKDOWN
        model, result, experiment = solve_case_v25(
            p, simulation_id; nodes=nodes, reltol=2e-5, abstol=2e-6, dtmax=60.0,
        )
        for j in eachindex(SENSOR_LABELS_v25_BREAKDOWN)
            sensor_losses[j] += signal_loss_v25(result.time, model[:, j], experiment[:, j])
            steady_errors[j] += model[end, j] - experiment[end, j]
        end
        ordering += ORDERING_WEIGHT_v25 * ordering_loss_v25(model, experiment)
    end
    sensor_losses ./= length(LOW_IRRADIANCE_KEYS_v25_BREAKDOWN)
    steady_errors ./= length(LOW_IRRADIANCE_KEYS_v25_BREAKDOWN)
    ordering /= length(LOW_IRRADIANCE_KEYS_v25_BREAKDOWN)
    println(join((
        factor,
        p[9],
        sensor_losses...,
        ordering,
        steady_errors...,
    ), ','))
end

println(join((
    "scale256_factor",
    "scale_256",
    ("loss_" * String(label) for label in SENSOR_LABELS_v25_BREAKDOWN)...,
    "ordering_loss",
    ("steady_err_" * String(label) for label in SENSOR_LABELS_v25_BREAKDOWN)...,
), ','))

for factor in (1.0, 1.1, 1.2, 1.35, 1.5)
    p = copy(FITTED_PARAMS_v25_LOW_IRRADIANCE_BREAKDOWN)
    p[9] *= factor
    print_breakdown_row_v25(factor, p)
end
