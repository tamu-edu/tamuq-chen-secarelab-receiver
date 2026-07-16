# ============================================================================
# 1D_v7.jl - 1D receiver model with measured T2 cavity boundary
# ============================================================================
# v7 keeps the reduced v5/v6 heat-transfer formulation, but replaces the free
# rear thermal mass and empirical side/rear losses with a geometry-based cavity
# resistance to the measured insulation temperature T2.
#
# Calibration stages:
#   1. Heating fits gas heat-transfer shape and three irradiance-level factors.
#   2. Cooling fits thermophysical / conduction scales and T3 lag.
#   3. Heating refits gas heat-transfer shape and irradiance-level factors.
#
# The T2 thermocouple is treated as a measured boundary condition, not as a
# fitted state. This intentionally removes the nearly adiabatic rear-mass degree
# of freedom that v6 used to store energy during cooling.
# ============================================================================

include("1D_v5.jl")

begin # v7 fixed constants
    const B_RE_FIXED_V7 = B_RE_FIXED_V5
    const C_PR_FIXED_V7 = C_PR_FIXED_V5
    const EPS_FRONT_FIXED_V7 = EPS_FRONT_FIXED_V5
    const ETA_ABS_FIXED_V7 = ETA_ABS_FIXED_V5
    const BETA_OPT_FIXED_V7 = BETA_OPT_FIXED_V5
    const FRONT_DEPOSITION_FIXED_V7 = FRONT_DEPOSITION_FIXED_V5
    const H_FRONT_FIXED_V7 = 10.0

    # Cavity geometry inherited from the older 0D_v3 resistance network.
    const ADAPTOR_RADIUS_V7 = 19.4e-3
    const ADAPTOR_LENGTH_V7 = 57.0e-3
    const ADAPTOR_TUBE_RADIUS_V7 = 6.5e-3
    const ADAPTOR_OVERLAP_LENGTH_V7 = 29.0e-3
    const ADAPTOR_FREE_LENGTH_V7 = 28.0e-3
    const INSULATION_OUTER_RADIUS_V7 = 75.0e-3
    const METAL_THICKNESS_V7 = 18.0e-3
    const METAL_OUTER_RADIUS_V7 = INSULATION_OUTER_RADIUS_V7 + METAL_THICKNESS_V7
    const METAL_LENGTH_V7 = L + ADAPTOR_FREE_LENGTH_V7
    const BACKPLATE_THICKNESS_V7 = METAL_THICKNESS_V7
    const INSULATION_MONOLITH_LENGTH_V7 = 108.0e-3
    const INSULATION_ADAPTOR_LENGTH_V7 = 57.0e-3
    const INSULATION_CONDUCTIVITY_V7 = 0.08

    # T2 is interpreted as a boundary inside the felt, about 40 mm radially
    # outside the equivalent receiver wall. The metal shell and outer felt are
    # downstream of this measured boundary and are therefore not solved as
    # additional hidden states in v7.
    const RECEIVER_EQ_RADIUS_V7 = sqrt(A_frt / pi)
    const T2_RADIAL_OFFSET_V7 = 40.0e-3
    const T2_BOUNDARY_RADIUS_V7 = min(
        RECEIVER_EQ_RADIUS_V7 + T2_RADIAL_OFFSET_V7,
        0.98 * INSULATION_OUTER_RADIUS_V7,
    )
    const RECEIVER_TO_T2_CONDUCTANCE_PER_LENGTH_V7 =
        2.0 * pi * INSULATION_CONDUCTIVITY_V7 /
        log(T2_BOUNDARY_RADIUS_V7 / RECEIVER_EQ_RADIUS_V7)
    const ADAPTOR_TO_T2_RESISTANCE_V7 =
        log(T2_BOUNDARY_RADIUS_V7 / ADAPTOR_RADIUS_V7) /
        (2.0 * pi * INSULATION_CONDUCTIVITY_V7 * ADAPTOR_LENGTH_V7)
    const ADAPTOR_CONTACT_RESISTANCE_V7 =
        1.0 / (100.0 * 4.0 * w_t * ADAPTOR_OVERLAP_LENGTH_V7)
    const RECEIVER_REAR_TO_T2_CONDUCTANCE_V7 =
        1.0 / (ADAPTOR_CONTACT_RESISTANCE_V7 + ADAPTOR_TO_T2_RESISTANCE_V7)
end

# Final v7 parameter vector p[1:10]
# Cooling / thermophysical:
#   p[1]  gamma_C       receiver solid heat-capacity multiplier
#   p[2]  k_scale       axial receiver conductivity multiplier
#   p[3]  k_ins_scale   felt/adaptor-to-T2 conductance multiplier
#   p[7]  tau_T3        outlet sensor/downstream lag [s]
# Heating:
#   p[4]  A_Nu          Nusselt multiplier exponent in 10^A_Nu
#   p[5]  h_floor       downstream exchange floor, 0<h_floor<=1
#   p[6]  L_h           axial exchange decay length [m]
#   p[8]  f_I_high      irradiance correction for high-Io group
#   p[9]  f_I_mid       irradiance correction for mid-Io group
#   p[10] f_I_low       irradiance correction for low-Io group
begin # v7 parameter values and bounds
    pnew_v7 = [
        1.0, 1.0, 1.0,
        -1.0, 0.25, 0.050,
        20.0,
        1.0, 1.0, 1.0,
    ]
    lb_full_v7 = [
        0.20, 0.20, 0.25,
        -4.00, 0.02, 0.005,
        0.50,
        0.60, 0.60, 0.60,
    ]
    ub_full_v7 = [
        6.00, 3.00, 4.00,
        2.00, 1.00, 0.500,
        500.0,
        1.60, 1.60, 1.60,
    ]

    fit_heat_transfer_indices_v7 = [4, 5, 6, 8, 9, 10]
    fit_cooling_thermal_indices_v7 = [1, 2, 3, 7]
end

irradiance_factor_v7(irradiance, p) =
    irradiance >= 3.80e5 ? p[8] :
    irradiance >= 2.80e5 ? p[9] : p[10]

function axial_exchange_shape_v7(z, p)
    floor_value = clamp(p[5], 0.0, 1.0)
    length_scale = max(p[6], 1e-6)
    return floor_value + (1.0 - floor_value) * exp(-max(0.0, z) / length_scale)
end

function gas_profile_v7!(Tg, Qgas, hcell, Ts, time, p, operating, z, dx)
    flow = max(0.0, operating.flow_lpm(time))
    Tin = operating.inlet_temperature(time)
    Tg[1] = Tin

    if flow <= 1e-12
        fill!(Tg, Tin)
        fill!(Qgas, 0.0)
        fill!(hcell, 0.0)
        return nothing
    end

    mdot = m_dot(flow, Tin)
    for i in eachindex(Ts)
        Tfilm = 0.5 * (Ts[i] + Tg[i])
        cp = cpf_f(Tfilm)
        mu = mu_f_f(Tfilm)
        kf = kf_f(Tfilm)
        reynolds = mdot * Dh / (A_flow * mu)
        prandtl = cp * mu / kf
        nusselt = 10.0^p[4] *
                  reynolds^B_RE_FIXED_V7 *
                  prandtl^(10.0^C_PR_FIXED_V7)
        hcell[i] = nusselt * kf / Dh * axial_exchange_shape_v7(z[i], p)
        UA = hcell[i] * P_exchange * dx
        effectiveness = -expm1(-UA / (mdot * cp))
        Tg[i + 1] = Tg[i] + effectiveness * (Ts[i] - Tg[i])
        Qgas[i] = mdot * cp * (Tg[i + 1] - Tg[i])
    end
    return nothing
end

function receiver_rhs_v7!(dTs, Ts, time, p, operating, z, dx, solar_weights,
                          t2_boundary, Tg, Qgas, hcell)
    nodes = length(z)
    ambient = operating.ambient_temperature(time)
    T2 = t2_boundary(time)

    gas_profile_v7!(Tg, Qgas, hcell, Ts, time, p, operating, z, dx)
    fill!(dTs, 0.0)

    for i in 1:(nodes - 1)
        ki = p[2] * ks_f(Ts[i])
        kj = p[2] * ks_f(Ts[i + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        Qcond = kface * A_solid * (Ts[i] - Ts[i + 1]) / dx
        dTs[i] -= Qcond
        dTs[i + 1] += Qcond
    end

    irradiance = max(0.0, operating.irradiance(time))
    Qsolar = ETA_ABS_FIXED_V7 * irradiance_factor_v7(irradiance, p) *
             irradiance * A_frt
    radial_conductance = p[3] * RECEIVER_TO_T2_CONDUCTANCE_PER_LENGTH_V7
    for i in 1:nodes
        Qside_T2 = radial_conductance * dx * (Ts[i] - T2)
        dTs[i] += Qsolar * solar_weights[i] - Qgas[i] - Qside_T2
    end

    Qfront = H_FRONT_FIXED_V7 * A_frt * (Ts[1] - ambient) +
             EPS_FRONT_FIXED_V7 * sigma_sb * A_frt * (Ts[1]^4 - ambient^4)
    Qrear_T2 = p[3] * RECEIVER_REAR_TO_T2_CONDUCTANCE_V7 * (Ts[end] - T2)

    dTs[1] -= Qfront
    dTs[end] -= Qrear_T2

    cell_volume = A_solid * dx
    for i in 1:nodes
        capacity = rho_s * p[1] * Cps_f(Ts[i]) * cell_volume
        dTs[i] /= capacity
    end
    return nothing
end

function receiver_ode_v7!(dTs, Ts, context, time)
    receiver_rhs_v7!(
        dTs, Ts, time, context.p, context.operating,
        context.z, context.dx, context.solar_weights, context.t2_boundary,
        context.Tg, context.Qgas, context.hcell,
    )
    return nothing
end

function simulate_v7(p, operating, save_times;
                     initial_temperature=Tamb, nodes=default_nodes,
                     t2_boundary=t -> Tamb,
                     solver=Rodas5P(autodiff=AutoFiniteDiff()),
                     reltol=1e-6, abstol=1e-7, dtmax=30.0)
    length(p) == 10 || throw(ArgumentError("1D_v7 expects 10 parameters"))
    all(isfinite, p) || throw(ArgumentError("all model parameters must be finite"))
    p[1] > 0.0 && p[2] > 0.0 && p[3] > 0.0 ||
        throw(ArgumentError("capacity, conductivity, and insulation scale must be positive"))
    nodes >= 3 || throw(ArgumentError("at least three axial nodes are required"))

    time = Float64.(save_times)
    isempty(time) && throw(ArgumentError("save_times cannot be empty"))
    (length(time) == 1 || all(diff(time) .> 0.0)) ||
        throw(ArgumentError("save_times must be strictly increasing"))

    dx = L / nodes
    z = collect(range(dx / 2, step=dx, length=nodes))
    zg = collect(range(0.0, L, length=nodes + 1))
    weights = solar_weights_v5(BETA_OPT_FIXED_V7, FRONT_DEPOSITION_FIXED_V7, nodes)
    Ts_initial = initial_profile_v3(initial_temperature, z)

    Tg_history = Matrix{Float64}(undef, nodes + 1, length(time))
    h_history = Matrix{Float64}(undef, nodes, length(time))
    t2_history = Vector{Float64}(undef, length(time))
    Tg = zeros(nodes + 1)
    Qgas = zeros(nodes)
    hcell = zeros(nodes)

    context = (
        p=p, operating=operating, z=z, dx=dx, solar_weights=weights,
        t2_boundary=t2_boundary, Tg=Tg, Qgas=Qgas, hcell=hcell,
    )
    problem = ODEProblem(receiver_ode_v7!, Ts_initial, (time[1], time[end]), context)
    solution = solve(
        problem, solver; saveat=time, save_everystep=false, dense=false,
        reltol=reltol, abstol=abstol, dtmax=dtmax,
    )
    successful_retcode(solution) ||
        error("1D_v7 ODE solve failed with retcode $(solution.retcode)")

    Ts_history = reduce(hcat, solution.u)
    for output_index in eachindex(time)
        t = time[output_index]
        gas_profile_v7!(Tg, Qgas, hcell, view(Ts_history, :, output_index),
                        t, p, operating, z, dx)
        Tg_history[:, output_index] .= Tg
        h_history[:, output_index] .= hcell
        t2_history[output_index] = t2_boundary(t)
    end

    return (
        time=time, z_solid=z, z_gas=zg,
        solid_temperature=Ts_history,
        boundary_temperature=t2_history,
        gas_temperature=Tg_history,
        heat_transfer_coefficient=h_history,
        ode_solution=solution,
    )
end

begin # v7 experimental interface
    paired_observation_v7(data, simulation_id, obs_a, obs_b) =
        0.5 .* (observation(data, simulation_id, obs_a) .+
                observation(data, simulation_id, obs_b))

    function measured_initial_profile_v7(data, simulation_id)
        z8 = sensor_positions[:T8]
        z9 = sensor_positions[:T9]
        z10 = sensor_positions[:T10]
        T8 = observation(data, simulation_id, "_T8")[1]
        T9 = paired_observation_v7(data, simulation_id, "_T9", "_T12")[1]
        T10 = paired_observation_v7(data, simulation_id, "_T10", "_T11")[1]
        return function (z)
            z <= z8 && return T8
            z >= z10 && return T10
            if z <= z9
                return T8 + (T9 - T8) * (z - z8) / (z9 - z8)
            end
            return T9 + (T10 - T9) * (z - z9) / (z10 - z9)
        end
    end

    function experimental_case_v7(simulation_id; is_cooling=false)
        data = is_cooling ? measurements_cooling : measurements
        conditions = is_cooling ? simulation_conditions_cooling : simulation_conditions
        time = observation_time(data, simulation_id)
        experiment = hcat(
            observation(data, simulation_id, "_T8"),
            paired_observation_v7(data, simulation_id, "_T9", "_T12"),
            paired_observation_v7(data, simulation_id, "_T10", "_T11"),
            observation(data, simulation_id, "_Tf"),
        )
        return conditions[simulation_id], time, experiment, data
    end

    function extract_outputs_v7(result, p, T3_initial)
        T8 = solid_at_v3(result, sensor_positions[:T8])
        T9_pair = solid_at_v3(result, sensor_positions[:T9])
        T10_pair = solid_at_v3(result, sensor_positions[:T10])
        Tf_true = vec(result.gas_temperature[end, :])
        Tf = apply_T3_lag_v3(Tf_true, T3_initial, result.time, p[7])
        h_mean = vec(mean(result.heat_transfer_coefficient, dims=1))
        return hcat(T8, T9_pair, T10_pair, Tf, h_mean)
    end

    function solve_case_v7(p, simulation_id; is_cooling=false,
                           nodes=default_nodes,
                           solver=Rodas5P(autodiff=AutoFiniteDiff()),
                           reltol=1e-6, abstol=1e-7, dtmax=30.0)
        conditions, time, experiment, data = experimental_case_v7(
            simulation_id; is_cooling=is_cooling
        )
        flow = observation(data, simulation_id, "_flow")
        Tin = observation(data, simulation_id, "_Tin")
        ambient = observation(data, simulation_id, "_Tamb")
        T2 = observation(data, simulation_id, "_T2")
        irradiance = fill(conditions[Io], length(time))
        operating = operating_condition_v3(
            irradiance=linear_history_v3(time, irradiance),
            flow_lpm=linear_history_v3(time, flow),
            inlet_temperature=linear_history_v3(time, Tin),
            ambient_temperature=linear_history_v3(time, ambient),
        )
        result = simulate_v7(
            p, operating, time;
            initial_temperature=measured_initial_profile_v7(data, simulation_id),
            t2_boundary=linear_history_v3(time, T2),
            nodes=nodes, solver=solver, reltol=reltol, abstol=abstol, dtmax=dtmax,
        )
        outputs = extract_outputs_v7(result, p, experiment[1, 4])
        return outputs, result, experiment
    end
end

begin # v7 calibration
    const TIMING_WEIGHT_V7 = 0.15
    const SLOPE_WEIGHT_V7 = 0.25

    function normalized_slope_mse_v7(model, experiment)
        length(model) < 2 && return 0.0
        scale = max(maximum(experiment) - minimum(experiment), 1.0)
        return mean(abs2, diff(model) .- diff(experiment)) / scale^2
    end

    function timing_penalty_v7(time, model, experiment)
        duration = max(time[end] - time[1], 1.0)
        return ((get_t90_v3(time, model) - get_t90_v3(time, experiment)) / duration)^2
    end

    function signal_loss_v7(time, model, experiment)
        return normalized_mse_v3(model, experiment) +
               SLOPE_WEIGHT_V7 * normalized_slope_mse_v7(model, experiment) +
               TIMING_WEIGHT_V7 * timing_penalty_v7(time, model, experiment)
    end

    function loss_cases_v7(p, keys; is_cooling=false, nodes=15)
        total = 0.0
        for simulation_id in keys
            model, result, experiment = solve_case_v7(
                p, simulation_id; is_cooling=is_cooling, nodes=nodes,
                reltol=2e-5, abstol=2e-6, dtmax=60.0,
            )
            all(isfinite, model) || return Inf
            total += mean(
                signal_loss_v7(result.time, model[:, j], experiment[:, j])
                for j in 1:4
            )
        end
        return total / length(keys)
    end

    function loss_heating_v7(p, keys=sim_key_heat; nodes=15)
        regularization = 0.003 * sum(log(p[i])^2 for i in 8:10)
        return loss_cases_v7(p, keys; is_cooling=false, nodes=nodes) + regularization
    end

    function loss_cooling_v7(p, keys=sim_key_cool; nodes=15)
        regularization = 0.005 * (log(p[1])^2 + log(p[2])^2 + log(p[3])^2) +
                         0.001 * log(p[7] / 20.0)^2
        return loss_cases_v7(p, keys; is_cooling=true, nodes=nodes) + regularization
    end
end

function optimize_parameter_subset_v7(objective, p_initial, indices, lower, upper;
                                      maximum_iterations=200,
                                      maximum_time_seconds=120.0,
                                      optimizer=NLopt.LN_NELDERMEAD(),
                                      label="subset-v7")
    base = Float64.(p_initial)
    selected = collect(indices)
    function subset_objective(theta)
        trial = copy(base)
        trial[selected] .= theta
        return objective(trial)
    end
    fit = optimize_with_nlopt_v3(
        subset_objective,
        base[selected],
        Float64.(lower[selected]),
        Float64.(upper[selected]);
        maximum_iterations=maximum_iterations,
        maximum_time_seconds=maximum_time_seconds,
        algorithm=optimizer,
        label=label,
    )
    parameters = copy(base)
    parameters[selected] .= fit.minimizer
    return merge(fit, (parameters=parameters, indices=selected))
end

function calibrate_v7(; cooling_keys=sim_key_cool,
                      heating_keys=sim_key_heat,
                      nodes=15,
                      heating_base_iterations=20,
                      cooling_iterations=20,
                      heating_refit_iterations=20,
                      cooling_time_seconds=120.0,
                      heating_base_time_seconds=120.0,
                      heating_refit_time_seconds=120.0,
                      optimizer=NLopt.LN_NELDERMEAD())
    p0 = copy(pnew_v7)

    println("\n=== 1D_v7 STAGE 1: heating heat-transfer + irradiance factors ===")
    heating_base_fit = optimize_parameter_subset_v7(
        p -> loss_heating_v7(p, heating_keys; nodes=nodes),
        p0, fit_heat_transfer_indices_v7, lb_full_v7, ub_full_v7;
        maximum_iterations=nlopt_evaluation_budget_v3(
            heating_base_iterations, length(fit_heat_transfer_indices_v7)
        ),
        maximum_time_seconds=heating_base_time_seconds,
        optimizer=optimizer,
        label="heating-base-v7",
    )
    p_after_heating_base = heating_base_fit.parameters

    println("\n=== 1D_v7 STAGE 2: cooling thermophysical / T2-boundary calibration ===")
    cooling_fit = optimize_parameter_subset_v7(
        p -> loss_cooling_v7(p, cooling_keys; nodes=nodes),
        p_after_heating_base, fit_cooling_thermal_indices_v7, lb_full_v7, ub_full_v7;
        maximum_iterations=nlopt_evaluation_budget_v3(
            cooling_iterations, length(fit_cooling_thermal_indices_v7)
        ),
        maximum_time_seconds=cooling_time_seconds,
        optimizer=optimizer,
        label="cooling-v7",
    )
    p_after_cooling = cooling_fit.parameters

    println("\n=== 1D_v7 STAGE 3: heating heat-transfer + irradiance refit ===")
    heating_refit = optimize_parameter_subset_v7(
        p -> loss_heating_v7(p, heating_keys; nodes=nodes),
        p_after_cooling, fit_heat_transfer_indices_v7, lb_full_v7, ub_full_v7;
        maximum_iterations=nlopt_evaluation_budget_v3(
            heating_refit_iterations, length(fit_heat_transfer_indices_v7)
        ),
        maximum_time_seconds=heating_refit_time_seconds,
        optimizer=optimizer,
        label="heating-refit-v7",
    )
    parameters = heating_refit.parameters

    global pnew_v7 = parameters
    global calibration_v7 = (
        heating_base=heating_base_fit,
        cooling=cooling_fit,
        heating_refit=heating_refit,
        fixed=(
            eta_abs=ETA_ABS_FIXED_V7,
            beta_opt=BETA_OPT_FIXED_V7,
            front_deposition_fraction=FRONT_DEPOSITION_FIXED_V7,
            h_front=H_FRONT_FIXED_V7,
            eps_front=EPS_FRONT_FIXED_V7,
            t2_boundary_radius=T2_BOUNDARY_RADIUS_V7,
            receiver_to_t2_conductance_per_length=RECEIVER_TO_T2_CONDUCTANCE_PER_LENGTH_V7,
            receiver_rear_to_t2_conductance=RECEIVER_REAR_TO_T2_CONDUCTANCE_V7,
        ),
        parameters=parameters,
    )
    println("\nFinal 1D_v7 parameters:")
    println(parameters)
    return calibration_v7
end

function quick_calibration_v7()
    return calibrate_v7(
        cooling_keys=["C69"], heating_keys=["E74"],
        nodes=11,
        heating_base_iterations=2,
        cooling_iterations=2,
        heating_refit_iterations=2,
    )
end

begin # v7 post-processing
    function transient_case_v7(simulation_id, p=pnew_v7; is_cooling=false,
                               nodes=default_nodes)
        model, result, experiment = solve_case_v7(
            p, simulation_id; is_cooling=is_cooling, nodes=nodes
        )
        return (
            time=result.time,
            T8_model=model[:, 1], T8_experiment=experiment[:, 1],
            T9_pair_model=model[:, 2], T9_pair_experiment=experiment[:, 2],
            T10_pair_model=model[:, 3], T10_pair_experiment=experiment[:, 3],
            T3_model=model[:, 4], T3_experiment=experiment[:, 4],
            T2_boundary=result.boundary_temperature,
            h_effective=model[:, 5],
            full_result=result,
        )
    end

    function plot_case_v7(simulation_id, params=pnew_v7; is_cooling=false,
                          nodes=default_nodes)
        @eval using StatsPlots
        data = transient_case_v7(simulation_id, params;
                                 is_cooling=is_cooling, nodes=nodes)
        sensors = (:T8, :T9_pair, :T10_pair, :T3)
        colors = (:blue, :red, :green, :purple)
        plot_object = plot(title="1D_v7 transient: $simulation_id",
                           xlabel="Time (s)", ylabel="Temperature (K)",
                           legend=:outerright)
        for (sensor, color) in zip(sensors, colors)
            plot!(plot_object, data.time, getproperty(data, Symbol(sensor, "_model"));
                  label="$(sensor) model", lw=2, color=color)
            scatter!(plot_object, data.time,
                     getproperty(data, Symbol(sensor, "_experiment"));
                     label="$(sensor) experiment", ms=2, markerstrokewidth=0,
                     color=color)
        end
        plot!(plot_object, data.time, data.T2_boundary;
              label="T2 boundary", color=:black, linestyle=:dash)
        return plot_object
    end

    function compute_metrics_v7(p=pnew_v7; heating_keys=sim_key_heat,
                                cooling_keys=sim_key_cool, nodes=default_nodes)
        metrics = NamedTuple[]
        for (keys, cooling) in ((heating_keys, false), (cooling_keys, true))
            for simulation_id in keys
                model, result, experiment = solve_case_v7(
                    p, simulation_id; is_cooling=cooling, nodes=nodes
                )
                for (sensor, column) in ((:T8, 1), (:T9_pair, 2),
                                         (:T10_pair, 3), (:T3, 4))
                    push!(metrics, (
                        simulation_id=simulation_id,
                        phase=cooling ? :cooling : :heating,
                        sensor=sensor,
                        rmse_K=sqrt(mean(abs2, model[:, column] .- experiment[:, column])),
                        steady_error_K=model[end, column] - experiment[end, column],
                        t90_error_s=get_t90_v3(result.time, model[:, column]) -
                                    get_t90_v3(result.time, experiment[:, column]),
                        shape_loss=normalized_slope_mse_v7(
                            model[:, column], experiment[:, column]
                        ),
                    ))
                end
            end
        end
        return metrics
    end
end

println("""
[1D_v7] Model loaded.
Run quick_calibration_v7() for a short check.
Run calibration_v7 = calibrate_v7() for the full staged calibration.
The calibrated vector is stored in pnew_v7.
""")
