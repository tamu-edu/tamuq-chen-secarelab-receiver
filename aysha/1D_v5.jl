# ============================================================================
# 1D_v5.jl - Reduced-parameter staged 1D receiver calibration model
# ============================================================================
# v5 keeps the conservative finite-volume core from v3 but reduces the fitted
# parameter set. Cooling runs identify thermophysical and gas-exchange behavior;
# heating runs then identify optical/front-source behavior with cooling
# parameters fixed.
#
# New modeling choices:
#   - Gas HTC uses Nu = 10^A * Re^B_fixed * Pr^(10^C_fixed).
#   - Axial exchange shape is simple and bounded:
#       s(z) = h_floor + (1 - h_floor) exp(-z / L_h)
#   - Solar deposition uses the v1 optical parameterization:
#       Beer-Lambert cell-integrated weights plus front deposition fraction.
#   - Irradiance/optical parameters are fixed during calibration so incoming
#     energy can be adjusted manually.
#   - T9/T12 and T10/T11 are averaged for fitting/comparison because the
#     model has no radial degree of freedom.
# ============================================================================

include("1D_v3.jl")

begin # v5 fixed modeling constants
    const B_RE_FIXED_V5 = 1.0
    const C_PR_FIXED_V5 = 0.5
    const EPS_FRONT_FIXED_V5 = 0.95
    const ETA_ABS_FIXED_V5 = 0.80
    const BETA_OPT_FIXED_V5 = 50.0
    const FRONT_DEPOSITION_FIXED_V5 = 1.0
end

# Final v5 parameter vector p[1:12]
# Cooling / thermophysical parameters:
#   p[1]  gamma_C     effective solid heat-capacity multiplier
#   p[2]  k_scale     effective axial-conductivity multiplier
#   p[3]  U_side      side-loss conductance per receiver length [W/m/K]
#   p[4]  U_rear      rear/adaptor conductance [W/K]
#   p[5]  A_Nu        Nusselt multiplier exponent in 10^A_Nu
#   p[6]  h_floor     downstream exchange floor, 0<h_floor<=1
#   p[7]  L_h         axial exchange decay length [m]
#   p[8]  tau_T3      outlet sensor/downstream lag [s]
# Fixed optical / manually adjusted irradiance parameters:
#   p[9]  eta_abs     absorbed solar fraction
#   p[10] beta_opt    exponential optical attenuation [1/m]
#   p[11] front_dep   front-deposition fraction added to first cell
# Heating-fitted parameter:
#   p[12] h_front     front convection coefficient [W/m2/K]
begin # v5 parameter values and bounds
    p_cool_init_v5 = [1.0, 1.0, 0.35, 0.10, -1.0, 0.20, 0.050, 20.0]
    lb_cool_v5 = [0.10, 0.10, 0.005, 0.001, -4.0, 0.02, 0.005, 0.50]
    ub_cool_v5 = [5.00, 5.00, 5.000, 5.000, 2.0, 1.00, 0.500, 500.0]

    p_heat_init_v5 = [ETA_ABS_FIXED_V5, BETA_OPT_FIXED_V5, FRONT_DEPOSITION_FIXED_V5, 10.0]
    lb_heat_v5 = [ETA_ABS_FIXED_V5, BETA_OPT_FIXED_V5, FRONT_DEPOSITION_FIXED_V5, 0.10]
    ub_heat_v5 = [ETA_ABS_FIXED_V5, BETA_OPT_FIXED_V5, FRONT_DEPOSITION_FIXED_V5, 100.0]

    assemble_parameters_v5(p_cool, p_heat) = vcat(Float64.(p_cool), Float64.(p_heat))
    pnew_v5 = assemble_parameters_v5(p_cool_init_v5, p_heat_init_v5)
    lb_full_v5 = assemble_parameters_v5(lb_cool_v5, lb_heat_v5)
    ub_full_v5 = assemble_parameters_v5(ub_cool_v5, ub_heat_v5)

    fit_nonirradiance_indices_v5 = [1, 2, 3, 4, 5, 6, 7, 8, 12]
    fit_cooling_indices_v5 = collect(1:8)
end

function axial_exchange_shape_v5(z, p)
    floor_value = clamp(p[6], 0.0, 1.0)
    length_scale = max(p[7], 1e-6)
    return floor_value + (1.0 - floor_value) * exp(-max(0.0, z) / length_scale)
end

function solar_weights_v5(beta, front_deposition_fraction, nodes)
    dx = L / nodes
    beta_value = max(0.0, Float64(beta))
    weights = zeros(nodes)
    if !isfinite(beta_value) || beta_value > 1e8
        weights[1] = 1.0
    elseif beta_value <= 1e-12
        fill!(weights, 1.0 / nodes)
    else
        for i in 1:nodes
            left = (i - 1) * dx
            right = i * dx
            weights[i] = exp(-beta_value * left) - exp(-beta_value * right)
        end
        weights ./= sum(weights)
    end
    front = clamp(Float64(front_deposition_fraction), 0.0, 1.0)
    weights .*= 1.0 - front
    weights[1] += front
    return Float64.(weights)
end

function gas_profile_v5!(Tg, Qgas, hcell, Ts, time, p, operating, z, dx)
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
        nusselt = 10.0^p[5] *
                  reynolds^B_RE_FIXED_V5 *
                  prandtl^(10.0^C_PR_FIXED_V5)
        hcell[i] = nusselt * kf / Dh * axial_exchange_shape_v5(z[i], p)
        UA = hcell[i] * P_exchange * dx
        effectiveness = -expm1(-UA / (mdot * cp))
        Tg[i + 1] = Tg[i] + effectiveness * (Ts[i] - Tg[i])
        Qgas[i] = mdot * cp * (Tg[i + 1] - Tg[i])
    end
    return nothing
end

function receiver_rhs_v5!(dTs, Ts, time, p, operating, z, dx, solar_weights,
                          Tg, Qgas, hcell)
    nodes = length(Ts)
    ambient = operating.ambient_temperature(time)
    gas_profile_v5!(Tg, Qgas, hcell, Ts, time, p, operating, z, dx)
    fill!(dTs, 0.0)

    for i in 1:(nodes - 1)
        ki = p[2] * ks_f(Ts[i])
        kj = p[2] * ks_f(Ts[i + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        Qcond = kface * A_solid * (Ts[i] - Ts[i + 1]) / dx
        dTs[i] -= Qcond
        dTs[i + 1] += Qcond
    end

    Qsolar = p[9] * max(0.0, operating.irradiance(time)) * A_frt
    for i in 1:nodes
        Qside = p[3] * dx * (Ts[i] - ambient)
        dTs[i] += Qsolar * solar_weights[i] - Qgas[i] - Qside
    end

    Qfront = p[12] * A_frt * (Ts[1] - ambient) +
             EPS_FRONT_FIXED_V5 * sigma_sb * A_frt * (Ts[1]^4 - ambient^4)
    Qrear = p[4] * (Ts[end] - ambient)
    dTs[1] -= Qfront
    dTs[end] -= Qrear

    cell_volume = A_solid * dx
    for i in 1:nodes
        capacity = rho_s * p[1] * Cps_f(Ts[i]) * cell_volume
        dTs[i] /= capacity
    end
    return nothing
end

function receiver_ode_v5!(dTs, Ts, context, time)
    receiver_rhs_v5!(
        dTs, Ts, time, context.p, context.operating,
        context.z, context.dx, context.solar_weights,
        context.Tg, context.Qgas, context.hcell,
    )
    return nothing
end

function simulate_v5(p, operating, save_times;
                     initial_temperature=Tamb, nodes=default_nodes,
                     solver=Rodas5P(autodiff=AutoFiniteDiff()),
                     reltol=1e-6, abstol=1e-7, dtmax=30.0)
    length(p) == 12 || throw(ArgumentError("1D_v5 expects 12 parameters"))
    all(isfinite, p) || throw(ArgumentError("all model parameters must be finite"))
    p[1] > 0.0 && p[2] > 0.0 ||
        throw(ArgumentError("capacity and conductivity parameters must be positive"))
    nodes >= 3 || throw(ArgumentError("at least three axial nodes are required"))

    time = Float64.(save_times)
    isempty(time) && throw(ArgumentError("save_times cannot be empty"))
    (length(time) == 1 || all(diff(time) .> 0.0)) ||
        throw(ArgumentError("save_times must be strictly increasing"))

    dx = L / nodes
    z = collect(range(dx / 2, step=dx, length=nodes))
    zg = collect(range(0.0, L, length=nodes + 1))
    weights = solar_weights_v5(p[10], p[11], nodes)
    Ts_initial = initial_profile_v3(initial_temperature, z)
    Tg_history = Matrix{Float64}(undef, nodes + 1, length(time))
    h_history = Matrix{Float64}(undef, nodes, length(time))
    Tg = zeros(nodes + 1)
    Qgas = zeros(nodes)
    hcell = zeros(nodes)

    context = (
        p=p, operating=operating, z=z, dx=dx, solar_weights=weights,
        Tg=Tg, Qgas=Qgas, hcell=hcell,
    )
    problem = ODEProblem(receiver_ode_v5!, Ts_initial, (time[1], time[end]), context)
    solution = solve(
        problem, solver; saveat=time, save_everystep=false, dense=false,
        reltol=reltol, abstol=abstol, dtmax=dtmax,
    )
    successful_retcode(solution) ||
        error("1D_v5 ODE solve failed with retcode $(solution.retcode)")

    Ts_history = reduce(hcat, solution.u)
    for output_index in eachindex(time)
        gas_profile_v5!(Tg, Qgas, hcell, view(Ts_history, :, output_index),
                        time[output_index], p, operating, z, dx)
        Tg_history[:, output_index] .= Tg
        h_history[:, output_index] .= hcell
    end

    return (
        time=time, z_solid=z, z_gas=zg,
        solid_temperature=Ts_history,
        gas_temperature=Tg_history,
        heat_transfer_coefficient=h_history,
        ode_solution=solution,
    )
end

begin # v5 experimental interface
    paired_observation_v5(data, simulation_id, obs_a, obs_b) =
        0.5 .* (observation(data, simulation_id, obs_a) .+
                observation(data, simulation_id, obs_b))

    function measured_initial_profile_v5(data, simulation_id)
        z8 = sensor_positions[:T8]
        z9 = sensor_positions[:T9]
        z10 = sensor_positions[:T10]
        T8 = observation(data, simulation_id, "_T8")[1]
        T9 = paired_observation_v5(data, simulation_id, "_T9", "_T12")[1]
        T10 = paired_observation_v5(data, simulation_id, "_T10", "_T11")[1]
        return function (z)
            z <= z8 && return T8
            z >= z10 && return T10
            if z <= z9
                return T8 + (T9 - T8) * (z - z8) / (z9 - z8)
            end
            return T9 + (T10 - T9) * (z - z9) / (z10 - z9)
        end
    end

    function experimental_case_v5(simulation_id; is_cooling=false)
        data = is_cooling ? measurements_cooling : measurements
        conditions = is_cooling ? simulation_conditions_cooling : simulation_conditions
        time = observation_time(data, simulation_id)
        experiment = hcat(
            observation(data, simulation_id, "_T8"),
            paired_observation_v5(data, simulation_id, "_T9", "_T12"),
            paired_observation_v5(data, simulation_id, "_T10", "_T11"),
            observation(data, simulation_id, "_Tf"),
        )
        return conditions[simulation_id], time, experiment, data
    end

    function extract_outputs_v5(result, p, T3_initial)
        T8 = solid_at_v3(result, sensor_positions[:T8])
        T9_pair = solid_at_v3(result, sensor_positions[:T9])
        T10_pair = solid_at_v3(result, sensor_positions[:T10])
        Tf_true = vec(result.gas_temperature[end, :])
        Tf = apply_T3_lag_v3(Tf_true, T3_initial, result.time, p[8])
        h_mean = vec(mean(result.heat_transfer_coefficient, dims=1))
        return hcat(T8, T9_pair, T10_pair, Tf, h_mean)
    end

    function solve_case_v5(p, simulation_id; is_cooling=false,
                           nodes=default_nodes,
                           solver=Rodas5P(autodiff=AutoFiniteDiff()),
                           reltol=1e-6, abstol=1e-7, dtmax=30.0)
        conditions, time, experiment, data = experimental_case_v5(
            simulation_id; is_cooling=is_cooling
        )
        flow = observation(data, simulation_id, "_flow")
        Tin = observation(data, simulation_id, "_Tin")
        ambient = observation(data, simulation_id, "_Tamb")
        irradiance = fill(conditions[Io], length(time))
        operating = operating_condition_v3(
            irradiance=linear_history_v3(time, irradiance),
            flow_lpm=linear_history_v3(time, flow),
            inlet_temperature=linear_history_v3(time, Tin),
            ambient_temperature=linear_history_v3(time, ambient),
        )
        result = simulate_v5(
            p, operating, time;
            initial_temperature=measured_initial_profile_v5(data, simulation_id),
            nodes=nodes, solver=solver, reltol=reltol, abstol=abstol, dtmax=dtmax,
        )
        outputs = extract_outputs_v5(result, p, experiment[1, 4])
        return outputs, result, experiment
    end
end

begin # v5 calibration
    function loss_cases_v5(p, keys; is_cooling=false, nodes=15)
        total = 0.0
        for simulation_id in keys
            model, _, experiment = solve_case_v5(
                p, simulation_id; is_cooling=is_cooling, nodes=nodes,
                reltol=2e-5, abstol=2e-6, dtmax=60.0,
            )
            all(isfinite, model) || return Inf
            total += mean(normalized_mse_v3(model[:, j], experiment[:, j]) for j in 1:4)
        end
        return total / length(keys)
    end

    function loss_cooling_v5(p, keys=sim_key_cool; nodes=15)
        regularization = 0.005 * (log(p[1])^2 + log(p[2])^2) +
                         0.001 * (log(p[6] / 0.20)^2 +
                                  log(p[7] / 0.050)^2)
        return loss_cases_v5(p, keys; is_cooling=true, nodes=nodes) + regularization
    end

    function loss_heating_v5(p, keys=sim_key_heat; nodes=15)
        return loss_cases_v5(p, keys; is_cooling=false, nodes=nodes)
    end
end

function optimize_parameter_subset_v5(objective, p_initial, indices, lower, upper;
                                      maximum_iterations=200,
                                      maximum_time_seconds=120.0,
                                      optimizer=NLopt.LN_NELDERMEAD(),
                                      label="subset")
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

function calibrate_v5(; cooling_keys=sim_key_cool,
                      heating_keys=sim_key_heat,
                      nodes=15,
                      heating_base_iterations=20,
                      cooling_iterations=20,
                      heating_refit_iterations=20,
                      cooling_time_seconds=120.0,
                      heating_base_time_seconds=120.0,
                      heating_refit_time_seconds=120.0,
                      optimizer=NLopt.LN_NELDERMEAD())
    p0 = copy(pnew_v5)

    println("\n=== 1D_v5 STAGE 1: heating base calibration, optical fixed ===")
    heating_base_fit = optimize_parameter_subset_v5(
        p -> loss_heating_v5(p, heating_keys; nodes=nodes),
        p0, fit_nonirradiance_indices_v5, lb_full_v5, ub_full_v5;
        maximum_iterations=nlopt_evaluation_budget_v3(
            heating_base_iterations, length(fit_nonirradiance_indices_v5)
        ),
        maximum_time_seconds=heating_base_time_seconds,
        optimizer=optimizer,
        label="heating-base-v5",
    )
    p_after_heating_base = heating_base_fit.parameters

    println("\n=== 1D_v5 STAGE 2: cooling thermophysical calibration ===")
    cooling_fit = optimize_parameter_subset_v5(
        p -> loss_cooling_v5(p, cooling_keys; nodes=nodes),
        p_after_heating_base, fit_cooling_indices_v5, lb_full_v5, ub_full_v5;
        maximum_iterations=nlopt_evaluation_budget_v3(
            cooling_iterations, length(fit_cooling_indices_v5)
        ),
        maximum_time_seconds=cooling_time_seconds,
        optimizer=optimizer,
        label="cooling-v5",
    )
    p_after_cooling = cooling_fit.parameters

    println("\n=== 1D_v5 STAGE 3: heating refit, optical fixed ===")
    heating_refit = optimize_parameter_subset_v5(
        p -> loss_heating_v5(p, heating_keys; nodes=nodes),
        p_after_cooling, fit_nonirradiance_indices_v5, lb_full_v5, ub_full_v5;
        maximum_iterations=nlopt_evaluation_budget_v3(
            heating_refit_iterations, length(fit_nonirradiance_indices_v5)
        ),
        maximum_time_seconds=heating_refit_time_seconds,
        optimizer=optimizer,
        label="heating-refit-v5",
    )
    parameters = heating_refit.parameters

    global pnew_v5 = parameters
    global calibration_v5 = (
        heating_base=heating_base_fit,
        cooling=cooling_fit,
        heating_refit=heating_refit,
        fixed_optical=(
            eta_abs=ETA_ABS_FIXED_V5,
            beta_opt=BETA_OPT_FIXED_V5,
            front_deposition_fraction=FRONT_DEPOSITION_FIXED_V5,
            eps_front=EPS_FRONT_FIXED_V5,
        ),
        parameters=parameters,
    )
    println("\nFinal 1D_v5 parameters:")
    println(parameters)
    return calibration_v5
end

function quick_calibration_v5()
    return calibrate_v5(
        cooling_keys=["C69"], heating_keys=["E74"],
        nodes=11,
        heating_base_iterations=2,
        cooling_iterations=2,
        heating_refit_iterations=2,
    )
end

begin # v5 post-processing
    function transient_case_v5(simulation_id, p=pnew_v5; is_cooling=false,
                               nodes=default_nodes)
        model, result, experiment = solve_case_v5(
            p, simulation_id; is_cooling=is_cooling, nodes=nodes
        )
        return (
            time=result.time,
            T8_model=model[:, 1], T8_experiment=experiment[:, 1],
            T9_pair_model=model[:, 2], T9_pair_experiment=experiment[:, 2],
            T10_pair_model=model[:, 3], T10_pair_experiment=experiment[:, 3],
            T3_model=model[:, 4], T3_experiment=experiment[:, 4],
            h_effective=model[:, 5],
            full_result=result,
        )
    end

    function plot_case_v5(simulation_id, params=pnew_v5; is_cooling=false,
                          nodes=default_nodes)
        @eval using StatsPlots
        data = transient_case_v5(simulation_id, params;
                                 is_cooling=is_cooling, nodes=nodes)
        sensors = (:T8, :T9_pair, :T10_pair, :T3)
        colors = (:blue, :red, :green, :purple)
        plot_object = plot(title="1D_v5 transient: $simulation_id",
                           xlabel="Time (s)", ylabel="Temperature (K)",
                           legend=:outerright)
        for (sensor, color) in zip(sensors, colors)
            plot!(plot_object, data.time, getproperty(data, Symbol(sensor, "_model")),
                  label="$(sensor) model", lw=2, color=color)
            scatter!(plot_object, data.time,
                     getproperty(data, Symbol(sensor, "_experiment")),
                     label="$(sensor) experiment", ms=2, markerstrokewidth=0,
                     color=color)
        end
        return plot_object
    end

    function compute_metrics_v5(p=pnew_v5; heating_keys=sim_key_heat,
                                cooling_keys=sim_key_cool, nodes=default_nodes)
        metrics = NamedTuple[]
        for (keys, cooling) in ((heating_keys, false), (cooling_keys, true))
            for simulation_id in keys
                model, result, experiment = solve_case_v5(
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
                    ))
                end
            end
        end
        return metrics
    end
end

println("""
[1D_v5] Model loaded.
Run quick_calibration_v5() for a short check.
Run calibration_v5 = calibrate_v5() for the full staged calibration.
The calibrated vector is stored in pnew_v5.
""")
