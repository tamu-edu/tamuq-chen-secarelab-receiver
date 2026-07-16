# ============================================================================
# 1D_v6.jl - Staged 1D receiver model with rear thermal mass
# ============================================================================
# v6 extends the reduced v5 model with one additional rear/adaptor thermal
# mass state. Calibration is staged by identifiability:
#   1. Heating fits gas heat-transfer shape.
#   2. Cooling fits thermophysical properties and hidden thermal masses.
#   3. Heating refits gas heat-transfer shape with thermophysics fixed.
#
# Optical/input-energy parameters and front convection are fixed so incoming
# energy can be adjusted manually and h_front remains literature-based.
# ============================================================================

include("1D_v5.jl")

begin # v6 fixed constants
    const B_RE_FIXED_V6 = B_RE_FIXED_V5
    const C_PR_FIXED_V6 = C_PR_FIXED_V5
    const EPS_FRONT_FIXED_V6 = EPS_FRONT_FIXED_V5
    const ETA_ABS_FIXED_V6 = ETA_ABS_FIXED_V5
    const BETA_OPT_FIXED_V6 = BETA_OPT_FIXED_V5
    const FRONT_DEPOSITION_FIXED_V6 = FRONT_DEPOSITION_FIXED_V5
    const H_FRONT_FIXED_V6 = 10.0
end

# Final v6 parameter vector p[1:11]
# Thermophysical / cooling-fitted parameters:
#   p[1]  gamma_C       effective receiver solid heat-capacity multiplier
#   p[2]  k_scale       effective axial-conductivity multiplier
#   p[3]  U_side        side-loss conductance per receiver length [W/m/K]
#   p[4]  U_rear        direct rear conductance from receiver end [W/K]
# Heat-transfer / heating-fitted parameters:
#   p[5]  A_Nu          Nusselt multiplier exponent in 10^A_Nu
#   p[6]  h_floor       downstream exchange floor, 0<h_floor<=1
#   p[7]  L_h           axial exchange decay length [m]
# Thermophysical / cooling-fitted parameters:
#   p[8]  tau_T3        outlet sensor/downstream lag [s]
#   p[9]  C_rear_scale  rear mass heat-capacity scale vs receiver solid
#   p[10] K_rear        receiver-end to rear-mass conductance [W/K]
#   p[11] U_rear_mass   rear-mass loss conductance to ambient [W/K]
begin # v6 parameter values and bounds
    pnew_v6 = [
        1.0, 1.0, 0.35, 0.10,
        -1.0, 0.20, 0.050,
        20.0, 0.50, 0.20, 0.05,
    ]
    lb_full_v6 = [
        0.10, 0.10, 0.005, 0.000,
        -4.0, 0.02, 0.005,
        0.50, 0.01, 0.001, 0.000,
    ]
    ub_full_v6 = [
        8.00, 5.00, 5.000, 5.000,
        2.00, 1.00, 0.500,
        500.0, 10.0, 10.0, 5.000,
    ]

    fit_heat_transfer_indices_v6 = [5, 6, 7]
    fit_cooling_thermal_indices_v6 = [1, 2, 3, 4, 8, 9, 10, 11]
end

rear_capacity_v6(Trear, p) =
    p[9] * rho_s * Cps_f(Trear) * A_solid * L

function gas_profile_v6!(Tg, Qgas, hcell, Ts, time, p, operating, z, dx)
    flow = max(0.0, operating.flow_lpm(time))
    Tin = operating.inlet_temperature(time)
    Tg[1] = Tin

    if flow <= 1e-12
        fill!(Qgas, 0.0)
        fill!(hcell, 0.0)
        Tg[1] = Tin
        for i in eachindex(Ts)
            Tg[i + 1] = Ts[i]
        end
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
                  reynolds^B_RE_FIXED_V6 *
                  prandtl^(10.0^C_PR_FIXED_V6)
        hcell[i] = nusselt * kf / Dh * axial_exchange_shape_v5(z[i], p)
        UA = hcell[i] * P_exchange * dx
        effectiveness = -expm1(-UA / (mdot * cp))
        Tg[i + 1] = Tg[i] + effectiveness * (Ts[i] - Tg[i])
        Qgas[i] = mdot * cp * (Tg[i + 1] - Tg[i])
    end
    return nothing
end

function receiver_rhs_v6!(du, u, time, p, operating, z, dx, solar_weights,
                          Tg, Qgas, hcell)
    nodes = length(z)
    Ts = view(u, 1:nodes)
    dTs = view(du, 1:nodes)
    Trear = u[nodes + 1]
    ambient = operating.ambient_temperature(time)

    gas_profile_v6!(Tg, Qgas, hcell, Ts, time, p, operating, z, dx)
    fill!(du, 0.0)

    for i in 1:(nodes - 1)
        ki = p[2] * ks_f(Ts[i])
        kj = p[2] * ks_f(Ts[i + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        Qcond = kface * A_solid * (Ts[i] - Ts[i + 1]) / dx
        dTs[i] -= Qcond
        dTs[i + 1] += Qcond
    end

    Qsolar = ETA_ABS_FIXED_V6 * max(0.0, operating.irradiance(time)) * A_frt
    for i in 1:nodes
        Qside = p[3] * dx * (Ts[i] - ambient)
        dTs[i] += Qsolar * solar_weights[i] - Qgas[i] - Qside
    end

    Qfront = H_FRONT_FIXED_V6 * A_frt * (Ts[1] - ambient) +
             EPS_FRONT_FIXED_V6 * sigma_sb * A_frt * (Ts[1]^4 - ambient^4)
    Qrear_direct = p[4] * (Ts[end] - ambient)
    Qrear_coupling = p[10] * (Ts[end] - Trear)
    Qrear_mass_loss = p[11] * (Trear - ambient)

    dTs[1] -= Qfront
    dTs[end] -= Qrear_direct + Qrear_coupling
    du[nodes + 1] = (Qrear_coupling - Qrear_mass_loss) /
                    max(rear_capacity_v6(Trear, p), eps(Float64))

    cell_volume = A_solid * dx
    for i in 1:nodes
        capacity = rho_s * p[1] * Cps_f(Ts[i]) * cell_volume
        dTs[i] /= capacity
    end
    return nothing
end

function receiver_ode_v6!(du, u, context, time)
    receiver_rhs_v6!(
        du, u, time, context.p, context.operating,
        context.z, context.dx, context.solar_weights,
        context.Tg, context.Qgas, context.hcell,
    )
    return nothing
end

function simulate_v6(p, operating, save_times;
                     initial_temperature=Tamb, nodes=default_nodes,
                     initial_rear_temperature=nothing,
                     solver=Rodas5P(autodiff=AutoFiniteDiff()),
                     reltol=1e-6, abstol=1e-7, dtmax=30.0)
    length(p) == 11 || throw(ArgumentError("1D_v6 expects 11 parameters"))
    all(isfinite, p) || throw(ArgumentError("all model parameters must be finite"))
    p[1] > 0.0 && p[2] > 0.0 && p[9] > 0.0 ||
        throw(ArgumentError("capacity, conductivity, and rear mass must be positive"))
    nodes >= 3 || throw(ArgumentError("at least three axial nodes are required"))

    time = Float64.(save_times)
    isempty(time) && throw(ArgumentError("save_times cannot be empty"))
    (length(time) == 1 || all(diff(time) .> 0.0)) ||
        throw(ArgumentError("save_times must be strictly increasing"))

    dx = L / nodes
    z = collect(range(dx / 2, step=dx, length=nodes))
    zg = collect(range(0.0, L, length=nodes + 1))
    weights = solar_weights_v5(BETA_OPT_FIXED_V6, FRONT_DEPOSITION_FIXED_V6, nodes)
    Ts_initial = initial_profile_v3(initial_temperature, z)
    Trear_initial = isnothing(initial_rear_temperature) ? Ts_initial[end] :
                    Float64(initial_rear_temperature)
    u_initial = vcat(Ts_initial, Trear_initial)

    Tg_history = Matrix{Float64}(undef, nodes + 1, length(time))
    h_history = Matrix{Float64}(undef, nodes, length(time))
    rear_history = Vector{Float64}(undef, length(time))
    Tg = zeros(nodes + 1)
    Qgas = zeros(nodes)
    hcell = zeros(nodes)

    context = (
        p=p, operating=operating, z=z, dx=dx, solar_weights=weights,
        Tg=Tg, Qgas=Qgas, hcell=hcell,
    )
    problem = ODEProblem(receiver_ode_v6!, u_initial, (time[1], time[end]), context)
    solution = solve(
        problem, solver; saveat=time, save_everystep=false, dense=false,
        reltol=reltol, abstol=abstol, dtmax=dtmax,
    )
    successful_retcode(solution) ||
        error("1D_v6 ODE solve failed with retcode $(solution.retcode)")

    state_history = reduce(hcat, solution.u)
    Ts_history = state_history[1:nodes, :]
    rear_history .= state_history[nodes + 1, :]
    for output_index in eachindex(time)
        gas_profile_v6!(Tg, Qgas, hcell, view(Ts_history, :, output_index),
                        time[output_index], p, operating, z, dx)
        Tg_history[:, output_index] .= Tg
        h_history[:, output_index] .= hcell
    end

    return (
        time=time, z_solid=z, z_gas=zg,
        solid_temperature=Ts_history,
        rear_temperature=rear_history,
        gas_temperature=Tg_history,
        heat_transfer_coefficient=h_history,
        ode_solution=solution,
    )
end

begin # v6 experimental interface
    paired_observation_v6(data, simulation_id, obs_a, obs_b) =
        0.5 .* (observation(data, simulation_id, obs_a) .+
                observation(data, simulation_id, obs_b))

    function measured_initial_profile_v6(data, simulation_id)
        z8 = sensor_positions[:T8]
        z9 = sensor_positions[:T9]
        z10 = sensor_positions[:T10]
        T8 = observation(data, simulation_id, "_T8")[1]
        T9 = paired_observation_v6(data, simulation_id, "_T9", "_T12")[1]
        T10 = paired_observation_v6(data, simulation_id, "_T10", "_T11")[1]
        return function (z)
            z <= z8 && return T8
            z >= z10 && return T10
            if z <= z9
                return T8 + (T9 - T8) * (z - z8) / (z9 - z8)
            end
            return T9 + (T10 - T9) * (z - z9) / (z10 - z9)
        end
    end

    function experimental_case_v6(simulation_id; is_cooling=false)
        data = is_cooling ? measurements_cooling : measurements
        conditions = is_cooling ? simulation_conditions_cooling : simulation_conditions
        time = observation_time(data, simulation_id)
        experiment = hcat(
            observation(data, simulation_id, "_T8"),
            paired_observation_v6(data, simulation_id, "_T9", "_T12"),
            paired_observation_v6(data, simulation_id, "_T10", "_T11"),
            observation(data, simulation_id, "_Tf"),
        )
        return conditions[simulation_id], time, experiment, data
    end

    initial_rear_temperature_v6(data, simulation_id) =
        observation(data, simulation_id, "_T2")[1]

    function extract_outputs_v6(result, p, T3_initial)
        T8 = solid_at_v3(result, sensor_positions[:T8])
        T9_pair = solid_at_v3(result, sensor_positions[:T9])
        T10_pair = solid_at_v3(result, sensor_positions[:T10])
        Tf_true = vec(result.gas_temperature[end, :])
        Tf = apply_T3_lag_v3(Tf_true, T3_initial, result.time, p[8])
        h_mean = vec(mean(result.heat_transfer_coefficient, dims=1))
        return hcat(T8, T9_pair, T10_pair, Tf, h_mean)
    end

    function solve_case_v6(p, simulation_id; is_cooling=false,
                           nodes=default_nodes,
                           solver=Rodas5P(autodiff=AutoFiniteDiff()),
                           reltol=1e-6, abstol=1e-7, dtmax=30.0)
        conditions, time, experiment, data = experimental_case_v6(
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
        result = simulate_v6(
            p, operating, time;
            initial_temperature=measured_initial_profile_v6(data, simulation_id),
            initial_rear_temperature=initial_rear_temperature_v6(data, simulation_id),
            nodes=nodes, solver=solver, reltol=reltol, abstol=abstol, dtmax=dtmax,
        )
        outputs = extract_outputs_v6(result, p, experiment[1, 4])
        return outputs, result, experiment
    end
end

begin # v6 calibration
    heat_shape_regularization_v6(p) =
        0.001 * (log(p[6] / 0.20)^2 + log(p[7] / 0.050)^2)

    function loss_cases_v6(p, keys; is_cooling=false, nodes=15)
        total = 0.0
        for simulation_id in keys
            try
                model, _, experiment = solve_case_v6(
                    p, simulation_id; is_cooling=is_cooling, nodes=nodes,
                    reltol=2e-5, abstol=2e-6, dtmax=60.0,
                )
                all(isfinite, model) || return Inf
                total += mean(normalized_mse_v3(model[:, j], experiment[:, j]) for j in 1:4)
            catch err
                err isa InterruptException && rethrow()
                return Inf
            end
        end
        return total / length(keys)
    end

    loss_heating_v6(p, keys=sim_key_heat; nodes=15) =
        loss_cases_v6(p, keys; is_cooling=false, nodes=nodes) +
        heat_shape_regularization_v6(p)

    function loss_cooling_v6(p, keys=sim_key_cool; nodes=15)
        regularization = 0.005 * (log(p[1])^2 + log(p[2])^2) +
                         0.001 * log(p[9] / 0.50)^2
        return loss_cases_v6(p, keys; is_cooling=true, nodes=nodes) + regularization
    end
end

function optimize_parameter_subset_v6(objective, p_initial, indices, lower, upper;
                                      maximum_iterations=200,
                                      maximum_time_seconds=120.0,
                                      optimizer=NLopt.LN_NELDERMEAD(),
                                      label="subset-v6")
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

function calibrate_v6(; cooling_keys=sim_key_cool,
                      heating_keys=sim_key_heat,
                      nodes=15,
                      heating_base_iterations=20,
                      cooling_iterations=20,
                      heating_refit_iterations=20,
                      cooling_time_seconds=120.0,
                      heating_base_time_seconds=120.0,
                      heating_refit_time_seconds=120.0,
                      optimizer=NLopt.LN_NELDERMEAD())
    p0 = copy(pnew_v6)

    println("\n=== 1D_v6 STAGE 1: heating heat-transfer calibration ===")
    heating_base_fit = optimize_parameter_subset_v6(
        p -> loss_heating_v6(p, heating_keys; nodes=nodes),
        p0, fit_heat_transfer_indices_v6, lb_full_v6, ub_full_v6;
        maximum_iterations=nlopt_evaluation_budget_v3(
            heating_base_iterations, length(fit_heat_transfer_indices_v6)
        ),
        maximum_time_seconds=heating_base_time_seconds,
        optimizer=optimizer,
        label="heating-base-v6",
    )
    p_after_heating_base = heating_base_fit.parameters

    println("\n=== 1D_v6 STAGE 2: cooling thermophysical/rear-mass calibration ===")
    cooling_fit = optimize_parameter_subset_v6(
        p -> loss_cooling_v6(p, cooling_keys; nodes=nodes),
        p_after_heating_base, fit_cooling_thermal_indices_v6, lb_full_v6, ub_full_v6;
        maximum_iterations=nlopt_evaluation_budget_v3(
            cooling_iterations, length(fit_cooling_thermal_indices_v6)
        ),
        maximum_time_seconds=cooling_time_seconds,
        optimizer=optimizer,
        label="cooling-v6",
    )
    p_after_cooling = cooling_fit.parameters

    println("\n=== 1D_v6 STAGE 3: heating heat-transfer refit ===")
    heating_refit = optimize_parameter_subset_v6(
        p -> loss_heating_v6(p, heating_keys; nodes=nodes),
        p_after_cooling, fit_heat_transfer_indices_v6, lb_full_v6, ub_full_v6;
        maximum_iterations=nlopt_evaluation_budget_v3(
            heating_refit_iterations, length(fit_heat_transfer_indices_v6)
        ),
        maximum_time_seconds=heating_refit_time_seconds,
        optimizer=optimizer,
        label="heating-refit-v6",
    )
    parameters = heating_refit.parameters

    global pnew_v6 = parameters
    global calibration_v6 = (
        heating_base=heating_base_fit,
        cooling=cooling_fit,
        heating_refit=heating_refit,
        fixed=(
            eta_abs=ETA_ABS_FIXED_V6,
            beta_opt=BETA_OPT_FIXED_V6,
            front_deposition_fraction=FRONT_DEPOSITION_FIXED_V6,
            h_front=H_FRONT_FIXED_V6,
            eps_front=EPS_FRONT_FIXED_V6,
        ),
        parameters=parameters,
    )
    println("\nFinal 1D_v6 parameters:")
    println(parameters)
    return calibration_v6
end

function quick_calibration_v6()
    return calibrate_v6(
        cooling_keys=["C69"], heating_keys=["E74"],
        nodes=11,
        heating_base_iterations=2,
        cooling_iterations=2,
        heating_refit_iterations=2,
    )
end

begin # v6 post-processing
    function transient_case_v6(simulation_id, p=pnew_v6; is_cooling=false,
                               nodes=default_nodes)
        model, result, experiment = solve_case_v6(
            p, simulation_id; is_cooling=is_cooling, nodes=nodes
        )
        return (
            time=result.time,
            T8_model=model[:, 1], T8_experiment=experiment[:, 1],
            T9_pair_model=model[:, 2], T9_pair_experiment=experiment[:, 2],
            T10_pair_model=model[:, 3], T10_pair_experiment=experiment[:, 3],
            T3_model=model[:, 4], T3_experiment=experiment[:, 4],
            T_rear=result.rear_temperature,
            h_effective=model[:, 5],
            full_result=result,
        )
    end

    function plot_case_v6(simulation_id, params=pnew_v6; is_cooling=false,
                          nodes=default_nodes)
        @eval using StatsPlots
        data = transient_case_v6(simulation_id, params;
                                 is_cooling=is_cooling, nodes=nodes)
        sensors = (:T8, :T9_pair, :T10_pair, :T3)
        colors = (:blue, :red, :green, :purple)
        plot_object = plot(title="1D_v6 transient: $simulation_id",
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
        plot!(plot_object, data.time, data.T_rear;
              label="rear mass", color=:black, linestyle=:dash)
        return plot_object
    end

    function compute_metrics_v6(p=pnew_v6; heating_keys=sim_key_heat,
                                cooling_keys=sim_key_cool, nodes=default_nodes)
        metrics = NamedTuple[]
        for (keys, cooling) in ((heating_keys, false), (cooling_keys, true))
            for simulation_id in keys
                model, result, experiment = solve_case_v6(
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
[1D_v6] Model loaded.
Run quick_calibration_v6() for a short check.
Run calibration_v6 = calibrate_v6() for the full staged calibration.
The calibrated vector is stored in pnew_v6.
""")
