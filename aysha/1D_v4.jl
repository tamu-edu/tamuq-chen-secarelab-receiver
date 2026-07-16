# ============================================================================
# 1D_v4.jl - 1D receiver model with axial heat-exchange shaping
# ============================================================================
# v4 extends the v3 finite-volume model with:
#   1. A Graetz-style axial multiplier on the A/B/C Nu(Re,Pr) gas HTC law.
#   2. Paired thermocouple comparison channels:
#        T9_pair  = mean(T9, T12)
#        T10_pair = mean(T10, T11)
#
# The radial averaging is a data-interface choice: the 1D model has no radial
# degree of freedom, so paired thermocouples at the same axial station are
# treated as one axial measurement.
# ============================================================================

include("1D_v3.jl")

# Final v4 parameter vector p[1:15]
# Cooling / flow-exchange parameters:
#   p[1]  gamma_C     effective solid heat-capacity multiplier
#   p[2]  k_scale     effective axial-conductivity multiplier
#   p[3]  U_side      side-loss conductance per receiver length [W/m/K]
#   p[4]  U_rear      rear/adaptor conductance [W/K]
#   p[5]  A_Nu        Nusselt multiplier exponent in 10^A_Nu
#   p[6]  B_Re        Reynolds exponent in Re^B_Re
#   p[7]  C_Pr        Prandtl double-exponent in Pr^(10^C_Pr)
#   p[8]  a_Gz        Graetz shape amplitude in 1 - a_Gz*Gz^n_Gz
#   p[9]  c_Gz        Graetz exponential damping in exp(-c_Gz/Gz)
#   p[10] n_Gz        Graetz power exponent
#   p[11] tau_T3      outlet sensor/downstream lag [s]
# Heating parameters:
#   p[12] eta_abs     absorbed solar fraction
#   p[13] beta_opt    Beer-Lambert extinction coefficient [1/m]
#   p[14] h_front     front convection coefficient [W/m2/K]
#   p[15] eps_front   effective front emissivity
begin # v4 parameter values and bounds
    p_cool_init_v4 = [1.0, 1.0, 0.35, 0.10, -1.0, 1.0, 0.5, 0.0, 0.0, 0.0, 20.0]
    lb_cool_v4 = [0.10, 0.10, 0.01, 0.001, -3.0, 0.50, -1.0, -1.0, 0.0, -1.0, 0.50]
    ub_cool_v4 = [5.00, 5.00, 5.00, 5.000, 1.00, 2.00, 2.00, 1.0, 50.0, 1.0, 500.0]

    p_heat_init_v4 = copy(p_heat_init)
    lb_heat_v4 = copy(lb_heat)
    ub_heat_v4 = copy(ub_heat)

    assemble_parameters_v4(p_cool, p_heat) = vcat(Float64.(p_cool), Float64.(p_heat))
    pnew_v4 = assemble_parameters_v4(p_cool_init_v4, p_heat_init_v4)
end

graetz_number_v4(reynolds, prandtl, z) =
    reynolds * prandtl * Dh / max(Float64(z), Dh)

function axial_exchange_shape_v4(graetz, p)
    gz = max(Float64(graetz), eps(Float64))
    raw = (1.0 - p[8] * gz^p[10]) * exp(-p[9] / gz)
    return clamp(raw, 0.05, 5.0)
end

function gas_profile_v4!(Tg, Qgas, hcell, Ts, time, p, operating, z, dx)
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
        graetz = graetz_number_v4(reynolds, prandtl, z[i])
        nusselt_base = 10.0^p[5] * reynolds^p[6] * prandtl^(10.0^p[7])
        nusselt = nusselt_base * axial_exchange_shape_v4(graetz, p)
        hcell[i] = nusselt * kf / Dh
        UA = hcell[i] * P_exchange * dx
        effectiveness = -expm1(-UA / (mdot * cp))
        Tg[i + 1] = Tg[i] + effectiveness * (Ts[i] - Tg[i])
        Qgas[i] = mdot * cp * (Tg[i + 1] - Tg[i])
    end
    return nothing
end

function receiver_rhs_v4!(dTs, Ts, time, p, operating, z, dx, solar_weights,
                          Tg, Qgas, hcell)
    nodes = length(Ts)
    ambient = operating.ambient_temperature(time)
    gas_profile_v4!(Tg, Qgas, hcell, Ts, time, p, operating, z, dx)
    fill!(dTs, 0.0)

    for i in 1:(nodes - 1)
        ki = p[2] * ks_f(Ts[i])
        kj = p[2] * ks_f(Ts[i + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        Qcond = kface * A_solid * (Ts[i] - Ts[i + 1]) / dx
        dTs[i] -= Qcond
        dTs[i + 1] += Qcond
    end

    Qsolar = p[12] * max(0.0, operating.irradiance(time)) * A_frt
    for i in 1:nodes
        Qside = p[3] * dx * (Ts[i] - ambient)
        dTs[i] += Qsolar * solar_weights[i] - Qgas[i] - Qside
    end

    Qfront = p[14] * A_frt * (Ts[1] - ambient) +
             p[15] * sigma_sb * A_frt * (Ts[1]^4 - ambient^4)
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

function receiver_ode_v4!(dTs, Ts, context, time)
    receiver_rhs_v4!(
        dTs, Ts, time, context.p, context.operating,
        context.z, context.dx, context.solar_weights,
        context.Tg, context.Qgas, context.hcell,
    )
    return nothing
end

function simulate_v4(p, operating, save_times;
                     initial_temperature=Tamb, nodes=default_nodes,
                     solver=Rodas5P(autodiff=AutoFiniteDiff()),
                     reltol=1e-6, abstol=1e-7, dtmax=30.0)
    length(p) == 15 || throw(ArgumentError("1D_v4 expects 15 parameters"))
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
    weights = solar_weights_v3(p[13], nodes)
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
    problem = ODEProblem(receiver_ode_v4!, Ts_initial, (time[1], time[end]), context)
    solution = solve(
        problem, solver; saveat=time, save_everystep=false, dense=false,
        reltol=reltol, abstol=abstol, dtmax=dtmax,
    )
    successful_retcode(solution) ||
        error("1D_v4 ODE solve failed with retcode $(solution.retcode)")

    Ts_history = reduce(hcat, solution.u)
    for output_index in eachindex(time)
        gas_profile_v4!(Tg, Qgas, hcell, view(Ts_history, :, output_index),
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

begin # v4 experimental interface
    function paired_observation_v4(data, simulation_id, obs_a, obs_b)
        return 0.5 .* (observation(data, simulation_id, obs_a) .+
                       observation(data, simulation_id, obs_b))
    end

    function measured_initial_profile_v4(data, simulation_id)
        z8 = sensor_positions[:T8]
        z9 = sensor_positions[:T9]
        z10 = sensor_positions[:T10]
        T8 = observation(data, simulation_id, "_T8")[1]
        T9 = paired_observation_v4(data, simulation_id, "_T9", "_T12")[1]
        T10 = paired_observation_v4(data, simulation_id, "_T10", "_T11")[1]
        return function (z)
            z <= z8 && return T8
            z >= z10 && return T10
            if z <= z9
                return T8 + (T9 - T8) * (z - z8) / (z9 - z8)
            end
            return T9 + (T10 - T9) * (z - z9) / (z10 - z9)
        end
    end

    function experimental_case_v4(simulation_id; is_cooling=false)
        data = is_cooling ? measurements_cooling : measurements
        conditions = is_cooling ? simulation_conditions_cooling : simulation_conditions
        time = observation_time(data, simulation_id)
        experiment = hcat(
            observation(data, simulation_id, "_T8"),
            paired_observation_v4(data, simulation_id, "_T9", "_T12"),
            paired_observation_v4(data, simulation_id, "_T10", "_T11"),
            observation(data, simulation_id, "_Tf"),
        )
        return conditions[simulation_id], time, experiment, data
    end

    function extract_outputs_v4(result, p, T3_initial)
        T8 = solid_at_v3(result, sensor_positions[:T8])
        T9_pair = solid_at_v3(result, sensor_positions[:T9])
        T10_pair = solid_at_v3(result, sensor_positions[:T10])
        Tf_true = vec(result.gas_temperature[end, :])
        Tf = apply_T3_lag_v3(Tf_true, T3_initial, result.time, p[11])
        h_mean = vec(mean(result.heat_transfer_coefficient, dims=1))
        return hcat(T8, T9_pair, T10_pair, Tf, h_mean)
    end

    function solve_case_v4(p, simulation_id; is_cooling=false,
                           nodes=default_nodes,
                           solver=Rodas5P(autodiff=AutoFiniteDiff()),
                           reltol=1e-6, abstol=1e-7, dtmax=30.0)
        conditions, time, experiment, data = experimental_case_v4(
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
        result = simulate_v4(
            p, operating, time;
            initial_temperature=measured_initial_profile_v4(data, simulation_id),
            nodes=nodes, solver=solver, reltol=reltol, abstol=abstol, dtmax=dtmax,
        )
        outputs = extract_outputs_v4(result, p, experiment[1, 4])
        return outputs, result, experiment
    end
end

begin # v4 calibration
    function loss_cases_v4(p, keys; is_cooling=false, nodes=15)
        total = 0.0
        for simulation_id in keys
            model, _, experiment = solve_case_v4(
                p, simulation_id; is_cooling=is_cooling, nodes=nodes,
                reltol=2e-5, abstol=2e-6, dtmax=60.0,
            )
            all(isfinite, model) || return Inf
            total += mean(normalized_mse_v3(model[:, j], experiment[:, j]) for j in 1:4)
        end
        return total / length(keys)
    end

    function loss_cooling_v4(p_cool, keys=sim_key_cool; nodes=15)
        p = assemble_parameters_v4(p_cool, p_heat_init_v4)
        regularization = 0.005 * (log(p_cool[1])^2 + log(p_cool[2])^2) +
                         0.001 * (p_cool[8]^2 + (p_cool[9] / 10.0)^2 + p_cool[10]^2)
        return loss_cases_v4(p, keys; is_cooling=true, nodes=nodes) + regularization
    end

    function loss_heating_v4(p_heat, keys=sim_key_heat;
                             p_cool=p_cool_init_v4, nodes=15)
        p = assemble_parameters_v4(p_cool, p_heat)
        return loss_cases_v4(p, keys; is_cooling=false, nodes=nodes)
    end
end

function calibrate_v4(; cooling_keys=sim_key_cool,
                      heating_keys=sim_key_heat,
                      nodes=15,
                      cooling_iterations=20,
                      heating_iterations=20,
                      cooling_time_seconds=120.0,
                      heating_time_seconds=120.0,
                      optimizer=NLopt.LN_NELDERMEAD())
    println("\n=== 1D_v4 STAGE 1: cooling calibration ===")
    cooling_fit = optimize_with_nlopt_v3(
        p -> loss_cooling_v4(p, cooling_keys; nodes=nodes),
        p_cool_init_v4, lb_cool_v4, ub_cool_v4;
        maximum_iterations=nlopt_evaluation_budget_v3(cooling_iterations, length(p_cool_init_v4)),
        maximum_time_seconds=cooling_time_seconds,
        algorithm=optimizer,
        label="cooling-v4",
    )
    p_cool = cooling_fit.minimizer
    println("Cooling v4 parameters: $p_cool")

    println("\n=== 1D_v4 STAGE 2: heating/optical calibration ===")
    heating_fit = optimize_with_nlopt_v3(
        p -> loss_heating_v4(p, heating_keys; p_cool=p_cool, nodes=nodes),
        p_heat_init_v4, lb_heat_v4, ub_heat_v4;
        maximum_iterations=nlopt_evaluation_budget_v3(heating_iterations, length(p_heat_init_v4)),
        maximum_time_seconds=heating_time_seconds,
        algorithm=optimizer,
        label="heating-v4",
    )
    p_heat = heating_fit.minimizer
    parameters = assemble_parameters_v4(p_cool, p_heat)

    global pnew_v4 = parameters
    global calibration_v4 = (
        cooling=cooling_fit, heating=heating_fit,
        p_cool=p_cool, p_heat=p_heat, parameters=parameters,
    )
    println("\nFinal 1D_v4 parameters:")
    println(parameters)
    return calibration_v4
end

function quick_calibration_v4()
    return calibrate_v4(
        cooling_keys=["C69"], heating_keys=["E74"],
        nodes=11, cooling_iterations=2, heating_iterations=2,
    )
end

begin # v4 post-processing
    function transient_case_v4(simulation_id, p=pnew_v4; is_cooling=false,
                               nodes=default_nodes)
        model, result, experiment = solve_case_v4(
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

    function plot_case_v4(simulation_id, params=pnew_v4; is_cooling=false,
                          nodes=default_nodes)
        @eval using StatsPlots
        data = transient_case_v4(simulation_id, params;
                                 is_cooling=is_cooling, nodes=nodes)
        sensors = (:T8, :T9_pair, :T10_pair, :T3)
        colors = (:blue, :red, :green, :purple)
        plot_object = plot(title="1D_v4 transient: $simulation_id",
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

    function compute_metrics_v4(p=pnew_v4; heating_keys=sim_key_heat,
                                cooling_keys=sim_key_cool, nodes=default_nodes)
        metrics = NamedTuple[]
        for (keys, cooling) in ((heating_keys, false), (cooling_keys, true))
            for simulation_id in keys
                model, result, experiment = solve_case_v4(
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
[1D_v4] Model loaded.
Run quick_calibration_v4() for a short check.
Run calibration_v4 = calibrate_v4() for the full cooling/heating calibration.
The calibrated vector is stored in pnew_v4.
""")
