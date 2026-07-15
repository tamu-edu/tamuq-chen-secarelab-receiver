# ============================================================================
# 1D_v3.jl - Conservative finite-volume solid/gas receiver model
# ============================================================================
# Clean research-script flow:
#   libraries -> geometry/properties -> model -> external data import
#   -> solver/remake interface -> staged calibration -> post-processing.
#
# Default execution only imports data and defines functions. Set environment
# variables RECEIVER1D_V3_RUN_QUICK_CALIBRATION, RECEIVER1D_V3_RUN_CALIBRATION,
# RECEIVER1D_V3_RUN_POSTPROCESS, RECEIVER1D_V3_RUN_PLOTS, or
# RECEIVER1D_V3_EXPORT_METRICS to "true" for optional work.
#
# Coordinate: z = 0 at the irradiated air inlet; z = L at the gas outlet.
# Dynamic states: N axial solid temperatures.
# Gas: quasi-steady plug flow, marched cell-by-cell with an exact NTU update.
# ============================================================================

begin # libraries
    using DifferentialEquations
    using LinearAlgebra
    using Optimization
    using OptimizationNLopt
    using SciMLBase
    using Statistics
end

begin # execution flags
    env_flag_v3(name) = lowercase(get(ENV, name, "false")) == "true"

    const RUN_QUICK_CALIBRATION_V3 = env_flag_v3("RECEIVER1D_V3_RUN_QUICK_CALIBRATION")
    const RUN_CALIBRATION_V3 = env_flag_v3("RECEIVER1D_V3_RUN_CALIBRATION")
    const RUN_POSTPROCESS_V3 = env_flag_v3("RECEIVER1D_V3_RUN_POSTPROCESS")
    const RUN_PLOTS_V3 = env_flag_v3("RECEIVER1D_V3_RUN_PLOTS")
    const EXPORT_METRICS_V3 = env_flag_v3("RECEIVER1D_V3_EXPORT_METRICS")
end

# ============================================================================
# SECTION A: GEOMETRY, CONSTANTS, AND PARAMETER ORDER
# ============================================================================
begin # fixed parameters
    Tamb = 22.448 + 273.15
    sigma_sb = 5.670374419e-8

    w_t = 19.0e-3
    A_frt = w_t^2
    L = 137.0e-3
    w_chnl = 1.5e-3
    t_chnl = 0.4e-3
    n_chnl = 100
    A_flow = n_chnl * w_chnl^2
    A_solid = A_frt - A_flow
    Dh = w_chnl
    P_exchange = 4.0 * w_chnl * n_chnl
    rho_s = 3200.0

    sensor_positions = Dict(:T8 => 0.011, :T9 => 0.058, :T10 => 0.107)
    default_nodes = 25

    # Keys used by the external importer and legacy-style condition dictionaries.
    if !isdefined(@__MODULE__, :Io)
        Io = :Io
    end
    if !isdefined(@__MODULE__, :qlpm)
        qlpm = :qlpm
    end
    if !isdefined(@__MODULE__, :Tinit)
        Tinit = :Tinit
    end
    if !isdefined(@__MODULE__, :T_in)
        T_in = :T_in
    end
    if !isdefined(@__MODULE__, :T_amb)
        T_amb = :T_amb
    end
end;

# Final parameter vector p[1:11]
# Cooling parameters:
#   p[1]  gamma_C     effective solid heat-capacity multiplier
#   p[2]  k_scale     effective axial-conductivity multiplier
#   p[3]  U_side      side-loss conductance per receiver length [W/m/K]
#   p[4]  U_rear      rear/adaptor conductance [W/K]
#   p[5]  A_Nu        Nusselt multiplier exponent in 10^A_Nu
#   p[6]  B_Re        Reynolds exponent in Re^B_Re
#   p[7]  C_Pr        Prandtl double-exponent in Pr^(10^C_Pr)
#   p[8]  tau_T3      outlet sensor/downstream lag [s]
# Heating parameters:
#   p[9]  eta_abs     absorbed solar fraction
#   p[10] beta_opt    Beer-Lambert extinction coefficient [1/m]
#   p[11] h_front     front convection coefficient [W/m2/K]
#   p[12] eps_front   effective front emissivity
begin # parameter values and bounds
    p_cool_init = [1.0, 1.0, 0.35, 0.10, -1.0, 1.0, 0.5, 20.0]
    lb_cool = [0.10, 0.10, 0.01, 0.001, -3.0, 0.50, -1.0, 0.50]
    ub_cool = [5.00, 5.00, 5.00, 5.000, 1.00, 2.00, 2.00, 500.0]

    p_heat_init = [0.85, 50.0, 10.0, 0.85]
    lb_heat = [0.10, 1.0, 0.10, 0.10]
    ub_heat = [0.99, 1000.0, 100.0, 1.00]

    assemble_parameters_v3(p_cool, p_heat) = vcat(Float64.(p_cool), Float64.(p_heat))
    pnew = assemble_parameters_v3(p_cool_init, p_heat_init)
end

# ============================================================================
# SECTION B: TEMPERATURE-DEPENDENT PROPERTIES
# ============================================================================
begin # material and fluid properties
    property_temperature(T) = clamp(Float64(T), 250.0, 1600.0)

    function rho_f_f(T)
        return 352.716 / property_temperature(T)
    end

    function mu_f_f(T)
        x = property_temperature(T)
        return max(5e-6, -8.38278e-7 + 8.35717342e-8x - 7.69429583e-11x^2 +
                        4.6437266e-14x^3 - 1.06585607e-17x^4)
    end

    function kf_f(T)
        x = property_temperature(T)
        return max(0.01, -0.00227583562 + 1.15480022e-4x - 7.90252856e-8x^2 +
                         4.11702505e-11x^3 - 7.43864331e-15x^4)
    end

    function cpf_f(T)
        x = property_temperature(T)
        return max(700.0, 1.93e-10x^4 - 8.0e-7x^3 + 1.14e-3x^2 - 0.449x + 1060.0)
    end

    function Cps_f(T)
        x = property_temperature(T)
        return max(400.0, 1110.0 + 0.15x - 425.0exp(-0.003x))
    end

    function ks_f(T)
        x = property_temperature(T)
        return max(2.0, 191.9216 - 0.3261784x + 2.739462e-4x^2 - 7.70926e-8x^3)
    end

    m_dot(flow_lpm, Tin=Tamb) = rho_f_f(Tin) * max(0.0, flow_lpm) / 60000.0
end;

# ============================================================================
# SECTION C: FINITE-VOLUME MODEL
# ============================================================================
begin # input histories and grid helpers
    constant_history_v3(value::Real) = (_ -> Float64(value))
    constant_history_v3(f) = f

    function linear_history_v3(time, values)
        t = Float64.(time)
        y = Float64.(values)
        length(t) == length(y) || throw(ArgumentError("history time/value lengths differ"))
        issorted(t) || throw(ArgumentError("history time must be sorted"))
        return function (x)
            x <= t[1] && return y[1]
            x >= t[end] && return y[end]
            i = searchsortedlast(t, x)
            fraction = (x - t[i]) / (t[i + 1] - t[i])
            return y[i] + fraction * (y[i + 1] - y[i])
        end
    end

    function operating_condition_v3(; irradiance=0.0, flow_lpm=0.0,
                                     inlet_temperature=Tamb,
                                     ambient_temperature=Tamb)
        return (
            irradiance=constant_history_v3(irradiance),
            flow_lpm=constant_history_v3(flow_lpm),
            inlet_temperature=constant_history_v3(inlet_temperature),
            ambient_temperature=constant_history_v3(ambient_temperature),
        )
    end

    function solar_weights_v3(beta, nodes)
        dx = L / nodes
        weights = zeros(nodes)
        if beta <= 1e-12
            fill!(weights, 1.0 / nodes)
        else
            for i in 1:nodes
                weights[i] = exp(-beta * (i - 1) * dx) - exp(-beta * i * dx)
            end
            weights ./= sum(weights)
        end
        return weights
    end

    function initial_profile_v3(initial_temperature, z)
        if initial_temperature isa Real
            return fill(Float64(initial_temperature), length(z))
        elseif initial_temperature isa AbstractVector
            length(initial_temperature) == length(z) ||
                throw(DimensionMismatch("initial profile length must equal node count"))
            return Float64.(initial_temperature)
        end
        return Float64[initial_temperature(location) for location in z]
    end
end

function gas_profile_v3!(Tg, Qgas, hcell, Ts, time, p, operating, z, dx)
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
        nusselt = 10.0^p[5] * reynolds^p[6] * prandtl^(10.0^p[7])
        hcell[i] = nusselt * kf / Dh
        UA = hcell[i] * P_exchange * dx
        effectiveness = -expm1(-UA / (mdot * cp))
        Tg[i + 1] = Tg[i] + effectiveness * (Ts[i] - Tg[i])
        Qgas[i] = mdot * cp * (Tg[i + 1] - Tg[i])
    end
    return nothing
end

function receiver_rhs_v3!(dTs, Ts, time, p, operating, z, dx, solar_weights,
                          Tg, Qgas, hcell)
    nodes = length(Ts)
    ambient = operating.ambient_temperature(time)
    gas_profile_v3!(Tg, Qgas, hcell, Ts, time, p, operating, z, dx)
    fill!(dTs, 0.0)

    # Conservative axial conduction: each interface flux is removed from one
    # cell and added to its neighbor.
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

    Qfront = p[11] * A_frt * (Ts[1] - ambient) +
             p[12] * sigma_sb * A_frt * (Ts[1]^4 - ambient^4)
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

function receiver_ode_v3!(dTs, Ts, context, time)
    receiver_rhs_v3!(
        dTs, Ts, time, context.p, context.operating,
        context.z, context.dx, context.solar_weights,
        context.Tg, context.Qgas, context.hcell,
    )
    return nothing
end

function simulate_v3(p, operating, save_times;
                     initial_temperature=Tamb, nodes=default_nodes,
                     solver=Rodas5P(autodiff=AutoFiniteDiff()),
                     reltol=1e-6, abstol=1e-7, dtmax=30.0)
    length(p) == 12 || throw(ArgumentError("1D_v3 expects 12 parameters"))
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
    weights = solar_weights_v3(p[10], nodes)
    Ts_initial = initial_profile_v3(initial_temperature, z)
    Tg_history = Matrix{Float64}(undef, nodes + 1, length(time))
    h_history = Matrix{Float64}(undef, nodes, length(time))
    Tg = zeros(nodes + 1)
    Qgas = zeros(nodes)
    hcell = zeros(nodes)

    context = (
        p=p,
        operating=operating,
        z=z,
        dx=dx,
        solar_weights=weights,
        Tg=Tg,
        Qgas=Qgas,
        hcell=hcell,
    )
    problem = ODEProblem(receiver_ode_v3!, Ts_initial, (time[1], time[end]), context)
    solution = solve(
        problem,
        solver;
        saveat=time,
        save_everystep=false,
        dense=false,
        reltol=reltol,
        abstol=abstol,
        dtmax=dtmax,
    )
    successful_retcode(solution) ||
        error("1D_v3 ODE solve failed with retcode $(solution.retcode)")
    length(solution.t) == length(time) ||
        error("1D_v3 solver returned $(length(solution.t)) points; expected $(length(time))")

    Ts_history = reduce(hcat, solution.u)
    for output_index in eachindex(time)
        gas_profile_v3!(Tg, Qgas, hcell, view(Ts_history, :, output_index),
                        time[output_index], p, operating, z, dx)
        Tg_history[:, output_index] .= Tg
        h_history[:, output_index] .= hcell
    end

    return (
        time=time,
        z_solid=z,
        z_gas=zg,
        solid_temperature=Ts_history,
        gas_temperature=Tg_history,
        heat_transfer_coefficient=h_history,
        ode_solution=solution,
    )
end

function solid_at_v3(result, location)
    z = result.z_solid
    x = clamp(Float64(location), z[1], z[end])
    i = searchsortedlast(z, x)
    i >= length(z) && return copy(result.solid_temperature[end, :])
    i < 1 && return copy(result.solid_temperature[1, :])
    fraction = (x - z[i]) / (z[i + 1] - z[i])
    return @. result.solid_temperature[i, :] + fraction *
              (result.solid_temperature[i + 1, :] - result.solid_temperature[i, :])
end

function apply_T3_lag_v3(Tf_true, Tf_initial, time, tau)
    tau <= 1e-12 && return copy(Tf_true)
    Tf_model = similar(Tf_true)
    Tf_model[1] = Tf_initial
    for i in 2:length(time)
        gain = -expm1(-(time[i] - time[i - 1]) / tau)
        Tf_model[i] = Tf_model[i - 1] + gain * (Tf_true[i] - Tf_model[i - 1])
    end
    return Tf_model
end

function energy_rates_v3(Ts, time, p, operating)
    nodes = length(Ts)
    dx = L / nodes
    z = collect(range(dx / 2, step=dx, length=nodes))
    weights = solar_weights_v3(p[10], nodes)
    Tg = zeros(nodes + 1)
    Qgas_cell = zeros(nodes)
    hcell = zeros(nodes)
    dTs = zeros(nodes)
    receiver_rhs_v3!(dTs, Ts, time, p, operating, z, dx, weights,
                     Tg, Qgas_cell, hcell)

    ambient = operating.ambient_temperature(time)
    absorbed = p[9] * max(0.0, operating.irradiance(time)) * A_frt
    gas = sum(Qgas_cell)
    side = p[3] * dx * sum(Ts .- ambient)
    front = p[11] * A_frt * (Ts[1] - ambient) +
            p[12] * sigma_sb * A_frt * (Ts[1]^4 - ambient^4)
    rear = p[4] * (Ts[end] - ambient)
    volume = A_solid * dx
    storage = sum(rho_s * p[1] * Cps_f(Ts[i]) * volume * dTs[i] for i in 1:nodes)
    return (
        absorbed=absorbed, gas=gas, side=side, front=front, rear=rear,
        storage=storage,
        residual=absorbed - gas - side - front - rear - storage,
    )
end

# ============================================================================
# SECTION D: EXTERNAL EXPERIMENTAL DATA IMPORT
# ============================================================================
begin # external import
    println("[1D_v3] Importing experimental data...")
    include("import_exp_1D_v2.jl")
    sim_key_heat = copy(IDs)
    sim_key_cool = copy(IDs_cooling)
    println("[1D_v3] Loaded $(length(sim_key_heat)) heating and $(length(sim_key_cool)) cooling runs.")
end

# ============================================================================
# SECTION E: MODEL/DATA INTERFACE
# ============================================================================
begin # experimental profiles
    function measured_initial_profile_v3(data, simulation_id)
        z8 = sensor_positions[:T8]
        z9 = sensor_positions[:T9]
        z10 = sensor_positions[:T10]
        T8 = observation(data, simulation_id, "_T8")[1]
        T9 = observation(data, simulation_id, "_T9")[1]
        T10 = observation(data, simulation_id, "_T10")[1]
        return function (z)
            z <= z8 && return T8
            z >= z10 && return T10
            if z <= z9
                return T8 + (T9 - T8) * (z - z8) / (z9 - z8)
            end
            return T9 + (T10 - T9) * (z - z9) / (z10 - z9)
        end
    end

    function experimental_case_v3(simulation_id; is_cooling=false)
        data = is_cooling ? measurements_cooling : measurements
        conditions = is_cooling ? simulation_conditions_cooling : simulation_conditions
        time = observation_time(data, simulation_id)
        experiment = hcat(
            observation(data, simulation_id, "_T8"),
            observation(data, simulation_id, "_T9"),
            observation(data, simulation_id, "_T10"),
            observation(data, simulation_id, "_Tf"),
        )
        return conditions[simulation_id], time, experiment, data
    end

    function extract_outputs_v3(result, p, T3_initial)
        T8 = solid_at_v3(result, sensor_positions[:T8])
        T9 = solid_at_v3(result, sensor_positions[:T9])
        T10 = solid_at_v3(result, sensor_positions[:T10])
        Tf_true = vec(result.gas_temperature[end, :])
        Tf = apply_T3_lag_v3(Tf_true, T3_initial, result.time, p[8])
        h_mean = vec(mean(result.heat_transfer_coefficient, dims=1))
        return hcat(T8, T9, T10, Tf, h_mean)
    end

    function solve_case_v3(p, simulation_id; is_cooling=false,
                           nodes=default_nodes, solver=Rodas5P(autodiff=AutoFiniteDiff()),
                           reltol=1e-6, abstol=1e-7, dtmax=30.0)
        conditions, time, experiment, data = experimental_case_v3(
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
        result = simulate_v3(
            p, operating, time;
            initial_temperature=measured_initial_profile_v3(data, simulation_id),
            nodes=nodes,
            solver=solver,
            reltol=reltol,
            abstol=abstol,
            dtmax=dtmax,
        )
        outputs = extract_outputs_v3(result, p, experiment[1, 4])
        return outputs, result, experiment
    end

    # Legacy-style constant-condition interface.
    function NLmodeloptim(tvalues, rmp, pguess, initial_values;
                          nodes=default_nodes,
                          solver=Rodas5P(autodiff=AutoFiniteDiff()),
                          reltol=1e-6, abstol=1e-7, dtmax=30.0)
        operating = operating_condition_v3(
            irradiance=rmp[Io],
            flow_lpm=rmp[qlpm],
            inlet_temperature=get(rmp, T_in, Tamb),
            ambient_temperature=get(rmp, T_amb, Tamb),
        )
        T8_init, T9_init, T10_init, Tf_init = initial_values
        profile = z -> begin
            z8, z9, z10 = sensor_positions[:T8], sensor_positions[:T9], sensor_positions[:T10]
            z <= z8 && return T8_init
            z >= z10 && return T10_init
            z <= z9 && return T8_init + (T9_init - T8_init) * (z - z8) / (z9 - z8)
            return T9_init + (T10_init - T9_init) * (z - z9) / (z10 - z9)
        end
        result = simulate_v3(pguess, operating, tvalues;
                             initial_temperature=profile, nodes=nodes,
                             solver=solver, reltol=reltol, abstol=abstol,
                             dtmax=dtmax)
        return extract_outputs_v3(result, pguess, Tf_init)
    end

    function remakeAysha(pguess, cond_k, time_opt, initial_values;
                         nodes=default_nodes,
                         solver=Rodas5P(autodiff=AutoFiniteDiff()),
                         reltol=1e-6, abstol=1e-7, dtmax=30.0)
        return NLmodeloptim(time_opt, cond_k, pguess, initial_values;
                            nodes=nodes, solver=solver, reltol=reltol,
                            abstol=abstol, dtmax=dtmax)
    end
end

# ============================================================================
# SECTION F: COOLING-FIRST / HEATING-SECOND CALIBRATION
# ============================================================================
begin # objective functions
    function normalized_mse_v3(model, experiment)
        scale = max(maximum(experiment) - minimum(experiment), 20.0)
        return mean(abs2, (model .- experiment) ./ scale)
    end

    function loss_cases_v3(p, keys; is_cooling=false, nodes=15)
        total = 0.0
        for simulation_id in keys
            model, _, experiment = solve_case_v3(
                p, simulation_id;
                is_cooling=is_cooling,
                nodes=nodes,
                reltol=2e-5,
                abstol=2e-6,
                dtmax=60.0,
            )
            all(isfinite, model) || return Inf
            case_loss = mean(normalized_mse_v3(model[:, j], experiment[:, j]) for j in 1:4)
            total += case_loss
        end
        return total / length(keys)
    end

    function loss_cooling_v3(p_cool, keys=sim_key_cool; nodes=15)
        p = assemble_parameters_v3(p_cool, p_heat_init)
        regularization = 0.005 * (log(p_cool[1])^2 + log(p_cool[2])^2)
        return loss_cases_v3(p, keys; is_cooling=true, nodes=nodes) + regularization
    end

    function loss_heating_v3(p_heat, keys=sim_key_heat;
                             p_cool=p_cool_init, nodes=15)
        p = assemble_parameters_v3(p_cool, p_heat)
        return loss_cases_v3(p, keys; is_cooling=false, nodes=nodes)
    end
end

function optimize_with_nlopt_v3(objective, initial, lower, upper;
                                maximum_iterations=200,
                                maximum_time_seconds=120.0,
                                algorithm=NLopt.LN_NELDERMEAD(),
                                label="calibration")
    x0 = clamp.(Float64.(initial), lower, upper)
    lb = Float64.(lower)
    ub = Float64.(upper)
    evaluations = Ref(0)
    best_objective = Ref(Inf)

    function counted_objective(x, _)
        evaluations[] += 1
        value = objective(x)
        if isfinite(value) && value < best_objective[]
            best_objective[] = value
            println("[$label] evaluation $(evaluations[]): objective = $value")
            flush(stdout)
        end
        return value
    end

    optf = OptimizationFunction(counted_objective, SciMLBase.NoAD())
    optprob = Optimization.OptimizationProblem(optf, x0, nothing; lb=lb, ub=ub)
    solution = solve(optprob, algorithm;
                     maxiters=maximum_iterations,
                     maxtime=maximum_time_seconds)
    retcode_text = string(solution.retcode)
    converged = retcode_text in ("Success", "Terminated", "MaxIters", "MaxTime")
    return (
        minimizer=Vector{Float64}(solution.u),
        minimum=Float64(solution.objective),
        iterations=maximum_iterations,
        evaluations=evaluations[],
        converged=converged,
        retcode=solution.retcode,
        optimizer=algorithm,
        raw_solution=solution,
    )
end

nlopt_evaluation_budget_v3(iterations, parameter_count) =
    max(1, 1 + 2 * parameter_count * Int(iterations))

function calibrate_v3(; cooling_keys=sim_key_cool,
                      heating_keys=sim_key_heat,
                      nodes=15,
                      cooling_iterations=20,
                      heating_iterations=20,
                      cooling_time_seconds=120.0,
                      heating_time_seconds=120.0,
                      optimizer=NLopt.LN_NELDERMEAD())
    println("\n=== 1D_v3 STAGE 1: cooling calibration ===")
    cooling_fit = optimize_with_nlopt_v3(
        p -> loss_cooling_v3(p, cooling_keys; nodes=nodes),
        p_cool_init, lb_cool, ub_cool;
        maximum_iterations=nlopt_evaluation_budget_v3(cooling_iterations, length(p_cool_init)),
        maximum_time_seconds=cooling_time_seconds,
        algorithm=optimizer,
        label="cooling",
    )
    p_cool = cooling_fit.minimizer
    println("Cooling parameters: $p_cool")

    println("\n=== 1D_v3 STAGE 2: heating/optical calibration ===")
    heating_fit = optimize_with_nlopt_v3(
        p -> loss_heating_v3(p, heating_keys; p_cool=p_cool, nodes=nodes),
        p_heat_init, lb_heat, ub_heat;
        maximum_iterations=nlopt_evaluation_budget_v3(heating_iterations, length(p_heat_init)),
        maximum_time_seconds=heating_time_seconds,
        algorithm=optimizer,
        label="heating",
    )
    p_heat = heating_fit.minimizer
    parameters = assemble_parameters_v3(p_cool, p_heat)

    global pnew = parameters
    global calibration_v3 = (
        cooling=cooling_fit,
        heating=heating_fit,
        p_cool=p_cool,
        p_heat=p_heat,
        parameters=parameters,
    )
    println("\nFinal 1D_v3 parameters:")
    println(parameters)
    return calibration_v3
end

function quick_calibration_v3()
    return calibrate_v3(
        cooling_keys=["C69"],
        heating_keys=["E74"],
        nodes=11,
        cooling_iterations=2,
        heating_iterations=2,
    )
end

begin # optional fit driver
    if RUN_QUICK_CALIBRATION_V3 && RUN_CALIBRATION_V3
        error("Choose only one of RECEIVER1D_V3_RUN_QUICK_CALIBRATION and RECEIVER1D_V3_RUN_CALIBRATION")
    elseif RUN_QUICK_CALIBRATION_V3
        calibration_v3 = quick_calibration_v3()
    elseif RUN_CALIBRATION_V3
        calibration_v3 = calibrate_v3()
    end
end

# ============================================================================
# SECTION G: POST-PROCESSING
# ============================================================================
begin # result helpers
    function build_steady_results_v3(params=pnew; keys=sim_key_heat,
                                     nodes=default_nodes)
        results = NamedTuple[]
        for simulation_id in keys
            model, _, experiment = solve_case_v3(params, simulation_id; nodes=nodes)
            conditions = simulation_conditions[simulation_id]
            push!(results, (
                simulation_id=simulation_id,
                flow_lpm=conditions[qlpm],
                irradiance=conditions[Io],
                T8_model=model[end, 1],
                T8_experiment=experiment[end, 1],
                T9_model=model[end, 2],
                T9_experiment=experiment[end, 2],
                T10_model=model[end, 3],
                T10_experiment=experiment[end, 3],
                T3_model=model[end, 4],
                T3_experiment=experiment[end, 4],
                h_effective=model[end, 5],
            ))
        end
        return results
    end

    function transient_case_v3(simulation_id, p=pnew; is_cooling=false,
                               nodes=default_nodes)
        model, result, experiment = solve_case_v3(
            p, simulation_id; is_cooling=is_cooling, nodes=nodes
        )
        return (
            time=result.time,
            T8_model=model[:, 1], T8_experiment=experiment[:, 1],
            T9_model=model[:, 2], T9_experiment=experiment[:, 2],
            T10_model=model[:, 3], T10_experiment=experiment[:, 3],
            T3_model=model[:, 4], T3_experiment=experiment[:, 4],
            h_effective=model[:, 5],
            full_result=result,
        )
    end

    function plot_case_v3(simulation_id, params=pnew; is_cooling=false,
                          nodes=default_nodes)
        @eval using StatsPlots
        data = transient_case_v3(simulation_id, params;
                                 is_cooling=is_cooling, nodes=nodes)
        sensors = (:T8, :T9, :T10, :T3)
        colors = (:blue, :red, :green, :purple)
        plot_object = plot(title="1D_v3 transient: $simulation_id",
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

    function get_t90_v3(time, temperature)
        target = temperature[1] + 0.9 * (temperature[end] - temperature[1])
        index = temperature[end] >= temperature[1] ?
                findfirst(>=(target), temperature) : findfirst(<=(target), temperature)
        return isnothing(index) ? time[end] : time[index]
    end

    function compute_metrics_v3(p=pnew; heating_keys=sim_key_heat,
                                cooling_keys=sim_key_cool, nodes=default_nodes)
        metrics = NamedTuple[]
        for (keys, cooling) in ((heating_keys, false), (cooling_keys, true))
            for simulation_id in keys
                model, result, experiment = solve_case_v3(
                    p, simulation_id; is_cooling=cooling, nodes=nodes
                )
                for (sensor, column) in ((:T8, 1), (:T9, 2), (:T10, 3), (:T3, 4))
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

    function write_metrics_v3(path, metrics=compute_metrics_v3())
        open(path, "w") do io
            println(io, "simulation_id,phase,sensor,rmse_K,steady_error_K,t90_error_s")
            for row in metrics
                println(io, join((row.simulation_id, row.phase, row.sensor,
                                  row.rmse_K, row.steady_error_K,
                                  row.t90_error_s), ','))
            end
        end
        return path
    end
end

begin # optional top-level outputs
    if RUN_POSTPROCESS_V3
        steady_results_v3 = build_steady_results_v3(pnew)
        metrics_v3 = compute_metrics_v3(pnew)
        display(steady_results_v3)
        display(metrics_v3)
    end

    if RUN_PLOTS_V3
        display(plot_case_v3(first(sim_key_heat), pnew))
        display(plot_case_v3(first(sim_key_cool), pnew; is_cooling=true))
    end

    if EXPORT_METRICS_V3
        output_path = joinpath(@__DIR__, "summaries", "analysis_results_1D_v3.csv")
        write_metrics_v3(output_path, compute_metrics_v3(pnew))
        println("Saved 1D_v3 metrics to $output_path")
    end
end

println("""
[1D_v3] Model loaded.
Run quick_calibration_v3() for a short check.
Run calibration_v3 = calibrate_v3() for the full cooling/heating calibration.
The calibrated vector is stored in pnew.
""")
