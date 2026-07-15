# ============================================================================
# 1D_v2.jl - Axial finite-volume solid/gas solar receiver model
# ============================================================================
# Refactor of 1D_v1 into the research-script flow used by 1D_v1.exp.All.jl and
# 0D_v4.jl. Numerical conservation remains in 1D_v1.jl; this file provides:
#   A. libraries and execution flags
#   B. geometry and parameter ordering
#   C. model construction
#   D. external experimental data import
#   E. solver/remake interface
#   F. cooling-first and heating-second fitting
#   G. post-processing and metrics
#
# Default execution only loads data and defines functions. Set environment
# variables RECEIVER1D_RUN_OPTIMIZATION, RECEIVER1D_RUN_POSTPROCESS,
# RECEIVER1D_RUN_PLOTS, or RECEIVER1D_EXPORT_METRICS to "true" to run those
# expensive/side-effecting sections.
# ============================================================================

begin # libraries
    using Statistics

    if !isdefined(@__MODULE__, :Receiver1D)
        include("1D_v1.jl")
    end
    using .Receiver1D
end

begin # execution controls
    const RUN_OPTIMIZATION = lowercase(get(ENV, "RECEIVER1D_RUN_OPTIMIZATION", "false")) == "true"
    const RUN_POSTPROCESS = lowercase(get(ENV, "RECEIVER1D_RUN_POSTPROCESS", "false")) == "true"
    const RUN_PLOTS = lowercase(get(ENV, "RECEIVER1D_RUN_PLOTS", "false")) == "true"
    const EXPORT_METRICS = lowercase(get(ENV, "RECEIVER1D_EXPORT_METRICS", "false")) == "true"
end

# ============================================================================
# SECTION A: GEOMETRY AND FIXED PARAMETERS
# ============================================================================
begin # fixed parameters
    Tamb = 22.448 + 273.15
    w_t = 19.0e-3
    A_frt = w_t^2
    L = 137e-3
    w_chnl = 1.5e-3
    t_chnl = 0.4e-3
    n_chnl = 10 * 10
    A_chnl_frt_all = n_chnl * w_chnl^2
    A_frt_solid = A_frt - A_chnl_frt_all
    Dh = w_chnl
    Pi_exchange = 4.0 * w_chnl * n_chnl
    Vs = A_frt_solid * L
    sigma_sb = 5.670374419e-8
    rho_s = 3200.0
    N_axial = 25

    geometry_v2 = Geometry(
        length=L,
        receiver_width=w_t,
        channel_width=w_chnl,
        wall_thickness=t_chnl,
        channel_count=n_chnl,
    )

    # Dictionary keys intentionally mirror the legacy import/model interface.
    const Io = :Io
    const qlpm = :qlpm
    const Tinit = :Tinit
    const T_in = :T_in
    const T_amb = :T_amb
end;

# ============================================================================
# SECTION B: FITTED PARAMETER ORDER
# ============================================================================
# Cooling parameters p[1:7]
#   p[1] = k_scale       effective SiC axial-conductivity multiplier
#   p[2] = cp_scale      effective solid heat-capacity multiplier
#   p[3] = U_side_prime  side-loss conductance per length (W/m/K)
#   p[4] = U_rear        rear/adaptor loss conductance (W/K)
#   p[5] = C_Nu          effective Nusselt correlation coefficient
#   p[6] = m_Re          Reynolds-number exponent
#   p[7] = tau_T3        outlet thermocouple/downstream lag (s)
# Heating parameters p[8:11]
#   p[8] = eta_abs       absorbed solar fraction
#   p[9] = beta_opt      optical extinction coefficient (1/m)
#   p[10] = h_front      front convection coefficient (W/m2/K)
#   p[11] = C_h_refine   multiplicative heating refinement of C_Nu
# ============================================================================
begin # parameter definitions
    p_cool_init = [1.0, 1.0, 0.35, 0.10, 6.1e-4, 1.3835, 20.0]
    lb_cool = [0.1, 0.3, 0.02, 0.005, 1e-5, 0.0, 0.5]
    ub_cool = [5.0, 3.0, 5.0, 5.0, 0.1, 3.0, 500.0]

    p_heat_init = [0.85, 50.0, 10.0, 1.0]
    lb_heat = [0.10, 1.0, 0.1, 0.25]
    ub_heat = [0.99, 2000.0, 100.0, 4.0]

    assemble_parameters(p_cool, p_heat) = vcat(Float64.(p_cool), Float64.(p_heat))
end

# ============================================================================
# SECTION C: 1D MODEL CONSTRUCTION
# ============================================================================
function build_model_v2(p; nodes=N_axial)
    length(p) == 11 || throw(ArgumentError("1D_v2 expects 11 parameters; received $(length(p))"))
    all(isfinite, p) || throw(ArgumentError("all model parameters must be finite"))

    solid = SolidProperties(
        density=rho_s,
        conductivity_scale=p[1],
        heat_capacity_scale=p[2],
    )
    htc = HeatTransferParameters(
        coefficient=p[5] * p[11],
        reynolds_exponent=p[6],
        prandtl_exponent=0.0,
        temperature_exponent=-0.4883,
        reference_temperature=600.0,
        development_coefficient=0.0,
        development_exponent=0.45,
        minimum_nusselt=1e-5,
    )
    losses = LossParameters(
        front_convection=p[10],
        front_emissivity=0.85,
        side_conductance_per_length=p[3],
        rear_conductance=p[4],
        rear_emissivity=0.0,
    )
    optics = OpticalParameters(
        absorbed_fraction=p[8],
        extinction_coefficient=p[9],
        front_deposition_fraction=0.0,
    )
    return ModelParameters(
        geometry=geometry_v2,
        solid=solid,
        heat_transfer=htc,
        losses=losses,
        optics=optics,
        nodes=nodes,
    )
end

# ============================================================================
# SECTION D: EXPERIMENTAL DATA IMPORT
# ============================================================================
begin # external import, analogous to import_exp_0D.jl
    include("import_exp_1D_v2.jl")
    sim_key_heat = copy(IDs)
    sim_key_cool = copy(IDs_cooling)
end

# ============================================================================
# SECTION E: MODEL SOLVER AND LEGACY-STYLE REMAKE INTERFACE
# ============================================================================
begin # solver helpers
    const SENSOR_POSITIONS_V2 = Dict("_T8" => 0.011, "_T9" => 0.058, "_T10" => 0.107)
    const MODEL_OUTPUT_COLUMNS = (:T8, :T9, :T10, :Tf, :h_mean)

    function initial_solid_profile(measurement_df, simulation_id)
        zdata = (0.011, 0.058, 0.107)
        tdata = (
            observation(measurement_df, simulation_id, "_T8")[1],
            observation(measurement_df, simulation_id, "_T9")[1],
            observation(measurement_df, simulation_id, "_T10")[1],
        )
        return function (z)
            z <= zdata[1] && return tdata[1]
            z >= zdata[3] && return tdata[3]
            i = z <= zdata[2] ? 1 : 2
            fraction = (z - zdata[i]) / (zdata[i + 1] - zdata[i])
            return tdata[i] + fraction * (tdata[i + 1] - tdata[i])
        end
    end

    function apply_T3_lag(target, measured_initial, time, tau)
        tau <= 1e-12 && return copy(target)
        filtered = similar(target)
        filtered[1] = measured_initial
        for i in 2:length(time)
            gain = -expm1(-(time[i] - time[i - 1]) / tau)
            filtered[i] = filtered[i - 1] + gain * (target[i] - filtered[i - 1])
        end
        return filtered
    end

    function extract_model_outputs(result, pguess; T3_initial=result.gas_temperature[end, 1])
        T8 = solid_at(result, SENSOR_POSITIONS_V2["_T8"])
        T9 = solid_at(result, SENSOR_POSITIONS_V2["_T9"])
        T10 = solid_at(result, SENSOR_POSITIONS_V2["_T10"])
        Tf_true = vec(result.gas_temperature[end, :])
        Tf = apply_T3_lag(Tf_true, T3_initial, result.time, pguess[7])
        h_mean = vec(mean(result.local_heat_transfer_coefficient, dims=1))
        return hcat(T8, T9, T10, Tf, h_mean)
    end

    function solve_model_v2(pguess, operating_condition, tvalues;
                            initial_temperature=operating_condition.ambient_temperature(first(tvalues)),
                            T3_initial=operating_condition.inlet_temperature(first(tvalues)),
                            nodes=N_axial, maximum_step=2.0)
        model = build_model_v2(pguess; nodes=nodes)
        result = simulate(model, operating_condition, tvalues;
                          initial_temperature=initial_temperature,
                          maximum_step=maximum_step)
        outputs = extract_model_outputs(result, pguess; T3_initial=T3_initial)
        return outputs, result
    end

    # Same constant-condition call shape used by 0D_v4.jl and 1D_v1.exp.All.jl.
    function solve_model(pguess, Io_val, qlpm_val, Tinit_val, tvalues;
                         Tin_val=Tamb, Tamb_val=Tamb, initial_profile=Tinit_val,
                         T3_initial=Tinit_val, tolr=1e-7, nodes=N_axial)
        operating_condition = OperatingCondition(
            irradiance=Io_val,
            flow_lpm=qlpm_val,
            inlet_temperature=Tin_val,
            ambient_temperature=Tamb_val,
        )
        outputs, _ = solve_model_v2(pguess, operating_condition, tvalues;
                                    initial_temperature=initial_profile,
                                    T3_initial=T3_initial, nodes=nodes)
        return outputs
    end

    function NLmodeloptim(tvalues, rmp, pguess, tolr=1e-7;
                          initial_profile=rmp[Tinit], T3_initial=rmp[Tinit],
                          nodes=N_axial)
        return solve_model(
            pguess, rmp[Io], rmp[qlpm], rmp[Tinit], tvalues;
            Tin_val=get(rmp, T_in, Tamb), Tamb_val=get(rmp, T_amb, Tamb),
            initial_profile=initial_profile, T3_initial=T3_initial,
            tolr=tolr, nodes=nodes,
        )
    end

    function remakeAysha(pguess, cond_k, time_opt, Tinit_val; tolr=1e-7,
                         initial_profile=Tinit_val, T3_initial=Tinit_val,
                         nodes=N_axial)
        rmp = copy(cond_k)
        rmp[Tinit] = Tinit_val
        return NLmodeloptim(time_opt, rmp, pguess, tolr;
                            initial_profile=initial_profile,
                            T3_initial=T3_initial, nodes=nodes)
    end

    function solve_case_v2(pguess, simulation_id; is_cooling=false,
                           nodes=N_axial, maximum_step=2.0)
        measurement_df = is_cooling ? measurements_cooling : measurements
        conditions = is_cooling ? simulation_conditions_cooling : simulation_conditions
        cond_k = conditions[simulation_id]
        time = observation_time(measurement_df, simulation_id)
        flow = observation(measurement_df, simulation_id, "_flow")
        inlet = observation(measurement_df, simulation_id, "_Tin")
        ambient = observation(measurement_df, simulation_id, "_Tamb")
        irradiance = fill(cond_k[Io], length(time))
        operating_condition = OperatingCondition(
            irradiance=linear_history(time, irradiance),
            flow_lpm=linear_history(time, flow),
            inlet_temperature=linear_history(time, inlet),
            ambient_temperature=linear_history(time, ambient),
        )
        profile = initial_solid_profile(measurement_df, simulation_id)
        T3_initial = observation(measurement_df, simulation_id, "_Tf")[1]
        return solve_model_v2(pguess, operating_condition, time;
                              initial_temperature=profile, T3_initial=T3_initial,
                              nodes=nodes, maximum_step=maximum_step)
    end
end

# ============================================================================
# SECTION F: TWO-STAGE OPTIMIZATION
# ============================================================================
begin # objective functions
    function normalized_mse(model_values, experimental_values)
        scale = max(maximum(experimental_values) - minimum(experimental_values), 20.0)
        return mean(abs2, (model_values .- experimental_values) ./ scale)
    end

    function loss_case_v2(p_full, simulation_id; is_cooling=false, nodes=15)
        measurement_df = is_cooling ? measurements_cooling : measurements
        outputs, _ = solve_case_v2(p_full, simulation_id;
                                   is_cooling=is_cooling, nodes=nodes,
                                   maximum_step=3.0)
        error = 0.0
        for (column, obs_id) in ((1, "_T8"), (2, "_T9"), (3, "_T10"))
            error += normalized_mse(outputs[:, column],
                                    observation(measurement_df, simulation_id, obs_id))
        end
        error += 0.5 * normalized_mse(outputs[:, 4],
                                      observation(measurement_df, simulation_id, "_Tf"))
        return error / 3.5
    end

    function loss_cooling(p_cool, keys=sim_key_cool; nodes=15)
        p_full = assemble_parameters(p_cool, p_heat_init)
        return mean(loss_case_v2(p_full, sm; is_cooling=true, nodes=nodes) for sm in keys)
    end

    function loss_heating(p_heat, keys=sim_key_heat; p_cool=p_cool_init, nodes=15)
        p_full = assemble_parameters(p_cool, p_heat)
        return mean(loss_case_v2(p_full, sm; is_cooling=false, nodes=nodes) for sm in keys)
    end

    function bounded_pattern_search(objective, initial, lower, upper;
                                    initial_step=0.20, maximum_iterations=30,
                                    tolerance=1e-3)
        x = clamp.(Float64.(initial), lower, upper)
        span = Float64.(upper .- lower)
        # Scale each search direction around its initial magnitude. Using only
        # the full bound span is ineffective for C_Nu, which is O(1e-4), while
        # beta_opt is O(10-100).
        step = initial_step .* max.(abs.(x), 0.02 .* span)
        fx = objective(x)
        evaluations = 1
        for iteration in 1:maximum_iterations
            improved = false
            for j in eachindex(x)
                for direction in (-1.0, 1.0)
                    trial = copy(x)
                    trial[j] = clamp(x[j] + direction * step[j], lower[j], upper[j])
                    ftrial = objective(trial)
                    evaluations += 1
                    if ftrial < fx
                        x, fx = trial, ftrial
                        improved = true
                    end
                end
            end
            if !improved
                step .*= 0.5
                maximum(step ./ max.(span, eps())) < tolerance &&
                    return (minimizer=x, minimum=fx, iterations=iteration,
                            evaluations=evaluations, converged=true)
            end
        end
        return (minimizer=x, minimum=fx, iterations=maximum_iterations,
                evaluations=evaluations, converged=false)
    end
end

begin # staged fit driver
    if RUN_OPTIMIZATION
        println("--- STAGE 1: 1D cooling optimization ---")
        cooling_fit_v2 = bounded_pattern_search(
            p -> loss_cooling(p, sim_key_cool),
            p_cool_init, lb_cool, ub_cool;
            maximum_iterations=30,
        )
        p_cool_opt = cooling_fit_v2.minimizer
        println("Cooling objective: ", cooling_fit_v2.minimum)

        println("--- STAGE 2: 1D heating optimization ---")
        heating_fit_v2 = bounded_pattern_search(
            p -> loss_heating(p, sim_key_heat; p_cool=p_cool_opt),
            p_heat_init, lb_heat, ub_heat;
            maximum_iterations=30,
        )
        p_heat_opt = heating_fit_v2.minimizer
        println("Heating objective: ", heating_fit_v2.minimum)
    else
        p_cool_opt = copy(p_cool_init)
        p_heat_opt = copy(p_heat_init)
    end
    pnew = assemble_parameters(p_cool_opt, p_heat_opt)
end

# ============================================================================
# SECTION G: POST-PROCESSING - STEADY STATE, TRANSIENTS, AND METRICS
# ============================================================================
begin # post-processing functions
    function build_steady_results_v2(params=pnew; keys=sim_key_heat, nodes=N_axial)
        results = NamedTuple[]
        for sm in keys
            outputs, _ = solve_case_v2(params, sm; nodes=nodes)
            cond = simulation_conditions[sm]
            push!(results, (
                simulation_id=sm, flow_lpm=cond[qlpm], irradiance=cond[Io],
                T8_model=outputs[end, 1], T8_experiment=observation(measurements, sm, "_T8")[end],
                T9_model=outputs[end, 2], T9_experiment=observation(measurements, sm, "_T9")[end],
                T10_model=outputs[end, 3], T10_experiment=observation(measurements, sm, "_T10")[end],
                T3_model=outputs[end, 4], T3_experiment=observation(measurements, sm, "_Tf")[end],
                h_effective=outputs[end, 5],
            ))
        end
        return results
    end

    function transient_case_v2(simulation_id, params=pnew; is_cooling=false,
                               nodes=N_axial)
        measurement_df = is_cooling ? measurements_cooling : measurements
        outputs, result = solve_case_v2(params, simulation_id;
                                        is_cooling=is_cooling, nodes=nodes)
        return (
            time=copy(result.time),
            T8_model=outputs[:, 1], T8_experiment=observation(measurement_df, simulation_id, "_T8"),
            T9_model=outputs[:, 2], T9_experiment=observation(measurement_df, simulation_id, "_T9"),
            T10_model=outputs[:, 3], T10_experiment=observation(measurement_df, simulation_id, "_T10"),
            T3_model=outputs[:, 4], T3_experiment=observation(measurement_df, simulation_id, "_Tf"),
            h_effective=outputs[:, 5]
        )
    end

    function plot_case_v2(simulation_id, params=pnew; is_cooling=false,
                          nodes=N_axial)
        @eval using StatsPlots
        data = transient_case_v2(simulation_id, params;
                                 is_cooling=is_cooling, nodes=nodes)
        plot_object = plot(title="1D_v2 transient: $simulation_id",
                           xlabel="Time (s)", ylabel="Temperature (K)",
                           legend=:outerright)
        for sensor in (:T8, :T9, :T10, :T3)
            plot!(plot_object, data.time, getproperty(data, Symbol(sensor, "_model")),
                  label="$(sensor) model", lw=2)
            scatter!(plot_object, data.time, getproperty(data, Symbol(sensor, "_experiment")),
                     label="$(sensor) experiment", ms=2, markerstrokewidth=0)
        end
        return plot_object
    end

    function get_t90(time, temperature)
        target = temperature[1] + 0.9 * (temperature[end] - temperature[1])
        index = temperature[end] >= temperature[1] ?
                findfirst(>=(target), temperature) : findfirst(<=(target), temperature)
        return isnothing(index) ? time[end] : time[index]
    end

    function compute_metrics_v2(params=pnew; heating_keys=sim_key_heat,
                                cooling_keys=sim_key_cool, nodes=N_axial)
        metrics = NamedTuple[]
        for (keys, cooling) in ((heating_keys, false), (cooling_keys, true))
            measurement_df = cooling ? measurements_cooling : measurements
            for sm in keys
                outputs, result = solve_case_v2(params, sm;
                                                is_cooling=cooling, nodes=nodes)
                for (sensor, column, obs_id) in (
                    (:T8, 1, "_T8"), (:T9, 2, "_T9"),
                    (:T10, 3, "_T10"), (:T3, 4, "_Tf"),
                )
                    experiment = observation(measurement_df, sm, obs_id)
                    model = outputs[:, column]
                    push!(metrics, (
                        simulation_id=sm,
                        phase=cooling ? :cooling : :heating,
                        sensor=sensor,
                        rmse_K=sqrt(mean(abs2, model .- experiment)),
                        steady_error_K=model[end] - experiment[end],
                        t90_error_s=get_t90(result.time, model) - get_t90(result.time, experiment),
                    ))
                end
            end
        end
        return metrics
    end
end

begin # optional top-level outputs
    if RUN_POSTPROCESS
        steady_results_v2 = build_steady_results_v2(pnew)
        metrics_v2 = compute_metrics_v2(pnew)
        display(steady_results_v2)
        display(metrics_v2)
    end

    if RUN_PLOTS
        display(plot_case_v2(first(sim_key_heat), pnew))
        display(plot_case_v2(first(sim_key_cool), pnew; is_cooling=true))
    end

    if EXPORT_METRICS
        output_path = joinpath(@__DIR__, "summaries", "analysis_results_1D_v2.csv")
        metrics_to_write = compute_metrics_v2(pnew)
        open(output_path, "w") do io
            println(io, "simulation_id,phase,sensor,rmse_K,steady_error_K,t90_error_s")
            for row in metrics_to_write
                println(io, join((row.simulation_id, row.phase, row.sensor,
                                  row.rmse_K, row.steady_error_K, row.t90_error_s), ','))
            end
        end
        println("Saved 1D_v2 metrics to $output_path")
    end
end
