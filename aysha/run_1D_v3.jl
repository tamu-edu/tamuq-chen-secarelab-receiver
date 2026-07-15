# ============================================================================
# run_1D_v3.jl - End-to-end 1D_v3 calibration and comparison workflow
# ============================================================================
# This runner mirrors the 0D scripts: load data/model, calibrate, build
# steady-state and transient comparison plots, and export summary metrics.
#
# Optional environment controls:
#   RECEIVER1D_V3_RUNNER_QUICK=true      use quick_calibration_v3()
#   RECEIVER1D_V3_RUNNER_NODES=15        calibration node count
#   RECEIVER1D_V3_RUNNER_PLOT_NODES=25   plotting/metrics node count
#   RECEIVER1D_V3_RUNNER_COOL_ITERS=20   cooling calibration iterations
#   RECEIVER1D_V3_RUNNER_HEAT_ITERS=20   heating calibration iterations
#   RECEIVER1D_V3_RUNNER_COOL_TIME=120   cooling max optimizer seconds
#   RECEIVER1D_V3_RUNNER_HEAT_TIME=120   heating max optimizer seconds
# ============================================================================

begin # libraries
    using Statistics
    using StatsPlots
end

include("1D_v3.jl")

begin # runner configuration
    runner_flag_v3(name, default=false) =
        lowercase(get(ENV, name, string(default))) == "true"

    runner_int_v3(name, default) = parse(Int, get(ENV, name, string(default)))
    runner_float_v3(name, default) = parse(Float64, get(ENV, name, string(default)))

    const RUNNER_QUICK_V3 = runner_flag_v3("RECEIVER1D_V3_RUNNER_QUICK", false)
    const RUNNER_NODES_V3 = runner_int_v3("RECEIVER1D_V3_RUNNER_NODES",
                                          RUNNER_QUICK_V3 ? 11 : 15)
    const RUNNER_PLOT_NODES_V3 = runner_int_v3("RECEIVER1D_V3_RUNNER_PLOT_NODES",
                                               default_nodes)
    const RUNNER_COOL_ITERS_V3 = runner_int_v3("RECEIVER1D_V3_RUNNER_COOL_ITERS",
                                               RUNNER_QUICK_V3 ? 2 : 20)
    const RUNNER_HEAT_ITERS_V3 = runner_int_v3("RECEIVER1D_V3_RUNNER_HEAT_ITERS",
                                               RUNNER_QUICK_V3 ? 2 : 20)
    const RUNNER_COOL_TIME_V3 = runner_float_v3("RECEIVER1D_V3_RUNNER_COOL_TIME", 120.0)
    const RUNNER_HEAT_TIME_V3 = runner_float_v3("RECEIVER1D_V3_RUNNER_HEAT_TIME", 120.0)

    const RUNNER_OUTPUT_DIR_V3 = joinpath(@__DIR__, "summaries", "1D_v3")
    const RUNNER_PLOT_DIR_V3 = joinpath(RUNNER_OUTPUT_DIR_V3, "plots")
    mkpath(RUNNER_OUTPUT_DIR_V3)
    mkpath(RUNNER_PLOT_DIR_V3)
end

begin # plot helpers
    const RUNNER_SENSOR_COLORS_V3 = Dict(
        :T8 => :blue,
        :T9 => :red,
        :T10 => :green,
        :T3 => :purple,
    )

    function save_runner_plot_v3(plot_object, filename)
        path = joinpath(RUNNER_PLOT_DIR_V3, filename)
        savefig(plot_object, path)
        println("[run_1D_v3] Saved plot: $path")
        return path
    end

    function steady_comparison_plot_v3(steady_results)
        sensors = (:T8, :T9, :T10, :T3)
        plot_object = plot(
            title="1D_v3 steady-state comparison",
            xlabel="Model temperature (K)",
            ylabel="Experiment temperature (K)",
            legend=:outerright,
            aspect_ratio=:equal,
        )

        all_values = Float64[]
        for sensor in sensors
            model_key = Symbol(sensor, "_model")
            experiment_key = Symbol(sensor, "_experiment")
            model_values = [getproperty(row, model_key) for row in steady_results]
            experiment_values = [getproperty(row, experiment_key) for row in steady_results]
            append!(all_values, model_values)
            append!(all_values, experiment_values)
            scatter!(plot_object, model_values, experiment_values;
                     label=String(sensor), ms=4, markerstrokewidth=0,
                     color=RUNNER_SENSOR_COLORS_V3[sensor])
        end

        lower = 50.0 * floor(minimum(all_values) / 50.0)
        upper = 50.0 * ceil(maximum(all_values) / 50.0)
        plot!(plot_object, [lower, upper], [lower, upper];
              label="1:1", color=:gray, linestyle=:dash)
        xlims!(plot_object, lower, upper)
        ylims!(plot_object, lower, upper)
        return plot_object
    end

    function steady_sensor_comparison_plot_v3(steady_results, sensor)
        model_key = Symbol(sensor, "_model")
        experiment_key = Symbol(sensor, "_experiment")
        model_values = [getproperty(row, model_key) for row in steady_results]
        experiment_values = [getproperty(row, experiment_key) for row in steady_results]
        all_values = vcat(model_values, experiment_values)
        lower = 50.0 * floor(minimum(all_values) / 50.0)
        upper = 50.0 * ceil(maximum(all_values) / 50.0)

        plot_object = scatter(
            model_values,
            experiment_values;
            title="1D_v3 steady-state comparison: $(sensor)",
            xlabel="$(sensor) model temperature (K)",
            ylabel="$(sensor) experiment temperature (K)",
            label=String(sensor),
            color=RUNNER_SENSOR_COLORS_V3[sensor],
            ms=5,
            markerstrokewidth=0,
            legend=:bottomright,
            aspect_ratio=:equal,
            xlims=(lower, upper),
            ylims=(lower, upper),
        )
        plot!(plot_object, [lower, upper], [lower, upper];
              label="1:1", color=:gray, linestyle=:dash)
        return plot_object
    end
end

begin # 0D-style metrics helpers
    tavg_v4_from_columns_v3(values) =
        0.248 .* values[:, 1] .+ 0.365 .* values[:, 2] .+ 0.387 .* values[:, 3]

    function leak_ratio_v3(Tgas, Tsolid, ambient; is_cooling=false)
        if is_cooling && abs(Tsolid - ambient) < 1.0
            return NaN
        end
        abs(Tsolid - ambient) < eps(Float64) && return NaN
        return (Tgas - ambient) / (Tsolid - ambient)
    end

    function effective_gas_metrics_v3(model, result, p, simulation_id; is_cooling=false)
        data = is_cooling ? measurements_cooling : measurements
        time = result.time
        Tin = observation(data, simulation_id, "_Tin")
        Tin_ss = Tin[end]
        Ts_avg_ss = tavg_v4_from_columns_v3(model)[end]
        Tf_true_ss = result.gas_temperature[end, end]
        h_eff = model[end, 5]

        denominator = Ts_avg_ss - Tin_ss
        epsilon = abs(denominator) < eps(Float64) ? NaN :
                  (Tf_true_ss - Tin_ss) / denominator

        mdot = m_dot(observation(data, simulation_id, "_flow")[end], Tin_ss)
        film = 0.5 * (Ts_avg_ss + Tin_ss)
        cp = cpf_f(film)
        ntu_from_h = mdot <= 0.0 ? NaN : h_eff * P_exchange * L / (mdot * cp)
        ntu_from_epsilon = isfinite(epsilon) && 0.0 <= epsilon < 1.0 ?
                           -log1p(-epsilon) : NaN
        ntu = isfinite(ntu_from_h) ? ntu_from_h : ntu_from_epsilon

        return (epsilon=epsilon, ntu=ntu, h_eff=h_eff)
    end

    function compute_0d_style_metrics_1D_v3(p=pnew; heating_keys=sim_key_heat,
                                            cooling_keys=sim_key_cool,
                                            nodes=RUNNER_PLOT_NODES_V3,
                                            run_id=1.0)
        rows = NamedTuple[]
        for (keys, cooling) in ((heating_keys, false), (cooling_keys, true))
            data = cooling ? measurements_cooling : measurements
            for simulation_id in keys
                model, result, experiment = solve_case_v3(
                    p, simulation_id; is_cooling=cooling, nodes=nodes
                )
                time = result.time
                Ts_mod = tavg_v4_from_columns_v3(model)
                Ts_exp = observation(data, simulation_id, "_Tavg_v4")
                Tf_mod = model[:, 4]
                Tf_exp = experiment[:, 4]

                Ts_ss_sim = Ts_mod[end]
                Ts_ss_exp = Ts_exp[end]
                Tf_ss_sim = Tf_mod[end]
                Tf_ss_exp = Tf_exp[end]

                dT_Ts = Ts_ss_sim - Ts_ss_exp
                dT_T3 = Tf_ss_sim - Tf_ss_exp
                R_leak_sim = leak_ratio_v3(Tf_ss_sim, Ts_ss_sim, Tamb;
                                           is_cooling=cooling)
                R_leak_exp = leak_ratio_v3(Tf_ss_exp, Ts_ss_exp, Tamb;
                                           is_cooling=cooling)

                t90_sim_Ts = get_t90_v3(time, Ts_mod)
                t90_exp_Ts = get_t90_v3(time, Ts_exp)
                t90_sim_T3 = get_t90_v3(time, Tf_mod)
                t90_exp_T3 = get_t90_v3(time, Tf_exp)

                gap_sim = Tf_ss_sim - Ts_ss_sim
                gap_exp = Tf_ss_exp - Ts_ss_exp
                gas_metrics = effective_gas_metrics_v3(
                    model, result, p, simulation_id; is_cooling=cooling
                )

                push!(rows, (
                    RunID=run_id,
                    Case=simulation_id,
                    Type=cooling ? "cooling" : "heating",
                    Ts_avg_SS_sim=Ts_ss_sim,
                    Ts_avg_SS_exp=Ts_ss_exp,
                    dT_Ts_avg=dT_Ts,
                    T3_SS_sim=Tf_ss_sim,
                    T3_SS_exp=Tf_ss_exp,
                    dT_T03=dT_T3,
                    R_leak_sim=R_leak_sim,
                    R_leak_exp=R_leak_exp,
                    t90_sim_Ts_s=t90_sim_Ts,
                    t90_exp_Ts_s=t90_exp_Ts,
                    dt90_Ts_s=t90_sim_Ts - t90_exp_Ts,
                    t90_sim_T3_s=t90_sim_T3,
                    t90_exp_T3_s=t90_exp_T3,
                    dt90_T3_s=t90_sim_T3 - t90_exp_T3,
                    Gap_ss_sim=gap_sim,
                    Gap_ss_exp=gap_exp,
                    dGap_ss=gap_sim - gap_exp,
                    epsilon_eff=gas_metrics.epsilon,
                    NTU_eff=gas_metrics.ntu,
                    h_eff_sim=gas_metrics.h_eff,
                ))
            end
        end
        return rows
    end

    function next_run_id_from_metrics_v3(path)
        isfile(path) || return 1.0
        values = Float64[]
        open(path, "r") do io
            first_line = true
            for line in eachline(io)
                if first_line
                    first_line = false
                    continue
                end
                isempty(strip(line)) && continue
                field = first(split(line, ','))
                try
                    push!(values, parse(Float64, field))
                catch
                end
            end
        end
        return isempty(values) ? 1.0 : maximum(values) + 1.0
    end

    function write_0d_style_metrics_1D_v3(path, rows; append=true)
        headers = (
            :RunID, :Case, :Type,
            :Ts_avg_SS_sim, :Ts_avg_SS_exp, :dT_Ts_avg,
            :T3_SS_sim, :T3_SS_exp, :dT_T03,
            :R_leak_sim, :R_leak_exp,
            :t90_sim_Ts_s, :t90_exp_Ts_s, :dt90_Ts_s,
            :t90_sim_T3_s, :t90_exp_T3_s, :dt90_T3_s,
            :Gap_ss_sim, :Gap_ss_exp, :dGap_ss,
            :epsilon_eff, :NTU_eff, :h_eff_sim,
        )
        mode = append && isfile(path) ? "a" : "w"
        open(path, mode) do io
            mode == "w" && println(io, join(string.(headers), ','))
            for row in rows
                println(io, join((getproperty(row, header) for header in headers), ','))
            end
        end
        return path
    end
end

begin # data and calibration
    println("[run_1D_v3] Loaded $(length(sim_key_heat)) heating runs and " *
            "$(length(sim_key_cool)) cooling runs.")

    if RUNNER_QUICK_V3
        println("[run_1D_v3] Running quick 1D_v3 calibration.")
        calibration_v3 = quick_calibration_v3()
    else
        println("[run_1D_v3] Running full 1D_v3 calibration.")
        calibration_v3 = calibrate_v3(
            nodes=RUNNER_NODES_V3,
            cooling_iterations=RUNNER_COOL_ITERS_V3,
            heating_iterations=RUNNER_HEAT_ITERS_V3,
            cooling_time_seconds=RUNNER_COOL_TIME_V3,
            heating_time_seconds=RUNNER_HEAT_TIME_V3,
        )
    end

    calibrated_parameters_v3 = calibration_v3.parameters
    parameter_path_v3 = joinpath(RUNNER_OUTPUT_DIR_V3, "parameters_1D_v3.csv")
    open(parameter_path_v3, "w") do io
        println(io, "index,name,value")
        names = (
            "gamma_C", "k_scale", "U_side", "U_rear",
            "A_Nu", "B_Re", "C_Pr", "tau_T3",
            "eta_abs", "beta_opt", "h_front", "eps_front",
        )
        for i in eachindex(calibrated_parameters_v3)
            println(io, join((i, names[i], calibrated_parameters_v3[i]), ','))
        end
    end
    println("[run_1D_v3] Saved calibrated parameters: $parameter_path_v3")
end

begin # post-processing and plots
    steady_results_v3 = build_steady_results_v3(pnew; nodes=RUNNER_PLOT_NODES_V3)
    metrics_v3 = compute_metrics_v3(pnew; nodes=RUNNER_PLOT_NODES_V3)

    metrics_path_v3 = joinpath(RUNNER_OUTPUT_DIR_V3, "analysis_results_1D_v3.csv")
    write_metrics_v3(metrics_path_v3, metrics_v3)
    println("[run_1D_v3] Saved metrics: $metrics_path_v3")

    metrics_0d_style_path_v3 = joinpath(
        RUNNER_OUTPUT_DIR_V3, "analysis_results_1D_v3_0D_style.csv"
    )
    next_run_id_v3 = next_run_id_from_metrics_v3(metrics_0d_style_path_v3)
    metrics_0d_style_v3 = compute_0d_style_metrics_1D_v3(
        pnew; nodes=RUNNER_PLOT_NODES_V3, run_id=next_run_id_v3
    )
    write_0d_style_metrics_1D_v3(metrics_0d_style_path_v3, metrics_0d_style_v3)
    println("[run_1D_v3] Saved 0D-style metrics: $metrics_0d_style_path_v3")

    steady_plot_v3 = steady_comparison_plot_v3(steady_results_v3)
    save_runner_plot_v3(steady_plot_v3, "steady_comparison_1D_v3.png")
    display(steady_plot_v3)

    for sensor in (:T8, :T9, :T10, :T3)
        sensor_plot = steady_sensor_comparison_plot_v3(steady_results_v3, sensor)
        save_runner_plot_v3(sensor_plot, "steady_comparison_$(sensor)_1D_v3.png")
        display(sensor_plot)
    end

    transient_plot_paths_v3 = String[]
    for simulation_id in sim_key_heat
        plot_object = plot_case_v3(simulation_id, pnew; nodes=RUNNER_PLOT_NODES_V3)
        push!(transient_plot_paths_v3,
              save_runner_plot_v3(plot_object, "transient_$(simulation_id)_1D_v3.png"))
        display(plot_object)
    end

    for simulation_id in sim_key_cool
        plot_object = plot_case_v3(simulation_id, pnew;
                                   is_cooling=true, nodes=RUNNER_PLOT_NODES_V3)
        push!(transient_plot_paths_v3,
              save_runner_plot_v3(plot_object, "transient_$(simulation_id)_1D_v3.png"))
        display(plot_object)
    end

    println("[run_1D_v3] Complete.")
    println("[run_1D_v3] Output directory: $RUNNER_OUTPUT_DIR_V3")
end
