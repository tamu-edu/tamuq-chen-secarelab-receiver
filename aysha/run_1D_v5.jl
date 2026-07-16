# ============================================================================
# run_1D_v5.jl - End-to-end reduced 1D_v5 workflow
# ============================================================================
# v5 uses a staged calibration:
#   1. Cooling cases fit thermophysical/gas-exchange parameters.
#   2. Heating cases fit optical/front-source parameters with cooling fixed.
#
# Optional environment controls:
#   RECEIVER1D_V5_RUNNER_QUICK=true      use quick_calibration_v5()
#   RECEIVER1D_V5_RUNNER_NODES=15        calibration node count
#   RECEIVER1D_V5_RUNNER_PLOT_NODES=25   plotting/metrics node count
#   RECEIVER1D_V5_RUNNER_HEAT_BASE_ITERS=20  first heating calibration iterations
#   RECEIVER1D_V5_RUNNER_COOL_ITERS=20       cooling calibration iterations
#   RECEIVER1D_V5_RUNNER_HEAT_REFIT_ITERS=20 final heating calibration iterations
# ============================================================================

begin # libraries
    using Statistics
    using StatsPlots
end

include("1D_v5.jl")

begin # runner configuration
    runner_flag_v5(name, default=false) =
        lowercase(get(ENV, name, string(default))) == "true"

    runner_int_v5(name, default) = parse(Int, get(ENV, name, string(default)))
    runner_float_v5(name, default) = parse(Float64, get(ENV, name, string(default)))

    const RUNNER_QUICK_V5 = runner_flag_v5("RECEIVER1D_V5_RUNNER_QUICK", false)
    const RUNNER_NODES_V5 = runner_int_v5("RECEIVER1D_V5_RUNNER_NODES",
                                          RUNNER_QUICK_V5 ? 11 : 15)
    const RUNNER_PLOT_NODES_V5 = runner_int_v5("RECEIVER1D_V5_RUNNER_PLOT_NODES",
                                               default_nodes)
    const RUNNER_HEAT_BASE_ITERS_V5 = runner_int_v5(
        "RECEIVER1D_V5_RUNNER_HEAT_BASE_ITERS", RUNNER_QUICK_V5 ? 2 : 20
    )
    const RUNNER_COOL_ITERS_V5 = runner_int_v5("RECEIVER1D_V5_RUNNER_COOL_ITERS",
                                               RUNNER_QUICK_V5 ? 2 : 20)
    const RUNNER_HEAT_REFIT_ITERS_V5 = runner_int_v5(
        "RECEIVER1D_V5_RUNNER_HEAT_REFIT_ITERS", RUNNER_QUICK_V5 ? 2 : 20
    )
    const RUNNER_COOL_TIME_V5 = runner_float_v5("RECEIVER1D_V5_RUNNER_COOL_TIME", 120.0)
    const RUNNER_HEAT_BASE_TIME_V5 = runner_float_v5(
        "RECEIVER1D_V5_RUNNER_HEAT_BASE_TIME", 120.0
    )
    const RUNNER_HEAT_REFIT_TIME_V5 = runner_float_v5(
        "RECEIVER1D_V5_RUNNER_HEAT_REFIT_TIME", 120.0
    )

    const RUNNER_OUTPUT_DIR_V5 = joinpath(@__DIR__, "summaries", "1D_v5")
    const RUNNER_PLOT_DIR_V5 = joinpath(RUNNER_OUTPUT_DIR_V5, "plots")
    mkpath(RUNNER_OUTPUT_DIR_V5)
    mkpath(RUNNER_PLOT_DIR_V5)
end

begin # runner helpers
    const RUNNER_SENSOR_COLORS_V5 = Dict(
        :T8 => :blue,
        :T9_pair => :red,
        :T10_pair => :green,
        :T3 => :purple,
    )

    function save_runner_plot_v5(plot_object, filename)
        path = joinpath(RUNNER_PLOT_DIR_V5, filename)
        savefig(plot_object, path)
        println("[run_1D_v5] Saved plot: $path")
        return path
    end

    function write_parameters_v5(path, parameters)
        names = (
            "gamma_C", "k_scale", "U_side", "U_rear",
            "A_Nu", "h_floor", "L_h", "tau_T3",
            "eta_abs", "beta_opt", "front_dep", "h_front",
        )
        open(path, "w") do io
            println(io, "index,name,value")
            for i in eachindex(parameters)
                println(io, join((i, names[i], parameters[i]), ','))
            end
            println(io, "# fixed,B_Re,$B_RE_FIXED_V5")
            println(io, "# fixed,C_Pr,$C_PR_FIXED_V5")
            println(io, "# fixed,eps_front,$EPS_FRONT_FIXED_V5")
            println(io, "# fixed_in_calibration,eta_abs,$ETA_ABS_FIXED_V5")
            println(io, "# fixed_in_calibration,beta_opt,$BETA_OPT_FIXED_V5")
            println(io, "# fixed_in_calibration,front_dep,$FRONT_DEPOSITION_FIXED_V5")
        end
        return path
    end

    function write_metrics_v5(path, metrics=compute_metrics_v5())
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

    function build_steady_results_v5(params=pnew_v5; keys=sim_key_heat,
                                     nodes=RUNNER_PLOT_NODES_V5)
        results = NamedTuple[]
        for simulation_id in keys
            model, _, experiment = solve_case_v5(params, simulation_id; nodes=nodes)
            conditions = simulation_conditions[simulation_id]
            push!(results, (
                simulation_id=simulation_id,
                flow_lpm=conditions[qlpm],
                irradiance=conditions[Io],
                T8_model=model[end, 1],
                T8_experiment=experiment[end, 1],
                T9_pair_model=model[end, 2],
                T9_pair_experiment=experiment[end, 2],
                T10_pair_model=model[end, 3],
                T10_pair_experiment=experiment[end, 3],
                T3_model=model[end, 4],
                T3_experiment=experiment[end, 4],
                h_effective=model[end, 5],
            ))
        end
        return results
    end

    function steady_comparison_plot_v5(steady_results)
        sensors = (:T8, :T9_pair, :T10_pair, :T3)
        plot_object = plot(
            title="1D_v5 steady-state comparison",
            xlabel="Model temperature (K)",
            ylabel="Experiment temperature (K)",
            legend=:outerright,
            aspect_ratio=:equal,
        )

        all_values = Float64[]
        for sensor in sensors
            model_values = [getproperty(row, Symbol(sensor, "_model")) for row in steady_results]
            experiment_values = [getproperty(row, Symbol(sensor, "_experiment")) for row in steady_results]
            append!(all_values, model_values)
            append!(all_values, experiment_values)
            scatter!(plot_object, model_values, experiment_values;
                     label=String(sensor), color=RUNNER_SENSOR_COLORS_V5[sensor],
                     ms=4, markerstrokewidth=0)
        end

        lower = 50.0 * floor(minimum(all_values) / 50.0)
        upper = 50.0 * ceil(maximum(all_values) / 50.0)
        plot!(plot_object, [lower, upper], [lower, upper];
              label="1:1", color=:gray, linestyle=:dash)
        xlims!(plot_object, lower, upper)
        ylims!(plot_object, lower, upper)
        return plot_object
    end
end

begin # calibration
    println("[run_1D_v5] Loaded $(length(sim_key_heat)) heating runs and " *
            "$(length(sim_key_cool)) cooling runs.")

    if RUNNER_QUICK_V5
        println("[run_1D_v5] Running quick 1D_v5 calibration.")
        calibration_v5 = quick_calibration_v5()
    else
        println("[run_1D_v5] Running full 1D_v5 calibration.")
        calibration_v5 = calibrate_v5(
            nodes=RUNNER_NODES_V5,
            heating_base_iterations=RUNNER_HEAT_BASE_ITERS_V5,
            cooling_iterations=RUNNER_COOL_ITERS_V5,
            heating_refit_iterations=RUNNER_HEAT_REFIT_ITERS_V5,
            heating_base_time_seconds=RUNNER_HEAT_BASE_TIME_V5,
            cooling_time_seconds=RUNNER_COOL_TIME_V5,
            heating_refit_time_seconds=RUNNER_HEAT_REFIT_TIME_V5,
        )
    end

    parameter_path_v5 = joinpath(RUNNER_OUTPUT_DIR_V5, "parameters_1D_v5.csv")
    write_parameters_v5(parameter_path_v5, calibration_v5.parameters)
    println("[run_1D_v5] Saved calibrated parameters: $parameter_path_v5")
end

begin # post-processing and plots
    metrics_v5 = compute_metrics_v5(pnew_v5; nodes=RUNNER_PLOT_NODES_V5)
    metrics_path_v5 = joinpath(RUNNER_OUTPUT_DIR_V5, "analysis_results_1D_v5.csv")
    write_metrics_v5(metrics_path_v5, metrics_v5)
    println("[run_1D_v5] Saved metrics: $metrics_path_v5")

    steady_results_v5 = build_steady_results_v5(pnew_v5; nodes=RUNNER_PLOT_NODES_V5)
    steady_plot_v5 = steady_comparison_plot_v5(steady_results_v5)
    save_runner_plot_v5(steady_plot_v5, "steady_comparison_1D_v5.png")
    display(steady_plot_v5)

    for simulation_id in sim_key_heat
        plot_object = plot_case_v5(simulation_id, pnew_v5; nodes=RUNNER_PLOT_NODES_V5)
        save_runner_plot_v5(plot_object, "transient_$(simulation_id)_1D_v5.png")
        display(plot_object)
    end

    for simulation_id in sim_key_cool
        plot_object = plot_case_v5(simulation_id, pnew_v5;
                                   is_cooling=true, nodes=RUNNER_PLOT_NODES_V5)
        save_runner_plot_v5(plot_object, "transient_$(simulation_id)_1D_v5.png")
        display(plot_object)
    end

    println("[run_1D_v5] Complete.")
    println("[run_1D_v5] Output directory: $RUNNER_OUTPUT_DIR_V5")
end
