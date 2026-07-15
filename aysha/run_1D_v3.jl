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
                     label=String(sensor), ms=4, markerstrokewidth=0)
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
            "gamma_C", "k_scale", "U_side", "U_rear", "h_ref", "n_exp",
            "tau_T3", "eta_abs", "beta_opt", "h_front", "eps_front",
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

    steady_plot_v3 = steady_comparison_plot_v3(steady_results_v3)
    save_runner_plot_v3(steady_plot_v3, "steady_comparison_1D_v3.png")
    display(steady_plot_v3)

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
