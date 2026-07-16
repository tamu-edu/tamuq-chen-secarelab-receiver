# ============================================================================
# run_1D_v7.jl - End-to-end 1D_v7 T2-boundary workflow
# ============================================================================
# Stages:
#   1. Heating fits gas heat-transfer shape and irradiance-level factors.
#   2. Cooling fits thermophysical / T2-boundary parameters.
#   3. Heating refits gas heat-transfer shape and irradiance-level factors.
#
# Optional environment controls:
#   RECEIVER1D_V7_RUNNER_QUICK=true
#   RECEIVER1D_V7_RUNNER_NODES=15
#   RECEIVER1D_V7_RUNNER_PLOT_NODES=25
#   RECEIVER1D_V7_RUNNER_HEAT_BASE_ITERS=20
#   RECEIVER1D_V7_RUNNER_COOL_ITERS=20
#   RECEIVER1D_V7_RUNNER_HEAT_REFIT_ITERS=20
# ============================================================================

begin # libraries
    using Statistics
    using StatsPlots
end

include("1D_v7.jl")

begin # runner configuration
    runner_flag_v7(name, default=false) =
        lowercase(get(ENV, name, string(default))) == "true"

    runner_int_v7(name, default) = parse(Int, get(ENV, name, string(default)))
    runner_float_v7(name, default) = parse(Float64, get(ENV, name, string(default)))

    const RUNNER_QUICK_V7 = runner_flag_v7("RECEIVER1D_V7_RUNNER_QUICK", false)
    const RUNNER_NODES_V7 = runner_int_v7("RECEIVER1D_V7_RUNNER_NODES",
                                          RUNNER_QUICK_V7 ? 11 : 15)
    const RUNNER_PLOT_NODES_V7 = runner_int_v7("RECEIVER1D_V7_RUNNER_PLOT_NODES",
                                               default_nodes)
    const RUNNER_HEAT_BASE_ITERS_V7 = runner_int_v7(
        "RECEIVER1D_V7_RUNNER_HEAT_BASE_ITERS", RUNNER_QUICK_V7 ? 2 : 20
    )
    const RUNNER_COOL_ITERS_V7 = runner_int_v7("RECEIVER1D_V7_RUNNER_COOL_ITERS",
                                               RUNNER_QUICK_V7 ? 2 : 20)
    const RUNNER_HEAT_REFIT_ITERS_V7 = runner_int_v7(
        "RECEIVER1D_V7_RUNNER_HEAT_REFIT_ITERS", RUNNER_QUICK_V7 ? 2 : 20
    )
    const RUNNER_COOL_TIME_V7 = runner_float_v7("RECEIVER1D_V7_RUNNER_COOL_TIME", 120.0)
    const RUNNER_HEAT_BASE_TIME_V7 = runner_float_v7(
        "RECEIVER1D_V7_RUNNER_HEAT_BASE_TIME", 120.0
    )
    const RUNNER_HEAT_REFIT_TIME_V7 = runner_float_v7(
        "RECEIVER1D_V7_RUNNER_HEAT_REFIT_TIME", 120.0
    )

    const RUNNER_OUTPUT_DIR_V7 = joinpath(@__DIR__, "summaries", "1D_v7")
    const RUNNER_PLOT_DIR_V7 = joinpath(RUNNER_OUTPUT_DIR_V7, "plots")
    mkpath(RUNNER_OUTPUT_DIR_V7)
    mkpath(RUNNER_PLOT_DIR_V7)
end

begin # runner helpers
    const RUNNER_SENSOR_COLORS_V7 = Dict(
        :T8 => :blue,
        :T9_pair => :red,
        :T10_pair => :green,
        :T3 => :purple,
    )

    function save_runner_plot_v7(plot_object, filename)
        path = joinpath(RUNNER_PLOT_DIR_V7, filename)
        savefig(plot_object, path)
        println("[run_1D_v7] Saved plot: $path")
        return path
    end

    function write_parameters_v7(path, parameters)
        names = (
            "gamma_C", "k_scale", "k_ins_scale",
            "A_Nu", "h_floor", "L_h",
            "tau_T3",
            "f_I_high", "f_I_mid", "f_I_low",
        )
        open(path, "w") do io
            println(io, "index,name,value")
            for i in eachindex(parameters)
                println(io, join((i, names[i], parameters[i]), ','))
            end
            println(io, "# fixed,B_Re,$B_RE_FIXED_V7")
            println(io, "# fixed,C_Pr,$C_PR_FIXED_V7")
            println(io, "# fixed,eta_abs,$ETA_ABS_FIXED_V7")
            println(io, "# fixed,beta_opt,$BETA_OPT_FIXED_V7")
            println(io, "# fixed,front_dep,$FRONT_DEPOSITION_FIXED_V7")
            println(io, "# fixed,h_front,$H_FRONT_FIXED_V7")
            println(io, "# fixed,eps_front,$EPS_FRONT_FIXED_V7")
            println(io, "# fixed,T2_boundary_radius_m,$T2_BOUNDARY_RADIUS_V7")
            println(io, "# fixed,G_receiver_to_T2_W_per_m_K,$RECEIVER_TO_T2_CONDUCTANCE_PER_LENGTH_V7")
            println(io, "# fixed,G_rear_to_T2_W_per_K,$RECEIVER_REAR_TO_T2_CONDUCTANCE_V7")
        end
        return path
    end

    function write_metrics_v7(path, metrics=compute_metrics_v7())
        open(path, "w") do io
            println(io, "simulation_id,phase,sensor,rmse_K,steady_error_K,t90_error_s,shape_loss")
            for row in metrics
                println(io, join((row.simulation_id, row.phase, row.sensor,
                                  row.rmse_K, row.steady_error_K,
                                  row.t90_error_s, row.shape_loss), ','))
            end
        end
        return path
    end

    function build_steady_results_v7(params=pnew_v7; keys=sim_key_heat,
                                     nodes=RUNNER_PLOT_NODES_V7)
        results = NamedTuple[]
        for simulation_id in keys
            model, result, experiment = solve_case_v7(params, simulation_id; nodes=nodes)
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
                T2_boundary=result.boundary_temperature[end],
                h_effective=model[end, 5],
            ))
        end
        return results
    end

    function steady_comparison_plot_v7(steady_results)
        sensors = (:T8, :T9_pair, :T10_pair, :T3)
        plot_object = plot(
            title="1D_v7 steady-state comparison",
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
                     label=String(sensor), color=RUNNER_SENSOR_COLORS_V7[sensor],
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
    println("[run_1D_v7] Loaded $(length(sim_key_heat)) heating runs and " *
            "$(length(sim_key_cool)) cooling runs.")

    if RUNNER_QUICK_V7
        println("[run_1D_v7] Running quick 1D_v7 calibration.")
        calibration_v7 = quick_calibration_v7()
    else
        println("[run_1D_v7] Running full 1D_v7 calibration.")
        calibration_v7 = calibrate_v7(
            nodes=RUNNER_NODES_V7,
            heating_base_iterations=RUNNER_HEAT_BASE_ITERS_V7,
            cooling_iterations=RUNNER_COOL_ITERS_V7,
            heating_refit_iterations=RUNNER_HEAT_REFIT_ITERS_V7,
            heating_base_time_seconds=RUNNER_HEAT_BASE_TIME_V7,
            cooling_time_seconds=RUNNER_COOL_TIME_V7,
            heating_refit_time_seconds=RUNNER_HEAT_REFIT_TIME_V7,
        )
    end

    parameter_path_v7 = joinpath(RUNNER_OUTPUT_DIR_V7, "parameters_1D_v7.csv")
    write_parameters_v7(parameter_path_v7, calibration_v7.parameters)
    println("[run_1D_v7] Saved calibrated parameters: $parameter_path_v7")
end

begin # post-processing and plots
    metrics_v7 = compute_metrics_v7(pnew_v7; nodes=RUNNER_PLOT_NODES_V7)
    metrics_path_v7 = joinpath(RUNNER_OUTPUT_DIR_V7, "analysis_results_1D_v7.csv")
    write_metrics_v7(metrics_path_v7, metrics_v7)
    println("[run_1D_v7] Saved metrics: $metrics_path_v7")

    steady_results_v7 = build_steady_results_v7(pnew_v7; nodes=RUNNER_PLOT_NODES_V7)
    steady_plot_v7 = steady_comparison_plot_v7(steady_results_v7)
    save_runner_plot_v7(steady_plot_v7, "steady_comparison_1D_v7.png")
    display(steady_plot_v7)

    for simulation_id in sim_key_heat
        plot_object = plot_case_v7(simulation_id, pnew_v7; nodes=RUNNER_PLOT_NODES_V7)
        save_runner_plot_v7(plot_object, "transient_$(simulation_id)_1D_v7.png")
        display(plot_object)
    end

    for simulation_id in sim_key_cool
        plot_object = plot_case_v7(simulation_id, pnew_v7;
                                   is_cooling=true, nodes=RUNNER_PLOT_NODES_V7)
        save_runner_plot_v7(plot_object, "transient_$(simulation_id)_1D_v7.png")
        display(plot_object)
    end

    println("[run_1D_v7] Complete.")
    println("[run_1D_v7] Output directory: $RUNNER_OUTPUT_DIR_V7")
end
