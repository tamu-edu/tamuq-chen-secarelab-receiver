# ============================================================================
# run_1D_v8.jl - End-to-end 1D_v8 rear-tube workflow
# ============================================================================
# Stages:
#   1. Heating fits gas heat-transfer shape and irradiance-level factors.
#   2. Cooling fits receiver thermophysical parameters with fixed rear extension.
#   3. Heating refits gas heat-transfer shape and irradiance-level factors.
#
# Optional environment controls:
#   RECEIVER1D_V8_RUNNER_QUICK=true
#   RECEIVER1D_V8_RUNNER_NODES=15
#   RECEIVER1D_V8_RUNNER_PLOT_NODES=25
#   RECEIVER1D_V8_RUNNER_HEAT_BASE_ITERS=20
#   RECEIVER1D_V8_RUNNER_COOL_ITERS=20
#   RECEIVER1D_V8_RUNNER_HEAT_REFIT_ITERS=20
# ============================================================================

begin # libraries
    using Statistics
    using StatsPlots
end

include("1D_v8.jl")

begin # runner configuration
    runner_flag_v8(name, default=false) =
        lowercase(get(ENV, name, string(default))) == "true"

    runner_int_v8(name, default) = parse(Int, get(ENV, name, string(default)))
    runner_float_v8(name, default) = parse(Float64, get(ENV, name, string(default)))

    const RUNNER_QUICK_v8 = runner_flag_v8("RECEIVER1D_V8_RUNNER_QUICK", false)
    const RUNNER_NODES_v8 = runner_int_v8("RECEIVER1D_V8_RUNNER_NODES",
                                          RUNNER_QUICK_v8 ? 11 : 15)
    const RUNNER_PLOT_NODES_v8 = runner_int_v8("RECEIVER1D_V8_RUNNER_PLOT_NODES",
                                               default_nodes)
    const RUNNER_HEAT_BASE_ITERS_v8 = runner_int_v8(
        "RECEIVER1D_V8_RUNNER_HEAT_BASE_ITERS", RUNNER_QUICK_v8 ? 2 : 20
    )
    const RUNNER_COOL_ITERS_v8 = runner_int_v8("RECEIVER1D_V8_RUNNER_COOL_ITERS",
                                               RUNNER_QUICK_v8 ? 2 : 20)
    const RUNNER_HEAT_REFIT_ITERS_v8 = runner_int_v8(
        "RECEIVER1D_V8_RUNNER_HEAT_REFIT_ITERS", RUNNER_QUICK_v8 ? 2 : 20
    )
    const RUNNER_COOL_TIME_v8 = runner_float_v8("RECEIVER1D_V8_RUNNER_COOL_TIME", 120.0)
    const RUNNER_HEAT_BASE_TIME_v8 = runner_float_v8(
        "RECEIVER1D_V8_RUNNER_HEAT_BASE_TIME", 120.0
    )
    const RUNNER_HEAT_REFIT_TIME_v8 = runner_float_v8(
        "RECEIVER1D_V8_RUNNER_HEAT_REFIT_TIME", 120.0
    )

    const RUNNER_OUTPUT_DIR_v8 = joinpath(@__DIR__, "summaries", "1D_v8")
    const RUNNER_PLOT_DIR_v8 = joinpath(RUNNER_OUTPUT_DIR_v8, "plots")
    mkpath(RUNNER_OUTPUT_DIR_v8)
    mkpath(RUNNER_PLOT_DIR_v8)
end

begin # runner helpers
    const RUNNER_SENSOR_COLORS_v8 = Dict(
        :T8 => :blue,
        :T9_pair => :red,
        :T10_pair => :green,
        :T3 => :purple,
    )

    function save_runner_plot_v8(plot_object, filename)
        path = joinpath(RUNNER_PLOT_DIR_v8, filename)
        savefig(plot_object, path)
        println("[run_1D_v8] Saved plot: $path")
        return path
    end

    function write_parameters_v8(path, parameters)
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
            println(io, "# fixed,B_Re,$B_RE_FIXED_v8")
            println(io, "# fixed,C_Pr,$C_PR_FIXED_v8")
            println(io, "# fixed,eta_abs,$ETA_ABS_FIXED_v8")
            println(io, "# fixed,beta_opt,$BETA_OPT_FIXED_v8")
            println(io, "# fixed,front_dep,$FRONT_DEPOSITION_FIXED_v8")
            println(io, "# fixed,h_front,$H_FRONT_FIXED_v8")
            println(io, "# fixed,eps_front,$EPS_FRONT_FIXED_v8")
            println(io, "# fixed,T2_boundary_radius_m,$T2_BOUNDARY_RADIUS_v8")
            println(io, "# fixed,G_receiver_to_T2_W_per_m_K,$RECEIVER_TO_T2_CONDUCTANCE_PER_LENGTH_v8")
            println(io, "# fixed,G_receiver_to_tube_W_per_K,$RECEIVER_TO_TUBE_CONDUCTANCE_v8")
            println(io, "# fixed,rear_tube_length_m,$REAR_TUBE_LENGTH_v8")
            println(io, "# fixed,rear_tube_cavity_length_m,$REAR_TUBE_CAVITY_LENGTH_v8")
            println(io, "# fixed,rear_tube_flange_length_m,$REAR_TUBE_FLANGE_LENGTH_v8")
            println(io, "# fixed,rear_tube_wall_thickness_m,$REAR_TUBE_WALL_THICKNESS_v8")
            println(io, "# fixed,water_flange_temperature_K,$WATER_FLANGE_TEMPERATURE_v8")
        end
        return path
    end

    function write_metrics_v8(path, metrics=compute_metrics_v8())
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

    function build_steady_results_v8(params=pnew_v8; keys=sim_key_heat,
                                     nodes=RUNNER_PLOT_NODES_v8)
        results = NamedTuple[]
        for simulation_id in keys
            model, result, experiment = solve_case_v8(params, simulation_id; nodes=nodes)
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
                rear_tube_exit=result.rear_tube_temperature[end, end],
                receiver_to_tube_heat=result.receiver_to_tube_heat[end],
                tube_to_t2_heat_loss=result.tube_to_t2_heat_loss[end],
                flange_heat_loss=result.flange_heat_loss[end],
                h_effective=model[end, 5],
            ))
        end
        return results
    end

    function steady_comparison_plot_v8(steady_results)
        sensors = (:T8, :T9_pair, :T10_pair, :T3)
        plot_object = plot(
            title="1D_v8 steady-state comparison",
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
                     label=String(sensor), color=RUNNER_SENSOR_COLORS_v8[sensor],
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
    println("[run_1D_v8] Loaded $(length(sim_key_heat)) heating runs and " *
            "$(length(sim_key_cool)) cooling runs.")

    if RUNNER_QUICK_v8
        println("[run_1D_v8] Running quick 1D_v8 calibration.")
        calibration_v8 = quick_calibration_v8()
    else
        println("[run_1D_v8] Running full 1D_v8 calibration.")
        calibration_v8 = calibrate_v8(
            nodes=RUNNER_NODES_v8,
            heating_base_iterations=RUNNER_HEAT_BASE_ITERS_v8,
            cooling_iterations=RUNNER_COOL_ITERS_v8,
            heating_refit_iterations=RUNNER_HEAT_REFIT_ITERS_v8,
            heating_base_time_seconds=RUNNER_HEAT_BASE_TIME_v8,
            cooling_time_seconds=RUNNER_COOL_TIME_v8,
            heating_refit_time_seconds=RUNNER_HEAT_REFIT_TIME_v8,
        )
    end

    parameter_path_v8 = joinpath(RUNNER_OUTPUT_DIR_v8, "parameters_1D_v8.csv")
    write_parameters_v8(parameter_path_v8, calibration_v8.parameters)
    println("[run_1D_v8] Saved calibrated parameters: $parameter_path_v8")
end

begin # post-processing and plots
    metrics_v8 = compute_metrics_v8(pnew_v8; nodes=RUNNER_PLOT_NODES_v8)
    metrics_path_v8 = joinpath(RUNNER_OUTPUT_DIR_v8, "analysis_results_1D_v8.csv")
    write_metrics_v8(metrics_path_v8, metrics_v8)
    println("[run_1D_v8] Saved metrics: $metrics_path_v8")

    steady_results_v8 = build_steady_results_v8(pnew_v8; nodes=RUNNER_PLOT_NODES_v8)
    steady_plot_v8 = steady_comparison_plot_v8(steady_results_v8)
    save_runner_plot_v8(steady_plot_v8, "steady_comparison_1D_v8.png")
    display(steady_plot_v8)

    for simulation_id in sim_key_heat
        plot_object = plot_case_v8(simulation_id, pnew_v8; nodes=RUNNER_PLOT_NODES_v8)
        save_runner_plot_v8(plot_object, "transient_$(simulation_id)_1D_v8.png")
        display(plot_object)
    end

    for simulation_id in sim_key_cool
        plot_object = plot_case_v8(simulation_id, pnew_v8;
                                   is_cooling=true, nodes=RUNNER_PLOT_NODES_v8)
        save_runner_plot_v8(plot_object, "transient_$(simulation_id)_1D_v8.png")
        display(plot_object)
    end

    println("[run_1D_v8] Complete.")
    println("[run_1D_v8] Output directory: $RUNNER_OUTPUT_DIR_v8")
end
