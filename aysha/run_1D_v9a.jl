# ============================================================================
# run_1D_v9a.jl - End-to-end 1D_v9a predicted-cavity workflow
# ============================================================================
# Stages:
#   1. Heating fits A/B/C heat-transfer, exchange area, and irradiance factors.
#   2. Cooling fits receiver/cavity conduction scales.
#   3. Heating refits A/B/C heat-transfer, exchange area, and irradiance factors.
#
# Optional environment controls:
#   RECEIVER1D_V9A_RUNNER_QUICK=true
#   RECEIVER1D_V9A_RUNNER_REUSE_PARAMETERS=true
#   RECEIVER1D_V9A_RUNNER_NODES=15
#   RECEIVER1D_V9A_RUNNER_PLOT_NODES=25
#   RECEIVER1D_V9A_RUNNER_HEAT_BASE_ITERS=20
#   RECEIVER1D_V9A_RUNNER_COOL_ITERS=20
#   RECEIVER1D_V9A_RUNNER_HEAT_REFIT_ITERS=20
# ============================================================================

begin # libraries
    using Statistics
    using StatsPlots
end

include("1D_v9a.jl")

begin # runner configuration
    runner_flag_v9a(name, default=false) =
        lowercase(get(ENV, name, string(default))) == "true"

    runner_int_v9a(name, default) = parse(Int, get(ENV, name, string(default)))
    runner_float_v9a(name, default) = parse(Float64, get(ENV, name, string(default)))

    const RUNNER_QUICK_v9a = runner_flag_v9a("RECEIVER1D_V9A_RUNNER_QUICK", false)
    const RUNNER_REUSE_PARAMETERS_v9a =
        runner_flag_v9a("RECEIVER1D_V9A_RUNNER_REUSE_PARAMETERS", false)
    const RUNNER_NODES_v9a = runner_int_v9a("RECEIVER1D_V9A_RUNNER_NODES",
                                          RUNNER_QUICK_v9a ? 11 : 15)
    const RUNNER_PLOT_NODES_v9a = runner_int_v9a("RECEIVER1D_V9A_RUNNER_PLOT_NODES",
                                               default_nodes)
    const RUNNER_HEAT_BASE_ITERS_v9a = runner_int_v9a(
        "RECEIVER1D_V9A_RUNNER_HEAT_BASE_ITERS", RUNNER_QUICK_v9a ? 2 : 20
    )
    const RUNNER_COOL_ITERS_v9a = runner_int_v9a("RECEIVER1D_V9A_RUNNER_COOL_ITERS",
                                               RUNNER_QUICK_v9a ? 2 : 20)
    const RUNNER_HEAT_REFIT_ITERS_v9a = runner_int_v9a(
        "RECEIVER1D_V9A_RUNNER_HEAT_REFIT_ITERS", RUNNER_QUICK_v9a ? 2 : 20
    )
    const RUNNER_COOL_TIME_v9a = runner_float_v9a("RECEIVER1D_V9A_RUNNER_COOL_TIME", 120.0)
    const RUNNER_HEAT_BASE_TIME_v9a = runner_float_v9a(
        "RECEIVER1D_V9A_RUNNER_HEAT_BASE_TIME", 120.0
    )
    const RUNNER_HEAT_REFIT_TIME_v9a = runner_float_v9a(
        "RECEIVER1D_V9A_RUNNER_HEAT_REFIT_TIME", 120.0
    )

    const RUNNER_OUTPUT_DIR_v9a = joinpath(@__DIR__, "summaries", "1D_v9a")
    const RUNNER_PLOT_DIR_v9a = joinpath(RUNNER_OUTPUT_DIR_v9a, "plots")
    mkpath(RUNNER_OUTPUT_DIR_v9a)
    mkpath(RUNNER_PLOT_DIR_v9a)
end

begin # runner helpers
    const RUNNER_SENSOR_COLORS_v9a = Dict(
        :T8 => :blue,
        :T9_pair => :red,
        :T10_pair => :green,
        :T3 => :purple,
        :T2 => :black,
    )

    function save_runner_plot_v9a(plot_object, filename)
        path = joinpath(RUNNER_PLOT_DIR_v9a, filename)
        savefig(plot_object, path)
        println("[run_1D_v9a] Saved plot: $path")
        return path
    end

    function write_parameters_v9a(path, parameters)
        names = (
            "k_scale", "k_ins_scale",
            "A_Nu", "B_Re", "C_Pr", "f_exchange",
            "f_I_high", "f_I_mid", "f_I_low",
        )
        open(path, "w") do io
            println(io, "index,name,value")
            for i in eachindex(parameters)
                println(io, join((i, names[i], parameters[i]), ','))
            end
            println(io, "# fixed,Nu_fd_receiver,$NU_FD_RECEIVER_v9a")
            println(io, "# fixed,eta_abs,$ETA_ABS_FIXED_v9a")
            println(io, "# fixed,beta_opt,$BETA_OPT_FIXED_v9a")
            println(io, "# fixed,front_dep,$FRONT_DEPOSITION_FIXED_v9a")
            println(io, "# fixed,h_front,$H_FRONT_FIXED_v9a")
            println(io, "# fixed,eps_front,$EPS_FRONT_FIXED_v9a")
            println(io, "# fixed,cavity_outer_radius_m,$CAVITY_OUTER_RADIUS_v9a")
            println(io, "# fixed,insulation_outer_radius_m,$INSULATION_OUTER_RADIUS_v9a")
            println(io, "# fixed,G_receiver_to_cavity_W_per_m_K,$RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v9a")
            println(io, "# fixed,G_receiver_to_tube_W_per_K,$RECEIVER_TO_TUBE_CONDUCTANCE_v9a")
            println(io, "# fixed,adaptor_diameter_m,$ADAPTOR_DIAMETER_v9a")
            println(io, "# fixed,rear_tube_length_m,$REAR_TUBE_LENGTH_v9a")
            println(io, "# fixed,rear_tube_cavity_length_m,$REAR_TUBE_CAVITY_LENGTH_v9a")
            println(io, "# fixed,rear_tube_flange_length_m,$REAR_TUBE_FLANGE_LENGTH_v9a")
            println(io, "# fixed,T3_sample_position_m,$T3_SAMPLE_POSITION_v9a")
            println(io, "# fixed,rear_tube_wall_thickness_m,$REAR_TUBE_WALL_THICKNESS_v9a")
            println(io, "# fixed,cavity_heat_capacity_J_per_K,$CAVITY_HEAT_CAPACITY_v9a")
            println(io, "# fixed,water_flange_temperature_K,$WATER_FLANGE_TEMPERATURE_v9a")
        end
        return path
    end

    function read_parameters_v9a(path)
        values = Float64[]
        open(path, "r") do io
            for line in eachline(io)
                startswith(line, "#") && continue
                startswith(line, "index,") && continue
                isempty(strip(line)) && continue
                fields = split(line, ',')
                length(fields) >= 3 || continue
                push!(values, parse(Float64, fields[3]))
            end
        end
        length(values) == length(pnew_v9a) ||
            error("Expected $(length(pnew_v9a)) parameters in $path, found $(length(values)).")
        return values
    end

    function write_metrics_v9a(path, metrics=compute_metrics_v9a())
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

    function build_steady_results_v9a(params=pnew_v9a; keys=sim_key_heat,
                                     nodes=RUNNER_PLOT_NODES_v9a)
        results = NamedTuple[]
        for simulation_id in keys
            model, result, experiment = solve_case_v9a(params, simulation_id; nodes=nodes)
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
                T2_model=model[end, 5],
                T2_experiment=experiment[end, 5],
                rear_tube_exit=result.rear_tube_temperature[end, end],
                receiver_to_tube_heat=result.receiver_to_tube_heat[end],
                receiver_to_cavity_heat=result.receiver_to_cavity_heat[end],
                tube_to_cavity_heat=result.tube_to_cavity_heat[end],
                cavity_ambient_heat_loss=result.cavity_ambient_heat_loss[end],
                flange_heat_loss=result.flange_heat_loss[end],
                h_effective=model[end, 6],
            ))
        end
        return results
    end

    function steady_comparison_plot_v9a(steady_results)
        sensors = (:T8, :T9_pair, :T10_pair, :T3, :T2)
        plot_object = plot(
            title="1D_v9a steady-state comparison",
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
                     label=String(sensor), color=RUNNER_SENSOR_COLORS_v9a[sensor],
                     ms=4, markerstrokewidth=0)
        end

        lower, upper = TEMPERATURE_AXIS_LIMITS_v9a
        plot!(plot_object, [lower, upper], [lower, upper];
              label="1:1", color=:gray, linestyle=:dash)
        xlims!(plot_object, lower, upper)
        ylims!(plot_object, lower, upper)
        return plot_object
    end

    function axial_profile_plot_v9a(simulation_id, params=pnew_v9a;
                                   is_cooling=false,
                                   nodes=RUNNER_PLOT_NODES_v9a)
        model, result, experiment = solve_case_v9a(
            params, simulation_id; is_cooling=is_cooling, nodes=nodes
        )
        conditions = is_cooling ? simulation_conditions_cooling[simulation_id] :
                                  simulation_conditions[simulation_id]
        phase = is_cooling ? "cooling" : "heating"

        solid_x_mm = vcat(result.z_solid, result.z_rear_tube) .* 1000.0
        solid_T = vcat(result.solid_temperature[:, end],
                       result.rear_tube_temperature[:, end])
        gas_x_mm = result.z_gas .* 1000.0
        gas_T = result.gas_temperature[:, end]

        measurement_data = is_cooling ? measurements_cooling : measurements
        solid_measurement_points = (
            (:T8, sensor_positions[:T8], "_T8", :square, :blue),
            (:T9, sensor_positions[:T9], "_T9", :circle, :blue),
            (:T12, sensor_positions[:T9], "_T12", :utriangle, :cyan),
            (:T10, sensor_positions[:T10], "_T10", :diamond, :blue),
            (:T11, sensor_positions[:T10], "_T11", :dtriangle, :cyan),
        )
        t3_x_mm = [T3_SAMPLE_POSITION_v9a * 1000.0]
        t3_experiment_T = [experiment[end, 4]]

        title_text = "1D_v9a axial profile: $simulation_id $phase"
        subtitle_text = "flow=$(round(conditions[qlpm]; digits=2)) L/min, " *
                        "Io=$(round(conditions[Io] / 1000.0; digits=1)) kW/m^2"
        plot_object = plot(
            title=title_text * "\n" * subtitle_text,
            xlabel="Length (mm)",
            ylabel="Temperature (K)",
            legend=:outerright,
            grid=true,
            xlims=(0.0, 1000.0 * (L + REAR_TUBE_LENGTH_v9a)),
            ylims=TEMPERATURE_AXIS_LIMITS_v9a,
        )
        plot!(plot_object, solid_x_mm, solid_T;
              label="solid/tube model", color=:blue, linestyle=:dash, lw=2)
        plot!(plot_object, gas_x_mm, gas_T;
              label="gas model", color=:green, lw=2)
        for (sensor, position, observation_id, marker, color) in solid_measurement_points
            scatter!(plot_object, [position * 1000.0],
                     [observation(measurement_data, simulation_id, observation_id)[end]];
                     label=String(sensor), color=color, markershape=marker,
                     markerstrokewidth=0, ms=5)
        end
        scatter!(plot_object, t3_x_mm, t3_experiment_T;
                 label="T3 experiment", color=:red, markershape=:star5,
                 markerstrokewidth=0, ms=7)
        scatter!(plot_object, [CAVITY_LENGTH_v9a * 1000.0], [experiment[end, 5]];
                 label="T2 experiment", color=:black, markershape=:diamond,
                 markerstrokewidth=0, ms=5)
        scatter!(plot_object, [CAVITY_LENGTH_v9a * 1000.0], [model[end, 5]];
                 label="T2 model", color=:black, markershape=:xcross,
                 markerstrokewidth=1, ms=5)
        vline!(plot_object, [L * 1000.0, (L + REAR_TUBE_CAVITY_LENGTH_v9a) * 1000.0];
               label=false, color=:gray, linestyle=:dot, lw=1)
        return plot_object
    end
end

begin # calibration
    println("[run_1D_v9a] Loaded $(length(sim_key_heat)) heating runs and " *
            "$(length(sim_key_cool)) cooling runs.")

    parameter_path_v9a = joinpath(RUNNER_OUTPUT_DIR_v9a, "parameters_1D_v9a.csv")
    if RUNNER_REUSE_PARAMETERS_v9a
        println("[run_1D_v9a] Reusing saved calibrated parameters.")
        global pnew_v9a = read_parameters_v9a(parameter_path_v9a)
        global calibration_v9a = (parameters=pnew_v9a,)
    elseif RUNNER_QUICK_v9a
        println("[run_1D_v9a] Running quick 1D_v9a calibration.")
        calibration_v9a = quick_calibration_v9a()
    else
        println("[run_1D_v9a] Running full 1D_v9a calibration.")
        calibration_v9a = calibrate_v9a(
            nodes=RUNNER_NODES_v9a,
            heating_base_iterations=RUNNER_HEAT_BASE_ITERS_v9a,
            cooling_iterations=RUNNER_COOL_ITERS_v9a,
            heating_refit_iterations=RUNNER_HEAT_REFIT_ITERS_v9a,
            heating_base_time_seconds=RUNNER_HEAT_BASE_TIME_v9a,
            cooling_time_seconds=RUNNER_COOL_TIME_v9a,
            heating_refit_time_seconds=RUNNER_HEAT_REFIT_TIME_v9a,
        )
    end

    if RUNNER_REUSE_PARAMETERS_v9a
        println("[run_1D_v9a] Kept calibrated parameters unchanged: $parameter_path_v9a")
    else
        write_parameters_v9a(parameter_path_v9a, calibration_v9a.parameters)
        println("[run_1D_v9a] Saved calibrated parameters: $parameter_path_v9a")
    end
end

begin # post-processing and plots
    metrics_v9a = compute_metrics_v9a(pnew_v9a; nodes=RUNNER_PLOT_NODES_v9a)
    metrics_path_v9a = joinpath(RUNNER_OUTPUT_DIR_v9a, "analysis_results_1D_v9a.csv")
    write_metrics_v9a(metrics_path_v9a, metrics_v9a)
    println("[run_1D_v9a] Saved metrics: $metrics_path_v9a")

    steady_results_v9a = build_steady_results_v9a(pnew_v9a; nodes=RUNNER_PLOT_NODES_v9a)
    steady_plot_v9a = steady_comparison_plot_v9a(steady_results_v9a)
    save_runner_plot_v9a(steady_plot_v9a, "steady_comparison_1D_v9a.png")
    display(steady_plot_v9a)

    for simulation_id in sim_key_heat
        plot_object = plot_case_v9a(simulation_id, pnew_v9a; nodes=RUNNER_PLOT_NODES_v9a)
        save_runner_plot_v9a(plot_object, "transient_$(simulation_id)_1D_v9a.png")
        display(plot_object)

        profile_plot_object = axial_profile_plot_v9a(
            simulation_id, pnew_v9a; nodes=RUNNER_PLOT_NODES_v9a
        )
        save_runner_plot_v9a(profile_plot_object,
                             "axial_profile_$(simulation_id)_1D_v9a.png")
        display(profile_plot_object)
    end

    for simulation_id in sim_key_cool
        plot_object = plot_case_v9a(simulation_id, pnew_v9a;
                                   is_cooling=true, nodes=RUNNER_PLOT_NODES_v9a)
        save_runner_plot_v9a(plot_object, "transient_$(simulation_id)_1D_v9a.png")
        display(plot_object)

        profile_plot_object = axial_profile_plot_v9a(
            simulation_id, pnew_v9a; is_cooling=true, nodes=RUNNER_PLOT_NODES_v9a
        )
        save_runner_plot_v9a(profile_plot_object,
                             "axial_profile_$(simulation_id)_1D_v9a.png")
        display(profile_plot_object)
    end

    println("[run_1D_v9a] Complete.")
    println("[run_1D_v9a] Output directory: $RUNNER_OUTPUT_DIR_v9a")
end
