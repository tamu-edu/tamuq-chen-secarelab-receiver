# ============================================================================
# run_1D_v8b.jl - End-to-end 1D_v8b predicted-cavity workflow
# ============================================================================
# Stages:
#   1. Heating fits gas heat-transfer shape and irradiance-level factors.
#   2. Cooling fits receiver thermophysical parameters with predicted T2/cavity.
#   3. Heating refits gas heat-transfer shape and irradiance-level factors.
#
# Optional environment controls:
#   RECEIVER1D_V8B_RUNNER_QUICK=true
#   RECEIVER1D_V8B_RUNNER_REUSE_PARAMETERS=true
#   RECEIVER1D_V8B_RUNNER_NODES=15
#   RECEIVER1D_V8B_RUNNER_PLOT_NODES=25
#   RECEIVER1D_V8B_RUNNER_HEAT_BASE_ITERS=20
#   RECEIVER1D_V8B_RUNNER_COOL_ITERS=20
#   RECEIVER1D_V8B_RUNNER_HEAT_REFIT_ITERS=20
# ============================================================================

begin # libraries
    using Statistics
    using StatsPlots
end

include("1D_v8b.jl")

begin # runner configuration
    runner_flag_v8b(name, default=false) =
        lowercase(get(ENV, name, string(default))) == "true"

    runner_int_v8b(name, default) = parse(Int, get(ENV, name, string(default)))
    runner_float_v8b(name, default) = parse(Float64, get(ENV, name, string(default)))

    const RUNNER_QUICK_v8b = runner_flag_v8b("RECEIVER1D_V8B_RUNNER_QUICK", false)
    const RUNNER_REUSE_PARAMETERS_v8b =
        runner_flag_v8b("RECEIVER1D_V8B_RUNNER_REUSE_PARAMETERS", false)
    const RUNNER_NODES_v8b = runner_int_v8b("RECEIVER1D_V8B_RUNNER_NODES",
                                          RUNNER_QUICK_v8b ? 11 : 15)
    const RUNNER_PLOT_NODES_v8b = runner_int_v8b("RECEIVER1D_V8B_RUNNER_PLOT_NODES",
                                               default_nodes)
    const RUNNER_HEAT_BASE_ITERS_v8b = runner_int_v8b(
        "RECEIVER1D_V8B_RUNNER_HEAT_BASE_ITERS", RUNNER_QUICK_v8b ? 2 : 20
    )
    const RUNNER_COOL_ITERS_v8b = runner_int_v8b("RECEIVER1D_V8B_RUNNER_COOL_ITERS",
                                               RUNNER_QUICK_v8b ? 2 : 20)
    const RUNNER_HEAT_REFIT_ITERS_v8b = runner_int_v8b(
        "RECEIVER1D_V8B_RUNNER_HEAT_REFIT_ITERS", RUNNER_QUICK_v8b ? 2 : 20
    )
    const RUNNER_COOL_TIME_v8b = runner_float_v8b("RECEIVER1D_V8B_RUNNER_COOL_TIME", 120.0)
    const RUNNER_HEAT_BASE_TIME_v8b = runner_float_v8b(
        "RECEIVER1D_V8B_RUNNER_HEAT_BASE_TIME", 120.0
    )
    const RUNNER_HEAT_REFIT_TIME_v8b = runner_float_v8b(
        "RECEIVER1D_V8B_RUNNER_HEAT_REFIT_TIME", 120.0
    )

    const RUNNER_OUTPUT_DIR_v8b = joinpath(@__DIR__, "summaries", "1D_v8b")
    const RUNNER_PLOT_DIR_v8b = joinpath(RUNNER_OUTPUT_DIR_v8b, "plots")
    mkpath(RUNNER_OUTPUT_DIR_v8b)
    mkpath(RUNNER_PLOT_DIR_v8b)
end

begin # runner helpers
    const RUNNER_SENSOR_COLORS_v8b = Dict(
        :T8 => :blue,
        :T9_pair => :red,
        :T10_pair => :green,
        :T3 => :purple,
        :T2 => :black,
    )

    function save_runner_plot_v8b(plot_object, filename)
        path = joinpath(RUNNER_PLOT_DIR_v8b, filename)
        savefig(plot_object, path)
        println("[run_1D_v8b] Saved plot: $path")
        return path
    end

    function write_parameters_v8b(path, parameters)
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
            println(io, "# fixed,B_Re,$B_RE_FIXED_v8b")
            println(io, "# fixed,C_Pr,$C_PR_FIXED_v8b")
            println(io, "# fixed,eta_abs,$ETA_ABS_FIXED_v8b")
            println(io, "# fixed,beta_opt,$BETA_OPT_FIXED_v8b")
            println(io, "# fixed,front_dep,$FRONT_DEPOSITION_FIXED_v8b")
            println(io, "# fixed,h_front,$H_FRONT_FIXED_v8b")
            println(io, "# fixed,eps_front,$EPS_FRONT_FIXED_v8b")
            println(io, "# fixed,cavity_outer_radius_m,$CAVITY_OUTER_RADIUS_v8b")
            println(io, "# fixed,insulation_outer_radius_m,$INSULATION_OUTER_RADIUS_v8b")
            println(io, "# fixed,G_receiver_to_cavity_W_per_m_K,$RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v8b")
            println(io, "# fixed,G_receiver_to_tube_W_per_K,$RECEIVER_TO_TUBE_CONDUCTANCE_v8b")
            println(io, "# fixed,adaptor_diameter_m,$ADAPTOR_DIAMETER_v8b")
            println(io, "# fixed,rear_tube_length_m,$REAR_TUBE_LENGTH_v8b")
            println(io, "# fixed,rear_tube_cavity_length_m,$REAR_TUBE_CAVITY_LENGTH_v8b")
            println(io, "# fixed,rear_tube_flange_length_m,$REAR_TUBE_FLANGE_LENGTH_v8b")
            println(io, "# fixed,T3_sample_position_m,$T3_SAMPLE_POSITION_v8b")
            println(io, "# fixed,rear_tube_wall_thickness_m,$REAR_TUBE_WALL_THICKNESS_v8b")
            println(io, "# fixed,cavity_heat_capacity_J_per_K,$CAVITY_HEAT_CAPACITY_v8b")
            println(io, "# fixed,water_flange_temperature_K,$WATER_FLANGE_TEMPERATURE_v8b")
        end
        return path
    end

    function read_parameters_v8b(path)
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
        length(values) == length(pnew_v8b) ||
            error("Expected $(length(pnew_v8b)) parameters in $path, found $(length(values)).")
        return values
    end

    function write_metrics_v8b(path, metrics=compute_metrics_v8b())
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

    function build_steady_results_v8b(params=pnew_v8b; keys=sim_key_heat,
                                     nodes=RUNNER_PLOT_NODES_v8b)
        results = NamedTuple[]
        for simulation_id in keys
            model, result, experiment = solve_case_v8b(params, simulation_id; nodes=nodes)
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

    function steady_comparison_plot_v8b(steady_results)
        sensors = (:T8, :T9_pair, :T10_pair, :T3, :T2)
        plot_object = plot(
            title="1D_v8b steady-state comparison",
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
                     label=String(sensor), color=RUNNER_SENSOR_COLORS_v8b[sensor],
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

    function axial_profile_plot_v8b(simulation_id, params=pnew_v8b;
                                   is_cooling=false,
                                   nodes=RUNNER_PLOT_NODES_v8b)
        model, result, experiment = solve_case_v8b(
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
        t3_x_mm = [T3_SAMPLE_POSITION_v8b * 1000.0]
        t3_experiment_T = [experiment[end, 4]]

        title_text = "1D_v8b axial profile: $simulation_id $phase"
        subtitle_text = "flow=$(round(conditions[qlpm]; digits=2)) L/min, " *
                        "Io=$(round(conditions[Io] / 1000.0; digits=1)) kW/m^2"
        plot_object = plot(
            title=title_text * "\n" * subtitle_text,
            xlabel="Length (mm)",
            ylabel="Temperature (K)",
            legend=:outerright,
            grid=true,
            xlims=(0.0, 1000.0 * (L + REAR_TUBE_LENGTH_v8b)),
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
        scatter!(plot_object, [CAVITY_LENGTH_v8b * 1000.0], [experiment[end, 5]];
                 label="T2 experiment", color=:black, markershape=:diamond,
                 markerstrokewidth=0, ms=5)
        scatter!(plot_object, [CAVITY_LENGTH_v8b * 1000.0], [model[end, 5]];
                 label="T2 model", color=:black, markershape=:xcross,
                 markerstrokewidth=1, ms=5)
        vline!(plot_object, [L * 1000.0, (L + REAR_TUBE_CAVITY_LENGTH_v8b) * 1000.0];
               label=false, color=:gray, linestyle=:dot, lw=1)
        return plot_object
    end
end

begin # calibration
    println("[run_1D_v8b] Loaded $(length(sim_key_heat)) heating runs and " *
            "$(length(sim_key_cool)) cooling runs.")

    parameter_path_v8b = joinpath(RUNNER_OUTPUT_DIR_v8b, "parameters_1D_v8b.csv")
    if RUNNER_REUSE_PARAMETERS_v8b
        println("[run_1D_v8b] Reusing saved calibrated parameters.")
        global pnew_v8b = read_parameters_v8b(parameter_path_v8b)
        global calibration_v8b = (parameters=pnew_v8b,)
    elseif RUNNER_QUICK_v8b
        println("[run_1D_v8b] Running quick 1D_v8b calibration.")
        calibration_v8b = quick_calibration_v8b()
    else
        println("[run_1D_v8b] Running full 1D_v8b calibration.")
        calibration_v8b = calibrate_v8b(
            nodes=RUNNER_NODES_v8b,
            heating_base_iterations=RUNNER_HEAT_BASE_ITERS_v8b,
            cooling_iterations=RUNNER_COOL_ITERS_v8b,
            heating_refit_iterations=RUNNER_HEAT_REFIT_ITERS_v8b,
            heating_base_time_seconds=RUNNER_HEAT_BASE_TIME_v8b,
            cooling_time_seconds=RUNNER_COOL_TIME_v8b,
            heating_refit_time_seconds=RUNNER_HEAT_REFIT_TIME_v8b,
        )
    end

    if RUNNER_REUSE_PARAMETERS_v8b
        println("[run_1D_v8b] Kept calibrated parameters unchanged: $parameter_path_v8b")
    else
        write_parameters_v8b(parameter_path_v8b, calibration_v8b.parameters)
        println("[run_1D_v8b] Saved calibrated parameters: $parameter_path_v8b")
    end
end

begin # post-processing and plots
    metrics_v8b = compute_metrics_v8b(pnew_v8b; nodes=RUNNER_PLOT_NODES_v8b)
    metrics_path_v8b = joinpath(RUNNER_OUTPUT_DIR_v8b, "analysis_results_1D_v8b.csv")
    write_metrics_v8b(metrics_path_v8b, metrics_v8b)
    println("[run_1D_v8b] Saved metrics: $metrics_path_v8b")

    steady_results_v8b = build_steady_results_v8b(pnew_v8b; nodes=RUNNER_PLOT_NODES_v8b)
    steady_plot_v8b = steady_comparison_plot_v8b(steady_results_v8b)
    save_runner_plot_v8b(steady_plot_v8b, "steady_comparison_1D_v8b.png")
    display(steady_plot_v8b)

    for simulation_id in sim_key_heat
        plot_object = plot_case_v8b(simulation_id, pnew_v8b; nodes=RUNNER_PLOT_NODES_v8b)
        save_runner_plot_v8b(plot_object, "transient_$(simulation_id)_1D_v8b.png")
        display(plot_object)

        profile_plot_object = axial_profile_plot_v8b(
            simulation_id, pnew_v8b; nodes=RUNNER_PLOT_NODES_v8b
        )
        save_runner_plot_v8b(profile_plot_object,
                             "axial_profile_$(simulation_id)_1D_v8b.png")
        display(profile_plot_object)
    end

    for simulation_id in sim_key_cool
        plot_object = plot_case_v8b(simulation_id, pnew_v8b;
                                   is_cooling=true, nodes=RUNNER_PLOT_NODES_v8b)
        save_runner_plot_v8b(plot_object, "transient_$(simulation_id)_1D_v8b.png")
        display(plot_object)

        profile_plot_object = axial_profile_plot_v8b(
            simulation_id, pnew_v8b; is_cooling=true, nodes=RUNNER_PLOT_NODES_v8b
        )
        save_runner_plot_v8b(profile_plot_object,
                             "axial_profile_$(simulation_id)_1D_v8b.png")
        display(profile_plot_object)
    end

    println("[run_1D_v8b] Complete.")
    println("[run_1D_v8b] Output directory: $RUNNER_OUTPUT_DIR_v8b")
end
