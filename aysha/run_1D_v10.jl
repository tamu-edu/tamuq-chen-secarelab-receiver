# ============================================================================
# run_1D_v10.jl - v10 ECM forward-comparison workflow
# ============================================================================
# Runs a small comparison matrix before any v10 coefficient fitting:
#   M1: monolith Nu + fixed properties, no Rosseland axial radiation
#   M2: M1 + Rosseland axial radiation at beta_tr = 1000 1/m
#   M3_low_beta: Rosseland sensitivity at beta_tr = 300 1/m
#   M3_high_beta: Rosseland sensitivity at beta_tr = 2700 1/m
#
# Optional environment controls:
#   RECEIVER1D_V10_RUNNER_PLOT_NODES=25
# ============================================================================

begin # libraries
    using Statistics
    using StatsPlots
end

include("1D_v10.jl")

begin # runner configuration
    runner_flag_v10(name, default=false) =
        lowercase(get(ENV, name, string(default))) == "true"

    runner_int_v10(name, default) = parse(Int, get(ENV, name, string(default)))
    runner_float_v10(name, default) = parse(Float64, get(ENV, name, string(default)))

    const RUNNER_PLOT_NODES_v10 = runner_int_v10("RECEIVER1D_V10_RUNNER_PLOT_NODES",
                                               default_nodes)

    const RUNNER_OUTPUT_DIR_v10 = joinpath(@__DIR__, "summaries", "1D_v10")
    const RUNNER_PLOT_DIR_v10 = joinpath(RUNNER_OUTPUT_DIR_v10, "plots")
    mkpath(RUNNER_OUTPUT_DIR_v10)
    mkpath(RUNNER_PLOT_DIR_v10)
end

begin # runner helpers
    const RUNNER_SENSOR_COLORS_v10 = Dict(
        :T8 => :blue,
        :T9_pair => :red,
        :T10_pair => :green,
        :T3 => :purple,
        :T2 => :black,
    )

    function save_runner_plot_v10(plot_object, filename)
        path = joinpath(RUNNER_PLOT_DIR_v10, filename)
        savefig(plot_object, path)
        println("[run_1D_v10] Saved plot: $path")
        return path
    end

    function write_parameters_v10(path, parameters)
        names = (
            "C_Nu_model",
            "ell_rad_m",
            "eta_opt",
            "front_dep",
            "beta_opt",
        )
        open(path, "w") do io
            println(io, "index,name,value")
            for i in eachindex(parameters)
                println(io, join((i, names[i], parameters[i]), ','))
            end
            println(io, "# fixed,Nu_fd_receiver,$NU_FD_RECEIVER_v10")
            println(io, "# fixed,eta_abs,$ETA_ABS_FIXED_v10")
            println(io, "# derived,beta_tr,$(rosseland_beta_v10(parameters))")
            println(io, "# derived,tau_receiver,$(rosseland_optical_thickness_v10(L, parameters))")
            println(io, "# derived,tau_channel,$(rosseland_optical_thickness_v10(Dh, parameters))")
            println(io, "# fixed,h_front,$H_FRONT_FIXED_v10")
            println(io, "# fixed,eps_front,$EPS_FRONT_FIXED_v10")
            println(io, "# fixed,cavity_outer_radius_m,$CAVITY_OUTER_RADIUS_v10")
            println(io, "# fixed,insulation_outer_radius_m,$INSULATION_OUTER_RADIUS_v10")
            println(io, "# fixed,G_receiver_to_cavity_W_per_m_K,$RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v10")
            println(io, "# fixed,G_receiver_to_tube_W_per_K,$RECEIVER_TO_TUBE_CONDUCTANCE_v10")
            println(io, "# fixed,adaptor_diameter_m,$ADAPTOR_DIAMETER_v10")
            println(io, "# fixed,rear_tube_length_m,$REAR_TUBE_LENGTH_v10")
            println(io, "# fixed,rear_tube_cavity_length_m,$REAR_TUBE_CAVITY_LENGTH_v10")
            println(io, "# fixed,rear_tube_flange_length_m,$REAR_TUBE_FLANGE_LENGTH_v10")
            println(io, "# fixed,T3_sample_position_m,$T3_SAMPLE_POSITION_v10")
            println(io, "# fixed,rear_tube_wall_thickness_m,$REAR_TUBE_WALL_THICKNESS_v10")
            println(io, "# fixed,cavity_heat_capacity_J_per_K,$CAVITY_HEAT_CAPACITY_v10")
            println(io, "# fixed,water_flange_temperature_K,$WATER_FLANGE_TEMPERATURE_v10")
        end
        return path
    end

    function read_parameters_v10(path)
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
        length(values) == length(pnew_v10) ||
            error("Expected $(length(pnew_v10)) parameters in $path, found $(length(values)).")
        return values
    end

    function write_metrics_v10(path, metrics=compute_metrics_v10())
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

    function build_steady_results_v10(params=pnew_v10; keys=sim_key_heat,
                                     nodes=RUNNER_PLOT_NODES_v10,
                                     use_rosseland_radiation=false)
        results = NamedTuple[]
        for simulation_id in keys
            model, result, experiment = solve_case_v10(
                params, simulation_id; nodes=nodes,
                use_rosseland_radiation=use_rosseland_radiation,
            )
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

    function steady_comparison_plot_v10(steady_results)
        sensors = (:T8, :T9_pair, :T10_pair, :T3, :T2)
        plot_object = plot(
            title="1D_v10 steady-state comparison",
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
                     label=String(sensor), color=RUNNER_SENSOR_COLORS_v10[sensor],
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

    function axial_profile_plot_v10(simulation_id, params=pnew_v10;
                                   is_cooling=false,
                                   nodes=RUNNER_PLOT_NODES_v10,
                                   use_rosseland_radiation=false,
                                   variant_name="v10")
        model, result, experiment = solve_case_v10(
            params, simulation_id; is_cooling=is_cooling, nodes=nodes,
            use_rosseland_radiation=use_rosseland_radiation,
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
        t3_x_mm = [T3_SAMPLE_POSITION_v10 * 1000.0]
        t3_experiment_T = [experiment[end, 4]]

        title_text = "1D_v10 $variant_name axial profile: $simulation_id $phase"
        subtitle_text = "flow=$(round(conditions[qlpm]; digits=2)) L/min, " *
                        "Io=$(round(conditions[Io] / 1000.0; digits=1)) kW/m^2"
        plot_object = plot(
            title=title_text * "\n" * subtitle_text,
            xlabel="Length (mm)",
            ylabel="Temperature (K)",
            legend=:outerright,
            grid=true,
            xlims=(0.0, 1000.0 * (L + REAR_TUBE_LENGTH_v10)),
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
        scatter!(plot_object, [CAVITY_LENGTH_v10 * 1000.0], [experiment[end, 5]];
                 label="T2 experiment", color=:black, markershape=:diamond,
                 markerstrokewidth=0, ms=5)
        scatter!(plot_object, [CAVITY_LENGTH_v10 * 1000.0], [model[end, 5]];
                 label="T2 model", color=:black, markershape=:xcross,
                 markerstrokewidth=1, ms=5)
        vline!(plot_object, [L * 1000.0, (L + REAR_TUBE_CAVITY_LENGTH_v10) * 1000.0];
               label=false, color=:gray, linestyle=:dot, lw=1)
        return plot_object
    end
end

begin # v10 comparison matrix
    println("[run_1D_v10] Loaded $(length(sim_key_heat)) heating runs and " *
            "$(length(sim_key_cool)) cooling runs.")

    variant_params_v10(beta_tr) = begin
        p = copy(pnew_v10)
        p[2] = 1.0 / beta_tr
        p
    end

    variants_v10 = (
        (name="M1_no_rad", params=copy(pnew_v10), use_rad=false),
        (name="M2_rad_beta1000", params=variant_params_v10(1000.0), use_rad=true),
        (name="M3_rad_beta300", params=variant_params_v10(300.0), use_rad=true),
        (name="M3_rad_beta2700", params=variant_params_v10(2700.0), use_rad=true),
    )

    representative_heat = ("E67", "E76", "E80")
    representative_cool = ("C80",)
    summary_rows = NamedTuple[]

    for variant in variants_v10
        println("[run_1D_v10] Running variant $(variant.name).")
        parameter_path = joinpath(RUNNER_OUTPUT_DIR_v10, "parameters_$(variant.name)_1D_v10.csv")
        write_parameters_v10(parameter_path, variant.params)

        metrics = compute_metrics_v10(
            variant.params; nodes=RUNNER_PLOT_NODES_v10,
            use_rosseland_radiation=variant.use_rad,
        )
        metrics_path = joinpath(RUNNER_OUTPUT_DIR_v10, "analysis_results_$(variant.name)_1D_v10.csv")
        write_metrics_v10(metrics_path, metrics)
        println("[run_1D_v10] Saved metrics: $metrics_path")

        for row in metrics
            push!(summary_rows, merge(row, (variant=variant.name,)))
        end

        steady_results = build_steady_results_v10(
            variant.params; nodes=RUNNER_PLOT_NODES_v10,
            use_rosseland_radiation=variant.use_rad,
        )
        steady_plot = steady_comparison_plot_v10(steady_results)
        plot!(steady_plot, title="1D_v10 $(variant.name) steady-state comparison")
        save_runner_plot_v10(steady_plot, "steady_comparison_$(variant.name)_1D_v10.png")

        for simulation_id in representative_heat
            profile_plot = axial_profile_plot_v10(
                simulation_id, variant.params; nodes=RUNNER_PLOT_NODES_v10,
                use_rosseland_radiation=variant.use_rad,
                variant_name=variant.name,
            )
            save_runner_plot_v10(profile_plot,
                                 "axial_profile_$(simulation_id)_$(variant.name)_1D_v10.png")
        end

        for simulation_id in representative_cool
            profile_plot = axial_profile_plot_v10(
                simulation_id, variant.params; is_cooling=true,
                nodes=RUNNER_PLOT_NODES_v10,
                use_rosseland_radiation=variant.use_rad,
                variant_name=variant.name,
            )
            save_runner_plot_v10(profile_plot,
                                 "axial_profile_$(simulation_id)_$(variant.name)_1D_v10.png")
        end
    end

    summary_path = joinpath(RUNNER_OUTPUT_DIR_v10, "analysis_results_all_variants_1D_v10.csv")
    open(summary_path, "w") do io
        println(io, "variant,simulation_id,phase,sensor,rmse_K,steady_error_K,t90_error_s,shape_loss")
        for row in summary_rows
            println(io, join((row.variant, row.simulation_id, row.phase, row.sensor,
                              row.rmse_K, row.steady_error_K,
                              row.t90_error_s, row.shape_loss), ','))
        end
    end
    println("[run_1D_v10] Saved combined variant metrics: $summary_path")

    println("[run_1D_v10] Complete.")
    println("[run_1D_v10] Output directory: $RUNNER_OUTPUT_DIR_v10")
end
