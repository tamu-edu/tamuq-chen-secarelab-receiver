# ============================================================================
# run_1D_v11.jl - v11 Stage 1 flow-shape isolation workflow
# ============================================================================
# v11 keeps eta_abs, eta_opt, beta_opt, and front_dep fixed, with all absorbed
# power placed in the front cell. The only optimized family is the receiver
# gas-solid entry/development Nu shape:
#
#   p = [A_Nu, B_Re, C_Pr, L_entry_m, ell_rad_m]
#
# Optional environment controls:
#   RECEIVER1D_V11_RUNNER_PLOT_NODES=25
#   RECEIVER1D_V11_FIT_NODES=15
#   RECEIVER1D_V11_FIT_ITERATIONS=80
#   RECEIVER1D_V11_FIT_SECONDS=240
#   RECEIVER1D_V11_SKIP_FIT=false
# ============================================================================

begin # libraries
    using Statistics
    using StatsPlots
end

include("1D_v11.jl")

begin # runner configuration
    runner_flag_v11(name, default=false) =
        lowercase(get(ENV, name, string(default))) == "true"

    runner_int_v11(name, default) = parse(Int, get(ENV, name, string(default)))
    runner_float_v11(name, default) = parse(Float64, get(ENV, name, string(default)))

    const RUNNER_PLOT_NODES_v11 = runner_int_v11(
        "RECEIVER1D_V11_RUNNER_PLOT_NODES", default_nodes,
    )
    const RUNNER_FIT_NODES_v11 = runner_int_v11("RECEIVER1D_V11_FIT_NODES", 15)
    const RUNNER_FIT_ITERATIONS_v11 = runner_int_v11(
        "RECEIVER1D_V11_FIT_ITERATIONS", 80,
    )
    const RUNNER_FIT_SECONDS_v11 = runner_float_v11(
        "RECEIVER1D_V11_FIT_SECONDS", 240.0,
    )
    const RUNNER_SKIP_FIT_v11 = runner_flag_v11("RECEIVER1D_V11_SKIP_FIT", false)

    const RUNNER_OUTPUT_DIR_v11 = joinpath(@__DIR__, "summaries", "1D_v11")
    const RUNNER_PLOT_DIR_v11 = joinpath(RUNNER_OUTPUT_DIR_v11, "plots")
    const RUNNER_DIAGNOSTIC_DIR_v11 = joinpath(RUNNER_OUTPUT_DIR_v11, "diagnostics")
    mkpath(RUNNER_OUTPUT_DIR_v11)
    mkpath(RUNNER_PLOT_DIR_v11)
    mkpath(RUNNER_DIAGNOSTIC_DIR_v11)
end

begin # runner helpers
    const RUNNER_SENSOR_COLORS_v11 = Dict(
        :T8 => :blue,
        :T9_pair => :red,
        :T10_pair => :green,
        :T3 => :purple,
        :T2 => :black,
    )

    function save_runner_plot_v11(plot_object, filename)
        path = joinpath(RUNNER_PLOT_DIR_v11, filename)
        savefig(plot_object, path)
        println("[run_1D_v11] Saved plot: $path")
        return path
    end

    function write_parameters_v11(path, parameters)
        names = (
            "A_Nu",
            "B_Re",
            "C_Pr",
            "L_entry_m",
            "ell_rad_m",
        )
        open(path, "w") do io
            println(io, "index,name,value")
            for i in eachindex(parameters)
                println(io, join((i, names[i], parameters[i]), ','))
            end
            println(io, "# fixed,Nu_fd_receiver,$NU_FD_RECEIVER_v11")
            println(io, "# fixed,entry_shape_exponent,$ENTRY_SHAPE_EXPONENT_FIXED_v11")
            println(io, "# fixed,eta_abs,$ETA_ABS_FIXED_v11")
            println(io, "# fixed,eta_opt,$ETA_OPT_FIXED_v11")
            println(io, "# fixed,beta_opt,$BETA_OPT_FIXED_v11")
            println(io, "# fixed,front_dep,$FRONT_DEPOSITION_STAGE1_v11")
            println(io, "# fixed,h_front,$H_FRONT_FIXED_v11")
            println(io, "# fixed,eps_front,$EPS_FRONT_FIXED_v11")
            println(io, "# derived,beta_tr,$(rosseland_beta_v11(parameters))")
            println(io, "# derived,tau_receiver,$(rosseland_optical_thickness_v11(L, parameters))")
            println(io, "# derived,tau_channel,$(rosseland_optical_thickness_v11(Dh, parameters))")
            println(io, "# fixed,cavity_outer_radius_m,$CAVITY_OUTER_RADIUS_v11")
            println(io, "# fixed,insulation_outer_radius_m,$INSULATION_OUTER_RADIUS_v11")
            println(io, "# fixed,G_receiver_to_cavity_W_per_m_K,$RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v11")
            println(io, "# fixed,G_receiver_to_tube_W_per_K,$RECEIVER_TO_TUBE_CONDUCTANCE_v11")
            println(io, "# fixed,adaptor_diameter_m,$ADAPTOR_DIAMETER_v11")
            println(io, "# fixed,rear_tube_length_m,$REAR_TUBE_LENGTH_v11")
            println(io, "# fixed,rear_tube_cavity_length_m,$REAR_TUBE_CAVITY_LENGTH_v11")
            println(io, "# fixed,rear_tube_flange_length_m,$REAR_TUBE_FLANGE_LENGTH_v11")
            println(io, "# fixed,T3_sample_position_m,$T3_SAMPLE_POSITION_v11")
            println(io, "# fixed,cavity_heat_capacity_J_per_K,$CAVITY_HEAT_CAPACITY_v11")
            println(io, "# fixed,water_flange_temperature_K,$WATER_FLANGE_TEMPERATURE_v11")
        end
        return path
    end

    function write_metrics_v11(path, metrics=compute_metrics_v11())
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

    function build_steady_results_v11(params=pnew_v11; keys=sim_key_heat,
                                     nodes=RUNNER_PLOT_NODES_v11)
        results = NamedTuple[]
        for simulation_id in keys
            model, result, experiment = solve_case_v11(params, simulation_id; nodes=nodes)
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
                T9_minus_T8_model=model[end, 2] - model[end, 1],
                T9_minus_T8_experiment=experiment[end, 2] - experiment[end, 1],
                T10_minus_T8_model=model[end, 3] - model[end, 1],
                T10_minus_T8_experiment=experiment[end, 3] - experiment[end, 1],
                rear_tube_exit=result.rear_tube_temperature[end, end],
                receiver_to_tube_heat=result.receiver_to_tube_heat[end],
                receiver_to_cavity_heat=result.receiver_to_cavity_heat[end],
                tube_to_cavity_heat=result.tube_to_cavity_heat[end],
                cavity_ambient_heat_loss=result.cavity_ambient_heat_loss[end],
                flange_heat_loss=result.flange_heat_loss[end],
                mean_nu=mean(result.receiver_nusselt[:, end]),
                max_effectiveness=maximum(result.receiver_effectiveness[:, end]),
                mean_effectiveness=mean(result.receiver_effectiveness[:, end]),
                h_effective=model[end, 6],
            ))
        end
        return results
    end

    function write_steady_results_v11(path, steady_results)
        fields = propertynames(first(steady_results))
        open(path, "w") do io
            println(io, join(fields, ','))
            for row in steady_results
                println(io, join((getproperty(row, field) for field in fields), ','))
            end
        end
        return path
    end

    function steady_comparison_plot_v11(steady_results; title_suffix="")
        sensors = (:T8, :T9_pair, :T10_pair, :T3, :T2)
        plot_object = plot(
            title="1D_v11 $title_suffix steady-state comparison",
            xlabel="Model temperature (K)",
            ylabel="Experiment temperature (K)",
            legend=:outerright,
            aspect_ratio=:equal,
            xlims=TEMPERATURE_AXIS_LIMITS_v11,
            ylims=TEMPERATURE_AXIS_LIMITS_v11,
        )
        for sensor in sensors
            model_values = [getproperty(row, Symbol(sensor, "_model")) for row in steady_results]
            experiment_values = [getproperty(row, Symbol(sensor, "_experiment")) for row in steady_results]
            scatter!(plot_object, model_values, experiment_values;
                     label=String(sensor), color=RUNNER_SENSOR_COLORS_v11[sensor],
                     ms=4, markerstrokewidth=0)
        end
        lower, upper = TEMPERATURE_AXIS_LIMITS_v11
        plot!(plot_object, [lower, upper], [lower, upper];
              label="1:1", color=:gray, linestyle=:dash)
        return plot_object
    end

    function axial_profile_plot_v11(simulation_id, params=pnew_v11;
                                   is_cooling=false,
                                   nodes=RUNNER_PLOT_NODES_v11,
                                   variant_name="v11")
        model, result, experiment = solve_case_v11(
            params, simulation_id; is_cooling=is_cooling, nodes=nodes,
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
            (:T9, sensor_positions[:T9], "_T9", :circle, :red),
            (:T12, sensor_positions[:T9], "_T12", :utriangle, :red),
            (:T10, sensor_positions[:T10], "_T10", :diamond, :green),
            (:T11, sensor_positions[:T10], "_T11", :dtriangle, :green),
        )

        title_text = "1D_v11 $variant_name axial profile: $simulation_id $phase"
        subtitle_text = "flow=$(round(conditions[qlpm]; digits=2)) L/min, " *
                        "Io=$(round(conditions[Io] / 1000.0; digits=1)) kW/m^2"
        plot_object = plot(
            title=title_text * "\n" * subtitle_text,
            xlabel="Length (mm)",
            ylabel="Temperature (K)",
            legend=:outerright,
            grid=true,
            xlims=(0.0, 1000.0 * (L + REAR_TUBE_LENGTH_v11)),
            ylims=TEMPERATURE_AXIS_LIMITS_v11,
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
        scatter!(plot_object, [T3_SAMPLE_POSITION_v11 * 1000.0], [experiment[end, 4]];
                 label="T3 experiment", color=:red, markershape=:star5,
                 markerstrokewidth=0, ms=7)
        scatter!(plot_object, [CAVITY_LENGTH_v11 * 1000.0], [experiment[end, 5]];
                 label="T2 experiment", color=:black, markershape=:diamond,
                 markerstrokewidth=0, ms=5)
        scatter!(plot_object, [CAVITY_LENGTH_v11 * 1000.0], [model[end, 5]];
                 label="T2 model", color=:black, markershape=:xcross,
                 markerstrokewidth=1, ms=5)
        vline!(plot_object, [L * 1000.0, (L + REAR_TUBE_CAVITY_LENGTH_v11) * 1000.0];
               label=false, color=:gray, linestyle=:dot, lw=1)
        return plot_object
    end

    function write_cell_diagnostics_v11(path, simulation_id, params; is_cooling=false,
                                       nodes=RUNNER_PLOT_NODES_v11,
                                       variant_name="v11")
        model, result, experiment = solve_case_v11(
            params, simulation_id; is_cooling=is_cooling, nodes=nodes,
        )
        phase = is_cooling ? "cooling" : "heating"
        open(path, "w") do io
            println(io, "variant,simulation_id,phase,z_mm,Tsolid_K,Tgas_in_K,Tgas_out_K,Nu,h_W_m2K,UA_over_mcp,effectiveness,Qgas_W")
            for i in eachindex(result.z_solid)
                println(io, join((
                    variant_name,
                    simulation_id,
                    phase,
                    result.z_solid[i] * 1000.0,
                    result.solid_temperature[i, end],
                    result.gas_temperature[i, end],
                    result.gas_temperature[i + 1, end],
                    result.receiver_nusselt[i, end],
                    result.heat_transfer_coefficient[i, end],
                    result.receiver_ua_over_mcp[i, end],
                    result.receiver_effectiveness[i, end],
                    result.receiver_gas_heat[i, end],
                ), ','))
            end
        end
        return path
    end
end

begin # v11 Stage 1 study
    println("[run_1D_v11] Loaded $(length(sim_key_heat)) heating runs and " *
            "$(length(sim_key_cool)) cooling runs.")
    println("[run_1D_v11] Fixed source: eta_abs=$ETA_ABS_FIXED_v11, " *
            "eta_opt=$ETA_OPT_FIXED_v11, front_dep=$FRONT_DEPOSITION_STAGE1_v11.")

    variants_v11 = NamedTuple[
        (name="baseline_front_source", params=copy(pnew_v11), fitted=false),
    ]

    if !RUNNER_SKIP_FIT_v11
        println("[run_1D_v11] Fitting Stage 1 Nu shape on heating runs.")
        fit = calibrate_v11(
            nodes=RUNNER_FIT_NODES_v11,
            maximum_iterations=RUNNER_FIT_ITERATIONS_v11,
            maximum_time_seconds=RUNNER_FIT_SECONDS_v11,
        )
        push!(variants_v11, (
            name="stage1_nu_fit",
            params=fit.parameters,
            fitted=true,
        ))
        open(joinpath(RUNNER_OUTPUT_DIR_v11, "optimization_summary_1D_v11.txt"), "w") do io
            println(io, "objective=$(fit.minimum)")
            println(io, "return_code=$(fit.retcode)")
            println(io, "iterations=$(get(fit, :iterations, missing))")
            println(io, "parameters=$(fit.parameters)")
        end
    end

    representative_heat = ("E67", "E76", "E80")
    representative_cool = ("C80",)
    summary_rows = NamedTuple[]

    for variant in variants_v11
        println("[run_1D_v11] Running variant $(variant.name).")
        parameter_path = joinpath(
            RUNNER_OUTPUT_DIR_v11, "parameters_$(variant.name)_1D_v11.csv",
        )
        write_parameters_v11(parameter_path, variant.params)

        metrics = compute_metrics_v11(variant.params; nodes=RUNNER_PLOT_NODES_v11)
        metrics_path = joinpath(
            RUNNER_OUTPUT_DIR_v11, "analysis_results_$(variant.name)_1D_v11.csv",
        )
        write_metrics_v11(metrics_path, metrics)
        println("[run_1D_v11] Saved metrics: $metrics_path")
        for row in metrics
            push!(summary_rows, merge(row, (variant=variant.name,)))
        end

        steady_results = build_steady_results_v11(
            variant.params; nodes=RUNNER_PLOT_NODES_v11,
        )
        write_steady_results_v11(
            joinpath(RUNNER_OUTPUT_DIR_v11, "steady_results_$(variant.name)_1D_v11.csv"),
            steady_results,
        )
        steady_plot = steady_comparison_plot_v11(
            steady_results; title_suffix=variant.name,
        )
        save_runner_plot_v11(
            steady_plot, "steady_comparison_$(variant.name)_1D_v11.png",
        )

        for simulation_id in representative_heat
            profile_plot = axial_profile_plot_v11(
                simulation_id, variant.params; nodes=RUNNER_PLOT_NODES_v11,
                variant_name=variant.name,
            )
            save_runner_plot_v11(
                profile_plot,
                "axial_profile_$(simulation_id)_$(variant.name)_1D_v11.png",
            )
            write_cell_diagnostics_v11(
                joinpath(
                    RUNNER_DIAGNOSTIC_DIR_v11,
                    "cell_diagnostics_$(simulation_id)_$(variant.name)_1D_v11.csv",
                ),
                simulation_id,
                variant.params;
                nodes=RUNNER_PLOT_NODES_v11,
                variant_name=variant.name,
            )
        end

        for simulation_id in representative_cool
            profile_plot = axial_profile_plot_v11(
                simulation_id, variant.params; is_cooling=true,
                nodes=RUNNER_PLOT_NODES_v11,
                variant_name=variant.name,
            )
            save_runner_plot_v11(
                profile_plot,
                "axial_profile_$(simulation_id)_$(variant.name)_1D_v11.png",
            )
            write_cell_diagnostics_v11(
                joinpath(
                    RUNNER_DIAGNOSTIC_DIR_v11,
                    "cell_diagnostics_$(simulation_id)_$(variant.name)_1D_v11.csv",
                ),
                simulation_id,
                variant.params;
                is_cooling=true,
                nodes=RUNNER_PLOT_NODES_v11,
                variant_name=variant.name,
            )
        end
    end

    summary_path = joinpath(RUNNER_OUTPUT_DIR_v11, "analysis_results_all_variants_1D_v11.csv")
    open(summary_path, "w") do io
        println(io, "variant,simulation_id,phase,sensor,rmse_K,steady_error_K,t90_error_s,shape_loss")
        for row in summary_rows
            println(io, join((row.variant, row.simulation_id, row.phase, row.sensor,
                              row.rmse_K, row.steady_error_K,
                              row.t90_error_s, row.shape_loss), ','))
        end
    end
    println("[run_1D_v11] Saved combined variant metrics: $summary_path")

    println("[run_1D_v11] Complete.")
    println("[run_1D_v11] Output directory: $RUNNER_OUTPUT_DIR_v11")
end
