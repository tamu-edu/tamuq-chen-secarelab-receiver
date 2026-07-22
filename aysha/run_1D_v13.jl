# ============================================================================
# run_1D_v13.jl - v13 axial-radiation and thermocouple diagnostics
# ============================================================================
# Forward comparison only. v13 keeps the v12a source distribution and v11 Nu
# shape fixed, then compares:
#   - baseline_v12a_source
#   - axial_view_rad_fixed
#   - axial_view_rad_upper
#   - tc_wire_fixed
#   - combined_fixed
#
# Optional environment controls:
#   RECEIVER1D_V13_RUNNER_PLOT_NODES=25
# ============================================================================

begin # libraries
    using Statistics
    using StatsPlots
end

include("1D_v13.jl")

begin # runner configuration
    runner_int_v13(name, default) = parse(Int, get(ENV, name, string(default)))

    const RUNNER_PLOT_NODES_v13 = runner_int_v13(
        "RECEIVER1D_V13_RUNNER_PLOT_NODES", default_nodes,
    )

    const RUNNER_OUTPUT_DIR_v13 = joinpath(@__DIR__, "summaries", "1D_v13")
    const RUNNER_PLOT_DIR_v13 = joinpath(RUNNER_OUTPUT_DIR_v13, "plots")
    const RUNNER_DIAGNOSTIC_DIR_v13 = joinpath(RUNNER_OUTPUT_DIR_v13, "diagnostics")
    mkpath(RUNNER_OUTPUT_DIR_v13)
    mkpath(RUNNER_PLOT_DIR_v13)
    mkpath(RUNNER_DIAGNOSTIC_DIR_v13)
end

begin # runner helpers
    const RUNNER_SENSOR_COLORS_v13 = Dict(
        :T8 => :blue,
        :T9_pair => :red,
        :T10_pair => :green,
        :T3 => :purple,
        :T2 => :black,
    )

    function save_runner_plot_v13(plot_object, filename)
        path = joinpath(RUNNER_PLOT_DIR_v13, filename)
        savefig(plot_object, path)
        println("[run_1D_v13] Saved plot: $path")
        return path
    end

    function variant_kwargs_v13(variant)
        return (
            use_axial_view_radiation=variant.use_axial_view_radiation,
            axial_view_strength=variant.axial_view_strength,
            use_tc_measurement_model=variant.use_tc_measurement_model,
        )
    end

    function solve_variant_case_v13(variant, simulation_id; is_cooling=false,
                                    nodes=RUNNER_PLOT_NODES_v13)
        return solve_case_v13(
            variant.params, simulation_id; is_cooling=is_cooling, nodes=nodes,
            variant_kwargs_v13(variant)...,
        )
    end

    function write_parameters_v13(path, variant)
        parameters = variant.params
        names = (
            "A_Nu",
            "B_Re",
            "C_Pr",
            "L_entry_m",
            "beta_opt",
            "front_dep",
        )
        open(path, "w") do io
            println(io, "index,name,value")
            for i in eachindex(parameters)
                println(io, join((i, names[i], parameters[i]), ','))
            end
            weights = solar_weights_v5(parameters[5], parameters[6], default_nodes)
            println(io, "# fixed,Nu_fd_receiver,$NU_FD_RECEIVER_v13")
            println(io, "# fixed,entry_shape_exponent,$ENTRY_SHAPE_EXPONENT_FIXED_v13")
            println(io, "# fixed,eta_abs,$ETA_ABS_FIXED_v13")
            println(io, "# fixed,eta_opt,$ETA_OPT_FIXED_v13")
            println(io, "# fixed,beta_opt,$(parameters[5])")
            println(io, "# fixed,front_dep,$(parameters[6])")
            println(io, "# derived,front_cell_source_fraction,$(weights[1])")
            println(io, "# derived,downstream_source_fraction,$(1.0 - weights[1])")
            println(io, "# fixed,axial_view_enabled,$(variant.use_axial_view_radiation)")
            println(io, "# fixed,axial_view_strength,$(variant.axial_view_strength)")
            println(io, "# fixed,axial_view_emissivity,$AXIAL_VIEW_RADIATION_EMISSIVITY_v13")
            println(io, "# fixed,axial_view_row_sum,$AXIAL_VIEW_RADIATION_ROW_SUM_v13")
            println(io, "# fixed,axial_view_length_m,$AXIAL_VIEW_RADIATION_LENGTH_v13")
            println(io, "# fixed,tc_model_enabled,$(variant.use_tc_measurement_model)")
            println(io, "# fixed,tc_wire_loss_front_fraction,$TC_WIRE_LOSS_FRONT_FRACTION_v13")
            println(io, "# fixed,tc_wire_loss_length_m,$TC_WIRE_LOSS_LENGTH_v13")
            println(io, "# fixed,tc_wire_sink_temperature_K,$TC_WIRE_SINK_TEMPERATURE_v13")
            println(io, "# fixed,h_front,$H_FRONT_FIXED_v13")
            println(io, "# fixed,eps_front,$EPS_FRONT_FIXED_v13")
            println(io, "# fixed,ell_rad_m,$DEFAULT_ROSSELAND_MEAN_PATH_v13")
            println(io, "# fixed,T3_sample_position_m,$T3_SAMPLE_POSITION_v13")
            println(io, "# fixed,cavity_heat_capacity_J_per_K,$CAVITY_HEAT_CAPACITY_v13")
            println(io, "# fixed,water_flange_temperature_K,$WATER_FLANGE_TEMPERATURE_v13")
        end
        return path
    end

    function write_metrics_v13(path, metrics)
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

    function build_steady_results_v13(variant; keys=sim_key_heat,
                                     nodes=RUNNER_PLOT_NODES_v13)
        results = NamedTuple[]
        weights = solar_weights_v5(variant.params[5], variant.params[6], nodes)
        for simulation_id in keys
            model, result, experiment = solve_variant_case_v13(
                variant, simulation_id; nodes=nodes,
            )
            conditions = simulation_conditions[simulation_id]
            raw_T8 = solid_at_v3(result, sensor_positions[:T8])[end]
            raw_T9 = solid_at_v3(result, sensor_positions[:T9])[end]
            raw_T10 = solid_at_v3(result, sensor_positions[:T10])[end]
            qview = result.axial_view_radiation_heat[:, end]
            push!(results, (
                simulation_id=simulation_id,
                flow_lpm=conditions[qlpm],
                irradiance=conditions[Io],
                beta_opt=variant.params[5],
                front_dep=variant.params[6],
                front_cell_source_fraction=weights[1],
                downstream_source_fraction=1.0 - weights[1],
                axial_view_strength=variant.axial_view_strength,
                tc_model=variant.use_tc_measurement_model,
                T8_model=model[end, 1],
                T8_raw_solid=raw_T8,
                T8_tc_shift=model[end, 1] - raw_T8,
                T8_experiment=experiment[end, 1],
                T9_pair_model=model[end, 2],
                T9_pair_raw_solid=raw_T9,
                T9_pair_tc_shift=model[end, 2] - raw_T9,
                T9_pair_experiment=experiment[end, 2],
                T10_pair_model=model[end, 3],
                T10_pair_raw_solid=raw_T10,
                T10_pair_tc_shift=model[end, 3] - raw_T10,
                T10_pair_experiment=experiment[end, 3],
                T3_model=model[end, 4],
                T3_experiment=experiment[end, 4],
                T2_model=model[end, 5],
                T2_experiment=experiment[end, 5],
                T9_minus_T8_model=model[end, 2] - model[end, 1],
                T9_minus_T8_experiment=experiment[end, 2] - experiment[end, 1],
                T10_minus_T8_model=model[end, 3] - model[end, 1],
                T10_minus_T8_experiment=experiment[end, 3] - experiment[end, 1],
                axial_view_abs_heat=sum(abs.(qview)),
                axial_view_front_loss=-min(qview[1], 0.0),
                axial_view_max_gain=maximum(qview),
                axial_view_max_loss=minimum(qview),
                mean_nu=mean(result.receiver_nusselt[:, end]),
                mean_effectiveness=mean(result.receiver_effectiveness[:, end]),
                receiver_to_tube_heat=result.receiver_to_tube_heat[end],
                receiver_to_cavity_heat=result.receiver_to_cavity_heat[end],
                flange_heat_loss=result.flange_heat_loss[end],
            ))
        end
        return results
    end

    function write_steady_results_v13(path, steady_results)
        fields = propertynames(first(steady_results))
        open(path, "w") do io
            println(io, join(fields, ','))
            for row in steady_results
                println(io, join((getproperty(row, field) for field in fields), ','))
            end
        end
        return path
    end

    function steady_comparison_plot_v13(steady_results; title_suffix="")
        sensors = (:T8, :T9_pair, :T10_pair, :T3, :T2)
        plot_object = plot(
            title="1D_v13 $title_suffix steady-state comparison",
            xlabel="Model temperature (K)",
            ylabel="Experiment temperature (K)",
            legend=:outerright,
            aspect_ratio=:equal,
            xlims=TEMPERATURE_AXIS_LIMITS_v13,
            ylims=TEMPERATURE_AXIS_LIMITS_v13,
        )
        for sensor in sensors
            model_values = [getproperty(row, Symbol(sensor, "_model")) for row in steady_results]
            experiment_values = [getproperty(row, Symbol(sensor, "_experiment")) for row in steady_results]
            scatter!(plot_object, model_values, experiment_values;
                     label=String(sensor), color=RUNNER_SENSOR_COLORS_v13[sensor],
                     ms=4, markerstrokewidth=0)
        end
        lower, upper = TEMPERATURE_AXIS_LIMITS_v13
        plot!(plot_object, [lower, upper], [lower, upper];
              label="1:1", color=:gray, linestyle=:dash)
        return plot_object
    end

    function axial_profile_plot_v13(simulation_id, variant;
                                   is_cooling=false,
                                   nodes=RUNNER_PLOT_NODES_v13)
        model, result, experiment = solve_variant_case_v13(
            variant, simulation_id; is_cooling=is_cooling, nodes=nodes,
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

        title_text = "1D_v13 $(variant.name): $simulation_id $phase"
        subtitle_text = "flow=$(round(conditions[qlpm]; digits=2)) L/min, " *
                        "Io=$(round(conditions[Io] / 1000.0; digits=1)) kW/m^2"
        plot_object = plot(
            title=title_text * "\n" * subtitle_text,
            xlabel="Length (mm)",
            ylabel="Temperature (K)",
            legend=:outerright,
            grid=true,
            xlims=(0.0, 1000.0 * (L + REAR_TUBE_LENGTH_v13)),
            ylims=TEMPERATURE_AXIS_LIMITS_v13,
        )
        plot!(plot_object, solid_x_mm, solid_T;
              label="bulk solid/tube", color=:blue, linestyle=:dash, lw=2)
        plot!(plot_object, gas_x_mm, gas_T;
              label="gas model", color=:green, lw=2)
        scatter!(plot_object, [sensor_positions[:T8] * 1000.0],
                 [model[end, 1]]; label="T8 model compare",
                 color=:blue, markershape=:xcross, ms=5)
        scatter!(plot_object, [sensor_positions[:T9] * 1000.0],
                 [model[end, 2]]; label="T9 pair model compare",
                 color=:red, markershape=:xcross, ms=5)
        scatter!(plot_object, [sensor_positions[:T10] * 1000.0],
                 [model[end, 3]]; label="T10 pair model compare",
                 color=:green, markershape=:xcross, ms=5)
        for (sensor, position, observation_id, marker, color) in solid_measurement_points
            scatter!(plot_object, [position * 1000.0],
                     [observation(measurement_data, simulation_id, observation_id)[end]];
                     label=String(sensor), color=color, markershape=marker,
                     markerstrokewidth=0, ms=5)
        end
        scatter!(plot_object, [T3_SAMPLE_POSITION_v13 * 1000.0], [experiment[end, 4]];
                 label="T3 experiment", color=:red, markershape=:star5,
                 markerstrokewidth=0, ms=7)
        scatter!(plot_object, [CAVITY_LENGTH_v13 * 1000.0], [experiment[end, 5]];
                 label="T2 experiment", color=:black, markershape=:diamond,
                 markerstrokewidth=0, ms=5)
        scatter!(plot_object, [CAVITY_LENGTH_v13 * 1000.0], [model[end, 5]];
                 label="T2 model", color=:black, markershape=:xcross,
                 markerstrokewidth=1, ms=5)
        vline!(plot_object, [L * 1000.0, (L + REAR_TUBE_CAVITY_LENGTH_v13) * 1000.0];
               label=false, color=:gray, linestyle=:dot, lw=1)
        return plot_object
    end

    function write_cell_diagnostics_v13(path, simulation_id, variant;
                                       is_cooling=false,
                                       nodes=RUNNER_PLOT_NODES_v13)
        model, result, experiment = solve_variant_case_v13(
            variant, simulation_id; is_cooling=is_cooling, nodes=nodes,
        )
        phase = is_cooling ? "cooling" : "heating"
        open(path, "w") do io
            println(io, "variant,simulation_id,phase,z_mm,Tsolid_K,Tgas_in_K,Tgas_out_K,Nu,h_W_m2K,UA_over_mcp,effectiveness,Qgas_W,Qaxial_view_W")
            for i in eachindex(result.z_solid)
                println(io, join((
                    variant.name,
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
                    result.axial_view_radiation_heat[i, end],
                ), ','))
            end
        end
        return path
    end
end

begin # v13 forward diagnostic matrix
    println("[run_1D_v13] Loaded $(length(sim_key_heat)) heating runs and " *
            "$(length(sim_key_cool)) cooling runs.")
    println("[run_1D_v13] No fitted parameters. Comparing fixed radiation/TC variants.")

    variants_v13 = NamedTuple[
        (
            name="baseline_v12a_source",
            params=copy(pnew_v13),
            use_axial_view_radiation=false,
            axial_view_strength=0.0,
            use_tc_measurement_model=false,
        ),
        (
            name="axial_view_rad_fixed",
            params=copy(pnew_v13),
            use_axial_view_radiation=true,
            axial_view_strength=1.0,
            use_tc_measurement_model=false,
        ),
        (
            name="axial_view_rad_upper",
            params=copy(pnew_v13),
            use_axial_view_radiation=true,
            axial_view_strength=3.0,
            use_tc_measurement_model=false,
        ),
        (
            name="tc_wire_fixed",
            params=copy(pnew_v13),
            use_axial_view_radiation=false,
            axial_view_strength=0.0,
            use_tc_measurement_model=true,
        ),
        (
            name="combined_fixed",
            params=copy(pnew_v13),
            use_axial_view_radiation=true,
            axial_view_strength=1.0,
            use_tc_measurement_model=true,
        ),
    ]

    representative_heat = ("E67", "E76", "E80")
    representative_cool = ("C80",)
    summary_rows = NamedTuple[]

    for variant in variants_v13
        println("[run_1D_v13] Running variant $(variant.name).")
        write_parameters_v13(
            joinpath(RUNNER_OUTPUT_DIR_v13, "parameters_$(variant.name)_1D_v13.csv"),
            variant,
        )

        metrics = compute_metrics_v13(
            variant.params; nodes=RUNNER_PLOT_NODES_v13,
            variant_kwargs_v13(variant)...,
        )
        metrics_path = joinpath(
            RUNNER_OUTPUT_DIR_v13, "analysis_results_$(variant.name)_1D_v13.csv",
        )
        write_metrics_v13(metrics_path, metrics)
        println("[run_1D_v13] Saved metrics: $metrics_path")
        for row in metrics
            push!(summary_rows, merge(row, (variant=variant.name,)))
        end

        steady_results = build_steady_results_v13(
            variant; nodes=RUNNER_PLOT_NODES_v13,
        )
        write_steady_results_v13(
            joinpath(RUNNER_OUTPUT_DIR_v13, "steady_results_$(variant.name)_1D_v13.csv"),
            steady_results,
        )
        save_runner_plot_v13(
            steady_comparison_plot_v13(steady_results; title_suffix=variant.name),
            "steady_comparison_$(variant.name)_1D_v13.png",
        )

        for simulation_id in representative_heat
            save_runner_plot_v13(
                axial_profile_plot_v13(
                    simulation_id, variant; nodes=RUNNER_PLOT_NODES_v13,
                ),
                "axial_profile_$(simulation_id)_$(variant.name)_1D_v13.png",
            )
            write_cell_diagnostics_v13(
                joinpath(
                    RUNNER_DIAGNOSTIC_DIR_v13,
                    "cell_diagnostics_$(simulation_id)_$(variant.name)_1D_v13.csv",
                ),
                simulation_id,
                variant;
                nodes=RUNNER_PLOT_NODES_v13,
            )
        end

        for simulation_id in representative_cool
            save_runner_plot_v13(
                axial_profile_plot_v13(
                    simulation_id, variant; is_cooling=true,
                    nodes=RUNNER_PLOT_NODES_v13,
                ),
                "axial_profile_$(simulation_id)_$(variant.name)_1D_v13.png",
            )
            write_cell_diagnostics_v13(
                joinpath(
                    RUNNER_DIAGNOSTIC_DIR_v13,
                    "cell_diagnostics_$(simulation_id)_$(variant.name)_1D_v13.csv",
                ),
                simulation_id,
                variant;
                is_cooling=true,
                nodes=RUNNER_PLOT_NODES_v13,
            )
        end
    end

    summary_path = joinpath(RUNNER_OUTPUT_DIR_v13, "analysis_results_all_variants_1D_v13.csv")
    open(summary_path, "w") do io
        println(io, "variant,simulation_id,phase,sensor,rmse_K,steady_error_K,t90_error_s,shape_loss")
        for row in summary_rows
            println(io, join((row.variant, row.simulation_id, row.phase, row.sensor,
                              row.rmse_K, row.steady_error_K,
                              row.t90_error_s, row.shape_loss), ','))
        end
    end
    println("[run_1D_v13] Saved combined variant metrics: $summary_path")

    println("[run_1D_v13] Complete.")
    println("[run_1D_v13] Output directory: $RUNNER_OUTPUT_DIR_v13")
end
