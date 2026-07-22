# ============================================================================
# run_1D_v14.jl - absorbed-power scale audit
# ============================================================================
# Diagnostic only. v14 keeps Nu and source distribution fixed, then asks:
#   1. What absorbed-power scale is needed per irradiance level?
#   2. If each heating case is solved separately, do required scales cluster by
#      irradiance or drift with flow?
#
# Optional environment controls:
#   RECEIVER1D_V14_RUNNER_PLOT_NODES=25
#   RECEIVER1D_V14_FIT_NODES=15
#   RECEIVER1D_V14_FIT_ITERATIONS=80
#   RECEIVER1D_V14_FIT_SECONDS=240
# ============================================================================

begin # libraries
    using Statistics
    using StatsPlots
end

include("1D_v14.jl")

begin # runner configuration
    runner_int_v14(name, default) = parse(Int, get(ENV, name, string(default)))
    runner_float_v14(name, default) = parse(Float64, get(ENV, name, string(default)))

    const RUNNER_PLOT_NODES_v14 = runner_int_v14(
        "RECEIVER1D_V14_RUNNER_PLOT_NODES", default_nodes,
    )
    const RUNNER_FIT_NODES_v14 = runner_int_v14("RECEIVER1D_V14_FIT_NODES", 15)
    const RUNNER_FIT_ITERATIONS_v14 = runner_int_v14(
        "RECEIVER1D_V14_FIT_ITERATIONS", 80,
    )
    const RUNNER_FIT_SECONDS_v14 = runner_float_v14(
        "RECEIVER1D_V14_FIT_SECONDS", 240.0,
    )

    const RUNNER_OUTPUT_DIR_v14 = joinpath(@__DIR__, "summaries", "1D_v14")
    const RUNNER_PLOT_DIR_v14 = joinpath(RUNNER_OUTPUT_DIR_v14, "plots")
    const RUNNER_DIAGNOSTIC_DIR_v14 = joinpath(RUNNER_OUTPUT_DIR_v14, "diagnostics")
    mkpath(RUNNER_OUTPUT_DIR_v14)
    mkpath(RUNNER_PLOT_DIR_v14)
    mkpath(RUNNER_DIAGNOSTIC_DIR_v14)
end

begin # runner helpers
    const RUNNER_SENSOR_COLORS_v14 = Dict(
        :T8 => :blue,
        :T9_pair => :red,
        :T10_pair => :green,
        :T3 => :purple,
        :T2 => :black,
    )

    function save_runner_plot_v14(plot_object, filename)
        path = joinpath(RUNNER_PLOT_DIR_v14, filename)
        savefig(plot_object, path)
        println("[run_1D_v14] Saved plot: $path")
        return path
    end

    absorbed_power_watts_v14(irradiance, scale) =
        ETA_ABS_FIXED_v14 * scale * irradiance * A_frt

    function write_parameters_v14(path, variant)
        parameters = variant.params
        names = (
            "A_Nu",
            "B_Re",
            "C_Pr",
            "L_entry_m",
            "beta_opt",
            "front_dep",
            "scale_456",
            "scale_304",
            "scale_256",
        )
        open(path, "w") do io
            println(io, "index,name,value")
            for i in eachindex(parameters)
                println(io, join((i, names[i], parameters[i]), ','))
            end
            weights = solar_weights_v5(parameters[5], parameters[6], default_nodes)
            println(io, "# fixed,Nu_fd_receiver,$NU_FD_RECEIVER_v14")
            println(io, "# fixed,entry_shape_exponent,$ENTRY_SHAPE_EXPONENT_FIXED_v14")
            println(io, "# fixed,eta_abs,$ETA_ABS_FIXED_v14")
            println(io, "# fixed,eta_opt,$ETA_OPT_FIXED_v14")
            println(io, "# fixed,beta_opt,$(parameters[5])")
            println(io, "# fixed,front_dep,$(parameters[6])")
            println(io, "# derived,front_cell_source_fraction,$(weights[1])")
            println(io, "# derived,downstream_source_fraction,$(1.0 - weights[1])")
            for (level_index, irradiance) in enumerate(IRRADIANCE_LEVELS_v14)
                scale = parameters[fit_power_scale_indices_v14[level_index]]
                nominal = absorbed_power_watts_v14(irradiance, 1.0)
                scaled = absorbed_power_watts_v14(irradiance, scale)
                println(io, "# power_audit,$irradiance,$scale,$nominal,$scaled")
            end
            println(io, "# fixed,h_front,$H_FRONT_FIXED_v14")
            println(io, "# fixed,eps_front,$EPS_FRONT_FIXED_v14")
            println(io, "# fixed,T3_sample_position_m,$T3_SAMPLE_POSITION_v14")
            println(io, "# fixed,cavity_heat_capacity_J_per_K,$CAVITY_HEAT_CAPACITY_v14")
            println(io, "# fixed,water_flange_temperature_K,$WATER_FLANGE_TEMPERATURE_v14")
        end
        return path
    end

    function write_metrics_v14(path, metrics)
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

    function build_steady_results_v14(params=pnew_v14; keys=sim_key_heat,
                                     nodes=RUNNER_PLOT_NODES_v14)
        results = NamedTuple[]
        weights = solar_weights_v5(params[5], params[6], nodes)
        for simulation_id in keys
            model, result, experiment = solve_case_v14(params, simulation_id; nodes=nodes)
            conditions = simulation_conditions[simulation_id]
            irradiance = conditions[Io]
            scale = absorbed_power_scale_v14(irradiance, params)
            nominal_absorbed = absorbed_power_watts_v14(irradiance, 1.0)
            scaled_absorbed = absorbed_power_watts_v14(irradiance, scale)
            push!(results, (
                simulation_id=simulation_id,
                flow_lpm=conditions[qlpm],
                irradiance=irradiance,
                power_scale=scale,
                nominal_absorbed_W=nominal_absorbed,
                scaled_absorbed_W=scaled_absorbed,
                added_absorbed_W=scaled_absorbed - nominal_absorbed,
                front_cell_source_fraction=weights[1],
                downstream_source_fraction=1.0 - weights[1],
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
                mean_nu=mean(result.receiver_nusselt[:, end]),
                mean_effectiveness=mean(result.receiver_effectiveness[:, end]),
                receiver_gas_heat=sum(result.receiver_gas_heat[:, end]),
                receiver_to_tube_heat=result.receiver_to_tube_heat[end],
                receiver_to_cavity_heat=result.receiver_to_cavity_heat[end],
                flange_heat_loss=result.flange_heat_loss[end],
            ))
        end
        return results
    end

    function write_steady_results_v14(path, steady_results)
        fields = propertynames(first(steady_results))
        open(path, "w") do io
            println(io, join(fields, ','))
            for row in steady_results
                println(io, join((getproperty(row, field) for field in fields), ','))
            end
        end
        return path
    end

    function slope_v14(xs, ys)
        xmean = mean(xs)
        ymean = mean(ys)
        denominator = sum((x .- xmean)^2 for x in xs)
        denominator <= eps(Float64) && return NaN
        return sum((xs[i] - xmean) * (ys[i] - ymean) for i in eachindex(xs)) /
               denominator
    end

    function write_flow_slopes_v14(path, steady_results)
        sensors = (:T8, :T9_pair, :T10_pair, :T3)
        irradiances = sort(unique(row.irradiance for row in steady_results); rev=true)
        open(path, "w") do io
            println(io, "irradiance,sensor,model_slope_K_per_Lmin,experiment_slope_K_per_Lmin")
            for irradiance in irradiances
                rows = [row for row in steady_results if row.irradiance == irradiance]
                flows = [row.flow_lpm for row in rows]
                for sensor in sensors
                    model_values = [getproperty(row, Symbol(sensor, "_model")) for row in rows]
                    experiment_values = [getproperty(row, Symbol(sensor, "_experiment")) for row in rows]
                    println(io, join((
                        irradiance,
                        sensor,
                        slope_v14(flows, model_values),
                        slope_v14(flows, experiment_values),
                    ), ','))
                end
            end
        end
        return path
    end

    function steady_comparison_plot_v14(steady_results; title_suffix="")
        sensors = (:T8, :T9_pair, :T10_pair, :T3, :T2)
        plot_object = plot(
            title="1D_v14 $title_suffix steady-state comparison",
            xlabel="Model temperature (K)",
            ylabel="Experiment temperature (K)",
            legend=:outerright,
            aspect_ratio=:equal,
            xlims=TEMPERATURE_AXIS_LIMITS_v14,
            ylims=TEMPERATURE_AXIS_LIMITS_v14,
        )
        for sensor in sensors
            model_values = [getproperty(row, Symbol(sensor, "_model")) for row in steady_results]
            experiment_values = [getproperty(row, Symbol(sensor, "_experiment")) for row in steady_results]
            scatter!(plot_object, model_values, experiment_values;
                     label=String(sensor), color=RUNNER_SENSOR_COLORS_v14[sensor],
                     ms=4, markerstrokewidth=0)
        end
        lower, upper = TEMPERATURE_AXIS_LIMITS_v14
        plot!(plot_object, [lower, upper], [lower, upper];
              label="1:1", color=:gray, linestyle=:dash)
        return plot_object
    end

    function axial_profile_plot_v14(simulation_id, params=pnew_v14;
                                   nodes=RUNNER_PLOT_NODES_v14,
                                   variant_name="v14")
        model, result, experiment = solve_case_v14(params, simulation_id; nodes=nodes)
        conditions = simulation_conditions[simulation_id]

        solid_x_mm = vcat(result.z_solid, result.z_rear_tube) .* 1000.0
        solid_T = vcat(result.solid_temperature[:, end],
                       result.rear_tube_temperature[:, end])
        gas_x_mm = result.z_gas .* 1000.0
        gas_T = result.gas_temperature[:, end]

        solid_measurement_points = (
            (:T8, sensor_positions[:T8], "_T8", :square, :blue),
            (:T9, sensor_positions[:T9], "_T9", :circle, :red),
            (:T12, sensor_positions[:T9], "_T12", :utriangle, :red),
            (:T10, sensor_positions[:T10], "_T10", :diamond, :green),
            (:T11, sensor_positions[:T10], "_T11", :dtriangle, :green),
        )

        scale = absorbed_power_scale_v14(conditions[Io], params)
        title_text = "1D_v14 $variant_name axial profile: $simulation_id"
        subtitle_text = "flow=$(round(conditions[qlpm]; digits=2)) L/min, " *
                        "Io=$(round(conditions[Io] / 1000.0; digits=1)) kW/m^2, " *
                        "scale=$(round(scale; digits=3))"
        plot_object = plot(
            title=title_text * "\n" * subtitle_text,
            xlabel="Length (mm)",
            ylabel="Temperature (K)",
            legend=:outerright,
            grid=true,
            xlims=(0.0, 1000.0 * (L + REAR_TUBE_LENGTH_v14)),
            ylims=TEMPERATURE_AXIS_LIMITS_v14,
        )
        plot!(plot_object, solid_x_mm, solid_T;
              label="solid/tube model", color=:blue, linestyle=:dash, lw=2)
        plot!(plot_object, gas_x_mm, gas_T;
              label="gas model", color=:green, lw=2)
        for (sensor, position, observation_id, marker, color) in solid_measurement_points
            scatter!(plot_object, [position * 1000.0],
                     [observation(measurements, simulation_id, observation_id)[end]];
                     label=String(sensor), color=color, markershape=marker,
                     markerstrokewidth=0, ms=5)
        end
        scatter!(plot_object, [T3_SAMPLE_POSITION_v14 * 1000.0], [experiment[end, 4]];
                 label="T3 experiment", color=:red, markershape=:star5,
                 markerstrokewidth=0, ms=7)
        scatter!(plot_object, [CAVITY_LENGTH_v14 * 1000.0], [experiment[end, 5]];
                 label="T2 experiment", color=:black, markershape=:diamond,
                 markerstrokewidth=0, ms=5)
        scatter!(plot_object, [CAVITY_LENGTH_v14 * 1000.0], [model[end, 5]];
                 label="T2 model", color=:black, markershape=:xcross,
                 markerstrokewidth=1, ms=5)
        vline!(plot_object, [L * 1000.0, (L + REAR_TUBE_CAVITY_LENGTH_v14) * 1000.0];
               label=false, color=:gray, linestyle=:dot, lw=1)
        return plot_object
    end

    function write_cell_diagnostics_v14(path, simulation_id, params;
                                       nodes=RUNNER_PLOT_NODES_v14,
                                       variant_name="v14")
        model, result, experiment = solve_case_v14(params, simulation_id; nodes=nodes)
        open(path, "w") do io
            println(io, "variant,simulation_id,z_mm,Tsolid_K,Tgas_in_K,Tgas_out_K,Nu,h_W_m2K,UA_over_mcp,effectiveness,Qgas_W")
            for i in eachindex(result.z_solid)
                println(io, join((
                    variant_name,
                    simulation_id,
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

    function fit_single_case_scale_v14(simulation_id; nodes=RUNNER_FIT_NODES_v14)
        conditions = simulation_conditions[simulation_id]
        scale_index = power_scale_index_v14(conditions[Io])
        scale_index > 0 || error("No power scale slot for $simulation_id")
        fit = optimize_parameter_subset_v14(
            p -> loss_heating_v14(p, [simulation_id]; nodes=nodes),
            pnew_v14,
            [scale_index],
            lb_full_v14,
            ub_full_v14;
            maximum_iterations=40,
            maximum_time_seconds=120.0,
            label="v14-case-$simulation_id",
        )
        model, result, experiment = solve_case_v14(
            fit.parameters, simulation_id; nodes=RUNNER_PLOT_NODES_v14,
        )
        scale = fit.parameters[scale_index]
        irradiance = conditions[Io]
        return (
            simulation_id=simulation_id,
            irradiance=irradiance,
            flow_lpm=conditions[qlpm],
            scale=scale,
            objective=fit.minimum,
            nominal_absorbed_W=absorbed_power_watts_v14(irradiance, 1.0),
            scaled_absorbed_W=absorbed_power_watts_v14(irradiance, scale),
            added_absorbed_W=absorbed_power_watts_v14(irradiance, scale) -
                             absorbed_power_watts_v14(irradiance, 1.0),
            T8_steady_error=model[end, 1] - experiment[end, 1],
            T9_pair_steady_error=model[end, 2] - experiment[end, 2],
            T10_pair_steady_error=model[end, 3] - experiment[end, 3],
            T3_steady_error=model[end, 4] - experiment[end, 4],
            T2_steady_error=model[end, 5] - experiment[end, 5],
            T9_minus_T8_error=(model[end, 2] - model[end, 1]) -
                              (experiment[end, 2] - experiment[end, 1]),
            T10_minus_T8_error=(model[end, 3] - model[end, 1]) -
                               (experiment[end, 3] - experiment[end, 1]),
        )
    end

    function write_namedtuple_rows_v14(path, rows)
        fields = propertynames(first(rows))
        open(path, "w") do io
            println(io, join(fields, ','))
            for row in rows
                println(io, join((getproperty(row, field) for field in fields), ','))
            end
        end
        return path
    end
end

begin # v14 power-scale audit
    println("[run_1D_v14] Loaded $(length(sim_key_heat)) heating runs and " *
            "$(length(sim_key_cool)) cooling runs.")
    println("[run_1D_v14] Fitting only per-irradiance absorbed-power scales.")

    fit = calibrate_v14(
        nodes=RUNNER_FIT_NODES_v14,
        maximum_iterations=RUNNER_FIT_ITERATIONS_v14,
        maximum_time_seconds=RUNNER_FIT_SECONDS_v14,
    )
    fitted_params = fit.parameters
    open(joinpath(RUNNER_OUTPUT_DIR_v14, "optimization_summary_1D_v14.txt"), "w") do io
        println(io, "objective=$(fit.minimum)")
        println(io, "return_code=$(fit.retcode)")
        println(io, "iterations=$(get(fit, :iterations, missing))")
        println(io, "parameters=$(fit.parameters)")
    end

    variants_v14 = NamedTuple[
        (name="baseline_scale1", params=copy(pnew_v14)),
        (name="per_irradiance_power_fit", params=fitted_params),
    ]

    summary_rows = NamedTuple[]
    representative_heat = ("E67", "E76", "E80")

    for variant in variants_v14
        println("[run_1D_v14] Running variant $(variant.name).")
        write_parameters_v14(
            joinpath(RUNNER_OUTPUT_DIR_v14, "parameters_$(variant.name)_1D_v14.csv"),
            variant,
        )
        metrics = compute_metrics_v14(variant.params; nodes=RUNNER_PLOT_NODES_v14)
        write_metrics_v14(
            joinpath(RUNNER_OUTPUT_DIR_v14, "analysis_results_$(variant.name)_1D_v14.csv"),
            metrics,
        )
        for row in metrics
            push!(summary_rows, merge(row, (variant=variant.name,)))
        end
        steady_results = build_steady_results_v14(
            variant.params; nodes=RUNNER_PLOT_NODES_v14,
        )
        write_steady_results_v14(
            joinpath(RUNNER_OUTPUT_DIR_v14, "steady_results_$(variant.name)_1D_v14.csv"),
            steady_results,
        )
        write_flow_slopes_v14(
            joinpath(RUNNER_OUTPUT_DIR_v14, "flow_slopes_$(variant.name)_1D_v14.csv"),
            steady_results,
        )
        save_runner_plot_v14(
            steady_comparison_plot_v14(steady_results; title_suffix=variant.name),
            "steady_comparison_$(variant.name)_1D_v14.png",
        )
        for simulation_id in representative_heat
            save_runner_plot_v14(
                axial_profile_plot_v14(
                    simulation_id, variant.params; variant_name=variant.name,
                ),
                "axial_profile_$(simulation_id)_$(variant.name)_1D_v14.png",
            )
            write_cell_diagnostics_v14(
                joinpath(
                    RUNNER_DIAGNOSTIC_DIR_v14,
                    "cell_diagnostics_$(simulation_id)_$(variant.name)_1D_v14.csv",
                ),
                simulation_id,
                variant.params;
                variant_name=variant.name,
            )
        end
    end

    summary_path = joinpath(RUNNER_OUTPUT_DIR_v14, "analysis_results_all_variants_1D_v14.csv")
    open(summary_path, "w") do io
        println(io, "variant,simulation_id,phase,sensor,rmse_K,steady_error_K,t90_error_s,shape_loss")
        for row in summary_rows
            println(io, join((row.variant, row.simulation_id, row.phase, row.sensor,
                              row.rmse_K, row.steady_error_K,
                              row.t90_error_s, row.shape_loss), ','))
        end
    end

    println("[run_1D_v14] Solving one power scale per heating case.")
    per_case_rows = [fit_single_case_scale_v14(simulation_id) for simulation_id in sim_key_heat]
    write_namedtuple_rows_v14(
        joinpath(RUNNER_OUTPUT_DIR_v14, "per_case_power_scale_audit_1D_v14.csv"),
        per_case_rows,
    )

    println("[run_1D_v14] Complete.")
    println("[run_1D_v14] Output directory: $RUNNER_OUTPUT_DIR_v14")
end
