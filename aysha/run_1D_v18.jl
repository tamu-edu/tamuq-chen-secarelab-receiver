# ============================================================================
# run_1D_v18.jl - T9/T10 gas-bias measurement diagnostic fit
# ============================================================================
# Diagnostic/refit runner. v18 keeps the v14 per-irradiance absorbed-power
# scales fixed, keeps the fitted v17 thermal topology fixed, and fits only a
# constrained T9/T10 gas-bias measurement layer.
#
# Optional environment controls:
#   RECEIVER1D_v18_RUNNER_PLOT_NODES=25
#   RECEIVER1D_v18_FIT_NODES=15
#   RECEIVER1D_v18_FIT_ITERATIONS=120
#   RECEIVER1D_v18_FIT_SECONDS=360
# ============================================================================

begin # libraries
    using Statistics
    using StatsPlots
end

include("1D_v18.jl")

begin # runner configuration
    runner_int_v18(name, default) = parse(Int, get(ENV, name, string(default)))
    runner_float_v18(name, default) = parse(Float64, get(ENV, name, string(default)))

    const RUNNER_PLOT_NODES_v18 = runner_int_v18(
        "RECEIVER1D_v18_RUNNER_PLOT_NODES", default_nodes,
    )
    const RUNNER_FIT_NODES_v18 = runner_int_v18("RECEIVER1D_v18_FIT_NODES", 15)
    const RUNNER_FIT_ITERATIONS_v18 = runner_int_v18(
        "RECEIVER1D_v18_FIT_ITERATIONS", 120,
    )
    const RUNNER_FIT_SECONDS_v18 = runner_float_v18(
        "RECEIVER1D_v18_FIT_SECONDS", 360.0,
    )

    const RUNNER_OUTPUT_DIR_v18 = joinpath(@__DIR__, "summaries", "1D_v18")
    const RUNNER_PLOT_DIR_v18 = joinpath(RUNNER_OUTPUT_DIR_v18, "plots")
    const RUNNER_TRANSIENT_DIR_v18 = joinpath(RUNNER_PLOT_DIR_v18, "transients")
    const RUNNER_AXIAL_DIR_v18 = joinpath(RUNNER_PLOT_DIR_v18, "axial_profiles")
    const RUNNER_DIAGNOSTIC_DIR_v18 = joinpath(RUNNER_OUTPUT_DIR_v18, "diagnostics")
    mkpath(RUNNER_OUTPUT_DIR_v18)
    mkpath(RUNNER_PLOT_DIR_v18)
    mkpath(RUNNER_TRANSIENT_DIR_v18)
    mkpath(RUNNER_AXIAL_DIR_v18)
    mkpath(RUNNER_DIAGNOSTIC_DIR_v18)
end

begin # runner helpers
    const RUNNER_SENSOR_COLORS_v18 = Dict(
        :T8 => :blue,
        :T12_wall => :red,
        :T11_wall => :green,
        :T9_active => :orange,
        :T10_active => :brown,
        :T3 => :purple,
        :T2 => :black,
    )

    function save_runner_plot_v18(plot_object, filename)
        path = joinpath(RUNNER_PLOT_DIR_v18, filename)
        savefig(plot_object, path)
        println("[run_1D_v18] Saved plot: $path")
        return path
    end

    absorbed_power_watts_v18(irradiance, scale) =
        ETA_ABS_FIXED_v18 * scale * irradiance * A_frt

    function write_parameters_v18(path, variant)
        parameters = variant.params
        names = (
            "A_Nu_app",
            "B_Re",
            "C_Pr_fixed",
            "reserved",
            "beta_opt",
            "front_dep",
            "scale_456",
            "scale_304",
            "scale_256",
            "f_side",
            "G_core_side_W_m_K",
            "C_side_eff_J_K",
            "k_side_ref_W_m_K",
            "alpha9_max",
            "alpha10_max",
            "Re50_alpha",
            "alpha_m",
        )
        open(path, "w") do io
            println(io, "index,name,value")
            for i in eachindex(parameters)
                println(io, join((i, names[i], parameters[i]), ','))
            end
            weights = solar_weights_v5(parameters[5], parameters[6], default_nodes)
            println(io, "# law,Nu_app,A*Re^B*Pr^(1/3)")
            println(io, "# fixed,Pr_exponent,$PRANDTL_EXPONENT_FIXED_v18")
            println(io, "# fixed,standard_air_density_kg_m3,$STANDARD_AIR_DENSITY_v18")
            println(io, "# fixed,eta_abs,$ETA_ABS_FIXED_v18")
            println(io, "# fixed,eta_opt,$ETA_OPT_FIXED_v18")
            println(io, "# fixed,beta_opt,$(parameters[5])")
            println(io, "# fixed,front_dep,$(parameters[6])")
            println(io, "# fitted,f_side,$(side_source_fraction_v18(parameters))")
            println(io, "# fitted,G_core_side_W_m_K,$(core_side_conductance_per_length_v18(parameters))")
            println(io, "# fitted,C_side_eff_J_K,$(side_heat_capacity_total_v18(parameters))")
            println(io, "# fitted,k_side_ref_W_m_K,$(side_axial_participation_conductivity_v18(900.0, parameters))")
            println(io, "# fitted,alpha9_max,$(parameters[14])")
            println(io, "# fitted,alpha10_max,$(parameters[15])")
            println(io, "# fitted,Re50_alpha,$(parameters[16])")
            println(io, "# fitted,alpha_m,$(parameters[17])")
            println(io, "# fixed_thermal,true")
            println(io, "# derived,C_active_eff_J_K,$(active_receiver_heat_capacity_v18(default_nodes))")
            println(io, "# derived,C_participating_eff_J_K,$(participating_assembly_heat_capacity_v18(parameters))")
            println(io, "# reference,measured_C_eff_J_K,$MEASURED_ASSEMBLY_CAPACITANCE_v18")
            println(io, "# derived,front_cell_source_fraction,$(weights[1])")
            println(io, "# derived,downstream_source_fraction,$(1.0 - weights[1])")
            for (level_index, irradiance) in enumerate(IRRADIANCE_LEVELS_v18)
                scale = parameters[power_scale_index_v18(irradiance)]
                nominal = absorbed_power_watts_v18(irradiance, 1.0)
                scaled = absorbed_power_watts_v18(irradiance, scale)
                println(io, "# fixed_power_scale,$irradiance,$scale,$nominal,$scaled")
            end
            println(io, "# fixed,h_front,$H_FRONT_FIXED_v18")
            println(io, "# fixed,eps_front,$EPS_FRONT_FIXED_v18")
            println(io, "# fixed,T3_sample_position_m,$T3_SAMPLE_POSITION_v18")
            println(io, "# fixed,cavity_heat_capacity_J_per_K,$CAVITY_HEAT_CAPACITY_v18")
            println(io, "# fixed,water_flange_temperature_K,$WATER_FLANGE_TEMPERATURE_v18")
        end
        return path
    end

    function write_metrics_v18(path, metrics)
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

    function build_steady_results_v18(params=pnew_v18; keys=sim_key_heat,
                                     nodes=RUNNER_PLOT_NODES_v18)
        results = NamedTuple[]
        weights = solar_weights_v5(params[5], params[6], nodes)
        for simulation_id in keys
            model, result, experiment = solve_case_v18(params, simulation_id; nodes=nodes)
            conditions = simulation_conditions[simulation_id]
            irradiance = conditions[Io]
            scale = absorbed_power_scale_v18(irradiance, params)
            nominal_absorbed = absorbed_power_watts_v18(irradiance, 1.0)
            scaled_absorbed = absorbed_power_watts_v18(irradiance, scale)
            T9_core = active_at_v18(result, sensor_positions[:T9])[end]
            T10_core = active_at_v18(result, sensor_positions[:T10])[end]
            T9_gas = gas_at_position_v18(result, sensor_positions[:T9])[end]
            T10_gas = gas_at_position_v18(result, sensor_positions[:T10])[end]
            Re9 = reynolds_at_v18(result, sensor_positions[:T9])[end]
            Re10 = reynolds_at_v18(result, sensor_positions[:T10])[end]
            alpha9 = alpha9_measurement_v18(Re9, params)
            alpha10 = alpha10_measurement_v18(Re10, params)
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
                T12_wall_model=model[end, 2],
                T12_wall_experiment=experiment[end, 2],
                T11_wall_model=model[end, 3],
                T11_wall_experiment=experiment[end, 3],
                T9_active_model=model[end, 4],
                T9_active_experiment=experiment[end, 4],
                T10_active_model=model[end, 5],
                T10_active_experiment=experiment[end, 5],
                T9_core_model=T9_core,
                T10_core_model=T10_core,
                T9_gas_model=T9_gas,
                T10_gas_model=T10_gas,
                T9_alpha=alpha9,
                T10_alpha=alpha10,
                T9_reynolds=Re9,
                T10_reynolds=Re10,
                T3_model=model[end, 6],
                T3_experiment=experiment[end, 6],
                T2_model=model[end, 7],
                T2_experiment=experiment[end, 7],
                T12_minus_T8_model=model[end, 2] - model[end, 1],
                T12_minus_T8_experiment=experiment[end, 2] - experiment[end, 1],
                T11_minus_T8_model=model[end, 3] - model[end, 1],
                T11_minus_T8_experiment=experiment[end, 3] - experiment[end, 1],
                T12_minus_T9_model=model[end, 2] - model[end, 4],
                T12_minus_T9_experiment=experiment[end, 2] - experiment[end, 4],
                T11_minus_T10_model=model[end, 3] - model[end, 5],
                T11_minus_T10_experiment=experiment[end, 3] - experiment[end, 5],
                Lambda58_experiment=(experiment[end, 2] -
                                     experiment[end, 4]) /
                                    max(experiment[end, 2] -
                                        observation(measurements, simulation_id, "_Tamb")[end],
                                        eps(Float64)),
                Lambda107_experiment=(experiment[end, 3] -
                                      experiment[end, 5]) /
                                     max(experiment[end, 3] -
                                         observation(measurements, simulation_id, "_Tamb")[end],
                                         eps(Float64)),
                mean_nu=mean(result.receiver_nusselt[:, end]),
                mean_effectiveness=mean(result.receiver_effectiveness[:, end]),
                receiver_gas_heat=sum(result.receiver_gas_heat[:, end]),
                receiver_to_tube_heat=result.receiver_to_tube_heat[end],
                receiver_to_cavity_heat=result.receiver_to_cavity_heat[end],
                active_wall_heat=result.active_wall_heat[end],
                flange_heat_loss=result.flange_heat_loss[end],
                f_side=side_source_fraction_v18(params),
                G_core_side_W_m_K=core_side_conductance_per_length_v18(params),
                C_side_eff_J_K=side_heat_capacity_total_v18(params),
                C_participating_eff_J_K=participating_assembly_heat_capacity_v18(params, nodes),
                k_side_ref_W_m_K=side_axial_participation_conductivity_v18(900.0, params),
            ))
        end
        return results
    end

    function write_steady_results_v18(path, steady_results)
        fields = propertynames(first(steady_results))
        open(path, "w") do io
            println(io, join(fields, ','))
            for row in steady_results
                println(io, join((getproperty(row, field) for field in fields), ','))
            end
        end
        return path
    end

    function slope_v18(xs, ys)
        xmean = mean(xs)
        ymean = mean(ys)
        denominator = sum((x .- xmean)^2 for x in xs)
        denominator <= eps(Float64) && return NaN
        return sum((xs[i] - xmean) * (ys[i] - ymean) for i in eachindex(xs)) /
               denominator
    end

    function write_flow_slopes_v18(path, steady_results)
        sensors = (:T8, :T12_wall, :T11_wall, :T9_active, :T10_active, :T3)
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
                        slope_v18(flows, model_values),
                        slope_v18(flows, experiment_values),
                    ), ','))
                end
            end
        end
        return path
    end

    function steady_comparison_plot_v18(steady_results; title_suffix="")
        sensors = (:T8, :T12_wall, :T11_wall, :T9_active, :T10_active, :T3, :T2)
        plot_object = plot(
            title="1D_v18 $title_suffix steady-state comparison",
            xlabel="Model temperature (K)",
            ylabel="Experiment temperature (K)",
            legend=:outerright,
            aspect_ratio=:equal,
            xlims=TEMPERATURE_AXIS_LIMITS_v18,
            ylims=TEMPERATURE_AXIS_LIMITS_v18,
        )
        for sensor in sensors
            model_values = [getproperty(row, Symbol(sensor, "_model")) for row in steady_results]
            experiment_values = [getproperty(row, Symbol(sensor, "_experiment")) for row in steady_results]
            scatter!(plot_object, model_values, experiment_values;
                     label=String(sensor), color=RUNNER_SENSOR_COLORS_v18[sensor],
                     ms=4, markerstrokewidth=0)
        end
        lower, upper = TEMPERATURE_AXIS_LIMITS_v18
        plot!(plot_object, [lower, upper], [lower, upper];
              label="1:1", color=:gray, linestyle=:dash)
        return plot_object
    end

    function axial_profile_plot_v18(simulation_id, params=pnew_v18;
                                   nodes=RUNNER_PLOT_NODES_v18,
                                   variant_name="v18")
        model, result, experiment = solve_case_v18(params, simulation_id; nodes=nodes)
        conditions = simulation_conditions[simulation_id]

        wall_x_mm = vcat(result.z_solid, result.z_rear_tube) .* 1000.0
        wall_T = vcat(result.wall_temperature[:, end],
                      result.rear_tube_temperature[:, end])
        active_x_mm = result.z_solid .* 1000.0
        active_T = result.active_temperature[:, end]
        gas_x_mm = result.z_gas .* 1000.0
        gas_T = result.gas_temperature[:, end]

        solid_measurement_points = (
            (:T8, sensor_positions[:T8], "_T8", :square, :blue),
            (:T9, sensor_positions[:T9], "_T9", :circle, :red),
            (:T12, sensor_positions[:T9], "_T12", :utriangle, :red),
            (:T10, sensor_positions[:T10], "_T10", :diamond, :green),
            (:T11, sensor_positions[:T10], "_T11", :dtriangle, :green),
        )

        scale = absorbed_power_scale_v18(conditions[Io], params)
        title_text = "1D_v18 $variant_name axial profile: $simulation_id"
        subtitle_text = "flow=$(round(conditions[qlpm]; digits=2)) L/min, " *
                        "Io=$(round(conditions[Io] / 1000.0; digits=1)) kW/m^2, " *
                        "scale=$(round(scale; digits=3))"
        plot_object = plot(
            title=title_text * "\n" * subtitle_text,
            xlabel="Length (mm)",
            ylabel="Temperature (K)",
            legend=:outerright,
            grid=true,
            xlims=(0.0, 1000.0 * (L + REAR_TUBE_LENGTH_v18)),
            ylims=TEMPERATURE_AXIS_LIMITS_v18,
        )
        plot!(plot_object, wall_x_mm, wall_T;
              label="wall/tube model", color=:blue, linestyle=:dash, lw=2)
        plot!(plot_object, active_x_mm, active_T;
              label="active solid model", color=:orange, linestyle=:dot, lw=2)
        plot!(plot_object, gas_x_mm, gas_T;
              label="gas model", color=:green, lw=2)
        for (sensor, position, observation_id, marker, color) in solid_measurement_points
            scatter!(plot_object, [position * 1000.0],
                     [observation(measurements, simulation_id, observation_id)[end]];
                     label=String(sensor), color=color, markershape=marker,
                     markerstrokewidth=0, ms=5)
        end
        scatter!(plot_object, [T3_SAMPLE_POSITION_v18 * 1000.0], [experiment[end, 6]];
                 label="T3 experiment", color=:red, markershape=:star5,
                 markerstrokewidth=0, ms=7)
        scatter!(plot_object, [CAVITY_LENGTH_v18 * 1000.0], [experiment[end, 7]];
                 label="T2 experiment", color=:black, markershape=:diamond,
                 markerstrokewidth=0, ms=5)
        scatter!(plot_object, [CAVITY_LENGTH_v18 * 1000.0], [model[end, 7]];
                 label="T2 model", color=:black, markershape=:xcross,
                 markerstrokewidth=1, ms=5)
        vline!(plot_object, [L * 1000.0, (L + REAR_TUBE_CAVITY_LENGTH_v18) * 1000.0];
               label=false, color=:gray, linestyle=:dot, lw=1)
        return plot_object
    end

    function write_cell_diagnostics_v18(path, simulation_id, params;
                                       nodes=RUNNER_PLOT_NODES_v18,
                                       variant_name="v18")
        model, result, experiment = solve_case_v18(params, simulation_id; nodes=nodes)
        open(path, "w") do io
            println(io, "variant,simulation_id,z_mm,Tside_K,Tcore_K,Tgas_in_K,Tgas_out_K,Nu,h_W_m2K,k_side_W_m_K,UA_over_mcp,effectiveness,Qgas_W,Qcore_side_W")
            for i in eachindex(result.z_solid)
                println(io, join((
                    variant_name,
                    simulation_id,
                    result.z_solid[i] * 1000.0,
                    result.wall_temperature[i, end],
                    result.active_temperature[i, end],
                    result.gas_temperature[i, end],
                    result.gas_temperature[i + 1, end],
                    result.receiver_nusselt[i, end],
                    result.heat_transfer_coefficient[i, end],
                    side_axial_participation_conductivity_v18(result.wall_temperature[i, end], params),
                    result.receiver_ua_over_mcp[i, end],
                    result.receiver_effectiveness[i, end],
                    result.receiver_gas_heat[i, end],
                    core_side_conductance_per_length_v18(params) *
                    (L / length(result.z_solid)) *
                    (result.active_temperature[i, end] - result.wall_temperature[i, end]),
                ), ','))
            end
        end
        return path
    end

    function write_namedtuple_rows_v18(path, rows)
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

begin # v18 heat-transfer refit
    println("[run_1D_v18] Loaded $(length(sim_key_heat)) heating runs and " *
            "$(length(sim_key_cool)) cooling runs.")
    println("[run_1D_v18] Fitting only T9/T10 measurement-bias parameters; " *
            "the v17 thermal topology is fixed.")

    fit = calibrate_v18(
        nodes=RUNNER_FIT_NODES_v18,
        maximum_iterations=RUNNER_FIT_ITERATIONS_v18,
        maximum_time_seconds=RUNNER_FIT_SECONDS_v18,
    )
    fitted_params = fit.parameters
    open(joinpath(RUNNER_OUTPUT_DIR_v18, "optimization_summary_1D_v18.txt"), "w") do io
        println(io, "objective=$(fit.minimum)")
        println(io, "return_code=$(fit.retcode)")
        println(io, "iterations=$(get(fit, :iterations, missing))")
        println(io, "parameters=$(fit.parameters)")
    end

    variants_v18 = NamedTuple[
        (name="initial_measurement_bias", params=copy(pnew_v18)),
        (name="fitted_measurement_bias", params=fitted_params),
    ]

    summary_rows = NamedTuple[]
    representative_heat = ("E67", "E76", "E80")

    for variant in variants_v18
        println("[run_1D_v18] Running variant $(variant.name).")
        write_parameters_v18(
            joinpath(RUNNER_OUTPUT_DIR_v18, "parameters_$(variant.name)_1D_v18.csv"),
            variant,
        )
        metrics = compute_metrics_v18(variant.params; nodes=RUNNER_PLOT_NODES_v18)
        write_metrics_v18(
            joinpath(RUNNER_OUTPUT_DIR_v18, "analysis_results_$(variant.name)_1D_v18.csv"),
            metrics,
        )
        for row in metrics
            push!(summary_rows, merge(row, (variant=variant.name,)))
        end
        steady_results = build_steady_results_v18(
            variant.params; nodes=RUNNER_PLOT_NODES_v18,
        )
        write_steady_results_v18(
            joinpath(RUNNER_OUTPUT_DIR_v18, "steady_results_$(variant.name)_1D_v18.csv"),
            steady_results,
        )
        write_flow_slopes_v18(
            joinpath(RUNNER_OUTPUT_DIR_v18, "flow_slopes_$(variant.name)_1D_v18.csv"),
            steady_results,
        )
        save_runner_plot_v18(
            steady_comparison_plot_v18(steady_results; title_suffix=variant.name),
            "steady_comparison_$(variant.name)_1D_v18.png",
        )
        axial_keys = variant.name == "fitted_measurement_bias" ? sim_key_heat : representative_heat
        transient_heat_keys = variant.name == "fitted_measurement_bias" ? sim_key_heat : representative_heat
        transient_cool_keys = variant.name == "fitted_measurement_bias" ? sim_key_cool : ()
        for simulation_id in axial_keys
            save_runner_plot_v18(
                axial_profile_plot_v18(
                    simulation_id, variant.params; variant_name=variant.name,
                ),
                joinpath(
                    "axial_profiles",
                    "axial_profile_$(simulation_id)_$(variant.name)_1D_v18.png",
                ),
            )
            write_cell_diagnostics_v18(
                joinpath(
                    RUNNER_DIAGNOSTIC_DIR_v18,
                    "cell_diagnostics_$(simulation_id)_$(variant.name)_1D_v18.csv",
                ),
                simulation_id,
                variant.params;
                variant_name=variant.name,
            )
        end
        for simulation_id in transient_heat_keys
            save_runner_plot_v18(
                plot_case_v18(
                    simulation_id, variant.params; nodes=RUNNER_PLOT_NODES_v18,
                ),
                joinpath(
                    "transients",
                    "transient_$(simulation_id)_$(variant.name)_1D_v18.png",
                ),
            )
        end
        for simulation_id in transient_cool_keys
            save_runner_plot_v18(
                plot_case_v18(
                    simulation_id, variant.params; is_cooling=true,
                    nodes=RUNNER_PLOT_NODES_v18,
                ),
                joinpath(
                    "transients",
                    "transient_$(simulation_id)_$(variant.name)_cooling_1D_v18.png",
                ),
            )
        end
    end

    summary_path = joinpath(RUNNER_OUTPUT_DIR_v18, "analysis_results_all_variants_1D_v18.csv")
    open(summary_path, "w") do io
        println(io, "variant,simulation_id,phase,sensor,rmse_K,steady_error_K,t90_error_s,shape_loss")
        for row in summary_rows
            println(io, join((row.variant, row.simulation_id, row.phase, row.sensor,
                              row.rmse_K, row.steady_error_K,
                              row.t90_error_s, row.shape_loss), ','))
        end
    end

    println("[run_1D_v18] Complete.")
    println("[run_1D_v18] Output directory: $RUNNER_OUTPUT_DIR_v18")
end



