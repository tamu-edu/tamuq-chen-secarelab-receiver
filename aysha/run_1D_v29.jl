# ============================================================================
# run_1D_v29.jl - Energy-accounting 2-Zone Core/Perimeter Runner
# ============================================================================
# Runner for v29: shape-only flux-prior two-zone macro model with corrected
# zero-irradiance cooling, flux-guided perimeter partition, a bounded rear
# reservoir, and a smooth gated hardware flange path. The flux-map absolute
# magnitudes are not used.
# Generates parameters, metrics, steady results, flow slopes, cell diagnostics,
# axial profile plots, and transient plots.
# ============================================================================

begin # libraries
    using Statistics
    using StatsPlots
end

include("1D_v29.jl")

begin # runner configuration
    runner_int_v29(name, default) = parse(Int, get(ENV, name, string(default)))
    runner_float_v29(name, default) = parse(Float64, get(ENV, name, string(default)))

    const RUNNER_PLOT_NODES_v29 = runner_int_v29("RECEIVER1D_v29_RUNNER_PLOT_NODES", default_nodes)
    const RUNNER_FIT_NODES_v29 = runner_int_v29("RECEIVER1D_v29_FIT_NODES", 15)
    const RUNNER_FIT_ITERATIONS_v29 = runner_int_v29("RECEIVER1D_v29_FIT_ITERATIONS", 100)
    const RUNNER_FIT_SECONDS_v29 = runner_float_v29("RECEIVER1D_v29_FIT_SECONDS", 360.0)

    const RUNNER_OUTPUT_DIR_v29 = joinpath(@__DIR__, "summaries", "1D_v29")
    const RUNNER_PLOT_DIR_v29 = joinpath(RUNNER_OUTPUT_DIR_v29, "plots")
    const RUNNER_TRANSIENT_DIR_v29 = joinpath(RUNNER_PLOT_DIR_v29, "transients")
    const RUNNER_AXIAL_DIR_v29 = joinpath(RUNNER_PLOT_DIR_v29, "axial_profiles")
    const RUNNER_DIAGNOSTIC_DIR_v29 = joinpath(RUNNER_OUTPUT_DIR_v29, "diagnostics")
    mkpath(RUNNER_OUTPUT_DIR_v29)
    mkpath(RUNNER_PLOT_DIR_v29)
    mkpath(RUNNER_TRANSIENT_DIR_v29)
    mkpath(RUNNER_AXIAL_DIR_v29)
    mkpath(RUNNER_DIAGNOSTIC_DIR_v29)
end

begin # runner helpers
    const RUNNER_SENSOR_COLORS_v29 = Dict(
        :T8 => :blue,
        :T12_perim => :red,
        :T11_perim => :green,
        :T9_core => :orange,
        :T10_core => :brown,
        :T3 => :purple,
        :T2 => :black,
    )

    const TEMPERATURE_AXIS_LIMITS_v29 = (280.0, 1400.0)

    function save_runner_plot_v29(plot_object, filename)
        path = joinpath(RUNNER_PLOT_DIR_v29, filename)
        savefig(plot_object, path)
        println("[run_1D_v29] Saved plot: $path")
        return path
    end

    scaled_incident_flux_w_m2_v29(irradiance, scale) =
        max(0.0, irradiance) * scale

    absorbed_receiver_power_watts_v29(irradiance, scale) =
        ETA_ABS_FIXED_v29 * scaled_incident_flux_w_m2_v29(irradiance, scale) * A_frt

    function write_parameters_v29(path, variant)
        parameters = variant.params
        names = (
            "A_Nu", "B_Re", "C_Pr_fixed", "phi_0", "m_rec",
            "front_dep", "scale_456", "scale_304", "scale_256",
            "G_core_perim_W_m_K", "C_perim_eff_J_K", "k_perim_ref_W_m_K", "beta_opt",
            "spill_capture", "beta_perim_1_m", "f_core_rear", "flange_scale",
            "flange_cool_gain", "flange_cool_tau_s", "k_core_axial_scale",
            "C_rear_eff_J_K", "G_receiver_rear_W_K", "G_rear_tube_W_K",
            "G_rear_cavity_W_K",
        )
        open(path, "w") do io
            println(io, "index,name,value")
            for i in eachindex(parameters)
                println(io, join((i, names[i], parameters[i]), ','))
            end
            println(io, "# fitted,phi_0,$(parameters[4])")
            println(io, "# fitted,m_rec,$(parameters[5])")
            println(io, "# fitted,G_core_perim_W_m_K,$(core_perimeter_conductance_per_length_v29(parameters))")
            println(io, "# fitted,C_perim_eff_J_K,$(perimeter_heat_capacity_total_v29(parameters))")
            println(io, "# fitted,k_perim_ref_W_m_K,$(perimeter_axial_conductivity_v29(900.0, parameters))")
            println(io, "# fitted,spill_capture,$(perimeter_spillover_capture_v29(parameters))")
            println(io, "# derived,f_perim_source,$(perimeter_source_fraction_v29(parameters))")
            println(io, "# derived,flux_receiver_fraction,$(flux_receiver_fraction_v29())")
            println(io, "# derived,flux_spillover_fraction,$(flux_spillover_fraction_v29())")
            println(io, "# fitted,beta_perim_1_m,$(perimeter_source_attenuation_v29(parameters))")
            println(io, "# fitted,f_core_rear,$(rear_core_fraction_v29(parameters))")
            println(io, "# fitted,flange_scale,$(flange_loss_scale_v29(parameters))")
            println(io, "# fitted,flange_cool_gain,$(flange_cooling_gain_v29(parameters))")
            println(io, "# fitted,flange_cool_tau_s,$(flange_cooling_time_constant_v29(parameters))")
            println(io, "# fitted,k_core_axial_scale,$(core_axial_conduction_scale_v29(parameters))")
            println(io, "# fitted,C_rear_eff_J_K,$(rear_reservoir_heat_capacity_v29(parameters))")
            println(io, "# fitted,G_receiver_rear_W_K,$(receiver_rear_conductance_v29(parameters))")
            println(io, "# fitted,G_rear_tube_W_K,$(rear_tube_conductance_v29(parameters))")
            println(io, "# fitted,G_rear_cavity_W_K,$(rear_cavity_conductance_v29(parameters))")
            println(io, "# derived,C_core_eff_J_K,$(core_receiver_heat_capacity_v29(default_nodes))")
            println(io, "# derived,C_receiver_participating_eff_J_K,$(participating_assembly_heat_capacity_v29(parameters))")
            println(io, "# derived,C_total_with_rear_eff_J_K,$(participating_total_heat_capacity_v29(parameters))")
            println(io, "# reference,measured_C_eff_J_K,$MEASURED_ASSEMBLY_CAPACITANCE_v29")
        end
        return path
    end

    function compute_metrics_v29(p=pnew_v29; heating_keys=sim_key_heat,
                                cooling_keys=sim_key_cool, nodes=RUNNER_PLOT_NODES_v29)
        metrics = NamedTuple[]
        for (keys, cooling) in ((heating_keys, false), (cooling_keys, true))
            for simulation_id in keys
                outputs, result, experiment = solve_case_v29(
                    p, simulation_id; is_cooling=cooling, nodes=nodes,
                )
                sensors = (:T8, :T12_perim, :T11_perim, :T9_core, :T10_core, :T3, :T2)
                for (j, sensor) in enumerate(sensors)
                    push!(metrics, (
                        simulation_id=simulation_id,
                        phase=cooling ? :cooling : :heating,
                        sensor=sensor,
                        rmse_K=sqrt(mean(abs2, outputs[:, j] .- experiment[:, j])),
                        steady_error_K=outputs[end, j] - experiment[end, j],
                        t90_error_s=get_t90_v3(result.time, outputs[:, j]) - get_t90_v3(result.time, experiment[:, j]),
                        shape_loss=normalized_slope_mse_v29(outputs[:, j], experiment[:, j]),
                    ))
                end
            end
        end
        return metrics
    end

    function write_metrics_v29(path, metrics)
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

    function build_steady_results_v29(params=pnew_v29; keys=sim_key_heat, nodes=RUNNER_PLOT_NODES_v29)
        results = NamedTuple[]
        core_weights = solar_weights_v5(params[13], params[6], nodes)
        perim_weights = axial_exponential_weights_v29(
            perimeter_source_attenuation_v29(params), nodes, L / nodes,
        )
        for simulation_id in keys
            outputs, result, experiment = solve_case_v29(params, simulation_id; nodes=nodes)
            conditions = simulation_conditions[simulation_id]
            irradiance = conditions[Io]
            scale = absorbed_power_scale_v29(irradiance, params)
            nominal_flux = max(0.0, irradiance)
            scaled_flux = scaled_incident_flux_w_m2_v29(irradiance, scale)
            nominal_absorbed = absorbed_receiver_power_watts_v29(irradiance, 1.0)
            scaled_absorbed = absorbed_receiver_power_watts_v29(irradiance, scale)
            f_perim_source = perimeter_source_fraction_v29(params)
            core_absorbed = core_absorbed_power_v29(irradiance, params)
            perim_absorbed = perimeter_absorbed_power_v29(irradiance, params)
            participating_absorbed = modeled_participating_absorbed_power_v29(irradiance, params)
            push!(results, (
                simulation_id=simulation_id,
                flow_lpm=conditions[qlpm],
                irradiance=irradiance,
                power_scale=scale,
                nominal_incident_flux_W_m2=nominal_flux,
                scaled_incident_flux_W_m2=scaled_flux,
                nominal_receiver_absorbed_W=nominal_absorbed,
                scaled_receiver_absorbed_W=scaled_absorbed,
                participating_absorbed_W=participating_absorbed,
                core_absorbed_W=core_absorbed,
                perim_absorbed_W=perim_absorbed,
                spill_capture=perimeter_spillover_capture_v29(params),
                flux_receiver_fraction=flux_receiver_fraction_v29(),
                flux_spillover_fraction=flux_spillover_fraction_v29(),
                f_perim_source=f_perim_source,
                beta_perim_1_m=perimeter_source_attenuation_v29(params),
                core_front_cell_source_fraction=core_weights[1],
                perim_front_cell_source_fraction=perim_weights[1],
                core_downstream_source_fraction=1.0 - core_weights[1],
                perim_downstream_source_fraction=1.0 - perim_weights[1],
                T8_model=outputs[end, 1],
                T8_experiment=experiment[end, 1],
                T12_perim_model=outputs[end, 2],
                T12_perim_experiment=experiment[end, 2],
                T11_perim_model=outputs[end, 3],
                T11_perim_experiment=experiment[end, 3],
                T9_core_model=outputs[end, 4],
                T9_core_experiment=experiment[end, 4],
                T10_core_model=outputs[end, 5],
                T10_core_experiment=experiment[end, 5],
                T3_model=outputs[end, 6],
                T3_experiment=experiment[end, 6],
                T2_model=outputs[end, 7],
                T2_experiment=experiment[end, 7],
                T12_minus_T8_model=outputs[end, 2] - outputs[end, 1],
                T12_minus_T8_experiment=experiment[end, 2] - experiment[end, 1],
                T11_minus_T8_model=outputs[end, 3] - outputs[end, 1],
                T11_minus_T8_experiment=experiment[end, 3] - experiment[end, 1],
                T12_minus_T9_model=outputs[end, 2] - outputs[end, 4],
                T12_minus_T9_experiment=experiment[end, 2] - experiment[end, 4],
                T11_minus_T10_model=outputs[end, 3] - outputs[end, 5],
                T11_minus_T10_experiment=experiment[end, 3] - experiment[end, 5],
                mean_nu=mean(result.receiver_nusselt[:, end]),
                mean_effectiveness=mean(result.receiver_effectiveness[:, end]),
                receiver_gas_heat=sum(result.receiver_gas_heat[:, end]),
                rear_tube_gas_heat=sum(result.rear_tube_gas_heat[:, end]),
                active_flow_fraction=result.active_flow_fraction[end],
                G_core_perim_W_m_K=core_perimeter_conductance_per_length_v29(params),
                C_perim_eff_J_K=perimeter_heat_capacity_total_v29(params),
                C_receiver_participating_eff_J_K=participating_assembly_heat_capacity_v29(params, nodes),
                C_rear_eff_J_K=rear_reservoir_heat_capacity_v29(params),
                C_total_with_rear_eff_J_K=participating_total_heat_capacity_v29(params, nodes),
                f_core_rear=rear_core_fraction_v29(params),
                G_receiver_rear_W_K=receiver_rear_conductance_v29(params),
                G_rear_tube_W_K=rear_tube_conductance_v29(params),
                G_rear_cavity_W_K=rear_cavity_conductance_v29(params),
                T_rear_reservoir=result.rear_reservoir_temperature[end],
                flange_scale=flange_loss_scale_v29(params),
                flange_cool_gain=flange_cooling_gain_v29(params),
                flange_cool_tau_s=flange_cooling_time_constant_v29(params),
                flange_gate_end=lamp_off_gate_v29(result.time[end], irradiance, params),
                effective_flange_scale_end=effective_flange_loss_scale_v29(params, result.time[end], irradiance),
                k_core_axial_scale=core_axial_conduction_scale_v29(params),
            ))
        end
        return results
    end

    function write_steady_results_v29(path, steady_results)
        fields = propertynames(first(steady_results))
        open(path, "w") do io
            println(io, join(fields, ','))
            for row in steady_results
                println(io, join((getproperty(row, field) for field in fields), ','))
            end
        end
        return path
    end

    function slope_v29(xs, ys)
        xmean = mean(xs)
        ymean = mean(ys)
        denominator = sum((x .- xmean)^2 for x in xs)
        denominator <= eps(Float64) && return NaN
        return sum((xs[i] - xmean) * (ys[i] - ymean) for i in eachindex(xs)) / denominator
    end

    function write_flow_slopes_v29(path, steady_results)
        sensors = (:T8, :T12_perim, :T11_perim, :T9_core, :T10_core, :T3)
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
                        slope_v29(flows, model_values),
                        slope_v29(flows, experiment_values),
                    ), ','))
                end
            end
        end
        return path
    end

    function steady_comparison_plot_v29(steady_results; title_suffix="")
        sensors = (:T8, :T12_perim, :T11_perim, :T9_core, :T10_core, :T3, :T2)
        plot_object = plot(
            title="1D_v29 $title_suffix steady-state comparison",
            xlabel="Model temperature (K)",
            ylabel="Experiment temperature (K)",
            legend=:outerright,
            aspect_ratio=:equal,
            xlims=TEMPERATURE_AXIS_LIMITS_v29,
            ylims=TEMPERATURE_AXIS_LIMITS_v29,
        )
        for sensor in sensors
            model_values = [getproperty(row, Symbol(sensor, "_model")) for row in steady_results]
            experiment_values = [getproperty(row, Symbol(sensor, "_experiment")) for row in steady_results]
            scatter!(plot_object, model_values, experiment_values;
                     label=String(sensor), color=RUNNER_SENSOR_COLORS_v29[sensor],
                     ms=4, markerstrokewidth=0)
        end
        lower, upper = TEMPERATURE_AXIS_LIMITS_v29
        plot!(plot_object, [lower, upper], [lower, upper];
              label="1:1", color=:gray, linestyle=:dash)
        return plot_object
    end

    function axial_profile_plot_v29(simulation_id, params=pnew_v29;
                                   nodes=RUNNER_PLOT_NODES_v29,
                                   variant_name="v29")
        outputs, result, experiment = solve_case_v29(params, simulation_id; nodes=nodes)
        conditions = simulation_conditions[simulation_id]

        perim_x_mm = vcat(result.z_solid, result.z_rear_tube) .* 1000.0
        perim_T = vcat(result.perim_temperature[:, end], result.rear_tube_temperature[:, end])
        core_x_mm = result.z_solid .* 1000.0
        core_T = result.core_temperature[:, end]
        gas_x_mm = result.z_gas .* 1000.0
        gas_T = result.gas_temperature[:, end]
        rear_reservoir_T = result.rear_reservoir_temperature[end]

        solid_measurement_points = (
            (:T8, sensor_positions[:T8], "_T8", :square, :blue),
            (:T9_core, sensor_positions[:T9], "_T9", :circle, :orange),
            (:T12_perim, sensor_positions[:T9], "_T12", :utriangle, :red),
            (:T10_core, sensor_positions[:T10], "_T10", :diamond, :brown),
            (:T11_perim, sensor_positions[:T10], "_T11", :dtriangle, :green),
        )

        scale = absorbed_power_scale_v29(conditions[Io], params)
        title_text = "1D_v29 $variant_name axial profile: $simulation_id"
        subtitle_text = "flow=$(round(conditions[qlpm]; digits=2)) L/min, " *
                        "Io=$(round(conditions[Io] / 1000.0; digits=1)) kW/m^2, " *
                        "scale=$(round(scale; digits=3))"
        plot_object = plot(
            title=title_text * "\n" * subtitle_text,
            xlabel="Length (mm)",
            ylabel="Temperature (K)",
            legend=:outerright,
            grid=true,
            xlims=(0.0, 1000.0 * (L + REAR_TUBE_LENGTH_v29)),
            ylims=TEMPERATURE_AXIS_LIMITS_v29,
        )
        plot!(plot_object, perim_x_mm, perim_T;
              label="perimeter/tube model", color=:red, linestyle=:dash, lw=2)
        plot!(plot_object, core_x_mm, core_T;
              label="core solid model", color=:orange, linestyle=:dot, lw=2)
        plot!(plot_object, gas_x_mm, gas_T;
              label="gas model", color=:green, lw=2)
        scatter!(plot_object, [L * 1000.0], [rear_reservoir_T];
                 label="rear reservoir model", color=:purple,
                 markershape=:diamond, markerstrokewidth=0, ms=6)
        for (sensor, position, observation_id, marker, color) in solid_measurement_points
            scatter!(plot_object, [position * 1000.0],
                     [observation(measurements, simulation_id, observation_id)[end]];
                     label=String(sensor), color=color, markershape=marker,
                     markerstrokewidth=0, ms=5)
        end
        scatter!(plot_object, [T3_SAMPLE_POSITION_v29 * 1000.0], [experiment[end, 6]];
                 label="T3 experiment", color=:purple, markershape=:star5,
                 markerstrokewidth=0, ms=7)
        scatter!(plot_object, [CAVITY_LENGTH_v29 * 1000.0], [experiment[end, 7]];
                 label="T2 experiment", color=:black, markershape=:diamond,
                 markerstrokewidth=0, ms=5)
        scatter!(plot_object, [CAVITY_LENGTH_v29 * 1000.0], [outputs[end, 7]];
                 label="T2 model", color=:black, markershape=:xcross,
                 markerstrokewidth=1, ms=5)
        vline!(plot_object, [L * 1000.0, (L + REAR_TUBE_CAVITY_LENGTH_v29) * 1000.0];
               label=false, color=:gray, linestyle=:dot, lw=1)
        return plot_object
    end

    function plot_case_v29(simulation_id, params=pnew_v29; is_cooling=false,
                          nodes=RUNNER_PLOT_NODES_v29, variant_name="v29")
        outputs, result, experiment = solve_case_v29(params, simulation_id; is_cooling=is_cooling, nodes=nodes)
        sensors = (:T8, :T12_perim, :T11_perim, :T9_core, :T10_core, :T3, :T2)
        colors = (:blue, :red, :green, :orange, :brown, :purple, :black)
        
        plot_object = plot(
            title="1D_v29 $variant_name transient: $simulation_id",
            xlabel="Time (s)", ylabel="Temperature (K)",
            legend=:outerright,
            ylims=TEMPERATURE_AXIS_LIMITS_v29,
        )
        for (j, (sensor, color)) in enumerate(zip(sensors, colors))
            plot!(plot_object, result.time, outputs[:, j];
                  label="$(sensor) model", lw=2, color=color)
            scatter!(plot_object, result.time, experiment[:, j];
                     label="$(sensor) experiment", ms=2, markerstrokewidth=0,
                     color=color)
        end
        return plot_object
    end

    function write_cell_diagnostics_v29(path, simulation_id, params;
                                       nodes=RUNNER_PLOT_NODES_v29,
                                       variant_name="v29")
        outputs, result, experiment = solve_case_v29(params, simulation_id; nodes=nodes)
        conditions = simulation_conditions[simulation_id]
        scaled_absorbed = absorbed_receiver_power_watts_v29(
            conditions[Io], absorbed_power_scale_v29(conditions[Io], params),
        )
        f_perim_source = perimeter_source_fraction_v29(params)
        core_absorbed = core_absorbed_power_v29(conditions[Io], params)
        perim_absorbed = perimeter_absorbed_power_v29(conditions[Io], params)
        dx = L / length(result.z_solid)
        f_core_rear = rear_core_fraction_v29(params)
        G_receiver_rear = receiver_rear_conductance_v29(params)
        Trear = result.rear_reservoir_temperature[end]
        Ttube_in = result.rear_tube_temperature[1, end]
        Tcavity = result.cavity_temperature[end]
        Qcore_rear_exit = G_receiver_rear * f_core_rear *
                          (result.core_temperature[end, end] - Trear)
        Qperim_rear_exit = G_receiver_rear * (1.0 - f_core_rear) *
                           (result.perim_temperature[end, end] - Trear)
        Qrear_tube_exit = rear_tube_conductance_v29(params) * (Trear - Ttube_in)
        Qrear_cavity_exit = rear_cavity_conductance_v29(params) * (Trear - Tcavity)
        open(path, "w") do io
            println(io, "variant,simulation_id,z_mm,Tperim_K,Tcore_K,Tgas_in_K,Tgas_out_K,Nu,h_W_m2K,k_perim_W_m_K,effectiveness,Qgas_W,Qcore_source_W,Qperim_source_W,Qcore_perim_W,Qcore_rear_W,Qperim_rear_W,Trear_K,Qrear_tube_W,Qrear_cavity_W,active_flow_fraction,flange_scale,flange_cool_gain,flange_cool_tau_s,k_core_axial_scale")
            for i in eachindex(result.z_solid)
                is_exit_cell = i == length(result.z_solid)
                println(io, join((
                    variant_name,
                    simulation_id,
                    result.z_solid[i] * 1000.0,
                    result.perim_temperature[i, end],
                    result.core_temperature[i, end],
                    result.gas_temperature[i, end],
                    result.gas_temperature[i + 1, end],
                    result.receiver_nusselt[i, end],
                    result.heat_transfer_coefficient[i, end],
                    perimeter_axial_conductivity_v29(result.perim_temperature[i, end], params),
                    result.receiver_effectiveness[i, end],
                    result.receiver_gas_heat[i, end],
                    core_absorbed * result.core_source_weights[i],
                    perim_absorbed * result.perim_source_weights[i],
                    core_perimeter_conductance_per_length_v29(params) *
                    dx *
                    (result.perim_temperature[i, end] - result.core_temperature[i, end]),
                    is_exit_cell ? Qcore_rear_exit : 0.0,
                    is_exit_cell ? Qperim_rear_exit : 0.0,
                    Trear,
                    is_exit_cell ? Qrear_tube_exit : 0.0,
                    is_exit_cell ? Qrear_cavity_exit : 0.0,
                    result.active_flow_fraction[end],
                    flange_loss_scale_v29(params),
                    flange_cooling_gain_v29(params),
                    flange_cooling_time_constant_v29(params),
                    core_axial_conduction_scale_v29(params),
                ), ','))
            end
        end
        return path
    end
end

begin # execution & output generation
    println("[run_1D_v29] Running energy-accounting 2-Zone Core/Perimeter Runner.")
    println("[run_1D_v29] Fitting power scales, source partition, core-perimeter coupling, rear reservoir, and flange path...")

    fit = calibrate_v29(
        nodes=RUNNER_FIT_NODES_v29,
        maximum_iterations=RUNNER_FIT_ITERATIONS_v29,
        maximum_time_seconds=RUNNER_FIT_SECONDS_v29,
    )
    fitted_params = fit.parameters

    open(joinpath(RUNNER_OUTPUT_DIR_v29, "optimization_summary_1D_v29.txt"), "w") do io
        println(io, "objective=$(fit.objective)")
        println(io, "return_code=$(fit.retcode)")
        println(io, "parameters=$(fit.parameters)")
    end

    variants_v29 = NamedTuple[
        (name="fitted_energy_accounting", params=fitted_params),
    ]

    summary_rows = NamedTuple[]
    for variant in variants_v29
        println("[run_1D_v29] Running variant $(variant.name).")
        write_parameters_v29(
            joinpath(RUNNER_OUTPUT_DIR_v29, "parameters_$(variant.name)_1D_v29.csv"),
            variant,
        )
        metrics = compute_metrics_v29(variant.params; nodes=RUNNER_PLOT_NODES_v29)
        write_metrics_v29(
            joinpath(RUNNER_OUTPUT_DIR_v29, "analysis_results_$(variant.name)_1D_v29.csv"),
            metrics,
        )
        for row in metrics
            push!(summary_rows, merge(row, (variant=variant.name,)))
        end
        steady_results = build_steady_results_v29(
            variant.params; nodes=RUNNER_PLOT_NODES_v29,
        )
        write_steady_results_v29(
            joinpath(RUNNER_OUTPUT_DIR_v29, "steady_results_$(variant.name)_1D_v29.csv"),
            steady_results,
        )
        write_flow_slopes_v29(
            joinpath(RUNNER_OUTPUT_DIR_v29, "flow_slopes_$(variant.name)_1D_v29.csv"),
            steady_results,
        )
        save_runner_plot_v29(
            steady_comparison_plot_v29(steady_results; title_suffix=variant.name),
            "steady_comparison_$(variant.name)_1D_v29.png",
        )
        
        axial_keys = sim_key_heat
        transient_heat_keys = sim_key_heat
        transient_cool_keys = sim_key_cool
        
        for simulation_id in axial_keys
            save_runner_plot_v29(
                axial_profile_plot_v29(
                    simulation_id, variant.params; variant_name=variant.name,
                    nodes=RUNNER_PLOT_NODES_v29,
                ),
                joinpath("axial_profiles", "axial_profile_$(simulation_id)_$(variant.name)_1D_v29.png"),
            )
            write_cell_diagnostics_v29(
                joinpath(RUNNER_DIAGNOSTIC_DIR_v29, "cell_diagnostics_$(simulation_id)_$(variant.name)_1D_v29.csv"),
                simulation_id, variant.params; variant_name=variant.name,
                nodes=RUNNER_PLOT_NODES_v29,
            )
        end
        
        for simulation_id in transient_heat_keys
            save_runner_plot_v29(
                plot_case_v29(
                    simulation_id, variant.params; is_cooling=false,
                    variant_name=variant.name, nodes=RUNNER_PLOT_NODES_v29,
                ),
                joinpath("transients", "transient_$(simulation_id)_$(variant.name)_1D_v29.png"),
            )
        end
        
        for simulation_id in transient_cool_keys
            save_runner_plot_v29(
                plot_case_v29(
                    simulation_id, variant.params; is_cooling=true,
                    variant_name=variant.name, nodes=RUNNER_PLOT_NODES_v29,
                ),
                joinpath("transients", "transient_$(simulation_id)_$(variant.name)_cooling_1D_v29.png"),
            )
        end
    end

    summary_path = joinpath(RUNNER_OUTPUT_DIR_v29, "analysis_results_all_variants_1D_v29.csv")
    open(summary_path, "w") do io
        println(io, "variant,simulation_id,phase,sensor,rmse_K,steady_error_K,t90_error_s,shape_loss")
        for row in summary_rows
            println(io, join((row.variant, row.simulation_id, row.phase, row.sensor,
                              row.rmse_K, row.steady_error_K,
                              row.t90_error_s, row.shape_loss), ','))
        end
    end

    println("[run_1D_v29] Optimization completed with objective: $(fit.objective)")
    println("[run_1D_v29] Fitted parameters: $fitted_params")
    println("[run_1D_v29] Complete.")
    println("[run_1D_v29] Output directory: $RUNNER_OUTPUT_DIR_v29")
end





