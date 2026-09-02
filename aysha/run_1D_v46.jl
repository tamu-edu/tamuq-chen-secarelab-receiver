# ============================================================================
# run_1D_v46.jl - Energy-accounting 2-Zone Core/Perimeter Runner
# ============================================================================
# Runner for v46 Entire Converter Model:
#   - Option A: Calibrate on 15 heating runs (E67-E81), validate out-of-sample
#     on 3 cooling runs (C69, C80, C81).
#   - Synchronized geometry and polynomial properties from 1D_v3.jl.
#   - Total mass flow through honeycomb core and rear tube (phi_act = 1.0).
#   - Coupled front suction preheating at z = 0.
#   - Perimeter spillage absorbed at z = 0 front face with aluminum casing conduction.
#   - Phase-invariant physical flange conductance (no lamp-off multipliers).
#   - Full energy accounting with exact conservation residual verification.
# ============================================================================

begin # libraries
    using Statistics
    using StatsPlots
    using LinearAlgebra
    using Optimization
    using OptimizationNLopt
end

include("1D_v46.jl")

begin # runner configuration
    runner_int_v46(name, default) = parse(Int, get(ENV, name, string(default)))
    runner_float_v46(name, default) = parse(Float64, get(ENV, name, string(default)))
    runner_bool_v46(name, default) =
        lowercase(get(ENV, name, default ? "true" : "false")) in ("1", "true", "yes", "on")

    const RUNNER_PLOT_NODES_v46 = runner_int_v46("RECEIVER1D_v46_RUNNER_PLOT_NODES", default_nodes)
    const RUNNER_FIT_NODES_v46 = runner_int_v46("RECEIVER1D_v46_FIT_NODES", 15)
    const RUNNER_FIT_ITERATIONS_v46 = runner_int_v46("RECEIVER1D_v46_FIT_ITERATIONS", 1000)
    const RUNNER_FIT_SECONDS_v46 = runner_float_v46("RECEIVER1D_v46_FIT_SECONDS", 1800.0)
    const RUNNER_FIT_STAGE_v46 = Symbol(get(ENV, "RECEIVER1D_v46_FIT_STAGE", "full"))
    const RUNNER_CALIBRATION_DATASET_v46 = Symbol(get(ENV, "RECEIVER1D_v46_CALIBRATION_DATASET", "heating"))
    const RUNNER_WRITE_PLOTS_v46 = runner_bool_v46("RECEIVER1D_v46_WRITE_PLOTS", true)

    const RUNNER_OUTPUT_DIR_v46 = joinpath(@__DIR__, "summaries", "1D_v46")
    const RUNNER_PLOT_DIR_v46 = joinpath(RUNNER_OUTPUT_DIR_v46, "plots")
    const RUNNER_TRANSIENT_DIR_v46 = joinpath(RUNNER_PLOT_DIR_v46, "transients")
    const RUNNER_AXIAL_DIR_v46 = joinpath(RUNNER_PLOT_DIR_v46, "axial_profiles")
    const RUNNER_DIAGNOSTIC_DIR_v46 = joinpath(RUNNER_OUTPUT_DIR_v46, "diagnostics")
    mkpath(RUNNER_OUTPUT_DIR_v46)
    mkpath(RUNNER_PLOT_DIR_v46)
    mkpath(RUNNER_TRANSIENT_DIR_v46)
    mkpath(RUNNER_AXIAL_DIR_v46)
    mkpath(RUNNER_DIAGNOSTIC_DIR_v46)
end

begin # runner helpers
    const RUNNER_SENSOR_COLORS_v46 = Dict(
        :T8 => :blue,
        :T12_perim => :red,
        :T11_perim => :green,
        :T9_core => :orange,
        :T10_core => :brown,
        :T3 => :purple,
        :T2 => :black,
    )

    const TEMPERATURE_AXIS_LIMITS_v46 = (280.0, 1400.0)

    function save_runner_plot_v46(plot_object, filename)
        path = joinpath(RUNNER_PLOT_DIR_v46, filename)
        savefig(plot_object, path)
        println("[run_1D_v46] Saved plot: $path")
        return path
    end

    scaled_incident_flux_w_m2_v46(irradiance, scale) =
        max(0.0, irradiance) * scale

    absorbed_receiver_power_watts_v46(irradiance, scale) =
        ETA_ABS_FIXED_v46 * scaled_incident_flux_w_m2_v46(irradiance, scale) * A_frt

    function write_parameters_v46(path, variant)
        parameters = variant.params
        names = (
            "A_Nu", "B_Re", "C_Pr_fixed", "front_dep",
            "scale_456", "scale_304", "scale_256",
            "G_core_perim_W_m_K", "C_perim_eff_J_K", "k_perim_ref_W_m_K", "beta_opt",
            "chi", "f_core_rear", "flange_scale",
            "k_core_axial_scale", "C_rear_eff_J_K", "G_receiver_rear_W_K", "G_rear_tube_W_K",
            "G_rear_cavity_W_K", "G_rear_axial_W_K", "delta_web", "C_z", "h_suction"
        )
        open(path, "w") do io
            println(io, "index,name,value")
            for i in eachindex(parameters)
                println(io, join((i, names[i], parameters[i]), ','))
            end
            println(io, "# derived,f_perim_source,$(perimeter_source_fraction_v46(parameters))")
            println(io, "# derived,flux_receiver_fraction,$(flux_receiver_fraction_v46())")
            println(io, "# derived,flux_spillover_fraction,$(flux_spillover_fraction_v46())")
            println(io, "# derived,C_total_with_rear_eff_J_K,$(participating_total_heat_capacity_v46(parameters))")
            println(io, "# reference,measured_C_eff_J_K,$MEASURED_ASSEMBLY_CAPACITANCE_v46")
        end
        return path
    end

    function compute_metrics_v46(p=pnew_v46; heating_keys=keys_heating,
                                cooling_keys=keys_cooling, nodes=RUNNER_PLOT_NODES_v46)
        metrics = NamedTuple[]
        for (keys, cooling) in ((heating_keys, false), (cooling_keys, true))
            for simulation_id in keys
                outputs, result, experiment = solve_case_v46(
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
                        shape_loss=normalized_slope_mse_v46(outputs[:, j], experiment[:, j]),
                    ))
                end
            end
        end
        return metrics
    end

    function write_metrics_v46(path, metrics)
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

    function build_steady_results_v46(params=pnew_v46; keys=keys_heating, nodes=RUNNER_PLOT_NODES_v46)
        results = NamedTuple[]
        for simulation_id in keys
            outputs, result, experiment = solve_case_v46(params, simulation_id; nodes=nodes)
            conditions, time_vec, _, data = experimental_case_v46(simulation_id; is_cooling=false)
            irradiance = conditions[Io]
            flow = observation(data, simulation_id, "_flow")
            Tin = observation(data, simulation_id, "_Tin")
            ambient = observation(data, simulation_id, "_Tamb")
            operating = operating_condition_v3(
                irradiance=linear_history_v3(time_vec, fill(irradiance, length(time_vec))),
                flow_lpm=linear_history_v3(time_vec, flow),
                inlet_temperature=linear_history_v3(time_vec, Tin),
                ambient_temperature=linear_history_v3(time_vec, ambient),
            )
            
            balance = compute_energy_balance_v46(result, operating, result.time[end], params; nodes=nodes)

            push!(results, (
                simulation_id=simulation_id,
                flow_lpm=conditions[qlpm],
                irradiance=irradiance,
                power_scale=balance.delivered_W / max(irradiance * A_frt, eps(Float64)),
                delivered_power_W=balance.delivered_W,
                core_solar_W=balance.core_solar_W,
                perim_solar_W=balance.perim_solar_W,
                gas_heat_W=balance.gas_heat_W,
                front_rad_loss_W=balance.front_rad_loss_W,
                cavity_amb_loss_W=balance.cavity_amb_loss_W,
                flange_loss_W=balance.flange_loss_W,
                residual_W=balance.residual_W,
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
                T9_minus_T10_model=outputs[end, 4] - outputs[end, 5],
                T9_minus_T10_experiment=experiment[end, 4] - experiment[end, 5],
                mean_htc=mean(result.heat_transfer_coefficient[:, end]),
            ))
        end
        return results
    end

    function write_steady_results_v46(path, steady_results)
        fields = propertynames(first(steady_results))
        open(path, "w") do io
            println(io, join(fields, ','))
            for row in steady_results
                println(io, join((getproperty(row, field) for field in fields), ','))
            end
        end
        return path
    end

    function compute_flow_slopes_v46(steady_results)
        rows = NamedTuple[]
        for target_irradiance in IRRADIANCE_LEVELS_v46
            subset = filter(r -> abs(r.irradiance - target_irradiance) < 1000.0, steady_results)
            length(subset) < 2 && continue
            flows = [r.flow_lpm for r in subset]
            signals = (
                (:T8, [r.T8_model for r in subset], [r.T8_experiment for r in subset]),
                (:T12_perim, [r.T12_perim_model for r in subset], [r.T12_perim_experiment for r in subset]),
                (:T11_perim, [r.T11_perim_model for r in subset], [r.T11_perim_experiment for r in subset]),
                (:T9_core, [r.T9_core_model for r in subset], [r.T9_core_experiment for r in subset]),
                (:T10_core, [r.T10_core_model for r in subset], [r.T10_core_experiment for r in subset]),
                (:T3, [r.T3_model for r in subset], [r.T3_experiment for r in subset]),
                (:T2, [r.T2_model for r in subset], [r.T2_experiment for r in subset]),
            )
            for (signal_name, model_vals, exp_vals) in signals
                slope_model = cov(flows, model_vals) / var(flows)
                slope_exp = cov(flows, exp_vals) / var(flows)
                push!(rows, (
                    irradiance=target_irradiance,
                    signal=signal_name,
                    slope_model_K_per_LPM=slope_model,
                    slope_experiment_K_per_LPM=slope_exp,
                    error_K_per_LPM=slope_model - slope_exp,
                ))
            end
        end
        return rows
    end

    function write_flow_slopes_v46(path, rows)
        open(path, "w") do io
            println(io, "irradiance,signal,slope_model_K_per_LPM,slope_experiment_K_per_LPM,error_K_per_LPM")
            for r in rows
                println(io, join((r.irradiance, r.signal, r.slope_model_K_per_LPM, r.slope_experiment_K_per_LPM, r.error_K_per_LPM), ','))
            end
        end
        return path
    end

    function compute_invariants_v46(p, steady_results)
        rows = NamedTuple[]
        C_tot = participating_total_heat_capacity_v46(p)
        push!(rows, (metric="C_total_with_rear", value=C_tot, target="301 +/- 23", status=abs(C_tot - 301.0) <= 23.0 ? "pass" : "diagnostic"))
        
        # Max steady conservation residual
        max_res = maximum(abs.(r.residual_W for r in steady_results))
        push!(rows, (metric="max_steady_residual_W", value=max_res, target="< 0.01 W", status=max_res < 0.01 ? "pass" : "diagnostic"))
        
        return rows
    end

    function write_invariants_v46(path, rows)
        open(path, "w") do io
            println(io, "metric,value,target,status")
            for r in rows
                println(io, join((r.metric, r.value, r.target, r.status), ','))
            end
        end
        return path
    end

    const TEMPERATURE_AXIS_LIMITS_v46 = (250.0, 1400.0)
    const RUNNER_SENSOR_COLORS_v46 = Dict(
        :T8 => :blue,
        :T12_perim => :red,
        :T11_perim => :green,
        :T9_core => :orange,
        :T10_core => :brown,
        :T3 => :purple,
        :T2 => :black,
    )

    function steady_comparison_plot_v46(steady_results; title_suffix="fitted")
        sensors = (:T8, :T12_perim, :T11_perim, :T9_core, :T10_core, :T3, :T2)
        plot_object = plot(
            title="1D_v46 $title_suffix steady-state comparison",
            xlabel="Model temperature [K]",
            ylabel="Experiment temperature [K]",
            legend=:outerright,
            aspect_ratio=:equal,
            xlims=TEMPERATURE_AXIS_LIMITS_v46,
            ylims=TEMPERATURE_AXIS_LIMITS_v46,
            size=(700, 600),
        )
        for sensor in sensors
            model_values = [getproperty(row, Symbol(sensor, "_model")) for row in steady_results]
            experiment_values = [getproperty(row, Symbol(sensor, "_experiment")) for row in steady_results]
            scatter!(plot_object, model_values, experiment_values;
                     label=String(sensor), color=RUNNER_SENSOR_COLORS_v46[sensor],
                     ms=4, markerstrokewidth=0)
        end
        lower, upper = TEMPERATURE_AXIS_LIMITS_v46
        plot!(plot_object, [lower, upper], [lower, upper];
              label="1:1 line", color=:gray, linestyle=:dash, lw=1.5)
        return plot_object
    end

    function axial_profile_plot_v46(simulation_id, params=pnew_v46;
                                   nodes=RUNNER_PLOT_NODES_v46,
                                   variant_name="v46")
        outputs, result, experiment = solve_case_v46(params, simulation_id; nodes=nodes)
        conditions = simulation_conditions[string(simulation_id)]

        perim_x_mm = vcat(result.z_solid, result.z_rear) .* 1000.0
        perim_T = vcat(result.perim_temperature[:, end], result.rear_temperature[:, end])
        core_x_mm = result.z_solid .* 1000.0
        core_T = result.core_temperature[:, end]
        gas_x_mm = result.z_gas .* 1000.0
        gas_T = result.gas_temperature[:, end]
        rear_x_mm = result.z_solid .* 1000.0
        rear_T = result.rear_reservoir_temperature[:, end]

        solid_measurement_points = (
            (:T8, sensor_positions[:T8], "_T8", :square, :blue),
            (:T9_core, sensor_positions[:T9], "_T9", :circle, :orange),
            (:T12_perim, sensor_positions[:T9], "_T12", :utriangle, :red),
            (:T10_core, sensor_positions[:T10], "_T10", :diamond, :brown),
            (:T11_perim, sensor_positions[:T10], "_T11", :dtriangle, :green),
        )

        scale = absorbed_power_scale_v46(conditions[Io], params)
        title_text = "1D_v46 $variant_name axial profile: $simulation_id"
        subtitle_text = "flow=$(round(conditions[qlpm]; digits=2)) LPM, " *
                        "Io=$(round(conditions[Io] / 1000.0; digits=1)) kW/m^2, " *
                        "scale=$(round(scale; digits=3))"
        plot_object = plot(
            title=title_text * "\n" * subtitle_text,
            xlabel="Axial Position [mm]",
            ylabel="Temperature [K]",
            legend=:outerright,
            grid=true,
            xlims=(0.0, 1000.0 * (L + REAR_TUBE_LENGTH_v46)),
            ylims=TEMPERATURE_AXIS_LIMITS_v46,
            size=(850, 500),
        )
        plot!(plot_object, perim_x_mm, perim_T;
              label="perimeter/tube model", color=:red, linestyle=:dash, lw=2)
        plot!(plot_object, core_x_mm, core_T;
              label="core solid model", color=:orange, linestyle=:dot, lw=2)
        plot!(plot_object, gas_x_mm, gas_T;
              label="gas model", color=:green, lw=2)
        plot!(plot_object, rear_x_mm, rear_T;
              label="rear rail model", color=:purple,
              linestyle=:dashdot, lw=2)
        
        data = measurements
        id_str = string(simulation_id)
        for (sensor, position, observation_id, marker, color) in solid_measurement_points
            scatter!(plot_object, [position * 1000.0],
                     [observation(data, id_str, observation_id)[end]];
                     label=String(sensor), color=color, markershape=marker,
                     markerstrokewidth=0, ms=6)
        end
        scatter!(plot_object, [T3_SAMPLE_POSITION_v46 * 1000.0], [experiment[end, 6]];
                 label="T3 experiment", color=:purple, markershape=:star5,
                 markerstrokewidth=0, ms=7)
        scatter!(plot_object, [CAVITY_LENGTH_v46 * 1000.0], [experiment[end, 7]];
                 label="T2 experiment", color=:black, markershape=:diamond,
                 markerstrokewidth=0, ms=6)
        scatter!(plot_object, [CAVITY_LENGTH_v46 * 1000.0], [outputs[end, 7]];
                 label="T2 model", color=:black, markershape=:xcross,
                 markerstrokewidth=1.5, ms=6)
        vline!(plot_object, [L * 1000.0, (L + REAR_TUBE_CAVITY_LENGTH_v46) * 1000.0];
               label=false, color=:gray, linestyle=:dot, lw=1)
        return plot_object
    end

    function plot_transient_cases_v46(params=pnew_v46; keys=vcat(collect(keys_heating), collect(keys_cooling)), nodes=RUNNER_PLOT_NODES_v46)
        RUNNER_WRITE_PLOTS_v46 || return String[]
        paths = String[]
        for simulation_id in keys
            is_cool = simulation_id in keys_cooling
            outputs, result, experiment = solve_case_v46(params, simulation_id; is_cooling=is_cool, nodes=nodes)
            time_m = result.time ./ 60.0
            
            p_plot = plot(
                title="1D_v46 $(is_cool ? "Cooling" : "Heating") $(simulation_id)",
                xlabel="Time [min]",
                ylabel="Temperature [K]",
                legend=:outerright,
                size=(900, 550),
            )
            
            signals = [
                (:T8, outputs[:, 1], experiment[:, 1], :blue),
                (:T12, outputs[:, 2], experiment[:, 2], :red),
                (:T11, outputs[:, 3], experiment[:, 3], :green),
                (:T9, outputs[:, 4], experiment[:, 4], :orange),
                (:T10, outputs[:, 5], experiment[:, 5], :brown),
                (:T3, outputs[:, 6], experiment[:, 6], :purple),
                (:T2, outputs[:, 7], experiment[:, 7], :black),
            ]
            
            for (name, model, exp_v, col) in signals
                plot!(p_plot, time_m, exp_v, label="$(name) Exp", color=col, linestyle=:dash, lw=1.5)
                plot!(p_plot, time_m, model, label="$(name) Model", color=col, linestyle=:solid, lw=2.0)
            end
            
            file_name = "transients_$(simulation_id)_1D_v46.png"
            path = joinpath(RUNNER_TRANSIENT_DIR_v46, file_name)
            savefig(p_plot, path)
            push!(paths, path)
        end
        return paths
    end

    function plot_axial_profiles_v46(params=pnew_v46; keys=keys_heating, nodes=RUNNER_PLOT_NODES_v46)
        RUNNER_WRITE_PLOTS_v46 || return String[]
        paths = String[]
        for simulation_id in keys
            p_plot = axial_profile_plot_v46(simulation_id, params; nodes=nodes)
            file_name = "axial_profile_$(simulation_id)_fitted_energy_accounting_1D_v46.png"
            path = joinpath(RUNNER_AXIAL_DIR_v46, file_name)
            savefig(p_plot, path)
            push!(paths, path)
        end
        return paths
    end

    function write_cell_diagnostics_v46(params=pnew_v46; keys=keys_heating, nodes=RUNNER_PLOT_NODES_v46)
        for simulation_id in keys
            id_str = string(simulation_id)
            outputs, result, experiment = solve_case_v46(params, id_str; nodes=nodes)
            conditions = simulation_conditions[id_str]
            irradiance = conditions[Io]
            M = absorbed_power_scale_v46(irradiance, params)
            Q_delivered = M * irradiance * A_frt
            chi = core_source_fraction_v46(params)
            core_absorbed = chi * Q_delivered
            perim_absorbed = (1.0 - chi) * Q_delivered
            dx = L / length(result.z_solid)
            f_core_rear = rear_core_fraction_v46(params)
            G_receiver_rear = params[17]
            Trear = result.rear_reservoir_temperature[:, end]
            Tcavity = result.cavity_temperature[end]

            path = joinpath(RUNNER_AXIAL_DIR_v46, "cell_diagnostics_$(id_str)_1D_v46.csv")
            open(path, "w") do io
                println(io, "simulation_id,z_mm,Tperim_K,Tcore_K,Trear_K,Tgas_in_K,Tgas_out_K,h_W_m2K,Qgas_W,Qcore_source_W,Qperim_source_W")
                for i in eachindex(result.z_solid)
                    println(io, join((
                        id_str,
                        result.z_solid[i] * 1000.0,
                        result.perim_temperature[i, end],
                        result.core_temperature[i, end],
                        Trear[i],
                        result.gas_temperature[i, end],
                        result.gas_temperature[i + 1, end],
                        result.heat_transfer_coefficient[i, end],
                        result.receiver_gas_heat[i, end],
                        core_absorbed * result.core_source_weights[i],
                        i == 1 ? perim_absorbed : 0.0,
                    ), ','))
                end
            end
        end
    end
end

begin # optimization driver
    function run_calibration_v46(; stage=RUNNER_FIT_STAGE_v46, max_time=RUNNER_FIT_SECONDS_v46,
                                nodes=RUNNER_FIT_NODES_v46, dataset=RUNNER_CALIBRATION_DATASET_v46)
        println("[run_1D_v46] Starting Option A calibration on heating runs (stage=$(stage), dataset=$(dataset), nodes=$(nodes), max_time=$(max_time) s)...")
        res = calibrate_v46(
            nodes=nodes,
            maximum_iterations=RUNNER_FIT_ITERATIONS_v46,
            maximum_time_seconds=max_time,
            stage=stage,
            dataset=dataset,
        )
        println("[run_1D_v46] Optimization completed with status: $(res.retcode)")
        println("[run_1D_v46] Final objective: $(res.objective)")
        return (params=res.parameters, objective=res.objective, retcode=res.retcode)
    end
end

# Main execution entry
function main_1D_v46()
    println("============================================================================")
    println("1D_v46 Calibration and Analysis Pipeline")
    println("============================================================================")

    run_cal = runner_bool_v46("RECEIVER1D_v46_RUN_CALIBRATION", false)
    
    variant = if run_cal
        run_calibration_v46()
    else
        println("[run_1D_v46] Using baseline v46 parameter vector (calibration skipped).")
        (params=pnew_v46, objective=calibration_loss_v46(pnew_v46), retcode=:Manual)
    end

    # Write summary
    sum_path = joinpath(RUNNER_OUTPUT_DIR_v46, "optimization_summary_1D_v46.txt")
    open(sum_path, "w") do io
        println(io, "objective=$(variant.objective)")
        println(io, "return_code=$(variant.retcode)")
        println(io, "stage=$(RUNNER_FIT_STAGE_v46)")
        println(io, "calibration_dataset=$(RUNNER_CALIBRATION_DATASET_v46)")
        println(io, "parameters=$(variant.params)")
    end

    # Write parameters
    param_path = joinpath(RUNNER_OUTPUT_DIR_v46, "parameters_fitted_energy_accounting_1D_v46.csv")
    write_parameters_v46(param_path, variant)

    # Compute & write metrics
    println("[run_1D_v46] Computing multi-signal metrics...")
    metrics = compute_metrics_v46(variant.params)
    write_metrics_v46(joinpath(RUNNER_OUTPUT_DIR_v46, "analysis_results_fitted_energy_accounting_1D_v46.csv"), metrics)

    # Build & write steady results
    println("[run_1D_v46] Computing steady-state results and conservation residuals...")
    steady_results = build_steady_results_v46(variant.params)
    write_steady_results_v46(joinpath(RUNNER_OUTPUT_DIR_v46, "steady_results_fitted_energy_accounting_1D_v46.csv"), steady_results)

    # Compute flow slopes
    println("[run_1D_v46] Computing flow slopes...")
    slopes = compute_flow_slopes_v46(steady_results)
    write_flow_slopes_v46(joinpath(RUNNER_OUTPUT_DIR_v46, "flow_slopes_fitted_energy_accounting_1D_v46.csv"), slopes)

    # Compute invariants
    println("[run_1D_v46] Computing invariants...")
    inv_rows = compute_invariants_v46(variant.params, steady_results)
    write_invariants_v46(joinpath(RUNNER_OUTPUT_DIR_v46, "invariant_summary_1D_v46.csv"), inv_rows)

    # Plotting & Diagnostics
    if RUNNER_WRITE_PLOTS_v46
        println("[run_1D_v46] Generating transient plots...")
        plot_transient_cases_v46(variant.params)

        println("[run_1D_v46] Generating axial temperature profile plots...")
        plot_axial_profiles_v46(variant.params)

        println("[run_1D_v46] Generating steady-state comparison parity plot...")
        p_steady = steady_comparison_plot_v46(steady_results; title_suffix="fitted")
        savefig(p_steady, joinpath(RUNNER_PLOT_DIR_v46, "steady_comparison_fitted_energy_accounting_1D_v46.png"))

        println("[run_1D_v46] Writing cell-by-cell axial diagnostics CSVs...")
        write_cell_diagnostics_v46(variant.params)
    end

    println("[run_1D_v46] Pipeline completed successfully. Output in: $RUNNER_OUTPUT_DIR_v46")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_1D_v46()
end
