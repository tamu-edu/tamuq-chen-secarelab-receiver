# ============================================================================
# run_1D_v47.jl - Calibration, Validation, and Post-Processing for 1D_v47
# ============================================================================

using DifferentialEquations
using Plots
using StatsPlots
using LinearAlgebra
using Statistics
using DelimitedFiles
using Optimization
using OptimizationNLopt

include("1D_v47.jl")

begin # configuration & environment options
    function runner_env_v47(key, default)
        haskey(ENV, key) ? ENV[key] : default
    end

    function runner_float_v47(key, default)
        val = runner_env_v47(key, "")
        isempty(val) ? default : parse(Float64, val)
    end

    function runner_int_v47(key, default)
        val = runner_env_v47(key, "")
        isempty(val) ? default : parse(Int, val)
    end

    function runner_bool_v47(key, default)
        val = lowercase(strip(runner_env_v47(key, "")))
        isempty(val) ? default : val in ("1", "true", "yes", "on")
    end

    const RUNNER_OUTPUT_DIR_v47 = normpath(joinpath(@__DIR__, "summaries", "1D_v47"))
    const RUNNER_PLOT_DIR_v47 = joinpath(RUNNER_OUTPUT_DIR_v47, "plots")
    const RUNNER_TRANSIENT_DIR_v47 = joinpath(RUNNER_PLOT_DIR_v47, "transients")
    const RUNNER_AXIAL_DIR_v47 = joinpath(RUNNER_PLOT_DIR_v47, "axial_profiles")

    const RUNNER_CALIBRATION_DATASET_v47 = Symbol(runner_env_v47("RECEIVER1D_v47_CALIBRATION_DATASET", "heating"))
    const RUNNER_FIT_STAGE_v47 = Symbol(runner_env_v47("RECEIVER1D_v47_FIT_STAGE", "full"))
    const RUNNER_FIT_NODES_v47 = runner_int_v47("RECEIVER1D_v47_FIT_NODES", 15)
    const RUNNER_FIT_SECONDS_v47 = runner_float_v47("RECEIVER1D_v47_FIT_SECONDS", 600.0)
    const RUNNER_FIT_ITERATIONS_v47 = runner_int_v47("RECEIVER1D_v47_FIT_ITERATIONS", 500)
    const RUNNER_PLOT_NODES_v47 = runner_int_v47("RECEIVER1D_v47_PLOT_NODES", 25)
    const RUNNER_WRITE_PLOTS_v47 = runner_bool_v47("RECEIVER1D_v47_WRITE_PLOTS", true)

    mkpath(RUNNER_OUTPUT_DIR_v47)
    mkpath(RUNNER_PLOT_DIR_v47)
    mkpath(RUNNER_TRANSIENT_DIR_v47)
    mkpath(RUNNER_AXIAL_DIR_v47)
end

begin # analysis metrics & reporting
    const SIGNAL_NAMES_v47 = [:T8, :T12_perim, :T11_perim, :T9_core, :T10_core, :T3, :T2]

    function calculate_t90_error_v47(time, model, exp_v)
        exp_range = maximum(exp_v) - minimum(exp_v)
        exp_range <= 5.0 && return 0.0
        t90_exp_val = minimum(exp_v) + 0.90 * exp_range
        idx_exp = findfirst(>=(t90_exp_val), exp_v)
        idx_exp === nothing && return 0.0
        t90_exp = time[idx_exp]

        mod_range = maximum(model) - minimum(model)
        mod_range <= 5.0 && return 0.0
        t90_mod_val = minimum(model) + 0.90 * mod_range
        idx_mod = findfirst(>=(t90_mod_val), model)
        idx_mod === nothing && return 0.0
        t90_mod = time[idx_mod]

        return t90_mod - t90_exp
    end

    function compute_metrics_v47(params; heating_keys=keys_heating, cooling_keys=keys_cooling, nodes=RUNNER_PLOT_NODES_v47)
        rows = NamedTuple[]
        data_heat = measurements
        data_cool = measurements_cooling

        for simulation_id in heating_keys
            id_str = string(simulation_id)
            outputs, result, experiment = solve_case_v47(params, id_str; is_cooling=false, nodes=nodes)
            time = result.time
            for (j, sensor) in enumerate(SIGNAL_NAMES_v47)
                m = outputs[:, j]
                e = experiment[:, j]
                rmse = sqrt(mean(abs2, m .- e))
                steady_err = m[end] - e[end]
                t90_err = calculate_t90_error_v47(time, m, e)
                shape_l = signal_loss_v47(time, m, e)
                push!(rows, (
                    simulation_id=id_str,
                    phase="heating",
                    sensor=String(sensor),
                    rmse_K=rmse,
                    steady_error_K=steady_err,
                    t90_error_s=t90_err,
                    shape_loss=shape_l,
                ))
            end
        end

        for simulation_id in cooling_keys
            id_str = string(simulation_id)
            outputs, result, experiment = solve_case_v47(params, id_str; is_cooling=true, nodes=nodes)
            time = result.time
            for (j, sensor) in enumerate(SIGNAL_NAMES_v47)
                m = outputs[:, j]
                e = experiment[:, j]
                rmse = sqrt(mean(abs2, m .- e))
                steady_err = m[end] - e[end]
                t90_err = calculate_t90_error_v47(time, m, e)
                shape_l = signal_loss_v47(time, m, e)
                push!(rows, (
                    simulation_id=id_str,
                    phase="cooling",
                    sensor=String(sensor),
                    rmse_K=rmse,
                    steady_error_K=steady_err,
                    t90_error_s=t90_err,
                    shape_loss=shape_l,
                ))
            end
        end

        return rows
    end

    function write_metrics_v47(path, rows)
        open(path, "w") do io
            println(io, "simulation_id,phase,sensor,rmse_K,steady_error_K,t90_error_s,shape_loss")
            for r in rows
                println(io, join((r.simulation_id, r.phase, r.sensor, r.rmse_K, r.steady_error_K, r.t90_error_s, r.shape_loss), ','))
            end
        end
        return path
    end

    function write_parameters_v47(path, variant)
        p = variant.params
        param_names = [
            "A_Nu", "B_Re", "C_Pr_fixed", "front_dep", "scale_456", "scale_304", "scale_256",
            "G_core_perim_W_m_K", "C_perim_eff_J_K", "k_perim_ref_W_m_K", "beta_opt",
            "chi", "f_core_rear", "flange_scale", "k_core_axial_scale", "C_rear_eff_J_K",
            "G_receiver_rear_W_K", "G_rear_tube_W_K", "G_rear_cavity_W_K", "G_rear_axial_W_K",
            "delta_web", "C_z", "h_suction"
        ]
        open(path, "w") do io
            println(io, "index,name,value")
            for i in eachindex(p)
                println(io, join((i, param_names[i], p[i]), ','))
            end
            println(io, "# derived,C_total_with_rear_eff_J_K,$(participating_total_heat_capacity_v47(p))")
            println(io, "# reference,measured_C_eff_J_K,$(MEASURED_ASSEMBLY_CAPACITANCE_v47)")
        end
        return path
    end

    function build_steady_results_v47(params; heating_keys=keys_heating, nodes=RUNNER_PLOT_NODES_v47)
        rows = NamedTuple[]
        data = measurements
        for simulation_id in heating_keys
            id_str = string(simulation_id)
            outputs, result, experiment = solve_case_v47(params, id_str; is_cooling=false, nodes=nodes)
            conditions = simulation_conditions[id_str]
            flow = conditions[qlpm]
            irradiance = conditions[Io]
            scale = absorbed_power_scale_v47(irradiance, params)

            operating = operating_condition_v3(
                irradiance=linear_history_v3(result.time, fill(irradiance, length(result.time))),
                flow_lpm=linear_history_v3(result.time, observation(data, id_str, "_flow")),
                inlet_temperature=linear_history_v3(result.time, observation(data, id_str, "_Tin")),
                ambient_temperature=linear_history_v3(result.time, observation(data, id_str, "_Tamb")),
            )
            balance = compute_energy_balance_v47(result, operating, result.time[end], params; nodes=nodes)

            # Area mean and heat-transfer-weighted average HTC
            h_local = result.heat_transfer_coefficient[:, end]
            T_core_end = result.core_temperature[:, end]
            T_gas_end = 0.5 .* (result.gas_temperature[1:nodes, end] .+ result.gas_temperature[2:nodes+1, end])
            dT_drive = abs.(T_core_end .- T_gas_end)
            sum_dT = sum(dT_drive)
            h_weighted = sum_dT > 1e-6 ? sum(h_local .* dT_drive) / sum_dT : mean(h_local)

            push!(rows, (
                simulation_id=id_str,
                flow_lpm=flow,
                irradiance=irradiance,
                power_scale=scale,
                delivered_power_W=balance.delivered_W,
                core_solar_W=balance.core_solar_W,
                perim_solar_W=balance.perim_solar_W,
                suction_gas_heat_W=balance.suction_gas_heat_W,
                core_gas_heat_W=balance.core_gas_heat_W,
                rear_tube_gas_heat_W=balance.rear_tube_gas_heat_W,
                gas_heat_W=balance.gas_heat_W,
                front_rad_loss_W=balance.front_rad_loss_W,
                cavity_amb_loss_W=balance.cavity_amb_loss_W,
                flange_loss_W=balance.flange_loss_W,
                dE_stored_dt_W=balance.dE_stored_dt_W,
                instantaneous_residual_W=balance.instantaneous_residual_W,
                steady_flux_gap_W=balance.steady_flux_gap_W,
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
                mean_htc=mean(h_local),
                weighted_mean_htc=h_weighted,
            ))
        end
        return rows
    end

    function write_steady_results_v47(path, rows)
        fields = propertynames(first(rows))
        open(path, "w") do io
            println(io, join(fields, ','))
            for r in rows
                println(io, join((getproperty(r, f) for f in fields), ','))
            end
        end
        return path
    end

    function compute_flow_slopes_v47(steady_results)
        irradiances = sort(unique(r.irradiance for r in steady_results), rev=true)
        rows = NamedTuple[]
        sensors = [
            (:T8, :T8_model, :T8_experiment),
            (:T12_perim, :T12_perim_model, :T12_perim_experiment),
            (:T11_perim, :T11_perim_model, :T11_perim_experiment),
            (:T9_core, :T9_core_model, :T9_core_experiment),
            (:T10_core, :T10_core_model, :T10_core_experiment),
            (:T3, :T3_model, :T3_experiment),
            (:T2, :T2_model, :T2_experiment),
        ]
        for irr in irradiances
            sub = filter(r -> r.irradiance == irr, steady_results)
            flows = [r.flow_lpm for r in sub]
            for (s_name, mod_f, exp_f) in sensors
                mod_vals = [getproperty(r, mod_f) for r in sub]
                exp_vals = [getproperty(r, exp_f) for r in sub]
                
                slope_mod = (cov(flows, mod_vals) / var(flows))
                slope_exp = (cov(flows, exp_vals) / var(flows))
                push!(rows, (
                    irradiance=irr,
                    signal=String(s_name),
                    slope_model_K_per_LPM=slope_mod,
                    slope_experiment_K_per_LPM=slope_exp,
                    error_K_per_LPM=slope_mod - slope_exp,
                ))
            end
        end
        return rows
    end

    function write_flow_slopes_v47(path, rows)
        open(path, "w") do io
            println(io, "irradiance,signal,slope_model_K_per_LPM,slope_experiment_K_per_LPM,error_K_per_LPM")
            for r in rows
                println(io, join((r.irradiance, r.signal, r.slope_model_K_per_LPM, r.slope_experiment_K_per_LPM, r.error_K_per_LPM), ','))
            end
        end
        return path
    end

    function compute_invariants_v47(p, steady_results)
        rows = NamedTuple[]
        C_tot = participating_total_heat_capacity_v47(p)
        push!(rows, (metric="C_total_with_rear", value=C_tot, target="301 +/- 23", status=abs(C_tot - 301.0) <= 23.0 ? "pass" : "diagnostic"))
        
        max_res = maximum(abs.(r.instantaneous_residual_W for r in steady_results))
        push!(rows, (metric="max_instantaneous_residual_W", value=max_res, target="< 0.0001 W", status=max_res < 0.0001 ? "pass" : "diagnostic"))
        
        return rows
    end

    function write_invariants_v47(path, rows)
        open(path, "w") do io
            println(io, "metric,value,target,status")
            for r in rows
                println(io, join((r.metric, r.value, r.target, r.status), ','))
            end
        end
        return path
    end
end

begin # plotting and visualization routines
    const TEMPERATURE_AXIS_LIMITS_v47 = (250.0, 1400.0)
    const RUNNER_SENSOR_COLORS_v47 = Dict(
        :T8 => :blue,
        :T12_perim => :red,
        :T11_perim => :green,
        :T9_core => :orange,
        :T10_core => :brown,
        :T3 => :purple,
        :T2 => :black,
    )

    function steady_comparison_plot_v47(steady_results; title_suffix="fitted")
        sensors = (:T8, :T12_perim, :T11_perim, :T9_core, :T10_core, :T3, :T2)
        plot_object = plot(
            title="1D_v47 $title_suffix steady-state comparison",
            xlabel="Model temperature [K]",
            ylabel="Experiment temperature [K]",
            legend=:outerright,
            aspect_ratio=:equal,
            xlims=TEMPERATURE_AXIS_LIMITS_v47,
            ylims=TEMPERATURE_AXIS_LIMITS_v47,
            size=(700, 600),
        )
        for sensor in sensors
            model_values = [getproperty(row, Symbol(sensor, "_model")) for row in steady_results]
            experiment_values = [getproperty(row, Symbol(sensor, "_experiment")) for row in steady_results]
            scatter!(plot_object, model_values, experiment_values;
                     label=String(sensor), color=RUNNER_SENSOR_COLORS_v47[sensor],
                     ms=4, markerstrokewidth=0)
        end
        lower, upper = TEMPERATURE_AXIS_LIMITS_v47
        plot!(plot_object, [lower, upper], [lower, upper];
              label="1:1 line", color=:gray, linestyle=:dash, lw=1.5)
        return plot_object
    end

    function axial_profile_plot_v47(simulation_id, params=pnew_v47;
                                   nodes=RUNNER_PLOT_NODES_v47,
                                   variant_name="v47")
        outputs, result, experiment = solve_case_v47(params, simulation_id; nodes=nodes)
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

        scale = absorbed_power_scale_v47(conditions[Io], params)
        title_text = "1D_v47 $variant_name axial profile: $simulation_id"
        subtitle_text = "flow=$(round(conditions[qlpm]; digits=2)) LPM, " *
                        "Io=$(round(conditions[Io] / 1000.0; digits=1)) kW/m^2, " *
                        "scale=$(round(scale; digits=3))"
        plot_object = plot(
            title=title_text * "\n" * subtitle_text,
            xlabel="Axial Position [mm]",
            ylabel="Temperature [K]",
            legend=:outerright,
            grid=true,
            xlims=(0.0, 1000.0 * (L + REAR_TUBE_LENGTH_v47)),
            ylims=TEMPERATURE_AXIS_LIMITS_v47,
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
        scatter!(plot_object, [T3_SAMPLE_POSITION_v47 * 1000.0], [experiment[end, 6]];
                 label="T3 experiment", color=:purple, markershape=:star5,
                 markerstrokewidth=0, ms=7)
        scatter!(plot_object, [CAVITY_LENGTH_v47 * 1000.0], [experiment[end, 7]];
                 label="T2 experiment", color=:black, markershape=:diamond,
                 markerstrokewidth=0, ms=6)
        scatter!(plot_object, [CAVITY_LENGTH_v47 * 1000.0], [outputs[end, 7]];
                 label="T2 model", color=:black, markershape=:xcross,
                 markerstrokewidth=1.5, ms=6)
        vline!(plot_object, [L * 1000.0, (L + REAR_TUBE_CAVITY_LENGTH_v47) * 1000.0];
               label=false, color=:gray, linestyle=:dot, lw=1)
        return plot_object
    end

    function plot_transient_cases_v47(params=pnew_v47; keys=vcat(collect(keys_heating), collect(keys_cooling)), nodes=RUNNER_PLOT_NODES_v47)
        RUNNER_WRITE_PLOTS_v47 || return String[]
        paths = String[]
        for simulation_id in keys
            is_cool = simulation_id in keys_cooling
            outputs, result, experiment = solve_case_v47(params, simulation_id; is_cooling=is_cool, nodes=nodes)
            time_m = result.time ./ 60.0
            
            p_plot = plot(
                title="1D_v47 $(is_cool ? "Cooling" : "Heating") $(simulation_id)",
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
            
            file_name = "transients_$(simulation_id)_1D_v47.png"
            path = joinpath(RUNNER_TRANSIENT_DIR_v47, file_name)
            savefig(p_plot, path)
            push!(paths, path)
        end
        return paths
    end

    function plot_axial_profiles_v47(params=pnew_v47; keys=keys_heating, nodes=RUNNER_PLOT_NODES_v47)
        RUNNER_WRITE_PLOTS_v47 || return String[]
        paths = String[]
        for simulation_id in keys
            p_plot = axial_profile_plot_v47(simulation_id, params; nodes=nodes)
            file_name = "axial_profile_$(simulation_id)_fitted_energy_accounting_1D_v47.png"
            path = joinpath(RUNNER_AXIAL_DIR_v47, file_name)
            savefig(p_plot, path)
            push!(paths, path)
        end
        return paths
    end

    function write_cell_diagnostics_v47(params=pnew_v47; keys=keys_heating, nodes=RUNNER_PLOT_NODES_v47)
        for simulation_id in keys
            id_str = string(simulation_id)
            outputs, result, experiment = solve_case_v47(params, id_str; nodes=nodes)
            conditions = simulation_conditions[id_str]
            irradiance = conditions[Io]
            M = absorbed_power_scale_v47(irradiance, params)
            Q_delivered = M * irradiance * A_frt
            chi = core_source_fraction_v47(params)
            core_absorbed = chi * Q_delivered
            perim_absorbed = (1.0 - chi) * Q_delivered
            dx = L / length(result.z_solid)
            Trear = result.rear_reservoir_temperature[:, end]

            path = joinpath(RUNNER_AXIAL_DIR_v47, "cell_diagnostics_$(id_str)_1D_v47.csv")
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
    function run_calibration_v47(; stage=RUNNER_FIT_STAGE_v47, max_time=RUNNER_FIT_SECONDS_v47,
                                nodes=RUNNER_FIT_NODES_v47, dataset=RUNNER_CALIBRATION_DATASET_v47)
        println("[run_1D_v47] Starting Option A calibration on heating runs (stage=$(stage), dataset=$(dataset), nodes=$(nodes), max_time=$(max_time) s)...")
        res = calibrate_v47(
            nodes=nodes,
            maximum_iterations=RUNNER_FIT_ITERATIONS_v47,
            maximum_time_seconds=max_time,
            stage=stage,
            dataset=dataset,
        )
        println("[run_1D_v47] Optimization completed with status: $(res.retcode)")
        println("[run_1D_v47] Final objective: $(res.objective)")
        return (params=res.parameters, objective=res.objective, retcode=res.retcode)
    end
end

# Main execution entry
function main_1D_v47()
    println("============================================================================")
    println("1D_v47 Calibration and Analysis Pipeline")
    println("============================================================================")

    run_cal = runner_bool_v47("RECEIVER1D_v47_RUN_CALIBRATION", true)
    
    variant = if run_cal
        run_calibration_v47()
    else
        println("[run_1D_v47] Using baseline v47 parameter vector.")
        (params=pnew_v47, objective=calibration_loss_v47(pnew_v47), retcode=:Manual)
    end

    # Write summary
    sum_path = joinpath(RUNNER_OUTPUT_DIR_v47, "optimization_summary_1D_v47.txt")
    open(sum_path, "w") do io
        println(io, "objective=$(variant.objective)")
        println(io, "return_code=$(variant.retcode)")
        println(io, "stage=$(RUNNER_FIT_STAGE_v47)")
        println(io, "calibration_dataset=$(RUNNER_CALIBRATION_DATASET_v47)")
        println(io, "parameters=$(variant.params)")
    end

    # Write parameters
    param_path = joinpath(RUNNER_OUTPUT_DIR_v47, "parameters_fitted_energy_accounting_1D_v47.csv")
    write_parameters_v47(param_path, variant)

    # Compute & write metrics
    println("[run_1D_v47] Computing multi-signal metrics...")
    metrics = compute_metrics_v47(variant.params)
    write_metrics_v47(joinpath(RUNNER_OUTPUT_DIR_v47, "analysis_results_fitted_energy_accounting_1D_v47.csv"), metrics)

    # Build & write steady results
    println("[run_1D_v47] Computing steady-state results and conservation residuals...")
    steady_results = build_steady_results_v47(variant.params)
    write_steady_results_v47(joinpath(RUNNER_OUTPUT_DIR_v47, "steady_results_fitted_energy_accounting_1D_v47.csv"), steady_results)

    # Compute flow slopes
    println("[run_1D_v47] Computing flow slopes...")
    slopes = compute_flow_slopes_v47(steady_results)
    write_flow_slopes_v47(joinpath(RUNNER_OUTPUT_DIR_v47, "flow_slopes_fitted_energy_accounting_1D_v47.csv"), slopes)

    # Compute invariants
    println("[run_1D_v47] Computing invariants...")
    inv_rows = compute_invariants_v47(variant.params, steady_results)
    write_invariants_v47(joinpath(RUNNER_OUTPUT_DIR_v47, "invariant_summary_1D_v47.csv"), inv_rows)

    # Plotting & Diagnostics
    if RUNNER_WRITE_PLOTS_v47
        println("[run_1D_v47] Generating transient plots...")
        plot_transient_cases_v47(variant.params)

        println("[run_1D_v47] Generating axial temperature profile plots...")
        plot_axial_profiles_v47(variant.params)

        println("[run_1D_v47] Generating steady-state comparison parity plot...")
        p_steady = steady_comparison_plot_v47(steady_results; title_suffix="fitted")
        savefig(p_steady, joinpath(RUNNER_PLOT_DIR_v47, "steady_comparison_fitted_energy_accounting_1D_v47.png"))

        println("[run_1D_v47] Writing cell-by-cell axial diagnostics CSVs...")
        write_cell_diagnostics_v47(variant.params)
    end

    println("[run_1D_v47] Pipeline completed successfully. Output in: $RUNNER_OUTPUT_DIR_v47")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_1D_v47()
end
