# run_1D_v49.jl
# Calibration, Verification, and Post-Processing Runner for Entire Converter Model 1D_v49.

begin # imports & directory setups
    using Plots
    using Statistics
    using DelimitedFiles
    using Dates

    include("1D_v49.jl")

    const RUNNER_OUTPUT_DIR_v49 = joinpath(@__DIR__, "summaries", "1D_v49")
    const RUNNER_PLOT_DIR_v49 = joinpath(RUNNER_OUTPUT_DIR_v49, "plots")
    const RUNNER_TRANSIENT_DIR_v49 = joinpath(RUNNER_PLOT_DIR_v49, "transients")
    const RUNNER_AXIAL_DIR_v49 = joinpath(RUNNER_PLOT_DIR_v49, "axial_profiles")

    mkpath(RUNNER_OUTPUT_DIR_v49)
    mkpath(RUNNER_PLOT_DIR_v49)
    mkpath(RUNNER_TRANSIENT_DIR_v49)
    mkpath(RUNNER_AXIAL_DIR_v49)

    function runner_env_v49(key, default)
        haskey(ENV, key) ? ENV[key] : default
    end

    function runner_bool_v49(key, default)
        val = runner_env_v49(key, string(default))
        return lowercase(strip(val)) in ("1", "true", "yes", "on")
    end

    function runner_float_v49(key, default)
        parse(Float64, runner_env_v49(key, string(default)))
    end

    function runner_int_v49(key, default)
        parse(Int, runner_env_v49(key, string(default)))
    end

    const RUNNER_FIT_STAGE_v49 = runner_env_v49("RECEIVER1D_v49_FIT_STAGE", "full")
    const RUNNER_FIT_SECONDS_v49 = runner_float_v49("RECEIVER1D_v49_FIT_SECONDS", 600.0)
    const RUNNER_FIT_ITERATIONS_v49 = runner_int_v49("RECEIVER1D_v49_FIT_ITERATIONS", 500)
    const RUNNER_FIT_NODES_v49 = runner_int_v49("RECEIVER1D_v49_FIT_NODES", 15)
    const RUNNER_CALIBRATION_DATASET_v49 = runner_env_v49("RECEIVER1D_v49_CALIBRATION_DATASET", "heating")
    const RUNNER_WRITE_PLOTS_v49 = runner_bool_v49("RECEIVER1D_v49_WRITE_PLOTS", true)
end

begin # metrics & analysis helpers
    const SENSOR_NAMES_v49 = ["T8", "T12_perim", "T11_perim", "T9_core", "T10_core", "T3", "T2"]

    function compute_t90_v49(time, trace; is_cooling=false)
        t0 = time[1]
        t_end = time[end]
        v0 = trace[1]
        v_end = trace[end]
        
        if is_cooling
            # For cooling (decaying signal), t90 is when signal drops 90% towards endpoint
            v_target = v0 - 0.90 * (v0 - v_end)
            for i in eachindex(trace)
                if trace[i] <= v_target
                    return time[i] - t0
                end
            end
        else
            # For heating (rising signal), t90 is when signal rises 90% towards endpoint
            v_target = v0 + 0.90 * (v_end - v0)
            for i in eachindex(trace)
                if trace[i] >= v_target
                    return time[i] - t0
                end
            end
        end
        return t_end - t0
    end

    function evaluate_case_metrics_v49(p, simulation_id; is_cooling=false, nodes=15)
        id_str = string(simulation_id)
        outputs, result, experiment = solve_case_v49(p, id_str; is_cooling=is_cooling, nodes=nodes)
        data = is_cooling ? measurements_cooling : measurements
        time = observation_time(data, id_str)
        
        metrics = []
        for j in 1:7
            m_trace = outputs[:, j]
            e_trace = experiment[:, j]
            
            rmse = sqrt(mean(abs2, m_trace .- e_trace))
            steady_err = m_trace[end] - e_trace[end]
            
            t90_m = compute_t90_v49(time, m_trace; is_cooling=is_cooling)
            t90_e = compute_t90_v49(time, e_trace; is_cooling=is_cooling)
            t90_err = t90_m - t90_e
            
            shape = signal_loss_v49(time, m_trace, e_trace)
            
            push!(metrics, (
                simulation_id=id_str,
                phase=is_cooling ? "cooling" : "heating",
                sensor=SENSOR_NAMES_v49[j],
                rmse_K=rmse,
                steady_error_K=steady_err,
                t90_error_s=t90_err,
                shape_loss=shape,
            ))
        end
        return metrics, outputs, result, experiment
    end

    function compute_metrics_v49(p; nodes=15)
        all_metrics = []
        for id in keys_heating_v49
            m, _, _, _ = evaluate_case_metrics_v49(p, id; is_cooling=false, nodes=nodes)
            append!(all_metrics, m)
        end
        for id in keys_cooling_v49
            m, _, _, _ = evaluate_case_metrics_v49(p, id; is_cooling=true, nodes=nodes)
            append!(all_metrics, m)
        end
        return all_metrics
    end
end

begin # steady state accounting & reporting
    function build_steady_results_v49(p; nodes=15)
        rows = []
        for id in keys_heating_v49
            id_str = string(id)
            conditions, time, experiment, data = experimental_case_v49(id_str)
            flow = conditions[qlpm]
            irradiance = conditions[Io]
            
            operating = operating_condition_v3(
                irradiance=linear_history_v3(time, fill(irradiance, length(time))),
                flow_lpm=linear_history_v3(time, observation(data, id_str, "_flow")),
                inlet_temperature=linear_history_v3(time, observation(data, id_str, "_Tin")),
                ambient_temperature=linear_history_v3(time, observation(data, id_str, "_Tamb")),
            )
            
            outputs, result, _ = solve_case_v49(p, id_str; is_cooling=false, nodes=nodes)
            t_end = time[end]
            balance = compute_energy_balance_v49(result, operating, t_end, p; nodes=nodes)
            
            push!(rows, (
                simulation_id=id_str,
                flow_lpm=flow,
                irradiance=irradiance,
                power_scale=absorbed_power_scale_v49(irradiance, p),
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
                macro_htc=balance.macro_htc,
            ))
        end
        return rows
    end
end

begin # flow slope calculation
    function compute_flow_slopes_v49(steady_rows)
        # Group by irradiance cluster
        clusters = Dict{Float64, Vector{Any}}()
        for r in steady_rows
            irr = r.irradiance
            if !haskey(clusters, irr)
                clusters[irr] = []
            end
            push!(clusters[irr], r)
        end

        slopes = []
        for (irr, group) in clusters
            flows = [g.flow_lpm for g in group]
            if length(flows) >= 2
                for (sig_idx, sig_name) in enumerate(SENSOR_NAMES_v49)
                    m_vals = [sig_idx == 1 ? g.T8_model :
                              sig_idx == 2 ? g.T12_perim_model :
                              sig_idx == 3 ? g.T11_perim_model :
                              sig_idx == 4 ? g.T9_core_model :
                              sig_idx == 5 ? g.T10_core_model :
                              sig_idx == 6 ? g.T3_model : g.T2_model for g in group]
                    e_vals = [sig_idx == 1 ? g.T8_experiment :
                              sig_idx == 2 ? g.T12_perim_experiment :
                              sig_idx == 3 ? g.T11_perim_experiment :
                              sig_idx == 4 ? g.T9_core_experiment :
                              sig_idx == 5 ? g.T10_core_experiment :
                              sig_idx == 6 ? g.T3_experiment : g.T2_experiment for g in group]
                    
                    # Linear regression slope: cov(flows, vals) / var(flows)
                    f_mean = mean(flows)
                    m_mean = mean(m_vals)
                    e_mean = mean(e_vals)
                    denom = sum((flows .- f_mean).^2)
                    
                    slope_m = sum((flows .- f_mean) .* (m_vals .- m_mean)) / denom
                    slope_e = sum((flows .- f_mean) .* (e_vals .- e_mean)) / denom
                    
                    push!(slopes, (
                        irradiance=irr,
                        signal=sig_name,
                        slope_model_K_per_LPM=slope_m,
                        slope_experiment_K_per_LPM=slope_e,
                        error_K_per_LPM=slope_m - slope_e,
                    ))
                end
            end
        end
        return slopes
    end
end

begin # invariant checking
    function compute_invariants_v49(p, steady_rows)
        invs = []
        C_tot = total_participating_heat_capacity_v49(p; nodes=15)
        push!(invs, (metric="C_total_with_rear", value=C_tot, target="301 +/- 23", status=abs(C_tot - 301.0) <= 23.0 ? "pass" : "fail"))
        
        max_res = maximum(abs(r.instantaneous_residual_W) for r in steady_rows)
        push!(invs, (metric="max_instantaneous_residual_W", value=max_res, target="< 0.0001 W", status=max_res < 1e-4 ? "pass" : "fail"))
        
        return invs
    end
end

begin # table export writers
    function write_parameters_v49(path, variant)
        p = variant.params
        names = [
            "C1_Nu", "C2_Nu", "Nu_inf", "front_dep",
            "scale_456", "scale_304", "scale_256",
            "G_core_perim_W_m_K", "C_perim_eff_J_K", "k_perim_ref_W_m_K",
            "beta_opt", "chi", "f_core_rear", "flange_scale",
            "k_core_axial_scale", "C_rear_eff_J_K", "G_receiver_rear_W_K",
            "G_rear_tube_W_K", "beta_rad", "kA_rear_eff_W_m_K",
            "delta_web", "F_LoS", "h_suction",
            "h_probe_ref", "w_probe_rad", "G_probe_stem",
            "w10_stem", "w11_stem"
        ]
        open(path, "w") do io
            println(io, "index,name,value")
            for i in eachindex(p)
                println(io, "$i,$(names[i]),$(p[i])")
            end
            C_tot = total_participating_heat_capacity_v49(p; nodes=15)
            println(io, "# derived,C_total_with_rear_eff_J_K,$C_tot")
            println(io, "# reference,measured_C_eff_J_K,$MEASURED_ASSEMBLY_CAPACITANCE_v49")
        end
    end

    function write_metrics_v49(path, metrics)
        open(path, "w") do io
            println(io, "simulation_id,phase,sensor,rmse_K,steady_error_K,t90_error_s,shape_loss")
            for m in metrics
                println(io, "$(m.simulation_id),$(m.phase),$(m.sensor),$(m.rmse_K),$(m.steady_error_K),$(m.t90_error_s),$(m.shape_loss)")
            end
        end
    end

    function write_steady_results_v49(path, rows)
        open(path, "w") do io
            header = join(keys(rows[1]), ',')
            println(io, header)
            for r in rows
                println(io, join(values(r), ','))
            end
        end
    end

    function write_flow_slopes_v49(path, slopes)
        open(path, "w") do io
            println(io, "irradiance,signal,slope_model_K_per_LPM,slope_experiment_K_per_LPM,error_K_per_LPM")
            for s in slopes
                println(io, "$(s.irradiance),$(s.signal),$(s.slope_model_K_per_LPM),$(s.slope_experiment_K_per_LPM),$(s.error_K_per_LPM)")
            end
        end
    end

    function write_invariants_v49(path, invs)
        open(path, "w") do io
            println(io, "metric,value,target,status")
            for iv in invs
                println(io, "$(iv.metric),$(iv.value),$(iv.target),$(iv.status)")
            end
        end
    end
end

begin # plotting generators
    function plot_transient_cases_v49(p; nodes=15)
        for id in vcat(keys_heating_v49, keys_cooling_v49)
            is_cool = id in keys_cooling_v49
            id_str = string(id)
            outputs, result, experiment = solve_case_v49(p, id_str; is_cooling=is_cool, nodes=nodes)
            data = is_cool ? measurements_cooling : measurements
            time = observation_time(data, id_str)
            
            p_plot = plot(title="1D_v49 Transient: Case $(id_str) ($(is_cool ? "Cooling" : "Heating"))",
                          xlabel="Time [s]", ylabel="Temperature [K]", legend=:outertopright, size=(900, 600))
            
            colors = [:crimson, :darkorange, :goldenrod, :forestgreen, :teal, :mediumblue, :purple]
            labels = ["T8 (Front wall)", "T12 (Mid wall)", "T11 (Rear wall)", "T9 (Mid core)", "T10 (Rear core)", "T3 (Exit gas)", "T2 (Cavity)"]
            
            for j in 1:7
                plot!(p_plot, time, experiment[:, j], label="$(labels[j]) Exp", color=colors[j], linestyle=:dash, linewidth=1.5)
                plot!(p_plot, time, outputs[:, j], label="$(labels[j]) Model", color=colors[j], linestyle=:solid, linewidth=2.0)
            end
            
            out_file = joinpath(RUNNER_TRANSIENT_DIR_v49, "transients_$(id_str)_1D_v49.png")
            savefig(p_plot, out_file)
        end
    end

    function plot_axial_profiles_v49(p; nodes=25)
        for id in keys_heating_v49
            id_str = string(id)
            outputs, result, experiment = solve_case_v49(p, id_str; is_cooling=false, nodes=nodes)
            
            p_ax = plot(title="1D_v49 Axial Temperature Profile: Case $(id_str)",
                        xlabel="Axial Position z [mm]", ylabel="Steady-State Temperature [K]",
                        legend=:outertopright, size=(900, 550))
            
            z_s_mm = result.z_solid .* 1000.0
            z_g_mm = result.z_gas .* 1000.0
            
            plot!(p_ax, z_s_mm, result.core_temperature[:, end], label="Core Solid Matrix T_core(z)", color=:red, linewidth=2.5)
            plot!(p_ax, z_s_mm, result.perim_temperature[:, end], label="Perimeter Housing T_perim(z)", color=:orange, linewidth=2.0)
            plot!(p_ax, z_s_mm, result.rear_reservoir_temperature[:, end], label="Rear Hardware Rail T_rear(z)", color=:purple, linewidth=2.0, linestyle=:dot)
            plot!(p_ax, z_g_mm, result.gas_temperature[:, end], label="Gas Stream T_gas(z)", color=:blue, linewidth=2.5, linestyle=:dash)
            
            # Scatter experimental points at their axial positions
            exp_z = [sensor_positions[:T8], WALL_CHAIN_POSITIONS_v49.T12, WALL_CHAIN_POSITIONS_v49.T11,
                     sensor_positions[:T9], sensor_positions[:T10]] .* 1000.0
            exp_t = [experiment[end, 1], experiment[end, 2], experiment[end, 3], experiment[end, 4], experiment[end, 5]]
            scatter!(p_ax, exp_z, exp_t, label="Experimental Thermocouples", color=:black, markersize=6, markerstrokewidth=1.5)
            
            out_file = joinpath(RUNNER_AXIAL_DIR_v49, "axial_profile_$(id_str)_fitted_energy_accounting_1D_v49.png")
            savefig(p_ax, out_file)
        end
    end

    function steady_comparison_plot_v49(steady_results; title_suffix="fitted")
        p = scatter(
            title="1D_v49 Steady Parity Comparison ($(title_suffix))",
            xlabel="Experimental Temperature [K]",
            ylabel="Model Predicted Temperature [K]",
            legend=:topleft,
            aspect_ratio=:equal,
            size=(750, 700)
        )
        
        all_exp = Float64[]
        all_mod = Float64[]
        
        colors = [:crimson, :darkorange, :goldenrod, :forestgreen, :teal, :mediumblue, :purple]
        names = ["T8", "T12_perim", "T11_perim", "T9_core", "T10_core", "T3", "T2"]
        
        for (j, name) in enumerate(names)
            exp_v = [j == 1 ? r.T8_experiment :
                     j == 2 ? r.T12_perim_experiment :
                     j == 3 ? r.T11_perim_experiment :
                     j == 4 ? r.T9_core_experiment :
                     j == 5 ? r.T10_core_experiment :
                     j == 6 ? r.T3_experiment : r.T2_experiment for r in steady_results]
            mod_v = [j == 1 ? r.T8_model :
                     j == 2 ? r.T12_perim_model :
                     j == 3 ? r.T11_perim_model :
                     j == 4 ? r.T9_core_model :
                     j == 5 ? r.T10_core_model :
                     j == 6 ? r.T3_model : r.T2_model for r in steady_results]
            
            append!(all_exp, exp_v)
            append!(all_mod, mod_v)
            
            scatter!(p, exp_v, mod_v, label=name, color=colors[j], markersize=5, alpha=0.85)
        end
        
        min_val = min(minimum(all_exp), minimum(all_mod)) - 20.0
        max_val = max(maximum(all_exp), maximum(all_mod)) + 20.0
        
        plot!(p, [min_val, max_val], [min_val, max_val], label="1:1 Parity", color=:black, linestyle=:dash, linewidth=1.5)
        plot!(p, [min_val, max_val], [min_val + 50, max_val + 50], label="+/- 50 K Bounds", color=:gray, linestyle=:dot)
        plot!(p, [min_val, max_val], [min_val - 50, max_val - 50], label=false, color=:gray, linestyle=:dot)
        
        xlims!(p, min_val, max_val)
        ylims!(p, min_val, max_val)
        
        return p
    end

    function write_cell_diagnostics_v49(p; nodes=25)
        for id in keys_heating_v49
            id_str = string(id)
            outputs, result, experiment = solve_case_v49(p, id_str; is_cooling=false, nodes=nodes)
            conditions, time, _, data = experimental_case_v49(id_str)
            irradiance = conditions[Io]
            delivered = absorbed_power_scale_v49(irradiance, p) * irradiance * A_frt_v49
            chi = core_source_fraction_v49(p)
            core_absorbed = chi * delivered
            perim_absorbed = (1.0 - chi) * delivered
            dx = L_v49 / length(result.z_solid)
            Trear = result.rear_reservoir_temperature[:, end]

            path = joinpath(RUNNER_AXIAL_DIR_v49, "cell_diagnostics_$(id_str)_1D_v49.csv")
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

begin # optimization runner
    function run_calibration_v49(; stage=RUNNER_FIT_STAGE_v49, max_time=RUNNER_FIT_SECONDS_v49,
                                nodes=RUNNER_FIT_NODES_v49, dataset=RUNNER_CALIBRATION_DATASET_v49)
        println("[run_1D_v49] Starting Option A calibration on heating runs (stage=$(stage), dataset=$(dataset), nodes=$(nodes), max_time=$(max_time) s)...")
        res = calibrate_v49(
            nodes=nodes,
            maximum_iterations=RUNNER_FIT_ITERATIONS_v49,
            maximum_time_seconds=max_time,
            stage=stage,
            dataset=dataset,
        )
        println("[run_1D_v49] Optimization completed with status: $(res.retcode)")
        println("[run_1D_v49] Final objective: $(res.objective)")
        return (params=res.parameters, objective=res.objective, retcode=res.retcode)
    end
end

function main_1D_v49()
    println("============================================================================")
    println("1D_v49 Calibration and Analysis Pipeline")
    println("============================================================================")

    run_cal = runner_bool_v49("RECEIVER1D_v49_RUN_CALIBRATION", true)
    
    variant = if run_cal
        run_calibration_v49()
    else
        println("[run_1D_v49] Using baseline v49 parameter vector.")
        (params=pnew_v49, objective=calibration_loss_v49(pnew_v49), retcode=:Manual)
    end

    # Write summary
    sum_path = joinpath(RUNNER_OUTPUT_DIR_v49, "optimization_summary_1D_v49.txt")
    open(sum_path, "w") do io
        println(io, "objective=$(variant.objective)")
        println(io, "return_code=$(variant.retcode)")
        println(io, "stage=$(RUNNER_FIT_STAGE_v49)")
        println(io, "calibration_dataset=$(RUNNER_CALIBRATION_DATASET_v49)")
        println(io, "parameters=$(variant.params)")
    end

    # Write parameters
    param_path = joinpath(RUNNER_OUTPUT_DIR_v49, "parameters_fitted_energy_accounting_1D_v49.csv")
    write_parameters_v49(param_path, variant)

    # Compute & write metrics
    println("[run_1D_v49] Computing multi-signal metrics...")
    metrics = compute_metrics_v49(variant.params)
    write_metrics_v49(joinpath(RUNNER_OUTPUT_DIR_v49, "analysis_results_fitted_energy_accounting_1D_v49.csv"), metrics)

    # Build & write steady results
    println("[run_1D_v49] Computing steady-state results and conservation residuals...")
    steady_results = build_steady_results_v49(variant.params)
    write_steady_results_v49(joinpath(RUNNER_OUTPUT_DIR_v49, "steady_results_fitted_energy_accounting_1D_v49.csv"), steady_results)

    # Compute flow slopes
    println("[run_1D_v49] Computing flow slopes...")
    slopes = compute_flow_slopes_v49(steady_results)
    write_flow_slopes_v49(joinpath(RUNNER_OUTPUT_DIR_v49, "flow_slopes_fitted_energy_accounting_1D_v49.csv"), slopes)

    # Compute invariants
    println("[run_1D_v49] Computing invariants...")
    inv_rows = compute_invariants_v49(variant.params, steady_results)
    write_invariants_v49(joinpath(RUNNER_OUTPUT_DIR_v49, "invariant_summary_1D_v49.csv"), inv_rows)

    # Plotting & Diagnostics
    if RUNNER_WRITE_PLOTS_v49
        println("[run_1D_v49] Generating transient plots...")
        plot_transient_cases_v49(variant.params)

        println("[run_1D_v49] Generating axial temperature profile plots...")
        plot_axial_profiles_v49(variant.params)

        println("[run_1D_v49] Generating steady-state comparison parity plot...")
        p_steady = steady_comparison_plot_v49(steady_results; title_suffix="fitted")
        savefig(p_steady, joinpath(RUNNER_PLOT_DIR_v49, "steady_comparison_fitted_energy_accounting_1D_v49.png"))

        println("[run_1D_v49] Writing cell-by-cell axial diagnostics CSVs...")
        write_cell_diagnostics_v49(variant.params)
    end

    println("[run_1D_v49] Pipeline completed successfully. Output in: $RUNNER_OUTPUT_DIR_v49")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_1D_v49()
end
