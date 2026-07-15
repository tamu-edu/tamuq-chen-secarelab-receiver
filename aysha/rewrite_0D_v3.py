import sys

with open(r'd:\kkakosim\github\tamuq-chen-secarelab-receiver\aysha\0D_v3.jl', 'r', encoding='utf-8') as f:
    lines = f.readlines()

new_content = "".join(lines[:335]) # up to line 335 (which is the blank line before "begin")

replacement = """
begin 
    sim_key_heat = ["E67","E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"]
    sim_key_cool = ["C69", "C80", "C81"]

    # ============================================================================
    # SECTION H: POST-PROCESSING — STEADY-STATE COMPARISON (GAS)
    # ============================================================================
    T_steady_gas = DataFrame(sim_id=String[], time=Vector[], T_mod=Vector[], T_exp=Vector[])

    for sm in sim_key_heat
        cond_k = simulation_conditions[sm]
        sel_meas = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tf"), :]
        expdata_tf = sel_meas[:, :temperatures][1]
        time_opt = sel_meas[:, :time][1]
        Ts_avg_init = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tavg_v4"), :temperatures][1][1]

        temp_T = remakeAysha(pnew, cond_k, time_opt, Ts_avg_init; tolr=1e-7)
        push!(T_steady_gas, (sm, time_opt, temp_T[:, 2], expdata_tf))
    end

    T_mod_gas_ss = [T_steady_gas[i, :T_mod][end] for i in 1:length(sim_key_heat)]
    T_exp_gas_ss = [T_steady_gas[i, :T_exp][end] for i in 1:length(sim_key_heat)]

    lims_gas = (300, 850)
    plot_gas_ss = scatter(T_mod_gas_ss, T_exp_gas_ss, label="",
        title="Gas SS Temperature (v3)", xlabel="T_mod (K)", ylabel="T_exp (K)",
        ylims=lims_gas, xlims=lims_gas, aspect_ratio=:equal,
        markersize=6, markerstrokewidth=0)
    plot!(collect(lims_gas), collect(lims_gas), label="1:1", ls=:dash, color=:gray)

    # ============================================================================
    # SECTION I: POST-PROCESSING — STEADY-STATE COMPARISON (SOLID AVG)
    # ============================================================================
    T_steady_solid = DataFrame(sim_id=String[], time=Vector[], T_mod=Vector[], T_exp=Vector[])

    for sm in sim_key_heat
        cond_k = simulation_conditions[sm]
        sel_meas = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tavg_v4"), :]
        expdata_ts = sel_meas[:, :temperatures][1]
        time_opt = sel_meas[:, :time][1]
        Ts_avg_init = expdata_ts[1]

        temp_T = remakeAysha(pnew, cond_k, time_opt, Ts_avg_init; tolr=1e-7)
        push!(T_steady_solid, (sm, time_opt, temp_T[:, 1], expdata_ts))
    end

    T_mod_solid_ss = [T_steady_solid[i, :T_mod][end] for i in 1:length(sim_key_heat)]
    T_exp_solid_ss = [T_steady_solid[i, :T_exp][end] for i in 1:length(sim_key_heat)]

    lims_solid = (300, 1150)
    plot_solid_ss = scatter(T_mod_solid_ss, T_exp_solid_ss, label="",
        title="Solid SS Temperature (v3)", xlabel="T_mod (K)", ylabel="T_exp (K)",
        ylims=lims_solid, xlims=lims_solid, aspect_ratio=:equal,
        markersize=6, markerstrokewidth=0)
    plot!(collect(lims_solid), collect(lims_solid), label="1:1", ls=:dash, color=:gray)

    display(plot(plot_gas_ss, plot_solid_ss, layout=(1, 2), size=(800, 400)))

    # ============================================================================
    # SECTION J: POST-PROCESSING — TRANSIENT PROFILES
    # ============================================================================
    function plot_transient_case(sm, params, is_cooling=false)
        if is_cooling
            cond_k = simulation_conditions_cooling[sm]
            meas_df = measurements_cooling
        else
            cond_k = simulation_conditions[sm]
            meas_df = measurements
        end

        time_exp = meas_df[(meas_df.simulation_id .== sm) .& (meas_df.obs_id .== "_Tavg_v4"), :time][1]
        Ts_exp = meas_df[(meas_df.simulation_id .== sm) .& (meas_df.obs_id .== "_Tavg_v4"), :temperatures][1]
        Tf_exp = meas_df[(meas_df.simulation_id .== sm) .& (meas_df.obs_id .== "_Tf"), :temperatures][1]
        Tinit_val = Ts_exp[1]

        temp_T = remakeAysha(params, cond_k, time_exp, Tinit_val)

        p = plot(title="Transient Profile: $sm (v3)", xlabel="Time (s)", ylabel="Temperature (K)",
                 legend=:outerright, ylims=(280, 1150))
        
        plot!(p, time_exp, temp_T[:, 1], label="T_s,avg (mod)", lw=2, color=colors[1])
        scatter!(p, time_exp, Ts_exp, label="T_s,avg (exp)", ms=3, color=colors[1], markerstrokewidth=0)
        
        plot!(p, time_exp, temp_T[:, 2], label="T_g,out (mod)", lw=2, color=colors[2])
        scatter!(p, time_exp, Tf_exp, label="T_g,out (exp)", ms=3, color=colors[2], markerstrokewidth=0)

        display(p)
    end

    println("\\nPlotting transient profiles for all experiments...")
    for sm in sim_key_heat
        plot_transient_case(sm, pnew, false)
    end
    for sm in sim_key_cool
        plot_transient_case(sm, pnew, true)
    end
end

A, B, C = pnew
T_test=range(300, 1000,100)
plot_flow = plot(xlabel="T (K)", ylabel="h_avg (W/m2.K)", title="h vs T", legend=:topright, color_palette=colors)
flow_rates= [1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0]
for fl in flow_rates
    plot!(T_test, h_avg_f6.(fl, T_test), label="$(fl) Lpm")
end
display(plot_flow)

# ==========================================
# ADDITIONAL METRICS COMPUTATION AND EXPORT
# ==========================================
begin
    println("\\n=== Computing Additional Metrics (0D v3 NTU) ===")
    
    function get_t90(time, temp)
        t_init = temp[1]
        t_ss = temp[end]
        target = t_init + 0.9 * (t_ss - t_init)
        idx = findfirst(t -> (t_ss > t_init ? t >= target : t <= target), temp)
        if idx !== nothing
            return time[idx]
        else
            return time[end]
        end
    end

    function run_analysis_0D(pguess_l, cond_k, time_opt, Tinit_val; tolr=1e-7)
        pguess_temp = Dict(k => pguess_l[i] for (i, (k, v)) in enumerate(p_opt))
        rmp = Dict(pguess_temp..., cond_k...)
        rmp[Tinit] = Tinit_val
        u0_map = Dict(Ts => rmp[Tinit], Ts2 => rmp[Tinit] + 2.0, Ts3 => rmp[Tinit] + 3.0, Ts4 => rmp[Tinit] + 4.0)
        rmp_clean = filter(pair -> !isequal(pair.first, Tinit), rmp)
        modeloptim = remake(prob, u0 = u0_map, p = rmp_clean, tspan = (0.0, time_opt[end]))
        sol = solve(modeloptim, FBDF(), saveat=time_opt, reltol=tolr)
        
        Ts_sol = float.(sol[Ts])
        Tf_sol = float.(sol[Tf])
        Ts2_sol = float.(sol[Ts2])
        Ts3_sol = float.(sol[Ts3])
        Ts4_sol = float.(sol[Ts4])
        
        return Ts_sol, Tf_sol, Ts2_sol, Ts3_sol, Ts4_sol
    end

    csv_path = "D:/kkakosim/sim_comsol/analysis_results_0D_v3.csv"
    global next_run_id = 1.0
    file_exists = isfile(csv_path)
    if file_exists
        try
            existing_df = CSV.read(csv_path, DataFrame)
            if !isempty(existing_df) && "RunID" in names(existing_df)
                global next_run_id = maximum(existing_df.RunID) + 1.0
            end
        catch err
            println("Warning: could not read existing CSV file, starting RunID at 1.0. Error: ", err)
        end
    end

    results_df = DataFrame(
        RunID = Float64[], Case = String[],
        T9_SS_sim = Float64[], T9_SS_exp = Float64[], dT_T09 = Float64[],
        T3_SS_sim = Float64[], T3_SS_exp = Float64[], dT_T03 = Float64[],
        T2_SS_sim = Float64[], T2_SS_exp = Float64[], dT_T02 = Float64[],
        R_leak_sim = Float64[], R_leak_exp = Float64[],
        t90_sim_T09_s = Float64[], t90_exp_T09_s = Float64[], dt90_T09_s = Float64[],
        t90_sim_T03_s = Float64[], t90_exp_T03_s = Float64[], dt90_T03_s = Float64[],
        Gap_ss_sim = Float64[], Gap_ss_exp = Float64[], dGap_ss = Float64[]
    )

    for sm in sim_key_heat
        cond_k = simulation_conditions[sm]
        time_exp = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tavg_v4"), :time][1]
        Tavg_v4_exp = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tavg_v4"), :temperatures][1]
        T3_exp = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tf"), :temperatures][1]
        T2_exp = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_T2"), :temperatures][1]
        Tinit_val = Tavg_v4_exp[1]
        
        Ts_sim, Tf_sim, Ts2_sim, Ts3_sim, Ts4_sim = run_analysis_0D(pnew, cond_k, time_exp, Tinit_val)
        
        Ts_ss = Ts_sim[end]; Tf_ss = Tf_sim[end]; Ts3_ss = Ts3_sim[end]
        T9_exp_ss = Tavg_v4_exp[end]; T3_exp_ss = T3_exp[end]; T2_exp_ss = T2_exp[end]
        
        dT_T09 = Ts_ss - T9_exp_ss; dT_T03 = Tf_ss - T3_exp_ss; dT_T02 = Ts3_ss - T2_exp_ss
        R_leak_sim = (Tf_ss - Tamb) / (Ts_ss - Tamb)
        R_leak_exp = (T3_exp_ss - Tamb) / (T9_exp_ss - Tamb)
        
        t90_sim_T09 = get_t90(time_exp, Ts_sim)
        t90_exp_T09 = get_t90(time_exp, Tavg_v4_exp)
        dt90_T09 = t90_sim_T09 - t90_exp_T09
        t90_sim_T03 = get_t90(time_exp, Tf_sim)
        t90_exp_T03 = get_t90(time_exp, T3_exp)
        dt90_T03 = t90_sim_T03 - t90_exp_T03
        
        gap_sim = Tf_ss - Ts_ss
        gap_exp = T3_exp_ss - T9_exp_ss
        dgap = gap_sim - gap_exp
        
        push!(results_df, (
            next_run_id, sm, Ts_ss, T9_exp_ss, dT_T09, Tf_ss, T3_exp_ss, dT_T03, Ts3_ss, T2_exp_ss, dT_T02,
            R_leak_sim, R_leak_exp, t90_sim_T09, t90_exp_T09, dt90_T09, t90_sim_T03, t90_exp_T03, dt90_T03, gap_sim, gap_exp, dgap
        ))
    end

    for sm in sim_key_cool
        cond_k = simulation_conditions_cooling[sm]
        time_exp = measurements_cooling[(measurements_cooling.simulation_id .== sm) .& (measurements_cooling.obs_id .== "_Tavg_v4"), :time][1]
        Tavg_v4_exp = measurements_cooling[(measurements_cooling.simulation_id .== sm) .& (measurements_cooling.obs_id .== "_Tavg_v4"), :temperatures][1]
        T3_exp = measurements_cooling[(measurements_cooling.simulation_id .== sm) .& (measurements_cooling.obs_id .== "_Tf"), :temperatures][1]
        T2_exp = measurements_cooling[(measurements_cooling.simulation_id .== sm) .& (measurements_cooling.obs_id .== "_T2"), :temperatures][1]
        Tinit_val = Tavg_v4_exp[1]
        
        Ts_sim, Tf_sim, Ts2_sim, Ts3_sim, Ts4_sim = run_analysis_0D(pnew, cond_k, time_exp, Tinit_val)
        
        Ts_ss = Ts_sim[end]; Tf_ss = Tf_sim[end]; Ts3_ss = Ts3_sim[end]
        T9_exp_ss = Tavg_v4_exp[end]; T3_exp_ss = T3_exp[end]; T2_exp_ss = T2_exp[end]
        
        dT_T09 = Ts_ss - T9_exp_ss; dT_T03 = Tf_ss - T3_exp_ss; dT_T02 = Ts3_ss - T2_exp_ss
        R_leak_sim = (Tf_ss - Tamb) / (Ts_ss - Tamb)
        R_leak_exp = (T3_exp_ss - Tamb) / (T9_exp_ss - Tamb)
        
        t90_sim_T09 = get_t90(time_exp, Ts_sim)
        t90_exp_T09 = get_t90(time_exp, Tavg_v4_exp)
        dt90_T09 = t90_sim_T09 - t90_exp_T09
        t90_sim_T03 = get_t90(time_exp, Tf_sim)
        t90_exp_T03 = get_t90(time_exp, T3_exp)
        dt90_T03 = t90_sim_T03 - t90_exp_T03
        
        gap_sim = Tf_ss - Ts_ss
        gap_exp = T3_exp_ss - T9_exp_ss
        dgap = gap_sim - gap_exp
        
        push!(results_df, (
            next_run_id, sm, Ts_ss, T9_exp_ss, dT_T09, Tf_ss, T3_exp_ss, dT_T03, Ts3_ss, T2_exp_ss, dT_T02,
            R_leak_sim, R_leak_exp, t90_sim_T09, t90_exp_T09, dt90_T09, t90_sim_T03, t90_exp_T03, dt90_T03, gap_sim, gap_exp, dgap
        ))
    end
    
    println(results_df)
    
    if file_exists
        CSV.write(csv_path, results_df, append=true)
        println("Appended additional metrics to $csv_path with RunID = $next_run_id")
    else
        CSV.write(csv_path, results_df)
        println("Created new CSV file $csv_path and saved metrics with RunID = $next_run_id")
    end
end
"""

with open(r'd:\kkakosim\github\tamuq-chen-secarelab-receiver\aysha\0D_v3.jl', 'w', encoding='utf-8') as f:
    f.write(new_content + replacement)
