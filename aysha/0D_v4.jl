# ============================================================================
# 0D_v4.jl — Two-State Lumped Parameter Solar Receiver Model with Hysteresis
# ============================================================================
# Model: Two dynamic states (T_s,avg and a*) with algebraic gas outlet (T_g,out)
# Reference: 0D_modelling_v4.md (Section D)
# Changes from initial v4:
#   - Added second state (a*) representing axial thermal distribution profile.
#   - Three-stage parameter estimation: 
#       1. Base Heating (7 params, no hysteresis)
#       2. Cooling Hysteresis (2 params, lock base)
#       3. Heating Hysteresis Refinement (8 params + tight bounds, lock cooling)
# ============================================================================

begin #libraries
    using ModelingToolkit
    using DifferentialEquations, DataFrames
    using CSV, Statistics, Interpolations
    using StatsPlots, ColorSchemes
    using Optimization, OptimizationNLopt
    using ForwardDiff
end

# ============================================================================
# SECTION A: GEOMETRY AND FIXED PARAMETERS
# ============================================================================
begin #FIXED Parameters
    Tamb = (22.448 + 273.15) #K
    w_t = 20.0e-3
    A_frt = w_t * w_t
    L = 137e-3
    w_chnl = 1.65e-3
    A_chnl_frt = w_chnl * w_chnl
    n_chnl = 10 * 10
    A_chnl_frt_all = A_chnl_frt * n_chnl
    Dh = 4 * A_chnl_frt / (4 * w_chnl)
    A_frt_solid = A_frt - A_chnl_frt_all
    A_exchange = 4 * w_chnl * L * n_chnl
    Pi = 4 * w_chnl * n_chnl
    Vs = A_frt_solid * L
    Vf = n_chnl * A_chnl_frt * L
    A_irr = A_frt
    A_rad = A_frt
    σ = 5.670374419e-8
    ρs = 3200.0
end;

# ============================================================================
# SECTION B: MATERIAL PROPERTIES (Temperature-dependent)
# ============================================================================
begin #Material Properties
    ρf_f(T) = 352.716 * T^(-1)
    μf_f(T) = -8.38278e-7 + 8.35717342e-8 * T - 7.69429583e-11 * T^2 + 4.6437266e-14 * T^3 - 1.06585607e-17 * T^4
    kf_f(T) = -0.00227583562 + 1.15480022e-4 * T - 7.90252856e-8 * T^2 + 4.11702505e-11 * T^3 - 7.43864331e-15 * T^4
    cpf_f(T) = 1.93e-10 * T^4 - 8.0e-7 * T^3 + 1.14e-3 * T^2 - 4.49e-1 * T + 1.06e3
    Cps_f(T) = 1110.0 + 0.15 * T - 425.0 * exp(-0.003 * T)
    ks_f(T) = 191.9216 - 0.3261784 * T + 2.739462e-4 * T^2 - 7.70926e-8 * T^3
    m_dot(qlpm) = ρf_f(Tamb) * (qlpm / 60.0 / 1000.0)
end;

# ============================================================================
# SECTION C: TWO-STATE ODE MODEL (Hysteresis)
# ============================================================================
# The model has 10 fitted parameters mapped as:
#   p[1] = γ_C     : effective capacitance correction factor
#   p[2] = K_lin   : linear loss conductance (W/K)
#   p[3] = χε_rad  : χ_rad * ε_rad product for radiation loss
#   p[4] = A       : Nusselt multiplier exponent (10^A)
#   p[5] = B       : Nusselt Reynolds exponent (Re^B)
#   p[6] = C       : Nusselt Prandtl double-exponent (Pr^(10^C))
#   p[7] = τ_a     : characteristic time for axial profile redistribution
#   p[8] = R_g_star: profile decay resistance due to gas heat removal
#   p[9] = η_eff   : effective absorbed solar fraction
#   p[10]= R_q_star: profile generation resistance due to solar flux
# ============================================================================

function compute_Tg_out(Ts_avg, a_star, qlpm, p, Tg_in)
    A = p[4]
    B = p[5]
    C = p[6]

    mdot = m_dot(qlpm)
    T_film = (Ts_avg + Tg_in) / 2.0
    
    cp_g = cpf_f(T_film)
    mu_g = μf_f(T_film)
    k_g = kf_f(T_film)

    Re_f = mdot * Dh / (A_chnl_frt_all * mu_g)
    Pr_f = cp_g * mu_g / k_g
    
    Nu_f = 10^A * (Re_f^B) * (Pr_f^(10^C))
    h_eff = Nu_f * k_g / Dh

    NTU_eff = h_eff * A_exchange / (mdot * cp_g)
    ε_eff = 1.0 - exp(-NTU_eff)

    # Hysteresis adjustment: gas weighted wall temp = Ts_avg + a_star
    Tg_out = Tg_in + ε_eff * (Ts_avg + a_star - Tg_in)

    return Tg_out, ε_eff, cp_g, mdot
end

function receiver_ode!(du, u, params, t)
    Ts_avg = u[1]
    a_star = u[2]

    p = params.p
    Io_val = params.Io
    qlpm_val = params.qlpm
    Tg_in = Tamb

    γ_C = p[1]
    K_lin = p[2]
    χε_rad = p[3]
    τ_a = p[7]
    R_g_star = p[8]
    η_eff = p[9]
    R_q_star = p[10]

    Tg_out, ε_eff, cp_g, mdot = compute_Tg_out(Ts_avg, a_star, qlpm_val, p, Tg_in)
    C_eff = γ_C * ρs * Cps_f(Ts_avg) * Vs

    Q_abs = η_eff * Io_val * A_irr
    Q_g = mdot * cp_g * (Tg_out - Tg_in)
    Q_loss_lin = K_lin * (Ts_avg - Tamb)
    Q_loss_rad = χε_rad * σ * A_rad * (Ts_avg^4 - Tamb^4)

    du[1] = (Q_abs - Q_g - Q_loss_lin - Q_loss_rad) / C_eff
    du[2] = (R_q_star * Q_abs - R_g_star * Q_g - a_star) / τ_a
end

# ============================================================================
# SECTION D: EXPERIMENTAL DATA IMPORT
# ============================================================================
begin
    # Define backward-compatible dictionary keys for the import script
    const Io = :Io
    const qlpm = :qlpm
    const Tinit = :Tinit

    include("import_exp_0D.jl")
    colors = theme_palette(:auto)
end

# ============================================================================
# SECTION E: MODEL SOLVER
# ============================================================================
begin 
    function solve_model(pguess, Io_val, qlpm_val, Tinit_val, tvalues; tolr=1e-7)
        u0 = [Tinit_val, 0.0]  # Start with uniform temperature profile (a* = 0)
        p_ode = (p=pguess, Io=Io_val, qlpm=qlpm_val)
        prob = ODEProblem(receiver_ode!, u0, (tvalues[1], tvalues[end]), p_ode)
        sol = solve(prob, Tsit5(), saveat=tvalues, reltol=tolr, abstol=tolr)
        
        Ts_sol = [u[1] for u in sol.u]
        a_sol = [u[2] for u in sol.u]
        Tf_sol = zeros(length(tvalues))
        for i in 1:length(tvalues)
            Tg_out, _, _, _ = compute_Tg_out(Ts_sol[i], a_sol[i], qlpm_val, pguess, Tamb)
            Tf_sol[i] = Tg_out
        end
        return Ts_sol, Tf_sol, a_sol
    end

    function remakeAysha(pguess, cond_k, time_opt, Tinit_val; tolr=1e-7)
        Io_val = cond_k[Io]
        qlpm_val = cond_k[qlpm]
        Ts_sol, Tf_sol, a_sol = solve_model(pguess, Io_val, qlpm_val, Tinit_val, time_opt; tolr=tolr)
        return hcat(Ts_sol, Tf_sol, a_sol)
    end
end

# ============================================================================
# SECTION F: THREE-STAGE OPTIMIZATION
# ============================================================================
begin
    sim_key_heat = ["E67","E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"]
    sim_key_cool = ["C69", "C80", "C81"]

    # Trigger to use both temperatures (false) or only gas temperature (true) for fitting
    fit_gas_only = true

    # ==========================================================
    # STAGE 1: BASE HEATING OPTIMIZATION (NO HYSTERESIS)
    # ==========================================================
    function loss_stage1(p_base, keys)
        lossr = 0.0
        # p_base = [γ_C, K_lin, χε_rad, a_h, n_exp, χ_L, η_eff]
        # full =   [γ_C, K_lin, χε_rad, a_h, n_exp, χ_L, τ_a=200, R_g=0, η_eff, R_q=0]
        p_full = [p_base[1], p_base[2], p_base[3], p_base[4], p_base[5], p_base[6], 200.0, 0.0, p_base[7], 0.0]
        
        for sm in keys
            cond_k = simulation_conditions[sm]
            sel_meas = measurements[(measurements.simulation_id .== sm), :]
            time_opt = sel_meas[(sel_meas.obs_id .== "_Tavg_v4"), :time][1]
            Ts_exp = sel_meas[(sel_meas.obs_id .== "_Tavg_v4"), :temperatures][1]
            Tf_exp = sel_meas[(sel_meas.obs_id .== "_Tf"), :temperatures][1]
            
            temp_T = remakeAysha(p_full, cond_k, time_opt, Ts_exp[1])
            N_samples = length(time_opt)
            
            error_s = sum((temp_T[:, 1] .- Ts_exp) .^ 2) / (maximum(Ts_exp) - minimum(Ts_exp))^2
            error_f = sum((temp_T[:, 2] .- Tf_exp) .^ 2) / (maximum(Tf_exp) - minimum(Tf_exp))^2
            lossr += fit_gas_only ? (error_f / N_samples) : ((error_s + error_f) / N_samples)
        end
        return lossr / length(keys)
    end

    println("\n--- STAGE 1: Heating Base Optimization (No Hysteresis) ---")
    p_s1_init = [1.0, 0.5, 1.0, -1.0, 1.0, 0.5, 0.8]
    lb_s1 = [0.1, 0.0, 0.0, -3.0, 0.5, -1.0, 0.1]
    ub_s1 = [5.0, 5.0, 5.0, 1.0, 2.0, 2.0, 1.0]

    optf_s1 = OptimizationFunction(loss_stage1, Optimization.AutoForwardDiff())
    optprob_s1 = Optimization.OptimizationProblem(optf_s1, p_s1_init, sim_key_heat, lb=lb_s1, ub=ub_s1)
    optsol_s1 = solve(optprob_s1, NLopt.GN_MLSL_LDS(), local_method=NLopt.LN_NELDERMEAD(), maxtime=60)
    p_base_opt = optsol_s1.u
    println("Stage 1 Loss: ", optsol_s1.objective)

    # ==========================================================
    # STAGE 2: COOLING HYSTERESIS FIT
    # ==========================================================
    function loss_stage2(p_cool, keys)
        lossr = 0.0
        # p_cool = [τ_a, R_g_star]
        p_full = [p_base_opt[1], p_base_opt[2], p_base_opt[3], p_base_opt[4], p_base_opt[5], p_base_opt[6], p_cool[1], p_cool[2], 0.0, 0.0]
        for sm in keys
            cond_k = simulation_conditions_cooling[sm]
            sel_meas = measurements_cooling[(measurements_cooling.simulation_id .== sm), :]
            time_opt = sel_meas[(sel_meas.obs_id .== "_Tavg_v4"), :time][1]
            Ts_exp = sel_meas[(sel_meas.obs_id .== "_Tavg_v4"), :temperatures][1]
            Tf_exp = sel_meas[(sel_meas.obs_id .== "_Tf"), :temperatures][1]
            
            temp_T = remakeAysha(p_full, cond_k, time_opt, Ts_exp[1])
            N_samples = length(time_opt)
            
            error_s = sum((temp_T[:, 1] .- Ts_exp) .^ 2) / (maximum(Ts_exp) - minimum(Ts_exp))^2
            error_f = sum((temp_T[:, 2] .- Tf_exp) .^ 2) / (maximum(Tf_exp) - minimum(Tf_exp))^2
            lossr += fit_gas_only ? (error_f / N_samples) : ((error_s + error_f) / N_samples)
        end
        return lossr / length(keys)
    end

    println("\n--- STAGE 2: Cooling Hysteresis Fit ---")
    p_s2_init = [200.0, 0.5]
    lb_s2 = [10.0, 0.0]
    ub_s2 = [1000.0, 5.0]

    optf_s2 = OptimizationFunction(loss_stage2, Optimization.AutoForwardDiff())
    optprob_s2 = Optimization.OptimizationProblem(optf_s2, p_s2_init, sim_key_cool, lb=lb_s2, ub=ub_s2)
    optsol_s2 = solve(optprob_s2, NLopt.GN_MLSL_LDS(), local_method=NLopt.LN_NELDERMEAD(), maxtime=60)
    p_cool_opt = optsol_s2.u
    println("Stage 2 Loss: ", optsol_s2.objective)

    # ==========================================================
    # STAGE 3: HEATING HYSTERESIS REFINEMENT
    # ==========================================================
    function loss_stage3(p_heat, keys)
        lossr = 0.0
        # p_heat = [γ_C, K_lin, χε_rad, A, B, C, η_eff, R_q_star]
        p_full = [p_heat[1], p_heat[2], p_heat[3], p_heat[4], p_heat[5], p_heat[6], p_cool_opt[1], p_cool_opt[2], p_heat[7], p_heat[8]]
        for sm in keys
            cond_k = simulation_conditions[sm]
            sel_meas = measurements[(measurements.simulation_id .== sm), :]
            time_opt = sel_meas[(sel_meas.obs_id .== "_Tavg_v4"), :time][1]
            Ts_exp = sel_meas[(sel_meas.obs_id .== "_Tavg_v4"), :temperatures][1]
            Tf_exp = sel_meas[(sel_meas.obs_id .== "_Tf"), :temperatures][1]
            
            temp_T = remakeAysha(p_full, cond_k, time_opt, Ts_exp[1])
            N_samples = length(time_opt)
            
            error_s = sum((temp_T[:, 1] .- Ts_exp) .^ 2) / (maximum(Ts_exp) - minimum(Ts_exp))^2
            error_f = sum((temp_T[:, 2] .- Tf_exp) .^ 2) / (maximum(Tf_exp) - minimum(Tf_exp))^2
            lossr += fit_gas_only ? (error_f / N_samples) : ((error_s + error_f) / N_samples)
        end
        return lossr / length(keys)
    end

    println("\n--- STAGE 3: Heating Hysteresis Refinement ---")
    p_s3_init = [p_base_opt[1], p_base_opt[2], p_base_opt[3], p_base_opt[4], p_base_opt[5], p_base_opt[6], p_base_opt[7], 0.5]
    # Tighten bounds slightly around base to prevent drift: +/- 20%
    lb_s3 = [max(0.1, p_base_opt[1]*0.8), max(0.0, p_base_opt[2]*0.8), max(0.0, p_base_opt[3]*0.8), max(-3.0, p_base_opt[4]-0.5), max(0.5, p_base_opt[5]*0.8), max(-1.0, p_base_opt[6]-0.5), max(0.1, p_base_opt[7]*0.8), 0.0]
    ub_s3 = [min(5.0, p_base_opt[1]*1.2), min(5.0, p_base_opt[2]*1.2), min(5.0, p_base_opt[3]*1.2), min(1.0, p_base_opt[4]+0.5), min(2.0, p_base_opt[5]*1.2), min(2.0, p_base_opt[6]+0.5), min(1.0, p_base_opt[7]*1.2), 5.0]

    optf_s3 = OptimizationFunction(loss_stage3, Optimization.AutoForwardDiff())
    optprob_s3 = Optimization.OptimizationProblem(optf_s3, p_s3_init, sim_key_heat, lb=lb_s3, ub=ub_s3)
    optsol_s3 = solve(optprob_s3, NLopt.GN_MLSL_LDS(), local_method=NLopt.LN_NELDERMEAD(), maxtime=60)
    p_heat_opt = optsol_s3.u
    println("Stage 3 Loss: ", optsol_s3.objective)

    # FINAL 10-PARAM ARRAY
    global pnew = [p_heat_opt[1], p_heat_opt[2], p_heat_opt[3], p_heat_opt[4], p_heat_opt[5], p_heat_opt[6], p_cool_opt[1], p_cool_opt[2], p_heat_opt[7], p_heat_opt[8]]
    println("\nFinal Tuned Parameters: ", pnew)
end

# ============================================================================
# SECTION G: POST-PROCESSING — STEADY-STATE COMPARISON
# ============================================================================
begin
    T_steady_gas = DataFrame(sim_id=String[], time=Vector[], T_mod=Vector[], T_exp=Vector[])
    T_steady_solid = DataFrame(sim_id=String[], time=Vector[], T_mod=Vector[], T_exp=Vector[])

    for sm in sim_key_heat
        cond_k = simulation_conditions[sm]
        Ts_exp = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tavg_v4"), :temperatures][1]
        Tf_exp = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tf"), :temperatures][1]
        time_opt = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tf"), :time][1]
        
        temp_T = remakeAysha(pnew, cond_k, time_opt, Ts_exp[1])
        push!(T_steady_gas, (sm, time_opt, temp_T[:, 2], Tf_exp))
        push!(T_steady_solid, (sm, time_opt, temp_T[:, 1], Ts_exp))
    end

    T_mod_gas_ss = [T_steady_gas[i, :T_mod][end] for i in 1:length(sim_key_heat)]
    T_exp_gas_ss = [T_steady_gas[i, :T_exp][end] for i in 1:length(sim_key_heat)]
    T_mod_solid_ss = [T_steady_solid[i, :T_mod][end] for i in 1:length(sim_key_heat)]
    T_exp_solid_ss = [T_steady_solid[i, :T_exp][end] for i in 1:length(sim_key_heat)]

    lims_gas = (300, 850)
    plot_gas_ss = scatter(T_mod_gas_ss, T_exp_gas_ss, label="", title="Gas SS Temp (v4 Hys)",
        xlabel="T_mod (K)", ylabel="T_exp (K)", ylims=lims_gas, xlims=lims_gas, aspect_ratio=:equal)
    plot!(collect(lims_gas), collect(lims_gas), label="1:1", ls=:dash, color=:gray)

    lims_solid = (300, 1150)
    plot_solid_ss = scatter(T_mod_solid_ss, T_exp_solid_ss, label="", title="Solid SS Temp (v4 Hys)",
        xlabel="T_mod (K)", ylabel="T_exp (K)", ylims=lims_solid, xlims=lims_solid, aspect_ratio=:equal)
    plot!(collect(lims_solid), collect(lims_solid), label="1:1", ls=:dash, color=:gray)

    display(plot(plot_gas_ss, plot_solid_ss, layout=(1, 2), size=(800, 400)))
end

# ============================================================================
# SECTION J: TRANSIENT COMPARISON PLOTS
# ============================================================================
begin
    function plt_case_v4(sm, params; is_cooling=false)
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
        
        temp_T = remakeAysha(params, cond_k, time_exp, Ts_exp[1])

        p = plot(title="Transient Profile: $sm (v4 Hysteresis)", xlabel="Time (s)", ylabel="Temperature (K)", legend=:outerright)
        plot!(p, time_exp, temp_T[:, 1], label="T_s,avg (mod)", lw=2, color=colors[1])
        scatter!(p, time_exp, Ts_exp, label="T_s,avg (exp)", ms=3, color=colors[1], markerstrokewidth=0)
        
        plot!(p, time_exp, temp_T[:, 2], label="T_g,out (mod)", lw=2, color=colors[2])
        scatter!(p, time_exp, Tf_exp, label="T_g,out (exp)", ms=3, color=colors[2], markerstrokewidth=0)

        display(p)
    end

    println("\nPlotting transient profiles...")
    for sm in sim_key_heat
        plt_case_v4(sm, pnew, is_cooling=false)
    end
    for sm in sim_key_cool
        plt_case_v4(sm, pnew, is_cooling=true)
    end
end

# ============================================================================
# SECTION K: EFFECTIVENESS SURFACE PLOTS
# ============================================================================
begin
    T_range = range(300, 1000, length=50)
    flow_rates = [1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0]

    plot_h = plot(xlabel="T_film (K)", ylabel="ε_eff", title="Gas effectiveness ε_eff vs T", legend=:topright)
    for fl in flow_rates
        mdot_fl = m_dot(fl)
        ε_vals = [begin
            T_film = (T + Tamb) / 2.0
            cp = cpf_f(T_film)
            NTU = (pnew[4] * Pi * L / cp) * mdot_fl^(pnew[5] - 1.0)
            pnew[6] * (1.0 - exp(-NTU))
        end for T in T_range]
        plot!(T_range, ε_vals, label="$(fl) L/min")
    end
    display(plot_h)

    fl_range = range(1, 30, length=30)
    plot_ntu = plot(xlabel="Flow rate (L/min)", ylabel="NTU_eff", title="NTU vs Flow Rate", legend=:topright)
    for T_test in [400.0, 600.0, 800.0, 1000.0]
        ntu_vals = [begin
            mdot_fl = m_dot(fl)
            T_film = (T_test + Tamb) / 2.0
            cp = cpf_f(T_film)
            (pnew[4] * Pi * L / cp) * mdot_fl^(pnew[5] - 1.0)
        end for fl in fl_range]
        plot!(fl_range, ntu_vals, label="T_s=$(T_test) K")
    end
    display(plot_ntu)
end

# ============================================================================
# SECTION L: METRICS COMPUTATION AND EXPORT
# ============================================================================
begin
    println("\n=== Computing Metrics (0D v4 Hysteresis) ===")

    function get_t90(time, temp)
        t_init = temp[1]
        t_ss = temp[end]
        target = t_init + 0.9 * (t_ss - t_init)
        if t_ss > t_init
            idx = findfirst(t -> t >= target, temp)
        else
            idx = findfirst(t -> t <= target, temp)
        end
        return idx !== nothing ? time[idx] : time[end]
    end

    csv_path = "D:/kkakosim/sim_comsol/analysis_results_0D_v4.csv"
    global next_run_id = 1.0
    file_exists = isfile(csv_path)

    if file_exists
        try
            existing_df = CSV.read(csv_path, DataFrame)
            if !isempty(existing_df) && "RunID" in names(existing_df)
                global next_run_id = maximum(existing_df.RunID) + 1.0
            end
        catch err
            println("Warning: could not read existing CSV. Error: ", err)
        end
    end

    results_df = DataFrame(
        RunID = Float64[], Case = String[], Type = String[],
        Ts_avg_SS_sim = Float64[], Ts_avg_SS_exp = Float64[], dT_Ts_avg = Float64[],
        T3_SS_sim = Float64[], T3_SS_exp = Float64[], dT_T03 = Float64[],
        R_leak_sim = Float64[], R_leak_exp = Float64[],
        t90_sim_Ts_s = Float64[], t90_exp_Ts_s = Float64[], dt90_Ts_s = Float64[],
        t90_sim_T3_s = Float64[], t90_exp_T3_s = Float64[], dt90_T3_s = Float64[],
        Gap_ss_sim = Float64[], Gap_ss_exp = Float64[], dGap_ss = Float64[],
        ε_eff = Float64[], NTU_eff = Float64[], h_eff_sim = Float64[]
    )

    function process_metrics(sm, is_cooling)
        cond_k = is_cooling ? simulation_conditions_cooling[sm] : simulation_conditions[sm]
        meas_df = is_cooling ? measurements_cooling : measurements
        
        Io_val = cond_k[Io]
        qlpm_val = cond_k[qlpm]
        
        time_exp = meas_df[(meas_df.simulation_id .== sm) .& (meas_df.obs_id .== "_Tavg_v4"), :time][1]
        Ts_exp = meas_df[(meas_df.simulation_id .== sm) .& (meas_df.obs_id .== "_Tavg_v4"), :temperatures][1]
        Tf_exp = meas_df[(meas_df.simulation_id .== sm) .& (meas_df.obs_id .== "_Tf"), :temperatures][1]
        
        Ts_mod, Tf_mod, a_mod = solve_model(pnew, Io_val, qlpm_val, Ts_exp[1], time_exp)

        Ts_ss_sim = Ts_mod[end]
        Tf_ss_sim = Tf_mod[end]
        Ts_ss_exp = Ts_exp[end]
        Tf_ss_exp = Tf_exp[end]

        dT_Ts = Ts_ss_sim - Ts_ss_exp
        dT_T3 = Tf_ss_sim - Tf_ss_exp

        if is_cooling && abs(Ts_ss_sim - Tamb) < 1.0
            R_leak_sim = NaN
            R_leak_exp = NaN
        else
            R_leak_sim = (Tf_ss_sim - Tamb) / (Ts_ss_sim - Tamb)
            R_leak_exp = (Tf_ss_exp - Tamb) / (Ts_ss_exp - Tamb)
        end

        t90_sim_Ts = get_t90(time_exp, Ts_mod)
        t90_exp_Ts = get_t90(time_exp, Ts_exp)
        dt90_Ts = t90_sim_Ts - t90_exp_Ts

        t90_sim_T3 = get_t90(time_exp, Tf_mod)
        t90_exp_T3 = get_t90(time_exp, Tf_exp)
        dt90_T3 = t90_sim_T3 - t90_exp_T3

        gap_sim = Tf_ss_sim - Ts_ss_sim
        gap_exp = Tf_ss_exp - Ts_ss_exp
        dgap = gap_sim - gap_exp

        _, ε_ss, cp_ss, _ = compute_Tg_out(Ts_ss_sim, a_mod[end], qlpm_val, pnew, Tamb)
        mdot_val = m_dot(qlpm_val)
        
        T_film_ss = (Ts_ss_sim + Tamb) / 2.0
        mu_g_ss = μf_f(T_film_ss)
        k_g_ss = kf_f(T_film_ss)
        
        Re_f_ss = mdot_val * Dh / (A_chnl_frt_all * mu_g_ss)
        Pr_f_ss = cp_ss * mu_g_ss / k_g_ss
        
        Nu_f_ss = 10^pnew[4] * (Re_f_ss^pnew[5]) * (Pr_f_ss^(10^pnew[6]))
        h_eff_val = Nu_f_ss * k_g_ss / Dh
        NTU_ss = h_eff_val * A_exchange / (mdot_val * cp_ss)

        push!(results_df, (
            next_run_id, sm, is_cooling ? "cooling" : "heating",
            Ts_ss_sim, Ts_ss_exp, dT_Ts,
            Tf_ss_sim, Tf_ss_exp, dT_T3,
            R_leak_sim, R_leak_exp,
            t90_sim_Ts, t90_exp_Ts, dt90_Ts,
            t90_sim_T3, t90_exp_T3, dt90_T3,
            gap_sim, gap_exp, dgap,
            ε_ss, NTU_ss, h_eff_val
        ))
    end

    for sm in sim_key_heat process_metrics(sm, false) end
    for sm in sim_key_cool process_metrics(sm, true) end

    println(results_df)

    if file_exists
        CSV.write(csv_path, results_df, append=true)
        println("Appended metrics to $csv_path with RunID = $next_run_id")
    else
        CSV.write(csv_path, results_df)
        println("Created $csv_path with RunID = $next_run_id")
    end
end
