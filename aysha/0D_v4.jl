# ============================================================================
# 0D_v4.jl — Single-Node Lumped Parameter Solar Receiver Model
# ============================================================================
# Model: One dynamic state (T_s,avg) with algebraic gas outlet (T_g,out)
# Reference: 0D_modelling_v4.md (Sections B.1–B.15)
# Changes from v3:
#   - Single effective thermal mass instead of 5-node network
#   - Length-corrected NTU effectiveness (χ_L correction)
#   - Nonlinear radiation loss retained
#   - Film-temperature gas property evaluation
#   - Channel hydraulic diameter Dh used consistently
#   - Joint heating + cooling data fitting
#   - Trapezoidal-weighted T_s,avg from TC positions
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

    Tamb = (22.448 + 273.15) #K (same for all exp)

    # Receiver dimensions (COMSOL D_tot = 20mm, L_channel = 137mm)
    w_t = 20.0e-3  # Width of the receiver in meters
    A_frt = w_t * w_t  # Total frontal area of the receiver (m^2)
    L = 137e-3  # Length (m)

    # Channel dimensions (COMSOL R_channel = 1.65mm, t_channel = 0.35mm)
    w_chnl = 1.65e-3  # Width of a single channel (m)
    A_chnl_frt = w_chnl * w_chnl  # Frontal area of a single channel (m^2)

    # Channel array
    n_chnl = 10 * 10  # Number of channels (10x10 array)
    A_chnl_frt_all = A_chnl_frt * n_chnl  # Total frontal area of all channels (m^2)

    # Hydraulic diameters
    Dh = 4 * A_chnl_frt / (4 * w_chnl) # Channel hydraulic diameter (m) = w_chnl for square
    # Note: v3 used Lc = w_t = 20mm in the Nu correlation. v4 uses Dh = w_chnl correctly.

    # Solid geometrical calculations
    A_frt_solid = A_frt - A_chnl_frt_all  # Frontal area of just solid (m^2)
    A_exchange = 4 * w_chnl * L * n_chnl  # Total internal wetted area (m^2)

    # Perimeter for NTU calculation
    Pi = 4 * w_chnl * n_chnl  # Total internal wetted perimeter (m)

    # Volume calculations
    Vs = A_frt_solid * L  # Solid volume (m^3)
    Vf = n_chnl * A_chnl_frt * L  # Fluid volume (m^3)

    # Irradiated area (front face of receiver)
    A_irr = A_frt  # (m^2)

    # Radiation loss area (front face — the only surface directly exposed)
    A_rad = A_frt  # (m^2)

    # Constants
    σ = 5.670374419e-8  # W/m^2·K^4 Stefan-Boltzmann constant

    # SiC Material properties (fixed reference values)
    ρs = 3200.0  # kg/m^3

end;

# ============================================================================
# SECTION B: MATERIAL PROPERTIES (Temperature-dependent)
# ============================================================================
begin #Material Properties

    ### Fluid (Air) Properties ###
    # Density (ideal gas approximation)
    ρf_f(T) = 352.716 * T^(-1) # kg/m^3 (dry air, COMSOL correlation)

    # Dynamic viscosity (Pa·s) — COMSOL polynomial
    μf_f(T) = -8.38278e-7 + 8.35717342e-8 * T - 7.69429583e-11 * T^2 +
              4.6437266e-14 * T^3 - 1.06585607e-17 * T^4

    # Thermal conductivity (W/m·K) — COMSOL polynomial
    kf_f(T) = -0.00227583562 + 1.15480022e-4 * T - 7.90252856e-8 * T^2 +
              4.11702505e-11 * T^3 - 7.43864331e-15 * T^4

    # Specific heat capacity (J/kg·K) — polynomial fit
    cpf_f(T) = 1.93e-10 * T^4 - 8.0e-7 * T^3 + 1.14e-3 * T^2 -
               4.49e-1 * T + 1.06e3

    ### Solid (SiC) Properties ###
    # Specific heat capacity (J/kg·K) — Munro 1997
    Cps_f(T) = 1110.0 + 0.15 * T - 425.0 * exp(-0.003 * T)

    # Thermal conductivity (W/m·K) — COMSOL SiC Alpha Polycrystalline
    ks_f(T) = 191.9216 - 0.3261784 * T + 2.739462e-4 * T^2 - 7.70926e-8 * T^3

    ### Mass flow rate from volumetric flow ###
    # qlpm: volumetric flow rate at flowmeter conditions (L/min)
    # Mass flow rate computed at ambient density
    m_dot(qlpm) = ρf_f(Tamb) * (qlpm / 60.0 / 1000.0) # kg/s

end;

# ============================================================================
# SECTION C: SINGLE-STATE ODE MODEL
# ============================================================================
# The model has 7 fitted parameters:
#   p[1] = γ_C     : effective capacitance correction factor
#   p[2] = η_eff   : effective absorbed solar fraction
#   p[3] = K_lin   : linear loss conductance (W/K)
#   p[4] = χε_rad  : χ_rad * ε_rad product for radiation loss
#   p[5] = a_h     : gas HTC coefficient in h_eff = a_h * ṁ^n
#   p[6] = n_exp   : gas HTC flow-rate exponent
#   p[7] = χ_L     : length shape correction factor
#
# Inputs per experiment:
#   Io    : incident solar flux (W/m^2)
#   qlpm  : volumetric flow rate (L/min)
#   Tinit : initial solid temperature (K)
# ============================================================================

"""
    compute_Tg_out(Ts_avg, qlpm, p, Tg_in)

Compute the gas outlet temperature algebraically via the length-corrected
NTU-effectiveness relation:

    T_g,out = T_g,in + ε_eff * (T_s,avg - T_g,in)

where:
    ε_eff = χ_L * [1 - exp(-NTU_eff)]
    NTU_eff = (a_h * P_i * L / cp_g) * ṁ^(n-1)

Gas properties are evaluated at the film temperature T_film = (Ts_avg + Tg_in)/2.
"""
function compute_Tg_out(Ts_avg, qlpm, p, Tg_in)
    a_h = p[5]
    n_exp = p[6]
    χ_L = p[7]

    mdot = m_dot(qlpm)

    # Film temperature for gas property evaluation
    T_film = (Ts_avg + Tg_in) / 2.0
    cp_g = cpf_f(T_film)

    # NTU_eff = (a_h * Pi * L / cp_g) * ṁ^(n-1)
    NTU_eff = (a_h * Pi * L / cp_g) * mdot^(n_exp - 1.0)

    # Effective gas-side effectiveness
    ε_eff = χ_L * (1.0 - exp(-NTU_eff))

    # Gas outlet temperature
    Tg_out = Tg_in + ε_eff * (Ts_avg - Tg_in)

    return Tg_out
end

"""
    receiver_ode!(du, u, params, t)

Single-state ODE for the receiver energy balance:

    C_eff * dT_s,avg/dt = Q_abs - Q_g - Q_loss_lin - Q_loss_rad

where:
    Q_abs = η_eff * Io * A_irr
    Q_g   = ṁ * cp_g * ε_eff * (T_s,avg - T_g,in)
    Q_loss_lin = K_lin * (T_s,avg - T_∞)
    Q_loss_rad = χε_rad * σ * A_rad * (T_s,avg^4 - T_sur^4)
"""
function receiver_ode!(du, u, params, t)
    Ts_avg = u[1]

    # Unpack parameters
    p = params.p       # fitted parameters [γ_C, η_eff, K_lin, χε_rad, a_h, n_exp, χ_L]
    Io_val = params.Io     # incident solar flux (W/m^2), 0 for cooling
    qlpm_val = params.qlpm # volumetric flow rate (L/min)
    Tg_in = Tamb       # gas inlet temperature = ambient

    γ_C = p[1]
    η_eff = p[2]
    K_lin = p[3]
    χε_rad = p[4]
    a_h = p[5]
    n_exp = p[6]
    χ_L = p[7]

    mdot = m_dot(qlpm_val)

    # Effective thermal capacitance
    C_eff = γ_C * ρs * Cps_f(Ts_avg) * Vs

    # Film temperature for gas property evaluation
    T_film = (Ts_avg + Tg_in) / 2.0
    cp_g = cpf_f(T_film)

    # NTU and effectiveness
    NTU_eff = (a_h * Pi * L / cp_g) * mdot^(n_exp - 1.0)
    ε_eff = χ_L * (1.0 - exp(-NTU_eff))

    # Heat fluxes
    Q_abs = η_eff * Io_val * A_irr                      # Absorbed solar power (W)
    Q_g = mdot * cp_g * ε_eff * (Ts_avg - Tg_in)       # Gas heat removal (W)
    Q_loss_lin = K_lin * (Ts_avg - Tamb)                 # Linear thermal losses (W)
    Q_loss_rad = χε_rad * σ * A_rad * (Ts_avg^4 - Tamb^4) # Radiation losses (W)

    # Energy balance
    du[1] = (Q_abs - Q_g - Q_loss_lin - Q_loss_rad) / C_eff
end

# ============================================================================
# SECTION D: EXPERIMENTAL DATA IMPORT
# ============================================================================
# Define the symbolic parameters that import_exp_0D.jl uses as Dict keys
# These must match the @parameters defined in v3's context
@parameters Io qlpm Tinit

include("import_exp_0D.jl") #import the experimental data
colors = ColorSchemes.tab10[1:4]

# ============================================================================
# SECTION E: MODEL SOLVER AND PARAMETER ESTIMATION
# ============================================================================
begin # define solver and loss functions

    """
        solve_model(pguess, Io, qlpm, Tinit, tvalues; tolr=1e-7)

    Solve the single-state ODE for given parameters and conditions.
    Returns (Ts_avg_vec, Tg_out_vec) at the requested time points.
    """
    function solve_model(pguess, Io_val, qlpm_val, Tinit_val, tvalues; tolr=1e-7)
        params = (p=pguess, Io=Io_val, qlpm=qlpm_val)

        u0 = [Tinit_val]
        tspan = (0.0, tvalues[end])

        prob = ODEProblem(receiver_ode!, u0, tspan, params)
        sol = solve(prob, FBDF(), saveat=tvalues, reltol=tolr, abstol=1e-10)

        if sol.retcode != :Success && sol.retcode != ReturnCode.Success
            return fill(NaN, length(tvalues)), fill(NaN, length(tvalues))
        end

        Ts_avg_sol = [sol.u[i][1] for i in 1:length(sol.t)]
        Tg_out_sol = [compute_Tg_out(sol.u[i][1], qlpm_val, pguess, Tamb) for i in 1:length(sol.t)]

        return Ts_avg_sol, Tg_out_sol
    end

    """
        lossAll_heating(pguess, sim_key)

    Compute normalized MSE across all heating experiments.
    Matches both T_s,avg (weighted) and T_g,out (T3).
    """
    function lossAll_heating(pguess, sim_key)
        lossr = 0.0
        n_valid = 0
        for sm in sim_key
            cond_k = simulation_conditions[sm]
            Io_val = cond_k[Io]
            qlpm_val = cond_k[qlpm]
            Tinit_val = cond_k[Tinit]

            # Experimental data: T_s,avg (weighted) and T_g,out
            time_exp = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tavg_v4"), :time][1]
            Ts_exp = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tavg_v4"), :temperatures][1]
            Tf_exp = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tf"), :temperatures][1]

            # Use initial T_s,avg from weighted data as the solid state IC
            Tinit_solid = Ts_exp[1]

            # Run model
            Ts_mod, Tf_mod = solve_model(pguess, Io_val, qlpm_val, Tinit_solid, time_exp)

            if any(isnan, Ts_mod) || length(Ts_exp) != length(Ts_mod)
                return Inf
            end

            N_samples = length(Ts_exp)

            # Normalization by total temperature rise
            delta_Ts = Ts_exp[end] - Ts_exp[1]
            delta_Tf = Tf_exp[end] - Tf_exp[1]
            delta_Ts = abs(delta_Ts) > 1e-3 ? delta_Ts : 1.0
            delta_Tf = abs(delta_Tf) > 1e-3 ? delta_Tf : 1.0

            # Normalized SSE
            error_s = sum(((Ts_mod .- Ts_exp) .^ 2)) / (delta_Ts^2)
            error_f = sum(((Tf_mod .- Tf_exp) .^ 2)) / (delta_Tf^2)

            lossr += (error_s + error_f) / N_samples
            n_valid += 1
        end
        return n_valid > 0 ? lossr / n_valid : Inf
    end

    """
        lossAll_cooling(pguess, sim_key_cool)

    Compute normalized MSE across cooling experiments.
    During cooling: Io = 0, flow continues, system decays.
    Matches both T_s,avg (weighted) and T_g,out (T3).
    """
    function lossAll_cooling(pguess, sim_key_cool)
        lossr = 0.0
        n_valid = 0
        for sm in sim_key_cool
            cond_k = simulation_conditions_cooling[sm]
            qlpm_val = cond_k[qlpm]
            Tinit_val = cond_k[Tinit]

            # Experimental data
            time_exp = measurements_cooling[(measurements_cooling.simulation_id .== sm) .& (measurements_cooling.obs_id .== "_Tavg_v4"), :time][1]
            Ts_exp = measurements_cooling[(measurements_cooling.simulation_id .== sm) .& (measurements_cooling.obs_id .== "_Tavg_v4"), :temperatures][1]
            Tf_exp = measurements_cooling[(measurements_cooling.simulation_id .== sm) .& (measurements_cooling.obs_id .== "_Tf"), :temperatures][1]

            # Use initial T_s,avg from weighted data as the solid state IC
            Tinit_solid = Ts_exp[1]

            # Run model with Io = 0 (no solar input during cooling)
            Ts_mod, Tf_mod = solve_model(pguess, 0.0, qlpm_val, Tinit_solid, time_exp)

            if any(isnan, Ts_mod) || length(Ts_exp) != length(Ts_mod)
                return Inf
            end

            N_samples = length(Ts_exp)

            # Normalization by total temperature drop
            delta_Ts = abs(Ts_exp[1] - Ts_exp[end])
            delta_Tf = abs(Tf_exp[1] - Tf_exp[end])
            delta_Ts = delta_Ts > 1e-3 ? delta_Ts : 1.0
            delta_Tf = delta_Tf > 1e-3 ? delta_Tf : 1.0

            error_s = sum(((Ts_mod .- Ts_exp) .^ 2)) / (delta_Ts^2)
            error_f = sum(((Tf_mod .- Tf_exp) .^ 2)) / (delta_Tf^2)

            lossr += (error_s + error_f) / N_samples
            n_valid += 1
        end
        return n_valid > 0 ? lossr / n_valid : Inf
    end

    """
        lossAll_combined(pguess, data)

    Joint loss function: weighted combination of heating and cooling losses.
    data = (sim_key_heat, sim_key_cool, w_cool)
    where w_cool is the relative weight given to cooling experiments.
    """
    function lossAll_combined(pguess, data)
        sim_key_heat, sim_key_cool, w_cool = data
        loss_heat = lossAll_heating(pguess, sim_key_heat)
        loss_cool = lossAll_cooling(pguess, sim_key_cool)

        # Weighted average: heating dominates (15 exps) but cooling is valuable
        n_heat = length(sim_key_heat)
        n_cool = length(sim_key_cool)

        if isinf(loss_heat) || isinf(loss_cool)
            return Inf
        end

        # Combine: normalize so that each experiment contributes equally,
        # then apply weight w_cool to boost cooling influence
        total = (loss_heat * n_heat + w_cool * loss_cool * n_cool) / (n_heat + w_cool * n_cool)
        return total
    end

end


# ============================================================================
# SECTION F: OPTIMIZATION
# ============================================================================
begin #Optimization

    # Experiment keys
    sim_key_heat = ["E67","E68", "E69", "E70", "E71",
                    "E72", "E73", "E74", "E75", "E76",
                    "E77", "E78", "E79", "E80", "E81"]

    sim_key_cool = ["C69", "C80", "C81"]

    # Initial parameter guesses
    #   γ_C, η_eff, K_lin, χε_rad, a_h, n_exp, χ_L
    p0 = [1.0, 0.70, 1.0, 0.50, 50.0, 0.5, 1.1]

    # Parameter bounds
    lb = [0.3,  0.30, 0.01, 0.01,  1.0, 0.0, 0.5]
    ub = [3.0,  1.00, 10.0, 2.00, 500.0, 1.5, 2.0]

    # Cooling weight (boost cooling influence relative to per-experiment average)
    w_cool = 2.0  # Give cooling 2x weight per experiment since they're rare and informative

    # Combined data tuple
    opt_data = (sim_key_heat, sim_key_cool, w_cool)

    # Compute initial error
    initial_error = lossAll_combined(p0, opt_data)
    println("Initial Error: $initial_error")

    # Setup optimization
    optf = OptimizationFunction(lossAll_combined, Optimization.AutoForwardDiff())
    optprob = Optimization.OptimizationProblem(optf, p0, opt_data, lb=lb, ub=ub)

    # Solve using MLSL (global) + Nelder-Mead (local)
    optsol = solve(optprob, NLopt.GN_MLSL_LDS(),
                   local_method=NLopt.LN_NELDERMEAD(),
                   maxtime=120, local_maxiters=10000)

    println(optsol.retcode)
    pnew = optsol.u
    println("\nFitted Parameters:")
    println("  γ_C      = $(pnew[1])")
    println("  η_eff    = $(pnew[2])")
    println("  K_lin    = $(pnew[3]) W/K")
    println("  χε_rad   = $(pnew[4])")
    println("  a_h      = $(pnew[5])")
    println("  n_exp    = $(pnew[6])")
    println("  χ_L      = $(pnew[7])")
    println("\nFinal Combined Error: $(optsol.objective)")
    println("  Heating-only Error: $(lossAll_heating(pnew, sim_key_heat))")
    println("  Cooling-only Error: $(lossAll_cooling(pnew, sim_key_cool))")

end

# ============================================================================
# SECTION G: DIAGNOSTICS — DERIVED QUANTITIES
# ============================================================================
begin
    println("\n=== Derived Physical Quantities ===")

    # Effective thermal capacitance at reference temperature
    T_ref = 600.0  # K — typical mid-range operating temperature
    C_eff_ref = pnew[1] * ρs * Cps_f(T_ref) * Vs
    println("  C_eff($(T_ref)K) = $(round(C_eff_ref, digits=2)) J/K")
    println("  τ_thermal ≈ C_eff/(K_g+K_loss) — compute per case")

    # Effective NTU and ε at different flow rates
    println("\n  Flow-dependent effectiveness (at T_s = $(T_ref) K):")
    for qlpm_test in [5.0, 10.0, 15.0, 20.0]
        mdot_test = m_dot(qlpm_test)
        T_film_test = (T_ref + Tamb) / 2.0
        cp_test = cpf_f(T_film_test)
        NTU_test = (pnew[5] * Pi * L / cp_test) * mdot_test^(pnew[6] - 1.0)
        ε_test = pnew[7] * (1.0 - exp(-NTU_test))
        K_g_test = mdot_test * cp_test * ε_test
        τ_test = C_eff_ref / (K_g_test + pnew[3])
        Tg_out_test = Tamb + ε_test * (T_ref - Tamb)
        println("    qlpm=$(qlpm_test): NTU=$(round(NTU_test,digits=3)), ε_eff=$(round(ε_test,digits=3)), " *
                "K_g=$(round(K_g_test,digits=3)) W/K, τ=$(round(τ_test,digits=0)) s, " *
                "T_g,out=$(round(Tg_out_test,digits=1)) K")
    end
end

# ============================================================================
# SECTION H: POST-PROCESSING — STEADY-STATE COMPARISON (GAS)
# ============================================================================
begin
    T_steady_gas = DataFrame(sim_id=String[], time=Vector[], T_mod=Vector[], T_exp=Vector[])

    for sm in sim_key_heat
        cond_k = simulation_conditions[sm]
        sel_meas = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tf"), :]
        expdata_tf = sel_meas[:, :temperatures][1]
        time_opt = sel_meas[:, :time][1]
        Ts_avg_init = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tavg_v4"), :temperatures][1][1]

        Ts_mod, Tf_mod = solve_model(pnew, cond_k[Io], cond_k[qlpm], Ts_avg_init, time_opt)
        push!(T_steady_gas, (sm, time_opt, Tf_mod, expdata_tf))
    end

    # Steady-state parity plot (gas)
    T_mod_gas_ss = [T_steady_gas[i, :T_mod][end] for i in 1:length(sim_key_heat)]
    T_exp_gas_ss = [T_steady_gas[i, :T_exp][end] for i in 1:length(sim_key_heat)]

    lims_gas = (300, 850)
    plot_gas_ss = scatter(T_mod_gas_ss, T_exp_gas_ss, label="",
        title="Gas SS Temperature (v4)", xlabel="T_mod (K)", ylabel="T_exp (K)",
        ylims=lims_gas, xlims=lims_gas, aspect_ratio=:equal,
        markersize=6, markerstrokewidth=0)
    plot!(collect(lims_gas), collect(lims_gas), label="1:1", ls=:dash, color=:gray)
    display(plot_gas_ss)
end

# ============================================================================
# SECTION I: POST-PROCESSING — STEADY-STATE COMPARISON (SOLID AVG)
# ============================================================================
begin
    T_steady_solid = DataFrame(sim_id=String[], time=Vector[], T_mod=Vector[], T_exp=Vector[])

    for sm in sim_key_heat
        cond_k = simulation_conditions[sm]
        sel_meas = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tavg_v4"), :]
        expdata_ts = sel_meas[:, :temperatures][1]
        time_opt = sel_meas[:, :time][1]
        Ts_avg_init = expdata_ts[1]

        Ts_mod, Tf_mod = solve_model(pnew, cond_k[Io], cond_k[qlpm], Ts_avg_init, time_opt)
        push!(T_steady_solid, (sm, time_opt, Ts_mod, expdata_ts))
    end

    # Steady-state parity plot (solid)
    T_mod_solid_ss = [T_steady_solid[i, :T_mod][end] for i in 1:length(sim_key_heat)]
    T_exp_solid_ss = [T_steady_solid[i, :T_exp][end] for i in 1:length(sim_key_heat)]

    lims_solid = (300, 1150)
    plot_solid_ss = scatter(T_mod_solid_ss, T_exp_solid_ss, label="",
        title="Solid Avg SS Temperature (v4)", xlabel="T_mod (K)", ylabel="T_exp (K)",
        ylims=lims_solid, xlims=lims_solid, aspect_ratio=:equal,
        markersize=6, markerstrokewidth=0)
    plot!(collect(lims_solid), collect(lims_solid), label="1:1", ls=:dash, color=:gray)
    display(plot_solid_ss)

    # Combined SS parity plot
    display(plot(plot_gas_ss, plot_solid_ss, layout=(1, 2), size=(1000, 450)))
end

# ============================================================================
# SECTION J: TRANSIENT COMPARISON PLOTS
# ============================================================================
begin
    """Plot transient comparison for a specific experiment."""
    function plt_case_v4(sm, params; show_cooling=false)
        if show_cooling
            cond_k = simulation_conditions_cooling[sm]
            Io_val = 0.0
            meas_df = measurements_cooling
        else
            cond_k = simulation_conditions[sm]
            Io_val = cond_k[Io]
            meas_df = measurements
        end
        qlpm_val = cond_k[qlpm]
        Tinit_val = cond_k[Tinit]

        time_exp = meas_df[(meas_df.simulation_id .== sm) .& (meas_df.obs_id .== "_Tavg_v4"), :time][1]
        Ts_exp = meas_df[(meas_df.simulation_id .== sm) .& (meas_df.obs_id .== "_Tavg_v4"), :temperatures][1]
        Tf_exp = meas_df[(meas_df.simulation_id .== sm) .& (meas_df.obs_id .== "_Tf"), :temperatures][1]
        Tinit_val = Ts_exp[1]  # Use initial T_s,avg as IC

        Ts_mod, Tf_mod = solve_model(params, Io_val, qlpm_val, Tinit_val, time_exp)

        mode_str = show_cooling ? "Cooling" : "Heating"
        plt_t = plot(title="$mode_str: $sm  (q=$(qlpm_val) L/min)",
                     xlabel="Time (s)", ylabel="Temperature (K)",
                     legend=:right, ylim=(250, 1100))
        scatter!(time_exp, Ts_exp, label="T_s,avg exp", ms=2, color=colors[1], alpha=0.5)
        scatter!(time_exp, Tf_exp, label="T_g,out exp", ms=2, color=colors[2], alpha=0.5)
        plot!(time_exp, Ts_mod, label="T_s,avg mod", lw=2, color=colors[1])
        plot!(time_exp, Tf_mod, label="T_g,out mod", lw=2, color=colors[2])
        return plt_t
    end

    # Plot all heating experiments
    println("\n=== Plotting Heating Experiments ===")
    map(x -> display(plt_case_v4(x, pnew)), sim_key_heat)

    # Plot cooling experiments
    println("\n=== Plotting Cooling Experiments ===")
    map(x -> display(plt_case_v4(x, pnew; show_cooling=true)), sim_key_cool)
end

# ============================================================================
# SECTION K: EFFECTIVENESS AND h_eff SURFACE PLOTS
# ============================================================================
begin
    T_range = range(300, 1000, length=50)
    flow_rates = [1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0]

    # h_eff = a_h * ṁ^n — this is the effective overall HTC
    plot_h = plot(xlabel="T_film (K)", ylabel="ε_eff", title="Gas effectiveness ε_eff vs T",
                  legend=:topright)
    for fl in flow_rates
        mdot_fl = m_dot(fl)
        ε_vals = [begin
            T_film = (T + Tamb) / 2.0
            cp = cpf_f(T_film)
            NTU = (pnew[5] * Pi * L / cp) * mdot_fl^(pnew[6] - 1.0)
            pnew[7] * (1.0 - exp(-NTU))
        end for T in T_range]
        plot!(T_range, ε_vals, label="$(fl) L/min")
    end
    display(plot_h)

    # NTU vs flow rate
    fl_range = range(1, 30, length=30)
    plot_ntu = plot(xlabel="Flow rate (L/min)", ylabel="NTU_eff",
                    title="NTU vs Flow Rate", legend=:topright)
    for T_test in [400.0, 600.0, 800.0, 1000.0]
        ntu_vals = [begin
            mdot_fl = m_dot(fl)
            T_film = (T_test + Tamb) / 2.0
            cp = cpf_f(T_film)
            (pnew[5] * Pi * L / cp) * mdot_fl^(pnew[6] - 1.0)
        end for fl in fl_range]
        plot!(fl_range, ntu_vals, label="T_s=$(T_test) K")
    end
    display(plot_ntu)
end

# ============================================================================
# SECTION L: METRICS COMPUTATION AND EXPORT
# ============================================================================
begin
    println("\n=== Computing Metrics (0D v4) ===")

    # Helper function for t90
    function get_t90(time, temp)
        t_init = temp[1]
        t_ss = temp[end]
        target = t_init + 0.9 * (t_ss - t_init)
        # For heating: target > t_init; for cooling: target < t_init
        if t_ss > t_init  # heating
            idx = findfirst(t -> t >= target, temp)
        else  # cooling
            idx = findfirst(t -> t <= target, temp)
        end
        return idx !== nothing ? time[idx] : time[end]
    end

    # Output CSV path
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
        RunID = Float64[],
        Case = String[],
        Type = String[],  # "heating" or "cooling"
        Ts_avg_SS_sim = Float64[],
        Ts_avg_SS_exp = Float64[],
        dT_Ts_avg = Float64[],
        T3_SS_sim = Float64[],
        T3_SS_exp = Float64[],
        dT_T03 = Float64[],
        R_leak_sim = Float64[],
        R_leak_exp = Float64[],
        t90_sim_Ts_s = Float64[],
        t90_exp_Ts_s = Float64[],
        dt90_Ts_s = Float64[],
        t90_sim_T3_s = Float64[],
        t90_exp_T3_s = Float64[],
        dt90_T3_s = Float64[],
        Gap_ss_sim = Float64[],
        Gap_ss_exp = Float64[],
        dGap_ss = Float64[],
        ε_eff = Float64[],
        NTU_eff = Float64[]
    )

    # Process heating experiments
    for sm in sim_key_heat
        cond_k = simulation_conditions[sm]
        Io_val = cond_k[Io]
        qlpm_val = cond_k[qlpm]
        Tinit_val = cond_k[Tinit]

        time_exp = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tavg_v4"), :time][1]
        Ts_exp = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tavg_v4"), :temperatures][1]
        Tf_exp = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tf"), :temperatures][1]
        Tinit_solid = Ts_exp[1]  # Use initial T_s,avg as IC

        Ts_mod, Tf_mod = solve_model(pnew, Io_val, qlpm_val, Tinit_solid, time_exp)

        # Steady-state values
        Ts_ss_sim = Ts_mod[end]
        Tf_ss_sim = Tf_mod[end]
        Ts_ss_exp = Ts_exp[end]
        Tf_ss_exp = Tf_exp[end]

        dT_Ts = Ts_ss_sim - Ts_ss_exp
        dT_T3 = Tf_ss_sim - Tf_ss_exp

        R_leak_sim = (Tf_ss_sim - Tamb) / (Ts_ss_sim - Tamb)
        R_leak_exp = (Tf_ss_exp - Tamb) / (Ts_ss_exp - Tamb)

        t90_sim_Ts = get_t90(time_exp, Ts_mod)
        t90_exp_Ts = get_t90(time_exp, Ts_exp)
        dt90_Ts = t90_sim_Ts - t90_exp_Ts

        t90_sim_T3 = get_t90(time_exp, Tf_mod)
        t90_exp_T3 = get_t90(time_exp, Tf_exp)
        dt90_T3 = t90_sim_T3 - t90_exp_T3

        gap_sim = Tf_ss_sim - Ts_ss_sim
        gap_exp = Tf_ss_exp - Ts_ss_exp
        dgap = gap_sim - gap_exp

        # Compute ε_eff and NTU at SS
        mdot_val = m_dot(qlpm_val)
        T_film_ss = (Ts_ss_sim + Tamb) / 2.0
        cp_ss = cpf_f(T_film_ss)
        NTU_ss = (pnew[5] * Pi * L / cp_ss) * mdot_val^(pnew[6] - 1.0)
        ε_ss = pnew[7] * (1.0 - exp(-NTU_ss))

        push!(results_df, (
            next_run_id, sm, "heating",
            Ts_ss_sim, Ts_ss_exp, dT_Ts,
            Tf_ss_sim, Tf_ss_exp, dT_T3,
            R_leak_sim, R_leak_exp,
            t90_sim_Ts, t90_exp_Ts, dt90_Ts,
            t90_sim_T3, t90_exp_T3, dt90_T3,
            gap_sim, gap_exp, dgap,
            ε_ss, NTU_ss
        ))
    end

    # Process cooling experiments
    for sm in sim_key_cool
        cond_k = simulation_conditions_cooling[sm]
        qlpm_val = cond_k[qlpm]
        Tinit_val = cond_k[Tinit]

        time_exp = measurements_cooling[(measurements_cooling.simulation_id .== sm) .& (measurements_cooling.obs_id .== "_Tavg_v4"), :time][1]
        Ts_exp = measurements_cooling[(measurements_cooling.simulation_id .== sm) .& (measurements_cooling.obs_id .== "_Tavg_v4"), :temperatures][1]
        Tf_exp = measurements_cooling[(measurements_cooling.simulation_id .== sm) .& (measurements_cooling.obs_id .== "_Tf"), :temperatures][1]
        Tinit_solid = Ts_exp[1]  # Use initial T_s,avg as IC

        Ts_mod, Tf_mod = solve_model(pnew, 0.0, qlpm_val, Tinit_solid, time_exp)

        Ts_ss_sim = Ts_mod[end]
        Tf_ss_sim = Tf_mod[end]
        Ts_ss_exp = Ts_exp[end]
        Tf_ss_exp = Tf_exp[end]

        dT_Ts = Ts_ss_sim - Ts_ss_exp
        dT_T3 = Tf_ss_sim - Tf_ss_exp

        R_leak_sim = abs(Ts_ss_sim - Tamb) > 1.0 ? (Tf_ss_sim - Tamb) / (Ts_ss_sim - Tamb) : NaN
        R_leak_exp = abs(Ts_ss_exp - Tamb) > 1.0 ? (Tf_ss_exp - Tamb) / (Ts_ss_exp - Tamb) : NaN

        t90_sim_Ts = get_t90(time_exp, Ts_mod)
        t90_exp_Ts = get_t90(time_exp, Ts_exp)
        dt90_Ts = t90_sim_Ts - t90_exp_Ts

        t90_sim_T3 = get_t90(time_exp, Tf_mod)
        t90_exp_T3 = get_t90(time_exp, Tf_exp)
        dt90_T3 = t90_sim_T3 - t90_exp_T3

        gap_sim = Tf_ss_sim - Ts_ss_sim
        gap_exp = Tf_ss_exp - Ts_ss_exp
        dgap = gap_sim - gap_exp

        mdot_val = m_dot(qlpm_val)
        T_film_ss = (Ts_ss_sim + Tamb) / 2.0
        cp_ss = cpf_f(T_film_ss)
        NTU_ss = (pnew[5] * Pi * L / cp_ss) * mdot_val^(pnew[6] - 1.0)
        ε_ss = pnew[7] * (1.0 - exp(-NTU_ss))

        push!(results_df, (
            next_run_id, sm, "cooling",
            Ts_ss_sim, Ts_ss_exp, dT_Ts,
            Tf_ss_sim, Tf_ss_exp, dT_T3,
            R_leak_sim, R_leak_exp,
            t90_sim_Ts, t90_exp_Ts, dt90_Ts,
            t90_sim_T3, t90_exp_T3, dt90_T3,
            gap_sim, gap_exp, dgap,
            ε_ss, NTU_ss
        ))
    end

    # Print results
    println(results_df)

    # Save to CSV
    if file_exists
        CSV.write(csv_path, results_df, append=true)
        println("Appended metrics to $csv_path with RunID = $next_run_id")
    else
        CSV.write(csv_path, results_df)
        println("Created $csv_path with RunID = $next_run_id")
    end
end

