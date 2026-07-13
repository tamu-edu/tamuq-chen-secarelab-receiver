# =============================================================
# insulation_inverse.jl  —  revised 2025‑07‑13
# -------------------------------------------------------------
# Inverse identification of the temperature‑dependent thermal
# conductivity  k(T)=k₀[1+a₁(T−Tf)+a₂(T−Tf)²]  and volumetric heat
# capacity ρc  of a ceramic‑fibre insulation layer, using the
# transient temperature recorded at
#   • x = 0          → interface with the hot SiC monolith  (Tin)
#   • x = x_probe    → thermocouple embedded in the felt (Tprobe)
# for 15 experiments that differ only in heater power.  The
# aluminium shell at x = L rejects heat to the laboratory air
# by natural convection + radiation.
#
# Compared with the first version:  (i) left boundary is now a
# time‑dependent Dirichlet condition (Tin(t)), not an imposed
# flux;  (ii) conductivity depends on temperature;  (iii) heat
# fluxes q″ are *not* parameters — they become model outputs.
# =============================================================
begin
    using CSV, DataFrames, Statistics, LinearAlgebra
    using DifferentialEquations, Interpolations
    using Optim
    using Plots, XLSX, Dates, DelimitedFiles
end
# -------------------------------------------------------------
# 1.  DATA HANDLING
# -------------------------------------------------------------
"""
    load_data(path) -> Vector{NamedTuple}

Input CSV must contain columns:
    expID   : 1‥15 (Int)
    time    : seconds (Float64)
    Tin     : K       (Float64)  — monolith / felt interface
    Tprobe  : K       (Float64)  — depth x_probe in felt
Returns a vector of 15 named tuples (t, Tin, Tp).
"""
function load_data(measurements)
    #df   = CSV.read(path, DataFrame)
    gexp = measurements.simulation_id
    out  = Vector{NamedTuple}(undef, length(gexp))
    for (i, g) in enumerate(gexp)
        t = measurements[(measurements.simulation_id .== g) .& (measurements.obs_id .== "_T2"), :time][1]
        Tin = measurements[(measurements.simulation_id .== g) .& (measurements.obs_id .== "_T9"), :temperatures][1]
        Tprobe = measurements[(measurements.simulation_id .== g) .& (measurements.obs_id .== "_T2"), :temperatures][1] 

        out[i] = (t     = t,
                   Tin   = Tin,
                   Tp    = Tprobe)
    end
    return out
end

# -------------------------------------------------------------
# 2.  THERMAL‑PROPERTY MODEL
# -------------------------------------------------------------
"""
    k_of_T(T, k0, a1, a2, Tf)
Piecewise‑polynomial positive‑definite conductivity.
"""
@inline k_of_T(T, k0, a1, a2, Tf) = k0 * (1 + a1*(T-Tf) + a2*(T-Tf)^2)

# -------------------------------------------------------------
# 3.  FORWARD SOLVER (1‑D IMPLICIT FD, VARIABLE k)
# -------------------------------------------------------------
function simulate_probe(temp_time, k0, a1, a2, ρc, x_probe, L;
                        Tf = 300.0, Tamb = 298.15, h = 10.0,
                        eps = 0.15, nx = 61)
    """Return T(x_probe,t) for a single experiment.
    temp_time : (t, Tin) named‑tuple
    """
    tvec   = temp_time.t
    Tinvec = temp_time.Tin
    dx     = L/(nx-1)
    idx_p  = round(Int, x_probe/dx) + 1

    # Tin(t)  — linear interpolation on provided knots
    Tin_itp = LinearInterpolation(tvec, Tinvec,extrapolation_bc = Interpolations.Flat())

    function heat_ode!(dT, T, p, t)
        # Update boundary node to prescribed Dirichlet value
        T_left = Tin_itp(t)
        T[1]   = T_left
        dT[1]  = 0.0

        # Interior nodes 2 … nx-1  (variable k)
        for i in 2:nx-1
            T_im, T_i, T_ip = T[i-1], T[i], T[i+1]
            k_im = k_of_T( (T_im+T_i)/2,  k0, a1, a2, Tf)
            k_ip = k_of_T( (T_i +T_ip)/2, k0, a1, a2, Tf)
            dT[i] = (k_ip*(T_ip-T_i) - k_im*(T_i-T_im)) / (ρc*dx^2)
        end

        # Right boundary, convection+rad to Tamb
        T_s  = T[nx]
        T_sm = T[nx-1]
        k_s  = k_of_T( (T_s+T_sm)/2, k0, a1, a2, Tf)
        flux_cond = -k_s*(T_s - T_sm)/dx
        flux_env  = h*(T_s - Tamb) + eps*5.670374419e-8*(T_s^4 - Tamb^4)
        dT[nx] = 2*(flux_cond - flux_env)/(ρc*dx)   # minus sign already in flux_cond
    end

    # ODE solve
    T0_vec = fill(Tinvec[1], nx)           # init with first Tin
    prob   = ODEProblem(heat_ode!, T0_vec, (tvec[1], tvec[end]))
    sol    = solve(prob, Tsit5(); saveat = tvec, reltol = 1e-6, abstol = 1e-6)
    return sol[idx_p, :]
end

# -------------------------------------------------------------
# 4.  GLOBAL OBJECTIVE (LEAST SQUARES OVER 15 EXPERIMENTS)
# -------------------------------------------------------------
function objective(p, data, x_probe, L, Tf, Tamb, h, eps)
    k0, a1, a2, ρc = p
    err = 0.0
    for ex in data
        Tp_sim = simulate_probe(ex, k0, a1, a2, ρc, x_probe, L;
                                Tf = Tf, Tamb = Tamb, h = h, eps = eps)
        err   += sum((Tp_sim .- ex.Tp) .^ 2)
    end
    return err
end

# -------------------------------------------------------------
# 5.  MAIN SCRIPT
# -------------------------------------------------------------
#if abspath(PROGRAM_FILE) == @__FILE__

    # ---------------- data ----------------------------------
    include("import_exp_0D.jl") #import the experimental data

    dataset = load_data(measurements)   

begin

    # ---------------- user inputs ---------------------------
    csv_path    = "transient_Tin_Tp.csv"   # edit accordingly
    L_felt      = 0.15    # m, thickness of insulation
    x_probe     = 0.040 - 0.018    # m, depth of thermocouple
    Tf_ref      = 450.0    # K, centre of fitted k(T) window
    Tamb        = 298.15   # K, lab air
    h_nat       = 10.0     # W m⁻² K⁻¹, est. natural conv.
    eps_shell   = 0.15     # emissivity of brushed Al shell

    # ---------------- optimisation --------------------------
    p0   = [0.10,   1.0e-3, 0.0, 9.6e4]   # k0, a1, a2, ρc
    low  = [0.05,  -5e-3,  -5e-6,  2e4]
    high = [0.25,   5e-3,   5e-6,  2e5]

    println("\nFitting k(T) and ρc …\n")
    cost(p) = objective(p, dataset, x_probe, L_felt,
                        Tf_ref, Tamb, h_nat, eps_shell)
    #res = optimize(cost, low, high, p0, Fminbox(LBFGS()); f_reltol = 1e-6, f_calls_limit=1000)
    #res = optimize(cost, p0, LBFGS(), Optim.Options(f_reltol = 1e-6, iterations = 1000))
    
    p̂ = Optim.minimizer(res)

    println("\n=========== BEST‑FIT PARAMETERS ===========")
    printf("k₀          = %.4f W/m/K\n", p̂[1])
    printf("a₁          = %.3e 1/K\n", p̂[2])
    printf("a₂          = %.3e 1/K²\n", p̂[3])
    printf("ρ·c         = %.0f J/m³/K\n", p̂[4])
    println("==========================================\n")
end
