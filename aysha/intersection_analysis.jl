# =============================================================
#  intersection_analysis.jl   (SiC monolith – T8/T9 method)
# -------------------------------------------------------------
#  Purpose :
#    • Read 15 experimental runs (3 irradiance levels × 5 flow‑rates)
#      containing time–temperature profiles at two wall locations
#      "T8" (upstream) and "T9" (downstream).
#    • Fit a single‑mode exponential tail for each trace → τ₈, τ₉.
#    • Optimise a common SiC conductivity  k_SiC  and a gas‑side
#      Nusselt correlation  h = k_g/D_h · C · Re^m · Pr^n  (with
#      n fixed, C & m fitted) by minimising the error between the
#      measured τ’s and the eigen‑value model of a cylindrical wall.
#
#  Geometry  (edit to match rig):
#      L_s8  = 0.0040   m   # T8 depth from hot face
#      L_s9  = 0.0160   m   # T9 depth from hot face
#      D_h   = 0.0017   m   # hydraulic diameter of a single channel
#
#  Input CSV columns (one file for all runs):
#      runID         : 1‥15
#      irradianceID  : 1‥3         (low / mid / high)
#      Re            : Reynolds number of the gas
#      Pr            : Prandtl number  (film temperature)
#      k_g           : Gas thermal conductivity (W/m/K)
#      time          : seconds
#      T8            : K
#      T9            : K
#
#  Output printed to REPL:
#      • fitted k_SiC (W/m/K)
#      • fitted C and m of the Nu law (n fixed at 0.40)
#      • table of predicted h for the 5 flow‑rates of each irradiance
#      • RMS error on τ’s
# =============================================================

using CSV, DataFrames, Statistics, LinearAlgebra
using Interpolations, Optim, Roots

# ------------------------- user editable constants ------------------------
const L_s8 = 0.011     # m
const L_s9 = 0.058     # m
const ρ_SiC  = 3100.0   # kg/m³
const c_pSiC = 900.0    # J/kg/K  (film temperature 700–800 °C)
const n_Pr   = 0.40     # fixed Pr exponent in Nu law
const λ_tol  = 1e-6     # eigen‑value solver tolerance
const Tamb = 298.     # K
const D_h = 0.0017   # m  (edit if necessary)


    # Receiver dimensions
    w_t = 20.e-3  # Width of the receiver in meters
    A_frt = w_t * w_t  # Total frontal area of the receiver (m^2)
    L = 137e-3 #m Length
    A_s_p = w_t * L * 4 #total area solid periphery m2

    # Channel dimensions
    w_chnl = 1.7e-3  # Width of a single channel (m)
    A_chnl_frt = w_chnl * w_chnl  # Frontal area of a single channel (m^2)
    A_chnl_p = w_chnl * L * 4  # Periphery area of a single channel (m^2)

    # Channel array
    n_chnl = 10 * 10  # Number of channels (10x10 array)
    A_chnl_frt_all = A_chnl_frt * n_chnl  # Total frontal area of all channels (m^2)
    Lc = 4 * (w_t * w_t) / (4 * w_t) #hydraulic receiver diameter

    # Solid geometrical calculations
    A_frt_solid = A_frt - A_chnl_frt_all  # Frontal area of solid (m^2)
    A_exchange = A_chnl_p * n_chnl  # Total contact area between fluid and solid (m^2)

    # Volume calculations
    Vs = A_frt_solid * L  # Corrected solid volume (m^3)
    Vf = n_chnl * A_chnl_frt * L  # Fluid volume in the channels (m^3)


Arqplm=[15.27, 12.50, 10.50, 9.10, 7.12,
    18.34, 13.16, 9.03, 6.95, 4.53,
    13.85, 10.02, 8.04, 6.62, 4.53]
RH = 0.4
ρf_f(T) = 352.716 * T^-1 * (1+RH) /(1+ 1.609 * RH) #COMSOL + https://www.engineeringtoolbox.com/density-air-d_680.html with some error for w vs RH
kf_f(T) = -0.00227583562 + 1.15480022e-4 * T^1 - 7.90252856e-8 * T^2 + 4.11702505e-11 * T^3 - 7.43864331e-15 * T^4 #COMSOL
#Cpw_f(T) = 1745.96354 + 0.185114553 * T^1 + 6.19448731e-4 * T^2 - 3.0267851e-7 * T^3 + 4.19053122e-11 * T^4 #COMSOL
cpf_f(T) = (1.93e-10 * T^4 - 8.0e-7 * T^3 + 1.14e-3 * T^2 - 4.49e-1 * T + 1.06e3) 
μf_f(T) = -8.38278E-7 + 8.35717342e-8 * T^1 - 7.69429583e-11 * T^2 + 4.6437266e-14 * T^3 - 1.06585607e-17 * T^4
m(qlpm) = ρf_f(Tamb) * (qlpm / 60.0 / 1000.0) #kg/s mass at the flowmeter conditions T=25oC
G(qlpm) = m(qlpm) / A_chnl_frt_all #kg/2/m2 mass at the flowmeter conditions T=25oC
Re_f(qlpm, T) = G(qlpm) * Lc / μf_f(T) #use of flux instead of u*ρf
Pr_f(T) = (cpf_f(T) * μf_f(T)) / kf_f(T)

# ---------------------- helper: first eigen‑value λ₁ ----------------------
"λ1(Bi) : solve λ tan λ + Bi = 0 for 0<λ<π/2"
function λ1(Bi)
    fλ(λ) = λ * tan(λ)+Bi
    # Bracket: λ→0  => f≈Bi>0 ; λ→π/2 => f→∞ ; root lives in (0,π/2)
    #return find_zero(f, (1e-6, 1.5607963267948966); atol=λ_tol)
    return find_zero(fλ, (0, 1.5607963267948966); atol=λ_tol)
end

# ---------------------- τ model for given h, k_SiC -----------------------
α_SiC(k) = k / (ρ_SiC * c_pSiC)
function τ_model(h, k)
    Bi8 = h*L_s8/k
    Bi9 = h*L_s9/k
    α   = α_SiC(k)
    
    τ8 = 0
    τ9 = 0
    if Bi8 > 0.1
        λ8  = λ1(Bi8)
        τ8  = L_s8^2 / (α*λ8^2)
    else
        #Vs8 = Vs * L_s8 / L
        #As8 = A_exchange * L_s8 / L
        τ8 = ρ_SiC * c_pSiC / h * (Vs/A_exchange * L_s8/L)
    end
    if Bi9 > 0.1
        λ9  = λ1(Bi9)
        τ9  = L_s9^2 / (α*λ9^2)
    else
        #Vs9 = Vs * L_s9 / L
        #As9 = A_exchange * L_s9 / L
        #τ9 = ρ_SiC * c_pSiC * Vs9 / h / As9
        τ9 = ρ_SiC * c_pSiC / h * (Vs/A_exchange * L_s9/L)
    end
    return τ8, τ9
end

# ---------------------- fit τ from one trace -----------------------------
function fit_tau(t::Vector{Float64}, Temp::Vector{Float64}; tail_frac = 0.4)
    N   = length(t)
    T_inf, Nt = findmax(Temp)
    Nt  = Int(round(tail_frac*Nt))
    idx = 10:Nt
    #T_inf = maximum(Temp[idx])
    ΔT    = T_inf .- Temp[idx]
    ΔT=ΔT[1:end-2]
    lnΔ   = log.(ΔT)
    tsub  = t[idx]
    tsub = tsub[1:end-2]
    # linear regression lnΔ = a + b t  →  τ = -1/b
    A = [sum(tsub.^2)  sum(tsub);  sum(tsub)  Nt-10]
    X = [sum(tsub .* lnΔ); sum(lnΔ)]
    Y = zero(X)
    b, a = LinearAlgebra.ldiv!(Y,qr(A),X)

    τ = -1 / b
    return τ
end

# ---------------------- load CSV & compute τ's ---------------------------
function load_data(measurements)
    #df   = CSV.read(path, DataFrame)
    gexp = unique(measurements.simulation_id)
    out  = Vector{NamedTuple}(undef, length(gexp))
    for (i, g) in enumerate(gexp)
        t = measurements[(measurements.simulation_id .== g) .& (measurements.obs_id .== "_T2"), :time][1]
        T9 = measurements[(measurements.simulation_id .== g) .& (measurements.obs_id .== "_T9"), :temperatures][1]
        #Tprobe = measurements[(measurements.simulation_id .== g) .& (measurements.obs_id .== "_T2"), :temperatures][1]
        T8 = measurements[(measurements.simulation_id .== g) .& (measurements.obs_id .== "_T8"), :temperatures][1]

        out[i] = (t = t,
                   T9 = T9,
                   T8 = T8)
    end
    return gexp, out
end
function load_and_fit()
    include("import_exp_0D.jl") #import the experimental data
    runs, out = load_data(measurements) 
    #df = CSV.read(csv_path, DataFrame)
    #runs = groupby(df, :runID)
    τ8  = Float64[]; τ9 = Float64[]; Re = Float64[]; Pr = Float64[]; k_g = Float64[];
    #limit = end
    for (i, g) in enumerate(out)
        #println(i)
        push!(Re, Re_f(Arqplm[i], 500))
        push!(Pr, Pr_f(500))
        push!(k_g, kf_f(500))
        push!(τ8, fit_tau(g.t, g.T8))
        push!(τ9, fit_tau(g.t, g.T9))
    end
    return (; τ8, τ9, Re, Pr, k_g)
end

# ---------------------- objective for Optim.jl ---------------------------
# params = [C, m, k_SiC]
function objective(params, data)
    oC, om, ok = params
    err = 0.0
    for (i, Re_i) in enumerate(data.Re)
        h_i = data.k_g[i] / D_h * oC * Re_i^om * data.Pr[i]^n_Pr
        τ8_pred, τ9_pred = τ_model(h_i, ok)
        err += (τ8_pred - data.τ8[i])^2 + (τ9_pred - data.τ9[i])^2
    end
    return err
end
function objective(params)
    return objective(params, data)
end


# ---------------------- main routine ------------------------------------
#if abspath(PROGRAM_FILE) == @__FILE__

    #csv_file = "T8_T9_profiles.csv"   # <- edit path
    data     = load_and_fit()

    # initial guesses: C 0.02, m 0.80, k_SiC 28 W/m/K
    p0   = [0.02, 0.80, 28.0]
    low  = [0.001, 0.50, 10.0]
    high = [0.20,  1.20, 240.0]

    println("Optimising ... this may take a minute")
    res = optimize(p -> objective(p, data), low, high, p0, Fminbox(LBFGS()), Optim.Options(f_reltol = 1e-8, f_abstol = 1e-8))
   #res = optimize(objective, p0, LBFGS(), Optim.Options(f_reltol = 1e-9, iterations = 100000))

    p̂  = Optim.minimizer(res)
    Ĉ, m̂, k̂ = p̂

    print("\n=========== fitted parameters ===========")
    print("C (Nu prefactor)      = $Ĉ")
    print("m (Re exponent)       = $m̂")
    print("k_SiC (700–800 °C)    = $k̂  W/m/K", )
    print("========================================\n")

    # table of h predictions
    print("Flow‑rate table (by run):")
    print("%4s %12s %12s\n", "run", "Re", "h [W/m²K]")
    for i in 1:length(data.Re)
        h_i = data.k_g[i] / D_h * Ĉ * data.Re[i]^m̂ * data.Pr[i]^n_Pr
        println("$i, $(data.Re[i]), $h_i")
    end

    rms = sqrt(objective(p̂, data) / (2*length(data.Re)))
    print("RMS error on τ's  = $rms")
#end
