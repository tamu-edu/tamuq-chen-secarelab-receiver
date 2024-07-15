begin #libraries
    using MethodOfLines
    using ModelingToolkit
    using DomainSets, OrdinaryDiffEq
    using NonlinearSolve, DifferentialEquations, DataFrames
    using Plots, XLSX, Statistics, Symbolics, Interpolations
    using LsqFit, LossFunctions
    using Optimization, OptimizationNLopt, Symbolics, ForwardDiff
end
begin #define parameters
    L = 137e-3 #m
    #α = keff/ρ*Cp #kW/m2.K 
    Tamb = (22.448 + 273.15) #K (same for all exp) 
    deltax = 0.002795918367346939 #discretization (m)
    w_t = 19 / 1000 #m
    A_t = w_t * w_t #m2 - for the whole receiver (19x19mm2)
    #A_st = 176e-6 #total area - solid m2
    #A_ft = 49e-6 # total area - fluid m2
    w_t = 19.e-3 #m
    A_frt = w_t * w_t #m2 - for the whole receiver (19x19mm2)
    A_s_p = w_t * L * 4 #total area solid periphery m2
    w_chnl = 1.5e-3 #m
    A_chnl_p = w_chnl * L * 4  #m2 channel periphery
	A_chnl_frt = w_chnl * w_chnl #m2 channel front
	n_chnl = 10*10
    A_exchange = A_chnl_p * n_chnl #m2 - contact area between fluid and solid
	A_chnl_frt_all = A_chnl_frt * n_chnl #m2 all frontal area of channels
    Vs =  A_frt * L #m3
	Vf = n_chnl * A_chnl_frt * L 
    # channel_w = 1.5 / 1000 #m
    # A_channel = channel_w * channel_w  #m2 (1.5x1.5mm2) 
    # n_channel = 100
    # A_exchange = 162.5e-6 #m2 - contact area between fluid and solid
    # Vs = A_st * L #* deltax
    # Vf = A_ft * L #* deltax
    qlpm = 7.12 #lpm
    q = qlpm / (1000 * 60) #m3/s
    ρ = 1.2 #kg/m3 - density of air at lab conditions
    m = q * ρ #kg/s
    #V= q/(A_channel*n_channel) #m/s - using the area of the whole receiver (18x18mm2)
    V = 0.57 #m/s (calculated from excel sheet and COMSOL)
    th_s = 0.4e-3 #m
    th_f = 0.7e-3 #m
    #Qv= I0*exp(-1000*x)  #kW/m2 - (K extinction coefficient taken from Howell and Hendrick paper pg.86, measure pore diameter and fix)
    ϵ = 0.8
    σ = 5.17e-8 #W/m2.K^4 Stefan-Boltzmann constant
    Lc = 4 * (w_t * w_t) / (4 * w_t)
    kf = 0.056 #W/m.K
    ρf = 0.5 #kg/m3
    Cpf = 1090 #J/kg.K
    mu = 2.0921e-5 #Pa.s
    e = 0.425
    Af = n_chnl * A_chnl_frt #m2
    #Gz = (w_t / L) * Re * Pr
    # A = 3.657
    # B = 0.5272
    # C = 55.1239
    # n = 0.3056
    #Nu = A*(1+(B*((Gz)^n)*exp(-C/Gz)))
    #Nu = 3.657 
    Lc = 4 * (w_t * w_t) / (4 * w_t)
    w = 1.5e-3 #width of channel (m)
    Vi = w * w * n_chnl * L #m3
    Av = 4 * (w * L) / (w^2 * L) #specific area (m-1)
    hext = 10 #W/m2.K
    kins = 0.078 #W/m*K
    r0 = 23 / 1000 #m
    r_ins = 42 / 1000 #m
    r_H = 4 * (w_t * w_t) / (4 * w_t) #hydraulic receiver diameter
    em = 0.8 #emissivity
    aCp = 4.4
end;
#for interpolations 
#1. extract T2 data
#Exp 71 
begin
    D3 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0071_T2.xlsx")["Sheet 1 - Data_FPT0071_T2"]["A3:D7087"]
    y1d3_data = D3[:, 2] .+ 273.0 #T2 (insulation)
    x1 = D3[:, 1]
    #2.create interpolation function
    Tins = LinearInterpolation(x1, y1d3_data, extrapolation_bc=y1d3_data[end])
    Tins_f(t) = Tins(t)
    @register_symbolic Tins_f(t)
end

begin 
    x11 = 0.0001:0.001383838383838384:0.137 #T2 (insulation)
    #Re_f(qlpm) = (ρf * (qlpm/60/1000/Af) * w_t) / mu
    Re = (ρf * V * w_t) / mu
    Pr = (Cpf * mu) / kf
    Gz = (1 ./ x11) * Re * Pr * w_t
    #2.create interpolation function
    Gz_ = LinearInterpolation(x11, Gz)
    Gz_f(x) = Gz_(x)
    @register_symbolic Gz_f(x)
end
begin
    # Parameters, variables, and derivatives for system 1
    @variables t x
    @parameters A B C #ks h_average A n
    @parameters Io qlpm
    @variables Ts(..) Tf(..)
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2
    ks = (48.78) * 1.97 * ((1 - e)^1.5) #W/m.K
    ρs = 3200  #kg/m3
    Cps = 1290  #J/kg*K
    #Nu = A*(1+(B*((Gz_f(x))^n)*exp(-C/Gz_f(x))))
    #Nu = A*(Re^B)*(Pr^C)
    #Cps(Ts) = (0.27+0.135e-4*(Ts)-9720*((Ts)^-2)+0.204e-7*((Ts)^2))/1000 #kJ/kg*K from manufacturer data
    #p_opt = [A => 2., B => 0.5, n=> 0.5, C=> 20.]
    p_opt = [A => 250., B => 0.5, C=> 5.]
    p_cond = [Io => 456000.0, qlpm => 15.27]
    p_math = vcat(p_opt, p_cond)
    #h_average = hfa * (qlpm^hfn)
    #h_average = (Nu * kf) / Lc
    Re_f(qlpm) = (ρf * (qlpm/60/1000/Af) * w) / mu
    Nu_f(qlpm) = A*(Re_f(qlpm)^B)*(Pr^C)
    h_avg_f(qlpm) = (Nu_f(qlpm) * kf) / Lc
    
 
    # MOL Discretization parameters for system 1
    x_max1 = L
    x_min1 = 0.0
    t_min = 0.0
    t_max = 7084.0

    nc1 = 100

    x_num1 = range(x_min1, x_max1, length=nc1)

    dx = (x_max1 - x_min1) / (nc1 - 1)


    # PDE equation for system 1

    # eq1 = [
    #     Vs * (ρs*Cps) * Dt(Ts(t, x)) ~ Vs * (ks) * Dxx(Ts(t, x)) - ((h_average) * Av * Vi * ((Ts(t, x)) - Tf(t, x))) .- (kins * (r / r0) .* (Ts(t, x) .- Tins_f(t)) * A_t / (r - r0)),
    #     Vf * ρf * Cpf * Dt(Tf(t, x)) ~ Vf * kf * Dxx(Tf(t, x)) - Vf * ρf * Cpf * V * Dx(Tf(t, x)) + (h_average) * Av * Vi * ((Ts(t, x) - Tf(t, x)))
    # ]
    eq1 = [
        (1-e) * (aCp * ρs * Cps) * Vs * Dt(Ts(t, x)) ~ A_frt * ks * Dxx(Ts(t, x)) - (h_avg_f(qlpm) * A_exchange * ((Ts(t, x)) - Tf(t, x))) .- (kins * (r_ins/r0).* (Ts(t, x) .- Tins_f(t)) * A_s_p/ (r_ins - r0)),
        e * ρf * Cpf * Vf * Dt(Tf(t, x)) ~ Af * e * kf * Dxx(Tf(t, x)) - m * Cpf * Dx(Tf(t, x)) + (h_avg_f(qlpm) * A_exchange * ((Ts(t, x)) - Tf(t, x)))
    ]
    bcs1 = [
        Ts(0.0, x) ~ Tamb, # initial
        Tf(0.0, x) ~ Tamb, # initial
        -A_frt * ks * Dx(Ts(t, x_max1)) ~ 0.0, # far right
        -A_frt * ks * Dx(Ts(t, x_min1)) ~ Io * A_chnl_frt_all - ϵ * σ * A_chnl_frt_all * (Ts(t, x_min1)^4 - Tamb^4) - hext * A_chnl_frt_all * (Ts(t, x_min1) - Tamb),  # far left
        -Af * kf * Dx(Tf(t, x_max1)) ~ 0.0, #-ρf * Cpf * V * A_ft * (Tf(t, x_max1) - Tamb), # exiting fluid
        -Af * kf * Dx(Tf(t, x_min1)) ~ m * Cpf * (Tf(t, x_min1) - Tamb) # entering fluid (upstream temperature)
    ]
    # Space and time domain for system 1
    domains1 = [t ∈ Interval(t_min, t_max),
        x ∈ Interval(x_min1, x_max1)]

    # ODE system for system 1
    @named pdesys = PDESystem(eq1, bcs1, domains1, [t, x], [Ts(t, x), Tf(t, x)], p_math)


    # MOL parameters for system 1

    order = 2
    discretization = MOLFiniteDifference([x => dx], t, approx_order=order)

    prob = discretize(pdesys, discretization)

end

sol1 = solve(prob, FBDF(), saveat=2)
begin
    Ts_front_t = sol1.u[(Ts(t, x))][:, 1]
    Tf_front_t = sol1.u[(Tf(t, x))][:, 2]
    Ts_back_t = sol1.u[(Ts(t, x))][:, end-1]
    Tf_back_t = sol1.u[(Tf(t, x))][:, end]
    x_domain = collect(sol1.ivdomain[2])
    plot(xlabel="time [s]")
    plot!(sol1.t, Ts_front_t, label="T_fr_s")
    plot!(sol1.t, Tf_front_t, label="T_fr_f")
    plot!(sol1.t, Ts_back_t, label="T_bck_s")
    plot!(sol1.t, Tf_back_t, label="T_bck_f")
end

x__domain = collect(sol1.ivdomain[2])

begin
    plot(title="Solid Temperature Profile T8")
    plot!(sol1.t,
        sol1.u[Ts(t, x)][:, 4], # around 5 mm in T8
        title="Solid Temperature Profile T8",
        label="Numerical",
        xlabel="Time (s)",
        ylabel="Temperature (K)")
    plot()
    plot!(sol1.t,
        sol1.u[Tf(t, x)][:, end-1], # around 136 mm in T3
        title="Gas Temperature Profile T3",
        label="Numerical",
        xlabel="Time (s)",
        ylabel="Temperature (K)")
    # scatter!(
    #     xd1_data,
    #     y1d1_data,
    #     label="Experimental",
    #     xlabel="Time (s)",
    #     ylabel="Temperature (K)")
end
begin
    # model results for different thermocouples
    sol1.u[Ts(t, x)][:, 40] #T9 model 58mm (external)
    sol1.u[Ts(t, x)][:, 77] #T10 model 107mm (external)
    sol1.u[Tf(t, x)][:, 77] #T11 model 107mm (internal-gas)
    sol1.u[Ts(t, x)][:, 77] #T11 model 107mm (internal-solid)
    sol1.u[Tf(t, x)][:, 40] #T12 model 58mm (internal-gas)
    sol1.u[Ts(t, x)][:, 40] #T12 model 58mm (internal-solid)
end

dec = 100
function decimate!(x, step)
    x = [x[i] for i in 1:step:length(x)]
end

#Exp. data to extract temp.
begin
    #Exp 67 - T3, T8
    Z = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0067_231125_161757.xlsx")["Sheet 1 - Data_FPT0067_231125_1"]["A3:C3932"]
    E67t = float.(Z[:, 1]) #xz_data
    E67Tf = Z[:, 2] .+ 273.15 #y1z_data
    E67T8 = Z[:, 3] .+ 273.15
 #decimate experimental data E67 for T3
    E67t = decimate!(E67t, dec)
    E67Tf = decimate!(E67Tf, dec)
    scatter(E67t, E67Tf)

    #Exp 67 - T9, T10, T11, T12
    Z1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0067_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0067_T9"]["A3:E3932"]
    E67T9 = Z1[:, 2] .+ 273.15
    E67T10 = Z1[:, 3] .+ 273.15
    E67T11 = Z1[:, 4] .+ 273.15
    E67T12 = Z1[:, 5] .+ 273.15
    #decimate experimental data E67 for T8
    E67T8 = decimate!(E67T8, dec)
    scatter(E67t, E67T8)
    #decimate experimental data E67 for T9
    E67T9 = decimate!(E67T9, dec)
    scatter(E67t, E67T9)
    #decimate experimental data E67 for T10
    E67T10 = decimate!(E67T10, dec)
    scatter(E67t, E67T10)
    #decimate experimental data E67 for T11
    E67T11 = decimate!(E67T11, dec)
    scatter(E67t, E67T11)
    #decimate experimental data E67 for T12
    E67T12 = decimate!(E67T12, dec)
    scatter(E67t, E67T12)
end 

begin 
#Exp 68 - T3, T8
    A1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0068_231126_115725.xlsx")["Sheet 1 - Data_FPT0068_231126_1"]["A3:C5365"]
    E68t = float.(A1[:, 1]) #xa1_data
    E68Tf = A1[:, 2] .+ 273.15 #y1a1_data
    E68T8 = A1[:, 3] .+ 273.15
 #decimate experimental data E68 for T3
    E68t = decimate!(E68t, dec)
    E68Tf = decimate!(E68Tf, dec)
    scatter(E68t, E68Tf)
    #decimate experimental data E68 for T8
    E68T8 = decimate!(E68T8, dec)
    scatter(E68t, E68T8)


    #Exp 68 - T9, T10, T11, T12
    A11 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0068_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0068_231126_1"]["A3:E5365"]
    E68T9 = A11[:, 2] .+ 273.15
    E68T10 = A11[:, 3] .+ 273.15
    E68T11 = A11[:, 4] .+ 273.15
    E68T12 = A11[:, 5] .+ 273.15
    #decimate experimental data E68 for T9
    E68T9 = decimate!(E68T9, dec)
    scatter(E68t, E68T9)
    #decimate experimental data E68 for T10
    E68T10 = decimate!(E68T10, dec)
    scatter(E68t, E68T10)
    #decimate experimental data E68 for T11
    E68T11 = decimate!(E68T11, dec)
    scatter(E68t, E68T11)
    #decimate experimental data E68 for T12
    E68T12 = decimate!(E68T12, dec)
    scatter(E68t, E68T12)
end

begin
    #Exp 69 - T3, T8
    B1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0069_231126_140153.xlsx")["Sheet 1 - Data_FPT0069_231126_1"]["A3:C5366"]
    E69t = float.(B1[:, 1])
    E69Tf = B1[:, 2] .+ 273.15
    E69T8 = B1[:, 3] .+ 273.15
#decimate experimental data E69 for T3
    E69t = decimate!(E69t, dec)
    E69Tf = decimate!(E69Tf, dec)
    scatter(E69t, E69Tf)
    #decimate experimental data E69 for T8
    E69T8 = decimate!(E69T8, dec)
    scatter(E69t, E69T8)

    #Exp 69 - T9, T10, T11, T12
    B11 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0069_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0069_231126_1"]["A3:E5366"]
    E69T9 = B11[:, 2] .+ 273.15
    E69T10 = B11[:, 3] .+ 273.15
    E69T11 = B11[:, 4] .+ 273.15
    E69T12 = B11[:, 5] .+ 273.15
    #decimate experimental data E69 for T9
    E69T9 = decimate!(E69T9, dec)
    scatter(E69t, E69T9)
    #decimate experimental data E69 for T10
    E69T10 = decimate!(E69T10, dec)
    scatter(E69t, E69T10)
    #decimate experimental data E69 for T11
    E69T11 = decimate!(E69T11, dec)
    scatter(E69t, E69T11)
    #decimate experimental data E69 for T12
    E69T12 = decimate!(E69T12, dec)
    scatter(E69t, E69T12)
end


begin
    #Exp 70 - T3, T8
    C1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0070_231127_090339.xlsx")["Sheet 1 - Data_FPT0070_231127_0"]["A3:C6705"]
    E70t = float.(C1[:, 1]) #xc1_data 
    E70Tf = C1[:, 2] .+ 273.15 #y1c1_data
    E70T8 = C1[:, 3] .+ 273.15
 #decimate experimental data E70 for T3
    E70t = decimate!(E70t, dec)
    E70Tf = decimate!(E70Tf, dec)
    scatter(E70t, E70Tf)
    #decimate experimental data E70 for T8
    E70T8 = decimate!(E70T8, dec)
    scatter(E70t, E70T8)

    #Exp 70 - T9, T10, T11, T12
    C11 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0070_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0070_231127_0"]["A3:E6705"]
    E70T9 = C11[:, 2] .+ 273.15
    E70T10 = C11[:, 3] .+ 273.15
    E70T11 = C11[:, 4] .+ 273.15
    E70T12 = C11[:, 5] .+ 273.15
    #decimate experimental data E70 for T9
    E70T9 = decimate!(E70T9, dec)
    scatter(E70t, E70T9)
    #decimate experimental data E70 for T10
    E70T10 = decimate!(E70T10, dec)
    scatter(E70t, E70T10)
    #decimate experimental data E70 for T11
    E70T11 = decimate!(E70T11, dec)
    scatter(E70t, E70T11)
    #decimate experimental data E70 for T12
    E70T12 = decimate!(E70T12, dec)
    scatter(E70t, E70T12)
end
begin
    #Exp 71 - T3, T8
    D1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0071_231128_102707.xlsx")["Sheet 1 - Data_FPT0071_231128_1"]["A3:C7087"]
    E71t = float.(D1[:, 1]) #xd1_data  
    E71Tf = D1[:, 2] .+ 273.15 #y1d1_data
    E71T8 = D1[:, 3] .+ 273.15
#decimate experimental data E71 for T3
    E71t = decimate!(E71t, dec)
    E71Tf = decimate!(E71Tf, dec)
    scatter(E71t, E71Tf)
    #decimate experimental data E71 for T8
    E71T8 = decimate!(E71T8, dec)
    scatter(E71t, E71T8)
    #Exp 71 - T9, T10, T11, T12
    D11 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0071_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0071_231128_1"]["A3:E7087"]
    E71T9 = D11[:, 2] .+ 273.15
    E71T10 = D11[:, 3] .+ 273.15
    E71T11 = D11[:, 4] .+ 273.15
    E71T12 = D11[:, 5] .+ 273.15
    #decimate experimental data E71 for T9
    E71T9 = decimate!(E71T9, dec)
    scatter(E71t, E71T9)
    #decimate experimental data E71 for T10
    E71T10 = decimate!(E71T10, dec)
    scatter(E71t, E71T10)
    #decimate experimental data E71 for T11
    E71T11 = decimate!(E71T11, dec)
    scatter(E71t, E71T11)
    #decimate experimental data E71 for T12
    E71T12 = decimate!(E71T12, dec)
    scatter(E71t, E71T12)
end

begin
#Exp 72 - T3, T8
    E1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0072_231129_104140.xlsx")["Sheet 1 - Data_FPT0072_231129_1"]["A3:C3217"]
    E72t = float.(E1[:, 1]) #xe1_data
    E72Tf = E1[:, 2] .+ 273.15 #y1e1_data
    E72T8 = E1[:, 3] .+ 273.15
 #decimate experimental data E72 for T3
    E72t = decimate!(E72t , dec)
    E72Tf = decimate!(E72Tf, dec)
    scatter(E72t, E72Tf)
    #decimate experimental data E72 for T8
    E72T8 = decimate!(E72T8, dec)
    scatter(E72t, E72T8)

    #Exp 72 - T9, T10, T11, T12
    E11 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0072_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0072_231129_1"]["A3:E3217"]
    E72T9 = E11[:, 2] .+ 273.15
    E72T10 = E11[:, 3] .+ 273.15
    E72T11 = E11[:, 4] .+ 273.15
    E72T12 = E11[:, 5] .+ 273.15
    #decimate experimental data E72 for T9
    E72T9 = decimate!(E72T9, dec)
    scatter(E72t, E72T9)
    #decimate experimental data E72 for T10
    E72T10 = decimate!(E72T10, dec)
    scatter(E72t, E72T10)
    #decimate experimental data E72 for T11
    E72T11 = decimate!(E72T11, dec)
    scatter(E72t, E72T11)
    #decimate experimental data E72 for T12
    E72T12 = decimate!(E72T12, dec)
    scatter(E72t, E72T12)
end

begin
#Exp 73 - T3, T8
    F1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0073_231129_132744.xlsx")["Sheet 1 - Data_FPT0073_231129_1"]["A3:C4575"]
    E73t = float.(F1[:, 1])#xf1_data
    E73Tf = F1[:, 2] .+ 273.15 #y1f1_data
    E73T8= F1[:, 3] .+ 273.15
 #decimate experimental data E73 for T3
    E73t = decimate!(E73t, dec)
    E73Tf = decimate!(E73Tf, dec)
    scatter(E73t, E73Tf)
    #decimate experimental data E73 for T8
    E73T8 = decimate!(E73T8, dec)
    scatter(E73t, E73T8)

    #Exp 73 - T9, T10, T11, T12
    F11 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0073_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0073_231129_1"]["A3:E4575"]
    E73T9 = F11[:, 2] .+ 273.15
    E73T10 = F11[:, 3] .+ 273.15
    E73T11 = F11[:, 4] .+ 273.15
    E73T12 = F11[:, 5] .+ 273.15
    #decimate experimental data E73 for T9
    E73T9 = decimate!(E73T9, dec)
    scatter(E73t, E73T9)
    #decimate experimental data E73 for T10
    E73T10 = decimate!(E73T10, dec)
    scatter(E73t, E73T10)
    #decimate experimental data E73 for T11
    E73T11 = decimate!(E73T11, dec)
    scatter(E73t, E73T11)
    #decimate experimental data E73 for T12
    E73T12 = decimate!(E73T12, dec)
    scatter(E73t, E73T12)

end

begin  
#Exp 74 - T3, T8
    G = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0074_231130_123228.xlsx")["Sheet 1 - Data_FPT0074_231130_1"]["A3:C6018"]
    E74t = float.(G[:, 1]) #xg_data 
    E74Tf = G[:, 2] .+ 273.15 #y1g_data 
    E74T8 = G[:, 3] .+ 273.15
 #decimate experimental data E74 for T3
    E74t = decimate!(E74t, dec)
    E74Tf = decimate!(E74Tf, dec)
    scatter(E74t, E74Tf)
    #decimate experimental data E74 for T8
    E74T8 = decimate!(E74T8, dec)
    scatter(E74t, E74T8)

    #Exp 74 - T9, T10, T11, T12
    G1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0074_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0074_231130_1"]["A3:E6018"]
    E74T9 = G1[:, 2] .+ 273.15
    E74T10 = G1[:, 3] .+ 273.15
    E74T11 = G1[:, 4] .+ 273.15
    E74T12 = G1[:, 5] .+ 273.15
    #decimate experimental data E74 for T9
    E74T9 = decimate!(E74T9, dec)
    scatter(E74t, E74T9)
    #decimate experimental data E74 for T10
    E74T10 = decimate!(E74T10, dec)
    scatter(E74t, E74T10)
    #decimate experimental data E74 for T11
    E74T11 = decimate!(E74T11, dec)
    scatter(E74t, E74T11)
    #decimate experimental data E74 for T12
    E74T12 = decimate!(E74T12, dec)
    scatter(E74t, E74T12)
end

begin
    #Exp 75 - T3, T8
    H = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0075_231201_162138.xlsx")["Sheet 1 - Data_FPT0075_231201_1"]["A3:C6354"]
    E75t= float.(H[:, 1]) #xh_data
    E75Tf = H[:, 2] .+ 273.15 #y1h_data
    E75T8 = H[:, 3] .+ 273.15
 #decimate experimental data E75 for T3
    E75t= decimate!(E75t, dec)
    E75Tf = decimate!(E75Tf, dec)
    scatter(E75t, E75Tf)
    #decimate experimental data E75 for T8
    E75T8 = decimate!(E75T8, dec)
    scatter(E75t, E75T8)


    #Exp 75 - T9, T10, T11, T12
    H1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0075_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0075_231201_1"]["A3:E6354"]
    E75T9 = H1[:, 2] .+ 273.15
    E75T10 = H1[:, 3] .+ 273.15
    E75T11 = H1[:, 4] .+ 273.15
    E75T12 = H1[:, 5] .+ 273.15
    #decimate experimental data E75 for T9
    E75T9 = decimate!(E75T9, dec)
    scatter(E75t, E75T9)
    #decimate experimental data E75 for T10
    E75T10 = decimate!(E75T10, dec)
    scatter(E75t, E75T10)
    #decimate experimental data E75 for T11
    E75T11 = decimate!(E75T11, dec)
    scatter(E75t, E75T11)
    #decimate experimental data E75 for T12
    E75T12 = decimate!(E75T12, dec)
    scatter(E75t, E75T12)

end

begin
#Exp 76 - T3, T8
    I = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0076_231203_120521.xlsx")["Sheet 1 - Data_FPT0076_231203_1"]["A3:C7147"]
    E76t = float.(I[:, 1]) #xi_data
    E76Tf = I[:, 2] .+ 273.15 #y1i_data
    E76T8 = I[:, 3] .+ 273.15
 #decimate experimental data E76 for T3
    E76t = decimate!(E76t, dec)
    E76Tf = decimate!(E76Tf, dec)
    scatter(E76t, E76Tf)
    #decimate experimental data E76 for T8
    E76T8 = decimate!(E76T8, dec)
    scatter(E76t, E76T8)

    #Exp 76 - T9, T10, T11, T12
    I1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0076_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0076_231203_1"]["A3:E7147"]
    E76T9 = I1[:, 2] .+ 273.15
    E76T10 = I1[:, 3] .+ 273.15
    E76T11 = I1[:, 4] .+ 273.15
    E76T12 = I1[:, 5] .+ 273.15
    #decimate experimental data E76 for T9
    E76T9 = decimate!(E76T9, dec)
    scatter(E76t, E76T9)
    #decimate experimental data E76 for T10
    E76T10 = decimate!(E76T10, dec)
    scatter(E76t, E76T10)
    #decimate experimental data E76 for T11
    E76T11 = decimate!(E76T11, dec)
    scatter(E76t, E76T11)
    #decimate experimental data E76 for T12
    E76T12 = decimate!(E76T12, dec)
    scatter(E76t, E76T12)

end

begin
    #Exp 77 - T3, T8
    J = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0077_231203_161315.xlsx")["Sheet 1 - Data_FPT0077_231203_1"]["A3:C3044"]
    E77t = float.(J[:, 1]) #xj_data
    E77Tf = J[:, 2] .+ 273.15 #y1j_data
    E77T8 = J[:, 3] .+ 273.15
 #decimate experimental data E77 for T3
    E77t = decimate!(E77t, dec)
    E77Tf = decimate!(E77Tf, dec)
    scatter(E77t, E77Tf)
    #decimate experimental data E77 for T8
    E77T8 = decimate!(E77T8, dec)
    scatter(E77t, E77T8)

    #Exp 77 - T9, T10, T11, T12
    J1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0077_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0077_231203_1"]["A3:E3044"]
    E77T9 = J1[:, 2] .+ 273.15
    E77T10 = J1[:, 3] .+ 273.15
    E77T11 = J1[:, 4] .+ 273.15
    E77T12 = J1[:, 5] .+ 273.15
    #decimate experimental data E77 for T9
    E77T9 = decimate!(E77T9, dec)
    scatter(E77t, E77T9)
    #decimate experimental data E77 for T10
    E77T10 = decimate!(E77T10, dec)
    scatter(E77t, E77T10)
    #decimate experimental data E77 for T11
    E77T11 = decimate!(E77T11, dec)
    scatter(E77t, E77T11)
    #decimate experimental data E77 for T12
    E77T12 = decimate!(E77T12, dec)
    scatter(E77t, E77T12)
end

begin
#Exp 78 - T3, T8
    K = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0078_231204_132252.xlsx")["Sheet 1 - Data_FPT0078_231204_1"]["A3:C5384"]
    E78t = float.(K[:, 1]) #xk_data
    E78Tf = K[:, 2] .+ 273.15 #y1k_data
    E78T8 = K[:, 3] .+ 273.15
 #decimate experimental data E78 for T3
    E78t = decimate!(E78t, dec)
    E78Tf = decimate!(E78Tf, dec)
    scatter(E78t, E78Tf)
    #decimate experimental data E78 for T8
    E78T8 = decimate!(E78T8, dec)
    scatter(E78t, E78T8)

    #Exp 78 - T9, T10, T11, T12
    K1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0078_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0078_231204_1"]["A3:E5384"]
    E78T9 = K1[:, 2] .+ 273.15
    E78T10 = K1[:, 3] .+ 273.15
    E78T11 = K1[:, 4] .+ 273.15
    E78T12 = K1[:, 5] .+ 273.15
    #decimate experimental data E78 for T9
    E78T9 = decimate!(E78T9, dec)
    scatter(E78t, E78T9)
    #decimate experimental data E78 for T10
    E78T10 = decimate!(E78T10, dec)
    scatter(E78t, E78T10)
    #decimate experimental data E78 for T11
    E78T11 = decimate!(E78T11, dec)
    scatter(E78t, E78T11)
    #decimate experimental data E78 for T12
    E78T12 = decimate!(E78T12, dec)
    scatter(E78t, E78T12)
end

begin
#Exp 79 - T3, T8
    L0 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0079_231204_172244.xlsx")["Sheet 1 - Data_FPT0079_231204_1"]["A3:C5233"]
    E79t = float.(L0[:, 1]) #xl_data 
    E79Tf = L0[:, 2] .+ 273.15 #y1l_data
    E79T8 = L0[:, 3] .+ 273.15  
 #decimate experimental data E79 for T3
    E79t= decimate!(E79t, dec)
    E79Tf = decimate!(E79Tf, dec)
    scatter(E79t, E79Tf)
    #decimate experimental data E79 for T8
    E79T8 = decimate!(E79T8, dec)
    scatter(E79t, E79T8)

    #Exp 79 - T9, T10, T11, T12
    L1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0079_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0079_231204_1"]["A3:E5233"]
    E79T9 = L1[:, 2] .+ 273.15
    E79T10 = L1[:, 3] .+ 273.15
    E79T11 = L1[:, 4] .+ 273.15
    E79T12 = L1[:, 5] .+ 273.15
    #decimate experimental data E79 for T9
    E79T9 = decimate!(E79T9, dec)
    scatter(E79t, E79T9)
    #decimate experimental data E79 for T10
    E79T10 = decimate!(E79T10, dec)
    scatter(E79t, E79T10)
    #decimate experimental data E79 for T11
    E79T11 = decimate!(E79T11, dec)
    scatter(E79t, E79T11)
    #decimate experimental data E79 for T12
    E79T12 = decimate!(E79T12, dec)
    scatter(E79t, E79T12)
end

begin
    #Exp 80 - T3, T8
    M = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0080_231205_095122.xlsx")["Sheet 1 - Data_FPT0080_231205_0"]["A3:C5814"]
    E80t = float.(M[:, 1]) # xm_data
    E80Tf = M[:, 2] .+ 273.15 #y1m_data
    E80T8 = M[:, 3] .+ 273.15   
 #decimate experimental data E80 for T3
    E80t = decimate!(E80t, dec)
    E80Tf = decimate!(E80Tf, dec)
    scatter(E80t, E80Tf)
 
  #decimate experimental data E80 for T8
    E80T8 = decimate!(E80T8, dec)
    scatter(E80t, E80T8)
end
begin
    #Exp 80 - T9, T10, T11, T12
    M1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0080_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0080_231205_0"]["A3:E5814"]
    E80T9 = M1[:, 2] .+ 273.15
    E80T10 = M1[:, 3] .+ 273.15
    E80T11 = M1[:, 4] .+ 273.15
    E80T12 = M1[:, 5] .+ 273.15
    #decimate experimental data E80 for T9
    E80T9 = decimate!(E80T9, dec)
    scatter(E80t, E80T9)
    #decimate experimental data E80 for T10
    E80T10 = decimate!(E80T10, dec)
    scatter(E80t, E80T10)
    #decimate experimental data E80 for T11
    E80T11 = decimate!(E80T11, dec)
    scatter(E80t, E80T11)
    #decimate experimental data E80 for T12
    E80T12 = decimate!(E80T12, dec)
    scatter(E80t, E80T12)
end

begin
    #Exp 81 - T3, T8
    N = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0081_231205_135354.xlsx")["Sheet 1 - Data_FPT0081_231205_1"]["A3:C5989"]
    E81t= float.(N[:, 1]) #xn_data 
    E81Tf = N[:, 2] .+ 273.15 #y1n_data 
    E81T8 = N[:, 3] .+ 273.15   
 #decimate experimental data E81 for T3
    E81t = decimate!(E81t, dec)
    E81Tf = decimate!(E81Tf, dec)
    scatter(E81t, E81Tf)
    #decimate experimental data E81 for T8
    E81T8 = decimate!(E81T8, dec)
    scatter(E81t, E81T8)

    #Exp 81 - T9, T10, T11, T12
    N1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0081_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0081_231205_1"]["A3:E5989"]
    E81T9 = N1[:, 2] .+ 273.15
    E81T10 = N1[:, 3] .+ 273.15
    E81T11 = N1[:, 4] .+ 273.15
    E81T12 = N1[:, 5] .+ 273.15
    #decimate experimental data E81 for T9
    E81T9 = decimate!(E81T9, dec)
    scatter(E81t, E81T9)
    #decimate experimental data E81 for T10
    E81T10 = decimate!(E81T10, dec)
    scatter(E81t, E81T10)
    #decimate experimental data E81 for T11
    E81T11 = decimate!(E81T11, dec)
    scatter(E81t, E81T11)
    #decimate experimental data E81 for T12
    E81T12 = decimate!(E81T12, dec)
    scatter(E81t, E81T12)

end


#measurements and conditions#Defining simulation conditions
begin
    condition_E67 = Dict(Io => 469680.0, qlpm => 15.27)
    condition_E68 = Dict(Io => 469680.0, qlpm => 12.50)
    condition_E69 = Dict(Io => 469680.0, qlpm => 10.50)
    condition_E70 = Dict(Io => 469680.0, qlpm => 9.10)
    condition_E71 = Dict(Io => 469680.0, qlpm => 7.12)
    condition_E72 = Dict(Io => 370880.0, qlpm => 18.34)
    condition_E73 = Dict(Io => 370880.0, qlpm => 13.16)
    condition_E74 = Dict(Io => 370880.0, qlpm => 9.03)
    condition_E75 = Dict(Io => 370880.0, qlpm => 6.95)
    condition_E76 = Dict(Io => 370880.0, qlpm => 4.53)
    condition_E77 = Dict(Io => 248320.0, qlpm => 13.85)
    condition_E78 = Dict(Io => 248320.0, qlpm => 10.02)
    condition_E79 = Dict(Io => 248320.0, qlpm => 8.04)
    condition_E80 = Dict(Io => 248320.0, qlpm => 6.62)
    condition_E81 = Dict(Io => 248320.0, qlpm => 4.53)


    simulation_conditions = Dict("E67" => condition_E67, "E68" => condition_E68,
        "E69" => condition_E69, "E70" => condition_E70,
        "E71" => condition_E71, "E72" => condition_E72,
        "E73" => condition_E73, "E74" => condition_E74,
        "E75" => condition_E75, "E76" => condition_E76,
        "E77" => condition_E77, "E78" => condition_E78,
        "E79" => condition_E79, "E80" => condition_E80,
        "E81" => condition_E81)
    #Defining measurement data
    # measurements = DataFrame(
    #     simulation_id=repeat(["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"], inner=6, outer=1),
    #     obs_id=repeat(["_T8", "_T9", "_T10", "_T11", "_T12", "_T3"], inner=1, outer=15),
    #     time=repeat([3929, 5362, 5363, 6702, 7084, 3214, 4572, 6015, 6351, 7144, 3041, 5381, 5230, 5811, 5986], inner=6, outer=1),
    #     temperatures=[965.407, 975.144, 825.592, 880.867, 1004.165, 763.859,
    #         1031.574, 1023.115, 850.099, 898.754, 1050.691, 773.207,
    #         1070.803, 1045.898, 852.837, 896.788, 1072.727, 769.76,
    #         1167.978, 1093.849, 871.496, 912.173, 1120.42, 779.53,
    #         1210.945, 1095.322, 847.476, 882.823, 1120.417, 753.56,
    #         742.125, 778.592, 684.246, 736.626, 807.125, 652.955,
    #         844.257, 870.26, 747.444, 791.958, 898.081, 694.626,
    #         962.113, 938.106, 767.803, 804.082, 965.702, 697.672,
    #         1015.081, 954.214, 757.393, 788.678, 979.795, 681.066,
    #         1069.567, 947.372, 726.308, 751.159, 970.498, 647.019,
    #         569.248, 604.984, 543.727, 574.299, 627.16, 525.356,
    #         634.731, 664.296, 583.704, 612.092, 686.936, 554.455,
    #         677.817, 694.156, 595.766, 622.325, 716.314, 560.033,
    #         711.537, 713.686, 601.77, 626.485, 735.984, 561.254,
    #         763.299, 729.461, 597.766, 618.975, 751.15, 550.499])
    # measurements = DataFrame( #T3 only
    #             simulation_id=repeat(["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"], inner=1, outer=1),
    #             obs_id=repeat(["_T3"], inner=1, outer=15),
    #             time=repeat([3929, 5362, 5363, 6702, 7084, 3214, 4572, 6015, 6351, 7144, 3041, 5381, 5230, 5811, 5986], inner=1, outer=1),
    #             temperatures=[763.859, 773.207, 769.76, 779.53, 753.56, 652.955, 694.626, 697.672, 681.066,
    #                 647.019, 525.356, 554.455, 560.033, 561.254, 550.499])
    measurements = DataFrame(
        simulation_id = repeat(["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"], inner=6, outer=1),
        obs_id=repeat(["_T8", "_T9", "_T10", "_T11", "_T12", "_T3"], inner=1, outer=15),
        time= repeat([E67t, E68t, E69t, E70t, E71t, E72t, E73t, E74t, E75t, E76t, E77t, E78t, E79t, E80t, E81t], inner=6, outer=1),
            temperatures= [E67T8, E67T9, E67T10, E67T11, E67T12, E67Tf,
            E68T8, E68T9, E68T10, E68T11, E68T12, E68Tf,
            E69T8, E69T9, E69T10, E69T11, E69T12, E69Tf,
            E70T8, E70T9, E70T10, E70T11, E70T12, E70Tf,
            E71T8, E71T9, E71T10, E71T11, E71T12, E71Tf,
            E72T8, E72T9, E72T10, E72T11, E72T12, E72Tf,
            E73T8, E73T9, E73T10, E73T11, E73T12, E73Tf,
            E74T8, E74T9, E74T10, E74T11, E74T12, E74Tf,
            E75T8, E75T9, E75T10, E75T11, E75T12, E75Tf,
            E76T8, E76T9, E76T10, E76T11, E76T12, E76Tf,
            E77T8, E77T9, E77T10, E77T11, E77T12, E77Tf,
            E78T8, E78T9, E78T10, E78T11, E78T12, E78Tf,
            E79T8, E79T9, E79T10, E79T11, E79T12, E79Tf,
            E80T8, E80T9, E80T10, E80T11, E80T12, E80Tf,
            E81T8, E81T9, E81T10, E81T11, E81T12, E81Tf])   
    end



#Optimization using NLOpt
rmp = ModelingToolkit.varmap_to_vars([Io => 469680.0, A => 250.1, B => 0.5, C=> 5., qlpm => 7.12], parameters(pdesys))
function NLmodeloptim(tvalues, rmp)

    #p = [hlocal => p_vary[1]]
    modeloptim = remake(prob, p=rmp, tspan=(0.0, tvalues[end]))
    modeloptim_sol = solve(modeloptim, FBDF(), saveat=tvalues)#, reltol=1e-12, abstol=1e-12)
    #time = modelfit_sol.t
    tempT8_op = float.(modeloptim_sol.u[Ts(t, x)][:, 8])
    tempT9_op = float.(modeloptim_sol.u[Ts(t, x)][:, 42])
    tempT10_op = float.(modeloptim_sol.u[Ts(t, x)][:, 78])
    tempT3_op = float.(modeloptim_sol.u[Tf(t, x)][:, end-1])
    # T12_modelmean = (modeloptim_sol.u[Tf(t, x)][end, 20] .+ modeloptim_sol.u[Ts(t, x)][end, 20]) ./ 2
    # T11_modelmean = (modeloptim_sol.u[Tf(t, x)][end, 59] .+ modeloptim_sol.u[Ts(t, x)][end, 59]) ./ 2
    tempT11_op = float.(modeloptim_sol.u[Ts(t, x)][:, 78])
    tempT12_op = float.(modeloptim_sol.u[Ts(t, x)][:, 42])
    return append!(tempT8_op, tempT9_op, tempT10_op, tempT3_op, tempT12_op, tempT11_op)
    #return ([tempT3_op])
end

function remakeAysha(pguess_l, cond, time_opt)
        #p_math_vec = copy(cond)
        # Initialize a new p_math vector
        rmp = Vector{Pair{Symbol, Float64}}(undef, length(pguess_l) + length(cond))
        # Fill the new_p_math with values from pguess_l
        for i in 1:length(pguess_l)
            rmp[i] = Symbol(p_math[i][1]) => pguess_l[i]
        end
        # Fill the new_p_math with values from cond
        for (i, (key, value)) in enumerate(cond)
            rmp[length(pguess_l) + i] = Symbol(key) => value
        end# pguess is the initial guess for the optimization
        #println(rmp)
        temp_T = NLmodeloptim(time_opt, rmp)
        return temp_T
end
#the first entry in the simulation_conditions
pguess = p_opt

function lossAll(pguess_l, _)

    #to place the conditions loop

    sim_key = collect(keys(simulation_conditions))
    lossr = zeros(length(sim_key))
   #Threads.@threads 
   for it in 1:length(sim_key)
    # Retrieve from measurements the experimental data for the current simulation condition
        sm = sim_key[it]
        #println(sm)
        cond = simulation_conditions[sm]
        expdata = reduce(vcat, measurements[measurements.simulation_id.==sm, :temperatures])
        time_opt = measurements[measurements.simulation_id.==sm, :time][1]

        # Ensure time_opt are vectors of floats
        
        #run selected simulation and get the steady temperature values
        #print(sm)
        temp_T = remakeAysha(pguess_l, cond, time_opt)
        
        temp_error = (temp_T .- expdata) .^ 2
        lossr[it] = sqrt(sum(temp_error))
    end
    return sum(lossr) #MSE
end
p0 = [x[2] for x in p_opt]
optf = OptimizationFunction(lossAll, Optimization.AutoForwardDiff())
#p_opt = [aCp => 1., hfa => 8., hfn =>0.66, aIo => 1.] 
lb = [70., 0.2, 5.]
ub = [5000., 2., 20.]
 
#pguess_opt = ModelingToolkit.varmap_to_vars([Io => 456000, h_average => 14., qlpm => 7.12], parameters(pdesys))
initialerror = (lossAll(p0, []))
println(initialerror)

optprob = Optimization.OptimizationProblem(optf, p0, [], lb=lb, ub=ub)

optsol = solve(optprob, NLopt.GN_MLSL_LDS(), local_method=NLopt.LN_NELDERMEAD(), maxtime=100, local_maxiters=100000)

println(optsol.retcode)
pnew = optsol.u
println(pnew)
res_error = lossAll(pnew, [])


display(res_error)

# pnew = [0.01, 0.6, 5.]
begin
    T_steady = DataFrame(sim_id=[], T_mod=[], T_exp=[])
    sim_key = collect(keys(simulation_conditions))
    lossr = zeros(length(sim_key))
   #Threads.@threads 
   for it in 1:length(sim_key)
    # Retrieve from measurements the experimental data for the current simulation condition
        sm = sim_key[it]
        #println(sm)
        cond = simulation_conditions[sm]
        expdata = measurements[measurements.simulation_id.==sm, :temperatures]
        time_opt = measurements[measurements.simulation_id.==sm, :time][1]
        
        #run selected simulation and get the steady temperature values
        temp_T = remakeAysha(pnew, cond, time_opt)
        temp_T = reshape(temp_T, length(expdata[1]), length(expdata))
        temp_T = [temp_T[:, i] for i in 1:length(expdata)]
        push!(T_steady, (sm, temp_T, expdata))
   end

end

begin #TI-8
    color_model = :blue
    color_exp = :orange
    ordered_sim_conditions = ["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"]
    order_indices = [findfirst(x -> x == condition, T_steady.sim_id) for condition in ordered_sim_conditions]
    plot1 = bar(
        [1, 2], [T_steady[order_indices[1], :T_mod][end][end], T_steady[order_indices[1], :T_exp][end][end]],
        title="TI-8 Steady State Temperature", ylabel="Temperature (K)", xlabel="Experimental Runs",
        bar_width=0.4, xticks=(1:2:30, ordered_sim_conditions),
       label="T_model", color=[color_model, color_exp], ylimit=(0, 1250))
 
    # Loop through the remaining data to add the bars without labels
    for i in 2:15
        bar!(plot1, [2i-1, 2i], [T_steady[order_indices[i], :T_mod][end][end], T_steady[order_indices[i], :T_exp][end][end]], bar_width=0.4, label=["" ""], color=[color_model, color_exp])
    end
    # Manually add a legend entry for T_exp and T_model
    plot1 = bar!(plot1, [0], [0], label="T_exp", color=color_exp)
    xlims!(plot1, (0, 31))
end

begin

    ordered_sim_conditions = ["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"]
    order_indices = [findfirst(x -> x == condition, T_steady.sim_id) for condition in ordered_sim_conditions]

    # Extract the data for the scatter plot
    model_temps = [T_steady[order_indices[i], :T_mod][end][end] for i in 1:15]
    exp_temps = [T_steady[order_indices[i], :T_exp][end][end] for i in 1:15]

    # Create the scatter plot
    plot1 = scatter(
        exp_temps, model_temps,
        title="TI-3 Steady State Temperature",
        xlabel="Experimental Temperature (K)",
        ylabel="Model Temperature (K)",
        label="Temperature data points",
        color=color_model
    )
    # Add the best fit line to the plot
plot!([500, 1250], [500, 1250], label="Ideal case (y=x)", color=:red)
end

new_78t = LinRange(sol1.t[1], sol1.t[end], 54) #E78
begin
plot1 = plot(
        new_78t,
        T_steady[1, :T_mod], # around 136 mm in T3
        title="Gas Temperature Profile T3 in E78",
        label="model",
        xlabel="Time (s)",
        ylabel="Temperature (K)", xlimit = (0, 4050))
        scatter!(E78t, E78Tf, label="experiment")
        display(plot1)
end 

new_69t = LinRange(sol1.t[1], sol1.t[end], 54)
begin
plot1 = plot(
        new_69t,
        T_steady[2, :T_mod], # around 136 mm in T3
        title="Gas Temperature Profile T3 in E69",
        label="model",
        xlabel="Time (s)",
        ylabel="Temperature (K)", xlimit = (0, 4050))
        scatter!(E69t, E69Tf, label="experiment")
        display(plot1)
end 

new_79t = LinRange(sol1.t[1], sol1.t[end], 53)
begin
    plot1 = plot(
            new_79t,
            T_steady[3, :T_mod], # around 136 mm in T3
            title="Gas Temperature Profile T3 in E79",
            label="model",
            xlabel="Time (s)",
            ylabel="Temperature (K)", xlimit = (0, 10050))
            scatter!(E79t, E79Tf, label="experiment")
            display(plot1)
end 

new_73t = LinRange(sol1.t[1], sol1.t[end], 46)
begin
    plot1 = plot(
            new_73t,
            T_steady[4, :T_mod], # around 136 mm in T3
            title="Gas Temperature Profile T3 in E73",
            label="model",
            xlabel="Time (s)",
            ylabel="Temperature (K)", xlimit = (0, 4050))
            scatter!(E73t, E73Tf, label="experiment")
            display(plot1)
end 
    
new_81t = LinRange(sol1.t[1], sol1.t[end], 60)
begin
    plot1 = plot(
            new_81t,
            T_steady[5, :T_mod], # around 136 mm in T3
            title="Gas Temperature Profile T3 in E81",
            label="model",
            xlabel="Time (s)",
            ylabel="Temperature (K)", xlimit = (0, 4050))
            scatter!(E81t, E81Tf, label="experiment")
            display(plot1)
end 

new_72t = LinRange(sol1.t[1], sol1.t[end], 33)
begin
    plot1 = plot(
            new_72t,
            T_steady[6, :T_mod], # around 136 mm in T3
            title="Gas Temperature Profile T3 in E72",
            label="model",
            xlabel="Time (s)",
            ylabel="Temperature (K)", xlimit = (0, 3150))
            scatter!(E72t, E72Tf, label="experiment")
            display(plot1)
end 

new_76t = LinRange(sol1.t[1], sol1.t[end], 72)
begin
    plot1 = plot(
            new_76t,
            T_steady[7, :T_mod], # around 136 mm in T3
            title="Gas Temperature Profile T3 in E76",
            label="model",
            xlabel="Time (s)",
            ylabel="Temperature (K)", xlimit = (0, 3600))
            scatter!(E76t, E76Tf, label="experiment")
            display(plot1)
end 

new_77t = LinRange(sol1.t[1], sol1.t[end], 31)
begin
    plot1 = plot(
            new_77t,
            T_steady[8, :T_mod], # around 136 mm in T3
            title="Gas Temperature Profile T3 in E77",
            label="model",
            xlabel="Time (s)",
            ylabel="Temperature (K)", xlimit = (0, 3050))
            scatter!(E77t, E77Tf, label="experiment")
            display(plot1)
end 

new_71t = LinRange(sol1.t[1], sol1.t[end], 71)
begin
    plot1 = plot(
            new_71t,
            T_steady[9, :T_mod], # around 136 mm in T3
            title="Gas Temperature Profile T3 in E71",
            label="model",
            xlabel="Time (s)",
            ylabel="Temperature (K)", xlimit = (0, 3050))
            scatter!(E71t, E71Tf, label="experiment")
            display(plot1)
end 

new_75t = LinRange(sol1.t[1], sol1.t[end], 64)
begin
    plot1 = plot(
            new_75t,
            T_steady[10, :T_mod], # around 136 mm in T3
            title="Gas Temperature Profile T3 in E75",
            label="model",
            xlabel="Time (s)",
            ylabel="Temperature (K)", xlimit = (0, 3050))
            scatter!(E75t, E75Tf, label="experiment")
            display(plot1)
end 

new_67t = LinRange(sol1.t[1], sol1.t[end], 40)
begin
    plot1 = plot(
            new_67t,
            T_steady[11, :T_mod], # around 136 mm in T3
            title="Gas Temperature Profile T3 in E67",
            label="model",
            xlabel="Time (s)",
            ylabel="Temperature (K)", xlimit = (0, 3050))
            scatter!(E67t, E67Tf, label="experiment")
            display(plot1)
end

new_68t = LinRange(sol1.t[1], sol1.t[end], 54)
begin
    plot1 = plot(
            new_68t,
            T_steady[12, :T_mod], # around 136 mm in T3
            title="Gas Temperature Profile T3 in E68",
            label="model",
            xlabel="Time (s)",
            ylabel="Temperature (K)", xlimit = (0, 3050))
            scatter!(E68t, E68Tf, label="experiment")
            display(plot1)
end 

new_70t = LinRange(sol1.t[1], sol1.t[end], 68)
begin
    plot1 = plot(
            new_70t,
            T_steady[13, :T_mod], # around 136 mm in T3
            title="Gas Temperature Profile T3 in E70",
            label="model",
            xlabel="Time (s)",
            ylabel="Temperature (K)", xlimit = (0, 3050))
            scatter!(E70t, E70Tf, label="experiment")
            display(plot1)
end 

plotly()

new_80t = LinRange(sol1.t[1], sol1.t[end], 59)
begin
    plot1 = plot(
            new_80t,
            T_steady[14, :T_mod], # around 136 mm in T3
            title="Gas Temperature Profile T3 in E80",
            label="model",
            xlabel="Time (s)",
            ylabel="Temperature (K)", xlimit = (0,3050))
            scatter!(E80t, E80Tf, label="experiment")
            display(plot1)
end 

new_74t = LinRange(sol1.t[1], sol1.t[end], 61)
begin
    plot1 = plot(
            new_74t,
            T_steady[15, :T_mod], # around 136 mm in T3
            title="Gas Temperature Profile T3 in E74",
            label="model",
            xlabel="Time (s)",
            ylabel="Temperature (K)", xlimit = (0, 3050))
            scatter!(E74t, E74Tf, label="experiment")
            display(plot1)
end 


begin #TI-9
    color_model = :blue
    color_exp = :orange
    ordered_sim_conditions = ["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"]
    order_indices = [findfirst(x -> x == condition, T_steady.sim_id) for condition in ordered_sim_conditions]
    plot1 = bar(
        [1, 2], [T_steady[order_indices[1], :T_mod][end][2], T_steady[order_indices[1], :T_exp][end][2]],
        title="TI-9 Steady State Temperature", ylabel="Temperature (K)", xlabel="Experimental Runs",
        bar_width=0.4, xticks=(1:2:30, ordered_sim_conditions),
       label="T_model", color=[color_model, color_exp],  ylimit=(0, 850))
 
    # Loop through the remaining data to add the bars without labels
    for i in 2:15
        bar!(plot1, [2i-1, 2i], [T_steady[order_indices[i], :T_mod][2], T_steady[order_indices[i], :T_exp][2]], bar_width=0.4, label=["" ""], color=[color_model, color_exp])
    end
    # Manually add a legend entry for T_exp and T_model
    plot1 = bar!(plot1, [0], [0], label="T_exp", color=color_exp)
    xlims!(plot1, (0, 31))
end

    # plot!(
    #     psol.t,
    #     psol.u[Tf(t, x)][:, 4], # around 5 mm in T8
    #     label="gas 5mm")

    # plot3 = plot(
    #     psol.t,
    #     psol.u[Ts(t, x)][:, 40], # around 55 mm in T9
    #     title="Solid Temperature Profile T9",
    #     label="model",
    #     xlabel="Time (s)",
    #     ylabel="Temperature (K)")
    begin #TI-10
        color_model = :blue
        color_exp = :orange
        ordered_sim_conditions = ["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"]
        order_indices = [findfirst(x -> x == condition, T_steady.sim_id) for condition in ordered_sim_conditions]
        plot1 = bar(
            [1, 2], [T_steady[order_indices[1], :T_mod][3], T_steady[order_indices[1], :T_exp][3]],
            title="TI-10 Steady State Temperature", ylabel="Temperature (K)", xlabel="Experimental Runs",
            bar_width=0.4, xticks=(1:2:30, ordered_sim_conditions),
           label="T_model", color=[color_model, color_exp], ylimit=(0, 850))
     
        # Loop through the remaining data to add the bars without labels
        for i in 2:15
            bar!(plot1, [2i-1, 2i], [T_steady[order_indices[i], :T_mod][3], T_steady[order_indices[i], :T_exp][3]], bar_width=0.4, label=["" ""], color=[color_model, color_exp])
        end
        # Manually add a legend entry for T_exp and T_model
        plot1 = bar!(plot1, [0], [0], label="T_exp", color=color_exp)
        xlims!(plot1, (0, 31))
    end
    # plot!(
    #     psol.t,
    #     psol.u[Tf(t, x)][:, 40],
    #     label="gas 55mm")
   
# plot4 = plot(
#         psol.t,
#         psol.u[Ts(t, x)][:, 77], # around 106 mm in T10
#         title="Solid Temperature Profile T10",
#         label="model",
#         xlabel="Time (s)",
#         ylabel="Temperature (K)")
begin #TI-11
    color_model = :blue
    color_exp = :orange
    ordered_sim_conditions = ["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"]
    order_indices = [findfirst(x -> x == condition, T_steady.sim_id) for condition in ordered_sim_conditions]
    plot1 = bar(
        [1, 2], [T_steady[order_indices[1], :T_mod][4], T_steady[order_indices[1], :T_exp][4]],
        title="TI-11 Steady State Temperature", ylabel="Temperature (K)", xlabel="Experimental Runs",
        bar_width=0.4, xticks=(1:2:30, ordered_sim_conditions),
       label="T_model", color=[color_model, color_exp], ylimit =(0, 1000))
 
    # Loop through the remaining data to add the bars without labels
    for i in 2:15
        bar!(plot1, [2i-1, 2i], [T_steady[order_indices[i], :T_mod][4], T_steady[order_indices[i], :T_exp][4]], bar_width=0.4, label=["" ""], color=[color_model, color_exp])
    end
    # Manually add a legend entry for T_exp and T_model
    plot1 = bar!(plot1, [0], [0], label="T_exp", color=color_exp)
    xlims!(plot1, (0, 31))
end
    
    # plot!(
    #     psol.t,
    #     psol.u[Tf(t, x)][:, 77],
    #     label="gas 106mm")
    # plot5 = plot(
    #     psol.t,
    #     psol.u[Ts(t, x)][:, 59], # around 82 mm in T11
    #     title="Solid Temperature Profile T11",
    #     label="model",
    #     xlabel="Time (s)",
    #     ylabel="Temperature (K)")
    begin #TI-12
        color_model = :blue
        color_exp = :orange
        ordered_sim_conditions = ["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"]
        order_indices = [findfirst(x -> x == condition, T_steady.sim_id) for condition in ordered_sim_conditions]
        plot1 = bar(
            [1, 2], [T_steady[order_indices[1], :T_mod][5], T_steady[order_indices[1], :T_exp][5]],
            title="TI-12 Steady State Temperature", ylabel="Temperature (K)", xlabel="Experimental Runs",
            bar_width=0.4, xticks=(1:2:30, ordered_sim_conditions),
           label="T_model", color=[color_model, color_exp], ylimit =(0, 1200))
     
        # Loop through the remaining data to add the bars without labels
        for i in 2:15
            bar!(plot1, [2i-1, 2i], [T_steady[order_indices[i], :T_mod][5], T_steady[order_indices[i], :T_exp][5]], bar_width=0.4, label=["" ""], color=[color_model, color_exp])
        end
        # Manually add a legend entry for T_exp and T_model
        plot1 = bar!(plot1, [0], [0], label="T_exp", color=color_exp)
        xlims!(plot1, (0, 31))
    end
    # plot!(
    #     psol.t,
    #     psol.u[Tf(t, x)][:, 59],
    #     label="gas 82mm")
    # plot6 = plot(
    #     psol.t,
    #     psol.u[Ts(t, x)][:, 20], # around 27 mm in T12
    #     title="Solid Temperature Profile T12",
    #     label="model",
    #     xlabel="Time (s)",
    #     ylabel="Temperature (K)")
    begin #TI-3
        color_model = :blue
        color_exp = :orange
        ordered_sim_conditions = ["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"]
        order_indices = [findfirst(x -> x == condition, T_steady.sim_id) for condition in ordered_sim_conditions]
        plot1 = bar(
            [1, 2], [T_steady[order_indices[1], :T_mod][end][end], T_steady[order_indices[2], :T_exp][end][end]],
            title="TI-3 Steady State Temperature", ylabel="Temperature (K)", xlabel="Experimental Runs",
            bar_width=0.4, xticks=(1:2:30, ordered_sim_conditions),
           label="T_model", color=[color_model, color_exp], ylimit =(0, 1000))
     
        # Loop through the remaining data to add the bars without labels
        for i in 2:15
            bar!(plot1, [2i-1, 2i], [T_steady[order_indices[i], :T_mod][6], T_steady[order_indices[i], :T_exp][6]], bar_width=0.4, label=["" ""], color=[color_model, color_exp])
        end
        # Manually add a legend entry for T_exp and T_model
        plot1 = bar!(plot1, [0], [0], label="T_exp", color=color_exp)
        xlims!(plot1, (0, 31))
    end
 
    # plot!(
    #     psol.t,
    #     psol.u[Tf(t, x)][:, 20],
    #     label="gas 27mm")
# Extract T_mod values
T_mod_values = [row.T_mod for row in eachrow(T_steady)]
begin
# Plot T_mod as a function of x
plot(
    x_num1, T_mod_values,
    label="T_mod",
    xlabel="Length (m)",
    ylabel="Temperature (K)",
    title="T_mod vs Length",
    legend=:topright
)
end


# Initialize the DataFrame
    T_steady = DataFrame(sim_id=[], T_mod=[], T_exp=[])
    sim_key = collect(keys(simulation_conditions))
    lossr = zeros(length(sim_key))
   #Threads.@threads 
   for it in 1:length(sim_key)
    # Retrieve from measurements the experimental data for the current simulation condition
        sm = sim_key[it]
        #println(sm)
        cond = simulation_conditions[sm]
        expdata = (measurements[measurements.simulation_id.==sm, :temperatures][1])
        time_opt = (measurements[measurements.simulation_id.==sm, :time][1])
        
        #run selected simulation and get the steady temperature values
        temp_T = remakeAysha(pnew, cond, time_opt)
        #push!(T_steady, [sm, temp_T, expdata])
        push!(T_steady, (sm, temp_T, expdata))
   end
# Extract T_mod and T_exp from T_steady (TI-8)
T_mod = [row[2][1] for row in eachrow(T_steady)]
T_exp = [row[3][1] for row in eachrow(T_steady)]

# Create scatter plot
scatter(T_exp, T_mod, label="Experimental and Model Temperatures", legend=:bottomright,
        xlabel="Experimental Temperature (K)", ylabel="Model Temperature (K)", title=" External Solid Temperature (TI-8)")

# Add the best fit line to the plot
plot!([490, 1250], [490, 1250], label="Best Fit", color=:red)

# Extract T_mod and T_exp from T_steady (TI-9)
T_mod = [row[2][2] for row in eachrow(T_steady)]
T_exp = [row[3][2] for row in eachrow(T_steady)]

# Create scatter plot
scatter(T_exp, T_mod, label="Experimental and Model Temperatures", legend=:bottomright,
        xlabel="Experimental Temperature (K)", ylabel="Model Temperature (K)", title=" External Solid Temperature (TI-9)")

# Add the best fit line to the plot
plot!([550, 1100], [550, 1100], label="Best Fit", color=:red)

# Extract T_mod and T_exp from T_steady (TI-10)
T_mod = [row[2][3] for row in eachrow(T_steady)]
T_exp = [row[3][3] for row in eachrow(T_steady)]

# Create scatter plot
scatter(T_exp, T_mod, label="Experimental and Model Temperatures", legend=:bottomright,
        xlabel="Experimental Temperature (K)", ylabel="Model Temperature (K)", title=" External Solid Temperature (TI-10)")

# Add the best fit line to the plot
plot!([550, 1000], [550, 1000], label="Best Fit", color=:red)

# Extract T_mod and T_exp from T_steady (TI-11)
T_mod = [row[2][4] for row in eachrow(T_steady)]
T_exp = [row[3][4] for row in eachrow(T_steady)]

# Create scatter plot
scatter(T_exp, T_mod, label="Experimental and Model Temperatures", legend=:bottomright,
        xlabel="Experimental Temperature (K)", ylabel="Model Temperature (K)", title=" Internal Solid Temperature (TI-11)")

# Add the best fit line to the plot
plot!([550, 1000], [550, 1000], label="Best Fit", color=:red)

# Extract T_mod and T_exp from T_steady (TI-12)
T_mod = [row[2][5] for row in eachrow(T_steady)]
T_exp = [row[3][5] for row in eachrow(T_steady)]

# Create scatter plot
scatter(T_exp, T_mod, label="Experimental and Model Temperatures", legend=:bottomright,
        xlabel="Experimental Temperature (K)", ylabel="Model Temperature (K)", title=" Internal Solid Temperature (TI-12)")

# Add the best fit line to the plot
plot!([550, 1150], [550, 1150], label="Best Fit", color=:red)

# Extract T_mod and T_exp from T_steady (TI-3)
T_mod = [row[2][6] for row in eachrow(T_steady)]
T_exp = [row[3][6] for row in eachrow(T_steady)]

# Create scatter plot
scatter(T_exp, T_mod, label="Experimental and Model Temperatures", legend=:bottomright,
        xlabel="Experimental Temperature (K)", ylabel="Model Temperature (K)", title=" Air Outlet Temperature (TI-3)")

# Add the best fit line to the plot
plot!([550, 1100], [550, 1100], label="Best Fit", color=:red)