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
    A_st = 176e-6 #total area - solid m2
    A_ft = 49e-6 # total area - fluid m2
    channel_w = 1.5 / 1000 #m
    A_channel = channel_w * channel_w  #m2 (1.5x1.5mm2) 
    n_channel = 100
    A_exchange = 162.5e-6 #m2 - contact area between fluid and solid
    Vs = A_st * L #* deltax
    Vf = A_ft * L #* deltax
    qlpm = 7.12 #lpm
    q = qlpm / (1000 * 60) #m3/s
    m = 3.2923e-4 #kg/s
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
    Cpf = 1090 #kJ/kg.K
    mu = 2.0921e-5 #Pa.s
    #Gz = (w_t / L) * Re * Pr
    # A = 3.657
    # B = 0.5272
    # C = 55.1239
    # n = 0.3056
    #Nu = A*(1+(B*((Gz)^n)*exp(-C/Gz)))
    #Nu = 3.657 
    Lc = 4 * (w_t * w_t) / (4 * w_t)
    w = 1.5e-3 #width of channel (m)
    Vi = w * w * n_channel * L #m3
    Av = 4 * (w * L) / (w^2 * L) #specific area (m-1)
    hext = 10 #W/m2.K
    kins = 0.078 #W/m*K
    r0 = 23 / 1000 #m
    r_ins = 42 / 1000 #m
    r_H = 4 * (w_t * w_t) / (4 * w_t) #hydraulic receiver diameter
end;
#for interpolations 
#1. extract T2 data
#Exp 71 
begin
    D3 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0071_T2.xlsx")["Sheet 1 - Data_FPT0071_T2"]["A3:D7087"]
    y1d3_data = D3[:, 2] .+ 273.0 #T2 (insulation)
    x1 = D3[:, 1]
    #2.create interpolation function
    Tins = linear_interpolation(x1, y1d3_data, extrapolation_bc=y1d3_data[end])
    Tins_f(t) = Tins(t)
    @register_symbolic Tins_f(t)
end

begin
    # Parameters, variables, and derivatives for system 1
    @variables t x
    @parameters ρsCps ks h_average #ks h_average A n
    @parameters I0 v
    @variables Ts(..) Tf(..)
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2
    e = 0.62
    #ks = (48.78 / 1000) * 1.97 * ((1 - e)^1.5) #kW/m.K
    ρs = 3200 * e #kg/m3
    Cps = 1290 / 1000  #kJ/kg*K
    Re = (ρf * v * w_t) / mu
    Pr = (Cpf * mu) / kf
    Gz = (1/L) * Re * Pr * w_t
    #Interpolation for Gz number
    # x11 = 0.0001:0.001383838383838384:0.137 #T2 (insulation)
    # Gz = (1 ./ x11) * Re * Pr * w_t
    # #2.create interpolation function
    # Gz_ = linear_interpolation(x11, Gz)
    # Gz_f(x) = Gz_(x)
    # @register_symbolic Gz_f(x)
    #Cps(Ts) = (0.27+0.135e-4*(Ts)-9720*((Ts)^-2)+0.204e-7*((Ts)^2))/1000 #kJ/kg*K from manufacturer data
    p_opt = [ρsCps => 706048., ks => 22.5, h_average => 50000.]#[A => 20., n => 0.005]
    p_cond = [I0 => 456000.0, v => 1.22]
    p_math = vcat(p_opt, p_cond)
    #Nu = A * (Gz_f(x)^n)
    #Nu = A * (Gz)^n

    #p_math_vec = collect(p_math)
    #Nu = A * (1 + (B * ((Gz_f(x))^n) * exp(-C / Gz_f(x))))
    #nu = 4.364 * (1 + (0.7 * ((Gz_f(0.134))^10) * exp(-40 / Gz_f(0.134))))
    #h_average = (Nu * kf) / Lc

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
        A_st * (ρsCps) * Dt(Ts(t, x)) ~ A_st * (ks) * Dxx(Ts(t, x)) - (((h_average)/ Av) * ((Ts(t, x)) - Tf(t, x))) .- (kins * (r_ins/r0).* (Ts(t, x) .- Tins_f(t)) * L / (r_ins - r0)),
        A_ft * ρf * Cpf * Dt(Tf(t, x)) ~ A_ft * kf * Dxx(Tf(t, x)) -  m * Cpf * Dx(Tf(t, x)) + (((h_average)/ Av) * ((Ts(t, x)) - Tf(t, x)))
    ]
    bcs1 = [
        Ts(0.0, x) ~ Tamb, # initial
        Tf(0.0, x) ~ Tamb, # initial
        -A_st * (ks) * Dx(Ts(t, x_max1)) ~ 0.0, # far right
        -A_st * (ks) * Dx(Ts(t, x_min1)) ~ I0 * A_st - ϵ * σ * A_st * (Ts(t, x_min1)^4 - Tamb^4) - hext * A_st * (Ts(t, x_min1) - Tamb),  # far left
        -A_ft * kf * Dx(Tf(t, x_max1)) ~ 0.0, #-ρf * Cpf * V * A_ft * (Tf(t, x_max1) - Tamb), # exiting fluid
        -A_ft * kf * Dx(Tf(t, x_min1)) ~ m * Cpf * (Tf(t, x_min1) - Tamb) # entering fluid (upstream temperature)
    ]
    # Space and time domain for system 1
    domains1 = [t ∈ Interval(t_min, t_max),
        x ∈ Interval(x_min1, x_max1)]

    # ODE system for system 1
    @named pdesys = PDESystem(eq1, bcs1, domains1, [t, x], [Ts(t, x), Tf(t, x)], p_math)

    #end
    #begin

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
    sol1.u[Ts(t, x)][:, 40] #T9 model 55mm (external)
    sol1.u[Ts(t, x)][:, 77] #T10 model 106mm (external)
    sol1.u[Tf(t, x)][:, 59] #T11 model 82mm (internal-gas)
    sol1.u[Ts(t, x)][:, 59] #T11 model 82mm (internal-solid)
    sol1.u[Tf(t, x)][:, 20] #T12 model 27mm (internal-gas)
    sol1.u[Ts(t, x)][:, 20] #T12 model 27mm (internal-solid)
end

#measurements and conditions#Defining simulation conditions
begin
    condition_E67 = Dict(:I0 => 456000.0, :v => 1.22)
    condition_E68 = Dict(:I0 => 456000.0, :v => 1.00)
    condition_E69 = Dict(:I0 => 456000.0, :v => 0.84)
    condition_E70 = Dict(:I0 => 456000.0, :v => 0.73)
    condition_E71 = Dict(:I0 => 456000.0, :v => 0.57)
    condition_E72 = Dict(:I0 => 304000.0, :v => 1.46)
    condition_E73 = Dict(:I0 => 304000.0, :v => 1.05)
    condition_E74 = Dict(:I0 => 304000.0, :v => 0.72)
    condition_E75 = Dict(:I0 => 304000.0, :v => 0.55)
    condition_E76 = Dict(:I0 => 304000.0, :v => 0.36)
    condition_E77 = Dict(:I0 => 256000.0, :v => 1.10)
    condition_E78 = Dict(:I0 => 256000.0, :v => 0.80)
    condition_E79 = Dict(:I0 => 256000.0, :v => 0.64)
    condition_E80 = Dict(:I0 => 256000.0, :v => 0.53)
    condition_E81 = Dict(:I0 => 256000.0, :v => 0.36)

    simulation_conditions = Dict("E67" => condition_E67, "E68" => condition_E68,
        "E69" => condition_E69, "E70" => condition_E70,
        "E71" => condition_E71, "E72" => condition_E72,
        "E73" => condition_E73, "E74" => condition_E74,
        "E75" => condition_E75, "E76" => condition_E76,
        "E77" => condition_E77, "E78" => condition_E78,
        "E79" => condition_E79, "E80" => condition_E80,
        "E81" => condition_E81)
    #Defining measurement data
    measurements = DataFrame(
        simulation_id=repeat(["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"], inner=6, outer=1),
        obs_id=repeat(["_T8", "_T9", "_T10", "_T11", "_T12", "_T3"], inner=1, outer=15),
        time=repeat([3929, 5362, 5363, 6702, 7084, 3214, 4572, 6015, 6351, 7144, 3041, 5381, 5230, 5811, 5986], inner=6, outer=1),
        temperatures=[965.407, 975.144, 825.592, 880.867, 1004.165, 763.859,
            1031.574, 1023.115, 850.099, 898.754, 1050.691, 773.207,
            1070.803, 1045.898, 852.837, 896.788, 1072.727, 769.76,
            1167.978, 1093.849, 871.496, 912.173, 1120.42, 779.53,
            1210.945, 1095.322, 847.476, 882.823, 1120.417, 753.56,
            742.125, 778.592, 684.246, 736.626, 807.125, 652.955,
            844.257, 870.26, 747.444, 791.958, 898.081, 694.626,
            962.113, 938.106, 767.803, 804.082, 965.702, 697.672,
            1015.081, 954.214, 757.393, 788.678, 979.795, 681.066,
            1069.567, 947.372, 726.308, 751.159, 970.498, 647.019,
            569.248, 604.984, 543.727, 574.299, 627.16, 525.356,
            634.731, 664.296, 583.704, 612.092, 686.936, 554.455,
            677.817, 694.156, 595.766, 622.325, 716.314, 560.033,
            711.537, 713.686, 601.77, 626.485, 735.984, 561.254,
            763.299, 729.461, 597.766, 618.975, 751.15, 550.499])
end

#Exp. data to extract temp.
begin
    #Exp 67 - T3, T8
    Z = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0067_231125_161757.xlsx")["Sheet 1 - Data_FPT0067_231125_1"]["A3:C3932"]
    xz_data = Z[:, 1]
    y1z_data = Z[:, 2] .+ 273.15
    y2z_data = Z[:, 3] .+ 273.15

    #Exp 67 - T9, T10, T11, T12
    Z1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0067_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0067_T9"]["A3:E3932"]
    y3z_data = Z1[:, 2] .+ 273.15
    y4z_data = Z1[:, 3] .+ 273.15
    y5z_data = Z1[:, 4] .+ 273.15
    y6z_data = Z1[:, 5] .+ 273.15

    #Exp 68 - T3, T8
    A1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0068_231126_115725.xlsx")["Sheet 1 - Data_FPT0068_231126_1"]["A3:C5365"]
    xa1_data = A1[:, 1]
    y1a1_data = A1[:, 2] .+ 273.15
    y2a1_data = A1[:, 3] .+ 273.15

    #Exp 68 - T9, T10, T11, T12
    A11 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0068_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0068_231126_1"]["A3:E5365"]
    y3a_data = A11[:, 2] .+ 273.15
    y4a_data = A11[:, 3] .+ 273.15
    y5a_data = A11[:, 4] .+ 273.15
    y6a_data = A11[:, 5] .+ 273.15

    #Exp 69 - T3, T8
    B1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0069_231126_140153.xlsx")["Sheet 1 - Data_FPT0069_231126_1"]["A3:C5366"]
    xb1_data = B1[:, 1]
    y1b1_data = B1[:, 2] .+ 273.15
    y2b1_data = B1[:, 3] .+ 273.15

    #Exp 69 - T9, T10, T11, T12
    B11 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0069_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0069_231126_1"]["A3:E5366"]
    y3b_data = B11[:, 2] .+ 273.15
    y4b_data = B11[:, 3] .+ 273.15
    y5b_data = B11[:, 4] .+ 273.15
    y6b_data = B11[:, 5] .+ 273.15

    #Exp 70 - T3, T8
    C1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0070_231127_090339.xlsx")["Sheet 1 - Data_FPT0070_231127_0"]["A3:C6705"]
    xc1_data = C1[:, 1]
    y1c1_data = C1[:, 2] .+ 273.15
    y2c1_data = C1[:, 3] .+ 273.15

    #Exp 70 - T9, T10, T11, T12

    C11 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0070_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0070_231127_0"]["A3:E6705"]
    y3c_data = C11[:, 2] .+ 273.15
    y4c_data = C11[:, 3] .+ 273.15
    y5c_data = C11[:, 4] .+ 273.15
    y6c_data = C11[:, 5] .+ 273.15

    #Exp 71 - T3, T8
    D1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0071_231128_102707.xlsx")["Sheet 1 - Data_FPT0071_231128_1"]["A3:C7087"]
    xd1_data = D1[:, 1]
    y1d1_data = D1[:, 2] .+ 273.15
    y2d1_data = D1[:, 3] .+ 273.15

    #Exp 71 - T9, T10, T11, T12

    D11 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0071_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0071_231128_1"]["A3:E7087"]
    y3d_data = D11[:, 2] .+ 273.15
    y4d_data = D11[:, 3] .+ 273.15
    y5d_data = D11[:, 4] .+ 273.15
    y6d_data = D11[:, 5] .+ 273.15

    #Exp 72 - T3, T8
    E1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0072_231129_104140.xlsx")["Sheet 1 - Data_FPT0072_231129_1"]["A3:C3217"]
    xe1_data = E1[:, 1]
    y1e1_data = E1[:, 2] .+ 273.15
    y2e1_data = E1[:, 3] .+ 273.15

    #Exp 72 - T9, T10, T11, T12

    E11 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0072_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0072_231129_1"]["A3:E3217"]
    y3e_data = E11[:, 2] .+ 273.15
    y4e_data = E11[:, 3] .+ 273.15
    y5e_data = E11[:, 4] .+ 273.15
    y6e_data = E11[:, 5] .+ 273.15

    #Exp 73 - T3, T8
    F1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0073_231129_132744.xlsx")["Sheet 1 - Data_FPT0073_231129_1"]["A3:C4575"]
    xf1_data = F1[:, 1]
    y1f1_data = F1[:, 2] .+ 273.15
    y2f1_data = F1[:, 3] .+ 273.15

    #Exp 73 - T9, T10, T11, T12

    F11 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0073_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0073_231129_1"]["A3:E4575"]
    y3f_data = F11[:, 2] .+ 273.15
    y4f_data = F11[:, 3] .+ 273.15
    y5f_data = F11[:, 4] .+ 273.15
    y6f_data = F11[:, 5] .+ 273.15
    
    #Exp 74 - T3, T8
    G = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0074_231130_123228.xlsx")["Sheet 1 - Data_FPT0074_231130_1"]["A3:C6018"]
    xg_data = G[:, 1]
    y1g_data = G[:, 2] .+ 273.15
    y2g_data = G[:, 3] .+ 273.15

    #Exp 74 - T9, T10, T11, T12

    G1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0074_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0074_231130_1"]["A3:E6018"]
    y3g1_data = G1[:, 2] .+ 273.15
    y4g1_data = G1[:, 3] .+ 273.15
    y5g1_data = G1[:, 4] .+ 273.15
    y6g1_data = G1[:, 5] .+ 273.15

    #Exp 75 - T3, T8
    H = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0075_231201_162138.xlsx")["Sheet 1 - Data_FPT0075_231201_1"]["A3:C6354"]
    xh_data = H[:, 1]
    y1h_data = H[:, 2] .+ 273.15
    y2h_data = H[:, 3] .+ 273.15

    #Exp 75 - T9, T10, T11, T12

    H1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0075_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0075_231201_1"]["A3:E6354"]
    y3h1_data = H1[:, 2] .+ 273.15
    y4h1_data = H1[:, 3] .+ 273.15
    y5h1_data = H1[:, 4] .+ 273.15
    y6h1_data = H1[:, 5] .+ 273.15

    #Exp 76 - T3, T8
    I = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0076_231203_120521.xlsx")["Sheet 1 - Data_FPT0076_231203_1"]["A3:C7147"]
    xi_data = I[:, 1]
    y1i_data = I[:, 2] .+ 273.15
    y2i_data = I[:, 3] .+ 273.15

    #Exp 76 - T9, T10, T11, T12

    I1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0076_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0076_231203_1"]["A3:E7147"]
    y3i1_data = I1[:, 2] .+ 273.15
    y4i1_data = I1[:, 3] .+ 273.15
    y5i1_data = I1[:, 4] .+ 273.15
    y6i1_data = I1[:, 5] .+ 273.15

    #Exp 77 - T3, T8
    J = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0077_231203_161315.xlsx")["Sheet 1 - Data_FPT0077_231203_1"]["A3:C3044"]
    xj_data = J[:, 1]
    y1j_data = J[:, 2] .+ 273.15
    y2j_data = J[:, 3] .+ 273.15

    #Exp 77 - T9, T10, T11, T12

    J1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0077_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0077_231203_1"]["A3:E3044"]
    y3j1_data = J1[:, 2] .+ 273.15
    y4j1_data = J1[:, 3] .+ 273.15
    y5j1_data = J1[:, 4] .+ 273.15
    y6j1_data = J1[:, 5] .+ 273.15

    #Exp 78 - T3, T8
    K = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0078_231204_132252.xlsx")["Sheet 1 - Data_FPT0078_231204_1"]["A3:C5384"]
    xk_data = K[:, 1]
    y1k_data = K[:, 2] .+ 273.15
    y2k_data = K[:, 3] .+ 273.15

    #Exp 78 - T9, T10, T11, T12

    K1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0078_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0078_231204_1"]["A3:E5384"]
    y3k1_data = K1[:, 2] .+ 273.15
    y4k1_data = K1[:, 3] .+ 273.15
    y5k1_data = K1[:, 4] .+ 273.15
    y6k1_data = K1[:, 5] .+ 273.15
    
    #Exp 79 - T3, T8
    L = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0079_231204_172244.xlsx")["Sheet 1 - Data_FPT0079_231204_1"]["A3:C5233"]
    xl_data = L[:, 1]
    y1l_data = L[:, 2] .+ 273.15
    y2l_data = L[:, 3] .+ 273.15  
    
    #Exp 79 - T9, T10, T11, T12

    L1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0079_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0079_231204_1"]["A3:E5233"]
    y3l1_data = L1[:, 2] .+ 273.15
    y4l1_data = L1[:, 3] .+ 273.15
    y5l1_data = L1[:, 4] .+ 273.15
    y6l1_data = L1[:, 5] .+ 273.15
     
      
    #Exp 80 - T3, T8
    M = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0080_231205_095122.xlsx")["Sheet 1 - Data_FPT0080_231205_0"]["A3:C5814"]
    xm_data = M[:, 1]
    y1m_data = M[:, 2] .+ 273.15
    y2m_data = M[:, 3] .+ 273.15   

    #Exp 80 - T9, T10, T11, T12

    M1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0080_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0080_231205_0"]["A3:E5814"]
    y3m1_data = M1[:, 2] .+ 273.15
    y4m1_data = M1[:, 3] .+ 273.15
    y5m1_data = M1[:, 4] .+ 273.15
    y6m1_data = M1[:, 5] .+ 273.15
     

    #Exp 81 - T3, T8
    N = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0081_231205_135354.xlsx")["Sheet 1 - Data_FPT0081_231205_1"]["A3:C5989"]
    xn_data = N[:, 1]
    y1n_data = N[:, 2] .+ 273.15
    y2n_data = N[:, 3] .+ 273.15   

    #Exp 81 - T9, T10, T11, T12

    N1 = XLSX.readxlsx("./SolarSimulator/EXCEL/Data_FPT0081_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0081_231205_1"]["A3:E5989"]
    y3n1_data = N1[:, 2] .+ 273.15
    y4n1_data = N1[:, 3] .+ 273.15
    y5n1_data = N1[:, 4] .+ 273.15
    y6n1_data = N1[:, 5] .+ 273.15
end
#Optimization using NLOpt
function NLmodeloptim(tvalues, p_math_vec)

    #p = [hlocal => p_vary[1]]
    modeloptim = remake(prob, p=p_math_vec, tspan=(1.0, tvalues[end]))
    modeloptim_sol = solve(modeloptim, FBDF(), saveat=tvalues[end])#, reltol=1e-12, abstol=1e-12)
    #time = modelfit_sol.t
    tempT8_op = modeloptim_sol.u[Ts(t, x)][end, 4]
    tempT9_op = modeloptim_sol.u[Ts(t, x)][end, 40]
    tempT10_op = modeloptim_sol.u[Ts(t, x)][end, 77]
    tempT3_op = modeloptim_sol.u[Tf(t, x)][end, end-1]
    T12_modelmean = (modeloptim_sol.u[Tf(t, x)][end, 20] .+ modeloptim_sol.u[Ts(t, x)][end, 20]) ./ 2
    T11_modelmean = (modeloptim_sol.u[Tf(t, x)][end, 59] .+ modeloptim_sol.u[Ts(t, x)][end, 59]) ./ 2

    #return append!(tempT8_op, tempT9_op, tempT10_op, tempT3_op, T12_modelmean, T11_modelmean)
    return [tempT8_op, tempT9_op, tempT10_op, tempT3_op, T12_modelmean, T11_modelmean]
end

#the first entry in the simulation_conditions
pguess = p_opt

function lossAll(pguess_l, _)

    #to place the conditions loop

    sim_key = collect(keys(simulation_conditions))
    lossr = zeros(length(sim_key))
   Threads.@threads for it in 1:length(sim_key)
        # Retrieve from measurements the experimental data for the current simulation condition
        sm = sim_key[it]
        #println(sm)
        cond = simulation_conditions[sm]
        expdata = (measurements[measurements.simulation_id.==sm, :temperatures])
        time_opt = (measurements[measurements.simulation_id.==sm, :time])
        p_math_vec = copy(cond)
        j = 1
        for (k, v) in p_opt #FIX DICTIONARY PARAMETERS OR VECTORS
            merge!(p_math_vec, Dict(Symbol(k) => pguess_l[j]))
            j += 1
        end  # pguess is the initial guess for the optimization
        temp_T = NLmodeloptim(time_opt, p_math_vec)
        temp_error = (temp_T .- expdata) .^ 2
        lossr[it] = sqrt(sum(temp_error))
    end

    return sum(lossr) #MSE
end

optf = OptimizationFunction(lossAll, Optimization.AutoForwardDiff())

lb = [0.0, 0.0, 0.0]
ub = [1000000., 50., 100000.]
pguess_opt = [x[2] for x in p_opt]
initialerror = (lossAll(pguess_opt, []))
println(initialerror)

optprob = Optimization.OptimizationProblem(optf, pguess_opt, [], lb=lb, ub=ub)
#optsol = solve(optprob, NLopt.GN_MLSL_LDS(), local_method=NLopt.LN_NELDERMEAD(), maxtime=10, local_maxiters=10)
optsol = solve(optprob, NLopt.LN_NELDERMEAD(), maxtime=100, local_maxiters=10)

println(optsol.retcode)
pnew = optsol.u
println(pnew)
res_error = lossAll(pnew, [])
display(res_error)

# begin

#     #p_opt = [A => 8.0, B => 0.5, C => 55.0, n => 0.2]
#     case = "E71" #"E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E78", "E79", "E80", "E81")
#     p_cond = simulation_conditions[case]
#     j=1
#     for (k, v) in p_opt #FIX DICTIONARY PARAMETERS OR VECTORS
#         merge!(p_cond, Dict(Symbol(k) => pnew[j]))
#         j += 1
#     end
#     modelfit = remake(prob, p=p_cond)
#     modelfit_sol = solve(modelfit, FBDF(), saveat=xd1_data)


#     psol = modelfit_sol
# end

# plot1 = plot(
#         psol.t,
#         psol.u[Tf(t, x)][:, end-1], # around 136 mm in T3
#         title="Gas Temperature Profile T3",
#         label="model",
#         xlabel="Time (s)",
#         ylabel="Temperature (K)")
begin
plot1 =  bar([psol.u[Tf(t, x)][end, end-1]],[measurements[6, :temperatures]], label="E67", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Tf(t, x)][end, end-1]],[measurements[12, :temperatures]], label="E68", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Tf(t, x)][end, end-1]],[measurements[18, :temperatures]], label="E69", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Tf(t, x)][end, end-1]],[measurements[24, :temperatures]], label="E70", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Tf(t, x)][end, end-1]],[measurements[30, :temperatures]], label="E71", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Tf(t, x)][end, end-1]],[measurements[36, :temperatures]], label="E72", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Tf(t, x)][end, end-1]],[measurements[42, :temperatures]], label="E73", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Tf(t, x)][end, end-1]],[measurements[48, :temperatures]], label="E74", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Tf(t, x)][end, end-1]],[measurements[54, :temperatures]], label="E75", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Tf(t, x)][end, end-1]],[measurements[60, :temperatures]], label="E76", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Tf(t, x)][end, end-1]],[measurements[66, :temperatures]], label="E77", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Tf(t, x)][end, end-1]],[measurements[72, :temperatures]], label="E78", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Tf(t, x)][end, end-1]],[measurements[78, :temperatures]], label="E79", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Tf(t, x)][end, end-1]],[measurements[84, :temperatures]], label="E80", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Tf(t, x)][end, end-1]],[measurements[90, :temperatures]], label="E81", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")

# plot2 = plot(
#         psol.t,
#         psol.u[Ts(t, x)][:, 4], # around 5 mm in T8
#         title="Solid Temperature Profile T8",
#         label="model",
#         xlabel="Time (s)",
#         ylabel="Temperature (K)")

plot2 = bar([psol.u[Ts(t, x)][end, 4]],[measurements[1, :temperatures]], label="E67", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Ts(t, x)][end, 4]],[measurements[7, :temperatures]], label="E68", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Ts(t, x)][end, 4]],[measurements[13, :temperatures]], label="E69", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Ts(t, x)][end, 4]],[measurements[19, :temperatures]], label="E70", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Ts(t, x)][end, 4]],[measurements[25, :temperatures]], label="E71", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Ts(t, x)][end, 4]],[measurements[31, :temperatures]], label="E72", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Ts(t, x)][end, 4]],[measurements[37, :temperatures]], label="E73", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Ts(t, x)][end, 4]],[measurements[43, :temperatures]], label="E74", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)") 
    bar!([psol.u[Ts(t, x)][end, 4]],[measurements[49, :temperatures]], label="E75", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Ts(t, x)][end, 4]],[measurements[55, :temperatures]], label="E76", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Ts(t, x)][end, 4]],[measurements[61, :temperatures]], label="E77", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Ts(t, x)][end, 4]],[measurements[67, :temperatures]], label="E78", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Ts(t, x)][end, 4]],[measurements[73, :temperatures]], label="E79", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Ts(t, x)][end, 4]],[measurements[79, :temperatures]], label="E80", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
    bar!([psol.u[Ts(t, x)][end, 4]],[measurements[85, :temperatures]], label="E81", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")

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
plot3 = bar([psol.u[Ts(t, x)][end, 40]],[measurements[2, :temperatures]], label="E67", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
        bar!([psol.u[Ts(t, x)][end, 40]],[measurements[8, :temperatures]], label="E68", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
        bar!([psol.u[Ts(t, x)][end, 40]],[measurements[14, :temperatures]], label="E69", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
        bar!([psol.u[Ts(t, x)][end, 40]],[measurements[20, :temperatures]], label="E70", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
        bar!([psol.u[Ts(t, x)][end, 40]],[measurements[26, :temperatures]], label="E71", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
        bar!([psol.u[Ts(t, x)][end, 40]],[measurements[32, :temperatures]], label="E72", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
        bar!([psol.u[Ts(t, x)][end, 40]],[measurements[38, :temperatures]], label="E73", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
        bar!([psol.u[Ts(t, x)][end, 40]],[measurements[44, :temperatures]], label="E74", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
        bar!([psol.u[Ts(t, x)][end, 40]],[measurements[50, :temperatures]], label="E75", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
        bar!([psol.u[Ts(t, x)][end, 40]],[measurements[56, :temperatures]], label="E76", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
        bar!([psol.u[Ts(t, x)][end, 40]],[measurements[62, :temperatures]], label="E77", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
        bar!([psol.u[Ts(t, x)][end, 40]],[measurements[68, :temperatures]], label="E78", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
        bar!([psol.u[Ts(t, x)][end, 40]],[measurements[74, :temperatures]], label="E79", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
        bar!([psol.u[Ts(t, x)][end, 40]],[measurements[80, :temperatures]], label="E80", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
        bar!([psol.u[Ts(t, x)][end, 40]],[measurements[86, :temperatures]], label="E81", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")

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
plot4 = bar([psol.u[Ts(t, x)][end, 77]],[measurements[3, :temperatures]], label="E67", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
        bar!([psol.u[Ts(t, x)][end, 77]],[measurements[9, :temperatures]], label="E68", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
        bar!([psol.u[Ts(t, x)][end, 77]],[measurements[15, :temperatures]], label="E69", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
        bar!([psol.u[Ts(t, x)][end, 77]],[measurements[21, :temperatures]], label="E70")
        bar!([psol.u[Ts(t, x)][end, 77]],[measurements[27, :temperatures]], label="E71")
        bar!([psol.u[Ts(t, x)][end, 77]],[measurements[33, :temperatures]], label="E72")
        bar!([psol.u[Ts(t, x)][end, 77]],[measurements[39, :temperatures]], label="E73")
        bar!([psol.u[Ts(t, x)][end, 77]],[measurements[45, :temperatures]], label="E74")
        bar!([psol.u[Ts(t, x)][end, 77]],[measurements[51, :temperatures]], label="E75")
        bar!([psol.u[Ts(t, x)][end, 77]],[measurements[57, :temperatures]], label="E76")
        bar!([psol.u[Ts(t, x)][end, 77]],[measurements[63, :temperatures]], label="E77")
        bar!([psol.u[Ts(t, x)][end, 77]],[measurements[69, :temperatures]], label="E78")
        bar!([psol.u[Ts(t, x)][end, 77]],[measurements[75, :temperatures]], label="E79")
        bar!([psol.u[Ts(t, x)][end, 77]],[measurements[81, :temperatures]], label="E80")
        bar!([psol.u[Ts(t, x)][end, 77]],[measurements[87, :temperatures]], label="E81")
    
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
plot5 = bar([psol.u[Ts(t, x)][end, 59]],[measurements[4, :temperatures]], label="E67", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
        bar!([psol.u[Ts(t, x)][end, 59]],[measurements[10, :temperatures]], label="E68")
        bar!([psol.u[Ts(t, x)][end, 59]],[measurements[16, :temperatures]], label="E69")
        bar!([psol.u[Ts(t, x)][end, 59]],[measurements[22, :temperatures]], label="E70")
        bar!([psol.u[Ts(t, x)][end, 59]],[measurements[28, :temperatures]], label="E71")
        bar!([psol.u[Ts(t, x)][end, 59]],[measurements[34, :temperatures]], label="E72")
        bar!([psol.u[Ts(t, x)][end, 59]],[measurements[40, :temperatures]], label="E73")
        bar!([psol.u[Ts(t, x)][end, 59]],[measurements[46, :temperatures]], label="E74")
        bar!([psol.u[Ts(t, x)][end, 59]],[measurements[52, :temperatures]], label="E75")
        bar!([psol.u[Ts(t, x)][end, 59]],[measurements[58, :temperatures]], label="E76")
        bar!([psol.u[Ts(t, x)][end, 59]],[measurements[64, :temperatures]], label="E77")
        bar!([psol.u[Ts(t, x)][end, 59]],[measurements[70, :temperatures]], label="E78")
        bar!([psol.u[Ts(t, x)][end, 59]],[measurements[76, :temperatures]], label="E79")
        bar!([psol.u[Ts(t, x)][end, 59]],[measurements[82, :temperatures]], label="E80")
        bar!([psol.u[Ts(t, x)][end, 59]],[measurements[88, :temperatures]], label="E81")
   
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
plot6 = bar([psol.u[Ts(t, x)][end, 20]],[measurements[5, :temperatures]], label="E67", ylabel="Experimental Temperature (K)", xlabel="Model Temperature (K)")
        bar!([psol.u[Ts(t, x)][end, 20]],[measurements[11, :temperatures]], label="E68")
        bar!([psol.u[Ts(t, x)][end, 20]],[measurements[17, :temperatures]], label="E69")
        bar!([psol.u[Ts(t, x)][end, 20]],[measurements[23, :temperatures]], label="E70")
        bar!([psol.u[Ts(t, x)][end, 20]],[measurements[29, :temperatures]], label="E71")
        bar!([psol.u[Ts(t, x)][end, 20]],[measurements[35, :temperatures]], label="E72")
        bar!([psol.u[Ts(t, x)][end, 20]],[measurements[41, :temperatures]], label="E73")
        bar!([psol.u[Ts(t, x)][end, 20]],[measurements[47, :temperatures]], label="E74")
        bar!([psol.u[Ts(t, x)][end, 20]],[measurements[53, :temperatures]], label="E75")
        bar!([psol.u[Ts(t, x)][end, 20]],[measurements[59, :temperatures]], label="E76")
        bar!([psol.u[Ts(t, x)][end, 20]],[measurements[65, :temperatures]], label="E77")
        bar!([psol.u[Ts(t, x)][end, 20]],[measurements[71, :temperatures]], label="E78")
        bar!([psol.u[Ts(t, x)][end, 20]],[measurements[77, :temperatures]], label="E79")
        bar!([psol.u[Ts(t, x)][end, 20]],[measurements[83, :temperatures]], label="E80")
        bar!([psol.u[Ts(t, x)][end, 20]],[measurements[89, :temperatures]], label="E81")
 
    # plot!(
    #     psol.t,
    #     psol.u[Tf(t, x)][:, 20],
    #     label="gas 27mm")
plot7 = plot(
        collect(x_num1),
        psol.u[Tf(t, x)][end, :],
        label="gas",
        xlabel="Length (m)",
        ylabel="Temperature (K)")
plot!(
        collect(x_num1),
        psol.u[Ts(t, x)][end, :],
        label="solid",
        xlabel="Length (m)",
        ylabel="Temperature (K)")

plot(plot1, plot2, plot3, plot4, plot5, plot6, plot7, layout=(10, 2), size=(1500, 1500))
end