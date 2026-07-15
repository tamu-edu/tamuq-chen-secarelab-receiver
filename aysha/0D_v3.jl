begin #libraries
    using MethodOfLines
    using ModelingToolkit
    using DomainSets, OrdinaryDiffEq
    using NonlinearSolve, DifferentialEquations, DataFrames
    using XLSX, Statistics, Symbolics, Interpolations
    using StatsPlots, ColorSchemes #, Plots
    using LsqFit, LossFunctions
    using Optimization, OptimizationNLopt, Symbolics, ForwardDiff, OptimizationMetaheuristics
    using SciMLStructures: replace!, Tunable
    using SymbolicIndexingInterface: parameter_values, state_values
    using Latexify

end

begin #FIXED Parameters

    Tamb = (22.448 + 273.15) #K (same for all exp) 

    # Receiver dimensions (COMSOL D_tot = 20mm, L_channel = 137mm)
    w_t = 20.e-3  # Width of the receiver in meters
    A_frt = w_t * w_t  # Total frontal area of the receiver (m^2)
    L = 137e-3 #m Length
    A_s_p = w_t * L * 4 #total area solid periphery m2

    # Channel dimensions (COMSOL R_channel = 1.65mm, t_channel = 0.35mm)
    w_chnl = 1.65e-3  # Width of a single channel (m)
    A_chnl_frt = w_chnl * w_chnl  # Frontal area of a single channel (m^2)
    A_chnl_p = w_chnl * L * 4  # Periphery area of a single channel (m^2)

    # Channel array
    n_chnl = 10 * 10  # Number of channels (10x10 array)
    A_chnl_frt_all = A_chnl_frt * n_chnl  # Total frontal area of all channels (m^2)
    Lc = 4 * (w_t * w_t) / (4 * w_t) #hydraulic receiver diameter (20mm)
    Dh = 4 * A_chnl_frt / (4 * w_chnl) #hydraulic diameter of the channel (m)

    # Solid geometrical calculations
    A_frt_solid = A_frt - A_chnl_frt_all  # Frontal area of just solid (m^2)
    A_exchange = A_chnl_p * n_chnl  # Total contact area between fluid and solid (m^2)

    # Volume calculations
    Vs = A_frt_solid * L  # Corrected solid volume (m^3)
    Vf = n_chnl * A_chnl_frt * L  # Fluid volume in the channels (m^3)

    # Alumina Adaptor dimensions (COMSOL Al_radius = 19.4mm, Al_L1 = 29mm, Al_L2 = 28mm, Al_tube = 13mm)
    rs2 = 19.4e-3 # radius of the second solid (m)
    Ls2 = 57.e-3 # total length of the adaptor (m)
    r_tube = 6.5e-3 # radius of the central gas exit tube (m)
    L_overlap = 29.e-3 # length of overlap with receiver (m)
    L_free = 28.e-3 # length of exit section (m)
    
    # Adaptor Volume: (cyl2 - receiver) in overlap + (cyl2 - tube) in free section
    Vs2 = (π * rs2^2 * L_overlap - A_frt * L_overlap) + π * (rs2^2 - r_tube^2) * L_free
    
    A_s1s2 = 4 * w_t * L_overlap # lateral contact area between receiver and adaptor (m^2)
    A_s_p2 = 2 * π * rs2 * Ls2 # periphery area of the second solid (m^2)

    # Constants
    σ = 5.670374419e-8 # W/m2.K^4 Stefan-Boltzmann constant (COMSOL value)
    hext = 10.0 # W/m2.K convective heat transfer coefficient for the front 

    # SiC Properties
    ϵ = 0.85 # emissivity
    α = 0.85 # absorptivity
    e = 0.425 # porosity

    # Insulation (3rd solid: felt) & Casing (4th solid: metal) parameters
    r_ins = 75.e-3 # outer radius of insulation felt (m, COMSOL R_ins = 75mm)
    t_metal = 18.e-3 # thickness of metal shell (m, COMSOL M_Length = 18mm)
    r_metal = r_ins + t_metal # outer radius of metal casing (m, 93mm)
    L_metal = L + L_free # total length of metal casing (m, 165mm)
    L_backplate = t_metal # backplate thickness (m, 18mm)

    # Volumes of insulation & metal
    # Insulation Volume: surrounds receiver (length 108mm) and adaptor (length 57mm)
    L_overlap_ins = 108.e-3 # monolith-only insulation length (m)
    L_adaptor_ins = 57.e-3 # adaptor insulation length (m)
    Vs_ins = (π * r_ins^2 - A_frt) * L_overlap_ins + π * (r_ins^2 - rs2^2) * L_adaptor_ins # m^3
    
    # Metal Volume: outer sleeve + backplate with tube hole
    Vs_metal = π * (r_metal^2 - r_ins^2) * L_metal + π * (r_metal^2 - r_tube^2) * L_backplate # m^3

    # Surface areas for loss calculations
    A_s_p3 = 2 * π * r_metal * L_metal # periphery area of metal casing (m^2)
    A_s_b3 = π * r_metal^2 # back area of metal casing (m^2)
    A_outer_casing = A_s_p3 + A_s_b3 # total outer area exposed to environment (m^2)

    # Material properties
    kins = 0.08 # W/m*K (COMSOL ins_k = 0.08)
    kmet = 201.0 # W/m*K (COMSOL M_k = 201)

    # Thermal resistance network for radial conduction through insulation
    r0_eq = sqrt(A_frt / π) # equivalent inner radius of square receiver (m)
    R_ins_total = log(r_ins / r0_eq) / (2 * π * kins * L_overlap_ins) # radial resistance monolith section
    R_out_monolith = R_ins_total / 2 # inner-midpoint resistance
    
    R_ins_adpt = log(r_ins / rs2) / (2 * π * kins * Ls2) # radial resistance adaptor section
    R_out_adpt = R_ins_adpt / 2 # inner-midpoint resistance
    
    # Parallel resistance for insulation midpoint (Ts3) to metal casing (Ts4)
    R_s3s4 = 1.0 / (1.0 / R_out_monolith + 1.0 / R_out_adpt)

end;

#Exp 71 One typical Insulation Temperature Profile
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
    # Parameters, variables, and derivatives for system 1
    #@independent_variables t 
    @parameters kfs_c cps_c I_loss A B C n #ks h_average A n
    @parameters Io qlpm Tinit
    @variables t Ts(t) Tf(t) Ts2(t) Ts3(t) Ts4(t)
    Dt = Differential(t)
    
    ### Fluid Properties
    RH = 0.0 #relative humidity
    @register_symbolic ρf_f(Tf)
    ρf_f(T) = 352.716 * T^-1 * (1+RH) /(1+ 1.609 * RH) #COMSOL + https://www.engineeringtoolbox.com/density-air-d_680.html with some error for w vs RH
  
    ρf = ρf_f(293.0) #kg/m3
    m(qlpm) = ρf_f(Tamb) * (qlpm / 60.0 / 1000.0) #kg/s mass at the flowmeter conditions T=25oC
    G(qlpm) = m(qlpm) / A_chnl_frt_all #kg/2/m2 mass at the flowmeter conditions T=25oC


    mu = 2.0921e-5 #Pa.s
    kf = 0.056 #W/m.K
    @register_symbolic kf_f(Tf)
    kf_f(T) = -0.00227583562 + 1.15480022e-4 * T^1 - 7.90252856e-8 * T^2 + 4.11702505e-11 * T^3 - 7.43864331e-15 * T^4 #COMSOL
    Cpf = 1090 #J/kg.K
    Cpw_f(T) = 1745.96354 + 0.185114553 * T^1 + 6.19448731e-4 * T^2 - 3.0267851e-7 * T^3 + 4.19053122e-11 * T^4 #COMSOL
    cpf_f(T) = (1.93e-10 * T^4 - 8.0e-7 * T^3 + 1.14e-3 * T^2 - 4.49e-1 * T + 1.06e3) + RH * Cpw_f(T)
    μf_f(T) = -8.38278E-7 + 8.35717342e-8 * T^1 - 7.69429583e-11 * T^2 + 4.6437266e-14 * T^3 - 1.06585607e-17 * T^4
    @register_symbolic cpf_f(Tf)

    #Solid Properties
    #ks = (500.0) * 1.97 * ((1 - e)^1.5)/10. #W/m.K
    ρs = 3200  #kg/m3
    ks_f(T) = (191.9216 - 0.3261784 * T^1 + 2.739462e-4 * T^2 - 7.70926e-8 * T^3) #COMSOL SiC Alpha Polycrystaline
    #ks_f(T) = 52000*exp(-1.24*10^(-5)*(T-273)) / ((T-273)+437) # Munro 1997
    @register_symbolic ks_f(Tf)

    #Cps = 1290  #J/kg*K
    #Cps(T) = (8.5e-5 * T^2 + 5.63e-2 * T - 4.05e7 * T^-2 + 1125.8) *2# * 4.4 #J/kg*K COMSOL
    Cps(T) = (1110 + 0.15 * T - 425 * exp(-0.003 * T)) *1 # Munro 1997
    @register_symbolic Cps(Ts)
    aCp = 1.0 #correction factor

    #alumina adaptor Properties
    ρs2 = 3900  #kg/m3
    ks2_f(T) = 5.5+34.5*exp(-0.0033*(T-273)) # k by Morrell 1987 (Auerkari 1996)
    @register_symbolic ks2_f(Tf)
    Cps2(T) = (1.00446+1.742e-4 *T - 2.796e4 * T^-2)*1000 #J/kg*K Cp by Touloukian 1970 (Auerkari 1996)
    @register_symbolic Cps2(Ts2)
    h_s1s2 = 100.0 #W/m2.K conduction/convective exchange with first solid

    #Insulation (felt) Properties
    ρs_ins = 140.0 # kg/m3 (COMSOL ins_rho = 0.140 g/cm3)
    Cps_ins = 1360.0 # J/kg*K (COMSOL ins_cp = 3500)
    
    #Metal Casing Properties
    ρs_metal = 2700.0 # kg/m3 (COMSOL M_rho = 2700)
    Cps_metal = 900.0 # J/kg*K (COMSOL M_Cp = 900)
    ϵ_metal = 0.2 # metal outer emissivity (COMSOL M_emis = 0.2)

    h_nat = 10.0 #W/m2.K natural convection heat transfer coefficient (COMSOL h_nat = 10)
    
    
    Re_f(qlpm, T) = G(qlpm) * Dh / μf_f(T) #use of flux instead of u*ρf
    Pr_f(T) = (cpf_f(T) * μf_f(T)) / kf_f(T)
    Nu_f6(qlpm, T) = 10^A * (Re_f(qlpm, T)^B) * (Pr_f(T)^(10^C))

    Gz(qlpm, T) = (1 / w_chnl) * Re_f(qlpm, T) * Pr_f(T) * Dh
    Nu_f5(qlpm, T) = A * (1 - B * Gz(qlpm, T)^n) * exp(-C / Gz(qlpm, T))

    #Gz_f(qlpm, x, T) = (1/x) * Re_f(qlpm) * Pr(T) * Dh
    #Nu_f5(qlpm, T) = A * [1 - B*Gz_f(qlpm, x, T)^Re_f(qlpm)^B) * (Pr(T)^C)

    h_avg_f5(qlpm, T) = (Nu_f5(qlpm, T) * kf_f(T)) / Dh
    h_avg_f6(qlpm, T) = (Nu_f6(qlpm, T) * kf_f(T)) / Dh
    #h_f = ha + hb / (w_chnl) #W/m2.K

    # h_avg_f_q(qlpm) = A * Re_f(qlpm)^B
    # h_avg_f_T(T) = (((cpf_f(T) * mu) / kf_f(T))^(C)) * kf_f(T) / Lc #9. instead of C because library error
    #@register_symbolic h_avg_f_T(x) #does not work if symbolic

    # MOL Discretization parameters for system 1
    t_min = 0.0
    t_max = 7084.0

    ## SETUP Parameters

    p_opt = [A => -1., B => 0.4, C => -5.]  #[A => 4., B => 0.06, C => 51., n => 0.5] 
    p_cond = [Io => 456000.0, qlpm => 12.5]
    p_math = Dict(vcat(p_opt, p_cond))
    #extract the parameters names from p_math
    p_names = [i for i in keys(p_math)]

    # ODE equations for 5-node system with algebraic NTU fluid node
    eq1 = [
        Vs * Dt(Ts) ~ ( # receiver monolith
            α * Io * A_frt # irradiance absorbed by the receiver
            - hext * A_frt_solid * (Ts - Tamb) # front convection losses
            - ϵ * σ * A_frt * (Ts^4 - Tamb^4) # radiation losses
            - (m(qlpm) * cpf_f((Ts + Tamb)/2) * (Tf - Tamb)) # enthalpy rise of gas (NTU model)
            - (h_s1s2 * A_s1s2 * (Ts - Ts2)) # contact conduction to alumina adaptor
            - ((Ts - Ts3) / (R_ins_total / 2)) # conduction loss to insulation felt
            ) / (ρs * Cps(Ts)),
            
        Tf ~ Ts - (Ts - Tamb) * exp(-(h_avg_f6(qlpm, (Ts + Tamb)/2) * A_exchange) / (m(qlpm) * cpf_f((Ts + Tamb)/2))), # algebraic NTU relation
            
        Vs2 * Dt(Ts2) ~ ( # alumina adaptor
            - h_s1s2 * A_s1s2 * (Ts2 - Ts) # conduction from receiver monolith
            - ((Ts2 - Ts3) / (R_ins_adpt / 2)) # conduction loss to insulation felt
            ) / (ρs2 * Cps2(Ts2)),
            
        Vs_ins * Dt(Ts3) ~ ( # insulation felt midpoint
            ((Ts - Ts3) / (R_ins_total / 2)) # heat from receiver monolith
            + ((Ts2 - Ts3) / (R_ins_adpt / 2)) # heat from adaptor
            - ((Ts3 - Ts4) / R_s3s4) # heat to metal casing
            ) / (ρs_ins * Cps_ins),
            
        Vs_metal * Dt(Ts4) ~ ( # aluminum casing outer shell
            ((Ts3 - Ts4) / R_s3s4) # heat from insulation felt
            - h_nat * A_outer_casing * (Ts4 - Tamb) # natural convection to environment
            - ϵ_metal * σ * A_outer_casing * (Ts4^4 - Tamb^4) # radiation to environment
            ) / (ρs_metal * Cps_metal)
    ]
  

    # Initial values for Ts, Ts2, Ts3, Ts4
    u0 = [Ts => 293.0, Ts2 => 295.0, Ts3 => 296.0, Ts4 => 297.0]

    # Time span for the solution
    tspan = (t_min, t_max)
    
    # ODE system for system 1
    @named odesys = ODESystem(eq1, t)
    simplified_odesys = structural_simplify(odesys)
    prob = ODEProblem(simplified_odesys, merge(Dict(u0), p_math), tspan)
end

sol_temp = solve(prob, FBDF(), reltol=1e-12, abstol=1e-12)

include("import_exp_0D.jl") #import the experimental data
colors = ColorSchemes.tab10[1:4]

begin # define functions
    #Optimization using NLOpt
    function NLmodeloptim(tvalues, rmp, tolr)
        u0_map = Dict(Ts => rmp[Tinit], Ts2 => rmp[Tinit] + 2.0, Ts3 => rmp[Tinit] + 3.0, Ts4 => rmp[Tinit] + 4.0)
        rmp_clean = filter(pair -> !isequal(pair.first, Tinit), rmp)
        modeloptim = remake(prob, u0 = u0_map, p = rmp_clean, tspan = (0.0, tvalues[end]))
        modeloptim_sol = solve(modeloptim, FBDF(), saveat=tvalues, reltol=tolr)
        tempT8_op = float.(modeloptim_sol[Ts])  # corresponds to Ts 
        tempT3_op = float.(modeloptim_sol[Tf])  # corresponds to Tf
        return ([tempT8_op tempT3_op])
    end
    function remakeAysha(pguess_l, cond_k, time_opt, Tinit_val; tolr=1e-7)
        # combine the pguess_l and cond into a single Dict
        pguess_temp = Dict(k => pguess_l[i] for (i, (k, v)) in enumerate(p_opt))
        rmp = Dict(pguess_temp..., cond_k...)
        rmp[Tinit] = Tinit_val
        #display(rmp)
        local temp_T = NLmodeloptim(time_opt, rmp, tolr)
        return temp_T
    end

    function lossAll(pguess_l, sim_key)
        lossr = 0.0
        for sm in sim_key
            local cond_k = simulation_conditions[sm]
            local expdata = reduce(hcat, measurements[(measurements.simulation_id .== sm), :temperatures])[:, 2:3]
            local time_exp = measurements[measurements.simulation_id .== sm, :time][1]
            local Tinit_val = expdata[1, 1]
            # Run model and get temperatures
            temp_T = remakeAysha(pguess_l, cond_k, time_exp, Tinit_val)[:, 1:2]
            
            if length(expdata[:, 1]) != length(temp_T[:, 1])
                return Inf
            end
            
            N_samples = length(expdata[:, 1])
            
            # Calculate total experimental temperature rise for normalization
            delta_Ts_exp = expdata[end, 1] - expdata[1, 1]
            delta_Tf_exp = expdata[end, 2] - expdata[1, 2]
            
            # Prevent division by zero
            delta_Ts_exp = abs(delta_Ts_exp) > 1e-3 ? delta_Ts_exp : 1.0
            delta_Tf_exp = abs(delta_Tf_exp) > 1e-3 ? delta_Tf_exp : 1.0
            
            # Normalized sum of squared errors
            error_s = sum(((temp_T[:, 1] .- expdata[:, 1]) .^ 2)) / (delta_Ts_exp^2)
            error_f = sum(((temp_T[:, 2] .- expdata[:, 2]) .^ 2)) / (delta_Tf_exp^2)
            
            # Average normalized error for this experiment
            temp_error = (error_s + error_f) / N_samples
            lossr += temp_error
        end
        return lossr / length(sim_key)
    end    
end
begin #Optimization
    #p_opt = [ha => 8.99, hb => 1.0, cps_c => 6.0] 
    p_opt = [A => -2., B => 1.61, C => -2.]  #[A => 4., B => 0.06, C => 51., n => 0.5]

    p0 = [x[2] for x in p_opt]
    optf = OptimizationFunction(lossAll, Optimization.AutoForwardDiff())
    lb = [-9., 0.01, -5.]
    ub = [4., 10., 2.] 

    sim_key = ["E67","E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"]
    #sim_key = ["E70"]#, "E74"]#["E67"]#, "E78", "E79", "E80", "E81"]
    initialerror = (lossAll(p0, sim_key))
    println("Initial Error $initialerror")
    optprob = Optimization.OptimizationProblem(optf, p0, sim_key, lb=lb, ub=ub)
    optsol = solve(optprob, NLopt.GN_MLSL_LDS(), local_method=NLopt.LN_NELDERMEAD(), maxtime=60, local_maxiters=10000)
    
    println(optsol.retcode)
    pnew = optsol.u
    display(pnew)
    #res_error = lossAll(pnew, [])
    println("Final Error $(optsol.objective)")
end



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

    println("\nPlotting transient profiles for all experiments...")
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
    println("\n=== Computing Additional Metrics (0D v3 NTU) ===")
    
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
