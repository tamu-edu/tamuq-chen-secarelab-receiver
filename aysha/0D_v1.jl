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
    using CSV

end

begin #FIXED Parameters

    Tamb = (22.448 + 273.15) #K (same for all exp) 

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
    A_frt_solid = A_frt - A_chnl_frt_all  # Frontal area of just solid (m^2)
    A_exchange = A_chnl_p * n_chnl  # Total contact area between fluid and solid (m^2)

    # Volume calculations
    Vs = A_frt_solid * L  # Corrected solid volume (m^3)
    Vf = n_chnl * A_chnl_frt * L  # Fluid volume in the channels (m^3)

    # Alumina Adaptor dimensions
    rs2 = 38.8E-3/2 # radius of the second solid
    Ls2 = 59.E-3 # length of the second solid
    Vs2 = π * rs2^2 * Ls2 # full volume of the second solid (m^3)
        - A_frt * Ls2 # minus the volume of the first solid in contact with the solid (omitted the tube)
    A_s1s2 = A_frt * Ls2 # area of the second solid in contact with the first solid (m^2)
    A_s_p2 = 2 * π * rs2 * Ls2 #periphery area of the second solid (m^2)

       #Constants
    σ = 5.670374419e-8 # W/m2.K^4 Stefan-Boltzmann constant (COMSOL value)
    hext = 10.0 #W/m2.K convective heat transfer coefficient for the front 

    #SiC Properties
    ϵ = 0.85 #emissivity
    α = 0.85 #absroptivity
    e = 0.425 #porosity

    #parameters for the insulation + metal (3rd solid)
    r0 = 23 /2 / 1000 #m
    r_ins = 150/2 / 1000 #m alumina felt
    r_metal = 180/2/1000 #m aluminum
    Vs3 = π * r_metal^2 * L # full volume of the third solid (m^3)
        + 2 * π * r_metal * (r_metal - r_ins) / 2  # plus the back volume of the metal
    Vs_ins = π * (r_ins^2 - r0^2) * L # full volume of the insulation (m^3)

    A_s1s3 = A_frt * (L- 0.029) # area of the third solid in contact with the first solid (m^2)
    A_s2s3 = A_s_p2 # area of the third solid in contact with the second solid (m^2)
    A_s_p3 = 2 * π * r_metal * L #periphery area of the third solid (m^2)
    A_s_b3 = π * r_metal^2 #back area of the third solid (m^2)
    

    kins = 0.078 #W/m*K -->assume very few losses through the insulation
    kmet = 201 #W/m*K 

    R1 = log(r_ins / r0) / (2 * π * kins * L) #thermal resistance of the insulation
    R2 = log(r_metal / r_ins) / (2 * π * kmet * L) #thermal resistance of the metal
    k_eff = log(r_metal / r0) / (2* π * L *(R1 + R2)) #effective thermal conductivity of the insulation

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
    @independent_variables t 
    @parameters kfs_c cps_c I_loss A B C n #ks h_average A n
    @parameters Io qlpm Tinit
    @variables Ts(t) Tf(t) Ts2(t) Ts3(t)
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

    #Insulation + metal Properties
    ρs3 = (160 * Vs_ins + 2700 * (Vs3 - Vs_ins)) / Vs3  #kg/m3
    Cps3 = (1360 * Vs_ins + 900 * (Vs3 - Vs_ins)) / Vs3 #J/kg*K
    h_s1s3 = h_s1s2 #W/m2.K conduction/convective exchange with first solid
    h_s2s3 = h_s1s2 #W/m2.K conduction/convective exchange with second solid
    h_nat = 10.0 #W/m2.K natural convection heat transfer coefficient
    
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
    # h_avg_f_T(T) = (((cpf_f(T) * mu) / kf_f(T))^(C)) * kf_f(T) / Dh #9. instead of C because library error
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

    # PDE equation for system 1
    eq1 = [
        Vs * Dt(Ts) ~ ( # receiver
            α * Io * A_frt
            - hext * A_frt * (Ts - Tamb)
            - ϵ * σ * A_frt * (Ts^4 - Tamb^4)
            - (h_avg_f6(qlpm, Tf) * A_exchange * (Ts - Tf))
            - (h_s1s2 * A_s1s2 * (Ts - Ts2)) # conduction/convective exchange with second solid
            - (h_s1s3 * A_s1s3 * (Ts - Ts3)) # conduction/convective exchange with third solid 
            ) / (ρs * 1. * Cps(Ts)),
        Vf * Dt(Tf) ~ ( # fluid
             - m(qlpm) * cpf_f(Tf) * (Tf - Tamb)
             + (h_avg_f6(qlpm, Tf)* A_exchange * (Ts - Tf))
            ) / (ρf_f(Tf) * cpf_f(Tf)),
        Vs2 * Dt(Ts2) ~ ( # adaptor
            - h_s1s2 * A_s1s2 * (Ts2 - Ts)                     # conduction/convective exchange with first solid
            - h_s2s3 * A_s2s3 * (Ts2 - Ts3)                 # conduction/convective exchange with third solid
            #- h_avg_f6(qlpm, Tf) * A_s2f * (Ts2 - Tf)          # heat exchange with fluid (if applicable)
            ) / (ρs2 * Cps2(Ts2)),
        Vs3 * Dt(Ts3) ~ ( # insulation + metal
            - h_s1s3 * A_s1s3 * (Ts3 - Ts)                     # conduction/convective exchange with first solid
            - h_s2s3 * A_s2s3 * (Ts3 - Ts2)                 # conduction/convective exchange with second solid
            - h_nat * (A_s_p3 + 2 * A_s_b3) * (Ts3 - (Tamb+50))
            - ϵ * σ * (A_s_p3 + 2 * A_s_b3)* (Ts3^4 - Tamb^4)  # insulation losses periphery, back and front areas
            ) / (ρs3 * Cps3)
    ]
  

    # Initial values for Ts and Tf
    u0 = [Ts => 293.0, Tf => 294.0, Ts2 => 295.0, Ts3 => 296.0]

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
        u0_map = Dict(Ts => rmp[Tinit], Tf => rmp[Tinit] + 1.0, Ts2 => rmp[Tinit] + 2.0, Ts3 => rmp[Tinit] + 3.0)
        rmp_clean = filter(pair -> !isequal(pair.first, Tinit), rmp)
        modeloptim = remake(prob, u0 = u0_map, p = rmp_clean, tspan = (0.0, tvalues[end]))
        modeloptim_sol = solve(modeloptim, FBDF(), saveat=tvalues, reltol=tolr)
        tempT8_op = float.(modeloptim_sol[Ts])  # corresponds to Ts 
        tempT3_op = float.(modeloptim_sol[Tf])  # corresponds to Tf
        return ([tempT8_op tempT3_op])
    end
    function remakeAysha(pguess_l, cond_k, time_opt; tolr=1e-7)
        # combine the pguess_l and cond into a single Dict
        pguess_temp = Dict(k => pguess_l[i] for (i, (k, v)) in enumerate(p_opt))
        rmp = Dict(pguess_temp..., cond_k...)
        #display(rmp)
        local temp_T = NLmodeloptim(time_opt, rmp, tolr)
        return temp_T
    end

    function lossAll(pguess_l, sim_key)
        lossr = 0
        for sm in sim_key
            local cond_k = simulation_conditions[sm]
            local expdata = reduce(hcat, measurements[(measurements.simulation_id .== sm), :temperatures])[:, 1:2]
            local time_exp = measurements[measurements.simulation_id .== sm, :time][1]
            # Run model and get temperatures
            temp_T = remakeAysha(pguess_l, cond_k, time_exp)[:, 1:2]
            #expdata_last = expdata[end, :]   # Last row of experimental data
            #temp_T_last = temp_T[end, :]     # Last row of simulated temperature data
            # Compute error
            if length(expdata) != length(temp_T)
                return temp_error = Inf
            end
            temp_error = sqrt(sum((temp_T .- expdata) .^ 2)) #/ length(expdata)
            #temp_error = sqrt(sum((temp_T_last .- expdata_last) .^ 2))
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
    #gas 
    T_steady = DataFrame(sim_id=[], time=[], T_mod=[], T_exp=[])
    #plotting the scatter plots for the different thermocouples for the steady state temperatures
    ordered_sim_conditions = sim_key#["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"]
    
    for sm in ordered_sim_conditions
        # Retrieve from measurements the experimental data for the current simulation condition
        #println(sm)
        local cond_k = simulation_conditions[sm]
        local sel_meas = measurements[(measurements.simulation_id.==sm).&(measurements.obs_id.=="_Tf"), :]
        local expdata = sel_meas[:, :temperatures]
        local time_opt = sel_meas[:, :time][1]

        #run selected simulation and get the steady temperature values
        local temp_T = remakeAysha(pnew, cond_k, time_opt; tolr=1e-7)
        push!(T_steady, (sm, time_opt, temp_T[:,2], expdata[1]))
    end
    color_model = :blue
    color_exp = :orange

    # Extract T_mod and T_exp values
    T_mod_values = [T_steady[i, :T_mod][end] for i in 1:length(ordered_sim_conditions)]
    T_exp_values = [T_steady[i, :T_exp][end] for i in 1:length(ordered_sim_conditions)]

    # Create the bar plot
    bar_data = hcat(T_mod_values, T_exp_values)
    plot0 = groupedbar(bar_data, label=["T_mod" "T_exp"], legend=:topright, title="T_steady Comparison",
        xlabel="Experimental Runs", ylabel="Temperature (K)",
        bar_width=0.4, xticks=(1:length(ordered_sim_conditions), ordered_sim_conditions),
        bar_position=:dodge, color=[color_model color_exp], ylimit=(0, 1250))

    lims = (300, 800)
    plot1 = scatter(T_mod_values, T_exp_values, label="", title="Gas SS Temp",
        xlabel="T_mod", ylabel="T_exp", ylims=lims, xlims=lims, aspect_ratio=:equal)
    plot1 = plot!(collect(lims), collect(lims))

    #plot(plot0, plot1)
    
    #Plotting the temperature profiles for the gas domanin and a few solid domain profiles across all experimental conditions

    function plt_own(it)
        plot2 = plot(
            T_steady[it, :time],
            T_steady[it, :T_mod],
            title="Gas Temperature Profile T3",
            label=permutedims(ordered_sim_conditions[it]),
            xlabel="Time (s)",
            ylabel="Temperature (K)", xlimit=(0, 4050),
            legend=:bottomright, color_palette=colors)
        scatter!(T_steady[it, :time],
            T_steady[it, :T_exp], label=permutedims(ordered_sim_conditions[it]),
            color_palette=colors
        )

        return plot2
    end
    function plt_case(sm, params)
        local cond_k = simulation_conditions[sm]
        local expdata = reduce(hcat, measurements[(measurements.simulation_id.==sm), :temperatures])[:, 1:2] 
        local time_opt = measurements[measurements.simulation_id.==sm, :time][1]
        local temp_T = remakeAysha(params, cond_k, time_opt)[:, 1:2] 

        labels = measurements[(measurements.simulation_id.==sm), :obs_id]
        plt_t = plot(title="Temperature Profiles $sm", xlabel="Time (s)", ylabel="Temperature (K)", ylim=(300, 1000),
            legend=:outerright, color_palette=colors)
        scatter!(time_opt, expdata, label=permutedims(labels),
            color_palette=colors)
        plot!(time_opt, temp_T, label=permutedims(labels), lw=3)
        return plt_t
    end
    #display(plot(plot1, plt(1:length(sim_key)), layout=(2, 1), size=(900, 600)))
    display(plt_own(1:length(sim_key)))
    map(x -> display(plt_case(x , pnew)), sim_key)
end

begin 
    #solid
    T_steady_solid = DataFrame(sim_id=[], time=[], T_mod=[], T_exp=[])
    #plotting the scatter plots for the different thermocouples for the steady state temperatures
    ordered_sim_conditions = sim_key#["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"]
    
    for sm in ordered_sim_conditions
        # Retrieve from measurements the experimental data for the current simulation condition
        #println(sm)
        local cond_k = simulation_conditions[sm]
        local sel_meas = measurements[(measurements.simulation_id.==sm).&(measurements.obs_id.=="_Tavg"), :]
        local expdata = sel_meas[:, :temperatures]
        local time_opt = sel_meas[:, :time][1]

        #run selected simulation and get the steady temperature values
        local temp_T = remakeAysha(pnew, cond_k, time_opt; tolr=1e-7)
        push!(T_steady_solid, (sm, time_opt, temp_T[:,1], expdata[1]))
    end
    color_model = :blue
    color_exp = :orange

    # Extract T_mod and T_exp values
    T_mod_values_solid = [T_steady_solid[i, :T_mod][end] for i in 1:length(ordered_sim_conditions)]
    T_exp_values_solid = [T_steady_solid[i, :T_exp][end] for i in 1:length(ordered_sim_conditions)]

    # Create the bar plot
    bar_data = hcat(T_mod_values_solid, T_exp_values_solid)
    plot3 = groupedbar(bar_data, label=["T_mod" "T_exp"], legend=:topright, title="T_steady Comparison",
        xlabel="Experimental Runs", ylabel="Temperature (K)",
        bar_width=0.4, xticks=(1:length(ordered_sim_conditions), ordered_sim_conditions),
        bar_position=:dodge, color=[color_model color_exp], ylimit=(0, 1250))

    lims = (300, 1150)
    plot4 = scatter(title="Solid SS Temp",T_mod_values_solid, T_exp_values_solid,
    xlabel="T_mod", ylabel="T_exp", ylims=lims, xlims=lims, aspect_ratio=:equal)
    plot4 = plot!(collect(lims), collect(lims)) 
end

#plot(plot1, plot4)
display(plot(plot1, plot4, layout=(2, 1), size=(600, 900)))

A, B, C = pnew

T_test=range(300, 1000,100)
plot_flow = plot(xlabel="T (K)", ylabel="h_avg (W/m2.K)", title="h vs T", legend=:topright, color_palette=colors)
flow_rates= [1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0]
for fl in flow_rates
    plot!(T_test, h_avg_f6.(fl, T_test), label="$(fl) Lpm")
end
display(plot_flow)

    
fl_test=range(1, 30, 30)
plot(fl_test, h_avg_f6.(fl_test, 300), label="T=300 K", xlabel="qlpm (Lpm)", 
    ylabel="h_avg (W/m2.K)", title="h vs qlpm", legend=:bottomright, color_palette=colors,
    xlims=(0, 30))
plot!(fl_test, h_avg_f6.(fl_test, 400), label="T=400 K")
plot!(fl_test, h_avg_f6.(fl_test, 500), label="T=500 K")
plot!(fl_test, h_avg_f6.(fl_test, 600), label="T=600 K")
plot!(fl_test, h_avg_f6.(fl_test, 700), label="T=700 K")
plot!(fl_test, h_avg_f6.(fl_test, 800), label="T=800 K")
plot!(fl_test, h_avg_f6.(fl_test, 900), label="T=900 K")
plot!(fl_test, h_avg_f6.(fl_test, 1000), label="T=1000 K")

# ==========================================
# ADDITIONAL METRICS COMPUTATION AND EXPORT
# ==========================================
begin
    println("\n=== Computing Additional Metrics (0D v1) ===")
    
    # Helper function for t90
    function get_t90(time, temp)
        t_init = temp[1]
        t_ss = temp[end]
        target = t_init + 0.9 * (t_ss - t_init)
        idx = findfirst(t -> t >= target, temp)
        if idx !== nothing
            return time[idx]
        else
            return time[end]
        end
    end

    # Define the 0D analysis function to extract all nodes
    function run_analysis_0D_v1(pguess_l, cond_k, time_opt; tolr=1e-7)
        pguess_temp = Dict(k => pguess_l[i] for (i, (k, v)) in enumerate(p_opt))
        rmp = Dict(pguess_temp..., cond_k...)
        u0_map = Dict(Ts => rmp[Tinit], Tf => rmp[Tinit] + 1.0, Ts2 => rmp[Tinit] + 2.0, Ts3 => rmp[Tinit] + 3.0)
        rmp_clean = filter(pair -> !isequal(pair.first, Tinit), rmp)
        modeloptim = remake(prob, u0 = u0_map, p = rmp_clean, tspan = (0.0, time_opt[end]))
        sol = solve(modeloptim, FBDF(), saveat=time_opt, reltol=tolr)
        
        Ts_sol = float.(sol[Ts])
        Tf_sol = float.(sol[Tf])
        Ts3_sol = float.(sol[Ts3])
        
        return Ts_sol, Tf_sol, Ts3_sol
    end

    # Determine the next RunID
    csv_path = "D:/kkakosim/sim_comsol/analysis_results_0D_v1.csv"
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
        RunID = Float64[],
        Case = String[],
        T9_SS_sim = Float64[],
        T9_SS_exp = Float64[],
        dT_T09 = Float64[],
        T3_SS_sim = Float64[],
        T3_SS_exp = Float64[],
        dT_T03 = Float64[],
        T2_SS_sim = Float64[],
        T2_SS_exp = Float64[],
        dT_T02 = Float64[],
        R_leak_sim = Float64[],
        R_leak_exp = Float64[],
        t90_sim_T09_s = Float64[],
        t90_exp_T09_s = Float64[],
        dt90_T09_s = Float64[],
        t90_sim_T03_s = Float64[],
        t90_exp_T03_s = Float64[],
        dt90_T03_s = Float64[],
        Gap_ss_sim = Float64[],
        Gap_ss_exp = Float64[],
        dGap_ss = Float64[]
    )

    for sm in sim_key
        local cond_k = simulation_conditions[sm]
        
        # Extract experimental profiles
        local time_exp = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tavg"), :time][1]
        local T9_exp = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_T9"), :temperatures][1]
        local T3_exp = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_Tf"), :temperatures][1] # Tf is T3
        local T2_exp = measurements[(measurements.simulation_id .== sm) .& (measurements.obs_id .== "_T2"), :temperatures][1]
        
        # Run simulation
        Ts_sim, Tf_sim, Ts3_sim = run_analysis_0D_v1(pnew, cond_k, time_exp)
        
        # Steady-state values (last element)
        Ts_ss = Ts_sim[end]
        Tf_ss = Tf_sim[end]
        Ts3_ss = Ts3_sim[end]
        
        T9_exp_ss = T9_exp[end]
        T3_exp_ss = T3_exp[end]
        T2_exp_ss = T2_exp[end]
        
        # Temperature differences
        dT_T09 = Ts_ss - T9_exp_ss
        dT_T03 = Tf_ss - T3_exp_ss
        dT_T02 = Ts3_ss - T2_exp_ss
        
        # Energy partitioning ratio (R_leak)
        R_leak_sim = (Tf_ss - Tamb) / (Ts_ss - Tamb)
        R_leak_exp = (T3_exp_ss - Tamb) / (T9_exp_ss - Tamb)
        
        # t90 times
        t90_sim_T09 = get_t90(time_exp, Ts_sim)
        t90_exp_T09 = get_t90(time_exp, T9_exp)
        dt90_T09 = t90_sim_T09 - t90_exp_T09
        
        t90_sim_T03 = get_t90(time_exp, Tf_sim)
        t90_exp_T03 = get_t90(time_exp, T3_exp)
        dt90_T03 = t90_sim_T03 - t90_exp_T03
        
        # Gas-to-solid gap
        gap_sim = Tf_ss - Ts_ss
        gap_exp = T3_exp_ss - T9_exp_ss
        dgap = gap_sim - gap_exp
        
        push!(results_df, (
            next_run_id, sm, Ts_ss, T9_exp_ss, dT_T09, Tf_ss, T3_exp_ss, dT_T03, Ts3_ss, T2_exp_ss, dT_T02,
            R_leak_sim, R_leak_exp, t90_sim_T09, t90_exp_T09, dt90_T09,
            t90_sim_T03, t90_exp_T03, dt90_T03, gap_sim, gap_exp, dgap
        ))
    end
    
    # Print the table to console
    println(results_df)
    
    # Save the dataframe to a CSV file (appending if file exists)
    if file_exists
        CSV.write(csv_path, results_df, append=true)
        println("Appended additional metrics to $csv_path with RunID = $next_run_id")
    else
        CSV.write(csv_path, results_df)
        println("Created new CSV file $csv_path and saved metrics with RunID = $next_run_id")
    end
end
