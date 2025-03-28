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

    # Receiver dimensions
    w_t = 19.e-3  # Width of the receiver in meters
    A_frt = w_t * w_t  # Total frontal area of the receiver (m^2)
    L = 137e-3 #m Length
    A_s_p = w_t * L * 4 #total area solid periphery m2

    # Channel dimensions
    w_chnl = 1.5e-3  # Width of a single channel (m)
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

    #Constants
    σ = 5.17e-8 #W/m2.K^4 Stefan-Boltzmann constant
    hext = 1.0 #W/m2.K convective heat transfer coefficient for the front 

    #SiC Properties
    ϵ = 0.75 #emissivity
    α = 0.625 #absroptivity
    e = 0.425 #porosity

    #parameters for the insulation
    r0 = 23 / 1000 #m
    r_ins = 42 / 1000 #m
    kins = 1. * 0.078 #W/m*K -->assume very few losses through the insulation
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
    @variables t Ts(t) Tf(t)
    Dt = Differential(t)
    
    ### Fluid Properties
    RH = 0.4 #relative humidity
    @register_symbolic ρf_f(Tf)
    ρf_f(T) = 352.716 * T^-1 * (1+RH) /(1+ 1.609 * RH) #COMSOL + https://www.engineeringtoolbox.com/density-air-d_680.html with some error for w vs RH
  
    ρf = ρf_f(1000.0) #kg/m3
    m(qlpm) = ρf * (qlpm / 60.0 / 1000.0) #mass at the flowmeter conditions kg/s

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
    @register_symbolic ks_f(Tf)
    #Cps = 1290  #J/kg*K
    Cps(T) = (8.5e-5 * T^2 + 5.63e-2 * T - 4.05e7 * T^-2 + 1125.8) #.* 4. #J/kg*K
    @register_symbolic Cps(Ts)
    aCp = 1.0 #correction factor


    Re_f(qlpm, T) = (m(qlpm) / A_chnl_frt_all) * Lc / μf_f(T) #use of flux instead of u*ρf
    Pr(T) = (cpf_f(T) * μf_f(T)) / kf_f(T)
    Nu_f6(qlpm, T) = A * (Re_f(qlpm, T)^B) * (Pr(T)^C)

    Gz(qlpm, T) = (1 / w_chnl) * Re_f(qlpm, T) * Pr(T) * Lc
    Nu_f5(qlpm, T) = A * (1 - B * Gz(qlpm, T)^n) * exp(-C / Gz(qlpm, T))

    #Gz_f(qlpm, x, T) = (1/x) * Re_f(qlpm) * Pr(T) * Lc
    #Nu_f5(qlpm, T) = A * [1 - B*Gz_f(qlpm, x, T)^Re_f(qlpm)^B) * (Pr(T)^C)

    h_avg_f5(qlpm, T) = (Nu_f5(qlpm, T) * kf_f(T)) / Lc
    h_avg_f6(qlpm, T) = (Nu_f6(qlpm, T) * kf_f(T)) / Lc
    #h_f = ha + hb / (w_chnl) #W/m2.K

    # h_avg_f_q(qlpm) = A * Re_f(qlpm)^B
    # h_avg_f_T(T) = (((cpf_f(T) * mu) / kf_f(T))^(C)) * kf_f(T) / Lc #9. instead of C because library error
    #@register_symbolic h_avg_f_T(x) #does not work if symbolic

    # MOL Discretization parameters for system 1
    t_min = 0.0
    t_max = 7084.0

    ## SETUP Parameters

    p_opt = [A => 100., B => 0.2, C => 0.6]  #[A => 4., B => 0.06, C => 51., n => 0.5] 
    p_cond = [Io => 456000.0, qlpm => 9.1, Tinit => 293.0]
    p_math = Dict(vcat(p_opt, p_cond))
    #extract the parameters names from p_math
    p_names = [i for i in keys(p_math)]

    # PDE equation for system 1
    eq1 = [
        Vs * Dt(Ts) ~ (
            α * Io * A_frt
            - hext * A_frt * (Ts - Tamb)
            - ϵ * σ * A_frt * (Ts^4 - Tamb^4)
            - (h_avg_f6(qlpm, Tf) * A_exchange * (Ts - Tf)) 
            - (kins * (r_ins / r0) * (Ts - Tins_f(t)) * A_s_p / (r_ins - r0))
        ) / (ρs * 1. * Cps(Ts)),
        Vf * Dt(Tf) ~ (
             - m(qlpm) * cpf_f(Tf) * (Tf - Tamb)
             + (h_avg_f6(qlpm, Tf)* A_exchange * (Ts - Tf))
        ) / (ρf_f(Tf) * cpf_f(Tf))
    ]
  

    # Initial values for Ts and Tf
    u0 = [Ts => Tinit, Tf => Tinit]

    # Time span for the solution
    tspan = (t_min, t_max)
    
    # ODE system for system 1
    @named odesys = ODESystem(eq1, t, [Ts, Tf], p_names; defaults=p_math)
    simplified_odesys = structural_simplify(odesys)
    prob = ODEProblem(simplified_odesys, u0, tspan)
end

include("import_exp_0D.jl") #import the experimental data
colors = ColorSchemes.tab10[1:4]

begin # define functions
    #Optimization using NLOpt
    function NLmodeloptim(tvalues, rmp, tolr)
        modeloptim = remake(prob, u0=ones(length(prob.u0)) * rmp[Tinit], tspan=(0.0, tvalues[end]))
        replace!(Tunable(), modeloptim.p, collect(values(rmp)))
        modeloptim_sol = solve(modeloptim, FBDF(), saveat=tvalues, reltol=tolr)#, reltol=1e-12, abstol=1e-12)
        #time = modelfit_sol.t
        tempT8_op = float.(modeloptim_sol'[:, 1])  # Index 1 corresponds to Ts 
        tempT3_op = float.(modeloptim_sol'[:, 2])  # Index 2 corresponds to Tf
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
            temp_error = sqrt(sum((temp_T .- expdata) .^ 2)) / length(expdata)
            #temp_error = sqrt(sum((temp_T_last .- expdata_last) .^ 2))
            lossr += temp_error
        end
        return lossr / length(sim_key)
    end    
end
begin #Optimization
    #p_opt = [ha => 8.99, hb => 1.0, cps_c => 6.0] 
    p_opt = [A => 9., B => 0.4, C => 10.]  #[A => 4., B => 0.06, C => 51., n => 0.5]

    p0 = [x[2] for x in p_opt]
    optf = OptimizationFunction(lossAll, Optimization.AutoForwardDiff())
    lb = [0.001, 0.001, 0.0001]
    ub = [10., 10., 20.] 

    sim_key = ["E67","E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"]
    #sim_key = ["E70"]#, "E74"]#["E67"]#, "E78", "E79", "E80", "E81"]
    initialerror = (lossAll(p0, sim_key))
    println("Initial Error $initialerror")
    optprob = Optimization.OptimizationProblem(optf, p0, sim_key, lb=lb, ub=ub)
    optsol = solve(optprob, NLopt.GN_MLSL_LDS(), local_method=NLopt.LN_NELDERMEAD(), maxtime=60, local_maxiters=100000)
    
    println(optsol.retcode)
    pnew = optsol.u
    println(pnew)
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

    lims = (500, 900)
    plot1 = scatter(T_mod_values, T_exp_values, label="",
        xlabel="T_mod", ylabel="T_exp", ylims=lims, xlims=lims, aspect_ratio=:equal)
    plot1 = plot!(collect(lims), collect(lims))

    #plot(plot0, plot1)
    
    #Plotting the temperature profiles for the gas domanin and a few solid domain profiles across all experimental conditions

    function plt(it)
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
        plt = plot(title="Temperature Profiles $sm", xlabel="Time (s)", ylabel="Temperature (K)", ylim=(300, 1600),
            legend=:outerright, color_palette=colors)
        scatter!(time_opt, expdata, label=permutedims(labels),
            color_palette=colors)
        plot!(time_opt, temp_T, label=permutedims(labels), lw=3)
        return plt
    end
    display(plot(plot1, plt(1:length(sim_key)), layout=(2, 1), size=(900, 600)))
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

    lims = (500, 1250)
    plot4 = scatter(title="Solid Steady State Temperature",T_mod_values_solid, T_exp_values_solid,
    xlabel="T_mod", ylabel="T_exp", ylims=lims, xlims=lims, aspect_ratio=:equal)
    plot4 = plot!(collect(lims), collect(lims)) 
end

plot(plot1, plot4)

