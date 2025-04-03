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

    #α = keff/ρ*Cp #kW/m2.K 
    Tamb = (22.448 + 273.15) #K (same for all exp) 
    #deltax = 0.002795918367346939 #discretization (m)
    # w_t = 19 / 1000 #m
    # A_t = w_t * w_t #m2 - for the whole receiver (19x19mm2)
    # #A_st = 176e-6 #total area - solid m2
    # #A_ft = 49e-6 # total area - fluid m2
    # w_t = 19.e-3 #m
    # A_frt = w_t * w_t #m2 - for the whole receiver (19x19mm2)
    # w_chnl = 1.5e-3 #m
    # A_chnl_p = w_chnl * L * 4  #m2 channel periphery
    # A_chnl_frt = w_chnl * w_chnl #m2 channel front
    # n_chnl = 10 * 10
    # A_exchange = A_chnl_p * n_chnl #m2 - contact area between fluid and solid
    # A_chnl_frt_all = A_chnl_frt * n_chnl #m2 all frontal area of channels
    # A_frt_solid = A_frt - A_chnl_frt_all #m2 - frontal area of solid
    # Vs = A_frt * L #m3
    # Vf = n_chnl * A_chnl_frt * L

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


    # #V= q/(A_channel*n_channel) #m/s - using the area of the whole receiver (18x18mm2)
    # V = 0.57 #m/s (calculated from excel sheet and COMSOL)
    # th_s = 0.4e-3 #m
    # th_f = 0.7e-3 #m
    #Qv= I0*exp(-1000*x)  #kW/m2 - (K extinction coefficient taken from Howell and Hendrick paper pg.86, measure pore diameter and fix)

    #Af = n_chnl * A_chnl_frt #m2
    # Lc = 4 * (w_t * w_t) / (4 * w_t)
    # w = 1.5e-3 #width of channel (m)
    # Vi = w * w * n_chnl * L #m3
    # Av = 4 * (w * L) / (w^2 * L) #specific area (m-1)

    #Constants
    σ = 5.17e-8 #W/m2.K^4 Stefan-Boltzmann constant
    hext = 1.0 #W/m2.K convective heat transfer coefficient for the front 

    #SiC Properties
    ϵ = 0.75 #emissivity
    α = 0.625 #absroptivity
    #e = 0.425 #porosity

    #parameters for the insulation
    r0 = 23 / 1000 #m
    r_ins = 42 / 1000 #m
    kins = 1. * 0.078 #W/m*K -->assume very few losses through the insulation
end;
# #for interpolations 
# #1. extract T2 data

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
# begin
#     x11 = 0.0001:0.001383838383838384:0.137 #T2 (insulation)
#     Re_f(qlpm) = (ρf * (qlpm / 60 / 1000 / Af) * w_t) / mu
#     Re = (ρ * V * w) / mu
#     Gz = (1 ./ x11) * Re_f(qlpm) * Pr * w
#     #2.create interpolation function
#     Gz_ = LinearInterpolation(x11, Gz)
#     Gz_f(x) = Gz_(x)
#     @register_symbolic Gz_f(x)
#end

begin
    # Parameters, variables, and derivatives for system 1
    @independent_variables t x
    @parameters ha hb cps_c kfs_c hext_c A B C n #ks h_average A n
    @parameters Io qlpm Tinit
    @variables t Ts(..) Tf(..)
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2

    ### Fluid Properties
    RH = 0.4 #relative humidity

    ρf_f(T) = 352.716 * T^-1 * (1+RH) /(1+ 1.609 * RH) #COMSOL + https://www.engineeringtoolbox.com/density-air-d_680.html with some error for w vs RH
    @register_symbolic ρf_f(x)
    ρf = ρf_f(1000.0) #kg/m3
    m(qlpm) = ρf * (qlpm / 60.0 / 1000.0) #mass at the flowmeter conditions kg/s

    mu = 2.0921e-5 #Pa.s
    kf = 0.056 #W/m.K
    kf_f(T) = -0.00227583562 + 1.15480022e-4 * T^1 - 7.90252856e-8 * T^2 + 4.11702505e-11 * T^3 - 7.43864331e-15 * T^4 #COMSOL
    @register_symbolic kf_f(x)
    Cpf = 1090 #J/kg.K
    Cpw_f(T) = 1745.96354 + 0.185114553 * T^1 + 6.19448731e-4 * T^2 - 3.0267851e-7 * T^3 + 4.19053122e-11 * T^4 #COMSOL
    cpf_f(T) = (1.93e-10 * T^4 - 8.0e-7 * T^3 + 1.14e-3 * T^2 - 4.49e-1 * T + 1.06e3) + RH * Cpw_f(T)
    μf_f(T) = -8.38278E-7 + 8.35717342e-8 * T^1 - 7.69429583e-11 * T^2 + 4.6437266e-14 * T^3 - 1.06585607e-17 * T^4
    @register_symbolic cpf_f(x)

    #Solid Properties
    #ks = (500.0) * 1.97 * ((1 - e)^1.5)/10. #W/m.K
    ρs = 3200  #kg/m3
    ks_f(T) = (191.9216 - 0.3261784 * T^1 + 2.739462e-4 * T^2 - 7.70926e-8 * T^3) #COMSOL SiC Alpha Polycrystaline
    @register_symbolic ks_f(X)
    #Cps = 1290  #J/kg*K
    Cps(T) = (8.5e-5 * T^2 + 5.63e-2 * T - 4.05e7 * T^-2 + 1125.8) #J/kg*K
    @register_symbolic Cps(x)
    aCp = 1.0 #correction factor

    #Nu = A*(1+(B*((Gz_f(x))^n)*exp(-C/Gz_f(x))))
    #Nu = A*(Re^B)*(Pr^C)
    #Cps(Ts) = (0.27+0.135e-4*(Ts)-9720*((Ts)^-2)+0.204e-7*((Ts)^2))/1000 #kJ/kg*K from manufacturer data

    Re_f(qlpm, T) = (m(qlpm) / A_chnl_frt_all) * Lc / μf_f(T) #use of flux instead of u*ρf
    Pr(T) = (cpf_f(T) * μf_f(T)) / kf_f(T)
    Nu_f6(qlpm, T) = A * (Re_f(qlpm, T)^B) * (Pr(T)^C)

    Gz(qlpm, x, T) = (1 / x) * Re_f(qlpm, T) * Pr(T) * Lc
    Nu_f5(qlpm, x, T) = A * (1 - B * Gz(qlpm, x, T)^n) * exp(-C / Gz(qlpm, x, T))

    #Gz_f(qlpm, x, T) = (1/x) * Re_f(qlpm) * Pr(T) * Lc
    #Nu_f5(qlpm, T) = A * [1 - B*Gz_f(qlpm, x, T)^Re_f(qlpm)^B) * (Pr(T)^C)

    h_avg_f5(qlpm, T) = (Nu_f5(qlpm, x, T) * kf_f(T)) / Lc
    h_avg_f6(qlpm, T) = (Nu_f6(qlpm, T) * kf_f(T)) / Lc
    h_local(x) = ha + hb / (1+x)

    # h_avg_f_q(qlpm) = A * Re_f(qlpm)^B
    # h_avg_f_T(T) = (((cpf_f(T) * mu) / kf_f(T))^(C)) * kf_f(T) / Lc #9. instead of C because library error
    #@register_symbolic h_avg_f_T(x) #does not work if symbolic

    # MOL Discretization parameters for system 1
    x_max1 = L
    x_min1 = 0.0
    t_min = 0.0
    t_max = 7084.0

    nc1 = 100
    x_num1 = range(x_min1, x_max1, length=nc1)
    dx = (x_max1 - x_min1) / (nc1 - 1)

    ## SETUP Parameters

    #p_opt = [A => 2., B => 0.5, n=> 0.5, C=> 20.]
    #p_opt = [A => 1000.0, B => 0.6, C => 7.0, aloss => 1.0]
    p_opt = [A => 9., B => 0.4, C => 10.] #, hext_c => 1.] 
    p_cond = [Io => 456000.0, qlpm => 9.1, Tinit => 293.0]
    p_math = Dict(vcat(p_opt, p_cond))
    #extract the parameters names from p_math
    p_names = [i for i in keys(p_math)]

    # PDE equation for system 1

    # eq1 = [
    #     Vs * (ρs*Cps) * Dt(Ts(t, x)) ~ Vs * (ks) * Dxx(Ts(t, x)) - ((h_average) * Av * Vi * ((Ts(t, x)) - Tf(t, x))) .- (kins * (r / r0) .* (Ts(t, x) .- Tins_f(t)) * A_t / (r - r0)),
    #     Vf * ρf * Cpf * Dt(Tf(t, x)) ~ Vf * kf * Dxx(Tf(t, x)) - Vf * ρf * Cpf * V * Dx(Tf(t, x)) + (h_average) * Av * Vi * ((Ts(t, x) - Tf(t, x)))
    # ]
    eq1 = [
        Vs * Dt(Ts(t, x)) ~ (
            A_frt_solid * ks_f(Ts(t, x)) * Dxx(Ts(t, x))
            - (h_avg_f6(qlpm, (Tf(t, x))) * A_exchange * ((Ts(t, x)) - Tf(t, x))) #h_avg_f(qlpm, x, Tf(t, x)) 
            - (kins * (r_ins / r0) * (Ts(t, x) - Tins_f(t)) * A_s_p / (r_ins - r0))
        ) / (ρs * Cps(Ts(t, x))),
        Vf * Dt(Tf(t, x)) ~ (
            A_chnl_frt_all * kf_f(Tf(t, x)) * Dxx(Tf(t, x))
            -
            m(qlpm) * cpf_f(Tf(t, x)) * Dx(Tf(t, x))
            +
            (h_avg_f6(qlpm, (Tf(t, x))) * A_exchange * ((Ts(t, x)) - Tf(t, x))) # h_avg_f(qlpm, x, Tf(t, x))
        ) / (ρf_f(Tf(t, x)) * cpf_f(Tf(t, x)))
    ]
    bcs1 = [
        Ts(0.0, x) ~ Tinit, # initial
        Tf(0.0, x) ~ Tinit, # initial
        -A_frt_solid *  ks_f(Ts(t, x)) * Dx(Ts(t, x_max1)) ~ 0.0, # far right
        -A_frt_solid *  ks_f(Ts(t, x)) * Dx(Ts(t, x_min1)) ~ α * Io * A_frt
                                                                    - ϵ * σ * A_frt * (Ts(t, x_min1)^4 - Tamb^4)
                                                                    - hext * A_frt * (Ts(t, x_min1) - Tamb),  # far left
        -A_chnl_frt_all * kf_f(Tf(t, x)) * Dx(Tf(t, x_max1)) ~ 0.0, #-ρf * Cpf * V * A_ft * (Tf(t, x_max1) - Tamb), # exiting fluid
        -A_chnl_frt_all * kf_f(Tf(t, x)) * Dx(Tf(t, x_min1)) ~ m(qlpm) * cpf_f(Tf(t, x)) * (Tf(t, x_min1) - Tamb) # entering fluid (upstream temperature)
    ]
    # Space and time domain for system 1
    domains1 = [t ∈ Interval(t_min, t_max),
        x ∈ Interval(x_min1, x_max1)]

    # ODE system for system 1
    @named pdesys = PDESystem(eq1, bcs1, domains1, [t, x], [Ts(t, x), Tf(t, x)], p_names; defaults=p_math)


    # MOL parameters for system 1

    order = 2
    discretization = MOLFiniteDifference([x => dx], t, approx_order=order)

    prob = discretize(pdesys, discretization)

end

include("import_exp2.jl") #import the experimental data
colors = ColorSchemes.tab10[1:4]

#### SETUP OPTIMIZATION ####
begin # define functions
    #Optimization using NLOpt
    #rmp = ModelingToolkit.varmap_to_vars([Io => 442320.0, A => 700.0, B => 0.4, C => 30.0, aloss => 1.0, qlpm => 7.12], parameters(pdesys))
    function NLmodeloptim(tvalues, rmp, tolr)
        #p = [hlocal => p_vary[1]]
        modeloptim = remake(prob, u0=ones(length(prob.u0)) * rmp[Tinit], tspan=(0.0, tvalues[end]))
        #equations(modeloptim.f.initialization_data.initializeprob.f.sys)
        replace!(Tunable(), modeloptim.p, collect(values(rmp)))
        modeloptim_sol = solve(modeloptim, FBDF(), saveat=tvalues, reltol=tolr)#, reltol=1e-12, abstol=1e-12)
        #time = modelfit_sol.t
        #tempT8_op = float.(modeloptim_sol.u[Ts(t, x)][:, 4])
        tempT9_op = float.(modeloptim_sol.u[Ts(t, x)][:, 43])
        tempT10_op = float.(modeloptim_sol.u[Ts(t, x)][:, 78])
        tempT3_op = float.(modeloptim_sol.u[Tf(t, x)][:, end-1])
        # T12_modelmean = (modeloptim_sol.u[Tf(t, x)][end, 20] .+ modeloptim_sol.u[Ts(t, x)][end, 20]) ./ 2
        # T11_modelmean = (modeloptim_sol.u[Tf(t, x)][end, 59] .+ modeloptim_sol.u[Ts(t, x)][end, 59]) ./ 2
        tempT11_op = float.(modeloptim_sol.u[Tf(t, x)][:, 78])
        tempT12_op = float.(modeloptim_sol.u[Tf(t, x)][:, 43])

        # T2 = (tempT9_op .+ tempT11_op) ./2
        # T3 = (tempT12_op .+ tempT10_op) ./2
        return ([tempT9_op tempT10_op tempT3_op tempT11_op tempT12_op])
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

    #to place the conditions loop

        lossr = 0 #zeros(length(sim_key))
        #Threads.@threads 
        for sm in sim_key
            # Retrieve from measurements the experimental data for the current simulation condition
            #sm = sim_key[it]
            #println(sm)
            local cond_k = simulation_conditions[sm]
            local expdata = reduce(hcat, measurements[(measurements.simulation_id.==sm), :temperatures])[:, 1:4] #ignore T8
            #local expdata = measurements[(measurements.simulation_id.==sm), :temperatures]
            local time_opt = measurements[measurements.simulation_id.==sm, :time][1]

            # Ensure time_opt are vectors of floats

            #run selected simulation and get the steady temperature values
            #print(sm)
            temp_T = remakeAysha(pguess_l, cond_k, time_opt)[:, 1:4] #ignore T8
            #temp_T = reshape(temp_T, length(expdata[1]), length(expdata))
            temp_error = sqrt(sum((temp_T .- expdata) .^ 2)) ./ length(expdata)
            lossr += temp_error
        end
        return lossr / length(sim_key) #MSE
    end
end

begin
    p_opt = [A => 0.5, B => 0.04, C => 15.]

    p0 = [x[2] for x in p_opt]
    optf = OptimizationFunction(lossAll, Optimization.AutoForwardDiff())
    lb = [0.001, 0.0001, 0.0001]
    ub = [10., 10., 20.] 

    #pguess_opt = ModelingToolkit.varmap_to_vars([Io => 456000, h_average => 14., qlpm => 7.12], parameters(pdesys))
    #sim_key = collect(keys(simulation_conditions))
    sim_key = ["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"]
    initialerror = (lossAll(p0, sim_key))
    println("Initial Error $initialerror")
    optprob = Optimization.OptimizationProblem(optf, p0, sim_key, lb=lb, ub=ub)
    optsol = solve(optprob, NLopt.GN_MLSL_LDS(), local_method=NLopt.LN_NELDERMEAD(), maxtime=60, local_maxiters=100000)
    #optsol = solve(optprob, ECA(), maxtime=20, maxiters=100000)
    
    println(optsol.retcode)
    pnew = optsol.u
    println(pnew)
    #res_error = lossAll(pnew, [])
    println("Final Error $(optsol.objective)")
end

#### FINAL RESULTS - FIX DATA EXTRACTION
#pnew = [636.0, 0.44, 8.488, 1.5]
begin #ONLY T3 gas temperature (as in the measurements dataframe)
    T_steady = DataFrame(sim_id=[], time=[], T_mod=[], T_exp=[])
    #sim_key = collect(keys(simulation_conditions))
    #plotting the scatter plots for the different thermocouples for the steady state temperatures
    ordered_sim_conditions = sim_key#["E67", "E68", "E69", "E70", "E71", "E72", "E73", "E74", "E75", "E76", "E77", "E78", "E79", "E80", "E81"]
    #order_indices = [findfirst(x -> x == condition, T_steady.sim_id) for condition in ordered_sim_conditions]
    #Threads.@threads 
    for sm in ordered_sim_conditions
        # Retrieve from measurements the experimental data for the current simulation condition
        #println(sm)
        local cond_k = simulation_conditions[sm]
        local sel_meas = measurements[(measurements.simulation_id.==sm).&(measurements.obs_id.=="_Tf"), :]
        local expdata = sel_meas[:, :temperatures]
        local time_opt = sel_meas[:, :time][1]

        #run selected simulation and get the steady temperature values
        local temp_T = remakeAysha(pnew, cond_k, time_opt; tolr=1e-7)
        #temp_T = reshape(temp_T, length(expdata[1]), length(expdata))
        #temp_T = reshape(temp_T, length(expdata), length(expdata))
        #temp_T = [temp_T[:, i] for i in 1:length(expdata)]
        push!(T_steady, (sm, time_opt, temp_T[:,4], expdata[1]))
    end

    #TI-8
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

    #plot1 = qqplot(T_mod_values, T_exp_values, label="",
    #xlabel="T_mod", ylabel="T_exp",ylims=(500, 900),xlims=(500, 900), aspect_ratio = :equal)
    lims = (500, 900)
    plot1 = scatter(T_mod_values, T_exp_values, label="",
        xlabel="T_mod", ylabel="T_exp", ylims=lims, xlims=lims, aspect_ratio=:equal)
    plot1 = plot!(collect(lims), collect(lims))

    plot(plot0, plot1)
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
        local expdata = reduce(hcat, measurements[(measurements.simulation_id.==sm), :temperatures])[:, 1:4] #ignore T8
        local time_opt = measurements[measurements.simulation_id.==sm, :time][1]
        local temp_T = remakeAysha(params, cond_k, time_opt)[:, 1:4] #ignore T8

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


