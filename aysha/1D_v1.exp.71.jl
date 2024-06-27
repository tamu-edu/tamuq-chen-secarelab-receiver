begin #libraries
    using MethodOfLines
    using ModelingToolkit
    using DomainSets, OrdinaryDiffEq
    using NonlinearSolve, DifferentialEquations
    using Plots, XLSX, Statistics, Symbolics, Interpolations
    using Optim, LsqFit, LossFunctions
    using Optimization, OptimizationNLopt, Symbolics, OptimizationOptimJL, ForwardDiff, OptimizationMOI
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
    th_s = 0.4e-3 #m
    th_f = 0.7e-3 #m
    w = 1.5e-3 #width of channel (m)
    n_chnl = 10*10
    A_frt = w * w * n_chnl #m2 - for the whole receiver (19x19mm2)
    A_s_p = w_t * L * 4 #total area solid periphery m2
    w_chnl = 1.5e-3 #m
    A_chnl_p = w_chnl * L * 4  #m2 channel periphery
	A_chnl_frt = w_chnl * w_chnl #m2 channel front
    #A_exchange = A_chnl_p * n_chnl #m2 - contact area between fluid and solid
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
    #Qv= I0*exp(-1000*x)  #kW/m2 - (K extinction coefficient taken from Howell and Hendrick paper pg.86, measure pore diameter and fix)
    ϵ = 0.8
    σ = 5.17e-8 #W/m2.K^4 Stefan-Boltzmann constant
    Lc = 4 * (w_t * w_t) / (4 * w_t)
    kf = 0.056 #W/m.K
    ρf = 0.5 #kg/m3
    Cpf = 1090 #J/kg.K
    mu = 2.0921e-5 #Pa.s
    e = 0.425
    Af = n_chnl * (w + th_s)^2 #m2
    Vi = w * w * n_chnl * L #m3
    #Av = 4 * (w * L) / (w^2 * L) #specific area (m-1)
    Av = (4 * w * L) / ((w + th_s)^2 * L) #1/m
    A_exchange = Av * Vi #m2 - contact area between fluid and solid
    hext = 10 #W/m2.K
    kins = 0.078 #W/m*K
    r0 = 23 / 1000 #m
    r_ins = 42 / 1000 #m
    r_H = 4 * (w_t * w_t) / (4 * w_t) #hydraulic receiver diameter
    em = 0.8 #emissivity
    aCp = 4.
    aIo = 2.
    al = 0.22
end;
    #for interpolations 
    #1. extract T2 data
    #Exp 71 
    begin
     D3 = XLSX.readxlsx("/Users/aishamelhim/Documents/GitHub/tamuq-chen-secarelab-receiver/aysha/SolarSimulator/EXCEL/Data_FPT0071_T2.xlsx")["Sheet 1 - Data_FPT0071_T2"]["A3:D7087"]
     y1d3_data  = D3[:,2].+ 273. #T2 (insulation)
     x1 = D3[:,1]
     #2.create interpolation function
     Tins = linear_interpolation(x1, y1d3_data)
     Tins_f(t) = Tins(t)
     @register_symbolic Tins_f(t)
    end

    begin 
        Re = (ρf * V * w_t) / mu
        Pr = (Cpf * mu) / kf
        x11  = 0.0001:0.001383838383838384:0.137 #T2 (insulation)
        Gz = (1 ./x11)*Re*Pr*w_t
        #2.create interpolation function
        Gz_ = linear_interpolation(x11, Gz)
        Gz_f(x) = Gz_(x)
        @register_symbolic Gz_f(x)
    end
    begin
        # Parameters, variables, and derivatives for system 1
        @variables t x
         @parameters h_average #A B ks h_average A n
        @parameters Io qlpm
        @variables Ts(..) Tf(..)
        Dt = Differential(t)
        Dx = Differential(x)
        Dxx = Differential(x)^2
        ks = (48.78) * 1.97 * ((1 - e)^1.5) #W/m.K
        ρs = 3200  #kg/m3
        Cps = 1290  #J/kg*K
        #Gz = (1/L) * Re * Pr * w_t
        #Cps(Ts) = (0.27+0.135e-4*(Ts)-9720*((Ts)^-2)+0.204e-7*((Ts)^2))/1000 #kJ/kg*K from manufacturer data
        # Nu = A*(1+(B*((Gz_f(x))^n)*exp(-C/Gz_f(x))))
        # Nu = A * (Re)^B
        p_opt = [h_average => 0.5]
        p_cond = [Io => 456000.0, qlpm => 7.12]
        p_math = vcat(p_opt, p_cond)
        # h_average = (Nu * kf) / Lc
     
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
            (1-e) * (aCp * ρs * Cps) * Vs * Dt(Ts(t, x)) ~ A_frt * ks * Dxx(Ts(t, x)) - (h_average * A_exchange * ((Ts(t, x)) - Tf(t, x))) .-  al * (kins * (r_ins/r0).* (Ts(t, x) .- Tins_f(t)) * A_s_p/ (r_ins - r0)),
            e * ρf * Cpf * Vf * Dt(Tf(t, x)) ~ Af * kf * e * Dxx(Tf(t, x)) -  m * Cpf * Dx(Tf(t, x)) + (h_average * A_exchange * ((Ts(t, x)) - Tf(t, x)))
        ]
        bcs1 = [
            Ts(0.0, x) ~ Tamb, # initial
            Tf(0.0, x) ~ Tamb, # initial
            -A_frt * ks * Dx(Ts(t, x_max1)) ~ 0.0, # far right
            -A_frt * ks * Dx(Ts(t, x_min1)) ~ aIo * Io * A_frt - ϵ * σ * A_frt * (Ts(t, x_min1)^4 - Tamb^4) - hext * A_frt * (Ts(t, x_min1) - Tamb),  # far left
            -Af * kf * Dx(Tf(t, x_max1)) ~ 0.0, #-ρf * Cpf * V * A_ft * (Tf(t, x_max1) - Tamb), # exiting fluid
            -Af * kf * Dx(Tf(t, x_min1)) ~ m * Cpf * (Tf(t, x_min1) - Tamb) # entering fluid (upstream temperature)
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
   
            
    sol1 = solve(prob, FBDF(), saveat=2, maxiters = 100)
        begin
            Ts_front_t = sol1.u[(Ts(t,x))][:,1]
            Tf_front_t = sol1.u[(Tf(t,x))][:,2]
            Ts_back_t = sol1.u[(Ts(t,x))][:,end-1]
            Tf_back_t = sol1.u[(Tf(t,x))][:,end]
            x_domain = collect(sol1.ivdomain[2])
            plot(xlabel="time [s]")
            plot!(sol1.t, Ts_front_t, label="T_fr_s")
            plot!(sol1.t, Tf_front_t, label="T_fr_f")
            plot!(sol1.t, Ts_back_t, label="T_bck_s")
            plot!(sol1.t, Tf_back_t, label="T_bck_f")
        end
        
        x__domain = collect(sol1.ivdomain[2])
        #Data extraction
        begin
            #Exp 71
            D1 = XLSX.readxlsx("/Users/aishamelhim/Documents/GitHub/tamuq-chen-secarelab-receiver/aysha/SolarSimulator/EXCEL/Data_FPT0071_231128_102707.xlsx")["Sheet 1 - Data_FPT0071_231128_1"]["A3:C7087"]
            xd1_data = D1[:,1] #time
            y1d1_data = D1[:,2] .+ 273. #T3
            y2d1_data = D1[:,3] .+ 273. #T8
        end
        begin
            #Exp 71 
            D2 = XLSX.readxlsx("/Users/aishamelhim/Documents/GitHub/tamuq-chen-secarelab-receiver/aysha/SolarSimulator/EXCEL/Data_FPT0071_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0071_231128_1"]["A3:E7087"]
            y1d2_data  = D2[:,2].+ 273. #T9 (external)
            y2d2_data  = D2[:,3].+ 273. #T10 (external)
            y3d2_data  = D2[:,4].+ 273. #T11 (internal)
            y4d2_data = D2[:,5].+ 273. #T12 (internal)
        end
        # Model solution plotted with experimental data
        begin
            plot(title = "Solid Temperature Profile T8")
            plot!(sol1.t, 
                sol1.u[Ts(t,x)][:,4], # around 5 mm in T8
                title = "Solid Temperature Profile T8", 
                label = "Numerical", 
                xlabel = "Time (s)", 
                ylabel = "Temperature (K)")
                scatter!(
                    xd1_data, 
                    y2d1_data,  
                    label = "Experimental",
                    xlabel = "Time (s)", 
                    ylabel = "Temperature (K)")
            
        end
        begin
            plot()
                plot!(sol1.t, 
                sol1.u[Tf(t,x)][:,end-1], # around 136 mm in T3
                title = "Gas Temperature Profile T3", 
                label = "Numerical", 
                xlabel = "Time (s)", 
                ylabel = "Temperature (K)")
                scatter!(
                    xd1_data, 
                    y1d1_data,  
                    label = "Experimental",
                    xlabel = "Time (s)", 
                    ylabel = "Temperature (K)")
        end
    begin
        # model results for different thermocouples
       sol1.u[Ts(t,x)][:, 4]
       sol1.u[Ts(t,x)][:,40] #T9 model 55mm (external)
       sol1.u[Ts(t,x)][:,77] #T10 model 106mm (external)
       sol1.u[Tf(t,x)][:,59] #T11 model 82mm (internal-gas)
       sol1.u[Ts(t,x)][:,59] #T11 model 82mm (internal-solid)
       sol1.u[Tf(t,x)][:,20] #T12 model 27mm (internal-gas)
       sol1.u[Ts(t,x)][:,20] #T12 model 27mm (internal-solid)
        # taking averages of internal thermocouples
       # T12_model = sol1.u[Tf(t,x)][:,20], sol1.u[Ts(t,x)][:,20]
       # T12_modelmean = mean(T12_model)
       # T11_model = sol1.u[Tf(t,x)][:,59], sol1.u[Ts(t,x)][:,59]
        #T11_modelmean = mean(T11_model)
    end
        p0 = [x[2] for x in p_opt]
        #expdata = append!(copy(y2d1_data), y1d2_data, y2d2_data, y1d1_data, y3d2_data, y4d2_data)
        expdata =  append!(copy(y2d1_data), y1d2_data, y2d2_data, y1d1_data, y3d2_data, y4d2_data)
        length(expdata)
        
        
        pguess = p0
        rmp = ModelingToolkit.varmap_to_vars([Io => 456000, ks => 22.5, h_average=> 20., qlpm => 7.12], parameters(pdesys))

            #Optimization using NLOpt
            function NLmodeloptim(xvalues, rmp)
                #p = [hlocal => p_vary[1]]
                modeloptim = remake(prob, p = rmp, tspan=(xvalues[1], xvalues[end]))
                modeloptim_sol = solve(modeloptim, FBDF(), saveat = xvalues)#, reltol=1e-12, abstol = 1e-12)
                #time = modelfit_sol.t
                tempT8_op = modeloptim_sol.u[Ts(t,x)][:, 4]
                tempT9_op = modeloptim_sol.u[Ts(t,x)][:,40]
                tempT10_op = modeloptim_sol.u[Ts(t,x)][:,77]
                tempT3_op = modeloptim_sol.u[Tf(t,x)][:, end-1]
                T12_modelmean = (modeloptim_sol.u[Tf(t,x)][:,20] .+  modeloptim_sol.u[Ts(t,x)][:,20]) ./2
                T11_modelmean = (modeloptim_sol.u[Tf(t,x)][:,59] .+ modeloptim_sol.u[Ts(t,x)][:,59]) ./2
                return append!(tempT8_op, tempT9_op, tempT10_op, tempT3_op, T12_modelmean, T11_modelmean) 
            end
        
            function loss(pguess, _)
                Temp = NLmodeloptim(xd1_data, pguess)
                lossr = (Temp .- expdata).^2
                return sqrt(sum(lossr) / length(expdata)) #MSE
            end

            #loss([1.], [])

        initialerror=(loss(pguess, []))

        optf = OptimizationFunction(loss, Optimization.AutoForwardDiff())
        
        lb = [0.01]
        ub = [100.]

        optprob = Optimization.OptimizationProblem(optf, p0, [], lb=lb, ub=ub)
        
        optsol = solve(optprob, NLopt.GN_MLSL_LDS(), local_method = NLopt.LN_NELDERMEAD(), maxtime = 500, local_maxiters = 1000)
        
        
        println(optsol.retcode)
        
        
        pnew = optsol.u
        println("Size of pnew: ", size(pnew))
        println("Contents of pnew: ", pnew)
        
        begin
            # res_error = loss(pnew, [])
            # display(res_error)
            modelfit = remake(prob, p = rmp)
            modelfit_sol = solve(modelfit, FBDF(), saveat = xd1_data, reltol=1e-12, abstol = 1e-12)
            psol = modelfit_sol
        end
begin
    plot1 =  plot(
        psol.t, 
        psol.u[Tf(t,x)][:,end-1], # around 136 mm in T3
        title = "Gas Temperature Profile T3 - exp. 71", 
        label = "optimized", 
        xlabel = "Time (s)", 
        ylabel = "Temperature (K)")
        scatter!(
        xd1_data, 
        y1d1_data,  
        label = "Experimental",
        xlabel = "Time (s)", 
        ylabel = "Temperature (K)")
        plot!(
        psol.t, 
        psol.u[Ts(t,x)][:,end-1], 
        label = "solid 136 mm")
            
    plot2 = plot(
        psol.t, 
        psol.u[Ts(t,x)][:,4], # around 5 mm in T8
        title = "Solid Temperature Profile T8 - exp. 71", 
        label = "Numerical", 
        xlabel = "Time (s)", 
        ylabel = "Temperature (K)")
        scatter!(
        xd1_data, 
        y2d1_data,  
        label = "Experimental",
        xlabel = "Time (s)", 
        ylabel = "Temperature (K)")
        plot!(
        psol.t, 
        psol.u[Tf(t,x)][:,4], # around 5 mm in T8
        label = "gas 5mm")
    plot3 = plot(
        psol.t, 
        psol.u[Ts(t,x)][:,40], # around 55 mm in T9
        title = "Solid Temperature Profile T9 - exp. 71", 
        label = "Numerical", 
        xlabel = "Time (s)", 
        ylabel = "Temperature (K)")
        scatter!(
        xd1_data, 
        y1d2_data,  
        label = "Experimental",
        xlabel = "Time (s)", 
        ylabel = "Temperature (K)")
        plot!(
        psol.t, 
        psol.u[Tf(t,x)][:,40], 
        label = "gas 55mm")
    plot4 = plot(
        psol.t, 
        psol.u[Ts(t,x)][:,77], # around 106 mm in T10
        title = "Solid Temperature Profile T10 - exp. 71", 
        label = "Numerical", 
        xlabel = "Time (s)", 
        ylabel = "Temperature (K)")
        scatter!(
        xd1_data, 
        y2d2_data,  
        label = "Experimental",
        xlabel = "Time (s)", 
        ylabel = "Temperature (K)")
        plot!(
        psol.t, 
        psol.u[Tf(t,x)][:,77], 
        label = "gas 106mm")
    plot5 = plot(
        psol.t, 
        psol.u[Ts(t,x)][:,20], # around 27 mm in T12
        title = "Solid Temperature Profile T12 - exp. 71", 
        label = "Numerical", 
        xlabel = "Time (s)", 
        ylabel = "Temperature (K)")
        scatter!(
        xd1_data, 
        y4d2_data,  
        label = "Experimental",
        xlabel = "Time (s)", 
        ylabel = "Temperature (K)")
        plot!(
        psol.t, 
        psol.u[Tf(t,x)][:,20], 
        label = "gas 27mm")
    plot6 = plot(
        psol.t, 
        psol.u[Ts(t,x)][:,59], # around 82 mm in T11
        title = "Solid Temperature Profile T11 - exp. 71", 
        label = "Numerical", 
        xlabel = "Time (s)", 
        ylabel = "Temperature (K)")
        scatter!(
        xd1_data, 
        y3d2_data,  
        label = "Experimental",
        xlabel = "Time (s)", 
        ylabel = "Temperature (K)")
        plot!(
        psol.t, 
        psol.u[Tf(t,x)][:,59], 
        label = "gas 82mm")
    plot7 = plot(
        collect(x_num1),
        psol.u[Tf(t,x)][end,:],
        label = "gas",
        xlabel = "Length (m)",
        ylabel = "Temperature (K)")
        plot!(
        collect(x_num1),
        psol.u[Ts(t,x)][end,:],
        label = "solid",
        xlabel = "Length (m)",
        ylabel = "Temperature (K)")
        
    plot(plot1, plot2, plot3, plot4, plot5, plot6, plot7, layout=(10,2), size=(1500, 1500))
end