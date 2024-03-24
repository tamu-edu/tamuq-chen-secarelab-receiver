begin #libraries
    using MethodOfLines
    using ModelingToolkit
    using DomainSets, OrdinaryDiffEq
    using NonlinearSolve, DifferentialEquations
    using Plots, XLSX, Statistics, Symbolics, Interpolations
    using Optim, LsqFit
    using Optimization, OptimizationNLopt, Symbolics, OptimizationOptimJL, ForwardDiff, OptimizationMOI
end
begin #define parameters

    #ks= (37/1000)*(1-e) #C SiC- #kW/m.K using equation plotted from (Ali et al.)
    #ks=(52000*exp(-1.24e-5*T)/(T+437))/1000 #kW/m.K (Ali et al.)
    #kf = (0.056/1000)*e #thermal conductivity of the fluid phase kW/m.K
    #kf=(1.52e-11*(T^3)-4.86e-8*(T^2)+1.02e-4*T-3.93e-3)/1000 #kW/m.K (Ali et al.)
    #keff=0.82*kf+(1-0.82)*ks #using "Thermal analysis and design of a volumetric solar absorber depending on the porosity" paper - note that I still need to measure the porosity of SiC 
    #ρf= 3.018*exp(-0.00574*T)+0.8063*exp(-0.0008381*T) #kg/m3 (Roldan et al.)
    #Cpf=(1090/1000)*e #kJ/kg.K 
    #Cps=1225/1000 #kJ/kg.K 
    #ρs= (3100)*(1-e) #kg/m3
    L= 137e-3 #m
    #α = keff/ρ*Cp #kW/m2.K 
    Tamb = (22.448 + 273.15) #K (same for all exp) 
    deltax = 0.002795918367346939 #discretization (m)
    A_t= 324e-6 #m2 - for the whole receiver (18x18mm2)
    
    A_st= 176e-6 #total area - solid m2
    A_ft= 49e-6 # total area - fluid m2
    A_channel = 2.25e-6 #m2 (1.5x1.5mm2) 
    n_channel = 100
    A_exchange = 162.5e-6 #m2 - contact area between fluid and solid
    Vs =  A_st * L #* deltax
    Vf =  A_ft * L #* deltax
    qlpm = 7.12 #lpm
    q= qlpm/(1000*60) #m3/s
    V= q/(A_channel*n_channel) #m/s - using the area of the whole receiver (18x18mm2)
    th_s = 0.4e-3 #m
    th_f = 0.7e-3 #m
    I0 = 456 #kW/m2
    #Qv= I0*exp(-1000*x)  #kW/m2 - (K extinction coefficient taken from Howell and Hendrick paper pg.86, measure pore diameter and fix)
    Q= I0 #kW/m2
    ϵ= 0.8
    σ= (5.17e-8)/1000 #kW/m2.K^4 Stefan-Boltzmann constant
    Nu = 3.657 
    #Nu = 0.00001
    w_t = 18/1000 #m
    Lc = 4*(w_t*w_t)/(4*w_t)
    #kf = (0.056/1000) #kW/m.K
    kf = (0.02/1000)
    w = 1.5e-3 #width of channel (m)
    Vi = w * w * n_channel * L #m3
    Av = 4*(w*L) / (w^2*L) #specific area (m-1)
    hext = 10/1000 #kW/m2.K
    kins = 0.078/1000 #kW/m*K
    r0 = 23/1000 #m
    r = 42/1000 #m
    Tins = 356 #K exp. 71
end;
    #for interpolations 
    #1. extract T2 data
    #Exp 71 
    begin
     D3 = XLSX.readxlsx("/Users/aishamelhim/Documents/GitHub/tamuq-chen-secarelab-receiver/aysha/SolarSimulator/EXCEL/Data_FPT0071_T2.xlsx")["Sheet 1 - Data_FPT0071_T2"]["A3:D7087"]
     y1d3_data  = D3[:,2].+ 273. #T2 (insulation)
     x1 = D3[:,1].+ 273.
     #2.create interpolation function
     itp = interpolate((x1,), y1d3_data, Gridded(Linear()))
     # Evaluate the interpolation function 
    Tins = itp(x1)
    #Tins_ = collect(Tins)
    end
    #  function Tins(t)
    #    return itp(x1)  # Assuming itp is an interpolation function of Tins
    #  end
     begin
        # Parameters, variables, and derivatives for system 1
        @variables t x 
        @parameters psCps h_average ks
        @variables Ts(..) Tf(..) 
        Dt = Differential(t) 
        Dx = Differential(x)
        Dxx = Differential(x)^2

        e = 0.62
        #ks = (37/1000)*(1-e) #kW/m.K
        #ρs = 2170 #kg/m3
        #Cps = (1225/1000) #kJ/kg
        ρf = 0.5 #kg/m3
        Cpf = (1090/1000) #kJ/kg.K

        p = [psCps => 900., ks=> 0.037, h_average => (Nu*kf/Lc)]
        #Cps(T)= (1110+0.15*(T-273)-425*exp(-0.003*(T-273)))/1000 #kJ/kg*K
        #Cpf(T)= (1.93e-10*(T^3)+1.14e-3*(T^2)-4.49e-1*T+1.06e3)/1000 #kJ/kg*K

        #psCps = ρs * Cps * (1-e)

        # MOL Discretization parameters for system 1
        x_max1 = L
        x_min1 = 0.
        t_min = 0.
        t_max = 7084.
        
        nc1 = 100
        
        x_num1 = range(x_min1, x_max1, length = nc1)
        
        
        dx = (x_max1 - x_min1) / (nc1 - 1)
        
        
        # PDE equation for system 1
        
        eq1 = [
               Vs * psCps * Dt(Ts(t,x)) ~ Vs * ks * Dxx(Ts(t,x)) - ((h_average) * Av * Vi * ((Ts(t,x)) - Tf(t,x))) .- (kins * (r/r0) .* (Ts(t,x) .- Tins)) * A_t / (r-r0)
               Vf* ρf * Cpf * Dt(Tf(t,x)) ~ Vf * kf * Dxx(Tf(t,x)) - Vf * ρf * Cpf * V * Dx(Tf(t,x)) + (h_average) * Av * Vi * ((Ts(t,x) - Tf(t,x)))
               ]
              
        bcs1 = [
            Ts(0., x) ~ Tamb, # initial
            Tf(0., x) ~ Tamb, # initial
            -A_st * ks * Dx(Ts(t, x_max1)) ~ 0.0, # far right
            -A_st * ks * Dx(Ts(t, x_min1)) ~ I0 * A_st - ϵ * σ * A_st * (Ts(t,x_min1)^4 - Tamb^4) - hext * A_st * (Ts(t, x_min1) - Tamb),  # far left
            -A_ft * kf * Dx(Tf(t, x_max1)) ~ -ρf * Cpf * V * A_ft * (Tf(t, x_max1) - Tamb),
            -A_ft * kf * Dx(Tf(t, x_min1)) ~ ρf * Cpf * V * A_ft * (Tf(t,x_min1) - Tamb) # fluid entering (fix upstream temperature)
              ] 
        # Space and time domain for system 1
        domains1 = [t ∈ Interval(t_min, t_max),
                    x ∈ Interval(x_min1, x_max1)]
        
        # ODE system for system 1
        @named pdesys = PDESystem(eq1, bcs1, domains1, [t, x], [Ts(t, x), Tf(t, x)], p)
        
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
            D2 = XLSX.readxlsx("/Users/aishamelhim/Documents/GitHub/tamuq-chen-secarelab-receiver/aysha/SolarSimulator/EXCEL/Data_FPT0071_T9,10,11,12.xlsx")["Sheet 1 - Data_FPT0071_231128_1"]["A3:D7087"]
            y1d2_data  = D2[:,1].+ 273. #T9 (external)
            y2d2_data  = D2[:,2].+ 273. #T10 (external)
            y3d2_data  = D2[:,3].+ 273. #T11 (internal)
            y4d2_data = D2[:,4].+ 273. #T12 (internal)
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
        p0 = [x[2] for x in p]
        #expdata = append!(copy(y2d1_data), y1d2_data, y2d2_data, y1d1_data, y3d2_data, y4d2_data)
        expdata =  append!(copy(y2d1_data), y1d2_data, y2d2_data, y1d1_data, y3d2_data, y4d2_data)
        length(expdata)
        pguess = p

            #Optimization using NLOpt
            function NLmodeloptim(xvalues, p_vary)
                #p = [hlocal => p_vary[1]]
                modeloptim = remake(prob, p = p_vary, tspan=(xvalues[1], xvalues[end]))
                modeloptim_sol = solve(modeloptim, FBDF(), saveat = xvalues, reltol=1e-12, abstol = 1e-12)
                #time = modelfit_sol.t
                tempT8_op = modeloptim_sol.u[Ts(t,x)][:, 4]
                tempT9_op = modeloptim_sol.u[Ts(t,x)][:,40]
                tempT10_op = modeloptim_sol.u[Ts(t,x)][:,77]
                tempT3_op = modeloptim_sol.u[Tf(t,x)][:, end-4]
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
    
        initialerror=(loss(p0, []))

        optf = OptimizationFunction(loss, Optimization.AutoForwardDiff())
        
        lb = [0.0, 0.0, 0.0]
        ub = [1.0e6, 1.0e2, 1.0e2]
        #lb = [100.]
        #ub = [1000.]
        optprob = Optimization.OptimizationProblem(optf, p0, [], lb=lb, ub=ub)
        optsol = solve(optprob, NLopt.GN_MLSL_LDS(), local_method = NLopt.LN_NELDERMEAD(), maxtime = 180, local_maxiters = 10000)
        
        
        println(optsol.retcode)
        
        
        pnew = optsol.u
begin
        
        #pnew = [0.3]
        res_error = loss(pnew, [])
        display(res_error)
        modelfit = remake(prob, p = pnew)
        modelfit_sol = solve(modelfit, FBDF(), saveat = xd1_data, reltol=1e-12, abstol = 1e-12)
        
    
        psol = modelfit_sol

    plot1 =  plot(
        psol.t, 
        psol.u[Tf(t,x)][:,end-4], # around 135 mm in T3
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
        psol.u[Ts(t,x)][:,end-4], 
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