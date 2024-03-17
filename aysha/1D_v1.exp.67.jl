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
    #solid
    #ks= (45/1000)*(1-0.62) #C SiC- #kW/m.K from the web (imetra.com))
    #ks=(52000*exp(-1.24e-5*T)/(T+437))/1000 #kW/m.K (Ali et al.)
    #k=ks*(1-0.88) # effective axial heat conductivty for monolithic structures pg.285 (structured cayalysts and reactors book), using apparent emissivity pg.240 (thermal radiation heat transfer book)
    #kf = (0.055/1000)*(0.62) #thermal conductivity of the fluid phase kW/m.K
    #kf=(1.52e-11*(T^3)-4.86e-8*(T^2)+1.02e-4*T-3.93e-3)/1000 #kW/m.K (Ali et al.)
    #keff=0.82*kf+(1-0.82)*ks #using "Thermal analysis and design of a volumetric solar absorber depending on the porosity" paper - note that I still need to measure the porosity of SiC 
    #ρf= 3.018*exp(-0.00574*T)+0.8063*exp(-0.0008381*T) #kg/m3 (Roldan et al.)
    #Cps=1110+0.15*T-425*exp(-0.003*T)
    #Cpf=1.93e-10*(T^3)+1.14e-3*(T^2)-4.49e-1*T+1.06e3
    #Cpf= (1090/1000)*0.62 #kJ/kg.K 
    #Cps=1225/1000 #kJ/kg.K 
    ρs= 3100*(1-0.62) #kg/m3
    ρf= 0.0155*0.62 #kg/m3
    L= 137e-3 #m
    #α = keff/ρ*Cp #kW/m2.K 
    Tamb= (23.409 + 273.15) #K (same for all exp) 
    deltax = 0.002795918367346939 #discretization (m)
    A_t= 324e-6 #m2 - for the whole receiver (18x18mm2)
    
    A_st= 176e-6 #total area - solid m2
    A_ft= 49e-6 # total area - fluid m2
    A_channel = 2.25e-6 #m2 (1.5x1.5mm2) 
    n_channel = 100
    A_exchange = 162.5e-6 #m2 - contact area between fluid and solid
    Vs =  A_st * deltax
    Vf =  A_ft * deltax
    qlpm = 15.27 #lpm
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
    w = 1.5e-3 #width of channel (m)
    Vi = w * w * deltax * n_channel #m3
    Av = 1662 #specific area (m-1)
    hext = 10/1000 #kW/m2.K
    kins = 0.078/1000 #kW/m*K
    r0 = 23/1000 #m
    r = 42/1000 #m
    Tins = 326.437 #K exp. 71
    #hf= 10/1000 #kW/m2.K - check value (free convection)
    #hs= 0.2/1000 #kW/m2.K 
    #I0 * exp(-2300 * T(t, x)) #attenuation term negligible
end;
    #begin
        # Parameters, variables, and derivatives for system 1
        @variables t x
        @parameters hlocal Cps 
        @variables Ts(..) Tf(..)
        Dt = Differential(t) 
        Dx = Differential(x)
        Dxx = Differential(x)^2
        
        p = [hlocal => (10000.), Cps => (1225.)]
        e = 0.3
        ks= (45/1000)*(1-e) #kW/m.K
        kf = (0.055/1000)*(e) #kW/m.K
        Cpf= (1090/1000)*e 
        Cpf= (1090/1000)*e
        ρs= 3100*(1-e) #kg/m3
        #p = [hlocal => (800.)]
        #Cps(T)= (1110+0.15*(T-273)-425*exp(-0.003*(T-273)))/1000 #kJ/kg*K
        #Cpf(T)= (1.93e-10*(T^3)+1.14e-3*(T^2)-4.49e-1*T+1.06e3)/1000 #kJ/kg*K

        
        # MOL Discretization parameters for system 1
        x_max1 = L
        x_min1 = 0.
        t_min = 0.
        t_max = 4000.
        nc1 = 50
        
        x_num1 = range(x_min1, x_max1, length = nc1)
        
        
        dx = (x_max1 - x_min1) / (nc1 - 1)
        
        
        # PDE equation for system 1
        
        eq1 = [
            Vs * ρs * Dt(Ts(t,x)) ~ (Vs * ks * Dxx(Ts(t,x))/(Cps/1000)) + (((hlocal/(1000)) * Av * Vi * ((Ts(t,x)) - Tf(t,x)))/(Cps/1000)) - (((kins * (r/r0) * (Ts(t,x)-Tins)) * A_t/ (r-r0))/(Cps/1000)),
            Vf* ρf * Dt(Tf(t, x)) ~ (- Vf * ρf * V * Dx(Tf(t,x))) + ((hlocal/(1000)) * Av * Vi * ((Ts(t,x) - Tf(t,x)))/(Cpf/1000))
            ]    
            
        bcs1 = [
            Ts(0., x) ~ Tamb, # initial
            Tf(0., x) ~ Tamb, # initial
            -A_st * ks * Dx(Ts(t, x_max1)) ~ 0, # far right
            -A_st * ks * Dx(Ts(t, x_min1)) ~ I0 * A_st - ϵ * σ * A_st * (Ts(t,x_min1)^4 - Tamb^4) - hext * A_st * (Ts(t, x_min1) - Tamb),  # far left
            -A_ft * kf * Dx(Tf(t, x_max1)) ~ -ρf * (Cpf/1000) * V * A_ft * (Tf(t, x_max1) - Tamb) # fluid exiting
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
            
            
            
    #end
    #begin
            
    sol1 = solve(prob, FBDF(), saveat=2, maxiters = 100)
        
     #   end
        begin
            Ts_front_t = sol1.u[(Ts(t,x))][:,1]
            Tf_front_t = sol1.u[(Tf(t,x))][:,1]
            Ts_back_t = sol1.u[(Ts(t,x))][:,end-1]
            Tf_back_t = sol1.u[(Tf(t,x))][:,end-1]
            x_domain = collect(sol1.ivdomain[2])
            plot(xlabel="time [s]")
            plot!(sol1.t, Ts_front_t, label="T_fr_s")
            plot!(sol1.t, Tf_front_t, label="T_fr_f")
            plot!(sol1.t, Ts_back_t, label="T_bck_s")
            plot!(sol1.t, Tf_back_t, label="T_bck_f")
        end
        x__domain = collect(sol1.ivdomain[2])
        
        begin
           #Exp 67
            Z = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0067_231125_161757.xlsx")["Sheet 1 - Data_FPT0067_231125_1"]["A3:C3932"]
            xz_data = Z[:,1]
            y1z_data = Z[:,2] .+ 273.
            y2z_data = Z[:,3] .+ 273.
        end
        
        begin
            plot(
                sol1.t, 
                sol1.u[Ts(t,x)][:,3], # around 5 mm in T8
                title = "Solid Temperature Profile T8", 
                label = "Numerical", 
                xlabel = "Time (s)", 
                ylabel = "Temperature (K)")
                scatter!(
                    xz_data, 
                    y2z_data,  
                    label = "Experimental",
                    xlabel = "Time (s)", 
                    ylabel = "Temperature (K)")
            
        end
        begin
            plot(
                sol1.t, 
                sol1.u[Tf(t,x)][:,49], # around 135 mm in T3
                title = "Gas Temperature Profile T3", 
                label = "Numerical", 
                xlabel = "Time (s)", 
                ylabel = "Temperature (K)")
                scatter!(
                    xz_data, 
                    y1z_data,  
                    label = "Experimental",
                    xlabel = "Time (s)", 
                    ylabel = "Temperature (K)")
        end
       
        p0 = [x[2] for x in p]
        #expdata = append!(copy(y2z_data), y1z_data)
        expdata =  y1z_data
        length(expdata)
        pguess = p
           
        # Set up initial guesses for parameters
            
            # Perform the optimization using Lsqfit
        
            #result = LsqFit.curve_fit(objective_function, xd1_data, expdata, p0)
            # Get the optimized parameters
            # optimized_parameters = result.param
            # result.converged, result.resid

            #Optimization using NLOpt
            function NLmodeloptim(xvalues, p_vary)
                #p = [hlocal => p_vary[1]]
                modeloptim = remake(prob, p = p_vary, tspan=(xvalues[1], xvalues[end]))
                modeloptim_sol = solve(modeloptim, FBDF(), saveat = xvalues, reltol=1e-12, abstol = 1e-12)
                #time = modelfit_sol.t
                #tempTs_op = modeloptim_sol.u[Ts(t,x)][:, 3]
                tempTf_op = modeloptim_sol.u[Tf(t,x)][:, end-1]
                #return  append!(tempTs_op, tempTf_op)
                return tempTf_op
            end
        
            function loss(pguess, _)
                Temp = NLmodeloptim(xz_data, pguess)
                lossr = (Temp .- expdata).^2
                return sqrt(sum(lossr) / length(expdata)) #MSE
            end
    #loss([1.], [])
    
    display(loss(p0, []))

        optf = OptimizationFunction(loss, Optimization.AutoForwardDiff())
        lb = [5000., 700.]
        ub = [20000., 2000.]
        #lb = [90.]
        #ub = [7000.]
        optprob = Optimization.OptimizationProblem(optf, p0, [], lb=lb, ub=ub)
        optsol = solve(optprob, NLopt.GN_MLSL_LDS(), local_method = NLopt.LN_NELDERMEAD(), maxtime = 180, local_maxiters = 10000)
        println(optsol.retcode)
        
        
         
        
        pnew = optsol.u
        
        res_error = loss(pnew, [])
        display(res_error)
        modelfit = remake(prob, p = pnew)
        modelfit_sol = solve(modelfit, Rodas4(), saveat = xz_data, reltol=1e-12, abstol = 1e-12)
        
    
        psol = modelfit_sol
    plot1 =  plot(
        psol.t, 
        psol.u[Tf(t,x)][:,end-1], # around 135 mm in T3
        title = "Gas Temperature Profile T3 - exp. 67", 
        label = "optimized", 
        xlabel = "Time (s)", 
        ylabel = "Temperature (K)")
        scatter!(
            xz_data, 
            y1z_data,  
            label = "Experimental",
            xlabel = "Time (s)", 
            ylabel = "Temperature (K)")
                
plot2 = plot(
   psol.t, 
   psol.u[Ts(t,x)][:,3], # around 5 mm in T8
   title = "Solid Temperature Profile T8 - exp. 67", 
   label = "Numerical", 
   xlabel = "Time (s)", 
   ylabel = "Temperature (K)")
   scatter!(
   xz_data, 
   y2z_data,  
   label = "Experimental",
   xlabel = "Time (s)", 
   ylabel = "Temperature (K)")
   plot(plot1, plot2, layout=(2,1), size=(700, 700))
