begin #libraries
    using MethodOfLines, PEtab, DataFrames
    using ModelingToolkit
    using DomainSets, OrdinaryDiffEq
    using NonlinearSolve, DifferentialEquations
    using Plots, XLSX, Statistics, Symbolics, Interpolations
    using Optim, LsqFit
    using Optimization, OptimizationNLopt, Symbolics, OptimizationOptimJL, ForwardDiff, OptimizationMOI
end
begin #define parameters
    L= 137e-3 #m
    #α = keff/ρ*Cp #kW/m2.K 
    Tamb = (22.448 + 273.15) #K (same for all exp) 
    deltax = 0.002795918367346939 #discretization (m)
    w_t = 19/1000 #m
    A_t= w_t*w_t #m2 - for the whole receiver (19x19mm2)
    A_st= 176e-6 #total area - solid m2
    A_ft= 49e-6 # total area - fluid m2
    channel_w = 1.5/1000 #m
    A_channel = channel_w * channel_w  #m2 (1.5x1.5mm2) 
    n_channel = 100
    A_exchange = 162.5e-6 #m2 - contact area between fluid and solid
    Vs =  A_st * L #* deltax
    Vf =  A_ft * L #* deltax
    qlpm = 7.12 #lpm
    q= qlpm/(1000*60) #m3/s
    #V= q/(A_channel*n_channel) #m/s - using the area of the whole receiver (18x18mm2)
    #V = 0.57 #m/s (calculated from excel sheet and COMSOL)
    th_s = 0.4e-3 #m
    th_f = 0.7e-3 #m
    I0 = 456 #kW/m2
    #Qv= I0*exp(-1000*x)  #kW/m2 - (K extinction coefficient taken from Howell and Hendrick paper pg.86, measure pore diameter and fix)
    #Q= I0 #kW/m2
    ϵ= 0.8
    σ= (5.17e-8)/1000 #kW/m2.K^4 Stefan-Boltzmann constant
    Lc = 4*(w_t*w_t)/(4*w_t)
    kf = (0.056/1000) #kW/m.K
    ρf = 0.5 #kg/m3
    Cpf = (1090/1000) #kJ/kg.K
    mu = 2.0921e-5 #Pa.s
    Re = (ρf*V*w_t)/mu
    Pr = (Cpf*mu)/kf
    Gz = (w_t/L)*Re*Pr
    # A = 3.657
    # B = 0.5272
    # C = 55.1239
    # n = 0.3056
    #Nu = A*(1+(B*((Gz)^n)*exp(-C/Gz)))
    #Nu = 3.657 
    Lc = 4*(w_t*w_t)/(4*w_t)
    w = 1.5e-3 #width of channel (m)
    Vi = w * w * n_channel * L #m3
    Av = 4*(w*L) / (w^2*L) #specific area (m-1)
    hext = 10/1000 #kW/m2.K
    kins = 0.078/1000 #kW/m*K
    r0 = 23/1000 #m
    r = 42/1000 #m
    #Tins = 356 #K exp. 71
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
        @parameters v I0 ρsCps A B n C #psCps ks h_average
        @variables Ts(..) Tf(..)
        Dt = Differential(t) 
        Dx = Differential(x)
        Dxx = Differential(x)^2
        e = 0.62
        ks = (37/1000)*(1-e) #kW/m.K
        #ρs = 3100*(1-e) #kg/m3
        #Cps = (1225/1000)*(1-e) #kJ/kg

        #p = [psCps => 90000., ks=> 37, h_average => (Nu*kf)/w_t]
        p = [ρsCps => 80000., A => 20., B =>0.7, C => 40., n => 0.1, V => 0.57, I0=>456]
        Nu = A*(1+(B*((Gz_f(x))^n)*exp(-C/Gz_f(x))))
        nu = 4.364*(1+(0.7*((Gz_f(0.134))^10)*exp(-40/Gz_f(0.134))))
        h_average = (Nu*kf)/Lc


        # MOL Discretization parameters for system 1
        x_max1 = L #T3
        x_min1 = 0.
        x1 = 5/1000 #m #T8
        x2 = 55/1000 #m #T9
        x3 = 106/1000 #m #T10
        x4 = 82/1000 #m #T11
        x5 = 27/1000 #m #T12
        t_min = 0.
        t_max = 7084.
        
        nc1 = 100
        
        x_num1 = range(x_min1, x_max1, length = nc1)
        
        
        dx = (x_max1 - x_min1) / (nc1 - 1)
        
        
        # PDE equation for system 1
        
        eq1 = [
               Vs * (ρsCps/1000) * Dt(Ts(t,x)) ~ Vs * (ks) * Dxx(Ts(t,x)) - ((h_average/1000) * Av * Vi * ((Ts(t,x)) - Tf(t,x))) .- (kins * (r/r0) .* (Ts(t,x) .- Tins_f(t)) * A_t / (r-r0)),
               Vf* ρf * Cpf * Dt(Tf(t,x)) ~ Vf * kf * Dxx(Tf(t,x)) - Vf * ρf * Cpf * V * Dx(Tf(t,x)) + (h_average/1000) * Av * Vi * ((Ts(t,x) - Tf(t,x)))
               ]
              
        bcs1 = [
            Ts(0., x) ~ Tamb, # initial
            Tf(0., x) ~ Tamb, # initial
            -A_st * (ks) * Dx(Ts(t, x_max1)) ~ 0.0, # far right
            -A_st * (ks) * Dx(Ts(t, x_min1)) ~ I0 * A_st - ϵ * σ * A_st * (Ts(t,x_min1)^4 - Tamb^4) - hext * A_st * (Ts(t, x_min1) - Tamb),  # far left
            -A_ft * kf * Dx(Tf(t, x_max1)) ~ 0.0, #-ρf * Cpf * V * A_ft * (Tf(t, x_max1) - Tamb), # exiting fluid
            -A_ft * kf * Dx(Tf(t, x_min1)) ~ ρf * Cpf * V * A_ft * (Tf(t,x_min1)- Tamb) # entering fluid (upstream temperature)
              ] 
        # Space and time domain for system 1
        domains1 = [t ∈ Interval(t_min, t_max),
                    x ∈ Interval(x_min1, x_max1)]
        
        # ODE system for system 1
        @named pdesys = PDESystem(eq1, bcs1, domains1, [t, x], [Ts(t, x), Tf(t,x)], p)
        
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
        #using PEtab
        obs_T8 = PEtabObservable(Ts(t, x1))
        obs_T9 = PEtabObservable(Ts(t, x2))
        obs_T10 = PEtabObservable(Ts(t, x3))
        obs_T11 = PEtabObservable((Tf(t, x4) + Ts(t, x4)) /2)
        obs_T12 = PEtabObservable((Tf(t, x5) + Ts(t, x5)) /2)
        obs_T3 = PEtabObservable(Tf(t, x_max1))
    
        observables = Dict("obs_T8" => obs_T8, "obs_T9" => obs_T9, "obs_T10" => obs_T10, 
        "obs_T11" => obs_T11, "obs_T12" => obs_T12, "obs_T3" => obs_T3)

            _ρsCps = PEtabParameter(ρsCps, lb=1e-3, ub=1e3, scale=:log10)
            _A = PEtabParameter(A, lb=1e-3, ub=1e3, scale=:log10)
            _B = PEtabParameter(B, lb=1e-3, ub=1e3, scale=:log10)
            _n = PEtabParameter(n, lb=1e-3, ub=1e3, scale=:log10)
            _C = PEtabParameter(C, lb=1e-3, ub=1e3, scale=:log10)
            
            parameters = [_ρsCps, _A, _B, _n, _C]
            
            #Defining simulation conditions
            begin
                condition_E67 = Dict(:I0 => 456, :v => 1.22)
                condition_E68 = Dict(:I0 => 456, :v => 1.00)
                condition_E69 = Dict(:I0 => 456, :v => 0.84)
                condition_E70 = Dict(:I0 => 456, :v => 0.73)
                condition_E71 = Dict(:I0 => 456, :v => 0.57)
                condition_E72 = Dict(:I0 => 304, :v => 1.46)
                condition_E73 = Dict(:I0 => 304, :v => 1.05)
                condition_E74 = Dict(:I0 => 304, :v => 0.72)
                condition_E75 = Dict(:I0 => 304, :v => 0.55)
                condition_E76 = Dict(:I0 => 304, :v => 0.36)
                condition_E77 = Dict(:I0 => 256, :v => 1.10)
                condition_E78 = Dict(:I0 => 256, :v => 0.80)
                condition_E79 = Dict(:I0 => 256, :v => 0.64)
                condition_E80 = Dict(:I0 => 256, :v => 0.53)
                condition_E81 = Dict(:I0 => 256, :v => 0.36)
            end
            
            begin
            simulation_conditions = Dict("E67" => condition_E67, "E68" => condition_E68, 
            "E69" => condition_E69, "E70" => condition_E70,
            "E71" => condition_E71,"E72" => condition_E72, 
            "E73" => condition_E73, "E74" => condition_E74, 
            "E75" => condition_E75, "E76" => condition_E76,
            "E77" => condition_E77,"E78" => condition_E78, 
            "E79" => condition_E79, "E80" => condition_E80,
            "E81" => condition_E81)
            end

            #Defining measurement data
            measurements = DataFrame(
            simulation_id = ["Ts1_T8", "Ts1_T9", "Ts1_T10", "Ts1_T11", "Ts1_T12", "Tf1_T3",
            "Ts2_T8", "Ts2_T9", "Ts2_T10", "Ts2_T11", "Ts2_T12", "Tf2_T3",
            "Ts3_T8", "Ts3_T9", "Ts3_T10", "Ts3_T11", "Ts3_T12", "Tf3_T3",
            "Ts4_T8", "Ts4_T9", "Ts4_T10", "Ts4_T11", "Ts4_T12", "Tf4_T3",
            "Ts5_T8", "Ts5_T9", "Ts5_T10", "Ts5_T11", "Ts5_T12", "Tf5_T3",
            "Ts6_T8", "Ts6_T9", "Ts6_T10", "Ts6_T11", "Ts6_T12", "Tf6_T3",
            "Ts7_T8", "Ts7_T9", "Ts7_T10", "Ts7_T11", "Ts7_T12", "Tf7_T3",
            "Ts8_T8", "Ts8_T9", "Ts8_T10", "Ts8_T11", "Ts8_T12", "Tf8_T3",
            "Ts9_T8", "Ts9_T9", "Ts9_T10", "Ts9_T11", "Ts9_T12", "Tf9_T3",
            "Ts10_T8", "Ts10_T9", "Ts10_T10", "Ts10_T11", "Ts10_T12", "Tf10_T3",
            "Ts11_T8", "Ts11_T9", "Ts11_T10", "Ts11_T11", "Ts11_T12", "Tf11_T3",
            "Ts12_T8", "Ts12_T9", "Ts12_T10", "Ts12_T11", "Ts12_T12", "Tf12_T3",
            "Ts13_T8", "Ts13_T9", "Ts13_T10", "Ts13_T11", "Ts13_T12", "Tf13_T3",
            "Ts14_T8", "Ts14_T9", "Ts14_T10", "Ts14_T11", "Ts14_T12", "Tf14_T3",
            "Ts15_T8", "Ts15_T9", "Ts15_T10", "Ts15_T11", "Ts15_T12", "Tf15_T3"],
            time = [3929, 3929, 3929, 3929, 3929, 3929,
            5362, 5362, 5362, 5362, 5362, 5362,
            5363, 5363, 5363, 5363, 5363, 5363,
            6702, 6702, 6702, 6702, 6702, 6702,
            7084, 7084, 7084, 7084, 7084, 7084,
            3214, 3214, 3214, 3214, 3214, 3214,
            4572, 4572, 4572, 4572, 4572, 4572,
            6015, 6015, 6015, 6015, 6015, 6015, 
            6351, 6351, 6351, 6351, 6351, 6351,
            7144, 7144, 7144, 7144, 7144, 7144, 
            3041, 3041, 3041, 3041, 3041, 3041,
            5381, 5381, 5381, 5381, 5381, 5381,
            5230, 5230, 5230, 5230, 5230, 5230,
            5811, 5811, 5811, 5811, 5811, 5811, 
            5986, 5986, 5986, 5986, 5986, 5986],
            temperatures = [965.407, 975.144, 825.592, 880.867, 1004.165, 763.859,
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


            petab_model = PEtabModel(pdesys, simulation_conditions, observables, measurements,
                         parameters, parameter_map = p, verbose=false)
            petab_problem = PEtabODEProblem(petab_model, verbose=false)

        #Optimization
        petab_problem_opt= PEtab.OptimizationProblem(petab_problem;
        interior_point_alg=true,
        box_constraints=true)

        sol = solve(petab_problem_opt, Optim.ParticleSwarm(); reltol=1e-8)


 