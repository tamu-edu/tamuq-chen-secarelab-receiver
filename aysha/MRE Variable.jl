using ModelingToolkit, MethodOfLines, LinearAlgebra, OrdinaryDiffEq, DomainSets, DifferentialEquations

#### CONSTANTS ###
ks= 120/1000 #C SiC- #kW/m.K from the web (imetra.com))
kf = 0.06763/1000 #thermal conductivity of the fluid phase kW/m.K
ρs=3100 #kg/m3
ρf=1.225 #kg/m3
L= 137e-3 #m
Tamb= (22.448 + 273.15) #K (same for all exp) 
deltax = 0.00223 #discretization (m)
A_t= 324e-6 #m2 - for the whole receiver (18x18mm2)
A_st= 176e-6 #total area - solid m2
A_ft= 49e-6 # total area - fluid m2
A_channel = 2.25e-6 #m2 (1.5x1.5mm2) 
n_channel = 100
A_exchange = 162.5e-6 #m2 - contact area between fluid and solid
Vs =  A_st * deltax
Vf =  A_ft * deltax
qlpm = 7.12 #lpm
q= qlpm/(1000*60) #m3/s
V= q/(A_channel*n_channel) #m/s - using the area of the whole receiver (18x18mm2)
th_s = 0.4e-3 #m
th_f = 0.7e-3 #m
I0 = 456 #kW/m2
Q= I0 #kW/m2
ϵ= 0.8
σ= (5.17e-8)/1000 #kW/m2.K^4 
Nu = 3.657 
w = 1.5e-3 #width of channel (m)
Vi = w * w * deltax * n_channel #m3
Av = 1662 #specific area (m-1)
hext = 10/1000 #kW/m2.K
kins = 0.078/1000 #kW/m*K
r0 = 23/1000 #m
r = 42/1000 #m
Tins = 356 #K exp. 71


#### PDE SYSTEM ####
@parameters t x
@parameters hlocal
@variables Ts(..) Tf(..)
Dt = Differential(t) 
Dx = Differential(x)
Dxx = Differential(x)^2

x_max1 = L
x_min1 = 0.
t_min = 0.
t_max = 6737.
nc1 = 11
order = 2
x_num1 = range(x_min1, x_max1, length = nc1)
dx = (x_max1 - x_min1) / (nc1 - 1)
p = [hlocal => 0.98/1000]

Cps(T)= 1110. + 0.15*T- 425. * exp(-0.003*T)
Cpf(T)= 1.93e-10 *(T^3)+1.14e-3*(T^2)-4.49e-1*T+1.06e3

eq1 = [
    Vs * ρs * Dt(Ts(t,x)) ~ 1/Cps(Ts(t, x)) * (Vs * ks * Dxx(Ts(t, x)) + hlocal * Av * Vi * ((Ts(t, x)) - Tf(t, x)) - (kins * (r/r0) * (Ts(t, x)-Tins)) * A_t/ (r-r0)) ,
    Vf * ρf * Dt(Tf(t, x)) ~ 1/Cpf(Tf(t, x)) * (-Vf * Cpf(Tf(t, x)) * ρf * V * Dx(Tf(t, x)) + hlocal * Av * Vi * ((Ts(t, x) - Tf(t, x))))
    ]      
bcs1 = [
    Ts(0., x) ~ Tamb, # initial
    Tf(0., x) ~ Tamb, # initial
    -A_st * ks * Dx(Ts(t, x_max1)) ~ 0, # far right
    -A_st * ks * Dx(Ts(t, x_min1)) ~ I0 * A_st - ϵ * σ * A_st * (Ts(t, x_min1)^4 - Tamb^4) - hext * A_st * (Ts(t, x_min1) - Tamb),  # far left
    -A_ft * kf * Dx(Tf(t, x_max1)) ~ -ρf * Cpf(Tf(t, x)) * V * A_ft * (Tf(t, x_max1) - Tamb) # fluid exiting
    ] 
domains1 = [
    t ∈ Interval(t_min, t_max),
    x ∈ Interval(x_min1, x_max1)
    ]
@named pdesys = PDESystem(eq1, bcs1, domains1, [t, x], [Ts(t, x), Tf(t, x)], p)
discretization = MOLFiniteDifference([x => dx], t, approx_order=order)
prob = discretize(pdesys, discretization)