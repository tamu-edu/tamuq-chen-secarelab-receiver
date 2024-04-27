using ModelingToolkit, DifferentialEquations, Plots
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters ks hf #to be fitted
@parameters qlpm Io Tins#varying conditions
@variables Tf(t) Ts(t)

L = 137e-3 #m
Tamb = (22.448 + 273.15) #K (same for all exp) 
w_t = 19.e-3 #m
A_frt = w_t * w_t #m2 - for the whole receiver (19x19mm2)
A_s_p = w_t * L * 4 #total area solid periphery m2
w_chnl = 1.5e-3 #m
A_chnl_p = w_chnl * L * 4  #m2 channel periphery
A_chnl_frt = w_chnl * w_chnl #m2 channel front	
n_chnl = 8 * 8
A_exchange = A_chnl_p * n_chnl #m2 - contact area between fluid and solid
Vs = A_frt * L #* deltax
Vf = n_chnl * A_chnl_frt * L #* deltax
ε = Vf / Vs #porosity
q = qlpm / (1000 * 60) #m3/s
ρf = 1.0 #kg/m3
mf = q * ρf #kg/s

th_s = 0.4e-3 #m
th_f = 0.7e-3 #m
em = 0.8
σ = 5.17e-8 #W/m2.K^4 Stefan-Boltzmann constant
kf = 0.056 #W/m.K
k_eff = 1.97 * (1 - ε)^1.5 * ks #Archie's Law
Cpf = 1090 #J/kg.K
ρCp_s_0 = 3200 * 1290  #J/kg.K * kg/m3
ks_0 = 0.1165 * 4184 #W/m.K
hf_0 = 500 #W/m2.K
μ = 2.0921e-5 #Pa.s

h_ext = 10 #W/m2.K natural convection
kins = 0.078 #W/m*K thermal conductivity of insulation
r_H_chnl = 4 * (w_chnl * w_chnl) / (4 * w_chnl) #hydraulic channel diameter
r_H = 4 * (w_t * w_t) / (4 * w_t) #hydraulic receiver diameter
r_ins = 42e-3 #m location of insulation thermocouple

ρCp_s = 3290.0 * 34.  * 4184. #J/kg.K kg/m3

#ρCp_sf(T) = 3290.0 * (0.27 + 0.135E-4 * T - 9720.0 * T^-2 + 0.204E-7 * T^2) * 4184.
#@register_symbolic ρCp_sf(T)   #J/kg.K kg/m3

eq1 = [D(Ts) ~ 1 / ((1 - ε) * ρCp_s * Vs) * (Io * A_frt - kins * (Ts - Tins) * A_s_p / (r_ins - r_H) - h_ext * A_frt * (Ts - Tamb) - em * σ * A_frt * (Ts^4 - Tamb^4) - hf * A_exchange * (Ts - Tf)),
    ε * ρf * Cpf * Vf * D(Tf) ~ hf * A_exchange * (Ts - Tf) - mf * (Tf - Tamb)
]

u0 = [Ts => Tamb, Tf => Tamb]

state_param = [qlpm => 15.0, Io => 300.0 * 1e3, Tins => (40.0 + 273.15)]
fit_param = [ks => ks_0, hf => h_ext]
p = vcat(state_param, fit_param)
tspan = (0, 3600.0)

@mtkbuild odes = ODESystem(eq1, t)
prob = ODEProblem(odes, u0, tspan, p)