### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 7f69353e-34af-11ee-075b-4961f0cbc6be
using MethodOfLines, ModelingToolkit, DomainSets, OrdinaryDiffEq, NonlinearSolve, PlutoUI, DifferentialEquations, Plots, XLSX, Statistics

# ╔═╡ 38b5f97c-eefe-48ba-936d-0c78096e2e94
begin 

	md" qlpm: $(@bind qlpm Slider(0.0:0.01:17.7, 4.52, true))"
	
end 

# ╔═╡ fd8f5434-7ea7-4b7e-80f7-eee48c308a56
begin 
	md" I0: $(@bind I0 Slider(0.0:0.01:700.0, 296.0, true))"
	
end 

# ╔═╡ 0f338119-c8cb-4d79-ab12-71119eaeb1b8
begin
#solid
#using exp. 028
ks= 120/1000 #C SiC- #kW/m.K from the web (imetra.com))
#ks=(52000*exp(-1.24e-5*T)/(T+437))/1000 #kW/m.K (Ali et al.)
#k=ks*(1-0.88) # effective axial heat conductivity for monolithic structures pg.285 (structured cayalysts and reactors book), using apparent emissivity pg.240 (thermal radiation heat transfer book)
kf = 0.06763/1000 #thermal conductivity of the fluid phase kW/m.K
#kf=(1.52e-11*(T^3)-4.86e-8*(T^2)+1.02e-4*T-3.93e-3)/1000 #kW/m.K (Ali et al.)
#keff=0.82*kf+(1-0.82)*ks #using "Thermal analysis and design of a volumetric solar absorber depending on the porosity" paper - note that I still need to measure the porosity of SiC 
#ρf= 3.018*exp(-0.00574*T)+0.8063*exp(-0.0008381*T) #kg/m3 (Roldan et al.)
#Cps= 750/1000 #KJ/kg.K SiC- from the web (imetra.com)
#Cps=1110+0.15*T-425*exp(-0.003*T)
#Cpf=1.93e-10*(T^3)+1.14e-3*(T^2)-4.49e-1*T+1.06e3
Cpf=700/1000 #kJ/kg.K - check
Cps=670/1000 #kJ/kg.K - check
ρs=3100 #kg/m3
ρf= 1.225 #kg/m3
L= 130e-3 #m
#α = keff/ρ*Cp #kW/m2.K 
Ts= (22.448 + 273.15) #K (same for all exp) 
T1= (127.931 + 273.15) #K - gas temp. initial for exp. 028
T8= (79.424 + 273.15) #K - solid temp.
A= 4e-4 #m2 - for the whole receiver
q= qlpm/(1000*60) #m3/s
V= q/(4e-6*81) #m/s - using the area of a single channel (2x2mm2)
#Qv= I0*exp(-1000*x)  #kW/m2 - (K extinction coefficient taken from Howell and Hendrick paper pg.86, measure pore diameter and fix)
Q= I0 #kW/m2
ϵ= 0.8
σ= (5.17e-8)/1000 #kW/m2.K^4 
h= 8/1000 #kW/m2.K - check value (free convection)
#(3.018 * exp(-0.00574 * T(t, x)) + 0.8063 * exp(-0.0008381 * T(t, x)))
end

# ╔═╡ 314b8647-4228-4deb-add4-080499726e23
V

# ╔═╡ a9f17892-f115-4bd3-9d90-8cb416d21602
begin
# Parameters, variables, and derivatives
@parameters t x
@variables T(..)

Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2

# Define your parameters as functions of time and space
param1(t, x) = 3.018 * exp(-0.00574 * T(t, x)) + 0.8063 * exp(-0.0008381 * T(t, x))
param2(t, x) = 1.93e-10 * (T(t, x)^3) + 1.14e-3 * (T(t, x)^2) - 4.49e-1 * T(t, x) + 1.06e3
param3(t, x) = (1.52e-11 * (T(t, x)^3) - 4.86e-8 * (T(t, x)^2) + 1.02e-4 * T(t, x) - 3.93e-3) / 1000

# MOL Discretization parameters
x_max = L
x_min = 0
t_min = 0
t_max = 6737 # 6737    # Experiment 28
nc = 20
x_num = range(x_min, x_max, length=(nc))
dx = (x_max - x_min) / (nc - 1)

# Define your equation using the parameters
eq = param1(t, x) * param2(t, x) * V * Dx(T(t, x)) ~ -param1(t, x) * param2(t, x) * Dt(T(t, x))

bcs = [T(t_min, x) ~ Ts,
       -(param3(t, x)* Dx(T(t, x_max))) ~ 0,
       -(param3(t, x) * Dx(T(t, x_min))) ~ -h * (T(t, x_min) - Ts) + Q]

# Space and time domain
domains = [t ∈ Interval(t_min, t_max),
           x ∈ Interval(x_min, x_max)]

# Define the PDE system
@named pdesys = PDESystem([eq], bcs, domains, [t, x], [T(t, x)])

# Solving using MOL
order = 2
discretization = MOLFiniteDifference([x => dx], t, approx_order=order)

# Convert the PDE problem into an ODE problem
prob = discretize(pdesys, discretization)

# Define initial conditions
u0 = [T(t_min, x) for x in x_num]

# Solve ODE problem
sol = solve(prob, saveat=3, u0=u0)

T_exp = sol.u[T(t, x)]


	
	md""" 
	## Numerical Solution
	"""
end 

# ╔═╡ 7664d53a-41be-4d01-aebd-0a8d4a5750e7
begin
#Exp 28
C = XLSX.readxlsx("/Users/aishamelhim/Documents/ResearchData/SolarSimulator/EXCEL/Data_FPT0028_230302_111716.xlsx")["Sheet 1 - Data_FPT0028_230302_1"]["A3:E6376"]
xc_data = C[:,1]
y1c_data = C[:,4]
y2c_data = C[:,5]
end

# ╔═╡ fbb27462-588d-4604-9f8c-44c49440e4f6
begin
plot(
	sol.t, 
	sol[T(t,x)][:,1], 
	title = " Temperature Profile", 
	label = "Numerical", 
	xlabel = "Time (s)", 
	ylabel = "Temperature (K)")
	#scatter!(
		#xc_data, 
		#y1c_data,  
		#label = "Experimental",
		#xlabel = "Time (s)", 
		#ylabel = "Temperature (K)")
end

# ╔═╡ 06661fef-3428-4607-82dc-f80792e3f6b0
begin
# Access the temperature solution at a specific time point 
y_min=400.0
y_max=1600.0
time_point = 2000
idx_time_point = argmin(abs.(sol.t .- time_point))  # Find the index closest to the specified time point

# Collect the temperature solution at all spatial points for the specified time point
T_exp_at_time_array = [sol[T(t,x)][idx_time_point, i] for i in 1:length(x_num)]

# Create a 1D plot of T(x) at the specified time point
plot(x_num, T_exp_at_time_array, xlabel="Position x (m)", ylabel="Temperature T (K)", label="Temperature Profile at t = $time_point s", linewidth=2, ylim=(y_min,y_max))
end

# ╔═╡ Cell order:
# ╠═7f69353e-34af-11ee-075b-4961f0cbc6be
# ╠═0f338119-c8cb-4d79-ab12-71119eaeb1b8
# ╠═314b8647-4228-4deb-add4-080499726e23
# ╠═38b5f97c-eefe-48ba-936d-0c78096e2e94
# ╠═fd8f5434-7ea7-4b7e-80f7-eee48c308a56
# ╠═a9f17892-f115-4bd3-9d90-8cb416d21602
# ╠═7664d53a-41be-4d01-aebd-0a8d4a5750e7
# ╠═fbb27462-588d-4604-9f8c-44c49440e4f6
# ╠═06661fef-3428-4607-82dc-f80792e3f6b0
