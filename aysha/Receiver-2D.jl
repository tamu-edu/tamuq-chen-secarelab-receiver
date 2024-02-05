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

# ╔═╡ 6c3b531e-cb51-4eb7-9223-694583fcc592
begin
import Pkg
Pkg.activate("/Users/aishamelhim/Documents/ResearchData")
end

# ╔═╡ 68934758-bb41-11ed-3036-c1066cced178
using MethodOfLines, ModelingToolkit, DomainSets, OrdinaryDiffEq, NonlinearSolve, PlutoUI, DifferentialEquations, Plots

# ╔═╡ 4d3f4033-6019-456f-b1c9-5c741de0e088
Pkg.status()

# ╔═╡ 7c12f7cb-5c00-4e36-a757-23ff0e771abe
begin 

	md" qlpm: $(@bind qlpm Slider(0.0:0.01:17.7, 4.52, true))"
	
end 

# ╔═╡ 203270e0-5f86-452f-9cba-94d6934bc350
begin 
	md" I0: $(@bind I0 Slider(0.0:0.01:700.0, 296.0, true))"
	
end 

# ╔═╡ c1b95018-c3cd-4b13-bedb-c5b7bdccc549
begin
#using exp. 028
k= 40/1000 #kW/m.K - check value
ρ= 1.293 #kg/m3
Cp= 700/1000 #kJ/kg.K
R1= 1.5e-3 #m - assumed
L= 130e-3 #m
α = k/ρ*Cp #W/m2.K - assumed value #35
Ts= (22.448 + 273.15) #K (same for all exp) 
T1= (127.931 + 273.15) #K - gas temp. initial for exp. 028
T8= (79.424 + 273.15) #K - solid temp.
A= 2*pi*R1*L+2*pi*R1^2  #m2
q= qlpm*1.67e-5 #m3/s
V= q/A #m/s
Qv= A*I0*exp(-2300*L)  #kW
Q= I0*A #kW
ϵ= 0.8
σ= (5.17e-8)/1000 #kW/m2.K^4 
h= 8/1000 #kW/m2.K - check value (free convection)
end

# ╔═╡ 6d05e7ef-7167-499f-a374-e6e93fc4707a
Q

# ╔═╡ d82bae0a-4368-4036-840b-4b0d5a281879
Qv

# ╔═╡ 24e53b97-8ce8-4c53-918f-2e9aa42661f8
begin
@parameters t, r, x
@variables T(..)
Dt= Differential(t)
Dr=Differential(r)
Drr=Differential(r)^2
Dx=Differential(x)
Dxx=Differential(x)^2


#Discretization parameters
x_max= L
x_min= 0
t_max= 6373
t_min= 0
r_max= R1
r_min= 0
nr1 = 30
r1= range(r_min, r_max, length=(nr1))
dr=(r_max-r_min)/(nr1-1)
x1= range(x_min, x_max, length=(nr1))
dx= (x_max-x_min)/(nr1-1)
#equation
#eq=[Dt(T(t,r,x)) ~ -V*Dx(T(t,r,x)) + (1/r)*Dr(α*r*Dr(T(t,r,x))) + α*Dxx(T(t,r,x)) + (Qv/ρ*Cp)]

eq=(1/r)*Dr(α*r*Dr(T(t,r,x))) + (α)*Dxx(T(t,r,x)) + (Qv/ρ*Cp) ~ Dt(T(t,r,x)) + V*Dx(T(t,r,x))

#BC
bcs=[T(t_min,r,x)~T1,
Dr(T(t,r_min,x))~ 0,
Dr(T(t,r_max,x))~0,
k*Dx(T(t,r,x_min)~ Q - h*A*((T(t,r,x_min))-Ts) - A*ϵ*σ*((T(t,r,x_min))^4-(Ts)^4)),
k*Dx(T(t,r,x_max))~ 0]


#space and time domain
domains=[t ∈ Interval(t_min, t_max),
r ∈ Interval(r_min, r_max),  
x ∈ Interval(x_min, x_max)]

@named pdesys = PDESystem([eq], bcs, domains, [t,r,x], [T(t,r,x)])

# Solving using MOL
order = 2
discretization = MOLFiniteDifference([r=> dr, x=>dx], t, approx_order = order)
# Convert the PDE problem into an ODE problem
prob = discretize(pdesys, discretization)
# Solve ODE problem
sol = solve(prob, saveat=3)
T_exp=sol.u[T(t,r,x)]

end

# ╔═╡ f2482038-c9b3-4dc9-8007-c96f8389c3db
begin

r_sol= sol[T(r)]
t_sol=sol[T(t)]
x_sol=sol[T(x)]
	
R = sol[r]
T_ = sol[t]
X = sol[x]
end 

# ╔═╡ Cell order:
# ╠═6c3b531e-cb51-4eb7-9223-694583fcc592
# ╠═4d3f4033-6019-456f-b1c9-5c741de0e088
# ╠═68934758-bb41-11ed-3036-c1066cced178
# ╠═c1b95018-c3cd-4b13-bedb-c5b7bdccc549
# ╠═6d05e7ef-7167-499f-a374-e6e93fc4707a
# ╠═d82bae0a-4368-4036-840b-4b0d5a281879
# ╠═7c12f7cb-5c00-4e36-a757-23ff0e771abe
# ╠═203270e0-5f86-452f-9cba-94d6934bc350
# ╠═24e53b97-8ce8-4c53-918f-2e9aa42661f8
# ╠═f2482038-c9b3-4dc9-8007-c96f8389c3db
