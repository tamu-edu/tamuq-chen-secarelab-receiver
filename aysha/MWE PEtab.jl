using ModelingToolkit, MethodOfLines, Plots, DomainSets, OrdinaryDiffEq, PEtab, DataFrames

pythonplot()

# Parameters, variables, and derivatives
@parameters k a To
@variables t x T(..)
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2

Tamb = 273 #K
To_test = 1000 #K
h = 30 #W/m^2K arbritrary value
L = 0.1 #m
ρ = 7500 #kg/m^3
cp = 420 #J/kgK
kin = 15 #W/m K
ain = kin/ρ/cp #m^2/s

param = [k => kin, a => ain, To => To_test]

# 1D PDE and boundary conditions
eq  = Dt(T(t, x)) ~ a * Dxx(T(t, x))
bcs = [T(0, x) ~ Tamb, #initial condition
        Dx(T(t, 0)) ~ -h/k *(T(t, 0) - Tamb), #cold end
        T(t, L) ~ To] #hot end

# Space and time domains
domains = [t ∈ Interval(0.0, 300.0),
           x ∈ Interval(0.0, L)]

# PDE system
@named pdesys = PDESystem(eq, bcs, domains, [t, x], [T(t, x)], param)

# Method of lines discretization
dx = L/100
order = 2
discretization = MOLFiniteDifference([x => dx], t)

# Convert the PDE problem into an ODE problem
prob = discretize(pdesys,discretization)


# Solve ODE problem
sol = solve(prob, Tsit5(), saveat=1.)

# Plot the solution
px = collect(sol[x])
pt = sol[t]
pT = sol[T(t,x)]
plt = plot(sol[t], sol[T(t,x)][:, 1], xlabel="Time", ylabel="Temperature", title="1D Heat Equation")
cnt = contourf(px, pt, pT, xlabel="Position", ylabel="Time", title="1D Heat Equation")
display(cnt)

#extract simulated experiments to create measurements for PEtab
idx =[3, 50, 97]
mx = px[idx]
mT = pT[:, idx]

#Create PEtab problem
odesys = symbolic_discretize(pdesys, discretization)
@unpack T, k, a, To = odesys[1]
obs_T3 = PEtabObservable(T[3], 0.5)
observables = Dict("obs_T3" => obs_T3)
_k = PEtabParameter(:k, lb=0.1, ub=100., scale=:lin)
_a = PEtabParameter(:a, lb=0.1, ub=100., scale=:lin)
params = [_k , _a]  
E0 = Dict[To => 1000.]
E1 = Dict[To => 900.]
sim_cond = Dict("c0" => E0, "c1" => E1)
measurements = DataFrame(
    sim_id = ["c0"],
    obs_id=["obs_T3"],
    time = [pt[end]],
    measurement=[400.])

petab_model = PEtabModel(odesys, simulation_conditions, observables, measurements,
    parameters, state_map=state_map, parameter_map=parameter_map,
    verbose=false)