using ModelingToolkit, MethodOfLines, Plots, DomainSets, OrdinaryDiffEq, DataFrames
using Optimization, OptimizationNLopt
using SciMLStructures: replace!, Tunable
using SymbolicIndexingInterface: parameter_values, state_values


# Parameters, variables, and derivatives
@parameters k aa
@independent_variables t x
@variables T(..)

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
aain = 0.3 #arbirtrary corrections value
af(x) = (x)/ρ/cp #m^2/s
af2(x) = (x/To_test)^aa
#@register_symbolic af2(x) #is not working??? but it is not needed

param = Dict([k => kin, aa => aain])
p_names = [i for i in keys(param)]

# 1D PDE and boundary conditions
eq  = Dt(T(t, x)) ~ (af(k) * af2(T(t,x)))Dxx(T(t, x))
bcs = [T(0, x) ~ Tamb, #initial condition
        Dx(T(t, 0)) ~ -h/k *(T(t, 0) - Tamb), #cold end
        T(t, L) ~ To_test] #hot end

# Space and time domains
domains = [t ∈ Interval(0.0, 300.0),
           x ∈ Interval(0.0, L)]

# PDE system
@named pdesys = PDESystem(eq, bcs, domains, [t, x], [T(t, x)], p_names; defaults = param)

# Method of lines discretization
dx = L/100
order = 2
discretization = MOLFiniteDifference([x => dx], t)

# Convert the PDE problem into an ODE problem
prob = discretize(pdesys,discretization)

# Solve ODE problem
timesteps = 0.0:1.:300.0
sol = solve(prob, FBDF(), saveat=timesteps, reltol=1e-12, abstol=1e-12)

# Plot the solution
px = collect(sol[x])
pt = sol[t]
pT = sol[T(t,x)]
plt = plot(sol[t], sol[T(t,x)][:, 1], xlabel="Time", ylabel="Temperature", title="1D Heat Equation")
cnt = contourf(px, pt, pT, xlabel="Position", ylabel="Time", title="1D Heat Equation")
display(cnt)

#extract simulated experiments to create measurements, pick 20 random points from px
n = 3
idx = rand(1:length(px), n)
mx = px[idx]
mT = pT[:, idx] #.* rand(size(mT)[1],size(mT)[2])

#Create Optimization problem
function loss(p_fit, p)
    odeprob = p[1] # ODEProblem stored as parameters to avoid using global variables
    newprob = remake(odeprob) # copy the problem with `remake`
    # update the parameter values in-place
    replace!(Tunable(), parameter_values(newprob), p_fit)
    timesteps = p[2]
    sol = solve(prob, FBDF(), saveat=timesteps, reltol=1e-12, abstol=1e-12)
    truth_id = p[3]
    truth_dt = p[4]
    res = sol[T(t,x)][:, truth_id]
    return sum((truth_dt .- res) .^ 2) / length(truth_dt)
end
p0 = [2., 0.5]
lb = [1.0, 0.1]
ub = [1000.0, 0.8]

#test loss function
loss([param[k]], (prob, timesteps, idx, mT))
loss(p0, (prob, timesteps, idx, mT))

optfn = OptimizationFunction(loss, Optimization.AutoForwardDiff())
# parameter object is a tuple, to store differently typed objects together
optprob = Optimization.OptimizationProblem(optfn, p0, (prob, timesteps, idx, mT), lb=lb, ub=ub)
optsol = solve(optprob, NLopt.GN_MLSL_LDS(), local_method=NLopt.LN_NELDERMEAD(), maxtime=30, local_maxiters=100000)
optsol.objective

newprob = remake(prob) # copy the problem with `remake`
replace!(Tunable(), parameter_values(newprob), optsol.minimizer)
scatter(sol[t], sol[T(t,x)][:, 1], xlabel="Time", ylabel="Temperature", title="1D Heat Equation", label="data")
sol = solve(prob, FBDF(), saveat=timesteps, reltol=1e-12, abstol=1e-12)
plot!(sol[t], sol[T(t,x)][:, 1], label="fit", lw=3)

