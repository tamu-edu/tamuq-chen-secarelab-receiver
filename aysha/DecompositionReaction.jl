using CSV, DataFrames, Plots
using Catalyst, DifferentialEquations, IfElse, DiffEqCallbacks
using StatsBase, Optimization, OptimizationOptimJL, OptimizationOptimisers
using OptimizationBBO, OptimizationNLopt
using LsqFit

#### Functions ####
function positive_domain()
    condition(u,t,integrator) = minimum(u) < 0
    affect!(integrator) = integrator.u .= integrator.uprev
    return DiscreteCallback(condition,affect!,save_positions = (false,false)) 
end
function derivf(time, f) 
	deriv = zeros(length(time))
	j=1
	for i in 2:length(time)-1
		xᵢ=time[i]
		xᵢ₋₁=time[i-1]
		xᵢ₊₁=time[i+1]
		deriv[j] = f[i-1]* (2*xᵢ - xᵢ - xᵢ₊₁)/(xᵢ₋₁ - xᵢ)/(xᵢ₋₁ - xᵢ₊₁) + f[i]* (2*xᵢ - xᵢ₋₁ - xᵢ₊₁)/(xᵢ - xᵢ₋₁)/(xᵢ - xᵢ₊₁) + f[i+1]* (2*xᵢ - xᵢ - xᵢ₋₁)/(xᵢ₊₁ - xᵢ₋₁)/(xᵢ₊₁ - xᵢ)
		j+=1
	end
	return deriv
end
function fwModel(pnew)
	rmprob = remake(oprob, p = pnew)
	rmsol  = solve(rmprob, Tsit5(), saveat = TG_time, reltol=1E-12, abstol=1E-12, 
		callback = PositiveDomain([1.0]; save = false, abstol = 1E-12));
	all_moles = Array(rmsol(rmsol.t)) #every specie per row
	
	mass = all_moles .* MW
	mass = mass .* solids
	mass = transpose(sum(mass, dims = 1))
	mass_init = mass[1] #mass at the first step, 1000ºC
	mass_end = mass[end] #mass at the end of isotherm, 1850ºC
	mass = (mass .- mass_init) #./ (mass_init)
	return rmsol, rmsol.t, mass[:,1]
end
function loss(pguess, p)
	_, _, mass = fwModel(pguess)
	lossr = sum(abs2, mass .- TG_mass) / abs(mean(TG_mass))
	return lossr
end
function lossd(pguess, p)
	_, time, mass = fwModel(pguess)
	dmass = derivf(time, mass) ./1000
	dTG_mass = derivf(TG_time, TG_mass) ./1000
	lossr = sum(abs2, dmass .- dTG_mass) / abs(mean(dTG_mass))
	return lossr
end
function resample(time, mass; span=0.1, filter=0.05)
	own_d = derivf(time, mass)	
	model = loess(time, own_d, span = span)
	new_d = predict(model, time)
	filter = maximum(abs.(new_d)) * 0.05
	mask_deriv = (abs.(new_d) .>  filter)
	ind_d = findall(mask_deriv)
	push!(ind_d, length(new_d))
	pushfirst!(ind_d, 1)
	d_Time = time[ind_d]
	d_TG = mass[ind_d]
	return [d_Time, d_TG]
end

#E062 - 5C - 9.23 mg @ 1000C
#E070 - 7C - 8.35 mg @ 1000C
#E061 - 10C - 9.01 mg @ 1000C
#E036c - 20C - 10.19 mg @ 1000C

#### Experimental Data ####
T_init = 1000 #oC
time_compl = 160. #min
mass_rxn = 10.19 #@ ~1000ºC
hr = 20 #oC/min
T_iso = 1850. #oC
area = 2 * pi * (3E-3)^2 #m2 indicative area of crucible
ΔHCa = 178E3 #J/mol sublimation enthalpy for Ca
frc = 0.8 #volume fraction of CaC2ᵥ vs CaO
impurities = mass_rxn * 0.

file="TGA/SHELL_CaC2/E036c.txt"
Exp = CSV.read(file, DataFrame, header=["time","TF","TG", "dTG"], skipto=14, ntasks = 1)

#select 1000ºC till end of reaction
mask = (Exp.TF .>= T_init) .& (Exp.time.<= time_compl)
indices = findall(mask)
data = Exp[indices[1]:indices[end],:]

#remove impurities
data.TG .+= impurities

mass_init = -data.TG[1] #in mg
mass_end = -data.TG[end]
data.a =  (data.TG .+ mass_init) ./ (mass_rxn - mass_init) 
data.mass = (data.TG .+ mass_init)

#reducing data points
res_rate = 60 #seconds
TG_mass = [data.mass[i] for i=1:res_rate:length(indices)]
TG_a = [data.a[i]  for i=1:res_rate:length(indices)]
TG_time = [data.time[i] for i=1:res_rate:length(indices)]  .- data.time[1]
TG_temp = [data.TF[i] for i=1:res_rate:length(indices)] 

CSV.write("TGA/E036c.csv", DataFrame(time=TG_time, mass=TG_mass, a=TG_a .+1, temp=TG_temp))


#### Reaction Model ####
R = 8.314 #J/K mol - Universal Gas Constant

MW = [64., 40., 12., 56., 28.] #g/mol or mg/mmol
solids = [1, 1, 1, 1, 0] #1 if solid, 0 if gas
@variables t θ(t)
@parameters Ad Ed Ac Ec Asub
@species CaC2(t) Ca(t) C(t) CO(t) CaO(t)
D = Differential(t)
Ed = 598.2 #kJ/mol #isoconversional activation energy - 598
Ec = 903.2 #kJ/mol #isoconversional activation energy - 903
# Ad = 2.86
# Ac = 2.55
time_step = (T_iso-T_init)/hr #minus the starting time
T(t) = IfElse.ifelse(t<time_step, hr*t+T_init+273, T_iso+273) #Kelvin
kd(t) = (Ad*1E14) * exp(-(Ed*1E3)/R/(T(t)))
kc(t) = (Ac*1E27) * exp(-(Ec*1E3)/R/(T(t)))
# ksub(t) = Asub
decomp = @reaction	kd(t), CaC2 --> Ca + 2C
calc = @reaction kc(t), CaC2 + 2CaO --> 3Ca + 2CO
sublim = @reaction $Asub, Ca --> 0
#srf = @reaction	ks(t), CaC2ₛ --> Ca + 2C
#eq = D(θ) ~ Kθ * (Cas * θ + 1)
@named rs = ReactionSystem([decomp, calc, sublim], t)
#osys = convert(ODESystem, rs)

p = (Ad => 2.86, Ac => 2.55, Asub => 0.095) #10 ºC/min, #(Ad => 0.63, Ac => 0.70, Asub => 0.23) - 20 ºC/min, #(Ad => 2.58, Ac => 2.54, Asub => 0.212) - 5 ºC/min

tspan = (0., TG_time[end]) #min	
u0 = [CaC2 => mass_rxn*frc/MW[1], Ca => 0.0, C => 0.0, CO => 0.0, CaO => mass_rxn*(1-frc)/MW[4]] #mmol
#osys = convert(ODESystem, rs)
#osyss = structural_simplify(osys)
oprob = ODEProblem(rs, u0, tspan, p)

#### Test Solve Data ####
osol = solve(oprob, Tsit5(), saveat = TG_time, reltol=1E-12, abstol=1E-12, 
	callback = PositiveDomain([0.0]; save = false, abstol = 1E-12))
time = osol.t

p0 = [v for (k,v) in p]
display(loss(p0, []))
_, time, mass = fwModel(p0)
plot(osol)
plot(time, mass)

#### Optimization
# optf = OptimizationFunction(loss, Optimization.AutoForwardDiff())
# lb = [0.01, 0.01, 0.01]
# ub = [100., 100., 10.]
# optprob = OptimizationProblem(optf, p0, [], lb=lb, ub=ub)
# optsol = solve(optprob, NLopt.GN_MLSL_LDS(), local_method = NLopt.LN_NELDERMEAD(), maxtime = 10.0,
# local_maxiters = 10000)
#  #optsol = solve(optprob, NLopt.LN_NELDERMEAD())
# println(optsol.retcode)
# pnew = optsol.u
# res_error = loss(pnew, [])
println("HR = 20ºC/min, Ed = $Ed, Ec = $Ec,  $(p[1]), $(p[2]), $(p[3])") #Ad = $Ad, Ac = $Ac,
# println("pnew = $pnew")
# println("Residual = $res_error")

# rmprob = remake(oprob, p = optsol.u)
# rmsol  = solve(rmprob, Tsit5(), saveat = TG_time, reltol=1E-12, abstol=1E-12, 
# callback = PositiveDomain([0.0]; abstol = 1E-12, save=false))

#### Plotting ####
psol = osol #change to rmsol for the optimized model
sol, time, mass = fwModel(p0) #change to pnew for the optimized model

plt = plot(psol, width=2, title="Reaction System - Commercial CaC2", xlabel="Time (min)", ylabel="Species (molar)")
plot!(twinx(), psol.t, T.(psol.t) .-273, legend=:bottomright, label="T (ºC)", color=:orange, ls=:dash)
#display("image/png", plt)

#scatter(TG_time, TG_mass, label="Experimental")
plt2 = plot(ylabel="Absolute Mass Change (mg)", xlabel="Time (min)", ylims=(minimum(TG_mass)*1.05,0.5))
plot!(time, mass, label="Modeled", lw=2)
scatter!(TG_time, TG_mass, label="Experimental", color=:pink, marker=:cross)
plot!(twinx(), time, T.(time) .-273, label="T", ylabel="T (ºC)", legend=:topright, color=:orange, ls=:dash)
hline!([0.375-1], ls=:dash, color=:gray, label="")
#display("image/png", plt2)

plt3 = 	plot(time, derivf(time, mass), label="Modeled", xlabel="Time (min)", ylabel="dM/dt (mg/min)", lw=2, legend=:bottomright)
plot!(TG_time, derivf(TG_time, TG_mass), legend=:bottomright, label="Experimental")
plot!(twinx(), TG_time, T.(TG_time), label="T", color=:orange, legend = :topleft)
#display("image/png", plt3)

plot(plt, plt2, plt3, layout=(3,1), size=(700, 800))

## Evaporation rate appears to be in the order of 5k mol/min pretty high to be a limiting factor
#VPCa(t) = IfElse.ifelse(T(t) < 1712., 10^(2.7843 - 3121.368 /(T(t) - 594.591)), 1. ) #in bar #exp(-ΔHCa/R/T(t)+6) 
#ks(t) = IfElse.ifelse(T(t) > 1254., 
#	area * 60. * VPCa(t)*0.98*1E5 / sqrt(2*pi*R*T(t) / (MW[2]*1E-3)) / (MW[2]*1E-3) *1000.,
#	0.) #from bar to Pa, eventually to mmol/ min

#Dif(t) = Ad*4E-2#) * exp(-(Ed*1E3)/R/(T(t)))