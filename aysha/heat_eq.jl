using OrdinaryDiffEq, ModelingToolkit, MethodOfLines, DomainSets
using Plots
using Interpolations

function load_data(measurements)
    #df   = CSV.read(path, DataFrame)
    gexp = unique(measurements.simulation_id)
    out  = Vector{NamedTuple}(undef, length(gexp))
    for (i, g) in enumerate(gexp)
        t = measurements[(measurements.simulation_id .== g) .& (measurements.obs_id .== "_T2"), :time][1]
        Tin = measurements[(measurements.simulation_id .== g) .& (measurements.obs_id .== "_T9"), :temperatures][1]
        Tprobe = measurements[(measurements.simulation_id .== g) .& (measurements.obs_id .== "_T2"), :temperatures][1]
        T8 = measurements[(measurements.simulation_id .== g) .& (measurements.obs_id .== "_T8"), :temperatures][1]


        out[i] = (t     = t,
                   Tin   = Tin,
                   Tp    = Tprobe,
                   T8 = T8)
    end
    return out
end

# Parameters, variables, and derivatives
@parameters t r
@variables T(..)
Dt = Differential(t)
Dr = Differential(r)

ρ = 0.016  #kg/m3 Density
cp = 1360  #J/kg/K Specific heat capacity
k_ins = 0.05  #W/m/K Thermal conductivity
k_air = 0.01
k = k_ins
res = 0.002 #m2 K / W
Tin = 1000  #K Internal temperature
Tout = 300  #K External temperature
Tinit = 300  #K Initial temperature
Rin = 0.02/2  # Inner radius
Rout = 0.150/2  # Outer radius
Rp = Rin + 0.06
hnat = 10.0  # W/m2/K Natural convection heat transfer coefficient
σ = 5.67e-8  # W/m2/K4 Stefan-Boltzmann constant
ε = 0.2  # Emissivity

include("import_exp_0D.jl") #import the experimental data
dataset = load_data(measurements)  

tvec   = dataset[1].t
Tinvec = dataset[1].Tin
Tpvec = dataset[1].Tp
Tin_itp = LinearInterpolation(tvec, Tinvec,extrapolation_bc = Interpolations.Flat())
f_Tin(t) = Tin_itp(t)
@register_symbolic f_Tin(t)
Tp_itp = LinearInterpolation(tvec, Tpvec,extrapolation_bc = Interpolations.Flat())
f_Tp(t) = Tp_itp(t)
@register_symbolic f_Tp(t)
plot(tvec,f_Tin.(tvec), label="Tin", color=:blue)
plot!(tvec,f_Tp.(tvec), label="Tp", color=:red)

# 1D PDE and boundary conditions
eq  = ρ * cp * Dt(T(t, r)) ~ k * 1/r * Dr(r * Dr(T(t, r)))
bcs = [T(0, r) ~ f_Tin(0),
        #T(t, Rin) ~ f_Tin(t),
        - k* Dr(T(t, Rin)) ~ (f_Tin(t) - T(t, Rin)) / res,
        #- k* Dr(T(t, Rout)) ~ - hnat * (Tout - T(t, Rout)) - σ * ε * (Tout^4 - T(t, Rout)^4)]
        T(t, Rout) ~ Tout]

# Space and time domains
domains = [t ∈ Interval(0.0, 3600.0),  # Time in seconds
           r ∈ Interval(Rin, Rout)]

# PDE system
@named pdesys = PDESystem(eq, bcs, domains, [t, r], [T(t, r)])

# Method of lines discretization
dr = (Rout - Rin)/100
order = 2
discretization = MOLFiniteDifference([r => dr], t, approx_oder=order)

# Convert the PDE problem into an ODE problem
prob = discretize(pdesys,discretization)

# Solve ODE problem
sol = solve(prob, saveat=10.)

# Plot results and compare with exact solution
discrete_r = sol[r]
discrete_t = sol[t]
solT = sol[T(t, r)]

plt1 = plot(xlabel="Radius (m)", ylabel="Temperature (K)")
time_samples = [1, 60, 120, 180, 360]
#plot!(tvec, Tinvec, label="Tin", color=:blue)
for ts in time_samples  # Time samples for plotting
    tt = discrete_t[ts]
    plot!(discrete_r, solT[ts, :], label="t=$tt")
end
display(plt1)

plt2 = plot(xlabel="time (s)", ylabel="Temperature (K)")
plot!(discrete_t, f_Tin.(discrete_t), label="Tin", color=:blue)
plot!(discrete_t, f_Tp.(discrete_t), label="Tp", color=:red, linestyle=:dash)
idx_p = round(Int, (Rp-Rin)/dr) + 1  # Index for the probe position
plot!(discrete_t, solT[:, idx_p], label="Simulated Tp", color=:green)
display(plt2)


#     plot!(discrete_x, solu[i, :], label="Numerical, t=$(discrete_t[i])")
#     scatter!(discrete_x, u_exact(discrete_x, discrete_t[i]), label="Exact, t=$(discrete_t[i])")
# end
# plt

# work with the cross section of temperatures
begin
    IoA = 456000.0 * 1.0
    IoB = 304000.0 * 1.0
    IoC = 256000.0 * 1.0 # arbitrary fix KK

    ArIo=[IoA, IoA, IoA, IoA, IoA,
        IoB, IoB, IoB, IoB, IoB,
        IoC, IoC, IoC, IoC, IoC] * 1.0
    Arqplm=[15.27, 12.50, 10.50, 9.10, 7.12,
         18.34, 13.16, 9.03, 6.95, 4.53,
         13.85, 10.02, 8.04, 6.62, 4.53]

    function get_cross(measurements)
        #df   = CSV.read(path, DataFrame)
        gexp = unique(measurements.simulation_id)
        out  = Vector{NamedTuple}(undef, 1)
        _t = zeros(length(gexp))
        _T8 = zeros(length(gexp))
        _dT = zeros(length(gexp))

        for (i, g) in enumerate(gexp)
            t = measurements[(measurements.simulation_id .== g) .& (measurements.obs_id .== "_T2"), :time][1]
            T8 = measurements[(measurements.simulation_id .== g) .& (measurements.obs_id .== "_T8"), :temperatures][1]
            T9 = measurements[(measurements.simulation_id .== g) .& (measurements.obs_id .== "_T9"), :temperatures][1] 

            delta_T = abs.(T8 .- T9)
            # f_T8= LinearInterpolation(t, T8,extrapolation_bc = Interpolations.Flat())
            # f_T9= LinearInterpolation(t, T9,extrapolation_bc = Interpolations.Flat())
            # deltaf(t) = abs(f_T8(t) - f_T9(t))
            skip_i = 3
            minv, idx = findmin(delta_T[(1+skip_i):end])
            idx += skip_i

            _t[i]     = t[idx]
            _T8[i]   = T8[idx]
            _dT[i]    = delta_T[idx]
            #sort!(temp, by=3)

        end
        out = (t=_t, T8= _T8, dT=_dT)
        return out
    end

    cross = get_cross(measurements)

    #plot in groups of 5
    plt_cross = plot(xlabel="time (s)", ylabel="Temperature (K)")
    for i in 0:2
        idx = i*5+1
        scatter!(cross.t[idx:(idx+1)], cross.T8[idx:(idx+1)], label="$(i)")
    end
    display(plt_cross)

     plt_cross = plot(xlabel="Flow rate", ylabel="Temperature (K)")
    for i in 0:2
        idx = i*5+1
        scatter!(Arqplm[idx:(idx+1)], cross.T8[idx:(idx+1)], label="$(i)")
    end
    display(plt_cross)

end