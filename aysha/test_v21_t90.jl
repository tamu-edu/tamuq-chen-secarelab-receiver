include("2D_v21.jl")
using .Receiver2D_v21

println("Running v20 (old SiC capacity, 400 J/kg/K) cooldown...")
op = Receiver2D_v21.V20.OperatingCondition2D(
    irradiance=0.0,
    flow_lpm=20.0,
    inlet_temperature=295.0,
    ambient_temperature=295.0
)
p_v20 = Receiver2D_v21.V20.default_parameters2D()
times = collect(0.0:10.0:2000.0)
res_v20 = Receiver2D_v21.V20.simulate2D(p_v20, op, times; initial_temperature=1000.0)

T_targ = 295.0 + 0.1 * (1000.0 - 295.0)
idx_v20 = findfirst(T -> T <= T_targ, res_v20.channel_temperature[1, 10, :])
t90_v20 = idx_v20 === nothing ? NaN : times[idx_v20]

println("Running v21 (new SiC capacity, ~1100 J/kg/K) cooldown...")
p_v21 = ModelParameters2D_v21()
res_v21 = simulate2D_v21(p_v21, op, times; initial_temperature=1000.0)
idx_v21 = findfirst(T -> T <= T_targ, res_v21.channel_temperature[1, 10, :])
t90_v21 = idx_v21 === nothing ? NaN : times[idx_v21]

println("\n--- T90 Results (Solid Node 1, 10) ---")
println("v20 t90: ", t90_v20, " s")
println("v21 t90: ", t90_v21, " s")
println("Factor increase: ", round(t90_v21 / t90_v20, digits=2))
