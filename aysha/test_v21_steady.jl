include("2D_v21.jl")
using .Receiver2D_v21

println("Setting up Phase 2 steady-state diagnostic...")
p = ModelParameters2D_v21()
# Spillage: 10% of 1500W = 150W deposited directly onto outer orbit
# Core preference: strong bias (0.5 means outer orbit gets half flow per channel compared to core)
p_mod = ModelParameters2D_v21(
    base=p.base,
    phase2=Phase2Parameters2D(
        spillage_power_W=150.0,
        core_preference=0.5,
        spillage_axial_spread_m=10.0e-3
    )
)

op = Receiver2D_v21.V20.OperatingCondition2D(
    irradiance=304e3, # ~ 1.5 kW total power
    flow_lpm=40.0,
    inlet_temperature=295.0,
    ambient_temperature=295.0
)

# Simulate 300s to reach roughly steady state
times = collect(0.0:10.0:300.0)
println("Simulating...")
res = simulate2D_v21(p_mod, op, times; initial_temperature=295.0)

println("\n--- Steady State Profile at 300s ---")
z_depths = [5, 10, 20, 35, 45] # Various depths along the 45 nodes
z_centers = res.z_solid
nr_rec = length(res.group_keys)

for j in z_depths
    z_m = z_centers[j]
    T_core = res.channel_temperature[1, j, end]
    T_perim = res.channel_temperature[nr_rec, j, end]
    println(round(z_m*1000, digits=1), " mm depth | Core: ", round(T_core, digits=1), " K | Perimeter: ", round(T_perim, digits=1), " K | Diff: ", round(T_perim - T_core, digits=1))
end

T_outlet_core = res.gas_temperature[1, end, end]
T_outlet_perim = res.gas_temperature[nr_rec, end, end]
println("\nGas Outlet Core: ", round(T_outlet_core, digits=1), " K")
println("Gas Outlet Perim: ", round(T_outlet_perim, digits=1), " K")
