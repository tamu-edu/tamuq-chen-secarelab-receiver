include("2D_v21.jl")
using .Receiver2D_v21

println("Building model parameters...")
p = ModelParameters2D_v21()

println("Setting up operating condition...")
op = Receiver2D_v21.V20.OperatingCondition2D(
    irradiance=304e3,
    flow_lpm=20.0,
    inlet_temperature=295.0,
    ambient_temperature=295.0
)

println("Simulating 5 seconds...")
res = simulate2D_v21(p, op, [0.0, 5.0]; initial_temperature=295.0)

println("Success! T3 at 5s = ", res.gas_temperature[1, end, end])
