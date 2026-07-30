include("../2D_v7.jl")
using .Receiver2D_v7

try
    p = default_parameters2D()
    grid = Receiver2D_v7.build_grid2D(p)
    t0_dict = Dict(:T8 => 760.0, :T12 => 750.0, :T11 => 620.0, :T9 => 730.0, :T10 => 600.0, :T3 => 550.0, :T2 => 320.0)
    u0_exp = build_initial_state_2D(grid, p, t0_dict, 295.0)
    times = [0.0, 50.0, 100.0]
    op = OperatingCondition2D(irradiance=304000.0, flow_lpm=10.0, inlet_temperature=295.0, ambient_temperature=295.0)
    res = simulate2D(p, op, times; initial_temperature=u0_exp)
    println("SUCCESS! Result length = ", length(res.time))
catch err
    println("EXACT ERROR: ", typeof(err))
    println(err)
end
