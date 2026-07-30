# ============================================================================
# smoke_2D_v6.jl - Unit smoke test suite for 2D_v6.jl (SciML ODE Solver)
# ============================================================================

using Test

include(joinpath(@__DIR__, "..", "2D_v6.jl"))
using .Receiver2D_v6

@testset "2D_v6 SciML ODE Solver Smoke Test" begin
    p = default_parameters2D()
    grid = Receiver2D_v6.build_grid2D(p)

    @test grid.nz == 30
    @test length(grid.dz_cells) == 30

    t0_dict = Dict(:T8 => 760.0, :T12 => 750.0, :T11 => 620.0, :T9 => 730.0, :T10 => 600.0, :T3 => 550.0, :T2 => 320.0)
    u0_exp = build_initial_state_2D(grid, p, t0_dict, 295.0)
    @test length(u0_exp) == grid.nr_total * grid.nz + grid.nz_rear + 1
    @test all(isfinite, u0_exp)

    times = [0.0, 50.0, 100.0]
    op = OperatingCondition2D(irradiance=304000.0, flow_lpm=10.0, inlet_temperature=295.0, ambient_temperature=295.0)

    result = simulate2D(p, op, times; initial_temperature=u0_exp)

    @test length(result.time) == 3
    @test size(result.solid_temperature) == (17, 30, 3)
    @test size(result.gas_temperature) == (10, 31, 3)
    @test all(isfinite, result.solid_temperature)

    preds = sensor_predictions2D(result)
    @test length(preds.T8) == 3
    @test length(preds.T12) == 3
    @test length(preds.T11) == 3
end
