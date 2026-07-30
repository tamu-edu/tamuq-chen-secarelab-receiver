using Test
using Statistics

include(joinpath(@__DIR__, "..", "run_2D_v17.jl"))

@testset "2D v17 verified casing/flange topology" begin
    low = parameters_2D_v17(
        casing_flange_conductance_W_K=0.0,
    )
    high = parameters_2D_v17(
        casing_flange_conductance_W_K=5.0,
    )
    @test low.casing_flange.conductance_W_K == 0.0
    @test high.casing_flange.conductance_W_K == 5.0
    @test V17.V11.SIDE_TC_FRONT_Z_2D == 11.0e-3

    grid = V17.build_network_grid2D(high)
    layout = V17.V15._state_layout(grid)
    state = fill(500.0, layout.total)
    state[layout.base.housing] = 700.0
    op = V17.OperatingCondition2D(
        irradiance=0.0, flow_lpm=0.0,
        inlet_temperature=500.0,
        ambient_temperature=500.0,
    )
    low_result = V17.simulate2D(
        low, op, [0.0, 20.0]; initial_temperature=state,
    )
    high_result = V17.simulate2D(
        high, op, [0.0, 20.0]; initial_temperature=state,
    )
    @test high_result.housing_extension_temperature[end] <
          low_result.housing_extension_temperature[end]

    ledger = V17.energy_rate_ledger2D(state, high, op, 0.0)
    @test ledger.casing_flange_loss ≈
          5.0 * (700.0 - V17.V11.WATER_FLANGE_TEMP)
    @test ledger.flange_loss ≈
          ledger.tube_flange_loss + ledger.casing_flange_loss
    @test abs(ledger.residual) < 1e-7
end

@testset "2D v17 equilibrium" begin
    p = parameters_2D_v17(
        casing_flange_conductance_W_K=5.0,
    )
    op = V17.OperatingCondition2D(
        irradiance=0.0, flow_lpm=0.0,
        inlet_temperature=V17.V11.WATER_FLANGE_TEMP,
        ambient_temperature=V17.V11.WATER_FLANGE_TEMP,
    )
    result = V17.simulate2D(
        p, op, [0.0, 20.0];
        initial_temperature=V17.V11.WATER_FLANGE_TEMP,
    )
    @test maximum(
        abs, result.skin_temperature .-
             V17.V11.WATER_FLANGE_TEMP,
    ) < 3e-7
end
