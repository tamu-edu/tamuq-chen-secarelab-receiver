using Test

include(joinpath(@__DIR__, "..", "2D_v13.jl"))
using .Receiver2D_v13

@testset "2D v13 receiver-population inventory" begin
    p = default_parameters2D()
    inventory = population_inventory2D(p)
    @test inventory.active_receiver_capacity_J_K +
          inventory.side_receiver_capacity_J_K ≈
          inventory.total_receiver_capacity_J_K rtol=1e-12
    @test inventory.side_solid_fraction == 0.30
end

@testset "2D v13 uniform equilibrium and observations" begin
    p = default_parameters2D()
    op = OperatingCondition2D(
        irradiance=0.0,
        flow_lpm=0.0,
        inlet_temperature=298.15,
        ambient_temperature=298.15,
    )
    result = simulate2D(
        p, op, [0.0, 10.0, 30.0];
        initial_temperature=298.15,
    )
    @test maximum(abs, result.solid_temperature .- 298.15) < 3e-7
    @test maximum(abs, result.side_temperature .- 298.15) < 3e-7
    @test maximum(abs, result.rear_tube_temperature .- 298.15) < 3e-7
    @test maximum(abs, result.adaptor_temperature .- 298.15) < 3e-7
    @test maximum(abs, result.housing_extension_temperature .- 298.15) < 3e-7
    predictions = sensor_predictions2D(result)
    for sensor in (:T8, :T12, :T11, :T9, :T10, :T3, :T2)
        @test maximum(abs, getproperty(predictions, sensor) .- 298.15) < 3e-7
    end
end

@testset "2D v13 heated response and hydraulics" begin
    p = default_parameters2D()
    op = OperatingCondition2D(
        irradiance=300000.0,
        flow_lpm=8.0,
        inlet_temperature=300.0,
        ambient_temperature=300.0,
    )
    times = collect(0.0:60.0:600.0)
    result = simulate2D(p, op, times; initial_temperature=300.0)
    predictions = sensor_predictions2D(result)
    @test all(isfinite, result.solid_temperature)
    @test all(isfinite, result.side_temperature)
    @test result.side_temperature[1, end] > 300.0
    @test predictions.T8[end] > 300.0
    @test maximum(abs, result.mass_flow_kg_s .-
        result.mass_flow_kg_s[1]) < 1e-12
    @test maximum(result.equal_pressure_relative_error) < 1e-4
end

@testset "2D v13 augmented energy closure" begin
    p = default_parameters2D()
    op = OperatingCondition2D(
        irradiance=250000.0,
        flow_lpm=9.0,
        inlet_temperature=300.0,
        ambient_temperature=300.0,
    )
    grid = Receiver2D_v13.V11.build_grid2D(p.base.base)
    nbase = grid.nr_total * grid.nz + grid.nz_rear
    n12 = nbase + 2
    state = fill(700.0, n12 + grid.nz)
    state[nbase + 1] = 650.0
    state[nbase + 2] = 500.0
    state[(n12 + 1):end] .= range(760.0, 680.0; length=grid.nz)
    ledger = energy_rate_ledger2D(state, p, op, 0.0)
    @test isfinite(ledger.residual)
    @test abs(ledger.residual) < 1e-7
end
