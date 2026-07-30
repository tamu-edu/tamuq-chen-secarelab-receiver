using Test

include(joinpath(@__DIR__, "..", "2D_v12.jl"))
using .Receiver2D_v12

@testset "2D v12 hardware inventory" begin
    p = default_parameters2D()
    inventory = hardware_inventory2D(p)
    @test inventory.receiver_mass_kg ≈ 0.0400588 rtol=1e-4
    @test inventory.adaptor_mass_kg ≈ 0.99604 rtol=1e-4
    @test inventory.felt_mass_kg ≈ 0.19258 rtol=1e-4
    @test inventory.aluminum_sleeve_mass_kg ≈ 3.3254 rtol=1e-4
    @test inventory.aluminum_backplate_mass_kg ≈ 0.85238 rtol=1e-4
    @test felt_conductivity_v12(773.15) ≈ 0.08 rtol=1e-8
    @test felt_conductivity_v12(1073.15) ≈ 0.11 rtol=1e-8
    @test felt_conductivity_v12(1373.15) ≈ 0.17 rtol=1e-8
end

@testset "2D v12 uniform equilibrium and observations" begin
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
    @test maximum(abs, result.rear_tube_temperature .- 298.15) < 3e-7
    @test maximum(abs, result.adaptor_temperature .- 298.15) < 3e-7
    @test maximum(abs, result.housing_extension_temperature .- 298.15) < 3e-7
    predictions = sensor_predictions2D(result)
    for sensor in (:T8, :T12, :T11, :T9, :T10, :T3, :T2)
        @test maximum(abs, getproperty(predictions, sensor) .- 298.15) < 3e-7
    end
end

@testset "2D v12 contact and observation sensitivity" begin
    base = default_parameters2D()
    op = OperatingCondition2D(
        irradiance=300000.0,
        flow_lpm=8.0,
        inlet_temperature=300.0,
        ambient_temperature=300.0,
    )
    times = collect(0.0:60.0:600.0)
    result = simulate2D(base, op, times; initial_temperature=300.0)
    predictions = sensor_predictions2D(result)
    @test all(isfinite, result.solid_temperature)
    @test result.adaptor_temperature[end] > 300.0
    @test result.housing_extension_temperature[end] >= 300.0
    @test predictions.T9[end] <=
        max(predictions.T9_wall[end], predictions.T9_gas[end]) + 5.0
    @test predictions.T10[end] <=
        max(predictions.T10_wall[end], predictions.T10_gas[end]) + 5.0
    @test predictions.T9[end] != predictions.T9_wall[end]
    @test predictions.T10[end] != predictions.T10_wall[end]
    @test maximum(abs, result.mass_flow_kg_s .-
        result.mass_flow_kg_s[1]) < 1e-12
    @test maximum(result.equal_pressure_relative_error) < 1e-4
end

@testset "2D v12 augmented energy closure" begin
    p = default_parameters2D()
    op = OperatingCondition2D(
        irradiance=250000.0,
        flow_lpm=9.0,
        inlet_temperature=300.0,
        ambient_temperature=300.0,
    )
    grid = Receiver2D_v12.Receiver2D_v11.build_grid2D(p.base)
    nbase = grid.nr_total * grid.nz + grid.nz_rear
    state = fill(700.0, nbase + 2)
    state[nbase + 1] = 650.0
    state[nbase + 2] = 500.0
    ledger = energy_rate_ledger2D(state, p, op, 0.0)
    @test isfinite(ledger.residual)
    @test abs(ledger.residual) < 1e-7
end
