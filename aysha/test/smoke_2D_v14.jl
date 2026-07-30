using Test

include(joinpath(@__DIR__, "..", "2D_v14.jl"))
using .Receiver2D_v14

@testset "2D v14 square-channel orbit inventory" begin
    p = default_parameters2D()
    grid = build_network_grid2D(p)
    inventory = network_inventory2D(p)
    @test grid.group_count == 15
    @test inventory.multiplicity_sum == 100
    @test inventory.boundary_face_sum == 40
    @test inventory.receiver_mass_kg ≈ 0.0400588 rtol=1e-8
    @test inventory.flow_area_m2 ≈ 100 * (1.5e-3)^2
    @test grid.group_keys[grid.core_group] == (1, 1)
    @test grid.group_keys[grid.side_group] == (9, 1)
    @test grid.group_keys[grid.corner_group] == (9, 9)
    @test grid.groove_overlap[grid.core_group] == 0.0
    @test grid.groove_overlap[grid.corner_group] == 1.0
end

@testset "2D v14 thermal-hydraulic feedback" begin
    diagnostic = hydraulic_feedback_diagnostic2D()
    @test diagnostic.heated_to_baseline_side_flow_ratio < 1.0
    @test diagnostic.heated_core_to_side_flow_ratio >
          diagnostic.baseline_core_to_side_flow_ratio
    @test diagnostic.baseline_equal_pressure_error < 1e-4
    @test diagnostic.heated_equal_pressure_error < 1e-4
    @test diagnostic.baseline_total_mass_flow ≈
          diagnostic.heated_total_mass_flow rtol=1e-12
    @test diagnostic.core_groove_overlap == 0.0
    @test diagnostic.corner_groove_overlap == 1.0
end

@testset "2D v14 uniform equilibrium and observations" begin
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
    @test maximum(abs, result.channel_temperature .- 298.15) <
          3e-7
    @test maximum(abs, result.outer_temperature .- 298.15) <
          3e-7
    @test maximum(abs, result.rear_tube_temperature .- 298.15) <
          3e-7
    @test maximum(abs, result.adaptor_temperature .- 298.15) <
          3e-7
    @test maximum(
        abs, result.housing_extension_temperature .- 298.15,
    ) < 3e-7
    predictions = sensor_predictions2D(result)
    for sensor in (:T8, :T12, :T11, :T9, :T10, :T3, :T2)
        @test maximum(
            abs, getproperty(predictions, sensor) .- 298.15,
        ) < 3e-7
    end
end

@testset "2D v14 heated response and mass conservation" begin
    p = default_parameters2D()
    op = OperatingCondition2D(
        irradiance=300000.0,
        flow_lpm=8.0,
        inlet_temperature=300.0,
        ambient_temperature=300.0,
    )
    result = simulate2D(
        p, op, collect(0.0:60.0:600.0);
        initial_temperature=300.0,
    )
    @test all(isfinite, result.channel_temperature)
    @test maximum(result.channel_temperature[:, :, end]) > 300.0
    @test maximum(result.equal_pressure_relative_error) < 1e-4
    @test maximum(
        abs,
        vec(sum(result.group_mass_flow_kg_s; dims=1)) .-
        result.mass_flow_kg_s,
    ) < 1e-12
end

@testset "2D v14 whole-assembly energy closure" begin
    p = default_parameters2D()
    grid = build_network_grid2D(p)
    layout = Receiver2D_v14._state_layout(grid)
    state = fill(700.0, layout.total)
    state[layout.adaptor] = 650.0
    state[layout.housing] = 500.0
    op = OperatingCondition2D(
        irradiance=250000.0,
        flow_lpm=9.0,
        inlet_temperature=300.0,
        ambient_temperature=300.0,
    )
    ledger = energy_rate_ledger2D(state, p, op, 0.0)
    @test isfinite(ledger.residual)
    @test abs(ledger.residual) < 1e-7
end
