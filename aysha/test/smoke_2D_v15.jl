using Test

include(joinpath(@__DIR__, "..", "2D_v15.jl"))
using .Receiver2D_v15

@testset "2D v15 conserved SiC partition" begin
    p = default_parameters2D()
    inventory = network_inventory2D(p)
    @test inventory.group_count == 15
    @test inventory.skin_area_m2 > 0.0
    @test inventory.residual_channel_area_m2 > 0.0
    @test inventory.receiver_solid_area_m2 ≈
          inventory.residual_channel_area_m2 +
          inventory.skin_area_m2 rtol=1e-14
    @test inventory.receiver_mass_kg ≈
          inventory.inherited_receiver_mass_kg rtol=1e-14
    @test inventory.receiver_capacity_J_K ≈
          inventory.inherited_receiver_capacity_J_K rtol=1e-14
end

@testset "2D v15 hydraulics inherited exactly" begin
    p = default_parameters2D()
    v15 = hydraulic_feedback_diagnostic2D(
        p; flow_lpm=9.0, base_temperature=700.0,
    )
    v14 = Receiver2D_v15.V14.hydraulic_feedback_diagnostic2D(
        p.base; flow_lpm=9.0, base_temperature=700.0,
    )
    @test v15.baseline_side_flow_kg_s ==
          v14.baseline_side_flow_kg_s
    @test v15.baseline_core_to_side_flow_ratio ==
          v14.baseline_core_to_side_flow_ratio
    @test v15.heated_core_to_side_flow_ratio ==
          v14.heated_core_to_side_flow_ratio
end

@testset "2D v15 uniform equilibrium and observations" begin
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
    @test maximum(abs, result.skin_temperature .- 298.15) < 3e-7
    @test maximum(abs, result.outer_temperature .- 298.15) < 3e-7
    @test maximum(abs, result.rear_tube_temperature .- 298.15) <
          3e-7
    @test maximum(abs, result.adaptor_temperature .- 298.15) <
          3e-7
    predictions = sensor_predictions2D(result)
    for sensor in (:T8, :T12, :T11, :T9, :T10, :T3, :T2)
        @test maximum(
            abs, getproperty(predictions, sensor) .- 298.15,
        ) < 3e-7
    end
end

@testset "2D v15 heated response and mass conservation" begin
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
    @test all(isfinite, result.skin_temperature)
    @test maximum(result.skin_temperature[:, end]) > 300.0
    @test maximum(result.equal_pressure_relative_error) < 1e-4
    @test maximum(
        abs,
        vec(sum(result.group_mass_flow_kg_s; dims=1)) .-
        result.mass_flow_kg_s,
    ) < 1e-12
end

@testset "2D v15 whole-assembly energy closure" begin
    p = default_parameters2D()
    grid = build_network_grid2D(p)
    layout = Receiver2D_v15._state_layout(grid)
    state = fill(700.0, layout.total)
    state[layout.skin] .= range(720.0, 780.0; length=length(layout.skin))
    state[layout.base.adaptor] = 650.0
    state[layout.base.housing] = 500.0
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
