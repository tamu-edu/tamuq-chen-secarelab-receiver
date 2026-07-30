using Test
using Statistics

include(joinpath(@__DIR__, "..", "2D_v16.jl"))
using .Receiver2D_v16

function with_contact_2D_v16(p, contact)
    skin = SkinParameters2D(
        thickness=p.skin.thickness,
        channel_conductance_scale=
            p.skin.channel_conductance_scale,
        felt_contact_scale=contact,
    )
    return ModelParameters2D(base=p.base, skin=skin)
end

@testset "2D v16 inherited conservation and equilibrium" begin
    p = default_parameters2D()
    inventory = network_inventory2D(p)
    @test inventory.receiver_mass_kg ≈
          inventory.inherited_receiver_mass_kg rtol=1e-14
    @test inventory.receiver_capacity_J_K ≈
          inventory.inherited_receiver_capacity_J_K rtol=1e-14
    op = OperatingCondition2D(
        irradiance=0.0, flow_lpm=0.0,
        inlet_temperature=298.15,
        ambient_temperature=298.15,
    )
    result = simulate2D(
        p, op, [0.0, 20.0]; initial_temperature=298.15,
    )
    @test maximum(abs, result.skin_temperature .- 298.15) < 3e-7
    @test maximum(abs, result.channel_temperature .- 298.15) <
          3e-7
end

@testset "2D v16 independent skin-felt path" begin
    base = default_parameters2D()
    low = with_contact_2D_v16(base, 0.10)
    high = with_contact_2D_v16(base, 10.0)
    grid = build_network_grid2D(base)
    layout = Receiver2D_v16.V15._state_layout(grid)
    state = fill(300.0, layout.total)
    state[layout.base.channel] .= 700.0
    state[layout.skin] .= 700.0
    op = OperatingCondition2D(
        irradiance=0.0, flow_lpm=0.0,
        inlet_temperature=300.0,
        ambient_temperature=300.0,
    )
    low_result = simulate2D(
        low, op, [0.0, 5.0]; initial_temperature=state,
    )
    high_result = simulate2D(
        high, op, [0.0, 5.0]; initial_temperature=state,
    )
    @test mean(high_result.skin_temperature[:, end]) <
          mean(low_result.skin_temperature[:, end])
    @test low.network.edge_felt_contact_scale ==
          high.network.edge_felt_contact_scale
    @test low.skin.felt_contact_scale !=
          high.skin.felt_contact_scale
end

@testset "2D v16 energy closure" begin
    p = with_contact_2D_v16(default_parameters2D(), 0.30)
    grid = build_network_grid2D(p)
    layout = Receiver2D_v16.V15._state_layout(grid)
    state = fill(700.0, layout.total)
    state[layout.skin] .= range(
        710.0, 780.0; length=length(layout.skin),
    )
    state[layout.base.adaptor] = 650.0
    state[layout.base.housing] = 500.0
    op = OperatingCondition2D(
        irradiance=250000.0, flow_lpm=9.0,
        inlet_temperature=300.0,
        ambient_temperature=300.0,
    )
    ledger = energy_rate_ledger2D(state, p, op, 0.0)
    @test isfinite(ledger.residual)
    @test abs(ledger.residual) < 1e-7
end
