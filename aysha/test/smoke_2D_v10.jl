using Test
using Statistics

include(joinpath(@__DIR__, "..", "2D_v10.jl"))
using .Receiver2D_v10

function with_front_coefficient(p::ModelParameters2D, coefficient::Real)
    h0 = p.heat_transfer
    heat = HeatTransferParameters2D(
        coefficient=h0.coefficient,
        reynolds_exponent=h0.reynolds_exponent,
        prandtl_exponent=h0.prandtl_exponent,
        minimum_nusselt=h0.minimum_nusselt,
        c_radial_flow=h0.c_radial_flow,
        front_coefficient=Float64(coefficient),
        front_reynolds_exponent=h0.front_reynolds_exponent,
    )
    return ModelParameters2D(
        p.geometry, p.solid, heat, p.losses, p.optics, p.hydraulics,
    )
end

@testset "2D_v10 conservative flow-sensitive front exchange" begin
    p0 = default_parameters2D()
    p1 = with_front_coefficient(p0, 0.5)
    grid = Receiver2D_v10.build_grid2D(p0)
    sensors = Dict(
        :T8 => 760.0, :T12 => 750.0, :T11 => 620.0,
        :T9 => 730.0, :T10 => 600.0, :T3 => 550.0, :T2 => 320.0,
    )
    u0 = build_initial_state_2D(grid, p0, sensors, 295.0)
    op = OperatingCondition2D(
        irradiance=0.0,
        flow_lpm=10.0,
        inlet_temperature=295.0,
        ambient_temperature=295.0,
    )

    no_front = simulate2D(p0, op, [0.0, 1.0e-6]; initial_temperature=u0)
    front = simulate2D(p1, op, [0.0, 1.0e-6]; initial_temperature=u0)

    @test all(iszero, no_front.front_heat_transfer_coefficient)
    @test all(iszero, no_front.front_gas_heat_transfer_W)
    @test all(front.front_heat_transfer_coefficient .> 0.0)
    @test all(front.front_gas_heat_transfer_W[:, 1] .> 0.0)
    @test all(front.gas_temperature[:, 1, 1] .> 295.0)
    @test all(front.gas_temperature[:, 1, 1] .<
              front.solid_temperature[1:grid.nr_rec, 1, 1])

    flow_lpm = op.flow_lpm(0.0)
    mdot_total = standard_mass_flow2D(flow_lpm, p1.hydraulics)
    radius = grid.r_faces[grid.nr_rec + 1]
    psi = [
        1.0 - p1.heat_transfer.c_radial_flow *
        (grid.r_centers[i] / radius)^2 for i in 1:grid.nr_rec
    ]
    psi_sum = sum(psi .* grid.area_flow[1:grid.nr_rec])
    cp_in = air_heat_capacity(op.inlet_temperature(0.0))
    reconstructed_front_q = [
        mdot_total * psi[i] * grid.area_flow[i] / psi_sum * cp_in *
        (front.gas_temperature[i, 1, 1] - op.inlet_temperature(0.0))
        for i in 1:grid.nr_rec
    ]
    @test reconstructed_front_q ≈ front.front_gas_heat_transfer_W[:, 1] rtol=1e-12
    @test front.mass_flow_kg_s == no_front.mass_flow_kg_s

    ledger = energy_rate_ledger2D(u0, p1, op, 0.0)
    @test ledger.gas_front ≈ sum(front.front_gas_heat_transfer_W[:, 1]) rtol=1e-12
    @test abs(ledger.residual) < 1e-8

    equilibrium_op = OperatingCondition2D(
        irradiance=0.0,
        flow_lpm=10.0,
        inlet_temperature=300.0,
        ambient_temperature=300.0,
    )
    equilibrium = simulate2D(
        p1, equilibrium_op, [0.0, 1.0e-6]; initial_temperature=300.0,
    )
    @test maximum(abs, equilibrium.front_gas_heat_transfer_W) < 1e-12
    @test maximum(abs, equilibrium.solid_temperature .- 300.0) < 1e-8
    equilibrium_ledger = energy_rate_ledger2D(
        fill(300.0, grid.nr_total * grid.nz + grid.nz_rear),
        p1,
        equilibrium_op,
        0.0,
    )
    @test abs(equilibrium_ledger.residual) < 1e-10

    zero_flow_op = OperatingCondition2D(
        irradiance=0.0,
        flow_lpm=0.0,
        inlet_temperature=295.0,
        ambient_temperature=295.0,
    )
    zero_flow = simulate2D(
        p1, zero_flow_op, [0.0, 1.0e-6]; initial_temperature=u0,
    )
    @test all(iszero, zero_flow.front_heat_transfer_coefficient)
    @test all(iszero, zero_flow.front_gas_heat_transfer_W)
    @test all(iszero, zero_flow.gas_velocity)
    @test all(iszero, zero_flow.cell_pressure_drop)

    @test size(front.front_heat_transfer_coefficient) == (grid.nr_rec, 2)
    @test size(front.front_gas_heat_transfer_W) == (grid.nr_rec, 2)
    @test length(pack_parameters2D(p1)) == 12
    @test unpack_parameters2D(pack_parameters2D(p1), p1).heat_transfer.front_coefficient ==
          p1.heat_transfer.front_coefficient
end
