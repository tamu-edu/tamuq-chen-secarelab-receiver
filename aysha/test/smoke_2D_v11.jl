using Test
using Statistics

include(joinpath(@__DIR__, "..", "2D_v11.jl"))
using .Receiver2D_v11

function v11_test_parameters(; groove_K=0.0, distribution=:equal_pressure)
    h = HydraulicParameters2D(
        flow_distribution_model=distribution,
        groove_free_diameter=13.0e-3,
        groove_loss_coefficient=groove_K,
        minor_loss_coefficient=0.0,
    )
    return ModelParameters2D(hydraulics=h)
end

function steady_micro_result(p, flow)
    op = OperatingCondition2D(
        irradiance=0.0,
        flow_lpm=flow,
        inlet_temperature=300.0,
        ambient_temperature=300.0,
    )
    return simulate2D(
        p, op, [0.0, 1.0e-6]; initial_temperature=300.0,
    )
end

@testset "2D_v11 local Graetz and measured groove geometry" begin
    p = v11_test_parameters()
    grid = Receiver2D_v11.build_grid2D(p)
    hp = p.heat_transfer

    low_gz = graetz_nusselt2D(1.0, 0.7, 1.0, grid.dh, hp)
    high_gz = graetz_nusselt2D(150.0, 0.7, 1.0e-3, grid.dh, hp)
    @test low_gz.nusselt > hp.fully_developed_nusselt
    @test low_gz.nusselt - hp.fully_developed_nusselt < 1.0e-3
    @test high_gz.nusselt > low_gz.nusselt
    @test high_gz.graetz > low_gz.graetz

    overlap = groove_overlap_fractions2D(grid, p)
    exposed_fraction = sum(
        overlap[i] * grid.area_frt[i] for i in 1:grid.nr_rec
    ) / grid.total_frt
    expected_exposed = 1.0 - pi * (6.5e-3)^2 / (19.0e-3)^2
    @test exposed_fraction ≈ expected_exposed rtol=1e-12
    @test minimum(overlap) == 0.0
    @test maximum(overlap) == 1.0
end

@testset "2D_v11 equal pressure and quadratic groove" begin
    p_open = v11_test_parameters(groove_K=0.0)
    p_groove = v11_test_parameters(groove_K=40.0)
    open_result = steady_micro_result(p_open, 18.0)
    groove_result = steady_micro_result(p_groove, 18.0)
    grid = Receiver2D_v11.build_grid2D(p_groove)

    @test sum(groove_result.ring_mass_flow_kg_s[:, end]) ≈
          groove_result.mass_flow_kg_s[end] rtol=1e-12
    @test groove_result.equal_pressure_relative_error[end] <
          p_groove.hydraulics.equal_pressure_tolerance
    @test maximum(groove_result.ring_pressure_drop_Pa[:, end]) -
          minimum(groove_result.ring_pressure_drop_Pa[:, end]) <
          1.0e-3
    @test all(iszero, open_result.groove_pressure_drop_Pa)
    @test maximum(groove_result.groove_pressure_drop_Pa[:, end]) > 0.0
    @test all(iszero, groove_result.front_gas_heat_transfer_W)

    mass_flux = groove_result.ring_mass_flow_kg_s[:, end] ./
                grid.area_flow[1:grid.nr_rec]
    core = mean(mass_flux[groove_result.groove_overlap_fraction .== 0.0])
    edge = mean(mass_flux[groove_result.groove_overlap_fraction .== 1.0])
    @test edge < core

    low_result = steady_micro_result(p_groove, 5.0)
    low_flux = low_result.ring_mass_flow_kg_s[:, end] ./
               grid.area_flow[1:grid.nr_rec]
    low_core = mean(low_flux[low_result.groove_overlap_fraction .== 0.0])
    low_edge = mean(low_flux[low_result.groove_overlap_fraction .== 1.0])
    @test low_core / low_edge < core / edge

    op = OperatingCondition2D(
        irradiance=0.0,
        flow_lpm=18.0,
        inlet_temperature=300.0,
        ambient_temperature=300.0,
    )
    ledger = energy_rate_ledger2D(
        fill(
            300.0,
            grid.nr_total * grid.nz + grid.nz_rear,
        ),
        p_groove,
        op,
        0.0,
    )
    @test abs(ledger.residual) < 1e-8
end

@testset "2D_v11 result diagnostics and zero flow" begin
    p = v11_test_parameters(groove_K=20.0)
    result = steady_micro_result(p, 10.0)
    grid = Receiver2D_v11.build_grid2D(p)
    @test size(result.gas_prandtl) == (grid.nr_rec, grid.nz, 2)
    @test size(result.gas_graetz) == (grid.nr_rec, grid.nz, 2)
    @test size(result.gas_nusselt) == (grid.nr_rec, grid.nz, 2)
    @test size(result.ring_mass_flow_kg_s) == (grid.nr_rec, 2)
    @test size(result.groove_pressure_drop_Pa) == (grid.nr_rec, 2)
    @test all(result.gas_nusselt .>=
              p.heat_transfer.fully_developed_nusselt)
    @test mean(result.gas_nusselt[:, 1, end]) >
          mean(result.gas_nusselt[:, end, end])

    zero_op = OperatingCondition2D(
        irradiance=0.0,
        flow_lpm=0.0,
        inlet_temperature=300.0,
        ambient_temperature=300.0,
    )
    zero = simulate2D(
        p, zero_op, [0.0, 1.0e-6]; initial_temperature=300.0,
    )
    @test all(iszero, zero.ring_mass_flow_kg_s)
    @test all(iszero, zero.gas_velocity)
    @test all(iszero, zero.gas_nusselt)
    @test all(iszero, zero.cell_pressure_drop)
end
