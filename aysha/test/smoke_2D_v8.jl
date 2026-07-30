# ============================================================================
# smoke_2D_v8.jl - Unit smoke test suite for 2D_v8.jl (SciML ODE Solver)
# ============================================================================

using Test
using Statistics

include(joinpath(@__DIR__, "..", "2D_v8.jl"))
using .Receiver2D_v8

@testset "2D_v8 SciML ODE Solver Smoke Test" begin
    p = default_parameters2D()
    grid = Receiver2D_v8.build_grid2D(p)

    @test grid.nz == 45
    @test length(grid.dz_cells) == 45
    invariants = geometry_invariants2D(p)
    @test invariants.receiver_area ≈ (19.0e-3)^2 rtol=1e-12
    @test invariants.equivalent_radius ≈ sqrt((19.0e-3)^2 / pi) rtol=1e-12
    @test invariants.porosity ≈ 0.6232686980609419 rtol=1e-12
    @test invariants.receiver_mass ≈ 0.0400588 rtol=1e-10
    @test 110.0 < sic_conductivity(300.0) < 125.0
    @test 55.0 < sic_conductivity(1200.0) < 70.0
    @test sic_conductivity(1200.0) < sic_conductivity(300.0)
    @test p.heat_transfer.minimum_nusselt == 0.0

    t0_dict = Dict(:T8 => 760.0, :T12 => 750.0, :T11 => 620.0, :T9 => 730.0, :T10 => 600.0, :T3 => 550.0, :T2 => 320.0)
    u0_exp = build_initial_state_2D(grid, p, t0_dict, 295.0)
    @test length(u0_exp) == grid.nr_total * grid.nz + grid.nz_rear
    @test all(isfinite, u0_exp)

    times = [0.0, 50.0, 100.0]
    op = OperatingCondition2D(irradiance=304000.0, flow_lpm=10.0, inlet_temperature=295.0, ambient_temperature=295.0)
    budget = solar_power_budget2D(p, 304000.0)
    @test budget.delivered ≈ 304000.0 * (19.0e-3)^2 rtol=1e-12
    @test budget.closure ≈ 0.0 atol=1e-12
    @test budget.transmitted > 0.0
    @test budget.deposited < budget.delivered

    ledger = energy_rate_ledger2D(u0_exp, p, op, 0.0)
    @test abs(ledger.residual) < 1e-8

    result = simulate2D(p, op, times; initial_temperature=u0_exp)

    @test length(result.time) == 3
    @test size(result.solid_temperature) == (24, 45, 3)
    @test size(result.gas_temperature) == (14, 46, 3)
    @test all(isfinite, result.solid_temperature)

    preds = sensor_predictions2D(result)
    @test length(preds.T8) == 3
    @test length(preds.T12) == 3
    @test length(preds.T11) == 3
    @test length(preds.T3) == 3
    @test get_t90_2D([0.0, 1.0, 2.0], [300.0, 390.0, 400.0]) == 1.0
    @test get_t90_2D([0.0, 1.0, 2.0], [400.0, 310.0, 300.0]) == 1.0

    loss_adiabatic = LossParameters2D(
        front_emissivity=p.losses.front_emissivity,
        casing_emissivity=p.losses.casing_emissivity,
        front_loss_scale=0.0,
        casing_loss_scale=0.0,
        rear_adaptor_conductance=p.losses.rear_adaptor_conductance,
        flange_conductance_scale=0.0,
    )
    p_adiabatic = ModelParameters2D(p.geometry, p.solid, p.heat_transfer, loss_adiabatic, p.optics)
    op_adiabatic = OperatingCondition2D(
        irradiance=0.0, flow_lpm=0.0,
        inlet_temperature=500.0, ambient_temperature=500.0,
    )
    result_adiabatic = simulate2D(p_adiabatic, op_adiabatic, [0.0, 10.0]; initial_temperature=500.0)
    @test maximum(abs.(result_adiabatic.solid_temperature .- 500.0)) < 1e-8
    @test maximum(abs.(result_adiabatic.rear_tube_temperature .- 500.0)) < 1e-8

    theta = pack_parameters2D(p)
    @test length(theta) == 12
    @test FIT_INDICES_2D == [1, 3, 4, 7, 8, 9, 10, 12]
    @test setdiff(1:12, FIT_INDICES_2D) == [2, 5, 6, 11]
    theta_high_nu = copy(theta)
    theta_high_nu[1] *= 2.0
    p_high_nu = unpack_parameters2D(theta_high_nu)
    op_exchange = OperatingCondition2D(
        irradiance=0.0, flow_lpm=10.0,
        inlet_temperature=295.0, ambient_temperature=600.0,
    )
    low_nu_result = simulate2D(p, op_exchange, [0.0, 1.0e-6]; initial_temperature=600.0)
    high_nu_result = simulate2D(p_high_nu, op_exchange, [0.0, 1.0e-6]; initial_temperature=600.0)
    @test mean(high_nu_result.gas_temperature[:, end, 1]) >
          mean(low_nu_result.gas_temperature[:, end, 1])
    @test pack_parameters2D(unpack_parameters2D(theta)) ≈ theta
end
