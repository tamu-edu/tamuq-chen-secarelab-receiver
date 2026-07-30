using Test

include(joinpath(@__DIR__, "..", "2D_v8.jl"))
using .Receiver2D_v8

function with_beta_opt(p, beta)
    opt = OpticalParameters2D(
        absorbed_fraction=p.optics.absorbed_fraction,
        extinction_coefficient=beta,
        beam_radius_sigma=p.optics.beam_radius_sigma,
        spillage_fraction=p.optics.spillage_fraction,
        front_deposition_fraction=p.optics.front_deposition_fraction,
        scale_456=p.optics.scale_456,
        scale_304=p.optics.scale_304,
        scale_256=p.optics.scale_256,
    )
    return ModelParameters2D(p.geometry, p.solid, p.heat_transfer, p.losses, opt)
end

@testset "2D_v8 physical invariants and diagnostics" begin
    p = default_parameters2D()
    grid = Receiver2D_v8.build_grid2D(p)
    geom = geometry_invariants2D(p)

    @test geom.receiver_area ≈ 3.61e-4 rtol=1e-12
    @test geom.receiver_mass ≈ 0.0400588 rtol=1e-10
    @test geom.porosity + geom.solid_fraction ≈ 1.0 atol=1e-14

    low_beta = solar_power_budget2D(with_beta_opt(p, 20.0), 456000.0)
    high_beta = solar_power_budget2D(with_beta_opt(p, 300.0), 456000.0)
    @test low_beta.transmitted > high_beta.transmitted
    @test low_beta.deposited < high_beta.deposited
    @test abs(low_beta.closure) < 1e-12
    @test abs(high_beta.closure) < 1e-12

    sensors = Dict(
        :T8 => 760.0, :T12 => 750.0, :T11 => 620.0,
        :T9 => 730.0, :T10 => 600.0, :T3 => 550.0, :T2 => 320.0,
    )
    u0 = build_initial_state_2D(grid, p, sensors, 295.0)
    op = OperatingCondition2D(
        irradiance=304000.0,
        flow_lpm=10.0,
        inlet_temperature=295.0,
        ambient_temperature=295.0,
    )
    ledger = energy_rate_ledger2D(u0, p, op, 0.0)
    scale = max(abs(ledger.expected_denergy_dt), 1.0)
    @test abs(ledger.residual) / scale < 1e-10

    result = simulate2D(p, op, [0.0, 1.0e-6]; initial_temperature=u0)
    predictions = sensor_predictions2D(result)
    @test predictions.T8[1] ≈ sensors[:T8] atol=1e-8
    @test predictions.T12[1] ≈ sensors[:T12] atol=1e-8
    @test predictions.T11[1] ≈ sensors[:T11] atol=1e-8
    @test predictions.T9[1] ≈ sensors[:T9] atol=1e-8
    @test predictions.T10[1] ≈ sensors[:T10] atol=1e-8

    z1 = result.z_rear_gas[1]
    z2 = result.z_rear_gas[2]
    w = (p.geometry.t3_distance_from_receiver - z1) / (z2 - z1)
    expected_t3 = (1.0 - w) * result.rear_gas_temperature[1, 1] +
                  w * result.rear_gas_temperature[2, 1]
    @test predictions.T3[1] ≈ expected_t3 atol=1e-10
end
