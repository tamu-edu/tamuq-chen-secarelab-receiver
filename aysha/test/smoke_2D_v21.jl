using Test

include(joinpath(@__DIR__, "..", "run_2D_v21.jl"))

@testset "2D v21 hot-front model" begin
    p = parameters_2D_v21(
        effective_front_area_cm2=10.0,
        power_scales=(1.19, 1.54, 0.94),
    )
    inventory = V21.network_inventory2D(p)
    @test inventory.hot_front_effective_area_m2 ≈ 10e-4
    @test inventory.hot_front_extra_area_m2 > 0.0
    @test inventory.receiver_web_radiating_area_m2 > 0.0

    times = [0.0, 10.0]
    op = V21.OperatingCondition2D(
        irradiance=V21.linear_history(times, [456e3, 456e3]),
        flow_lpm=V21.linear_history(times, [10.0, 10.0]),
        inlet_temperature=V21.linear_history(times, [300.0, 300.0]),
        ambient_temperature=V21.linear_history(times, [300.0, 300.0]),
    )
    result = V21.simulate2D(
        p, op, times; reltol=1e-5, abstol=1e-6,
    )
    @test V21.SciMLBase.successful_retcode(
        result.ode_solution.retcode,
    )
    diagnostic = V21.hot_front_radiative_diagnostic2D(
        result.ode_solution.u[end], p, op, times[end],
    )
    @test diagnostic.extra_radiative_loss_W >= 0.0
    ledger = V21.energy_rate_ledger2D(
        result.ode_solution.u[end], p, op, times[end],
    )
    @test abs(ledger.residual) < 1e-7

    # Setting the effective area equal to the represented web makes the v21
    # correction exactly zero, which is the no-double-counting invariant.
    pzero = parameters_2D_v21(
        effective_front_area_cm2=
            1e4 * inventory.receiver_web_radiating_area_m2,
    )
    zero = V21.hot_front_radiative_diagnostic2D(
        result.ode_solution.u[end], pzero, op, times[end],
    )
    @test zero.extra_area_m2 ≈ 0.0 atol=1e-18
    @test zero.extra_radiative_loss_W ≈ 0.0 atol=1e-12
end
