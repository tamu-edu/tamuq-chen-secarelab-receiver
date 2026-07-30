using Test

include(joinpath(@__DIR__, "..", "run_2D_v20.jl"))

@testset "2D v20 enthalpy thermodynamics" begin
    temperatures = (100.0, 199.9, 200.0, 300.0, 800.0, 1500.0, 3500.0)
    enthalpies = V20.air_specific_enthalpy2D.(temperatures)
    @test all(diff(collect(enthalpies)) .> 0.0)
    for T in temperatures
        recovered = V20.air_temperature_from_enthalpy2D(
            V20.air_specific_enthalpy2D(T);
            lower=100.0, upper=3500.0,
        )
        @test recovered ≈ T atol=1e-8
    end
    for T in (250.0, 300.0, 800.0, 1500.0, 3000.0)
        step = 1e-3
        derivative = (
            V20.air_specific_enthalpy2D(T + step) -
            V20.air_specific_enthalpy2D(T - step)
        ) / (2step)
        @test derivative ≈ V20.V11.air_heat_capacity(T) rtol=1e-8
    end

    for (Tin, Twall) in ((300.0, 1000.0), (1000.0, 300.0))
        cell = V20._enthalpy_cell_exchange2D(
            Tin, Twall, 0.04, 1e-4,
        )
        @test min(Tin, Twall) <= cell.temperature <= max(Tin, Twall)
        @test cell.heat_W ≈ 1e-4 * (
            V20.air_specific_enthalpy2D(cell.temperature) -
            V20.air_specific_enthalpy2D(Tin)
        ) rtol=1e-13
        integral = V20._air_exchange_primitive2D(
            cell.temperature, Twall,
        ) - V20._air_exchange_primitive2D(Tin, Twall)
        @test integral ≈ 0.04 / 1e-4 rtol=1e-10
    end
    @test V20._enthalpy_cell_exchange2D(
        300.0, 1000.0, 0.0, 1e-4,
    ).temperature == 300.0
    @test V20._enthalpy_cell_exchange2D(
        300.0, 1000.0, 0.04, 0.0,
    ).heat_W == 0.0
end

@testset "2D v20 architecture and short solve" begin
    p = parameters_2D_v20(
            mesh=:screen,
            source_model=:near_deep,
            deep_fraction=0.90,
            deep_length_m=0.12,
            nu_prefactor=3.693432e-4,
            reynolds_exponent=1.54346,
            power_scales=(1.05, 1.23, 0.84),
            felt_conductivity_scale=0.70,
            felt_heat_capacity_scale=0.75,
            probe_capacity_areal_J_m2_K=3000.0,
            probe_stem_conductance_areal_W_m2_K=0.0,
            t3_location=:receiver_136,
    )
    inventory = V20.network_inventory2D(p)
    @test inventory.multiplicity_sum == 100
    @test inventory.gas_thermodynamics == :integral_enthalpy
    @test V20.source_power_invariant2D(p).absolute_error < 1e-12

    times = [0.0, 1.0, 2.0]
    op = V20.OperatingCondition2D(
        irradiance=V20.linear_history(times, fill(0.0, 3)),
        flow_lpm=V20.linear_history(times, fill(10.0, 3)),
        inlet_temperature=V20.linear_history(times, fill(300.0, 3)),
        ambient_temperature=V20.linear_history(times, fill(300.0, 3)),
    )
    result = V20.simulate2D(
        p, op, times; reltol=1e-4, abstol=1e-5, dtmax=1.0,
    )
    predictions = V20.sensor_predictions2D(result)
    @test length(predictions.T3) == length(times)
    @test predictions.T3_global_z_m == 0.136
    audit = V20.enthalpy_transport_audit2D(result)
    @test audit.maximum_absolute_residual_W < 1e-8
    @test all(result.flow_solver_converged)
end
