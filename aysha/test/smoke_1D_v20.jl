using Test

include(joinpath(@__DIR__, "..", "1D_v20.jl"))

@testset "1D_v20 energy-accounting 2-Zone Core/Perimeter smoke test" begin
    @test length(pnew_v20) == 16
    @test H_FRONT_FIXED_v20 == 10.0
    @test fit_heat_transfer_indices_v20 == [1, 2, 4, 5, 7, 8, 9, 10, 11, 12, 14, 15, 16]
    @test fit_source_indices_v20 == [14, 15, 16]
    @test fit_power_scale_indices_v20 == [7, 8, 9]
    @test isempty(fit_radiation_indices_v20)
    @test pnew_v20[1] > 0.0
    @test pnew_v20[3] == PRANDTL_EXPONENT_FIXED_v20
    @test 0.1 <= pnew_v20[4] <= 1.0
    @test pnew_v20[5] >= 0.0
    @test pnew_v20[6] == 1.0
    @test all(pnew_v20[7:9] .> 0.0)
    @test pnew_v20[10] > 0.0
    @test pnew_v20[11] > 0.0
    @test pnew_v20[12] >= 0.0
    @test pnew_v20[13] > 0.0
    @test 0.0 <= perimeter_source_fraction_v20(pnew_v20) <= 0.80
    @test perimeter_source_attenuation_v20(pnew_v20) >= 0.0
    @test 0.0 <= core_tube_fraction_v20(pnew_v20) <= 1.0
    
    @test irradiance_level_index_v20(456000.0) == 1
    @test irradiance_level_index_v20(304000.0) == 2
    @test irradiance_level_index_v20(256000.0) == 3
    @test power_scale_index_v20(456000.0) == 7
    @test power_scale_index_v20(304000.0) == 8
    @test power_scale_index_v20(256000.0) == 9
    @test absorbed_power_scale_v20(456000.0, pnew_v20) == pnew_v20[7]
    @test isapprox(STANDARD_AIR_DENSITY_v20, 1.199; atol=1e-3)
    
    @test active_flow_fraction_v20(50.0, pnew_v20) == pnew_v20[4]
    @test active_flow_fraction_v20(100.0, pnew_v20) >= active_flow_fraction_v20(50.0, pnew_v20)
    @test 0.1 <= receiver_active_flow_fraction_v20(10.0, 295.0, pnew_v20) <= 1.0
    @test receiver_channel_nusselt_v20(0.01, 50.0, 0.7, pnew_v20) >= NU_FLOOR_v20
    @test isapprox(sum(axial_exponential_weights_v20(40.0, 11, L / 11)), 1.0; atol=1e-12)

    @test isapprox(CAVITY_OUTER_RADIUS_v20, 0.075)
    @test isapprox(INSULATION_OUTER_RADIUS_v20, 0.057)
    @test isapprox(ADAPTOR_DIAMETER_v20, 0.0776)
    @test RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v20 > 0.0
    @test CAVITY_HEAT_CAPACITY_v20 > 0.0
    @test isapprox(REAR_TUBE_LENGTH_v20, 0.150)
    @test isapprox(REAR_TUBE_CAVITY_LENGTH_v20, 0.046)
    @test isapprox(T3_SAMPLE_POSITION_v20, 0.140)
    @test RECEIVER_TO_TUBE_CONDUCTANCE_v20 > 0.0
    
    @test core_perimeter_conductance_per_length_v20(pnew_v20) == pnew_v20[10]
    @test perimeter_heat_capacity_total_v20(pnew_v20) == pnew_v20[11]
    @test perimeter_axial_conductivity_v20(900.0, pnew_v20) == pnew_v20[12]
    @test participating_assembly_heat_capacity_v20(pnew_v20, 11) > pnew_v20[11]

    # Simulation smoke test on E74
    times = [0.0, 100.0, 500.0]
    operating = (
        flow_lpm = t -> 10.0,
        irradiance = t -> 304000.0,
        inlet_temperature = t -> 295.0,
        ambient_temperature = t -> 295.0,
    )
    result = simulate_v20(pnew_v20, operating, times; nodes=11)
    
    @test size(result.core_temperature) == (11, 3)
    @test size(result.perim_temperature) == (11, 3)
    @test size(result.rear_tube_temperature) == (REAR_TUBE_DEFAULT_NODES_v20, 3)
    @test length(result.cavity_temperature) == 3
    @test all(isfinite, result.core_temperature)
    @test all(isfinite, result.perim_temperature)
    @test all(isfinite, result.gas_temperature)
    @test all(isfinite, result.cavity_temperature)
    @test all(isfinite, result.rear_tube_gas_heat)
    @test all((0.1 .<= result.active_flow_fraction) .& (result.active_flow_fraction .<= 1.0))
    @test isapprox(sum(result.core_source_weights), 1.0; atol=1e-12)
    @test isapprox(sum(result.perim_source_weights), 1.0; atol=1e-12)
    @test successful_retcode(result.ode_solution)
end
