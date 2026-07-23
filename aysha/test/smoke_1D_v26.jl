using Test

include(joinpath(@__DIR__, "..", "1D_v26.jl"))

@testset "1D_v26 energy-accounting rear-loss smoke test" begin
    @test length(pnew_v26) == 21
    @test H_FRONT_FIXED_v26 == 10.0
    @test fit_heat_transfer_indices_v26 == [1, 2, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20]
    @test fit_source_indices_v26 == [14, 15, 16, 17, 18, 19, 20]
    @test fit_power_scale_indices_v26 == [7, 8, 9]
    @test isempty(fit_radiation_indices_v26)
    @test pnew_v26[1] > 0.0
    @test pnew_v26[3] == PRANDTL_EXPONENT_FIXED_v26
    @test 0.1 <= pnew_v26[4] <= 1.0
    @test pnew_v26[5] >= 0.0
    @test pnew_v26[6] == 1.0
    @test all(pnew_v26[7:9] .> 0.0)
    @test pnew_v26[10] > 0.0
    @test pnew_v26[11] > 0.0
    @test pnew_v26[12] >= 0.0
    @test pnew_v26[13] > 0.0
    @test 0.0 <= perimeter_spillover_capture_v26(pnew_v26) <= 1.0
    @test 0.0 <= perimeter_source_fraction_v26(pnew_v26) <= 0.80
    @test 0.0 < flux_receiver_fraction_v26() < 1.0
    @test 0.0 < flux_spillover_fraction_v26() < 1.0
    @test isapprox(
        flux_receiver_fraction_v26() + flux_spillover_fraction_v26(),
        1.0; atol=1e-12,
    )
    @test perimeter_source_fraction_v26(pnew_v26) <= flux_spillover_fraction_v26()
    @test isapprox(
        modeled_absorbed_receiver_power_v26(456000.0, pnew_v26),
        ETA_ABS_FIXED_v26 * pnew_v26[7] * 456000.0 * A_frt;
        atol=1e-10,
    )
    @test isapprox(
        modeled_participating_absorbed_power_v26(456000.0, pnew_v26),
        core_absorbed_power_v26(456000.0, pnew_v26) +
        perimeter_absorbed_power_v26(456000.0, pnew_v26);
        atol=1e-10,
    )
    @test perimeter_source_attenuation_v26(pnew_v26) >= 0.0
    @test 0.0 <= core_tube_fraction_v26(pnew_v26) <= 1.0
    @test 0.10 <= flange_loss_scale_v26(pnew_v26) <= 20.0
    @test rear_core_loss_per_length_v26(pnew_v26) >= 0.0
    @test rear_perim_loss_per_length_v26(pnew_v26) >= 0.0
    @test 0.25 <= rear_sink_shape_v26(pnew_v26) <= 1.0
    @test 0.0 <= core_axial_conduction_scale_v26(pnew_v26) <= 0.125
    @test pnew_v26[21] == 0.0
    
    @test irradiance_level_index_v26(456000.0) == 1
    @test irradiance_level_index_v26(304000.0) == 2
    @test irradiance_level_index_v26(256000.0) == 3
    @test power_scale_index_v26(456000.0) == 7
    @test power_scale_index_v26(304000.0) == 8
    @test power_scale_index_v26(256000.0) == 9
    @test absorbed_power_scale_v26(456000.0, pnew_v26) == pnew_v26[7]
    @test isapprox(STANDARD_AIR_DENSITY_v26, 1.199; atol=1e-3)
    
    @test active_flow_fraction_v26(50.0, pnew_v26) == pnew_v26[4]
    @test active_flow_fraction_v26(100.0, pnew_v26) >= active_flow_fraction_v26(50.0, pnew_v26)
    @test 0.1 <= receiver_active_flow_fraction_v26(10.0, 295.0, pnew_v26) <= 1.0
    @test receiver_channel_nusselt_v26(0.01, 50.0, 0.7, pnew_v26) >= NU_FLOOR_v26
    @test isapprox(sum(axial_exponential_weights_v26(40.0, 11, L / 11)), 1.0; atol=1e-12)
    @test rear_receiver_sink_weight_v26(sensor_positions[:T8], pnew_v26) == 0.0
    @test rear_receiver_sink_weight_v26(sensor_positions[:T10], pnew_v26) >
          rear_receiver_sink_weight_v26(sensor_positions[:T9], pnew_v26)

    @test isapprox(CAVITY_OUTER_RADIUS_v26, 0.075)
    @test isapprox(INSULATION_OUTER_RADIUS_v26, 0.057)
    @test isapprox(ADAPTOR_DIAMETER_v26, 0.0776)
    @test RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v26 > 0.0
    @test CAVITY_HEAT_CAPACITY_v26 > 0.0
    @test isapprox(REAR_TUBE_LENGTH_v26, 0.150)
    @test isapprox(REAR_TUBE_CAVITY_LENGTH_v26, 0.046)
    @test isapprox(T3_SAMPLE_POSITION_v26, L)
    @test isapprox(COOLING_ROOM_TEMPERATURE_v26, 290.15)
    @test RECEIVER_TO_TUBE_CONDUCTANCE_v26 > 0.0
    
    @test core_perimeter_conductance_per_length_v26(pnew_v26) == pnew_v26[10]
    @test perimeter_heat_capacity_total_v26(pnew_v26) == pnew_v26[11]
    @test perimeter_axial_conductivity_v26(900.0, pnew_v26) == pnew_v26[12]
    @test participating_assembly_heat_capacity_v26(pnew_v26, 11) > pnew_v26[11]

    # Simulation smoke test on E74
    times = [0.0, 100.0, 500.0]
    operating = (
        flow_lpm = t -> 10.0,
        irradiance = t -> 304000.0,
        inlet_temperature = t -> 295.0,
        ambient_temperature = t -> 295.0,
    )
    result = simulate_v26(pnew_v26, operating, times; nodes=11)
    
    @test size(result.core_temperature) == (11, 3)
    @test size(result.perim_temperature) == (11, 3)
    @test size(result.rear_tube_temperature) == (REAR_TUBE_DEFAULT_NODES_v26, 3)
    @test length(result.cavity_temperature) == 3
    @test all(isfinite, result.core_temperature)
    @test all(isfinite, result.perim_temperature)
    @test all(isfinite, result.gas_temperature)
    @test all(isfinite, result.cavity_temperature)
    @test all(isfinite, result.rear_tube_gas_heat)
    @test all((0.1 .<= result.active_flow_fraction) .& (result.active_flow_fraction .<= 1.0))
    @test all(result.active_flow_fraction .== 1.0)
    @test isapprox(sum(result.core_source_weights), 1.0; atol=1e-12)
    @test isapprox(sum(result.perim_source_weights), 1.0; atol=1e-12)
    @test successful_retcode(result.ode_solution)

    outputs_cooling, _, experiment_cooling = solve_case_v26(
        pnew_v26, first(sim_key_cool); is_cooling=true, nodes=11,
    )
    @test outputs_cooling[1, 1:7] == experiment_cooling[1, 1:7]
end



