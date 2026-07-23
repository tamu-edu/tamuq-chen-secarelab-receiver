using Test

include(joinpath(@__DIR__, "..", "1D_v29.jl"))

@testset "1D_v29 rear-reservoir energy-accounting smoke test" begin
    @test length(pnew_v29) == 24
    @test H_FRONT_FIXED_v29 == 10.0
    @test fit_heat_transfer_indices_v29 == [1, 2, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
    @test fit_source_indices_v29 == [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
    @test fit_power_scale_indices_v29 == [7, 8, 9]
    @test isempty(fit_radiation_indices_v29)
    @test pnew_v29[1] > 0.0
    @test pnew_v29[3] == PRANDTL_EXPONENT_FIXED_v29
    @test 0.1 <= pnew_v29[4] <= 1.0
    @test pnew_v29[5] >= 0.0
    @test pnew_v29[6] == 1.0
    @test all(pnew_v29[7:9] .> 0.0)
    @test pnew_v29[10] > 0.0
    @test pnew_v29[11] > 0.0
    @test pnew_v29[12] >= 0.0
    @test pnew_v29[13] > 0.0
    @test 0.0 <= perimeter_spillover_capture_v29(pnew_v29) <= 1.0
    @test 0.0 <= perimeter_source_fraction_v29(pnew_v29) <= 0.80
    @test 0.0 < flux_receiver_fraction_v29() < 1.0
    @test 0.0 < flux_spillover_fraction_v29() < 1.0
    @test isapprox(
        flux_receiver_fraction_v29() + flux_spillover_fraction_v29(),
        1.0; atol=1e-12,
    )
    @test perimeter_source_fraction_v29(pnew_v29) <= flux_spillover_fraction_v29()
    @test isapprox(
        modeled_absorbed_receiver_power_v29(456000.0, pnew_v29),
        ETA_ABS_FIXED_v29 * pnew_v29[7] * 456000.0 * A_frt;
        atol=1e-10,
    )
    @test isapprox(
        modeled_participating_absorbed_power_v29(456000.0, pnew_v29),
        core_absorbed_power_v29(456000.0, pnew_v29) +
        perimeter_absorbed_power_v29(456000.0, pnew_v29);
        atol=1e-10,
    )
    @test perimeter_source_attenuation_v29(pnew_v29) >= 0.0
    @test 0.0 <= rear_core_fraction_v29(pnew_v29) <= 1.0
    @test core_tube_fraction_v29(pnew_v29) == rear_core_fraction_v29(pnew_v29)
    @test 0.10 <= flange_loss_scale_v29(pnew_v29) <= 20.0
    @test 0.0 <= flange_cooling_gain_v29(pnew_v29) <= 50.0
    @test 1.0 <= flange_cooling_time_constant_v29(pnew_v29) <= 2000.0
    @test 0.0 <= core_axial_conduction_scale_v29(pnew_v29) <= 0.50
    @test 50.0 <= rear_reservoir_heat_capacity_v29(pnew_v29) <= 500.0
    @test 0.0 <= receiver_rear_conductance_v29(pnew_v29) <= 50.0
    @test 0.0 <= rear_tube_conductance_v29(pnew_v29) <= 50.0
    @test 0.0 <= rear_cavity_conductance_v29(pnew_v29) <= 10.0
    @test pnew_v29[20] >= 0.0
    @test lamp_off_gate_v29(0.0, 0.0, pnew_v29) == 0.0
    @test lamp_off_gate_v29(1000.0, 0.0, pnew_v29) > lamp_off_gate_v29(10.0, 0.0, pnew_v29)
    @test lamp_off_gate_v29(1000.0, 456000.0, pnew_v29) < 1e-8
    @test effective_flange_loss_scale_v29(pnew_v29, 1000.0, 0.0) > flange_loss_scale_v29(pnew_v29)
    
    @test irradiance_level_index_v29(456000.0) == 1
    @test irradiance_level_index_v29(304000.0) == 2
    @test irradiance_level_index_v29(256000.0) == 3
    @test power_scale_index_v29(456000.0) == 7
    @test power_scale_index_v29(304000.0) == 8
    @test power_scale_index_v29(256000.0) == 9
    @test absorbed_power_scale_v29(456000.0, pnew_v29) == pnew_v29[7]
    @test isapprox(STANDARD_AIR_DENSITY_v29, 1.199; atol=1e-3)
    
    @test active_flow_fraction_v29(50.0, pnew_v29) == pnew_v29[4]
    @test active_flow_fraction_v29(100.0, pnew_v29) >= active_flow_fraction_v29(50.0, pnew_v29)
    @test 0.1 <= receiver_active_flow_fraction_v29(10.0, 295.0, pnew_v29) <= 1.0
    @test receiver_channel_nusselt_v29(0.01, 50.0, 0.7, pnew_v29) >= NU_FLOOR_v29
    @test isapprox(sum(axial_exponential_weights_v29(40.0, 11, L / 11)), 1.0; atol=1e-12)

    @test isapprox(CAVITY_OUTER_RADIUS_v29, 0.075)
    @test isapprox(INSULATION_OUTER_RADIUS_v29, 0.057)
    @test isapprox(ADAPTOR_DIAMETER_v29, 0.0776)
    @test RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v29 > 0.0
    @test CAVITY_HEAT_CAPACITY_v29 > 0.0
    @test isapprox(REAR_TUBE_LENGTH_v29, 0.150)
    @test isapprox(REAR_TUBE_CAVITY_LENGTH_v29, 0.046)
    @test isapprox(T3_SAMPLE_POSITION_v29, 0.140)
    @test isapprox(COOLING_ROOM_TEMPERATURE_v29, 290.15)
    @test ADAPTOR_CONTACT_CONDUCTANCE_REFERENCE_v29 > 0.0
    
    @test core_perimeter_conductance_per_length_v29(pnew_v29) == pnew_v29[10]
    @test perimeter_heat_capacity_total_v29(pnew_v29) == pnew_v29[11]
    @test perimeter_axial_conductivity_v29(900.0, pnew_v29) == pnew_v29[12]
    @test participating_assembly_heat_capacity_v29(pnew_v29, 11) > pnew_v29[11]
    @test participating_total_heat_capacity_v29(pnew_v29, 11) >
          participating_assembly_heat_capacity_v29(pnew_v29, 11)

    # Simulation smoke test on E74
    times = [0.0, 100.0, 500.0]
    operating = (
        flow_lpm = t -> 10.0,
        irradiance = t -> 304000.0,
        inlet_temperature = t -> 295.0,
        ambient_temperature = t -> 295.0,
    )
    result = simulate_v29(pnew_v29, operating, times; nodes=11)
    
    @test size(result.core_temperature) == (11, 3)
    @test size(result.perim_temperature) == (11, 3)
    @test size(result.rear_tube_temperature) == (REAR_TUBE_DEFAULT_NODES_v29, 3)
    @test length(result.cavity_temperature) == 3
    @test length(result.rear_reservoir_temperature) == 3
    @test all(isfinite, result.core_temperature)
    @test all(isfinite, result.perim_temperature)
    @test all(isfinite, result.gas_temperature)
    @test all(isfinite, result.cavity_temperature)
    @test all(isfinite, result.rear_reservoir_temperature)
    @test all(isfinite, result.rear_tube_gas_heat)
    @test all((0.1 .<= result.active_flow_fraction) .& (result.active_flow_fraction .<= 1.0))
    @test all(result.active_flow_fraction .== 1.0)
    @test isapprox(sum(result.core_source_weights), 1.0; atol=1e-12)
    @test isapprox(sum(result.perim_source_weights), 1.0; atol=1e-12)
    @test successful_retcode(result.ode_solution)

    outputs_cooling, _, experiment_cooling = solve_case_v29(
        pnew_v29, first(sim_key_cool); is_cooling=true, nodes=11,
    )
    @test outputs_cooling[1, 1:7] == experiment_cooling[1, 1:7]
end





