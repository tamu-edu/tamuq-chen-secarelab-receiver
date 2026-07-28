using Test

include(joinpath(@__DIR__, "..", "1D_v31.jl"))

@testset "1D_v31 distributed rear-rail energy-accounting smoke test" begin
    @test length(pnew_v31) == 25
    @test H_FRONT_FIXED_v31 == 10.0
    @test fit_rear_stage_indices_v31 == [17, 18, 19, 22, 23, 24, 25]
    @test fit_transport_stage_indices_v31 == [10, 12, 17, 18, 19, 20, 22, 23, 24, 25]
    @test fit_full_stage_indices_v31 == [1, 2, 7, 8, 9, 10, 12, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25]
    @test fit_heat_transfer_indices_v31 == fit_rear_stage_indices_v31
    @test fit_source_indices_v31 == [14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25]
    @test fit_power_scale_indices_v31 == [7, 8, 9]
    @test isempty(fit_radiation_indices_v31)
    @test fit_indices_for_stage_v31(:rear) == fit_rear_stage_indices_v31
    @test fit_indices_for_stage_v31(:transport) == fit_transport_stage_indices_v31
    @test fit_indices_for_stage_v31(:full) == fit_full_stage_indices_v31
    @test pnew_v31[1] > 0.0
    @test pnew_v31[3] == PRANDTL_EXPONENT_FIXED_v31
    @test 0.1 <= pnew_v31[4] <= 1.0
    @test pnew_v31[5] >= 0.0
    @test pnew_v31[6] == 1.0
    @test all(pnew_v31[7:9] .> 0.0)
    @test pnew_v31[10] > 0.0
    @test pnew_v31[11] > 0.0
    @test pnew_v31[12] >= 0.0
    @test pnew_v31[13] > 0.0
    @test 0.0 <= perimeter_spillover_capture_v31(pnew_v31) <= 1.0
    @test 0.0 <= perimeter_source_fraction_v31(pnew_v31) <= 0.80
    @test 0.0 < flux_receiver_fraction_v31() < 1.0
    @test 0.0 < flux_spillover_fraction_v31() < 1.0
    @test isapprox(
        flux_receiver_fraction_v31() + flux_spillover_fraction_v31(),
        1.0; atol=1e-12,
    )
    @test perimeter_source_fraction_v31(pnew_v31) <= flux_spillover_fraction_v31()
    @test isapprox(
        modeled_absorbed_receiver_power_v31(456000.0, pnew_v31),
        ETA_ABS_FIXED_v31 * pnew_v31[7] * 456000.0 * A_frt;
        atol=1e-10,
    )
    @test isapprox(
        modeled_participating_absorbed_power_v31(456000.0, pnew_v31),
        core_absorbed_power_v31(456000.0, pnew_v31) +
        perimeter_absorbed_power_v31(456000.0, pnew_v31);
        atol=1e-10,
    )
    @test perimeter_source_attenuation_v31(pnew_v31) >= 0.0
    @test 0.0 <= rear_core_fraction_v31(pnew_v31) <= 1.0
    @test core_tube_fraction_v31(pnew_v31) == rear_core_fraction_v31(pnew_v31)
    @test 0.10 <= flange_loss_scale_v31(pnew_v31) <= 20.0
    @test 0.0 <= flange_cooling_gain_v31(pnew_v31) <= 50.0
    @test 1.0 <= flange_cooling_time_constant_v31(pnew_v31) <= 2000.0
    @test 0.0 <= core_axial_conduction_scale_v31(pnew_v31) <= 0.50
    @test 80.0 <= rear_reservoir_heat_capacity_v31(pnew_v31) <= 150.0
    @test 0.0 <= receiver_rear_conductance_v31(pnew_v31) <= 50.0
    @test 0.0 <= rear_tube_conductance_v31(pnew_v31) <= 50.0
    @test 0.0 <= rear_cavity_conductance_v31(pnew_v31) <= 10.0
    @test 0.0 <= rear_axial_conductance_v31(pnew_v31) <= 100.0
    @test pnew_v31[20] >= 0.0
    @test lamp_off_gate_v31(0.0, 0.0, pnew_v31) == 0.0
    @test lamp_off_gate_v31(1000.0, 0.0, pnew_v31) > lamp_off_gate_v31(10.0, 0.0, pnew_v31)
    @test lamp_off_gate_v31(1000.0, 456000.0, pnew_v31) < 1e-8
    @test effective_flange_loss_scale_v31(pnew_v31, 1000.0, 0.0) > flange_loss_scale_v31(pnew_v31)
    
    @test irradiance_level_index_v31(456000.0) == 1
    @test irradiance_level_index_v31(304000.0) == 2
    @test irradiance_level_index_v31(256000.0) == 3
    @test power_scale_index_v31(456000.0) == 7
    @test power_scale_index_v31(304000.0) == 8
    @test power_scale_index_v31(256000.0) == 9
    @test absorbed_power_scale_v31(456000.0, pnew_v31) == pnew_v31[7]
    @test isapprox(STANDARD_AIR_DENSITY_v31, 1.199; atol=1e-3)
    
    @test active_flow_fraction_v31(50.0, pnew_v31) == pnew_v31[4]
    @test active_flow_fraction_v31(100.0, pnew_v31) >= active_flow_fraction_v31(50.0, pnew_v31)
    @test 0.1 <= receiver_active_flow_fraction_v31(10.0, 295.0, pnew_v31) <= 1.0
    @test receiver_channel_nusselt_v31(0.01, 50.0, 0.7, pnew_v31) >= NU_FLOOR_v31
    @test isapprox(sum(axial_exponential_weights_v31(40.0, 11, L / 11)), 1.0; atol=1e-12)
    rear_weights = rear_contact_weights_v31(collect(range(L / 22, step=L / 11, length=11)))
    @test length(rear_weights) == 11
    @test isapprox(sum(rear_weights), 1.0; atol=1e-12)
    @test all(rear_weights .> 0.0)
    @test rear_weights[end] > rear_weights[1]

    @test isapprox(CAVITY_OUTER_RADIUS_v31, 0.075)
    @test isapprox(INSULATION_OUTER_RADIUS_v31, 0.057)
    @test isapprox(ADAPTOR_DIAMETER_v31, 0.0776)
    @test RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v31 > 0.0
    @test CAVITY_HEAT_CAPACITY_v31 > 0.0
    @test isapprox(REAR_TUBE_LENGTH_v31, 0.150)
    @test isapprox(REAR_TUBE_CAVITY_LENGTH_v31, 0.046)
    @test isapprox(T3_SAMPLE_POSITION_v31, 0.140)
    @test isapprox(COOLING_ROOM_TEMPERATURE_v31, 290.15)
    @test ADAPTOR_CONTACT_CONDUCTANCE_REFERENCE_v31 > 0.0
    
    @test core_perimeter_conductance_per_length_v31(pnew_v31) == pnew_v31[10]
    @test 150.0 <= perimeter_heat_capacity_total_v31(pnew_v31) <= 230.0
    @test perimeter_heat_capacity_total_v31(pnew_v31) == pnew_v31[11]
    @test perimeter_axial_conductivity_v31(900.0, pnew_v31) == pnew_v31[12]
    @test participating_assembly_heat_capacity_v31(pnew_v31, 11) > pnew_v31[11]
    @test participating_total_heat_capacity_v31(pnew_v31, 11) >
          participating_assembly_heat_capacity_v31(pnew_v31, 11)
    @test 278.0 <= participating_total_heat_capacity_v31(pnew_v31, 11) <= 324.0

    # Simulation smoke test on E74
    times = [0.0, 100.0, 500.0]
    operating = (
        flow_lpm = t -> 10.0,
        irradiance = t -> 304000.0,
        inlet_temperature = t -> 295.0,
        ambient_temperature = t -> 295.0,
    )
    result = simulate_v31(pnew_v31, operating, times; nodes=11)
    
    @test size(result.core_temperature) == (11, 3)
    @test size(result.perim_temperature) == (11, 3)
    @test size(result.rear_tube_temperature) == (REAR_TUBE_DEFAULT_NODES_v31, 3)
    @test length(result.cavity_temperature) == 3
    @test size(result.rear_reservoir_temperature) == (11, 3)
    @test all(isfinite, result.core_temperature)
    @test all(isfinite, result.perim_temperature)
    @test all(isfinite, result.gas_temperature)
    @test all(isfinite, result.cavity_temperature)
    @test all(isfinite, result.rear_reservoir_temperature)
    @test isapprox(sum(result.rear_contact_weights), 1.0; atol=1e-12)
    @test all(isfinite, result.rear_tube_gas_heat)
    @test all((0.1 .<= result.active_flow_fraction) .& (result.active_flow_fraction .<= 1.0))
    @test all(result.active_flow_fraction .== 1.0)
    @test isapprox(sum(result.core_source_weights), 1.0; atol=1e-12)
    @test isapprox(sum(result.perim_source_weights), 1.0; atol=1e-12)
    @test successful_retcode(result.ode_solution)

    outputs_cooling, _, experiment_cooling = solve_case_v31(
        pnew_v31, first(sim_key_cool); is_cooling=true, nodes=11,
    )
    @test outputs_cooling[1, 1:7] == experiment_cooling[1, 1:7]
end







