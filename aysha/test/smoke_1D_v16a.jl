using Test

include(joinpath(@__DIR__, "..", "1D_v16a.jl"))

@testset "1D_v16a predicted-cavity smoke test" begin
    @test length(pnew_v16a) == 12
    @test H_FRONT_FIXED_v16a == 10.0
    @test fit_heat_transfer_indices_v16a == [1, 2, 10, 11, 12]
    @test isempty(fit_source_indices_v16a)
    @test isempty(fit_power_scale_indices_v16a)
    @test isempty(fit_radiation_indices_v16a)
    @test pnew_v16a[1] > 0.0
    @test pnew_v16a[3] == PRANDTL_EXPONENT_FIXED_v16a
    @test pnew_v16a[4] > 0.0
    @test BETA_OPT_MIN_v16a < pnew_v16a[5] < BETA_OPT_MAX_v16a
    @test 0.0 <= pnew_v16a[6] <= 1.0
    @test all(pnew_v16a[7:9] .> 0.0)
    @test 0.0 <= pnew_v16a[10] <= 0.95
    @test pnew_v16a[11] > 0.0
    @test pnew_v16a[12] > 0.0
    @test ETA_OPT_FIXED_v16a == 1.0
    @test FRONT_DEPOSITION_INITIAL_v16a == 1.0
    @test ETA_ABS_FIXED_v16a == ETA_ABS_FIXED_V5
    @test BETA_OPT_FIXED_v16a == BETA_OPT_FIXED_V5
    @test irradiance_level_index_v16a(456000.0) == 1
    @test irradiance_level_index_v16a(304000.0) == 2
    @test irradiance_level_index_v16a(256000.0) == 3
    @test power_scale_index_v16a(456000.0) == 7
    @test power_scale_index_v16a(304000.0) == 8
    @test power_scale_index_v16a(256000.0) == 9
    @test absorbed_power_scale_v16a(456000.0, pnew_v16a) == pnew_v16a[7]
    @test isapprox(STANDARD_AIR_DENSITY_v16a, 1.199; atol=1e-3)
    @test isapprox(
        receiver_nusselt_v16a(0.01, 80.0, 0.7, pnew_v16a),
        pnew_v16a[1] * 80.0^pnew_v16a[2] * 0.7^PRANDTL_EXPONENT_FIXED_v16a,
    )
    @test isapprox(
        receiver_nusselt_v16a(0.002, 80.0, 0.7, pnew_v16a),
        receiver_nusselt_v16a(0.100, 80.0, 0.7, pnew_v16a),
    )
    @test receiver_nusselt_v16a(0.010, 80.0, 0.7, pnew_v16a) >
          receiver_nusselt_v16a(0.010, 20.0, 0.7, pnew_v16a)
    @test rosseland_beta_v16a(pnew_v16a) > 0.0
    @test rosseland_conductivity_v16a(1000.0, pnew_v16a) > 0.0
    @test rosseland_optical_thickness_v16a(L, pnew_v16a) > 0.0
    @test isapprox(CAVITY_OUTER_RADIUS_v16a, 0.075)
    @test isapprox(INSULATION_OUTER_RADIUS_v16a, 0.057)
    @test isapprox(ADAPTOR_DIAMETER_v16a, 0.0776)
    @test RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v16a > 0.0
    @test CAVITY_HEAT_CAPACITY_v16a > 0.0
    @test isapprox(REAR_TUBE_LENGTH_v16a, 0.150)
    @test isapprox(REAR_TUBE_CAVITY_LENGTH_v16a, 0.046)
    @test isapprox(REAR_TUBE_FLANGE_LENGTH_v16a, 0.104)
    @test isapprox(T3_SAMPLE_POSITION_v16a, 0.140)
    @test RECEIVER_TO_TUBE_CONDUCTANCE_v16a > 0.0
    @test rear_tube_flange_conductance_per_length_v16a(0.100, 500.0) > 0.0
    @test isapprox(sum(solar_weights_v5(pnew_v16a[5], pnew_v16a[6], 11)), 1.0)
    @test isapprox(solar_weights_v5(pnew_v16a[5], pnew_v16a[6], 11)[1], 1.0)
    @test wall_power_fraction_v16a(pnew_v16a) == pnew_v16a[10]
    @test active_wall_conductance_per_length_v16a(pnew_v16a) == pnew_v16a[11]
    @test wall_heat_capacity_total_v16a(pnew_v16a) == pnew_v16a[12]
    @test participating_assembly_heat_capacity_v16a(pnew_v16a, 11) > pnew_v16a[12]
    @test AXIAL_VIEW_RADIATION_ROW_SUM_v16a > 0.0
    @test AXIAL_VIEW_RADIATION_LENGTH_v16a > 0.0
    @test TC_WIRE_LOSS_FRONT_FRACTION_v16a > 0.0
    view_matrix = axial_view_matrix_v16a(collect(range(L / 22, step=L / 11, length=11)))
    @test size(view_matrix) == (11, 11)
    @test all(diag(view_matrix) .== 0.0)
    @test maximum(sum(view_matrix, dims=2)) <= AXIAL_VIEW_RADIATION_ROW_SUM_v16a * (1.0 + 1e-10)

    model, result, experiment = solve_case_v16a(pnew_v16a, "E74"; nodes=11)
    @test size(model) == (length(result.time), 6)
    @test size(experiment) == (length(result.time), 5)
    @test length(result.boundary_temperature) == length(result.time)
    @test size(result.gas_temperature, 1) == 11 + REAR_TUBE_DEFAULT_NODES_v16a + 1
    @test size(result.active_temperature) == (11, length(result.time))
    @test size(result.wall_temperature) == (11, length(result.time))
    @test result.solid_temperature === result.wall_temperature
    @test size(result.rear_tube_temperature) == (REAR_TUBE_DEFAULT_NODES_v16a, length(result.time))
    @test result.cavity_temperature[1] == observation(measurements, "E74", "_T2")[1]
    @test all(isfinite, model)
    @test all(isfinite, result.boundary_temperature)
    @test all(isfinite, result.cavity_temperature)
    @test all(isfinite, result.rear_tube_temperature)
    @test all(isfinite, result.flange_heat_loss)
    @test all(isfinite, result.receiver_to_cavity_heat)
    @test all(isfinite, result.active_wall_heat)
    @test all(isfinite, result.tube_to_cavity_heat)
    @test all(isfinite, result.cavity_ambient_heat_loss)
    @test all(isfinite, result.receiver_nusselt)
    @test all(isfinite, result.receiver_effectiveness)
    @test all(isfinite, result.receiver_ua_over_mcp)
    @test all(isfinite, result.receiver_gas_heat)
    @test all(isfinite, result.axial_view_radiation_heat)
    @test all(result.axial_view_radiation_heat .== 0.0)
    @test size(result.receiver_nusselt) == (11, length(result.time))
    @test all(model[:, 6] .> 0.0)
    @test successful_retcode(result.ode_solution)

    model_rad, result_rad, experiment_rad = solve_case_v16a(
        pnew_v16a, "E74"; nodes=11, use_rosseland_radiation=true,
    )
    @test size(model_rad) == size(model)
    @test size(experiment_rad) == size(experiment)
    @test all(isfinite, model_rad)
    @test successful_retcode(result_rad.ode_solution)

    model_profile, result_profile, _ = solve_case_v16a(
        pnew_v16a, "E74"; nodes=11, use_rosseland_radiation=true,
        rosseland_profile=:front_strong,
    )
    @test size(model_profile) == size(model)
    @test all(isfinite, model_profile)
    @test successful_retcode(result_profile.ode_solution)

    model_view, result_view, experiment_view = solve_case_v16a(
        pnew_v16a, "E74"; nodes=11, use_axial_view_radiation=true,
    )
    @test size(model_view) == size(model)
    @test size(experiment_view) == size(experiment)
    @test all(isfinite, model_view)
    @test any(abs.(result_view.axial_view_radiation_heat) .> 0.0)
    @test successful_retcode(result_view.ode_solution)

    model_tc, result_tc, experiment_tc = solve_case_v16a(
        pnew_v16a, "E74"; nodes=11, use_tc_measurement_model=true,
    )
    @test size(model_tc) == size(model)
    @test size(experiment_tc) == size(experiment)
    @test all(isfinite, model_tc)
    @test model_tc[end, 1] < solid_at_v3(result_tc, sensor_positions[:T8])[end]
    @test successful_retcode(result_tc.ode_solution)

    @test isfinite(loss_heating_v16a(pnew_v16a, ["E74"]; nodes=11))
    @test isfinite(loss_cooling_v16a(pnew_v16a, ["C69"]; nodes=11))

    trial = optimize_parameter_subset_v16a(
        p -> sum(abs2, p[fit_heat_transfer_indices_v16a]),
        pnew_v16a,
        fit_heat_transfer_indices_v16a,
        lb_full_v16a,
        ub_full_v16a;
        maximum_iterations=5,
        maximum_time_seconds=5.0,
        label="test-v16a",
    )
    @test length(trial.parameters) == length(pnew_v16a)
end






