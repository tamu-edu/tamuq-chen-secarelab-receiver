using Test

include(joinpath(@__DIR__, "..", "1D_v12a.jl"))

@testset "1D_v12a predicted-cavity smoke test" begin
    @test length(pnew_v12a) == 6
    @test H_FRONT_FIXED_v12a == 10.0
    @test isempty(fit_heat_transfer_indices_v12a)
    @test fit_source_indices_v12a == [5, 6]
    @test isempty(fit_radiation_indices_v12a)
    @test pnew_v12a[1] > 0.0
    @test 0.2 <= pnew_v12a[3] <= 0.5
    @test pnew_v12a[4] > 0.0
    @test BETA_OPT_MIN_v12a < pnew_v12a[5] < BETA_OPT_MAX_v12a
    @test 0.0 <= pnew_v12a[6] <= 1.0
    @test ETA_OPT_FIXED_v12a == 1.0
    @test FRONT_DEPOSITION_INITIAL_v12a == 1.0
    @test ETA_ABS_FIXED_v12a == ETA_ABS_FIXED_V5
    @test BETA_OPT_FIXED_v12a == BETA_OPT_FIXED_V5
    @test NU_FD_RECEIVER_v12a > 0.0
    @test receiver_nusselt_v12a(0.01, 40.0, 0.7, pnew_v12a) >= NU_FD_RECEIVER_v12a
    @test receiver_nusselt_v12a(0.002, 80.0, 0.7, pnew_v12a) >
          receiver_nusselt_v12a(0.100, 80.0, 0.7, pnew_v12a)
    @test receiver_nusselt_v12a(0.010, 80.0, 0.7, pnew_v12a) >
          receiver_nusselt_v12a(0.010, 20.0, 0.7, pnew_v12a)
    @test rosseland_beta_v12a(pnew_v12a) > 0.0
    @test rosseland_conductivity_v12a(1000.0, pnew_v12a) > 0.0
    @test rosseland_optical_thickness_v12a(L, pnew_v12a) > 0.0
    @test isapprox(CAVITY_OUTER_RADIUS_v12a, 0.075)
    @test isapprox(INSULATION_OUTER_RADIUS_v12a, 0.057)
    @test isapprox(ADAPTOR_DIAMETER_v12a, 0.0776)
    @test RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v12a > 0.0
    @test CAVITY_HEAT_CAPACITY_v12a > 0.0
    @test isapprox(REAR_TUBE_LENGTH_v12a, 0.150)
    @test isapprox(REAR_TUBE_CAVITY_LENGTH_v12a, 0.046)
    @test isapprox(REAR_TUBE_FLANGE_LENGTH_v12a, 0.104)
    @test isapprox(T3_SAMPLE_POSITION_v12a, 0.140)
    @test RECEIVER_TO_TUBE_CONDUCTANCE_v12a > 0.0
    @test rear_tube_flange_conductance_per_length_v12a(0.100, 500.0) > 0.0
    @test isapprox(sum(solar_weights_v5(pnew_v12a[5], pnew_v12a[6], 11)), 1.0)
    @test solar_weights_v5(pnew_v12a[5], pnew_v12a[6], 11)[1] == 1.0

    model, result, experiment = solve_case_v12a(pnew_v12a, "E74"; nodes=11)
    @test size(model) == (length(result.time), 6)
    @test size(experiment) == (length(result.time), 5)
    @test length(result.boundary_temperature) == length(result.time)
    @test size(result.gas_temperature, 1) == 11 + REAR_TUBE_DEFAULT_NODES_v12a + 1
    @test size(result.rear_tube_temperature) == (REAR_TUBE_DEFAULT_NODES_v12a, length(result.time))
    @test result.cavity_temperature[1] == observation(measurements, "E74", "_T2")[1]
    @test all(isfinite, model)
    @test all(isfinite, result.boundary_temperature)
    @test all(isfinite, result.cavity_temperature)
    @test all(isfinite, result.rear_tube_temperature)
    @test all(isfinite, result.flange_heat_loss)
    @test all(isfinite, result.receiver_to_cavity_heat)
    @test all(isfinite, result.tube_to_cavity_heat)
    @test all(isfinite, result.cavity_ambient_heat_loss)
    @test all(isfinite, result.receiver_nusselt)
    @test all(isfinite, result.receiver_effectiveness)
    @test all(isfinite, result.receiver_ua_over_mcp)
    @test all(isfinite, result.receiver_gas_heat)
    @test size(result.receiver_nusselt) == (11, length(result.time))
    @test all(model[:, 6] .> 0.0)
    @test successful_retcode(result.ode_solution)

    model_rad, result_rad, experiment_rad = solve_case_v12a(
        pnew_v12a, "E74"; nodes=11, use_rosseland_radiation=true,
    )
    @test size(model_rad) == size(model)
    @test size(experiment_rad) == size(experiment)
    @test all(isfinite, model_rad)
    @test successful_retcode(result_rad.ode_solution)

    model_profile, result_profile, _ = solve_case_v12a(
        pnew_v12a, "E74"; nodes=11, use_rosseland_radiation=true,
        rosseland_profile=:front_strong,
    )
    @test size(model_profile) == size(model)
    @test all(isfinite, model_profile)
    @test successful_retcode(result_profile.ode_solution)

    @test isfinite(loss_heating_v12a(pnew_v12a, ["E74"]; nodes=11))
    @test isfinite(loss_cooling_v12a(pnew_v12a, ["C69"]; nodes=11))

    trial = optimize_parameter_subset_v12a(
        p -> sum(abs2, p[fit_source_indices_v12a]),
        pnew_v12a,
        fit_source_indices_v12a,
        lb_full_v12a,
        ub_full_v12a;
        maximum_iterations=5,
        maximum_time_seconds=5.0,
        label="test-v12a",
    )
    @test length(trial.parameters) == length(pnew_v12a)
end

