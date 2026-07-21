using Test

include(joinpath(@__DIR__, "..", "1D_v11.jl"))

@testset "1D_v11 predicted-cavity smoke test" begin
    @test length(pnew_v11) == 5
    @test H_FRONT_FIXED_v11 == 10.0
    @test fit_heat_transfer_indices_v11 == [1, 2, 3, 4]
    @test isempty(fit_radiation_indices_v11)
    @test pnew_v11[1] > 0.0
    @test 0.2 <= pnew_v11[3] <= 0.5
    @test pnew_v11[4] > 0.0
    @test ETA_OPT_FIXED_v11 == 1.0
    @test FRONT_DEPOSITION_STAGE1_v11 == 1.0
    @test ETA_ABS_FIXED_v11 == ETA_ABS_FIXED_V5
    @test BETA_OPT_FIXED_v11 == BETA_OPT_FIXED_V5
    @test NU_FD_RECEIVER_v11 > 0.0
    @test receiver_nusselt_v11(0.01, 40.0, 0.7, pnew_v11) >= NU_FD_RECEIVER_v11
    @test receiver_nusselt_v11(0.002, 80.0, 0.7, pnew_v11) >
          receiver_nusselt_v11(0.100, 80.0, 0.7, pnew_v11)
    @test receiver_nusselt_v11(0.010, 80.0, 0.7, pnew_v11) >
          receiver_nusselt_v11(0.010, 20.0, 0.7, pnew_v11)
    @test rosseland_beta_v11(pnew_v11) > 0.0
    @test rosseland_conductivity_v11(1000.0, pnew_v11) > 0.0
    @test rosseland_optical_thickness_v11(L, pnew_v11) > 0.0
    @test isapprox(CAVITY_OUTER_RADIUS_v11, 0.075)
    @test isapprox(INSULATION_OUTER_RADIUS_v11, 0.057)
    @test isapprox(ADAPTOR_DIAMETER_v11, 0.0776)
    @test RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v11 > 0.0
    @test CAVITY_HEAT_CAPACITY_v11 > 0.0
    @test isapprox(REAR_TUBE_LENGTH_v11, 0.150)
    @test isapprox(REAR_TUBE_CAVITY_LENGTH_v11, 0.046)
    @test isapprox(REAR_TUBE_FLANGE_LENGTH_v11, 0.104)
    @test isapprox(T3_SAMPLE_POSITION_v11, 0.140)
    @test RECEIVER_TO_TUBE_CONDUCTANCE_v11 > 0.0
    @test rear_tube_flange_conductance_per_length_v11(0.100, 500.0) > 0.0
    @test isapprox(sum(solar_weights_v5(BETA_OPT_FIXED_v11, FRONT_DEPOSITION_STAGE1_v11, 11)), 1.0)
    @test solar_weights_v5(BETA_OPT_FIXED_v11, FRONT_DEPOSITION_STAGE1_v11, 11)[1] == 1.0

    model, result, experiment = solve_case_v11(pnew_v11, "E74"; nodes=11)
    @test size(model) == (length(result.time), 6)
    @test size(experiment) == (length(result.time), 5)
    @test length(result.boundary_temperature) == length(result.time)
    @test size(result.gas_temperature, 1) == 11 + REAR_TUBE_DEFAULT_NODES_v11 + 1
    @test size(result.rear_tube_temperature) == (REAR_TUBE_DEFAULT_NODES_v11, length(result.time))
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

    model_rad, result_rad, experiment_rad = solve_case_v11(
        pnew_v11, "E74"; nodes=11, use_rosseland_radiation=true,
    )
    @test size(model_rad) == size(model)
    @test size(experiment_rad) == size(experiment)
    @test all(isfinite, model_rad)
    @test successful_retcode(result_rad.ode_solution)

    model_profile, result_profile, _ = solve_case_v11(
        pnew_v11, "E74"; nodes=11, use_rosseland_radiation=true,
        rosseland_profile=:front_strong,
    )
    @test size(model_profile) == size(model)
    @test all(isfinite, model_profile)
    @test successful_retcode(result_profile.ode_solution)

    @test isfinite(loss_heating_v11(pnew_v11, ["E74"]; nodes=11))
    @test isfinite(loss_cooling_v11(pnew_v11, ["C69"]; nodes=11))

    trial = optimize_parameter_subset_v11(
        p -> sum(abs2, p[fit_heat_transfer_indices_v11]),
        pnew_v11,
        fit_heat_transfer_indices_v11,
        lb_full_v11,
        ub_full_v11;
        maximum_iterations=5,
        maximum_time_seconds=5.0,
        label="test-v11",
    )
    @test length(trial.parameters) == length(pnew_v11)
end
