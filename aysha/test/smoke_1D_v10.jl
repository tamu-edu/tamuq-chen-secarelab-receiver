using Test

include(joinpath(@__DIR__, "..", "1D_v10.jl"))

@testset "1D_v10 predicted-cavity smoke test" begin
    @test length(pnew_v10) == 5
    @test H_FRONT_FIXED_v10 == 10.0
    @test length(fit_heat_transfer_indices_v10) == 1
    @test length(fit_radiation_indices_v10) == 1
    @test pnew_v10[3] == 1.0
    @test pnew_v10[4] == 1.0
    @test lb_full_v10[3] == ub_full_v10[3] == 1.0
    @test lb_full_v10[4] == ub_full_v10[4] == 1.0
    @test NU_FD_RECEIVER_v10 > 0.0
    @test receiver_nusselt_v10(0.01, 40.0, 0.7, pnew_v10) >= NU_FD_RECEIVER_v10
    @test rosseland_beta_v10(pnew_v10) > 0.0
    @test rosseland_conductivity_v10(1000.0, pnew_v10) > 0.0
    @test rosseland_optical_thickness_v10(L, pnew_v10) > 0.0
    @test isapprox(CAVITY_OUTER_RADIUS_v10, 0.075)
    @test isapprox(INSULATION_OUTER_RADIUS_v10, 0.057)
    @test isapprox(ADAPTOR_DIAMETER_v10, 0.0776)
    @test RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v10 > 0.0
    @test CAVITY_HEAT_CAPACITY_v10 > 0.0
    @test isapprox(REAR_TUBE_LENGTH_v10, 0.150)
    @test isapprox(REAR_TUBE_CAVITY_LENGTH_v10, 0.046)
    @test isapprox(REAR_TUBE_FLANGE_LENGTH_v10, 0.104)
    @test isapprox(T3_SAMPLE_POSITION_v10, 0.140)
    @test RECEIVER_TO_TUBE_CONDUCTANCE_v10 > 0.0
    @test rear_tube_flange_conductance_per_length_v10(0.100, 500.0) > 0.0
    @test isapprox(sum(solar_weights_v5(pnew_v10[5], pnew_v10[4], 11)), 1.0)
    @test solar_weights_v5(pnew_v10[5], pnew_v10[4], 11)[1] == 1.0

    model, result, experiment = solve_case_v10(pnew_v10, "E74"; nodes=11)
    @test size(model) == (length(result.time), 6)
    @test size(experiment) == (length(result.time), 5)
    @test length(result.boundary_temperature) == length(result.time)
    @test size(result.gas_temperature, 1) == 11 + REAR_TUBE_DEFAULT_NODES_v10 + 1
    @test size(result.rear_tube_temperature) == (REAR_TUBE_DEFAULT_NODES_v10, length(result.time))
    @test result.cavity_temperature[1] == observation(measurements, "E74", "_T2")[1]
    @test all(isfinite, model)
    @test all(isfinite, result.boundary_temperature)
    @test all(isfinite, result.cavity_temperature)
    @test all(isfinite, result.rear_tube_temperature)
    @test all(isfinite, result.flange_heat_loss)
    @test all(isfinite, result.receiver_to_cavity_heat)
    @test all(isfinite, result.tube_to_cavity_heat)
    @test all(isfinite, result.cavity_ambient_heat_loss)
    @test all(model[:, 6] .> 0.0)
    @test successful_retcode(result.ode_solution)

    model_rad, result_rad, experiment_rad = solve_case_v10(
        pnew_v10, "E74"; nodes=11, use_rosseland_radiation=true,
    )
    @test size(model_rad) == size(model)
    @test size(experiment_rad) == size(experiment)
    @test all(isfinite, model_rad)
    @test successful_retcode(result_rad.ode_solution)

    model_profile, result_profile, _ = solve_case_v10(
        pnew_v10, "E74"; nodes=11, use_rosseland_radiation=true,
        rosseland_profile=:front_strong,
    )
    @test size(model_profile) == size(model)
    @test all(isfinite, model_profile)
    @test successful_retcode(result_profile.ode_solution)

    @test isfinite(loss_heating_v10(pnew_v10, ["E74"]; nodes=11))
    @test isfinite(loss_cooling_v10(pnew_v10, ["C69"]; nodes=11))

    trial = optimize_parameter_subset_v10(
        p -> sum(abs2, p[fit_heat_transfer_indices_v10]),
        pnew_v10,
        fit_heat_transfer_indices_v10,
        lb_full_v10,
        ub_full_v10;
        maximum_iterations=5,
        maximum_time_seconds=5.0,
        label="test-v10",
    )
    @test length(trial.parameters) == length(pnew_v10)
end
