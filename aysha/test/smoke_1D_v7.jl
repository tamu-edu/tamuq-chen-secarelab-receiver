using Test

include(joinpath(@__DIR__, "..", "1D_v7.jl"))

@testset "1D_v7 T2-boundary smoke test" begin
    @test length(pnew_v7) == 10
    @test H_FRONT_FIXED_V7 == 10.0
    @test length(fit_heat_transfer_indices_v7) == 6
    @test length(fit_cooling_thermal_indices_v7) == 4
    @test T2_BOUNDARY_RADIUS_V7 > RECEIVER_EQ_RADIUS_V7
    @test T2_BOUNDARY_RADIUS_V7 < INSULATION_OUTER_RADIUS_V7
    @test RECEIVER_TO_T2_CONDUCTANCE_PER_LENGTH_V7 > 0.0
    @test RECEIVER_REAR_TO_T2_CONDUCTANCE_V7 > 0.0
    @test isapprox(sum(solar_weights_v5(BETA_OPT_FIXED_V7, FRONT_DEPOSITION_FIXED_V7, 11)), 1.0)

    model, result, experiment = solve_case_v7(pnew_v7, "E74"; nodes=11)
    @test size(model) == (length(result.time), 5)
    @test size(experiment) == (length(result.time), 4)
    @test length(result.boundary_temperature) == length(result.time)
    @test result.boundary_temperature[1] == observation(measurements, "E74", "_T2")[1]
    @test all(isfinite, model)
    @test all(isfinite, result.boundary_temperature)
    @test all(model[:, 5] .> 0.0)
    @test successful_retcode(result.ode_solution)

    @test isfinite(loss_heating_v7(pnew_v7, ["E74"]; nodes=11))
    @test isfinite(loss_cooling_v7(pnew_v7, ["C69"]; nodes=11))

    trial = optimize_parameter_subset_v7(
        p -> sum(abs2, p[fit_heat_transfer_indices_v7]),
        pnew_v7,
        fit_heat_transfer_indices_v7,
        lb_full_v7,
        ub_full_v7;
        maximum_iterations=5,
        maximum_time_seconds=5.0,
        label="test-v7",
    )
    @test length(trial.parameters) == length(pnew_v7)
end
