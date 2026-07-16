using Test

include(joinpath(@__DIR__, "..", "1D_v6.jl"))

@testset "1D_v6 rear thermal mass smoke test" begin
    @test length(pnew_v6) == 11
    @test H_FRONT_FIXED_V6 == 10.0
    @test length(fit_heat_transfer_indices_v6) == 3
    @test length(fit_cooling_thermal_indices_v6) == 8
    @test sum(solar_weights_v5(BETA_OPT_FIXED_V6, FRONT_DEPOSITION_FIXED_V6, 11)) ≈ 1.0

    model, result, experiment = solve_case_v6(pnew_v6, "E74"; nodes=11)
    @test size(model) == (length(result.time), 5)
    @test size(experiment) == (length(result.time), 4)
    @test length(result.rear_temperature) == length(result.time)
    @test result.rear_temperature[1] == observation(measurements, "E74", "_T2")[1]
    @test all(isfinite, model)
    @test all(isfinite, result.rear_temperature)
    @test all(model[:, 5] .> 0.0)
    @test successful_retcode(result.ode_solution)

    @test isfinite(loss_heating_v6(pnew_v6, ["E74"]; nodes=11))
    @test isfinite(loss_cooling_v6(pnew_v6, ["C69"]; nodes=11))

    trial = optimize_parameter_subset_v6(
        p -> sum(abs2, p[fit_heat_transfer_indices_v6]),
        pnew_v6,
        fit_heat_transfer_indices_v6,
        lb_full_v6,
        ub_full_v6;
        maximum_iterations=5,
        maximum_time_seconds=5.0,
        label="test-v6",
    )
    @test length(trial.parameters) == length(pnew_v6)
end
