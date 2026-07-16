using Test

include(joinpath(@__DIR__, "..", "1D_v5.jl"))

@testset "1D_v5 reduced staged model smoke test" begin
    @test length(pnew_v5) == 12
    @test length(p_cool_init_v5) == 8
    @test length(p_heat_init_v5) == 4

    @test axial_exchange_shape_v5(0.0, pnew_v5) ≈ 1.0
    @test axial_exchange_shape_v5(L, pnew_v5) >= pnew_v5[6]
    @test length(solar_weights_v5(pnew_v5[10], pnew_v5[11], 11)) == 11
    @test sum(solar_weights_v5(pnew_v5[10], pnew_v5[11], 11)) ≈ 1.0

    model, result, experiment = solve_case_v5(pnew_v5, "E74"; nodes=11)
    @test size(model) == (length(result.time), 5)
    @test size(experiment) == (length(result.time), 4)
    @test all(isfinite, model)
    @test all(model[:, 5] .> 0.0)
    @test successful_retcode(result.ode_solution)

    paired_T9 = 0.5 .* (observation(measurements, "E74", "_T9") .+
                        observation(measurements, "E74", "_T12"))
    paired_T10 = 0.5 .* (observation(measurements, "E74", "_T10") .+
                         observation(measurements, "E74", "_T11"))
    @test experiment[:, 2] == paired_T9
    @test experiment[:, 3] == paired_T10

    @test isfinite(loss_cooling_v5(pnew_v5, ["C69"]; nodes=11))
    @test isfinite(loss_heating_v5(pnew_v5, ["E74"]; nodes=11))

    trial = optimize_parameter_subset_v5(
        p -> sum(abs2, p[fit_cooling_indices_v5]),
        pnew_v5,
        fit_cooling_indices_v5,
        lb_full_v5,
        ub_full_v5;
        maximum_iterations=5,
        maximum_time_seconds=5.0,
        label="test-v5",
    )
    @test length(trial.parameters) == length(pnew_v5)
end
