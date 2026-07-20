using Test

include(joinpath(@__DIR__, "..", "1D_v8b.jl"))

@testset "1D_v8b predicted-cavity smoke test" begin
    @test length(pnew_v8b) == 10
    @test H_FRONT_FIXED_v8b == 10.0
    @test length(fit_heat_transfer_indices_v8b) == 6
    @test length(fit_cooling_thermal_indices_v8b) == 3
    @test isapprox(CAVITY_OUTER_RADIUS_v8b, 0.075)
    @test isapprox(INSULATION_OUTER_RADIUS_v8b, 0.057)
    @test isapprox(ADAPTOR_DIAMETER_v8b, 0.0776)
    @test RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v8b > 0.0
    @test CAVITY_HEAT_CAPACITY_v8b > 0.0
    @test isapprox(REAR_TUBE_LENGTH_v8b, 0.150)
    @test isapprox(REAR_TUBE_CAVITY_LENGTH_v8b, 0.046)
    @test isapprox(REAR_TUBE_FLANGE_LENGTH_v8b, 0.104)
    @test isapprox(T3_SAMPLE_POSITION_v8b, 0.140)
    @test RECEIVER_TO_TUBE_CONDUCTANCE_v8b > 0.0
    @test rear_tube_flange_conductance_per_length_v8b(0.100, 500.0) > 0.0
    @test isapprox(sum(solar_weights_v5(BETA_OPT_FIXED_v8b, FRONT_DEPOSITION_FIXED_v8b, 11)), 1.0)

    model, result, experiment = solve_case_v8b(pnew_v8b, "E74"; nodes=11)
    @test size(model) == (length(result.time), 6)
    @test size(experiment) == (length(result.time), 5)
    @test length(result.boundary_temperature) == length(result.time)
    @test size(result.gas_temperature, 1) == 11 + REAR_TUBE_DEFAULT_NODES_v8b + 1
    @test size(result.rear_tube_temperature) == (REAR_TUBE_DEFAULT_NODES_v8b, length(result.time))
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

    @test isfinite(loss_heating_v8b(pnew_v8b, ["E74"]; nodes=11))
    @test isfinite(loss_cooling_v8b(pnew_v8b, ["C69"]; nodes=11))

    trial = optimize_parameter_subset_v8b(
        p -> sum(abs2, p[fit_heat_transfer_indices_v8b]),
        pnew_v8b,
        fit_heat_transfer_indices_v8b,
        lb_full_v8b,
        ub_full_v8b;
        maximum_iterations=5,
        maximum_time_seconds=5.0,
        label="test-v8b",
    )
    @test length(trial.parameters) == length(pnew_v8b)
end
