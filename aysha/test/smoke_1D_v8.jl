using Test

include(joinpath(@__DIR__, "..", "1D_v8.jl"))

@testset "1D_v8 rear-tube smoke test" begin
    @test length(pnew_v8) == 10
    @test H_FRONT_FIXED_v8 == 10.0
    @test length(fit_heat_transfer_indices_v8) == 6
    @test length(fit_cooling_thermal_indices_v8) == 4
    @test T2_BOUNDARY_RADIUS_v8 > RECEIVER_EQ_RADIUS_v8
    @test T2_BOUNDARY_RADIUS_v8 < INSULATION_OUTER_RADIUS_v8
    @test RECEIVER_TO_T2_CONDUCTANCE_PER_LENGTH_v8 > 0.0
    @test REAR_TUBE_LENGTH_v8 == 0.200
    @test REAR_TUBE_CAVITY_LENGTH_v8 == 0.100
    @test REAR_TUBE_FLANGE_LENGTH_v8 == 0.100
    @test RECEIVER_TO_TUBE_CONDUCTANCE_v8 > 0.0
    @test rear_tube_flange_conductance_per_length_v8(0.150, 500.0) > 0.0
    @test isapprox(sum(solar_weights_v5(BETA_OPT_FIXED_v8, FRONT_DEPOSITION_FIXED_v8, 11)), 1.0)

    model, result, experiment = solve_case_v8(pnew_v8, "E74"; nodes=11)
    @test size(model) == (length(result.time), 5)
    @test size(experiment) == (length(result.time), 4)
    @test length(result.boundary_temperature) == length(result.time)
    @test size(result.gas_temperature, 1) == 11 + REAR_TUBE_DEFAULT_NODES_v8 + 1
    @test size(result.rear_tube_temperature) == (REAR_TUBE_DEFAULT_NODES_v8, length(result.time))
    @test result.boundary_temperature[1] == observation(measurements, "E74", "_T2")[1]
    @test all(isfinite, model)
    @test all(isfinite, result.boundary_temperature)
    @test all(isfinite, result.rear_tube_temperature)
    @test all(isfinite, result.flange_heat_loss)
    @test all(model[:, 5] .> 0.0)
    @test successful_retcode(result.ode_solution)

    @test isfinite(loss_heating_v8(pnew_v8, ["E74"]; nodes=11))
    @test isfinite(loss_cooling_v8(pnew_v8, ["C69"]; nodes=11))

    trial = optimize_parameter_subset_v8(
        p -> sum(abs2, p[fit_heat_transfer_indices_v8]),
        pnew_v8,
        fit_heat_transfer_indices_v8,
        lb_full_v8,
        ub_full_v8;
        maximum_iterations=5,
        maximum_time_seconds=5.0,
        label="test-v8",
    )
    @test length(trial.parameters) == length(pnew_v8)
end
