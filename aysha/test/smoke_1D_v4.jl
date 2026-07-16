using Test

include(joinpath(@__DIR__, "..", "1D_v4.jl"))

@testset "1D_v4 axial-exchange smoke test" begin
    @test length(pnew_v4) == 15
    @test length(p_cool_init_v4) == 11
    @test length(p_heat_init_v4) == 4

    model, result, experiment = solve_case_v4(pnew_v4, "E74"; nodes=11)
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

    @test isfinite(loss_cooling_v4(p_cool_init_v4, ["C69"]; nodes=11))
    @test isfinite(loss_heating_v4(p_heat_init_v4, ["E74"];
                                   p_cool=p_cool_init_v4, nodes=11))
end
