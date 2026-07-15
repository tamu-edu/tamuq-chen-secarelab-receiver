using Test

include(joinpath(@__DIR__, "..", "1D_v3.jl"))

@testset "1D_v3 finite-volume smoke test" begin
    @test length(pnew) == 11
    @test length(sim_key_heat) == 15
    @test length(sim_key_cool) == 3

    model, result, experiment = solve_case_v3(pnew, "E74"; nodes=11)
    @test size(model) == (length(result.time), 5)
    @test size(experiment) == (length(result.time), 4)
    @test all(isfinite, model)
    @test all(model[:, 5] .> 0.0)
    @test successful_retcode(result.ode_solution)
    @test result.ode_solution.stats.naccept > 0

    cond = simulation_conditions["E74"]
    operating = operating_condition_v3(
        irradiance=cond[Io], flow_lpm=cond[qlpm],
        inlet_temperature=cond[T_in], ambient_temperature=cond[T_amb],
    )
    rates = energy_rates_v3(result.solid_temperature[:, end], result.time[end],
                            pnew, operating)
    @test abs(rates.residual) < 1e-10

    @test isfinite(loss_cooling_v3(p_cool_init, ["C69"]; nodes=11))
    @test isfinite(loss_heating_v3(p_heat_init, ["E74"];
                                   p_cool=p_cool_init, nodes=11))
end
