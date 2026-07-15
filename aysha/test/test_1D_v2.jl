using Test

include(joinpath(@__DIR__, "..", "1D_v2.jl"))

@testset "1D_v2 legacy-style research flow" begin
    @test length(sim_key_heat) == 15
    @test length(sim_key_cool) == 3
    @test length(measurements) == 15 * 13
    @test length(measurements_cooling) == 3 * 13

    @testset "external import contract" begin
        @test haskey(simulation_conditions, "E67")
        @test haskey(simulation_conditions_cooling, "C69")
        @test 15.2 < simulation_conditions["E67"][qlpm] < 15.4
        @test simulation_conditions["E67"][Io] == 456000.0
        @test simulation_conditions["E72"][Io] == 304000.0
        @test simulation_conditions["E77"][Io] == 256000.0
        @test simulation_conditions_cooling["C69"][Io] == 0.0
        @test length(observation(measurements, "E67", "_T8")) ==
              length(observation_time(measurements, "E67"))
    end

    @testset "model and remake interfaces" begin
        case = "E74"
        outputs, result = solve_case_v2(pnew, case; nodes=11)
        @test size(outputs, 1) == length(observation_time(measurements, case))
        @test size(outputs, 2) == length(MODEL_OUTPUT_COLUMNS)
        @test all(isfinite, outputs)
        @test all(outputs[:, 5] .> 0.0)
        @test result.gas_temperature[end, end] > result.gas_temperature[1, end]

        cond = simulation_conditions[case]
        time = observation_time(measurements, case)
        remade = remakeAysha(pnew, cond, time, cond[Tinit]; nodes=11)
        @test size(remade) == size(outputs)
        @test all(isfinite, remade)
    end

    @testset "objectives and post-processing" begin
        @test isfinite(loss_cooling(p_cool_init, ["C69"]; nodes=11))
        @test isfinite(loss_heating(p_heat_init, ["E74"];
                                    p_cool=p_cool_init, nodes=11))
        steady = build_steady_results_v2(pnew; keys=["E74"], nodes=11)
        @test length(steady) == 1
        @test steady[1].h_effective > 0.0
        metrics = compute_metrics_v2(pnew; heating_keys=["E74"],
                                     cooling_keys=String[], nodes=11)
        @test length(metrics) == 4
        @test all(row -> row.rmse_K >= 0.0, metrics)
    end
end
