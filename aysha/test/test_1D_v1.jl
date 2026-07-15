using Test

include(joinpath(@__DIR__, "..", "1D_v1.jl"))
using .Receiver1D

@testset "Receiver1D conservative finite-volume model" begin
    ambient = 295.6

    @testset "isothermal equilibrium" begin
        p = ModelParameters(nodes=11)
        op = OperatingCondition(irradiance=0.0, flow_lpm=10.0,
                                inlet_temperature=ambient,
                                ambient_temperature=ambient)
        result = simulate(p, op, [0.0, 20.0, 100.0];
                          initial_temperature=ambient)
        @test maximum(abs.(result.solid_temperature .- ambient)) < 1e-10
        @test maximum(abs.(result.gas_temperature .- ambient)) < 1e-10
    end

    @testset "local gas exchange is physical" begin
        p = ModelParameters(nodes=17)
        op = OperatingCondition(flow_lpm=9.0, inlet_temperature=ambient,
                                ambient_temperature=ambient)
        profile = collect(range(900.0, 650.0, length=p.nodes))
        gas = gas_profile(profile, 0.0, p, op)
        @test all(gas.h .> 0.0)
        @test all(gas.heat_to_gas .> 0.0)
        @test gas.temperature[1] == ambient
        @test ambient < gas.temperature[end] < maximum(profile)
        @test all(diff(gas.temperature) .> 0.0)
    end

    @testset "power balance" begin
        p = ModelParameters(nodes=19)
        op = OperatingCondition(irradiance=304e3, flow_lpm=9.0,
                                inlet_temperature=ambient,
                                ambient_temperature=ambient)
        temperature = collect(range(800.0, 550.0, length=p.nodes))
        rates = energy_rates(temperature, 0.0, p, op)
        @test abs(rates.residual) < 1e-10 * max(rates.absorbed, 1.0)
        @test rates.absorbed > 0.0
        @test rates.gas > 0.0
    end

    @testset "heating produces axial and gas profiles" begin
        p = ModelParameters(nodes=15)
        op = OperatingCondition(irradiance=256e3, flow_lpm=6.6,
                                inlet_temperature=ambient,
                                ambient_temperature=ambient)
        result = simulate(p, op, 0.0:30.0:300.0;
                          initial_temperature=ambient)
        @test result.solid_temperature[1, end] > result.solid_temperature[end, end]
        @test result.gas_temperature[end, end] > ambient
        @test heat_transfer_summary(result).mean_h > 0.0
        @test length(solid_at(result, 0.058)) == length(result.time)
    end
end
