# ============================================================================
# smoke_2D_v1.jl - Unit smoke test suite for 2D_v1.jl (SciML ODE Solver)
# ============================================================================

using Test

include(joinpath(@__DIR__, "..", "2D_v1.jl"))
using .Receiver2D

@testset "2D_v1 SciML ODE Solver Smoke Test" begin
    p = default_parameters2D()
    
    @test p.geometry.length == 137e-3
    @test p.geometry.nodes_r_rec == 10
    @test p.geometry.nodes_r_felt == 5
    @test p.geometry.nodes_r_case == 2
    @test p.geometry.nodes_z == 25
    @test p.solid.radial_conductivity_scale == 0.03
    @test p.solid.axial_conductivity_scale == 1.00

    # Material property evaluations
    @test sic_conductivity(900.0) > 0.0
    @test sic_heat_capacity(900.0) > 0.0
    @test air_conductivity(300.0) > 0.0
    @test air_heat_capacity(300.0) > 0.0
    @test air_viscosity(300.0) > 0.0

    # Operating condition & simulation setup
    times = [0.0, 50.0, 100.0]
    op = OperatingCondition2D(irradiance=304000.0, flow_lpm=10.0, inlet_temperature=295.0, ambient_temperature=295.0)

    result = simulate2D(p, op, times)

    @test length(result.time) == 3
    @test size(result.solid_temperature) == (17, 25, 3) # 10 rec + 5 felt + 2 case = 17 total radial nodes
    @test size(result.gas_temperature) == (10, 26, 3)   # 10 receiver channel rings
    @test length(result.cavity_temperature) == 3
    @test all(isfinite, result.solid_temperature)
    @test all(isfinite, result.gas_temperature)
    @test all(isfinite, result.cavity_temperature)

    # Thermocouple prediction extraction
    preds = sensor_predictions2D(result)
    @test length(preds.T8) == 3
    @test length(preds.T12) == 3
    @test length(preds.T11) == 3
    @test length(preds.T9) == 3
    @test length(preds.T10) == 3
    @test length(preds.T3) == 3
    @test length(preds.T2) == 3

    # Check physics: irradiated front core temperature > ambient
    @test result.solid_temperature[1, 1, end] > 295.0
end
