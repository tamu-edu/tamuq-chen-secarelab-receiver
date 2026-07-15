using Test

include(joinpath(@__DIR__, "..", "1D_v1.jl"))
using .Receiver1D

# Reproduces the interactive-session order that previously made Geometry
# ambiguous in Main.
include(joinpath(@__DIR__, "..", "1D_v2.jl"))

@testset "v1 and v2 coexist without exported-name ambiguity" begin
    @test Receiver1D.Geometry !== Receiver1DV2.Geometry
    @test parentmodule(Receiver1D.Geometry) === Receiver1D
    @test parentmodule(Receiver1DV2.Geometry) === Receiver1DV2

    outputs, result = solve_case_v2(pnew, "E74"; nodes=11)
    @test size(outputs, 2) == 5
    @test all(isfinite, outputs)
    @test result.gas_temperature[end, end] > result.gas_temperature[1, end]
end
