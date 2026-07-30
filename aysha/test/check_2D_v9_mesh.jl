using Test

include(joinpath(@__DIR__, "..", "2D_v9.jl"))
using .Receiver2D_v9

function with_axial_nodes_2D_v9(p, nodes_z)
    names = fieldnames(Geometry2D)
    values = Tuple(getfield(p.geometry, name) for name in names)
    geometry_values = NamedTuple{names}(values)
    geometry = Geometry2D(; merge(geometry_values, (nodes_z=nodes_z,))...)
    return ModelParameters2D(
        geometry, p.solid, p.heat_transfer, p.losses, p.optics, p.hydraulics,
    )
end

@testset "2D_v9 hydraulic axial-mesh invariance" begin
    base = default_parameters2D()
    operating = OperatingCondition2D(
        irradiance=0.0,
        flow_lpm=10.0,
        inlet_temperature=500.0,
        ambient_temperature=500.0,
    )

    pressure_drops = Float64[]
    velocities = Float64[]
    for nodes_z in (25, 45, 75)
        p = with_axial_nodes_2D_v9(base, nodes_z)
        result = simulate2D(
            p, operating, [0.0, 1.0e-6]; initial_temperature=500.0,
        )
        push!(pressure_drops, result.receiver_pressure_drop_mbar[1])
        push!(velocities, sum(result.gas_velocity[:, :, 1]) /
                          length(result.gas_velocity[:, :, 1]))
    end

    @test all(isapprox.(pressure_drops, pressure_drops[2]; rtol=1e-12))
    @test all(isapprox.(velocities, velocities[2]; rtol=1e-12))
end
