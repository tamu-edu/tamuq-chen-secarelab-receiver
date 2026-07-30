using Test

include(joinpath(@__DIR__, "..", "2D_v8.jl"))
using .Receiver2D_v8

function with_mesh(p; nr_rec, nr_felt, nr_case, nz, nz_rear)
    g0 = p.geometry
    g = Geometry2D(
        length=g0.length,
        receiver_width=g0.receiver_width,
        receiver_radius=g0.receiver_radius,
        insulation_outer_radius=g0.insulation_outer_radius,
        casing_outer_radius=g0.casing_outer_radius,
        channel_width=g0.channel_width,
        wall_thickness=g0.wall_thickness,
        channel_count=g0.channel_count,
        nodes_r_rec=nr_rec,
        nodes_r_felt=nr_felt,
        nodes_r_case=nr_case,
        nodes_z=nz,
        rear_tube_length=g0.rear_tube_length,
        rear_tube_nodes=nz_rear,
        rear_tube_inner_radius=g0.rear_tube_inner_radius,
        rear_tube_outer_radius=g0.rear_tube_outer_radius,
        t3_distance_from_receiver=g0.t3_distance_from_receiver,
    )
    return ModelParameters2D(g, p.solid, p.heat_transfer, p.losses, p.optics)
end

function final_sensor_vector(p, op)
    result = simulate2D(p, op, [0.0, 300.0, 600.0])
    pred = sensor_predictions2D(result)
    return [pred.T8[end], pred.T12[end], pred.T11[end], pred.T9[end],
            pred.T10[end], pred.T3[end], pred.T2[end]]
end

@testset "2D_v8 mesh-sensitivity diagnostic" begin
    p0 = default_parameters2D()
    coarse = with_mesh(p0; nr_rec=10, nr_felt=5, nr_case=2, nz=30, nz_rear=20)
    nominal = p0
    fine = with_mesh(p0; nr_rec=20, nr_felt=10, nr_case=4, nz=60, nz_rear=40)
    op = OperatingCondition2D(
        irradiance=p0.optics.scale_304 * 304000.0,
        flow_lpm=10.0,
        inlet_temperature=295.0,
        ambient_temperature=295.0,
    )

    y_coarse = final_sensor_vector(coarse, op)
    y_nominal = final_sensor_vector(nominal, op)
    y_fine = final_sensor_vector(fine, op)
    coarse_step = maximum(abs.(y_nominal .- y_coarse))
    fine_step = maximum(abs.(y_fine .- y_nominal))

    @info "2D_v8 mesh comparison at 600 s" y_coarse y_nominal y_fine coarse_step fine_step
    println("coarse=", y_coarse)
    println("nominal=", y_nominal)
    println("fine=", y_fine)
    println("nominal_minus_coarse=", y_nominal .- y_coarse)
    println("fine_minus_nominal=", y_fine .- y_nominal)
    @test all(isfinite, y_coarse)
    @test all(isfinite, y_nominal)
    @test all(isfinite, y_fine)
    @test fine_step < coarse_step
    @test fine_step < 20.0
end
