include(joinpath(@__DIR__, "..", "run_2D_v8.jl"))

function read_fitted_theta_mesh_2D(path)
    lines = readlines(path)
    return [parse(Float64, split(line, ',')[3]) for line in lines[2:end]]
end

function fine_mesh_parameters_2D(p)
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
        nodes_r_rec=14,
        nodes_r_felt=7,
        nodes_r_case=3,
        nodes_z=45,
        rear_tube_length=g0.rear_tube_length,
        rear_tube_nodes=30,
        rear_tube_inner_radius=g0.rear_tube_inner_radius,
        rear_tube_outer_radius=g0.rear_tube_outer_radius,
        t3_distance_from_receiver=g0.t3_distance_from_receiver,
    )
    return ModelParameters2D(g, p.solid, p.heat_transfer, p.losses, p.optics)
end

function fitted_mesh_response_2D(p, op)
    result = simulate2D(p, op, [0.0, 300.0, 600.0])
    pred = sensor_predictions2D(result)
    return [pred.T8[end], pred.T12[end], pred.T11[end], pred.T9[end],
            pred.T10[end], pred.T3[end], pred.T2[end]]
end

output_dir = joinpath(@__DIR__, "..", "summaries", "2D_v8")
theta = read_fitted_theta_mesh_2D(joinpath(output_dir, "parameters_fitted_2D_v8.csv"))
p = unpack_parameters2D(theta)
op = OperatingCondition2D(
    irradiance=p.optics.scale_304 * 304000.0,
    flow_lpm=10.0,
    inlet_temperature=295.0,
    ambient_temperature=295.0,
)
coarse = fitted_mesh_response_2D(calibration_mesh_parameters_2D(p), op)
nominal = fitted_mesh_response_2D(p, op)
fine = fitted_mesh_response_2D(fine_mesh_parameters_2D(p), op)
coarse_to_nominal = maximum(abs.(nominal .- coarse))
nominal_to_fine = maximum(abs.(fine .- nominal))
accepted = nominal_to_fine < 20.0 && coarse_to_nominal < 20.0

open(joinpath(output_dir, "fitted_mesh_summary_2D_v8.txt"), "w") do io
    println(io, "coarse=$(coarse)")
    println(io, "nominal=$(nominal)")
    println(io, "fine=$(fine)")
    println(io, "coarse_to_nominal_max_K=$(coarse_to_nominal)")
    println(io, "nominal_to_fine_max_K=$(nominal_to_fine)")
    println(io, "acceptance_threshold_K=20.0")
    println(io, "accepted=$(accepted)")
end

println(
    "2D_v8 fitted mesh gate: accepted=", accepted,
    " coarse_to_nominal_max_K=", coarse_to_nominal,
    " nominal_to_fine_max_K=", nominal_to_fine,
)
