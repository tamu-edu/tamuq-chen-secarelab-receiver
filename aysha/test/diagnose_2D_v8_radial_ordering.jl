include(joinpath(@__DIR__, "..", "run_2D_v8.jl"))

function read_pilot_theta_radial_2D(path)
    lines = readlines(path)
    return [parse(Float64, split(line, ',')[3]) for line in lines[2:end]]
end

function with_radial_flow_2D(p, c_radial_flow)
    ht = HeatTransferParameters2D(
        coefficient=p.heat_transfer.coefficient,
        reynolds_exponent=p.heat_transfer.reynolds_exponent,
        prandtl_exponent=p.heat_transfer.prandtl_exponent,
        minimum_nusselt=p.heat_transfer.minimum_nusselt,
        c_radial_flow=c_radial_flow,
    )
    return ModelParameters2D(p.geometry, p.solid, ht, p.losses, p.optics)
end

parameter_path = joinpath(
    @__DIR__, "..", "summaries", "2D_v8",
    "parameters_pilot_apparent_nu_2D_v8.csv",
)
p0 = unpack_parameters2D(read_pilot_theta_radial_2D(parameter_path))
for c in (0.10, 0.50, 0.90)
    p = with_radial_flow_2D(p0, c)
    case_data = run_single_case_2D("E81", p; is_cooling=false)
    model = case_data.model[end, :]
    experiment = case_data.experiment[end, :]
    println(
        "c_radial_flow=", c,
        " mid_model=", model[2] - model[4],
        " mid_experiment=", experiment[2] - experiment[4],
        " deep_model=", model[3] - model[5],
        " deep_experiment=", experiment[3] - experiment[5],
        " rmse=", sqrt(mean(abs2, case_data.model .- case_data.experiment)),
    )
end
