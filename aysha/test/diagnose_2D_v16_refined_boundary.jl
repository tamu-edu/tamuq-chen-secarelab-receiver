using Statistics

include(joinpath(
    @__DIR__, "..", "refine_2D_v16_cooling_boundary.jl",
))

function _read_refined_2D_v16()
    values = Dict{String,Float64}()
    path = joinpath(
        OUTPUT_DIR_2D_v16,
        "parameters_cooling_refined_2D_v16.txt",
    )
    for line in eachline(path)
        occursin("=", line) || continue
        key, raw = split(line, "="; limit=2)
        value = tryparse(Float64, raw)
        value === nothing || (values[key] = value)
    end
    return values
end

function diagnose_refined_boundary_2D_v16()
    values = _read_refined_2D_v16()
    p = parameters_2D_v16(
        nominal_mesh=false, screen_mesh=true,
        skin_felt_contact_scale=
            values["skin_felt_contact_scale"],
        felt_conductivity_scale=
            values["felt_conductivity_scale"],
        felt_heat_capacity_scale=
            values["felt_heat_capacity_scale"],
    )
    rows = NamedTuple[]
    for id in ("E67", "E72", "E77")
        input = case_inputs_2D_v16(id; max_points=61)
        case = simulate_case_2D_v16(
            input, p;
            reltol=3e-4, abstol=1e-4, dtmax=120.0,
        )
        m = case.model[end, :]
        e = case.observed[end, :]
        push!(rows, (
            simulation_id=id,
            sensor_rmse_K=sqrt(mean(abs2, case.model .- case.observed)),
            model_mid_gap_K=m[2] - m[4],
            observed_mid_gap_K=e[2] - e[4],
            model_deep_gap_K=m[3] - m[5],
            observed_deep_gap_K=e[3] - e[5],
            mid_sign_correct=
                sign(m[2] - m[4]) == sign(e[2] - e[4]),
            deep_sign_correct=
                sign(m[3] - m[5]) == sign(e[3] - e[5]),
        ))
    end
    _write_namedtuples_csv_2D_v16(
        joinpath(
            OUTPUT_DIR_2D_v16,
            "refined_boundary_heating_diagnostic_2D_v16.csv",
        ), rows,
    )
    foreach(println, rows)
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    diagnose_refined_boundary_2D_v16()
end
