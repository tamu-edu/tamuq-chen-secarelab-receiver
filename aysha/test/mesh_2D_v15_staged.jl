using Test
using Statistics

include(joinpath(@__DIR__, "..", "run_2D_v15.jl"))

@testset "2D v15 selected screen-to-nominal transfer" begin
    selected_path = joinpath(
        OUTPUT_DIR_2D_v15,
        "parameters_skin_selected_2D_v15.txt",
    )
    values = Dict{String,Float64}()
    for line in eachline(selected_path)
        occursin("=", line) || continue
        key, raw = split(line, "="; limit=2)
        value = tryparse(Float64, raw)
        value === nothing || (values[key] = value)
    end
    thickness = values["skin_thickness_mm"] * 1e-3
    conductance = values["channel_conductance_scale"]
    nominal = selected_parameters_2D_v15(
        nominal_mesh=true, skin_thickness=thickness,
        skin_conductance_scale=conductance,
    )
    screen = selected_parameters_2D_v15(
        nominal_mesh=false, screen_mesh=true,
        skin_thickness=thickness,
        skin_conductance_scale=conductance,
    )
    rows = NamedTuple[]
    for id in ("E67", "E72", "E77")
        inputs = case_inputs_2D_v15(id; max_points=61)
        coarse = simulate_case_2D_v15(inputs, screen)
        fine = simulate_case_2D_v15(inputs, nominal)
        changes = abs.(coarse.model[end, :] .- fine.model[end, :])
        push!(rows, (
            simulation_id=id,
            maximum_final_sensor_change_K=maximum(changes),
            mean_final_sensor_change_K=mean(changes),
        ))
    end
    _write_namedtuples_csv_2D_v15(
        joinpath(
            OUTPUT_DIR_2D_v15,
            "staged_screen_to_nominal_mesh_2D_v15.csv",
        ), rows,
    )
    foreach(println, rows)
    @test all(
        isfinite(row.maximum_final_sensor_change_K)
        for row in rows
    )
    @test maximum(
        row.maximum_final_sensor_change_K for row in rows
    ) < 20.0
end
