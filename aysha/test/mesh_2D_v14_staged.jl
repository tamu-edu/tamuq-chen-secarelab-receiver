using Test
using Statistics

include(joinpath(@__DIR__, "..", "validate_2D_v14_staged.jl"))

@testset "2D v14 staged screen-to-nominal transfer" begin
    nominal = selected_parameters_2D_v14(; nominal_mesh=true)
    screen = selected_parameters_2D_v14(
        nominal_mesh=false, screen_mesh=true,
    )
    rows = NamedTuple[]
    for id in ("E67", "E72", "E77")
        inputs = case_inputs_2D_v14(id; max_points=61)
        coarse = simulate_case_2D_v14(inputs, screen)
        fine = simulate_case_2D_v14(inputs, nominal)
        changes = abs.(coarse.model[end, :] .- fine.model[end, :])
        push!(rows, (
            simulation_id=id,
            maximum_final_sensor_change_K=maximum(changes),
            mean_final_sensor_change_K=mean(changes),
        ))
    end
    _write_namedtuples_csv_2D_v14(
        joinpath(
            OUTPUT_DIR_2D_v14,
            "staged_screen_to_nominal_mesh_2D_v14.csv",
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
