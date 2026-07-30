using Test
using Statistics

include(joinpath(
    @__DIR__, "..", "validate_2D_v17_staged.jl",
))

@testset "2D v17 selected screen-to-nominal transfer" begin
    nominal = selected_parameters_2D_v17(
        nominal_mesh=true,
    )
    screen = selected_parameters_2D_v17(
        nominal_mesh=false,
    )
    rows = NamedTuple[]
    for id in ("C69", "C81", "E72")
        cooling = startswith(id, "C")
        inputs = case_inputs_2D_v17(
            id; cooling=cooling, max_points=61,
        )
        coarse = simulate_case_2D_v17(inputs, screen)
        fine = simulate_case_2D_v17(inputs, nominal)
        changes = abs.(
            coarse.model[end, :] .- fine.model[end, :]
        )
        push!(rows, (
            simulation_id=id,
            maximum_final_sensor_change_K=maximum(changes),
            mean_final_sensor_change_K=mean(changes),
        ))
    end
    _write_namedtuples_csv_2D_v17(
        joinpath(
            OUTPUT_DIR_2D_v17,
            "staged_screen_to_nominal_mesh_2D_v17.csv",
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
