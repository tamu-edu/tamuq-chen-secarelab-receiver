using Test
using Statistics

include(joinpath(@__DIR__, "..", "validate_2D_v12_staged.jl"))

@testset "2D v12 staged screen-to-nominal transfer" begin
    nominal = selected_parameters_2D_v12(; nominal_mesh=true)
    screen_base = cooling_selected_parameters_2D_v12(
        nominal_mesh=false, screen_mesh=true,
    )
    screen = rebuild_parameters_2D_v12(
        screen_base;
        beta_opt=150.0,
        heat_scale=1.20,
        power_scales=(1.25, 1.05, 0.77),
    )
    rows = NamedTuple[]
    for id in ("E67", "E72", "E77")
        inputs = case_inputs_2D_v12(id; max_points=61)
        coarse_case = simulate_case_2D_v12(inputs, screen)
        nominal_case = simulate_case_2D_v12(inputs, nominal)
        changes = abs.(coarse_case.model[end, :] .-
                       nominal_case.model[end, :])
        push!(rows, (
            simulation_id=id,
            maximum_final_sensor_change_K=maximum(changes),
            mean_final_sensor_change_K=mean(changes),
        ))
    end
    _write_namedtuples_csv_2D_v12(
        joinpath(
            OUTPUT_DIR_2D_v12,
            "staged_screen_to_nominal_mesh_2D_v12.csv",
        ),
        rows,
    )
    foreach(println, rows)
    @test all(isfinite(row.maximum_final_sensor_change_K) for row in rows)
    # Diagnostic gate.  A failure means the screening mesh cannot validate
    # coefficients even though the nominal confirmation remains usable.
    @test maximum(
        row.maximum_final_sensor_change_K for row in rows
    ) < 20.0
end
