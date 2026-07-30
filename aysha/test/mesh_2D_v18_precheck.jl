using Test
using Statistics

include(joinpath(@__DIR__, "..", "run_2D_v18.jl"))

const PRIMARY_INDICES_2D_v18 = [1, 2, 3, 6, 7]

@testset "2D v18 C/M/F pre-fit mesh check" begin
    specs = (
        (id="E72", cooling=false),
        (id="C69", cooling=true),
        (id="C81", cooling=true),
    )
    rows = NamedTuple[]
    cases = Dict{Tuple{String,Symbol},Any}()
    for spec in specs
        inputs = case_inputs_2D_v18(
            spec.id; cooling=spec.cooling, max_points=61,
        )
        for mesh in (:screen, :nominal, :refined)
            p = parameters_2D_v18(
                mesh=mesh,
                source_model=:beer_lambert,
                exchange_multiplier=1.0,
                power_scales=(1.20, 1.40, 0.975),
            )
            cases[(spec.id, mesh)] = simulate_case_2D_v18(
                inputs, p; reltol=1e-6, abstol=1e-7,
                dtmax=30.0,
            )
        end
        coarse = cases[(spec.id, :screen)]
        medium = cases[(spec.id, :nominal)]
        fine = cases[(spec.id, :refined)]
        for (label, left, right) in (
            ("C_to_M", coarse, medium),
            ("M_to_F", medium, fine),
        )
            delta = left.model[:, PRIMARY_INDICES_2D_v18] .-
                    right.model[:, PRIMARY_INDICES_2D_v18]
            final_delta = abs.(
                left.model[end, PRIMARY_INDICES_2D_v18] .-
                right.model[end, PRIMARY_INDICES_2D_v18]
            )
            push!(rows, (
                simulation_id=spec.id,
                comparison=label,
                history_rms_K=sqrt(mean(abs2, delta)),
                history_max_K=maximum(abs, delta),
                final_rms_K=sqrt(mean(abs2, final_delta)),
                final_max_K=maximum(final_delta),
            ))
        end
    end
    mkpath(OUTPUT_DIR_2D_v18)
    _write_namedtuples_csv_2D_v18(
        joinpath(
            OUTPUT_DIR_2D_v18,
            "prefit_mesh_CMF_2D_v18.csv",
        ), rows,
    )
    foreach(println, rows)
    mf = filter(row -> row.comparison == "M_to_F", rows)
    @test maximum(row.history_rms_K for row in mf) < 10.0
    @test maximum(row.history_max_K for row in mf) < 20.0
end
