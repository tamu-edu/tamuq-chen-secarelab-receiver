# ============================================================================
# Gate-free power-scale profile at the selected v20 structural stress point.
#
# Each irradiance group is independent in the forward model, so its scale is
# profiled on one representative experiment without applying plant gates or
# boundary penalties. T3/T9/T10 are diagnostics only.
# ============================================================================

using Statistics

include("stress_screen_2D_v20_nuisance.jl")

const STRESS_POWER_PROFILES_2D_v20 = (
    (
        name="456", id="E69", index=1,
        levels=(0.70, 0.85, 0.95, 1.05, 1.15, 1.25, 1.40),
    ),
    (
        name="304", id="E74", index=2,
        levels=(0.70, 0.90, 1.05, 1.15, 1.23, 1.35, 1.50),
    ),
    (
        name="256", id="E79", index=3,
        levels=(0.40, 0.55, 0.65, 0.70, 0.75, 0.84, 0.95, 1.10),
    ),
)

function stress_refine_2D_v20_power(; max_points=31)
    mkpath(OUTPUT_DIR_2D_v20)
    jobs = NamedTuple[]
    for group in STRESS_POWER_PROFILES_2D_v20
        for scale in group.levels
            push!(jobs, (
                name=group.name, id=group.id,
                index=group.index, scale=scale,
            ))
        end
    end

    rows = Vector{NamedTuple}(undef, length(jobs))
    Threads.@threads for job_index in eachindex(jobs)
        job = jobs[job_index]
        powers = [1.05, 1.23, 0.70]
        powers[job.index] = job.scale
        p = _stress_parameters_2D_v20(
            0.50, 2.00;
            felt_k=0.40, felt_cp=0.55,
            powers=Tuple(powers),
        )
        input = case_inputs_2D_v20(
            job.id; max_points=max_points,
        )
        case = simulate_case_2D_v20(
            input, p;
            reltol=5e-4, abstol=1e-4, dtmax=120.0,
        )
        score = _plant_case_score_2D_v20(case)
        rows[job_index] = merge((
            power_group=job.name,
            simulation_id=job.id,
            power_scale=job.scale,
            flow_converged=all(
                case.result.flow_solver_converged,
            ),
        ), score)
        println(
            "v20 gate-free power ", job.name,
            " scale=", job.scale, " complete",
        )
        flush(stdout)
    end

    sort!(rows; by=row -> (
        parse(Int, row.power_group), row.power_scale,
    ))
    output = joinpath(
        OUTPUT_DIR_2D_v20,
        "stress_power_refined_2D_v20.csv",
    )
    _write_namedtuples_csv_2D_v20(output, rows)

    selected = NamedTuple[]
    for group in STRESS_POWER_PROFILES_2D_v20
        candidates = filter(
            row -> row.power_group == group.name &&
                   row.flow_converged,
            rows,
        )
        push!(selected, first(sort(
            candidates; by=row -> row.objective,
        )))
    end
    _write_namedtuples_csv_2D_v20(
        joinpath(
            OUTPUT_DIR_2D_v20,
            "stress_power_refined_selected_2D_v20.csv",
        ),
        selected,
    )
    println("v20 gate-free refined powers: ", selected)
    return (rows=rows, selected=selected)
end

if abspath(PROGRAM_FILE) == @__FILE__
    points = isempty(ARGS) ? 31 : parse(Int, ARGS[1])
    stress_refine_2D_v20_power(; max_points=points)
end
