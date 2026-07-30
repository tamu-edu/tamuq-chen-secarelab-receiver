# ============================================================================
# Second gate-free felt-range extension at the best structural stress point.
# Cooling side/T2 only; no T3 objective.
# ============================================================================

using Statistics

include("stress_screen_2D_v20_nuisance.jl")

const STRESS_REFINED_FELT_K_2D_v20 = (0.05, 0.15, 0.30, 0.40)
const STRESS_REFINED_FELT_CP_2D_v20 = (0.30, 0.45, 0.55, 0.75)

function stress_refine_2D_v20_felt(
    ; k_indices=1:length(STRESS_REFINED_FELT_K_2D_v20),
    cp_indices=1:length(STRESS_REFINED_FELT_CP_2D_v20),
    max_points=31, suffix="",
)
    exponent = 0.50
    ratio = 2.00
    inputs = [
        case_inputs_2D_v20(
            id; cooling=true, max_points=max_points,
        ) for id in ("C69", "C80")
    ]
    rows = NamedTuple[]
    output = joinpath(
        OUTPUT_DIR_2D_v20,
        "stress_felt_refined$(suffix)_2D_v20.csv",
    )
    for index in k_indices
        felt_k = STRESS_REFINED_FELT_K_2D_v20[index]
        for cp_index in cp_indices
            felt_cp = STRESS_REFINED_FELT_CP_2D_v20[cp_index]
            println(
                "v20 refined felt k=$felt_k cp=$felt_cp",
            )
            flush(stdout)
            p = _stress_parameters_2D_v20(
                exponent, ratio;
                felt_k=felt_k, felt_cp=felt_cp,
                powers=(1.05, 1.23, 0.70),
            )
            cases = [
                simulate_case_2D_v20(
                    input, p; initialization=:side_T2_only,
                    reltol=5e-4, abstol=1e-4, dtmax=120.0,
                ) for input in inputs
            ]
            score = _aggregate_plant_score_2D_v20(
                _stress_cooling_score_2D_v20.(cases),
            )
            push!(rows, merge((
                felt_conductivity_scale=felt_k,
                felt_heat_capacity_scale=felt_cp,
                flow_converged=all(
                    all(case.result.flow_solver_converged)
                    for case in cases
                ),
            ), score))
            _write_namedtuples_csv_2D_v20(output, rows)
        end
    end
    sort!(rows; by=row -> row.objective)
    _write_namedtuples_csv_2D_v20(output, rows)
    println("v20 refined felt best: ", first(rows))
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    first_index = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
    last_index = length(ARGS) >= 2 ? parse(Int, ARGS[2]) :
                 length(STRESS_REFINED_FELT_K_2D_v20)
    points = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 31
    suffix = length(ARGS) >= 4 ? ARGS[4] : ""
    first_cp_index = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 1
    last_cp_index = length(ARGS) >= 6 ? parse(Int, ARGS[6]) :
                    length(STRESS_REFINED_FELT_CP_2D_v20)
    stress_refine_2D_v20_felt(
        ; k_indices=first_index:last_index,
        cp_indices=first_cp_index:last_cp_index,
        max_points=points, suffix=suffix,
    )
end
