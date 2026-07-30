# ============================================================================
# Gate-free v20 power profile over every heating training case in one
# irradiance group. No gate or boundary penalty enters the ranking.
# ============================================================================

using Statistics

include("stress_screen_2D_v20_nuisance.jl")

const STRESS_TRAIN_POWER_GROUPS_2D_v20 = (
    (
        name="456", index=1,
        train=("E67", "E69", "E71"),
        levels=(0.70, 0.85, 0.95, 1.05, 1.15),
    ),
    (
        name="304", index=2,
        train=("E72", "E74", "E76"),
        levels=(0.70, 0.90, 1.05, 1.15, 1.23),
    ),
    (
        name="256", index=3,
        train=("E77", "E79", "E81"),
        levels=(0.55, 0.65, 0.70, 0.75, 0.84),
    ),
)

function stress_refine_2D_v20_power_group(
    group_index; max_points=31,
)
    group = STRESS_TRAIN_POWER_GROUPS_2D_v20[group_index]
    inputs = [
        case_inputs_2D_v20(id; max_points=max_points)
        for id in group.train
    ]
    rows = NamedTuple[]
    output = joinpath(
        OUTPUT_DIR_2D_v20,
        "stress_power_group_$(group.name)_2D_v20.csv",
    )
    for scale in group.levels
        powers = [1.05, 1.23, 0.75]
        powers[group.index] = scale
        p = _stress_parameters_2D_v20(
            0.50, 2.00;
            felt_k=0.15, felt_cp=0.55,
            powers=Tuple(powers),
        )
        cases = [
            simulate_case_2D_v20(
                input, p;
                reltol=5e-4, abstol=1e-4, dtmax=120.0,
            ) for input in inputs
        ]
        score = _aggregate_plant_score_2D_v20(
            _plant_case_score_2D_v20.(cases),
        )
        push!(rows, merge((
            power_group=group.name,
            power_scale=scale,
            case_count=length(cases),
            flow_converged=all(
                all(case.result.flow_solver_converged)
                for case in cases
            ),
        ), score))
        _write_namedtuples_csv_2D_v20(output, rows)
        println(
            "v20 all-training power ", group.name,
            " scale=", scale, " J=", score.objective,
        )
        flush(stdout)
    end
    sort!(rows; by=row -> row.objective)
    _write_namedtuples_csv_2D_v20(output, rows)
    println("v20 all-training power winner: ", first(rows))
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error(
        "usage: group_index [max_points]",
    )
    group_index = parse(Int, ARGS[1])
    points = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 31
    stress_refine_2D_v20_power_group(
        group_index; max_points=points,
    )
end
