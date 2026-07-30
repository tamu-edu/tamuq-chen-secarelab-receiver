using Statistics

include(joinpath(@__DIR__, "..", "calibrate_2D_v13_staged.jl"))

base = inherited_parameters_2D_v13(
    nominal_mesh=false, screen_mesh=true,
)
heating_inputs = [
    case_inputs_2D_v13(id; max_points=FIT_POINTS_2D_v13)
    for id in TRAIN_HEATING_2D_v13
]
cooling_inputs = [
    case_inputs_2D_v13(
        id; cooling=true, max_points=FIT_POINTS_2D_v13,
    )
    for id in COOLING_IDS_2D_v13
]
candidates = (
    ("collapsed", (0.025, 0.0125, 100.0, 5.0)),
    ("first_selected", (0.075, 0.0225, 30.0, 5.0)),
    ("moderate", (0.30, 0.045, 2.0, 5.0)),
    ("separated", (0.30, 0.0, 0.2, 5.0)),
    ("large_side", (0.45, 0.0675, 2.0, 5.0)),
)
rows = NamedTuple[]
for (label, theta) in candidates
    result = evaluate_population_2D_v13(
        base, theta, heating_inputs, cooling_inputs,
    )
    hc = result.heating_cases
    cc = result.cooling_cases
    mid_sign = count(
        case -> sign(case.model[end, 2] - case.model[end, 4]) ==
                sign(case.observed[end, 2] - case.observed[end, 4]),
        hc,
    )
    deep_sign = count(
        case -> sign(case.model[end, 3] - case.model[end, 5]) ==
                sign(case.observed[end, 3] - case.observed[end, 5]),
        hc,
    )
    push!(rows, (
        label=label,
        objective=result.loss.total,
        heating_rmse_K=mean(
            sqrt(mean(abs2, case.model .- case.observed))
            for case in hc
        ),
        cooling_rmse_K=mean(
            sqrt(mean(abs2, case.model .- case.observed))
            for case in cc
        ),
        mid_sign_correct=mid_sign,
        deep_sign_correct=deep_sign,
        mid_gap_rmse_K=sqrt(mean(
            (
                (case.model[end, 2] - case.model[end, 4]) -
                (case.observed[end, 2] - case.observed[end, 4])
            )^2 for case in hc
        )),
        deep_gap_rmse_K=sqrt(mean(
            (
                (case.model[end, 3] - case.model[end, 5]) -
                (case.observed[end, 3] - case.observed[end, 5])
            )^2 for case in hc
        )),
    ))
end
foreach(println, rows)
_write_namedtuples_csv_2D_v13(
    joinpath(
        OUTPUT_DIR_2D_v13,
        "population_candidate_diagnostics_2D_v13.csv",
    ),
    rows,
)
