using Statistics

include(joinpath(@__DIR__, "..", "calibrate_2D_v14_staged.jl"))

base = inherited_parameters_2D_v14(
    nominal_mesh=false, screen_mesh=true,
)
heating_inputs = [
    case_inputs_2D_v14(id; max_points=FIT_POINTS_2D_v14)
    for id in TRAIN_HEATING_2D_v14
]
cooling_inputs = [
    case_inputs_2D_v14(
        id; cooling=true, max_points=FIT_POINTS_2D_v14,
    )
    for id in COOLING_IDS_2D_v14
]
candidates = (
    ("selected", 1.0, 3.0, 100e-3),
    ("nominal", 1.0, 1.0, 14e-3),
    ("moderate_isolation", 0.20, 0.30, 100e-3),
    ("strong_isolation", 0.05, 0.30, 30e-3),
)
rows = NamedTuple[]
for (label, lateral, felt, sigma) in candidates
    network = rebuild_network_2D_v14(
        base; lateral_scale=lateral, edge_felt_scale=felt,
    )
    p = rebuild_heating_2D_v14(
        network; beam_sigma=sigma,
    )
    heating = _simulate_set_2D_v14(heating_inputs, p)
    cooling = _simulate_set_2D_v14(cooling_inputs, p)
    loss = network_loss_2D_v14(heating, cooling)
    push!(rows, (
        label=label,
        objective=loss.total,
        heating_rmse_K=mean(
            sqrt(mean(abs2, case.model .- case.observed))
            for case in heating
        ),
        cooling_rmse_K=mean(
            sqrt(mean(abs2, case.model .- case.observed))
            for case in cooling
        ),
        mid_sign_correct=count(
            case -> sign(
                case.model[end, 2] - case.model[end, 4]
            ) == sign(
                case.observed[end, 2] - case.observed[end, 4]
            ), heating,
        ),
        deep_sign_correct=count(
            case -> sign(
                case.model[end, 3] - case.model[end, 5]
            ) == sign(
                case.observed[end, 3] - case.observed[end, 5]
            ), heating,
        ),
        mid_gap_rmse_K=sqrt(mean(
            (
                (case.model[end, 2] - case.model[end, 4]) -
                (case.observed[end, 2] - case.observed[end, 4])
            )^2 for case in heating
        )),
        deep_gap_rmse_K=sqrt(mean(
            (
                (case.model[end, 3] - case.model[end, 5]) -
                (case.observed[end, 3] - case.observed[end, 5])
            )^2 for case in heating
        )),
    ))
end
foreach(println, rows)
_write_namedtuples_csv_2D_v14(
    joinpath(
        OUTPUT_DIR_2D_v14,
        "network_candidate_diagnostics_2D_v14.csv",
    ), rows,
)
