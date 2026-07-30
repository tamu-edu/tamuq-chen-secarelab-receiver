# ============================================================================
# diagnose_2D_v17_flange_boundary.jl - high-G profile after boundary selection
# ============================================================================

include("calibrate_2D_v17_casing_flange.jl")

function diagnose_flange_boundary_2D_v17()
    inputs = [
        case_inputs_2D_v17(id; cooling=true, max_points=61)
        for id in COOLING_TRAIN_2D_v17
    ]
    holdout_inputs = [
        case_inputs_2D_v17(id; cooling=true, max_points=61)
        for id in COOLING_HOLDOUT_2D_v17
    ]
    rows = NamedTuple[]
    for (index, conductance) in enumerate(
        (10.0, 20.0, 50.0, 100.0, 200.0, 500.0),
    )
        candidate = _evaluate_cooling_candidate_2D_v17(
            inputs, conductance, 1.0; nominal_mesh=true,
        )
        holdout = cooling_loss_2D_v17(
            _simulate_cooling_set_2D_v17(
                holdout_inputs, candidate.parameters,
            ),
        )
        push!(rows, (
            evaluation=index,
            casing_flange_conductance_W_K=conductance,
            training_objective=candidate.loss.total,
            training_primary_rmse_K=
                candidate.loss.primary_rmse_K,
            training_t2_rmse_K=candidate.loss.t2_rmse_K,
            holdout_primary_rmse_K=holdout.primary_rmse_K,
            holdout_t2_rmse_K=holdout.t2_rmse_K,
        ))
        println(last(rows))
        flush(stdout)
    end
    _write_namedtuples_csv_2D_v17(
        joinpath(
            OUTPUT_DIR_2D_v17,
            "casing_flange_high_G_profile_2D_v17.csv",
        ), rows,
    )
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    diagnose_flange_boundary_2D_v17()
end
