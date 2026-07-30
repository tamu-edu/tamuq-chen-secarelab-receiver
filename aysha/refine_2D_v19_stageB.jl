# ============================================================================
# refine_2D_v19_stageB.jl
# Minimal constrained refinement at the Stage-B acceptance corner.
# This classifies the coarse miss; it does not reopen or widen the bounds.
# ============================================================================

include("calibrate_2D_v19_stageB.jl")

const REFINED_UA_CANDIDATES_2D_v19 = (
    (prefactor=3.693432e-4, exponent=1.54346),
    (prefactor=3.650000e-4, exponent=1.54346),
    (prefactor=3.693432e-4, exponent=1.53000),
    (prefactor=3.650000e-4, exponent=1.53000),
)

function refine_2D_v19_stageB()
    source = selected_source_2D_v19()
    measured = CSV.read(joinpath(
        @__DIR__, "analysis", "exp_analysis",
        "delivered_power_check.csv",
    ), DataFrame)
    anchor_inputs = [
        case_inputs_2D_v19(id; max_points=61)
        for id in UA_ANCHORS_2D_v19
    ]
    rows = NamedTuple[]
    anchor_results = NamedTuple[]
    for candidate in REFINED_UA_CANDIDATES_2D_v19
        cases = _simulate_set_2D_v19(
            anchor_inputs,
            _ua_parameters_2D_v19(
                source, candidate.prefactor,
                candidate.exponent,
            ),
        )
        score = ua_set_2D_v19(cases, measured)
        push!(anchor_results, (
            candidate=candidate, score=score,
        ))
        push!(rows, merge((
            stage="anchor_C",
            nu_prefactor=candidate.prefactor,
            reynolds_exponent=candidate.exponent,
        ), score))
        println("v19-B-refine anchor ", candidate,
                " J=", round(score.objective; digits=4))
        flush(stdout)
    end

    remaining_ids = [
        id for id in TRAIN_HEATING_2D_v19
        if !(id in UA_ANCHORS_2D_v19)
    ]
    remaining_inputs = [
        case_inputs_2D_v19(id; max_points=61)
        for id in remaining_ids
    ]
    completed = NamedTuple[]
    for item in first(sort(
        anchor_results; by=x -> x.score.objective,
    ), 2)
        remaining_score = ua_set_2D_v19(
            _simulate_set_2D_v19(
                remaining_inputs,
                _ua_parameters_2D_v19(
                    source, item.candidate.prefactor,
                    item.candidate.exponent,
                ),
            ), measured,
        )
        combined = (
            objective=(item.score.objective +
                       2.0 * remaining_score.objective) / 3.0,
            effectiveness_rmse=sqrt((
                item.score.effectiveness_rmse^2 +
                2.0 * remaining_score.effectiveness_rmse^2
            ) / 3.0),
            log_Nu_rmse=sqrt((
                item.score.log_Nu_rmse^2 +
                2.0 * remaining_score.log_Nu_rmse^2
            ) / 3.0),
            effectiveness_min=min(
                item.score.effectiveness_min,
                remaining_score.effectiveness_min,
            ),
            effectiveness_max=max(
                item.score.effectiveness_max,
                remaining_score.effectiveness_max,
            ),
        )
        result = (candidate=item.candidate, score=combined)
        push!(completed, result)
        push!(rows, merge((
            stage="training_complete_C",
            nu_prefactor=item.candidate.prefactor,
            reynolds_exponent=item.candidate.exponent,
        ), combined))
    end
    winner = first(sort(completed; by=x -> x.score.objective))
    train_signature_pass =
        winner.score.log_Nu_rmse <= 0.10 &&
        winner.score.effectiveness_rmse <= 0.05 &&
        winner.score.effectiveness_max <= 0.85

    holdout = nothing
    holdout_pass = false
    nominal = nothing
    nominal_pass = false
    if train_signature_pass
        holdout_inputs = [
            case_inputs_2D_v19(id; max_points=61)
            for id in VALID_HEATING_2D_v19
        ]
        holdout = ua_set_2D_v19(
            _simulate_set_2D_v19(
                holdout_inputs,
                _ua_parameters_2D_v19(
                    source, winner.candidate.prefactor,
                    winner.candidate.exponent,
                ),
            ), measured,
        )
        holdout_pass =
            holdout.log_Nu_rmse <= 0.10 &&
            holdout.effectiveness_rmse <= 0.05 &&
            holdout.effectiveness_max <= 0.85
        push!(rows, merge((
            stage="holdout_C",
            nu_prefactor=winner.candidate.prefactor,
            reynolds_exponent=winner.candidate.exponent,
        ), holdout))
    end
    if holdout_pass
        nominal = ua_set_2D_v19(
            _simulate_set_2D_v19(
                anchor_inputs,
                _ua_parameters_2D_v19(
                    source, winner.candidate.prefactor,
                    winner.candidate.exponent;
                    mesh=:nominal,
                ),
            ), measured,
        )
        nominal_pass =
            nominal.log_Nu_rmse <= 0.10 &&
            nominal.effectiveness_rmse <= 0.05 &&
            nominal.effectiveness_max <= 0.85
        push!(rows, merge((
            stage="anchor_M",
            nu_prefactor=winner.candidate.prefactor,
            reynolds_exponent=winner.candidate.exponent,
        ), nominal))
    end

    upper_A = 1.2 * 3.07786e-4
    upper_n = 1.44346 + 0.10
    interior = winner.candidate.prefactor < upper_A &&
               winner.candidate.exponent < upper_n
    stage_pass = train_signature_pass && holdout_pass &&
                 nominal_pass && interior
    _write_namedtuples_csv_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19,
        "stageB_refined_trace_2D_v19.csv",
    ), rows)
    open(joinpath(
        OUTPUT_DIR_2D_v19,
        "stageB_refined_selected_2D_v19.txt",
    ), "w") do io
        println(io, "nu_prefactor=", winner.candidate.prefactor)
        println(io, "reynolds_exponent=",
                winner.candidate.exponent)
        println(io, "training_effectiveness_rmse=",
                winner.score.effectiveness_rmse)
        println(io, "training_log_Nu_rmse=",
                winner.score.log_Nu_rmse)
        println(io, "training_signature_pass=",
                train_signature_pass)
        println(io, "holdout_run=", holdout !== nothing)
        holdout !== nothing && println(
            io, "holdout_effectiveness_rmse=",
            holdout.effectiveness_rmse,
        )
        holdout !== nothing && println(
            io, "holdout_log_Nu_rmse=",
            holdout.log_Nu_rmse,
        )
        println(io, "nominal_confirmation_run=",
                nominal !== nothing)
        println(io, "parameter_interior=", interior)
        println(io, "stageB_refined_pass=", stage_pass)
    end
    println(
        "v19-B-refine winner=", winner,
        " holdout=", holdout,
        " nominal=", nominal,
        " interior=", interior,
        " pass=", stage_pass,
    )
    return (
        winner=winner, holdout=holdout, nominal=nominal,
        interior=interior, pass=stage_pass,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    refine_2D_v19_stageB()
end
