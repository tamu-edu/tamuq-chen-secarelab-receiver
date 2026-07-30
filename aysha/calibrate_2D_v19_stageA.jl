# ============================================================================
# calibrate_2D_v19_stageA.jl
# Stage A: centered optical shape from normalized side-wall profiles only.
# ============================================================================

using Statistics

include("run_2D_v19.jl")

rho_2D_v19(x) = 2.0 * (sqrt(1.0 + x * x) - 1.0)

const OPTICAL_ANCHORS_2D_v19 = (
    "E67", "E71", "E72", "E76", "E77", "E81",
)
const OPTICAL_COMPLETION_2D_v19 = ("E69", "E74", "E79")

function _simulate_set_2D_v19(inputs, p; tight=false)
    cases = Vector{Any}(undef, length(inputs))
    Threads.@threads for index in eachindex(inputs)
        cases[index] = simulate_case_2D_v19(
            inputs[index], p;
            reltol=tight ? 2e-5 : 5e-4,
            abstol=tight ? 2e-6 : 1e-4,
            dtmax=tight ? 60.0 : 120.0,
        )
    end
    return cases
end

function side_shape_case_2D_v19(case)
    n = size(case.model, 1)
    final = max(2, floor(Int, 0.9 * n)):n
    ambient = mean(case.inputs.ambient[final])
    model = [
        mean(case.model[final, index]) - ambient
        for index in 1:3
    ]
    observed = [
        mean(case.observed[final, index]) - ambient
        for index in 1:3
    ]
    pmodel = model ./ max(sum(model), 1.0)
    pobserved = observed ./ max(sum(observed), 1.0)
    residual = pmodel .- pobserved
    loss = mean(rho_2D_v19(value / 0.05) for value in residual)
    rmse = sqrt(mean(abs2, residual))
    model_inversion = (model[2] - model[1]) /
        max(mean(model), 1.0)
    observed_inversion = (observed[2] - observed[1]) /
        max(mean(observed), 1.0)
    model_rear_fall = (model[3] - model[2]) /
        max(mean(model), 1.0)
    observed_rear_fall = (observed[3] - observed[2]) /
        max(mean(observed), 1.0)
    return (
        loss=loss, rmse=rmse,
        model_peak=model[2] > model[1],
        observed_peak=observed[2] > observed[1],
        model_inversion=model_inversion,
        observed_inversion=observed_inversion,
        model_rear_fall=model_rear_fall,
        observed_rear_fall=observed_rear_fall,
    )
end

function side_shape_set_2D_v19(cases)
    rows = side_shape_case_2D_v19.(cases)
    true_peaks = count(row -> row.observed_peak, rows)
    recovered = count(
        row -> row.observed_peak && row.model_peak, rows,
    )
    false_peaks = count(
        row -> !row.observed_peak && row.model_peak, rows,
    )
    return (
        objective=mean(row.loss for row in rows),
        normalized_rmse=mean(row.rmse for row in rows),
        observed_peaks=true_peaks,
        recovered_peaks=recovered,
        false_peaks=false_peaks,
        inversion_rmse=sqrt(mean(
            (row.model_inversion - row.observed_inversion)^2
            for row in rows
        )),
        rear_fall_rmse=sqrt(mean(
            (row.model_rear_fall - row.observed_rear_fall)^2
            for row in rows
        )),
    )
end

function _source_label_2D_v19(candidate)
    return candidate.model === :beer_lambert ?
        "beer_lambert" :
        "near_deep_f$(candidate.fraction)_L$(Int(round(
            1e3 * candidate.length_m
        )))mm"
end

function _source_parameters_2D_v19(candidate; mesh=:screen)
    return parameters_2D_v19(
        mesh=mesh,
        source_model=candidate.model,
        deep_fraction=candidate.fraction,
        deep_length_m=candidate.length_m,
        nu_prefactor=3.07786e-4,
        reynolds_exponent=1.44346,
        power_scales=(1.195, 1.405, 0.975),
        felt_conductivity_scale=1.0,
        felt_heat_capacity_scale=1.0,
        rear_tube_flange_contact_h_W_m2_K=0.0,
    )
end

function calibrate_2D_v19_stageA()
    mkpath(OUTPUT_DIR_2D_v19)
    candidates = [(
        model=:beer_lambert, fraction=0.0, length_m=0.08,
    )]
    append!(candidates, [(
        model=:near_deep, fraction=f, length_m=L,
    ) for f in (0.40, 0.60, 0.80, 0.90)
      for L in (0.030, 0.060, 0.120, 0.200)])
    anchor_inputs = [
        case_inputs_2D_v19(id; max_points=61)
        for id in OPTICAL_ANCHORS_2D_v19
    ]
    trace = NamedTuple[]
    evaluated = NamedTuple[]
    for (index, candidate) in enumerate(candidates)
        cases = _simulate_set_2D_v19(
            anchor_inputs, _source_parameters_2D_v19(candidate),
        )
        score = side_shape_set_2D_v19(cases)
        push!(evaluated, (
            candidate=candidate, score=score,
        ))
        push!(trace, merge((
            stage="anchor",
            candidate=_source_label_2D_v19(candidate),
        ), score))
        println(
            "v19-A $index/$(length(candidates)) ",
            _source_label_2D_v19(candidate),
            " J=", round(score.objective; digits=4),
            " peaks=", score.recovered_peaks, "/",
            score.observed_peaks, " false=", score.false_peaks,
        )
        flush(stdout)
    end
    top = first(sort(
        evaluated;
        by=x -> (
            -x.score.recovered_peaks,
            x.score.false_peaks,
            x.score.objective,
        ),
    ), 2)
    completion_inputs = [
        case_inputs_2D_v19(id; max_points=61)
        for id in OPTICAL_COMPLETION_2D_v19
    ]
    completed = NamedTuple[]
    for item in top
        cases = _simulate_set_2D_v19(
            completion_inputs,
            _source_parameters_2D_v19(item.candidate),
        )
        score = side_shape_set_2D_v19(cases)
        combined = (
            objective=(6.0 * item.score.objective +
                       3.0 * score.objective) / 9.0,
            normalized_rmse=(
                6.0 * item.score.normalized_rmse +
                3.0 * score.normalized_rmse
            ) / 9.0,
            observed_peaks=item.score.observed_peaks +
                score.observed_peaks,
            recovered_peaks=item.score.recovered_peaks +
                score.recovered_peaks,
            false_peaks=item.score.false_peaks +
                score.false_peaks,
            inversion_rmse=sqrt((
                6.0 * item.score.inversion_rmse^2 +
                3.0 * score.inversion_rmse^2
            ) / 9.0),
            rear_fall_rmse=sqrt((
                6.0 * item.score.rear_fall_rmse^2 +
                3.0 * score.rear_fall_rmse^2
            ) / 9.0),
        )
        push!(completed, (
            candidate=item.candidate, score=combined,
        ))
        push!(trace, merge((
            stage="training_complete",
            candidate=_source_label_2D_v19(item.candidate),
        ), combined))
    end
    winner = first(sort(
        completed;
        by=x -> (
            -x.score.recovered_peaks,
            x.score.false_peaks,
            x.score.objective,
        ),
    ))
    holdout_inputs = [
        case_inputs_2D_v19(id; max_points=61)
        for id in VALID_HEATING_2D_v19
    ]
    holdout = side_shape_set_2D_v19(_simulate_set_2D_v19(
        holdout_inputs,
        _source_parameters_2D_v19(winner.candidate),
    ))
    push!(trace, merge((
        stage="heating_holdout",
        candidate=_source_label_2D_v19(winner.candidate),
    ), holdout))
    _write_namedtuples_csv_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19, "stageA_optical_trace_2D_v19.csv",
    ), trace)
    combined_recovered =
        winner.score.recovered_peaks + holdout.recovered_peaks
    combined_false =
        winner.score.false_peaks + holdout.false_peaks
    pass = winner.score.normalized_rmse <= 0.05 &&
        winner.score.recovered_peaks >= 5 &&
        winner.score.false_peaks <= 1 &&
        holdout.recovered_peaks >= 3 &&
        holdout.false_peaks <= 1 &&
        combined_recovered >= 8 &&
        combined_false <= 2
    open(joinpath(
        OUTPUT_DIR_2D_v19, "stageA_selected_2D_v19.txt",
    ), "w") do io
        println(io, "source_model=", winner.candidate.model)
        println(io, "deep_fraction=", winner.candidate.fraction)
        println(io, "deep_length_m=", winner.candidate.length_m)
        println(io, "training_normalized_rmse=",
                winner.score.normalized_rmse)
        println(io, "training_recovered_peaks=",
                winner.score.recovered_peaks)
        println(io, "training_false_peaks=",
                winner.score.false_peaks)
        println(io, "holdout_recovered_peaks=",
                holdout.recovered_peaks)
        println(io, "holdout_false_peaks=",
                holdout.false_peaks)
        println(io, "stageA_pass=", pass)
    end
    println("v19-A winner=", winner, " holdout=", holdout,
            " pass=", pass)
    return (winner=winner, holdout=holdout, pass=pass)
end

if abspath(PROGRAM_FILE) == @__FILE__
    calibrate_2D_v19_stageA()
end
