# ============================================================================
# extend_2D_v18_full_joint.jl
#
# Boundary-triggered audit of every near/deep combination at weak exchange.
# This prevents the marginal-source pruning from hiding a strong interaction.
# ============================================================================

using Statistics

include("calibrate_2D_v18_staged.jl")

function extend_2D_v18_full_joint()
    inputs = [
        case_inputs_2D_v18(id; max_points=61)
        for id in ("E69", "E74", "E79")
    ]
    candidates = [(
        source_model=:beer_lambert, deep_fraction=0.0,
        deep_length_m=0.08, exchange_multiplier=0.05,
    )]
    append!(candidates, [(
        source_model=:near_deep, deep_fraction=f,
        deep_length_m=L, exchange_multiplier=x,
    ) for x in (0.05, 0.10)
      for f in (0.20, 0.40, 0.60, 0.80)
      for L in (0.02, 0.05, 0.10, 0.20)])
    rows = NamedTuple[]
    results = NamedTuple[]
    for (index, candidate) in enumerate(candidates)
        result = _evaluate_candidate_2D_v18(
            candidate, inputs;
            power_scales=(1.05, 1.23, 0.84),
            felt_conductivity_scale=1.30,
            felt_heat_capacity_scale=1.30,
        )
        model_peaks = count(
            case -> case.model[end, 2] > case.model[end, 1],
            result.cases,
        )
        observed_peaks = count(
            case -> case.observed[end, 2] > case.observed[end, 1],
            result.cases,
        )
        mean_inversion_error = mean(
            (case.model[end, 2] - case.model[end, 1]) -
            (case.observed[end, 2] - case.observed[end, 1])
            for case in result.cases
        )
        value = (
            candidate=candidate, loss=result.loss.total,
            model_peaks=model_peaks,
            mean_inversion_error=mean_inversion_error,
        )
        push!(results, value)
        push!(rows, (
            candidate=_candidate_label_2D_v18(candidate),
            objective=result.loss.total,
            curve=result.loss.curve,
            level=result.loss.level,
            side_shape=result.loss.side_shape,
            effectiveness=result.loss.effectiveness,
            model_middle_peaks=model_peaks,
            observed_middle_peaks=observed_peaks,
            mean_inversion_error_K=mean_inversion_error,
        ))
        println(
            "v18 full joint $index/$(length(candidates)) ",
            _candidate_label_2D_v18(candidate),
            " J=", round(result.loss.total; digits=4),
            " peaks=$model_peaks/$observed_peaks",
        )
        flush(stdout)
    end
    _write_namedtuples_csv_2D_v18(joinpath(
        OUTPUT_DIR_2D_v18, "full_joint_audit_2D_v18.csv",
    ), rows)
    best = first(sort(results; by=x -> x.loss))
    with_peaks = filter(x -> x.model_peaks > 0, results)
    best_shape = isempty(with_peaks) ? nothing : first(sort(
        with_peaks; by=x -> abs(x.mean_inversion_error),
    ))
    println("v18 full-joint best=", best)
    println("v18 full-joint best-with-peak=", best_shape)
    return (best=best, best_shape=best_shape)
end

if abspath(PROGRAM_FILE) == @__FILE__
    extend_2D_v18_full_joint()
end
