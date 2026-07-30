# ============================================================================
# diagnose_2D_v18_peak_candidate.jl
#
# Generalization check for the full-joint candidate that reproduces the
# representative middle peaks, irrespective of its worse composite objective.
# ============================================================================

using Statistics

include("calibrate_2D_v18_staged.jl")

function diagnose_2D_v18_peak_candidate()
    candidate = (
        source_model=:near_deep, deep_fraction=0.60,
        deep_length_m=0.050, exchange_multiplier=0.10,
    )
    p = _parameter_tuple_2D_v18(
        candidate; mesh=:screen,
        power_scales=(1.05, 1.23, 0.84),
        felt_conductivity_scale=1.30,
        felt_heat_capacity_scale=1.30,
    )
    inputs = [
        case_inputs_2D_v18(id; max_points=61)
        for id in HEATING_IDS_2D_v18
    ]
    cases = _simulate_set_2D_v18(inputs, p)
    rows = NamedTuple[]
    for case in cases
        model_delta = case.model[end, 2] - case.model[end, 1]
        observed_delta = case.observed[end, 2] -
                         case.observed[end, 1]
        push!(rows, (
            simulation_id=case.inputs.id,
            phase=case.inputs.id in TRAIN_HEATING_2D_v18 ?
                "heating_training" : "heating_validation",
            model_inversion_K=model_delta,
            observed_inversion_K=observed_delta,
            model_middle_peak=model_delta > 0.0,
            observed_middle_peak=observed_delta > 0.0,
            inversion_error_K=model_delta - observed_delta,
            side_rmse_K=sqrt(mean(abs2,
                case.model[:, 1:3] .- case.observed[:, 1:3],
            )),
            air_rmse_K=sqrt(mean(abs2,
                case.model[:, 6] .- case.observed[:, 6],
            )),
            felt_rmse_K=sqrt(mean(abs2,
                case.model[:, 7] .- case.observed[:, 7],
            )),
            primary_rmse_K=sqrt(mean(abs2,
                case.model[:, [1, 2, 3, 6, 7]] .-
                case.observed[:, [1, 2, 3, 6, 7]],
            )),
        ))
        println("v18 peak diagnostic ", case.inputs.id, " complete")
        flush(stdout)
    end
    _write_namedtuples_csv_2D_v18(joinpath(
        OUTPUT_DIR_2D_v18,
        "peak_candidate_screen_2D_v18.csv",
    ), rows)
    for phase in ("heating_training", "heating_validation")
        selected = filter(row -> row.phase == phase, rows)
        observed_peaks = count(row -> row.observed_middle_peak, selected)
        correct_peaks = count(
            row -> row.observed_middle_peak && row.model_middle_peak,
            selected,
        )
        false_peaks = count(
            row -> !row.observed_middle_peak && row.model_middle_peak,
            selected,
        )
        println((
            phase=phase,
            primary_rmse_K=mean(row.primary_rmse_K for row in selected),
            side_rmse_K=mean(row.side_rmse_K for row in selected),
            air_rmse_K=mean(row.air_rmse_K for row in selected),
            felt_rmse_K=mean(row.felt_rmse_K for row in selected),
            observed_peaks=observed_peaks,
            correct_peaks=correct_peaks,
            false_peaks=false_peaks,
            axial_rmse_K=sqrt(mean(
                row.inversion_error_K^2 for row in selected
            )),
        ))
    end
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    diagnose_2D_v18_peak_candidate()
end
