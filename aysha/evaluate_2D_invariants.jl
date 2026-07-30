# ============================================================================
# evaluate_2D_invariants.jl
#
# Julia port of summaries/invariant_evaluator_v0.py.  Applies the manuscript
# steady-state reduction identically to measurements and model predictions.
# Usage: julia --project=. evaluate_2D_invariants.jl 2D_v17
# ============================================================================

using CSV
using DataFrames
using Statistics
using Printf

const WALL_WEIGHTS = (0.251825, 0.350365, 0.397810)
const TTAB = [
    300.0, 400.0, 500.0, 600.0, 700.0,
    800.0, 900.0, 1000.0, 1100.0, 1200.0,
]
const K_AIR = [
    0.0263, 0.0338, 0.0407, 0.0469, 0.0524,
    0.0573, 0.0620, 0.0667, 0.0715, 0.0763,
]
const CP_AIR = [
    1005.0, 1014.0, 1030.0, 1051.0, 1075.0,
    1099.0, 1121.0, 1141.0, 1159.0, 1175.0,
]

function _linear_table_2D_v17(T, values)
    T <= first(TTAB) && return first(values)
    T >= last(TTAB) && return last(values)
    i = searchsortedlast(TTAB, T)
    w = (T - TTAB[i]) / (TTAB[i + 1] - TTAB[i])
    return (1.0 - w) * values[i] + w * values[i + 1]
end

function _invariants_2D_v17(
    T8, T12, T11, T9, T10, T3, Tamb, mdot_gs,
)
    Tw = sum(
        WALL_WEIGHTS[i] * (T8, T12, T11)[i]
        for i in 1:3
    )
    effectiveness = (T3 - Tamb) / (Tw - Tamb)
    Tgm = 0.5 * (Tamb + T3)
    Nu = if effectiveness >= 1.0
        NaN
    else
        NTU = -log(1.0 - effectiveness)
        mch = mdot_gs * 1e-3 / 100.0
        h = NTU * mch *
            _linear_table_2D_v17(Tgm, CP_AIR) /
            (4.0 * 1.5e-3 * 0.137)
        h * 1.5e-3 / _linear_table_2D_v17(
            0.5 * (Tw + Tgm), K_AIR,
        )
    end
    return (
        Tw=Tw, effectiveness=effectiveness, Nu=Nu,
        axial_inversion_K=T12 - T8,
        lambda_58=(T12 - T9) / (T12 - Tamb),
        lambda_107=(T11 - T10) / (T11 - Tamb),
    )
end

function _regression_2D_v17(x, y)
    mx = mean(x)
    my = mean(y)
    slope = sum((x .- mx) .* (y .- my)) /
            sum(abs2, x .- mx)
    intercept = my - slope * mx
    fitted = intercept .+ slope .* x
    r2 = 1.0 - sum(abs2, y .- fitted) /
         sum(abs2, y .- my)
    return (intercept=intercept, slope=slope, r2=r2)
end

function evaluate_invariants_2D(version="2D_v17")
    measured = CSV.read(
        joinpath(
            @__DIR__, "analysis", "exp_analysis",
            "delivered_power_check.csv",
        ), DataFrame,
    )
    metrics = CSV.read(
        joinpath(
            @__DIR__, "summaries", version,
            "staged_sensor_metrics_$(version).csv",
        ), DataFrame,
    )
    metrics = filter(
        row -> startswith(row.phase, "heating"), metrics,
    )
    rows = NamedTuple[]
    for measured_row in eachrow(measured)
        id = measured_row.ID
        case_metrics = filter(
            row -> row.simulation_id == id, metrics,
        )
        nrow(case_metrics) == 7 || continue
        errors = Dict(
            row.sensor => row.steady_error_K
            for row in eachrow(case_metrics)
        )
        measured_inv = _invariants_2D_v17(
            measured_row.T8_ss, measured_row.T12_ss,
            measured_row.T11_ss, measured_row.T9_ss,
            measured_row.T10_ss, measured_row.T3_ss,
            measured_row.Tamb, measured_row.mdot_gs,
        )
        model_inv = _invariants_2D_v17(
            measured_row.T8_ss + errors["T8"],
            measured_row.T12_ss + errors["T12"],
            measured_row.T11_ss + errors["T11"],
            measured_row.T9_ss + errors["T9"],
            measured_row.T10_ss + errors["T10"],
            measured_row.T3_ss + errors["T3"],
            measured_row.Tamb, measured_row.mdot_gs,
        )
        push!(rows, (
            simulation_id=id, Re=measured_row.Re,
            measured_Nu=measured_inv.Nu,
            model_Nu=model_inv.Nu,
            measured_effectiveness=
                measured_inv.effectiveness,
            model_effectiveness=model_inv.effectiveness,
            measured_axial_inversion_K=
                measured_inv.axial_inversion_K,
            model_axial_inversion_K=
                model_inv.axial_inversion_K,
            measured_lambda_58=measured_inv.lambda_58,
            model_lambda_58=model_inv.lambda_58,
            measured_lambda_107=measured_inv.lambda_107,
            model_lambda_107=model_inv.lambda_107,
        ))
    end
    output = DataFrame(rows)
    mkpath(joinpath(@__DIR__, "summaries", version))
    CSV.write(
        joinpath(
            @__DIR__, "summaries", version,
            "manuscript_invariants_$(version).csv",
        ), output,
    )

    summaries = NamedTuple[]
    for source in ("measured", "model")
        nu = getproperty(output, Symbol("$(source)_Nu"))
        finite = isfinite.(nu) .& (nu .> 0.0)
        nu_fit = _regression_2D_v17(
            log.(output.Re[finite]), log.(nu[finite]),
        )
        l58 = getproperty(
            output, Symbol("$(source)_lambda_58"),
        )
        l107 = getproperty(
            output, Symbol("$(source)_lambda_107"),
        )
        fit58 = _regression_2D_v17(output.Re, l58)
        fit107 = _regression_2D_v17(output.Re, l107)
        effectiveness = getproperty(
            output, Symbol("$(source)_effectiveness"),
        )
        inversion = getproperty(
            output, Symbol("$(source)_axial_inversion_K"),
        )
        push!(summaries, (
            source=source,
            Nu_prefactor=exp(nu_fit.intercept),
            Nu_exponent=nu_fit.slope,
            Nu_r2=nu_fit.r2,
            effectiveness_min=minimum(effectiveness),
            effectiveness_max=maximum(effectiveness),
            middle_peak_count=count(>(0.0), inversion),
            lambda_58_intercept=fit58.intercept,
            lambda_58_slope=fit58.slope,
            lambda_58_r2=fit58.r2,
            lambda_58_positive_count=count(>(0.0), l58),
            lambda_107_intercept=fit107.intercept,
            lambda_107_slope=fit107.slope,
            lambda_107_r2=fit107.r2,
            lambda_107_positive_count=count(>(0.0), l107),
        ))
    end
    summary_frame = DataFrame(summaries)
    CSV.write(
        joinpath(
            @__DIR__, "summaries", version,
            "manuscript_invariant_summary_$(version).csv",
        ), summary_frame,
    )
    foreach(println, summaries)
    return (cases=output, summary=summary_frame)
end

if abspath(PROGRAM_FILE) == @__FILE__
    version = isempty(ARGS) ? "2D_v17" : first(ARGS)
    evaluate_invariants_2D(version)
end
