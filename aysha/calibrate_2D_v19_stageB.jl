# ============================================================================
# calibrate_2D_v19_stageB.jl
# Stage B: integrated UA from effectiveness and apparent Nu only.
# ============================================================================

using CSV
using DataFrames
using Statistics

include("calibrate_2D_v19_stageA.jl")

const UA_ANCHORS_2D_v19 = ("E72", "E69", "E81")

function _read_keyvalues_2D_v19(path)
    values = Dict{String,String}()
    for line in eachline(path)
        occursin("=", line) || continue
        key, value = split(line, "="; limit=2)
        values[key] = value
    end
    return values
end

function selected_source_2D_v19()
    values = _read_keyvalues_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19, "stageA_selected_2D_v19.txt",
    ))
    return (
        model=Symbol(values["source_model"]),
        fraction=parse(Float64, values["deep_fraction"]),
        length_m=parse(Float64, values["deep_length_m"]),
    )
end

function _linear_air_table_2D_v19(T, values)
    knots = [
        300.0, 400.0, 500.0, 600.0, 700.0,
        800.0, 900.0, 1000.0, 1100.0, 1200.0,
    ]
    T <= first(knots) && return first(values)
    T >= last(knots) && return last(values)
    index = searchsortedlast(knots, T)
    weight = (T - knots[index]) /
             (knots[index + 1] - knots[index])
    return (1.0 - weight) * values[index] +
           weight * values[index + 1]
end

const K_AIR_2D_v19 = [
    0.0263, 0.0338, 0.0407, 0.0469, 0.0524,
    0.0573, 0.0620, 0.0667, 0.0715, 0.0763,
]
const CP_AIR_2D_v19 = [
    1005.0, 1014.0, 1030.0, 1051.0, 1075.0,
    1099.0, 1121.0, 1141.0, 1159.0, 1175.0,
]
const WALL_WEIGHTS_2D_v19 = (0.251825, 0.350365, 0.397810)

function _apparent_invariant_2D_v19(
    side, outlet, ambient, mdot_gs,
)
    wall = sum(
        WALL_WEIGHTS_2D_v19[i] * side[i] for i in 1:3
    )
    effectiveness = (outlet - ambient) /
        max(wall - ambient, 1.0)
    gas_mean = 0.5 * (ambient + outlet)
    nusselt = if !(0.0 < effectiveness < 1.0)
        NaN
    else
        ntu = -log1p(-effectiveness)
        mdot_channel = mdot_gs * 1e-3 / 100.0
        h = ntu * mdot_channel *
            _linear_air_table_2D_v19(
                gas_mean, CP_AIR_2D_v19,
            ) / (4.0 * 1.5e-3 * 0.137)
        h * 1.5e-3 / _linear_air_table_2D_v19(
            0.5 * (wall + gas_mean), K_AIR_2D_v19,
        )
    end
    return (wall=wall, effectiveness=effectiveness, Nu=nusselt)
end

function _steady_window_2D_v19(values)
    n = size(values, 1)
    return max(2, floor(Int, 0.9 * n)):n
end

function ua_case_2D_v19(case, measured_row)
    final = _steady_window_2D_v19(case.model)
    model_side = [
        mean(case.model[final, index]) for index in 1:3
    ]
    receiver_outlet = mean(
        case.result.rear_gas_temperature[1, final],
    )
    model = _apparent_invariant_2D_v19(
        model_side, receiver_outlet,
        Float64(measured_row.Tamb),
        Float64(measured_row.mdot_gs),
    )
    observed = _apparent_invariant_2D_v19(
        (
            Float64(measured_row.T8_ss),
            Float64(measured_row.T12_ss),
            Float64(measured_row.T11_ss),
        ),
        Float64(measured_row.T3_ss),
        Float64(measured_row.Tamb),
        Float64(measured_row.mdot_gs),
    )
    log_error = isfinite(model.Nu) && model.Nu > 0.0 ?
        log(model.Nu / observed.Nu) : 10.0
    objective = 0.5 * rho_2D_v19(
        (model.effectiveness - observed.effectiveness) / 0.04,
    ) + 0.5 * rho_2D_v19(log_error / 0.10)
    return (
        objective=objective,
        model_effectiveness=model.effectiveness,
        observed_effectiveness=observed.effectiveness,
        model_Nu=model.Nu, observed_Nu=observed.Nu,
        log_Nu_error=log_error,
    )
end

function ua_set_2D_v19(cases, measured)
    rows = [
        ua_case_2D_v19(
            case,
            measured[measured.ID .== case.inputs.id, :][1, :],
        ) for case in cases
    ]
    return (
        objective=mean(row.objective for row in rows),
        effectiveness_rmse=sqrt(mean(
            (row.model_effectiveness -
             row.observed_effectiveness)^2 for row in rows
        )),
        log_Nu_rmse=sqrt(mean(
            row.log_Nu_error^2 for row in rows
        )),
        effectiveness_min=minimum(
            row.model_effectiveness for row in rows
        ),
        effectiveness_max=maximum(
            row.model_effectiveness for row in rows
        ),
    )
end

function _ua_parameters_2D_v19(
    source, prefactor, exponent; mesh=:screen,
)
    return parameters_2D_v19(
        mesh=mesh,
        source_model=source.model,
        deep_fraction=source.fraction,
        deep_length_m=source.length_m,
        nu_prefactor=prefactor,
        reynolds_exponent=exponent,
        power_scales=(1.195, 1.405, 0.975),
        felt_conductivity_scale=1.0,
        felt_heat_capacity_scale=1.0,
        rear_tube_flange_contact_h_W_m2_K=0.0,
    )
end

function calibrate_2D_v19_stageB()
    source = selected_source_2D_v19()
    measured = CSV.read(joinpath(
        @__DIR__, "analysis", "exp_analysis",
        "delivered_power_check.csv",
    ), DataFrame)
    anchors = [
        case_inputs_2D_v19(id; max_points=61)
        for id in UA_ANCHORS_2D_v19
    ]
    candidates = [
        (prefactor=A, exponent=n)
        for A in (2.0e-4, 2.5e-4, 3.10e-4, 3.8e-4, 4.5e-4)
        for n in (1.30, 1.445, 1.54, 1.58, 1.65)
    ]
    trace = NamedTuple[]
    evaluated = NamedTuple[]
    for (index, candidate) in enumerate(candidates)
        cases = _simulate_set_2D_v19(
            anchors, _ua_parameters_2D_v19(
                source, candidate.prefactor,
                candidate.exponent,
            ),
        )
        score = ua_set_2D_v19(cases, measured)
        push!(evaluated, (
            candidate=candidate, score=score,
        ))
        push!(trace, merge((
            stage="anchor",
            nu_prefactor=candidate.prefactor,
            reynolds_exponent=candidate.exponent,
        ), score))
        println(
            "v19-B $index/$(length(candidates)) A=",
            candidate.prefactor, " n=", candidate.exponent,
            " J=", round(score.objective; digits=4),
        )
        flush(stdout)
    end
    top = first(sort(
        evaluated; by=x -> x.score.objective,
    ), 2)
    remaining_ids = [
        id for id in TRAIN_HEATING_2D_v19
        if !(id in UA_ANCHORS_2D_v19)
    ]
    remaining = [
        case_inputs_2D_v19(id; max_points=61)
        for id in remaining_ids
    ]
    completed = NamedTuple[]
    for item in top
        score = ua_set_2D_v19(
            _simulate_set_2D_v19(
                remaining, _ua_parameters_2D_v19(
                    source, item.candidate.prefactor,
                    item.candidate.exponent,
                ),
            ), measured,
        )
        combined = (
            objective=(item.score.objective + 2.0 * score.objective) / 3.0,
            effectiveness_rmse=sqrt((
                item.score.effectiveness_rmse^2 +
                2.0 * score.effectiveness_rmse^2
            ) / 3.0),
            log_Nu_rmse=sqrt((
                item.score.log_Nu_rmse^2 +
                2.0 * score.log_Nu_rmse^2
            ) / 3.0),
            effectiveness_min=min(
                item.score.effectiveness_min,
                score.effectiveness_min,
            ),
            effectiveness_max=max(
                item.score.effectiveness_max,
                score.effectiveness_max,
            ),
        )
        push!(completed, (
            candidate=item.candidate, score=combined,
        ))
        push!(trace, merge((
            stage="training_complete",
            nu_prefactor=item.candidate.prefactor,
            reynolds_exponent=item.candidate.exponent,
        ), combined))
    end
    winner = first(sort(
        completed; by=x -> x.score.objective,
    ))
    holdout_inputs = [
        case_inputs_2D_v19(id; max_points=61)
        for id in VALID_HEATING_2D_v19
    ]
    holdout = ua_set_2D_v19(
        _simulate_set_2D_v19(
            holdout_inputs, _ua_parameters_2D_v19(
                source, winner.candidate.prefactor,
                winner.candidate.exponent,
            ),
        ), measured,
    )
    push!(trace, merge((
        stage="heating_holdout",
        nu_prefactor=winner.candidate.prefactor,
        reynolds_exponent=winner.candidate.exponent,
    ), holdout))
    _write_namedtuples_csv_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19, "stageB_ua_trace_2D_v19.csv",
    ), trace)
    pass = 0.8 * 3.07786e-4 <=
            winner.candidate.prefactor <=
            1.2 * 3.07786e-4 &&
        abs(winner.candidate.exponent - 1.44346) <= 0.10 &&
        winner.score.log_Nu_rmse <= 0.10 &&
        holdout.log_Nu_rmse <= 0.10 &&
        winner.score.effectiveness_rmse <= 0.05 &&
        holdout.effectiveness_rmse <= 0.05 &&
        winner.score.effectiveness_max <= 0.85
    open(joinpath(
        OUTPUT_DIR_2D_v19, "stageB_selected_2D_v19.txt",
    ), "w") do io
        println(io, "nu_prefactor=", winner.candidate.prefactor)
        println(io, "reynolds_exponent=",
                winner.candidate.exponent)
        println(io, "training_effectiveness_rmse=",
                winner.score.effectiveness_rmse)
        println(io, "training_log_Nu_rmse=",
                winner.score.log_Nu_rmse)
        println(io, "holdout_effectiveness_rmse=",
                holdout.effectiveness_rmse)
        println(io, "holdout_log_Nu_rmse=",
                holdout.log_Nu_rmse)
        println(io, "stageB_pass=", pass)
    end
    println("v19-B winner=", winner, " holdout=", holdout,
            " pass=", pass)
    return (winner=winner, holdout=holdout, pass=pass)
end

if abspath(PROGRAM_FILE) == @__FILE__
    calibrate_2D_v19_stageB()
end
