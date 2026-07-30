# ============================================================================
# calibrate_2D_v19_stageD.jl
# Stage D: bounded felt-property screen followed by groupwise power fitting.
#
# This stage is intentionally downstream of the structural gates.  It cannot
# rescue a rejected optical, integrated-UA, or outlet-observation model.
# T9/T10 and the core-to-gas balance are not used in the fitted objective.
# ============================================================================

using Statistics

include("calibrate_2D_v19_stageC.jl")

const POWER_GROUPS_2D_v19 = (
    (name="456", train=("E67", "E69", "E71"),
     valid=("E68", "E70"),
     values=(1.05, 1.1225, 1.195, 1.2675, 1.34)),
    (name="304", train=("E72", "E74", "E76"),
     valid=("E73", "E75"),
     values=(1.23, 1.3175, 1.405, 1.4925, 1.58)),
    (name="256", train=("E77", "E79", "E81"),
     valid=("E78", "E80"),
     values=(0.84, 0.9075, 0.975, 1.0425, 1.11)),
)

function selected_stageC_2D_v19()
    values = _read_keyvalues_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19, "stageC_selected_2D_v19.txt",
    ))
    return (
        tube_h=parse(
            Float64,
            values["rear_tube_flange_contact_h_W_m2_K"],
        ),
        probe_capacity=parse(
            Float64,
            values["probe_capacity_areal_J_m2_K"],
        ),
        probe_stem=parse(
            Float64,
            values["probe_stem_conductance_areal_W_m2_K"],
        ),
    )
end

function _stageD_parameters_2D_v19(
    source, ua, outlet;
    felt_k=1.0, felt_cp=1.0,
    power_scales=(1.195, 1.405, 0.975),
    mesh=:screen,
)
    return parameters_2D_v19(
        mesh=mesh,
        source_model=source.model,
        deep_fraction=source.fraction,
        deep_length_m=source.length_m,
        nu_prefactor=ua.prefactor,
        reynolds_exponent=ua.exponent,
        probe_capacity_areal_J_m2_K=outlet.probe_capacity,
        probe_stem_conductance_areal_W_m2_K=
            outlet.probe_stem,
        rear_tube_flange_contact_h_W_m2_K=outlet.tube_h,
        felt_conductivity_scale=felt_k,
        felt_heat_capacity_scale=felt_cp,
        power_scales=power_scales,
    )
end

function fit_case_score_2D_v19(case; mode=:heating)
    rows = 2:size(case.model, 1)
    final = _steady_window_2D_v19(case.model)
    side_errors = case.model[rows, 1:3] .-
                  case.observed[rows, 1:3]
    air_errors = case.model[rows, 6] .-
                 case.observed[rows, 6]
    felt_errors = case.model[rows, 7] .-
                  case.observed[rows, 7]
    side_curve = mean(rho_2D_v19(
        value / 35.0
    ) for value in side_errors)
    air_curve = mean(
        rho_2D_v19(value / 25.0) for value in air_errors
    )
    felt_curve = mean(
        rho_2D_v19(value / 15.0) for value in felt_errors
    )
    side_level = mean(
        rho_2D_v19((
            mean(case.model[final, index]) -
            mean(case.observed[final, index])
        ) / 30.0) for index in 1:3
    )
    air_level = rho_2D_v19((
        mean(case.model[final, 6]) -
        mean(case.observed[final, 6])
    ) / 20.0)
    felt_level = rho_2D_v19((
        mean(case.model[final, 7]) -
        mean(case.observed[final, 7])
    ) / 10.0)
    if mode === :cooling
        curve = 0.40 * side_curve + 0.30 * air_curve +
                0.30 * felt_curve
        level = 0.40 * side_level + 0.30 * air_level +
                0.30 * felt_level
        objective = 0.70 * curve + 0.30 * level
    elseif mode === :heating
        curve = 0.50 * side_curve + 0.30 * air_curve +
                0.20 * felt_curve
        level = 0.45 * side_level + 0.35 * air_level +
                0.20 * felt_level
        objective = 0.40 * curve + 0.60 * level
    else
        error("unknown fit-score mode: $mode")
    end
    return (
        objective=objective,
        side_rmse_K=sqrt(mean(abs2, side_errors)),
        air_rmse_K=sqrt(mean(abs2, air_errors)),
        felt_rmse_K=sqrt(mean(abs2, felt_errors)),
        side_final_bias_K=mean(
            case.model[final, 1:3] .-
            case.observed[final, 1:3],
        ),
        air_final_bias_K=mean(
            case.model[final, 6] .-
            case.observed[final, 6],
        ),
        felt_final_bias_K=mean(
            case.model[final, 7] .-
            case.observed[final, 7],
        ),
    )
end

function fit_set_score_2D_v19(cases; mode=:heating)
    scores = [
        fit_case_score_2D_v19(case; mode=mode)
        for case in cases
    ]
    return NamedTuple(
        name => mean(getproperty(score, name) for score in scores)
        for name in propertynames(first(scores))
    )
end

function _fit_felt_2D_v19(source, ua, outlet)
    training_inputs = [
        case_inputs_2D_v19(
            id; cooling=true, max_points=61,
        ) for id in ("C69", "C80")
    ]
    conductivities = (0.70, 1.00, 1.30)
    capacities = (0.75, 1.00, 1.50)
    rows = NamedTuple[]
    candidates = NamedTuple[]
    for conductivity in conductivities, capacity in capacities
        p = _stageD_parameters_2D_v19(
            source, ua, outlet;
            felt_k=conductivity, felt_cp=capacity,
        )
        score = fit_set_score_2D_v19(
            _simulate_set_2D_v19(training_inputs, p);
            mode=:cooling,
        )
        candidate = (
            conductivity=conductivity,
            capacity=capacity,
            score=score,
        )
        push!(candidates, candidate)
        push!(rows, merge((
            stage="screen",
            felt_conductivity_scale=conductivity,
            felt_heat_capacity_scale=capacity,
        ), score))
        println(
            "v19-D felt k=", conductivity,
            " cp=", capacity,
            " J=", round(score.objective; digits=4),
        )
        flush(stdout)
    end
    # Screen-mesh ranking is only a preselection.  The best three pairs are
    # reranked on M so mesh error cannot choose a fitted thermal mass.
    nominal_candidates = NamedTuple[]
    for candidate in first(sort(
        candidates; by=x -> x.score.objective,
    ), 3)
        p = _stageD_parameters_2D_v19(
            source, ua, outlet;
            felt_k=candidate.conductivity,
            felt_cp=candidate.capacity,
            mesh=:nominal,
        )
        score = fit_set_score_2D_v19(
            _simulate_set_2D_v19(training_inputs, p);
            mode=:cooling,
        )
        nominal = (
            conductivity=candidate.conductivity,
            capacity=candidate.capacity,
            score=score,
        )
        push!(nominal_candidates, nominal)
        push!(rows, merge((
            stage="nominal_rerank",
            felt_conductivity_scale=candidate.conductivity,
            felt_heat_capacity_scale=candidate.capacity,
        ), score))
    end
    selected = first(sort(
        nominal_candidates; by=x -> x.score.objective,
    ))
    holdout_input = case_inputs_2D_v19(
        "C81"; cooling=true, max_points=61,
    )
    holdout = fit_set_score_2D_v19([
        simulate_case_2D_v19(
            holdout_input,
            _stageD_parameters_2D_v19(
                source, ua, outlet;
                felt_k=selected.conductivity,
                felt_cp=selected.capacity,
                mesh=:nominal,
            );
            reltol=5e-4, abstol=1e-4, dtmax=120.0,
        ),
    ]; mode=:cooling)
    _write_namedtuples_csv_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19, "stageD_felt_trace_2D_v19.csv",
    ), rows)
    interior =
        selected.conductivity != first(conductivities) &&
        selected.conductivity != last(conductivities) &&
        selected.capacity != first(capacities) &&
        selected.capacity != last(capacities)
    pass = interior &&
        selected.score.side_rmse_K <= 35.0 &&
        selected.score.air_rmse_K <= 30.0 &&
        selected.score.felt_rmse_K <= 15.0 &&
        holdout.side_rmse_K <= 35.0 &&
        holdout.air_rmse_K <= 30.0 &&
        holdout.felt_rmse_K <= 12.0 &&
        abs(holdout.air_final_bias_K) <= 20.0 &&
        abs(holdout.felt_final_bias_K) <= 8.0 &&
        (3.0 * holdout.side_rmse_K +
         holdout.air_rmse_K + holdout.felt_rmse_K) / 5.0 <= 27.0
    return (
        selected=selected, holdout=holdout,
        interior=interior, pass=pass,
    )
end

function _fit_power_2D_v19(source, ua, outlet, felt)
    selected_scales = [1.195, 1.405, 0.975]
    rows = NamedTuple[]
    group_results = NamedTuple[]
    for (group_index, group) in enumerate(POWER_GROUPS_2D_v19)
        training_inputs = [
            case_inputs_2D_v19(id; max_points=61)
            for id in group.train
        ]
        candidates = NamedTuple[]
        for value in group.values
            scales = copy(selected_scales)
            scales[group_index] = value
            p = _stageD_parameters_2D_v19(
                source, ua, outlet;
                felt_k=felt.conductivity,
                felt_cp=felt.capacity,
                power_scales=Tuple(scales),
                mesh=:nominal,
            )
            score = fit_set_score_2D_v19(
                _simulate_set_2D_v19(training_inputs, p),
            )
            push!(candidates, (value=value, score=score))
            push!(rows, merge((
                stage="training",
                power_group=group.name,
                power_scale=value,
            ), score))
            println(
                "v19-D power ", group.name,
                " scale=", value,
                " J=", round(score.objective; digits=4),
            )
            flush(stdout)
        end
        selected = first(sort(
            candidates; by=x -> x.score.objective,
        ))
        selected_scales[group_index] = selected.value
        holdout_inputs = [
            case_inputs_2D_v19(id; max_points=61)
            for id in group.valid
        ]
        holdout = fit_set_score_2D_v19(_simulate_set_2D_v19(
            holdout_inputs,
            _stageD_parameters_2D_v19(
                source, ua, outlet;
                felt_k=felt.conductivity,
                felt_cp=felt.capacity,
                power_scales=Tuple(selected_scales),
                mesh=:nominal,
            ),
        ))
        push!(rows, merge((
            stage="holdout",
            power_group=group.name,
            power_scale=selected.value,
        ), holdout))
        interior = selected.value != first(group.values) &&
                   selected.value != last(group.values)
        push!(group_results, (
            name=group.name, selected=selected,
            holdout=holdout, interior=interior,
        ))
    end
    _write_namedtuples_csv_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19, "stageD_power_trace_2D_v19.csv",
    ), rows)
    training_pass = all(
        result.selected.score.side_rmse_K <= 75.0 &&
        result.selected.score.air_rmse_K <= 45.0 &&
        result.selected.score.felt_rmse_K <= 30.0
        for result in group_results
    )
    holdout_pass = all(
        result.holdout.side_rmse_K <= 75.0 &&
        result.holdout.air_rmse_K <= 45.0 &&
        result.holdout.felt_rmse_K <= 30.0
        for result in group_results
    )
    pass = training_pass && holdout_pass &&
           all(result.interior for result in group_results)
    return (
        scales=Tuple(selected_scales),
        groups=group_results,
        pass=pass,
    )
end

function calibrate_2D_v19_stageD()
    source = selected_source_2D_v19()
    ua = selected_ua_2D_v19()
    outlet = selected_stageC_2D_v19()
    felt = _fit_felt_2D_v19(source, ua, outlet)
    power = _fit_power_2D_v19(
        source, ua, outlet, felt.selected,
    )
    upstream = (
        stageA=parse(
            Bool, _read_keyvalues_2D_v19(joinpath(
                OUTPUT_DIR_2D_v19,
                "stageA_selected_2D_v19.txt",
            ))["stageA_pass"],
        ),
        stageB=parse(
            Bool, _read_keyvalues_2D_v19(joinpath(
                OUTPUT_DIR_2D_v19,
                "stageB_selected_2D_v19.txt",
            ))["stageB_pass"],
        ),
        stageC=parse(
            Bool, _read_keyvalues_2D_v19(joinpath(
                OUTPUT_DIR_2D_v19,
                "stageC_selected_2D_v19.txt",
            ))["stageC_pass"],
        ),
    )
    prevalidation_pass = upstream.stageA && upstream.stageB &&
                         upstream.stageC &&
                         felt.pass && power.pass
    open(joinpath(
        OUTPUT_DIR_2D_v19, "parameters_selected_2D_v19.txt",
    ), "w") do io
        println(io, "source_model=", source.model)
        println(io, "deep_fraction=", source.fraction)
        println(io, "deep_length_m=", source.length_m)
        println(io, "nu_prefactor=", ua.prefactor)
        println(io, "reynolds_exponent=", ua.exponent)
        println(io, "rear_tube_flange_contact_h_W_m2_K=",
                outlet.tube_h)
        println(io, "probe_capacity_areal_J_m2_K=",
                outlet.probe_capacity)
        println(io, "probe_stem_conductance_areal_W_m2_K=",
                outlet.probe_stem)
        println(io, "felt_conductivity_scale=",
                felt.selected.conductivity)
        println(io, "felt_heat_capacity_scale=",
                felt.selected.capacity)
        println(io, "felt_contact_scale=0.30")
        println(io, "power_scale_456=", power.scales[1])
        println(io, "power_scale_304=", power.scales[2])
        println(io, "power_scale_256=", power.scales[3])
        println(io, "stageD_felt_pass=", felt.pass)
        println(io, "stageD_power_pass=", power.pass)
        println(io, "v19_prevalidation_pass=", prevalidation_pass)
        println(io, "selection_mesh=mixed_C_screen_M_rerank")
        println(io, "T9_T10_in_objective=false")
    end
    println(
        "v19-D felt=", felt,
        " power=", power,
        " upstream=", upstream,
        " prevalidation_pass=", prevalidation_pass,
    )
    return (
        felt=felt, power=power, upstream=upstream,
        prevalidation_pass=prevalidation_pass,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    calibrate_2D_v19_stageD()
end
