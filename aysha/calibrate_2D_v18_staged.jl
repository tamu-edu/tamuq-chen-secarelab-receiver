# ============================================================================
# calibrate_2D_v18_staged.jl
#
# Staged, bounded source/exchange/felt screen.  Selection uses the side-wall
# thermocouples, outlet air, and felt temperature only.  T9/T10 are excluded.
# ============================================================================

using Statistics

include("run_2D_v18.jl")

const SIDE_2D_v18 = (1, 2, 3)
const AIR_2D_v18 = 6
const FELT_2D_v18 = 7
const WALL_WEIGHTS_2D_v18 = (0.251825, 0.350365, 0.397810)

const POWER_GROUPS_2D_v18 = (
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

const TRAIN_COOLING_2D_v18 = ("C69", "C80")
const VALID_COOLING_2D_v18 = ("C81",)

rho_2D_v18(x) = 2.0 * (sqrt(1.0 + x * x) - 1.0)

function _robust_mean_2D_v18(values, scale)
    return mean(rho_2D_v18(value / scale) for value in values)
end

function _final_indices_2D_v18(n)
    first_index = max(2, floor(Int, 0.9 * n))
    return first_index:n
end

function case_objective_2D_v18(case)
    # t=0 is an initialization constraint, not an independent fitted datum.
    rows = 2:size(case.model, 1)
    final_rows = _final_indices_2D_v18(size(case.model, 1))
    side_curve = mean(
        _robust_mean_2D_v18(
            case.model[rows, index] .-
            case.observed[rows, index], 35.0,
        ) for index in SIDE_2D_v18
    )
    air_curve = _robust_mean_2D_v18(
        case.model[rows, AIR_2D_v18] .-
        case.observed[rows, AIR_2D_v18], 25.0,
    )
    felt_curve = _robust_mean_2D_v18(
        case.model[rows, FELT_2D_v18] .-
        case.observed[rows, FELT_2D_v18], 15.0,
    )
    curve = 0.50 * side_curve + 0.30 * air_curve +
            0.20 * felt_curve

    side_level = mean(
        rho_2D_v18((
            mean(case.model[final_rows, index]) -
            mean(case.observed[final_rows, index])
        ) / 35.0) for index in SIDE_2D_v18
    )
    air_level = rho_2D_v18((
        mean(case.model[final_rows, AIR_2D_v18]) -
        mean(case.observed[final_rows, AIR_2D_v18])
    ) / 25.0)
    felt_level = rho_2D_v18((
        mean(case.model[final_rows, FELT_2D_v18]) -
        mean(case.observed[final_rows, FELT_2D_v18])
    ) / 15.0)
    level = 0.45 * side_level + 0.35 * air_level +
            0.20 * felt_level

    model_side = [
        mean(case.model[final_rows, index]) -
        mean(case.inputs.ambient[final_rows])
        for index in SIDE_2D_v18
    ]
    observed_side = [
        mean(case.observed[final_rows, index]) -
        mean(case.inputs.ambient[final_rows])
        for index in SIDE_2D_v18
    ]
    model_shape = model_side ./ max(sum(model_side), 1.0)
    observed_shape = observed_side ./ max(sum(observed_side), 1.0)
    shape = _robust_mean_2D_v18(
        model_shape .- observed_shape, 0.05,
    )

    model_wall = sum(
        WALL_WEIGHTS_2D_v18[i] * model_side[i] for i in 1:3
    )
    observed_wall = sum(
        WALL_WEIGHTS_2D_v18[i] * observed_side[i] for i in 1:3
    )
    model_air = mean(case.model[final_rows, AIR_2D_v18])
    observed_air = mean(case.observed[final_rows, AIR_2D_v18])
    inlet = mean(case.inputs.inlet[final_rows])
    model_effectiveness = (model_air - inlet) /
                          max(model_wall, 1.0)
    observed_effectiveness = (observed_air - inlet) /
                             max(observed_wall, 1.0)
    effectiveness = rho_2D_v18(
        (model_effectiveness - observed_effectiveness) / 0.04,
    )
    total = 0.35 * curve + 0.25 * level +
            0.20 * shape + 0.20 * effectiveness
    return (
        total=total, curve=curve, level=level,
        side_shape=shape, effectiveness=effectiveness,
        model_effectiveness=model_effectiveness,
        observed_effectiveness=observed_effectiveness,
    )
end

function set_objective_2D_v18(cases)
    losses = case_objective_2D_v18.(cases)
    return (
        total=mean(loss.total for loss in losses),
        curve=mean(loss.curve for loss in losses),
        level=mean(loss.level for loss in losses),
        side_shape=mean(loss.side_shape for loss in losses),
        effectiveness=mean(loss.effectiveness for loss in losses),
    )
end

function cooling_set_objective_2D_v18(cases)
    case_losses = map(cases) do case
        rows = 2:size(case.model, 1)
        side = mean(
            _robust_mean_2D_v18(
                case.model[rows, index] .-
                case.observed[rows, index], 35.0,
            ) for index in SIDE_2D_v18
        )
        air = _robust_mean_2D_v18(
            case.model[rows, AIR_2D_v18] .-
            case.observed[rows, AIR_2D_v18], 25.0,
        )
        felt = _robust_mean_2D_v18(
            case.model[rows, FELT_2D_v18] .-
            case.observed[rows, FELT_2D_v18], 15.0,
        )
        (total=0.40 * side + 0.30 * air + 0.30 * felt,
         side=side, air=air, felt=felt)
    end
    return (
        total=mean(x.total for x in case_losses),
        side=mean(x.side for x in case_losses),
        air=mean(x.air for x in case_losses),
        felt=mean(x.felt for x in case_losses),
    )
end

function _simulate_set_2D_v18(inputs, p; tight=false)
    cases = Vector{Any}(undef, length(inputs))
    Threads.@threads for index in eachindex(inputs)
        cases[index] = simulate_case_2D_v18(
            inputs[index], p;
            reltol=tight ? 1e-6 : 5e-4,
            abstol=tight ? 1e-7 : 1e-4,
            dtmax=tight ? 30.0 : 120.0,
        )
    end
    return cases
end

function _parameter_tuple_2D_v18(candidate; kwargs...)
    return parameters_2D_v18(;
        mesh=get(kwargs, :mesh, :screen),
        source_model=candidate.source_model,
        deep_fraction=candidate.deep_fraction,
        deep_length_m=candidate.deep_length_m,
        exchange_multiplier=candidate.exchange_multiplier,
        felt_conductivity_scale=get(
            kwargs, :felt_conductivity_scale, 1.0,
        ),
        felt_heat_capacity_scale=get(
            kwargs, :felt_heat_capacity_scale, 1.0,
        ),
        felt_contact_scale=get(kwargs, :felt_contact_scale, 0.30),
        power_scales=get(
            kwargs, :power_scales, (1.195, 1.405, 0.975),
        ),
    )
end

function _candidate_label_2D_v18(candidate)
    return string(
        candidate.source_model, "_f",
        candidate.deep_fraction, "_L",
        round(1e3 * candidate.deep_length_m; digits=0),
        "_x", candidate.exchange_multiplier,
    )
end

function _evaluate_candidate_2D_v18(
    candidate, inputs; power_scales=(1.195, 1.405, 0.975),
    mesh=:screen, felt_conductivity_scale=1.0,
    felt_heat_capacity_scale=1.0,
)
    p = _parameter_tuple_2D_v18(
        candidate; mesh=mesh, power_scales=power_scales,
        felt_conductivity_scale=felt_conductivity_scale,
        felt_heat_capacity_scale=felt_heat_capacity_scale,
    )
    cases = _simulate_set_2D_v18(inputs, p)
    return (loss=set_objective_2D_v18(cases), cases=cases, p=p)
end

function _screen_candidates_2D_v18()
    # One interleaved training case per irradiance group keeps the broad
    # structural screen inexpensive.  All nine train cases enter nested power.
    representative = [
        case_inputs_2D_v18(id; max_points=61)
        for id in ("E69", "E74", "E79")
    ]
    rows = NamedTuple[]
    baseline = (
        source_model=:beer_lambert, deep_fraction=0.0,
        deep_length_m=0.08, exchange_multiplier=1.0,
    )
    candidates = [baseline]
    append!(candidates, [(
        source_model=:beer_lambert, deep_fraction=0.0,
        deep_length_m=0.08, exchange_multiplier=x,
    ) for x in (0.05, 0.10, 0.15, 0.25, 0.35, 0.45, 0.60, 0.75)])
    append!(candidates, [(
        source_model=:near_deep, deep_fraction=f,
        deep_length_m=L, exchange_multiplier=1.0,
    ) for f in (0.20, 0.40, 0.60, 0.80)
      for L in (0.02, 0.05, 0.10, 0.20)])

    evaluated = NamedTuple[]
    for (index, candidate) in enumerate(candidates)
        result = _evaluate_candidate_2D_v18(
            candidate, representative,
        )
        push!(evaluated, merge(candidate, (loss=result.loss.total,)))
        push!(rows, (
            stage="marginal", candidate=_candidate_label_2D_v18(candidate),
            objective=result.loss.total, curve=result.loss.curve,
            level=result.loss.level, side_shape=result.loss.side_shape,
            effectiveness=result.loss.effectiveness,
        ))
        println(
            "v18 marginal $index/$(length(candidates)) ",
            _candidate_label_2D_v18(candidate),
            " J=", round(result.loss.total; digits=4),
        )
        flush(stdout)
    end
    best_exchange = first(sort(
        filter(x -> x.source_model === :beer_lambert &&
                    x.exchange_multiplier < 1.0, evaluated);
        by=x -> x.loss,
    ))
    best_sources = first(sort(
        filter(x -> x.source_model === :near_deep, evaluated);
        by=x -> x.loss,
    ), 4)
    exchange_values = unique((
        best_exchange.exchange_multiplier,
        max(0.05, best_exchange.exchange_multiplier - 0.05),
        min(0.75, best_exchange.exchange_multiplier + 0.10),
    ))
    joint = [(
        source_model=:near_deep,
        deep_fraction=source.deep_fraction,
        deep_length_m=source.deep_length_m,
        exchange_multiplier=x,
    ) for source in best_sources for x in exchange_values]
    joint_evaluated = NamedTuple[]
    for (index, candidate) in enumerate(joint)
        result = _evaluate_candidate_2D_v18(
            candidate, representative,
        )
        push!(joint_evaluated, merge(
            candidate, (loss=result.loss.total,),
        ))
        push!(rows, (
            stage="joint", candidate=_candidate_label_2D_v18(candidate),
            objective=result.loss.total, curve=result.loss.curve,
            level=result.loss.level, side_shape=result.loss.side_shape,
            effectiveness=result.loss.effectiveness,
        ))
        println(
            "v18 joint $index/$(length(joint)) ",
            _candidate_label_2D_v18(candidate),
            " J=", round(result.loss.total; digits=4),
        )
        flush(stdout)
    end
    structural = (
        baseline,
        (
            source_model=:beer_lambert, deep_fraction=0.0,
            deep_length_m=0.08,
            exchange_multiplier=best_exchange.exchange_multiplier,
        ),
        (
            source_model=:near_deep,
            deep_fraction=best_sources[1].deep_fraction,
            deep_length_m=best_sources[1].deep_length_m,
            exchange_multiplier=1.0,
        ),
        let best_joint = first(sort(joint_evaluated; by=x -> x.loss))
            (
                source_model=best_joint.source_model,
                deep_fraction=best_joint.deep_fraction,
                deep_length_m=best_joint.deep_length_m,
                exchange_multiplier=best_joint.exchange_multiplier,
            )
        end,
    )
    return (structural=structural, rows=rows)
end

function _profile_power_2D_v18(candidate)
    selected = [1.195, 1.405, 0.975]
    rows = NamedTuple[]
    for (group_index, group) in enumerate(POWER_GROUPS_2D_v18)
        inputs = [
            case_inputs_2D_v18(id; max_points=61)
            for id in group.train
        ]
        results = NamedTuple[]
        for value in group.values
            scales = copy(selected)
            scales[group_index] = value
            result = _evaluate_candidate_2D_v18(
                candidate, inputs; power_scales=Tuple(scales),
            )
            push!(results, (value=value, loss=result.loss.total))
            push!(rows, (
                candidate=_candidate_label_2D_v18(candidate),
                power_group=group.name, power_scale=value,
                objective=result.loss.total,
            ))
        end
        selected[group_index] = first(sort(results; by=x -> x.loss)).value
    end
    all_inputs = [
        case_inputs_2D_v18(id; max_points=61)
        for group in POWER_GROUPS_2D_v18 for id in group.train
    ]
    result = _evaluate_candidate_2D_v18(
        candidate, all_inputs; power_scales=Tuple(selected),
    )
    return (
        candidate=candidate, power_scales=Tuple(selected),
        loss=result.loss.total, rows=rows,
    )
end

function _profile_felt_2D_v18(profile)
    cooling_inputs = [
        case_inputs_2D_v18(id; cooling=true, max_points=61)
        for id in TRAIN_COOLING_2D_v18
    ]
    rows = NamedTuple[]
    results = NamedTuple[]
    # Contact stays fixed: it is structurally confounded with conductivity.
    for conductivity in (0.70, 1.00, 1.30)
        for capacity in (0.70, 1.00, 1.30)
            result = _evaluate_candidate_2D_v18(
                profile.candidate, cooling_inputs;
                power_scales=profile.power_scales,
                felt_conductivity_scale=conductivity,
                felt_heat_capacity_scale=capacity,
            )
            cooling_loss = cooling_set_objective_2D_v18(result.cases)
            push!(results, (
                conductivity=conductivity, capacity=capacity,
                loss=cooling_loss.total,
            ))
            push!(rows, (
                felt_conductivity_scale=conductivity,
                felt_heat_capacity_scale=capacity,
                objective=cooling_loss.total,
                side_objective=cooling_loss.side,
                air_objective=cooling_loss.air,
                felt_objective=cooling_loss.felt,
            ))
        end
    end
    best = first(sort(results; by=x -> x.loss))
    return (best=best, rows=rows)
end

function calibrate_2D_v18_staged()
    mkpath(OUTPUT_DIR_2D_v18)
    screen = _screen_candidates_2D_v18()
    _write_namedtuples_csv_2D_v18(
        joinpath(OUTPUT_DIR_2D_v18, "structural_screen_2D_v18.csv"),
        screen.rows,
    )
    profiles = [_profile_power_2D_v18(c) for c in screen.structural]
    power_rows = reduce(vcat, (profile.rows for profile in profiles))
    _write_namedtuples_csv_2D_v18(
        joinpath(OUTPUT_DIR_2D_v18, "nested_power_trace_2D_v18.csv"),
        power_rows,
    )
    comparison = [(
        candidate=_candidate_label_2D_v18(profile.candidate),
        objective=profile.loss,
        power_scale_456=profile.power_scales[1],
        power_scale_304=profile.power_scales[2],
        power_scale_256=profile.power_scales[3],
    ) for profile in profiles]
    _write_namedtuples_csv_2D_v18(
        joinpath(OUTPUT_DIR_2D_v18, "structural_profiled_2D_v18.csv"),
        comparison,
    )
    winner = first(sort(profiles; by=x -> x.loss))
    felt = _profile_felt_2D_v18(winner)
    _write_namedtuples_csv_2D_v18(
        joinpath(OUTPUT_DIR_2D_v18, "felt_profile_2D_v18.csv"),
        felt.rows,
    )
    open(joinpath(
        OUTPUT_DIR_2D_v18, "parameters_selected_2D_v18.txt",
    ), "w") do io
        println(io, "source_model=", winner.candidate.source_model)
        println(io, "deep_fraction=", winner.candidate.deep_fraction)
        println(io, "deep_length_m=", winner.candidate.deep_length_m)
        println(io, "exchange_multiplier=",
                winner.candidate.exchange_multiplier)
        println(io, "power_scale_456=", winner.power_scales[1])
        println(io, "power_scale_304=", winner.power_scales[2])
        println(io, "power_scale_256=", winner.power_scales[3])
        println(io, "felt_conductivity_scale=", felt.best.conductivity)
        println(io, "felt_heat_capacity_scale=", felt.best.capacity)
        println(io, "felt_contact_scale=0.30")
        println(io, "selection_mesh=screen")
        println(io, "T9_T10_in_objective=false")
    end
    println("v18 selected: ", winner, " felt=", felt.best)
    return (winner=winner, felt=felt.best)
end

if abspath(PROGRAM_FILE) == @__FILE__
    calibrate_2D_v18_staged()
end
