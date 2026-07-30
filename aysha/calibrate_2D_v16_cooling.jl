# ============================================================================
# calibrate_2D_v16_cooling.jl - fit felt path on C69/C80; hold out C81
# ============================================================================

using Statistics

include("run_2D_v16.jl")

const CONTACT_GRID_2D_v16 = (0.10, 0.30, 1.0, 3.0, 10.0)
const FELT_K_GRID_2D_v16 = (0.60, 1.20, 2.40)
const FELT_CP_GRID_2D_v16 = (0.75, 1.50, 3.00)
const COOLING_TRAIN_2D_v16 = ("C69", "C80")
const COOLING_HOLDOUT_2D_v16 = ("C81",)

function cooling_loss_2D_v16(cases)
    primary_side = (1, 2, 3)
    guard = (4, 5, 6)
    side_transient = mean(
        mean(abs2, (
            case.model[:, index] .- case.observed[:, index]
        ) ./ 35.0)
        for case in cases for index in primary_side
    )
    felt_transient = mean(
        mean(abs2, (
            case.model[:, 7] .- case.observed[:, 7]
        ) ./ 15.0)
        for case in cases
    )
    guard_transient = mean(
        mean(abs2, (
            case.model[:, index] .- case.observed[:, index]
        ) ./ 50.0)
        for case in cases for index in guard
    )
    side_final = mean(
        (
            (case.model[end, index] -
             case.observed[end, index]) / 35.0
        )^2
        for case in cases for index in primary_side
    )
    felt_final = mean(
        (
            (case.model[end, 7] -
             case.observed[end, 7]) / 15.0
        )^2
        for case in cases
    )
    total = side_transient + felt_transient +
            0.20 * guard_transient +
            0.30 * side_final + 0.20 * felt_final
    primary_rmse = sqrt(mean(
        mean(abs2, (
            case.model[:, index] - case.observed[:, index]
        ))
        for case in cases for index in (1, 2, 3, 7)
    ))
    t2_rmse = sqrt(mean(
        mean(abs2, (
            case.model[:, 7] - case.observed[:, 7]
        ))
        for case in cases
    ))
    return (
        total=total,
        side_transient=side_transient,
        felt_transient=felt_transient,
        guard_transient=guard_transient,
        side_final=side_final,
        felt_final=felt_final,
        primary_rmse_K=primary_rmse,
        t2_rmse_K=t2_rmse,
    )
end

function _simulate_cooling_set_2D_v16(inputs, p)
    cases = Vector{Any}(undef, length(inputs))
    Threads.@threads for index in eachindex(inputs)
        cases[index] = simulate_case_2D_v16(
            inputs[index], p;
            reltol=3e-4, abstol=1e-4, dtmax=120.0,
        )
    end
    return cases
end

function fit_cooling_stage_2D_v16()
    mkpath(OUTPUT_DIR_2D_v16)
    training_inputs = [
        case_inputs_2D_v16(
            id; cooling=true, max_points=61,
        ) for id in COOLING_TRAIN_2D_v16
    ]
    holdout_inputs = [
        case_inputs_2D_v16(
            id; cooling=true, max_points=61,
        ) for id in COOLING_HOLDOUT_2D_v16
    ]
    trace = NamedTuple[]
    candidates = NamedTuple[]
    evaluation = 0
    total = length(CONTACT_GRID_2D_v16) *
            length(FELT_K_GRID_2D_v16) *
            length(FELT_CP_GRID_2D_v16)
    for contact in CONTACT_GRID_2D_v16,
        felt_k in FELT_K_GRID_2D_v16,
        felt_cp in FELT_CP_GRID_2D_v16
        evaluation += 1
        p = parameters_2D_v16(
            nominal_mesh=false, screen_mesh=true,
            skin_felt_contact_scale=contact,
            felt_conductivity_scale=felt_k,
            felt_heat_capacity_scale=felt_cp,
        )
        cases = _simulate_cooling_set_2D_v16(
            training_inputs, p,
        )
        loss = cooling_loss_2D_v16(cases)
        candidate = (
            contact=contact, felt_k=felt_k, felt_cp=felt_cp,
            loss=loss, parameters=p,
        )
        push!(candidates, candidate)
        push!(trace, (
            evaluation=evaluation, objective=loss.total,
            primary_rmse_K=loss.primary_rmse_K,
            t2_rmse_K=loss.t2_rmse_K,
            side_transient_loss=loss.side_transient,
            felt_transient_loss=loss.felt_transient,
            guard_transient_loss=loss.guard_transient,
            side_final_loss=loss.side_final,
            felt_final_loss=loss.felt_final,
            skin_felt_contact_scale=contact,
            felt_conductivity_scale=felt_k,
            felt_heat_capacity_scale=felt_cp,
        ))
        println(
            "v16 cooling $evaluation/$total " *
            "loss=$(round(loss.total,digits=5)) " *
            "primary=$(round(loss.primary_rmse_K,digits=2))K " *
            "contact=$contact k=$felt_k Cp=$felt_cp",
        )
        flush(stdout)
    end
    sort!(candidates; by=candidate -> candidate.loss.total)
    selected = first(candidates)
    holdout_cases = _simulate_cooling_set_2D_v16(
        holdout_inputs, selected.parameters,
    )
    holdout_loss = cooling_loss_2D_v16(holdout_cases)
    _write_namedtuples_csv_2D_v16(
        joinpath(
            OUTPUT_DIR_2D_v16,
            "cooling_calibration_trace_2D_v16.csv",
        ), trace,
    )
    open(joinpath(
        OUTPUT_DIR_2D_v16,
        "parameters_cooling_selected_2D_v16.txt",
    ), "w") do io
        println(io, "training_objective=", selected.loss.total)
        println(
            io, "training_primary_rmse_K=",
            selected.loss.primary_rmse_K,
        )
        println(
            io, "holdout_primary_rmse_K=",
            holdout_loss.primary_rmse_K,
        )
        println(
            io, "skin_felt_contact_scale=", selected.contact,
        )
        println(
            io, "felt_conductivity_scale=", selected.felt_k,
        )
        println(
            io, "felt_heat_capacity_scale=", selected.felt_cp,
        )
        println(io, "heating_parameters_frozen=true")
        println(io, "cooling_C81_held_out=true")
    end
    println("v16 cooling selected = ", selected)
    println("v16 C81 holdout = ", holdout_loss)
    return (selected=selected, holdout=holdout_loss)
end

if abspath(PROGRAM_FILE) == @__FILE__
    fit_cooling_stage_2D_v16()
end
