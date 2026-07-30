# ============================================================================
# calibrate_2D_v17_casing_flange.jl
#
# Cooling identification of the verified aluminum-casing/flange link.
# Felt k and skin/felt contact are fixed to avoid the v16 series-path nullspace.
# ============================================================================

using Statistics

include("run_2D_v17.jl")

const CASING_FLANGE_G_GRID_2D_v17 = (
    0.0, 0.02, 0.05, 0.10, 0.20, 0.50,
    1.0, 2.0, 5.0, 10.0, 20.0,
)
const FELT_CP_GRID_2D_v17 = (0.75, 1.00, 1.50, 2.00)
const COOLING_TRAIN_2D_v17 = ("C69", "C80")
const COOLING_HOLDOUT_2D_v17 = ("C81",)

function cooling_loss_2D_v17(cases)
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

function _simulate_cooling_set_2D_v17(inputs, p)
    cases = Vector{Any}(undef, length(inputs))
    Threads.@threads for index in eachindex(inputs)
        cases[index] = simulate_case_2D_v17(
            inputs[index], p;
            reltol=3e-4, abstol=1e-4, dtmax=120.0,
        )
    end
    return cases
end

function _evaluate_cooling_candidate_2D_v17(
    inputs, conductance, felt_cp; nominal_mesh,
)
    p = parameters_2D_v17(
        nominal_mesh=nominal_mesh,
        screen_mesh=!nominal_mesh,
        skin_felt_contact_scale=0.30,
        felt_conductivity_scale=1.00,
        felt_heat_capacity_scale=felt_cp,
        casing_flange_conductance_W_K=conductance,
    )
    cases = _simulate_cooling_set_2D_v17(inputs, p)
    return (
        conductance=conductance, felt_cp=felt_cp,
        loss=cooling_loss_2D_v17(cases), parameters=p,
    )
end

function fit_casing_flange_stage_2D_v17()
    mkpath(OUTPUT_DIR_2D_v17)
    screen_inputs = [
        case_inputs_2D_v17(id; cooling=true, max_points=41)
        for id in COOLING_TRAIN_2D_v17
    ]
    nominal_inputs = [
        case_inputs_2D_v17(id; cooling=true, max_points=61)
        for id in COOLING_TRAIN_2D_v17
    ]
    holdout_inputs = [
        case_inputs_2D_v17(id; cooling=true, max_points=61)
        for id in COOLING_HOLDOUT_2D_v17
    ]

    screen = NamedTuple[]
    total = length(CASING_FLANGE_G_GRID_2D_v17) *
            length(FELT_CP_GRID_2D_v17)
    evaluation = 0
    for conductance in CASING_FLANGE_G_GRID_2D_v17,
        felt_cp in FELT_CP_GRID_2D_v17
        evaluation += 1
        candidate = _evaluate_cooling_candidate_2D_v17(
            screen_inputs, conductance, felt_cp;
            nominal_mesh=false,
        )
        push!(screen, candidate)
        println(
            "v17 flange screen $evaluation/$total " *
            "loss=$(round(candidate.loss.total,digits=5)) " *
            "primary=$(round(candidate.loss.primary_rmse_K,digits=2))K " *
            "T2=$(round(candidate.loss.t2_rmse_K,digits=2))K " *
            "G=$conductance Cp=$felt_cp",
        )
        flush(stdout)
    end
    sort!(screen; by=candidate -> candidate.loss.total)

    # Nominal-mesh reranking prevents the screen mesh from selecting the
    # parameter.  Include the best 12 plus every zero-G baseline.
    keys_to_confirm = Set(
        (candidate.conductance, candidate.felt_cp)
        for candidate in screen[1:min(12, length(screen))]
    )
    foreach(
        cp -> push!(keys_to_confirm, (0.0, cp)),
        FELT_CP_GRID_2D_v17,
    )
    nominal = NamedTuple[]
    for (index, key) in enumerate(sort!(collect(keys_to_confirm)))
        candidate = _evaluate_cooling_candidate_2D_v17(
            nominal_inputs, key[1], key[2]; nominal_mesh=true,
        )
        push!(nominal, candidate)
        println(
            "v17 flange nominal $index/$(length(keys_to_confirm)) " *
            "loss=$(round(candidate.loss.total,digits=5)) " *
            "primary=$(round(candidate.loss.primary_rmse_K,digits=2))K " *
            "T2=$(round(candidate.loss.t2_rmse_K,digits=2))K " *
            "G=$(key[1]) Cp=$(key[2])",
        )
        flush(stdout)
    end
    sort!(nominal; by=candidate -> candidate.loss.total)
    selected = first(nominal)
    holdout_cases = _simulate_cooling_set_2D_v17(
        holdout_inputs, selected.parameters,
    )
    holdout = cooling_loss_2D_v17(holdout_cases)

    screen_rows = [
        (
            rank=index,
            objective=candidate.loss.total,
            primary_rmse_K=candidate.loss.primary_rmse_K,
            t2_rmse_K=candidate.loss.t2_rmse_K,
            casing_flange_conductance_W_K=
                candidate.conductance,
            felt_heat_capacity_scale=candidate.felt_cp,
        )
        for (index, candidate) in enumerate(screen)
    ]
    nominal_rows = [
        (
            rank=index,
            objective=candidate.loss.total,
            primary_rmse_K=candidate.loss.primary_rmse_K,
            t2_rmse_K=candidate.loss.t2_rmse_K,
            casing_flange_conductance_W_K=
                candidate.conductance,
            felt_heat_capacity_scale=candidate.felt_cp,
        )
        for (index, candidate) in enumerate(nominal)
    ]
    _write_namedtuples_csv_2D_v17(
        joinpath(
            OUTPUT_DIR_2D_v17,
            "casing_flange_screen_2D_v17.csv",
        ), screen_rows,
    )
    _write_namedtuples_csv_2D_v17(
        joinpath(
            OUTPUT_DIR_2D_v17,
            "casing_flange_nominal_confirmation_2D_v17.csv",
        ), nominal_rows,
    )
    open(joinpath(
        OUTPUT_DIR_2D_v17,
        "parameters_casing_flange_selected_2D_v17.txt",
    ), "w") do io
        println(io, "training_objective=", selected.loss.total)
        println(
            io, "training_primary_rmse_K=",
            selected.loss.primary_rmse_K,
        )
        println(
            io, "training_t2_rmse_K=",
            selected.loss.t2_rmse_K,
        )
        println(
            io, "holdout_primary_rmse_K=",
            holdout.primary_rmse_K,
        )
        println(io, "holdout_t2_rmse_K=", holdout.t2_rmse_K)
        println(
            io, "casing_flange_conductance_W_K=",
            selected.conductance,
        )
        println(
            io, "felt_heat_capacity_scale=",
            selected.felt_cp,
        )
        println(io, "skin_felt_contact_scale_fixed=0.30")
        println(io, "felt_conductivity_scale_fixed=1.00")
        println(io, "selected_on_nominal_mesh=true")
        println(io, "cooling_C81_held_out=true")
    end
    println("v17 casing/flange selected = ", selected)
    println("v17 C81 holdout = ", holdout)
    return (selected=selected, holdout=holdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    fit_casing_flange_stage_2D_v17()
end
