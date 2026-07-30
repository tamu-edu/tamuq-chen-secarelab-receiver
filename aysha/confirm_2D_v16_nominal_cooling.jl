# ============================================================================
# confirm_2D_v16_nominal_cooling.jl
#
# Re-rank the best screen candidates directly on the nominal mesh after the
# C69 screen-to-nominal change exceeded the declared transfer threshold.
# ============================================================================

include("calibrate_2D_v16_cooling.jl")

const NOMINAL_CANDIDATES_2D_v16 = (
    (1.0, 2.4, 1.50),
    (3.0, 1.2, 0.75),
    (1.0, 2.4, 0.75),
    (10.0, 0.6, 0.75),
    (1.0, 1.2, 0.75),
    (10.0, 1.2, 0.75),
    (3.0, 2.4, 0.75),
    (3.0, 0.6, 0.75),
    (3.0, 1.2, 1.50),
    (0.3, 7.2, 1.50),
)

function confirm_nominal_cooling_2D_v16()
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
    rows = NamedTuple[]
    best = nothing
    for (index, candidate) in enumerate(NOMINAL_CANDIDATES_2D_v16)
        contact, felt_k, felt_cp = candidate
        p = parameters_2D_v16(
            nominal_mesh=true,
            skin_felt_contact_scale=contact,
            felt_conductivity_scale=felt_k,
            felt_heat_capacity_scale=felt_cp,
        )
        cases = _simulate_cooling_set_2D_v16(
            training_inputs, p,
        )
        loss = cooling_loss_2D_v16(cases)
        result = (
            contact=contact, felt_k=felt_k, felt_cp=felt_cp,
            loss=loss, parameters=p,
        )
        (best === nothing || loss.total < best.loss.total) &&
            (best = result)
        push!(rows, (
            candidate=index, objective=loss.total,
            primary_rmse_K=loss.primary_rmse_K,
            t2_rmse_K=loss.t2_rmse_K,
            skin_felt_contact_scale=contact,
            felt_conductivity_scale=felt_k,
            felt_heat_capacity_scale=felt_cp,
        ))
        println(
            "v16 nominal cooling $index/" *
            "$(length(NOMINAL_CANDIDATES_2D_v16)) " *
            "loss=$(round(loss.total,digits=5)) " *
            "primary=$(round(loss.primary_rmse_K,digits=2))K " *
            "parameters=$candidate",
        )
        flush(stdout)
    end
    holdout = cooling_loss_2D_v16(
        _simulate_cooling_set_2D_v16(
            holdout_inputs, best.parameters,
        ),
    )
    _write_namedtuples_csv_2D_v16(
        joinpath(
            OUTPUT_DIR_2D_v16,
            "cooling_nominal_confirmation_2D_v16.csv",
        ), rows,
    )
    open(joinpath(
        OUTPUT_DIR_2D_v16,
        "parameters_cooling_nominal_selected_2D_v16.txt",
    ), "w") do io
        println(io, "training_objective=", best.loss.total)
        println(
            io, "training_primary_rmse_K=",
            best.loss.primary_rmse_K,
        )
        println(
            io, "holdout_primary_rmse_K=",
            holdout.primary_rmse_K,
        )
        println(io, "skin_felt_contact_scale=", best.contact)
        println(io, "felt_conductivity_scale=", best.felt_k)
        println(io, "felt_heat_capacity_scale=", best.felt_cp)
        println(io, "selected_on_nominal_mesh=true")
    end
    println("v16 nominal selected = ", best)
    println("v16 nominal holdout = ", holdout)
    return (selected=best, holdout=holdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    confirm_nominal_cooling_2D_v16()
end
