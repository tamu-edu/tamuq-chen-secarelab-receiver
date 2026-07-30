# ============================================================================
# refine_2D_v16_cooling_boundary.jl
#
# The first cooling grid selected its upper felt-conductivity bound.  Extend
# that dimension while retaining neighboring contact and capacity values.
# ============================================================================

include("calibrate_2D_v16_cooling.jl")

const CONTACT_REFINE_2D_v16 = (0.30, 1.0, 3.0)
const FELT_K_REFINE_2D_v16 = (2.40, 3.60, 4.80, 7.20)
const FELT_CP_REFINE_2D_v16 = (0.75, 1.50)

function refine_cooling_boundary_2D_v16()
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
    evaluation = 0
    total = length(CONTACT_REFINE_2D_v16) *
            length(FELT_K_REFINE_2D_v16) *
            length(FELT_CP_REFINE_2D_v16)
    for contact in CONTACT_REFINE_2D_v16,
        felt_k in FELT_K_REFINE_2D_v16,
        felt_cp in FELT_CP_REFINE_2D_v16
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
        (best === nothing || loss.total < best.loss.total) &&
            (best = candidate)
        push!(rows, (
            evaluation=evaluation, objective=loss.total,
            primary_rmse_K=loss.primary_rmse_K,
            t2_rmse_K=loss.t2_rmse_K,
            skin_felt_contact_scale=contact,
            felt_conductivity_scale=felt_k,
            felt_heat_capacity_scale=felt_cp,
        ))
        println(
            "v16 refine $evaluation/$total " *
            "loss=$(round(loss.total,digits=5)) " *
            "contact=$contact k=$felt_k Cp=$felt_cp",
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
            "cooling_boundary_refinement_2D_v16.csv",
        ), rows,
    )
    open(joinpath(
        OUTPUT_DIR_2D_v16,
        "parameters_cooling_refined_2D_v16.txt",
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
    end
    println("v16 refined selected = ", best)
    println("v16 refined holdout = ", holdout)
    return (selected=best, holdout=holdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    refine_cooling_boundary_2D_v16()
end
