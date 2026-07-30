# ============================================================================
# confirm_2D_v13_boundary.jl - bound-expansion check after staged v13 screen
# ============================================================================

include("calibrate_2D_v13_staged.jl")

function confirm_boundary_2D_v13()
    mkpath(OUTPUT_DIR_2D_v13)
    base = inherited_parameters_2D_v13(
        nominal_mesh=false, screen_mesh=true,
    )
    heating_inputs = [
        case_inputs_2D_v13(id; max_points=FIT_POINTS_2D_v13)
        for id in TRAIN_HEATING_2D_v13
    ]
    cooling_inputs = [
        case_inputs_2D_v13(
            id; cooling=true, max_points=FIT_POINTS_2D_v13,
        )
        for id in COOLING_IDS_2D_v13
    ]
    candidates = [
        (f, ratio * f, G, 5.0)
        for f in (0.025, 0.050, 0.075, 0.100)
        for ratio in (0.30, 0.50)
        for G in (30.0, 100.0)
    ]
    results = Vector{Any}(undef, length(candidates))
    rows = NamedTuple[]
    for (index, theta) in enumerate(candidates)
        result = evaluate_population_2D_v13(
            base, theta, heating_inputs, cooling_inputs,
        )
        results[index] = result
        push!(rows, (
            evaluation=index,
            objective=result.loss.total,
            transient_loss=result.loss.transient,
            steady_loss=result.loss.steady,
            radial_loss=result.loss.radial,
            axial_loss=result.loss.axial,
            cooling_loss=result.loss.cooling,
            side_solid_fraction=theta[1],
            side_gas_heat_fraction=theta[2],
            active_side_conductance_W_K_m=theta[3],
            side_felt_contact_h_W_m2K=theta[4],
        ))
        println(
            "v13 boundary $index/$(length(candidates)) " *
            "loss=$(round(result.loss.total,digits=5)) theta=$theta",
        )
        flush(stdout)
    end
    best = results[argmin(result.loss.total for result in results)]
    _write_namedtuples_csv_2D_v13(
        joinpath(
            OUTPUT_DIR_2D_v13,
            "population_boundary_confirmation_2D_v13.csv",
        ),
        rows,
    )
    open(joinpath(
        OUTPUT_DIR_2D_v13,
        "parameters_boundary_selected_2D_v13.txt",
    ), "w") do io
        println(io, "objective=", best.loss.total)
        println(io, "side_solid_fraction=", best.theta[1])
        println(io, "side_gas_heat_fraction=", best.theta[2])
        println(io, "active_side_conductance_W_K_m=", best.theta[3])
        println(io, "side_felt_contact_h_W_m2K=", best.theta[4])
    end
    println("v13 boundary selected theta = ", best.theta)
    println("v13 boundary selected loss = ", best.loss)
    return best
end

if abspath(PROGRAM_FILE) == @__FILE__
    confirm_boundary_2D_v13()
end

