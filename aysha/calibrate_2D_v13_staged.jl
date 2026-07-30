# ============================================================================
# calibrate_2D_v13_staged.jl
#
# Fits only the four new, shared receiver-population quantities.  All v12
# assembly, observation, optical, Graetz, power, and hydraulic values remain
# fixed at their definitive v12 values.
# ============================================================================

using Statistics

include("run_2D_v13.jl")

const FIT_POINTS_2D_v13 = 61
const SCREEN_CANDIDATES_2D_v13 = Tuple(
    (
        side_fraction,
        gas_ratio * side_fraction,
        active_side_G,
        5.0,
    )
    for side_fraction in (0.15, 0.30, 0.45)
    for gas_ratio in (0.0, 0.15)
    for active_side_G in (0.2, 2.0, 10.0)
)

function population_loss_2D_v13(
    heating_cases,
    cooling_cases,
)
    transient = mean(
        mean(abs2, (case.model .- case.observed) ./ 100.0)
        for case in heating_cases
    )
    steady = mean(
        mean(abs2, (case.model[end, :] .-
                    case.observed[end, :]) ./ 75.0)
        for case in heating_cases
    )
    radial = mean(
        (
            (
                (case.model[end, 2] - case.model[end, 4]) -
                (case.observed[end, 2] - case.observed[end, 4])
            ) / 30.0
        )^2 +
        (
            (
                (case.model[end, 3] - case.model[end, 5]) -
                (case.observed[end, 3] - case.observed[end, 5])
            ) / 30.0
        )^2
        for case in heating_cases
    ) / 2.0
    axial = mean(
        (
            (
                (case.model[end, 2] - case.model[end, 1]) -
                (case.observed[end, 2] - case.observed[end, 1])
            ) / 50.0
        )^2
        for case in heating_cases
    )
    cooling = mean(
        mean(abs2, (case.model .- case.observed) ./ 50.0)
        for case in cooling_cases
    )
    total = transient + 0.30 * steady + 1.50 * radial +
            0.30 * axial + 0.50 * cooling
    return (
        total=total, transient=transient, steady=steady,
        radial=radial, axial=axial, cooling=cooling,
    )
end

function evaluate_population_2D_v13(
    base,
    theta,
    heating_inputs,
    cooling_inputs,
)
    side_fraction, gas_fraction, active_side_G, side_felt_h = theta
    p = rebuild_population_2D_v13(
        base;
        side_fraction=side_fraction,
        side_gas_fraction=gas_fraction,
        active_side_G_per_m=active_side_G,
        side_felt_h=side_felt_h,
    )
    heating_cases = Vector{Any}(undef, length(heating_inputs))
    Threads.@threads for index in eachindex(heating_inputs)
        heating_cases[index] = simulate_case_2D_v13(
            heating_inputs[index], p;
            reltol=3e-4, abstol=1e-4, dtmax=120.0,
        )
    end
    cooling_cases = Vector{Any}(undef, length(cooling_inputs))
    Threads.@threads for index in eachindex(cooling_inputs)
        cooling_cases[index] = simulate_case_2D_v13(
            cooling_inputs[index], p;
            reltol=3e-4, abstol=1e-4, dtmax=120.0,
        )
    end
    loss = population_loss_2D_v13(heating_cases, cooling_cases)
    return (
        parameters=p, heating_cases=heating_cases,
        cooling_cases=cooling_cases, loss=loss, theta=theta,
    )
end

function refinement_candidates_2D_v13(theta)
    f0, g0, G0, _ = theta
    ratio0 = g0 / f0
    raw = Tuple{Float64,Float64,Float64,Float64}[]
    for f in unique(clamp.((f0 - 0.075, f0, f0 + 0.075),
                           0.075, 0.60)),
        ratio in unique(clamp.((ratio0, 0.075, 0.30), 0.0, 0.50)),
        G in unique((max(0.05, G0 / 3.0), G0, min(30.0, 3.0 * G0))),
        hf in (1.0, 10.0)
        push!(raw, (Float64(f), Float64(ratio * f),
                    Float64(G), Float64(hf)))
    end
    # Rank-one variations plus a compact crossed subset keep the refinement
    # bounded enough for a full staged run.
    candidates = [
        (f0, g0, G0, h) for h in (0.0, 1.0, 5.0, 10.0, 30.0)
    ]
    append!(candidates, [
        (f, ratio * f, G, 5.0)
        for f in unique(clamp.((f0 - 0.075, f0, f0 + 0.075),
                               0.075, 0.60))
        for ratio in unique(clamp.((ratio0, 0.075, 0.30),
                                   0.0, 0.50))
        for G in (max(0.05, G0 / 3.0), min(30.0, 3.0 * G0))
    ])
    append!(candidates, raw[1:min(6, length(raw))])
    return unique(candidates)
end

function fit_population_2D_v13(;
    candidates=SCREEN_CANDIDATES_2D_v13,
)
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
    trace = NamedTuple[]
    screen = Vector{Any}(undef, length(candidates))
    for (index, theta) in enumerate(candidates)
        selected = evaluate_population_2D_v13(
            base, theta, heating_inputs, cooling_inputs,
        )
        screen[index] = selected
        push!(trace, (
            stage="screen", evaluation=index,
            objective=selected.loss.total,
            transient_loss=selected.loss.transient,
            steady_loss=selected.loss.steady,
            radial_loss=selected.loss.radial,
            axial_loss=selected.loss.axial,
            cooling_loss=selected.loss.cooling,
            side_solid_fraction=theta[1],
            side_gas_heat_fraction=theta[2],
            active_side_conductance_W_K_m=theta[3],
            side_felt_contact_h_W_m2K=theta[4],
        ))
        println(
            "v13 screen $index/$(length(candidates)) " *
            "loss=$(round(selected.loss.total,digits=5)) theta=$theta",
        )
        flush(stdout)
    end
    best_screen = screen[
        argmin(item.loss.total for item in screen)
    ]

    refinement = refinement_candidates_2D_v13(best_screen.theta)
    refined = Vector{Any}(undef, length(refinement))
    for (index, theta) in enumerate(refinement)
        selected = evaluate_population_2D_v13(
            base, theta, heating_inputs, cooling_inputs,
        )
        refined[index] = selected
        push!(trace, (
            stage="refinement", evaluation=index,
            objective=selected.loss.total,
            transient_loss=selected.loss.transient,
            steady_loss=selected.loss.steady,
            radial_loss=selected.loss.radial,
            axial_loss=selected.loss.axial,
            cooling_loss=selected.loss.cooling,
            side_solid_fraction=theta[1],
            side_gas_heat_fraction=theta[2],
            active_side_conductance_W_K_m=theta[3],
            side_felt_contact_h_W_m2K=theta[4],
        ))
        println(
            "v13 refine $index/$(length(refinement)) " *
            "loss=$(round(selected.loss.total,digits=5)) theta=$theta",
        )
        flush(stdout)
    end
    all_results = vcat(screen, refined)
    best = all_results[
        argmin(item.loss.total for item in all_results)
    ]
    _write_namedtuples_csv_2D_v13(
        joinpath(
            OUTPUT_DIR_2D_v13,
            "population_calibration_trace_2D_v13.csv",
        ),
        trace,
    )
    open(joinpath(
        OUTPUT_DIR_2D_v13,
        "parameters_population_selected_2D_v13.txt",
    ), "w") do io
        println(io, "objective=", best.loss.total)
        println(io, "side_solid_fraction=", best.theta[1])
        println(io, "side_gas_heat_fraction=", best.theta[2])
        println(io, "active_side_conductance_W_K_m=", best.theta[3])
        println(io, "side_felt_contact_h_W_m2K=", best.theta[4])
        println(io, "inherited_v12_parameters_frozen=true")
    end
    println("v13 selected theta = ", best.theta)
    println("v13 selected loss = ", best.loss)
    return (
        parameters=best.parameters, theta=best.theta,
        loss=best.loss, trace=trace,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    fit_population_2D_v13()
end
