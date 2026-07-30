# ============================================================================
# calibrate_2D_v14_staged.jl
#
# Stage 1 fits only square-network/profile quantities with all v12 heating
# coefficients frozen.  C69/C80/C81 penalize loss of cooling performance.
# ============================================================================

using Statistics

include("run_2D_v14.jl")

const FIT_POINTS_2D_v14 = 61
const LATERAL_SCALE_GRID_2D_v14 = (0.05, 0.20, 1.00)
const EDGE_FELT_SCALE_GRID_2D_v14 = (0.30, 1.00, 3.00)
const BEAM_SIGMA_GRID_2D_v14 = (14e-3, 30e-3, 100e-3)

function network_loss_2D_v14(heating_cases, cooling_cases)
    transient = mean(
        mean(abs2, (case.model .- case.observed) ./ 100.0)
        for case in heating_cases
    )
    steady = mean(
        mean(abs2, (
            case.model[end, :] .- case.observed[end, :]
        ) ./ 75.0)
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
    total = transient + 0.30 * steady + 2.00 * radial +
            0.30 * axial + 0.50 * cooling
    return (
        total=total, transient=transient, steady=steady,
        radial=radial, axial=axial, cooling=cooling,
    )
end

function _simulate_set_2D_v14(inputs, p)
    cases = Vector{Any}(undef, length(inputs))
    Threads.@threads for index in eachindex(inputs)
        cases[index] = simulate_case_2D_v14(
            inputs[index], p;
            reltol=3e-4, abstol=1e-4, dtmax=120.0,
        )
    end
    return cases
end

function fit_network_stage_2D_v14()
    mkpath(OUTPUT_DIR_2D_v14)
    base = inherited_parameters_2D_v14(
        nominal_mesh=false, screen_mesh=true,
    )
    heating_inputs = [
        case_inputs_2D_v14(id; max_points=FIT_POINTS_2D_v14)
        for id in TRAIN_HEATING_2D_v14
    ]
    cooling_inputs = [
        case_inputs_2D_v14(
            id; cooling=true, max_points=FIT_POINTS_2D_v14,
        )
        for id in COOLING_IDS_2D_v14
    ]
    trace = NamedTuple[]
    best = nothing
    evaluation = 0
    total = length(LATERAL_SCALE_GRID_2D_v14) *
            length(EDGE_FELT_SCALE_GRID_2D_v14) *
            length(BEAM_SIGMA_GRID_2D_v14)
    for lateral in LATERAL_SCALE_GRID_2D_v14,
        felt in EDGE_FELT_SCALE_GRID_2D_v14,
        sigma in BEAM_SIGMA_GRID_2D_v14
        evaluation += 1
        pnetwork = rebuild_network_2D_v14(
            base; lateral_scale=lateral,
            edge_felt_scale=felt,
        )
        p = rebuild_heating_2D_v14(
            pnetwork; beam_sigma=sigma,
        )
        heating_cases = _simulate_set_2D_v14(
            heating_inputs, p,
        )
        cooling_cases = _simulate_set_2D_v14(
            cooling_inputs, p,
        )
        loss = network_loss_2D_v14(
            heating_cases, cooling_cases,
        )
        candidate = (
            parameters=p, lateral=lateral,
            edge_felt=felt, beam_sigma=sigma, loss=loss,
        )
        if best === nothing || loss.total < best.loss.total
            best = candidate
        end
        push!(trace, (
            evaluation=evaluation,
            objective=loss.total,
            transient_loss=loss.transient,
            steady_loss=loss.steady,
            radial_loss=loss.radial,
            axial_loss=loss.axial,
            cooling_loss=loss.cooling,
            lateral_conductivity_scale=lateral,
            edge_felt_contact_scale=felt,
            beam_sigma_mm=1000.0 * sigma,
        ))
        println(
            "v14 network $evaluation/$total " *
            "loss=$(round(loss.total,digits=5)) " *
            "lateral=$lateral felt=$felt sigma=$(1000sigma)mm",
        )
        flush(stdout)
    end
    _write_namedtuples_csv_2D_v14(
        joinpath(
            OUTPUT_DIR_2D_v14,
            "network_calibration_trace_2D_v14.csv",
        ), trace,
    )
    open(joinpath(
        OUTPUT_DIR_2D_v14,
        "parameters_network_selected_2D_v14.txt",
    ), "w") do io
        println(io, "objective=", best.loss.total)
        println(io, "lateral_conductivity_scale=", best.lateral)
        println(io, "edge_felt_contact_scale=", best.edge_felt)
        println(io, "beam_sigma_mm=", 1000.0 * best.beam_sigma)
        println(io, "v12_heating_parameters_frozen=true")
    end
    println("v14 network selected = ", best)
    return best
end

if abspath(PROGRAM_FILE) == @__FILE__
    fit_network_stage_2D_v14()
end

