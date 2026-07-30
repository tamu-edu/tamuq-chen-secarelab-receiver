# ============================================================================
# calibrate_2D_v15_staged.jl
#
# Stage A screens the radial mechanism on three representative heating cases.
# Stage B evaluates the four best structures on all training and cooling cases.
# All v14 hydraulic, optical, heat-transfer and hardware parameters are frozen.
# ============================================================================

using Statistics

include("run_2D_v15.jl")

const SKIN_THICKNESS_GRID_2D_v15 = (0.05e-3, 0.15e-3, 0.30e-3)
const SKIN_CONDUCTANCE_GRID_2D_v15 = (0.01, 0.03, 0.10, 0.30, 1.00)
const MECHANISM_IDS_2D_v15 = ("E67", "E72", "E77")

function structural_loss_2D_v15(heating_cases, cooling_cases=Any[])
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
    cooling = isempty(cooling_cases) ? 0.0 : mean(
        mean(abs2, (case.model .- case.observed) ./ 50.0)
        for case in cooling_cases
    )
    signs_mid = count(
        sign(case.model[end, 2] - case.model[end, 4]) ==
        sign(case.observed[end, 2] - case.observed[end, 4])
        for case in heating_cases
    )
    signs_deep = count(
        sign(case.model[end, 3] - case.model[end, 5]) ==
        sign(case.observed[end, 3] - case.observed[end, 5])
        for case in heating_cases
    )
    total = transient + 0.30 * steady + 2.00 * radial +
            0.30 * axial + 0.50 * cooling
    return (
        total=total, transient=transient, steady=steady,
        radial=radial, axial=axial, cooling=cooling,
        signs_mid=signs_mid, signs_deep=signs_deep,
    )
end

function _simulate_set_2D_v15(inputs, p)
    cases = Vector{Any}(undef, length(inputs))
    Threads.@threads for index in eachindex(inputs)
        cases[index] = simulate_case_2D_v15(
            inputs[index], p;
            reltol=3e-4, abstol=1e-4, dtmax=120.0,
        )
    end
    return cases
end

function fit_skin_stage_2D_v15()
    mkpath(OUTPUT_DIR_2D_v15)
    mechanism_inputs = [
        case_inputs_2D_v15(id; max_points=41)
        for id in MECHANISM_IDS_2D_v15
    ]
    screen_rows = NamedTuple[]
    candidates = NamedTuple[]
    evaluation = 0
    total = length(SKIN_THICKNESS_GRID_2D_v15) *
            length(SKIN_CONDUCTANCE_GRID_2D_v15)
    for thickness in SKIN_THICKNESS_GRID_2D_v15,
        conductance in SKIN_CONDUCTANCE_GRID_2D_v15
        evaluation += 1
        p = selected_parameters_2D_v15(
            nominal_mesh=false, screen_mesh=true,
            skin_thickness=thickness,
            skin_conductance_scale=conductance,
        )
        cases = _simulate_set_2D_v15(mechanism_inputs, p)
        loss = structural_loss_2D_v15(cases)
        candidate = (
            thickness=thickness, conductance=conductance,
            mechanism_loss=loss, parameters=p,
        )
        push!(candidates, candidate)
        push!(screen_rows, (
            evaluation=evaluation,
            objective=loss.total,
            transient_loss=loss.transient,
            steady_loss=loss.steady,
            radial_loss=loss.radial,
            axial_loss=loss.axial,
            mid_sign_correct=loss.signs_mid,
            deep_sign_correct=loss.signs_deep,
            skin_thickness_mm=1000.0 * thickness,
            channel_conductance_scale=conductance,
        ))
        println(
            "v15 mechanism $evaluation/$total " *
            "loss=$(round(loss.total,digits=5)) " *
            "signs=$(loss.signs_mid)/$(loss.signs_deep) " *
            "delta=$(1000thickness)mm Gscale=$conductance",
        )
        flush(stdout)
    end
    sort!(candidates; by=candidate -> candidate.mechanism_loss.total)
    shortlist = candidates[1:min(4, length(candidates))]

    heating_inputs = [
        case_inputs_2D_v15(id; max_points=61)
        for id in TRAIN_HEATING_2D_v15
    ]
    cooling_inputs = [
        case_inputs_2D_v15(
            id; cooling=true, max_points=61,
        )
        for id in COOLING_IDS_2D_v15
    ]
    final_rows = NamedTuple[]
    best = nothing
    for (index, candidate) in enumerate(shortlist)
        p = selected_parameters_2D_v15(
            nominal_mesh=false, screen_mesh=true,
            skin_thickness=candidate.thickness,
            skin_conductance_scale=candidate.conductance,
        )
        heating = _simulate_set_2D_v15(heating_inputs, p)
        cooling = _simulate_set_2D_v15(cooling_inputs, p)
        loss = structural_loss_2D_v15(heating, cooling)
        result = (
            thickness=candidate.thickness,
            conductance=candidate.conductance,
            loss=loss, parameters=p,
        )
        (best === nothing || loss.total < best.loss.total) &&
            (best = result)
        push!(final_rows, (
            rank=index, objective=loss.total,
            transient_loss=loss.transient,
            steady_loss=loss.steady,
            radial_loss=loss.radial,
            axial_loss=loss.axial,
            cooling_loss=loss.cooling,
            mid_sign_correct=loss.signs_mid,
            deep_sign_correct=loss.signs_deep,
            skin_thickness_mm=1000.0 * candidate.thickness,
            channel_conductance_scale=candidate.conductance,
        ))
        println(
            "v15 full $index/$(length(shortlist)) " *
            "loss=$(round(loss.total,digits=5)) " *
            "signs=$(loss.signs_mid)/$(loss.signs_deep)",
        )
        flush(stdout)
    end
    _write_namedtuples_csv_2D_v15(
        joinpath(
            OUTPUT_DIR_2D_v15,
            "skin_mechanism_screen_2D_v15.csv",
        ), screen_rows,
    )
    _write_namedtuples_csv_2D_v15(
        joinpath(
            OUTPUT_DIR_2D_v15,
            "skin_shortlist_calibration_2D_v15.csv",
        ), final_rows,
    )
    open(joinpath(
        OUTPUT_DIR_2D_v15,
        "parameters_skin_selected_2D_v15.txt",
    ), "w") do io
        println(io, "objective=", best.loss.total)
        println(io, "skin_thickness_mm=", 1000.0 * best.thickness)
        println(
            io, "channel_conductance_scale=", best.conductance,
        )
        println(io, "mid_sign_correct_training=",
                best.loss.signs_mid, "/9")
        println(io, "deep_sign_correct_training=",
                best.loss.signs_deep, "/9")
        println(io, "all_v14_parameters_frozen=true")
    end
    println("v15 selected = ", best)
    return best
end

if abspath(PROGRAM_FILE) == @__FILE__
    fit_skin_stage_2D_v15()
end
