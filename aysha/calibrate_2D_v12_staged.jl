# ============================================================================
# calibrate_2D_v12_staged.jl
#
# Stage 1A: bounded physical felt/contact candidates using C69/C80/C81 only.
# Stage 1B: shared thermocouple observation parameters, evaluated without
#           rerunning the thermal model.
# Stage 2 is deliberately separate and is entered only after cooling selects
# a credible assembly point.
# ============================================================================

using Statistics

include("run_2D_v12.jl")

const COOLING_FIT_POINTS_2D_v12 = 81
const PHYSICAL_CANDIDATES_2D_v12 = (
    (1.00, 1.00, 25.0),
    (0.60, 1.00, 25.0),
    (1.40, 1.00, 25.0),
    (1.00, 0.70, 25.0),
    (1.00, 1.30, 25.0),
    (1.00, 1.00, 5.0),
    (1.00, 1.00, 80.0),
    (0.60, 0.70, 5.0),
    (0.60, 0.70, 80.0),
    (0.60, 1.30, 5.0),
    (0.60, 1.30, 80.0),
    (1.40, 0.70, 5.0),
    (1.40, 0.70, 80.0),
    (1.40, 1.30, 5.0),
    (1.40, 1.30, 80.0),
)

const PHYSICAL_REFINEMENT_2D_v12 = (
    (1.00, 1.30, 25.0),
    (0.80, 1.30, 25.0),
    (1.20, 1.30, 25.0),
    (1.00, 1.10, 25.0),
    (1.00, 1.50, 25.0),
    (1.00, 1.30, 15.0),
    (1.00, 1.30, 40.0),
    (0.80, 1.50, 15.0),
    (1.20, 1.50, 15.0),
    (1.00, 1.50, 40.0),
)

const SIDE_TAU_GRID_2D_v12 = (0.0, 30.0, 90.0, 180.0)
const INTERIOR_TAU_GRID_2D_v12 = (0.0, 10.0, 30.0, 60.0)
const OUTLET_TAU_GRID_2D_v12 = (0.0, 30.0, 90.0, 180.0)
const WALL_FRACTION_GRID_2D_v12 = (0.20, 0.50, 0.80)

function observation_case_2D_v12(case, p_observation)
    result = SimulationResult2D(
        case.result.base_result,
        case.result.adaptor_temperature,
        case.result.housing_extension_temperature,
        p_observation,
    )
    predictions = sensor_predictions2D(result)
    model = hcat((
        getproperty(predictions, sensor)
        for sensor in SENSOR_NAMES_2D_v12
    )...)
    return (
        inputs=case.inputs, result=result, predictions=predictions,
        model=model, observed=case.observed,
    )
end

function cooling_residual_loss_2D_v12(cases)
    time_loss = 0.0
    steady_loss = 0.0
    t90_loss = 0.0
    count = 0
    for case in cases
        for index in eachindex(SENSOR_NAMES_2D_v12)
            model = case.model[:, index]
            observed = case.observed[:, index]
            time_loss += mean(abs2, (model .- observed) ./ 100.0)
            steady_loss += ((model[end] - observed[end]) / 75.0)^2
            t90_loss += (
                (
                    get_t90_2D(case.inputs.times, model) -
                    get_t90_2D(case.inputs.times, observed)
                ) / 1000.0
            )^2
            count += 1
        end
    end
    return (
        total=(time_loss + 0.30 * steady_loss + 0.35 * t90_loss) /
              max(count, 1),
        transient=time_loss / max(count, 1),
        steady=steady_loss / max(count, 1),
        t90=t90_loss / max(count, 1),
    )
end

function best_observation_for_physics_2D_v12(raw_cases, p_physical)
    best = nothing
    for side_tau in SIDE_TAU_GRID_2D_v12,
        interior_tau in INTERIOR_TAU_GRID_2D_v12,
        outlet_tau in OUTLET_TAU_GRID_2D_v12,
        wall_fraction in WALL_FRACTION_GRID_2D_v12

        candidate = rebuild_parameters_2D_v12(
            p_physical;
            side_tau=side_tau,
            interior_tau=interior_tau,
            outlet_tau=outlet_tau,
            interior_wall_fraction=wall_fraction,
        )
        observed_cases = [
            observation_case_2D_v12(case, candidate)
            for case in raw_cases
        ]
        loss = cooling_residual_loss_2D_v12(observed_cases)
        if best === nothing || loss.total < best.loss.total
            best = (
                parameters=candidate,
                cases=observed_cases,
                loss=loss,
                observation=(
                    side_tau, interior_tau, outlet_tau, wall_fraction,
                ),
            )
        end
    end
    return best
end

function fit_cooling_stage_2D_v12(;
    physical_candidates=PHYSICAL_CANDIDATES_2D_v12,
    output_tag="cooling",
)
    mkpath(OUTPUT_DIR_2D_v12)
    base = base_parameters_2D_v12(; nominal_mesh=false)
    inputs = [
        case_inputs_2D_v12(
            id; cooling=true,
            max_points=COOLING_FIT_POINTS_2D_v12,
        )
        for id in COOLING_IDS_2D_v12
    ]
    trace = NamedTuple[]
    candidates = Vector{Any}(undef, length(physical_candidates))
    for (candidate_index, theta) in
        enumerate(physical_candidates)

        felt_k, felt_cp, contact_h = theta
        physical = rebuild_parameters_2D_v12(
            base;
            felt_k_scale=felt_k,
            felt_cp_scale=felt_cp,
            receiver_adaptor_h=contact_h,
            side_tau=0.0,
            interior_tau=0.0,
            outlet_tau=0.0,
            interior_wall_fraction=0.5,
        )
        raw_cases = Vector{Any}(undef, length(inputs))
        Threads.@threads for index in eachindex(inputs)
            raw_cases[index] = simulate_case_2D_v12(
                inputs[index], physical,
            )
        end
        selected = best_observation_for_physics_2D_v12(
            raw_cases, physical,
        )
        candidates[candidate_index] = selected
        push!(trace, (
            evaluation=candidate_index,
            objective=selected.loss.total,
            transient_loss=selected.loss.transient,
            steady_loss=selected.loss.steady,
            t90_loss=selected.loss.t90,
            felt_k_scale=felt_k,
            felt_cp_scale=felt_cp,
            receiver_adaptor_h_W_m2K=contact_h,
            side_tau_s=selected.observation[1],
            interior_tau_s=selected.observation[2],
            outlet_tau_s=selected.observation[3],
            interior_wall_fraction=selected.observation[4],
        ))
        println(
            "cooling candidate $candidate_index/" *
            "$(length(physical_candidates)): " *
            "loss=$(round(selected.loss.total,digits=6)) " *
            "physical=$theta observation=$(selected.observation)",
        )
        flush(stdout)
    end
    best_index = argmin(candidate.loss.total for candidate in candidates)
    best = candidates[best_index]
    _write_namedtuples_csv_2D_v12(
        joinpath(
            OUTPUT_DIR_2D_v12,
            "$(output_tag)_calibration_trace_2D_v12.csv",
        ),
        trace,
    )
    p = best.parameters
    open(joinpath(
        OUTPUT_DIR_2D_v12,
        "parameters_$(output_tag)_selected_2D_v12.txt",
    ), "w") do io
        println(io, "objective=", best.loss.total)
        println(io, "felt_conductivity_scale=",
            p.assembly.felt_conductivity_scale)
        println(io, "felt_heat_capacity_scale=",
            p.assembly.felt_heat_capacity_scale)
        println(io, "receiver_adaptor_contact_h_W_m2K=",
            p.assembly.receiver_adaptor_contact_h)
        println(io, "side_time_constant_s=",
            p.observation.side_time_constant_s)
        println(io, "interior_time_constant_s=",
            p.observation.interior_time_constant_s)
        println(io, "outlet_time_constant_s=",
            p.observation.outlet_time_constant_s)
        println(io, "interior_wall_fraction=",
            p.observation.interior_wall_fraction)
    end
    return (
        parameters=best.parameters, cases=best.cases,
        loss=best.loss, trace=trace, selected_index=best_index,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    result = fit_cooling_stage_2D_v12()
    println("selected cooling index = ", result.selected_index)
    println("selected cooling loss = ", result.loss)
end
