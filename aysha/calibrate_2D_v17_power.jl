# ============================================================================
# calibrate_2D_v17_power.jl - group power profile after flange correction
# ============================================================================

using Statistics

include("run_2D_v17.jl")

const POWER_GROUPS_2D_v17 = (
    (
        name="456", train=("E67", "E69", "E71"),
        valid=("E68", "E70"),
        candidates=(0.90, 1.05, 1.20, 1.35, 1.50, 1.65),
        closure_low=1.05, closure_high=1.34,
    ),
    (
        name="304", train=("E72", "E74", "E76"),
        valid=("E73", "E75"),
        candidates=(1.00, 1.20, 1.40, 1.60, 1.80),
        closure_low=1.23, closure_high=1.58,
    ),
    (
        name="256", train=("E77", "E79", "E81"),
        valid=("E78", "E80"),
        candidates=(0.75, 0.90, 1.05, 1.20, 1.35),
        closure_low=0.84, closure_high=1.11,
    ),
)

function _selected_hardware_2D_v17()
    values = Dict{String,Float64}()
    path = joinpath(
        OUTPUT_DIR_2D_v17,
        "parameters_casing_flange_selected_2D_v17.txt",
    )
    for line in eachline(path)
        occursin("=", line) || continue
        key, raw = split(line, "="; limit=2)
        value = tryparse(Float64, raw)
        value === nothing || (values[key] = value)
    end
    return (
        conductance=values["casing_flange_conductance_W_K"],
        felt_cp=values["felt_heat_capacity_scale"],
    )
end

function power_level_loss_2D_v17(cases)
    receiver_transient = mean(
        mean(abs2, (
            case.model[:, index] .- case.observed[:, index]
        ) ./ 100.0)
        for case in cases for index in 1:6
    )
    receiver_final = mean(
        (
            (case.model[end, index] -
             case.observed[end, index]) / 75.0
        )^2
        for case in cases for index in 1:6
    )
    felt_guard = mean(
        mean(abs2, (
            case.model[:, 7] .- case.observed[:, 7]
        ) ./ 30.0)
        for case in cases
    )
    total = receiver_transient + 0.50 * receiver_final +
            0.10 * felt_guard
    rmse = sqrt(mean(
        mean(abs2, (
            case.model[:, 1:6] .- case.observed[:, 1:6]
        )) for case in cases
    ))
    return (
        total=total,
        receiver_transient=receiver_transient,
        receiver_final=receiver_final,
        felt_guard=felt_guard,
        receiver_rmse_K=rmse,
    )
end

function _simulate_heating_set_2D_v17(inputs, p)
    cases = Vector{Any}(undef, length(inputs))
    Threads.@threads for index in eachindex(inputs)
        cases[index] = simulate_case_2D_v17(
            inputs[index], p;
            reltol=3e-4, abstol=1e-4, dtmax=120.0,
        )
    end
    return cases
end

function fit_power_stage_2D_v17()
    mkpath(OUTPUT_DIR_2D_v17)
    hardware = _selected_hardware_2D_v17()
    selected_scales = [1.20, 1.40, 1.05]
    rows = NamedTuple[]
    group_results = NamedTuple[]
    for (group_index, group) in enumerate(POWER_GROUPS_2D_v17)
        training_inputs = [
            case_inputs_2D_v17(id; max_points=61)
            for id in group.train
        ]
        best = nothing
        for (evaluation, candidate) in enumerate(group.candidates)
            scales = copy(selected_scales)
            scales[group_index] = candidate
            p = parameters_2D_v17(
                nominal_mesh=true,
                felt_heat_capacity_scale=hardware.felt_cp,
                casing_flange_conductance_W_K=
                    hardware.conductance,
                power_scales=Tuple(scales),
            )
            cases = _simulate_heating_set_2D_v17(
                training_inputs, p,
            )
            loss = power_level_loss_2D_v17(cases)
            result = (
                scale=candidate, loss=loss, parameters=p,
            )
            (best === nothing || loss.total < best.loss.total) &&
                (best = result)
            inside = group.closure_low <= candidate <=
                     group.closure_high
            push!(rows, (
                irradiance_group=group.name,
                evaluation=evaluation,
                objective=loss.total,
                receiver_rmse_K=loss.receiver_rmse_K,
                power_scale=candidate,
                inside_independent_closure=inside,
            ))
            println(
                "v17 power $(group.name) $evaluation/" *
                "$(length(group.candidates)) " *
                "loss=$(round(loss.total,digits=5)) " *
                "rmse=$(round(loss.receiver_rmse_K,digits=1))K " *
                "scale=$candidate closure=$inside",
            )
            flush(stdout)
        end
        selected_scales[group_index] = best.scale
        validation_inputs = [
            case_inputs_2D_v17(id; max_points=61)
            for id in group.valid
        ]
        final_parameters = parameters_2D_v17(
            nominal_mesh=true,
            felt_heat_capacity_scale=hardware.felt_cp,
            casing_flange_conductance_W_K=
                hardware.conductance,
            power_scales=Tuple(selected_scales),
        )
        validation_loss = power_level_loss_2D_v17(
            _simulate_heating_set_2D_v17(
                validation_inputs, final_parameters,
            ),
        )
        push!(group_results, (
            irradiance_group=group.name,
            selected_power_scale=best.scale,
            inside_independent_closure=(
                group.closure_low <= best.scale <=
                group.closure_high
            ),
            training_objective=best.loss.total,
            training_receiver_rmse_K=
                best.loss.receiver_rmse_K,
            heldout_objective=validation_loss.total,
            heldout_receiver_rmse_K=
                validation_loss.receiver_rmse_K,
        ))
        println(
            "v17 power $(group.name) selected=$(best.scale), " *
            "heldout RMSE=" *
            "$(round(validation_loss.receiver_rmse_K,digits=1))K",
        )
    end
    _write_namedtuples_csv_2D_v17(
        joinpath(
            OUTPUT_DIR_2D_v17,
            "power_calibration_trace_2D_v17.csv",
        ), rows,
    )
    _write_namedtuples_csv_2D_v17(
        joinpath(
            OUTPUT_DIR_2D_v17,
            "power_group_validation_2D_v17.csv",
        ), group_results,
    )
    open(joinpath(
        OUTPUT_DIR_2D_v17,
        "parameters_power_selected_2D_v17.txt",
    ), "w") do io
        println(io, "power_scale_456=", selected_scales[1])
        println(io, "power_scale_304=", selected_scales[2])
        println(io, "power_scale_256=", selected_scales[3])
        println(io, "selected_on_nominal_mesh=true")
        println(io, "hardware_parameters_frozen=true")
    end
    println(
        "v17 selected power scales = ",
        Tuple(selected_scales),
    )
    return (
        power_scales=Tuple(selected_scales),
        groups=group_results,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    fit_power_stage_2D_v17()
end
