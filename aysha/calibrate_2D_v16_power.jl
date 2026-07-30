# ============================================================================
# calibrate_2D_v16_power.jl
#
# Corrective nominal-mesh power fit after cooling-path identification.
# Spatial, hydraulic, heat-transfer, and cooling parameters remain frozen.
# ============================================================================

using Statistics

include("run_2D_v16.jl")

const POWER_GROUPS_2D_v16 = (
    (
        name="456", train=("E67", "E69", "E71"),
        valid=("E68", "E70"),
        candidates=(1.50, 1.80, 2.10, 2.40, 2.70),
    ),
    (
        name="304", train=("E72", "E74", "E76"),
        valid=("E73", "E75"),
        candidates=(1.50, 1.90, 2.30, 2.70, 3.10),
    ),
    (
        name="256", train=("E77", "E79", "E81"),
        valid=("E78", "E80"),
        candidates=(1.10, 1.40, 1.70, 2.00, 2.30),
    ),
)

function _read_nominal_cooling_2D_v16()
    values = Dict{String,Float64}()
    path = joinpath(
        OUTPUT_DIR_2D_v16,
        "parameters_cooling_nominal_selected_2D_v16.txt",
    )
    for line in eachline(path)
        occursin("=", line) || continue
        key, raw = split(line, "="; limit=2)
        value = tryparse(Float64, raw)
        value === nothing || (values[key] = value)
    end
    return (
        contact=values["skin_felt_contact_scale"],
        felt_k=values["felt_conductivity_scale"],
        felt_cp=values["felt_heat_capacity_scale"],
    )
end

function power_level_loss_2D_v16(cases)
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

function _simulate_heating_set_2D_v16(inputs, p)
    cases = Vector{Any}(undef, length(inputs))
    Threads.@threads for index in eachindex(inputs)
        cases[index] = simulate_case_2D_v16(
            inputs[index], p;
            reltol=3e-4, abstol=1e-4, dtmax=120.0,
        )
    end
    return cases
end

function fit_power_stage_2D_v16()
    mkpath(OUTPUT_DIR_2D_v16)
    cooling = _read_nominal_cooling_2D_v16()
    selected_scales = [1.25, 1.2075, 0.8855]
    rows = NamedTuple[]
    group_results = NamedTuple[]
    for (group_index, group) in enumerate(POWER_GROUPS_2D_v16)
        training_inputs = [
            case_inputs_2D_v16(id; max_points=61)
            for id in group.train
        ]
        best = nothing
        for (evaluation, candidate) in enumerate(group.candidates)
            scales = copy(selected_scales)
            scales[group_index] = candidate
            p = parameters_2D_v16(
                nominal_mesh=true,
                skin_felt_contact_scale=cooling.contact,
                felt_conductivity_scale=cooling.felt_k,
                felt_heat_capacity_scale=cooling.felt_cp,
                power_scales=Tuple(scales),
            )
            cases = _simulate_heating_set_2D_v16(
                training_inputs, p,
            )
            loss = power_level_loss_2D_v16(cases)
            result = (
                scale=candidate, loss=loss, parameters=p,
            )
            (best === nothing || loss.total < best.loss.total) &&
                (best = result)
            push!(rows, (
                irradiance_group=group.name,
                evaluation=evaluation,
                objective=loss.total,
                receiver_rmse_K=loss.receiver_rmse_K,
                receiver_transient_loss=loss.receiver_transient,
                receiver_final_loss=loss.receiver_final,
                felt_guard_loss=loss.felt_guard,
                power_scale=candidate,
            ))
            println(
                "v16 power $(group.name) $evaluation/" *
                "$(length(group.candidates)) " *
                "loss=$(round(loss.total,digits=5)) " *
                "rmse=$(round(loss.receiver_rmse_K,digits=1))K " *
                "scale=$candidate",
            )
            flush(stdout)
        end
        selected_scales[group_index] = best.scale
        validation_inputs = [
            case_inputs_2D_v16(id; max_points=61)
            for id in group.valid
        ]
        final_parameters = parameters_2D_v16(
            nominal_mesh=true,
            skin_felt_contact_scale=cooling.contact,
            felt_conductivity_scale=cooling.felt_k,
            felt_heat_capacity_scale=cooling.felt_cp,
            power_scales=Tuple(selected_scales),
        )
        validation_loss = power_level_loss_2D_v16(
            _simulate_heating_set_2D_v16(
                validation_inputs, final_parameters,
            ),
        )
        push!(group_results, (
            irradiance_group=group.name,
            selected_power_scale=best.scale,
            training_objective=best.loss.total,
            training_receiver_rmse_K=best.loss.receiver_rmse_K,
            heldout_objective=validation_loss.total,
            heldout_receiver_rmse_K=
                validation_loss.receiver_rmse_K,
        ))
        println(
            "v16 power $(group.name) selected=$(best.scale), " *
            "heldout RMSE=$(round(validation_loss.receiver_rmse_K,digits=1))K",
        )
    end
    _write_namedtuples_csv_2D_v16(
        joinpath(
            OUTPUT_DIR_2D_v16,
            "power_calibration_trace_2D_v16.csv",
        ), rows,
    )
    _write_namedtuples_csv_2D_v16(
        joinpath(
            OUTPUT_DIR_2D_v16,
            "power_group_validation_2D_v16.csv",
        ), group_results,
    )
    open(joinpath(
        OUTPUT_DIR_2D_v16,
        "parameters_power_selected_2D_v16.txt",
    ), "w") do io
        println(io, "power_scale_456=", selected_scales[1])
        println(io, "power_scale_304=", selected_scales[2])
        println(io, "power_scale_256=", selected_scales[3])
        println(io, "selected_on_nominal_mesh=true")
        println(io, "cooling_parameters_frozen=true")
    end
    println("v16 selected power scales = ", Tuple(selected_scales))
    return (
        power_scales=Tuple(selected_scales),
        groups=group_results,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    fit_power_stage_2D_v16()
end
