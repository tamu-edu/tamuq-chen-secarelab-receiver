# ============================================================================
# refine_2D_v16_power.jl - local nominal refinement of group power minima
# ============================================================================

include("calibrate_2D_v16_power.jl")

const POWER_REFINE_GROUPS_2D_v16 = (
    (
        name="456", train=("E67", "E69", "E71"),
        valid=("E68", "E70"),
        candidates=(1.65, 1.75, 1.85, 1.95),
    ),
    (
        name="304", train=("E72", "E74", "E76"),
        valid=("E73", "E75"),
        candidates=(1.70, 1.80, 1.90, 2.00, 2.10),
    ),
    (
        name="256", train=("E77", "E79", "E81"),
        valid=("E78", "E80"),
        candidates=(0.85, 0.95, 1.05, 1.15, 1.25),
    ),
)

function refine_power_stage_2D_v16()
    cooling = _read_nominal_cooling_2D_v16()
    selected_scales = [1.8, 1.9, 1.1]
    rows = NamedTuple[]
    group_rows = NamedTuple[]
    for (group_index, group) in enumerate(
        POWER_REFINE_GROUPS_2D_v16,
    )
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
            loss = power_level_loss_2D_v16(
                _simulate_heating_set_2D_v16(
                    training_inputs, p,
                ),
            )
            result = (scale=candidate, loss=loss)
            (best === nothing || loss.total < best.loss.total) &&
                (best = result)
            push!(rows, (
                irradiance_group=group.name,
                evaluation=evaluation,
                objective=loss.total,
                receiver_rmse_K=loss.receiver_rmse_K,
                power_scale=candidate,
            ))
            println(
                "v16 power refine $(group.name) $evaluation/" *
                "$(length(group.candidates)) " *
                "loss=$(round(loss.total,digits=5)) " *
                "scale=$candidate",
            )
            flush(stdout)
        end
        selected_scales[group_index] = best.scale
        validation_inputs = [
            case_inputs_2D_v16(id; max_points=61)
            for id in group.valid
        ]
        pbest = parameters_2D_v16(
            nominal_mesh=true,
            skin_felt_contact_scale=cooling.contact,
            felt_conductivity_scale=cooling.felt_k,
            felt_heat_capacity_scale=cooling.felt_cp,
            power_scales=Tuple(selected_scales),
        )
        valid_loss = power_level_loss_2D_v16(
            _simulate_heating_set_2D_v16(
                validation_inputs, pbest,
            ),
        )
        push!(group_rows, (
            irradiance_group=group.name,
            selected_power_scale=best.scale,
            training_objective=best.loss.total,
            training_receiver_rmse_K=best.loss.receiver_rmse_K,
            heldout_objective=valid_loss.total,
            heldout_receiver_rmse_K=valid_loss.receiver_rmse_K,
        ))
        println(
            "v16 refined $(group.name)=$(best.scale), " *
            "heldout=$(round(valid_loss.receiver_rmse_K,digits=1))K",
        )
    end
    _write_namedtuples_csv_2D_v16(
        joinpath(
            OUTPUT_DIR_2D_v16,
            "power_refinement_trace_2D_v16.csv",
        ), rows,
    )
    _write_namedtuples_csv_2D_v16(
        joinpath(
            OUTPUT_DIR_2D_v16,
            "power_refined_group_validation_2D_v16.csv",
        ), group_rows,
    )
    open(joinpath(
        OUTPUT_DIR_2D_v16,
        "parameters_power_selected_2D_v16.txt",
    ), "w") do io
        println(io, "power_scale_456=", selected_scales[1])
        println(io, "power_scale_304=", selected_scales[2])
        println(io, "power_scale_256=", selected_scales[3])
        println(io, "selected_on_nominal_mesh=true")
        println(io, "locally_refined=true")
        println(io, "cooling_parameters_frozen=true")
    end
    println("v16 refined power scales = ", Tuple(selected_scales))
    return (
        power_scales=Tuple(selected_scales),
        groups=group_rows,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    refine_power_stage_2D_v16()
end
