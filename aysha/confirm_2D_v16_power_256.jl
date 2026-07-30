# ============================================================================
# confirm_2D_v16_power_256.jl - close the remaining 256-group edge
# ============================================================================

include("calibrate_2D_v16_power.jl")

function confirm_power_256_2D_v16()
    cooling = _read_nominal_cooling_2D_v16()
    inputs = [
        case_inputs_2D_v16(id; max_points=61)
        for id in ("E77", "E79", "E81")
    ]
    candidates = (1.20, 1.25, 1.30, 1.35)
    rows = NamedTuple[]
    best = nothing
    for (index, scale) in enumerate(candidates)
        p = parameters_2D_v16(
            nominal_mesh=true,
            skin_felt_contact_scale=cooling.contact,
            felt_conductivity_scale=cooling.felt_k,
            felt_heat_capacity_scale=cooling.felt_cp,
            power_scales=(1.65, 1.80, scale),
        )
        loss = power_level_loss_2D_v16(
            _simulate_heating_set_2D_v16(inputs, p),
        )
        result = (scale=scale, loss=loss)
        (best === nothing || loss.total < best.loss.total) &&
            (best = result)
        push!(rows, (
            evaluation=index, objective=loss.total,
            receiver_rmse_K=loss.receiver_rmse_K,
            power_scale_256=scale,
        ))
        println(
            "v16 256 confirm $index/$(length(candidates)) " *
            "loss=$(round(loss.total,digits=5)) scale=$scale",
        )
        flush(stdout)
    end
    _write_namedtuples_csv_2D_v16(
        joinpath(
            OUTPUT_DIR_2D_v16,
            "power_256_confirmation_2D_v16.csv",
        ), rows,
    )
    open(joinpath(
        OUTPUT_DIR_2D_v16,
        "parameters_power_selected_2D_v16.txt",
    ), "w") do io
        println(io, "power_scale_456=1.65")
        println(io, "power_scale_304=1.8")
        println(io, "power_scale_256=", best.scale)
        println(io, "selected_on_nominal_mesh=true")
        println(io, "locally_refined=true")
        println(io, "cooling_parameters_frozen=true")
    end
    println("v16 confirmed power scales = ", (1.65, 1.8, best.scale))
    return best
end

if abspath(PROGRAM_FILE) == @__FILE__
    confirm_power_256_2D_v16()
end
