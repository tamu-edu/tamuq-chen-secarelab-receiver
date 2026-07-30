# ============================================================================
# calibrate_2D_v14_hydraulics.jl
#
# Cold/t0-only square-network identification.  The physical DP1 offset is
# retained; channel resistance and rear-groove strength are re-identified
# because annular v11 overlap weights do not transfer to square channels.
# ============================================================================

using Statistics

include("run_2D_v14.jl")

function t0_hydraulic_inputs_2D_v14()
    return [
        case_inputs_2D_v14(id; max_points=nothing)
        for id in HEATING_IDS_2D_v14
    ]
end

function t0_dp_prediction_2D_v14(inputs, p)
    grid = build_network_grid2D(p)
    bg = grid.base_grid
    t0 = Dict(
        sensor => Float64(observation(
            inputs.data, inputs.id, "_$(sensor)",
        )[1])
        for sensor in SENSOR_NAMES_2D_v14
    )
    base_initial = build_initial_state_2D(
        bg, p.base.base, t0, inputs.ambient[1],
    )
    state = Receiver2D_v14._initial_network_state2D(
        base_initial, p, grid,
    )
    layout = Receiver2D_v14._state_layout(grid)
    Tch = reshape(
        view(state, layout.channel), grid.group_count, bg.nz,
    )
    Ttube = view(state, layout.tube)
    op = OperatingCondition2D(
        irradiance=0.0,
        flow_lpm=inputs.flow[1],
        inlet_temperature=inputs.inlet[1],
        ambient_temperature=inputs.ambient[1],
    )
    work = Receiver2D_v14.NetworkWorkspace2D(grid)
    Receiver2D_v14._gas_profile_network2D!(
        work, Tch, Ttube, inputs.times[1], p, op, grid,
    )
    return p.hydraulics.dp1_zero_offset_mbar +
           work.common_dp / 100.0
end

function hydraulic_objective_2D_v14(inputs, p)
    predictions = [
        t0_dp_prediction_2D_v14(input, p) for input in inputs
    ]
    observed = [input.dp1[1] for input in inputs]
    return (
        rmse=sqrt(mean(abs2, predictions .- observed)),
        predictions=predictions, observed=observed,
    )
end

function fit_hydraulics_2D_v14()
    mkpath(OUTPUT_DIR_2D_v14)
    base = inherited_parameters_2D_v14(
        nominal_mesh=false, screen_mesh=true,
    )
    inputs = t0_hydraulic_inputs_2D_v14()
    coarse_scales = (0.20, 0.50, 0.80, 1.10, 1.50, 2.00)
    coarse_grooves = (0.0, 25.0, 50.0, 100.0, 184.0, 250.0, 350.0)
    trace = NamedTuple[]
    best = nothing
    evaluation = 0
    for scale in coarse_scales, groove in coarse_grooves
        evaluation += 1
        p = rebuild_hydraulics_2D_v14(
            base; resistance_scale=scale,
            groove_coefficient=groove,
        )
        fit = hydraulic_objective_2D_v14(inputs, p)
        candidate = (
            parameters=p, scale=scale, groove=groove, fit=fit,
        )
        if best === nothing || fit.rmse < best.fit.rmse
            best = candidate
        end
        push!(trace, (
            stage="coarse", evaluation=evaluation,
            rmse_mbar=fit.rmse,
            hydraulic_resistance_scale=scale,
            groove_loss_coefficient=groove,
        ))
    end
    scale_grid = range(
        max(0.05, best.scale - 0.35),
        best.scale + 0.35; length=15,
    )
    groove_grid = range(
        max(0.0, best.groove - 75.0),
        best.groove + 75.0; length=16,
    )
    for scale in scale_grid, groove in groove_grid
        evaluation += 1
        p = rebuild_hydraulics_2D_v14(
            base; resistance_scale=scale,
            groove_coefficient=groove,
        )
        fit = hydraulic_objective_2D_v14(inputs, p)
        candidate = (
            parameters=p, scale=scale, groove=groove, fit=fit,
        )
        if fit.rmse < best.fit.rmse
            best = candidate
        end
        push!(trace, (
            stage="refinement", evaluation=evaluation,
            rmse_mbar=fit.rmse,
            hydraulic_resistance_scale=scale,
            groove_loss_coefficient=groove,
        ))
    end

    # Nested no-groove comparison.
    no_groove = nothing
    for scale in range(0.05, 3.0; length=120)
        p = rebuild_hydraulics_2D_v14(
            base; resistance_scale=scale,
            groove_coefficient=0.0,
        )
        fit = hydraulic_objective_2D_v14(inputs, p)
        candidate = (
            parameters=p, scale=scale, groove=0.0, fit=fit,
        )
        if no_groove === nothing ||
           fit.rmse < no_groove.fit.rmse
            no_groove = candidate
        end
    end

    rows = NamedTuple[]
    for (index, input) in enumerate(inputs)
        push!(rows, (
            simulation_id=input.id,
            flow_lpm=input.flow[1],
            dp1_observed_mbar=best.fit.observed[index],
            dp1_square_network_mbar=
                best.fit.predictions[index],
            residual_mbar=
                best.fit.predictions[index] -
                best.fit.observed[index],
        ))
    end
    _write_namedtuples_csv_2D_v14(
        joinpath(
            OUTPUT_DIR_2D_v14,
            "hydraulic_t0_calibration_trace_2D_v14.csv",
        ), trace,
    )
    _write_namedtuples_csv_2D_v14(
        joinpath(
            OUTPUT_DIR_2D_v14,
            "hydraulic_t0_predictions_2D_v14.csv",
        ), rows,
    )
    open(joinpath(
        OUTPUT_DIR_2D_v14,
        "parameters_hydraulic_selected_2D_v14.txt",
    ), "w") do io
        println(io, "t0_rmse_mbar=", best.fit.rmse)
        println(io, "hydraulic_resistance_scale=", best.scale)
        println(io, "groove_loss_coefficient=", best.groove)
        println(io, "dp1_zero_offset_mbar=",
                best.parameters.hydraulics.dp1_zero_offset_mbar)
        println(io, "no_groove_best_rmse_mbar=",
                no_groove.fit.rmse)
        println(io, "no_groove_resistance_scale=",
                no_groove.scale)
    end
    println("v14 hydraulic selected = ", (
        rmse=best.fit.rmse, scale=best.scale,
        groove=best.groove,
        no_groove_rmse=no_groove.fit.rmse,
        no_groove_scale=no_groove.scale,
    ))
    return (
        parameters=best.parameters, scale=best.scale,
        groove=best.groove, fit=best.fit,
        no_groove=no_groove,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    fit_hydraulics_2D_v14()
end

