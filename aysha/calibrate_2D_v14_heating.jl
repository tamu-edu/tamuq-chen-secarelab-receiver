# ============================================================================
# calibrate_2D_v14_heating.jl
#
# Independent heating refit across error-optimal and profile-producing square
# network candidates.  Only shared beta, Graetz scale, and group powers move.
# ============================================================================

include("calibrate_2D_v14_staged.jl")

const NETWORK_HEATING_CANDIDATES_2D_v14 = (
    (1.00, 3.00, 14e-3),
    (1.00, 1.00, 14e-3),
    (1.00, 1.00, 100e-3),
    (1.00, 3.00, 100e-3),
    (0.20, 0.30, 100e-3),
    (0.20, 1.00, 100e-3),
    (0.05, 0.30, 30e-3),
    (0.05, 1.00, 100e-3),
)
const BETA_GRID_2D_v14 = (50.0, 100.0, 150.0, 250.0, 400.0)
const HEAT_GRID_2D_v14 = (0.75, 1.20, 1.80)
const POWER_MULTIPLIERS_2D_v14 = (0.85, 1.00, 1.15)

function fit_heating_stage_2D_v14()
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
    total_shape = length(NETWORK_HEATING_CANDIDATES_2D_v14) *
                  length(BETA_GRID_2D_v14) *
                  length(HEAT_GRID_2D_v14)
    for config in NETWORK_HEATING_CANDIDATES_2D_v14
        lateral, felt, sigma = config
        network = rebuild_network_2D_v14(
            base; lateral_scale=lateral,
            edge_felt_scale=felt,
        )
        shape_base = rebuild_heating_2D_v14(
            network; beam_sigma=sigma,
        )
        cooling_cases = _simulate_set_2D_v14(
            cooling_inputs, shape_base,
        )
        for beta in BETA_GRID_2D_v14,
            heat_scale in HEAT_GRID_2D_v14
            evaluation += 1
            p = rebuild_heating_2D_v14(
                shape_base; beta_opt=beta,
                beam_sigma=sigma, heat_scale=heat_scale,
            )
            heating_cases = _simulate_set_2D_v14(
                heating_inputs, p,
            )
            loss = network_loss_2D_v14(
                heating_cases, cooling_cases,
            )
            candidate = (
                parameters=p, lateral=lateral,
                edge_felt=felt, beam_sigma=sigma,
                beta=beta, heat_scale=heat_scale,
                power_scales=(
                    p.optics.scale_456, p.optics.scale_304,
                    p.optics.scale_256,
                ), loss=loss,
            )
            if best === nothing || loss.total < best.loss.total
                best = candidate
            end
            push!(trace, (
                stage="shape", evaluation=evaluation,
                objective=loss.total,
                transient_loss=loss.transient,
                steady_loss=loss.steady,
                radial_loss=loss.radial,
                axial_loss=loss.axial,
                cooling_loss=loss.cooling,
                lateral_conductivity_scale=lateral,
                edge_felt_contact_scale=felt,
                beam_sigma_mm=1000.0 * sigma,
                beta_opt_1_m=beta,
                heat_transfer_scale=heat_scale,
                power_scale_456=p.optics.scale_456,
                power_scale_304=p.optics.scale_304,
                power_scale_256=p.optics.scale_256,
            ))
            println(
                "v14 shape $evaluation/$total_shape " *
                "loss=$(round(loss.total,digits=5)) " *
                "config=$config beta=$beta Sh=$heat_scale",
            )
            flush(stdout)
        end
    end

    network = rebuild_network_2D_v14(
        base; lateral_scale=best.lateral,
        edge_felt_scale=best.edge_felt,
    )
    cooling_base = rebuild_heating_2D_v14(
        network; beta_opt=best.beta,
        beam_sigma=best.beam_sigma,
        heat_scale=best.heat_scale,
    )
    cooling_cases = _simulate_set_2D_v14(
        cooling_inputs, cooling_base,
    )
    base_scales = (1.25, 1.05, 0.77)
    best_power = best
    power_index = 0
    for m456 in POWER_MULTIPLIERS_2D_v14,
        m304 in POWER_MULTIPLIERS_2D_v14,
        m256 in POWER_MULTIPLIERS_2D_v14
        power_index += 1
        scales = (
            base_scales[1] * m456,
            base_scales[2] * m304,
            base_scales[3] * m256,
        )
        p = rebuild_heating_2D_v14(
            network; beta_opt=best.beta,
            beam_sigma=best.beam_sigma,
            heat_scale=best.heat_scale,
            power_scales=scales,
        )
        heating_cases = _simulate_set_2D_v14(
            heating_inputs, p,
        )
        loss = network_loss_2D_v14(
            heating_cases, cooling_cases,
        )
        candidate = (
            parameters=p, lateral=best.lateral,
            edge_felt=best.edge_felt,
            beam_sigma=best.beam_sigma, beta=best.beta,
            heat_scale=best.heat_scale,
            power_scales=scales, loss=loss,
        )
        loss.total < best_power.loss.total &&
            (best_power = candidate)
        push!(trace, (
            stage="power", evaluation=power_index,
            objective=loss.total,
            transient_loss=loss.transient,
            steady_loss=loss.steady,
            radial_loss=loss.radial,
            axial_loss=loss.axial,
            cooling_loss=loss.cooling,
            lateral_conductivity_scale=best.lateral,
            edge_felt_contact_scale=best.edge_felt,
            beam_sigma_mm=1000.0 * best.beam_sigma,
            beta_opt_1_m=best.beta,
            heat_transfer_scale=best.heat_scale,
            power_scale_456=scales[1],
            power_scale_304=scales[2],
            power_scale_256=scales[3],
        ))
        println(
            "v14 power $power_index/27 " *
            "loss=$(round(loss.total,digits=5)) scales=$scales",
        )
        flush(stdout)
    end
    _write_namedtuples_csv_2D_v14(
        joinpath(
            OUTPUT_DIR_2D_v14,
            "heating_calibration_trace_2D_v14.csv",
        ), trace,
    )
    open(joinpath(
        OUTPUT_DIR_2D_v14,
        "parameters_heating_selected_2D_v14.txt",
    ), "w") do io
        println(io, "objective=", best_power.loss.total)
        println(io, "lateral_conductivity_scale=",
                best_power.lateral)
        println(io, "edge_felt_contact_scale=",
                best_power.edge_felt)
        println(io, "beam_sigma_mm=",
                1000.0 * best_power.beam_sigma)
        println(io, "beta_opt_1_m=", best_power.beta)
        println(io, "heat_transfer_scale=",
                best_power.heat_scale)
        println(io, "power_scale_456=",
                best_power.power_scales[1])
        println(io, "power_scale_304=",
                best_power.power_scales[2])
        println(io, "power_scale_256=",
                best_power.power_scales[3])
    end
    println("v14 heating selected = ", best_power)
    return best_power
end

if abspath(PROGRAM_FILE) == @__FILE__
    fit_heating_stage_2D_v14()
end
