# ============================================================================
# calibrate_2D_v13_heating.jl
#
# Independent heating stage.  It prevents a v13 topology from being rejected
# merely because it inherited the optical/Graetz point fitted to v12.
# ============================================================================

include("calibrate_2D_v13_staged.jl")

const POPULATION_HEATING_CANDIDATES_2D_v13 = (
    (0.025, 0.0125, 100.0, 5.0),
    (0.075, 0.0225, 30.0, 5.0),
    (0.15, 0.0225, 2.0, 5.0),
    (0.15, 0.0225, 10.0, 5.0),
    (0.30, 0.045, 2.0, 5.0),
    (0.30, 0.045, 10.0, 5.0),
    (0.45, 0.0675, 2.0, 5.0),
    (0.45, 0.0675, 10.0, 5.0),
    (0.30, 0.0, 0.2, 5.0),
)
const BETA_GRID_2D_v13 = (75.0, 150.0, 250.0, 400.0)
const HEAT_SCALE_GRID_2D_v13 = (0.75, 1.20, 1.80)
const POWER_MULTIPLIERS_2D_v13 = (0.85, 1.00, 1.15)

function simulate_heating_set_2D_v13(inputs, p)
    cases = Vector{Any}(undef, length(inputs))
    Threads.@threads for index in eachindex(inputs)
        cases[index] = simulate_case_2D_v13(
            inputs[index], p;
            reltol=3e-4, abstol=1e-4, dtmax=120.0,
        )
    end
    return cases
end

function simulate_cooling_set_2D_v13(inputs, p)
    cases = Vector{Any}(undef, length(inputs))
    Threads.@threads for index in eachindex(inputs)
        cases[index] = simulate_case_2D_v13(
            inputs[index], p;
            reltol=3e-4, abstol=1e-4, dtmax=120.0,
        )
    end
    return cases
end

function fit_heating_stage_2D_v13()
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
    best = nothing
    evaluation = 0
    total_shape = length(POPULATION_HEATING_CANDIDATES_2D_v13) *
                  length(BETA_GRID_2D_v13) *
                  length(HEAT_SCALE_GRID_2D_v13)
    for theta in POPULATION_HEATING_CANDIDATES_2D_v13
        ppop = rebuild_population_2D_v13(
            base;
            side_fraction=theta[1],
            side_gas_fraction=theta[2],
            active_side_G_per_m=theta[3],
            side_felt_h=theta[4],
        )
        cooling_cases = simulate_cooling_set_2D_v13(
            cooling_inputs, ppop,
        )
        for beta in BETA_GRID_2D_v13,
            heat_scale in HEAT_SCALE_GRID_2D_v13
            evaluation += 1
            p = rebuild_heating_2D_v13(
                ppop; beta_opt=beta, heat_scale=heat_scale,
            )
            heating_cases = simulate_heating_set_2D_v13(
                heating_inputs, p,
            )
            loss = population_loss_2D_v13(
                heating_cases, cooling_cases,
            )
            candidate = (
                parameters=p, theta=theta, beta=beta,
                heat_scale=heat_scale, power_scales=(
                    p.optics.scale_456, p.optics.scale_304,
                    p.optics.scale_256,
                ), loss=loss,
            )
            if best === nothing || loss.total < best.loss.total
                best = candidate
            end
            push!(trace, (
                stage="heating_shape", evaluation=evaluation,
                objective=loss.total,
                transient_loss=loss.transient,
                steady_loss=loss.steady,
                radial_loss=loss.radial,
                axial_loss=loss.axial,
                cooling_loss=loss.cooling,
                side_solid_fraction=theta[1],
                side_gas_heat_fraction=theta[2],
                active_side_conductance_W_K_m=theta[3],
                side_felt_contact_h_W_m2K=theta[4],
                beta_opt_1_m=beta,
                heat_transfer_scale=heat_scale,
                power_scale_456=p.optics.scale_456,
                power_scale_304=p.optics.scale_304,
                power_scale_256=p.optics.scale_256,
            ))
            println(
                "v13 heating $evaluation/$total_shape " *
                "loss=$(round(loss.total,digits=5)) " *
                "pop=$theta beta=$beta Sh=$heat_scale",
            )
            flush(stdout)
        end
    end

    # Power scales are refitted only after the population/shape point is
    # selected, using shared factors for the three established irradiance
    # groups.
    theta = best.theta
    ppop = rebuild_population_2D_v13(
        base;
        side_fraction=theta[1],
        side_gas_fraction=theta[2],
        active_side_G_per_m=theta[3],
        side_felt_h=theta[4],
    )
    cooling_cases = simulate_cooling_set_2D_v13(
        cooling_inputs, ppop,
    )
    base_scales = (1.25, 1.05, 0.77)
    best_power = best
    power_index = 0
    for m456 in POWER_MULTIPLIERS_2D_v13,
        m304 in POWER_MULTIPLIERS_2D_v13,
        m256 in POWER_MULTIPLIERS_2D_v13
        power_index += 1
        scales = (
            base_scales[1] * m456,
            base_scales[2] * m304,
            base_scales[3] * m256,
        )
        p = rebuild_heating_2D_v13(
            ppop; beta_opt=best.beta,
            heat_scale=best.heat_scale, power_scales=scales,
        )
        heating_cases = simulate_heating_set_2D_v13(
            heating_inputs, p,
        )
        loss = population_loss_2D_v13(
            heating_cases, cooling_cases,
        )
        candidate = (
            parameters=p, theta=theta, beta=best.beta,
            heat_scale=best.heat_scale, power_scales=scales,
            loss=loss,
        )
        if loss.total < best_power.loss.total
            best_power = candidate
        end
        push!(trace, (
            stage="power", evaluation=power_index,
            objective=loss.total,
            transient_loss=loss.transient,
            steady_loss=loss.steady,
            radial_loss=loss.radial,
            axial_loss=loss.axial,
            cooling_loss=loss.cooling,
            side_solid_fraction=theta[1],
            side_gas_heat_fraction=theta[2],
            active_side_conductance_W_K_m=theta[3],
            side_felt_contact_h_W_m2K=theta[4],
            beta_opt_1_m=best.beta,
            heat_transfer_scale=best.heat_scale,
            power_scale_456=scales[1],
            power_scale_304=scales[2],
            power_scale_256=scales[3],
        ))
        println(
            "v13 power $power_index/27 " *
            "loss=$(round(loss.total,digits=5)) scales=$scales",
        )
        flush(stdout)
    end
    _write_namedtuples_csv_2D_v13(
        joinpath(
            OUTPUT_DIR_2D_v13,
            "heating_calibration_trace_2D_v13.csv",
        ),
        trace,
    )
    open(joinpath(
        OUTPUT_DIR_2D_v13,
        "parameters_heating_selected_2D_v13.txt",
    ), "w") do io
        println(io, "objective=", best_power.loss.total)
        println(io, "side_solid_fraction=", best_power.theta[1])
        println(io, "side_gas_heat_fraction=", best_power.theta[2])
        println(io, "active_side_conductance_W_K_m=",
                best_power.theta[3])
        println(io, "side_felt_contact_h_W_m2K=",
                best_power.theta[4])
        println(io, "beta_opt_1_m=", best_power.beta)
        println(io, "heat_transfer_scale=", best_power.heat_scale)
        println(io, "power_scale_456=", best_power.power_scales[1])
        println(io, "power_scale_304=", best_power.power_scales[2])
        println(io, "power_scale_256=", best_power.power_scales[3])
    end
    println("v13 heating selected = ", best_power)
    return best_power
end

if abspath(PROGRAM_FILE) == @__FILE__
    fit_heating_stage_2D_v13()
end
