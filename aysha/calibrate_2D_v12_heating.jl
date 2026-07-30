# ============================================================================
# calibrate_2D_v12_heating.jl
#
# Cooling-selected assembly/observation parameters are frozen.  Stage 2A
# profiles beta_opt and the one order-unity Graetz scale on the nine heating
# training cases.  Stage 2B profiles each irradiance-group power scale.
# SiC k(T), Cp(T), density and the cooling-selected felt/contact terms do not
# move in this stage.
# ============================================================================

using Statistics

include("run_2D_v12.jl")

const HEATING_FIT_POINTS_2D_v12 = 31
const PROFILE_CANDIDATES_HEAT_2D_v12 = (
    (80.0, 0.80), (80.0, 1.20),
    (150.0, 0.80), (150.0, 1.20),
    (220.0, 0.80), (220.0, 1.20),
)
const POWER_CANDIDATES_HEAT_2D_v12 = (
    (1.25, 1.60, 1.90, 2.20),
    (1.05, 1.40, 1.75, 2.10),
    (0.77, 1.00, 1.25, 1.50),
)
const HEATING_GROUP_IDS_2D_v12 = (
    ("E67", "E69", "E71"),
    ("E72", "E74", "E76"),
    ("E77", "E79", "E81"),
)

function cooling_selected_parameters_2D_v12(;
    nominal_mesh=false, screen_mesh=false,
)
    base = base_parameters_2D_v12(
        ; nominal_mesh=nominal_mesh, screen_mesh=screen_mesh,
    )
    return rebuild_parameters_2D_v12(
        base;
        felt_k_scale=1.20,
        felt_cp_scale=1.50,
        receiver_adaptor_h=15.0,
        side_tau=180.0,
        interior_tau=60.0,
        outlet_tau=180.0,
        interior_wall_fraction=0.80,
    )
end

function heating_profile_loss_2D_v12(cases)
    residuals = Float64[]
    group_axial_model = Dict{String, Float64}()
    group_axial_observed = Dict{String, Float64}()
    for case in cases
        m = case.model[end, :]
        e = case.observed[end, :]
        append!(residuals, (
            ((m[2] - m[1]) - (e[2] - e[1])) / 50.0,
            ((m[3] - m[2]) - (e[3] - e[2])) / 40.0,
            ((m[5] - m[4]) - (e[5] - e[4])) / 40.0,
            (
                (m[6] - 0.5 * (m[3] + m[5])) -
                (e[6] - 0.5 * (e[3] + e[5]))
            ) / 60.0,
        ))
        group_axial_model[case.inputs.id] = m[2] - m[1]
        group_axial_observed[case.inputs.id] = e[2] - e[1]
    end
    for ids in HEATING_GROUP_IDS_2D_v12
        available = [id for id in ids if haskey(group_axial_model, id)]
        length(available) >= 2 || continue
        model_mean = mean(group_axial_model[id] for id in available)
        observed_mean = mean(group_axial_observed[id] for id in available)
        for id in available
            push!(
                residuals,
                (
                    (group_axial_model[id] - model_mean) -
                    (group_axial_observed[id] - observed_mean)
                ) / 25.0,
            )
        end
    end
    return mean(abs2, residuals)
end

function heating_level_loss_2D_v12(cases)
    residuals = Float64[]
    for case in cases
        m = case.model
        e = case.observed
        for sensor in eachindex(SENSOR_NAMES_2D_v12)
            append!(residuals, (m[:, sensor] .- e[:, sensor]) ./ 100.0)
            push!(residuals, (m[end, sensor] - e[end, sensor]) / 60.0)
        end
    end
    return mean(abs2, residuals)
end

function simulate_training_set_2D_v12(p; ids=TRAIN_HEATING_2D_v12)
    inputs = [
        case_inputs_2D_v12(
            id; max_points=HEATING_FIT_POINTS_2D_v12,
        )
        for id in ids
    ]
    cases = Vector{Any}(undef, length(inputs))
    Threads.@threads for index in eachindex(inputs)
        cases[index] = simulate_case_2D_v12(inputs[index], p)
    end
    return cases
end

function fit_heating_stage_2D_v12()
    mkpath(OUTPUT_DIR_2D_v12)
    cooling = cooling_selected_parameters_2D_v12(
        nominal_mesh=false, screen_mesh=true,
    )
    profile_trace = NamedTuple[]
    profile_cases = Vector{Any}(
        undef, length(PROFILE_CANDIDATES_HEAT_2D_v12),
    )
    for (index, theta) in enumerate(PROFILE_CANDIDATES_HEAT_2D_v12)
        beta_opt, heat_scale = theta
        p = rebuild_parameters_2D_v12(
            cooling; beta_opt=beta_opt, heat_scale=heat_scale,
        )
        cases = simulate_training_set_2D_v12(p)
        profile_cases[index] = cases
        profile_loss = heating_profile_loss_2D_v12(cases)
        level_loss = heating_level_loss_2D_v12(cases)
        push!(profile_trace, (
            evaluation=index,
            beta_opt_1_m=beta_opt,
            heat_transfer_scale=heat_scale,
            profile_loss=profile_loss,
            level_loss=level_loss,
        ))
        println(
            "heating profile $index/" *
            "$(length(PROFILE_CANDIDATES_HEAT_2D_v12)): " *
            "profile=$(round(profile_loss,digits=6)) " *
            "level=$(round(level_loss,digits=6)) theta=$theta",
        )
        flush(stdout)
    end
    profile_index = argmin(row.profile_loss for row in profile_trace)
    beta_selected, heat_selected =
        PROFILE_CANDIDATES_HEAT_2D_v12[profile_index]

    power = [
        cooling.optics.scale_456,
        cooling.optics.scale_304,
        cooling.optics.scale_256,
    ]
    power_trace = NamedTuple[]
    group_losses = Float64[]
    for group_index in 1:3
        losses = Float64[]
        for (candidate_index, scale) in
            enumerate(POWER_CANDIDATES_HEAT_2D_v12[group_index])

            trial = copy(power)
            trial[group_index] = scale
            p = rebuild_parameters_2D_v12(
                cooling;
                beta_opt=beta_selected,
                heat_scale=heat_selected,
                power_scales=Tuple(trial),
            )
            cases = simulate_training_set_2D_v12(
                p; ids=HEATING_GROUP_IDS_2D_v12[group_index],
            )
            loss = heating_level_loss_2D_v12(cases) +
                   0.20 * heating_profile_loss_2D_v12(cases)
            push!(losses, loss)
            push!(power_trace, (
                group=group_index,
                evaluation=candidate_index,
                tested_scale=scale,
                objective=loss,
            ))
            println(
                "power group $group_index candidate $candidate_index: " *
                "loss=$(round(loss,digits=6)) scale=$scale",
            )
            flush(stdout)
        end
        selected = argmin(losses)
        power[group_index] =
            POWER_CANDIDATES_HEAT_2D_v12[group_index][selected]
        push!(group_losses, losses[selected])
    end

    selected = rebuild_parameters_2D_v12(
        cooling;
        beta_opt=beta_selected,
        heat_scale=heat_selected,
        power_scales=Tuple(power),
    )
    _write_namedtuples_csv_2D_v12(
        joinpath(
            OUTPUT_DIR_2D_v12,
            "heating_profile_trace_2D_v12.csv",
        ),
        profile_trace,
    )
    _write_namedtuples_csv_2D_v12(
        joinpath(
            OUTPUT_DIR_2D_v12,
            "heating_power_trace_2D_v12.csv",
        ),
        power_trace,
    )
    open(joinpath(
        OUTPUT_DIR_2D_v12,
        "parameters_staged_selected_2D_v12.txt",
    ), "w") do io
        println(io, "felt_conductivity_scale=",
            selected.assembly.felt_conductivity_scale)
        println(io, "felt_heat_capacity_scale=",
            selected.assembly.felt_heat_capacity_scale)
        println(io, "receiver_adaptor_contact_h_W_m2K=",
            selected.assembly.receiver_adaptor_contact_h)
        println(io, "side_time_constant_s=",
            selected.observation.side_time_constant_s)
        println(io, "interior_time_constant_s=",
            selected.observation.interior_time_constant_s)
        println(io, "outlet_time_constant_s=",
            selected.observation.outlet_time_constant_s)
        println(io, "interior_wall_fraction=",
            selected.observation.interior_wall_fraction)
        println(io, "beta_opt_1_m=", beta_selected)
        println(io, "heat_transfer_scale=", heat_selected)
        println(io, "power_scale_456=", power[1])
        println(io, "power_scale_304=", power[2])
        println(io, "power_scale_256=", power[3])
    end
    return (
        parameters=selected,
        profile_trace=profile_trace,
        power_trace=power_trace,
        selected_profile_index=profile_index,
        group_losses=group_losses,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    result = fit_heating_stage_2D_v12()
    println("selected profile index = ",
        result.selected_profile_index)
    println("selected parameters written")
end
