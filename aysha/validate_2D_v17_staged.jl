# ============================================================================
# validate_2D_v17_staged.jl - complete flange-corrected validation
# ============================================================================

using Statistics
using Plots

include("run_2D_v17.jl")

function _read_selected_2D_v17(path)
    values = Dict{String,Float64}()
    for line in eachline(path)
        occursin("=", line) || continue
        key, raw = split(line, "="; limit=2)
        value = tryparse(Float64, raw)
        value === nothing || (values[key] = value)
    end
    return values
end

function selected_parameters_2D_v17(; nominal_mesh=true)
    hardware = _read_selected_2D_v17(joinpath(
        OUTPUT_DIR_2D_v17,
        "parameters_casing_flange_selected_2D_v17.txt",
    ))
    power = _read_selected_2D_v17(joinpath(
        OUTPUT_DIR_2D_v17,
        "parameters_power_selected_2D_v17.txt",
    ))
    return parameters_2D_v17(
        nominal_mesh=nominal_mesh,
        screen_mesh=!nominal_mesh,
        felt_heat_capacity_scale=
            hardware["felt_heat_capacity_scale"],
        casing_flange_conductance_W_K=
            hardware["casing_flange_conductance_W_K"],
        power_scales=(
            power["power_scale_456"],
            power["power_scale_304"],
            power["power_scale_256"],
        ),
    )
end

function _phase_2D_v17(id, cooling)
    cooling && return "cooling_validation"
    return id in TRAIN_HEATING_2D_v17 ?
        "heating_training" : "heating_validation"
end

function run_validation_2D_v17(; max_points=121)
    mkpath(OUTPUT_DIR_2D_v17)
    plot_root = joinpath(OUTPUT_DIR_2D_v17, "plots")
    transient_dir = joinpath(plot_root, "transients")
    axial_dir = joinpath(plot_root, "axial_profiles")
    foreach(mkpath, (plot_root, transient_dir, axial_dir))

    p = selected_parameters_2D_v17()
    specs = vcat(
        [(id=id, cooling=false) for id in HEATING_IDS_2D_v17],
        [(id=id, cooling=true) for id in COOLING_IDS_2D_v17],
    )
    cases = Vector{Any}(undef, length(specs))
    Threads.@threads for index in eachindex(specs)
        spec = specs[index]
        input = case_inputs_2D_v17(
            spec.id; cooling=spec.cooling,
            max_points=max_points,
        )
        cases[index] = simulate_case_2D_v17(input, p)
    end

    sensors = NamedTuple[]
    finals = NamedTuple[]
    ledgers = NamedTuple[]
    for case in cases
        phase = _phase_2D_v17(
            case.inputs.id, case.inputs.cooling,
        )
        for (index, sensor) in enumerate(SENSOR_NAMES_2D_v17)
            model = case.model[:, index]
            observed = case.observed[:, index]
            push!(sensors, (
                phase=phase,
                simulation_id=case.inputs.id,
                sensor=String(sensor),
                rmse_K=sqrt(mean(abs2, model .- observed)),
                steady_error_K=model[end] - observed[end],
                t90_error_s=V17.get_t90_2D(
                    case.inputs.times, model,
                ) - V17.get_t90_2D(
                    case.inputs.times, observed,
                ),
            ))
        end
        m = case.model[end, :]
        e = case.observed[end, :]
        result = case.result
        k = length(result.time)
        core_flow = result.channel_mass_flow_kg_s[
            result.core_group, k,
        ]
        side_flow = result.channel_mass_flow_kg_s[
            result.side_group, k,
        ]
        corner_flow = result.channel_mass_flow_kg_s[
            result.corner_group, k,
        ]
        push!(finals, (
            phase=phase, simulation_id=case.inputs.id,
            flow_lpm=mean(case.inputs.flow),
            model_T8_K=m[1], observed_T8_K=e[1],
            model_T12_K=m[2], observed_T12_K=e[2],
            model_T11_K=m[3], observed_T11_K=e[3],
            model_T9_K=m[4], observed_T9_K=e[4],
            model_T10_K=m[5], observed_T10_K=e[5],
            model_T3_K=m[6], observed_T3_K=e[6],
            model_T2_K=m[7], observed_T2_K=e[7],
            model_T12_minus_T8_K=m[2] - m[1],
            observed_T12_minus_T8_K=e[2] - e[1],
            model_T12_minus_T9_K=m[2] - m[4],
            observed_T12_minus_T9_K=e[2] - e[4],
            model_T11_minus_T10_K=m[3] - m[5],
            observed_T11_minus_T10_K=e[3] - e[5],
            core_to_side_flow_ratio=core_flow / side_flow,
            core_to_corner_flow_ratio=core_flow / corner_flow,
            dp1_model_mbar=result.dp1_prediction_mbar[end],
            dp1_observed_mbar=case.inputs.dp1[end],
            equal_pressure_error=
                result.equal_pressure_relative_error[end],
        ))
        op = result.operating_condition
        ledger = V17.energy_rate_ledger2D(
            result.ode_solution.u[end], p, op,
            result.time[end],
        )
        push!(ledgers, (
            phase=phase, simulation_id=case.inputs.id,
            residual_W=ledger.residual,
            solar_deposited_W=ledger.solar_deposited,
            tube_flange_loss_W=ledger.tube_flange_loss,
            casing_flange_loss_W=
                ledger.casing_flange_loss,
            total_flange_loss_W=ledger.flange_loss,
            housing_ambient_loss_W=ledger.housing_loss,
        ))
    end

    aggregates = NamedTuple[]
    for phase in unique(row.phase for row in sensors)
        phase_sensors = filter(
            row -> row.phase == phase, sensors,
        )
        phase_finals = filter(
            row -> row.phase == phase, finals,
        )
        push!(aggregates, (
            phase=phase,
            mean_sensor_rmse_K=mean(
                row.rmse_K for row in phase_sensors
            ),
            steady_mae_K=mean(
                abs(row.steady_error_K)
                for row in phase_sensors
            ),
            t90_mae_s=mean(
                abs(row.t90_error_s)
                for row in phase_sensors
            ),
            axial_rmse_K=sqrt(mean(
                (
                    row.model_T12_minus_T8_K -
                    row.observed_T12_minus_T8_K
                )^2 for row in phase_finals
            )),
            mid_radial_rmse_K=sqrt(mean(
                (
                    row.model_T12_minus_T9_K -
                    row.observed_T12_minus_T9_K
                )^2 for row in phase_finals
            )),
            deep_radial_rmse_K=sqrt(mean(
                (
                    row.model_T11_minus_T10_K -
                    row.observed_T11_minus_T10_K
                )^2 for row in phase_finals
            )),
            dp1_rmse_mbar=sqrt(mean(
                (
                    row.dp1_model_mbar -
                    row.dp1_observed_mbar
                )^2 for row in phase_finals
            )),
        ))
    end

    # The independent assessment proposed pure gas re-observation, but the
    # actual model already uses an 80% wall / 20% gas target.  Quantify the
    # complete steady-state blend sensitivity without re-solving.
    heating_cases = filter(
        case -> !case.inputs.cooling, cases,
    )
    observation_rows = NamedTuple[]
    for fraction in (0.0, 0.25, 0.55, 0.75, 0.80, 1.0)
        mid_errors = Float64[]
        deep_errors = Float64[]
        mid_sign = 0
        deep_sign = 0
        for case in heating_cases
            pred = case.predictions
            t9 = fraction * pred.T9_wall[end] +
                 (1.0 - fraction) * pred.T9_gas[end]
            t10 = fraction * pred.T10_wall[end] +
                  (1.0 - fraction) * pred.T10_gas[end]
            model_mid = pred.T12_skin[end] - t9
            model_deep = pred.T11_skin[end] - t10
            observed_mid =
                case.observed[end, 2] - case.observed[end, 4]
            observed_deep =
                case.observed[end, 3] - case.observed[end, 5]
            push!(mid_errors, model_mid - observed_mid)
            push!(deep_errors, model_deep - observed_deep)
            mid_sign += sign(model_mid) == sign(observed_mid)
            deep_sign += sign(model_deep) == sign(observed_deep)
        end
        push!(observation_rows, (
            wall_fraction=fraction,
            gas_fraction=1.0 - fraction,
            mid_radial_rmse_K=sqrt(mean(abs2, mid_errors)),
            deep_radial_rmse_K=sqrt(mean(abs2, deep_errors)),
            mid_sign_correct=mid_sign,
            deep_sign_correct=deep_sign,
        ))
    end

    _write_namedtuples_csv_2D_v17(
        joinpath(
            OUTPUT_DIR_2D_v17,
            "staged_sensor_metrics_2D_v17.csv",
        ), sensors,
    )
    _write_namedtuples_csv_2D_v17(
        joinpath(
            OUTPUT_DIR_2D_v17,
            "staged_final_profiles_2D_v17.csv",
        ), finals,
    )
    _write_namedtuples_csv_2D_v17(
        joinpath(
            OUTPUT_DIR_2D_v17,
            "staged_aggregate_2D_v17.csv",
        ), aggregates,
    )
    _write_namedtuples_csv_2D_v17(
        joinpath(
            OUTPUT_DIR_2D_v17,
            "energy_ledgers_2D_v17.csv",
        ), ledgers,
    )
    _write_namedtuples_csv_2D_v17(
        joinpath(
            OUTPUT_DIR_2D_v17,
            "observation_operator_sensitivity_2D_v17.csv",
        ), observation_rows,
    )

    phase_colors = Dict(
        "heating_training" => :royalblue,
        "heating_validation" => :darkorange,
        "cooling_validation" => :seagreen,
    )
    parity = plot(
        xlabel="Observed final temperature (K)",
        ylabel="Modeled final temperature (K)",
        title="2D v17 casing/flange-corrected parity",
        legend=:topleft, size=(850, 720),
    )
    all_values = Float64[]
    for phase in keys(phase_colors)
        xs = Float64[]
        ys = Float64[]
        for row in filter(r -> r.phase == phase, finals)
            for sensor in (
                "T8", "T12", "T11", "T9",
                "T10", "T3", "T2",
            )
                push!(xs, getproperty(
                    row, Symbol("observed_$(sensor)_K"),
                ))
                push!(ys, getproperty(
                    row, Symbol("model_$(sensor)_K"),
                ))
            end
        end
        append!(all_values, xs)
        append!(all_values, ys)
        scatter!(
            parity, xs, ys; label=replace(phase, "_" => " "),
            color=phase_colors[phase], alpha=0.8,
        )
    end
    lo, hi = extrema(all_values)
    plot!(
        parity, [lo, hi], [lo, hi]; label="1:1",
        color=:black, linestyle=:dash,
    )
    savefig(
        parity,
        joinpath(plot_root, "parity_staged_2D_v17.png"),
    )

    colors = (
        :red, :orange, :purple, :blue,
        :cyan, :green, :black,
    )
    for case in cases
        panels = Any[]
        for (index, sensor) in enumerate(SENSOR_NAMES_2D_v17)
            panel = plot(
                case.inputs.times, case.observed[:, index];
                label="measured", color=:black, linewidth=2,
                ylabel="$(sensor) (K)",
            )
            plot!(
                panel, case.inputs.times, case.model[:, index];
                label="v17", color=colors[index],
                linewidth=2, linestyle=:dash,
            )
            push!(panels, panel)
        end
        transient_plot = plot(
            panels...; layout=(4, 2), size=(1100, 1000),
            plot_title="2D v17 $(case.inputs.id)",
            xlabel="Time (s)",
        )
        savefig(
            transient_plot,
            joinpath(
                transient_dir,
                "transient_$(case.inputs.id)_2D_v17.png",
            ),
        )
        case.inputs.cooling && continue
        result = case.result
        zmm = 1000.0 .* result.z_solid
        axial = plot(
            zmm, result.skin_temperature[:, end];
            label="exterior SiC skin", linewidth=3,
            color=:red, xlabel="Axial position (mm)",
            ylabel="Final temperature (K)",
            title="2D v17 axial profile $(case.inputs.id)",
        )
        plot!(
            axial, zmm,
            result.channel_temperature[
                result.core_group, :, end,
            ]; label="central channel orbit",
            color=:blue, linewidth=2,
        )
        scatter!(
            axial, [11.0, 58.0, 107.0],
            case.observed[end, 1:3];
            label="measured side TCs", color=:black,
        )
        savefig(
            axial,
            joinpath(
                axial_dir,
                "axial_profile_$(case.inputs.id)_2D_v17.png",
            ),
        )
    end

    hardware = p.casing_flange
    power = (
        p.optics.scale_456,
        p.optics.scale_304,
        p.optics.scale_256,
    )
    heating_t2 = filter(
        row -> row.phase != "cooling_validation" &&
               row.sensor == "T2",
        sensors,
    )
    t2_bias = mean(row.steady_error_K for row in heating_t2)
    t2_t90 = mean(abs(row.t90_error_s) for row in heating_t2)
    max_energy_residual =
        maximum(abs(row.residual_W) for row in ledgers)
    open(joinpath(
        OUTPUT_DIR_2D_v17,
        "staged_acceptance_status_2D_v17.txt",
    ), "w") do io
        println(
            io, "casing_flange_conductance_W_K=",
            hardware.conductance_W_K,
        )
        println(
            io, "felt_conductivity_scale_fixed=",
            p.assembly.felt_conductivity_scale,
        )
        println(
            io, "felt_heat_capacity_scale=",
            p.assembly.felt_heat_capacity_scale,
        )
        println(io, "power_scale_456=", power[1])
        println(io, "power_scale_304=", power[2])
        println(io, "power_scale_256=", power[3])
        println(io, "heating_T2_mean_steady_error_K=", t2_bias)
        println(io, "heating_T2_t90_MAE_s=", t2_t90)
        println(io, "max_energy_residual_W=", max_energy_residual)
        foreach(row -> println(io, row), aggregates)
    end
    foreach(println, aggregates)
    foreach(println, observation_rows)
    println(
        "v17 heating T2 bias=$(round(t2_bias,digits=2)) K; " *
        "T2 t90 MAE=$(round(t2_t90,digits=1)) s; " *
        "max energy residual=$(max_energy_residual) W",
    )
    return (
        cases=cases, sensors=sensors, finals=finals,
        aggregates=aggregates, ledgers=ledgers,
        observation_rows=observation_rows,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_validation_2D_v17()
end
