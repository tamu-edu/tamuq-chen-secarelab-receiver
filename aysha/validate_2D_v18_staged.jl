# ============================================================================
# validate_2D_v18_staged.jl - untouched heating/cooling validation and plots
# ============================================================================

using Statistics
using Plots

include("run_2D_v18.jl")

function _read_selected_2D_v18(path)
    values = Dict{String,String}()
    for line in eachline(path)
        occursin("=", line) || continue
        key, value = split(line, "="; limit=2)
        values[key] = value
    end
    return values
end

function selected_parameters_2D_v18(; mesh=:nominal)
    values = _read_selected_2D_v18(joinpath(
        OUTPUT_DIR_2D_v18, "parameters_selected_2D_v18.txt",
    ))
    number(key) = parse(Float64, values[key])
    return parameters_2D_v18(
        mesh=mesh,
        source_model=Symbol(values["source_model"]),
        deep_fraction=number("deep_fraction"),
        deep_length_m=number("deep_length_m"),
        exchange_multiplier=number("exchange_multiplier"),
        power_scales=(
            number("power_scale_456"),
            number("power_scale_304"),
            number("power_scale_256"),
        ),
        felt_conductivity_scale=
            number("felt_conductivity_scale"),
        felt_heat_capacity_scale=
            number("felt_heat_capacity_scale"),
        felt_contact_scale=number("felt_contact_scale"),
    )
end

function _phase_2D_v18(id, cooling)
    cooling && return id == "C81" ?
        "cooling_validation" : "cooling_training"
    return id in TRAIN_HEATING_2D_v18 ?
        "heating_training" : "heating_validation"
end

function _rmse_2D_v18(model, observed)
    return sqrt(mean(abs2, model .- observed))
end

function run_validation_2D_v18(; max_points=121)
    mkpath(OUTPUT_DIR_2D_v18)
    plot_root = joinpath(OUTPUT_DIR_2D_v18, "plots")
    transient_dir = joinpath(plot_root, "transients")
    axial_dir = joinpath(plot_root, "axial_profiles")
    foreach(mkpath, (plot_root, transient_dir, axial_dir))
    p = selected_parameters_2D_v18()
    specs = vcat(
        [(id=id, cooling=false) for id in HEATING_IDS_2D_v18],
        [(id=id, cooling=true) for id in COOLING_IDS_2D_v18],
    )
    cases = Vector{Any}(undef, length(specs))
    Threads.@threads for index in eachindex(specs)
        spec = specs[index]
        inputs = case_inputs_2D_v18(
            spec.id; cooling=spec.cooling,
            max_points=max_points,
        )
        cases[index] = simulate_case_2D_v18(
            inputs, p; reltol=2e-5, abstol=2e-6, dtmax=60.0,
        )
        println("v18 validation ", spec.id, " complete")
        flush(stdout)
    end

    sensor_rows = NamedTuple[]
    final_rows = NamedTuple[]
    ledger_rows = NamedTuple[]
    for case in cases
        phase = _phase_2D_v18(
            case.inputs.id, case.inputs.cooling,
        )
        for (index, sensor) in enumerate(SENSOR_NAMES_2D_v18)
            model = case.model[:, index]
            observed = case.observed[:, index]
            push!(sensor_rows, (
                phase=phase, simulation_id=case.inputs.id,
                sensor=String(sensor),
                rmse_K=_rmse_2D_v18(model, observed),
                steady_error_K=model[end] - observed[end],
                t90_error_s=V18.get_t90_2D(
                    case.inputs.times, model,
                ) - V18.get_t90_2D(
                    case.inputs.times, observed,
                ),
            ))
        end
        m = case.model[end, :]
        o = case.observed[end, :]
        result = case.result
        k = length(result.time)
        push!(final_rows, (
            phase=phase, simulation_id=case.inputs.id,
            flow_lpm=mean(case.inputs.flow),
            model_T8_K=m[1], observed_T8_K=o[1],
            model_T12_K=m[2], observed_T12_K=o[2],
            model_T11_K=m[3], observed_T11_K=o[3],
            model_T9_K=m[4], observed_T9_K=o[4],
            model_T10_K=m[5], observed_T10_K=o[5],
            model_T3_K=m[6], observed_T3_K=o[6],
            model_T2_K=m[7], observed_T2_K=o[7],
            model_T12_minus_T8_K=m[2] - m[1],
            observed_T12_minus_T8_K=o[2] - o[1],
            model_T11_minus_T12_K=m[3] - m[2],
            observed_T11_minus_T12_K=o[3] - o[2],
            model_T12_minus_T9_K=m[2] - m[4],
            observed_T12_minus_T9_K=o[2] - o[4],
            model_T11_minus_T10_K=m[3] - m[5],
            observed_T11_minus_T10_K=o[3] - o[5],
            core_to_side_flow_ratio=
                result.channel_mass_flow_kg_s[result.core_group, k] /
                result.channel_mass_flow_kg_s[result.side_group, k],
            dp1_model_mbar=result.dp1_prediction_mbar[end],
            dp1_observed_mbar=case.inputs.dp1[end],
        ))
        ledger = V18.energy_rate_ledger2D(
            result.ode_solution.u[end], p,
            result.operating_condition, result.time[end],
        )
        push!(ledger_rows, (
            phase=phase, simulation_id=case.inputs.id,
            residual_W=ledger.residual,
            solar_deposited_W=ledger.solar_deposited,
            casing_flange_loss_W=ledger.casing_flange_loss,
            source_power_error=ledger.source_power_error,
        ))

        primary_indices = (1, 2, 3, 6, 7)
        transient = plot(
            xlabel="Time (s)", ylabel="Temperature (K)",
            title="2D v18 $(case.inputs.id): primary observables",
            size=(920, 650), legend=:outerright,
        )
        colors = (:red, :orange, :purple, :blue, :brown)
        for (local_index, sensor_index) in enumerate(primary_indices)
            sensor = SENSOR_NAMES_2D_v18[sensor_index]
            plot!(
                transient, case.inputs.times,
                case.observed[:, sensor_index];
                color=colors[local_index], linestyle=:dash,
                label="$(sensor) measured",
            )
            plot!(
                transient, case.inputs.times,
                case.model[:, sensor_index];
                color=colors[local_index], linewidth=2,
                label="$(sensor) model",
            )
        end
        savefig(transient, joinpath(
            transient_dir, "$(case.inputs.id)_transient.png",
        ))
        if !case.inputs.cooling
            axial = plot(
                [11.0, 58.0, 107.0], o[1:3];
                marker=:circle, linestyle=:dash,
                label="measured", xlabel="Axial position (mm)",
                ylabel="Side-wall temperature (K)",
                title="2D v18 $(case.inputs.id): side axial profile",
                size=(720, 520),
            )
            plot!(
                axial, [11.0, 58.0, 107.0], m[1:3];
                marker=:square, linewidth=2, label="model",
            )
            savefig(axial, joinpath(
                axial_dir, "$(case.inputs.id)_axial.png",
            ))
        end
    end

    aggregate_rows = NamedTuple[]
    for phase in unique(row.phase for row in sensor_rows)
        selected = filter(row -> row.phase == phase, sensor_rows)
        primary = filter(
            row -> row.sensor in ("T8", "T12", "T11", "T3", "T2"),
            selected,
        )
        side = filter(
            row -> row.sensor in ("T8", "T12", "T11"), selected,
        )
        air = filter(row -> row.sensor == "T3", selected)
        felt = filter(row -> row.sensor == "T2", selected)
        phase_final = filter(row -> row.phase == phase, final_rows)
        push!(aggregate_rows, (
            phase=phase,
            primary_mean_sensor_rmse_K=
                mean(row.rmse_K for row in primary),
            primary_pooled_rmse_K=sqrt(
                mean(row.rmse_K^2 for row in primary),
            ),
            side_mean_rmse_K=mean(row.rmse_K for row in side),
            air_mean_rmse_K=mean(row.rmse_K for row in air),
            felt_mean_rmse_K=mean(row.rmse_K for row in felt),
            axial_inversion_rmse_K=sqrt(mean(
                (
                    row.model_T12_minus_T8_K -
                    row.observed_T12_minus_T8_K
                )^2 for row in phase_final
            )),
        ))
    end
    _write_namedtuples_csv_2D_v18(joinpath(
        OUTPUT_DIR_2D_v18, "staged_sensor_metrics_2D_v18.csv",
    ), sensor_rows)
    _write_namedtuples_csv_2D_v18(joinpath(
        OUTPUT_DIR_2D_v18, "staged_final_profiles_2D_v18.csv",
    ), final_rows)
    _write_namedtuples_csv_2D_v18(joinpath(
        OUTPUT_DIR_2D_v18, "staged_aggregate_2D_v18.csv",
    ), aggregate_rows)
    _write_namedtuples_csv_2D_v18(joinpath(
        OUTPUT_DIR_2D_v18, "energy_ledgers_2D_v18.csv",
    ), ledger_rows)

    parity = plot(
        xlabel="Measured temperature (K)",
        ylabel="Modeled temperature (K)",
        title="2D v18 primary-observable parity",
        size=(760, 680), legend=:topleft,
    )
    phase_colors = Dict(
        "heating_training" => :royalblue,
        "heating_validation" => :darkorange,
        "cooling_training" => :seagreen,
        "cooling_validation" => :purple,
    )
    all_values = Float64[]
    for phase in keys(phase_colors)
        xs = Float64[]
        ys = Float64[]
        for row in filter(r -> r.phase == phase, final_rows)
            for sensor in ("T8", "T12", "T11", "T3", "T2")
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
            color=phase_colors[phase], alpha=0.80,
        )
    end
    bounds = extrema(all_values)
    plot!(parity, collect(bounds), collect(bounds);
          color=:black, linestyle=:dash, label="1:1")
    savefig(parity, joinpath(plot_root, "parity_primary_2D_v18.png"))
    foreach(println, aggregate_rows)
    return (
        cases=cases, sensors=sensor_rows, finals=final_rows,
        aggregates=aggregate_rows, ledgers=ledger_rows,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_validation_2D_v18()
end
