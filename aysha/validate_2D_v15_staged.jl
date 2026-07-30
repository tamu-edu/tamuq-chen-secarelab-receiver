# ============================================================================
# validate_2D_v15_staged.jl - nominal validation and profile plots
# ============================================================================

using Statistics
using Plots

include("run_2D_v15.jl")

function _read_selected_skin_2D_v15()
    path = joinpath(
        OUTPUT_DIR_2D_v15,
        "parameters_skin_selected_2D_v15.txt",
    )
    values = Dict{String,Float64}()
    for line in eachline(path)
        occursin("=", line) || continue
        key, raw = split(line, "="; limit=2)
        value = tryparse(Float64, raw)
        value === nothing || (values[key] = value)
    end
    return (
        thickness=values["skin_thickness_mm"] * 1e-3,
        conductance=values["channel_conductance_scale"],
    )
end

function validation_phase_2D_v15(id, cooling)
    cooling && return "cooling_validation"
    return id in TRAIN_HEATING_2D_v15 ?
        "heating_training" : "heating_validation"
end

function _channel_map_2D_v15(result, values)
    grid = V15.build_network_grid2D(result.parameters)
    return [
        values[grid.group_map[i, j]]
        for j in 1:grid.channel_side_count,
            i in 1:grid.channel_side_count
    ]
end

function run_validation_2D_v15(; max_points=121)
    mkpath(OUTPUT_DIR_2D_v15)
    plot_root = joinpath(OUTPUT_DIR_2D_v15, "plots")
    transient_dir = joinpath(plot_root, "transients")
    axial_dir = joinpath(plot_root, "axial_profiles")
    map_dir = joinpath(plot_root, "channel_maps")
    mkpath(transient_dir)
    mkpath(axial_dir)
    mkpath(map_dir)
    selected = _read_selected_skin_2D_v15()
    p = selected_parameters_2D_v15(
        nominal_mesh=true,
        skin_thickness=selected.thickness,
        skin_conductance_scale=selected.conductance,
    )
    specs = vcat(
        [(id=id, cooling=false) for id in HEATING_IDS_2D_v15],
        [(id=id, cooling=true) for id in COOLING_IDS_2D_v15],
    )
    cases = Vector{Any}(undef, length(specs))
    Threads.@threads for index in eachindex(specs)
        spec = specs[index]
        inputs = case_inputs_2D_v15(
            spec.id; cooling=spec.cooling, max_points=max_points,
        )
        cases[index] = simulate_case_2D_v15(inputs, p)
    end

    sensor_rows = NamedTuple[]
    final_rows = NamedTuple[]
    transient_rows = NamedTuple[]
    for case in cases
        phase = validation_phase_2D_v15(
            case.inputs.id, case.inputs.cooling,
        )
        for (index, sensor) in enumerate(SENSOR_NAMES_2D_v15)
            model = case.model[:, index]
            observed = case.observed[:, index]
            push!(sensor_rows, (
                phase=phase, simulation_id=case.inputs.id,
                sensor=String(sensor),
                rmse_K=sqrt(mean(abs2, model .- observed)),
                steady_error_K=model[end] - observed[end],
                t90_error_s=V15.get_t90_2D(
                    case.inputs.times, model,
                ) - V15.get_t90_2D(
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
        push!(final_rows, (
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
            final_skin_max_minus_side_channel_K=
                maximum(result.skin_temperature[:, end]) -
                maximum(result.channel_temperature[
                    result.side_group, :, end,
                ]),
        ))
        for index in eachindex(case.inputs.times)
            push!(transient_rows, (
                phase=phase, simulation_id=case.inputs.id,
                time_s=case.inputs.times[index],
                model_T8_K=case.model[index, 1],
                observed_T8_K=case.observed[index, 1],
                model_T12_K=case.model[index, 2],
                observed_T12_K=case.observed[index, 2],
                model_T11_K=case.model[index, 3],
                observed_T11_K=case.observed[index, 3],
                model_T9_K=case.model[index, 4],
                observed_T9_K=case.observed[index, 4],
                model_T10_K=case.model[index, 5],
                observed_T10_K=case.observed[index, 5],
                model_T3_K=case.model[index, 6],
                observed_T3_K=case.observed[index, 6],
                model_T2_K=case.model[index, 7],
                observed_T2_K=case.observed[index, 7],
            ))
        end
    end

    aggregate_rows = NamedTuple[]
    for phase in unique(row.phase for row in sensor_rows)
        sensors = filter(row -> row.phase == phase, sensor_rows)
        finals = filter(row -> row.phase == phase, final_rows)
        push!(aggregate_rows, (
            phase=phase,
            mean_sensor_rmse_K=mean(row.rmse_K for row in sensors),
            steady_mae_K=mean(
                abs(row.steady_error_K) for row in sensors
            ),
            t90_mae_s=mean(abs(row.t90_error_s) for row in sensors),
            axial_rmse_K=sqrt(mean(
                (row.model_T12_minus_T8_K -
                 row.observed_T12_minus_T8_K)^2
                for row in finals
            )),
            mid_radial_rmse_K=sqrt(mean(
                (row.model_T12_minus_T9_K -
                 row.observed_T12_minus_T9_K)^2
                for row in finals
            )),
            deep_radial_rmse_K=sqrt(mean(
                (row.model_T11_minus_T10_K -
                 row.observed_T11_minus_T10_K)^2
                for row in finals
            )),
            dp1_rmse_mbar=sqrt(mean(
                (row.dp1_model_mbar - row.dp1_observed_mbar)^2
                for row in finals
            )),
            mean_core_to_side_flow_ratio=mean(
                row.core_to_side_flow_ratio for row in finals
            ),
            mean_core_to_corner_flow_ratio=mean(
                row.core_to_corner_flow_ratio for row in finals
            ),
        ))
    end
    _write_namedtuples_csv_2D_v15(
        joinpath(OUTPUT_DIR_2D_v15,
                 "staged_sensor_metrics_2D_v15.csv"),
        sensor_rows,
    )
    _write_namedtuples_csv_2D_v15(
        joinpath(OUTPUT_DIR_2D_v15,
                 "staged_final_profiles_2D_v15.csv"),
        final_rows,
    )
    _write_namedtuples_csv_2D_v15(
        joinpath(OUTPUT_DIR_2D_v15,
                 "staged_aggregate_2D_v15.csv"),
        aggregate_rows,
    )
    _write_namedtuples_csv_2D_v15(
        joinpath(OUTPUT_DIR_2D_v15,
                 "staged_transients_2D_v15.csv"),
        transient_rows,
    )

    heating_finals = filter(
        row -> row.phase != "cooling_validation", final_rows,
    )
    mid_sign = count(
        sign(row.model_T12_minus_T9_K) ==
        sign(row.observed_T12_minus_T9_K)
        for row in heating_finals
    )
    deep_sign = count(
        sign(row.model_T11_minus_T10_K) ==
        sign(row.observed_T11_minus_T10_K)
        for row in heating_finals
    )
    validation = only(filter(
        row -> row.phase == "heating_validation", aggregate_rows,
    ))
    cooling = only(filter(
        row -> row.phase == "cooling_validation", aggregate_rows,
    ))
    accepted =
        mid_sign >= 12 && deep_sign >= 12 &&
        validation.mid_radial_rmse_K < 20.0 &&
        validation.deep_radial_rmse_K < 20.0 &&
        validation.mean_sensor_rmse_K <= 81.57 &&
        cooling.mean_sensor_rmse_K <= 37.73
    open(joinpath(
        OUTPUT_DIR_2D_v15,
        "staged_acceptance_status_2D_v15.txt",
    ), "w") do io
        println(io, "mid_radial_sign_correct=$mid_sign/15")
        println(io, "deep_radial_sign_correct=$deep_sign/15")
        foreach(row -> println(io, row), aggregate_rows)
        println(
            io, "status=",
            accepted ? "ACCEPTED" :
            "REJECTED_FOR_COEFFICIENT_EXTRACTION",
        )
    end

    phase_colors = Dict(
        "heating_training" => :royalblue,
        "heating_validation" => :darkorange,
        "cooling_validation" => :seagreen,
    )
    parity = plot(
        xlabel="Observed final temperature (K)",
        ylabel="Modeled final temperature (K)",
        title="2D v15 exterior-wall parity",
        legend=:topleft, size=(850, 720),
    )
    all_values = Float64[]
    for phase in keys(phase_colors)
        xs = Float64[]
        ys = Float64[]
        for row in filter(r -> r.phase == phase, final_rows)
            for sensor in ("T8", "T12", "T11", "T9", "T10", "T3", "T2")
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
    savefig(parity, joinpath(plot_root, "parity_staged_2D_v15.png"))

    colors = (:red, :orange, :purple, :blue, :cyan, :green, :black)
    for case in cases
        panels = Any[]
        for (index, sensor) in enumerate(SENSOR_NAMES_2D_v15)
            panel = plot(
                case.inputs.times, case.observed[:, index];
                label="measured", color=:black, linewidth=2,
                ylabel="$(sensor) (K)",
            )
            plot!(
                panel, case.inputs.times, case.model[:, index];
                label="v15", color=colors[index], linewidth=2,
                linestyle=:dash,
            )
            push!(panels, panel)
        end
        transient_plot = plot(
            panels...; layout=(4, 2), size=(1100, 1000),
            plot_title="2D v15 $(case.inputs.id)",
            xlabel="Time (s)",
        )
        savefig(
            transient_plot,
            joinpath(
                transient_dir,
                "transient_$(case.inputs.id)_2D_v15.png",
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
            title="2D v15 axial profile $(case.inputs.id)",
        )
        plot!(
            axial, zmm,
            result.channel_temperature[
                result.side_group, :, end,
            ]; label="outer channel orbit", linewidth=2,
            color=:orange,
        )
        plot!(
            axial, zmm,
            result.channel_temperature[
                result.core_group, :, end,
            ]; label="central channel orbit", linewidth=2,
            color=:blue,
        )
        scatter!(
            axial, [5.0, 58.0, 107.0],
            case.observed[end, 1:3];
            label="measured side TCs", color=:black,
            markersize=5,
        )
        savefig(
            axial,
            joinpath(
                axial_dir,
                "axial_profile_$(case.inputs.id)_2D_v15.png",
            ),
        )

        jmid = argmin(abs.(result.z_solid .- 58e-3))
        map_values = result.channel_temperature[:, jmid, end]
        channel_map = heatmap(
            _channel_map_2D_v15(result, map_values);
            xlabel="Channel index", ylabel="Channel index",
            colorbar_title="K", aspect_ratio=:equal,
            title=(
                "2D v15 $(case.inputs.id), z=58 mm\n" *
                "skin=$(round(result.skin_temperature[jmid,end],digits=1)) K"
            ),
        )
        savefig(
            channel_map,
            joinpath(
                map_dir,
                "channel_map_$(case.inputs.id)_2D_v15.png",
            ),
        )
    end
    foreach(println, aggregate_rows)
    println(
        "v15 radial signs: mid=$mid_sign/15, " *
        "deep=$deep_sign/15",
    )
    println("v15 accepted = $accepted")
    return (
        cases=cases, sensor_rows=sensor_rows,
        final_rows=final_rows, aggregate_rows=aggregate_rows,
        accepted=accepted, mid_sign=mid_sign, deep_sign=deep_sign,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_validation_2D_v15()
end
