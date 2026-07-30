# ============================================================================
# validate_2D_v16_staged.jl - frozen cooling fit, full nominal validation
# ============================================================================

using Statistics
using Plots

include("run_2D_v16.jl")

function _selected_cooling_2D_v16()
    nominal_path = joinpath(
        OUTPUT_DIR_2D_v16,
        "parameters_cooling_nominal_selected_2D_v16.txt",
    )
    screen_path = joinpath(
        OUTPUT_DIR_2D_v16,
        "parameters_cooling_selected_2D_v16.txt",
    )
    path = isfile(nominal_path) ? nominal_path : screen_path
    values = Dict{String,Float64}()
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

function _selected_power_2D_v16()
    path = joinpath(
        OUTPUT_DIR_2D_v16,
        "parameters_power_selected_2D_v16.txt",
    )
    isfile(path) || return (1.25, 1.2075, 0.8855)
    values = Dict{String,Float64}()
    for line in eachline(path)
        occursin("=", line) || continue
        key, raw = split(line, "="; limit=2)
        value = tryparse(Float64, raw)
        value === nothing || (values[key] = value)
    end
    return (
        values["power_scale_456"],
        values["power_scale_304"],
        values["power_scale_256"],
    )
end

function _phase_2D_v16(id, cooling)
    cooling && return "cooling_validation"
    return id in TRAIN_HEATING_2D_v16 ?
        "heating_training" : "heating_validation"
end

function _map_2D_v16(result, values)
    grid = V16.build_network_grid2D(result.parameters)
    return [
        values[grid.group_map[i, j]]
        for j in 1:grid.channel_side_count,
            i in 1:grid.channel_side_count
    ]
end

function run_validation_2D_v16(; max_points=121)
    mkpath(OUTPUT_DIR_2D_v16)
    plot_root = joinpath(OUTPUT_DIR_2D_v16, "plots")
    transient_dir = joinpath(plot_root, "transients")
    axial_dir = joinpath(plot_root, "axial_profiles")
    map_dir = joinpath(plot_root, "channel_maps")
    foreach(mkpath, (plot_root, transient_dir, axial_dir, map_dir))
    selected = _selected_cooling_2D_v16()
    power_scales = _selected_power_2D_v16()
    p = parameters_2D_v16(
        nominal_mesh=true,
        skin_felt_contact_scale=selected.contact,
        felt_conductivity_scale=selected.felt_k,
        felt_heat_capacity_scale=selected.felt_cp,
        power_scales=power_scales,
    )
    specs = vcat(
        [(id=id, cooling=false) for id in HEATING_IDS_2D_v16],
        [(id=id, cooling=true) for id in COOLING_IDS_2D_v16],
    )
    cases = Vector{Any}(undef, length(specs))
    Threads.@threads for index in eachindex(specs)
        spec = specs[index]
        input = case_inputs_2D_v16(
            spec.id; cooling=spec.cooling, max_points=max_points,
        )
        cases[index] = simulate_case_2D_v16(input, p)
    end

    sensors = NamedTuple[]
    finals = NamedTuple[]
    transients = NamedTuple[]
    for case in cases
        phase = _phase_2D_v16(case.inputs.id, case.inputs.cooling)
        for (index, sensor) in enumerate(SENSOR_NAMES_2D_v16)
            model = case.model[:, index]
            observed = case.observed[:, index]
            push!(sensors, (
                phase=phase, simulation_id=case.inputs.id,
                sensor=String(sensor),
                rmse_K=sqrt(mean(abs2, model .- observed)),
                steady_error_K=model[end] - observed[end],
                t90_error_s=V16.get_t90_2D(
                    case.inputs.times, model,
                ) - V16.get_t90_2D(
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
        for index in eachindex(case.inputs.times)
            push!(transients, (
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

    aggregates = NamedTuple[]
    for phase in unique(row.phase for row in sensors)
        phase_sensors = filter(row -> row.phase == phase, sensors)
        phase_finals = filter(row -> row.phase == phase, finals)
        push!(aggregates, (
            phase=phase,
            mean_sensor_rmse_K=mean(
                row.rmse_K for row in phase_sensors
            ),
            steady_mae_K=mean(
                abs(row.steady_error_K) for row in phase_sensors
            ),
            t90_mae_s=mean(
                abs(row.t90_error_s) for row in phase_sensors
            ),
            axial_rmse_K=sqrt(mean(
                (row.model_T12_minus_T8_K -
                 row.observed_T12_minus_T8_K)^2
                for row in phase_finals
            )),
            mid_radial_rmse_K=sqrt(mean(
                (row.model_T12_minus_T9_K -
                 row.observed_T12_minus_T9_K)^2
                for row in phase_finals
            )),
            deep_radial_rmse_K=sqrt(mean(
                (row.model_T11_minus_T10_K -
                 row.observed_T11_minus_T10_K)^2
                for row in phase_finals
            )),
            dp1_rmse_mbar=sqrt(mean(
                (row.dp1_model_mbar - row.dp1_observed_mbar)^2
                for row in phase_finals
            )),
            mean_core_to_side_flow_ratio=mean(
                row.core_to_side_flow_ratio for row in phase_finals
            ),
        ))
    end
    _write_namedtuples_csv_2D_v16(
        joinpath(OUTPUT_DIR_2D_v16,
                 "staged_sensor_metrics_2D_v16.csv"), sensors,
    )
    _write_namedtuples_csv_2D_v16(
        joinpath(OUTPUT_DIR_2D_v16,
                 "staged_final_profiles_2D_v16.csv"), finals,
    )
    _write_namedtuples_csv_2D_v16(
        joinpath(OUTPUT_DIR_2D_v16,
                 "staged_aggregate_2D_v16.csv"), aggregates,
    )
    _write_namedtuples_csv_2D_v16(
        joinpath(OUTPUT_DIR_2D_v16,
                 "staged_transients_2D_v16.csv"), transients,
    )

    heating_finals = filter(
        row -> row.phase != "cooling_validation", finals,
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
        row -> row.phase == "heating_validation", aggregates,
    ))
    cooling = only(filter(
        row -> row.phase == "cooling_validation", aggregates,
    ))
    cooling_sensor_rows = filter(
        row -> row.phase == "cooling_validation", sensors,
    )
    primary_cooling_rmse = mean(
        row.rmse_K for row in cooling_sensor_rows
        if row.sensor in ("T8", "T12", "T11", "T2")
    )
    t2_cooling_rmse = mean(
        row.rmse_K for row in cooling_sensor_rows
        if row.sensor == "T2"
    )
    cooling_fit_useful =
        primary_cooling_rmse <= 0.90 * 26.87 &&
        t2_cooling_rmse < 20.0
    parameter_identifiable = selected.felt_k < 7.20
    accepted =
        cooling_fit_useful &&
        mid_sign >= 12 && deep_sign >= 12 &&
        validation.mid_radial_rmse_K < 20.0 &&
        validation.deep_radial_rmse_K < 20.0 &&
        validation.mean_sensor_rmse_K <= 81.57 &&
        cooling.mean_sensor_rmse_K <= 37.73
    open(joinpath(
        OUTPUT_DIR_2D_v16,
        "staged_acceptance_status_2D_v16.txt",
    ), "w") do io
        println(io, "primary_cooling_rmse_K=", primary_cooling_rmse)
        println(io, "t2_cooling_rmse_K=", t2_cooling_rmse)
        println(io, "cooling_fit_useful=", cooling_fit_useful)
        println(
            io, "felt_conductivity_at_extended_bound=",
            !parameter_identifiable,
        )
        println(
            io, "parameter_identifiable=",
            parameter_identifiable,
        )
        println(io, "power_scale_456=", power_scales[1])
        println(io, "power_scale_304=", power_scales[2])
        println(io, "power_scale_256=", power_scales[3])
        println(io, "mid_radial_sign_correct=$mid_sign/15")
        println(io, "deep_radial_sign_correct=$deep_sign/15")
        foreach(row -> println(io, row), aggregates)
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
        title="2D v16 cooling + power-refitted parity",
        legend=:topleft, size=(850, 720),
    )
    all_values = Float64[]
    phase_markers = Dict(
        "heating_training" => :circle,
        "heating_validation" => :diamond,
        "cooling_validation" => :utriangle,
    )
    for phase in keys(phase_colors)
        xs = Float64[]
        ys = Float64[]
        for row in filter(r -> r.phase == phase, finals)
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
            marker=phase_markers[phase],
        )
    end
    lo, hi = extrema(all_values)
    plot!(
        parity, [lo, hi], [lo, hi]; label="1:1",
        color=:black, linestyle=:dash,
    )
    savefig(parity, joinpath(plot_root, "parity_staged_2D_v16.png"))

    colors = (:red, :orange, :purple, :blue, :cyan, :green, :black)
    for case in cases
        panels = Any[]
        for (index, sensor) in enumerate(SENSOR_NAMES_2D_v16)
            panel = plot(
                case.inputs.times, case.observed[:, index];
                label="measured", color=:black, linewidth=2,
                ylabel="$(sensor) (K)",
            )
            plot!(
                panel, case.inputs.times, case.model[:, index];
                label="v16", color=colors[index], linewidth=2,
                linestyle=:dash,
            )
            push!(panels, panel)
        end
        transient_plot = plot(
            panels...; layout=(4, 2), size=(1100, 1000),
            plot_title="2D v16 $(case.inputs.id)",
            xlabel="Time (s)",
        )
        savefig(
            transient_plot,
            joinpath(
                transient_dir,
                "transient_$(case.inputs.id)_2D_v16.png",
            ),
        )
        case.inputs.cooling && continue
        result = case.result
        zmm = 1000.0 .* result.z_solid
        axial = plot(
            zmm, result.skin_temperature[:, end];
            label="exterior SiC skin", linewidth=3, color=:red,
            xlabel="Axial position (mm)",
            ylabel="Final temperature (K)",
            title="2D v16 axial profile $(case.inputs.id)",
        )
        plot!(
            axial, zmm,
            result.channel_temperature[
                result.side_group, :, end,
            ]; label="outer channel orbit", color=:orange,
            linewidth=2,
        )
        plot!(
            axial, zmm,
            result.channel_temperature[
                result.core_group, :, end,
            ]; label="central channel orbit", color=:blue,
            linewidth=2,
        )
        scatter!(
            axial, [5.0, 58.0, 107.0],
            case.observed[end, 1:3];
            label="measured side TCs", color=:black,
        )
        savefig(
            axial,
            joinpath(
                axial_dir,
                "axial_profile_$(case.inputs.id)_2D_v16.png",
            ),
        )
        jmid = argmin(abs.(result.z_solid .- 58e-3))
        channel_map = heatmap(
            _map_2D_v16(
                result,
                result.channel_temperature[:, jmid, end],
            );
            xlabel="Channel index", ylabel="Channel index",
            colorbar_title="K", aspect_ratio=:equal,
            title=(
                "2D v16 $(case.inputs.id), z=58 mm\n" *
                "skin=$(round(result.skin_temperature[jmid,end],digits=1)) K"
            ),
        )
        savefig(
            channel_map,
            joinpath(
                map_dir,
                "channel_map_$(case.inputs.id)_2D_v16.png",
            ),
        )
    end
    foreach(println, aggregates)
    println(
        "v16 signs: mid=$mid_sign/15 deep=$deep_sign/15; " *
        "cooling useful=$cooling_fit_useful accepted=$accepted",
    )
    return (
        cases=cases, sensors=sensors, finals=finals,
        aggregates=aggregates, accepted=accepted,
        cooling_fit_useful=cooling_fit_useful,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_validation_2D_v16()
end
