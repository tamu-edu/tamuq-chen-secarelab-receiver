# ============================================================================
# validate_2D_v13_staged.jl - nominal-mesh held-out confirmation and plots
# ============================================================================

using Statistics
using Plots

include("run_2D_v13.jl")

function selected_parameters_2D_v13(;
    nominal_mesh=true,
    screen_mesh=false,
)
    base = inherited_parameters_2D_v13(
        nominal_mesh=nominal_mesh, screen_mesh=screen_mesh,
    )
    population = rebuild_population_2D_v13(
        base;
        side_fraction=0.025,
        side_gas_fraction=0.0125,
        active_side_G_per_m=100.0,
        side_felt_h=5.0,
    )
    return rebuild_heating_2D_v13(
        population;
        beta_opt=150.0,
        heat_scale=1.80,
        power_scales=(1.25, 1.2075, 0.77),
    )
end

function validation_phase_2D_v13(id, cooling)
    cooling && return "cooling_validation"
    return id in TRAIN_HEATING_2D_v13 ?
        "heating_training" : "heating_validation"
end

function run_validation_2D_v13(; max_points=121)
    mkpath(OUTPUT_DIR_2D_v13)
    plot_root = joinpath(OUTPUT_DIR_2D_v13, "plots")
    temporal_dir = joinpath(plot_root, "transients")
    axial_dir = joinpath(plot_root, "axial_profiles")
    mkpath(temporal_dir)
    mkpath(axial_dir)

    p = selected_parameters_2D_v13(; nominal_mesh=true)
    specs = vcat(
        [(id=id, cooling=false) for id in HEATING_IDS_2D_v13],
        [(id=id, cooling=true) for id in COOLING_IDS_2D_v13],
    )
    cases = Vector{Any}(undef, length(specs))
    Threads.@threads for index in eachindex(specs)
        spec = specs[index]
        inputs = case_inputs_2D_v13(
            spec.id; cooling=spec.cooling, max_points=max_points,
        )
        cases[index] = simulate_case_2D_v13(inputs, p)
    end

    sensor_rows = NamedTuple[]
    final_rows = NamedTuple[]
    transient_rows = NamedTuple[]
    for case in cases
        phase = validation_phase_2D_v13(
            case.inputs.id, case.inputs.cooling,
        )
        for (index, sensor) in enumerate(SENSOR_NAMES_2D_v13)
            model = case.model[:, index]
            observed = case.observed[:, index]
            push!(sensor_rows, (
                phase=phase,
                simulation_id=case.inputs.id,
                sensor=String(sensor),
                rmse_K=sqrt(mean(abs2, model .- observed)),
                steady_error_K=model[end] - observed[end],
                t90_error_s=get_t90_2D(case.inputs.times, model) -
                    get_t90_2D(case.inputs.times, observed),
            ))
        end
        m = case.model[end, :]
        e = case.observed[end, :]
        push!(final_rows, (
            phase=phase,
            simulation_id=case.inputs.id,
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
            dp1_model_mbar=case.result.dp1_prediction_mbar[end],
            dp1_observed_mbar=case.inputs.dp1[end],
            adaptor_T_K=case.result.adaptor_temperature[end],
            rear_housing_T_K=
                case.result.housing_extension_temperature[end],
            side_front_T_K=case.result.side_temperature[1, end],
            side_mid_T_K=case.result.side_temperature[
                argmin(abs.(case.result.z_solid .- 58e-3)), end
            ],
            side_deep_T_K=case.result.side_temperature[
                argmin(abs.(case.result.z_solid .- 107e-3)), end
            ],
        ))
        for index in eachindex(case.inputs.times)
            push!(transient_rows, (
                phase=phase,
                simulation_id=case.inputs.id,
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
            steady_mae_K=mean(abs(row.steady_error_K)
                              for row in sensors),
            t90_mae_s=mean(abs(row.t90_error_s)
                           for row in sensors),
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
        ))
    end

    _write_namedtuples_csv_2D_v13(
        joinpath(
            OUTPUT_DIR_2D_v13,
            "staged_sensor_metrics_2D_v13.csv",
        ), sensor_rows,
    )
    _write_namedtuples_csv_2D_v13(
        joinpath(
            OUTPUT_DIR_2D_v13,
            "staged_final_profiles_2D_v13.csv",
        ), final_rows,
    )
    _write_namedtuples_csv_2D_v13(
        joinpath(
            OUTPUT_DIR_2D_v13,
            "staged_aggregate_2D_v13.csv",
        ), aggregate_rows,
    )
    _write_namedtuples_csv_2D_v13(
        joinpath(
            OUTPUT_DIR_2D_v13,
            "staged_transients_2D_v13.csv",
        ), transient_rows,
    )

    colors = Dict(
        "heating_training" => :royalblue,
        "heating_validation" => :darkorange,
        "cooling_validation" => :seagreen,
    )
    parity = plot(
        xlabel="Observed final temperature (K)",
        ylabel="Modeled final temperature (K)",
        title="2D v13 staged nominal-mesh parity",
        legend=:topleft, size=(850, 720),
    )
    all_observed = Float64[]
    all_modeled = Float64[]
    for phase in keys(colors)
        xs = Float64[]
        ys = Float64[]
        for row in filter(r -> r.phase == phase, final_rows)
            for sensor in (
                "T8", "T12", "T11", "T9", "T10", "T3", "T2",
            )
                push!(xs, getproperty(
                    row, Symbol("observed_$(sensor)_K"),
                ))
                push!(ys, getproperty(
                    row, Symbol("model_$(sensor)_K"),
                ))
            end
        end
        append!(all_observed, xs)
        append!(all_modeled, ys)
        scatter!(
            parity, xs, ys; label=replace(phase, "_" => " "),
            color=colors[phase], alpha=0.8,
        )
    end
    lo = min(250.0, minimum(vcat(all_observed, all_modeled)))
    hi = max(1200.0, maximum(vcat(all_observed, all_modeled)))
    plot!(
        parity, [lo, hi], [lo, hi]; label="1:1",
        color=:black, linestyle=:dash,
    )
    savefig(
        parity,
        joinpath(plot_root, "parity_staged_2D_v13.png"),
    )

    sensor_colors = (
        :red, :orange, :purple, :blue, :cyan, :green, :black,
    )
    for case in cases
        panels = Any[]
        for (index, sensor) in enumerate(SENSOR_NAMES_2D_v13)
            panel = plot(
                case.inputs.times, case.observed[:, index];
                label="experiment", color=:black, linewidth=2,
                ylabel="$(sensor) (K)",
            )
            plot!(
                panel, case.inputs.times, case.model[:, index];
                label="v13", color=sensor_colors[index],
                linewidth=2, linestyle=:dash,
            )
            push!(panels, panel)
        end
        combined = plot(
            panels...; layout=(4, 2), size=(1100, 1200),
            plot_title="2D v13 $(case.inputs.id)",
        )
        savefig(combined, joinpath(
            temporal_dir,
            "transient_$(case.inputs.id)_2D_v13.png",
        ))
    end

    for case in filter(c -> !c.inputs.cooling, cases)
        m = case.model[end, :]
        e = case.observed[end, :]
        axial = plot(
            [5.0, 58.0, 107.0], m[1:3];
            marker=:circle, linewidth=2, label="v13 side",
            xlabel="Axial position from front (mm)",
            ylabel="Final temperature (K)",
            title="2D v13 axial profile $(case.inputs.id)",
        )
        plot!(
            axial, [5.0, 58.0, 107.0], e[1:3];
            marker=:diamond, linewidth=2, linestyle=:dash,
            label="experiment side",
        )
        plot!(
            axial, [58.0, 107.0], m[4:5];
            marker=:circle, linewidth=2,
            label="v13 interior TC",
        )
        plot!(
            axial, [58.0, 107.0], e[4:5];
            marker=:diamond, linewidth=2, linestyle=:dash,
            label="experiment interior",
        )
        savefig(axial, joinpath(
            axial_dir,
            "axial_profile_$(case.inputs.id)_2D_v13.png",
        ))
    end

    heating_finals = filter(
        row -> startswith(row.phase, "heating"), final_rows,
    )
    mid_sign = count(
        row -> sign(row.model_T12_minus_T9_K) ==
               sign(row.observed_T12_minus_T9_K),
        heating_finals,
    )
    deep_sign = count(
        row -> sign(row.model_T11_minus_T10_K) ==
               sign(row.observed_T11_minus_T10_K),
        heating_finals,
    )
    status = (
        mid_sign == 15 && deep_sign == 15 &&
        all(row.mean_sensor_rmse_K < 75.0
            for row in aggregate_rows if
            startswith(row.phase, "heating"))
    ) ? "CANDIDATE_FOR_COEFFICIENT_EXTRACTION" :
        "REJECTED_FOR_COEFFICIENT_EXTRACTION"
    open(joinpath(
        OUTPUT_DIR_2D_v13,
        "staged_acceptance_status_2D_v13.txt",
    ), "w") do io
        println(io, "mid_radial_sign_correct=$mid_sign/15")
        println(io, "deep_radial_sign_correct=$deep_sign/15")
        for row in aggregate_rows
            println(io, row)
        end
        println(io, "status=$status")
    end
    println("mid radial sign = $mid_sign/15")
    println("deep radial sign = $deep_sign/15")
    foreach(println, aggregate_rows)
    return (
        parameters=p, cases=cases, sensor_rows=sensor_rows,
        final_rows=final_rows, aggregate_rows=aggregate_rows,
        mid_sign=mid_sign, deep_sign=deep_sign, status=status,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_validation_2D_v13()
end

