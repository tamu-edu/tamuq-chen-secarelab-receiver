ENV["GKSwstype"] = "100"

using CSV
using DataFrames
using Plots
using Statistics

const STAGED_DIR_2D_v11 = joinpath(
    @__DIR__, "summaries", "2D_v11",
)
const STAGED_PLOT_DIR_2D_v11 = joinpath(
    STAGED_DIR_2D_v11, "plots", "staged",
)
const STAGED_GROUPS_2D_v11 = (
    (456.0, ("E67", "E68", "E69", "E70", "E71")),
    (304.0, ("E72", "E73", "E74", "E75", "E76")),
    (256.0, ("E77", "E78", "E79", "E80", "E81")),
)

function staged_linear_fit_2D_v11(x, y)
    xmean = mean(x)
    ymean = mean(y)
    slope = sum(
        (x[i] - xmean) * (y[i] - ymean)
        for i in eachindex(x)
    ) / sum((value - xmean)^2 for value in x)
    return (slope=slope, correlation=cor(x, y))
end

function staged_flow_slopes_2D_v11(cases)
    rows = NamedTuple[]
    for (irradiance, ids) in STAGED_GROUPS_2D_v11
        for phase in (
            "heating_training",
            "heating_validation",
            "all_heating",
        )
            selected = filter(
                row -> row.simulation_id in ids &&
                       (
                           phase == "all_heating" ||
                           row.phase == phase
                       ),
                cases,
            )
            nrow(selected) >= 2 || continue
            model = staged_linear_fit_2D_v11(
                selected.mean_flow_lpm,
                selected.model_T12_minus_T8_K,
            )
            observed = staged_linear_fit_2D_v11(
                selected.mean_flow_lpm,
                selected.observed_T12_minus_T8_K,
            )
            push!(rows, (
                phase=phase,
                irradiance_kW_m2=irradiance,
                n=nrow(selected),
                model_slope_K_per_lpm=model.slope,
                observed_slope_K_per_lpm=observed.slope,
                model_correlation=model.correlation,
                observed_correlation=observed.correlation,
            ))
        end
    end
    return DataFrame(rows)
end

function staged_mesh_comparison_2D_v11(mesh)
    sensors = (:T8_K, :T12_K, :T11_K, :T9_K, :T10_K, :T3_K, :T2_K)
    rows = NamedTuple[]
    for id in unique(mesh.simulation_id)
        coarse = only(filter(
            row -> row.mesh == "sensitivity" &&
                   row.simulation_id == id,
            mesh,
        ))
        nominal = only(filter(
            row -> row.mesh == "nominal" &&
                   row.simulation_id == id,
            mesh,
        ))
        deltas = [
            nominal[sensor] - coarse[sensor] for sensor in sensors
        ]
        push!(rows, (
            simulation_id=id,
            max_abs_sensor_delta_K=maximum(abs, deltas),
            rms_sensor_delta_K=sqrt(mean(abs2, deltas)),
            axial_offset_delta_K=
                (nominal.T12_K - nominal.T8_K) -
                (coarse.T12_K - coarse.T8_K),
            mid_radial_offset_delta_K=
                (nominal.T12_K - nominal.T9_K) -
                (coarse.T12_K - coarse.T9_K),
            dp1_delta_mbar=nominal.dp1_mbar - coarse.dp1_mbar,
        ))
    end
    return DataFrame(rows)
end

function plot_staged_traces_2D_v11(profile, power)
    p1 = scatter(
        profile.evaluation,
        profile.objective;
        xlabel="True-model evaluation",
        ylabel="Profile objective",
        title="Transport response-surface evaluations",
        marker_z=profile.beta_opt,
        colorbar_title="beta_opt",
        legend=false,
    )
    p2 = scatter(
        power.tested_scale,
        power.objective;
        group=power.group,
        xlabel="Power scale",
        ylabel="Power-stage objective",
        title="Separate irradiance-group profiles",
        legend=:topright,
    )
    fig = plot(
        p1, p2;
        layout=(1, 2),
        size=(1200, 430),
        bottom_margin=6 * Plots.mm,
        left_margin=4 * Plots.mm,
    )
    savefig(
        fig,
        joinpath(STAGED_PLOT_DIR_2D_v11, "optimization_traces_2D_v11.png"),
    )
end

function parity_rows_staged_2D_v11(cases)
    rows = NamedTuple[]
    for row in eachrow(cases), sensor in (
        "T8", "T12", "T11", "T9", "T10", "T3", "T2",
    )
        push!(rows, (
            phase=row.phase,
            simulation_id=row.simulation_id,
            sensor=sensor,
            observed=Float64(
                row[Symbol("observed_$(sensor)_K")],
            ),
            model=Float64(row[Symbol("model_$(sensor)_K")]),
        ))
    end
    return DataFrame(rows)
end

function plot_staged_parity_2D_v11(cases)
    parity = parity_rows_staged_2D_v11(cases)
    panels = Any[]
    titles = Dict(
        "heating_training" => "Heating training",
        "heating_validation" => "Held-out heating",
        "cooling_validation" => "Cooling validation",
    )
    for phase in (
        "heating_training",
        "heating_validation",
        "cooling_validation",
    )
        sub = filter(:phase => ==(phase), parity)
        lower = min(minimum(sub.observed), minimum(sub.model))
        upper = max(maximum(sub.observed), maximum(sub.model))
        p = scatter(
            sub.observed,
            sub.model;
            group=sub.sensor,
            xlabel="Observed temperature (K)",
            ylabel="Model temperature (K)",
            title=titles[phase],
            legend=phase == "heating_training" ? :topleft : false,
            markersize=4,
        )
        plot!(
            p, [lower, upper], [lower, upper];
            color=:gray, linestyle=:dash, label="",
        )
        push!(panels, p)
    end
    fig = plot(panels...; layout=(1, 3), size=(1450, 440))
    savefig(
        fig,
        joinpath(STAGED_PLOT_DIR_2D_v11, "validation_parity_2D_v11.png"),
    )
end

function plot_staged_flow_profiles_2D_v11(cases)
    panels = Any[]
    for (irradiance, ids) in STAGED_GROUPS_2D_v11
        selected = filter(
            row -> row.simulation_id in ids,
            cases,
        )
        sort!(selected, :mean_flow_lpm)
        p = plot(
            xlabel="Mean flow (standard L/min)",
            ylabel="T12 - T8 (K)",
            title="$(Int(irradiance)) kW/m^2",
            legend=:topleft,
        )
        plot!(
            p,
            selected.mean_flow_lpm,
            selected.observed_T12_minus_T8_K;
            marker=:circle,
            color=:black,
            label="observed",
        )
        plot!(
            p,
            selected.mean_flow_lpm,
            selected.model_T12_minus_T8_K;
            marker=:circle,
            color=:purple,
            label="staged v11",
        )
        training = filter(
            :phase => ==("heating_training"),
            selected,
        )
        validation = filter(
            :phase => ==("heating_validation"),
            selected,
        )
        scatter!(
            p,
            training.mean_flow_lpm,
            training.model_T12_minus_T8_K;
            marker=:circle,
            color=:purple,
            label="training model",
        )
        scatter!(
            p,
            validation.mean_flow_lpm,
            validation.model_T12_minus_T8_K;
            marker=:diamond,
            color=:orange,
            label="held-out model",
        )
        push!(panels, p)
    end
    fig = plot(panels...; layout=(1, 3), size=(1450, 430))
    savefig(
        fig,
        joinpath(STAGED_PLOT_DIR_2D_v11, "axial_flow_profiles_2D_v11.png"),
    )
end

function plot_staged_radial_2D_v11(cases)
    heating = filter(
        row -> row.phase != "cooling_validation",
        cases,
    )
    p1 = scatter(
        heating.observed_T12_minus_T9_K,
        heating.model_T12_minus_T9_K;
        group=heating.phase,
        xlabel="Observed T12-T9 (K)",
        ylabel="Model T12-T9 (K)",
        title="Mid-depth radial offset",
    )
    p2 = scatter(
        heating.observed_T11_minus_T10_K,
        heating.model_T11_minus_T10_K;
        group=heating.phase,
        xlabel="Observed T11-T10 (K)",
        ylabel="Model T11-T10 (K)",
        title="Deep radial offset",
    )
    for p in (p1, p2)
        plot!(
            p, [-10.0, 80.0], [-10.0, 80.0];
            color=:gray, linestyle=:dash, label="",
        )
    end
    fig = plot(p1, p2; layout=(1, 2), size=(1050, 430))
    savefig(
        fig,
        joinpath(STAGED_PLOT_DIR_2D_v11, "radial_validation_2D_v11.png"),
    )
end

function plot_staged_transients_2D_v11(transients)
    for id in unique(transients.simulation_id)
        selected = filter(:simulation_id => ==(id), transients)
        panels = Any[]
        for sensor in ("T8", "T12", "T10", "T3")
            p = plot(
                selected.time_s,
                selected[!, Symbol("observed_$(sensor)_K")];
                label="observed",
                color=:black,
                linewidth=2,
                xlabel="Time (s)",
                ylabel="$sensor (K)",
                title="$id $sensor",
            )
            plot!(
                p,
                selected.time_s,
                selected[!, Symbol("model_$(sensor)_K")];
                label="staged v11",
                color=:purple,
                linewidth=2,
            )
            push!(panels, p)
        end
        fig = plot(panels...; layout=(2, 2), size=(1050, 780))
        savefig(
            fig,
            joinpath(
                STAGED_PLOT_DIR_2D_v11,
                "transient_$(id)_2D_v11.png",
            ),
        )
    end
end

function plot_staged_sensor_rmse_2D_v11(sensor_metrics)
    aggregated = combine(
        groupby(sensor_metrics, [:phase, :sensor]),
        :rmse_K => mean => :mean_rmse_K,
    )
    sensors = ["T8", "T12", "T11", "T9", "T10", "T3", "T2"]
    p = plot(
        xlabel="Sensor",
        ylabel="Mean transient RMSE (K)",
        title="Staged v11 error by phase and sensor",
        xticks=(1:length(sensors), sensors),
    )
    for (phase, color) in (
        ("heating_training", :blue),
        ("heating_validation", :orange),
        ("cooling_validation", :green),
    )
        sub = filter(:phase => ==(phase), aggregated)
        values = [
            only(sub.mean_rmse_K[sub.sensor .== sensor])
            for sensor in sensors
        ]
        plot!(
            p, 1:length(sensors), values;
            marker=:circle, linewidth=2,
            label=replace(phase, "_" => " "),
            color=color,
        )
    end
    savefig(
        p,
        joinpath(STAGED_PLOT_DIR_2D_v11, "sensor_rmse_2D_v11.png"),
    )
end

function main_plot_staged_2D_v11()
    mkpath(STAGED_PLOT_DIR_2D_v11)
    profile = CSV.read(
        joinpath(STAGED_DIR_2D_v11, "staged_profile_trace_2D_v11.csv"),
        DataFrame,
    )
    power = CSV.read(
        joinpath(STAGED_DIR_2D_v11, "staged_power_trace_2D_v11.csv"),
        DataFrame,
    )
    cases = CSV.read(
        joinpath(STAGED_DIR_2D_v11, "staged_case_metrics_2D_v11.csv"),
        DataFrame,
    )
    sensors = CSV.read(
        joinpath(STAGED_DIR_2D_v11, "staged_sensor_metrics_2D_v11.csv"),
        DataFrame,
    )
    transients = CSV.read(
        joinpath(STAGED_DIR_2D_v11, "staged_transients_2D_v11.csv"),
        DataFrame,
    )
    mesh = CSV.read(
        joinpath(STAGED_DIR_2D_v11, "staged_mesh_cases_2D_v11.csv"),
        DataFrame,
    )

    slopes = staged_flow_slopes_2D_v11(cases)
    CSV.write(
        joinpath(STAGED_DIR_2D_v11, "staged_flow_slopes_2D_v11.csv"),
        slopes,
    )
    mesh_comparison = staged_mesh_comparison_2D_v11(mesh)
    CSV.write(
        joinpath(
            STAGED_DIR_2D_v11,
            "staged_mesh_comparison_2D_v11.csv",
        ),
        mesh_comparison,
    )

    plot_staged_traces_2D_v11(profile, power)
    plot_staged_parity_2D_v11(cases)
    plot_staged_flow_profiles_2D_v11(cases)
    plot_staged_radial_2D_v11(cases)
    plot_staged_transients_2D_v11(transients)
    plot_staged_sensor_rmse_2D_v11(sensors)
    println("Staged v11 plots and derived tables written.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_plot_staged_2D_v11()
end
