using CSV
using DataFrames
using Plots
using Statistics

const ROOT_2D_v11 = @__DIR__
const INPUT_2D_v11 = joinpath(ROOT_2D_v11, "summaries", "2D_v11")
const PLOT_DIR_2D_v11 = joinpath(INPUT_2D_v11, "plots")
const FOCUS_VARIANTS_2D_v11 = (
    "apparent_prescribed_v9",
    "graetz_fd361_prescribed",
    "graetz_fd361_equal_open",
    "graetz_fd361_equal_groove",
)
const VARIANT_LABELS_2D_v11 = Dict(
    "apparent_prescribed_v9" => "v9 apparent",
    "graetz_fd361_prescribed" => "Graetz, fixed flow",
    "graetz_fd298_prescribed" => "Graetz Nu∞=2.98",
    "graetz_capuano366_prescribed" => "Capuano",
    "graetz_fd361_equal_open" => "Graetz, equal Δp",
    "graetz_fd361_equal_groove" => "Graetz, equal Δp + groove",
)

label_for_2D_v11(name) = get(VARIANT_LABELS_2D_v11, name, name)

function plot_cold_dp_2D_v11()
    cases = CSV.read(
        joinpath(INPUT_2D_v11, "cold_dp1_cases_2D_v11.csv"),
        DataFrame,
    )
    profile = CSV.read(
        joinpath(INPUT_2D_v11, "cold_groove_profile_2D_v11.csv"),
        DataFrame,
    )
    p1 = plot(
        xlabel="Corrected MFC flow (standard L/min)",
        ylabel="DP1 (mbar)",
        title="Cold t0 hydraulic decomposition",
        legend=:topleft,
    )
    observed = unique(
        select(cases, :simulation_id, :flow_lpm, :observed_dp1_mbar),
    )
    scatter!(
        p1, observed.flow_lpm, observed.observed_dp1_mbar;
        label="observed", color=:black, markersize=5,
    )
    for (model, color) in (("no_groove", :blue), ("groove", :red))
        sub = sort(filter(:model => ==(model), cases), :flow_lpm)
        plot!(
            p1, sub.flow_lpm, sub.predicted_dp1_mbar;
            label=replace(model, "_" => " "), color=color, linewidth=2,
        )
    end
    p2 = plot(
        profile.groove_K, profile.rmse_mbar;
        marker=:circle,
        xlabel="Combined groove coefficient K",
        ylabel="Cold DP1 RMSE (mbar)",
        title="Groove profile likelihood",
        legend=false,
    )
    fig = plot(
        p1, p2;
        layout=(1, 2),
        size=(1100, 430),
        left_margin=7 * Plots.mm,
        bottom_margin=7 * Plots.mm,
    )
    savefig(fig, joinpath(PLOT_DIR_2D_v11, "cold_dp1_groove_fit_2D_v11.png"))
end

function parity_data_2D_v11(cases)
    sensors = ("T8", "T12", "T11", "T9", "T10", "T3", "T2")
    rows = NamedTuple[]
    for row in eachrow(cases), sensor in sensors
        push!(rows, (
            variant=row.variant,
            sensor=sensor,
            observed=Float64(row[Symbol("observed_$(sensor)_K")]),
            model=Float64(row[Symbol("model_$(sensor)_K")]),
        ))
    end
    return DataFrame(rows)
end

function plot_parity_2D_v11(cases)
    parity = parity_data_2D_v11(cases)
    panels = Any[]
    colors = [:black, :blue, :orange, :red]
    for (variant, color) in zip(FOCUS_VARIANTS_2D_v11, colors)
        sub = filter(:variant => ==(variant), parity)
        lo = min(minimum(sub.observed), minimum(sub.model))
        hi = max(maximum(sub.observed), maximum(sub.model))
        p = scatter(
            sub.observed, sub.model;
            group=sub.sensor,
            xlabel="Observed temperature (K)",
            ylabel="Model temperature (K)",
            title=label_for_2D_v11(variant),
            markersize=3,
            legend=:topleft,
        )
        plot!(p, [lo, hi], [lo, hi]; color=:gray, linestyle=:dash, label="")
        push!(panels, p)
    end
    fig = plot(panels...; layout=(2, 2), size=(1050, 850))
    savefig(fig, joinpath(PLOT_DIR_2D_v11, "steady_parity_2D_v11.png"))
end

function plot_flow_offsets_2D_v11(cases)
    panels = Any[]
    for (irradiance, ids) in ((456.0, 67:71), (304.0, 72:76), (256.0, 77:81))
        p = plot(
            xlabel="Mean flow (standard L/min)",
            ylabel="T12 - T8 (K)",
            title="$(Int(irradiance)) kW/m²",
            legend=:topright,
        )
        observed = filter(:irradiance_kW_m2 => ==(irradiance), cases)
        observed = unique(
            select(
                observed,
                :simulation_id,
                :mean_flow_lpm,
                :observed_T12_minus_T8_K,
            ),
        )
        sort!(observed, :mean_flow_lpm)
        scatter!(
            p,
            observed.mean_flow_lpm,
            observed.observed_T12_minus_T8_K;
            label="observed", color=:black, markersize=5,
        )
        for (variant, color) in zip(
            FOCUS_VARIANTS_2D_v11, [:gray, :blue, :orange, :red],
        )
            sub = filter(
                row -> row.variant == variant &&
                       row.irradiance_kW_m2 == irradiance,
                cases,
            )
            sort!(sub, :mean_flow_lpm)
            plot!(
                p, sub.mean_flow_lpm, sub.model_T12_minus_T8_K;
                marker=:circle,
                color=color,
                linewidth=2,
                label=label_for_2D_v11(variant),
            )
        end
        push!(panels, p)
    end
    fig = plot(panels...; layout=(1, 3), size=(1450, 430))
    savefig(fig, joinpath(PLOT_DIR_2D_v11, "axial_offset_flow_2D_v11.png"))
end

function plot_radial_and_flow_2D_v11(cases)
    groove = filter(
        :variant => ==("graetz_fd361_equal_groove"),
        cases,
    )
    p1 = scatter(
        groove.observed_T12_minus_T9_K,
        groove.model_T12_minus_T9_K;
        xlabel="Observed T12-T9 (K)",
        ylabel="Model T12-T9 (K)",
        title="Mid-depth radial offset",
        legend=false,
    )
    p2 = scatter(
        groove.observed_T11_minus_T10_K,
        groove.model_T11_minus_T10_K;
        xlabel="Observed T11-T10 (K)",
        ylabel="Model T11-T10 (K)",
        title="Deep radial offset",
        legend=false,
    )
    p3 = scatter(
        groove.mean_flow_lpm,
        groove.core_to_edge_mass_flux_ratio;
        group=groove.irradiance_kW_m2,
        xlabel="Mean flow (standard L/min)",
        ylabel="Core/edge mass-flux ratio",
        title="Groove-induced maldistribution",
        legend=:topleft,
    )
    for p in (p1, p2)
        lo = -20.0
        hi = 80.0
        plot!(p, [lo, hi], [lo, hi]; color=:gray, linestyle=:dash)
    end
    fig = plot(p1, p2, p3; layout=(1, 3), size=(1400, 420))
    savefig(fig, joinpath(PLOT_DIR_2D_v11, "radial_offsets_and_flow_2D_v11.png"))
end

function plot_axial_profiles_2D_v11(axial, cases)
    for id in ("E67", "E72", "E77")
        panels = Any[]
        case_rows = filter(:simulation_id => ==(id), cases)
        observed = first(case_rows)
        for variant in FOCUS_VARIANTS_2D_v11
            sub = filter(
                row -> row.variant == variant && row.simulation_id == id,
                axial,
            )
            p = plot(
                sub.z_mm, sub.core_solid_K;
                label="solid core", color=:red, linewidth=2,
                xlabel="Depth (mm)", ylabel="Temperature (K)",
                title=label_for_2D_v11(variant),
            )
            plot!(p, sub.z_mm, sub.edge_solid_K;
                  label="solid edge", color=:orange, linewidth=2)
            plot!(p, sub.z_mm, sub.core_gas_K;
                  label="gas core", color=:blue, linestyle=:dash)
            plot!(p, sub.z_mm, sub.edge_gas_K;
                  label="gas edge", color=:cyan, linestyle=:dash)
            scatter!(
                p, [5.0, 58.0, 107.0],
                [observed.observed_T8_K,
                 observed.observed_T9_K,
                 observed.observed_T10_K];
                label="observed core", color=:black, marker=:circle,
            )
            scatter!(
                p, [58.0, 107.0],
                [observed.observed_T12_K,
                 observed.observed_T11_K];
                label="observed side", color=:black, marker=:diamond,
            )
            push!(panels, p)
        end
        fig = plot(panels...; layout=(2, 2), size=(1150, 850))
        savefig(
            fig,
            joinpath(PLOT_DIR_2D_v11, "axial_profile_$(id)_2D_v11.png"),
        )
    end
end

function plot_transients_2D_v11(transients)
    sensors = ("T8", "T12", "T9", "T10", "T3")
    for id in ("E67", "E72", "E77")
        panels = Any[]
        for sensor in sensors
            p = plot(
                xlabel="Time (s)",
                ylabel="$(sensor) (K)",
                title="$id $sensor",
                legend=:bottomright,
            )
            observed = filter(
                row -> row.simulation_id == id &&
                       row.variant == first(FOCUS_VARIANTS_2D_v11),
                transients,
            )
            plot!(
                p, observed.time_s,
                observed[!, Symbol("observed_$(sensor)_K")];
                label="observed", color=:black, linewidth=2,
            )
            for (variant, color) in zip(
                FOCUS_VARIANTS_2D_v11, [:gray, :blue, :orange, :red],
            )
                sub = filter(
                    row -> row.simulation_id == id &&
                           row.variant == variant,
                    transients,
                )
                plot!(
                    p, sub.time_s,
                    sub[!, Symbol("model_$(sensor)_K")];
                    label=label_for_2D_v11(variant),
                    color=color,
                    linewidth=1.5,
                )
            end
            push!(panels, p)
        end
        fig = plot(panels...; layout=(3, 2), size=(1150, 1100))
        savefig(
            fig,
            joinpath(PLOT_DIR_2D_v11, "transient_$(id)_2D_v11.png"),
        )
    end
end

function plot_ring_profiles_2D_v11(rings)
    panels = Any[]
    for id in ("E67", "E72", "E77")
        p = plot(
            xlabel="Radius (mm)",
            ylabel="Ring mass flux (kg/m²/s)",
            title=id,
            legend=:topright,
        )
        for (variant, color) in (
            ("graetz_fd361_equal_open", :orange),
            ("graetz_fd361_equal_groove", :red),
        )
            sub = filter(
                row -> row.variant == variant &&
                       row.simulation_id == id,
                rings,
            )
            plot!(
                p, sub.radius_mm, sub.ring_mass_flux_kg_m2s;
                marker=:circle, color=color,
                label=label_for_2D_v11(variant),
            )
        end
        vline!(p, [6.5]; color=:black, linestyle=:dash,
               label="free radius")
        push!(panels, p)
    end
    fig = plot(panels...; layout=(1, 3), size=(1400, 420))
    savefig(fig, joinpath(PLOT_DIR_2D_v11, "ring_flow_profiles_2D_v11.png"))
end

function plot_inheritance_sensitivity_2D_v11()
    path = joinpath(
        INPUT_2D_v11,
        "inheritance_sensitivity_summary_2D_v11.csv",
    )
    isfile(path) || return
    summary = CSV.read(path, DataFrame)
    x = collect(1:nrow(summary))
    labels = summary.variant
    p1 = scatter(
        x, summary.axial_rmse_K;
        marker=:circle,
        xlabel="Sensitivity variant",
        ylabel="Endpoint axial RMSE (K)",
        title="Inherited-parameter axial response",
        xticks=(x, labels),
        xrotation=45,
        legend=false,
    )
    p2 = plot(
        x, summary.slope_456_K_per_lpm;
        marker=:circle,
        label="456",
        xlabel="Sensitivity variant",
        ylabel="Endpoint slope (K/(L/min))",
        title="Axial flow-slope signs",
        xticks=(x, labels),
        xrotation=45,
    )
    plot!(
        p2, x, summary.slope_304_K_per_lpm;
        marker=:circle, label="304",
    )
    plot!(
        p2, x, summary.slope_256_K_per_lpm;
        marker=:circle, label="256",
    )
    hline!(p2, [0.0]; color=:black, linestyle=:dash, label="")
    fig = plot(
        p1, p2;
        layout=(2, 1),
        size=(1250, 850),
        bottom_margin=14 * Plots.mm,
        left_margin=5 * Plots.mm,
    )
    savefig(
        fig,
        joinpath(
            PLOT_DIR_2D_v11,
            "inheritance_sensitivity_2D_v11.png",
        ),
    )
end

function plot_betaopt_confirmation_2D_v11(cases)
    path = joinpath(
        INPUT_2D_v11,
        "betaopt110_confirmation_cases_2D_v11.csv",
    )
    isfile(path) || return
    confirmation = CSV.read(path, DataFrame)
    panels = Any[]
    for irradiance in (456.0, 304.0, 256.0)
        p = plot(
            xlabel="Mean flow (standard L/min)",
            ylabel="T12 - T8 (K)",
            title="$(Int(irradiance)) kW/m^2",
            legend=:topright,
        )
        original = filter(
            row -> row.variant == "graetz_fd361_equal_groove" &&
                   row.irradiance_kW_m2 == irradiance,
            cases,
        )
        revised = filter(
            :irradiance_kW_m2 => ==(irradiance),
            confirmation,
        )
        sort!(original, :mean_flow_lpm)
        sort!(revised, :mean_flow_lpm)
        scatter!(
            p,
            original.mean_flow_lpm,
            original.observed_T12_minus_T8_K;
            color=:black,
            label="observed",
        )
        plot!(
            p,
            original.mean_flow_lpm,
            original.model_T12_minus_T8_K;
            marker=:circle,
            color=:red,
            label="v9-fitted beta_opt=21.3",
        )
        plot!(
            p,
            revised.mean_flow_lpm,
            revised.model_T12_minus_T8_K;
            marker=:circle,
            color=:purple,
            label="beta_opt=110",
        )
        push!(panels, p)
    end
    fig = plot(panels...; layout=(1, 3), size=(1450, 430))
    savefig(
        fig,
        joinpath(
            PLOT_DIR_2D_v11,
            "betaopt110_flow_confirmation_2D_v11.png",
        ),
    )
end

function main_plot_2D_v11()
    mkpath(PLOT_DIR_2D_v11)
    cases = CSV.read(
        joinpath(INPUT_2D_v11, "model_form_cases_2D_v11.csv"),
        DataFrame,
    )
    axial = CSV.read(
        joinpath(INPUT_2D_v11, "axial_profiles_2D_v11.csv"),
        DataFrame,
    )
    transients = CSV.read(
        joinpath(INPUT_2D_v11, "representative_transients_2D_v11.csv"),
        DataFrame,
    )
    rings = CSV.read(
        joinpath(INPUT_2D_v11, "ring_profiles_2D_v11.csv"),
        DataFrame,
    )
    plot_cold_dp_2D_v11()
    plot_parity_2D_v11(cases)
    plot_flow_offsets_2D_v11(cases)
    plot_radial_and_flow_2D_v11(cases)
    plot_axial_profiles_2D_v11(axial, cases)
    plot_transients_2D_v11(transients)
    plot_ring_profiles_2D_v11(rings)
    plot_inheritance_sensitivity_2D_v11()
    plot_betaopt_confirmation_2D_v11(cases)
    println("v11 plots written to $PLOT_DIR_2D_v11")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_plot_2D_v11()
end
