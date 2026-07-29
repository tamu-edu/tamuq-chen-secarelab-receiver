# ============================================================================
# plot_2D_v9_full.jl - Publication-oriented plots for the completed v9 test
# ============================================================================

using CSV
using DataFrames
using Statistics
using StatsPlots

include("run_2D_v9_full.jl")

const PLOT_ROOT_2D_v9 = joinpath(
    @__DIR__, "summaries", "2D_v9", "plots",
)
const PARITY_PLOT_DIR_2D_v9 = joinpath(PLOT_ROOT_2D_v9, "parity")
const TRANSIENT_PLOT_DIR_2D_v9 = joinpath(PLOT_ROOT_2D_v9, "transients")
const AXIAL_PLOT_DIR_2D_v9 = joinpath(PLOT_ROOT_2D_v9, "axial_profiles")
const FIELD_PLOT_DIR_2D_v9 = joinpath(PLOT_ROOT_2D_v9, "2d_profiles")
const DIAGNOSTIC_PLOT_DIR_2D_v9 = joinpath(PLOT_ROOT_2D_v9, "diagnostics")

for directory in (
    PARITY_PLOT_DIR_2D_v9,
    TRANSIENT_PLOT_DIR_2D_v9,
    AXIAL_PLOT_DIR_2D_v9,
    FIELD_PLOT_DIR_2D_v9,
    DIAGNOSTIC_PLOT_DIR_2D_v9,
)
    mkpath(directory)
end

const SENSOR_PLOTS_2D_v9 = (
    ("T8", :T8_model_K, :T8_experiment_K, :royalblue),
    ("T12 perimeter", :T12_model_K, :T12_experiment_K, :firebrick),
    ("T11 perimeter", :T11_model_K, :T11_experiment_K, :darkgreen),
    ("T9 core", :T9_model_K, :T9_experiment_K, :darkorange),
    ("T10 core", :T10_model_K, :T10_experiment_K, :saddlebrown),
    ("T3 gas", :T3_model_K, :T3_experiment_K, :purple),
    ("T2 felt", :T2_model_K, :T2_experiment_K, :black),
)

default(
    size=(1400, 900),
    dpi=140,
    linewidth=2.2,
    framestyle=:box,
    gridalpha=0.22,
    legendfontsize=8,
    guidefontsize=11,
    tickfontsize=9,
    titlefontsize=12,
)

function read_fitted_plot_parameters_2D_v9()
    parameter_path = joinpath(
        @__DIR__, "summaries", "2D_v9", "parameters_fitted_2D_v9.csv",
    )
    parameter_rows = CSV.read(parameter_path, DataFrame)
    theta = Float64.(parameter_rows.value)

    summary_lines = readlines(joinpath(
        @__DIR__, "summaries", "2D_v9", "full_test_summary_2D_v9.txt",
    ))
    prefix = "hot_excess_loss_coefficient="
    line = only(filter(row -> startswith(row, prefix), summary_lines))
    hot_loss = parse(Float64, split(line, '=')[2])

    cold = cold_t0_dp1_calibration_2D_v9()
    p0 = with_t0_hydraulics_2D_v9(default_parameters2D(), cold)
    return with_minor_loss_full_2D_v9(
        unpack_parameters2D(theta, p0), hot_loss,
    )
end

function final_heating_rows_2D_v9(predictions)
    rows = DataFrame()
    for id in IDs
        case_rows = filter(
            row -> row.simulation_id == id, predictions,
        )
        push!(rows, case_rows[argmax(case_rows.time_s), :])
    end
    return rows
end

function plot_temperature_parity_2D_v9(predictions)
    final_rows = final_heating_rows_2D_v9(predictions)
    all_values = Float64[]
    for (_, model_column, experiment_column, _) in SENSOR_PLOTS_2D_v9
        append!(all_values, Float64.(final_rows[!, model_column]))
        append!(all_values, Float64.(final_rows[!, experiment_column]))
    end
    lower = floor(minimum(all_values) / 50.0) * 50.0
    upper = ceil(maximum(all_values) / 50.0) * 50.0

    parity = plot(
        xlabel="Model final temperature (K)",
        ylabel="Experimental final temperature (K)",
        title="2D_v9 steady heating parity — circles: training, diamonds: held-out",
        xlims=(lower, upper),
        ylims=(lower, upper),
        aspect_ratio=:equal,
        legend=:outerright,
        size=(1300, 900),
    )
    plot!(
        parity, [lower, upper], [lower, upper];
        color=:gray45, linestyle=:dash, label="1:1",
    )
    marker_shapes = [
        row.phase == "heating_train" ? :circle : :diamond
        for row in eachrow(final_rows)
    ]
    for (sensor, model_column, experiment_column, color) in SENSOR_PLOTS_2D_v9
        scatter!(
            parity,
            final_rows[!, model_column],
            final_rows[!, experiment_column];
            label=sensor,
            color=color,
            markershape=marker_shapes,
            markersize=6,
            markerstrokewidth=0.5,
        )
    end
    path = joinpath(
        PARITY_PLOT_DIR_2D_v9, "steady_temperature_parity_2D_v9.png",
    )
    savefig(parity, path)
    return path
end

function plot_dp1_parity_2D_v9(predictions)
    limits_data = vcat(
        Float64.(predictions.dp1_observed_mbar),
        Float64.(predictions.dp1_base_mbar),
        Float64.(predictions.dp1_augmented_mbar),
    )
    lower = floor(minimum(limits_data) * 2.0) / 2.0
    upper = ceil(maximum(limits_data) * 2.0) / 2.0
    parity = plot(
        xlabel="Model DP1 (mbar)",
        ylabel="Observed DP1 (mbar)",
        title="2D_v9 full-transient DP1 parity",
        xlims=(lower, upper),
        ylims=(lower, upper),
        aspect_ratio=:equal,
        legend=:bottomright,
        size=(1050, 900),
    )
    plot!(
        parity, [lower, upper], [lower, upper];
        color=:gray45, linestyle=:dash, label="1:1",
    )
    scatter!(
        parity,
        predictions.dp1_base_mbar,
        predictions.dp1_observed_mbar;
        label="Cold-linear closure",
        color=:gray55,
        markersize=2.4,
        markerstrokewidth=0,
        alpha=0.40,
    )
    for (phase, label, color, shape) in (
        ("heating_train", "Hot closure: training", :royalblue, :circle),
        ("heating_validation", "Hot closure: held-out heating", :darkorange, :diamond),
        ("cooling_validation", "Hot closure: cooling", :darkgreen, :utriangle),
    )
        rows = filter(row -> row.phase == phase, predictions)
        scatter!(
            parity,
            rows.dp1_augmented_mbar,
            rows.dp1_observed_mbar;
            label=label,
            color=color,
            markershape=shape,
            markersize=3.3,
            markerstrokewidth=0,
            alpha=0.72,
        )
    end
    path = joinpath(PARITY_PLOT_DIR_2D_v9, "dp1_parity_2D_v9.png")
    savefig(parity, path)
    return path
end

function add_temperature_pair_2D_v9!(
    panel, rows, model_column, experiment_column, label, color,
)
    plot!(
        panel,
        rows.time_s,
        rows[!, model_column];
        color=color,
        label="$label model",
    )
    scatter!(
        panel,
        rows.time_s,
        rows[!, experiment_column];
        color=color,
        label="$label experiment",
        markersize=2.4,
        markerstrokewidth=0,
        alpha=0.72,
    )
    return panel
end

function plot_transient_case_2D_v9(rows, simulation_id)
    phase = first(rows.phase)
    p_mid = plot(
        ylabel="Temperature (K)",
        title="Front and mid-depth receiver",
        legend=:outerright,
    )
    for (label, model, experiment, color) in SENSOR_PLOTS_2D_v9[[1, 2, 4]]
        add_temperature_pair_2D_v9!(
            p_mid, rows, model, experiment, label, color,
        )
    end

    p_deep = plot(
        ylabel="Temperature (K)",
        title="Deep receiver",
        legend=:outerright,
    )
    for (label, model, experiment, color) in SENSOR_PLOTS_2D_v9[[3, 5]]
        add_temperature_pair_2D_v9!(
            p_deep, rows, model, experiment, label, color,
        )
    end

    p_aux = plot(
        xlabel="Time (s)",
        ylabel="Temperature (K)",
        title="Gas outlet and insulation",
        legend=:outerright,
    )
    for (label, model, experiment, color) in SENSOR_PLOTS_2D_v9[[6, 7]]
        add_temperature_pair_2D_v9!(
            p_aux, rows, model, experiment, label, color,
        )
    end

    p_dp = plot(
        rows.time_s,
        rows.dp1_observed_mbar;
        xlabel="Time (s)",
        ylabel="DP1 (mbar)",
        title="Flush-tap pressure",
        label="Observed",
        color=:black,
        linewidth=1.6,
        legend=:outerright,
    )
    plot!(
        p_dp,
        rows.time_s,
        rows.dp1_base_mbar;
        label="Cold-linear model",
        color=:gray55,
        linestyle=:dash,
    )
    plot!(
        p_dp,
        rows.time_s,
        rows.dp1_augmented_mbar;
        label="Hot-excess model",
        color=:royalblue,
    )

    combined = plot(
        p_mid, p_deep, p_aux, p_dp;
        layout=(2, 2),
        size=(1800, 1100),
        plot_title="2D_v9 transient: $simulation_id ($phase)",
        margin=5Plots.mm,
    )
    path = joinpath(
        TRANSIENT_PLOT_DIR_2D_v9,
        "transient_$(simulation_id)_2D_v9.png",
    )
    savefig(combined, path)
    return path
end

function plot_all_transients_2D_v9(predictions)
    paths = String[]
    for id in vcat(IDs, collect(VALIDATION_COOLING_2D_v9))
        rows = filter(row -> row.simulation_id == id, predictions)
        push!(paths, plot_transient_case_2D_v9(rows, id))
    end
    return paths
end

function plot_axial_and_field_case_2D_v9(id, p)
    case = run_full_case_2D_v9(id, p; is_cooling=false)
    result = case.result
    experiment = case.experiment
    flow = mean(
        result.mass_flow_kg_s ./
        standard_mass_flow2D(1.0, p.hydraulics)
    )
    irradiance = simulation_conditions[id][Io] / 1000.0

    profile = plot(
        result.z_solid .* 1000.0,
        result.solid_temperature[1, :, end];
        xlabel="Axial position from illuminated face (mm)",
        ylabel="Temperature (K)",
        title="$id final axial profile — $(round(flow, digits=2)) L/min, $(round(irradiance, digits=0)) kW/m²",
        label="Core solid",
        color=:darkorange,
        legend=:outerright,
        size=(1450, 850),
        left_margin=10Plots.mm,
        bottom_margin=8Plots.mm,
    )
    plot!(
        profile,
        result.z_solid .* 1000.0,
        result.solid_temperature[result.nr_rec, :, end];
        label="Perimeter solid",
        color=:firebrick,
        linestyle=:dash,
    )
    plot!(
        profile,
        result.z_solid .* 1000.0,
        result.solid_temperature[end, :, end];
        label="Outer casing",
        color=:black,
        linestyle=:dashdot,
    )
    plot!(
        profile,
        result.z_gas .* 1000.0,
        vec(mean(result.gas_temperature[:, :, end], dims=1));
        label="Mean channel gas",
        color=:darkgreen,
    )
    scatter!(
        profile,
        [5.0, 58.0, 107.0],
        [experiment[end, 1], experiment[end, 2], experiment[end, 3]];
        label="Perimeter measurements",
        color=:firebrick,
        markershape=:diamond,
        markersize=7,
    )
    scatter!(
        profile,
        [58.0, 107.0],
        [experiment[end, 4], experiment[end, 5]];
        label="Core measurements",
        color=:darkorange,
        markershape=:circle,
        markersize=7,
    )
    profile_path = joinpath(
        AXIAL_PLOT_DIR_2D_v9, "axial_profile_$(id)_2D_v9.png",
    )
    savefig(profile, profile_path)

    field = heatmap(
        result.z_solid .* 1000.0,
        result.r_solid .* 1000.0,
        result.solid_temperature[:, :, end];
        xlabel="Axial position (mm)",
        ylabel="Equivalent radius (mm)",
        title="$id final 2D solid temperature (K)",
        color=:inferno,
        size=(1400, 850),
        left_margin=10Plots.mm,
        bottom_margin=8Plots.mm,
    )
    field_path = joinpath(
        FIELD_PLOT_DIR_2D_v9, "temperature_field_$(id)_2D_v9.png",
    )
    savefig(field, field_path)
    return (profile_path=profile_path, field_path=field_path)
end

function plot_identifiability_2D_v9()
    data = CSV.read(joinpath(
        @__DIR__, "summaries", "2D_v9",
        "identifiability_correlation_2D_v9.csv",
    ), DataFrame)
    names = String.(data.parameter)
    matrix = Matrix{Float64}(data[:, 2:end])
    correlation = heatmap(
        names,
        names,
        matrix;
        xlabel="Parameter",
        ylabel="Parameter",
        title="2D_v9 fitted transient-sensitivity correlation",
        color=:balance,
        clims=(-1.0, 1.0),
        xrotation=35,
        size=(1100, 950),
        colorbar_title="Correlation",
    )
    path = joinpath(
        DIAGNOSTIC_PLOT_DIR_2D_v9,
        "identifiability_correlation_2D_v9.png",
    )
    savefig(correlation, path)
    return path
end

function main_plot_2D_v9()
    predictions = CSV.read(joinpath(
        @__DIR__, "summaries", "2D_v9",
        "transient_predictions_2D_v9.csv",
    ), DataFrame)
    p = read_fitted_plot_parameters_2D_v9()

    temperature_parity = plot_temperature_parity_2D_v9(predictions)
    dp1_parity = plot_dp1_parity_2D_v9(predictions)
    transient_paths = plot_all_transients_2D_v9(predictions)
    axial_and_fields = [
        plot_axial_and_field_case_2D_v9(id, p) for id in IDs
    ]
    identifiability = plot_identifiability_2D_v9()

    manifest_path = joinpath(PLOT_ROOT_2D_v9, "plot_manifest_2D_v9.txt")
    open(manifest_path, "w") do io
        println(io, "temperature_parity=$temperature_parity")
        println(io, "dp1_parity=$dp1_parity")
        println(io, "transient_count=$(length(transient_paths))")
        println(io, "axial_profile_count=$(length(axial_and_fields))")
        println(io, "temperature_field_count=$(length(axial_and_fields))")
        println(io, "identifiability=$identifiability")
    end
    println("2D_v9 plots complete: $manifest_path")
    return manifest_path
end

function main_profile_plots_2D_v9()
    p = read_fitted_plot_parameters_2D_v9()
    paths = [plot_axial_and_field_case_2D_v9(id, p) for id in IDs]
    println(
        "2D_v9 axial and field plots complete: $(length(paths)) cases",
    )
    return paths
end

if abspath(PROGRAM_FILE) == @__FILE__
    "--profiles-only" in ARGS ?
        main_profile_plots_2D_v9() :
        main_plot_2D_v9()
end
