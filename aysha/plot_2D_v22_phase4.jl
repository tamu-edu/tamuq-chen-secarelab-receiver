# ============================================================================
# plot_2D_v22_phase4.jl
#
# Generates validation and parity plots using the optimal Phase 4 parameters
# for ALL cases (E67 to E81)
# ============================================================================

using Statistics
using Plots
using CSV
using DataFrames

include("calibrate_2D_v22_phase4.jl")

const V22_PRIMARY_INDICES = (1, 2, 3, 6, 7)
const V22_PRIMARY_NAMES = ("T8", "T12", "T11", "T3", "T2")

function _phase_2D_v22(id, cooling)
    if cooling
        return id == "C81" ?
            "cooling_holdout" : "cooling_training"
    end
    return "heating_training"
end

function _final_rows_phase4(times)
    threshold = first(times) + 0.90 * (last(times) - first(times))
    return something(findfirst(t -> t >= threshold, times), length(times)):length(times)
end

function _continuous_axial_profile_2D_v22(case)
    result = case.result
    final = _final_rows_phase4(result.time)
    tau = result.parameters.base.base.base.base.base.base.observation.side_time_constant_s
    model = [
        mean(v22.V20.V12._filter_observation(
            result.time,
            vec(result.skin_temperature[index, :]),
            tau,
        )[final])
        for index in axes(result.skin_temperature, 1)
    ]
    return (
        id=case.inputs.id,
        z_mm=1e3 .* result.z_solid,
        model=model,
        observed=[
            mean(case.observed[final, sensor]) for sensor in 1:3
        ],
    )
end

function _common_limits_2D_v22(values; step=25.0)
    lower, upper = extrema(values)
    padding = max(0.05 * (upper - lower), 0.5 * step)
    return (
        step * floor((lower - padding) / step),
        step * ceil((upper + padding) / step),
    )
end

function _temporal_plot_2D_v22(case, path)
    colors = (:red, :darkorange, :purple, :blue, :brown)
    figure = plot(
        xlabel="Time (s)", ylabel="Temperature (K)",
        title="Phase 4 $(case.inputs.id) Transient",
        size=(980, 660), legend=:outerright,
        left_margin=12Plots.mm, top_margin=5Plots.mm,
    )
    for (local_index, sensor_index) in enumerate(V22_PRIMARY_INDICES)
        sensor = V22_PRIMARY_NAMES[local_index]
        plot!(
            figure, case.inputs.times,
            case.observed[:, sensor_index];
            color=colors[local_index], linestyle=:dash,
            linewidth=1.5, label="$sensor measured",
        )
        plot!(
            figure, case.inputs.times,
            case.model[:, sensor_index];
            color=colors[local_index], linewidth=2.2,
            label="$sensor model",
        )
    end
    savefig(figure, path)
end

function _axial_plot_2D_v22(profile, limits, path)
    figure = plot(
        profile.z_mm, profile.model;
        linewidth=2.6, color=:royalblue,
        label="continuous model profile",
        xlabel="Axial position (mm)",
        ylabel="Side-wall temperature (K)",
        title="Phase 4 $(profile.id): final-window axial profile",
        xlims=(0.0, 137.0), ylims=limits,
        size=(760, 540), legend=:topright,
        left_margin=7Plots.mm, top_margin=4Plots.mm,
    )
    scatter!(
        figure, [11.0, 58.0, 107.0], profile.observed;
        marker=:circle, markersize=6, color=:darkorange,
        markerstrokecolor=:black,
        label="measured thermocouples",
    )
    savefig(figure, path)
end

function _axial_overview_2D_v22(profiles, limits, path)
    panels = Any[]
    for (index, profile) in enumerate(profiles)
        panel = plot(
            profile.z_mm, profile.model;
            linewidth=2.0, color=:royalblue,
            label=index == 1 ? "continuous model" : "",
            title=profile.id, titlefontsize=9,
            xlims=(0.0, 137.0), ylims=limits,
            xticks=(0:50:137),
            legend=index == 1 ? :bottomleft : false,
        )
        scatter!(
            panel, [11.0, 58.0, 107.0], profile.observed;
            marker=:circle, markersize=3.5,
            color=:darkorange, markerstrokecolor=:black,
            label=index == 1 ? "measured TCs" : "",
        )
        push!(panels, panel)
    end
    overview = plot(
        panels...; layout=(3, ceil(Int, length(profiles)/3)), size=(1500, 900),
        xlabel="Axial position (mm)",
        ylabel="Side-wall temperature (K)",
        plot_title=(
            "Phase 4 Optimal Candidate — common axes and continuous model curves"
        ),
        plot_titlefontsize=15, left_margin=5Plots.mm,
        top_margin=3Plots.mm, bottom_margin=3Plots.mm,
    )
    savefig(overview, path)
end

function _parity_plots_2D_v22(final_rows, plot_root)
    phase_order = ("heating_training", "cooling_training", "cooling_holdout")
    colors = Dict(
        "heating_training" => :royalblue,
        "cooling_training" => :seagreen,
        "cooling_holdout" => :purple,
    )
    markers = Dict(
        "heating_training" => :circle,
        "cooling_training" => :utriangle,
        "cooling_holdout" => :rect,
    )
    values = reduce(vcat, (vcat(row.observed, row.model) for row in final_rows))
    limits = _common_limits_2D_v22(values; step=50.0)
    parity = plot(
        xlabel="Measured final-window temperature (K)",
        ylabel="Modeled final-window temperature (K)",
        title="Phase 4 Optimal Candidate Primary Parity",
        size=(820, 720), legend=:topleft,
        xlims=limits, ylims=limits, aspect_ratio=:equal,
        left_margin=8Plots.mm, top_margin=5Plots.mm,
    )
    for phase in phase_order
        selected = filter(row -> row.phase == phase, final_rows)
        if isempty(selected) continue end
        scatter!(
            parity,
            reduce(vcat, (row.observed for row in selected)),
            reduce(vcat, (row.model for row in selected));
            color=colors[phase], marker=markers[phase],
            alpha=0.82, markersize=6,
            label=replace(phase, "_" => " "),
        )
    end
    plot!(
        parity, collect(limits), collect(limits);
        color=:black, linestyle=:dash, linewidth=1.5,
        label="1:1",
    )
    savefig(parity, joinpath(plot_root, "parity_primary_final_phase4.png"))

    panels = Any[]
    for (sensor_index, sensor) in enumerate(V22_PRIMARY_NAMES)
        panel = plot(
            title=sensor, titlefontsize=10,
            xlims=limits, ylims=limits, aspect_ratio=:equal,
            legend=sensor_index == 1 ? :topleft : false,
        )
        for phase in phase_order
            selected = filter(row -> row.phase == phase, final_rows)
            if isempty(selected) continue end
            scatter!(
                panel,
                [row.observed[sensor_index] for row in selected],
                [row.model[sensor_index] for row in selected];
                color=colors[phase], marker=markers[phase],
                markersize=5, alpha=0.85,
                label=sensor_index == 1 ? replace(phase, "_" => " ") : "",
            )
        end
        plot!(
            panel, collect(limits), collect(limits);
            color=:black, linestyle=:dash, linewidth=1.2,
            label=sensor_index == 1 ? "1:1" : "",
        )
        push!(panels, panel)
    end
    figure = plot(
        panels...; layout=(2, 3), size=(1350, 850),
        xlabel="Measured final-window temperature (K)",
        ylabel="Modeled final-window temperature (K)",
        plot_title="Phase 4 Parity by sensor",
        plot_titlefontsize=15, left_margin=5Plots.mm,
        top_margin=4Plots.mm, bottom_margin=4Plots.mm,
    )
    savefig(figure, joinpath(plot_root, "parity_by_sensor_final_phase4.png"))
end

function plot_phase4()
    # Load optimal parameters from Phase 4 summary
    outdir = joinpath(@__DIR__, "summaries", "2D_v22")
    grid_path = joinpath(outdir, "phase4_calibration_grid.csv")
    if !isfile(grid_path)
        error("Calibration grid not found! Run calibrate_2D_v22_phase4.jl first.")
    end
    df = CSV.read(grid_path, DataFrame)
    best_row = df[argmin(df.objective), :]
    
    best_nu = best_row.nu_multiplier
    best_ross = best_row.ross_multiplier
    println("Plotting with Best Phase 4 Parameters: Nu_Mult=$best_nu, Ross_Mult=$best_ross")
    
    # build_base_parameters sets GLOBAL_NU_MULTIPLIER and GLOBAL_ROSS_MULTIPLIER
    # We call build_base_parameters from calibrate_2D_v22_phase4.jl
    p_mod = build_base_parameters(best_nu, best_ross)

    plot_root = joinpath(outdir, "plots")
    transient_dir = joinpath(plot_root, "transients")
    axial_dir = joinpath(plot_root, "axial_profiles")
    foreach(mkpath, (plot_root, transient_dir, axial_dir))

    ids = ["E$(66+i)" for i in 1:15] # E67 to E81
    
    cases = []
    
    for id in ids
        println("Simulating $id (heating)...")
        inputs_heat = Main.case_inputs_2D_v22(id; cooling=false)
        push!(cases, simulate_case_phase4(inputs_heat, p_mod; full_initial_data=false))
        
        c_id = replace(id, "E" => "C")
        println("Simulating $c_id (cooling)...")
        try
            inputs_cool = Main.case_inputs_2D_v22(c_id; cooling=true)
            push!(cases, simulate_case_phase4(inputs_cool, p_mod; full_initial_data=false))
        catch e
            println("Failed to simulate $c_id: $e")
        end
    end
    
    final_rows = NamedTuple[]
    axial_profiles = NamedTuple[]
    for case in cases
        final = _final_rows_phase4(case.inputs.times)
        push!(final_rows, (
            phase=_phase_2D_v22(case.inputs.id, case.inputs.cooling),
            observed=[mean(case.observed[final, sensor]) for sensor in (1, 2, 3, 6, 7)],
            model=[mean(case.model[final, sensor]) for sensor in (1, 2, 3, 6, 7)],
        ))
        _temporal_plot_2D_v22(
            case,
            joinpath(transient_dir, "$(case.inputs.id)_$(case.inputs.cooling ? "cool" : "heat")_transient.png"),
        )
        if !case.inputs.cooling
            push!(axial_profiles, _continuous_axial_profile_2D_v22(case))
        end
    end
    
    axial_values = reduce(vcat, (vcat(profile.model, profile.observed) for profile in axial_profiles))
    axial_limits = _common_limits_2D_v22(axial_values; step=25.0)
    
    for profile in axial_profiles
        _axial_plot_2D_v22(profile, axial_limits, joinpath(axial_dir, "$(profile.id)_axial.png"))
    end
    _axial_overview_2D_v22(axial_profiles, axial_limits, joinpath(plot_root, "axial_profiles_common_scale_phase4.png"))
    _parity_plots_2D_v22(final_rows, plot_root)
    println("Plotting complete! Saved to $plot_root")
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_phase4()
end
