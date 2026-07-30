# ============================================================================
# Complete plotting suite for the final v20 gate-free diagnostic stress point.
#
# Produces:
#   * final-window primary-observable parity;
#   * final-window parity faceted by sensor;
#   * one temporal-history plot for every heating and cooling experiment;
#   * one continuous side-wall axial profile for every heating experiment;
#   * a common-scale axial-profile overview.
#
# All stress parameters are diagnostic and are not validated coefficients.
# ============================================================================

using Statistics
using Plots

include("validate_2D_v20_gatefree.jl")

const GATEFREE_PLOT_EXPONENT_2D_v20 = 0.50
const GATEFREE_PLOT_NU50_RATIO_2D_v20 = 2.00
const GATEFREE_PLOT_FELT_K_2D_v20 = 0.15
const GATEFREE_PLOT_FELT_CP_2D_v20 = 0.55
const GATEFREE_PLOT_POWERS_2D_v20 = (1.05, 1.05, 0.70)
const GATEFREE_PLOT_T3_LOCATION_2D_v20 = :receiver_exit
const GATEFREE_PLOT_T3_CAPACITY_2D_v20 = 30000.0
const GATEFREE_PLOT_T3_STEM_2D_v20 = 0.0
const GATEFREE_PLOT_PRIMARY_INDICES_2D_v20 = (1, 2, 3, 6, 7)
const GATEFREE_PLOT_PRIMARY_NAMES_2D_v20 =
    ("T8", "T12", "T11", "T3", "T2")

function _gatefree_plot_parameters_2D_v20()
    exponent = GATEFREE_PLOT_EXPONENT_2D_v20
    nu50 = GATEFREE_PLOT_NU50_RATIO_2D_v20 *
           GATEFREE_NU50_REFERENCE_2D_v20
    return parameters_2D_v20(
        mesh=:nominal,
        source_model=:near_deep,
        deep_fraction=0.90,
        deep_length_m=0.12,
        nu_prefactor=nu50 / 50.0^exponent,
        reynolds_exponent=exponent,
        distributed_tube_flange_h_W_m2_K=0.0,
        probe_capacity_areal_J_m2_K=
            GATEFREE_PLOT_T3_CAPACITY_2D_v20,
        probe_stem_conductance_areal_W_m2_K=
            GATEFREE_PLOT_T3_STEM_2D_v20,
        felt_conductivity_scale=GATEFREE_PLOT_FELT_K_2D_v20,
        felt_heat_capacity_scale=GATEFREE_PLOT_FELT_CP_2D_v20,
        felt_contact_scale=0.30,
        power_scales=GATEFREE_PLOT_POWERS_2D_v20,
        t3_location=GATEFREE_PLOT_T3_LOCATION_2D_v20,
    )
end

function _gatefree_plot_model_2D_v20(case)
    model = copy(case.model)
    observer = _observer_prediction_2D_v20(
        case,
        GATEFREE_PLOT_T3_LOCATION_2D_v20,
        GATEFREE_PLOT_T3_CAPACITY_2D_v20,
        GATEFREE_PLOT_T3_STEM_2D_v20,
    )
    model[:, 6] .= observer.T3
    return model
end

function _gatefree_axial_profile_2D_v20(case)
    result = case.result
    tau = result.parameters.observation.side_time_constant_s
    final = _gatefree_final_rows_2D_v20(result.time)
    profile = [
        mean(V20.V12._filter_observation(
            result.time,
            vec(result.skin_temperature[index, :]),
            tau,
        )[final])
        for index in axes(result.skin_temperature, 1)
    ]
    observed = [
        mean(case.observed[final, sensor])
        for sensor in 1:3
    ]
    return (
        id=case.inputs.id,
        z_mm=1e3 .* result.z_solid,
        model=profile,
        observed=observed,
    )
end

function _gatefree_common_limits_2D_v20(values; step=25.0)
    lower, upper = extrema(values)
    span = max(upper - lower, step)
    padding = max(0.05span, 0.5step)
    return (
        step * floor((lower - padding) / step),
        step * ceil((upper + padding) / step),
    )
end

function _gatefree_temporal_plot_2D_v20(case, model, path)
    colors = (:red, :darkorange, :purple, :blue, :brown)
    plot_case = plot(
        xlabel="Time (s)",
        ylabel="Temperature (K)",
        title="2D v20 $(case.inputs.id): gate-free diagnostic histories",
        size=(980, 660),
        legend=:outerright,
        left_margin=12Plots.mm,
        top_margin=5Plots.mm,
    )
    for (local_index, sensor_index) in enumerate(
        GATEFREE_PLOT_PRIMARY_INDICES_2D_v20,
    )
        sensor = GATEFREE_PLOT_PRIMARY_NAMES_2D_v20[local_index]
        plot!(
            plot_case,
            case.inputs.times,
            case.observed[:, sensor_index];
            color=colors[local_index],
            linestyle=:dash,
            linewidth=1.5,
            label="$(sensor) measured",
        )
        plot!(
            plot_case,
            case.inputs.times,
            model[:, sensor_index];
            color=colors[local_index],
            linewidth=2.2,
            label="$(sensor) model",
        )
    end
    savefig(plot_case, path)
end

function _gatefree_axial_plot_2D_v20(profile, limits, path)
    axial = plot(
        profile.z_mm,
        profile.model;
        linewidth=2.6,
        color=:royalblue,
        label="continuous model profile",
        xlabel="Axial position (mm)",
        ylabel="Side-wall temperature (K)",
        title="2D v20 $(profile.id): final-window side axial profile",
        xlims=(0.0, 137.0),
        ylims=limits,
        size=(760, 540),
        legend=:topright,
        left_margin=7Plots.mm,
        top_margin=4Plots.mm,
    )
    scatter!(
        axial,
        [11.0, 58.0, 107.0],
        profile.observed;
        marker=:circle,
        markersize=6,
        color=:darkorange,
        markerstrokecolor=:black,
        label="measured thermocouples",
    )
    savefig(axial, path)
end

function _gatefree_axial_overview_2D_v20(
    profiles, limits, path,
)
    panels = Any[]
    for (index, profile) in enumerate(profiles)
        panel = plot(
            profile.z_mm,
            profile.model;
            linewidth=2.0,
            color=:royalblue,
            label=index == 1 ? "continuous model" : "",
            title=profile.id,
            titlefontsize=9,
            xlims=(0.0, 137.0),
            ylims=limits,
            xticks=(0:50:137),
            legend=index == 1 ? :bottomleft : false,
        )
        scatter!(
            panel,
            [11.0, 58.0, 107.0],
            profile.observed;
            marker=:circle,
            markersize=3.5,
            color=:darkorange,
            markerstrokecolor=:black,
            label=index == 1 ? "measured TCs" : "",
        )
        push!(panels, panel)
    end
    overview = plot(
        panels...;
        layout=(3, 5),
        size=(1500, 900),
        xlabel="Axial position (mm)",
        ylabel="Side-wall temperature (K)",
        plot_title=(
            "2D v20 gate-free diagnostic axial profiles — " *
            "common axes and continuous model curves"
        ),
        plot_titlefontsize=15,
        left_margin=5Plots.mm,
        top_margin=3Plots.mm,
        bottom_margin=3Plots.mm,
    )
    savefig(overview, path)
end

function _gatefree_parity_plots_2D_v20(
    final_rows, plot_root,
)
    phase_order = (
        "heating_training", "heating_holdout",
        "cooling_training", "cooling_holdout",
    )
    phase_colors = Dict(
        "heating_training" => :royalblue,
        "heating_holdout" => :darkorange,
        "cooling_training" => :seagreen,
        "cooling_holdout" => :purple,
    )
    phase_markers = Dict(
        "heating_training" => :circle,
        "heating_holdout" => :diamond,
        "cooling_training" => :utriangle,
        "cooling_holdout" => :rect,
    )
    all_values = Float64[]
    for row in final_rows
        append!(all_values, row.observed)
        append!(all_values, row.model)
    end
    limits = _gatefree_common_limits_2D_v20(
        all_values; step=50.0,
    )

    parity = plot(
        xlabel="Measured final-window temperature (K)",
        ylabel="Modeled final-window temperature (K)",
        title="2D v20 gate-free diagnostic primary parity",
        size=(820, 720),
        legend=:topleft,
        xlims=limits,
        ylims=limits,
        aspect_ratio=:equal,
        left_margin=8Plots.mm,
        top_margin=5Plots.mm,
    )
    for phase in phase_order
        selected = filter(row -> row.phase == phase, final_rows)
        xs = reduce(vcat, (row.observed for row in selected))
        ys = reduce(vcat, (row.model for row in selected))
        scatter!(
            parity,
            xs,
            ys;
            color=phase_colors[phase],
            marker=phase_markers[phase],
            alpha=0.82,
            markersize=6,
            label=replace(phase, "_" => " "),
        )
    end
    plot!(
        parity,
        collect(limits),
        collect(limits);
        color=:black,
        linestyle=:dash,
        linewidth=1.5,
        label="1:1",
    )
    savefig(
        parity,
        joinpath(
            plot_root,
            "parity_primary_final_2D_v20_gatefree.png",
        ),
    )

    sensor_panels = Any[]
    for (sensor_position, sensor) in enumerate(
        GATEFREE_PLOT_PRIMARY_NAMES_2D_v20,
    )
        panel = plot(
            title=sensor,
            titlefontsize=10,
            xlims=limits,
            ylims=limits,
            aspect_ratio=:equal,
            legend=sensor_position == 1 ? :topleft : false,
        )
        for phase in phase_order
            selected = filter(
                row -> row.phase == phase,
                final_rows,
            )
            scatter!(
                panel,
                [row.observed[sensor_position] for row in selected],
                [row.model[sensor_position] for row in selected];
                color=phase_colors[phase],
                marker=phase_markers[phase],
                markersize=5,
                alpha=0.85,
                label=sensor_position == 1 ?
                    replace(phase, "_" => " ") : "",
            )
        end
        plot!(
            panel,
            collect(limits),
            collect(limits);
            color=:black,
            linestyle=:dash,
            linewidth=1.2,
            label=sensor_position == 1 ? "1:1" : "",
        )
        push!(sensor_panels, panel)
    end
    faceted = plot(
        sensor_panels...;
        layout=(2, 3),
        size=(1350, 850),
        xlabel="Measured final-window temperature (K)",
        ylabel="Modeled final-window temperature (K)",
        plot_title="2D v20 gate-free diagnostic parity by sensor",
        plot_titlefontsize=15,
        left_margin=5Plots.mm,
        top_margin=4Plots.mm,
        bottom_margin=4Plots.mm,
    )
    savefig(
        faceted,
        joinpath(
            plot_root,
            "parity_by_sensor_final_2D_v20_gatefree.png",
        ),
    )
end

function plot_2D_v20_gatefree_full(; max_points=121)
    plot_root = joinpath(
        OUTPUT_DIR_2D_v20, "plots", "gatefree_full",
    )
    transient_dir = joinpath(plot_root, "transients")
    axial_dir = joinpath(plot_root, "axial_profiles")
    foreach(mkpath, (plot_root, transient_dir, axial_dir))

    p = _gatefree_plot_parameters_2D_v20()
    specs = vcat(
        [(id=id, cooling=false) for id in HEATING_IDS_2D_v20],
        [(id=id, cooling=true)
         for id in GATEFREE_COOLING_IDS_2D_v20],
    )
    cases = Vector{Any}(undef, length(specs))
    Threads.@threads for index in eachindex(specs)
        spec = specs[index]
        inputs = case_inputs_2D_v20(
            spec.id;
            cooling=spec.cooling,
            max_points=max_points,
        )
        cases[index] = simulate_case_2D_v20(
            inputs,
            p;
            initialization=spec.cooling ?
                :side_T2_only : :ambient,
            reltol=5e-4,
            abstol=1e-4,
            dtmax=120.0,
        )
        println(
            "v20 full plotting simulation ",
            index, "/", length(specs), ": ", spec.id,
        )
        flush(stdout)
    end

    models = [_gatefree_plot_model_2D_v20(case) for case in cases]
    final_rows = NamedTuple[]
    axial_profiles = NamedTuple[]
    for (case, model) in zip(cases, models)
        final = _gatefree_final_rows_2D_v20(case.inputs.times)
        push!(final_rows, (
            simulation_id=case.inputs.id,
            phase=_gatefree_phase_2D_v20(
                case.inputs.id, case.inputs.cooling,
            ),
            observed=[
                mean(case.observed[final, sensor])
                for sensor in GATEFREE_PLOT_PRIMARY_INDICES_2D_v20
            ],
            model=[
                mean(model[final, sensor])
                for sensor in GATEFREE_PLOT_PRIMARY_INDICES_2D_v20
            ],
        ))
        _gatefree_temporal_plot_2D_v20(
            case,
            model,
            joinpath(
                transient_dir,
                "$(case.inputs.id)_transient.png",
            ),
        )
        if !case.inputs.cooling
            push!(
                axial_profiles,
                _gatefree_axial_profile_2D_v20(case),
            )
        end
    end

    axial_values = Float64[]
    for profile in axial_profiles
        append!(axial_values, profile.model)
        append!(axial_values, profile.observed)
    end
    axial_limits = _gatefree_common_limits_2D_v20(
        axial_values; step=25.0,
    )
    for profile in axial_profiles
        _gatefree_axial_plot_2D_v20(
            profile,
            axial_limits,
            joinpath(
                axial_dir,
                "$(profile.id)_axial.png",
            ),
        )
    end
    _gatefree_axial_overview_2D_v20(
        axial_profiles,
        axial_limits,
        joinpath(
            plot_root,
            "axial_profiles_common_scale_2D_v20_gatefree.png",
        ),
    )
    _gatefree_parity_plots_2D_v20(final_rows, plot_root)

    open(joinpath(plot_root, "plot_manifest_2D_v20_gatefree.txt"), "w") do io
        println(io, "V20 GATE-FREE FULL PLOT MANIFEST")
        println(io, "Diagnostic stress parameters; not validated coefficients.")
        println(io, "heating_cases=", length(axial_profiles))
        println(io, "all_cases=", length(cases))
        println(io, "transient_plots=", length(cases))
        println(io, "axial_plots=", length(axial_profiles))
        println(io, "axial_x_limits_mm=0,137")
        println(io, "axial_y_limits_K=", first(axial_limits), ",",
                last(axial_limits))
        println(io, "axial_profiles_use_final_actual_time_10_percent=true")
        println(io, "axial_model_curves_use_all_axial_mesh_centers=true")
        println(io, "parity_uses_final_actual_time_10_percent=true")
        println(io, "reynolds_exponent=",
                GATEFREE_PLOT_EXPONENT_2D_v20)
        println(io, "Nu50_ratio=",
                GATEFREE_PLOT_NU50_RATIO_2D_v20)
        println(io, "felt_k_scale=",
                GATEFREE_PLOT_FELT_K_2D_v20)
        println(io, "felt_Cp_scale=",
                GATEFREE_PLOT_FELT_CP_2D_v20)
        println(io, "power_scales=",
                join(GATEFREE_PLOT_POWERS_2D_v20, ","))
        println(io, "T3_location=",
                GATEFREE_PLOT_T3_LOCATION_2D_v20)
        println(io, "T3_capacity_areal_J_m2_K=",
                GATEFREE_PLOT_T3_CAPACITY_2D_v20)
        println(io, "T3_stem_conductance_areal_W_m2_K=",
                GATEFREE_PLOT_T3_STEM_2D_v20)
    end
    return (
        cases=cases,
        models=models,
        finals=final_rows,
        axial_profiles=axial_profiles,
        axial_limits=axial_limits,
        plot_root=plot_root,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    points = isempty(ARGS) ? 121 : parse(Int, ARGS[1])
    plot_2D_v20_gatefree_full(; max_points=points)
end
