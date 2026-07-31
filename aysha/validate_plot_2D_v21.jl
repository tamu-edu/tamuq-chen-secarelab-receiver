# ============================================================================
# validate_plot_2D_v21.jl
#
# Full v21 validation of the pre-declared hot-front hypothesis.  The selected
# upper-bound candidate is plotted; the 20 cm2 interior side-optimal candidate
# is retained as a full-heating comparator.
# ============================================================================

using Statistics
using Plots

include("run_2D_v21.jl")

const V21_BOUNDARY = (
    name="boundary40",
    area_cm2=40.0,
    powers=(1.65, 1.80, 1.11),
)
const V21_INTERIOR = (
    name="interior20",
    area_cm2=20.0,
    powers=(1.50, 1.58, 0.94),
)
const V21_PRIMARY_INDICES = (1, 2, 3, 6, 7)
const V21_PRIMARY_NAMES = ("T8", "T12", "T11", "T3", "T2")
const V21_WALL_WEIGHTS = (0.251825, 0.350365, 0.397810)

function _final_rows_2D_v21(times)
    length(times) >= 2 || error(
        "v21 validation requires at least two samples",
    )
    threshold = first(times) +
        0.90 * (last(times) - first(times))
    start = something(
        findfirst(time -> time >= threshold, times),
        length(times),
    )
    return start:length(times)
end

function _phase_2D_v21(id, cooling)
    if cooling
        return id == "C81" ?
            "cooling_holdout" : "cooling_training"
    end
    return id in TRAIN_HEATING_2D_v21 ?
        "heating_training" : "heating_holdout"
end

function _continuous_axial_profile_2D_v21(case)
    result = case.result
    final = _final_rows_2D_v21(result.time)
    tau = result.parameters.observation.side_time_constant_s
    model = [
        mean(V21.V12._filter_observation(
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

function _case_row_2D_v21(case, candidate)
    final = _final_rows_2D_v21(case.inputs.times)
    observed = [
        mean(case.observed[final, sensor])
        for sensor in V21_PRIMARY_INDICES
    ]
    model = [
        mean(case.model[final, sensor])
        for sensor in V21_PRIMARY_INDICES
    ]
    side_observed = observed[1:3]
    side_model = model[1:3]
    ambient = mean(case.inputs.ambient[final])
    inlet = mean(case.inputs.inlet[final])
    observed_wall = sum(
        V21_WALL_WEIGHTS[i] *
        (side_observed[i] - ambient) for i in 1:3
    )
    model_wall = sum(
        V21_WALL_WEIGHTS[i] *
        (side_model[i] - ambient) for i in 1:3
    )
    observed_effectiveness =
        (observed[4] - inlet) / max(observed_wall, 1.0)
    model_effectiveness =
        (model[4] - inlet) / max(model_wall, 1.0)
    diagnostic = V21.hot_front_radiative_diagnostic2D(
        case.result.ode_solution.u[end],
        case.result.parameters,
        case.operating_condition,
        case.result.time[end],
    )
    ledger = V21.energy_rate_ledger2D(
        case.result.ode_solution.u[end],
        case.result.parameters,
        case.operating_condition,
        case.result.time[end],
    )
    scale = case.inputs.nominal >= 400000.0 ?
        case.result.parameters.optics.scale_456 :
        (case.inputs.nominal >= 280000.0 ?
         case.result.parameters.optics.scale_304 :
         case.result.parameters.optics.scale_256)
    absorbed = case.inputs.cooling ? 0.0 :
        case.result.parameters.optics.absorbed_fraction *
        scale * case.inputs.nominal *
        case.result.parameters.geometry.receiver_width^2
    side_residual = side_model .- side_observed
    return (
        candidate=candidate.name,
        simulation_id=case.inputs.id,
        phase=_phase_2D_v21(
            case.inputs.id, case.inputs.cooling,
        ),
        nominal_kW_m2=case.inputs.nominal / 1000.0,
        flow_L_min=mean(case.inputs.flow[final]),
        side_final_bias_K=mean(side_residual),
        side_final_rmse_K=sqrt(mean(abs2, side_residual)),
        T8_bias_K=side_residual[1],
        T12_bias_K=side_residual[2],
        T11_bias_K=side_residual[3],
        T3_bias_K=model[4] - observed[4],
        T2_bias_K=model[5] - observed[5],
        axial_rmse_K=sqrt(mean(abs2, side_residual)),
        measured_middle_peak=(
            side_observed[2] > side_observed[1] &&
            side_observed[2] > side_observed[3]
        ),
        modeled_middle_peak=(
            side_model[2] > side_model[1] &&
            side_model[2] > side_model[3]
        ),
        measured_effectiveness=observed_effectiveness,
        model_effectiveness=model_effectiveness,
        hot_front_extra_loss_W=
            diagnostic.extra_radiative_loss_W,
        hot_front_extra_share=absorbed > 0.0 ?
            diagnostic.extra_radiative_loss_W / absorbed : 0.0,
        receiver_front_loss_W=ledger.receiver_front_loss,
        outer_front_loss_W=ledger.outer_front_loss,
        total_front_boundary_loss_W=
            ledger.receiver_front_loss + ledger.outer_front_loss,
        total_front_boundary_share=absorbed > 0.0 ?
            (ledger.receiver_front_loss + ledger.outer_front_loss) /
            absorbed : 0.0,
        energy_residual_W=ledger.residual,
        ode_success=V21.SciMLBase.successful_retcode(
            case.result.ode_solution.retcode,
        ),
    )
end

function _aggregate_2D_v21(rows, phase)
    selected = filter(row -> row.phase == phase, rows)
    return (
        candidate=first(selected).candidate,
        phase=phase,
        cases=length(selected),
        side_final_rmse_K=sqrt(mean(
            row.side_final_rmse_K^2 for row in selected
        )),
        side_final_bias_K=mean(
            row.side_final_bias_K for row in selected
        ),
        T3_final_rmse_K=sqrt(mean(
            row.T3_bias_K^2 for row in selected
        )),
        T3_final_bias_K=mean(
            row.T3_bias_K for row in selected
        ),
        T2_final_rmse_K=sqrt(mean(
            row.T2_bias_K^2 for row in selected
        )),
        T2_final_bias_K=mean(
            row.T2_bias_K for row in selected
        ),
        axial_rmse_K=sqrt(mean(
            row.axial_rmse_K^2 for row in selected
        )),
        modeled_middle_peaks=count(
            row.modeled_middle_peak for row in selected
        ),
        measured_middle_peaks=count(
            row.measured_middle_peak for row in selected
        ),
        effectiveness_rmse=sqrt(mean(
            (row.model_effectiveness -
             row.measured_effectiveness)^2
            for row in selected
        )),
        model_effectiveness_min=minimum(
            row.model_effectiveness for row in selected
        ),
        model_effectiveness_max=maximum(
            row.model_effectiveness for row in selected
        ),
        max_energy_residual_W=maximum(
            abs(row.energy_residual_W) for row in selected
        ),
    )
end

function _common_limits_2D_v21(values; step=25.0)
    lower, upper = extrema(values)
    padding = max(0.05 * (upper - lower), 0.5 * step)
    return (
        step * floor((lower - padding) / step),
        step * ceil((upper + padding) / step),
    )
end

function _temporal_plot_2D_v21(case, path)
    colors = (:red, :darkorange, :purple, :blue, :brown)
    figure = plot(
        xlabel="Time (s)", ylabel="Temperature (K)",
        title="2D v21 $(case.inputs.id): hot-front boundary candidate",
        size=(980, 660), legend=:outerright,
        left_margin=12Plots.mm, top_margin=5Plots.mm,
    )
    for (local_index, sensor_index) in enumerate(
        V21_PRIMARY_INDICES,
    )
        sensor = V21_PRIMARY_NAMES[local_index]
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

function _axial_plot_2D_v21(profile, limits, path)
    figure = plot(
        profile.z_mm, profile.model;
        linewidth=2.6, color=:royalblue,
        label="continuous model profile",
        xlabel="Axial position (mm)",
        ylabel="Side-wall temperature (K)",
        title="2D v21 $(profile.id): final-window side axial profile",
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

function _axial_overview_2D_v21(profiles, limits, path)
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
        panels...; layout=(3, 5), size=(1500, 900),
        xlabel="Axial position (mm)",
        ylabel="Side-wall temperature (K)",
        plot_title=(
            "2D v21 hot-front boundary candidate — " *
            "common axes and continuous model curves"
        ),
        plot_titlefontsize=15, left_margin=5Plots.mm,
        top_margin=3Plots.mm, bottom_margin=3Plots.mm,
    )
    savefig(overview, path)
end

function _parity_plots_2D_v21(final_rows, plot_root)
    phase_order = (
        "heating_training", "heating_holdout",
        "cooling_training", "cooling_holdout",
    )
    colors = Dict(
        "heating_training" => :royalblue,
        "heating_holdout" => :darkorange,
        "cooling_training" => :seagreen,
        "cooling_holdout" => :purple,
    )
    markers = Dict(
        "heating_training" => :circle,
        "heating_holdout" => :diamond,
        "cooling_training" => :utriangle,
        "cooling_holdout" => :rect,
    )
    values = reduce(vcat, (
        vcat(row.observed, row.model) for row in final_rows
    ))
    limits = _common_limits_2D_v21(values; step=50.0)
    parity = plot(
        xlabel="Measured final-window temperature (K)",
        ylabel="Modeled final-window temperature (K)",
        title="2D v21 hot-front boundary candidate primary parity",
        size=(820, 720), legend=:topleft,
        xlims=limits, ylims=limits, aspect_ratio=:equal,
        left_margin=8Plots.mm, top_margin=5Plots.mm,
    )
    for phase in phase_order
        selected = filter(row -> row.phase == phase, final_rows)
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
    savefig(
        parity,
        joinpath(plot_root, "parity_primary_final_2D_v21.png"),
    )

    panels = Any[]
    for (sensor_index, sensor) in enumerate(V21_PRIMARY_NAMES)
        panel = plot(
            title=sensor, titlefontsize=10,
            xlims=limits, ylims=limits, aspect_ratio=:equal,
            legend=sensor_index == 1 ? :topleft : false,
        )
        for phase in phase_order
            selected = filter(
                row -> row.phase == phase, final_rows,
            )
            scatter!(
                panel,
                [row.observed[sensor_index] for row in selected],
                [row.model[sensor_index] for row in selected];
                color=colors[phase], marker=markers[phase],
                markersize=5, alpha=0.85,
                label=sensor_index == 1 ?
                    replace(phase, "_" => " ") : "",
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
        plot_title="2D v21 hot-front boundary parity by sensor",
        plot_titlefontsize=15, left_margin=5Plots.mm,
        top_margin=4Plots.mm, bottom_margin=4Plots.mm,
    )
    savefig(
        figure,
        joinpath(plot_root, "parity_by_sensor_final_2D_v21.png"),
    )
end

function _simulate_candidate_2D_v21(
    candidate, specs; max_points=121,
)
    p = parameters_2D_v21(
        effective_front_area_cm2=candidate.area_cm2,
        power_scales=candidate.powers,
    )
    cases = Vector{Any}(undef, length(specs))
    Threads.@threads for index in eachindex(specs)
        spec = specs[index]
        inputs = case_inputs_2D_v21(
            spec.id; cooling=spec.cooling,
            max_points=max_points,
        )
        cases[index] = simulate_case_2D_v21(
            inputs, p; full_initial_data=spec.cooling,
            reltol=5e-4, abstol=1e-4, dtmax=120.0,
        )
        println(
            "v21 ", candidate.name, " ", index, "/",
            length(specs), ": ", spec.id,
        )
        flush(stdout)
    end
    return cases
end

function validate_plot_2D_v21(; max_points=121)
    plot_root = joinpath(
        OUTPUT_DIR_2D_v21, "plots", "boundary40_full",
    )
    transient_dir = joinpath(plot_root, "transients")
    axial_dir = joinpath(plot_root, "axial_profiles")
    foreach(mkpath, (plot_root, transient_dir, axial_dir))

    full_specs = vcat(
        [(id=id, cooling=false) for id in HEATING_IDS_2D_v21],
        [(id=id, cooling=true) for id in COOLING_IDS_2D_v21],
    )
    heating_specs = [
        (id=id, cooling=false) for id in HEATING_IDS_2D_v21
    ]
    boundary_cases = _simulate_candidate_2D_v21(
        V21_BOUNDARY, full_specs; max_points=max_points,
    )
    interior_cases = _simulate_candidate_2D_v21(
        V21_INTERIOR, heating_specs; max_points=max_points,
    )
    boundary_rows = [
        _case_row_2D_v21(case, V21_BOUNDARY)
        for case in boundary_cases
    ]
    interior_rows = [
        _case_row_2D_v21(case, V21_INTERIOR)
        for case in interior_cases
    ]
    all_rows = vcat(boundary_rows, interior_rows)
    _write_namedtuples_csv_2D_v21(
        joinpath(OUTPUT_DIR_2D_v21, "validation_cases_2D_v21.csv"),
        all_rows,
    )
    phase_order = (
        "heating_training", "heating_holdout",
        "cooling_training", "cooling_holdout",
    )
    phase_rows = NamedTuple[]
    for phase in phase_order
        push!(
            phase_rows,
            _aggregate_2D_v21(boundary_rows, phase),
        )
    end
    for phase in ("heating_training", "heating_holdout")
        push!(
            phase_rows,
            _aggregate_2D_v21(interior_rows, phase),
        )
    end
    _write_namedtuples_csv_2D_v21(
        joinpath(OUTPUT_DIR_2D_v21, "validation_phases_2D_v21.csv"),
        phase_rows,
    )

    final_rows = NamedTuple[]
    axial_profiles = NamedTuple[]
    for case in boundary_cases
        final = _final_rows_2D_v21(case.inputs.times)
        push!(final_rows, (
            phase=_phase_2D_v21(
                case.inputs.id, case.inputs.cooling,
            ),
            observed=[
                mean(case.observed[final, sensor])
                for sensor in V21_PRIMARY_INDICES
            ],
            model=[
                mean(case.model[final, sensor])
                for sensor in V21_PRIMARY_INDICES
            ],
        ))
        _temporal_plot_2D_v21(
            case,
            joinpath(
                transient_dir,
                "$(case.inputs.id)_transient.png",
            ),
        )
        if !case.inputs.cooling
            push!(
                axial_profiles,
                _continuous_axial_profile_2D_v21(case),
            )
        end
    end
    axial_values = reduce(vcat, (
        vcat(profile.model, profile.observed)
        for profile in axial_profiles
    ))
    axial_limits = _common_limits_2D_v21(
        axial_values; step=25.0,
    )
    for profile in axial_profiles
        _axial_plot_2D_v21(
            profile, axial_limits,
            joinpath(
                axial_dir, "$(profile.id)_axial.png",
            ),
        )
    end
    _axial_overview_2D_v21(
        axial_profiles, axial_limits,
        joinpath(
            plot_root,
            "axial_profiles_common_scale_2D_v21.png",
        ),
    )
    _parity_plots_2D_v21(final_rows, plot_root)
    open(
        joinpath(plot_root, "plot_manifest_2D_v21.txt"), "w",
    ) do io
        println(io, "2D V21 FULL PLOT MANIFEST")
        println(io, "candidate=boundary40")
        println(io, "effective_front_area_cm2=40.0")
        println(io, "power_scales=1.65,1.80,1.11")
        println(io, "all_cases=18")
        println(io, "heating_axial_profiles=15")
        println(io, "continuous_axial_model_curves=true")
        println(io, "common_axial_y_limits_K=",
                first(axial_limits), ",", last(axial_limits))
        println(io, "final_window=last_10_percent_actual_time")
    end
    println("wrote v21 validation and plots to $OUTPUT_DIR_2D_v21")
    return (
        boundary_rows=boundary_rows,
        interior_rows=interior_rows,
        phase_rows=phase_rows,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    points = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 121
    validate_plot_2D_v21(; max_points=points)
end
