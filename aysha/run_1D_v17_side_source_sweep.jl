# ============================================================================
# run_1D_v17_side_source_sweep.jl - diagnostic forced side-source sweep
# ============================================================================
# Holds the fitted v17 parameters fixed except for f_side. This is not a
# calibration. It tests whether plausible side/corner source leakage can flip
# the measured side-core signs T12-T9 and T11-T10 without damaging T2/T3.
# ============================================================================

begin # libraries
    using Statistics
    using StatsPlots
end

include("1D_v17.jl")

begin # configuration
    sweep_int_v17(name, default) = parse(Int, get(ENV, name, string(default)))

    const SWEEP_NODES_v17 = sweep_int_v17(
        "RECEIVER1D_V17_SIDE_SWEEP_NODES", default_nodes,
    )
    const SWEEP_OUTPUT_DIR_v17 = joinpath(
        @__DIR__, "summaries", "1D_v17", "side_source_sweep",
    )
    const SWEEP_PLOT_DIR_v17 = joinpath(SWEEP_OUTPUT_DIR_v17, "plots")
    mkpath(SWEEP_OUTPUT_DIR_v17)
    mkpath(SWEEP_PLOT_DIR_v17)

    const FITTED_SIDE_TOPOLOGY_PARAMS_v17 = [
        0.001502107632562677,
        2.533521266762132,
        PRANDTL_EXPONENT_FIXED_v17,
        1.0,
        184.67243519237965,
        1.0,
        1.669539299381132,
        1.7170694424514947,
        0.9827026344016692,
        0.00031146399009212027,
        32.63508525236098,
        157.4455455398344,
        0.01488263433836151,
    ]
    const SIDE_SOURCE_SWEEP_v17 = (0.0, 0.02, 0.05, 0.10)
end

begin # helpers
    function write_rows_v17(path, rows)
        fields = propertynames(first(rows))
        open(path, "w") do io
            println(io, join(fields, ','))
            for row in rows
                println(io, join((getproperty(row, field) for field in fields), ','))
            end
        end
        return path
    end

    function mean_abs_v17(values)
        return mean(abs.(Float64.(values)))
    end

    function slope_v17(xs, ys)
        xmean = mean(xs)
        ymean = mean(ys)
        denominator = sum((x .- xmean)^2 for x in xs)
        denominator <= eps(Float64) && return NaN
        return sum((xs[i] - xmean) * (ys[i] - ymean) for i in eachindex(xs)) /
               denominator
    end

    function sweep_parameters_v17(f_side)
        params = copy(FITTED_SIDE_TOPOLOGY_PARAMS_v17)
        params[10] = f_side
        return params
    end
end

begin # run sweep
    println("[run_1D_v17_side_source_sweep] Running forced f_side sweep.")
    detailed_rows = NamedTuple[]

    for f_side in SIDE_SOURCE_SWEEP_v17
        params = sweep_parameters_v17(f_side)
        for simulation_id in sim_key_heat
            model, result, experiment = solve_case_v17(
                params, simulation_id; nodes=SWEEP_NODES_v17,
            )
            conditions = simulation_conditions[simulation_id]
            push!(detailed_rows, (
                f_side=f_side,
                simulation_id=simulation_id,
                flow_lpm=conditions[qlpm],
                irradiance=conditions[Io],
                T8_model=model[end, 1],
                T8_experiment=experiment[end, 1],
                T12_side_model=model[end, 2],
                T12_side_experiment=experiment[end, 2],
                T11_side_model=model[end, 3],
                T11_side_experiment=experiment[end, 3],
                T9_core_model=model[end, 4],
                T9_core_experiment=experiment[end, 4],
                T10_core_model=model[end, 5],
                T10_core_experiment=experiment[end, 5],
                T3_model=model[end, 6],
                T3_experiment=experiment[end, 6],
                T2_model=model[end, 7],
                T2_experiment=experiment[end, 7],
                T12_minus_T9_model=model[end, 2] - model[end, 4],
                T12_minus_T9_experiment=experiment[end, 2] - experiment[end, 4],
                T11_minus_T10_model=model[end, 3] - model[end, 5],
                T11_minus_T10_experiment=experiment[end, 3] - experiment[end, 5],
                mean_nu=mean(result.receiver_nusselt[:, end]),
                mean_effectiveness=mean(result.receiver_effectiveness[:, end]),
                receiver_gas_heat=sum(result.receiver_gas_heat[:, end]),
                receiver_to_cavity_heat=result.receiver_to_cavity_heat[end],
                core_side_heat=result.active_wall_heat[end],
                flange_heat_loss=result.flange_heat_loss[end],
            ))
        end
    end

    write_rows_v17(
        joinpath(SWEEP_OUTPUT_DIR_v17, "side_source_sweep_detailed_1D_v17.csv"),
        detailed_rows,
    )

    summary_rows = NamedTuple[]
    for f_side in SIDE_SOURCE_SWEEP_v17
        rows = [row for row in detailed_rows if row.f_side == f_side]
        gap12_errors = [
            row.T12_minus_T9_model - row.T12_minus_T9_experiment for row in rows
        ]
        gap11_errors = [
            row.T11_minus_T10_model - row.T11_minus_T10_experiment for row in rows
        ]
        push!(summary_rows, (
            f_side=f_side,
            mean_T12_minus_T9_model=mean(row.T12_minus_T9_model for row in rows),
            mean_T12_minus_T9_experiment=mean(row.T12_minus_T9_experiment for row in rows),
            T12_minus_T9_positive_model_fraction=mean(row.T12_minus_T9_model > 0.0 for row in rows),
            mean_abs_T12_minus_T9_error=mean_abs_v17(gap12_errors),
            mean_T11_minus_T10_model=mean(row.T11_minus_T10_model for row in rows),
            mean_T11_minus_T10_experiment=mean(row.T11_minus_T10_experiment for row in rows),
            T11_minus_T10_positive_model_fraction=mean(row.T11_minus_T10_model > 0.0 for row in rows),
            mean_abs_T11_minus_T10_error=mean_abs_v17(gap11_errors),
            mean_abs_T2_error=mean_abs_v17(row.T2_model - row.T2_experiment for row in rows),
            mean_abs_T3_error=mean_abs_v17(row.T3_model - row.T3_experiment for row in rows),
            mean_abs_side_error=mean_abs_v17(
                vcat(
                    [row.T8_model - row.T8_experiment for row in rows],
                    [row.T12_side_model - row.T12_side_experiment for row in rows],
                    [row.T11_side_model - row.T11_side_experiment for row in rows],
                ),
            ),
            mean_abs_core_error=mean_abs_v17(
                vcat(
                    [row.T9_core_model - row.T9_core_experiment for row in rows],
                    [row.T10_core_model - row.T10_core_experiment for row in rows],
                ),
            ),
        ))
    end

    write_rows_v17(
        joinpath(SWEEP_OUTPUT_DIR_v17, "side_source_sweep_summary_1D_v17.csv"),
        summary_rows,
    )

    slope_rows = NamedTuple[]
    for f_side in SIDE_SOURCE_SWEEP_v17
        frows = [row for row in detailed_rows if row.f_side == f_side]
        for irradiance in sort(unique(row.irradiance for row in frows); rev=true)
            rows = [row for row in frows if row.irradiance == irradiance]
            flows = [row.flow_lpm for row in rows]
            for sensor in (:T8, :T12_side, :T11_side, :T9_core, :T10_core, :T3)
                model_values = [getproperty(row, Symbol(sensor, "_model")) for row in rows]
                experiment_values = [getproperty(row, Symbol(sensor, "_experiment")) for row in rows]
                push!(slope_rows, (
                    f_side=f_side,
                    irradiance=irradiance,
                    sensor=sensor,
                    model_slope_K_per_Lmin=slope_v17(flows, model_values),
                    experiment_slope_K_per_Lmin=slope_v17(flows, experiment_values),
                ))
            end
        end
    end

    write_rows_v17(
        joinpath(SWEEP_OUTPUT_DIR_v17, "side_source_sweep_flow_slopes_1D_v17.csv"),
        slope_rows,
    )
end

begin # plots
    f_values = [row.f_side for row in summary_rows]
    gap_plot = plot(
        title="1D_v17 forced side-source gap sweep",
        xlabel="Forced side-source fraction",
        ylabel="Mean side-core gap (K)",
        legend=:outerright,
        ylims=(-40.0, 60.0),
        grid=true,
    )
    plot!(
        gap_plot, f_values,
        [row.mean_T12_minus_T9_model for row in summary_rows];
        label="T12-T9 model", marker=:circle, lw=2,
    )
    plot!(
        gap_plot, f_values,
        [row.mean_T12_minus_T9_experiment for row in summary_rows];
        label="T12-T9 experiment", linestyle=:dash, lw=2,
    )
    plot!(
        gap_plot, f_values,
        [row.mean_T11_minus_T10_model for row in summary_rows];
        label="T11-T10 model", marker=:square, lw=2,
    )
    plot!(
        gap_plot, f_values,
        [row.mean_T11_minus_T10_experiment for row in summary_rows];
        label="T11-T10 experiment", linestyle=:dashdot, lw=2,
    )
    savefig(
        gap_plot,
        joinpath(SWEEP_PLOT_DIR_v17, "side_source_gap_sweep_1D_v17.png"),
    )

    error_plot = plot(
        title="1D_v17 forced side-source error sweep",
        xlabel="Forced side-source fraction",
        ylabel="Mean absolute steady error (K)",
        legend=:outerright,
        grid=true,
    )
    for (field, label) in (
        (:mean_abs_T12_minus_T9_error, "T12-T9 gap"),
        (:mean_abs_T11_minus_T10_error, "T11-T10 gap"),
        (:mean_abs_side_error, "side sensors"),
        (:mean_abs_core_error, "core sensors"),
        (:mean_abs_T3_error, "T3"),
        (:mean_abs_T2_error, "T2"),
    )
        plot!(
            error_plot, f_values, [getproperty(row, field) for row in summary_rows];
            label=label, marker=:circle, lw=2,
        )
    end
    savefig(
        error_plot,
        joinpath(SWEEP_PLOT_DIR_v17, "side_source_error_sweep_1D_v17.png"),
    )

    println("[run_1D_v17_side_source_sweep] Complete.")
    println("[run_1D_v17_side_source_sweep] Output directory: $SWEEP_OUTPUT_DIR_v17")
end
