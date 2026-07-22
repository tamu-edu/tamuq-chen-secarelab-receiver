# ============================================================================
# run_1D_v19.jl - 2-Zone Core/Perimeter Macro-ECM Runner
# ============================================================================
# Runner for v19: 2-Zone Core/Perimeter Macro-ECM model.
#
# Environment variables:
#   RECEIVER1D_v19_RUNNER_PLOT_NODES=25
#   RECEIVER1D_v19_FIT_NODES=15
#   RECEIVER1D_v19_FIT_ITERATIONS=100
#   RECEIVER1D_v19_FIT_SECONDS=360
# ============================================================================

begin # libraries
    using Statistics
    using StatsPlots
end

include("1D_v19.jl")

begin # runner configuration
    runner_int_v19(name, default) = parse(Int, get(ENV, name, string(default)))
    runner_float_v19(name, default) = parse(Float64, get(ENV, name, string(default)))

    const RUNNER_PLOT_NODES_v19 = runner_int_v19("RECEIVER1D_v19_RUNNER_PLOT_NODES", default_nodes)
    const RUNNER_FIT_NODES_v19 = runner_int_v19("RECEIVER1D_v19_FIT_NODES", 15)
    const RUNNER_FIT_ITERATIONS_v19 = runner_int_v19("RECEIVER1D_v19_FIT_ITERATIONS", 100)
    const RUNNER_FIT_SECONDS_v19 = runner_float_v19("RECEIVER1D_v19_FIT_SECONDS", 360.0)

    const RUNNER_OUTPUT_DIR_v19 = joinpath(@__DIR__, "summaries", "1D_v19")
    const RUNNER_PLOT_DIR_v19 = joinpath(RUNNER_OUTPUT_DIR_v19, "plots")
    const RUNNER_TRANSIENT_DIR_v19 = joinpath(RUNNER_PLOT_DIR_v19, "transients")
    const RUNNER_AXIAL_DIR_v19 = joinpath(RUNNER_PLOT_DIR_v19, "axial_profiles")
    const RUNNER_DIAGNOSTIC_DIR_v19 = joinpath(RUNNER_OUTPUT_DIR_v19, "diagnostics")
    mkpath(RUNNER_OUTPUT_DIR_v19)
    mkpath(RUNNER_PLOT_DIR_v19)
    mkpath(RUNNER_TRANSIENT_DIR_v19)
    mkpath(RUNNER_AXIAL_DIR_v19)
    mkpath(RUNNER_DIAGNOSTIC_DIR_v19)
end

begin # runner helpers
    const RUNNER_SENSOR_COLORS_v19 = Dict(
        :T8 => :blue,
        :T12_perim => :red,
        :T11_perim => :green,
        :T9_core => :orange,
        :T10_core => :brown,
        :T3 => :purple,
        :T2 => :black,
    )

    function save_runner_plot_v19(plot_object, filename)
        path = joinpath(RUNNER_PLOT_DIR_v19, filename)
        savefig(plot_object, path)
        println("[run_1D_v19] Saved plot: $path")
        return path
    end

    absorbed_power_watts_v19(irradiance, scale) =
        ETA_ABS_FIXED_v19 * scale * irradiance * A_frt

    function write_parameters_v19(path, variant)
        parameters = variant.params
        names = (
            "A_Nu", "B_Re", "C_Pr_fixed", "phi_0", "m_rec",
            "front_dep", "scale_456", "scale_304", "scale_256",
            "G_core_perim_W_m_K", "C_perim_eff_J_K", "k_perim_ref_W_m_K", "beta_opt",
        )
        open(path, "w") do io
            println(io, "index,name,value")
            for i in eachindex(parameters)
                println(io, join((i, names[i], parameters[i]), ','))
            end
            println(io, "# fitted,phi_0,$(parameters[4])")
            println(io, "# fitted,m_rec,$(parameters[5])")
            println(io, "# fitted,G_core_perim_W_m_K,$(core_perimeter_conductance_per_length_v19(parameters))")
            println(io, "# fitted,C_perim_eff_J_K,$(perimeter_heat_capacity_total_v19(parameters))")
            println(io, "# fitted,k_perim_ref_W_m_K,$(perimeter_axial_conductivity_v19(900.0, parameters))")
            println(io, "# derived,C_core_eff_J_K,$(core_receiver_heat_capacity_v19(default_nodes))")
            println(io, "# derived,C_participating_eff_J_K,$(participating_assembly_heat_capacity_v19(parameters))")
            println(io, "# reference,measured_C_eff_J_K,$MEASURED_ASSEMBLY_CAPACITANCE_v19")
        end
        return path
    end

    function build_steady_results_v19(params=pnew_v19; keys=sim_key_heat, nodes=RUNNER_PLOT_NODES_v19)
        results = NamedTuple[]
        weights = solar_weights_v5(params[13], params[6], nodes)
        for simulation_id in keys
            outputs, result, experiment = solve_case_v19(params, simulation_id; nodes=nodes)
            conditions = simulation_conditions[simulation_id]
            irradiance = conditions[Io]
            scale = absorbed_power_scale_v19(irradiance, params)
            nominal_absorbed = absorbed_power_watts_v19(irradiance, 1.0)
            scaled_absorbed = absorbed_power_watts_v19(irradiance, scale)
            push!(results, (
                simulation_id=simulation_id,
                flow_lpm=conditions[qlpm],
                irradiance=irradiance,
                power_scale=scale,
                nominal_absorbed_W=nominal_absorbed,
                scaled_absorbed_W=scaled_absorbed,
                added_absorbed_W=scaled_absorbed - nominal_absorbed,
                front_cell_source_fraction=weights[1],
                downstream_source_fraction=1.0 - weights[1],
                T8_model=outputs[end, 1],
                T8_experiment=experiment[end, 1],
                T12_perim_model=outputs[end, 2],
                T12_perim_experiment=experiment[end, 2],
                T11_perim_model=outputs[end, 3],
                T11_perim_experiment=experiment[end, 3],
                T9_core_model=outputs[end, 4],
                T9_core_experiment=experiment[end, 4],
                T10_core_model=outputs[end, 5],
                T10_core_experiment=experiment[end, 5],
                T3_model=outputs[end, 6],
                T3_experiment=experiment[end, 6],
                T2_model=outputs[end, 7],
                T2_experiment=experiment[end, 7],
                T12_minus_T8_model=outputs[end, 2] - outputs[end, 1],
                T12_minus_T8_experiment=experiment[end, 2] - experiment[end, 1],
                T11_minus_T8_model=outputs[end, 3] - outputs[end, 1],
                T11_minus_T8_experiment=experiment[end, 3] - experiment[end, 1],
                T12_minus_T9_model=outputs[end, 2] - outputs[end, 4],
                T12_minus_T9_experiment=experiment[end, 2] - experiment[end, 4],
                T11_minus_T10_model=outputs[end, 3] - outputs[end, 5],
                T11_minus_T10_experiment=experiment[end, 3] - experiment[end, 5],
                mean_nu=mean(result.receiver_nusselt[:, end]),
                mean_effectiveness=mean(result.receiver_effectiveness[:, end]),
                receiver_gas_heat=sum(result.receiver_gas_heat[:, end]),
                G_core_perim_W_m_K=core_perimeter_conductance_per_length_v19(params),
                C_perim_eff_J_K=perimeter_heat_capacity_total_v19(params),
                C_participating_eff_J_K=participating_assembly_heat_capacity_v19(params, nodes),
            ))
        end
        return results
    end

    function write_steady_results_v19(path, steady_results)
        fields = propertynames(first(steady_results))
        open(path, "w") do io
            println(io, join(fields, ','))
            for row in steady_results
                println(io, join((getproperty(row, field) for field in fields), ','))
            end
        end
        return path
    end

    function steady_comparison_plot_v19(steady_results; title_suffix="")
        sensors = (:T8, :T12_perim, :T11_perim, :T9_core, :T10_core, :T3, :T2)
        plot_object = plot(
            title="1D_v19 $title_suffix steady-state comparison",
            xlabel="Model temperature (K)",
            ylabel="Experiment temperature (K)",
            legend=:outerright,
            aspect_ratio=:equal,
            xlims=TEMPERATURE_AXIS_LIMITS_v5,
            ylims=TEMPERATURE_AXIS_LIMITS_v5,
        )
        for sensor in sensors
            model_values = [getproperty(row, Symbol(sensor, "_model")) for row in steady_results]
            experiment_values = [getproperty(row, Symbol(sensor, "_experiment")) for row in steady_results]
            scatter!(plot_object, model_values, experiment_values;
                     label=String(sensor), color=RUNNER_SENSOR_COLORS_v19[sensor],
                     ms=4, markerstrokewidth=0)
        end
        lower, upper = TEMPERATURE_AXIS_LIMITS_v5
        plot!(plot_object, [lower, upper], [lower, upper];
              label="1:1", color=:gray, linestyle=:dash)
        return plot_object
    end
end

begin # execution & output generation
    println("[run_1D_v19] Running 2-Zone Core/Perimeter Macro-ECM Runner.")
    println("[run_1D_v19] Fitting heat-transfer, active flow participation, and core-perimeter parameters...")

    fit = calibrate_v19(
        nodes=RUNNER_FIT_NODES_v19,
        maximum_iterations=RUNNER_FIT_ITERATIONS_v19,
        maximum_time_seconds=RUNNER_FIT_SECONDS_v19,
    )
    fitted_params = fit.parameters

    open(joinpath(RUNNER_OUTPUT_DIR_v19, "optimization_summary_1D_v19.txt"), "w") do io
        println(io, "objective=$(fit.objective)")
        println(io, "return_code=$(fit.retcode)")
        println(io, "parameters=$(fit.parameters)")
    end

    variants_v19 = NamedTuple[
        (name="initial_2zone_macro_ecm", params=copy(pnew_v19)),
        (name="fitted_2zone_macro_ecm", params=fitted_params),
    ]

    for variant in variants_v19
        println("[run_1D_v19] Running variant $(variant.name).")
        write_parameters_v19(
            joinpath(RUNNER_OUTPUT_DIR_v19, "parameters_$(variant.name)_1D_v19.csv"),
            variant,
        )
        steady_results = build_steady_results_v19(
            variant.params; nodes=RUNNER_PLOT_NODES_v19,
        )
        write_steady_results_v19(
            joinpath(RUNNER_OUTPUT_DIR_v19, "steady_results_$(variant.name)_1D_v19.csv"),
            steady_results,
        )
        save_runner_plot_v19(
            steady_comparison_plot_v19(steady_results; title_suffix=variant.name),
            "steady_comparison_$(variant.name)_1D_v19.png",
        )
    end

    println("[run_1D_v19] Optimization completed with objective: $(fit.objective)")
    println("[run_1D_v19] Fitted parameters: $fitted_params")
    println("[run_1D_v19] Complete.")
    println("[run_1D_v19] Output directory: $RUNNER_OUTPUT_DIR_v19")
end
