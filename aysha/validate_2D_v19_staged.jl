# ============================================================================
# validate_2D_v19_staged.jl
# Full untouched v19 validation, plots, closure, mesh, and null-branch checks.
# ============================================================================

using Statistics
using Plots

include("calibrate_2D_v19_stageD.jl")

function selected_parameters_2D_v19(; mesh=:nominal, null_stageC=false)
    values = _read_keyvalues_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19, "parameters_selected_2D_v19.txt",
    ))
    number(key) = parse(Float64, values[key])
    return parameters_2D_v19(
        mesh=mesh,
        source_model=Symbol(values["source_model"]),
        deep_fraction=number("deep_fraction"),
        deep_length_m=number("deep_length_m"),
        nu_prefactor=number("nu_prefactor"),
        reynolds_exponent=number("reynolds_exponent"),
        rear_tube_flange_contact_h_W_m2_K=null_stageC ?
            0.0 :
            number("rear_tube_flange_contact_h_W_m2_K"),
        terminal_tube_flange_scale=null_stageC ? 0.0 : nothing,
        probe_capacity_areal_J_m2_K=null_stageC ?
            1000.0 :
            number("probe_capacity_areal_J_m2_K"),
        probe_stem_conductance_areal_W_m2_K=null_stageC ?
            20.0 :
            number("probe_stem_conductance_areal_W_m2_K"),
        felt_conductivity_scale=
            number("felt_conductivity_scale"),
        felt_heat_capacity_scale=
            number("felt_heat_capacity_scale"),
        felt_contact_scale=number("felt_contact_scale"),
        power_scales=(
            number("power_scale_456"),
            number("power_scale_304"),
            number("power_scale_256"),
        ),
    )
end

function _phase_2D_v19(id, cooling)
    cooling && return id == "C81" ?
        "cooling_validation" : "cooling_training"
    return id in TRAIN_HEATING_2D_v19 ?
        "heating_training" : "heating_validation"
end

_rmse_2D_v19(model, observed) =
    sqrt(mean(abs2, model .- observed))

function _aggregate_2D_v19(sensor_rows, final_rows)
    aggregate_rows = NamedTuple[]
    for phase in unique(row.phase for row in sensor_rows)
        selected = filter(row -> row.phase == phase, sensor_rows)
        primary = filter(
            row -> row.sensor in ("T8", "T12", "T11", "T3", "T2"),
            selected,
        )
        side = filter(
            row -> row.sensor in ("T8", "T12", "T11"), selected,
        )
        air = filter(row -> row.sensor == "T3", selected)
        felt = filter(row -> row.sensor == "T2", selected)
        phase_final = filter(row -> row.phase == phase, final_rows)
        push!(aggregate_rows, (
            phase=phase,
            primary_mean_sensor_rmse_K=
                mean(row.rmse_K for row in primary),
            primary_pooled_rmse_K=sqrt(mean(
                row.rmse_K^2 for row in primary
            )),
            side_mean_rmse_K=mean(row.rmse_K for row in side),
            air_mean_rmse_K=mean(row.rmse_K for row in air),
            felt_mean_rmse_K=mean(row.rmse_K for row in felt),
            axial_inversion_rmse_K=sqrt(mean(
                (
                    row.model_T12_minus_T8_K -
                    row.observed_T12_minus_T8_K
                )^2 for row in phase_final
            )),
        ))
    end
    return aggregate_rows
end

function _plot_case_2D_v19(case, transient_dir, axial_dir)
    primary_indices = (1, 2, 3, 6, 7)
    colors = (:red, :orange, :purple, :blue, :brown)
    transient = plot(
        xlabel="Time (s)", ylabel="Temperature (K)",
        title="2D v19 $(case.inputs.id): primary observables",
        size=(920, 650), legend=:outerright,
    )
    for (local_index, sensor_index) in enumerate(primary_indices)
        sensor = SENSOR_NAMES_2D_v19[sensor_index]
        plot!(
            transient, case.inputs.times,
            case.observed[:, sensor_index];
            color=colors[local_index], linestyle=:dash,
            label="$(sensor) measured",
        )
        plot!(
            transient, case.inputs.times,
            case.model[:, sensor_index];
            color=colors[local_index], linewidth=2,
            label="$(sensor) model",
        )
    end
    savefig(transient, joinpath(
        transient_dir, "$(case.inputs.id)_transient.png",
    ))
    case.inputs.cooling && return
    observed = case.observed[end, 1:3]
    result = case.result
    side_tau = result.parameters.observation.side_time_constant_s
    model_profile = [
        V19.V12._filter_observation(
            result.time,
            vec(result.skin_temperature[index, :]),
            side_tau,
        )[end]
        for index in axes(result.skin_temperature, 1)
    ]
    axial = plot(
        1e3 .* result.z_solid, model_profile;
        linewidth=2.5, label="model",
        xlabel="Axial position (mm)",
        ylabel="Side-wall temperature (K)",
        title="2D v19 $(case.inputs.id): side axial profile",
        size=(720, 520),
    )
    scatter!(
        axial, [11.0, 58.0, 107.0], observed;
        marker=:circle, markersize=6,
        label="measured thermocouples",
    )
    savefig(axial, joinpath(
        axial_dir, "$(case.inputs.id)_axial.png",
    ))
end

function _mesh_confirmations_2D_v19(max_points)
    specs = (
        (id="E67", cooling=false),
        (id="E71", cooling=false),
        (id="E75", cooling=false),
        (id="E80", cooling=false),
        (id="C69", cooling=true),
        (id="C81", cooling=true),
    )
    nominal = selected_parameters_2D_v19(mesh=:nominal)
    refined = selected_parameters_2D_v19(mesh=:refined)
    rows = NamedTuple[]
    for spec in specs
        input = case_inputs_2D_v19(
            spec.id; cooling=spec.cooling,
            max_points=max_points,
        )
        m = simulate_case_2D_v19(
            input, nominal;
            reltol=5e-4, abstol=1e-4, dtmax=120.0,
        )
        f = simulate_case_2D_v19(
            input, refined;
            reltol=5e-4, abstol=1e-4, dtmax=120.0,
        )
        primary = [1, 2, 3, 6, 7]
        difference = f.model[:, primary] .-
                     m.model[:, primary]
        audit_m = V19.exchange_energy_audit2D(
            m.result.ode_solution.u[end], nominal,
            m.result.operating_condition, m.result.time[end],
        )
        audit_f = V19.exchange_energy_audit2D(
            f.result.ode_solution.u[end], refined,
            f.result.operating_condition, f.result.time[end],
        )
        push!(rows, (
            simulation_id=spec.id,
            primary_history_rms_K=sqrt(mean(abs2, difference)),
            primary_history_max_abs_K=maximum(abs, difference),
            nominal_exchange_energy_relative_residual=
                audit_m.relative_residual,
            refined_exchange_energy_relative_residual=
                audit_f.relative_residual,
            nominal_flow_converged=
                all(m.result.flow_solver_converged),
            refined_flow_converged=
                all(f.result.flow_solver_converged),
        ))
        println("v19 mesh confirmation ", spec.id, " complete")
        flush(stdout)
    end
    return rows
end

function _probe_sampling_check_2D_v19()
    p = selected_parameters_2D_v19(mesh=:nominal)
    rows = NamedTuple[]
    for point_count in (61, 121, 241)
        input = case_inputs_2D_v19(
            "C81"; cooling=true, max_points=point_count,
        )
        case = simulate_case_2D_v19(
            input, p;
            reltol=5e-4, abstol=1e-4, dtmax=120.0,
        )
        push!(rows, (
            requested_points=point_count,
            actual_points=length(input.times),
            final_T3_K=case.model[end, 6],
            T3_t90_s=V19.get_t90_2D(
                input.times, case.model[:, 6],
            ),
            T3_rmse_K=_rmse_2D_v19(
                case.model[:, 6], case.observed[:, 6],
            ),
        ))
    end
    return rows
end

function _solver_tolerance_check_2D_v19()
    p = selected_parameters_2D_v19(mesh=:nominal)
    input = case_inputs_2D_v19("E75"; max_points=61)
    loose = simulate_case_2D_v19(
        input, p;
        reltol=5e-4, abstol=1e-4, dtmax=120.0,
    )
    tight = simulate_case_2D_v19(
        input, p;
        reltol=2e-5, abstol=2e-6, dtmax=60.0,
    )
    primary = [1, 2, 3, 6, 7]
    difference = tight.model[:, primary] .-
                 loose.model[:, primary]
    return (
        simulation_id="E75",
        primary_history_rms_K=sqrt(mean(abs2, difference)),
        primary_history_max_abs_K=maximum(abs, difference),
    )
end

function _final_acceptance_2D_v19(
    aggregate_rows, ledger_rows, mesh_rows,
    null_rows, sampling_rows, tolerance_row, cases,
)
    stageA = parse(
        Bool, _read_keyvalues_2D_v19(joinpath(
            OUTPUT_DIR_2D_v19,
            "stageA_selected_2D_v19.txt",
        ))["stageA_pass"],
    )
    stageB = parse(
        Bool, _read_keyvalues_2D_v19(joinpath(
            OUTPUT_DIR_2D_v19,
            "stageB_selected_2D_v19.txt",
        ))["stageB_pass"],
    )
    stageC = parse(
        Bool, _read_keyvalues_2D_v19(joinpath(
            OUTPUT_DIR_2D_v19,
            "stageC_selected_2D_v19.txt",
        ))["stageC_pass"],
    )
    selected = _read_keyvalues_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19,
        "parameters_selected_2D_v19.txt",
    ))
    felt_pass = parse(Bool, selected["stageD_felt_pass"])
    power_pass = parse(Bool, selected["stageD_power_pass"])
    mesh_pass = all(
        row.primary_history_rms_K < 10.0 &&
        row.primary_history_max_abs_K < 20.0 &&
        row.nominal_exchange_energy_relative_residual <= 1e-10 &&
        row.refined_exchange_energy_relative_residual <= 1e-10 &&
        row.nominal_flow_converged &&
        row.refined_flow_converged
        for row in mesh_rows
    )
    ledger_identity_pass = all(
        row.relative_residual <= 1e-8 &&
        abs(row.source_power_error) <= 1e-10
        for row in ledger_rows
    )
    independent_energy_audit_pass = all(
        row.independent_exchange_relative_residual <= 1e-10
        for row in ledger_rows
    )
    flow_solver_pass = all(
        all(case.result.flow_solver_converged) &&
        maximum(
            case.result.gas_reference_relative_error
        ) < 1e-4 &&
        maximum(
            case.result.equal_pressure_relative_error
        ) < 1e-4
        for case in cases
    )
    null_sensitive = maximum(
        row.primary_history_rms_K for row in null_rows
    ) > 10.0
    sampling_final_span = maximum(
        row.final_T3_K for row in sampling_rows
    ) - minimum(row.final_T3_K for row in sampling_rows)
    sampling_t90_span = maximum(
        row.T3_t90_s for row in sampling_rows
    ) - minimum(row.T3_t90_s for row in sampling_rows)
    probe_sampling_pass = sampling_final_span <= 2.0 &&
                          sampling_t90_span <= 120.0
    solver_tolerance_pass =
        tolerance_row.primary_history_rms_K <= 2.0 &&
        tolerance_row.primary_history_max_abs_K <= 5.0
    accepted = stageA && stageB && stageC &&
        felt_pass && power_pass && mesh_pass &&
        ledger_identity_pass &&
        independent_energy_audit_pass &&
        flow_solver_pass && probe_sampling_pass &&
        solver_tolerance_pass &&
        !null_sensitive
    path = joinpath(
        OUTPUT_DIR_2D_v19, "final_acceptance_2D_v19.txt",
    )
    open(path, "w") do io
        println(io, "stageA_pass=", stageA)
        println(io, "stageB_pass=", stageB)
        println(io, "stageC_pass=", stageC)
        println(io, "stageD_felt_pass=", felt_pass)
        println(io, "stageD_power_pass=", power_pass)
        println(io, "mesh_pass=", mesh_pass)
        println(io, "ledger_identity_pass=", ledger_identity_pass)
        println(io, "independent_energy_audit_pass=",
                independent_energy_audit_pass)
        println(io, "flow_solver_pass=", flow_solver_pass)
        println(io, "stageC_null_sensitive=", null_sensitive)
        println(io, "probe_sampling_pass=", probe_sampling_pass)
        println(io, "solver_tolerance_pass=",
                solver_tolerance_pass)
        println(io, "solver_tolerance_rms_K=",
                tolerance_row.primary_history_rms_K)
        println(io, "solver_tolerance_max_abs_K=",
                tolerance_row.primary_history_max_abs_K)
        println(io, "probe_sampling_final_span_K=",
                sampling_final_span)
        println(io, "probe_sampling_t90_span_s=",
                sampling_t90_span)
        println(io, "v19_accepted=", accepted)
    end
    return (
        stageA=stageA, stageB=stageB, stageC=stageC,
        felt=felt_pass, power=power_pass, mesh=mesh_pass,
        ledger_identity=ledger_identity_pass,
        independent_energy=independent_energy_audit_pass,
        flow_solver=flow_solver_pass,
        null_sensitive=null_sensitive,
        probe_sampling=probe_sampling_pass,
        solver_tolerance=solver_tolerance_pass,
        accepted=accepted,
    )
end

function _null_stageC_sensitivity_2D_v19(max_points)
    specs = (
        (id="E69", cooling=false),
        (id="E74", cooling=false),
        (id="E79", cooling=false),
        (id="C69", cooling=true),
        (id="C81", cooling=true),
    )
    selected = selected_parameters_2D_v19(mesh=:nominal)
    null = selected_parameters_2D_v19(
        mesh=:nominal, null_stageC=true,
    )
    rows = NamedTuple[]
    for spec in specs
        input = case_inputs_2D_v19(
            spec.id; cooling=spec.cooling,
            max_points=max_points,
        )
        base_case = simulate_case_2D_v19(
            input, selected;
            reltol=5e-4, abstol=1e-4, dtmax=120.0,
        )
        null_case = simulate_case_2D_v19(
            input, null;
            reltol=5e-4, abstol=1e-4, dtmax=120.0,
        )
        primary = [1, 2, 3, 6, 7]
        difference = null_case.model[:, primary] .-
                     base_case.model[:, primary]
        push!(rows, (
            simulation_id=spec.id,
            primary_history_rms_K=sqrt(mean(abs2, difference)),
            primary_history_max_abs_K=maximum(abs, difference),
            T3_history_rms_K=_rmse_2D_v19(
                null_case.model[:, 6],
                base_case.model[:, 6],
            ),
        ))
    end
    return rows
end

function run_validation_2D_v19(; max_points=121)
    mkpath(OUTPUT_DIR_2D_v19)
    plot_root = joinpath(OUTPUT_DIR_2D_v19, "plots")
    transient_dir = joinpath(plot_root, "transients")
    axial_dir = joinpath(plot_root, "axial_profiles")
    foreach(mkpath, (plot_root, transient_dir, axial_dir))
    p = selected_parameters_2D_v19(mesh=:nominal)
    specs = vcat(
        [(id=id, cooling=false) for id in HEATING_IDS_2D_v19],
        [(id=id, cooling=true) for id in COOLING_IDS_2D_v19],
    )
    cases = Vector{Any}(undef, length(specs))
    Threads.@threads for index in eachindex(specs)
        spec = specs[index]
        inputs = case_inputs_2D_v19(
            spec.id; cooling=spec.cooling,
            max_points=max_points,
        )
        cases[index] = simulate_case_2D_v19(
            inputs, p;
            reltol=5e-4, abstol=1e-4, dtmax=120.0,
        )
        println("v19 validation ", spec.id, " complete")
        flush(stdout)
    end

    sensor_rows = NamedTuple[]
    final_rows = NamedTuple[]
    ledger_rows = NamedTuple[]
    for case in cases
        phase = _phase_2D_v19(
            case.inputs.id, case.inputs.cooling,
        )
        for (index, sensor) in enumerate(SENSOR_NAMES_2D_v19)
            model = case.model[:, index]
            observed = case.observed[:, index]
            push!(sensor_rows, (
                phase=phase, simulation_id=case.inputs.id,
                sensor=String(sensor),
                rmse_K=_rmse_2D_v19(model, observed),
                steady_error_K=model[end] - observed[end],
                t90_error_s=V19.get_t90_2D(
                    case.inputs.times, model,
                ) - V19.get_t90_2D(
                    case.inputs.times, observed,
                ),
            ))
        end
        m = case.model[end, :]
        o = case.observed[end, :]
        result = case.result
        k = length(result.time)
        push!(final_rows, (
            phase=phase, simulation_id=case.inputs.id,
            flow_lpm=mean(case.inputs.flow),
            model_T8_K=m[1], observed_T8_K=o[1],
            model_T12_K=m[2], observed_T12_K=o[2],
            model_T11_K=m[3], observed_T11_K=o[3],
            model_T9_K=m[4], observed_T9_K=o[4],
            model_T10_K=m[5], observed_T10_K=o[5],
            model_T3_K=m[6], observed_T3_K=o[6],
            model_T2_K=m[7], observed_T2_K=o[7],
            model_T3_gas_raw_K=case.predictions.T3_gas_raw[end],
            model_T3_tube_K=case.predictions.T3_tube[end],
            model_T12_minus_T8_K=m[2] - m[1],
            observed_T12_minus_T8_K=o[2] - o[1],
            model_T11_minus_T12_K=m[3] - m[2],
            observed_T11_minus_T12_K=o[3] - o[2],
            core_to_side_flow_ratio=
                result.channel_mass_flow_kg_s[
                    result.core_group, k
                ] /
                result.channel_mass_flow_kg_s[
                    result.side_group, k
                ],
            dp1_model_mbar=result.dp1_prediction_mbar[end],
            dp1_observed_mbar=case.inputs.dp1[end],
        ))
        ledger = V19.energy_rate_ledger2D(
            result.ode_solution.u[end], p,
            result.operating_condition, result.time[end],
        )
        exchange_audit = V19.exchange_energy_audit2D(
            result.ode_solution.u[end], p,
            result.operating_condition, result.time[end],
        )
        push!(ledger_rows, (
            phase=phase, simulation_id=case.inputs.id,
            residual_W=ledger.residual,
            relative_residual=abs(ledger.residual) /
                max(abs(ledger.solar_deposited), 1.0),
            solar_deposited_W=ledger.solar_deposited,
            integrated_gas_removal_W=
                ledger.integrated_gas_removal_W,
            distributed_tube_flange_loss_W=
                ledger.distributed_tube_flange_loss_W,
            source_power_error=ledger.source_power_error,
            independent_exchange_residual_W=
                exchange_audit.residual_W,
            independent_exchange_relative_residual=
                exchange_audit.relative_residual,
        ))
        _plot_case_2D_v19(
            case, transient_dir, axial_dir,
        )
    end

    aggregate_rows = _aggregate_2D_v19(
        sensor_rows, final_rows,
    )
    _write_namedtuples_csv_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19,
        "staged_sensor_metrics_2D_v19.csv",
    ), sensor_rows)
    _write_namedtuples_csv_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19,
        "staged_final_profiles_2D_v19.csv",
    ), final_rows)
    _write_namedtuples_csv_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19,
        "staged_aggregate_2D_v19.csv",
    ), aggregate_rows)
    _write_namedtuples_csv_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19, "energy_ledgers_2D_v19.csv",
    ), ledger_rows)

    parity = plot(
        xlabel="Measured temperature (K)",
        ylabel="Modeled temperature (K)",
        title="2D v19 primary-observable parity",
        size=(760, 680), legend=:topleft,
    )
    phase_colors = Dict(
        "heating_training" => :royalblue,
        "heating_validation" => :darkorange,
        "cooling_training" => :seagreen,
        "cooling_validation" => :purple,
    )
    all_values = Float64[]
    for phase in keys(phase_colors)
        xs = Float64[]
        ys = Float64[]
        for row in filter(r -> r.phase == phase, final_rows)
            for sensor in ("T8", "T12", "T11", "T3", "T2")
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
            parity, xs, ys;
            label=replace(phase, "_" => " "),
            color=phase_colors[phase], alpha=0.80,
        )
    end
    bounds = extrema(all_values)
    plot!(
        parity, collect(bounds), collect(bounds);
        color=:black, linestyle=:dash, label="1:1",
    )
    savefig(parity, joinpath(
        plot_root, "parity_primary_2D_v19.png",
    ))

    mesh_rows = _mesh_confirmations_2D_v19(
        min(max_points, 61),
    )
    null_rows = _null_stageC_sensitivity_2D_v19(
        min(max_points, 61),
    )
    sampling_rows = _probe_sampling_check_2D_v19()
    tolerance_row = _solver_tolerance_check_2D_v19()
    _write_namedtuples_csv_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19,
        "mesh_confirmation_2D_v19.csv",
    ), mesh_rows)
    _write_namedtuples_csv_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19,
        "stageC_null_sensitivity_2D_v19.csv",
    ), null_rows)
    _write_namedtuples_csv_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19,
        "probe_sampling_check_2D_v19.csv",
    ), sampling_rows)
    _write_namedtuples_csv_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19,
        "solver_tolerance_check_2D_v19.csv",
    ), [tolerance_row])
    acceptance = _final_acceptance_2D_v19(
        aggregate_rows, ledger_rows, mesh_rows,
        null_rows, sampling_rows, tolerance_row, cases,
    )
    foreach(println, aggregate_rows)
    return (
        cases=cases, sensors=sensor_rows, finals=final_rows,
        aggregates=aggregate_rows, ledgers=ledger_rows,
        mesh=mesh_rows, null_stageC=null_rows,
        probe_sampling=sampling_rows,
        solver_tolerance=tolerance_row,
        acceptance=acceptance,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_validation_2D_v19()
end
