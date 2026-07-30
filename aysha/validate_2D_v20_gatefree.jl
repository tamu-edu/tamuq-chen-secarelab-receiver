# ============================================================================
# Full v20 validation of one gate-free plant and one T3 observer.
#
# This file evaluates a supplied point; it does not rank, fit, or reject
# candidates while solving. Historical thresholds are reported only after
# all train/holdout metrics have been calculated.
# ============================================================================

using Statistics
using SciMLBase: successful_retcode

include("profile_2D_v20_t3_observer.jl")

const GATEFREE_COOLING_IDS_2D_v20 = ("C69", "C80", "C81")
const GATEFREE_NU50_REFERENCE_2D_v20 =
    3.07786e-4 * 50.0^1.44346

function _gatefree_phase_2D_v20(id, cooling)
    cooling && return id == "C81" ?
        "cooling_holdout" : "cooling_training"
    return id in TRAIN_HEATING_2D_v20 ?
        "heating_training" : "heating_holdout"
end

function _gatefree_final_rows_2D_v20(times)
    length(times) >= 2 || error(
        "gate-free validation requires at least two samples",
    )
    threshold = first(times) +
        0.90 * (last(times) - first(times))
    start = something(
        findfirst(time -> time >= threshold, times),
        length(times),
    )
    return start:length(times)
end

function _gatefree_case_metrics_2D_v20(
    case, observer_location, observer_capacity, observer_stem,
)
    size(case.model, 1) == length(case.inputs.times) || error(
        "incomplete solution for $(case.inputs.id)",
    )
    rows = 2:size(case.model, 1)
    final = _gatefree_final_rows_2D_v20(case.inputs.times)
    observer = _observer_prediction_2D_v20(
        case, observer_location, observer_capacity, observer_stem,
    )

    side_error = vec(
        case.model[rows, 1:3] .-
        case.observed[rows, 1:3],
    )
    T2_error = vec(
        case.model[rows, 7] .-
        case.observed[rows, 7],
    )
    T3_error = vec(
        observer.T3[rows] .-
        case.observed[rows, 6],
    )
    side_final_error = vec(
        case.model[final, 1:3] .-
        case.observed[final, 1:3],
    )
    T2_final_error = vec(
        case.model[final, 7] .-
        case.observed[final, 7],
    )
    T3_final_error = vec(
        observer.T3[final] .-
        case.observed[final, 6],
    )
    total_power_scale = max(
        maximum(abs, case.result.total_gas_power_W), 1.0,
    )

    return (
        simulation_id=case.inputs.id,
        phase=_gatefree_phase_2D_v20(
            case.inputs.id, case.inputs.cooling,
        ),
        sample_count=length(rows),
        final_sample_count=length(final),
        final_window_start_s=case.inputs.times[first(final)],
        final_window_end_s=case.inputs.times[last(final)],
        side_value_count=length(side_error),
        T2_value_count=length(T2_error),
        T3_value_count=length(T3_error),
        side_final_value_count=length(side_final_error),
        T2_final_value_count=length(T2_final_error),
        T3_final_value_count=length(T3_final_error),
        side_pooled_rmse_K=sqrt(mean(abs2, side_error)),
        T2_rmse_K=sqrt(mean(abs2, T2_error)),
        T3_observer_rmse_K=sqrt(mean(abs2, T3_error)),
        side_final_bias_K=mean(side_final_error),
        T2_final_bias_K=mean(T2_final_error),
        T3_observer_final_bias_K=mean(T3_final_error),
        T3_observer_t90_error_s=
            V20.get_t90_2D(
                case.inputs.times, observer.T3,
            ) -
            V20.get_t90_2D(
                case.inputs.times, case.observed[:, 6],
            ),
        ode_success=successful_retcode(
            case.result.ode_solution,
        ),
        flow_converged=all(
            case.result.flow_solver_converged,
        ),
        max_equal_pressure_relative_error=maximum(
            case.result.equal_pressure_relative_error,
        ),
        max_gas_reference_relative_error=maximum(
            case.result.gas_reference_relative_error,
        ),
        max_total_enthalpy_residual_W=maximum(
            abs, case.result.total_gas_enthalpy_residual_W,
        ),
        max_total_enthalpy_relative_residual=maximum(
            abs, case.result.total_gas_enthalpy_residual_W,
        ) / total_power_scale,
    )
end

function _gatefree_aggregate_2D_v20(case_rows, phase)
    selected = filter(row -> row.phase == phase, case_rows)
    isempty(selected) && error("no cases found for phase $phase")

    side_count = sum(row.side_value_count for row in selected)
    T2_count = sum(row.T2_value_count for row in selected)
    T3_count = sum(row.T3_value_count for row in selected)
    side_final_count = sum(
        row.side_final_value_count for row in selected
    )
    T2_final_count = sum(
        row.T2_final_value_count for row in selected
    )
    T3_final_count = sum(
        row.T3_final_value_count for row in selected
    )
    side_rmse = sqrt(sum(
        row.side_value_count * row.side_pooled_rmse_K^2
        for row in selected
    ) / side_count)
    T2_rmse = sqrt(sum(
        row.T2_value_count * row.T2_rmse_K^2
        for row in selected
    ) / T2_count)
    T3_rmse = sqrt(sum(
        row.T3_value_count * row.T3_observer_rmse_K^2
        for row in selected
    ) / T3_count)
    side_bias = sum(
        row.side_final_value_count * row.side_final_bias_K
        for row in selected
    ) / side_final_count
    T2_bias = sum(
        row.T2_final_value_count * row.T2_final_bias_K
        for row in selected
    ) / T2_final_count
    T3_bias = sum(
        row.T3_final_value_count *
        row.T3_observer_final_bias_K
        for row in selected
    ) / T3_final_count
    heating = startswith(phase, "heating")
    original_side_limit = heating ? 60.0 : 35.0
    original_T2_limit = heating ? 15.0 : 12.0

    return (
        phase=phase,
        case_count=length(selected),
        side_value_count=side_count,
        T2_value_count=T2_count,
        T3_value_count=T3_count,
        side_pooled_rmse_K=side_rmse,
        T2_pooled_rmse_K=T2_rmse,
        T3_observer_pooled_rmse_K=T3_rmse,
        side_final_bias_K=side_bias,
        T2_final_bias_K=T2_bias,
        T3_observer_final_bias_K=T3_bias,
        T3_observer_t90_mae_s=mean(
            abs(row.T3_observer_t90_error_s)
            for row in selected
        ),
        ode_success=all(row.ode_success for row in selected),
        flow_converged=all(
            row.flow_converged for row in selected
        ),
        max_equal_pressure_relative_error=maximum(
            row.max_equal_pressure_relative_error
            for row in selected
        ),
        max_gas_reference_relative_error=maximum(
            row.max_gas_reference_relative_error
            for row in selected
        ),
        max_total_enthalpy_residual_W=maximum(
            row.max_total_enthalpy_residual_W
            for row in selected
        ),
        max_total_enthalpy_relative_residual=maximum(
            row.max_total_enthalpy_relative_residual
            for row in selected
        ),
        original_side_limit_K=original_side_limit,
        original_T2_limit_K=original_T2_limit,
        original_plant_threshold_pass=(
            side_rmse < original_side_limit &&
            T2_rmse < original_T2_limit
        ),
        original_T3_rmse_limit_K=30.0,
        original_T3_bias_limit_K=20.0,
        original_T3_t90_limit_s=500.0,
        original_T3_performance_threshold_pass=(
            T3_rmse <= 30.0 &&
            abs(T3_bias) <= 20.0 &&
            mean(
                abs(row.T3_observer_t90_error_s)
                for row in selected
            ) <= 500.0
        ),
    )
end

function validate_2D_v20_gatefree(
    exponent, ratio, felt_k, felt_cp,
    power456, power304, power256,
    observer_location, observer_capacity, observer_stem;
    mesh=:nominal, felt_contact=0.30,
    max_points=61, suffix="",
    reltol=5e-4, abstol=1e-4, dtmax=120.0,
)
    mkpath(OUTPUT_DIR_2D_v20)
    location = Symbol(observer_location)
    location in T3_LOCATIONS_2D_v20 || error(
        "unknown T3 location $location; expected one of " *
        "$(T3_LOCATIONS_2D_v20)",
    )
    nu50 = ratio * GATEFREE_NU50_REFERENCE_2D_v20
    p = parameters_2D_v20(
        mesh=mesh,
        source_model=:near_deep,
        deep_fraction=0.90,
        deep_length_m=0.12,
        nu_prefactor=nu50 / 50.0^exponent,
        reynolds_exponent=exponent,
        distributed_tube_flange_h_W_m2_K=0.0,
        probe_capacity_areal_J_m2_K=observer_capacity,
        probe_stem_conductance_areal_W_m2_K=observer_stem,
        felt_conductivity_scale=felt_k,
        felt_heat_capacity_scale=felt_cp,
        felt_contact_scale=felt_contact,
        power_scales=(power456, power304, power256),
        t3_location=location,
    )

    specs = vcat(
        [(id=id, cooling=false) for id in HEATING_IDS_2D_v20],
        [(id=id, cooling=true)
         for id in GATEFREE_COOLING_IDS_2D_v20],
    )
    case_rows = Vector{NamedTuple}(undef, length(specs))
    Threads.@threads for index in eachindex(specs)
        spec = specs[index]
        println(
            "v20 gate-free validation $index/$(length(specs)): " *
            spec.id,
        )
        flush(stdout)
        inputs = case_inputs_2D_v20(
            spec.id; cooling=spec.cooling,
            max_points=max_points,
        )
        case = simulate_case_2D_v20(
            inputs, p;
            initialization=spec.cooling ?
                :side_T2_only : :ambient,
            reltol=reltol, abstol=abstol, dtmax=dtmax,
        )
        case_rows[index] = _gatefree_case_metrics_2D_v20(
            case, location, observer_capacity, observer_stem,
        )
    end

    phase_order = (
        "heating_training", "heating_holdout",
        "cooling_training", "cooling_holdout",
    )
    phase_rows = [
        _gatefree_aggregate_2D_v20(case_rows, phase)
        for phase in phase_order
    ]
    by_phase = Dict(row.phase => row for row in phase_rows)
    observer_inside_original_bracket = (
        200.0 < observer_capacity < 3000.0 &&
        0.0 < observer_stem < 200.0
    )
    cooling_train = by_phase["cooling_training"]
    cooling_holdout = by_phase["cooling_holdout"]
    original_observer_gate = (
        cooling_train.T3_observer_pooled_rmse_K <= 30.0 &&
        cooling_holdout.T3_observer_pooled_rmse_K <= 30.0 &&
        abs(cooling_holdout.T3_observer_final_bias_K) <= 20.0 &&
        cooling_holdout.T3_observer_t90_mae_s <= 500.0 &&
        observer_inside_original_bracket
    )
    summary = (
        mesh=String(mesh),
        reynolds_exponent=exponent,
        Nu50_ratio=ratio,
        Nu50=nu50,
        nu_prefactor=nu50 / 50.0^exponent,
        felt_conductivity_scale=felt_k,
        felt_heat_capacity_scale=felt_cp,
        felt_contact_scale=felt_contact,
        power_scale_456=power456,
        power_scale_304=power304,
        power_scale_256=power256,
        T3_location=String(location),
        T3_capacity_areal_J_m2_K=observer_capacity,
        T3_stem_conductance_areal_W_m2_K=observer_stem,
        all_ode_success=all(row.ode_success for row in case_rows),
        all_flow_converged=all(
            row.flow_converged for row in case_rows
        ),
        original_heating_training_gate=
            by_phase["heating_training"].
                original_plant_threshold_pass,
        original_heating_holdout_gate=
            by_phase["heating_holdout"].
                original_plant_threshold_pass,
        original_cooling_training_gate=
            cooling_train.original_plant_threshold_pass,
        original_cooling_holdout_gate=
            cooling_holdout.original_plant_threshold_pass,
        observer_inside_original_bracket=
            observer_inside_original_bracket,
        original_cooling_observer_gate=original_observer_gate,
        original_all_posthoc_thresholds_pass=(
            all(
                row.original_plant_threshold_pass
                for row in phase_rows
            ) &&
            original_observer_gate &&
            all(row.ode_success for row in case_rows) &&
            all(row.flow_converged for row in case_rows)
        ),
    )

    stem = isempty(suffix) ? "" : suffix
    _write_namedtuples_csv_2D_v20(joinpath(
        OUTPUT_DIR_2D_v20,
        "gatefree_validation_cases$(stem)_2D_v20.csv",
    ), case_rows)
    _write_namedtuples_csv_2D_v20(joinpath(
        OUTPUT_DIR_2D_v20,
        "gatefree_validation_phases$(stem)_2D_v20.csv",
    ), phase_rows)
    _write_namedtuples_csv_2D_v20(joinpath(
        OUTPUT_DIR_2D_v20,
        "gatefree_validation_summary$(stem)_2D_v20.csv",
    ), [summary])
    println("v20 gate-free validation summary: ", summary)
    return (
        summary=summary, phases=phase_rows, cases=case_rows,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 10 || error(
        "usage: exponent ratio felt_k felt_cp p456 p304 p256 " *
        "T3_location T3_capacity T3_stem [points] [suffix]",
    )
    exponent = parse(Float64, ARGS[1])
    ratio = parse(Float64, ARGS[2])
    felt_k = parse(Float64, ARGS[3])
    felt_cp = parse(Float64, ARGS[4])
    p456 = parse(Float64, ARGS[5])
    p304 = parse(Float64, ARGS[6])
    p256 = parse(Float64, ARGS[7])
    location = Symbol(ARGS[8])
    capacity = parse(Float64, ARGS[9])
    stem = parse(Float64, ARGS[10])
    points = length(ARGS) >= 11 ? parse(Int, ARGS[11]) : 61
    suffix = length(ARGS) >= 12 ? ARGS[12] : ""
    validate_2D_v20_gatefree(
        exponent, ratio, felt_k, felt_cp,
        p456, p304, p256,
        location, capacity, stem;
        max_points=points, suffix=suffix,
    )
end
