# ============================================================================
# diagnose_2D_v16_t8_11mm.jl
#
# No-refit sensitivity of the completed v16 solution to sampling T8 at the
# corrected 11 mm position instead of the 5 mm position used by v9--v16.
# ============================================================================

using Statistics

include(joinpath(@__DIR__, "..", "run_2D_v16.jl"))

function selected_values_2D_v16_t8()
    cooling_path = joinpath(
        OUTPUT_DIR_2D_v16,
        "parameters_cooling_nominal_selected_2D_v16.txt",
    )
    power_path = joinpath(
        OUTPUT_DIR_2D_v16,
        "parameters_power_selected_2D_v16.txt",
    )
    values = Dict{String,Float64}()
    for path in (cooling_path, power_path), line in eachline(path)
        occursin("=", line) || continue
        key, raw = split(line, "="; limit=2)
        value = tryparse(Float64, raw)
        value === nothing || (values[key] = value)
    end
    return (
        contact=values["skin_felt_contact_scale"],
        felt_k=values["felt_conductivity_scale"],
        felt_cp=values["felt_heat_capacity_scale"],
        powers=(
            values["power_scale_456"],
            values["power_scale_304"],
            values["power_scale_256"],
        ),
    )
end

function phase_2D_v16_t8(id, cooling)
    cooling && return "cooling_validation"
    return id in TRAIN_HEATING_2D_v16 ?
        "heating_training" : "heating_validation"
end

function filtered_side_2D_v16_t8(result, z)
    raw = V16.V15._sample_skin2D(result, z)
    return V16.V12._filter_observation(
        result.time, raw,
        result.parameters.observation.side_time_constant_s,
    )
end

function aggregate_2D_v16_t8(rows, phase, key)
    selected = filter(row -> row.phase == phase, rows)
    sensor_rmse = mean(
        sqrt(mean(abs2, getproperty(row, key)[:, index] .-
                       row.observed[:, index]))
        for row in selected, index in 1:7
    )
    steady_mae = mean(
        abs(getproperty(row, key)[end, index] -
            row.observed[end, index])
        for row in selected, index in 1:7
    )
    axial_rmse = sqrt(mean(
        (
            getproperty(row, key)[end, 2] -
            getproperty(row, key)[end, 1] -
            (row.observed[end, 2] - row.observed[end, 1])
        )^2 for row in selected
    ))
    t8_rmse = mean(
        sqrt(mean(abs2, getproperty(row, key)[:, 1] .-
                       row.observed[:, 1]))
        for row in selected
    )
    modeled_mid_peak = count(
        getproperty(row, key)[end, 2] >
        getproperty(row, key)[end, 1]
        for row in selected if phase != "cooling_validation"
    )
    return (
        phase=phase,
        observation=key == :model5 ? "5_mm" : "11_mm",
        mean_sensor_rmse_K=sensor_rmse,
        steady_mae_K=steady_mae,
        axial_rmse_K=axial_rmse,
        mean_T8_rmse_K=t8_rmse,
        modeled_mid_peak_count=modeled_mid_peak,
        case_count=length(selected),
    )
end

function diagnose_2D_v16_t8_11mm(; max_points=121)
    selected = selected_values_2D_v16_t8()
    p = parameters_2D_v16(
        nominal_mesh=true,
        skin_felt_contact_scale=selected.contact,
        felt_conductivity_scale=selected.felt_k,
        felt_heat_capacity_scale=selected.felt_cp,
        power_scales=selected.powers,
    )
    specs = vcat(
        [(id=id, cooling=false) for id in HEATING_IDS_2D_v16],
        [(id=id, cooling=true) for id in COOLING_IDS_2D_v16],
    )
    cases = Vector{Any}(undef, length(specs))
    Threads.@threads for index in eachindex(specs)
        spec = specs[index]
        inputs = case_inputs_2D_v16(
            spec.id; cooling=spec.cooling, max_points=max_points,
        )
        cases[index] = simulate_case_2D_v16(inputs, p)
    end

    rows = NamedTuple[]
    details = NamedTuple[]
    for case in cases
        model5 = copy(case.model)
        model11 = copy(case.model)
        model5[:, 1] .= filtered_side_2D_v16_t8(
            case.result, 5e-3,
        )
        model11[:, 1] .= filtered_side_2D_v16_t8(
            case.result, 11e-3,
        )
        phase = phase_2D_v16_t8(
            case.inputs.id, case.inputs.cooling,
        )
        push!(rows, (
            phase=phase, id=case.inputs.id,
            model5=model5, model11=model11,
            observed=case.observed,
        ))
        push!(details, (
            phase=phase,
            simulation_id=case.inputs.id,
            model_T8_5mm_K=model5[end, 1],
            model_T8_11mm_K=model11[end, 1],
            delta_11_minus_5_K=model11[end, 1] - model5[end, 1],
            observed_T8_K=case.observed[end, 1],
            model_gap_5mm_K=model5[end, 2] - model5[end, 1],
            model_gap_11mm_K=model11[end, 2] - model11[end, 1],
            observed_gap_K=case.observed[end, 2] -
                           case.observed[end, 1],
        ))
    end

    aggregates = NamedTuple[]
    for phase in (
        "heating_training",
        "heating_validation",
        "cooling_validation",
    ), key in (:model5, :model11)
        push!(aggregates, aggregate_2D_v16_t8(rows, phase, key))
    end

    detail_path = joinpath(
        OUTPUT_DIR_2D_v16,
        "t8_position_sensitivity_cases_2D_v16.csv",
    )
    aggregate_path = joinpath(
        OUTPUT_DIR_2D_v16,
        "t8_position_sensitivity_aggregate_2D_v16.csv",
    )
    _write_namedtuples_csv_2D_v16(detail_path, details)
    _write_namedtuples_csv_2D_v16(aggregate_path, aggregates)
    foreach(println, aggregates)
    return (details=details, aggregates=aggregates)
end

if abspath(PROGRAM_FILE) == @__FILE__
    diagnose_2D_v16_t8_11mm()
end
