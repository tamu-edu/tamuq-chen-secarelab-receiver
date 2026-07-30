# ============================================================================
# Frozen-parameter v19/v20 comparison.
#
# This isolates the variable-cp enthalpy correction. All plant parameters,
# including the rejected v19 distributed rear-tube sink, are identical.
# ============================================================================

using Statistics

include("run_2D_v20.jl")

const COMPARISON_IDS_2D_v20 = (
    ("E67", false), ("E72", false), ("E76", false),
    ("E80", false), ("C69", true), ("C81", true),
)

function _frozen_parameter_kwargs_2D_v20(mesh)
    return (
        mesh=mesh,
        source_model=:near_deep,
        deep_fraction=0.90,
        deep_length_m=0.12,
        nu_prefactor=3.80e-4,
        reynolds_exponent=1.54,
        probe_capacity_areal_J_m2_K=3000.0,
        probe_stem_conductance_areal_W_m2_K=0.0,
        felt_conductivity_scale=0.70,
        felt_heat_capacity_scale=0.75,
        felt_contact_scale=0.30,
        power_scales=(1.05, 1.23, 0.84),
    )
end

function _rmse_2D_v20(values)
    return sqrt(mean(abs2, values))
end

function _result_at_location_2D_v20(result, location)
    parameters = V20.ModelParameters2D(
        base=result.parameters.base,
        t3_location=V20.T3LocationParameters2D(
            location=location,
        ),
    )
    return V20.SimulationResult2D(
        result.base_result, parameters,
        result.diagnostics, result.ode_solution,
    )
end

function compare_2D_v19_v20_enthalpy(; mesh=:nominal, max_points=61)
    mkpath(OUTPUT_DIR_2D_v20)
    kwargs = _frozen_parameter_kwargs_2D_v20(mesh)
    p19 = parameters_2D_v19(
        ; kwargs...,
        rear_tube_flange_contact_h_W_m2_K=400.0,
    )
    p20 = parameters_2D_v20(
        ; kwargs...,
        distributed_tube_flange_h_W_m2_K=400.0,
        t3_location=:rear_003,
    )
    rows = NamedTuple[]
    location_rows = NamedTuple[]
    output = joinpath(
        OUTPUT_DIR_2D_v20,
        "frozen_v19_v20_enthalpy_comparison.csv",
    )
    location_output = joinpath(
        OUTPUT_DIR_2D_v20,
        "frozen_t3_location_comparison.csv",
    )
    for (index, (id, cooling)) in enumerate(COMPARISON_IDS_2D_v20)
        println("v20 frozen comparison $index/$(length(COMPARISON_IDS_2D_v20)): $id")
        flush(stdout)
        inputs = case_inputs_2D_v20(
            id; cooling=cooling, max_points=max_points,
        )
        case19 = simulate_case_2D_v19(
            inputs, p19; reltol=5e-4, abstol=1e-4, dtmax=120.0,
        )
        case20 = simulate_case_2D_v20(
            inputs, p20; reltol=5e-4, abstol=1e-4, dtmax=120.0,
        )
        delta = case20.model .- case19.model
        audit = V20.enthalpy_transport_audit2D(case20.result)
        push!(rows, (
            ID=id, cooling=cooling, mesh=String(mesh),
            model_all_rms_delta_K=_rmse_2D_v20(delta),
            model_all_max_delta_K=maximum(abs, delta),
            side_rms_delta_K=_rmse_2D_v20(delta[:, 1:3]),
            T3_rms_delta_K=_rmse_2D_v20(delta[:, 6]),
            T2_rms_delta_K=_rmse_2D_v20(delta[:, 7]),
            receiver_exit_rms_delta_K=_rmse_2D_v20(
                case20.result.rear_gas_temperature[1, :] .-
                case19.result.rear_gas_temperature[1, :],
            ),
            receiver_exit_final_delta_K=
                case20.result.rear_gas_temperature[1, end] -
                case19.result.rear_gas_temperature[1, end],
            enthalpy_residual_max_W=
                audit.maximum_absolute_residual_W,
            enthalpy_residual_relative=
                audit.maximum_relative_residual,
            flow_converged=all(
                case20.result.flow_solver_converged,
            ),
        ))
        for location in (
            :receiver_136, :receiver_exit, :rear_003,
        )
            located = _result_at_location_2D_v20(
                case20.result, location,
            )
            initial_T3 = cooling ?
                Float64(observation(
                    inputs.data, id, "_T3",
                )[1]) : nothing
            prediction = V20.sensor_predictions2D(
                located; initial_T3=initial_T3,
            )
            observed = case20.observed[:, 6]
            error = prediction.T3 .- observed
            raw_error = prediction.T3_gas_raw .- observed
            final = max(2, floor(Int, 0.9length(error))):length(error)
            push!(location_rows, (
                ID=id, cooling=cooling,
                location=String(location),
                global_z_mm=1000prediction.T3_global_z_m,
                probe_rmse_K=_rmse_2D_v20(error[2:end]),
                probe_final_bias_K=mean(error[final]),
                raw_gas_rmse_K=_rmse_2D_v20(
                    raw_error[2:end],
                ),
                raw_gas_final_bias_K=mean(raw_error[final]),
            ))
        end
        _write_namedtuples_csv_2D_v20(output, rows)
        _write_namedtuples_csv_2D_v20(
            location_output, location_rows,
        )
    end
    return (comparison=rows, locations=location_rows)
end

if abspath(PROGRAM_FILE) == @__FILE__
    mesh = length(ARGS) >= 1 ? Symbol(ARGS[1]) : :nominal
    points = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 61
    compare_2D_v19_v20_enthalpy(
        ; mesh=mesh, max_points=points,
    )
end
