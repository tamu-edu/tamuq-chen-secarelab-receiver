# Regenerate v19 axial figures with the full continuous model skin profile.

include("validate_2D_v19_staged.jl")

function regenerate_2D_v19_axial_plots()
    axial_dir = joinpath(
        OUTPUT_DIR_2D_v19, "plots", "axial_profiles",
    )
    mkpath(axial_dir)
    p = selected_parameters_2D_v19(mesh=:nominal)
    flow_rows = NamedTuple[]
    for id in HEATING_IDS_2D_v19
        input = case_inputs_2D_v19(id; max_points=121)
        case = simulate_case_2D_v19(
            input, p;
            reltol=5e-4, abstol=1e-4, dtmax=120.0,
        )
        result = case.result
        side_tau =
            result.parameters.observation.side_time_constant_s
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
            title="2D v19 $id: side axial profile",
            size=(720, 520),
        )
        scatter!(
            axial, [11.0, 58.0, 107.0],
            case.observed[end, 1:3];
            marker=:circle, markersize=6,
            label="measured thermocouples",
        )
        savefig(axial, joinpath(axial_dir, "$(id)_axial.png"))
        push!(flow_rows, (
            simulation_id=id,
            phase="heating",
            nonconverged_points=count(
                !, result.flow_solver_converged,
            ),
            max_pressure_relative_error=maximum(
                result.equal_pressure_relative_error,
            ),
            max_gas_reference_relative_error=maximum(
                result.gas_reference_relative_error,
            ),
            max_iterations=maximum(
                result.flow_solver_iterations,
            ),
        ))
        println("v19 continuous axial plot $id complete")
        flush(stdout)
    end
    for id in COOLING_IDS_2D_v19
        input = case_inputs_2D_v19(
            id; cooling=true, max_points=121,
        )
        case = simulate_case_2D_v19(
            input, p;
            reltol=5e-4, abstol=1e-4, dtmax=120.0,
        )
        result = case.result
        push!(flow_rows, (
            simulation_id=id,
            phase="cooling",
            nonconverged_points=count(
                !, result.flow_solver_converged,
            ),
            max_pressure_relative_error=maximum(
                result.equal_pressure_relative_error,
            ),
            max_gas_reference_relative_error=maximum(
                result.gas_reference_relative_error,
            ),
            max_iterations=maximum(
                result.flow_solver_iterations,
            ),
        ))
        println("v19 flow-convergence audit $id complete")
        flush(stdout)
    end
    _write_namedtuples_csv_2D_v19(joinpath(
        OUTPUT_DIR_2D_v19,
        "flow_convergence_audit_2D_v19.csv",
    ), flow_rows)
end

if abspath(PROGRAM_FILE) == @__FILE__
    regenerate_2D_v19_axial_plots()
end
