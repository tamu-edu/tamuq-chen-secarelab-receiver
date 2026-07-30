# ============================================================================
# Full 15-case confirmation of the inheritance-screen optical-extinction
# result. This is a no-fit confirmation, not a calibrated v11 result.
# ============================================================================

include("run_2D_v11_inheritance_sensitivity.jl")

function main_betaopt_confirmation_2D_v11()
    fitted_v9 = sensitivity_mesh_2D_v11(
        unpack_parameters2D(
            fitted_v9_theta_2D_v11(), default_parameters2D(),
        ),
    )
    graetz = with_heat_model_2D_v11(
        fitted_v9;
        model=:graetz,
        nu_fd=3.61,
        temperature_exponent=0.0,
        scale=1.0,
    )
    groove = with_hydraulics_2D_v11(
        graetz;
        offset=-0.5428447496201336,
        resistance_scale=0.9701974617588522,
        old_hot_K=0.0,
        distribution=:equal_pressure,
        groove_diameter=13.0e-3,
        groove_K=184.15800344228472,
        max_iterations=16,
        tolerance=1.0e-5,
    )
    p = with_optical_sensitivity_2D_v11(
        groove; extinction=110.0,
    )
    label = "graetz_equal_groove_betaopt110"
    current_parameters_2D_v11[] = p

    rows = NamedTuple[]
    for id in IDs
        println("Optical-extinction confirmation: $id")
        flush(stdout)
        case = operating_case_2D_v11(id, p)
        result = simulate2D(
            p,
            case.op,
            case.times;
            initial_temperature=case.ambient[1],
        )
        push!(
            rows,
            final_case_row_2D_v11(label, id, case, result),
        )
    end
    slopes = slope_rows_2D_v11(rows)
    summary = aggregate_rows_2D_v11(rows, slopes)
    write_namedtuple_csv_2D_v11(
        joinpath(
            OUTPUT_DIR_2D_v11,
            "betaopt110_confirmation_cases_2D_v11.csv",
        ),
        rows,
    )
    write_namedtuple_csv_2D_v11(
        joinpath(
            OUTPUT_DIR_2D_v11,
            "betaopt110_confirmation_slopes_2D_v11.csv",
        ),
        slopes,
    )
    write_namedtuple_csv_2D_v11(
        joinpath(
            OUTPUT_DIR_2D_v11,
            "betaopt110_confirmation_summary_2D_v11.csv",
        ),
        summary,
    )
    println("Optical-extinction confirmation summary:")
    foreach(println, summary)
    println("Optical-extinction confirmation slopes:")
    foreach(println, slopes)
    return summary
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_betaopt_confirmation_2D_v11()
end
