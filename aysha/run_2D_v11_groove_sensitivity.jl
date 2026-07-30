# ============================================================================
# Rear-groove free-diameter sensitivity for the accepted v11 test matrix.
# ============================================================================

include("run_2D_v11.jl")

const GROOVE_DIAMETERS_MM_2D_v11 = (12.0, 13.0, 14.0)

function main_groove_sensitivity_2D_v11()
    mkpath(OUTPUT_DIR_2D_v11)
    base = sensitivity_mesh_2D_v11(
        unpack_parameters2D(
            fitted_v9_theta_2D_v11(), default_parameters2D(),
        ),
    )

    case_rows = NamedTuple[]
    fit_rows = NamedTuple[]
    fit_by_label = Dict{String, Any}()

    for diameter_mm in GROOVE_DIAMETERS_MM_2D_v11
        label = "graetz_equal_groove_d$(Int(diameter_mm))"
        diameter_m = diameter_mm * 1.0e-3
        println("Fitting cold hydraulics for free diameter $diameter_mm mm")
        flush(stdout)
        fitted = fit_cold_hydraulics_2D_v11(
            base;
            allow_groove=true,
            groove_diameter=diameter_m,
        )
        fit_by_label[label] = fitted
        push!(fit_rows, (
            variant=label,
            free_diameter_mm=diameter_mm,
            offset_mbar=fitted.theta[1],
            resistance_scale=fitted.theta[2],
            groove_K=fitted.theta[3],
            cold_rmse_mbar=fitted.rmse_mbar,
            cold_aicc=fitted.aicc,
        ))

        p = with_heat_model_2D_v11(
            fitted.parameters;
            model=:graetz,
            nu_fd=3.61,
            temperature_exponent=0.0,
            scale=1.0,
        )
        current_parameters_2D_v11[] = p
        for id in IDs
            println("  $label $id")
            flush(stdout)
            case = operating_case_2D_v11(id, p)
            result = simulate2D(
                p,
                case.op,
                case.times;
                initial_temperature=case.ambient[1],
            )
            push!(
                case_rows,
                final_case_row_2D_v11(label, id, case, result),
            )
        end
    end

    slopes = slope_rows_2D_v11(case_rows)
    aggregate = aggregate_rows_2D_v11(case_rows, slopes)
    summary_rows = NamedTuple[]
    for row in aggregate
        fitted = fit_by_label[row.variant]
        push!(summary_rows, merge((
            free_diameter_mm=
                fitted.parameters.hydraulics.groove_free_diameter * 1000.0,
            offset_mbar=fitted.theta[1],
            resistance_scale=fitted.theta[2],
            groove_K=fitted.theta[3],
            cold_rmse_mbar=fitted.rmse_mbar,
            cold_aicc=fitted.aicc,
        ), row))
    end

    write_namedtuple_csv_2D_v11(
        joinpath(
            OUTPUT_DIR_2D_v11,
            "groove_geometry_fits_2D_v11.csv",
        ),
        fit_rows,
    )
    write_namedtuple_csv_2D_v11(
        joinpath(
            OUTPUT_DIR_2D_v11,
            "groove_geometry_cases_2D_v11.csv",
        ),
        case_rows,
    )
    write_namedtuple_csv_2D_v11(
        joinpath(
            OUTPUT_DIR_2D_v11,
            "groove_geometry_slopes_2D_v11.csv",
        ),
        slopes,
    )
    write_namedtuple_csv_2D_v11(
        joinpath(
            OUTPUT_DIR_2D_v11,
            "groove_geometry_summary_2D_v11.csv",
        ),
        summary_rows,
    )

    println("Groove-geometry summary:")
    foreach(println, summary_rows)
    return summary_rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_groove_sensitivity_2D_v11()
end
