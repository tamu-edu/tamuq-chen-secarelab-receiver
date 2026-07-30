using Test
using Statistics

include(joinpath(@__DIR__, "..", "run_2D_v11.jl"))

const MESH_CASES_2D_v11 = ("E67", "E72", "E77")
const MESH_SENSORS_2D_v11 = (
    :T8_K, :T12_K, :T11_K, :T9_K, :T10_K, :T3_K, :T2_K,
)

function groove_parameters_2D_v11(base)
    thermal = with_heat_model_2D_v11(
        base;
        model=:graetz,
        nu_fd=3.61,
        temperature_exponent=0.0,
        scale=1.0,
    )
    return with_hydraulics_2D_v11(
        thermal;
        offset=-0.5428447496201336,
        resistance_scale=0.9701974617588522,
        old_hot_K=0.0,
        distribution=:equal_pressure,
        groove_diameter=13.0e-3,
        groove_K=184.15800344228472,
        max_iterations=16,
        tolerance=1.0e-5,
    )
end

function main_mesh_2D_v11()
    nominal_base = unpack_parameters2D(
        fitted_v9_theta_2D_v11(), default_parameters2D(),
    )
    parameter_sets = (
        ("sensitivity", groove_parameters_2D_v11(
            sensitivity_mesh_2D_v11(nominal_base),
        )),
        ("nominal", groove_parameters_2D_v11(nominal_base)),
    )

    rows = NamedTuple[]
    for (mesh, p) in parameter_sets
        for id in MESH_CASES_2D_v11
            println("Running mesh check: $mesh $id")
            flush(stdout)
            case = operating_case_2D_v11(id, p)
            result = simulate2D(
                p,
                case.op,
                case.times;
                initial_temperature=case.ambient[1],
            )
            pred = sensor_predictions2D(result)
            push!(rows, (
                mesh=mesh,
                simulation_id=id,
                T8_K=pred.T8[end],
                T12_K=pred.T12[end],
                T11_K=pred.T11[end],
                T9_K=pred.T9[end],
                T10_K=pred.T10[end],
                T3_K=pred.T3[end],
                T2_K=pred.T2[end],
                dp1_mbar=result.dp1_prediction_mbar[end],
                pressure_equalization_error=
                    result.equal_pressure_relative_error[end],
            ))
        end
    end

    comparison = NamedTuple[]
    for id in MESH_CASES_2D_v11
        coarse = only(filter(
            row -> row.mesh == "sensitivity" &&
                   row.simulation_id == id,
            rows,
        ))
        fine = only(filter(
            row -> row.mesh == "nominal" &&
                   row.simulation_id == id,
            rows,
        ))
        sensor_deltas = [
            getproperty(fine, sensor) - getproperty(coarse, sensor)
            for sensor in MESH_SENSORS_2D_v11
        ]
        push!(comparison, (
            simulation_id=id,
            max_abs_sensor_delta_K=maximum(abs, sensor_deltas),
            rms_sensor_delta_K=sqrt(mean(abs2, sensor_deltas)),
            axial_delta_difference_K=
                (fine.T12_K - fine.T8_K) -
                (coarse.T12_K - coarse.T8_K),
            mid_radial_delta_difference_K=
                (fine.T12_K - fine.T9_K) -
                (coarse.T12_K - coarse.T9_K),
            dp1_difference_mbar=fine.dp1_mbar - coarse.dp1_mbar,
        ))
    end

    mkpath(OUTPUT_DIR_2D_v11)
    write_namedtuple_csv_2D_v11(
        joinpath(OUTPUT_DIR_2D_v11, "mesh_cases_2D_v11.csv"),
        rows,
    )
    write_namedtuple_csv_2D_v11(
        joinpath(OUTPUT_DIR_2D_v11, "mesh_comparison_2D_v11.csv"),
        comparison,
    )

    @testset "v11 representative nominal-mesh confirmation" begin
        @test length(rows) == 6
        @test all(
            isfinite(getproperty(row, sensor))
            for row in rows for sensor in MESH_SENSORS_2D_v11
        )
        @test maximum(
            row.pressure_equalization_error for row in rows
        ) < 2.0e-5
        @test all(isfinite(row.dp1_mbar) for row in rows)
    end

    println("Mesh comparison:")
    foreach(println, comparison)
    return comparison
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_mesh_2D_v11()
end
