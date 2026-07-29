using LinearAlgebra
using Statistics
using Test

include(joinpath(@__DIR__, "..", "run_2D_v9_full.jl"))

const DIAGNOSTIC_OUTPUT_DIR_2D_v9 = joinpath(
    @__DIR__, "..", "summaries", "2D_v9",
)

function read_fitted_parameters_full_2D_v9()
    path = joinpath(
        DIAGNOSTIC_OUTPUT_DIR_2D_v9, "parameters_fitted_2D_v9.csv",
    )
    lines = readlines(path)
    length(lines) == 13 ||
        error("Expected header plus 12 parameter rows in $path")
    theta = [parse(Float64, split(line, ',')[3]) for line in lines[2:end]]

    summary_lines = readlines(joinpath(
        DIAGNOSTIC_OUTPUT_DIR_2D_v9, "full_test_summary_2D_v9.txt",
    ))
    prefix = "hot_excess_loss_coefficient="
    line = only(filter(row -> startswith(row, prefix), summary_lines))
    coefficient = parse(Float64, split(line, '=')[2])

    calibration = cold_t0_dp1_calibration_2D_v9()
    p0 = with_t0_hydraulics_2D_v9(default_parameters2D(), calibration)
    return with_minor_loss_full_2D_v9(
        unpack_parameters2D(theta, p0), coefficient,
    )
end

function transient_response_full_2D_v9(theta, p_base)
    p = unpack_parameters2D(theta, p_base)
    values = Float64[]
    for id in TRAIN_HEATING_2D_v9
        case = run_full_case_2D_v9(
            id, calibration_mesh_full_2D_v9(p); is_cooling=false,
        )
        for sensor_index in axes(case.model, 2)
            span = max(
                maximum(case.experiment[:, sensor_index]) -
                minimum(case.experiment[:, sensor_index]),
                1.0,
            )
            append!(values, case.model[:, sensor_index] ./ span)
        end
    end
    return values
end

function identifiability_full_2D_v9(p)
    theta = pack_parameters2D(p)
    response0 = transient_response_full_2D_v9(theta, p)
    jacobian = zeros(length(response0), length(FIT_INDICES_2D))

    for (column, parameter_index) in enumerate(FIT_INDICES_2D)
        range_i = UB_2D[parameter_index] - LB_2D[parameter_index]
        step = 0.01 * range_i
        theta_low = copy(theta)
        theta_high = copy(theta)
        theta_low[parameter_index] = max(
            LB_2D[parameter_index], theta[parameter_index] - step,
        )
        theta_high[parameter_index] = min(
            UB_2D[parameter_index], theta[parameter_index] + step,
        )
        response_low = transient_response_full_2D_v9(theta_low, p)
        response_high = transient_response_full_2D_v9(theta_high, p)
        delta = theta_high[parameter_index] - theta_low[parameter_index]
        jacobian[:, column] .= (
            (response_high .- response_low) ./ delta
        ) .* range_i
    end

    singular_values = svdvals(jacobian)
    relative_cutoff = 1.0e-3
    numerical_rank = count(
        value -> value > singular_values[1] * relative_cutoff,
        singular_values,
    )
    condition_number = singular_values[1] /
                       max(singular_values[end], eps(Float64))
    column_norms = [
        norm(view(jacobian, :, j)) for j in axes(jacobian, 2)
    ]
    correlation = zeros(length(FIT_INDICES_2D), length(FIT_INDICES_2D))
    for i in axes(correlation, 1), j in axes(correlation, 2)
        correlation[i, j] = dot(
            view(jacobian, :, i), view(jacobian, :, j),
        ) / max(column_norms[i] * column_norms[j], eps(Float64))
    end
    maximum_correlation = maximum(
        abs(correlation[i, j])
        for i in axes(correlation, 1), j in axes(correlation, 2)
        if i != j
    )
    return (
        singular_values=singular_values,
        relative_cutoff=relative_cutoff,
        numerical_rank=numerical_rank,
        condition_number=condition_number,
        column_norms=column_norms,
        correlation=correlation,
        maximum_correlation=maximum_correlation,
    )
end

function fitted_mesh_full_2D_v9(p)
    coarse = with_mesh_full_2D_v9(
        p; nr_rec=10, nr_felt=5, nr_case=2, nz=30, nz_rear=20,
    )
    fine = with_mesh_full_2D_v9(
        p; nr_rec=20, nr_felt=10, nr_case=4, nz=60, nz_rear=40,
    )
    cases = [
        run_full_case_2D_v9("E73", candidate; is_cooling=false)
        for candidate in (coarse, p, fine)
    ]
    final_sensors = [case.model[end, :] for case in cases]
    final_dp1 = [
        case.result.dp1_prediction_mbar[end] for case in cases
    ]
    coarse_to_nominal = maximum(abs.(
        final_sensors[2] .- final_sensors[1]
    ))
    nominal_to_fine = maximum(abs.(
        final_sensors[3] .- final_sensors[2]
    ))
    dp_coarse_to_nominal = abs(final_dp1[2] - final_dp1[1])
    dp_nominal_to_fine = abs(final_dp1[3] - final_dp1[2])
    return (
        final_sensors=final_sensors,
        final_dp1=final_dp1,
        coarse_to_nominal_max_K=coarse_to_nominal,
        nominal_to_fine_max_K=nominal_to_fine,
        dp_coarse_to_nominal_mbar=dp_coarse_to_nominal,
        dp_nominal_to_fine_mbar=dp_nominal_to_fine,
        accepted=nominal_to_fine < coarse_to_nominal &&
                 nominal_to_fine < 20.0,
        nominal_case=cases[2],
    )
end

p_fitted = read_fitted_parameters_full_2D_v9()
mesh = fitted_mesh_full_2D_v9(p_fitted)
identifiability = identifiability_full_2D_v9(p_fitted)

case = mesh.nominal_case
times = case.times
index = max(1, length(times) ÷ 2)
operating = OperatingCondition2D(
    irradiance=p_fitted.optics.scale_304 * 304000.0,
    flow_lpm=10.0,
    inlet_temperature=295.0,
    ambient_temperature=295.0,
)
ledger = energy_rate_ledger2D(
    case.result.ode_solution.u[index],
    p_fitted,
    operating,
    times[index],
)

open(joinpath(
    DIAGNOSTIC_OUTPUT_DIR_2D_v9,
    "fitted_mesh_summary_2D_v9.txt",
), "w") do io
    println(io, "case=E73")
    println(io, "coarse_mesh=(10,5,2,30,20)")
    println(io, "nominal_mesh=(14,7,3,45,30)")
    println(io, "fine_mesh=(20,10,4,60,40)")
    println(io, "final_sensor_vectors=$(mesh.final_sensors)")
    println(io, "final_dp1_mbar=$(mesh.final_dp1)")
    println(
        io,
        "coarse_to_nominal_max_K=$(mesh.coarse_to_nominal_max_K)",
    )
    println(
        io,
        "nominal_to_fine_max_K=$(mesh.nominal_to_fine_max_K)",
    )
    println(
        io,
        "dp_coarse_to_nominal_mbar=$(mesh.dp_coarse_to_nominal_mbar)",
    )
    println(
        io,
        "dp_nominal_to_fine_mbar=$(mesh.dp_nominal_to_fine_mbar)",
    )
    println(io, "accepted=$(mesh.accepted)")
end

open(joinpath(
    DIAGNOSTIC_OUTPUT_DIR_2D_v9,
    "identifiability_summary_2D_v9.txt",
), "w") do io
    println(io, "training_cases=$(collect(TRAIN_HEATING_2D_v9))")
    println(io, "fitted_indices=$(FIT_INDICES_2D)")
    println(
        io,
        "singular_values=$(identifiability.singular_values)",
    )
    println(
        io,
        "relative_rank_cutoff=$(identifiability.relative_cutoff)",
    )
    println(io, "numerical_rank=$(identifiability.numerical_rank)")
    println(io, "condition_number=$(identifiability.condition_number)")
    println(
        io,
        "maximum_absolute_column_correlation=$(identifiability.maximum_correlation)",
    )
    println(io, "column_norms=$(identifiability.column_norms)")
end

names = PARAMETER_NAMES_FULL_2D_v9[FIT_INDICES_2D]
open(joinpath(
    DIAGNOSTIC_OUTPUT_DIR_2D_v9,
    "identifiability_correlation_2D_v9.csv",
), "w") do io
    println(io, join(vcat("parameter", collect(names)), ','))
    for (i, name) in enumerate(names)
        println(io, join(
            vcat(name, string.(identifiability.correlation[i, :])), ',',
        ))
    end
end

open(joinpath(
    DIAGNOSTIC_OUTPUT_DIR_2D_v9,
    "fitted_energy_ledger_2D_v9.txt",
), "w") do io
    for name in propertynames(ledger)
        println(io, "$name=$(getproperty(ledger, name))")
    end
end

@testset "2D_v9 fitted full-test diagnostics" begin
    @test mesh.accepted
    @test all(isfinite, identifiability.singular_values)
    @test identifiability.numerical_rank <= length(FIT_INDICES_2D)
    @test isfinite(identifiability.condition_number)
    @test abs(ledger.residual) < 1.0e-7
end

println(
    "2D_v9 fitted diagnostics: mesh accepted=", mesh.accepted,
    ", rank=", identifiability.numerical_rank, "/",
    length(FIT_INDICES_2D),
    ", condition=", identifiability.condition_number,
    ", max_correlation=", identifiability.maximum_correlation,
    ", energy_residual_W=", ledger.residual,
)
