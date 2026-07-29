using LinearAlgebra
using Statistics

include(joinpath(@__DIR__, "..", "run_2D_v8.jl"))

const FITTED_PARAMETER_NAMES_2D = [
    "A_Nu", "B_Re", "k_scale_r", "k_scale_z", "sigma_beam_mm",
    "spillage_fraction", "scale_456", "scale_304", "scale_256",
    "beta_rad", "receiver_felt_contact_resistance_m2K_W", "beta_opt",
]
const IDENTIFIABILITY_CASES_2D = ("E67", "E71", "E72", "E76", "E77", "E81")

function read_fitted_theta_2D(path)
    lines = readlines(path)
    length(lines) == 13 || error("Expected header plus 12 fitted-parameter rows in $path")
    return [parse(Float64, split(line, ',')[3]) for line in lines[2:end]]
end

function calibration_response_2D(theta)
    p = unpack_parameters2D(theta)
    values = Float64[]
    for id in IDENTIFIABILITY_CASES_2D
        case_data = run_single_case_2D(
            id, calibration_mesh_parameters_2D(p);
            is_cooling=false,
        )
        for sensor_index in axes(case_data.model, 2)
            scale = max(
                maximum(case_data.experiment[:, sensor_index]) -
                minimum(case_data.experiment[:, sensor_index]),
                1.0,
            )
            append!(values, case_data.model[:, sensor_index] ./ scale)
        end
    end
    return values
end

parameter_path = joinpath(
    @__DIR__, "..", "summaries", "2D_v8",
    "parameters_pilot_apparent_nu_2D_v8.csv",
)
theta = read_fitted_theta_2D(parameter_path)
response0 = calibration_response_2D(theta)
jacobian = zeros(length(response0), length(FIT_INDICES_2D))

for (column, parameter_index) in enumerate(FIT_INDICES_2D)
    range_i = UB_2D[parameter_index] - LB_2D[parameter_index]
    step = 0.01 * range_i
    theta_low = copy(theta)
    theta_high = copy(theta)
    theta_low[parameter_index] = max(LB_2D[parameter_index], theta[parameter_index] - step)
    theta_high[parameter_index] = min(UB_2D[parameter_index], theta[parameter_index] + step)
    response_low = calibration_response_2D(theta_low)
    response_high = calibration_response_2D(theta_high)
    delta = theta_high[parameter_index] - theta_low[parameter_index]
    # Sensitivity of span-normalized transient sensors with respect to one
    # unit of normalized parameter range.
    jacobian[:, column] .= ((response_high .- response_low) ./ delta) .* range_i
end

singular_values = svdvals(jacobian)
relative_cutoff = 1.0e-3
numerical_rank = count(value -> value > singular_values[1] * relative_cutoff, singular_values)
condition_number = singular_values[1] / max(singular_values[end], eps(Float64))
column_norms = [norm(view(jacobian, :, j)) for j in axes(jacobian, 2)]
correlation = zeros(length(FIT_INDICES_2D), length(FIT_INDICES_2D))
for i in axes(correlation, 1), j in axes(correlation, 2)
    correlation[i, j] = dot(view(jacobian, :, i), view(jacobian, :, j)) /
                        max(column_norms[i] * column_norms[j], eps(Float64))
end
off_diagonal = [
    abs(correlation[i, j])
    for i in axes(correlation, 1), j in axes(correlation, 2) if i != j
]
maximum_correlation = maximum(off_diagonal)

output_dir = joinpath(@__DIR__, "..", "summaries", "2D_v8")
open(joinpath(output_dir, "identifiability_pilot_apparent_nu_summary_2D_v8.txt"), "w") do io
    println(io, "cases=$(collect(IDENTIFIABILITY_CASES_2D))")
    println(io, "fitted_indices=$(FIT_INDICES_2D)")
    println(io, "singular_values=$(singular_values)")
    println(io, "relative_rank_cutoff=$(relative_cutoff)")
    println(io, "numerical_rank=$(numerical_rank)")
    println(io, "condition_number=$(condition_number)")
    println(io, "maximum_absolute_column_correlation=$(maximum_correlation)")
    println(io, "column_norms=$(column_norms)")
end

open(joinpath(output_dir, "identifiability_pilot_apparent_nu_correlation_2D_v8.csv"), "w") do io
    fitted_names = FITTED_PARAMETER_NAMES_2D[FIT_INDICES_2D]
    println(io, join(vcat("parameter", fitted_names), ','))
    for (i, name) in enumerate(fitted_names)
        println(io, join(vcat(name, string.(correlation[i, :])), ','))
    end
end

println(
    "2D_v8 local transient-response identifiability: rank=", numerical_rank,
    "/", length(FIT_INDICES_2D),
    " condition_number=", condition_number,
    " maximum_abs_correlation=", maximum_correlation,
)
