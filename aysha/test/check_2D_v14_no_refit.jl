using Statistics

include(joinpath(@__DIR__, "..", "run_2D_v14.jl"))

mkpath(OUTPUT_DIR_2D_v14)
p = inherited_parameters_2D_v14(
    nominal_mesh=false, screen_mesh=true,
)
cases = Vector{Any}(undef, length(TRAIN_HEATING_2D_v14))
Threads.@threads for index in eachindex(TRAIN_HEATING_2D_v14)
    id = TRAIN_HEATING_2D_v14[index]
    inputs = case_inputs_2D_v14(id; max_points=61)
    cases[index] = simulate_case_2D_v14(
        inputs, p; reltol=3e-4, abstol=1e-4, dtmax=120.0,
    )
end
rows = NamedTuple[]
for case in cases
    m = case.model[end, :]
    e = case.observed[end, :]
    result = case.result
    k = length(result.time)
    core_flow = result.channel_mass_flow_kg_s[
        result.core_group, k,
    ]
    side_flow = result.channel_mass_flow_kg_s[
        result.side_group, k,
    ]
    corner_flow = result.channel_mass_flow_kg_s[
        result.corner_group, k,
    ]
    push!(rows, (
        simulation_id=case.inputs.id,
        rmse_K=sqrt(mean(abs2, case.model .- case.observed)),
        model_mid_gap_K=m[2] - m[4],
        observed_mid_gap_K=e[2] - e[4],
        model_deep_gap_K=m[3] - m[5],
        observed_deep_gap_K=e[3] - e[5],
        model_axial_gap_K=m[2] - m[1],
        observed_axial_gap_K=e[2] - e[1],
        core_to_side_flow_ratio=core_flow / side_flow,
        core_to_corner_flow_ratio=core_flow / corner_flow,
        equal_pressure_error=
            result.equal_pressure_relative_error[end],
        dp1_model_mbar=result.dp1_prediction_mbar[end],
        dp1_observed_mbar=case.inputs.dp1[end],
    ))
end
_write_namedtuples_csv_2D_v14(
    joinpath(
        OUTPUT_DIR_2D_v14,
        "no_refit_training_diagnostics_2D_v14.csv",
    ), rows,
)
println(
    "no-refit mean RMSE = ",
    mean(row.rmse_K for row in rows),
)
println(
    "mid radial sign = ",
    count(row -> sign(row.model_mid_gap_K) ==
                 sign(row.observed_mid_gap_K), rows),
    "/9",
)
println(
    "deep radial sign = ",
    count(row -> sign(row.model_deep_gap_K) ==
                 sign(row.observed_deep_gap_K), rows),
    "/9",
)
println(
    "mean core/side flow ratio = ",
    mean(row.core_to_side_flow_ratio for row in rows),
)
println(
    "mean core/corner flow ratio = ",
    mean(row.core_to_corner_flow_ratio for row in rows),
)
foreach(println, rows)
