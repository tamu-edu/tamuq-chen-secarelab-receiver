include(joinpath(@__DIR__, "..", "run_2D_v8.jl"))

const PILOT_CASES_2D = ("E67", "E71", "E72", "E76", "E77", "E81")

function run_pilot_case_2D(simulation_id, p=default_parameters2D(); is_cooling=false)
    return run_single_case_2D(
        simulation_id, calibration_mesh_parameters_2D(p);
        is_cooling=is_cooling,
    )
end

function phase_summary_2D(rows, phase)
    selected = filter(row -> row.phase == phase, rows)
    return (
        rmse_mean=mean(row.rmse_K for row in selected),
        steady_mae=mean(abs(row.steady_error_K) for row in selected),
        t90_mae=mean(abs(row.t90_error_s) for row in selected),
    )
end

calibration = calibrate2D(
    run_pilot_case_2D;
    heating_cases=PILOT_CASES_2D,
    max_iters=60,
    max_time=1200.0,
    p_init=default_parameters2D(),
)
p_fitted = calibration.parameters
println("pilot_objective=$(calibration.objective)")
println("pilot_return_code=$(calibration.retcode)")
println("pilot_evaluations=$(calibration.evaluations)")
println("pilot_parameters=$(pack_parameters2D(p_fitted))")
flush(stdout)

metrics = compute_all_metrics_2D(p_fitted)
heating = phase_summary_2D(metrics, :heating_calibrated)
cooling = phase_summary_2D(metrics, :cooling_benchmark)
steady = build_steady_results_2D(p_fitted)
mid_sign_correct = count(
    row -> row.T12_minus_T9_model * row.T12_minus_T9_experiment > 0.0,
    steady,
)
deep_sign_correct = count(
    row -> row.T11_minus_T10_model * row.T11_minus_T10_experiment > 0.0,
    steady,
)

output_dir = joinpath(@__DIR__, "..", "summaries", "2D_v8")
write_parameters_2D(
    joinpath(output_dir, "parameters_pilot_apparent_nu_2D_v8.csv"),
    p_fitted,
)
open(joinpath(output_dir, "pilot_apparent_nu_summary_2D_v8.txt"), "w") do io
    println(io, "objective=$(calibration.objective)")
    println(io, "return_code=$(calibration.retcode)")
    println(io, "evaluations=$(calibration.evaluations)")
    println(io, "calibration_cases=$(collect(PILOT_CASES_2D))")
    println(io, "calibration_mesh=(10,5,2,30,20)")
    println(io, "parameters=$(pack_parameters2D(p_fitted))")
    println(io, "heating=$(heating)")
    println(io, "cooling=$(cooling)")
    println(io, "mid_sign_correct=$(mid_sign_correct)/$(length(steady))")
    println(io, "deep_sign_correct=$(deep_sign_correct)/$(length(steady))")
end
println("pilot_heating=$(heating)")
println("pilot_cooling=$(cooling)")
println("pilot_mid_sign_correct=$(mid_sign_correct)/$(length(steady))")
println("pilot_deep_sign_correct=$(deep_sign_correct)/$(length(steady))")
