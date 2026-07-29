include(joinpath(@__DIR__, "..", "run_2D_v8.jl"))

function summarize_phase(rows, phase)
    selected = filter(row -> row.phase == phase, rows)
    return (
        rmse_mean=mean(row.rmse_K for row in selected),
        rmse_median=median(row.rmse_K for row in selected),
        steady_mae=mean(abs(row.steady_error_K) for row in selected),
        t90_mae=mean(abs(row.t90_error_s) for row in selected),
    )
end

p = default_parameters2D()
for (id, cooling) in (("E81", false), ("C81", true))
    data = run_single_case_2D(id, p; is_cooling=cooling)
    error = data.model - data.experiment
    println(
        id,
        " rmse=", round(sqrt(mean(abs2, error)), digits=3),
        " final_mae=", round(mean(abs.(error[end, :])), digits=3),
        " final_model=", round.(data.model[end, :], digits=1),
        " final_experiment=", round.(data.experiment[end, :], digits=1),
    )
end

if "--all" in ARGS
    metrics = compute_all_metrics_2D(p)
    println("heating=", summarize_phase(metrics, :heating_calibrated))
    println("cooling=", summarize_phase(metrics, :cooling_benchmark))
end
