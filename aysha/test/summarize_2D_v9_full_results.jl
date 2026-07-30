using CSV
using DataFrames
using Statistics

const RESULT_DIR_2D_v9 = joinpath(
    @__DIR__, "..", "summaries", "2D_v9",
)

function slope_summary_2D_v9(x, y)
    xmean = mean(x)
    ymean = mean(y)
    denominator = sum((value - xmean)^2 for value in x)
    return sum(
        (x[i] - xmean) * (y[i] - ymean) for i in eachindex(x)
    ) / denominator
end

predictions = CSV.read(
    joinpath(RESULT_DIR_2D_v9, "transient_predictions_2D_v9.csv"),
    DataFrame,
)
heating = filter(row -> row.phase != "cooling_validation", predictions)
final_indices = combine(
    groupby(heating, :simulation_id),
    :time_s => argmax => :row_index,
)

final_rows = DataFrame()
for id in final_indices.simulation_id
    rows = filter(row -> row.simulation_id == id, heating)
    push!(final_rows, rows[argmax(rows.time_s), :])
end

function irradiance_group_2D_v9(id)
    number = parse(Int, id[2:end])
    number <= 71 && return 456000.0
    number <= 76 && return 304000.0
    return 256000.0
end
final_rows.irradiance_W_m2 = irradiance_group_2D_v9.(
    final_rows.simulation_id,
)

sensors = (
    ("T8", :T8_model_K, :T8_experiment_K),
    ("T12_perim", :T12_model_K, :T12_experiment_K),
    ("T11_perim", :T11_model_K, :T11_experiment_K),
    ("T9_core", :T9_model_K, :T9_experiment_K),
    ("T10_core", :T10_model_K, :T10_experiment_K),
    ("T3", :T3_model_K, :T3_experiment_K),
    ("T2", :T2_model_K, :T2_experiment_K),
)

rows = NamedTuple[]
for irradiance in (456000.0, 304000.0, 256000.0)
    group = filter(
        row -> row.irradiance_W_m2 == irradiance, final_rows,
    )
    flows = Float64.(group.flow_lpm)
    for (sensor, model_column, experiment_column) in sensors
        model_slope = slope_summary_2D_v9(
            flows, Float64.(group[!, model_column]),
        )
        experiment_slope = slope_summary_2D_v9(
            flows, Float64.(group[!, experiment_column]),
        )
        push!(rows, (
            irradiance_W_m2=irradiance,
            sensor=sensor,
            model_slope_K_per_Lmin=model_slope,
            experiment_slope_K_per_Lmin=experiment_slope,
            sign_correct=model_slope * experiment_slope > 0.0,
        ))
    end
end

CSV.write(joinpath(
    RESULT_DIR_2D_v9, "flow_slopes_2D_v9.csv",
), DataFrame(rows))

sign_correct = count(row -> row.sign_correct, rows)
open(joinpath(
    RESULT_DIR_2D_v9, "flow_slope_summary_2D_v9.txt",
), "w") do io
    println(io, "sign_correct=$sign_correct/$(length(rows))")
    println(
        io,
        "incorrect=$(collect((row.irradiance_W_m2, row.sensor) for row in rows if !row.sign_correct))",
    )
end

println("2D_v9 fitted flow-slope signs correct=$sign_correct/$(length(rows))")
