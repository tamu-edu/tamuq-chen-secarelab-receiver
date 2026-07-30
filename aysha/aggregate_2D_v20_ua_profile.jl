using CSV
using DataFrames

const directory_2D_v20 = joinpath(
    @__DIR__, "summaries", "2D_v20",
)
const parts_2D_v20 = (
    "ua_plant_profile_corrected_part12_2D_v20.csv",
    "ua_plant_profile_corrected_part345_2D_v20.csv",
)

function aggregate_2D_v20_ua_profile()
    frames = [
        CSV.read(joinpath(directory_2D_v20, part), DataFrame)
        for part in parts_2D_v20
    ]
    combined = vcat(frames...)
    sort!(combined, :objective)
    output = joinpath(
        directory_2D_v20,
        "ua_plant_profile_2D_v20.csv",
    )
    CSV.write(output, combined)
    open(joinpath(
        directory_2D_v20,
        "ua_plant_profile_summary_2D_v20.txt",
    ), "w") do io
        best = combined[1, :]
        println(io, "candidate_count=", nrow(combined))
        println(io, "best_reynolds_exponent=",
                best.reynolds_exponent)
        println(io, "best_Nu50_ratio=", best.Nu50_ratio)
        println(io, "best_nu_prefactor=", best.nu_prefactor)
        println(io, "best_objective=", best.objective)
        println(io, "best_side_rmse_K=", best.side_rmse_K)
        println(io, "best_T2_rmse_K=", best.T2_rmse_K)
        println(io, "best_T3_diagnostic_rmse_K=",
                best.T3_diagnostic_rmse_K)
        println(io, "best_exponent_on_boundary=",
                best.reynolds_exponent in (1.25, 1.65))
        println(io, "best_Nu50_ratio_on_boundary=",
                best.Nu50_ratio in (0.75, 1.95))
        near = combined[
            combined.objective .<= 1.10best.objective, :,
        ]
        println(io, "within_10pct_count=", nrow(near))
        println(io, "within_10pct_n_min=",
                minimum(near.reynolds_exponent))
        println(io, "within_10pct_n_max=",
                maximum(near.reynolds_exponent))
        println(io, "within_10pct_Nu50_ratio_min=",
                minimum(near.Nu50_ratio))
        println(io, "within_10pct_Nu50_ratio_max=",
                maximum(near.Nu50_ratio))
        println(io, "ua_identified=false")
    end
    return combined
end

if abspath(PROGRAM_FILE) == @__FILE__
    aggregate_2D_v20_ua_profile()
end
