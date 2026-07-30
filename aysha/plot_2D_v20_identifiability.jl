using CSV
using DataFrames
using Plots

const PLOT_DIR_2D_v20 = joinpath(
    @__DIR__, "summaries", "2D_v20", "plots",
)

function plot_2D_v20_identifiability()
    mkpath(PLOT_DIR_2D_v20)
    default(
        fontfamily="Computer Modern",
        linewidth=2.2,
        markersize=5,
        legendfontsize=8,
        guidefontsize=10,
        tickfontsize=8,
        titlefontsize=11,
        framestyle=:box,
        gridalpha=0.20,
        dpi=180,
    )

    ua = CSV.read(joinpath(
        @__DIR__, "summaries", "2D_v20",
        "ua_plant_profile_2D_v20.csv",
    ), DataFrame)
    p_objective = plot(
        xlabel="Nu50 / measured Nu50",
        ylabel="T3-free plant objective",
        title="V20 integrated-UA identifiability profile",
        xlims=(0.68, 2.02),
    )
    palette = cgrad(:viridis, length(unique(ua.reynolds_exponent));
                    categorical=true)
    for (index, exponent) in enumerate(sort(unique(
        ua.reynolds_exponent,
    )))
        subset = sort(
            ua[ua.reynolds_exponent .== exponent, :],
            :Nu50_ratio,
        )
        plot!(
            p_objective, subset.Nu50_ratio, subset.objective,
            marker=:circle, color=palette[index],
            label="n=$(round(exponent; digits=2))",
        )
    end
    hline!(
        p_objective, [1.10minimum(ua.objective)],
        color=:black, linestyle=:dash, linewidth=1.3,
        label="10% profile band",
    )
    savefig(
        p_objective,
        joinpath(PLOT_DIR_2D_v20, "ua_plant_profile_2D_v20.png"),
    )

    observer = CSV.read(joinpath(
        @__DIR__, "summaries", "2D_v20",
        "t3_observer_profile_2D_v20.csv",
    ), DataFrame)
    locations = ["receiver_136", "receiver_exit", "rear_003"]
    labels = ["136 mm", "137 mm", "140 mm"]
    train = Float64[]
    valid = Float64[]
    for location in locations
        subset = observer[observer.location .== location, :]
        winner = subset[argmin(subset.training_objective), :]
        push!(train, winner.training_T3_rmse_K)
        push!(valid, winner.validation_T3_rmse_K)
    end
    xlocations = collect(1:length(labels))
    p_observer = bar(
        xlocations .- 0.18, train,
        label="C69/C80 training",
        ylabel="T3 RMSE (K)",
        title="V20 discrete T3 observer profile minima",
        color=:steelblue,
        bar_width=0.34,
        xticks=(xlocations, labels),
        ylims=(0.0, 1.12maximum(vcat(train, valid))),
        legend=:topright,
    )
    bar!(
        p_observer, xlocations .+ 0.18, valid,
        label="C81 validation", color=:darkorange,
        bar_width=0.34,
    )
    hline!(
        p_observer, [30.0],
        color=:black, linestyle=:dash, linewidth=1.3,
        label="30 K gate",
    )
    savefig(
        p_observer,
        joinpath(PLOT_DIR_2D_v20, "t3_location_profile_2D_v20.png"),
    )
    return (
        ua_plot=p_objective,
        observer_plot=p_observer,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_2D_v20_identifiability()
end
