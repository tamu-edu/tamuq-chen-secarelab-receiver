using CSV
using DataFrames
using Plots
using Plots.PlotMeasures

const EXPONENT_INPUT_DIR = joinpath(@__DIR__, "summaries", "2D_v10")
const EXPONENT_PLOT_DIR = joinpath(EXPONENT_INPUT_DIR, "plots")
mkpath(EXPONENT_PLOT_DIR)

summary = CSV.read(
    joinpath(
        EXPONENT_INPUT_DIR, "front_exponent_summary_2D_v10.csv",
    ),
    DataFrame,
)

default(
    fontfamily="Computer Modern",
    linewidth=2.2,
    framestyle=:box,
    gridalpha=0.20,
    legendfontsize=8,
    guidefontsize=10,
    tickfontsize=8,
    titlefontsize=11,
    dpi=180,
)

irradiances = (456, 304, 256)
colors = (:dodgerblue, :seagreen3, :darkgoldenrod2)
observed_slopes = Dict(456 => 17.33862389158247,
                       304 => 10.476819976890445,
                       256 => 6.898088829268151)

slope_panels = Any[]
for strength in sort(unique(summary.equivalent_strength))
    group = summary[summary.equivalent_strength .== strength, :]
    sort!(group, :exponent)
    panel = plot(
        xlabel="Front Reynolds exponent, mf",
        ylabel="Slope [K/(L/min)]",
        title="Equivalent strength = $(strength)",
        legend=:topleft,
    )
    for (irradiance, color) in zip(irradiances, colors)
        column = Symbol("slope_$(irradiance)_K_per_lpm")
        plot!(
            panel, group.exponent, group[!, column];
            color=color,
            marker=:circle,
            label="Model $(irradiance) kW/m2",
        )
        hline!(
            panel, [observed_slopes[irradiance]];
            color=color,
            linestyle=:dash,
            label="Observed $(irradiance)",
        )
    end
    push!(slope_panels, panel)
end
p1 = plot(
    slope_panels...;
    layout=(1, length(slope_panels)),
    size=(1400, 600),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    plot_title="V10 normalized exponent test: axial flow slopes",
)
savefig(
    p1,
    joinpath(
        EXPONENT_PLOT_DIR, "front_exponent_flow_slopes_2D_v10.png",
    ),
)

p2a = plot(
    xlabel="Front Reynolds exponent, mf",
    ylabel="Mean T12-T8 [K]",
    title="Mean axial offset",
    legend=:bottomright,
)
p2b = plot(
    xlabel="Front Reynolds exponent, mf",
    ylabel="Error [K]",
    title="Axial and T3 error",
    legend=:topright,
)
for strength in sort(unique(summary.equivalent_strength))
    group = summary[summary.equivalent_strength .== strength, :]
    sort!(group, :exponent)
    plot!(
        p2a, group.exponent, group.axial_mean_K;
        marker=:circle, label="Ceq=$(strength)",
    )
    plot!(
        p2b, group.exponent, group.axial_rmse_K;
        marker=:circle, label="Axial RMSE, Ceq=$(strength)",
    )
    plot!(
        p2b, group.exponent, group.T3_mae_K;
        marker=:square, linestyle=:dash, label="T3 MAE, Ceq=$(strength)",
    )
end
hline!(
    p2a, [10.8754];
    color=:black, linestyle=:dash, label="Observed mean",
)
p2 = plot(
    p2a, p2b;
    layout=(1, 2),
    size=(1400, 600),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    plot_title="V10 normalized exponent test: profile tradeoff",
)
savefig(
    p2,
    joinpath(
        EXPONENT_PLOT_DIR, "front_exponent_profile_tradeoff_2D_v10.png",
    ),
)

p3a = plot(
    xlabel="Front Reynolds exponent, mf",
    ylabel="Mean radial difference [K]",
    title="Radial offsets",
    legend=:bottomleft,
)
p3b = plot(
    xlabel="Front Reynolds exponent, mf",
    ylabel="Mean front exchange",
    title="Exchange magnitude",
    legend=:topleft,
)
for strength in sort(unique(summary.equivalent_strength))
    group = summary[summary.equivalent_strength .== strength, :]
    sort!(group, :exponent)
    plot!(
        p3a, group.exponent, group.mid_radial_mean_K;
        marker=:circle, label="T12-T9, Ceq=$(strength)",
    )
    plot!(
        p3a, group.exponent, group.deep_radial_mean_K;
        marker=:square, label="T11-T10, Ceq=$(strength)",
    )
    plot!(
        p3b, group.exponent, group.mean_front_h_W_m2K;
        marker=:circle, label="h_front, Ceq=$(strength)",
    )
    plot!(
        p3b, group.exponent, 10 .* group.mean_front_to_gas_W;
        marker=:square, linestyle=:dash,
        label="10 x Q_front [W], Ceq=$(strength)",
    )
end
hline!(p3a, [24.4803]; color=:purple, linestyle=:dash,
       label="Observed mean T12-T9")
hline!(p3a, [35.9729]; color=:teal, linestyle=:dash,
       label="Observed mean T11-T10")
p3 = plot(
    p3a, p3b;
    layout=(1, 2),
    size=(1400, 600),
    left_margin=10mm,
    bottom_margin=10mm,
    top_margin=5mm,
    plot_title="V10 normalized exponent test: radial and exchange response",
)
savefig(
    p3,
    joinpath(
        EXPONENT_PLOT_DIR, "front_exponent_radial_exchange_2D_v10.png",
    ),
)

println("V10 exponent-sensitivity plots written to $EXPONENT_PLOT_DIR")
