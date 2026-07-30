using CSV
using DataFrames
using Plots
using Statistics

const INPUT_DIR = joinpath(@__DIR__, "summaries", "2D_v10")
const PLOT_DIR = joinpath(INPUT_DIR, "plots")
mkpath(PLOT_DIR)

cases = CSV.read(
    joinpath(INPUT_DIR, "front_sensitivity_cases_2D_v10.csv"), DataFrame,
)
slopes = CSV.read(
    joinpath(INPUT_DIR, "front_sensitivity_slopes_2D_v10.csv"), DataFrame,
)

coefficients = sort(unique(cases.front_coefficient))
xvalue(c) = log10(1.0 + c)
xvalues = xvalue.(coefficients)
tick_coefficients = [0.0, 0.1, 0.5, 1.0, 2.0, 4.0, 8.0]
tick_values = xvalue.(tick_coefficients)
tick_labels = string.(tick_coefficients)

function grouped_mean(column)
    return [
        mean(cases[cases.front_coefficient .== coefficient, column])
        for coefficient in coefficients
    ]
end

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

p1 = plot(
    xlabel="Front coefficient, Cf",
    ylabel="Final temperature difference [K]",
    title="V10 front exchange: mean spatial-profile response",
    xticks=(tick_values, tick_labels),
    legend=:bottomright,
)
plot!(
    p1, xvalues, grouped_mean(:model_T12_minus_T8_K);
    marker=:circle, label="Model T12-T8",
)
hline!(
    p1, [mean(cases.observed_T12_minus_T8_K)];
    linestyle=:dash, label="Observed mean T12-T8",
)
plot!(
    p1, xvalues, grouped_mean(:model_T12_minus_T9_K);
    marker=:square, label="Model T12-T9",
)
hline!(
    p1, [mean(cases.observed_T12_minus_T9_K)];
    linestyle=:dash, label="Observed mean T12-T9",
)
plot!(
    p1, xvalues, grouped_mean(:model_T11_minus_T10_K);
    marker=:diamond, label="Model T11-T10",
)
hline!(
    p1, [mean(cases.observed_T11_minus_T10_K)];
    linestyle=:dash, label="Observed mean T11-T10",
)
savefig(p1, joinpath(PLOT_DIR, "front_sensitivity_profile_offsets_2D_v10.png"))

p2 = plot(
    xlabel="Front coefficient, Cf",
    ylabel="Slope of T12-T8 vs flow [K/(L/min)]",
    title="V10 front exchange: flow sensitivity by irradiance",
    xticks=(tick_values, tick_labels),
    legend=:topright,
)
for irradiance in sort(unique(slopes.irradiance_kW_m2); rev=true)
    group = slopes[slopes.irradiance_kW_m2 .== irradiance, :]
    sort!(group, :front_coefficient)
    plot!(
        p2,
        xvalue.(group.front_coefficient),
        group.model_slope_K_per_lpm;
        marker=:circle,
        label="Model $(Int(round(irradiance))) kW/m2",
    )
    hline!(
        p2, [first(group.observed_slope_K_per_lpm)];
        linestyle=:dash,
        label="Observed $(Int(round(irradiance))) kW/m2",
    )
end
savefig(p2, joinpath(PLOT_DIR, "front_sensitivity_flow_slopes_2D_v10.png"))

p3 = plot(
    xlabel="Front coefficient, Cf",
    ylabel="Front exchange",
    title="V10 conservative front-exchange magnitude",
    xticks=(tick_values, tick_labels),
    legend=:topleft,
)
plot!(
    p3, xvalues, grouped_mean(:front_to_gas_W);
    marker=:circle, label="Heat to inlet gas [W]",
)
plot!(
    p3, xvalues, grouped_mean(:mean_gas_preheat_K);
    marker=:square, label="Mean gas preheat [K]",
)
plot!(
    p3, xvalues, grouped_mean(:mean_front_h_W_m2K) ./ 10.0;
    marker=:diamond, label="Mean h_front / 10 [W/m2/K]",
)
savefig(p3, joinpath(PLOT_DIR, "front_sensitivity_exchange_2D_v10.png"))

println("V10 sensitivity plots written to $PLOT_DIR")
