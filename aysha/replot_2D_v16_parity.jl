using DelimitedFiles
using Plots

const ROOT_2D_v16_REPLOT = joinpath(@__DIR__, "summaries", "2D_v16")
const FINAL_2D_v16_REPLOT = joinpath(
    ROOT_2D_v16_REPLOT, "staged_final_profiles_2D_v16.csv",
)

function replot_parity_2D_v16()
    raw, header_raw = readdlm(
        FINAL_2D_v16_REPLOT, ',', Any, '\n'; header=true,
    )
    header = vec(String.(header_raw))
    column = Dict(name => index for (index, name) in enumerate(header))
    phase_colors = Dict(
        "heating_training" => :royalblue,
        "heating_validation" => :darkorange,
        "cooling_validation" => :seagreen,
    )
    phase_markers = Dict(
        "heating_training" => :circle,
        "heating_validation" => :diamond,
        "cooling_validation" => :utriangle,
    )
    parity = plot(
        xlabel="Observed final temperature (K)",
        ylabel="Modeled final temperature (K)",
        title="2D v16 cooling + power-refitted parity",
        legend=:topleft, size=(850, 720),
    )
    all_values = Float64[]
    sensors = ("T8", "T12", "T11", "T9", "T10", "T3", "T2")
    for phase in keys(phase_colors)
        xs = Float64[]
        ys = Float64[]
        for row in axes(raw, 1)
            String(raw[row, column["phase"]]) == phase || continue
            for sensor in sensors
                push!(
                    xs,
                    Float64(raw[
                        row, column["observed_$(sensor)_K"],
                    ]),
                )
                push!(
                    ys,
                    Float64(raw[
                        row, column["model_$(sensor)_K"],
                    ]),
                )
            end
        end
        append!(all_values, xs)
        append!(all_values, ys)
        scatter!(
            parity, xs, ys;
            label=replace(phase, "_" => " "),
            color=phase_colors[phase],
            marker=phase_markers[phase], alpha=0.8,
        )
    end
    lo, hi = extrema(all_values)
    plot!(
        parity, [lo, hi], [lo, hi];
        label="1:1", color=:black, linestyle=:dash,
    )
    savefig(
        parity,
        joinpath(
            ROOT_2D_v16_REPLOT, "plots",
            "parity_staged_2D_v16.png",
        ),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    replot_parity_2D_v16()
end
