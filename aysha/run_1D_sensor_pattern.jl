# ============================================================================
# run_1D_sensor_pattern.jl - experimental thermocouple pattern diagnostic
# ============================================================================
# Summarizes the steady solid thermocouple ordering by irradiance and flow.
# ============================================================================

include("1D_v10.jl")
using StatsPlots

const SENSOR_PATTERN_DIR = joinpath(@__DIR__, "summaries", "1D_sensor_pattern")
const SENSOR_PATTERN_PLOT_DIR = joinpath(SENSOR_PATTERN_DIR, "plots")
mkpath(SENSOR_PATTERN_DIR)
mkpath(SENSOR_PATTERN_PLOT_DIR)

function write_sensor_pattern()
    rows = NamedTuple[]
    for simulation_id in sim_key_heat
        conditions = simulation_conditions[simulation_id]
        T8 = observation(measurements, simulation_id, "_T8")[end]
        T9 = observation(measurements, simulation_id, "_T9")[end]
        T10 = observation(measurements, simulation_id, "_T10")[end]
        T11 = observation(measurements, simulation_id, "_T11")[end]
        T12 = observation(measurements, simulation_id, "_T12")[end]
        T9_pair = 0.5 * (T9 + T12)
        T10_pair = 0.5 * (T10 + T11)
        push!(rows, (
            simulation_id=simulation_id,
            irradiance_kW_m2=conditions[Io] / 1000.0,
            flow_lpm=conditions[qlpm],
            T8_K=T8,
            T9_K=T9,
            T12_K=T12,
            T9_pair_K=T9_pair,
            T10_K=T10,
            T11_K=T11,
            T10_pair_K=T10_pair,
            T9_minus_T8_K=T9 - T8,
            T12_minus_T8_K=T12 - T8,
            T9_pair_minus_T8_K=T9_pair - T8,
            T10_pair_minus_T8_K=T10_pair - T8,
        ))
    end

    path = joinpath(SENSOR_PATTERN_DIR, "steady_thermocouple_pattern.csv")
    open(path, "w") do io
        header = propertynames(first(rows))
        println(io, join(header, ','))
        for row in rows
            println(io, join((getproperty(row, name) for name in header), ','))
        end
    end

    println("[sensor-pattern] Saved: $path")
    println("[sensor-pattern] Key grouped trend:")
    for irradiance in sort(unique(row.irradiance_kW_m2 for row in rows); rev=true)
        group = sort(filter(row -> row.irradiance_kW_m2 == irradiance, rows),
                     by=row -> row.flow_lpm, rev=true)
        println("I=$(round(irradiance; digits=1)) kW/m2")
        for row in group
            println("  $(row.simulation_id): flow=$(round(row.flow_lpm; digits=2)) L/min, " *
                    "T8=$(round(row.T8_K; digits=1)) K, " *
                    "T9pair-T8=$(round(row.T9_pair_minus_T8_K; digits=1)) K, " *
                    "T12-T8=$(round(row.T12_minus_T8_K; digits=1)) K, " *
                    "T10pair-T8=$(round(row.T10_pair_minus_T8_K; digits=1)) K")
        end

        flow = [row.flow_lpm for row in group]
        plot_temperature = plot(
            title="Experimental solid thermocouples, I=$(round(irradiance; digits=1)) kW/m2",
            xlabel="Flow rate (L/min)",
            ylabel="Temperature (K)",
            ylims=TEMPERATURE_AXIS_LIMITS_v10,
            legend=:outerright,
            xflip=true,
        )
        plot!(plot_temperature, flow, [row.T8_K for row in group];
              label="T8", marker=:square, lw=2)
        plot!(plot_temperature, flow, [row.T9_pair_K for row in group];
              label="(T9+T12)/2", marker=:circle, lw=2)
        plot!(plot_temperature, flow, [row.T12_K for row in group];
              label="T12", marker=:utriangle, lw=2)
        plot!(plot_temperature, flow, [row.T10_pair_K for row in group];
              label="(T10+T11)/2", marker=:diamond, lw=2)
        plot!(plot_temperature, flow, [row.T10_K for row in group];
              label="T10", marker=:dtriangle, lw=2, linestyle=:dash)
        plot!(plot_temperature, flow, [row.T11_K for row in group];
              label="T11", marker=:hexagon, lw=2, linestyle=:dash)
        savefig(
            plot_temperature,
            joinpath(SENSOR_PATTERN_PLOT_DIR,
                     "temperature_vs_flow_I$(round(Int, irradiance))_1D.png"),
        )

        plot_difference = plot(
            title="Internal-minus-front difference, I=$(round(irradiance; digits=1)) kW/m2",
            xlabel="Flow rate (L/min)",
            ylabel="Temperature difference (K)",
            legend=:outerright,
            xflip=true,
        )
        plot!(plot_difference, flow, [row.T9_pair_minus_T8_K for row in group];
              label="(T9+T12)/2 - T8", marker=:circle, lw=2)
        plot!(plot_difference, flow, [row.T12_minus_T8_K for row in group];
              label="T12 - T8", marker=:utriangle, lw=2)
        plot!(plot_difference, flow, [row.T10_pair_minus_T8_K for row in group];
              label="(T10+T11)/2 - T8", marker=:diamond, lw=2)
        plot!(plot_difference, flow, [row.T10_K - row.T8_K for row in group];
              label="T10 - T8", marker=:dtriangle, lw=2, linestyle=:dash)
        plot!(plot_difference, flow, [row.T11_K - row.T8_K for row in group];
              label="T11 - T8", marker=:hexagon, lw=2, linestyle=:dash)
        hline!(plot_difference, [0.0]; label=false, color=:gray, linestyle=:dash)
        savefig(
            plot_difference,
            joinpath(SENSOR_PATTERN_PLOT_DIR,
                     "internal_minus_front_I$(round(Int, irradiance))_1D.png"),
        )
    end
    return path
end

write_sensor_pattern()
