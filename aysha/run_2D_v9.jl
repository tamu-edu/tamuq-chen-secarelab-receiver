# ============================================================================
# run_2D_v9.jl - Pressure-informed 2D_v9 hydraulic audit and model runner
# ============================================================================

using Statistics

begin
    const Io = :Io
    const qlpm = :qlpm
    const inlet_temperature = :inlet_temperature
    const ambient_temperature = :ambient_temperature
    const time_data = :time_data
    const Tinit = :Tinit
    const T_in = :T_in
    const T_amb = :T_amb
end

include("2D_v9.jl")
include("import_exp_1D_v2.jl")

using .Receiver2D_v9

const OUTPUT_DIR_2D_v9 = joinpath(@__DIR__, "summaries", "2D_v9")
const COLD_T0_WINDOW_2D_v9 = 10
const COLD_T0_TOLERANCE_K_2D_v9 = 5.0

function mean_first(values, count=COLD_T0_WINDOW_2D_v9)
    n = min(length(values), count)
    return mean(view(values, 1:n))
end

function linear_fit_2D_v9(x, y)
    length(x) == length(y) || error("linear-fit vectors must have equal length")
    length(x) >= 3 || error("at least three points are required")
    xmean = mean(x)
    ymean = mean(y)
    sxx = sum((xi - xmean)^2 for xi in x)
    sxx > eps(Float64) || error("linear-fit flow span is zero")
    slope = sum((x[i] - xmean) * (y[i] - ymean) for i in eachindex(x)) / sxx
    intercept = ymean - slope * xmean
    residuals = y .- (intercept .+ slope .* x)
    syy = sum((yi - ymean)^2 for yi in y)
    r2 = 1.0 - sum(abs2, residuals) / syy
    rmse = sqrt(mean(abs2, residuals))
    return (intercept=intercept, slope=slope, r2=r2, rmse=rmse)
end

function cold_t0_dp1_calibration_2D_v9()
    rows = NamedTuple[]
    for meta in filter(row -> row.phase == :heating, hydraulic_t0_metadata)
        is_cold = abs(meta.T3_K - meta.Tamb_K) <= COLD_T0_TOLERANCE_K_2D_v9 &&
                  abs(meta.Tsolid_mean_K - meta.Tamb_K) <= COLD_T0_TOLERANCE_K_2D_v9
        push!(rows, (
            simulation_id=meta.simulation_id,
            flow_lpm=meta.flow_lpm,
            dp1_mbar=meta.dp1_mbar,
            T3_K=meta.T3_K,
            Tsolid_mean_K=meta.Tsolid_mean_K,
            Tamb_K=meta.Tamb_K,
            selected=is_cold,
        ))
    end

    selected = filter(row -> row.selected, rows)
    fit = linear_fit_2D_v9(
        [row.flow_lpm for row in selected],
        [row.dp1_mbar for row in selected],
    )
    p_ref = default_parameters2D()
    ideal_slope = ideal_square_channel_dp_slope2D(p_ref)
    return (
        rows=rows,
        selected=selected,
        intercept_mbar=fit.intercept,
        slope_mbar_per_lpm=fit.slope,
        r2=fit.r2,
        rmse_mbar=fit.rmse,
        ideal_slope_mbar_per_lpm=ideal_slope,
        hydraulic_resistance_scale=fit.slope / ideal_slope,
    )
end

function with_t0_hydraulics_2D_v9(p::ModelParameters2D, calibration)
    h0 = p.hydraulics
    hydraulics = HydraulicParameters2D(
        standard_pressure=h0.standard_pressure,
        standard_temperature=h0.standard_temperature,
        atmospheric_pressure=h0.atmospheric_pressure,
        mass_flow_scale=1.0,
        dp1_zero_offset_mbar=calibration.intercept_mbar,
        hydraulic_resistance_scale=calibration.hydraulic_resistance_scale,
        minor_loss_coefficient=0.0,
    )
    return ModelParameters2D(
        p.geometry, p.solid, p.heat_transfer, p.losses, p.optics, hydraulics,
    )
end

function experimental_initial_state_2D_v9(data, id, grid, p)
    sensors = Dict(
        :T8 => Float64(observation(data, id, "_T8")[1]),
        :T12 => Float64(observation(data, id, "_T12")[1]),
        :T11 => Float64(observation(data, id, "_T11")[1]),
        :T9 => Float64(observation(data, id, "_T9")[1]),
        :T10 => Float64(observation(data, id, "_T10")[1]),
        :T3 => Float64(observation(data, id, "_T3")[1]),
        :T2 => Float64(observation(data, id, "_T2")[1]),
    )
    ambient = Float64(observation(data, id, "_Tamb")[1])
    return build_initial_state_2D(grid, p, sensors, ambient)
end

function run_hydraulic_case_2D_v9(id, p; cooling=false)
    data = cooling ? measurements_cooling : measurements
    times = Float64.(observation_time(data, id))
    flow = Float64.(observation(data, id, "_flow"))
    Tin = Float64.(observation(data, id, "_Tin"))
    Tamb = Float64.(observation(data, id, "_Tamb"))
    nominal_irradiance = cooling ? 0.0 :
                         Float64(simulation_conditions[id][Io])
    scale_factor = nominal_irradiance >= 400000.0 ? p.optics.scale_456 :
                   (nominal_irradiance >= 280000.0 ? p.optics.scale_304 :
                    p.optics.scale_256)
    irradiance = cooling ? zeros(length(times)) :
                 fill(scale_factor * nominal_irradiance, length(times))
    op = OperatingCondition2D(
        irradiance=linear_history(times, irradiance),
        flow_lpm=linear_history(times, flow),
        inlet_temperature=linear_history(times, Tin),
        ambient_temperature=linear_history(times, Tamb),
    )
    grid = Receiver2D_v9.build_grid2D(p)
    initial = cooling ?
        experimental_initial_state_2D_v9(data, id, grid, p) :
        Float64(Tamb[1])
    result = simulate2D(p, op, times; initial_temperature=initial)
    observed_dp1 = Float64.(observation(data, id, "_DP1"))
    return (times=times, flow=flow, observed_dp1=observed_dp1, result=result)
end

function write_t0_calibration_2D_v9(calibration)
    mkpath(OUTPUT_DIR_2D_v9)
    path = joinpath(OUTPUT_DIR_2D_v9, "dp1_cold_t0_calibration_2D_v9.csv")
    open(path, "w") do io
        println(io, "simulation_id,flow_lpm,dp1_mbar,T3_K,Tsolid_mean_K,Tamb_K,selected")
        for row in calibration.rows
            println(io, join(values(row), ','))
        end
        println(io, "# intercept_mbar,$(calibration.intercept_mbar)")
        println(io, "# slope_mbar_per_lpm,$(calibration.slope_mbar_per_lpm)")
        println(io, "# r2,$(calibration.r2)")
        println(io, "# rmse_mbar,$(calibration.rmse_mbar)")
        println(io, "# ideal_slope_mbar_per_lpm,$(calibration.ideal_slope_mbar_per_lpm)")
        println(io, "# hydraulic_resistance_scale,$(calibration.hydraulic_resistance_scale)")
    end
    return path
end

function write_dp1_summary_2D_v9(p)
    mkpath(OUTPUT_DIR_2D_v9)
    path = joinpath(OUTPUT_DIR_2D_v9, "dp1_summary_2D_v9.csv")
    open(path, "w") do io
        println(io, "simulation_id,phase,flow_t0_lpm,dp1_t0_observed_mbar,dp1_t0_model_mbar,dp1_final_observed_mbar,dp1_final_model_mbar,velocity_t0_mean_m_s,velocity_final_mean_m_s")
        for (ids, cooling) in ((IDs, false), (IDs_cooling, true))
            for id in ids
                case = run_hydraulic_case_2D_v9(id, p; cooling=cooling)
                result = case.result
                println(io, join((
                    id,
                    cooling ? "cooling" : "heating",
                    case.flow[1],
                    case.observed_dp1[1],
                    result.dp1_prediction_mbar[1],
                    case.observed_dp1[end],
                    result.dp1_prediction_mbar[end],
                    mean(result.gas_velocity[:, :, 1]),
                    mean(result.gas_velocity[:, :, end]),
                ), ','))
            end
        end
    end
    return path
end

function main_2D_v9()
    calibration = cold_t0_dp1_calibration_2D_v9()
    p = with_t0_hydraulics_2D_v9(default_parameters2D(), calibration)
    calibration_path = write_t0_calibration_2D_v9(calibration)
    summary_path = write_dp1_summary_2D_v9(p)
    println("2D_v9 pressure-informed hydraulic audit complete.")
    println("  selected cold t0 cases = $(length(calibration.selected))")
    println("  DP1 offset = $(calibration.intercept_mbar) mbar")
    println("  cold slope = $(calibration.slope_mbar_per_lpm) mbar/(standard L/min)")
    println("  hydraulic scale = $(calibration.hydraulic_resistance_scale)")
    println("  calibration artifact = $calibration_path")
    println("  DP1 summary artifact = $summary_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_2D_v9()
end
