# ============================================================================
# run_2D_v7.jl - 2D_v7 Heating-Only Runner with Cooling Validation Benchmarks
# ============================================================================

using Statistics
using StatsPlots

begin # Dictionary Key Symbols for Importer
    Io = :Io
    qlpm = :qlpm
    inlet_temperature = :inlet_temperature
    ambient_temperature = :ambient_temperature
    time_data = :time_data
    Tinit = :Tinit
    T_in = :T_in
    T_amb = :T_amb
end

include("2D_v7.jl")
include("import_exp_1D_v2.jl")

using .Receiver2D_v7

begin # Directory setup
    const OUTPUT_DIR_2D_v7 = joinpath(@__DIR__, "summaries", "2D_v7")
    const PLOT_DIR_2D_v7 = joinpath(OUTPUT_DIR_2D_v7, "plots")
    const TRANSIENT_DIR_2D_v7 = joinpath(PLOT_DIR_2D_v7, "transients")
    const PROFILES_DIR_2D_v7 = joinpath(PLOT_DIR_2D_v7, "2d_profiles")
    const AXIAL_DIR_2D_v7 = joinpath(PLOT_DIR_2D_v7, "axial_profiles")
    const DIAGNOSTICS_DIR_2D_v7 = joinpath(OUTPUT_DIR_2D_v7, "diagnostics")
    
    mkpath(OUTPUT_DIR_2D_v7)
    mkpath(PLOT_DIR_2D_v7)
    mkpath(TRANSIENT_DIR_2D_v7)
    mkpath(PROFILES_DIR_2D_v7)
    mkpath(AXIAL_DIR_2D_v7)
    mkpath(DIAGNOSTICS_DIR_2D_v7)

    const TEMPERATURE_AXIS_LIMITS_2D = (280.0, 1400.0)
    const SENSOR_NAMES_2D = (:T8, :T12_perim, :T11_perim, :T9_core, :T10_core, :T3, :T2)
    const SENSOR_COLORS_2D = Dict(
        :T8 => :blue, :T12_perim => :red, :T11_perim => :darkgreen,
        :T9_core => :orange, :T10_core => :brown, :T3 => :purple, :T2 => :black
    )
    const SENSOR_STYLES_2D = Dict(
        :T8 => :solid, :T12_perim => :dash, :T11_perim => :dot,
        :T9_core => :solid, :T10_core => :solid, :T3 => :solid, :T2 => :dashdot
    )
end

function get_experimental_t0_dict(sim_id; is_cooling=false)
    data = is_cooling ? measurements_cooling : measurements
    Dict(
        :T8 => Float64(observation(data, sim_id, "_T8")[1]),
        :T12 => Float64(observation(data, sim_id, "_T12")[1]),
        :T11 => Float64(observation(data, sim_id, "_T11")[1]),
        :T9 => Float64(observation(data, sim_id, "_T9")[1]),
        :T10 => Float64(observation(data, sim_id, "_T10")[1]),
        :T3 => Float64(observation(data, sim_id, "_T3")[1]),
        :T2 => Float64(observation(data, sim_id, "_T2")[1]),
    )
end

function run_single_case_2D(simulation_id, p=default_parameters2D(); is_cooling=false)
    data = is_cooling ? measurements_cooling : measurements
    conds = is_cooling ? simulation_conditions_cooling : simulation_conditions
    
    conditions = conds[simulation_id]
    save_times = Float64.(observation_time(data, simulation_id))
    
    flow = Float64.(observation(data, simulation_id, "_flow"))
    Tin = Float64.(observation(data, simulation_id, "_Tin"))
    Tamb = Float64.(observation(data, simulation_id, "_Tamb"))
    
    irr_val = Float64(conditions[Io])
    scale_factor = irr_val >= 400000.0 ? p.optics.scale_456 : (irr_val >= 280000.0 ? p.optics.scale_304 : p.optics.scale_256)
    irradiance = is_cooling ? zeros(length(save_times)) : fill(scale_factor * irr_val, length(save_times))
    
    op = OperatingCondition2D(
        irradiance = linear_history(save_times, irradiance),
        flow_lpm = linear_history(save_times, flow),
        inlet_temperature = linear_history(save_times, Tin),
        ambient_temperature = linear_history(save_times, Tamb),
    )
    
    grid = Receiver2D_v7.build_grid2D(p)

    if is_cooling
        t0_sensors = get_experimental_t0_dict(simulation_id; is_cooling=true)
        u0_exp = build_initial_state_2D(grid, p, t0_sensors, Tamb[1])
        result = simulate2D(p, op, save_times; initial_temperature=u0_exp)
    else
        result = simulate2D(p, op, save_times)
    end

    preds = sensor_predictions2D(result)
    
    obs_T8 = Float64.(observation(data, simulation_id, "_T8"))
    obs_T12 = Float64.(observation(data, simulation_id, "_T12"))
    obs_T11 = Float64.(observation(data, simulation_id, "_T11"))
    obs_T9 = Float64.(observation(data, simulation_id, "_T9"))
    obs_T10 = Float64.(observation(data, simulation_id, "_T10"))
    obs_T3 = Float64.(observation(data, simulation_id, "_T3"))
    obs_T2 = Float64.(observation(data, simulation_id, "_T2"))
    
    experiment = hcat(obs_T8, obs_T12, obs_T11, obs_T9, obs_T10, obs_T3, obs_T2)
    model = hcat(preds.T8, preds.T12, preds.T11, preds.T9, preds.T10, preds.T3, preds.T2)
    
    return (result=result, model=model, experiment=experiment, times=save_times)
end

function compute_all_metrics_2D(p=default_parameters2D(); heating_keys=IDs, cooling_keys=IDs_cooling)
    metrics = NamedTuple[]
    for (keys, is_cooling) in ((heating_keys, false), (cooling_keys, true))
        for sim_id in keys
            case_data = run_single_case_2D(sim_id, p; is_cooling=is_cooling)
            model = case_data.model
            exp = case_data.experiment
            times = case_data.times
            
            for (j, sensor) in enumerate(SENSOR_NAMES_2D)
                rmse = sqrt(mean(abs2, model[:, j] .- exp[:, j]))
                steady_err = model[end, j] - exp[end, j]
                t90_err = get_t90_2D(times, model[:, j]) - get_t90_2D(times, exp[:, j])
                shape_err = normalized_slope_mse_2D(model[:, j], exp[:, j])
                
                push!(metrics, (
                    simulation_id=sim_id,
                    phase=is_cooling ? :cooling_benchmark : :heating_calibrated,
                    sensor=sensor,
                    rmse_K=rmse,
                    steady_error_K=steady_err,
                    t90_error_s=t90_err,
                    shape_loss=shape_err,
                ))
            end
        end
    end
    return metrics
end

function slope_2D(xs, ys)
    xmean = mean(xs)
    ymean = mean(ys)
    denom = sum((x - xmean)^2 for x in xs)
    denom <= eps(Float64) && return NaN
    return sum((xs[i] - xmean) * (ys[i] - ymean) for i in eachindex(xs)) / denom
end

function build_steady_results_2D(p=default_parameters2D(); keys=IDs)
    steady_rows = NamedTuple[]
    for sim_id in keys
        case_data = run_single_case_2D(sim_id, p; is_cooling=false)
        conds = simulation_conditions[sim_id]
        model = case_data.model
        exp = case_data.experiment
        
        push!(steady_rows, (
            simulation_id=sim_id,
            flow_lpm=conds[qlpm],
            irradiance=conds[Io],
            T8_model=model[end, 1], T8_experiment=exp[end, 1],
            T12_perim_model=model[end, 2], T12_perim_experiment=exp[end, 2],
            T11_perim_model=model[end, 3], T11_perim_experiment=exp[end, 3],
            T9_core_model=model[end, 4], T9_core_experiment=exp[end, 4],
            T10_core_model=model[end, 5], T10_core_experiment=exp[end, 5],
            T3_model=model[end, 6], T3_experiment=exp[end, 6],
            T2_model=model[end, 7], T2_experiment=exp[end, 7],
            T12_minus_T9_model=model[end, 2] - model[end, 4],
            T12_minus_T9_experiment=exp[end, 2] - exp[end, 4],
            T11_minus_T10_model=model[end, 3] - model[end, 5],
            T11_minus_T10_experiment=exp[end, 3] - exp[end, 5],
        ))
    end
    return steady_rows
end

function plot_steady_parity_2D(steady_results)
    p_plot = plot(
        title="2D_v7 Heating-Calibrated Steady-State Parity Plot",
        xlabel="Model Temperature (K)",
        ylabel="Experiment Temperature (K)",
        legend=:outerright,
        aspect_ratio=:equal,
        xlims=TEMPERATURE_AXIS_LIMITS_2D,
        ylims=TEMPERATURE_AXIS_LIMITS_2D,
    )
    for sensor in SENSOR_NAMES_2D
        model_vals = [getproperty(row, Symbol(sensor, "_model")) for row in steady_results]
        exp_vals = [getproperty(row, Symbol(sensor, "_experiment")) for row in steady_results]
        scatter!(p_plot, model_vals, exp_vals;
                 label=String(sensor), color=SENSOR_COLORS_2D[sensor],
                 ms=5, markerstrokewidth=0)
    end
    lower, upper = TEMPERATURE_AXIS_LIMITS_2D
    plot!(p_plot, [lower, upper], [lower, upper]; label="1:1 Parity", color=:gray, linestyle=:dash)
    path = joinpath(PLOT_DIR_2D_v7, "steady_comparison_2D_v7.png")
    savefig(p_plot, path)
    return path
end

function plot_axial_profile_2D(sim_id, case_data)
    result = case_data.result
    exp = case_data.experiment
    conds = simulation_conditions[sim_id]

    core_T = result.solid_temperature[1, :, end]
    perim_T = result.solid_temperature[min(10, size(result.solid_temperature, 1)), :, end]
    casing_T = result.solid_temperature[end, :, end]
    gas_T = mean(result.gas_temperature[:, :, end], dims=1)

    z_solid_mm = result.z_solid .* 1000.0
    z_gas_mm = result.z_gas .* 1000.0

    p_plot = plot(
        title="2D_v7 Multi-Domain Axial Profiles: $sim_id\nflow=$(conds[qlpm]) L/min, Io=$(round(conds[Io]/1000, digits=1)) kW/m^2",
        xlabel="Axial Position z (mm)", ylabel="Temperature (K)",
        legend=:outerright, grid=true,
        xlims=(0.0, 150.0), ylims=TEMPERATURE_AXIS_LIMITS_2D
    )

    plot!(p_plot, z_solid_mm, core_T, label="Core Solid (r=0)", color=:orange, lw=2.5, linestyle=:dot)
    plot!(p_plot, z_solid_mm, perim_T, label="Perimeter SiC (r=R_rec)", color=:red, lw=2.5, linestyle=:dash)
    plot!(p_plot, z_solid_mm, casing_T, label="Outer Casing (r=R_case)", color=:black, lw=1.5, linestyle=:dashdot)
    plot!(p_plot, z_gas_mm, vec(gas_T), label="Mixed Gas Profile", color=:green, lw=2.5)

    scatter!(p_plot, [11.0], [exp[end, 1]], label="T8 (Perim 11mm)", color=:blue, markershape=:square, ms=6)
    scatter!(p_plot, [58.0], [exp[end, 4]], label="T9 (Core 58mm)", color=:orange, markershape=:circle, ms=6)
    scatter!(p_plot, [58.0], [exp[end, 2]], label="T12 (Perim 58mm)", color=:red, markershape=:utriangle, ms=6)
    scatter!(p_plot, [107.0], [exp[end, 5]], label="T10 (Core 107mm)", color=:brown, markershape=:diamond, ms=6)
    scatter!(p_plot, [107.0], [exp[end, 3]], label="T11 (Perim 107mm)", color=:darkgreen, markershape=:dtriangle, ms=6)
    scatter!(p_plot, [140.0], [exp[end, 6]], label="T3 Gas (Exit)", color=:purple, markershape=:star5, ms=8)
    scatter!(p_plot, [58.0], [exp[end, 7]], label="T2 Felt (Insulation)", color=:black, markershape=:pentagon, ms=6)

    path = joinpath(AXIAL_DIR_2D_v7, "axial_profile_$(sim_id)_2D_v7.png")
    savefig(p_plot, path)
    return path
end

function plot_transient_2D(sim_id, case_data; is_cooling=false)
    times = case_data.times
    model = case_data.model
    exp = case_data.experiment
    
    p_plot = plot(
        title="2D_v7 Fitted Transient: $sim_id $(is_cooling ? "(Cooling Benchmark)" : "(Heating Calibrated)")",
        xlabel="Time (s)", ylabel="Temperature (K)",
        legend=:outerright, grid=true, ylims=TEMPERATURE_AXIS_LIMITS_2D
    )
    
    for (j, sensor) in enumerate(SENSOR_NAMES_2D)
        col = SENSOR_COLORS_2D[sensor]
        style = SENSOR_STYLES_2D[sensor]
        plot!(p_plot, times, model[:, j], label="$(sensor) Model", lw=3.0, color=col, linestyle=style)
        scatter!(p_plot, times, exp[:, j], label="$(sensor) Exp", ms=2.5, color=col, markerstrokewidth=0)
    end
    
    suffix = is_cooling ? "_cooling_benchmark" : ""
    path = joinpath(TRANSIENT_DIR_2D_v7, "transient_$(sim_id)$(suffix)_2D_v7.png")
    savefig(p_plot, path)
    return path
end

function plot_2d_heatmap(sim_id, result)
    r_mm = result.r_solid .* 1000.0
    z_mm = result.z_solid .* 1000.0
    T_final = result.solid_temperature[:, :, end]
    
    p_heat = heatmap(
        z_mm, r_mm, T_final,
        title="2D_v7 Multi-Domain Temperature Field T_s(r, z) [K] - $sim_id",
        xlabel="Axial Position z (mm)", ylabel="Radial Radius r (mm)",
        c=:inferno, aspect_ratio=:equal, clim=(minimum(T_final), maximum(T_final))
    )
    
    path = joinpath(PROFILES_DIR_2D_v7, "heatmap_2D_$(sim_id)_2D_v7.png")
    savefig(p_heat, path)
    return path
end

function write_cell_diagnostics_2D(path, sim_id, result)
    open(path, "w") do io
        println(io, "simulation_id,r_mm,z_mm,T_solid_K,T_gas_K,h_W_m2K")
        nr_rec = size(result.gas_temperature, 1)
        for i in 1:nr_rec, j in 1:length(result.z_solid)
            println(io, join((
                sim_id,
                result.r_solid[i] * 1000.0,
                result.z_solid[j] * 1000.0,
                result.solid_temperature[i, j, end],
                result.gas_temperature[i, j, end],
                result.heat_transfer_coefficient[i, j, end],
            ), ','))
        end
    end
    return path
end

function write_parameters_2D(path, p::ModelParameters2D)
    names = (
        "A_Nu", "B_Re", "k_scale_r", "k_scale_z", "sigma_beam_mm",
        "spillage_fraction", "C_cavity_eff_J_K", "scale_456", "scale_304", "scale_256",
        "A_entry", "sic_k_temp_scale", "beta_rad", "perim_gap_resistance", "beta_opt"
    )
    vals = pack_parameters2D(p)
    open(path, "w") do io
        println(io, "index,name,value")
        for i in eachindex(vals)
            println(io, join((i, names[i], vals[i]), ','))
        end
    end
    return path
end

begin # Main Execution Loop
    println("==========================================================================")
    println("Running 2D_v7 Multi-Domain HEATING-ONLY Calibration with Radiative Transport...")
    println("==========================================================================")
    
    p_initial = default_parameters2D()
    
    calib = calibrate2D(
        run_single_case_2D;
        heating_cases=("E67", "E76", "E80"),
        max_iters=150,
        max_time=1800.0,
        p_init=p_initial,
    )
    
    p_fitted = calib.parameters
    
    open(joinpath(OUTPUT_DIR_2D_v7, "optimization_summary_2D_v7.txt"), "w") do io
        println(io, "objective=$(calib.objective)")
        println(io, "return_code=$(calib.retcode)")
        println(io, "parameters=$(pack_parameters2D(p_fitted))")
    end
    
    write_parameters_2D(joinpath(OUTPUT_DIR_2D_v7, "parameters_fitted_2D_v7.csv"), p_fitted)
    
    println("[run_2D_v7] Computing metrics and generating figures for fitted 2D_v7 model...")
    metrics = compute_all_metrics_2D(p_fitted)
    open(joinpath(OUTPUT_DIR_2D_v7, "analysis_results_2D_v7.csv"), "w") do io
        println(io, "simulation_id,phase,sensor,rmse_K,steady_error_K,t90_error_s,shape_loss")
        for row in metrics
            println(io, join((row.simulation_id, row.phase, row.sensor,
                              row.rmse_K, row.steady_error_K,
                              row.t90_error_s, row.shape_loss), ','))
        end
    end

    steady_results = build_steady_results_2D(p_fitted)
    steady_path = joinpath(OUTPUT_DIR_2D_v7, "steady_results_2D_v7.csv")
    fields = propertynames(first(steady_results))
    open(steady_path, "w") do io
        println(io, join(fields, ','))
        for row in steady_results
            println(io, join((getproperty(row, f) for f in fields), ','))
        end
    end

    plot_steady_parity_2D(steady_results)

    slopes_path = joinpath(OUTPUT_DIR_2D_v7, "flow_slopes_2D_v7.csv")
    irradiances = sort(unique(row.irradiance for row in steady_results); rev=true)
    open(slopes_path, "w") do io
        println(io, "irradiance,sensor,model_slope_K_per_Lmin,experiment_slope_K_per_Lmin")
        for irr in irradiances
            rows = [row for row in steady_results if row.irradiance == irr]
            flows = [row.flow_lpm for row in rows]
            for sensor in SENSOR_NAMES_2D
                model_vals = [getproperty(row, Symbol(sensor, "_model")) for row in rows]
                exp_vals = [getproperty(row, Symbol(sensor, "_experiment")) for row in rows]
                println(io, join((irr, sensor, slope_2D(flows, model_vals), slope_2D(flows, exp_vals)), ','))
            end
        end
    end

    for sim_id in IDs
        case_data = run_single_case_2D(sim_id, p_fitted; is_cooling=false)
        plot_transient_2D(sim_id, case_data; is_cooling=false)
        plot_axial_profile_2D(sim_id, case_data)
        plot_2d_heatmap(sim_id, case_data.result)
        write_cell_diagnostics_2D(
            joinpath(DIAGNOSTICS_DIR_2D_v7, "cell_diagnostics_$(sim_id)_2D_v7.csv"),
            sim_id, case_data.result
        )
    end

    for sim_id in IDs_cooling
        case_data = run_single_case_2D(sim_id, p_fitted; is_cooling=true)
        plot_transient_2D(sim_id, case_data; is_cooling=true)
    end

    println("==========================================================================")
    println("Multi-Domain 2D_v7 Calibration and Runner Execution Complete.")
    println("Final Objective Loss: $(calib.objective)")
    println("Fitted Parameter Vector: $(pack_parameters2D(p_fitted))")
    println("==========================================================================")
end
