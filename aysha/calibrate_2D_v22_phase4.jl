# ============================================================================
# calibrate_2D_v22_phase4.jl
#
# Phase 4: Microscopic Parameter Calibration
# Sweeps over Nu_multiplier and Rosseland_multiplier using fixed macroscopic bounds.
# ============================================================================

using Statistics
using Plots
using CSV
using DataFrames

include("run_2D_v22.jl")

const SENSOR_NAMES_PHASE4 = ("T8", "T12", "T11", "T9", "T10", "T3", "T2")
const SENSOR_WEIGHTS_PHASE4 = [1.0, 1.0, 1.0, 2.0, 2.0, 0.5, 0.5] # Bias towards core/side fit

const GROUP_IDS_PHASE4 = Dict(
    456 => ("E67", "E69", "E71")
)

function build_base_parameters(nu_mult, ross_mult)
    # Set globals BEFORE constructing V22 parameters
    Main.Receiver2D_v22.Receiver2D_v22_Properties.GLOBAL_NU_MULTIPLIER[] = nu_mult
    Main.Receiver2D_v22.Receiver2D_v22_Properties.GLOBAL_ROSS_MULTIPLIER[] = ross_mult
    
    # Get the correctly constructed parameters with experimental power scales etc.
    p_v22_defaults = Main.parameters_2D_v22()
    
    p = Main.Receiver2D_v22.ModelParameters2D_v22(
        base = p_v22_defaults,
        phase2 = Main.Receiver2D_v22.Phase2Parameters2D(
            spillage_power_W = 0.05 * 10500.0, # 5% baseline spillage approximation
            core_preference = 0.6,             # Physically realistic core preference
            spillage_axial_spread_m = 10.0e-3
        ),
        phase4 = Main.Receiver2D_v22.Phase4Parameters2D(
            nu_multiplier = nu_mult,
            ross_multiplier = ross_mult
        )
    )
    return p
end

function simulate_case_phase4(inputs, p; full_initial_data=false)
    return Main.simulate_case_2D_v22(inputs, p; full_initial_data=full_initial_data)
end

function _calibration_metrics_phase4(case)
    times = case.inputs.times
    threshold = first(times) + 0.90 * (last(times) - first(times))
    final = something(findfirst(t -> t >= threshold, times), length(times)):length(times)
    
    # 1:3 are side (T8, T12, T11)
    # 4:5 are core (T9, T10)
    # 7 is rear (T2)
    side_final_sse = sum(abs2, (case.model[final, 1:3] .- case.observed[final, 1:3]) .* SENSOR_WEIGHTS_PHASE4[1:3]')
    core_final_sse = sum(abs2, (case.model[final, 4:5] .- case.observed[final, 4:5]) .* SENSOR_WEIGHTS_PHASE4[4:5]')
    T2_final_sse = sum(abs2, (case.model[final, 7] .- case.observed[final, 7]) .* SENSOR_WEIGHTS_PHASE4[7])
    
    return (
        side_final_sse = side_final_sse,
        core_final_sse = core_final_sse,
        T2_final_sse = T2_final_sse,
        objective = side_final_sse + core_final_sse + T2_final_sse
    )
end

function _profile_candidate_phase4(nu_mult, ross_mult; max_points=21)
    p_mod = build_base_parameters(nu_mult, ross_mult)
    
    metrics = NamedTuple[]
    for id in GROUP_IDS_PHASE4[456]
        inputs_heat = Main.case_inputs_2D_v22(id; cooling=false, max_points=max_points)
        case_heat = simulate_case_phase4(inputs_heat, p_mod; full_initial_data=false)
        push!(metrics, _calibration_metrics_phase4(case_heat))
    end
    
    return (
        nu_multiplier = nu_mult,
        ross_multiplier = ross_mult,
        objective = mean(row.objective for row in metrics),
        side_final_rmse_K = sqrt(mean(row.side_final_sse for row in metrics) / 3),
        core_final_rmse_K = sqrt(mean(row.core_final_sse for row in metrics) / 2),
        T2_final_rmse_K = sqrt(mean(row.T2_final_sse for row in metrics)),
    )
end

function calibrate_phase4(; max_points=21)
    Main.Receiver2D_v22.apply_v22_property_fixes!()
    
    outdir = joinpath(@__DIR__, "summaries", "2D_v22")
    mkpath(outdir)
    
    nu_grid = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
    ross_grid = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    
    results = NamedTuple[]
    total = length(nu_grid) * length(ross_grid)
    idx = 1
    
    for nu in nu_grid
        for ross in ross_grid
            println("Phase 4 profile $idx/$total: Nu_Mult=$nu Ross_Mult=$ross")
            try
                res = _profile_candidate_phase4(nu, ross; max_points=max_points)
                push!(results, res)
            catch e
                println("Failed on Nu=$nu, Ross=$ross: $e")
            end
            idx += 1
        end
    end
    
    df = DataFrame(results)
    CSV.write(joinpath(outdir, "phase4_calibration_grid.csv"), df)
    
    best_idx = argmin(df.objective)
    best = df[best_idx, :]
    
    println("--- BEST PHASE 4 PARAMETERS ---")
    println("Nu Multiplier:    ", best.nu_multiplier)
    println("Ross Multiplier:  ", best.ross_multiplier)
    println("Objective:        ", best.objective)
    println("Side Final RMSE:  ", best.side_final_rmse_K)
    println("Core Final RMSE:  ", best.core_final_rmse_K)
    
    return best
end

if abspath(PROGRAM_FILE) == @__FILE__
    calibrate_phase4()
end

