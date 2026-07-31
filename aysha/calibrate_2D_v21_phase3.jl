# ============================================================================
# calibrate_2D_v21_phase3.jl
#
# Phase 3 Formal Calibration: Grid search over macroscopic reduced-order
# physics introduced in Phase 2 (spillage_fraction, core_preference).
# ============================================================================

using Statistics

# Include run_2D_v20.jl to get case inputs (data loading, etc)
include("run_2D_v20.jl")

include("2D_v21.jl")
using .Receiver2D_v21
const V21 = Receiver2D_v21

const SPILLAGE_GRID = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
const CORE_PREF_GRID = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
const BASE_POWER_SCALES = (1.19, 1.54, 0.94)

const GROUP_IDS_PHASE3 = Dict(
    456 => ("E67", "E69", "E71"),
)
const SENSOR_NAMES_PHASE3 = ("T8", "T12", "T11", "T9", "T10", "T3", "T2")

function build_base_parameters(power_scales)
    # Start from V21.V20.default_parameters2D()
    def = V21.V20.default_parameters2D()
    
    # Unpack all the way down to V11 where OpticalParameters2D lives
    v19_old = def.base
    v18_old = v19_old.base
    v17_old = v18_old.base
    v15_old = v17_old.base
    v14_old = v15_old.base
    v12_old = v14_old.base
    v11_old = v12_old.base
    
    # Rebuild V11 with new optical scales
    v11_new = V21.V20.V11.ModelParameters2D(
        v11_old.geometry, v11_old.solid, v11_old.heat_transfer, v11_old.losses,
        V21.V20.V11.OpticalParameters2D(
            absorbed_fraction=v11_old.optics.absorbed_fraction,
            extinction_coefficient=v11_old.optics.extinction_coefficient,
            beam_radius_sigma=v11_old.optics.beam_radius_sigma,
            spillage_fraction=v11_old.optics.spillage_fraction,
            front_deposition_fraction=v11_old.optics.front_deposition_fraction,
            scale_456=power_scales[1],
            scale_304=power_scales[2],
            scale_256=power_scales[3]
        ),
        v11_old.hydraulics
    )
    
    # Rebuild the rest of the stack
    v12_new = V21.V20.V12.ModelParameters2D(base=v11_new, assembly=v12_old.assembly, observation=v12_old.observation)
    v14_new = V21.V20.V14.ModelParameters2D(base=v12_new, network=v14_old.network)
    v15_new = V21.V20.V15.ModelParameters2D(base=v14_new, skin=v15_old.skin)
    v17_new = V21.V20.V17.ModelParameters2D(base=v15_new, casing_flange=v17_old.casing_flange)
    v18_new = V21.V20.V18.ModelParameters2D(base=v17_new, source=v18_old.source)
    v19_new = V21.V20.V19.ModelParameters2D(base=v18_new, integrated_exchange=v19_old.integrated_exchange, outlet_observation=v19_old.outlet_observation, rear_tube_flange=v19_old.rear_tube_flange)
    
    # Final V20 base params
    v20_new = V21.V20.ModelParameters2D(base=v19_new, t3_location=def.t3_location)
    
    return v20_new
end


function _final_rows_phase3(times)
    length(times) >= 2 || error("calibration requires at least two samples")
    threshold = first(times) + 0.90 * (last(times) - first(times))
    start = something(findfirst(time -> time >= threshold, times), length(times))
    return start:length(times)
end

function _calibration_metrics_phase3(case)
    final = _final_rows_phase3(case.inputs.times)
    residual = case.model .- case.observed
    side_final = [mean(residual[final, sensor]) for sensor in 1:3]
    T2_final = mean(residual[final, 7])
    
    # Exclude T3 (gas outlet) from objective entirely in Phase 3
    final_objective = mean((side_final ./ 50.0).^2) / 2.0 + (T2_final / 35.0)^2 / 2.0
    
    transient_objective = mean((residual[:, 1:3] ./ 75.0).^2) / 2.0 + mean((residual[:, 7] ./ 50.0).^2) / 2.0
    
    return (
        objective = 0.8 * final_objective + 0.2 * transient_objective,
        side_final_sse = sum(abs2, side_final),
        T2_final_sse = T2_final^2,
        side_final_bias_K = mean(side_final),
        side_final_rmse_K = sqrt(mean(abs2, side_final)),
        T2_final_bias_K = T2_final,
        transient_side_sse = sum(abs2, residual[:, 1:3]),
        transient_side_n = length(residual[:, 1:3]),
        transient_T2_sse = sum(abs2, residual[:, 7]),
        transient_T2_n = length(residual[:, 7])
    )
end

function simulate_case_phase3(inputs, p_mod::V21.ModelParameters2D_v21; full_initial_data=true, reltol=5e-4, abstol=1e-4, dtmax=120.0)
    scale = inputs.nominal >= 400000.0 ? p_mod.base.base.base.base.base.base.base.optics.scale_456 :
            (inputs.nominal >= 280000.0 ? p_mod.base.base.base.base.base.base.base.optics.scale_304 :
             p_mod.base.base.base.base.base.base.base.optics.scale_256)
    
    irradiance = inputs.cooling ? zeros(length(inputs.times)) :
                 fill(scale * inputs.nominal, length(inputs.times))
                 
    op = V21.V20.OperatingCondition2D(
        irradiance = V21.V20.linear_history(inputs.times, irradiance),
        flow_lpm = V21.V20.linear_history(inputs.times, inputs.flow),
        inlet_temperature = V21.V20.linear_history(inputs.times, inputs.inlet),
        ambient_temperature = V21.V20.linear_history(inputs.times, inputs.ambient),
    )
    
    initial = if inputs.cooling && full_initial_data
        t0 = Dict(
            sensor => Float64(Main.observation(inputs.data, inputs.id, "_$(sensor)")[1])
            for sensor in SENSOR_NAMES_PHASE3
        )
        grid = V21.V20.V11.build_grid2D(p_mod.base.base.base.base.base.base.base.base)
        V21.V20.build_initial_state_2D(grid, p_mod.base.base.base.base.base.base.base.base, t0, inputs.ambient[1])
    else
        inputs.ambient[1]
    end
    
    result = V21.simulate2D_v21(
        p_mod, op, inputs.times;
        initial_temperature=initial,
        reltol=reltol, abstol=abstol, dtmax=dtmax
    )
    predictions = V21.V20.sensor_predictions2D(result)
    model_matrix = hcat((predictions[Symbol(sensor)] for sensor in SENSOR_NAMES_PHASE3)...)
    return (
        inputs=inputs,
        parameters=p_mod,
        result=result,
        model=model_matrix,
        observed=inputs.observed
    )
end

function _profile_candidate_phase3(spillage_fraction, core_pref; max_points=21)
    base_p = build_base_parameters(BASE_POWER_SCALES)
    
    spillage_power_W = 1500.0 * spillage_fraction
    
    p_mod = V21.ModelParameters2D_v21(
        base = base_p,
        phase2 = V21.Phase2Parameters2D(
            spillage_power_W = spillage_power_W,
            core_preference = core_pref,
            spillage_axial_spread_m = 10.0e-3
        )
    )
    
    metrics = NamedTuple[]
    for id in GROUP_IDS_PHASE3[456]
        inputs_heat = Main.case_inputs_2D_v20(id; cooling=false, max_points=max_points)
        case_heat = simulate_case_phase3(inputs_heat, p_mod; full_initial_data=false)
        push!(metrics, _calibration_metrics_phase3(case_heat))
    end
    
    return (
        spillage_fraction = spillage_fraction,
        core_preference = core_pref,
        objective = mean(row.objective for row in metrics),
        side_final_rmse_K = sqrt(sum(row.side_final_sse for row in metrics) / (3 * length(metrics))),
        T2_final_rmse_K = sqrt(mean(row.T2_final_sse for row in metrics)),
        side_final_bias_K = mean(row.side_final_bias_K for row in metrics),
        T2_final_bias_K = mean(row.T2_final_bias_K for row in metrics),
        transient_side_rmse_K = sqrt(sum(row.transient_side_sse for row in metrics) / sum(row.transient_side_n for row in metrics)),
        transient_T2_rmse_K = sqrt(sum(row.transient_T2_sse for row in metrics) / sum(row.transient_T2_n for row in metrics)),
    )
end

function calibrate_phase3(; max_points=21)
    V21.apply_v21_property_fixes!()
    
    outdir = joinpath(@__DIR__, "summaries", "2D_v21")
    mkpath(outdir)
    candidates = [
        (spillage_fraction=sf, core_pref=cp)
        for sf in SPILLAGE_GRID
        for cp in CORE_PREF_GRID
    ]
    
    rows = Vector{NamedTuple}(undef, length(candidates))
    for index in eachindex(candidates)
        candidate = candidates[index]
        println("Phase 3 profile $index/$(length(candidates)): Spillage=", candidate.spillage_fraction, " CorePref=", candidate.core_pref)
        flush(stdout)
        rows[index] = _profile_candidate_phase3(candidate.spillage_fraction, candidate.core_pref; max_points=max_points)
    end
    
    open(joinpath(outdir, "phase3_calibration_grid.csv"), "w") do io
        keys_tuple = keys(rows[1])
        println(io, join(String.(keys_tuple), ","))
        for row in rows
            println(io, join([string(row[k]) for k in keys_tuple], ","))
        end
    end
    
    best = rows[argmin(row.objective for row in rows)]
    println("--- BEST PHASE 3 PARAMETERS ---")
    println("Spillage Fraction: ", best.spillage_fraction)
    println("Core Preference:   ", best.core_preference)
    println("Objective:         ", best.objective)
    println("Side Final RMSE:   ", best.side_final_rmse_K)
    println("T2 Final RMSE:     ", best.T2_final_rmse_K)
end

if abspath(PROGRAM_FILE) == @__FILE__
    calibrate_phase3(max_points=21)
end
