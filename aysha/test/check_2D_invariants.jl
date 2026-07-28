# ============================================================================
# check_2D_invariants.jl - Direct Emergent Invariant Extraction (I1-I6)
# ============================================================================

using Statistics

begin # Dictionary Key Symbols
    Io = :Io
    qlpm = :qlpm
    inlet_temperature = :inlet_temperature
    ambient_temperature = :ambient_temperature
    time_data = :time_data
    Tinit = :Tinit
    T_in = :T_in
    T_amb = :T_amb
end

include(joinpath(@__DIR__, "..", "2D_v3.jl"))
include(joinpath(@__DIR__, "..", "import_exp_1D_v2.jl"))

using .Receiver2D_v3

println("==========================================================================")
println("Extracting Emergent Manuscript Invariants (I1-I6) from 2D_v3 Model...")
println("==========================================================================")

p_fitted = unpack_parameters2D([2.5, 0.333, 0.030, 1.0, 14.0, 0.2, 301.0, 0.55, 0.55, 0.45, 1.0, 0.0, 0.2, 1.675, 0.1])

const COOLING_HEATING_MAP = Dict("C69" => "E69", "C80" => "E80", "C81" => "E81")

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

function run_case(sim_id, p; is_cooling=false)
    data = is_cooling ? measurements_cooling : measurements
    conds = is_cooling ? simulation_conditions_cooling : simulation_conditions
    conditions = conds[sim_id]
    save_times = Float64.(observation_time(data, sim_id))
    flow = Float64.(observation(data, sim_id, "_flow"))
    Tin = Float64.(observation(data, sim_id, "_Tin"))
    Tamb = Float64.(observation(data, sim_id, "_Tamb"))
    
    irr_val = Float64(conditions[Io])
    scale_factor = irr_val >= 400000.0 ? p.optics.scale_456 : (irr_val >= 280000.0 ? p.optics.scale_304 : p.optics.scale_256)
    irradiance = is_cooling ? zeros(length(save_times)) : fill(scale_factor * irr_val, length(save_times))
    
    op = OperatingCondition2D(
        irradiance = linear_history(save_times, irradiance),
        flow_lpm = linear_history(save_times, flow),
        inlet_temperature = linear_history(save_times, Tin),
        ambient_temperature = linear_history(save_times, Tamb),
    )
    
    grid = Receiver2D_v3.build_grid2D(p)
    if is_cooling
        t0_sensors = get_experimental_t0_dict(sim_id; is_cooling=true)
        u0_exp = build_initial_state_2D(grid, p, t0_sensors, Tamb[1])
        result = simulate2D(p, op, save_times; initial_temperature=u0_exp)
    else
        result = simulate2D(p, op, save_times)
    end
    return result
end

# Invariant I3: Deep Nonequilibrium Deficit T_core(107) - T_perim(107) vs Re
println("\n--- Invariant I3: Deep Nonequilibrium Deficit (107 mm) ---")
for sim_id in ["E67", "E72", "E76", "E77", "E80", "E81"]
    res = run_case(sim_id, p_fitted; is_cooling=false)
    conds = simulation_conditions[sim_id]
    flow_lpm = Float64(conds[qlpm])
    
    T10_core = res.solid_temperature[1, 20, end]
    T11_perim = res.solid_temperature[10, 20, end]
    deficit = T10_core - T11_perim
    println("  $sim_id: Flow = $(flow_lpm) L/min -> T10_core - T11_perim = $(round(deficit, digits=1)) K")
end

# Invariant I5: Total Parasitic Loss Conductance K_loss [W/K]
println("\n--- Invariant I5: Total Parasitic Loss Conductance K_loss ---")
for sim_id in ["E67", "E76", "E80"]
    res = run_case(sim_id, p_fitted; is_cooling=false)
    T_s = res.solid_temperature[:, :, end]
    Tamb = 295.0
    
    # Casing loss
    dz = res.z_solid[2] - res.z_solid[1]
    Q_casing = sum(10.0 * (2.0 * pi * 75.0e-3 * dz) * (T_s[end, j] - Tamb) for j in 1:size(T_s, 2))
    dT_excess = mean(T_s[1:10, :, end]) - Tamb
    K_loss = Q_casing / max(dT_excess, 1.0)
    println("  $sim_id: K_loss = $(round(K_loss, digits=3)) W/K  (Manuscript Target: 0.10 - 0.16 W/K)")
end

println("\n==========================================================================")
println("Invariant Extraction Complete.")
println("==========================================================================")
