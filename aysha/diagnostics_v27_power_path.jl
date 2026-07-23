include("1D_v27.jl")

const FITTED_PARAMS_v27_POWER_PATH = copy(pnew_v27)

const PARAMETER_NAMES_v27_POWER_PATH = (
    "A_Nu", "B_Re", "C_Pr_fixed", "phi_0", "m_rec", "front_dep",
    "scale_456", "scale_304", "scale_256", "G_core_perim_W_m_K",
    "C_perim_eff_J_K", "k_perim_ref_W_m_K", "beta_opt", "spill_capture",
    "beta_perim_1_m", "f_core_tube", "flange_scale", "G_rear_core_heat_W_m_K",
    "G_rear_perim_W_m_K", "rear_sink_shape", "k_core_axial_scale",
    "G_rear_core_cool_W_m_K", "f_T3_wall",
)

function print_parameter_limiter_report_v27(p)
    println("parameter,value,lower,upper,near_lower,near_upper")
    for i in eachindex(p)
        lower = lb_full_v27[i]
        upper = ub_full_v27[i]
        span = max(upper - lower, eps(Float64))
        near_lower = (p[i] - lower) / span < 0.02
        near_upper = (upper - p[i]) / span < 0.02
        println(join((PARAMETER_NAMES_v27_POWER_PATH[i], p[i], lower, upper, near_lower, near_upper), ','))
    end
end

function receiver_rear_sink_totals_v27(result, p, irradiance)
    z = result.z_solid
    dx = L / length(z)
    Tcore = result.core_temperature[:, end]
    Tperim = result.perim_temperature[:, end]
    core = 0.0
    perim = 0.0
    for i in eachindex(z)
        weight = rear_receiver_sink_weight_v27(z[i], p)
        core += rear_core_loss_per_length_v27(p, irradiance) * dx * weight *
                (Tcore[i] - WATER_FLANGE_TEMPERATURE_v27)
        perim += rear_perim_loss_per_length_v27(p) * dx * weight *
                 (Tperim[i] - WATER_FLANGE_TEMPERATURE_v27)
    end
    return core, perim
end

function tube_loss_totals_v27(result, p)
    z_local = result.z_rear_tube .- L
    dx = REAR_TUBE_LENGTH_v27 / length(z_local)
    Ttube = result.rear_tube_temperature[:, end]
    Tcavity = result.cavity_temperature[end]
    cavity = 0.0
    flange = 0.0
    for j in eachindex(z_local)
        G_cavity = rear_tube_cavity_conductance_per_length_v27(z_local[j], p)
        G_flange = rear_tube_flange_conductance_per_length_v27(z_local[j], Ttube[j], p)
        cavity += G_cavity * dx * (Ttube[j] - Tcavity)
        flange += G_flange * dx * (Ttube[j] - WATER_FLANGE_TEMPERATURE_v27)
    end
    return cavity, flange
end

function receiver_to_tube_heat_v27(result, p)
    Tcore_exit = result.core_temperature[end, end]
    Tperim_exit = result.perim_temperature[end, end]
    Ttube_in = result.rear_tube_temperature[1, end]
    f_core = core_tube_fraction_v27(p)
    core = RECEIVER_TO_TUBE_CONDUCTANCE_v27 * f_core * (Tcore_exit - Ttube_in)
    perim = RECEIVER_TO_TUBE_CONDUCTANCE_v27 * (1.0 - f_core) * (Tperim_exit - Ttube_in)
    return core, perim
end

function print_case_power_path_v27(simulation_id, p; nodes=15)
    outputs, result, experiment = solve_case_v27(p, simulation_id; nodes=nodes)
    conditions = simulation_conditions[simulation_id]
    irradiance = conditions[Io]
    _, _, _, case_data = experimental_case_v27(simulation_id)
    ambient = observation(case_data, simulation_id, "_Tamb")[end]
    scale = absorbed_power_scale_v27(irradiance, p)
    Qtotal = modeled_absorbed_receiver_power_v27(irradiance, p)
    Qcore = core_absorbed_power_v27(irradiance, p)
    Qperim = perimeter_absorbed_power_v27(irradiance, p)
    Qgas = sum(result.receiver_gas_heat[:, end])
    Qrear_gas = sum(result.rear_tube_gas_heat[:, end])
    Tfront = result.perim_temperature[1, end]
    Qfront = H_FRONT_FIXED_v27 * A_frt * (Tfront - ambient) +
             EPS_FRONT_FIXED_v27 * sigma_sb * A_frt * (Tfront^4 - ambient^4)
    Qcavity_ambient = cavity_loss_to_ambient_v27(result.cavity_temperature[end], ambient)
    Qrear_core, Qrear_perim = receiver_rear_sink_totals_v27(result, p, irradiance)
    Qtube_cavity, Qtube_flange = tube_loss_totals_v27(result, p)
    Qreceiver_tube_core, Qreceiver_tube_perim = receiver_to_tube_heat_v27(result, p)
    println(join((
        simulation_id,
        conditions[qlpm],
        irradiance,
        scale,
        Qtotal,
        Qcore,
        Qperim,
        Qgas,
        Qrear_gas,
        Qfront,
        Qcavity_ambient,
        Qrear_core,
        Qrear_perim,
        Qtube_cavity,
        Qtube_flange,
        Qreceiver_tube_core,
        Qreceiver_tube_perim,
        outputs[end, 1] - experiment[end, 1],
        outputs[end, 2] - experiment[end, 2],
        outputs[end, 3] - experiment[end, 3],
        outputs[end, 4] - experiment[end, 4],
        outputs[end, 5] - experiment[end, 5],
        outputs[end, 6] - experiment[end, 6],
        outputs[end, 7] - experiment[end, 7],
    ), ','))
end

println("# parameter bounds")
print_parameter_limiter_report_v27(FITTED_PARAMS_v27_POWER_PATH)

println("# power path")
println("case,flow_lpm,irradiance,power_scale,Qtotal_abs_W,Qcore_abs_W,Qperim_abs_W,Qreceiver_gas_W,Qrear_gas_W,Qfront_loss_W,Qcavity_ambient_W,Qrear_core_sink_W,Qrear_perim_sink_W,Qtube_cavity_W,Qtube_flange_W,Qreceiver_tube_core_W,Qreceiver_tube_perim_W,T8_err_K,T12_err_K,T11_err_K,T9_err_K,T10_err_K,T3_err_K,T2_err_K")
for simulation_id in ("E67", "E69", "E71", "E72", "E76", "E77", "E80", "E81")
    print_case_power_path_v27(simulation_id, FITTED_PARAMS_v27_POWER_PATH)
end
