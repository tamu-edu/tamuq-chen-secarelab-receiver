include("1D_v30.jl")

const FITTED_PARAMS_v30_POWER_PATH = copy(pnew_v30)

const PARAMETER_NAMES_v30_POWER_PATH = (
    "A_Nu", "B_Re", "C_Pr_fixed", "phi_0", "m_rec", "front_dep",
    "scale_456", "scale_304", "scale_256", "G_core_perim_W_m_K",
    "C_perim_eff_J_K", "k_perim_ref_W_m_K", "beta_opt", "spill_capture",
    "beta_perim_1_m", "f_core_rear", "flange_scale", "flange_cool_gain",
    "flange_cool_tau_s", "k_core_axial_scale", "C_rear_eff_J_K",
    "G_receiver_rear_W_K", "G_rear_tube_W_K", "G_rear_cavity_W_K",
    "G_rear_axial_W_K",
)

function print_parameter_limiter_report_v30(p)
    println("parameter,value,lower,upper,near_lower,near_upper")
    for i in eachindex(p)
        lower = lb_full_v30[i]
        upper = ub_full_v30[i]
        span = max(upper - lower, eps(Float64))
        near_lower = (p[i] - lower) / span < 0.02
        near_upper = (upper - p[i]) / span < 0.02
        println(join((PARAMETER_NAMES_v30_POWER_PATH[i], p[i], lower, upper, near_lower, near_upper), ','))
    end
end

function tube_loss_totals_v30(result, p, time, irradiance)
    z_local = result.z_rear_tube .- L
    dx = REAR_TUBE_LENGTH_v30 / length(z_local)
    Ttube = result.rear_tube_temperature[:, end]
    Tcavity = result.cavity_temperature[end]
    cavity = 0.0
    flange = 0.0
    for j in eachindex(z_local)
        G_cavity = rear_tube_cavity_conductance_per_length_v30(z_local[j], p)
        G_flange = rear_tube_flange_conductance_per_length_v30(z_local[j], Ttube[j], p, time, irradiance)
        cavity += G_cavity * dx * (Ttube[j] - Tcavity)
        flange += G_flange * dx * (Ttube[j] - WATER_FLANGE_TEMPERATURE_v30)
    end
    return cavity, flange
end

function rear_reservoir_heat_path_v30(result, p)
    Tcore = result.core_temperature[:, end]
    Tperim = result.perim_temperature[:, end]
    Trear = result.rear_reservoir_temperature[:, end]
    Ttube_in = result.rear_tube_temperature[1, end]
    Tcavity = result.cavity_temperature[end]
    f_core = rear_core_fraction_v30(p)
    G_receiver_rear = receiver_rear_conductance_v30(p)
    core = sum(
        G_receiver_rear * result.rear_contact_weights[i] * f_core *
        (Tcore[i] - Trear[i])
        for i in eachindex(Tcore)
    )
    perim = sum(
        G_receiver_rear * result.rear_contact_weights[i] * (1.0 - f_core) *
        (Tperim[i] - Trear[i])
        for i in eachindex(Tperim)
    )
    tube = rear_tube_conductance_v30(p) * (Trear[end] - Ttube_in)
    cavity = sum(
        rear_cavity_conductance_v30(p) * result.rear_contact_weights[i] *
        (Trear[i] - Tcavity)
        for i in eachindex(Trear)
    )
    return core, perim, tube, cavity
end

function print_case_power_path_v30(simulation_id, p; nodes=15)
    outputs, result, experiment = solve_case_v30(p, simulation_id; nodes=nodes)
    conditions = simulation_conditions[simulation_id]
    irradiance = conditions[Io]
    _, _, _, case_data = experimental_case_v30(simulation_id)
    ambient = observation(case_data, simulation_id, "_Tamb")[end]
    scale = absorbed_power_scale_v30(irradiance, p)
    Qtotal = modeled_absorbed_receiver_power_v30(irradiance, p)
    Qcore = core_absorbed_power_v30(irradiance, p)
    Qperim = perimeter_absorbed_power_v30(irradiance, p)
    Qgas = sum(result.receiver_gas_heat[:, end])
    Qrear_gas = sum(result.rear_tube_gas_heat[:, end])
    Tfront = result.perim_temperature[1, end]
    Qfront = H_FRONT_FIXED_v30 * A_frt * (Tfront - ambient) +
             EPS_FRONT_FIXED_v30 * sigma_sb * A_frt * (Tfront^4 - ambient^4)
    Qcavity_ambient = cavity_loss_to_ambient_v30(result.cavity_temperature[end], ambient)
    Qtube_cavity, Qtube_flange = tube_loss_totals_v30(result, p, result.time[end], irradiance)
    Qcore_rear, Qperim_rear, Qrear_tube, Qrear_cavity = rear_reservoir_heat_path_v30(result, p)
    gate = lamp_off_gate_v30(result.time[end], irradiance, p)
    effective_flange_scale = effective_flange_loss_scale_v30(p, result.time[end], irradiance)
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
        Qtube_cavity,
        Qtube_flange,
        Qcore_rear,
        Qperim_rear,
        Qrear_tube,
        Qrear_cavity,
        result.rear_reservoir_temperature[1, end],
        result.rear_reservoir_temperature[cld(size(result.rear_reservoir_temperature, 1), 2), end],
        result.rear_reservoir_temperature[end, end],
        gate,
        effective_flange_scale,
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
print_parameter_limiter_report_v30(FITTED_PARAMS_v30_POWER_PATH)

println("# power path")
println("case,flow_lpm,irradiance,power_scale,Qtotal_abs_W,Qcore_abs_W,Qperim_abs_W,Qreceiver_gas_W,Qrear_gas_W,Qfront_loss_W,Qcavity_ambient_W,Qtube_cavity_W,Qtube_flange_W,Qcore_rear_W,Qperim_rear_W,Qrear_tube_W,Qrear_cavity_W,Trear_inlet_K,Trear_mid_K,Trear_exit_K,flange_gate,effective_flange_scale,T8_err_K,T12_err_K,T11_err_K,T9_err_K,T10_err_K,T3_err_K,T2_err_K")
for simulation_id in ("E67", "E69", "E71", "E72", "E76", "E77", "E80", "E81")
    print_case_power_path_v30(simulation_id, FITTED_PARAMS_v30_POWER_PATH)
end



