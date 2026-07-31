content = read("1D_v38.jl", String)

# Replace gas_profile_v38! signature and body
old_gas_func = r"function gas_profile_v38!\(Tg, Qgas, hcell, Qgas_rear, hrear,(.*?)\nend\n"s
new_gas_func = \"\"\"function gas_profile_v38!(Tg, Qgas, hcell, Qgas_rear, hrear,
                          Tcore, Tperim, Ttube, time, p, operating, z, dx, z_rear, dx_rear, Tgas_perim, Qgas_perim)
    flow = max(0.0, operating.flow_lpm(time))
    Tin = operating.inlet_temperature(time)
    receiver_nodes = length(Tcore)
    Tg[1] = Tin
    Tgas_perim[1] = Tin

    if flow <= 1e-12
        fill!(Qgas, 0.0)
        fill!(hcell, 0.0)
        fill!(Qgas_rear, 0.0)
        fill!(hrear, 0.0)
        fill!(Qgas_perim, 0.0)
        for i in eachindex(Tcore)
            Tg[i + 1] = Tcore[i]
            Tgas_perim[i + 1] = Tperim[i]
        end
        for j in eachindex(Ttube)
            Tg[receiver_nodes + j + 1] = Ttube[j]
        end
        return nothing
    end

    mdot_total = m_dot_standard_v38(flow)
    phi_act = clamp(p[4], 0.05, 1.0)
    mdot_core = phi_act * mdot_total
    mdot_perim = (1.0 - phi_act) * mdot_total

    Tactive_core = Tin
    Tactive_perim = Tin

    for i in eachindex(Tcore)
        # Core gas
        Tfilm = 0.5 * (Tcore[i] + Tactive_core)
        cp = cpf_f(Tfilm)
        mu = mu_f_f(Tfilm)
        kf = kf_f(Tfilm)
        reynolds_core = mdot_core * Dh / (A_flow * mu)
        prandtl = cp * mu / kf
        
        nusselt = receiver_channel_nusselt_v38(z[i], reynolds_core, prandtl, p)
        hcell[i] = nusselt * kf / Dh
        
        if mdot_core > 1e-12
            UA = hcell[i] * P_exchange * dx
            effectiveness = -expm1(-UA / (mdot_core * cp))
            Tactive_next = Tactive_core + effectiveness * (Tcore[i] - Tactive_core)
            Qgas[i] = mdot_core * cp * (Tactive_next - Tactive_core)
            Tactive_core = Tactive_next
        else
            Qgas[i] = 0.0
            Tactive_core = Tcore[i]
        end
        Tg[i + 1] = Tactive_core

        # Perimeter gas
        Tfilm_p = 0.5 * (Tperim[i] + Tactive_perim)
        cp_p = cpf_f(Tfilm_p)
        if mdot_perim > 1e-12
            h_perim = 20.0
            A_perim = 0.16 * dx
            UA_p = h_perim * A_perim
            effectiveness_p = -expm1(-UA_p / (mdot_perim * cp_p))
            Tactive_next_p = Tactive_perim + effectiveness_p * (Tperim[i] - Tactive_perim)
            Qgas_perim[i] = mdot_perim * cp_p * (Tactive_next_p - Tactive_perim)
            Tactive_perim = Tactive_next_p
        else
            Qgas_perim[i] = 0.0
            Tactive_perim = Tperim[i]
        end
        Tgas_perim[i + 1] = Tactive_perim
    end

    for j in eachindex(Ttube)
        inlet_index = receiver_nodes + j
        Tfilm = 0.5 * (Ttube[j] + Tg[inlet_index])
        cp = cpf_f(Tfilm)
        hrear[j] = rear_tube_htc_v38(Ttube[j], Tg[inlet_index], flow * phi_act, Tin)
        
        if mdot_core > 1e-12
            UA = hrear[j] * REAR_TUBE_EXCHANGE_PERIMETER_v38 * dx_rear
            effectiveness = -expm1(-UA / (mdot_core * cp))
            Tg[inlet_index + 1] = Tg[inlet_index] + effectiveness * (Ttube[j] - Tg[inlet_index])
            Qgas_rear[j] = mdot_core * cp * (Tg[inlet_index + 1] - Tg[inlet_index])
        else
            Tg[inlet_index + 1] = Tg[inlet_index]
            Qgas_rear[j] = 0.0
        end
    end
    return nothing
end
\"\"\"
content = replace(content, old_gas_func => new_gas_func)

# Replace receiver_rhs_v38!
old_rhs = r"function receiver_rhs_v38!\(.*?\nend\n"s
new_rhs = \"\"\"function receiver_rhs_v38!(du, u, time, p, operating, z, dx, core_weights,
                          perim_weights, rear_weights,
                          z_rear, dx_rear,
                          Tg, Qgas, hcell, Qgas_rear, hrear, zg, Tgas_perim, Qgas_perim)
    nodes = length(z)
    rear_nodes = length(z_rear)
    Tcore = view(u, 1:nodes)
    Tperim = view(u, (nodes + 1):(2 * nodes))
    tube_start = 2 * nodes + 1
    tube_stop = 2 * nodes + rear_nodes
    Ttube = view(u, tube_start:tube_stop)
    cavity_index = 2 * nodes + rear_nodes + 1
    Tcavity = u[cavity_index]
    rear_start = cavity_index + 1
    rear_stop = cavity_index + nodes
    Trear = view(u, rear_start:rear_stop)
    T3_sensor_index = rear_stop + 1
    T3_sensor = u[T3_sensor_index]
    
    dTcore = view(du, 1:nodes)
    dTperim = view(du, (nodes + 1):(2 * nodes))
    dTtube = view(du, tube_start:tube_stop)
    dTrear = view(du, rear_start:rear_stop)
    ambient = operating.ambient_temperature(time)

    gas_profile_v38!(Tg, Qgas, hcell, Qgas_rear, hrear,
                    Tcore, Tperim, Ttube, time, p, operating, z, dx, z_rear, dx_rear, Tgas_perim, Qgas_perim)
    fill!(du, 0.0)

    # Core axial solid conduction
    for i in 1:(nodes - 1)
        ki = ks_f(Tcore[i])
        kj = ks_f(Tcore[i + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        Qcond_core =
            core_axial_conduction_scale_v38(p) * kface * A_solid *
            (Tcore[i] - Tcore[i + 1]) / dx
        dTcore[i] -= Qcond_core
        dTcore[i + 1] += Qcond_core

        kwi = ks_f(Tperim[i])
        kwj = ks_f(Tperim[i + 1])
        kwface = 2.0 * kwi * kwj / (kwi + kwj)
        Tface_perim = 0.5 * (Tperim[i] + Tperim[i + 1])
        k_perim = perimeter_axial_conductivity_v38(Tface_perim, p)
        Qcond_perim = (kwface * A_solid * 0.20 + k_perim * A_frt) * (Tperim[i] - Tperim[i + 1]) / dx
        dTperim[i] -= Qcond_perim
        dTperim[i + 1] += Qcond_perim
    end

    irradiance = max(0.0, operating.irradiance(time))
    Q_aperture = irradiance * A_frt
    
    f_recv = max(flux_receiver_fraction_v38(), 1e-6)
    f_spill = flux_spillover_fraction_v38()
    Q_spill_total = Q_aperture * (f_spill / f_recv)
    
    scale = absorbed_power_scale_v38(irradiance, p)
    Qcore_solar = ETA_ABS_FIXED_v38 * scale * Q_aperture
    Qperim_solar = scale * Q_spill_total * perimeter_spillover_capture_v38(p)

    G_cp = core_perimeter_conductance_per_length_v38(p)
    radial_cavity_conductance = RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v38

    for i in 1:nodes
        Qrad_in = G_cp * dx * (Tperim[i] - Tcore[i])
        Qperim_cavity = radial_cavity_conductance * dx * (Tperim[i] - Tcavity)
        
        Qcell_core = Qcore_solar * core_weights[i]
        Qcell_perim = Qperim_solar * perim_weights[i]
        dTcore[i] += Qcell_core - Qgas[i] + Qrad_in
        dTperim[i] += Qcell_perim - Qgas_perim[i] - Qrad_in - Qperim_cavity
        du[cavity_index] += Qperim_cavity
    end

    Qfront_perim = H_FRONT_FIXED_v38 * A_frt * (Tperim[1] - ambient) +
                   EPS_FRONT_FIXED_v38 * sigma_sb * A_frt * (Tperim[1]^4 - ambient^4)
    
    f_core_rear = rear_core_fraction_v38(p)
    G_receiver_rear = receiver_rear_conductance_v38(p)
    G_rear_cavity = rear_cavity_conductance_v38(p)

    dTperim[1] -= Qfront_perim
    for i in 1:nodes
        receiver_contact = G_receiver_rear * rear_weights[i]
        cavity_contact = G_rear_cavity * rear_weights[i]
        Qcore_rear = receiver_contact * f_core_rear * (Tcore[i] - Trear[i])
        Qperim_rear = receiver_contact * (1.0 - f_core_rear) * (Tperim[i] - Trear[i])
        Qrear_cavity = cavity_contact * (Trear[i] - Tcavity)
        dTcore[i] -= Qcore_rear
        dTperim[i] -= Qperim_rear
        dTrear[i] += Qcore_rear + Qperim_rear - Qrear_cavity
        du[cavity_index] += Qrear_cavity
    end

    G_rear_axial = rear_axial_conductance_v38(p)
    for i in 1:(nodes - 1)
        Qrear_axial = G_rear_axial * (Trear[i] - Trear[i + 1])
        dTrear[i] -= Qrear_axial
        dTrear[i + 1] += Qrear_axial
    end

    Qrear_tube = rear_tube_conductance_v38(p) * (Trear[end] - Ttube[1])
    dTrear[end] -= Qrear_tube
    dTtube[1] += Qrear_tube

    for j in 1:(rear_nodes - 1)
        ki = alumina_conductivity_v38(Ttube[j])
        kj = alumina_conductivity_v38(Ttube[j + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        Aface = 0.5 * (rear_tube_solid_area_v38(z_rear[j]) + rear_tube_solid_area_v38(z_rear[j + 1]))
        Qcond = kface * Aface * (Ttube[j] - Ttube[j + 1]) / dx_rear
        dTtube[j] -= Qcond
        dTtube[j + 1] += Qcond
    end

    for j in 1:rear_nodes
        G_cavity = rear_tube_cavity_conductance_per_length_v38(z_rear[j], p)
        G_flange = rear_tube_flange_conductance_per_length_v38(z_rear[j], Ttube[j], p, time, irradiance)
        Q_cavity = G_cavity * dx_rear * (Ttube[j] - Tcavity)
        Q_flange = G_flange * dx_rear * (Ttube[j] - WATER_FLANGE_TEMPERATURE_v38)
        dTtube[j] -= Qgas_rear[j] + Q_cavity + Q_flange
        du[cavity_index] += Q_cavity
    end

    Qcavity_ambient = cavity_loss_to_ambient_v38(Tcavity, ambient)
    du[cavity_index] -= Qcavity_ambient

    cell_volume = A_solid * dx
    perim_capacity_cell = perimeter_heat_capacity_total_v38(p) / nodes
    for i in 1:nodes
        core_capacity = rho_s * Cps_f(Tcore[i]) * cell_volume
        dTcore[i] /= max(core_capacity, eps(Float64))
        dTperim[i] /= max(perim_capacity_cell, eps(Float64))
    end
    for j in 1:rear_nodes
        dTtube[j] /= max(rear_tube_capacity_v38(Ttube[j], z_rear[j], dx_rear), eps(Float64))
    end
    du[cavity_index] /= max(CAVITY_HEAT_CAPACITY_v38, eps(Float64))
    C_rear_total = rear_reservoir_heat_capacity_v38(p)
    for i in 1:nodes
        dTrear[i] /= max(C_rear_total * rear_weights[i], eps(Float64))
    end
    
    Tgas_sensor_core = interpolate_state_at_position_v38(Tg, zg, T3_SAMPLE_POSITION_v38)
    
    flow = max(0.0, operating.flow_lpm(time))
    mdot_total = m_dot_standard_v38(flow)
    phi_act = clamp(p[4], 0.05, 1.0)
    mdot_core = phi_act * mdot_total
    mdot_perim = (1.0 - phi_act) * mdot_total

    T_mixed = mdot_total > 1e-12 ? (mdot_core * Tgas_sensor_core + mdot_perim * Tgas_perim[end]) / mdot_total : Tgas_sensor_core
    
    Twall_sensor = interpolate_state_at_position_v38(Ttube, z_rear, T3_SAMPLE_POSITION_v38 - L)
    
    Q_conv_sensor = T3_SENSOR_H_v38 * T3_SENSOR_AREA_v38 * (T_mixed - T3_sensor)
    Q_rad_sensor = T3_SENSOR_EMISSIVITY_v38 * sigma_sb * T3_SENSOR_AREA_v38 * (Twall_sensor^4 - T3_sensor^4)
    du[T3_sensor_index] = (Q_conv_sensor + Q_rad_sensor) / max(T3_SENSOR_CAPACITY_v38, eps(Float64))

    return nothing
end
\"\"\"
content = replace(content, old_rhs => new_rhs)


old_ode = r"function receiver_ode_v38!\(dTs, Ts, context, time\)\n    receiver_rhs_v38!\(\n        dTs, Ts, time, context\.p, context\.operating,\n        context\.z, context\.dx, context\.core_weights, context\.perim_weights,\n        context\.rear_weights,\n        context\.z_rear, context\.dx_rear,\n        context\.Tg, context\.Qgas, context\.hcell, context\.Qgas_rear, context\.hrear, context\.zg\n    \)\n    return nothing\nend"s
new_ode = \"\"\"function receiver_ode_v38!(dTs, Ts, context, time)
    receiver_rhs_v38!(
        dTs, Ts, time, context.p, context.operating,
        context.z, context.dx, context.core_weights, context.perim_weights,
        context.rear_weights,
        context.z_rear, context.dx_rear,
        context.Tg, context.Qgas, context.hcell, context.Qgas_rear, context.hrear, context.zg, context.Tgas_perim, context.Qgas_perim
    )
    return nothing
end\"\"\"
content = replace(content, old_ode => new_ode)

# In simulate_v38, we need to add Tgas_perim, Qgas_perim
old_simulate_alloc = "    Tg = zeros(nodes + rear_nodes + 1)\n    Qgas = zeros(nodes)\n    hcell = zeros(nodes)\n    Qgas_rear = zeros(rear_nodes)\n    hrear = zeros(rear_nodes)\n\n    context = (\n        p=p, operating=operating, z=z, dx=dx,\n        core_weights=core_weights, perim_weights=perim_weights,\n        rear_weights=rear_weights,\n        z_rear=z_rear, dx_rear=dx_rear,\n        Tg=Tg, Qgas=Qgas, hcell=hcell, Qgas_rear=Qgas_rear, hrear=hrear, zg=zg,\n    )"

new_simulate_alloc = \"\"\"    Tg = zeros(nodes + rear_nodes + 1)
    Tgas_perim = zeros(nodes + 1)
    Qgas = zeros(nodes)
    Qgas_perim = zeros(nodes)
    hcell = zeros(nodes)
    Qgas_rear = zeros(rear_nodes)
    hrear = zeros(rear_nodes)

    context = (
        p=p, operating=operating, z=z, dx=dx,
        core_weights=core_weights, perim_weights=perim_weights,
        rear_weights=rear_weights,
        z_rear=z_rear, dx_rear=dx_rear,
        Tg=Tg, Qgas=Qgas, hcell=hcell, Qgas_rear=Qgas_rear, hrear=hrear, zg=zg,
        Tgas_perim=Tgas_perim, Qgas_perim=Qgas_perim
    )\"\"\"
content = replace(content, old_simulate_alloc => new_simulate_alloc)

# inside simulate_v38 loop
old_sim_loop = r"        Tcore_now = view\(core_temperature_history, :, output_index\)\n        Ttube_now = view\(rear_temperature_history, :, output_index\)\n        gas_profile_v38!\(Tg, Qgas, hcell, Qgas_rear, hrear, Tcore_now, Ttube_now,\n                         t, p, operating, z, dx, z_rear, dx_rear\)"s
new_sim_loop = \"\"\"        Tcore_now = view(core_temperature_history, :, output_index)
        Tperim_now = view(perim_temperature_history, :, output_index)
        Ttube_now = view(rear_temperature_history, :, output_index)
        gas_profile_v38!(Tg, Qgas, hcell, Qgas_rear, hrear, Tcore_now, Tperim_now, Ttube_now,
                         t, p, operating, z, dx, z_rear, dx_rear, context.Tgas_perim, context.Qgas_perim)\"\"\"
content = replace(content, old_sim_loop => new_sim_loop)

write("1D_v38.jl", content)
