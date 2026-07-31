import re

with open('1D_v38.jl', 'r') as f:
    content = f.read()

# Make changes to gas_profile_v38!
new_gas_profile = '''function gas_profile_v38!(Tg, Qgas, hcell, Qgas_rear, hrear,
                          Tcore, Tperim, Ttube, time, p, operating, z, dx, z_rear, dx_rear, Tgas_perim)
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
        Tg[1] = Tin
        Tgas_perim[1] = Tin
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
        h_perim = 20.0
        Tfilm_p = 0.5 * (Tperim[i] + Tactive_perim)
        cp_p = cpf_f(Tfilm_p)
        if mdot_perim > 1e-12
            # assume some perimeter exchange area, maybe matching cavity or some fixed value
            # Actually, user said: "integrate Tgas_perim across the Tperim nodes. The perimeter gas is heated by the perimeter solid. Use a simple constant h_perim = 20.0 W/m2K"
            # Let's use P_exchange or similar area? The instruction didn't specify the area. 
            # I will use CAVITY_OUTER_AREA_v38/nodes, or P_exchange? Let's just use a reasonable area.
            # Actually let's use the outer perimeter of the core: 4 * w_t = 0.16 m?
            # Or just use the same P_exchange as the core? No, core P_exchange is large.
            # "perimeter gas exchange" -> A_perim = 4 * w_t * dx = 4 * 0.04 * dx = 0.16 * dx
            A_perim = 0.16 * dx
            UA_p = h_perim * A_perim
            effectiveness_p = -expm1(-UA_p / (mdot_perim * cp_p))
            Tactive_next_p = Tactive_perim + effectiveness_p * (Tperim[i] - Tactive_perim)
            Qgas_perim_cell = mdot_perim * cp_p * (Tactive_next_p - Tactive_perim)
            Tactive_perim = Tactive_next_p
        else
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
'''

# Find the old gas_profile_v38!
content = re.sub(r'function gas_profile_v38!\(.*?\nend\n', new_gas_profile, content, flags=re.DOTALL)

# Now modify receiver_rhs_v38!
# Add Tgas_perim to context and change gas_profile_v38! call
# also we need to extract Qgas_perim_cell to subtract from dTperim.
# But wait, Qgas_perim_cell is not returned. I should make gas_profile_v38! return or populate Qgas_perim.
'''

with open('rewrite.py', 'w') as f:
    f.write()
