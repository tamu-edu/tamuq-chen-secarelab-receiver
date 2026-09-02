# ============================================================================
# 1D_v46.jl - Energy-Accounting Conservative 2-Zone Core/Perimeter Macro Model
# ============================================================================
# v46 addresses all critique points from the theory audit:
#   1. Geometry and polynomial thermophysical properties synchronized with 1D_v3.jl.
#   2. Complete removal of bypass parameter phi_act (phi_act = 1.0) ensuring strict
#      mass and enthalpy conservation throughout receiver channels and rear tube.
#   3. Front aperture suction heat is coupled directly into the incoming gas
#      stream at z = 0, preheating air before channel entry.
#   4. Perimeter optical spillage (1 - chi)*Q_delivered is absorbed at the front
#      face (z = 0, cell 1) with axial heat spreading via aluminum casing conduction.
#   5. Removal of lamp-off flange cooling multipliers (g_cool, tau_cool) for
#      phase-invariant physical conductances.
#   6. Exact energy conservation audit with explicit residual calculation.
#
# Parameter Vector p[1:23] in v46:
#   p[1]  A_Nu          laminar developing Nu prefactor
#   p[2]  B_Re          laminar developing Reynolds exponent (bounded 0.0 to 0.6)
#   p[3]  C_Pr          fixed Prandtl exponent (= 1/3)
#   p[4]  front_dep     front-cell optical deposition fraction in core
#   p[5]  scale_456     delivered power scale M at 456 kW/m2
#   p[6]  scale_304     delivered power scale M at 304 kW/m2
#   p[7]  scale_256     delivered power scale M at 256 kW/m2
#   p[8]  G_core_perim  radial core-perimeter conductance per length [W/m/K]
#   p[9]  C_perim_eff   perimeter/housing participating heat capacity [J/K]
#   p[10] k_perim_ref   perimeter axial conductivity at 900 K [W/m/K]
#   p[11] beta_opt      core Beer-Lambert optical extinction coefficient [1/m]
#   p[12] chi           core solar source fraction (spillover to perimeter = 1 - chi)
#   p[13] f_core_rear   fraction of receiver-to-rear coupling from core
#   p[14] flange_scale  multiplier for rear-tube-to-water-flange base conductance
#   p[15] k_core_axial_scale effective core axial solid-conduction multiplier
#   p[16] C_rear_eff    rear/adaptor rail participating heat capacity [J/K]
#   p[17] G_recv_rear   receiver-to-rear-rail conductance [W/K]
#   p[18] G_rear_tube   rear-rail exit to tube-inlet conductance [W/K]
#   p[19] G_rear_cavity rear-rail to cavity conductance [W/K]
#   p[20] G_rear_axial  rear-rail axial conductance per link [W/K]
#   p[21] delta_web     intra-webbing solid conduction resistance length [m]
#   p[22] C_z           axial spatial decay exponent for Nu [-]
#   p[23] h_suction     aperture suction convection scaling parameter [W/m2/K]
# ============================================================================

include("1D_v5.jl")

begin # v46 fixed constants and synchronized geometry
    const EPS_FRONT_FIXED_v46 = 0.95
    const ETA_ABS_FIXED_v46 = 1.0
    const PRANDTL_EXPONENT_FIXED_v46 = 1.0 / 3.0
    const STANDARD_AIR_DENSITY_v46 = 101325.0 / (287.05 * 294.25)
    const REYNOLDS_REFERENCE_v46 = 50.0
    const NU_FLOOR_v46 = 0.0
    const MEASURED_ASSEMBLY_CAPACITANCE_v46 = 301.0
    const IRRADIANCE_LEVELS_v46 = (456000.0, 304000.0, 256000.0)
    const T3_SAMPLE_POSITION_v46 = 140.0e-3
    const COOLING_ROOM_TEMPERATURE_v46 = 17.0 + 273.15
    const CAPACITANCE_REGULARIZATION_WEIGHT_v46 = 0.00
    const FLUX_SIGMA_FIXED_v46 = 30.0e-3
    const FLUX_INTEGRATION_POINTS_v46 = 161

    # Cavity and housing geometry (synchronized with CAD & 1D_v3.jl)
    const CAVITY_OUTER_RADIUS_v46 = 75.0e-3
    const METAL_THICKNESS_v46 = 18.0e-3
    const INSULATION_OUTER_RADIUS_v46 = CAVITY_OUTER_RADIUS_v46 - METAL_THICKNESS_v46
    const METAL_OUTER_RADIUS_v46 = CAVITY_OUTER_RADIUS_v46
    const CAVITY_LENGTH_v46 = 165.0e-3

    const ADAPTOR_DIAMETER_v46 = 77.6e-3
    const ADAPTOR_RADIUS_v46 = ADAPTOR_DIAMETER_v46 / 2.0
    const ADAPTOR_LENGTH_v46 = 57.0e-3
    const ADAPTOR_TUBE_RADIUS_v46 = 6.5e-3
    const ADAPTOR_RECEIVER_OVERLAP_LENGTH_v46 = 28.0e-3
    const ADAPTOR_TUBE_OVERLAP_LENGTH_v46 = ADAPTOR_LENGTH_v46 - ADAPTOR_RECEIVER_OVERLAP_LENGTH_v46
    const ADAPTOR_OVERLAP_LENGTH_v46 = ADAPTOR_RECEIVER_OVERLAP_LENGTH_v46

    const RECEIVER_EQ_RADIUS_v46 = sqrt(A_frt / pi)
    const RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v46 =
        4.0 * pi * 0.08 / log(INSULATION_OUTER_RADIUS_v46 / RECEIVER_EQ_RADIUS_v46)
    const ADAPTOR_CONTACT_RESISTANCE_v46 =
        1.0 / (100.0 * 4.0 * w_t * ADAPTOR_OVERLAP_LENGTH_v46)
    const ADAPTOR_CONTACT_CONDUCTANCE_REFERENCE_v46 = 1.0 / ADAPTOR_CONTACT_RESISTANCE_v46

    # Rear tube geometry
    const REAR_TUBE_LENGTH_v46 = 150.0e-3
    const REAR_TUBE_CAVITY_LENGTH_v46 = 46.0e-3
    const REAR_TUBE_DEFAULT_NODES_v46 = 12
    const REAR_TUBE_GAS_RADIUS_v46 = ADAPTOR_TUBE_RADIUS_v46
    const REAR_TUBE_WALL_THICKNESS_v46 = 1.5e-3
    const REAR_TUBE_OUTER_RADIUS_v46 = REAR_TUBE_GAS_RADIUS_v46 + REAR_TUBE_WALL_THICKNESS_v46
    const REAR_TUBE_FLOW_AREA_v46 = pi * REAR_TUBE_GAS_RADIUS_v46^2
    const REAR_TUBE_HYDRAULIC_DIAMETER_v46 = 2.0 * REAR_TUBE_GAS_RADIUS_v46
    const REAR_TUBE_EXCHANGE_PERIMETER_v46 = 2.0 * pi * REAR_TUBE_GAS_RADIUS_v46

    # T3 Sensor properties
    const T3_SENSOR_CAPACITY_v46 = 0.05 # J/K (small thermocouple bead)
    const T3_SENSOR_AREA_v46 = 4.0 * pi * (0.0015)^2 # 1.5mm radius bead
    const T3_SENSOR_EMISSIVITY_v46 = 0.80
    const T3_SENSOR_H_v46 = 150.0 # W/m2/K (cross-flow cylinder)

    const ALUMINA_DENSITY_v46 = 3900.0
    const ALUMINUM_DENSITY_v46 = 2700.0
    const ALUMINUM_HEAT_CAPACITY_v46 = 900.0
    const ALUMINUM_CONDUCTIVITY_v46 = 205.0
    const METAL_EMISSIVITY_v46 = 0.20
    const H_NAT_CAVITY_v46 = 10.0
    const INSULATION_DENSITY_v46 = 140.0
    const INSULATION_HEAT_CAPACITY_v46 = 1360.0
    const WATER_FLANGE_TEMPERATURE_v46 = 25.0 + 273.15

    const CAVITY_FELT_VOLUME_v46 = max(
        0.0,
        pi * INSULATION_OUTER_RADIUS_v46^2 * CAVITY_LENGTH_v46 -
        A_frt * L - pi * ADAPTOR_RADIUS_v46^2 * ADAPTOR_LENGTH_v46,
    )
    const CAVITY_METAL_VOLUME_v46 =
        pi * (METAL_OUTER_RADIUS_v46^2 - INSULATION_OUTER_RADIUS_v46^2) * CAVITY_LENGTH_v46 +
        pi * METAL_OUTER_RADIUS_v46^2 * METAL_THICKNESS_v46
    const CAVITY_OUTER_AREA_v46 =
        2.0 * pi * METAL_OUTER_RADIUS_v46 * CAVITY_LENGTH_v46 + pi * METAL_OUTER_RADIUS_v46^2
    const CAVITY_HEAT_CAPACITY_v46 =
        INSULATION_DENSITY_v46 * INSULATION_HEAT_CAPACITY_v46 * CAVITY_FELT_VOLUME_v46 +
        ALUMINUM_DENSITY_v46 * ALUMINUM_HEAT_CAPACITY_v46 * CAVITY_METAL_VOLUME_v46
end

begin # v46 flux spatial partition
    transverse_flux_shape_v46(x, y; sigma=FLUX_SIGMA_FIXED_v46) =
        exp(-0.5 * (x^2 + y^2) / sigma^2)

    function flux_zone_fractions_v46(; radius=CAVITY_OUTER_RADIUS_v46,
                                     receiver_half_width=w_t / 2.0,
                                     points=FLUX_INTEGRATION_POINTS_v46)
        dx_flux = 2.0 * radius / (points - 1)
        receiver_power = 0.0
        aperture_power = 0.0
        for ix in 1:points
            x = -radius + (ix - 1) * dx_flux
            for iy in 1:points
                y = -radius + (iy - 1) * dx_flux
                r2 = x^2 + y^2
                r2 > radius^2 && continue
                flux = transverse_flux_shape_v46(x, y)
                aperture_power += flux
                if abs(x) <= receiver_half_width && abs(y) <= receiver_half_width
                    receiver_power += flux
                end
            end
        end
        f_rec = receiver_power / max(aperture_power, eps(Float64))
        return (receiver_fraction=f_rec, spillover_fraction=max(0.0, 1.0 - f_rec))
    end

    const FLUX_SPLIT_v46 = flux_zone_fractions_v46()
    flux_receiver_fraction_v46() = FLUX_SPLIT_v46.receiver_fraction
    flux_spillover_fraction_v46() = FLUX_SPLIT_v46.spillover_fraction
end

# Default v46 parameter vector p[1:23]
begin # v46 parameter values and bounds
    pnew_v46 = [
        2.50,               # p[1]  A_Nu (laminar developing prefactor)
        0.50,               # p[2]  B_Re (Reynolds exponent bounded 0.0 to 0.6)
        PRANDTL_EXPONENT_FIXED_v46, # p[3]  C_Pr (= 1/3)
        0.35,               # p[4]  front_dep (front-cell core deposition)
        1.34,               # p[5]  scale_456
        1.58,               # p[6]  scale_304
        1.11,               # p[7]  scale_256
        6.00,               # p[8]  G_core_perim (radial core-perimeter conductance [W/m/K])
        108.0,              # p[9]  C_perim_eff (perimeter participating capacity [J/K])
        15.0,               # p[10] k_perim_ref (perimeter axial conductivity [W/m/K])
        220.0,              # p[11] beta_opt [1/m]
        0.45,               # p[12] chi (core source fraction; 1 - chi is perimeter spill)
        0.99,               # p[13] f_core_rear (receiver-rear coupling core fraction)
        0.10,               # p[14] flange_scale (rear/flange base conductance multiplier)
        0.15,               # p[15] k_core_axial_scale
        90.0,               # p[16] C_rear_eff (rear/adaptor participating capacity [J/K])
        0.15,               # p[17] G_recv_rear (receiver-to-rear conductance [W/K])
        3.50,               # p[18] G_rear_tube (rear rail to tube conductance [W/K])
        0.15,               # p[19] G_rear_cavity (rear rail to cavity conductance [W/K])
        9.00,               # p[20] G_rear_axial (rear rail axial conductance per link [W/K])
        6.0e-5,             # p[21] delta_web [m]
        0.80,               # p[22] C_z (Nusselt axial decay exponent)
        250.0,              # p[23] h_suction (front-face suction convection scale [W/m2/K])
    ]

    lb_full_v46 = [
        0.01,   # p[1] A_Nu
        0.0,    # p[2] B_Re
        pnew_v46[3], # p[3] C_Pr
        0.0,    # p[4] front_dep
        0.5,    # p[5] scale_456
        0.5,    # p[6] scale_304
        0.5,    # p[7] scale_256
        2.0,    # p[8] G_core_perim
        105.0,  # p[9] C_perim_eff
        0.0,    # p[10] k_perim_ref
        20.0,   # p[11] beta_opt
        0.0,    # p[12] chi
        0.0,    # p[13] f_core_rear
        0.05,   # p[14] flange_scale
        0.0,    # p[15] k_core_axial_scale
        50.0,   # p[16] C_rear_eff
        0.0,    # p[17] G_recv_rear
        0.0,    # p[18] G_rear_tube
        0.0,    # p[19] G_rear_cavity
        0.0,    # p[20] G_rear_axial
        0.0,    # p[21] delta_web
        0.0,    # p[22] C_z
        0.0,    # p[23] h_suction
    ]

    ub_full_v46 = [
        25.00,  # p[1] A_Nu
        0.60,   # p[2] B_Re
        pnew_v46[3], # p[3] C_Pr
        1.0,    # p[4] front_dep
        2.0,    # p[5] scale_456
        2.0,    # p[6] scale_304
        2.0,    # p[7] scale_256
        80.0,   # p[8] G_core_perim
        230.0,  # p[9] C_perim_eff
        2000.0, # p[10] k_perim_ref
        500.0,  # p[11] beta_opt
        1.0,    # p[12] chi
        1.0,    # p[13] f_core_rear
        20.0,   # p[14] flange_scale
        0.50,   # p[15] k_core_axial_scale
        150.0,  # p[16] C_rear_eff
        50.0,   # p[17] G_recv_rear
        50.0,   # p[18] G_rear_tube
        10.0,   # p[19] G_rear_cavity
        100.0,  # p[20] G_rear_axial
        2.0e-3, # p[21] delta_web
        1.0,    # p[22] C_z
        1000.0, # p[23] h_suction
    ]

    fit_rear_stage_indices_v46 = [14, 17, 18, 19, 20]
    fit_transport_stage_indices_v46 = [8, 10, 14, 15, 17, 18, 19, 20, 21, 22, 23]
    fit_full_stage_indices_v46 = [1, 2, 4, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
    fit_heat_transfer_indices_v46 = [1, 2, 21, 22, 23]
    fit_source_indices_v46 = [4, 11, 12, 13, 14, 15, 17, 18, 19, 20]
    fit_power_scale_indices_v46 = [5, 6, 7]
end

function irradiance_level_index_v46(irradiance)
    irradiance <= 0.0 && return 0
    differences = abs.(collect(IRRADIANCE_LEVELS_v46) .- Float64(irradiance))
    return argmin(differences)
end

power_scale_index_v46(irradiance) = begin
    level_index = irradiance_level_index_v46(irradiance)
    level_index == 0 ? 0 : 4 + level_index
end

function absorbed_power_scale_v46(irradiance, p)
    index = power_scale_index_v46(irradiance)
    return index == 0 ? 1.0 : p[index]
end

core_source_fraction_v46(p) = clamp(p[12], 0.0, 1.0)
perimeter_source_fraction_v46(p) = 1.0 - core_source_fraction_v46(p)

core_perimeter_conductance_per_length_v46(p) = max(p[8], 0.0)
perimeter_heat_capacity_total_v46(p) = max(p[9], 105.0)

perimeter_axial_conductivity_v46(T, p) = max(p[10], 0.0)
core_axial_conduction_scale_v46(p) = clamp(p[15], 0.0, 1.0)

rear_core_fraction_v46(p) = clamp(p[13], 0.0, 1.0)
rear_reservoir_heat_capacity_v46(p) = max(p[16], 50.0)
receiver_rear_conductance_v46(p) = max(p[17], 0.0)
rear_tube_conductance_v46(p) = max(p[18], 0.0)
rear_cavity_conductance_v46(p) = max(p[19], 0.0)
rear_axial_conductance_v46(p) = max(p[20], 0.0)
flange_loss_scale_v46(p) = max(p[14], 0.0)

# Rear contact spatial weights (normalized affine profile)
function rear_contact_weights_v46(z; L_recv=L)
    s = z ./ L_recv
    raw = 0.02 .+ 0.98 .* s
    return raw ./ sum(raw)
end

function receiver_channel_nusselt_v46(z, reynolds, prandtl, p)
    A_nu = max(p[1], 0.0)
    B_re = clamp(p[2], 0.0, 0.6)
    C_pr = p[3]
    C_z = clamp(p[22], 0.0, 1.0)
    entry_factor = (Dh / max(z, Dh / 2.0))^C_z
    Nu_dev = A_nu * max(reynolds, eps(Float64))^B_re * max(prandtl, eps(Float64))^C_pr * entry_factor
    return max(NU_FLOOR_v46, Nu_dev)
end

function front_heat_transfer_coefficient_v46(flow_lpm, Tin, p)
    flow_lpm <= 1e-12 && return 10.0
    mdot_total = m_dot_standard_v46(flow_lpm)
    mu_in = mu_f_f(Tin)
    reynolds_inlet = mdot_total * Dh / (A_flow * mu_in)
    h_suction_scale = max(p[23], 0.0)
    h_front = 10.0 + h_suction_scale * sqrt(max(reynolds_inlet, 0.0) / REYNOLDS_REFERENCE_v46)
    return h_front
end

m_dot_standard_v46(flow_lpm) =
    STANDARD_AIR_DENSITY_v46 * max(0.0, flow_lpm) / 60000.0

alumina_conductivity_v46(T) =
    5.5 + 34.5 * exp(-0.0033 * (property_temperature(T) - 273.15))

function rear_tube_solid_area_v46(z_rear)
    if z_rear <= ADAPTOR_TUBE_OVERLAP_LENGTH_v46
        return pi * (ADAPTOR_RADIUS_v46^2 - REAR_TUBE_GAS_RADIUS_v46^2)
    end
    return pi * (REAR_TUBE_OUTER_RADIUS_v46^2 - REAR_TUBE_GAS_RADIUS_v46^2)
end

function rear_tube_capacity_v46(T, z_rear, dx_rear)
    cp = max(500.0, (1.00446 + 1.742e-4 * property_temperature(T) -
                     2.796e4 * property_temperature(T)^(-2)) * 1000.0)
    return ALUMINA_DENSITY_v46 * cp * rear_tube_solid_area_v46(z_rear) * dx_rear
end

function rear_tube_cavity_conductance_per_length_v46(z_rear, p)
    z_rear > REAR_TUBE_CAVITY_LENGTH_v46 && return 0.0
    k_insulation = 0.08
    r_outer = INSULATION_OUTER_RADIUS_v46
    r_inner = z_rear <= ADAPTOR_TUBE_OVERLAP_LENGTH_v46 ? ADAPTOR_RADIUS_v46 : REAR_TUBE_OUTER_RADIUS_v46
    return 2.0 * pi * k_insulation / max(log(r_outer / r_inner), 0.1)
end

function rear_tube_flange_conductance_per_length_v46(z_rear, T, p)
    z_rear < REAR_TUBE_CAVITY_LENGTH_v46 && return 0.0
    scale = flange_loss_scale_v46(p)
    k_alumina = alumina_conductivity_v46(T)
    L_to_flange = max(REAR_TUBE_LENGTH_v46 - z_rear, 1e-3)
    A_solid = rear_tube_solid_area_v46(z_rear)
    G_solid = k_alumina * A_solid / L_to_flange
    h_flange = 100.0 * scale
    A_contact = 2.0 * pi * REAR_TUBE_OUTER_RADIUS_v46
    G_contact = h_flange * A_contact
    return 1.0 / (1.0 / max(G_solid, 1e-4) + 1.0 / max(G_contact, 1e-4))
end

function cavity_loss_to_ambient_v46(Tcavity, Tambient)
    Qconv = H_NAT_CAVITY_v46 * CAVITY_OUTER_AREA_v46 * (Tcavity - Tambient)
    Qrad = METAL_EMISSIVITY_v46 * sigma_sb * CAVITY_OUTER_AREA_v46 * (Tcavity^4 - Tambient^4)
    return Qconv + Qrad
end

function rear_tube_htc_v46(Twall, Tgas, flow, Tin)
    flow <= 1e-12 && return 0.0
    mdot = m_dot_standard_v46(flow)
    Tfilm = 0.5 * (Twall + Tgas)
    mu = mu_f_f(Tfilm)
    cp = cpf_f(Tfilm)
    kf = kf_f(Tfilm)
    re = mdot * REAR_TUBE_HYDRAULIC_DIAMETER_v46 / (REAR_TUBE_FLOW_AREA_v46 * mu)
    pr = cp * mu / kf
    nu = 3.66 + 0.0668 * (REAR_TUBE_HYDRAULIC_DIAMETER_v46 / REAR_TUBE_LENGTH_v46) * re * pr /
                (1.0 + 0.04 * ((REAR_TUBE_HYDRAULIC_DIAMETER_v46 / REAR_TUBE_LENGTH_v46) * re * pr)^(2.0 / 3.0))
    return nu * kf / REAR_TUBE_HYDRAULIC_DIAMETER_v46
end

function interpolate_state_at_position_v46(profile, grid, position)
    position <= grid[1] && return profile[1]
    position >= grid[end] && return profile[end]
    idx = searchsortedlast(grid, position)
    f = (position - grid[idx]) / (grid[idx + 1] - grid[idx])
    return (1.0 - f) * profile[idx] + f * profile[idx + 1]
end

# ============================================================================
# GAS MARCHING & ENTHALPY BALANCE (v46: 100% conservative, phi_act = 1.0)
# ============================================================================
function gas_profile_v46!(Tg, Qgas, hcell, Qgas_rear, hrear,
                          Tcore, Tperim, Ttube, time, p, operating,
                          z, dx, z_rear, dx_rear)
    receiver_nodes = length(Tcore)
    flow = max(0.0, operating.flow_lpm(time))
    Tin = operating.inlet_temperature(time)
    Tg[1] = Tin

    if flow <= 1e-12
        fill!(Qgas, 0.0)
        fill!(hcell, 0.0)
        fill!(Qgas_rear, 0.0)
        fill!(hrear, 0.0)
        Tg[1] = Tin
        for i in eachindex(Tcore)
            Tg[i + 1] = Tcore[i]
        end
        for j in eachindex(Ttube)
            Tg[receiver_nodes + j + 1] = Ttube[j]
        end
        return nothing
    end

    mdot_total = m_dot_standard_v46(flow)
    
    # 1. Front suction preheating across full frontal face
    h_front = front_heat_transfer_coefficient_v46(flow, Tin, p)
    cp_in = cpf_f(Tin)
    effectiveness_suction = -expm1(-h_front * A_frt / (mdot_total * cp_in))
    Tg[1] = Tin + effectiveness_suction * (Tperim[1] - Tin)

    # 2. Receiver channel flow (total mass flow through 100 honeycomb channels)
    Tgas_curr = Tg[1]
    for i in eachindex(Tcore)
        Tfilm = 0.5 * (Tcore[i] + Tgas_curr)
        cp = cpf_f(Tfilm)
        mu = mu_f_f(Tfilm)
        kf = kf_f(Tfilm)
        reynolds = mdot_total * Dh / (A_flow * mu)
        prandtl = cp * mu / kf

        nusselt = receiver_channel_nusselt_v46(z[i], reynolds, prandtl, p)
        hcell[i] = nusselt * kf / Dh
        delta_web = max(p[21], 0.0)
        k_solid = max(ks_f(Tcore[i]), eps(Float64))
        U_eff = 1.0 / (1.0 / max(hcell[i], eps(Float64)) + delta_web / k_solid)
        UA = U_eff * P_exchange * dx
        effectiveness = -expm1(-UA / (mdot_total * cp))

        Tgas_next = Tgas_curr + effectiveness * (Tcore[i] - Tgas_curr)
        Qgas[i] = mdot_total * cp * (Tgas_next - Tgas_curr)
        Tgas_curr = Tgas_next
        Tg[i + 1] = Tgas_curr
    end

    # 3. Rear exit tube flow (conservative total mass flow entry at Tgas(L))
    for j in eachindex(Ttube)
        inlet_index = receiver_nodes + j
        Tfilm = 0.5 * (Ttube[j] + Tg[inlet_index])
        cp = cpf_f(Tfilm)
        hrear[j] = rear_tube_htc_v46(Ttube[j], Tg[inlet_index], flow, Tin)
        UA = hrear[j] * REAR_TUBE_EXCHANGE_PERIMETER_v46 * dx_rear
        effectiveness = -expm1(-UA / (mdot_total * cp))
        Tg[inlet_index + 1] = Tg[inlet_index] + effectiveness * (Ttube[j] - Tg[inlet_index])
        Qgas_rear[j] = mdot_total * cp * (Tg[inlet_index + 1] - Tg[inlet_index])
    end
    return nothing
end

function receiver_rhs_v46!(du, u, time, p, operating, z, dx, core_weights,
                          rear_weights, z_rear, dx_rear,
                          Tg, Qgas, hcell, Qgas_rear, hrear, zg)
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

    fill!(du, 0.0)

    # 1. Evaluate gas profile
    gas_profile_v46!(
        Tg, Qgas, hcell, Qgas_rear, hrear,
        Tcore, Tperim, Ttube, time, p, operating,
        z, dx, z_rear, dx_rear
    )

    # 2. Solid axial conduction inside core
    k_core_scale = core_axial_conduction_scale_v46(p)
    for i in 1:(nodes - 1)
        ki = ks_f(Tcore[i]) * k_core_scale
        kj = ks_f(Tcore[i + 1]) * k_core_scale
        kface = 2.0 * ki * kj / (ki + kj)
        Qcond = kface * A_solid * (Tcore[i] - Tcore[i + 1]) / dx
        dTcore[i] -= Qcond
        dTcore[i + 1] += Qcond
    end

    # 3. Solid axial conduction along perimeter housing (aluminum casing)
    for i in 1:(nodes - 1)
        k_perim = perimeter_axial_conductivity_v46(0.5 * (Tperim[i] + Tperim[i + 1]), p)
        ki = (0.20 * ks_f(Tperim[i]) * A_solid + k_perim * A_frt) / A_frt
        kj = (0.20 * ks_f(Tperim[i + 1]) * A_solid + k_perim * A_frt) / A_frt
        kface = 2.0 * ki * kj / (ki + kj)
        Qcond = kface * A_frt * (Tperim[i] - Tperim[i + 1]) / dx
        dTperim[i] -= Qcond
        dTperim[i + 1] += Qcond
    end

    # 4. Solar radiative source distribution
    irradiance = operating.irradiance(time)
    ambient = operating.ambient_temperature(time)
    Q_aperture = max(0.0, irradiance) * A_frt
    
    M = absorbed_power_scale_v46(irradiance, p)
    Q_delivered = M * Q_aperture
    chi = core_source_fraction_v46(p)
    Qcore_solar = chi * Q_delivered
    Qperim_solar = (1.0 - chi) * Q_delivered

    G_cp = core_perimeter_conductance_per_length_v46(p)
    radial_cavity_conductance = RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v46

    for i in 1:nodes
        Qrad_in = G_cp * dx * (Tperim[i] - Tcore[i])
        Qperim_cavity = radial_cavity_conductance * dx * (Tperim[i] - Tcavity)
        
        Qcell_core = Qcore_solar * core_weights[i]
        # Perimeter spillage is absorbed at the front face (cell 1)
        Qcell_perim = (i == 1) ? Qperim_solar : 0.0

        dTcore[i] += Qcell_core - Qgas[i] + Qrad_in
        dTperim[i] += Qcell_perim - Qrad_in - Qperim_cavity
        du[cavity_index] += Qperim_cavity
    end

    # 5. Front face suction convection & radiative loss
    flow = max(0.0, operating.flow_lpm(time))
    Tin = operating.inlet_temperature(time)
    mdot_total = m_dot_standard_v46(flow)
    Q_suction = mdot_total * cpf_f(Tin) * (Tg[1] - Tin) # Enthalpy transferred into gas
    Qfront_rad = EPS_FRONT_FIXED_v46 * sigma_sb * A_frt * (Tperim[1]^4 - ambient^4)
    dTperim[1] -= (Q_suction + Qfront_rad)

    # 6. Rear contact rail and cavity couplings
    f_core_rear = rear_core_fraction_v46(p)
    G_receiver_rear = receiver_rear_conductance_v46(p)
    G_rear_cavity = rear_cavity_conductance_v46(p)

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

    # 7. Rear rail internal axial conduction
    G_rear_axial = rear_axial_conductance_v46(p)
    for i in 1:(nodes - 1)
        Qrear_axial = G_rear_axial * (Trear[i] - Trear[i + 1])
        dTrear[i] -= Qrear_axial
        dTrear[i + 1] += Qrear_axial
    end

    # 8. Rear rail to exit tube conduction link
    Qrear_tube = rear_tube_conductance_v46(p) * (Trear[end] - Ttube[1])
    dTrear[end] -= Qrear_tube
    dTtube[1] += Qrear_tube

    # 9. Rear tube solid axial conduction
    for j in 1:(rear_nodes - 1)
        ki = alumina_conductivity_v46(Ttube[j])
        kj = alumina_conductivity_v46(Ttube[j + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        Aface = 0.5 * (rear_tube_solid_area_v46(z_rear[j]) + rear_tube_solid_area_v46(z_rear[j + 1]))
        Qcond = kface * Aface * (Ttube[j] - Ttube[j + 1]) / dx_rear
        dTtube[j] -= Qcond
        dTtube[j + 1] += Qcond
    end

    # 10. Rear tube heat loss to cavity and water-cooled mounting flange
    for j in 1:rear_nodes
        G_cavity = rear_tube_cavity_conductance_per_length_v46(z_rear[j], p)
        G_flange = rear_tube_flange_conductance_per_length_v46(z_rear[j], Ttube[j], p)
        Q_cavity = G_cavity * dx_rear * (Ttube[j] - Tcavity)
        Q_flange = G_flange * dx_rear * (Ttube[j] - WATER_FLANGE_TEMPERATURE_v46)
        dTtube[j] -= Qgas_rear[j] + Q_cavity + Q_flange
        du[cavity_index] += Q_cavity
    end

    # 11. Cavity shell loss to ambient
    Qcavity_ambient = cavity_loss_to_ambient_v46(Tcavity, ambient)
    du[cavity_index] -= Qcavity_ambient

    # 12. Convert net heat rates to temperature derivatives dT/dt
    cell_volume = A_solid * dx
    perim_capacity_cell = perimeter_heat_capacity_total_v46(p) / nodes
    for i in 1:nodes
        core_capacity = rho_s * Cps_f(Tcore[i]) * cell_volume
        dTcore[i] /= max(core_capacity, eps(Float64))
        dTperim[i] /= max(perim_capacity_cell, eps(Float64))
    end
    for j in 1:rear_nodes
        dTtube[j] /= max(rear_tube_capacity_v46(Ttube[j], z_rear[j], dx_rear), eps(Float64))
    end
    du[cavity_index] /= max(CAVITY_HEAT_CAPACITY_v46, eps(Float64))
    C_rear_total = rear_reservoir_heat_capacity_v46(p)
    for i in 1:nodes
        dTrear[i] /= max(C_rear_total * rear_weights[i], eps(Float64))
    end
    
    # 13. Centerline T3 Thermocouple Observer at 140 mm
    Tgas_sensor = interpolate_state_at_position_v46(Tg, zg, T3_SAMPLE_POSITION_v46)
    Twall_sensor = interpolate_state_at_position_v46(Ttube, z_rear, T3_SAMPLE_POSITION_v46 - L)
    
    Q_conv_sensor = T3_SENSOR_H_v46 * T3_SENSOR_AREA_v46 * (Tgas_sensor - T3_sensor)
    Q_rad_sensor = T3_SENSOR_EMISSIVITY_v46 * sigma_sb * T3_SENSOR_AREA_v46 * (Twall_sensor^4 - T3_sensor^4)
    du[T3_sensor_index] = (Q_conv_sensor + Q_rad_sensor) / max(T3_SENSOR_CAPACITY_v46, eps(Float64))

    return nothing
end

function receiver_ode_v46!(dTs, Ts, context, time)
    receiver_rhs_v46!(
        dTs, Ts, time, context.p, context.operating,
        context.z, context.dx, context.core_weights,
        context.rear_weights,
        context.z_rear, context.dx_rear,
        context.Tg, context.Qgas, context.hcell, context.Qgas_rear, context.hrear, context.zg
    )
    return nothing
end

function simulate_v46(p, operating, save_times;
                      initial_temperature=Tamb, nodes=default_nodes,
                      rear_nodes=REAR_TUBE_DEFAULT_NODES_v46,
                      initial_perim_temperature=nothing,
                      initial_rear_temperature=nothing,
                      initial_rear_reservoir_temperature=nothing,
                      initial_cavity_temperature=Tamb,
                      solver=Rodas5P(autodiff=AutoFiniteDiff()),
                      reltol=1e-6, abstol=1e-7, dtmax=30.0)
    length(p) == 23 || throw(ArgumentError("1D_v46 expects 23 parameters"))
    all(isfinite, p) || throw(ArgumentError("all model parameters must be finite"))
    all(lb_full_v46[i] <= p[i] <= ub_full_v46[i] for i in eachindex(p)) ||
        throw(ArgumentError("1D_v46 parameters are outside declared bounds"))

    time = Float64.(save_times)
    dx = L / nodes
    dx_rear = REAR_TUBE_LENGTH_v46 / rear_nodes
    z = collect(range(dx / 2, step=dx, length=nodes))
    z_rear = collect(range(dx_rear / 2, step=dx_rear, length=rear_nodes))
    zg = vcat(
        collect(range(0.0, L, length=nodes + 1)),
        L .+ collect(range(dx_rear, step=dx_rear, length=rear_nodes)),
    )
    core_weights = solar_weights_v5(p[11], p[4], nodes)
    rear_weights = rear_contact_weights_v46(z)

    Tcore_initial = initial_profile_v3(initial_temperature, z)
    Tperim_initial = isnothing(initial_perim_temperature) ?
                     copy(Tcore_initial) :
                     initial_profile_v3(initial_perim_temperature, z)
    Ttube_initial = if isnothing(initial_rear_temperature)
        fill(Tcore_initial[end], rear_nodes)
    elseif initial_rear_temperature isa Function
        Float64.([initial_rear_temperature(zr) for zr in z_rear])
    else
        fill(Float64(initial_rear_temperature), rear_nodes)
    end
    Trear_initial = if isnothing(initial_rear_reservoir_temperature)
        rear_core_fraction_v46(p) .* Tcore_initial .+
        (1.0 - rear_core_fraction_v46(p)) .* Tperim_initial
    elseif initial_rear_reservoir_temperature isa Function
        Float64.([initial_rear_reservoir_temperature(zi) for zi in z])
    else
        fill(Float64(initial_rear_reservoir_temperature), nodes)
    end
    u_initial = vcat(
        Tcore_initial, Tperim_initial, Ttube_initial,
        Float64(initial_cavity_temperature), Trear_initial, Tcore_initial[end],
    )

    Tg_history = Matrix{Float64}(undef, nodes + rear_nodes + 1, length(time))
    qgas_history = Matrix{Float64}(undef, nodes, length(time))
    hcell_history = Matrix{Float64}(undef, nodes, length(time))
    qgas_rear_history = Matrix{Float64}(undef, rear_nodes, length(time))
    hrear_history = Matrix{Float64}(undef, rear_nodes, length(time))

    context = (
        p=p,
        operating=operating,
        z=z,
        dx=dx,
        core_weights=core_weights,
        rear_weights=rear_weights,
        z_rear=z_rear,
        dx_rear=dx_rear,
        Tg=zeros(nodes + rear_nodes + 1),
        Qgas=zeros(nodes),
        hcell=zeros(nodes),
        Qgas_rear=zeros(rear_nodes),
        hrear=zeros(rear_nodes),
        zg=zg,
    )

    prob = ODEProblem(receiver_ode_v46!, u_initial, (time[1], time[end]), context)
    solution = solve(
        prob,
        solver;
        saveat=time,
        reltol=reltol,
        abstol=abstol,
        dtmax=dtmax,
        maxiters=1000000,
    )

    solution.t == time || error("solution times do not match requested grid")

    for k in eachindex(time)
        u_k = solution.u[k]
        Tcore_k = view(u_k, 1:nodes)
        Tperim_k = view(u_k, (nodes + 1):(2 * nodes))
        Ttube_k = view(u_k, (2 * nodes + 1):(2 * nodes + rear_nodes))
        gas_profile_v46!(
            context.Tg, context.Qgas, context.hcell,
            context.Qgas_rear, context.hrear,
            Tcore_k, Tperim_k, Ttube_k, time[k], p, operating,
            z, dx, z_rear, dx_rear
        )
        Tg_history[:, k] .= context.Tg
        qgas_history[:, k] .= context.Qgas
        hcell_history[:, k] .= context.hcell
        qgas_rear_history[:, k] .= context.Qgas_rear
        hrear_history[:, k] .= context.hrear
    end

    cavity_index = 2 * nodes + rear_nodes + 1
    rear_start = cavity_index + 1
    rear_stop = cavity_index + nodes
    T3_sensor_index = rear_stop + 1

    return (
        time=time,
        z_solid=z,
        z_gas=zg,
        z_rear=z_rear,
        core_temperature=Array(solution[(1:nodes), :]),
        perim_temperature=Array(solution[((nodes + 1):(2 * nodes)), :]),
        rear_temperature=Array(solution[((2 * nodes + 1):(2 * nodes + rear_nodes)), :]),
        cavity_temperature=vec(solution[cavity_index, :]),
        rear_reservoir_temperature=Array(solution[(rear_start:rear_stop), :]),
        T3_sensor=vec(solution[T3_sensor_index, :]),
        gas_temperature=Tg_history,
        heat_transfer_coefficient=hcell_history,
        rear_heat_transfer_coefficient=hrear_history,
        receiver_gas_heat=qgas_history,
        rear_tube_gas_heat=qgas_rear_history,
        core_source_weights=core_weights,
        rear_contact_weights=rear_weights,
        ode_solution=solution,
    )
end

begin # v46 experimental interface & sensor extraction
    const WALL_CHAIN_POSITIONS_v46 = (
        T8=sensor_positions[:T8],   # 0.011 m (11 mm)
        T12=sensor_positions[:T9],  # 0.058 m (58 mm)
        T11=sensor_positions[:T10], # 0.107 m (107 mm)
    )

    function measured_initial_perim_profile_v46(data, simulation_id)
        id_str = string(simulation_id)
        z8 = WALL_CHAIN_POSITIONS_v46.T8
        z12 = WALL_CHAIN_POSITIONS_v46.T12
        z11 = WALL_CHAIN_POSITIONS_v46.T11
        T8 = observation(data, id_str, "_T8")[1]
        T12 = observation(data, id_str, "_T12")[1]
        T11 = observation(data, id_str, "_T11")[1]
        return function (z)
            z <= z8 && return T8
            z >= z11 && return T11
            if z <= z12
                return T8 + (T12 - T8) * (z - z8) / (z12 - z8)
            end
            return T12 + (T11 - T12) * (z - z12) / (z11 - z12)
        end
    end

    function measured_initial_core_profile_v46(data, simulation_id)
        id_str = string(simulation_id)
        z8 = sensor_positions[:T8]
        z9 = sensor_positions[:T9]
        z10 = sensor_positions[:T10]
        T8 = observation(data, id_str, "_T8")[1]
        T9 = observation(data, id_str, "_T9")[1]
        T10 = observation(data, id_str, "_T10")[1]
        return function (z)
            z <= z8 && return T8
            z >= z10 && return T10
            if z <= z9
                return T8 + (T9 - T8) * (z - z8) / (z9 - z8)
            end
            return T9 + (T10 - T9) * (z - z9) / (z10 - z9)
        end
    end

    function measured_initial_rear_profile_v46(data, simulation_id, p)
        id_str = string(simulation_id)
        core_profile = measured_initial_core_profile_v46(data, id_str)
        perim_profile = measured_initial_perim_profile_v46(data, id_str)
        f_core_rear = rear_core_fraction_v46(p)
        T_receiver_exit = f_core_rear * core_profile(L) + (1.0 - f_core_rear) * perim_profile(L)
        T_tube_exit = observation(data, id_str, "_Tf")[1]
        return function (z_rear)
            fraction = clamp(z_rear / REAR_TUBE_LENGTH_v46, 0.0, 1.0)
            return T_receiver_exit + fraction * (T_tube_exit - T_receiver_exit)
        end
    end

    function measured_initial_rear_reservoir_profile_v46(data, simulation_id, p)
        core_profile = measured_initial_core_profile_v46(data, simulation_id)
        perim_profile = measured_initial_perim_profile_v46(data, simulation_id)
        f_core = rear_core_fraction_v46(p)
        return function (z)
            return f_core * core_profile(z) + (1.0 - f_core) * perim_profile(z)
        end
    end

    function experimental_case_v46(simulation_id; is_cooling=false)
        id_str = string(simulation_id)
        data = is_cooling ? measurements_cooling : measurements
        conditions = is_cooling ? simulation_conditions_cooling : simulation_conditions
        time = observation_time(data, id_str)
        experiment = hcat(
            observation(data, id_str, "_T8"),
            observation(data, id_str, "_T12"),
            observation(data, id_str, "_T11"),
            observation(data, id_str, "_T9"),
            observation(data, id_str, "_T10"),
            observation(data, id_str, "_Tf"),
            observation(data, id_str, "_T2"),
        )
        cond_dict = haskey(conditions, id_str) ? conditions[id_str] : conditions[simulation_id]
        return cond_dict, time, experiment, data
    end

    function core_at_v46(result, position)
        z = result.z_solid
        values = result.core_temperature
        position <= z[1] && return vec(values[1, :])
        position >= z[end] && return vec(values[end, :])
        i = searchsortedlast(z, position)
        fraction = (position - z[i]) / (z[i + 1] - z[i])
        return vec((1.0 - fraction) .* values[i, :] .+ fraction .* values[i + 1, :])
    end

    function perim_at_v46(result, position)
        z = result.z_solid
        values = result.perim_temperature
        position <= z[1] && return vec(values[1, :])
        position >= z[end] && return vec(values[end, :])
        i = searchsortedlast(z, position)
        fraction = (position - z[i]) / (z[i + 1] - z[i])
        return vec((1.0 - fraction) .* values[i, :] .+ fraction .* values[i + 1, :])
    end

    function extract_outputs_v46(result, p)
        T8 = perim_at_v46(result, sensor_positions[:T8])
        T12_perim = perim_at_v46(result, WALL_CHAIN_POSITIONS_v46.T12)
        T11_perim = perim_at_v46(result, WALL_CHAIN_POSITIONS_v46.T11)
        T9_core = core_at_v46(result, sensor_positions[:T9])
        T10_core = core_at_v46(result, sensor_positions[:T10])
        Tf = result.T3_sensor
        T2 = result.cavity_temperature
        h_mean = vec(mean(result.heat_transfer_coefficient, dims=1))
        return hcat(T8, T12_perim, T11_perim, T9_core, T10_core, Tf, T2, h_mean)
    end

    function align_cooling_initial_outputs_v46!(outputs, experiment)
        size(outputs, 1) >= 1 || return outputs
        outputs[1, 1:7] .= experiment[1, 1:7]
        return outputs
    end

    function solve_case_v46(p, simulation_id; is_cooling=false,
                           nodes=default_nodes,
                           rear_nodes=REAR_TUBE_DEFAULT_NODES_v46,
                           solver=Rodas5P(autodiff=AutoFiniteDiff()),
                           reltol=1e-6, abstol=1e-7, dtmax=30.0)
        id_str = string(simulation_id)
        conditions, time, experiment, data = experimental_case_v46(
            id_str; is_cooling=is_cooling
        )
        flow = observation(data, id_str, "_flow")
        Tin = observation(data, id_str, "_Tin")
        ambient = observation(data, id_str, "_Tamb")
        if is_cooling
            Tin = fill(COOLING_ROOM_TEMPERATURE_v46, length(time))
            ambient = fill(COOLING_ROOM_TEMPERATURE_v46, length(time))
        end
        T2 = observation(data, id_str, "_T2")
        irradiance = is_cooling ? zeros(length(time)) : fill(conditions[Io], length(time))
        operating = operating_condition_v3(
            irradiance=linear_history_v3(time, irradiance),
            flow_lpm=linear_history_v3(time, flow),
            inlet_temperature=linear_history_v3(time, Tin),
            ambient_temperature=linear_history_v3(time, ambient),
        )
        result = simulate_v46(
            p, operating, time;
            nodes=nodes,
            rear_nodes=rear_nodes,
            initial_temperature=measured_initial_core_profile_v46(data, id_str),
            initial_perim_temperature=measured_initial_perim_profile_v46(data, id_str),
            initial_rear_temperature=measured_initial_rear_profile_v46(data, id_str, p),
            initial_rear_reservoir_temperature=measured_initial_rear_reservoir_profile_v46(data, id_str, p),
            initial_cavity_temperature=T2[1],
            solver=solver,
            reltol=reltol,
            abstol=abstol,
            dtmax=dtmax,
        )
        outputs = extract_outputs_v46(result, p)
        is_cooling && align_cooling_initial_outputs_v46!(outputs, experiment)
        return outputs, result, experiment
    end
end

# ============================================================================
# OBJECTIVE FUNCTION & MULTI-SIGNAL LOSS (v46)
# ============================================================================
begin # v46 objective helpers
    const ORDERING_SCALE_v46 = 100.0
    const ORDERING_WEIGHT_v46 = 1.00
    const TIMING_WEIGHT_v46 = 0.15
    const SLOPE_WEIGHT_v46 = 0.25
    const COOLING_UPTURN_WEIGHT_v46 = 10000.0

    function normalized_slope_mse_v46(model, experiment)
        length(model) < 2 && return 0.0
        scale = max(maximum(experiment) - minimum(experiment), 1.0)
        return mean(abs2, diff(model) .- diff(experiment)) / scale^2
    end

    function timing_penalty_v46(time, model, experiment)
        duration = max(time[end] - time[1], 1.0)
        return ((get_t90_v3(time, model) - get_t90_v3(time, experiment)) / duration)^2
    end

    function signal_loss_v46(time, model, experiment)
        return normalized_mse_v3(model, experiment) +
               SLOPE_WEIGHT_v46 * normalized_slope_mse_v46(model, experiment) +
               TIMING_WEIGHT_v46 * timing_penalty_v46(time, model, experiment)
    end

    function cooling_upturn_loss_v46(model)
        target_columns = (3, 5, 6) # T11, T10, T3
        total = 0.0
        for j in target_columns
            scale = max(maximum(model[:, j]) - minimum(model[:, j]), 20.0)
            upturn = max.(diff(model[:, j]), 0.0) ./ scale
            total += maximum(upturn)^2 + mean(abs2, upturn)
        end
        return total / length(target_columns)
    end

    function ordering_loss_v46(model, experiment)
        model_offsets = (
            model[end, 2] - model[end, 1],
            model[end, 3] - model[end, 1],
            model[end, 2] - model[end, 4],
            model[end, 3] - model[end, 5],
        )
        experiment_offsets = (
            experiment[end, 2] - experiment[end, 1],
            experiment[end, 3] - experiment[end, 1],
            experiment[end, 2] - experiment[end, 4],
            experiment[end, 3] - experiment[end, 5],
        )
        return mean(abs2, collect(model_offsets) .- collect(experiment_offsets)) / ORDERING_SCALE_v46^2
    end

    function participating_total_heat_capacity_v46(p, nodes=15)
        dx = L / nodes
        C_core = rho_s * Cps_f(900.0) * A_solid * L
        C_perim = perimeter_heat_capacity_total_v46(p)
        C_rear = rear_reservoir_heat_capacity_v46(p)
        return C_core + C_perim + C_rear
    end

    function capacitance_regularization_v46(p; nodes=15)
        capacitance = participating_total_heat_capacity_v46(p, nodes)
        return CAPACITANCE_REGULARIZATION_WEIGHT_v46 *
               ((capacitance - MEASURED_ASSEMBLY_CAPACITANCE_v46) / MEASURED_ASSEMBLY_CAPACITANCE_v46)^2
    end
    
    function power_scale_regularization_v46(p)
        min_upper = min(p[5], p[6])
        if p[7] < min_upper - 0.1
            return 100.0 * (min_upper - 0.1 - p[7])^2
        end
        return 0.0
    end

    const keys_heating = sim_key_heat
    const keys_cooling = sim_key_cool

    function loss_cases_v46(p, keys; is_cooling=false, nodes=15)
        total_loss = 0.0
        data = is_cooling ? measurements_cooling : measurements
        for key in keys
            try
                outputs, _, experiment = solve_case_v46(
                    p, key;
                    is_cooling=is_cooling,
                    nodes=nodes,
                    solver=Rodas5P(autodiff=AutoFiniteDiff()),
                    reltol=1e-4,
                    abstol=1e-5,
                    dtmax=60.0,
                )
                case_loss = 0.0
                time = observation_time(data, key)
                for j in 1:7
                    case_loss += signal_loss_v46(time, outputs[:, j], experiment[:, j])
                end
                case_loss /= 7.0
                if is_cooling
                    case_loss += COOLING_UPTURN_WEIGHT_v46 * cooling_upturn_loss_v46(outputs)
                else
                    case_loss += ORDERING_WEIGHT_v46 * ordering_loss_v46(outputs, experiment)
                end
                total_loss += case_loss
            catch e
                return 1e6
            end
        end
        return total_loss / length(keys)
    end

    function calibration_loss_v46(p; stage=:full, nodes=15)
        for i in eachindex(p)
            (p[i] < lb_full_v46[i] || p[i] > ub_full_v46[i]) && return 1e6 + 1e3 * (abs(p[i] - clamp(p[i], lb_full_v46[i], ub_full_v46[i])))
        end
        
        heating_loss = loss_cases_v46(p, keys_heating; is_cooling=false, nodes=nodes)
        reg = capacitance_regularization_v46(p; nodes=nodes) + power_scale_regularization_v46(p)
        return heating_loss + reg
    end
end

# ============================================================================
# EXACT ENERGY CONSERVATION AUDIT FUNCTION (v46)
# ============================================================================
begin # v46 conservation audit
    function compute_energy_balance_v46(result, operating, time_val, p; nodes=default_nodes, rear_nodes=REAR_TUBE_DEFAULT_NODES_v46)
        dx = L / nodes
        dx_rear = REAR_TUBE_LENGTH_v46 / rear_nodes
        irradiance = operating.irradiance(time_val)
        ambient = operating.ambient_temperature(time_val)
        flow = max(0.0, operating.flow_lpm(time_val))
        Tin = operating.inlet_temperature(time_val)

        Q_aperture = max(0.0, irradiance) * A_frt
        M = absorbed_power_scale_v46(irradiance, p)
        Q_delivered = M * Q_aperture
        chi = core_source_fraction_v46(p)
        Qcore_solar = chi * Q_delivered
        Qperim_solar = (1.0 - chi) * Q_delivered

        # Temperatures at this time
        idx = argmin(abs.(result.time .- time_val))
        Tcore = result.core_temperature[:, idx]
        Tperim = result.perim_temperature[:, idx]
        Ttube = result.rear_temperature[:, idx]
        Tcavity = result.cavity_temperature[idx]
        Tg = result.gas_temperature[:, idx]

        mdot = m_dot_standard_v46(flow)
        
        # Suction heat preheating the gas at z = 0
        Q_suction = mdot * cpf_f(Tin) * (Tg[1] - Tin)
        Q_gas_core = sum(result.receiver_gas_heat[:, idx])
        Q_gas_rear = sum(result.rear_tube_gas_heat[:, idx])
        Q_gas_total = Q_suction + Q_gas_core + Q_gas_rear
        
        Q_front_rad = EPS_FRONT_FIXED_v46 * sigma_sb * A_frt * (Tperim[1]^4 - ambient^4)
        Q_cavity_amb = cavity_loss_to_ambient_v46(Tcavity, ambient)
        
        Q_flange = 0.0
        for j in 1:rear_nodes
            G_fl = rear_tube_flange_conductance_per_length_v46(result.z_rear[j], Ttube[j], p)
            Q_flange += G_fl * dx_rear * (Ttube[j] - WATER_FLANGE_TEMPERATURE_v46)
        end

        # Calculate instantaneous sensible storage rate dE_stored/dt for the solid network
        u_k = result.ode_solution.u[idx]
        du_k = similar(u_k)
        core_weights = solar_weights_v5(p[11], p[4], nodes)
        rear_weights = rear_contact_weights_v46(result.z_solid)
        context = (
            p=p, operating=operating, z=result.z_solid, dx=dx,
            core_weights=core_weights, rear_weights=rear_weights,
            z_rear=result.z_rear, dx_rear=dx_rear,
            Tg=zeros(nodes + rear_nodes + 1), Qgas=zeros(nodes), hcell=zeros(nodes),
            Qgas_rear=zeros(rear_nodes), hrear=zeros(rear_nodes), zg=result.z_gas,
        )
        receiver_rhs_v46!(
            du_k, u_k, time_val, p, operating,
            result.z_solid, dx, core_weights, rear_weights,
            result.z_rear, dx_rear,
            context.Tg, context.Qgas, context.hcell, context.Qgas_rear, context.hrear, context.zg
        )
        
        dE_stored = 0.0
        cell_volume = A_solid * dx
        perim_cap_cell = perimeter_heat_capacity_total_v46(p) / nodes
        for i in 1:nodes
            core_cap = rho_s * Cps_f(Tcore[i]) * cell_volume
            dE_stored += core_cap * du_k[i]
            dE_stored += perim_cap_cell * du_k[nodes + i]
        end
        for j in 1:rear_nodes
            tube_cap = rear_tube_capacity_v46(Ttube[j], result.z_rear[j], dx_rear)
            dE_stored += tube_cap * du_k[2 * nodes + j]
        end
        cavity_index = 2 * nodes + rear_nodes + 1
        dE_stored += CAVITY_HEAT_CAPACITY_v46 * du_k[cavity_index]
        C_rear_tot = rear_reservoir_heat_capacity_v46(p)
        for i in 1:nodes
            dE_stored += (C_rear_tot * rear_weights[i]) * du_k[cavity_index + i]
        end

        instantaneous_residual = Q_delivered - (Q_gas_total + Q_front_rad + Q_cavity_amb + Q_flange + dE_stored)
        steady_flux_residual = Q_delivered - (Q_gas_total + Q_front_rad + Q_cavity_amb + Q_flange)

        return (
            delivered_W=Q_delivered,
            core_solar_W=Qcore_solar,
            perim_solar_W=Qperim_solar,
            gas_heat_W=Q_gas_total,
            front_rad_loss_W=Q_front_rad,
            cavity_amb_loss_W=Q_cavity_amb,
            flange_loss_W=Q_flange,
            dE_stored_dt_W=dE_stored,
            instantaneous_residual_W=instantaneous_residual,
            residual_W=instantaneous_residual,
            steady_flux_residual_W=steady_flux_residual,
        )
    end

    function fit_indices_for_stage_v46(stage)
        stage_symbol = Symbol(stage)
        stage_symbol == :rear && return fit_rear_stage_indices_v46
        stage_symbol == :transport && return fit_transport_stage_indices_v46
        stage_symbol == :full && return fit_full_stage_indices_v46
        stage_symbol == :all && return collect(1:23)
        throw(ArgumentError("unknown v46 calibration stage: $stage"))
    end

    function calibrate_v46(; nodes=15, maximum_iterations=500,
                           maximum_time_seconds=1800.0, stage=:full, dataset=:heating)
        p0 = copy(pnew_v46)
        idx_fit = fit_indices_for_stage_v46(stage)
        function obj_free(theta)
            p_full = copy(p0)
            p_full[idx_fit] .= theta
            if dataset == :heating
                return loss_cases_v46(p_full, keys_heating; is_cooling=false, nodes=nodes) +
                       capacitance_regularization_v46(p_full; nodes=nodes) + power_scale_regularization_v46(p_full)
            elseif dataset == :cooling
                return loss_cases_v46(p_full, keys_cooling; is_cooling=true, nodes=nodes) +
                       capacitance_regularization_v46(p_full; nodes=nodes)
            else
                return loss_cases_v46(p_full, keys_heating; is_cooling=false, nodes=nodes) +
                       loss_cases_v46(p_full, keys_cooling; is_cooling=true, nodes=nodes) +
                       capacitance_regularization_v46(p_full; nodes=nodes) + power_scale_regularization_v46(p_full)
            end
        end
        res = optimize_with_nlopt_v3(
            obj_free,
            p0[idx_fit],
            lb_full_v46[idx_fit],
            ub_full_v46[idx_fit];
            maximum_iterations=maximum_iterations,
            maximum_time_seconds=maximum_time_seconds,
            algorithm=NLopt.LN_BOBYQA(),
            label="fit-v46-$(Symbol(stage))",
        )
        p_opt = copy(p0)
        p_opt[idx_fit] .= res.minimizer
        return (objective=res.minimum, parameters=p_opt, minimizer=res.minimizer, retcode=res.retcode)
    end
end
