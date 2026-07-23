# ============================================================================
# 1D_v24.jl - Energy-accounting 2-Zone Core/Perimeter Macro Model
# ============================================================================
# v24 keeps the v19 two-zone idea but makes the energy accounting explicit:
# aperture-average irradiance, flux-derived core/perimeter power partition, corrected
# zero-irradiance cooling cases, and a smooth downstream rear/flange heat-sink
# influence that can reach T10 without a discontinuous sink onset. The rear
# tube initial state is tied to the same core/perimeter exit blend used by the
# receiver-to-tube conductance, avoiding artificial heat injection into T10 at
# the beginning of cooling.
#
# Key Physical Features of v24:
#   1. Core Channel Zone (T_core): 100-channel matrix. Receives direct central
#      solar beam fraction. Exchanges heat with the active gas stream using
#      bounded laminar developing-flow convection. Compares to T9 and T10.
#   2. Perimeter Housing Zone (T_perim): Side walls, alumina felt insulation,
#      and front housing. Receives a captured fraction of off-receiver flux
#      spillover, with its own axial deposition shape. Compares to T8,
#      T12, T11 and drives the T2/cavity state.
#   3. Rear Heat Sink: the water-cooled flange conductance can be scaled to
#      test whether a stronger downstream loss fixes T3 and rear thermocouples.
#
# Parameter Vector p[1:21] in v24:
#   p[1]  A_Nu          laminar developing Nu prefactor
#   p[2]  B_Re          laminar developing Reynolds exponent (bounded 0.0 to 0.6)
#   p[3]  C_Pr          fixed Prandtl exponent (= 1/3)
#   p[4]  phi_0         active flow fraction at Re = 50
#   p[5]  m_rec         flow recruitment exponent
#   p[6]  front_dep     frozen v12a front-cell deposition fraction (= 1.0 or 0.0)
#   p[7]  scale_456     total absorbed power scale at 456 kW/m2
#   p[8]  scale_304     total absorbed power scale at 304 kW/m2
#   p[9]  scale_256     total absorbed power scale at 256 kW/m2
#   p[10] G_core_perim  radial core-perimeter conductance per length [W/m/K]
#   p[11] C_perim_eff   perimeter/housing participating heat capacity [J/K]
#   p[12] k_perim_ref   perimeter axial conductivity at 900 K [W/m/K]
#   p[13] beta_opt      frozen v12a optical attenuation [1/m]
#   p[14] spill_capture captured fraction of off-receiver flux spillover
#   p[15] beta_perim    axial attenuation of perimeter/rim source [1/m]
#   p[16] f_core_tube   fraction of receiver-exit solid heat leak from core
#   p[17] flange_scale  multiplier for rear-tube-to-water-flange conductance
#   p[18] G_rear_core   downstream core-to-flange loss per length [W/m/K]
#   p[19] G_rear_perim  downstream perimeter-to-flange loss per length [W/m/K]
#   p[20] rear_sink_shape exponent controlling how abruptly rear sink grows
#   p[21] k_core_axial_scale fixed effective axial solid-conduction multiplier
# ============================================================================

include("1D_v5.jl")

begin # v24 fixed constants
    const EPS_FRONT_FIXED_v24 = EPS_FRONT_FIXED_V5
    const ETA_ABS_FIXED_v24 = ETA_ABS_FIXED_V5
    const BETA_OPT_FIXED_v24 = BETA_OPT_FIXED_V5
    const H_FRONT_FIXED_v24 = 10.0
    const PRANDTL_EXPONENT_FIXED_v24 = 1.0 / 3.0
    const STANDARD_AIR_DENSITY_v24 = 101325.0 / (287.05 * 294.25)
    const REYNOLDS_REFERENCE_v24 = 50.0
    const NU_FLOOR_v24 = 3.61
    const MEASURED_ASSEMBLY_CAPACITANCE_v24 = 301.0
    const IRRADIANCE_LEVELS_v24 = (456000.0, 304000.0, 256000.0)
    const T3_SAMPLE_POSITION_v24 = 140.0e-3
    const CAPACITANCE_REGULARIZATION_WEIGHT_v24 = 0.10
    const FLUX_SIGMA_FIXED_v24 = 30.0e-3
    const FLUX_INTEGRATION_POINTS_v24 = 161

    # Cavity and housing geometry
    const CAVITY_OUTER_RADIUS_v24 = 75.0e-3
    const APERTURE_AREA_v24 = pi * CAVITY_OUTER_RADIUS_v24^2
    const METAL_THICKNESS_v24 = 18.0e-3
    const INSULATION_OUTER_RADIUS_v24 = CAVITY_OUTER_RADIUS_v24 - METAL_THICKNESS_v24
    const METAL_OUTER_RADIUS_v24 = CAVITY_OUTER_RADIUS_v24
    const CAVITY_LENGTH_v24 = 165.0e-3

    const ADAPTOR_DIAMETER_v24 = 77.6e-3
    const ADAPTOR_RADIUS_v24 = ADAPTOR_DIAMETER_v24 / 2.0
    const ADAPTOR_LENGTH_v24 = 57.0e-3
    const ADAPTOR_TUBE_RADIUS_v24 = 6.5e-3
    const ADAPTOR_RECEIVER_OVERLAP_LENGTH_v24 = 28.0e-3
    const ADAPTOR_TUBE_OVERLAP_LENGTH_v24 = ADAPTOR_LENGTH_v24 - ADAPTOR_RECEIVER_OVERLAP_LENGTH_v24
    const ADAPTOR_OVERLAP_LENGTH_v24 = ADAPTOR_RECEIVER_OVERLAP_LENGTH_v24

    const RECEIVER_EQ_RADIUS_v24 = sqrt(A_frt / pi)
    const RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v24 =
        4.0 * pi * 0.08 / log(INSULATION_OUTER_RADIUS_v24 / RECEIVER_EQ_RADIUS_v24)
    const ADAPTOR_CONTACT_RESISTANCE_v24 =
        1.0 / (100.0 * 4.0 * w_t * ADAPTOR_OVERLAP_LENGTH_v24)
    const RECEIVER_TO_TUBE_CONDUCTANCE_v24 = 1.0 / ADAPTOR_CONTACT_RESISTANCE_v24

    # Rear tube geometry
    const REAR_TUBE_LENGTH_v24 = 150.0e-3
    const REAR_TUBE_CAVITY_LENGTH_v24 = 46.0e-3
    const REAR_TUBE_DEFAULT_NODES_v24 = 12
    const REAR_TUBE_GAS_RADIUS_v24 = ADAPTOR_TUBE_RADIUS_v24
    const REAR_TUBE_WALL_THICKNESS_v24 = 1.5e-3
    const REAR_TUBE_OUTER_RADIUS_v24 = REAR_TUBE_GAS_RADIUS_v24 + REAR_TUBE_WALL_THICKNESS_v24
    const REAR_TUBE_FLOW_AREA_v24 = pi * REAR_TUBE_GAS_RADIUS_v24^2
    const REAR_TUBE_HYDRAULIC_DIAMETER_v24 = 2.0 * REAR_TUBE_GAS_RADIUS_v24
    const REAR_TUBE_EXCHANGE_PERIMETER_v24 = 2.0 * pi * REAR_TUBE_GAS_RADIUS_v24

    const ALUMINA_DENSITY_v24 = 3900.0
    const ALUMINUM_DENSITY_v24 = 2700.0
    const ALUMINUM_HEAT_CAPACITY_v24 = 900.0
    const ALUMINUM_CONDUCTIVITY_v24 = 205.0
    const METAL_EMISSIVITY_v24 = 0.2
    const H_NAT_CAVITY_v24 = 10.0
    const INSULATION_DENSITY_v24 = 140.0
    const INSULATION_HEAT_CAPACITY_v24 = 1360.0
    const WATER_FLANGE_TEMPERATURE_v24 = 25.0 + 273.15

    const CAVITY_FELT_VOLUME_v24 = max(
        0.0,
        pi * INSULATION_OUTER_RADIUS_v24^2 * CAVITY_LENGTH_v24 -
        A_frt * L - pi * ADAPTOR_RADIUS_v24^2 * ADAPTOR_LENGTH_v24,
    )
    const CAVITY_METAL_VOLUME_v24 =
        pi * (METAL_OUTER_RADIUS_v24^2 - INSULATION_OUTER_RADIUS_v24^2) * CAVITY_LENGTH_v24 +
        pi * METAL_OUTER_RADIUS_v24^2 * METAL_THICKNESS_v24
    const CAVITY_OUTER_AREA_v24 =
        2.0 * pi * METAL_OUTER_RADIUS_v24 * CAVITY_LENGTH_v24 + pi * METAL_OUTER_RADIUS_v24^2
    const CAVITY_HEAT_CAPACITY_v24 =
        INSULATION_DENSITY_v24 * INSULATION_HEAT_CAPACITY_v24 * CAVITY_FELT_VOLUME_v24 +
        ALUMINUM_DENSITY_v24 * ALUMINUM_HEAT_CAPACITY_v24 * CAVITY_METAL_VOLUME_v24
end

begin # v24 non-uniform front-flux partition
    transverse_flux_shape_v24(x, y; sigma=FLUX_SIGMA_FIXED_v24) =
        exp(-0.5 * (x^2 + y^2) / sigma^2)

    function flux_zone_fractions_v24(; radius=CAVITY_OUTER_RADIUS_v24,
                                     receiver_half_width=w_t / 2.0,
                                     points=FLUX_INTEGRATION_POINTS_v24)
        dx_flux = 2.0 * radius / (points - 1)
        receiver_power = 0.0
        aperture_power = 0.0
        for ix in 1:points
            x = -radius + (ix - 1) * dx_flux
            for iy in 1:points
                y = -radius + (iy - 1) * dx_flux
                x^2 + y^2 > radius^2 && continue
                weight = transverse_flux_shape_v24(x, y) * dx_flux^2
                aperture_power += weight
                if abs(x) <= receiver_half_width && abs(y) <= receiver_half_width
                    receiver_power += weight
                end
            end
        end
        receiver_fraction = receiver_power / max(aperture_power, eps(Float64))
        return (
            receiver=receiver_fraction,
            spillover=clamp(1.0 - receiver_fraction, 0.0, 1.0),
        )
    end

    const FLUX_ZONE_FRACTIONS_FIXED_v24 = flux_zone_fractions_v24()
end

# Default v24 parameter vector p[1:21]
begin # v24 parameter values and bounds
    pnew_v24 = [
        4.924815049912292,  # p[1]  A_Nu (laminar developing prefactor)
        0.5162678835987395, # p[2]  B_Re (Reynolds exponent bounded 0.0 to 0.6)
        PRANDTL_EXPONENT_FIXED_v24, # p[3]  C_Pr (= 1/3)
        1.0,                # p[4]  phi_0 (fixed full active flow fraction)
        0.0,                # p[5]  m_rec (flow recruitment exponent)
        1.0,                # p[6]  front_dep (frozen v12a front-cell deposition)
        0.70,               # p[7]  scale_456
        0.70,               # p[8]  scale_304
        0.70,               # p[9]  scale_256
        16.044477788186573, # p[10] G_core_perim (radial core-perimeter conductance [W/m/K])
        214.29323236142696, # p[11] C_perim_eff (perimeter participating capacity [J/K])
        7.591183315111345,  # p[12] k_perim_ref (perimeter axial conductivity at 900K [W/m/K])
        184.67243519237965, # p[13] beta_opt [1/m]
        0.02,               # p[14] spill_capture (captured spillover fraction)
        1.8019923585149422, # p[15] beta_perim (perimeter source attenuation [1/m])
        0.9993097154917425, # p[16] f_core_tube (rear-tube heat leak core fraction)
        0.1013817714666468, # p[17] flange_scale (rear/flange heat-removal multiplier)
        3.336402359922935,  # p[18] G_rear_core (rear core sink [W/m/K])
        0.0,                # p[19] G_rear_perim (rear perimeter sink [W/m/K])
        0.9997676098325551, # p[20] rear_sink_shape (lower => earlier rear-sink influence)
        0.0,                # p[21] k_core_axial_scale (fixed off for monotonic cooling)
    ]
    lb_full_v24 = [
        0.01,
        0.0,
        pnew_v24[3],
        pnew_v24[4],
        pnew_v24[5],
        pnew_v24[6],
        0.02,
        0.02,
        0.02,
        0.5,
        50.0,
        0.0,
        pnew_v24[13],
        0.0,
        0.0,
        0.0,
        0.10,
        0.0,
        0.0,
        0.25,
        0.0,
    ]
    ub_full_v24 = [
        5.00,
        0.60,
        pnew_v24[3],
        pnew_v24[4],
        pnew_v24[5],
        pnew_v24[6],
        2.50,
        2.50,
        2.50,
        100.0,
        450.0,
        50.0,
        pnew_v24[13],
        1.00,
        300.0,
        1.0,
        20.0,
        300.0,
        300.0,
        1.0,
        0.125,
    ]

    fit_heat_transfer_indices_v24 = [1, 2, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20]
    fit_source_indices_v24 = [14, 15, 16, 17, 18, 19, 20]
    fit_power_scale_indices_v24 = [7, 8, 9]
    fit_radiation_indices_v24 = Int[]
end

function irradiance_level_index_v24(irradiance)
    irradiance <= 0.0 && return 0
    differences = abs.(collect(IRRADIANCE_LEVELS_v24) .- Float64(irradiance))
    return argmin(differences)
end

power_scale_index_v24(irradiance) = begin
    level_index = irradiance_level_index_v24(irradiance)
    level_index == 0 ? 0 : 6 + level_index
end

function absorbed_power_scale_v24(irradiance, p)
    index = power_scale_index_v24(irradiance)
    return index == 0 ? 1.0 : p[index]
end

core_perimeter_conductance_per_length_v24(p) = max(p[10], eps(Float64))

perimeter_heat_capacity_total_v24(p) = max(p[11], eps(Float64))

core_receiver_heat_capacity_v24(nodes) =
    sum(rho_s * Cps_f(900.0) * A_solid * (L / nodes) for _ in 1:nodes)

participating_assembly_heat_capacity_v24(p, nodes=default_nodes) =
    core_receiver_heat_capacity_v24(nodes) + perimeter_heat_capacity_total_v24(p)

perimeter_axial_conductivity_v24(T, p) =
    max(p[12], 0.0) * (property_temperature(T) / 900.0)^3

flux_receiver_fraction_v24() = FLUX_ZONE_FRACTIONS_FIXED_v24.receiver

flux_spillover_fraction_v24() = FLUX_ZONE_FRACTIONS_FIXED_v24.spillover

perimeter_spillover_capture_v24(p) = clamp(p[14], 0.0, 1.0)

perimeter_source_fraction_v24(p) =
    clamp(perimeter_spillover_capture_v24(p) * flux_spillover_fraction_v24(), 0.0, 0.80)

modeled_absorbed_aperture_power_v24(irradiance, p) =
    ETA_ABS_FIXED_v24 * absorbed_power_scale_v24(irradiance, p) *
    max(0.0, irradiance) * APERTURE_AREA_v24

core_absorbed_power_v24(irradiance, p) =
    modeled_absorbed_aperture_power_v24(irradiance, p) * flux_receiver_fraction_v24()

perimeter_absorbed_power_v24(irradiance, p) =
    modeled_absorbed_aperture_power_v24(irradiance, p) * perimeter_source_fraction_v24(p)

modeled_participating_absorbed_power_v24(irradiance, p) =
    core_absorbed_power_v24(irradiance, p) + perimeter_absorbed_power_v24(irradiance, p)

perimeter_source_attenuation_v24(p) = clamp(p[15], 0.0, 300.0)

core_tube_fraction_v24(p) = clamp(p[16], 0.0, 1.0)

flange_loss_scale_v24(p) = clamp(p[17], 0.10, 20.0)

rear_core_loss_per_length_v24(p) = clamp(p[18], 0.0, 300.0)

rear_perim_loss_per_length_v24(p) = clamp(p[19], 0.0, 300.0)

rear_sink_shape_v24(p) = clamp(p[20], 0.25, 1.0)

core_axial_conduction_scale_v24(p) = clamp(p[21], 0.0, 0.125)

function rear_receiver_sink_weight_v24(z, p=pnew_v24)
    z0 = sensor_positions[:T8]
    z <= z0 && return 0.0
    s = clamp((z - z0) / max(L - z0, eps(Float64)), 0.0, 1.0)
    return s^rear_sink_shape_v24(p)
end

function axial_exponential_weights_v24(beta, nodes, dx)
    beta_eff = max(beta, 0.0)
    if beta_eff <= 1e-9
        return fill(1.0 / nodes, nodes)
    end
    raw = [
        exp(-beta_eff * ((i - 1) * dx)) - exp(-beta_eff * (i * dx))
        for i in 1:nodes
    ]
    total = sum(raw)
    total <= eps(Float64) && return fill(1.0 / nodes, nodes)
    return raw ./ total
end

function active_flow_fraction_v24(reynolds, p)
    phi0 = clamp(p[4], 0.1, 1.0)
    m = max(p[5], 0.0)
    ratio = max(reynolds, eps(Float64)) / REYNOLDS_REFERENCE_v24
    return clamp(phi0 * ratio^m, 0.1, 1.0)
end

function receiver_active_flow_fraction_v24(flow_lpm, Tin, p)
    flow_lpm <= 1e-12 && return 1.0
    mdot = m_dot_standard_v24(flow_lpm)
    mu = mu_f_f(Tin)
    reynolds = mdot * Dh / (A_flow * mu)
    return active_flow_fraction_v24(reynolds, p)
end

function receiver_channel_nusselt_v24(z, reynolds, prandtl, p)
    A_nu = max(p[1], 0.0)
    B_re = clamp(p[2], 0.0, 0.6)
    C_pr = p[3]
    entry_factor = (Dh / max(z, Dh / 2.0))^(1.0 / 3.0)
    Nu_dev = A_nu * max(reynolds, eps(Float64))^B_re * max(prandtl, eps(Float64))^C_pr * entry_factor
    return max(NU_FLOOR_v24, Nu_dev)
end

m_dot_standard_v24(flow_lpm) =
    STANDARD_AIR_DENSITY_v24 * max(0.0, flow_lpm) / 60000.0

alumina_conductivity_v24(T) =
    5.5 + 34.5 * exp(-0.0033 * (property_temperature(T) - 273.15))

alumina_heat_capacity_v24(T) =
    max(500.0, (1.00446 + 1.742e-4 * property_temperature(T) -
                2.796e4 * property_temperature(T)^-2) * 1000.0)

function rear_tube_outer_radius_v24(z_rear)
    return z_rear <= ADAPTOR_TUBE_OVERLAP_LENGTH_v24 ?
           ADAPTOR_RADIUS_v24 : REAR_TUBE_OUTER_RADIUS_v24
end

function rear_tube_solid_area_v24(z_rear)
    outer_radius = rear_tube_outer_radius_v24(z_rear)
    return pi * (outer_radius^2 - REAR_TUBE_GAS_RADIUS_v24^2)
end

function rear_tube_cavity_conductance_per_length_v24(z_rear, p)
    z_rear > REAR_TUBE_CAVITY_LENGTH_v24 && return 0.0
    outer_radius = min(rear_tube_outer_radius_v24(z_rear), 0.98 * INSULATION_OUTER_RADIUS_v24)
    return 2.0 * pi * 0.08 / log(INSULATION_OUTER_RADIUS_v24 / outer_radius)
end

function rear_tube_flange_conductance_per_length_v24(z_rear, Ttube, p)
    z_rear <= REAR_TUBE_CAVITY_LENGTH_v24 && return 0.0
    tube_wall_resistance =
        log(REAR_TUBE_OUTER_RADIUS_v24 / REAR_TUBE_GAS_RADIUS_v24) /
        (2.0 * pi * alumina_conductivity_v24(Ttube))
    aluminum_resistance =
        log(METAL_OUTER_RADIUS_v24 / REAR_TUBE_OUTER_RADIUS_v24) /
        (2.0 * pi * ALUMINUM_CONDUCTIVITY_v24)
    return flange_loss_scale_v24(p) / (tube_wall_resistance + aluminum_resistance)
end

function cavity_loss_to_ambient_v24(Tcavity, ambient)
    return H_NAT_CAVITY_v24 * CAVITY_OUTER_AREA_v24 * (Tcavity - ambient) +
           METAL_EMISSIVITY_v24 * sigma_sb * CAVITY_OUTER_AREA_v24 *
           (Tcavity^4 - ambient^4)
end

function rear_tube_capacity_v24(Ttube, z_rear, dx_rear)
    return ALUMINA_DENSITY_v24 * alumina_heat_capacity_v24(Ttube) *
           rear_tube_solid_area_v24(z_rear) * dx_rear
end

function rear_tube_htc_v24(Twall, Tgas, flow_lpm, Tin)
    flow_lpm <= 1e-12 && return 0.0
    mdot = m_dot_standard_v24(flow_lpm)
    Tfilm = 0.5 * (Twall + Tgas)
    cp = cpf_f(Tfilm)
    mu = mu_f_f(Tfilm)
    kf = kf_f(Tfilm)
    reynolds = mdot * REAR_TUBE_HYDRAULIC_DIAMETER_v24 / (REAR_TUBE_FLOW_AREA_v24 * mu)
    prandtl = cp * mu / kf
    nusselt = reynolds < 2300.0 ? 3.66 : 0.023 * reynolds^0.8 * prandtl^0.4
    return nusselt * kf / REAR_TUBE_HYDRAULIC_DIAMETER_v24
end

function gas_profile_v24!(Tg, Qgas, hcell, Qgas_rear, hrear,
                          Tcore, Ttube, time, p, operating, z, dx, z_rear, dx_rear)
    flow = max(0.0, operating.flow_lpm(time))
    Tin = operating.inlet_temperature(time)
    receiver_nodes = length(Tcore)
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

    mdot_total = m_dot_standard_v24(flow)
    phi_act = receiver_active_flow_fraction_v24(flow, Tin, p)
    mdot_active = phi_act * mdot_total
    Tactive = Tin
    for i in eachindex(Tcore)
        Tfilm = 0.5 * (Tcore[i] + Tactive)
        cp = cpf_f(Tfilm)
        mu = mu_f_f(Tfilm)
        kf = kf_f(Tfilm)
        reynolds_total = mdot_total * Dh / (A_flow * mu)
        prandtl = cp * mu / kf

        reynolds_active = reynolds_total * phi_act
        
        nusselt = receiver_channel_nusselt_v24(z[i], reynolds_active, prandtl, p)
        hcell[i] = nusselt * kf / Dh
        UA = hcell[i] * P_exchange * dx
        effectiveness = -expm1(-UA / (mdot_active * cp))

        Tactive_next = Tactive + effectiveness * (Tcore[i] - Tactive)
        Qgas[i] = mdot_active * cp * (Tactive_next - Tactive)
        Tactive = Tactive_next
        Tg[i + 1] = phi_act * Tactive + (1.0 - phi_act) * Tin
    end

    for j in eachindex(Ttube)
        inlet_index = receiver_nodes + j
        Tfilm = 0.5 * (Ttube[j] + Tg[inlet_index])
        cp = cpf_f(Tfilm)
        hrear[j] = rear_tube_htc_v24(Ttube[j], Tg[inlet_index], flow, Tin)
        UA = hrear[j] * REAR_TUBE_EXCHANGE_PERIMETER_v24 * dx_rear
        effectiveness = -expm1(-UA / (mdot_total * cp))
        Tg[inlet_index + 1] = Tg[inlet_index] + effectiveness * (Ttube[j] - Tg[inlet_index])
        Qgas_rear[j] = mdot_total * cp * (Tg[inlet_index + 1] - Tg[inlet_index])
    end
    return nothing
end

function receiver_rhs_v24!(du, u, time, p, operating, z, dx, core_weights,
                          perim_weights,
                          z_rear, dx_rear,
                          Tg, Qgas, hcell, Qgas_rear, hrear)
    nodes = length(z)
    rear_nodes = length(z_rear)
    Tcore = view(u, 1:nodes)
    Tperim = view(u, (nodes + 1):(2 * nodes))
    tube_start = 2 * nodes + 1
    tube_stop = 2 * nodes + rear_nodes
    Ttube = view(u, tube_start:tube_stop)
    cavity_index = 2 * nodes + rear_nodes + 1
    Tcavity = u[cavity_index]
    
    dTcore = view(du, 1:nodes)
    dTperim = view(du, (nodes + 1):(2 * nodes))
    dTtube = view(du, tube_start:tube_stop)
    ambient = operating.ambient_temperature(time)

    gas_profile_v24!(Tg, Qgas, hcell, Qgas_rear, hrear,
                    Tcore, Ttube, time, p, operating, z, dx, z_rear, dx_rear)
    fill!(du, 0.0)

    # Core axial solid conduction
    for i in 1:(nodes - 1)
        ki = ks_f(Tcore[i])
        kj = ks_f(Tcore[i + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        Qcond_core =
            core_axial_conduction_scale_v24(p) * kface * A_solid *
            (Tcore[i] - Tcore[i + 1]) / dx
        dTcore[i] -= Qcond_core
        dTcore[i + 1] += Qcond_core

        kwi = ks_f(Tperim[i])
        kwj = ks_f(Tperim[i + 1])
        kwface = 2.0 * kwi * kwj / (kwi + kwj)
        Tface_perim = 0.5 * (Tperim[i] + Tperim[i + 1])
        k_perim = perimeter_axial_conductivity_v24(Tface_perim, p)
        Qcond_perim = (kwface * A_solid * 0.20 + k_perim * A_frt) * (Tperim[i] - Tperim[i + 1]) / dx
        dTperim[i] -= Qcond_perim
        dTperim[i + 1] += Qcond_perim
    end

    irradiance = max(0.0, operating.irradiance(time))
    scale = absorbed_power_scale_v24(irradiance, p)
    Qaperture_absorbed = ETA_ABS_FIXED_v24 * scale * irradiance * APERTURE_AREA_v24
    f_perim_source = perimeter_source_fraction_v24(p)
    Qcore_solar = Qaperture_absorbed * flux_receiver_fraction_v24()
    Qperim_solar = Qaperture_absorbed * f_perim_source

    G_cp = core_perimeter_conductance_per_length_v24(p)
    radial_cavity_conductance = RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v24

    for i in 1:nodes
        Qrad_in = G_cp * dx * (Tperim[i] - Tcore[i])
        Qperim_cavity = radial_cavity_conductance * dx * (Tperim[i] - Tcavity)
        
        Qcell_core = Qcore_solar * core_weights[i]
        Qcell_perim = Qperim_solar * perim_weights[i]
        rear_sink_weight = rear_receiver_sink_weight_v24(z[i], p)
        Qrear_core_sink =
            rear_core_loss_per_length_v24(p) * dx * rear_sink_weight *
            (Tcore[i] - WATER_FLANGE_TEMPERATURE_v24)
        Qrear_perim_sink =
            rear_perim_loss_per_length_v24(p) * dx * rear_sink_weight *
            (Tperim[i] - WATER_FLANGE_TEMPERATURE_v24)

        dTcore[i] += Qcell_core - Qgas[i] + Qrad_in - Qrear_core_sink
        dTperim[i] += Qcell_perim - Qrad_in - Qperim_cavity - Qrear_perim_sink
        du[cavity_index] += Qperim_cavity
    end

    Qfront_perim = H_FRONT_FIXED_v24 * A_frt * (Tperim[1] - ambient) +
                   EPS_FRONT_FIXED_v24 * sigma_sb * A_frt * (Tperim[1]^4 - ambient^4)
    
    f_core_tube = core_tube_fraction_v24(p)
    Qreceiver_tube_core = RECEIVER_TO_TUBE_CONDUCTANCE_v24 * f_core_tube * (Tcore[end] - Ttube[1])
    Qreceiver_tube_perim = RECEIVER_TO_TUBE_CONDUCTANCE_v24 * (1.0 - f_core_tube) * (Tperim[end] - Ttube[1])
    Qreceiver_tube = Qreceiver_tube_core + Qreceiver_tube_perim

    dTperim[1] -= Qfront_perim
    dTcore[end] -= Qreceiver_tube_core
    dTperim[end] -= Qreceiver_tube_perim
    dTtube[1] += Qreceiver_tube

    for j in 1:(rear_nodes - 1)
        ki = alumina_conductivity_v24(Ttube[j])
        kj = alumina_conductivity_v24(Ttube[j + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        Aface = 0.5 * (rear_tube_solid_area_v24(z_rear[j]) + rear_tube_solid_area_v24(z_rear[j + 1]))
        Qcond = kface * Aface * (Ttube[j] - Ttube[j + 1]) / dx_rear
        dTtube[j] -= Qcond
        dTtube[j + 1] += Qcond
    end

    for j in 1:rear_nodes
        G_cavity = rear_tube_cavity_conductance_per_length_v24(z_rear[j], p)
        G_flange = rear_tube_flange_conductance_per_length_v24(z_rear[j], Ttube[j], p)
        Q_cavity = G_cavity * dx_rear * (Ttube[j] - Tcavity)
        Q_flange = G_flange * dx_rear * (Ttube[j] - WATER_FLANGE_TEMPERATURE_v24)
        dTtube[j] -= Qgas_rear[j] + Q_cavity + Q_flange
        du[cavity_index] += Q_cavity
    end

    Qcavity_ambient = cavity_loss_to_ambient_v24(Tcavity, ambient)
    du[cavity_index] -= Qcavity_ambient

    cell_volume = A_solid * dx
    perim_capacity_cell = perimeter_heat_capacity_total_v24(p) / nodes
    for i in 1:nodes
        core_capacity = rho_s * Cps_f(Tcore[i]) * cell_volume
        dTcore[i] /= max(core_capacity, eps(Float64))
        dTperim[i] /= max(perim_capacity_cell, eps(Float64))
    end
    for j in 1:rear_nodes
        dTtube[j] /= max(rear_tube_capacity_v24(Ttube[j], z_rear[j], dx_rear), eps(Float64))
    end
    du[cavity_index] /= max(CAVITY_HEAT_CAPACITY_v24, eps(Float64))
    return nothing
end

function receiver_ode_v24!(dTs, Ts, context, time)
    receiver_rhs_v24!(
        dTs, Ts, time, context.p, context.operating,
        context.z, context.dx, context.core_weights, context.perim_weights,
        context.z_rear, context.dx_rear,
        context.Tg, context.Qgas, context.hcell, context.Qgas_rear, context.hrear,
    )
    return nothing
end

function simulate_v24(p, operating, save_times;
                      initial_temperature=Tamb, nodes=default_nodes,
                      rear_nodes=REAR_TUBE_DEFAULT_NODES_v24,
                      initial_perim_temperature=nothing,
                      initial_rear_temperature=nothing,
                      initial_cavity_temperature=Tamb,
                      solver=Rodas5P(autodiff=AutoFiniteDiff()),
                      reltol=1e-6, abstol=1e-7, dtmax=30.0)
    length(p) == 21 || throw(ArgumentError("1D_v24 expects 21 parameters"))
    all(isfinite, p) || throw(ArgumentError("all model parameters must be finite"))
    all(lb_full_v24[i] <= p[i] <= ub_full_v24[i] for i in eachindex(p)) ||
        throw(ArgumentError("1D_v24 parameters are outside declared bounds"))

    time = Float64.(save_times)
    dx = L / nodes
    dx_rear = REAR_TUBE_LENGTH_v24 / rear_nodes
    z = collect(range(dx / 2, step=dx, length=nodes))
    z_rear = collect(range(dx_rear / 2, step=dx_rear, length=rear_nodes))
    zg = vcat(
        collect(range(0.0, L, length=nodes + 1)),
        L .+ collect(range(dx_rear, step=dx_rear, length=rear_nodes)),
    )
    core_weights = solar_weights_v5(p[13], p[6], nodes)
    perim_weights = axial_exponential_weights_v24(perimeter_source_attenuation_v24(p), nodes, dx)

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
    u_initial = vcat(Tcore_initial, Tperim_initial, Ttube_initial, Float64(initial_cavity_temperature))

    Tg_history = Matrix{Float64}(undef, nodes + rear_nodes + 1, length(time))
    h_history = Matrix{Float64}(undef, nodes, length(time))
    nu_history = Matrix{Float64}(undef, nodes, length(time))
    effectiveness_history = Matrix{Float64}(undef, nodes, length(time))
    qgas_history = Matrix{Float64}(undef, nodes, length(time))
    qgas_rear_history = Matrix{Float64}(undef, rear_nodes, length(time))
    active_fraction_history = Vector{Float64}(undef, length(time))
    core_temperature_history = Matrix{Float64}(undef, nodes, length(time))
    perim_temperature_history = Matrix{Float64}(undef, nodes, length(time))
    rear_temperature_history = Matrix{Float64}(undef, rear_nodes, length(time))
    cavity_temperature_history = Vector{Float64}(undef, length(time))

    Tg = zeros(nodes + rear_nodes + 1)
    Qgas = zeros(nodes)
    hcell = zeros(nodes)
    Qgas_rear = zeros(rear_nodes)
    hrear = zeros(rear_nodes)

    context = (
        p=p, operating=operating, z=z, dx=dx,
        core_weights=core_weights, perim_weights=perim_weights,
        z_rear=z_rear, dx_rear=dx_rear,
        Tg=Tg, Qgas=Qgas, hcell=hcell, Qgas_rear=Qgas_rear, hrear=hrear,
    )
    problem = ODEProblem(receiver_ode_v24!, u_initial, (time[1], time[end]), context)
    solution = solve(
        problem, solver; saveat=time, save_everystep=false, dense=false,
        reltol=reltol, abstol=abstol, dtmax=dtmax, tstops=time,
    )
    successful_retcode(solution) || error("1D_v24 ODE solve failed with retcode $(solution.retcode)")

    state_history = reduce(hcat, solution.u)
    core_temperature_history .= state_history[1:nodes, :]
    perim_temperature_history .= state_history[(nodes + 1):(2 * nodes), :]
    tube_start = 2 * nodes + 1
    tube_stop = 2 * nodes + rear_nodes
    cavity_index = 2 * nodes + rear_nodes + 1
    rear_temperature_history .= state_history[tube_start:tube_stop, :]
    cavity_temperature_history .= state_history[cavity_index, :]

    nu_history = Matrix{Float64}(undef, nodes, length(time))
    effectiveness_history = Matrix{Float64}(undef, nodes, length(time))

    for output_index in eachindex(time)
        t = time[output_index]
        Tcore_now = view(core_temperature_history, :, output_index)
        Ttube_now = view(rear_temperature_history, :, output_index)
        gas_profile_v24!(Tg, Qgas, hcell, Qgas_rear, hrear, Tcore_now, Ttube_now,
                         t, p, operating, z, dx, z_rear, dx_rear)
        Tg_history[:, output_index] .= Tg
        h_history[:, output_index] .= hcell
        qgas_history[:, output_index] .= Qgas
        qgas_rear_history[:, output_index] .= Qgas_rear

        flow = max(0.0, operating.flow_lpm(t))
        Tin = operating.inlet_temperature(t)
        mdot_total = m_dot_standard_v24(flow)
        phi_act = receiver_active_flow_fraction_v24(flow, Tin, p)
        active_fraction_history[output_index] = phi_act
        for i in 1:nodes
            Tfilm = 0.5 * (Tcore_now[i] + Tg_history[i, output_index])
            cp = cpf_f(Tfilm)
            kf = kf_f(Tfilm)
            nu_history[i, output_index] = h_history[i, output_index] * Dh / max(kf, eps(Float64))
            if flow > 1e-12
                mdot_active = phi_act * mdot_total
                UA = h_history[i, output_index] * P_exchange * dx
                effectiveness_history[i, output_index] = -expm1(-UA / (mdot_active * cp))
            else
                effectiveness_history[i, output_index] = 1.0
            end
        end
    end

    return (
        time=time, z_solid=z, z_rear_tube=L .+ z_rear, z_gas=zg,
        core_temperature=core_temperature_history,
        perim_temperature=perim_temperature_history,
        active_temperature=core_temperature_history,
        wall_temperature=perim_temperature_history,
        solid_temperature=perim_temperature_history,
        rear_tube_temperature=rear_temperature_history,
        cavity_temperature=cavity_temperature_history,
        boundary_temperature=vec(perim_temperature_history[1, :]),
        gas_temperature=Tg_history,
        heat_transfer_coefficient=h_history,
        receiver_nusselt=nu_history,
        receiver_effectiveness=effectiveness_history,
        receiver_gas_heat=qgas_history,
        rear_tube_gas_heat=qgas_rear_history,
        active_flow_fraction=active_fraction_history,
        core_source_weights=core_weights,
        perim_source_weights=perim_weights,
        ode_solution=solution,
    )
end

begin # v24 experimental interface
    const WALL_CHAIN_POSITIONS_v24 = (
        T8=sensor_positions[:T8],
        T12=sensor_positions[:T9],
        T11=sensor_positions[:T10],
    )

    function measured_initial_perim_profile_v24(data, simulation_id)
        z8 = WALL_CHAIN_POSITIONS_v24.T8
        z12 = WALL_CHAIN_POSITIONS_v24.T12
        z11 = WALL_CHAIN_POSITIONS_v24.T11
        T8 = observation(data, simulation_id, "_T8")[1]
        T12 = observation(data, simulation_id, "_T12")[1]
        T11 = observation(data, simulation_id, "_T11")[1]
        return function (z)
            z <= z8 && return T8
            z >= z11 && return T11
            if z <= z12
                return T8 + (T12 - T8) * (z - z8) / (z12 - z8)
            end
            return T12 + (T11 - T12) * (z - z12) / (z11 - z12)
        end
    end

    function measured_initial_core_profile_v24(data, simulation_id)
        z8 = sensor_positions[:T8]
        z9 = sensor_positions[:T9]
        z10 = sensor_positions[:T10]
        T8 = observation(data, simulation_id, "_T8")[1]
        T9 = observation(data, simulation_id, "_T9")[1]
        T10 = observation(data, simulation_id, "_T10")[1]
        return function (z)
            z <= z8 && return T8
            z >= z10 && return T10
            if z <= z9
                return T8 + (T9 - T8) * (z - z8) / (z9 - z8)
            end
            return T9 + (T10 - T9) * (z - z9) / (z10 - z9)
        end
    end

    function measured_initial_rear_profile_v24(data, simulation_id, p)
        core_profile = measured_initial_core_profile_v24(data, simulation_id)
        perim_profile = measured_initial_perim_profile_v24(data, simulation_id)
        f_core_tube = core_tube_fraction_v24(p)
        T_receiver_exit =
            f_core_tube * core_profile(L) + (1.0 - f_core_tube) * perim_profile(L)
        T_tube_exit = observation(data, simulation_id, "_Tf")[1]
        return function (z_rear)
            fraction = clamp(z_rear / REAR_TUBE_LENGTH_v24, 0.0, 1.0)
            return T_receiver_exit + fraction * (T_tube_exit - T_receiver_exit)
        end
    end

    function experimental_case_v24(simulation_id; is_cooling=false)
        data = is_cooling ? measurements_cooling : measurements
        conditions = is_cooling ? simulation_conditions_cooling : simulation_conditions
        time = observation_time(data, simulation_id)
        experiment = hcat(
            observation(data, simulation_id, "_T8"),
            observation(data, simulation_id, "_T12"),
            observation(data, simulation_id, "_T11"),
            observation(data, simulation_id, "_T9"),
            observation(data, simulation_id, "_T10"),
            observation(data, simulation_id, "_Tf"),
            observation(data, simulation_id, "_T2"),
        )
        return conditions[simulation_id], time, experiment, data
    end

    function gas_at_position_v24(result, position)
        z = result.z_gas
        values = result.gas_temperature
        position <= z[1] && return vec(values[1, :])
        position >= z[end] && return vec(values[end, :])
        i = searchsortedlast(z, position)
        fraction = (position - z[i]) / (z[i + 1] - z[i])
        return vec((1.0 - fraction) .* values[i, :] .+ fraction .* values[i + 1, :])
    end

    function core_at_v24(result, position)
        z = result.z_solid
        values = result.core_temperature
        position <= z[1] && return vec(values[1, :])
        position >= z[end] && return vec(values[end, :])
        i = searchsortedlast(z, position)
        fraction = (position - z[i]) / (z[i + 1] - z[i])
        return vec((1.0 - fraction) .* values[i, :] .+ fraction .* values[i + 1, :])
    end

    function perim_at_v24(result, position)
        z = result.z_solid
        values = result.perim_temperature
        position <= z[1] && return vec(values[1, :])
        position >= z[end] && return vec(values[end, :])
        i = searchsortedlast(z, position)
        fraction = (position - z[i]) / (z[i + 1] - z[i])
        return vec((1.0 - fraction) .* values[i, :] .+ fraction .* values[i + 1, :])
    end

    function extract_outputs_v24(result, p)
        T8 = perim_at_v24(result, sensor_positions[:T8])
        T12_perim = perim_at_v24(result, WALL_CHAIN_POSITIONS_v24.T12)
        T11_perim = perim_at_v24(result, WALL_CHAIN_POSITIONS_v24.T11)
        T9_core = core_at_v24(result, sensor_positions[:T9])
        T10_core = core_at_v24(result, sensor_positions[:T10])
        Tf = gas_at_position_v24(result, T3_SAMPLE_POSITION_v24)
        T2 = result.cavity_temperature
        h_mean = vec(mean(result.heat_transfer_coefficient, dims=1))
        return hcat(T8, T12_perim, T11_perim, T9_core, T10_core, Tf, T2, h_mean)
    end

    function solve_case_v24(p, simulation_id; is_cooling=false,
                           nodes=default_nodes,
                           rear_nodes=REAR_TUBE_DEFAULT_NODES_v24,
                           solver=Rodas5P(autodiff=AutoFiniteDiff()),
                           reltol=1e-6, abstol=1e-7, dtmax=30.0)
        conditions, time, experiment, data = experimental_case_v24(
            simulation_id; is_cooling=is_cooling
        )
        flow = observation(data, simulation_id, "_flow")
        Tin = observation(data, simulation_id, "_Tin")
        ambient = observation(data, simulation_id, "_Tamb")
        T2 = observation(data, simulation_id, "_T2")
        irradiance = is_cooling ? zeros(length(time)) : fill(conditions[Io], length(time))
        operating = operating_condition_v3(
            irradiance=linear_history_v3(time, irradiance),
            flow_lpm=linear_history_v3(time, flow),
            inlet_temperature=linear_history_v3(time, Tin),
            ambient_temperature=linear_history_v3(time, ambient),
        )
        result = simulate_v24(
            p, operating, time;
            initial_temperature=measured_initial_core_profile_v24(data, simulation_id),
            initial_perim_temperature=measured_initial_perim_profile_v24(data, simulation_id),
            initial_rear_temperature=measured_initial_rear_profile_v24(data, simulation_id, p),
            initial_cavity_temperature=T2[1],
            nodes=nodes, rear_nodes=rear_nodes, solver=solver,
            reltol=reltol, abstol=abstol, dtmax=dtmax,
        )
        outputs = extract_outputs_v24(result, p)
        return outputs, result, experiment
    end
end

begin # v24 objective helpers
    const ORDERING_SCALE_v24 = 100.0
    const ORDERING_WEIGHT_v24 = 1.00
    const TIMING_WEIGHT_v24 = 0.15
    const SLOPE_WEIGHT_v24 = 0.25
    const COOLING_UPTURN_WEIGHT_v24 = 10000.0

    function normalized_slope_mse_v24(model, experiment)
        length(model) < 2 && return 0.0
        scale = max(maximum(experiment) - minimum(experiment), 1.0)
        return mean(abs2, diff(model) .- diff(experiment)) / scale^2
    end

    function timing_penalty_v24(time, model, experiment)
        duration = max(time[end] - time[1], 1.0)
        return ((get_t90_v3(time, model) - get_t90_v3(time, experiment)) / duration)^2
    end

    function signal_loss_v24(time, model, experiment)
        return normalized_mse_v3(model, experiment) +
               SLOPE_WEIGHT_v24 * normalized_slope_mse_v24(model, experiment) +
               TIMING_WEIGHT_v24 * timing_penalty_v24(time, model, experiment)
    end

    function cooling_upturn_loss_v24(model)
        target_columns = (3, 5, 6) # T11, T10, T3
        total = 0.0
        for j in target_columns
            scale = max(maximum(model[:, j]) - minimum(model[:, j]), 20.0)
            upturn = max.(diff(model[:, j]), 0.0) ./ scale
            total += maximum(upturn)^2 + mean(abs2, upturn)
        end
        return total / length(target_columns)
    end

    function ordering_loss_v24(model, experiment)
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
        return mean(abs2, collect(model_offsets) .- collect(experiment_offsets)) / ORDERING_SCALE_v24^2
    end

    function capacitance_regularization_v24(p; nodes=15)
        capacitance = participating_assembly_heat_capacity_v24(p, nodes)
        return CAPACITANCE_REGULARIZATION_WEIGHT_v24 *
               ((capacitance - MEASURED_ASSEMBLY_CAPACITANCE_v24) / MEASURED_ASSEMBLY_CAPACITANCE_v24)^2
    end

    function loss_cases_v24(p, keys; is_cooling=false, nodes=15)
        total = 0.0
        for simulation_id in keys
            try
                model, result, experiment = solve_case_v24(
                    p, simulation_id; is_cooling=is_cooling, nodes=nodes,
                    reltol=2e-5, abstol=2e-6, dtmax=60.0,
                )
                all(isfinite, model) || return Inf
                total += mean(
                    signal_loss_v24(result.time, model[:, j], experiment[:, j])
                    for j in 1:7
                )
                if !is_cooling
                    total += ORDERING_WEIGHT_v24 * ordering_loss_v24(model, experiment)
                else
                    total += COOLING_UPTURN_WEIGHT_v24 * cooling_upturn_loss_v24(model)
                end
            catch err
                err isa InterruptException && rethrow()
                return Inf
            end
        end
        return total / length(keys)
    end

    function loss_heating_v24(p, keys=sim_key_heat; nodes=15)
        regularization = capacitance_regularization_v24(p; nodes=nodes)
        return loss_cases_v24(p, keys; is_cooling=false, nodes=nodes) + regularization
    end

    function loss_cooling_v24(p, keys=sim_key_cool; nodes=15)
        regularization = capacitance_regularization_v24(p; nodes=nodes)
        return loss_cases_v24(p, keys; is_cooling=true, nodes=nodes) + regularization
    end

    function calibrate_v24(; nodes=15, maximum_iterations=100, maximum_time_seconds=360.0)
        p0 = copy(pnew_v24)
        idx_fit = fit_heat_transfer_indices_v24
        function obj_free(theta)
            p_full = copy(p0)
            p_full[idx_fit] .= theta
            return loss_heating_v24(p_full; nodes=nodes) + loss_cooling_v24(p_full; nodes=nodes)
        end
        res = optimize_with_nlopt_v3(
            obj_free,
            pnew_v24[idx_fit],
            lb_full_v24[idx_fit],
            ub_full_v24[idx_fit];
            maximum_iterations=maximum_iterations,
            maximum_time_seconds=maximum_time_seconds,
            algorithm=NLopt.LN_BOBYQA(),
            label="fit-v24",
        )
        p_opt = copy(p0)
        p_opt[idx_fit] .= res.minimizer
        return (objective=res.minimum, parameters=p_opt, minimizer=res.minimizer, retcode=res.retcode)
    end
end




