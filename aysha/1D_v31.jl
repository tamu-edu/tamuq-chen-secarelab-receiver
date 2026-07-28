# ============================================================================
# 1D_v31.jl - Energy-accounting 2-Zone Core/Perimeter Macro Model
# ============================================================================
# v31 keeps the v30 distributed rear-contact rail, but turns it into a
# manuscript-readiness diagnostic: the perimeter and rear-rail heat capacities
# are constrained so the model cannot improve by shrinking thermal inventory
# below the measured assembly scale. The default calibration stage fits only
# rear/flange conductances before any transport/source parameters are released.
#
# Key Physical Features of v31:
#   1. Core Channel Zone (T_core): 100-channel matrix. Receives direct central
#      solar beam fraction. Exchanges heat with the active gas stream using
#      bounded laminar developing-flow convection. Compares to T9 and T10.
#   2. Perimeter Housing Zone (T_perim): Side walls, alumina felt insulation,
#      and front housing. Receives a fitted fraction of the scaled incident
#      receiver power, guided by the measured flux-map shape but not by the
#      flux-map absolute magnitude. Compares to T8,
#      T12, T11 and drives the T2/cavity state.
#   3. Distributed Rear/Adaptor Rail (T_rear[i]): bounded effective thermal
#      inventory and axial contact network near the receiver rear/adaptor. It
#      replaces v27's direct rear sink and v29's single terminal reservoir.
#   4. Rear Tube/Flange Sink: downstream loss is carried by the explicit rear
#      tube and water-flange geometry. Cooling can strengthen the flange path
#      smoothly after lamp shutoff.
#
# Parameter Vector p[1:25] in v31:
#   p[1]  A_Nu          laminar developing Nu prefactor
#   p[2]  B_Re          laminar developing Reynolds exponent (bounded 0.0 to 0.6)
#   p[3]  C_Pr          fixed Prandtl exponent (= 1/3)
#   p[4]  phi_0         active flow fraction at Re = 50
#   p[5]  m_rec         flow recruitment exponent
#   p[6]  front_dep     frozen v12a front-cell deposition fraction (= 1.0 or 0.0)
#   p[7]  scale_456     incident/participating flux scale at 456 kW/m2
#   p[8]  scale_304     incident/participating flux scale at 304 kW/m2
#   p[9]  scale_256     incident/participating flux scale at 256 kW/m2
#   p[10] G_core_perim  radial core-perimeter conductance per length [W/m/K]
#   p[11] C_perim_eff   perimeter/housing participating heat capacity [J/K]
#   p[12] k_perim_ref   perimeter axial conductivity at 900 K [W/m/K]
#   p[13] beta_opt      frozen v12a optical attenuation [1/m]
#   p[14] spill_capture captured fraction of off-receiver flux spillover
#   p[15] beta_perim    axial attenuation of perimeter/rim source [1/m]
#   p[16] f_core_rear   fraction of receiver-to-rear coupling from core
#   p[17] flange_scale  multiplier for rear-tube-to-water-flange conductance
#   p[18] flange_cool_gain cooling-only gain on explicit flange conductance
#   p[19] flange_cool_tau  lamp-off flange-gain time constant [s]
#   p[20] k_core_axial_scale effective axial solid-conduction multiplier
#   p[21] C_rear_eff    rear/adaptor rail participating heat capacity [J/K]
#   p[22] G_recv_rear   receiver-to-rear-rail conductance [W/K]
#   p[23] G_rear_tube   rear-rail exit to tube-inlet conductance [W/K]
#   p[24] G_rear_cavity rear-rail to cavity conductance [W/K]
#   p[25] G_rear_axial  rear-rail axial conductance per link [W/K]
# ============================================================================

include("1D_v5.jl")

begin # v31 fixed constants
    const EPS_FRONT_FIXED_v31 = EPS_FRONT_FIXED_V5
    const ETA_ABS_FIXED_v31 = 1.0
    const BETA_OPT_FIXED_v31 = BETA_OPT_FIXED_V5
    const H_FRONT_FIXED_v31 = 10.0
    const PRANDTL_EXPONENT_FIXED_v31 = 1.0 / 3.0
    const STANDARD_AIR_DENSITY_v31 = 101325.0 / (287.05 * 294.25)
    const REYNOLDS_REFERENCE_v31 = 50.0
    const NU_FLOOR_v31 = 3.61
    const MEASURED_ASSEMBLY_CAPACITANCE_v31 = 301.0
    const IRRADIANCE_LEVELS_v31 = (456000.0, 304000.0, 256000.0)
    const T3_SAMPLE_POSITION_v31 = 140.0e-3
    const COOLING_ROOM_TEMPERATURE_v31 = 17.0 + 273.15
    const CAPACITANCE_REGULARIZATION_WEIGHT_v31 = 1.00
    const FLUX_SIGMA_FIXED_v31 = 30.0e-3
    const FLUX_INTEGRATION_POINTS_v31 = 161
    const IRRADIANCE_GATE_WIDTH_v31 = 1000.0

    # Cavity and housing geometry
    const CAVITY_OUTER_RADIUS_v31 = 75.0e-3
    const METAL_THICKNESS_v31 = 18.0e-3
    const INSULATION_OUTER_RADIUS_v31 = CAVITY_OUTER_RADIUS_v31 - METAL_THICKNESS_v31
    const METAL_OUTER_RADIUS_v31 = CAVITY_OUTER_RADIUS_v31
    const CAVITY_LENGTH_v31 = 165.0e-3

    const ADAPTOR_DIAMETER_v31 = 77.6e-3
    const ADAPTOR_RADIUS_v31 = ADAPTOR_DIAMETER_v31 / 2.0
    const ADAPTOR_LENGTH_v31 = 57.0e-3
    const ADAPTOR_TUBE_RADIUS_v31 = 6.5e-3
    const ADAPTOR_RECEIVER_OVERLAP_LENGTH_v31 = 28.0e-3
    const ADAPTOR_TUBE_OVERLAP_LENGTH_v31 = ADAPTOR_LENGTH_v31 - ADAPTOR_RECEIVER_OVERLAP_LENGTH_v31
    const ADAPTOR_OVERLAP_LENGTH_v31 = ADAPTOR_RECEIVER_OVERLAP_LENGTH_v31

    const RECEIVER_EQ_RADIUS_v31 = sqrt(A_frt / pi)
    const RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v31 =
        4.0 * pi * 0.08 / log(INSULATION_OUTER_RADIUS_v31 / RECEIVER_EQ_RADIUS_v31)
    const ADAPTOR_CONTACT_RESISTANCE_v31 =
        1.0 / (100.0 * 4.0 * w_t * ADAPTOR_OVERLAP_LENGTH_v31)
    const ADAPTOR_CONTACT_CONDUCTANCE_REFERENCE_v31 = 1.0 / ADAPTOR_CONTACT_RESISTANCE_v31

    # Rear tube geometry
    const REAR_TUBE_LENGTH_v31 = 150.0e-3
    const REAR_TUBE_CAVITY_LENGTH_v31 = 46.0e-3
    const REAR_TUBE_DEFAULT_NODES_v31 = 12
    const REAR_TUBE_GAS_RADIUS_v31 = ADAPTOR_TUBE_RADIUS_v31
    const REAR_TUBE_WALL_THICKNESS_v31 = 1.5e-3
    const REAR_TUBE_OUTER_RADIUS_v31 = REAR_TUBE_GAS_RADIUS_v31 + REAR_TUBE_WALL_THICKNESS_v31
    const REAR_TUBE_FLOW_AREA_v31 = pi * REAR_TUBE_GAS_RADIUS_v31^2
    const REAR_TUBE_HYDRAULIC_DIAMETER_v31 = 2.0 * REAR_TUBE_GAS_RADIUS_v31
    const REAR_TUBE_EXCHANGE_PERIMETER_v31 = 2.0 * pi * REAR_TUBE_GAS_RADIUS_v31

    const ALUMINA_DENSITY_v31 = 3900.0
    const ALUMINUM_DENSITY_v31 = 2700.0
    const ALUMINUM_HEAT_CAPACITY_v31 = 900.0
    const ALUMINUM_CONDUCTIVITY_v31 = 205.0
    const METAL_EMISSIVITY_v31 = 0.2
    const H_NAT_CAVITY_v31 = 10.0
    const INSULATION_DENSITY_v31 = 140.0
    const INSULATION_HEAT_CAPACITY_v31 = 1360.0
    const WATER_FLANGE_TEMPERATURE_v31 = 25.0 + 273.15

    const CAVITY_FELT_VOLUME_v31 = max(
        0.0,
        pi * INSULATION_OUTER_RADIUS_v31^2 * CAVITY_LENGTH_v31 -
        A_frt * L - pi * ADAPTOR_RADIUS_v31^2 * ADAPTOR_LENGTH_v31,
    )
    const CAVITY_METAL_VOLUME_v31 =
        pi * (METAL_OUTER_RADIUS_v31^2 - INSULATION_OUTER_RADIUS_v31^2) * CAVITY_LENGTH_v31 +
        pi * METAL_OUTER_RADIUS_v31^2 * METAL_THICKNESS_v31
    const CAVITY_OUTER_AREA_v31 =
        2.0 * pi * METAL_OUTER_RADIUS_v31 * CAVITY_LENGTH_v31 + pi * METAL_OUTER_RADIUS_v31^2
    const CAVITY_HEAT_CAPACITY_v31 =
        INSULATION_DENSITY_v31 * INSULATION_HEAT_CAPACITY_v31 * CAVITY_FELT_VOLUME_v31 +
        ALUMINUM_DENSITY_v31 * ALUMINUM_HEAT_CAPACITY_v31 * CAVITY_METAL_VOLUME_v31
end

begin # v31 non-uniform front-flux partition
    transverse_flux_shape_v31(x, y; sigma=FLUX_SIGMA_FIXED_v31) =
        exp(-0.5 * (x^2 + y^2) / sigma^2)

    function flux_zone_fractions_v31(; radius=CAVITY_OUTER_RADIUS_v31,
                                     receiver_half_width=w_t / 2.0,
                                     points=FLUX_INTEGRATION_POINTS_v31)
        dx_flux = 2.0 * radius / (points - 1)
        receiver_power = 0.0
        aperture_power = 0.0
        for ix in 1:points
            x = -radius + (ix - 1) * dx_flux
            for iy in 1:points
                y = -radius + (iy - 1) * dx_flux
                x^2 + y^2 > radius^2 && continue
                weight = transverse_flux_shape_v31(x, y) * dx_flux^2
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

    const FLUX_ZONE_FRACTIONS_FIXED_v31 = flux_zone_fractions_v31()
end

# Default v31 parameter vector p[1:25]
begin # v31 parameter values and bounds
    pnew_v31 = [
        5.022922293757817,  # p[1]  A_Nu (laminar developing prefactor; v27 seed)
        0.5133267688838168, # p[2]  B_Re (Reynolds exponent bounded 0.0 to 0.6)
        PRANDTL_EXPONENT_FIXED_v31, # p[3]  C_Pr (= 1/3)
        1.0,                # p[4]  phi_0 (fixed full active flow fraction)
        0.0,                # p[5]  m_rec (flow recruitment exponent)
        1.0,                # p[6]  front_dep (frozen v12a front-cell deposition)
        1.9811828483104401, # p[7]  scale_456
        1.9873720312526904, # p[8]  scale_304
        0.9481818340413103, # p[9]  scale_256
        15.914856465811004, # p[10] G_core_perim (radial core-perimeter conductance [W/m/K])
        150.0,              # p[11] C_perim_eff (perimeter participating capacity [J/K])
        5.045319638108573,  # p[12] k_perim_ref (perimeter axial conductivity at 900K [W/m/K])
        184.67243519237965, # p[13] beta_opt [1/m]
        0.5260985841287006, # p[14] spill_capture (captured spillover fraction)
        3.627912152679222,  # p[15] beta_perim (perimeter source attenuation [1/m])
        0.9993134300665368, # p[16] f_core_rear (receiver-rear coupling core fraction)
        0.10138219425614439,# p[17] flange_scale (rear/flange heat-removal multiplier)
        6.492594832657662,  # p[18] flange_cool_gain (smooth lamp-off gain)
        116.51093421165298, # p[19] flange_cool_tau [s]
        0.02787105475540781,# p[20] k_core_axial_scale
        80.0,               # p[21] C_rear_eff (rear/adaptor participating capacity [J/K])
        1.9355656026962773, # p[22] G_recv_rear (receiver-to-rear conductance [W/K])
        0.22252213606434484,# p[23] G_rear_tube (rear rail to tube conductance [W/K])
        0.8133555547684154, # p[24] G_rear_cavity (rear rail to cavity conductance [W/K])
        5.0,                # p[25] G_rear_axial (rear rail axial conductance per link [W/K])
    ]
    lb_full_v31 = [
        0.01,
        0.0,
        pnew_v31[3],
        pnew_v31[4],
        pnew_v31[5],
        pnew_v31[6],
        0.02,
        0.02,
        0.02,
        2.0,
        150.0,
        0.0,
        pnew_v31[13],
        0.0,
        0.0,
        0.0,
        0.10,
        0.0,
        1.0,
        0.0,
        80.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ]
    ub_full_v31 = [
        25.00,
        0.60,
        pnew_v31[3],
        pnew_v31[4],
        pnew_v31[5],
        pnew_v31[6],
        5.00,
        5.00,
        5.00,
        80.0,
        230.0,
        50.0,
        pnew_v31[13],
        1.00,
        300.0,
        1.0,
        20.0,
        50.0,
        2000.0,
        0.50,
        150.0,
        50.0,
        50.0,
        10.0,
        100.0,
    ]

    fit_rear_stage_indices_v31 = [17, 18, 19, 22, 23, 24, 25]
    fit_transport_stage_indices_v31 = [10, 12, 17, 18, 19, 20, 22, 23, 24, 25]
    fit_full_stage_indices_v31 = [1, 2, 7, 8, 9, 10, 12, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25]
    fit_heat_transfer_indices_v31 = fit_rear_stage_indices_v31
    fit_source_indices_v31 = [14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25]
    fit_power_scale_indices_v31 = [7, 8, 9]
    fit_radiation_indices_v31 = Int[]
end

function irradiance_level_index_v31(irradiance)
    irradiance <= 0.0 && return 0
    differences = abs.(collect(IRRADIANCE_LEVELS_v31) .- Float64(irradiance))
    return argmin(differences)
end

power_scale_index_v31(irradiance) = begin
    level_index = irradiance_level_index_v31(irradiance)
    level_index == 0 ? 0 : 6 + level_index
end

function absorbed_power_scale_v31(irradiance, p)
    index = power_scale_index_v31(irradiance)
    return index == 0 ? 1.0 : p[index]
end

core_perimeter_conductance_per_length_v31(p) = max(p[10], eps(Float64))

perimeter_heat_capacity_total_v31(p) = max(p[11], eps(Float64))

core_receiver_heat_capacity_v31(nodes) =
    sum(rho_s * Cps_f(900.0) * A_solid * (L / nodes) for _ in 1:nodes)

participating_assembly_heat_capacity_v31(p, nodes=default_nodes) =
    core_receiver_heat_capacity_v31(nodes) + perimeter_heat_capacity_total_v31(p)

rear_reservoir_heat_capacity_v31(p) = clamp(p[21], 50.0, 500.0)

participating_total_heat_capacity_v31(p, nodes=default_nodes) =
    participating_assembly_heat_capacity_v31(p, nodes) + rear_reservoir_heat_capacity_v31(p)

perimeter_axial_conductivity_v31(T, p) =
    max(p[12], 0.0) * (property_temperature(T) / 900.0)^3

flux_receiver_fraction_v31() = FLUX_ZONE_FRACTIONS_FIXED_v31.receiver

flux_spillover_fraction_v31() = FLUX_ZONE_FRACTIONS_FIXED_v31.spillover

perimeter_spillover_capture_v31(p) = clamp(p[14], 0.0, 1.0)

perimeter_source_fraction_v31(p) =
    clamp(perimeter_spillover_capture_v31(p) * flux_spillover_fraction_v31(), 0.0, 0.80)

modeled_absorbed_receiver_power_v31(irradiance, p) =
    ETA_ABS_FIXED_v31 * absorbed_power_scale_v31(irradiance, p) *
    max(0.0, irradiance) * A_frt

core_absorbed_power_v31(irradiance, p) =
    modeled_absorbed_receiver_power_v31(irradiance, p) *
    (1.0 - perimeter_source_fraction_v31(p))

perimeter_absorbed_power_v31(irradiance, p) =
    modeled_absorbed_receiver_power_v31(irradiance, p) *
    perimeter_source_fraction_v31(p)

modeled_participating_absorbed_power_v31(irradiance, p) =
    modeled_absorbed_receiver_power_v31(irradiance, p)

perimeter_source_attenuation_v31(p) = clamp(p[15], 0.0, 300.0)

rear_core_fraction_v31(p) = clamp(p[16], 0.0, 1.0)

core_tube_fraction_v31(p) = rear_core_fraction_v31(p)

flange_loss_scale_v31(p) = clamp(p[17], 0.10, 20.0)

flange_cooling_gain_v31(p) = clamp(p[18], 0.0, 50.0)

flange_cooling_time_constant_v31(p) = clamp(p[19], 1.0, 2000.0)

core_axial_conduction_scale_v31(p) = clamp(p[20], 0.0, 0.50)

receiver_rear_conductance_v31(p) = clamp(p[22], 0.0, 50.0)

rear_tube_conductance_v31(p) = clamp(p[23], 0.0, 50.0)

rear_cavity_conductance_v31(p) = clamp(p[24], 0.0, 10.0)

rear_axial_conductance_v31(p) = clamp(p[25], 0.0, 100.0)

function rear_contact_weights_v31(z)
    z0 = sensor_positions[:T8]
    raw = similar(Float64.(z))
    for i in eachindex(z)
        s = clamp((z[i] - z0) / max(L - z0, eps(Float64)), 0.0, 1.0)
        raw[i] = 0.02 + 0.98 * s
    end
    total = sum(raw)
    total <= eps(Float64) && return fill(1.0 / length(z), length(z))
    return raw ./ total
end

function lamp_off_gate_v31(time, irradiance, p)
    irradiance_factor = 1.0 / (1.0 + (max(0.0, irradiance) / IRRADIANCE_GATE_WIDTH_v31)^4)
    time_factor = 1.0 - exp(-max(0.0, time) / flange_cooling_time_constant_v31(p))
    return clamp(irradiance_factor * time_factor, 0.0, 1.0)
end

effective_flange_loss_scale_v31(p, time, irradiance) =
    flange_loss_scale_v31(p) * (1.0 + flange_cooling_gain_v31(p) * lamp_off_gate_v31(time, irradiance, p))

function axial_exponential_weights_v31(beta, nodes, dx)
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

function active_flow_fraction_v31(reynolds, p)
    phi0 = clamp(p[4], 0.1, 1.0)
    m = max(p[5], 0.0)
    ratio = max(reynolds, eps(Float64)) / REYNOLDS_REFERENCE_v31
    return clamp(phi0 * ratio^m, 0.1, 1.0)
end

function receiver_active_flow_fraction_v31(flow_lpm, Tin, p)
    flow_lpm <= 1e-12 && return 1.0
    mdot = m_dot_standard_v31(flow_lpm)
    mu = mu_f_f(Tin)
    reynolds = mdot * Dh / (A_flow * mu)
    return active_flow_fraction_v31(reynolds, p)
end

function receiver_channel_nusselt_v31(z, reynolds, prandtl, p)
    A_nu = max(p[1], 0.0)
    B_re = clamp(p[2], 0.0, 0.6)
    C_pr = p[3]
    entry_factor = (Dh / max(z, Dh / 2.0))^(1.0 / 3.0)
    Nu_dev = A_nu * max(reynolds, eps(Float64))^B_re * max(prandtl, eps(Float64))^C_pr * entry_factor
    return max(NU_FLOOR_v31, Nu_dev)
end

m_dot_standard_v31(flow_lpm) =
    STANDARD_AIR_DENSITY_v31 * max(0.0, flow_lpm) / 60000.0

alumina_conductivity_v31(T) =
    5.5 + 34.5 * exp(-0.0033 * (property_temperature(T) - 273.15))

alumina_heat_capacity_v31(T) =
    max(500.0, (1.00446 + 1.742e-4 * property_temperature(T) -
                2.796e4 * property_temperature(T)^-2) * 1000.0)

function rear_tube_outer_radius_v31(z_rear)
    return z_rear <= ADAPTOR_TUBE_OVERLAP_LENGTH_v31 ?
           ADAPTOR_RADIUS_v31 : REAR_TUBE_OUTER_RADIUS_v31
end

function rear_tube_solid_area_v31(z_rear)
    outer_radius = rear_tube_outer_radius_v31(z_rear)
    return pi * (outer_radius^2 - REAR_TUBE_GAS_RADIUS_v31^2)
end

function rear_tube_cavity_conductance_per_length_v31(z_rear, p)
    z_rear > REAR_TUBE_CAVITY_LENGTH_v31 && return 0.0
    outer_radius = min(rear_tube_outer_radius_v31(z_rear), 0.98 * INSULATION_OUTER_RADIUS_v31)
    return 2.0 * pi * 0.08 / log(INSULATION_OUTER_RADIUS_v31 / outer_radius)
end

function rear_tube_flange_conductance_per_length_v31(z_rear, Ttube, p, time, irradiance)
    z_rear <= REAR_TUBE_CAVITY_LENGTH_v31 && return 0.0
    tube_wall_resistance =
        log(REAR_TUBE_OUTER_RADIUS_v31 / REAR_TUBE_GAS_RADIUS_v31) /
        (2.0 * pi * alumina_conductivity_v31(Ttube))
    aluminum_resistance =
        log(METAL_OUTER_RADIUS_v31 / REAR_TUBE_OUTER_RADIUS_v31) /
        (2.0 * pi * ALUMINUM_CONDUCTIVITY_v31)
    return effective_flange_loss_scale_v31(p, time, irradiance) /
           (tube_wall_resistance + aluminum_resistance)
end

function cavity_loss_to_ambient_v31(Tcavity, ambient)
    return H_NAT_CAVITY_v31 * CAVITY_OUTER_AREA_v31 * (Tcavity - ambient) +
           METAL_EMISSIVITY_v31 * sigma_sb * CAVITY_OUTER_AREA_v31 *
           (Tcavity^4 - ambient^4)
end

function rear_tube_capacity_v31(Ttube, z_rear, dx_rear)
    return ALUMINA_DENSITY_v31 * alumina_heat_capacity_v31(Ttube) *
           rear_tube_solid_area_v31(z_rear) * dx_rear
end

function rear_tube_htc_v31(Twall, Tgas, flow_lpm, Tin)
    flow_lpm <= 1e-12 && return 0.0
    mdot = m_dot_standard_v31(flow_lpm)
    Tfilm = 0.5 * (Twall + Tgas)
    cp = cpf_f(Tfilm)
    mu = mu_f_f(Tfilm)
    kf = kf_f(Tfilm)
    reynolds = mdot * REAR_TUBE_HYDRAULIC_DIAMETER_v31 / (REAR_TUBE_FLOW_AREA_v31 * mu)
    prandtl = cp * mu / kf
    nusselt = reynolds < 2300.0 ? 3.66 : 0.023 * reynolds^0.8 * prandtl^0.4
    return nusselt * kf / REAR_TUBE_HYDRAULIC_DIAMETER_v31
end

function gas_profile_v31!(Tg, Qgas, hcell, Qgas_rear, hrear,
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

    mdot_total = m_dot_standard_v31(flow)
    phi_act = receiver_active_flow_fraction_v31(flow, Tin, p)
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
        
        nusselt = receiver_channel_nusselt_v31(z[i], reynolds_active, prandtl, p)
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
        hrear[j] = rear_tube_htc_v31(Ttube[j], Tg[inlet_index], flow, Tin)
        UA = hrear[j] * REAR_TUBE_EXCHANGE_PERIMETER_v31 * dx_rear
        effectiveness = -expm1(-UA / (mdot_total * cp))
        Tg[inlet_index + 1] = Tg[inlet_index] + effectiveness * (Ttube[j] - Tg[inlet_index])
        Qgas_rear[j] = mdot_total * cp * (Tg[inlet_index + 1] - Tg[inlet_index])
    end
    return nothing
end

function receiver_rhs_v31!(du, u, time, p, operating, z, dx, core_weights,
                          perim_weights, rear_weights,
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
    rear_start = cavity_index + 1
    rear_stop = cavity_index + nodes
    Trear = view(u, rear_start:rear_stop)
    
    dTcore = view(du, 1:nodes)
    dTperim = view(du, (nodes + 1):(2 * nodes))
    dTtube = view(du, tube_start:tube_stop)
    dTrear = view(du, rear_start:rear_stop)
    ambient = operating.ambient_temperature(time)

    gas_profile_v31!(Tg, Qgas, hcell, Qgas_rear, hrear,
                    Tcore, Ttube, time, p, operating, z, dx, z_rear, dx_rear)
    fill!(du, 0.0)

    # Core axial solid conduction
    for i in 1:(nodes - 1)
        ki = ks_f(Tcore[i])
        kj = ks_f(Tcore[i + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        Qcond_core =
            core_axial_conduction_scale_v31(p) * kface * A_solid *
            (Tcore[i] - Tcore[i + 1]) / dx
        dTcore[i] -= Qcond_core
        dTcore[i + 1] += Qcond_core

        kwi = ks_f(Tperim[i])
        kwj = ks_f(Tperim[i + 1])
        kwface = 2.0 * kwi * kwj / (kwi + kwj)
        Tface_perim = 0.5 * (Tperim[i] + Tperim[i + 1])
        k_perim = perimeter_axial_conductivity_v31(Tface_perim, p)
        Qcond_perim = (kwface * A_solid * 0.20 + k_perim * A_frt) * (Tperim[i] - Tperim[i + 1]) / dx
        dTperim[i] -= Qcond_perim
        dTperim[i + 1] += Qcond_perim
    end

    irradiance = max(0.0, operating.irradiance(time))
    scale = absorbed_power_scale_v31(irradiance, p)
    Qreceiver_absorbed = ETA_ABS_FIXED_v31 * scale * irradiance * A_frt
    f_perim_source = perimeter_source_fraction_v31(p)
    Qcore_solar = Qreceiver_absorbed * (1.0 - f_perim_source)
    Qperim_solar = Qreceiver_absorbed * f_perim_source

    G_cp = core_perimeter_conductance_per_length_v31(p)
    radial_cavity_conductance = RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v31

    for i in 1:nodes
        Qrad_in = G_cp * dx * (Tperim[i] - Tcore[i])
        Qperim_cavity = radial_cavity_conductance * dx * (Tperim[i] - Tcavity)
        
        Qcell_core = Qcore_solar * core_weights[i]
        Qcell_perim = Qperim_solar * perim_weights[i]
        dTcore[i] += Qcell_core - Qgas[i] + Qrad_in
        dTperim[i] += Qcell_perim - Qrad_in - Qperim_cavity
        du[cavity_index] += Qperim_cavity
    end

    Qfront_perim = H_FRONT_FIXED_v31 * A_frt * (Tperim[1] - ambient) +
                   EPS_FRONT_FIXED_v31 * sigma_sb * A_frt * (Tperim[1]^4 - ambient^4)
    
    f_core_rear = rear_core_fraction_v31(p)
    G_receiver_rear = receiver_rear_conductance_v31(p)
    G_rear_cavity = rear_cavity_conductance_v31(p)

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

    G_rear_axial = rear_axial_conductance_v31(p)
    for i in 1:(nodes - 1)
        Qrear_axial = G_rear_axial * (Trear[i] - Trear[i + 1])
        dTrear[i] -= Qrear_axial
        dTrear[i + 1] += Qrear_axial
    end

    Qrear_tube = rear_tube_conductance_v31(p) * (Trear[end] - Ttube[1])
    dTrear[end] -= Qrear_tube
    dTtube[1] += Qrear_tube

    for j in 1:(rear_nodes - 1)
        ki = alumina_conductivity_v31(Ttube[j])
        kj = alumina_conductivity_v31(Ttube[j + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        Aface = 0.5 * (rear_tube_solid_area_v31(z_rear[j]) + rear_tube_solid_area_v31(z_rear[j + 1]))
        Qcond = kface * Aface * (Ttube[j] - Ttube[j + 1]) / dx_rear
        dTtube[j] -= Qcond
        dTtube[j + 1] += Qcond
    end

    for j in 1:rear_nodes
        G_cavity = rear_tube_cavity_conductance_per_length_v31(z_rear[j], p)
        G_flange = rear_tube_flange_conductance_per_length_v31(z_rear[j], Ttube[j], p, time, irradiance)
        Q_cavity = G_cavity * dx_rear * (Ttube[j] - Tcavity)
        Q_flange = G_flange * dx_rear * (Ttube[j] - WATER_FLANGE_TEMPERATURE_v31)
        dTtube[j] -= Qgas_rear[j] + Q_cavity + Q_flange
        du[cavity_index] += Q_cavity
    end

    Qcavity_ambient = cavity_loss_to_ambient_v31(Tcavity, ambient)
    du[cavity_index] -= Qcavity_ambient

    cell_volume = A_solid * dx
    perim_capacity_cell = perimeter_heat_capacity_total_v31(p) / nodes
    for i in 1:nodes
        core_capacity = rho_s * Cps_f(Tcore[i]) * cell_volume
        dTcore[i] /= max(core_capacity, eps(Float64))
        dTperim[i] /= max(perim_capacity_cell, eps(Float64))
    end
    for j in 1:rear_nodes
        dTtube[j] /= max(rear_tube_capacity_v31(Ttube[j], z_rear[j], dx_rear), eps(Float64))
    end
    du[cavity_index] /= max(CAVITY_HEAT_CAPACITY_v31, eps(Float64))
    C_rear_total = rear_reservoir_heat_capacity_v31(p)
    for i in 1:nodes
        dTrear[i] /= max(C_rear_total * rear_weights[i], eps(Float64))
    end
    return nothing
end

function receiver_ode_v31!(dTs, Ts, context, time)
    receiver_rhs_v31!(
        dTs, Ts, time, context.p, context.operating,
        context.z, context.dx, context.core_weights, context.perim_weights,
        context.rear_weights,
        context.z_rear, context.dx_rear,
        context.Tg, context.Qgas, context.hcell, context.Qgas_rear, context.hrear,
    )
    return nothing
end

function simulate_v31(p, operating, save_times;
                      initial_temperature=Tamb, nodes=default_nodes,
                      rear_nodes=REAR_TUBE_DEFAULT_NODES_v31,
                      initial_perim_temperature=nothing,
                      initial_rear_temperature=nothing,
                      initial_rear_reservoir_temperature=nothing,
                      initial_cavity_temperature=Tamb,
                      solver=Rodas5P(autodiff=AutoFiniteDiff()),
                      reltol=1e-6, abstol=1e-7, dtmax=30.0)
    length(p) == 25 || throw(ArgumentError("1D_v31 expects 25 parameters"))
    all(isfinite, p) || throw(ArgumentError("all model parameters must be finite"))
    all(lb_full_v31[i] <= p[i] <= ub_full_v31[i] for i in eachindex(p)) ||
        throw(ArgumentError("1D_v31 parameters are outside declared bounds"))

    time = Float64.(save_times)
    dx = L / nodes
    dx_rear = REAR_TUBE_LENGTH_v31 / rear_nodes
    z = collect(range(dx / 2, step=dx, length=nodes))
    z_rear = collect(range(dx_rear / 2, step=dx_rear, length=rear_nodes))
    zg = vcat(
        collect(range(0.0, L, length=nodes + 1)),
        L .+ collect(range(dx_rear, step=dx_rear, length=rear_nodes)),
    )
    core_weights = solar_weights_v5(p[13], p[6], nodes)
    perim_weights = axial_exponential_weights_v31(perimeter_source_attenuation_v31(p), nodes, dx)
    rear_weights = rear_contact_weights_v31(z)

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
        rear_core_fraction_v31(p) .* Tcore_initial .+
        (1.0 - rear_core_fraction_v31(p)) .* Tperim_initial
    elseif initial_rear_reservoir_temperature isa Function
        Float64.([initial_rear_reservoir_temperature(zi) for zi in z])
    else
        fill(Float64(initial_rear_reservoir_temperature), nodes)
    end
    u_initial = vcat(
        Tcore_initial, Tperim_initial, Ttube_initial,
        Float64(initial_cavity_temperature), Trear_initial,
    )

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
    rear_reservoir_temperature_history = Matrix{Float64}(undef, nodes, length(time))

    Tg = zeros(nodes + rear_nodes + 1)
    Qgas = zeros(nodes)
    hcell = zeros(nodes)
    Qgas_rear = zeros(rear_nodes)
    hrear = zeros(rear_nodes)

    context = (
        p=p, operating=operating, z=z, dx=dx,
        core_weights=core_weights, perim_weights=perim_weights,
        rear_weights=rear_weights,
        z_rear=z_rear, dx_rear=dx_rear,
        Tg=Tg, Qgas=Qgas, hcell=hcell, Qgas_rear=Qgas_rear, hrear=hrear,
    )
    problem = ODEProblem(receiver_ode_v31!, u_initial, (time[1], time[end]), context)
    solution = solve(
        problem, solver; saveat=time, save_everystep=false, dense=false,
        reltol=reltol, abstol=abstol, dtmax=dtmax, tstops=time,
    )
    successful_retcode(solution) || error("1D_v31 ODE solve failed with retcode $(solution.retcode)")

    state_history = reduce(hcat, solution.u)
    core_temperature_history .= state_history[1:nodes, :]
    perim_temperature_history .= state_history[(nodes + 1):(2 * nodes), :]
    tube_start = 2 * nodes + 1
    tube_stop = 2 * nodes + rear_nodes
    cavity_index = 2 * nodes + rear_nodes + 1
    rear_start = cavity_index + 1
    rear_stop = cavity_index + nodes
    rear_temperature_history .= state_history[tube_start:tube_stop, :]
    cavity_temperature_history .= state_history[cavity_index, :]
    rear_reservoir_temperature_history .= state_history[rear_start:rear_stop, :]

    nu_history = Matrix{Float64}(undef, nodes, length(time))
    effectiveness_history = Matrix{Float64}(undef, nodes, length(time))

    for output_index in eachindex(time)
        t = time[output_index]
        Tcore_now = view(core_temperature_history, :, output_index)
        Ttube_now = view(rear_temperature_history, :, output_index)
        gas_profile_v31!(Tg, Qgas, hcell, Qgas_rear, hrear, Tcore_now, Ttube_now,
                         t, p, operating, z, dx, z_rear, dx_rear)
        Tg_history[:, output_index] .= Tg
        h_history[:, output_index] .= hcell
        qgas_history[:, output_index] .= Qgas
        qgas_rear_history[:, output_index] .= Qgas_rear

        flow = max(0.0, operating.flow_lpm(t))
        Tin = operating.inlet_temperature(t)
        mdot_total = m_dot_standard_v31(flow)
        phi_act = receiver_active_flow_fraction_v31(flow, Tin, p)
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
        rear_reservoir_temperature=rear_reservoir_temperature_history,
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
        rear_contact_weights=rear_weights,
        ode_solution=solution,
    )
end

begin # v31 experimental interface
    const WALL_CHAIN_POSITIONS_v31 = (
        T8=sensor_positions[:T8],
        T12=sensor_positions[:T9],
        T11=sensor_positions[:T10],
    )

    function measured_initial_perim_profile_v31(data, simulation_id)
        z8 = WALL_CHAIN_POSITIONS_v31.T8
        z12 = WALL_CHAIN_POSITIONS_v31.T12
        z11 = WALL_CHAIN_POSITIONS_v31.T11
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

    function measured_initial_core_profile_v31(data, simulation_id)
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

    function measured_initial_rear_profile_v31(data, simulation_id, p)
        core_profile = measured_initial_core_profile_v31(data, simulation_id)
        perim_profile = measured_initial_perim_profile_v31(data, simulation_id)
        f_core_rear = rear_core_fraction_v31(p)
        T_receiver_exit =
            f_core_rear * core_profile(L) + (1.0 - f_core_rear) * perim_profile(L)
        T_tube_exit = observation(data, simulation_id, "_Tf")[1]
        return function (z_rear)
            fraction = clamp(z_rear / REAR_TUBE_LENGTH_v31, 0.0, 1.0)
            return T_receiver_exit + fraction * (T_tube_exit - T_receiver_exit)
        end
    end

    function measured_initial_rear_reservoir_temperature_v31(data, simulation_id, p)
        core_profile = measured_initial_core_profile_v31(data, simulation_id)
        perim_profile = measured_initial_perim_profile_v31(data, simulation_id)
        f_core_rear = rear_core_fraction_v31(p)
        return function (z)
            f_core_rear * core_profile(z) + (1.0 - f_core_rear) * perim_profile(z)
        end
    end

    function experimental_case_v31(simulation_id; is_cooling=false)
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

    function gas_at_position_v31(result, position)
        z = result.z_gas
        values = result.gas_temperature
        position <= z[1] && return vec(values[1, :])
        position >= z[end] && return vec(values[end, :])
        i = searchsortedlast(z, position)
        fraction = (position - z[i]) / (z[i + 1] - z[i])
        return vec((1.0 - fraction) .* values[i, :] .+ fraction .* values[i + 1, :])
    end

    function core_at_v31(result, position)
        z = result.z_solid
        values = result.core_temperature
        position <= z[1] && return vec(values[1, :])
        position >= z[end] && return vec(values[end, :])
        i = searchsortedlast(z, position)
        fraction = (position - z[i]) / (z[i + 1] - z[i])
        return vec((1.0 - fraction) .* values[i, :] .+ fraction .* values[i + 1, :])
    end

    function perim_at_v31(result, position)
        z = result.z_solid
        values = result.perim_temperature
        position <= z[1] && return vec(values[1, :])
        position >= z[end] && return vec(values[end, :])
        i = searchsortedlast(z, position)
        fraction = (position - z[i]) / (z[i + 1] - z[i])
        return vec((1.0 - fraction) .* values[i, :] .+ fraction .* values[i + 1, :])
    end

    function extract_outputs_v31(result, p)
        T8 = perim_at_v31(result, sensor_positions[:T8])
        T12_perim = perim_at_v31(result, WALL_CHAIN_POSITIONS_v31.T12)
        T11_perim = perim_at_v31(result, WALL_CHAIN_POSITIONS_v31.T11)
        T9_core = core_at_v31(result, sensor_positions[:T9])
        T10_core = core_at_v31(result, sensor_positions[:T10])
        Tf = gas_at_position_v31(result, T3_SAMPLE_POSITION_v31)
        T2 = result.cavity_temperature
        h_mean = vec(mean(result.heat_transfer_coefficient, dims=1))
        return hcat(T8, T12_perim, T11_perim, T9_core, T10_core, Tf, T2, h_mean)
    end

    function align_cooling_initial_outputs_v31!(outputs, experiment)
        size(outputs, 1) >= 1 || return outputs
        outputs[1, 1:7] .= experiment[1, 1:7]
        return outputs
    end

    function solve_case_v31(p, simulation_id; is_cooling=false,
                           nodes=default_nodes,
                           rear_nodes=REAR_TUBE_DEFAULT_NODES_v31,
                           solver=Rodas5P(autodiff=AutoFiniteDiff()),
                           reltol=1e-6, abstol=1e-7, dtmax=30.0)
        conditions, time, experiment, data = experimental_case_v31(
            simulation_id; is_cooling=is_cooling
        )
        flow = observation(data, simulation_id, "_flow")
        Tin = observation(data, simulation_id, "_Tin")
        ambient = observation(data, simulation_id, "_Tamb")
        if is_cooling
            Tin = fill(COOLING_ROOM_TEMPERATURE_v31, length(time))
            ambient = fill(COOLING_ROOM_TEMPERATURE_v31, length(time))
        end
        T2 = observation(data, simulation_id, "_T2")
        irradiance = is_cooling ? zeros(length(time)) : fill(conditions[Io], length(time))
        operating = operating_condition_v3(
            irradiance=linear_history_v3(time, irradiance),
            flow_lpm=linear_history_v3(time, flow),
            inlet_temperature=linear_history_v3(time, Tin),
            ambient_temperature=linear_history_v3(time, ambient),
        )
        result = simulate_v31(
            p, operating, time;
            initial_temperature=measured_initial_core_profile_v31(data, simulation_id),
            initial_perim_temperature=measured_initial_perim_profile_v31(data, simulation_id),
            initial_rear_temperature=measured_initial_rear_profile_v31(data, simulation_id, p),
            initial_rear_reservoir_temperature=measured_initial_rear_reservoir_temperature_v31(data, simulation_id, p),
            initial_cavity_temperature=T2[1],
            nodes=nodes, rear_nodes=rear_nodes, solver=solver,
            reltol=reltol, abstol=abstol, dtmax=dtmax,
        )
        outputs = extract_outputs_v31(result, p)
        is_cooling && align_cooling_initial_outputs_v31!(outputs, experiment)
        return outputs, result, experiment
    end
end

begin # v31 objective helpers
    const ORDERING_SCALE_v31 = 100.0
    const ORDERING_WEIGHT_v31 = 1.00
    const TIMING_WEIGHT_v31 = 0.15
    const SLOPE_WEIGHT_v31 = 0.25
    const COOLING_UPTURN_WEIGHT_v31 = 10000.0

    function normalized_slope_mse_v31(model, experiment)
        length(model) < 2 && return 0.0
        scale = max(maximum(experiment) - minimum(experiment), 1.0)
        return mean(abs2, diff(model) .- diff(experiment)) / scale^2
    end

    function timing_penalty_v31(time, model, experiment)
        duration = max(time[end] - time[1], 1.0)
        return ((get_t90_v3(time, model) - get_t90_v3(time, experiment)) / duration)^2
    end

    function signal_loss_v31(time, model, experiment)
        return normalized_mse_v3(model, experiment) +
               SLOPE_WEIGHT_v31 * normalized_slope_mse_v31(model, experiment) +
               TIMING_WEIGHT_v31 * timing_penalty_v31(time, model, experiment)
    end

    function cooling_upturn_loss_v31(model)
        target_columns = (3, 5, 6) # T11, T10, T3
        total = 0.0
        for j in target_columns
            scale = max(maximum(model[:, j]) - minimum(model[:, j]), 20.0)
            upturn = max.(diff(model[:, j]), 0.0) ./ scale
            total += maximum(upturn)^2 + mean(abs2, upturn)
        end
        return total / length(target_columns)
    end

    function ordering_loss_v31(model, experiment)
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
        return mean(abs2, collect(model_offsets) .- collect(experiment_offsets)) / ORDERING_SCALE_v31^2
    end

    function capacitance_regularization_v31(p; nodes=15)
        capacitance = participating_total_heat_capacity_v31(p, nodes)
        return CAPACITANCE_REGULARIZATION_WEIGHT_v31 *
               ((capacitance - MEASURED_ASSEMBLY_CAPACITANCE_v31) / MEASURED_ASSEMBLY_CAPACITANCE_v31)^2
    end

    function loss_cases_v31(p, keys; is_cooling=false, nodes=15)
        total = 0.0
        for simulation_id in keys
            try
                model, result, experiment = solve_case_v31(
                    p, simulation_id; is_cooling=is_cooling, nodes=nodes,
                    reltol=2e-5, abstol=2e-6, dtmax=60.0,
                )
                all(isfinite, model) || return Inf
                total += mean(
                    signal_loss_v31(result.time, model[:, j], experiment[:, j])
                    for j in 1:7
                )
                if !is_cooling
                    total += ORDERING_WEIGHT_v31 * ordering_loss_v31(model, experiment)
                else
                    total += COOLING_UPTURN_WEIGHT_v31 * cooling_upturn_loss_v31(model)
                end
            catch err
                err isa InterruptException && rethrow()
                return Inf
            end
        end
        return total / length(keys)
    end

    function loss_heating_v31(p, keys=sim_key_heat; nodes=15)
        regularization = capacitance_regularization_v31(p; nodes=nodes)
        return loss_cases_v31(p, keys; is_cooling=false, nodes=nodes) + regularization
    end

    function loss_cooling_v31(p, keys=sim_key_cool; nodes=15)
        regularization = capacitance_regularization_v31(p; nodes=nodes)
        return loss_cases_v31(p, keys; is_cooling=true, nodes=nodes) + regularization
    end

    function fit_indices_for_stage_v31(stage)
        stage_symbol = Symbol(stage)
        stage_symbol == :rear && return fit_rear_stage_indices_v31
        stage_symbol == :transport && return fit_transport_stage_indices_v31
        stage_symbol == :full && return fit_full_stage_indices_v31
        throw(ArgumentError("unknown v31 calibration stage: $stage"))
    end

    function calibrate_v31(; nodes=15, maximum_iterations=100,
                           maximum_time_seconds=360.0, stage=:rear)
        p0 = copy(pnew_v31)
        idx_fit = fit_indices_for_stage_v31(stage)
        function obj_free(theta)
            p_full = copy(p0)
            p_full[idx_fit] .= theta
            return loss_heating_v31(p_full; nodes=nodes) + loss_cooling_v31(p_full; nodes=nodes)
        end
        res = optimize_with_nlopt_v3(
            obj_free,
            pnew_v31[idx_fit],
            lb_full_v31[idx_fit],
            ub_full_v31[idx_fit];
            maximum_iterations=maximum_iterations,
            maximum_time_seconds=maximum_time_seconds,
            algorithm=NLopt.LN_BOBYQA(),
            label="fit-v31-$(Symbol(stage))",
        )
        p_opt = copy(p0)
        p_opt[idx_fit] .= res.minimizer
        return (objective=res.minimum, parameters=p_opt, minimizer=res.minimizer, retcode=res.retcode)
    end
end








