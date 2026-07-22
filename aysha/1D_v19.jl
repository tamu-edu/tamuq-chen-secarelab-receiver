# ============================================================================
# 1D_v19.jl - 2-Zone Core/Perimeter Macro-ECM Receiver Model
# ============================================================================
# v19 represents the architectural transition from 1D single-channel fitting
# to a 2-Zone Core/Perimeter Macro-ECM formulation.
#
# Key Physical Features of v19:
#   1. Core Channel Zone (T_core): 100-channel matrix. Receives direct central
#      solar beam. Exchanges heat with gas using bounded laminar developing-flow
#      convection (Nu floor = 3.61). Compares to interior probes T9 and T10.
#   2. Perimeter Housing Zone (T_perim): Side walls, alumina felt insulation,
#      and front housing. Receives front-rim spillage power from the v14 power
#      scales (P_spill = (scale - 1.0) * Q_solar). Holds the ~301 J/K housing
#      thermal flywheel. Conducts heat radially into the core zone. Compares to
#      side-wall probes T8, T12, T11 and insulation probe T2.
#   3. Active Flow Participation Fraction:
#      phi_active(Re) = clamp(phi_0 * (Re / Re_ref)^m_rec, phi_min, 1.0)
#      Gas exchanges heat only with the active fraction phi_active(Re) of core flow.
#
# Parameter Vector p[1:13] in v19:
#   p[1]  A_Nu          laminar developing Nu prefactor
#   p[2]  B_Re          laminar developing Reynolds exponent (bounded 0.0 to 0.6)
#   p[3]  C_Pr          fixed Prandtl exponent (= 1/3)
#   p[4]  phi_0         active flow fraction at Re = 50
#   p[5]  m_rec         flow recruitment exponent
#   p[6]  front_dep     frozen v12a front-cell deposition fraction (= 1.0 or 0.0)
#   p[7]  scale_456     fixed v14 power scale at 456 kW/m2 (= 1.6695)
#   p[8]  scale_304     fixed v14 power scale at 304 kW/m2 (= 1.7171)
#   p[9]  scale_256     fixed v14 power scale at 256 kW/m2 (= 0.9827)
#   p[10] G_core_perim  radial core-perimeter conductance per length [W/m/K]
#   p[11] C_perim_eff   perimeter/housing participating heat capacity [J/K]
#   p[12] k_perim_ref   perimeter axial conductivity at 900 K [W/m/K]
#   p[13] beta_opt      frozen v12a optical attenuation [1/m]
# ============================================================================

include("1D_v5.jl")

begin # v19 fixed constants
    const EPS_FRONT_FIXED_v19 = EPS_FRONT_FIXED_V5
    const ETA_ABS_FIXED_v19 = ETA_ABS_FIXED_V5
    const BETA_OPT_FIXED_v19 = BETA_OPT_FIXED_V5
    const H_FRONT_FIXED_v19 = 10.0
    const PRANDTL_EXPONENT_FIXED_v19 = 1.0 / 3.0
    const STANDARD_AIR_DENSITY_v19 = 101325.0 / (287.05 * 294.25)
    const REYNOLDS_REFERENCE_v19 = 50.0
    const NU_FLOOR_v19 = 3.61
    const MEASURED_ASSEMBLY_CAPACITANCE_v19 = 301.0
    const IRRADIANCE_LEVELS_v19 = (456000.0, 304000.0, 256000.0)
    const T3_SAMPLE_POSITION_v19 = 140.0e-3
    const CAPACITANCE_REGULARIZATION_WEIGHT_v19 = 0.05

    # Cavity and housing geometry
    const CAVITY_OUTER_RADIUS_v19 = 75.0e-3
    const METAL_THICKNESS_v19 = 18.0e-3
    const INSULATION_OUTER_RADIUS_v19 = CAVITY_OUTER_RADIUS_v19 - METAL_THICKNESS_v19
    const METAL_OUTER_RADIUS_v19 = CAVITY_OUTER_RADIUS_v19
    const CAVITY_LENGTH_v19 = 165.0e-3

    const ADAPTOR_DIAMETER_v19 = 77.6e-3
    const ADAPTOR_RADIUS_v19 = ADAPTOR_DIAMETER_v19 / 2.0
    const ADAPTOR_LENGTH_v19 = 57.0e-3
    const ADAPTOR_TUBE_RADIUS_v19 = 6.5e-3
    const ADAPTOR_RECEIVER_OVERLAP_LENGTH_v19 = 28.0e-3
    const ADAPTOR_TUBE_OVERLAP_LENGTH_v19 = ADAPTOR_LENGTH_v19 - ADAPTOR_RECEIVER_OVERLAP_LENGTH_v19
    const ADAPTOR_OVERLAP_LENGTH_v19 = ADAPTOR_RECEIVER_OVERLAP_LENGTH_v19

    const RECEIVER_EQ_RADIUS_v19 = sqrt(A_frt / pi)
    const RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v19 =
        4.0 * pi * 0.08 / log(INSULATION_OUTER_RADIUS_v19 / RECEIVER_EQ_RADIUS_v19)
    const ADAPTOR_CONTACT_RESISTANCE_v19 =
        1.0 / (100.0 * 4.0 * w_t * ADAPTOR_OVERLAP_LENGTH_v19)
    const RECEIVER_TO_TUBE_CONDUCTANCE_v19 = 1.0 / ADAPTOR_CONTACT_RESISTANCE_v19

    # Rear tube geometry
    const REAR_TUBE_LENGTH_v19 = 150.0e-3
    const REAR_TUBE_CAVITY_LENGTH_v19 = 46.0e-3
    const REAR_TUBE_DEFAULT_NODES_v19 = 12
    const REAR_TUBE_GAS_RADIUS_v19 = ADAPTOR_TUBE_RADIUS_v19
    const REAR_TUBE_WALL_THICKNESS_v19 = 1.5e-3
    const REAR_TUBE_OUTER_RADIUS_v19 = REAR_TUBE_GAS_RADIUS_v19 + REAR_TUBE_WALL_THICKNESS_v19
    const REAR_TUBE_FLOW_AREA_v19 = pi * REAR_TUBE_GAS_RADIUS_v19^2
    const REAR_TUBE_HYDRAULIC_DIAMETER_v19 = 2.0 * REAR_TUBE_GAS_RADIUS_v19
    const REAR_TUBE_EXCHANGE_PERIMETER_v19 = 2.0 * pi * REAR_TUBE_GAS_RADIUS_v19

    const ALUMINA_DENSITY_v19 = 3900.0
    const ALUMINUM_DENSITY_v19 = 2700.0
    const ALUMINUM_HEAT_CAPACITY_v19 = 900.0
    const ALUMINUM_CONDUCTIVITY_v19 = 205.0
    const METAL_EMISSIVITY_v19 = 0.2
    const H_NAT_CAVITY_v19 = 10.0
    const INSULATION_DENSITY_v19 = 140.0
    const INSULATION_HEAT_CAPACITY_v19 = 1360.0
    const WATER_FLANGE_TEMPERATURE_v19 = 25.0 + 273.15

    const CAVITY_FELT_VOLUME_v19 = max(
        0.0,
        pi * INSULATION_OUTER_RADIUS_v19^2 * CAVITY_LENGTH_v19 -
        A_frt * L - pi * ADAPTOR_RADIUS_v19^2 * ADAPTOR_LENGTH_v19,
    )
    const CAVITY_METAL_VOLUME_v19 =
        pi * (METAL_OUTER_RADIUS_v19^2 - INSULATION_OUTER_RADIUS_v19^2) * CAVITY_LENGTH_v19 +
        pi * METAL_OUTER_RADIUS_v19^2 * METAL_THICKNESS_v19
    const CAVITY_OUTER_AREA_v19 =
        2.0 * pi * METAL_OUTER_RADIUS_v19 * CAVITY_LENGTH_v19 + pi * METAL_OUTER_RADIUS_v19^2
    const CAVITY_HEAT_CAPACITY_v19 =
        INSULATION_DENSITY_v19 * INSULATION_HEAT_CAPACITY_v19 * CAVITY_FELT_VOLUME_v19 +
        ALUMINUM_DENSITY_v19 * ALUMINUM_HEAT_CAPACITY_v19 * CAVITY_METAL_VOLUME_v19
end

# Default v19 parameter vector p[1:13]
begin # v19 parameter values and bounds
    pnew_v19 = [
        1.0996,             # p[1]  A_Nu (laminar developing prefactor)
        0.3333,             # p[2]  B_Re (Reynolds exponent bounded 0.0 to 0.6)
        PRANDTL_EXPONENT_FIXED_v19, # p[3]  C_Pr (= 1/3)
        0.65,               # p[4]  phi_0 (active flow fraction at Re = 50)
        0.50,               # p[5]  m_rec (flow recruitment exponent)
        1.0,                # p[6]  front_dep (frozen v12a front-cell deposition)
        1.669539299381132,  # p[7]  scale_456
        1.7170694424514947, # p[8]  scale_304
        0.9827026344016692, # p[9]  scale_256
        15.0,               # p[10] G_core_perim (radial core-perimeter conductance [W/m/K])
        250.0,              # p[11] C_perim_eff (perimeter participating capacity [J/K])
        2.0,                # p[12] k_perim_ref (perimeter axial conductivity at 900K [W/m/K])
        184.67243519237965, # p[13] beta_opt [1/m]
    ]
    lb_full_v19 = [
        0.01,
        0.0,
        pnew_v19[3],
        0.20,
        0.0,
        pnew_v19[6],
        0.50,
        0.50,
        0.50,
        0.5,
        50.0,
        0.0,
        pnew_v19[13],
    ]
    ub_full_v19 = [
        5.00,
        0.60,
        pnew_v19[3],
        1.00,
        1.50,
        pnew_v19[6],
        2.50,
        2.50,
        2.50,
        100.0,
        450.0,
        50.0,
        pnew_v19[13],
    ]

    fit_heat_transfer_indices_v19 = [1, 2, 4, 5, 7, 8, 9, 10, 11, 12]
    fit_source_indices_v19 = Int[]
    fit_power_scale_indices_v19 = [7, 8, 9]
    fit_radiation_indices_v19 = Int[]
end

function irradiance_level_index_v19(irradiance)
    irradiance <= 0.0 && return 0
    differences = abs.(collect(IRRADIANCE_LEVELS_v19) .- Float64(irradiance))
    return argmin(differences)
end

power_scale_index_v19(irradiance) = begin
    level_index = irradiance_level_index_v19(irradiance)
    level_index == 0 ? 0 : 6 + level_index
end

function absorbed_power_scale_v19(irradiance, p)
    index = power_scale_index_v19(irradiance)
    return index == 0 ? 1.0 : p[index]
end

core_perimeter_conductance_per_length_v19(p) = max(p[10], eps(Float64))

perimeter_heat_capacity_total_v19(p) = max(p[11], eps(Float64))

core_receiver_heat_capacity_v19(nodes) =
    sum(rho_s * Cps_f(900.0) * A_solid * (L / nodes) for _ in 1:nodes)

participating_assembly_heat_capacity_v19(p, nodes=default_nodes) =
    core_receiver_heat_capacity_v19(nodes) + perimeter_heat_capacity_total_v19(p)

perimeter_axial_conductivity_v19(T, p) =
    max(p[12], 0.0) * (property_temperature(T) / 900.0)^3

function active_flow_fraction_v19(reynolds, p)
    phi0 = clamp(p[4], 0.1, 1.0)
    m = max(p[5], 0.0)
    ratio = max(reynolds, eps(Float64)) / REYNOLDS_REFERENCE_v19
    return clamp(phi0 * ratio^m, 0.1, 1.0)
end

function receiver_channel_nusselt_v19(z, reynolds, prandtl, p)
    A_nu = max(p[1], 0.0)
    B_re = clamp(p[2], 0.0, 0.6)
    C_pr = p[3]
    entry_factor = (Dh / max(z, Dh / 2.0))^(1.0 / 3.0)
    Nu_dev = A_nu * max(reynolds, eps(Float64))^B_re * max(prandtl, eps(Float64))^C_pr * entry_factor
    return max(NU_FLOOR_v19, Nu_dev)
end

m_dot_standard_v19(flow_lpm) =
    STANDARD_AIR_DENSITY_v19 * max(0.0, flow_lpm) / 60000.0

alumina_conductivity_v19(T) =
    5.5 + 34.5 * exp(-0.0033 * (property_temperature(T) - 273.15))

alumina_heat_capacity_v19(T) =
    max(500.0, (1.00446 + 1.742e-4 * property_temperature(T) -
                2.796e4 * property_temperature(T)^-2) * 1000.0)

function rear_tube_outer_radius_v19(z_rear)
    return z_rear <= ADAPTOR_TUBE_OVERLAP_LENGTH_v19 ?
           ADAPTOR_RADIUS_v19 : REAR_TUBE_OUTER_RADIUS_v19
end

function rear_tube_solid_area_v19(z_rear)
    outer_radius = rear_tube_outer_radius_v19(z_rear)
    return pi * (outer_radius^2 - REAR_TUBE_GAS_RADIUS_v19^2)
end

function rear_tube_cavity_conductance_per_length_v19(z_rear, p)
    z_rear > REAR_TUBE_CAVITY_LENGTH_v19 && return 0.0
    outer_radius = min(rear_tube_outer_radius_v19(z_rear), 0.98 * INSULATION_OUTER_RADIUS_v19)
    return 2.0 * pi * 0.08 / log(INSULATION_OUTER_RADIUS_v19 / outer_radius)
end

function rear_tube_flange_conductance_per_length_v19(z_rear, Ttube)
    z_rear <= REAR_TUBE_CAVITY_LENGTH_v19 && return 0.0
    tube_wall_resistance =
        log(REAR_TUBE_OUTER_RADIUS_v19 / REAR_TUBE_GAS_RADIUS_v19) /
        (2.0 * pi * alumina_conductivity_v19(Ttube))
    aluminum_resistance =
        log(METAL_OUTER_RADIUS_v19 / REAR_TUBE_OUTER_RADIUS_v19) /
        (2.0 * pi * ALUMINUM_CONDUCTIVITY_v19)
    return 1.0 / (tube_wall_resistance + aluminum_resistance)
end

function cavity_loss_to_ambient_v19(Tcavity, ambient)
    return H_NAT_CAVITY_v19 * CAVITY_OUTER_AREA_v19 * (Tcavity - ambient) +
           METAL_EMISSIVITY_v19 * sigma_sb * CAVITY_OUTER_AREA_v19 *
           (Tcavity^4 - ambient^4)
end

function rear_tube_capacity_v19(Ttube, z_rear, dx_rear)
    return ALUMINA_DENSITY_v19 * alumina_heat_capacity_v19(Ttube) *
           rear_tube_solid_area_v19(z_rear) * dx_rear
end

function rear_tube_htc_v19(Twall, Tgas, flow_lpm, Tin)
    flow_lpm <= 1e-12 && return 0.0
    mdot = m_dot_standard_v19(flow_lpm)
    Tfilm = 0.5 * (Twall + Tgas)
    cp = cpf_f(Tfilm)
    mu = mu_f_f(Tfilm)
    kf = kf_f(Tfilm)
    reynolds = mdot * REAR_TUBE_HYDRAULIC_DIAMETER_v19 / (REAR_TUBE_FLOW_AREA_v19 * mu)
    prandtl = cp * mu / kf
    nusselt = reynolds < 2300.0 ? 3.66 : 0.023 * reynolds^0.8 * prandtl^0.4
    return nusselt * kf / REAR_TUBE_HYDRAULIC_DIAMETER_v19
end

function gas_profile_v19!(Tg, Qgas, hcell, Qgas_rear, hrear,
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

    mdot_total = m_dot_standard_v19(flow)
    for i in eachindex(Tcore)
        Tfilm = 0.5 * (Tcore[i] + Tg[i])
        cp = cpf_f(Tfilm)
        mu = mu_f_f(Tfilm)
        kf = kf_f(Tfilm)
        reynolds_total = mdot_total * Dh / (A_flow * mu)
        prandtl = cp * mu / kf
        
        phi_act = active_flow_fraction_v19(reynolds_total, p)
        mdot_active = phi_act * mdot_total
        reynolds_active = reynolds_total * phi_act
        
        nusselt = receiver_channel_nusselt_v19(z[i], reynolds_active, prandtl, p)
        hcell[i] = nusselt * kf / Dh
        UA = hcell[i] * P_exchange * dx
        effectiveness = -expm1(-UA / (mdot_active * cp))
        
        Tg[i + 1] = Tg[i] + effectiveness * (Tcore[i] - Tg[i])
        Qgas[i] = mdot_active * cp * (Tg[i + 1] - Tg[i])
    end

    for j in eachindex(Ttube)
        inlet_index = receiver_nodes + j
        Tfilm = 0.5 * (Ttube[j] + Tg[inlet_index])
        cp = cpf_f(Tfilm)
        hrear[j] = rear_tube_htc_v19(Ttube[j], Tg[inlet_index], flow, Tin)
        UA = hrear[j] * REAR_TUBE_EXCHANGE_PERIMETER_v19 * dx_rear
        effectiveness = -expm1(-UA / (mdot_total * cp))
        Tg[inlet_index + 1] = Tg[inlet_index] + effectiveness * (Ttube[j] - Tg[inlet_index])
        Qgas_rear[j] = mdot_total * cp * (Tg[inlet_index + 1] - Tg[inlet_index])
    end
    return nothing
end

function receiver_rhs_v19!(du, u, time, p, operating, z, dx, solar_weights,
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

    gas_profile_v19!(Tg, Qgas, hcell, Qgas_rear, hrear,
                    Tcore, Ttube, time, p, operating, z, dx, z_rear, dx_rear)
    fill!(du, 0.0)

    # Core axial solid conduction
    for i in 1:(nodes - 1)
        ki = ks_f(Tcore[i])
        kj = ks_f(Tcore[i + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        Qcond_core = kface * A_solid * (Tcore[i] - Tcore[i + 1]) / dx
        dTcore[i] -= Qcond_core
        dTcore[i + 1] += Qcond_core

        kwi = ks_f(Tperim[i])
        kwj = ks_f(Tperim[i + 1])
        kwface = 2.0 * kwi * kwj / (kwi + kwj)
        Tface_perim = 0.5 * (Tperim[i] + Tperim[i + 1])
        k_perim = perimeter_axial_conductivity_v19(Tface_perim, p)
        Qcond_perim = (kwface * A_solid * 0.20 + k_perim * A_frt) * (Tperim[i] - Tperim[i + 1]) / dx
        dTperim[i] -= Qcond_perim
        dTperim[i + 1] += Qcond_perim
    end

    irradiance = max(0.0, operating.irradiance(time))
    scale = absorbed_power_scale_v19(irradiance, p)
    Qsolar_total = ETA_ABS_FIXED_v19 * scale * irradiance * A_frt

    # Solar power partition: nominal 1.0 scale goes to Core, overspill (scale - 1.0) goes to Perimeter
    if scale > 1.0
        Qcore_solar = ETA_ABS_FIXED_v19 * 1.0 * irradiance * A_frt
        Qperim_solar = ETA_ABS_FIXED_v19 * (scale - 1.0) * irradiance * A_frt
    else
        Qcore_solar = Qsolar_total
        Qperim_solar = 0.0
    end

    G_cp = core_perimeter_conductance_per_length_v19(p)
    radial_cavity_conductance = RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v19

    for i in 1:nodes
        Qrad_in = G_cp * dx * (Tperim[i] - Tcore[i])
        Qperim_cavity = radial_cavity_conductance * dx * (Tperim[i] - Tcavity)
        
        Qcell_core = Qcore_solar * solar_weights[i]
        Qcell_perim = Qperim_solar * solar_weights[i]

        dTcore[i] += Qcell_core - Qgas[i] + Qrad_in
        dTperim[i] += Qcell_perim - Qrad_in - Qperim_cavity
        du[cavity_index] += Qperim_cavity
    end

    Qfront_perim = H_FRONT_FIXED_v19 * A_frt * (Tperim[1] - ambient) +
                   EPS_FRONT_FIXED_v19 * sigma_sb * A_frt * (Tperim[1]^4 - ambient^4)
    
    Qreceiver_tube_core = RECEIVER_TO_TUBE_CONDUCTANCE_v19 * 0.70 * (Tcore[end] - Ttube[1])
    Qreceiver_tube_perim = RECEIVER_TO_TUBE_CONDUCTANCE_v19 * 0.30 * (Tperim[end] - Ttube[1])
    Qreceiver_tube = Qreceiver_tube_core + Qreceiver_tube_perim

    dTperim[1] -= Qfront_perim
    dTcore[end] -= Qreceiver_tube_core
    dTperim[end] -= Qreceiver_tube_perim
    dTtube[1] += Qreceiver_tube

    for j in 1:(rear_nodes - 1)
        ki = alumina_conductivity_v19(Ttube[j])
        kj = alumina_conductivity_v19(Ttube[j + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        Aface = 0.5 * (rear_tube_solid_area_v19(z_rear[j]) + rear_tube_solid_area_v19(z_rear[j + 1]))
        Qcond = kface * Aface * (Ttube[j] - Ttube[j + 1]) / dx_rear
        dTtube[j] -= Qcond
        dTtube[j + 1] += Qcond
    end

    for j in 1:rear_nodes
        G_cavity = rear_tube_cavity_conductance_per_length_v19(z_rear[j], p)
        G_flange = rear_tube_flange_conductance_per_length_v19(z_rear[j], Ttube[j])
        Q_cavity = G_cavity * dx_rear * (Ttube[j] - Tcavity)
        Q_flange = G_flange * dx_rear * (Ttube[j] - WATER_FLANGE_TEMPERATURE_v19)
        dTtube[j] -= Qgas_rear[j] + Q_cavity + Q_flange
        du[cavity_index] += Q_cavity
    end

    Qcavity_ambient = cavity_loss_to_ambient_v19(Tcavity, ambient)
    du[cavity_index] -= Qcavity_ambient

    cell_volume = A_solid * dx
    perim_capacity_cell = perimeter_heat_capacity_total_v19(p) / nodes
    for i in 1:nodes
        core_capacity = rho_s * Cps_f(Tcore[i]) * cell_volume
        dTcore[i] /= max(core_capacity, eps(Float64))
        dTperim[i] /= max(perim_capacity_cell, eps(Float64))
    end
    for j in 1:rear_nodes
        dTtube[j] /= max(rear_tube_capacity_v19(Ttube[j], z_rear[j], dx_rear), eps(Float64))
    end
    du[cavity_index] /= max(CAVITY_HEAT_CAPACITY_v19, eps(Float64))
    return nothing
end

function receiver_ode_v19!(dTs, Ts, context, time)
    receiver_rhs_v19!(
        dTs, Ts, time, context.p, context.operating,
        context.z, context.dx, context.solar_weights,
        context.z_rear, context.dx_rear,
        context.Tg, context.Qgas, context.hcell, context.Qgas_rear, context.hrear,
    )
    return nothing
end

function simulate_v19(p, operating, save_times;
                      initial_temperature=Tamb, nodes=default_nodes,
                      rear_nodes=REAR_TUBE_DEFAULT_NODES_v19,
                      initial_perim_temperature=nothing,
                      initial_rear_temperature=nothing,
                      initial_cavity_temperature=Tamb,
                      solver=Rodas5P(autodiff=AutoFiniteDiff()),
                      reltol=1e-6, abstol=1e-7, dtmax=30.0)
    length(p) == 13 || throw(ArgumentError("1D_v19 expects 13 parameters"))
    all(isfinite, p) || throw(ArgumentError("all model parameters must be finite"))

    time = Float64.(save_times)
    dx = L / nodes
    dx_rear = REAR_TUBE_LENGTH_v19 / rear_nodes
    z = collect(range(dx / 2, step=dx, length=nodes))
    z_rear = collect(range(dx_rear / 2, step=dx_rear, length=rear_nodes))
    zg = vcat(
        collect(range(0.0, L, length=nodes + 1)),
        L .+ collect(range(dx_rear, step=dx_rear, length=rear_nodes)),
    )
    weights = solar_weights_v5(p[13], p[6], nodes)

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
        p=p, operating=operating, z=z, dx=dx, solar_weights=weights,
        z_rear=z_rear, dx_rear=dx_rear,
        Tg=Tg, Qgas=Qgas, hcell=hcell, Qgas_rear=Qgas_rear, hrear=hrear,
    )
    problem = ODEProblem(receiver_ode_v19!, u_initial, (time[1], time[end]), context)
    solution = solve(
        problem, solver; saveat=time, save_everystep=false, dense=false,
        reltol=reltol, abstol=abstol, dtmax=dtmax, tstops=time,
    )
    successful_retcode(solution) || error("1D_v19 ODE solve failed with retcode $(solution.retcode)")

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
        gas_profile_v19!(Tg, Qgas, hcell, Qgas_rear, hrear, Tcore_now, Ttube_now,
                         t, p, operating, z, dx, z_rear, dx_rear)
        Tg_history[:, output_index] .= Tg
        h_history[:, output_index] .= hcell
        qgas_history[:, output_index] .= Qgas

        flow = max(0.0, operating.flow_lpm(t))
        mdot_total = m_dot_standard_v19(flow)
        for i in 1:nodes
            Tfilm = 0.5 * (Tcore_now[i] + Tg_history[i, output_index])
            cp = cpf_f(Tfilm)
            kf = kf_f(Tfilm)
            nu_history[i, output_index] = h_history[i, output_index] * Dh / max(kf, eps(Float64))
            if flow > 1e-12
                phi_act = active_flow_fraction_v19(mdot_total * Dh / (A_flow * mu_f_f(Tfilm)), p)
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
        ode_solution=solution,
    )
end

begin # v19 experimental interface
    const WALL_CHAIN_POSITIONS_v19 = (
        T8=sensor_positions[:T8],
        T12=sensor_positions[:T9],
        T11=sensor_positions[:T10],
    )

    function measured_initial_perim_profile_v19(data, simulation_id)
        z8 = WALL_CHAIN_POSITIONS_v19.T8
        z12 = WALL_CHAIN_POSITIONS_v19.T12
        z11 = WALL_CHAIN_POSITIONS_v19.T11
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

    function measured_initial_core_profile_v19(data, simulation_id)
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

    function measured_initial_rear_profile_v19(data, simulation_id)
        perim_profile = measured_initial_perim_profile_v19(data, simulation_id)
        T_receiver_exit = perim_profile(L)
        T_tube_exit = observation(data, simulation_id, "_Tf")[1]
        return function (z_rear)
            fraction = clamp(z_rear / REAR_TUBE_LENGTH_v19, 0.0, 1.0)
            return T_receiver_exit + fraction * (T_tube_exit - T_receiver_exit)
        end
    end

    function experimental_case_v19(simulation_id; is_cooling=false)
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

    function gas_at_position_v19(result, position)
        z = result.z_gas
        values = result.gas_temperature
        position <= z[1] && return vec(values[1, :])
        position >= z[end] && return vec(values[end, :])
        i = searchsortedlast(z, position)
        fraction = (position - z[i]) / (z[i + 1] - z[i])
        return vec((1.0 - fraction) .* values[i, :] .+ fraction .* values[i + 1, :])
    end

    function core_at_v19(result, position)
        z = result.z_solid
        values = result.core_temperature
        position <= z[1] && return vec(values[1, :])
        position >= z[end] && return vec(values[end, :])
        i = searchsortedlast(z, position)
        fraction = (position - z[i]) / (z[i + 1] - z[i])
        return vec((1.0 - fraction) .* values[i, :] .+ fraction .* values[i + 1, :])
    end

    function perim_at_v19(result, position)
        z = result.z_solid
        values = result.perim_temperature
        position <= z[1] && return vec(values[1, :])
        position >= z[end] && return vec(values[end, :])
        i = searchsortedlast(z, position)
        fraction = (position - z[i]) / (z[i + 1] - z[i])
        return vec((1.0 - fraction) .* values[i, :] .+ fraction .* values[i + 1, :])
    end

    function extract_outputs_v19(result, p)
        T8 = perim_at_v19(result, sensor_positions[:T8])
        T12_perim = perim_at_v19(result, WALL_CHAIN_POSITIONS_v19.T12)
        T11_perim = perim_at_v19(result, WALL_CHAIN_POSITIONS_v19.T11)
        T9_core = core_at_v19(result, sensor_positions[:T9])
        T10_core = core_at_v19(result, sensor_positions[:T10])
        Tf = gas_at_position_v19(result, T3_SAMPLE_POSITION_v19)
        T2 = result.cavity_temperature
        h_mean = vec(mean(result.heat_transfer_coefficient, dims=1))
        return hcat(T8, T12_perim, T11_perim, T9_core, T10_core, Tf, T2, h_mean)
    end

    function solve_case_v19(p, simulation_id; is_cooling=false,
                           nodes=default_nodes,
                           rear_nodes=REAR_TUBE_DEFAULT_NODES_v19,
                           solver=Rodas5P(autodiff=AutoFiniteDiff()),
                           reltol=1e-6, abstol=1e-7, dtmax=30.0)
        conditions, time, experiment, data = experimental_case_v19(
            simulation_id; is_cooling=is_cooling
        )
        flow = observation(data, simulation_id, "_flow")
        Tin = observation(data, simulation_id, "_Tin")
        ambient = observation(data, simulation_id, "_Tamb")
        T2 = observation(data, simulation_id, "_T2")
        irradiance = fill(conditions[Io], length(time))
        operating = operating_condition_v3(
            irradiance=linear_history_v3(time, irradiance),
            flow_lpm=linear_history_v3(time, flow),
            inlet_temperature=linear_history_v3(time, Tin),
            ambient_temperature=linear_history_v3(time, ambient),
        )
        result = simulate_v19(
            p, operating, time;
            initial_temperature=measured_initial_core_profile_v19(data, simulation_id),
            initial_perim_temperature=measured_initial_perim_profile_v19(data, simulation_id),
            initial_rear_temperature=measured_initial_rear_profile_v19(data, simulation_id),
            initial_cavity_temperature=T2[1],
            nodes=nodes, rear_nodes=rear_nodes, solver=solver,
            reltol=reltol, abstol=abstol, dtmax=dtmax,
        )
        outputs = extract_outputs_v19(result, p)
        return outputs, result, experiment
    end
end

begin # v19 objective helpers
    const ORDERING_SCALE_v19 = 100.0
    const ORDERING_WEIGHT_v19 = 0.25
    const TIMING_WEIGHT_v19 = 0.15
    const SLOPE_WEIGHT_v19 = 0.25

    function normalized_slope_mse_v19(model, experiment)
        length(model) < 2 && return 0.0
        scale = max(maximum(experiment) - minimum(experiment), 1.0)
        return mean(abs2, diff(model) .- diff(experiment)) / scale^2
    end

    function timing_penalty_v19(time, model, experiment)
        duration = max(time[end] - time[1], 1.0)
        return ((get_t90_v3(time, model) - get_t90_v3(time, experiment)) / duration)^2
    end

    function signal_loss_v19(time, model, experiment)
        return normalized_mse_v3(model, experiment) +
               SLOPE_WEIGHT_v19 * normalized_slope_mse_v19(model, experiment) +
               TIMING_WEIGHT_v19 * timing_penalty_v19(time, model, experiment)
    end

    function ordering_loss_v19(model, experiment)
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
        return mean(abs2, collect(model_offsets) .- collect(experiment_offsets)) / ORDERING_SCALE_v19^2
    end

    function capacitance_regularization_v19(p; nodes=15)
        capacitance = participating_assembly_heat_capacity_v19(p, nodes)
        return CAPACITANCE_REGULARIZATION_WEIGHT_v19 *
               ((capacitance - MEASURED_ASSEMBLY_CAPACITANCE_v19) / MEASURED_ASSEMBLY_CAPACITANCE_v19)^2
    end

    function loss_cases_v19(p, keys; is_cooling=false, nodes=15)
        total = 0.0
        for simulation_id in keys
            try
                model, result, experiment = solve_case_v19(
                    p, simulation_id; is_cooling=is_cooling, nodes=nodes,
                    reltol=2e-5, abstol=2e-6, dtmax=60.0,
                )
                all(isfinite, model) || return Inf
                total += mean(
                    signal_loss_v19(result.time, model[:, j], experiment[:, j])
                    for j in 1:7
                )
                if !is_cooling
                    total += ORDERING_WEIGHT_v19 * ordering_loss_v19(model, experiment)
                end
            catch err
                err isa InterruptException && rethrow()
                return Inf
            end
        end
        return total / length(keys)
    end

    function loss_heating_v19(p, keys=sim_key_heat; nodes=15)
        regularization = capacitance_regularization_v19(p; nodes=nodes)
        return loss_cases_v19(p, keys; is_cooling=false, nodes=nodes) + regularization
    end

    function loss_cooling_v19(p, keys=sim_key_cool; nodes=15)
        regularization = capacitance_regularization_v19(p; nodes=nodes)
        return loss_cases_v19(p, keys; is_cooling=true, nodes=nodes) + regularization
    end

    function calibrate_v19(; nodes=15, maximum_iterations=100, maximum_time_seconds=360.0)
        p0 = copy(pnew_v19)
        idx_fit = fit_heat_transfer_indices_v19
        function obj_free(theta)
            p_full = copy(p0)
            p_full[idx_fit] .= theta
            return loss_heating_v19(p_full; nodes=nodes) + loss_cooling_v19(p_full; nodes=nodes)
        end
        res = optimize_with_nlopt_v3(
            obj_free,
            pnew_v19[idx_fit],
            lb_full_v19[idx_fit],
            ub_full_v19[idx_fit];
            maximum_iterations=maximum_iterations,
            maximum_time_seconds=maximum_time_seconds,
            algorithm=NLopt.LN_BOBYQA(),
            label="fit-v19",
        )
        p_opt = copy(p0)
        p_opt[idx_fit] .= res.minimizer
        return (objective=res.minimum, parameters=p_opt, minimizer=res.minimizer, retcode=res.retcode)
    end
end

