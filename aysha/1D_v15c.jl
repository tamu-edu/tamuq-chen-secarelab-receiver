# ============================================================================
# 1D_v15c.jl - apparent heat transfer with fitted axial redistribution
# ============================================================================
# v15c keeps the v8b rear-tube/flange, predicted cavity/T2 structure, v12a
# source-distribution fit, and v14 per-irradiance absorbed-power scales. It
# refits the apparent assembly heat-transfer law and adds a bounded,
# energy-coupled axial redistribution term against the corrected thermocouple
# topology: T8/T12/T11 are side-wall targets; T9/T10 are interior flow-exposed
# diagnostics.
#
# Geometry updates supplied for v15c:
#   - Alumina tube length downstream of receiver: 150 mm.
#   - First 46 mm of tube is inside the cavity.
#   - Remaining 104 mm of tube is inside the water-cooled flange.
#   - T3 comparison uses the modeled gas temperature at z = 140 mm.
#   - Alumina adaptor diameter is interpreted as 77.6 mm.
#
# Rosseland, axial view-radiation, and the thermocouple measurement model remain
# available as optional diagnostics, but the default v15c audit keeps them off.
#
# The state vector is [receiver solid cells, rear alumina tube cells, T_cavity].
# ============================================================================

include("1D_v5.jl")

begin # v15c fixed constants
    const EPS_FRONT_FIXED_v15c = EPS_FRONT_FIXED_V5
    const ETA_ABS_FIXED_v15c = ETA_ABS_FIXED_V5
    const BETA_OPT_FIXED_v15c = BETA_OPT_FIXED_V5
    const BETA_OPT_MIN_v15c = 0.10
    const BETA_OPT_MAX_v15c = 300.0
    const H_FRONT_FIXED_v15c = 10.0
    const PRANDTL_EXPONENT_FIXED_v15c = 1.0 / 3.0
    const STANDARD_AIR_DENSITY_v15c = 101325.0 / (287.05 * 294.25)
    const DEFAULT_ROSSELAND_MEAN_PATH_v15c = 1.0e-3
    const RADIATION_AREA_v15c = A_frt
    const AXIAL_VIEW_RADIATION_EMISSIVITY_v15c = 0.85
    const AXIAL_VIEW_RADIATION_ROW_SUM_v15c = 0.06
    const AXIAL_VIEW_RADIATION_LENGTH_v15c = 0.050
    const AXIAL_VIEW_RADIATION_AREA_FACTOR_v15c = 1.0
    const TC_WIRE_LOSS_FRONT_FRACTION_v15c = 0.25
    const TC_WIRE_LOSS_LENGTH_v15c = 0.025
    const TC_WIRE_SINK_TEMPERATURE_v15c = 25.0 + 273.15
    const TEMPERATURE_AXIS_LIMITS_v15c = (250.0, 1250.0)
    const ETA_OPT_FIXED_v15c = 1.0
    const FRONT_DEPOSITION_INITIAL_v15c = 1.0
    const POWER_SCALE_MIN_v15c = 0.25
    const POWER_SCALE_MAX_v15c = 4.00
    const IRRADIANCE_LEVELS_v15c = (456000.0, 304000.0, 256000.0)
    const ORDERING_WEIGHT_v15c = 0.25
    const ORDERING_SCALE_v15c = 100.0

    # Cavity geometry from the annotated v15c schematic.
    const CAVITY_OUTER_RADIUS_v15c = 75.0e-3
    const METAL_THICKNESS_v15c = 18.0e-3
    const INSULATION_OUTER_RADIUS_v15c = CAVITY_OUTER_RADIUS_v15c - METAL_THICKNESS_v15c
    const METAL_OUTER_RADIUS_v15c = CAVITY_OUTER_RADIUS_v15c
    const CAVITY_LENGTH_v15c = 165.0e-3
    const T3_SAMPLE_POSITION_v15c = 140.0e-3

    const ADAPTOR_DIAMETER_v15c = 77.6e-3
    const ADAPTOR_RADIUS_v15c = ADAPTOR_DIAMETER_v15c / 2.0
    const ADAPTOR_LENGTH_v15c = 57.0e-3
    const ADAPTOR_TUBE_RADIUS_v15c = 6.5e-3
    const ADAPTOR_RECEIVER_OVERLAP_LENGTH_v15c = 28.0e-3
    const ADAPTOR_TUBE_OVERLAP_LENGTH_v15c = ADAPTOR_LENGTH_v15c - ADAPTOR_RECEIVER_OVERLAP_LENGTH_v15c
    const ADAPTOR_OVERLAP_LENGTH_v15c = ADAPTOR_RECEIVER_OVERLAP_LENGTH_v15c
    const ADAPTOR_FREE_LENGTH_v15c = ADAPTOR_TUBE_OVERLAP_LENGTH_v15c
    const METAL_LENGTH_v15c = CAVITY_LENGTH_v15c
    const BACKPLATE_THICKNESS_v15c = METAL_THICKNESS_v15c
    const INSULATION_MONOLITH_LENGTH_v15c = L
    const INSULATION_ADAPTOR_LENGTH_v15c = ADAPTOR_LENGTH_v15c
    const INSULATION_CONDUCTIVITY_v15c = 0.08

    const RECEIVER_EQ_RADIUS_v15c = sqrt(A_frt / pi)
    const RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v15c =
        4.0 * pi * INSULATION_CONDUCTIVITY_v15c /
        log(INSULATION_OUTER_RADIUS_v15c / RECEIVER_EQ_RADIUS_v15c)
    const ADAPTOR_TO_CAVITY_RESISTANCE_v15c =
        log(INSULATION_OUTER_RADIUS_v15c / ADAPTOR_RADIUS_v15c) /
        (4.0 * pi * INSULATION_CONDUCTIVITY_v15c * ADAPTOR_LENGTH_v15c)
    const ADAPTOR_CONTACT_RESISTANCE_v15c =
        1.0 / (100.0 * 4.0 * w_t * ADAPTOR_OVERLAP_LENGTH_v15c)

    # Rear extension geometry supplied for v15c.
    const REAR_TUBE_LENGTH_v15c = 150.0e-3
    const REAR_TUBE_CAVITY_LENGTH_v15c = 46.0e-3
    const REAR_TUBE_FLANGE_LENGTH_v15c = REAR_TUBE_LENGTH_v15c - REAR_TUBE_CAVITY_LENGTH_v15c
    const REAR_TUBE_DEFAULT_NODES_v15c = 12
    const REAR_TUBE_GAS_RADIUS_v15c = ADAPTOR_TUBE_RADIUS_v15c
    const REAR_TUBE_WALL_THICKNESS_v15c = 1.5e-3
    const REAR_TUBE_OUTER_RADIUS_v15c = REAR_TUBE_GAS_RADIUS_v15c + REAR_TUBE_WALL_THICKNESS_v15c
    const REAR_TUBE_FLOW_AREA_v15c = pi * REAR_TUBE_GAS_RADIUS_v15c^2
    const REAR_TUBE_HYDRAULIC_DIAMETER_v15c = 2.0 * REAR_TUBE_GAS_RADIUS_v15c
    const REAR_TUBE_EXCHANGE_PERIMETER_v15c = 2.0 * pi * REAR_TUBE_GAS_RADIUS_v15c

    const RECEIVER_TO_TUBE_CONDUCTANCE_v15c =
        1.0 / ADAPTOR_CONTACT_RESISTANCE_v15c

    # Alumina and aluminum properties follow the older 0D material set. The
    # flange is treated as an isothermal water-cooled sink in v15c.
    const ALUMINA_DENSITY_v15c = 3900.0
    const ALUMINUM_DENSITY_v15c = 2700.0
    const ALUMINUM_HEAT_CAPACITY_v15c = 900.0
    const ALUMINUM_CONDUCTIVITY_v15c = 205.0
    const METAL_EMISSIVITY_v15c = 0.2
    const H_NAT_CAVITY_v15c = 10.0
    const INSULATION_DENSITY_v15c = 140.0
    const INSULATION_HEAT_CAPACITY_v15c = 1360.0
    const WATER_FLANGE_TEMPERATURE_v15c = 25.0 + 273.15

    const CAVITY_FELT_VOLUME_v15c = max(
        0.0,
        pi * INSULATION_OUTER_RADIUS_v15c^2 * CAVITY_LENGTH_v15c -
        A_frt * L -
        pi * ADAPTOR_RADIUS_v15c^2 * ADAPTOR_LENGTH_v15c,
    )
    const CAVITY_METAL_VOLUME_v15c =
        pi * (METAL_OUTER_RADIUS_v15c^2 - INSULATION_OUTER_RADIUS_v15c^2) *
        CAVITY_LENGTH_v15c +
        pi * METAL_OUTER_RADIUS_v15c^2 * BACKPLATE_THICKNESS_v15c
    const CAVITY_OUTER_AREA_v15c =
        2.0 * pi * METAL_OUTER_RADIUS_v15c * CAVITY_LENGTH_v15c +
        pi * METAL_OUTER_RADIUS_v15c^2
    const CAVITY_HEAT_CAPACITY_v15c =
        INSULATION_DENSITY_v15c * INSULATION_HEAT_CAPACITY_v15c *
        CAVITY_FELT_VOLUME_v15c +
        ALUMINUM_DENSITY_v15c * ALUMINUM_HEAT_CAPACITY_v15c *
        CAVITY_METAL_VOLUME_v15c
end

# Final v15c diagnostic parameter vector p[1:12]
#   p[1]  A_Nu_app      apparent assembly Nu prefactor
#   p[2]  B_Re          apparent assembly Reynolds exponent
#   p[3]  C_Pr          fixed Prandtl exponent (= 1/3)
#   p[4]  reserved      fixed placeholder for vector compatibility
#   p[5]  beta_opt      frozen v12a fitted optical attenuation [1/m]
#   p[6]  front_dep     frozen v12a fitted front-cell deposition fraction
#   p[7]  scale_456     fixed v14 absorbed-power scale at 456 kW/m2
#   p[8]  scale_304     fixed v14 absorbed-power scale at 304 kW/m2
#   p[9]  scale_256     fixed v14 absorbed-power scale at 256 kW/m2
#   p[10] k_axial_ref   fitted axial redistribution conductivity scale [W/m/K]
#   p[11] axial_exp     axial shape exponent for redistribution
#   p[12] Re_exp        Reynolds exponent for redistribution
begin # v15c parameter values and bounds
    pnew_v15c = [
        5.296511755883694e-4,
        2.40,
        PRANDTL_EXPONENT_FIXED_v15c,
        1.0,
        184.67243519237965,
        0.0,
        1.669539299381132,
        1.7170694424514947,
        0.9827026344016692,
        10.0,
        1.00,
        0.00,
    ]
    lb_full_v15c = [
        1.0e-5,
        0.20,
        pnew_v15c[3],
        pnew_v15c[4],
        pnew_v15c[5],
        pnew_v15c[6],
        pnew_v15c[7],
        pnew_v15c[8],
        pnew_v15c[9],
        0.00,
        0.20,
        -1.00,
    ]
    ub_full_v15c = [
        1.0e-2,
        3.00,
        pnew_v15c[3],
        pnew_v15c[4],
        pnew_v15c[5],
        pnew_v15c[6],
        pnew_v15c[7],
        pnew_v15c[8],
        pnew_v15c[9],
        200.00,
        4.00,
        1.00,
    ]

    fit_heat_transfer_indices_v15c = [1, 2, 10, 11, 12]
    fit_source_indices_v15c = Int[]
    fit_power_scale_indices_v15c = Int[]
    fit_radiation_indices_v15c = Int[]
end

function irradiance_level_index_v15c(irradiance)
    irradiance <= 0.0 && return 0
    differences = abs.(collect(IRRADIANCE_LEVELS_v15c) .- Float64(irradiance))
    return argmin(differences)
end

power_scale_index_v15c(irradiance) = begin
    level_index = irradiance_level_index_v15c(irradiance)
    level_index == 0 ? 0 : 6 + level_index
end

function absorbed_power_scale_v15c(irradiance, p)
    index = power_scale_index_v15c(irradiance)
    return index == 0 ? 1.0 : p[index]
end

irradiance_factor_v15c(irradiance, p) =
    ETA_OPT_FIXED_v15c * absorbed_power_scale_v15c(irradiance, p)

function receiver_nusselt_v15c(z, reynolds, prandtl, p)
    return max(0.0, p[1]) *
           max(reynolds, eps(Float64))^p[2] *
           max(prandtl, eps(Float64))^p[3]
end

m_dot_standard_v15c(flow_lpm) =
    STANDARD_AIR_DENSITY_v15c * max(0.0, flow_lpm) / 60000.0

smoothstep_v15c(x) = begin
    xi = clamp(x, 0.0, 1.0)
    xi^2 * (3.0 - 2.0 * xi)
end

function rosseland_mean_path_v15c(z, p; profile=:uniform)
    ell_uniform = DEFAULT_ROSSELAND_MEAN_PATH_v15c
    xi = smoothstep_v15c(z / L)
    if profile == :uniform
        return ell_uniform
    elseif profile == :front_strong
        ell_front = 3.0 * ell_uniform
        ell_rear = ell_uniform / 3.0
        return ell_front + (ell_rear - ell_front) * xi
    elseif profile == :rear_strong
        ell_front = ell_uniform / 3.0
        ell_rear = 3.0 * ell_uniform
        return ell_front + (ell_rear - ell_front) * xi
    elseif profile == :weak_gradient
        ell_front = ell_uniform
        ell_rear = ell_uniform / 3.0
        return ell_front + (ell_rear - ell_front) * xi
    else
        error("Unknown Rosseland axial profile: $profile")
    end
end

rosseland_beta_v15c(p; z=0.0, profile=:uniform) =
    1.0 / rosseland_mean_path_v15c(z, p; profile=profile)

rosseland_conductivity_v15c(T, p; z=0.0, profile=:uniform) =
    16.0 * sigma_sb * property_temperature(T)^3 /
    (3.0 * rosseland_beta_v15c(p; z=z, profile=profile))

rosseland_optical_thickness_v15c(length_scale, p) =
    rosseland_beta_v15c(p) * length_scale

function axial_redistribution_shape_v15c(zface, p)
    xi = clamp(zface / L, 0.0, 1.0)
    return xi^max(p[11], 1e-6)
end

function axial_redistribution_conductivity_v15c(Tface, zface, reynolds, p)
    temperature_factor = (property_temperature(Tface) / 900.0)^3
    reynolds_factor = max(reynolds, eps(Float64))^p[12] / 50.0^p[12]
    return max(0.0, p[10]) *
           temperature_factor *
           reynolds_factor *
           axial_redistribution_shape_v15c(zface, p)
end

function axial_view_matrix_v15c(z)
    nodes = length(z)
    raw = zeros(nodes, nodes)
    for i in 1:nodes, j in 1:nodes
        i == j && continue
        separation = abs(z[i] - z[j])
        geometry_decay = exp(-separation / AXIAL_VIEW_RADIATION_LENGTH_v15c)
        angular_decay = 1.0 / (1.0 + (separation / max(Dh, eps(Float64)))^2)
        raw[i, j] = geometry_decay * angular_decay
    end
    row_sum = maximum(sum(raw, dims=2))
    row_sum <= eps(Float64) && return raw
    return AXIAL_VIEW_RADIATION_ROW_SUM_v15c .* raw ./ row_sum
end

function axial_view_radiation_net_v15c(Ts, view_matrix, dx; strength=1.0)
    nodes = length(Ts)
    qnet = zeros(nodes)
    area = AXIAL_VIEW_RADIATION_AREA_FACTOR_v15c * P_exchange * dx
    coefficient = strength * AXIAL_VIEW_RADIATION_EMISSIVITY_v15c * sigma_sb * area
    for i in 1:(nodes - 1), j in (i + 1):nodes
        view_factor = view_matrix[i, j]
        view_factor <= 0.0 && continue
        Qij = coefficient * view_factor * (Ts[i]^4 - Ts[j]^4)
        qnet[i] -= Qij
        qnet[j] += Qij
    end
    return qnet
end

alumina_conductivity_v15c(T) =
    5.5 + 34.5 * exp(-0.0033 * (property_temperature(T) - 273.15))

alumina_heat_capacity_v15c(T) =
    max(500.0, (1.00446 + 1.742e-4 * property_temperature(T) -
                2.796e4 * property_temperature(T)^-2) * 1000.0)

function rear_tube_outer_radius_v15c(z_rear)
    return z_rear <= ADAPTOR_TUBE_OVERLAP_LENGTH_v15c ?
           ADAPTOR_RADIUS_v15c : REAR_TUBE_OUTER_RADIUS_v15c
end

function rear_tube_solid_area_v15c(z_rear)
    outer_radius = rear_tube_outer_radius_v15c(z_rear)
    return pi * (outer_radius^2 - REAR_TUBE_GAS_RADIUS_v15c^2)
end

function rear_tube_cavity_conductance_per_length_v15c(z_rear, p)
    z_rear > REAR_TUBE_CAVITY_LENGTH_v15c && return 0.0
    outer_radius = min(rear_tube_outer_radius_v15c(z_rear), 0.98 * INSULATION_OUTER_RADIUS_v15c)
    return 2.0 * pi * INSULATION_CONDUCTIVITY_v15c /
           log(INSULATION_OUTER_RADIUS_v15c / outer_radius)
end

function rear_tube_flange_conductance_per_length_v15c(z_rear, Ttube)
    z_rear <= REAR_TUBE_CAVITY_LENGTH_v15c && return 0.0
    tube_wall_resistance =
        log(REAR_TUBE_OUTER_RADIUS_v15c / REAR_TUBE_GAS_RADIUS_v15c) /
        (2.0 * pi * alumina_conductivity_v15c(Ttube))
    aluminum_resistance =
        log(METAL_OUTER_RADIUS_v15c / REAR_TUBE_OUTER_RADIUS_v15c) /
        (2.0 * pi * ALUMINUM_CONDUCTIVITY_v15c)
    return 1.0 / (tube_wall_resistance + aluminum_resistance)
end

function cavity_loss_to_ambient_v15c(Tcavity, ambient)
    return H_NAT_CAVITY_v15c * CAVITY_OUTER_AREA_v15c * (Tcavity - ambient) +
           METAL_EMISSIVITY_v15c * sigma_sb * CAVITY_OUTER_AREA_v15c *
           (Tcavity^4 - ambient^4)
end

function rear_tube_capacity_v15c(Ttube, z_rear, dx_rear)
    return ALUMINA_DENSITY_v15c * alumina_heat_capacity_v15c(Ttube) *
           rear_tube_solid_area_v15c(z_rear) * dx_rear
end

function rear_tube_htc_v15c(Twall, Tgas, flow_lpm, Tin)
    flow_lpm <= 1e-12 && return 0.0
    mdot = m_dot_standard_v15c(flow_lpm)
    Tfilm = 0.5 * (Twall + Tgas)
    cp = cpf_f(Tfilm)
    mu = mu_f_f(Tfilm)
    kf = kf_f(Tfilm)
    reynolds = mdot * REAR_TUBE_HYDRAULIC_DIAMETER_v15c /
               (REAR_TUBE_FLOW_AREA_v15c * mu)
    prandtl = cp * mu / kf
    nusselt = reynolds < 2300.0 ?
              3.66 :
              0.023 * reynolds^0.8 * prandtl^0.4
    return nusselt * kf / REAR_TUBE_HYDRAULIC_DIAMETER_v15c
end

function gas_profile_v15c!(Tg, Qgas, hcell, Qgas_rear, hrear,
                         Ts, Ttube, time, p, operating, z, dx, z_rear, dx_rear)
    flow = max(0.0, operating.flow_lpm(time))
    Tin = operating.inlet_temperature(time)
    receiver_nodes = length(Ts)
    Tg[1] = Tin

    if flow <= 1e-12
        fill!(Qgas, 0.0)
        fill!(hcell, 0.0)
        fill!(Qgas_rear, 0.0)
        fill!(hrear, 0.0)
        Tg[1] = Tin
        for i in eachindex(Ts)
            Tg[i + 1] = Ts[i]
        end
        for j in eachindex(Ttube)
            Tg[receiver_nodes + j + 1] = Ttube[j]
        end
        return nothing
    end

    mdot = m_dot_standard_v15c(flow)
    for i in eachindex(Ts)
        Tfilm = 0.5 * (Ts[i] + Tg[i])
        cp = cpf_f(Tfilm)
        mu = mu_f_f(Tfilm)
        kf = kf_f(Tfilm)
        reynolds = mdot * Dh / (A_flow * mu)
        prandtl = cp * mu / kf
        nusselt = receiver_nusselt_v15c(z[i], reynolds, prandtl, p)
        hcell[i] = nusselt * kf / Dh
        UA = hcell[i] * P_exchange * dx
        effectiveness = -expm1(-UA / (mdot * cp))
        Tg[i + 1] = Tg[i] + effectiveness * (Ts[i] - Tg[i])
        Qgas[i] = mdot * cp * (Tg[i + 1] - Tg[i])
    end

    for j in eachindex(Ttube)
        inlet_index = receiver_nodes + j
        Tfilm = 0.5 * (Ttube[j] + Tg[inlet_index])
        cp = cpf_f(Tfilm)
        hrear[j] = rear_tube_htc_v15c(Ttube[j], Tg[inlet_index], flow, Tin)
        UA = hrear[j] * REAR_TUBE_EXCHANGE_PERIMETER_v15c * dx_rear
        effectiveness = -expm1(-UA / (mdot * cp))
        Tg[inlet_index + 1] = Tg[inlet_index] + effectiveness *
                              (Ttube[j] - Tg[inlet_index])
        Qgas_rear[j] = mdot * cp * (Tg[inlet_index + 1] - Tg[inlet_index])
    end
    return nothing
end

function receiver_rhs_v15c!(du, u, time, p, operating, z, dx, solar_weights,
                          z_rear, dx_rear,
                          Tg, Qgas, hcell, Qgas_rear, hrear,
                          use_rosseland_radiation, rosseland_profile,
                          use_axial_view_radiation, axial_view_matrix,
                          axial_view_strength)
    nodes = length(z)
    rear_nodes = length(z_rear)
    Ts = view(u, 1:nodes)
    Ttube = view(u, (nodes + 1):(nodes + rear_nodes))
    Tcavity = u[nodes + rear_nodes + 1]
    dTs = view(du, 1:nodes)
    dTtube = view(du, (nodes + 1):(nodes + rear_nodes))
    ambient = operating.ambient_temperature(time)

    gas_profile_v15c!(Tg, Qgas, hcell, Qgas_rear, hrear,
                    Ts, Ttube, time, p, operating, z, dx, z_rear, dx_rear)
    fill!(du, 0.0)
    flow_now = max(0.0, operating.flow_lpm(time))
    mdot_now = m_dot_standard_v15c(flow_now)

    for i in 1:(nodes - 1)
        ki = ks_f(Ts[i])
        kj = ks_f(Ts[i + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        Qcond = kface * A_solid * (Ts[i] - Ts[i + 1]) / dx
        dTs[i] -= Qcond
        dTs[i + 1] += Qcond
        Tface = 0.5 * (Ts[i] + Ts[i + 1])
        zface = 0.5 * (z[i] + z[i + 1])
        if flow_now > 1e-12
            reynolds_face = mdot_now * Dh / (A_flow * mu_f_f(Tface))
            kredist = axial_redistribution_conductivity_v15c(
                Tface, zface, reynolds_face, p,
            )
            Qredist = kredist * RADIATION_AREA_v15c * (Ts[i] - Ts[i + 1]) / dx
            dTs[i] -= Qredist
            dTs[i + 1] += Qredist
        end
        if use_rosseland_radiation
            Qrad = rosseland_conductivity_v15c(
                Tface, p; z=zface, profile=rosseland_profile,
            ) * RADIATION_AREA_v15c *
                   (Ts[i] - Ts[i + 1]) / dx
            dTs[i] -= Qrad
            dTs[i + 1] += Qrad
        end
    end

    if use_axial_view_radiation
        Qview = axial_view_radiation_net_v15c(
            Ts, axial_view_matrix, dx; strength=axial_view_strength,
        )
        for i in 1:nodes
            dTs[i] += Qview[i]
        end
    end

    irradiance = max(0.0, operating.irradiance(time))
    Qsolar = ETA_ABS_FIXED_v15c * irradiance_factor_v15c(irradiance, p) *
             irradiance * A_frt
    radial_conductance = RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v15c
    for i in 1:nodes
        Qside_cavity = radial_conductance * dx * (Ts[i] - Tcavity)
        dTs[i] += Qsolar * solar_weights[i] - Qgas[i] - Qside_cavity
        du[nodes + rear_nodes + 1] += Qside_cavity
    end

    Qfront = H_FRONT_FIXED_v15c * A_frt * (Ts[1] - ambient) +
             EPS_FRONT_FIXED_v15c * sigma_sb * A_frt * (Ts[1]^4 - ambient^4)
    Qreceiver_tube = RECEIVER_TO_TUBE_CONDUCTANCE_v15c * (Ts[end] - Ttube[1])

    dTs[1] -= Qfront
    dTs[end] -= Qreceiver_tube
    dTtube[1] += Qreceiver_tube

    for j in 1:(rear_nodes - 1)
        ki = alumina_conductivity_v15c(Ttube[j])
        kj = alumina_conductivity_v15c(Ttube[j + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        Aface = 0.5 * (rear_tube_solid_area_v15c(z_rear[j]) +
                       rear_tube_solid_area_v15c(z_rear[j + 1]))
        Qcond = kface * Aface * (Ttube[j] - Ttube[j + 1]) / dx_rear
        dTtube[j] -= Qcond
        dTtube[j + 1] += Qcond
    end

    for j in 1:rear_nodes
        G_cavity = rear_tube_cavity_conductance_per_length_v15c(z_rear[j], p)
        G_flange = rear_tube_flange_conductance_per_length_v15c(z_rear[j], Ttube[j])
        Q_cavity = G_cavity * dx_rear * (Ttube[j] - Tcavity)
        Q_flange = G_flange * dx_rear * (Ttube[j] - WATER_FLANGE_TEMPERATURE_v15c)
        dTtube[j] -= Qgas_rear[j] + Q_cavity + Q_flange
        du[nodes + rear_nodes + 1] += Q_cavity
    end

    Qcavity_ambient = cavity_loss_to_ambient_v15c(Tcavity, ambient)
    du[nodes + rear_nodes + 1] -= Qcavity_ambient

    cell_volume = A_solid * dx
    for i in 1:nodes
        capacity = rho_s * Cps_f(Ts[i]) * cell_volume
        dTs[i] /= capacity
    end
    for j in 1:rear_nodes
        dTtube[j] /= max(rear_tube_capacity_v15c(Ttube[j], z_rear[j], dx_rear),
                         eps(Float64))
    end
    du[nodes + rear_nodes + 1] /= max(CAVITY_HEAT_CAPACITY_v15c, eps(Float64))
    return nothing
end

function receiver_ode_v15c!(dTs, Ts, context, time)
    receiver_rhs_v15c!(
        dTs, Ts, time, context.p, context.operating,
        context.z, context.dx, context.solar_weights,
        context.z_rear, context.dx_rear,
        context.Tg, context.Qgas, context.hcell, context.Qgas_rear, context.hrear,
        context.use_rosseland_radiation, context.rosseland_profile,
        context.use_axial_view_radiation, context.axial_view_matrix,
        context.axial_view_strength,
    )
    return nothing
end

function simulate_v15c(p, operating, save_times;
                     initial_temperature=Tamb, nodes=default_nodes,
                     rear_nodes=REAR_TUBE_DEFAULT_NODES_v15c,
                     initial_rear_temperature=nothing,
                     initial_cavity_temperature=Tamb,
                     use_rosseland_radiation=false,
                     rosseland_profile=:uniform,
                     use_axial_view_radiation=false,
                     axial_view_strength=1.0,
                     solver=Rodas5P(autodiff=AutoFiniteDiff()),
                     reltol=1e-6, abstol=1e-7, dtmax=30.0)
    length(p) == 12 || throw(ArgumentError("1D_v15c expects 12 parameters"))
    all(isfinite, p) || throw(ArgumentError("all model parameters must be finite"))
    all(>(0.0), p[1:5]) && 0.0 <= p[6] <= 1.0 && all(>(0.0), p[7:9]) &&
        p[10] >= 0.0 && p[11] > 0.0 ||
        throw(ArgumentError("Nu coefficients, beta_opt, power scales, and axial redistribution parameters must be finite and within bounds"))
    nodes >= 3 || throw(ArgumentError("at least three axial nodes are required"))
    rear_nodes >= 2 || throw(ArgumentError("at least two rear tube nodes are required"))

    time = Float64.(save_times)
    isempty(time) && throw(ArgumentError("save_times cannot be empty"))
    (length(time) == 1 || all(diff(time) .> 0.0)) ||
        throw(ArgumentError("save_times must be strictly increasing"))

    dx = L / nodes
    dx_rear = REAR_TUBE_LENGTH_v15c / rear_nodes
    z = collect(range(dx / 2, step=dx, length=nodes))
    z_rear = collect(range(dx_rear / 2, step=dx_rear, length=rear_nodes))
    zg = vcat(
        collect(range(0.0, L, length=nodes + 1)),
        L .+ collect(range(dx_rear, step=dx_rear, length=rear_nodes)),
    )
    weights = solar_weights_v5(
        p[5], p[6], nodes,
    )
    axial_view_matrix = axial_view_matrix_v15c(z)
    Ts_initial = initial_profile_v3(initial_temperature, z)
    Ttube_initial = if isnothing(initial_rear_temperature)
        fill(Ts_initial[end], rear_nodes)
    elseif initial_rear_temperature isa Function
        Float64.([initial_rear_temperature(zr) for zr in z_rear])
    else
        fill(Float64(initial_rear_temperature), rear_nodes)
    end
    u_initial = vcat(Ts_initial, Ttube_initial, Float64(initial_cavity_temperature))

    Tg_history = Matrix{Float64}(undef, nodes + rear_nodes + 1, length(time))
    h_history = Matrix{Float64}(undef, nodes, length(time))
    nu_history = Matrix{Float64}(undef, nodes, length(time))
    effectiveness_history = Matrix{Float64}(undef, nodes, length(time))
    ua_over_mcp_history = Matrix{Float64}(undef, nodes, length(time))
    qgas_history = Matrix{Float64}(undef, nodes, length(time))
    axial_view_radiation_history = Matrix{Float64}(undef, nodes, length(time))
    h_rear_history = Matrix{Float64}(undef, rear_nodes, length(time))
    rear_temperature_history = Matrix{Float64}(undef, rear_nodes, length(time))
    cavity_temperature_history = Vector{Float64}(undef, length(time))
    receiver_to_tube_heat_history = Vector{Float64}(undef, length(time))
    receiver_to_cavity_heat_history = Vector{Float64}(undef, length(time))
    tube_to_cavity_heat_history = Vector{Float64}(undef, length(time))
    cavity_ambient_heat_history = Vector{Float64}(undef, length(time))
    flange_heat_history = Vector{Float64}(undef, length(time))
    Tg = zeros(nodes + rear_nodes + 1)
    Qgas = zeros(nodes)
    hcell = zeros(nodes)
    Qgas_rear = zeros(rear_nodes)
    hrear = zeros(rear_nodes)

    context = (
        p=p, operating=operating, z=z, dx=dx, solar_weights=weights,
        z_rear=z_rear, dx_rear=dx_rear,
        Tg=Tg, Qgas=Qgas, hcell=hcell, Qgas_rear=Qgas_rear, hrear=hrear,
        use_rosseland_radiation=use_rosseland_radiation,
        rosseland_profile=rosseland_profile,
        use_axial_view_radiation=use_axial_view_radiation,
        axial_view_matrix=axial_view_matrix,
        axial_view_strength=axial_view_strength,
    )
    problem = ODEProblem(receiver_ode_v15c!, u_initial, (time[1], time[end]), context)
    solution = solve(
        problem, solver; saveat=time, save_everystep=false, dense=false,
        reltol=reltol, abstol=abstol, dtmax=dtmax, tstops=time,
    )
    successful_retcode(solution) ||
        error("1D_v15c ODE solve failed with retcode $(solution.retcode)")

    state_history = reduce(hcat, solution.u)
    Ts_history = state_history[1:nodes, :]
    rear_temperature_history .= state_history[(nodes + 1):(nodes + rear_nodes), :]
    cavity_temperature_history .= state_history[nodes + rear_nodes + 1, :]
    for output_index in eachindex(time)
        t = time[output_index]
        Ts_now = view(Ts_history, :, output_index)
        Ttube_now = view(rear_temperature_history, :, output_index)
        Tcavity_now = cavity_temperature_history[output_index]
        gas_profile_v15c!(Tg, Qgas, hcell, Qgas_rear, hrear, Ts_now, Ttube_now,
                        t, p, operating, z, dx, z_rear, dx_rear)
        Tg_history[:, output_index] .= Tg
        h_history[:, output_index] .= hcell
        qgas_history[:, output_index] .= Qgas
        axial_view_radiation_history[:, output_index] .=
            use_axial_view_radiation ?
            axial_view_radiation_net_v15c(
                Ts_now, axial_view_matrix, dx; strength=axial_view_strength,
            ) :
            zeros(nodes)
        flow_now = max(0.0, operating.flow_lpm(t))
        Tin_now = operating.inlet_temperature(t)
        if flow_now <= 1e-12
            fill!(view(nu_history, :, output_index), 0.0)
            fill!(view(effectiveness_history, :, output_index), 0.0)
            fill!(view(ua_over_mcp_history, :, output_index), 0.0)
        else
            mdot_now = m_dot_standard_v15c(flow_now)
            for i in 1:nodes
                Tfilm = 0.5 * (Ts_now[i] + Tg[i])
                cp = cpf_f(Tfilm)
                kf = kf_f(Tfilm)
                ratio = hcell[i] * P_exchange * dx / (mdot_now * cp)
                nu_history[i, output_index] = hcell[i] * Dh / kf
                ua_over_mcp_history[i, output_index] = ratio
                effectiveness_history[i, output_index] = -expm1(-ratio)
            end
        end
        h_rear_history[:, output_index] .= hrear
        receiver_to_tube_heat_history[output_index] =
            RECEIVER_TO_TUBE_CONDUCTANCE_v15c * (Ts_now[end] - Ttube_now[1])
        receiver_to_cavity_heat_history[output_index] = sum(
            RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v15c * dx *
            (Ts_now[i] - Tcavity_now)
            for i in 1:nodes
        )
        tube_to_cavity_heat_history[output_index] = sum(
            rear_tube_cavity_conductance_per_length_v15c(z_rear[j], p) * dx_rear *
            (Ttube_now[j] - Tcavity_now)
            for j in 1:rear_nodes
        )
        cavity_ambient_heat_history[output_index] =
            cavity_loss_to_ambient_v15c(Tcavity_now, operating.ambient_temperature(t))
        flange_heat_history[output_index] = sum(
            rear_tube_flange_conductance_per_length_v15c(z_rear[j], Ttube_now[j]) *
            dx_rear * (Ttube_now[j] - WATER_FLANGE_TEMPERATURE_v15c)
            for j in 1:rear_nodes
        )
    end

    return (
        time=time, z_solid=z, z_rear_tube=L .+ z_rear, z_gas=zg,
        solid_temperature=Ts_history,
        rear_tube_temperature=rear_temperature_history,
        cavity_temperature=cavity_temperature_history,
        boundary_temperature=cavity_temperature_history,
        gas_temperature=Tg_history,
        heat_transfer_coefficient=h_history,
        receiver_nusselt=nu_history,
        receiver_effectiveness=effectiveness_history,
        receiver_ua_over_mcp=ua_over_mcp_history,
        receiver_gas_heat=qgas_history,
        axial_view_radiation_heat=axial_view_radiation_history,
        axial_view_matrix=axial_view_matrix,
        rear_tube_heat_transfer_coefficient=h_rear_history,
        receiver_to_tube_heat=receiver_to_tube_heat_history,
        receiver_to_cavity_heat=receiver_to_cavity_heat_history,
        tube_to_cavity_heat=tube_to_cavity_heat_history,
        cavity_ambient_heat_loss=cavity_ambient_heat_history,
        flange_heat_loss=flange_heat_history,
        ode_solution=solution,
    )
end

begin # v15c experimental interface
    const WALL_CHAIN_POSITIONS_v15c = (
        T8=sensor_positions[:T8],
        T12=sensor_positions[:T9],
        T11=sensor_positions[:T10],
    )

    function measured_initial_profile_v15c(data, simulation_id)
        z8 = WALL_CHAIN_POSITIONS_v15c.T8
        z12 = WALL_CHAIN_POSITIONS_v15c.T12
        z11 = WALL_CHAIN_POSITIONS_v15c.T11
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

    function measured_initial_rear_profile_v15c(data, simulation_id)
        receiver_profile = measured_initial_profile_v15c(data, simulation_id)
        T_receiver_exit = receiver_profile(L)
        T_tube_exit = observation(data, simulation_id, "_Tf")[1]
        return function (z_rear)
            fraction = clamp(z_rear / REAR_TUBE_LENGTH_v15c, 0.0, 1.0)
            return T_receiver_exit + fraction * (T_tube_exit - T_receiver_exit)
        end
    end

    function experimental_case_v15c(simulation_id; is_cooling=false)
        data = is_cooling ? measurements_cooling : measurements
        conditions = is_cooling ? simulation_conditions_cooling : simulation_conditions
        time = observation_time(data, simulation_id)
        experiment = hcat(
            observation(data, simulation_id, "_T8"),
            observation(data, simulation_id, "_T12"),
            observation(data, simulation_id, "_T11"),
            observation(data, simulation_id, "_Tf"),
            observation(data, simulation_id, "_T2"),
        )
        return conditions[simulation_id], time, experiment, data
    end

    function gas_at_position_v15c(result, position)
        z = result.z_gas
        values = result.gas_temperature
        position <= z[1] && return vec(values[1, :])
        position >= z[end] && return vec(values[end, :])
        i = searchsortedlast(z, position)
        fraction = (position - z[i]) / (z[i + 1] - z[i])
        return vec((1.0 - fraction) .* values[i, :] .+ fraction .* values[i + 1, :])
    end

    tc_wire_loss_fraction_v15c(position) =
        TC_WIRE_LOSS_FRONT_FRACTION_v15c *
        exp(-max(position, 0.0) / TC_WIRE_LOSS_LENGTH_v15c)

    function thermocouple_at_v15c(result, position; use_tc_measurement_model=false)
        solid = solid_at_v3(result, position)
        use_tc_measurement_model || return solid
        fraction = tc_wire_loss_fraction_v15c(position)
        return solid .- fraction .* (solid .- TC_WIRE_SINK_TEMPERATURE_v15c)
    end

    function extract_outputs_v15c(result, p, T3_initial;
                                use_tc_measurement_model=false)
        T8 = thermocouple_at_v15c(
            result, sensor_positions[:T8];
            use_tc_measurement_model=use_tc_measurement_model,
        )
        T12_wall = thermocouple_at_v15c(
            result, WALL_CHAIN_POSITIONS_v15c.T12;
            use_tc_measurement_model=use_tc_measurement_model,
        )
        T11_wall = thermocouple_at_v15c(
            result, WALL_CHAIN_POSITIONS_v15c.T11;
            use_tc_measurement_model=use_tc_measurement_model,
        )
        Tf = gas_at_position_v15c(result, T3_SAMPLE_POSITION_v15c)
        T2 = result.cavity_temperature
        h_mean = vec(mean(result.heat_transfer_coefficient, dims=1))
        return hcat(T8, T12_wall, T11_wall, Tf, T2, h_mean)
    end

    function solve_case_v15c(p, simulation_id; is_cooling=false,
                           nodes=default_nodes,
                           rear_nodes=REAR_TUBE_DEFAULT_NODES_v15c,
                           use_rosseland_radiation=false,
                           rosseland_profile=:uniform,
                           use_axial_view_radiation=false,
                           axial_view_strength=1.0,
                           use_tc_measurement_model=false,
                           solver=Rodas5P(autodiff=AutoFiniteDiff()),
                           reltol=1e-6, abstol=1e-7, dtmax=30.0)
        conditions, time, experiment, data = experimental_case_v15c(
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
        result = simulate_v15c(
            p, operating, time;
            initial_temperature=measured_initial_profile_v15c(data, simulation_id),
            initial_rear_temperature=measured_initial_rear_profile_v15c(data, simulation_id),
            initial_cavity_temperature=T2[1],
            use_rosseland_radiation=use_rosseland_radiation,
            rosseland_profile=rosseland_profile,
            use_axial_view_radiation=use_axial_view_radiation,
            axial_view_strength=axial_view_strength,
            nodes=nodes, rear_nodes=rear_nodes, solver=solver,
            reltol=reltol, abstol=abstol, dtmax=dtmax,
        )
        outputs = extract_outputs_v15c(
            result, p, experiment[1, 4];
            use_tc_measurement_model=use_tc_measurement_model,
        )
        return outputs, result, experiment
    end
end

begin # v15c objective helpers
    const TIMING_WEIGHT_v15c = 0.15
    const SLOPE_WEIGHT_v15c = 0.25

    function normalized_slope_mse_v15c(model, experiment)
        length(model) < 2 && return 0.0
        scale = max(maximum(experiment) - minimum(experiment), 1.0)
        return mean(abs2, diff(model) .- diff(experiment)) / scale^2
    end

    function timing_penalty_v15c(time, model, experiment)
        duration = max(time[end] - time[1], 1.0)
        return ((get_t90_v3(time, model) - get_t90_v3(time, experiment)) / duration)^2
    end

    function signal_loss_v15c(time, model, experiment)
        return normalized_mse_v3(model, experiment) +
               SLOPE_WEIGHT_v15c * normalized_slope_mse_v15c(model, experiment) +
               TIMING_WEIGHT_v15c * timing_penalty_v15c(time, model, experiment)
    end

    function ordering_loss_v15c(model, experiment)
        model_offsets = (
            model[end, 2] - model[end, 1],
            model[end, 3] - model[end, 1],
        )
        experiment_offsets = (
            experiment[end, 2] - experiment[end, 1],
            experiment[end, 3] - experiment[end, 1],
        )
        return mean(abs2, collect(model_offsets) .- collect(experiment_offsets)) /
               ORDERING_SCALE_v15c^2
    end

    function loss_cases_v15c(p, keys; is_cooling=false, nodes=15,
                            use_rosseland_radiation=false,
                            rosseland_profile=:uniform,
                            use_axial_view_radiation=false,
                            axial_view_strength=1.0,
                            use_tc_measurement_model=false)
        total = 0.0
        for simulation_id in keys
            try
                model, result, experiment = solve_case_v15c(
                    p, simulation_id; is_cooling=is_cooling, nodes=nodes,
                    use_rosseland_radiation=use_rosseland_radiation,
                    rosseland_profile=rosseland_profile,
                    use_axial_view_radiation=use_axial_view_radiation,
                    axial_view_strength=axial_view_strength,
                    use_tc_measurement_model=use_tc_measurement_model,
                    reltol=2e-5, abstol=2e-6, dtmax=60.0,
                )
                all(isfinite, model) || return Inf
                total += mean(
                    signal_loss_v15c(result.time, model[:, j], experiment[:, j])
                    for j in 1:5
                )
                if !is_cooling
                    total += ORDERING_WEIGHT_v15c * ordering_loss_v15c(model, experiment)
                end
            catch err
                err isa InterruptException && rethrow()
                return Inf
            end
        end
        return total / length(keys)
    end

    function loss_heating_v15c(p, keys=sim_key_heat; nodes=15)
        regularization = 0.0
        return loss_cases_v15c(p, keys; is_cooling=false, nodes=nodes) + regularization
    end

    function loss_cooling_v15c(p, keys=sim_key_cool; nodes=15)
        regularization = 0.0
        return loss_cases_v15c(p, keys; is_cooling=true, nodes=nodes) + regularization
    end
end

function optimize_parameter_subset_v15c(objective, p_initial, indices, lower, upper;
                                      maximum_iterations=200,
                                      maximum_time_seconds=120.0,
                                      optimizer=NLopt.LN_NELDERMEAD(),
                                      label="subset-v15c")
    base = Float64.(p_initial)
    selected = collect(indices)
    function subset_objective(theta)
        trial = copy(base)
        trial[selected] .= theta
        return objective(trial)
    end
    fit = optimize_with_nlopt_v3(
        subset_objective,
        base[selected],
        Float64.(lower[selected]),
        Float64.(upper[selected]);
        maximum_iterations=maximum_iterations,
        maximum_time_seconds=maximum_time_seconds,
        algorithm=optimizer,
        label=label,
    )
    parameters = copy(base)
    parameters[selected] .= fit.minimizer
    return merge(fit, (parameters=parameters, indices=selected))
end

function calibrate_v15c(; heating_keys=sim_key_heat,
                      nodes=15,
                      maximum_iterations=80,
                      maximum_time_seconds=240.0,
                      optimizer=NLopt.LN_NELDERMEAD())
    objective = p -> loss_heating_v15c(p, heating_keys; nodes=nodes)
    return optimize_parameter_subset_v15c(
        objective, pnew_v15c, fit_heat_transfer_indices_v15c,
        lb_full_v15c, ub_full_v15c;
        maximum_iterations=maximum_iterations,
        maximum_time_seconds=maximum_time_seconds,
        optimizer=optimizer,
        label="v15c-apparent-heat-transfer",
    )
end

function quick_calibration_v15c()
    return calibrate_v15c(maximum_iterations=20, maximum_time_seconds=60.0)
end

begin # v15c post-processing
    function transient_case_v15c(simulation_id, p=pnew_v15c; is_cooling=false,
                               nodes=default_nodes, use_rosseland_radiation=false,
                               rosseland_profile=:uniform,
                               use_axial_view_radiation=false,
                               axial_view_strength=1.0,
                               use_tc_measurement_model=false)
        model, result, experiment = solve_case_v15c(
            p, simulation_id; is_cooling=is_cooling, nodes=nodes,
            use_rosseland_radiation=use_rosseland_radiation,
            rosseland_profile=rosseland_profile,
            use_axial_view_radiation=use_axial_view_radiation,
            axial_view_strength=axial_view_strength,
            use_tc_measurement_model=use_tc_measurement_model,
        )
        return (
            time=result.time,
            T8_model=model[:, 1], T8_experiment=experiment[:, 1],
            T12_wall_model=model[:, 2], T12_wall_experiment=experiment[:, 2],
            T11_wall_model=model[:, 3], T11_wall_experiment=experiment[:, 3],
            T3_model=model[:, 4], T3_experiment=experiment[:, 4],
            T2_model=model[:, 5], T2_experiment=experiment[:, 5],
            T2_cavity=result.cavity_temperature,
            receiver_exit_gas=vec(result.gas_temperature[length(result.z_solid) + 1, :]),
            T3_sample_gas=gas_at_position_v15c(result, T3_SAMPLE_POSITION_v15c),
            rear_tube_inlet=vec(result.rear_tube_temperature[1, :]),
            rear_tube_exit=vec(result.rear_tube_temperature[end, :]),
            receiver_to_tube_heat=result.receiver_to_tube_heat,
            receiver_to_cavity_heat=result.receiver_to_cavity_heat,
            tube_to_cavity_heat=result.tube_to_cavity_heat,
            cavity_ambient_heat_loss=result.cavity_ambient_heat_loss,
            flange_heat_loss=result.flange_heat_loss,
            h_effective=model[:, 6],
            receiver_nusselt=result.receiver_nusselt,
            receiver_effectiveness=result.receiver_effectiveness,
            receiver_ua_over_mcp=result.receiver_ua_over_mcp,
            receiver_gas_heat=result.receiver_gas_heat,
            axial_view_radiation_heat=result.axial_view_radiation_heat,
            full_result=result,
        )
    end

    function plot_case_v15c(simulation_id, params=pnew_v15c; is_cooling=false,
                          nodes=default_nodes, use_rosseland_radiation=false,
                          rosseland_profile=:uniform,
                          use_axial_view_radiation=false,
                          axial_view_strength=1.0,
                          use_tc_measurement_model=false)
        @eval using StatsPlots
        data = transient_case_v15c(simulation_id, params;
                                 is_cooling=is_cooling, nodes=nodes,
                                 use_rosseland_radiation=use_rosseland_radiation,
                                 rosseland_profile=rosseland_profile,
                                 use_axial_view_radiation=use_axial_view_radiation,
                                 axial_view_strength=axial_view_strength,
                                 use_tc_measurement_model=use_tc_measurement_model)
        sensors = (:T8, :T12_wall, :T11_wall, :T3, :T2)
        colors = (:blue, :red, :green, :purple, :black)
        plot_object = plot(title="1D_v15c transient: $simulation_id",
                           xlabel="Time (s)", ylabel="Temperature (K)",
                           legend=:outerright,
                           ylims=TEMPERATURE_AXIS_LIMITS_v15c)
        for (sensor, color) in zip(sensors, colors)
            plot!(plot_object, data.time, getproperty(data, Symbol(sensor, "_model"));
                  label="$(sensor) model", lw=2, color=color)
            scatter!(plot_object, data.time,
                     getproperty(data, Symbol(sensor, "_experiment"));
                     label="$(sensor) experiment", ms=2, markerstrokewidth=0,
                     color=color)
        end
        plot!(plot_object, data.time, data.rear_tube_exit;
              label="rear tube wall exit", color=:orange, linestyle=:dot)
        return plot_object
    end

    function compute_metrics_v15c(p=pnew_v15c; heating_keys=sim_key_heat,
                                cooling_keys=sim_key_cool, nodes=default_nodes,
                                use_rosseland_radiation=false,
                                rosseland_profile=:uniform,
                                use_axial_view_radiation=false,
                                axial_view_strength=1.0,
                                use_tc_measurement_model=false)
        metrics = NamedTuple[]
        for (keys, cooling) in ((heating_keys, false), (cooling_keys, true))
            for simulation_id in keys
                model, result, experiment = solve_case_v15c(
                    p, simulation_id; is_cooling=cooling, nodes=nodes,
                    use_rosseland_radiation=use_rosseland_radiation,
                    rosseland_profile=rosseland_profile,
                    use_axial_view_radiation=use_axial_view_radiation,
                    axial_view_strength=axial_view_strength,
                    use_tc_measurement_model=use_tc_measurement_model,
                )
                for (sensor, column) in ((:T8, 1), (:T12_wall, 2),
                                         (:T11_wall, 3), (:T3, 4), (:T2, 5))
                    push!(metrics, (
                        simulation_id=simulation_id,
                        phase=cooling ? :cooling : :heating,
                        sensor=sensor,
                        rmse_K=sqrt(mean(abs2, model[:, column] .- experiment[:, column])),
                        steady_error_K=model[end, column] - experiment[end, column],
                        t90_error_s=get_t90_v3(result.time, model[:, column]) -
                                    get_t90_v3(result.time, experiment[:, column]),
                        shape_loss=normalized_slope_mse_v15c(
                            model[:, column], experiment[:, column]
                        ),
                    ))
                end
                for (sensor, model_delta, experiment_delta) in (
                    (:T12_minus_T8, model[:, 2] .- model[:, 1],
                     experiment[:, 2] .- experiment[:, 1]),
                    (:T11_minus_T8, model[:, 3] .- model[:, 1],
                     experiment[:, 3] .- experiment[:, 1]),
                )
                    push!(metrics, (
                        simulation_id=simulation_id,
                        phase=cooling ? :cooling : :heating,
                        sensor=sensor,
                        rmse_K=sqrt(mean(abs2, model_delta .- experiment_delta)),
                        steady_error_K=model_delta[end] - experiment_delta[end],
                        t90_error_s=0.0,
                        shape_loss=normalized_slope_mse_v15c(model_delta, experiment_delta),
                    ))
                end
            end
        end
        return metrics
    end
end

println("""
[1D_v15c] Model loaded.
Use run_1D_v15c.jl for the heat-transfer and axial redistribution refit.
The default coefficient vector is stored in pnew_v15c.
""")







