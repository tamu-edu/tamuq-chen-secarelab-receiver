# ============================================================================
# 1D_v12a.jl - ECM-style 1D receiver model with source-distribution diagnostic
# ============================================================================
# v12a keeps the v8b rear-tube/flange, predicted cavity/T2 structure, and the
# v11 receiver gas-solid diagnostics. It is the Stage 2 source-distribution
# test: eta_abs and eta_opt stay fixed, the v11 Stage 1 Nu shape is frozen, and
# only beta_opt plus front_dep are fitted.
#
# Geometry updates supplied for v12a:
#   - Alumina tube length downstream of receiver: 150 mm.
#   - First 46 mm of tube is inside the cavity.
#   - Remaining 104 mm of tube is inside the water-cooled flange.
#   - T3 comparison uses the modeled gas temperature at z = 140 mm.
#   - Alumina adaptor diameter is interpreted as 77.6 mm.
#
# Rosseland remains available as an optional fixed diagnostic, but it is not
# part of the v12a source-distribution fit.
#
# The state vector is [receiver solid cells, rear alumina tube cells, T_cavity].
# ============================================================================

include("1D_v5.jl")

begin # v12a fixed constants
    const EPS_FRONT_FIXED_v12a = EPS_FRONT_FIXED_V5
    const ETA_ABS_FIXED_v12a = ETA_ABS_FIXED_V5
    const BETA_OPT_FIXED_v12a = BETA_OPT_FIXED_V5
    const BETA_OPT_MIN_v12a = 0.10
    const BETA_OPT_MAX_v12a = 300.0
    const H_FRONT_FIXED_v12a = 10.0
    const NU_FD_RECEIVER_v12a = 3.61
    const DEFAULT_ROSSELAND_MEAN_PATH_v12a = 1.0e-3
    const RADIATION_AREA_v12a = A_frt
    const TEMPERATURE_AXIS_LIMITS_v12a = (250.0, 1250.0)
    const ETA_OPT_FIXED_v12a = 1.0
    const ENTRY_SHAPE_EXPONENT_FIXED_v12a = 1.0 / 3.0
    const FRONT_DEPOSITION_INITIAL_v12a = 1.0
    const ORDERING_WEIGHT_v12a = 0.25
    const ORDERING_SCALE_v12a = 100.0

    # Cavity geometry from the annotated v12a schematic.
    const CAVITY_OUTER_RADIUS_v12a = 75.0e-3
    const METAL_THICKNESS_v12a = 18.0e-3
    const INSULATION_OUTER_RADIUS_v12a = CAVITY_OUTER_RADIUS_v12a - METAL_THICKNESS_v12a
    const METAL_OUTER_RADIUS_v12a = CAVITY_OUTER_RADIUS_v12a
    const CAVITY_LENGTH_v12a = 165.0e-3
    const T3_SAMPLE_POSITION_v12a = 140.0e-3

    const ADAPTOR_DIAMETER_v12a = 77.6e-3
    const ADAPTOR_RADIUS_v12a = ADAPTOR_DIAMETER_v12a / 2.0
    const ADAPTOR_LENGTH_v12a = 57.0e-3
    const ADAPTOR_TUBE_RADIUS_v12a = 6.5e-3
    const ADAPTOR_RECEIVER_OVERLAP_LENGTH_v12a = 28.0e-3
    const ADAPTOR_TUBE_OVERLAP_LENGTH_v12a = ADAPTOR_LENGTH_v12a - ADAPTOR_RECEIVER_OVERLAP_LENGTH_v12a
    const ADAPTOR_OVERLAP_LENGTH_v12a = ADAPTOR_RECEIVER_OVERLAP_LENGTH_v12a
    const ADAPTOR_FREE_LENGTH_v12a = ADAPTOR_TUBE_OVERLAP_LENGTH_v12a
    const METAL_LENGTH_v12a = CAVITY_LENGTH_v12a
    const BACKPLATE_THICKNESS_v12a = METAL_THICKNESS_v12a
    const INSULATION_MONOLITH_LENGTH_v12a = L
    const INSULATION_ADAPTOR_LENGTH_v12a = ADAPTOR_LENGTH_v12a
    const INSULATION_CONDUCTIVITY_v12a = 0.08

    const RECEIVER_EQ_RADIUS_v12a = sqrt(A_frt / pi)
    const RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v12a =
        4.0 * pi * INSULATION_CONDUCTIVITY_v12a /
        log(INSULATION_OUTER_RADIUS_v12a / RECEIVER_EQ_RADIUS_v12a)
    const ADAPTOR_TO_CAVITY_RESISTANCE_v12a =
        log(INSULATION_OUTER_RADIUS_v12a / ADAPTOR_RADIUS_v12a) /
        (4.0 * pi * INSULATION_CONDUCTIVITY_v12a * ADAPTOR_LENGTH_v12a)
    const ADAPTOR_CONTACT_RESISTANCE_v12a =
        1.0 / (100.0 * 4.0 * w_t * ADAPTOR_OVERLAP_LENGTH_v12a)

    # Rear extension geometry supplied for v12a.
    const REAR_TUBE_LENGTH_v12a = 150.0e-3
    const REAR_TUBE_CAVITY_LENGTH_v12a = 46.0e-3
    const REAR_TUBE_FLANGE_LENGTH_v12a = REAR_TUBE_LENGTH_v12a - REAR_TUBE_CAVITY_LENGTH_v12a
    const REAR_TUBE_DEFAULT_NODES_v12a = 12
    const REAR_TUBE_GAS_RADIUS_v12a = ADAPTOR_TUBE_RADIUS_v12a
    const REAR_TUBE_WALL_THICKNESS_v12a = 1.5e-3
    const REAR_TUBE_OUTER_RADIUS_v12a = REAR_TUBE_GAS_RADIUS_v12a + REAR_TUBE_WALL_THICKNESS_v12a
    const REAR_TUBE_FLOW_AREA_v12a = pi * REAR_TUBE_GAS_RADIUS_v12a^2
    const REAR_TUBE_HYDRAULIC_DIAMETER_v12a = 2.0 * REAR_TUBE_GAS_RADIUS_v12a
    const REAR_TUBE_EXCHANGE_PERIMETER_v12a = 2.0 * pi * REAR_TUBE_GAS_RADIUS_v12a

    const RECEIVER_TO_TUBE_CONDUCTANCE_v12a =
        1.0 / ADAPTOR_CONTACT_RESISTANCE_v12a

    # Alumina and aluminum properties follow the older 0D material set. The
    # flange is treated as an isothermal water-cooled sink in v12a.
    const ALUMINA_DENSITY_v12a = 3900.0
    const ALUMINUM_DENSITY_v12a = 2700.0
    const ALUMINUM_HEAT_CAPACITY_v12a = 900.0
    const ALUMINUM_CONDUCTIVITY_v12a = 205.0
    const METAL_EMISSIVITY_v12a = 0.2
    const H_NAT_CAVITY_v12a = 10.0
    const INSULATION_DENSITY_v12a = 140.0
    const INSULATION_HEAT_CAPACITY_v12a = 1360.0
    const WATER_FLANGE_TEMPERATURE_v12a = 25.0 + 273.15

    const CAVITY_FELT_VOLUME_v12a = max(
        0.0,
        pi * INSULATION_OUTER_RADIUS_v12a^2 * CAVITY_LENGTH_v12a -
        A_frt * L -
        pi * ADAPTOR_RADIUS_v12a^2 * ADAPTOR_LENGTH_v12a,
    )
    const CAVITY_METAL_VOLUME_v12a =
        pi * (METAL_OUTER_RADIUS_v12a^2 - INSULATION_OUTER_RADIUS_v12a^2) *
        CAVITY_LENGTH_v12a +
        pi * METAL_OUTER_RADIUS_v12a^2 * BACKPLATE_THICKNESS_v12a
    const CAVITY_OUTER_AREA_v12a =
        2.0 * pi * METAL_OUTER_RADIUS_v12a * CAVITY_LENGTH_v12a +
        pi * METAL_OUTER_RADIUS_v12a^2
    const CAVITY_HEAT_CAPACITY_v12a =
        INSULATION_DENSITY_v12a * INSULATION_HEAT_CAPACITY_v12a *
        CAVITY_FELT_VOLUME_v12a +
        ALUMINUM_DENSITY_v12a * ALUMINUM_HEAT_CAPACITY_v12a *
        CAVITY_METAL_VOLUME_v12a
end

# Final v12a Stage 2 parameter vector p[1:6]
#   p[1]  A_Nu          frozen v11 entry/development Nu prefactor
#   p[2]  B_Re          frozen v11 Reynolds exponent
#   p[3]  C_Pr          frozen v11 Prandtl exponent
#   p[4]  L_entry_m     frozen v11 axial decay length [m]
#   p[5]  beta_opt      fitted optical attenuation [1/m]
#   p[6]  front_dep     fitted front-cell deposition fraction
begin # v12a parameter values and bounds
    pnew_v12a = [
        1.099639225439012,
        1.0,
        0.31833744247285517,
        0.03790557942990763,
        BETA_OPT_FIXED_v12a,
        FRONT_DEPOSITION_INITIAL_v12a,
    ]
    lb_full_v12a = [
        pnew_v12a[1],
        pnew_v12a[2],
        pnew_v12a[3],
        pnew_v12a[4],
        BETA_OPT_MIN_v12a,
        0.0,
    ]
    ub_full_v12a = [
        pnew_v12a[1],
        pnew_v12a[2],
        pnew_v12a[3],
        pnew_v12a[4],
        BETA_OPT_MAX_v12a,
        1.0,
    ]

    fit_heat_transfer_indices_v12a = Int[]
    fit_source_indices_v12a = [5, 6]
    fit_radiation_indices_v12a = Int[]
end

irradiance_factor_v12a(_irradiance, _p) = ETA_OPT_FIXED_v12a

function receiver_nusselt_v12a(z, reynolds, prandtl, p)
    z_eff = max(z, 0.5e-3)
    entry_decay = exp(-z_eff / max(p[4], 1e-6))
    entry_shape = (Dh / z_eff)^ENTRY_SHAPE_EXPONENT_FIXED_v12a * entry_decay
    enhancement = p[1] * max(reynolds, eps(Float64))^p[2] *
                  max(prandtl, eps(Float64))^p[3] * entry_shape
    return NU_FD_RECEIVER_v12a + max(0.0, enhancement)
end

smoothstep_v12a(x) = begin
    xi = clamp(x, 0.0, 1.0)
    xi^2 * (3.0 - 2.0 * xi)
end

function rosseland_mean_path_v12a(z, p; profile=:uniform)
    ell_uniform = DEFAULT_ROSSELAND_MEAN_PATH_v12a
    xi = smoothstep_v12a(z / L)
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

rosseland_beta_v12a(p; z=0.0, profile=:uniform) =
    1.0 / rosseland_mean_path_v12a(z, p; profile=profile)

rosseland_conductivity_v12a(T, p; z=0.0, profile=:uniform) =
    16.0 * sigma_sb * property_temperature(T)^3 /
    (3.0 * rosseland_beta_v12a(p; z=z, profile=profile))

rosseland_optical_thickness_v12a(length_scale, p) =
    rosseland_beta_v12a(p) * length_scale

alumina_conductivity_v12a(T) =
    5.5 + 34.5 * exp(-0.0033 * (property_temperature(T) - 273.15))

alumina_heat_capacity_v12a(T) =
    max(500.0, (1.00446 + 1.742e-4 * property_temperature(T) -
                2.796e4 * property_temperature(T)^-2) * 1000.0)

function rear_tube_outer_radius_v12a(z_rear)
    return z_rear <= ADAPTOR_TUBE_OVERLAP_LENGTH_v12a ?
           ADAPTOR_RADIUS_v12a : REAR_TUBE_OUTER_RADIUS_v12a
end

function rear_tube_solid_area_v12a(z_rear)
    outer_radius = rear_tube_outer_radius_v12a(z_rear)
    return pi * (outer_radius^2 - REAR_TUBE_GAS_RADIUS_v12a^2)
end

function rear_tube_cavity_conductance_per_length_v12a(z_rear, p)
    z_rear > REAR_TUBE_CAVITY_LENGTH_v12a && return 0.0
    outer_radius = min(rear_tube_outer_radius_v12a(z_rear), 0.98 * INSULATION_OUTER_RADIUS_v12a)
    return 2.0 * pi * INSULATION_CONDUCTIVITY_v12a /
           log(INSULATION_OUTER_RADIUS_v12a / outer_radius)
end

function rear_tube_flange_conductance_per_length_v12a(z_rear, Ttube)
    z_rear <= REAR_TUBE_CAVITY_LENGTH_v12a && return 0.0
    tube_wall_resistance =
        log(REAR_TUBE_OUTER_RADIUS_v12a / REAR_TUBE_GAS_RADIUS_v12a) /
        (2.0 * pi * alumina_conductivity_v12a(Ttube))
    aluminum_resistance =
        log(METAL_OUTER_RADIUS_v12a / REAR_TUBE_OUTER_RADIUS_v12a) /
        (2.0 * pi * ALUMINUM_CONDUCTIVITY_v12a)
    return 1.0 / (tube_wall_resistance + aluminum_resistance)
end

function cavity_loss_to_ambient_v12a(Tcavity, ambient)
    return H_NAT_CAVITY_v12a * CAVITY_OUTER_AREA_v12a * (Tcavity - ambient) +
           METAL_EMISSIVITY_v12a * sigma_sb * CAVITY_OUTER_AREA_v12a *
           (Tcavity^4 - ambient^4)
end

function rear_tube_capacity_v12a(Ttube, z_rear, dx_rear)
    return ALUMINA_DENSITY_v12a * alumina_heat_capacity_v12a(Ttube) *
           rear_tube_solid_area_v12a(z_rear) * dx_rear
end

function rear_tube_htc_v12a(Twall, Tgas, flow_lpm, Tin)
    flow_lpm <= 1e-12 && return 0.0
    mdot = m_dot(flow_lpm, Tin)
    Tfilm = 0.5 * (Twall + Tgas)
    cp = cpf_f(Tfilm)
    mu = mu_f_f(Tfilm)
    kf = kf_f(Tfilm)
    reynolds = mdot * REAR_TUBE_HYDRAULIC_DIAMETER_v12a /
               (REAR_TUBE_FLOW_AREA_v12a * mu)
    prandtl = cp * mu / kf
    nusselt = reynolds < 2300.0 ?
              3.66 :
              0.023 * reynolds^0.8 * prandtl^0.4
    return nusselt * kf / REAR_TUBE_HYDRAULIC_DIAMETER_v12a
end

function gas_profile_v12a!(Tg, Qgas, hcell, Qgas_rear, hrear,
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

    mdot = m_dot(flow, Tin)
    for i in eachindex(Ts)
        Tfilm = 0.5 * (Ts[i] + Tg[i])
        cp = cpf_f(Tfilm)
        mu = mu_f_f(Tfilm)
        kf = kf_f(Tfilm)
        reynolds = mdot * Dh / (A_flow * mu)
        prandtl = cp * mu / kf
        nusselt = receiver_nusselt_v12a(z[i], reynolds, prandtl, p)
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
        hrear[j] = rear_tube_htc_v12a(Ttube[j], Tg[inlet_index], flow, Tin)
        UA = hrear[j] * REAR_TUBE_EXCHANGE_PERIMETER_v12a * dx_rear
        effectiveness = -expm1(-UA / (mdot * cp))
        Tg[inlet_index + 1] = Tg[inlet_index] + effectiveness *
                              (Ttube[j] - Tg[inlet_index])
        Qgas_rear[j] = mdot * cp * (Tg[inlet_index + 1] - Tg[inlet_index])
    end
    return nothing
end

function receiver_rhs_v12a!(du, u, time, p, operating, z, dx, solar_weights,
                          z_rear, dx_rear,
                          Tg, Qgas, hcell, Qgas_rear, hrear,
                          use_rosseland_radiation, rosseland_profile)
    nodes = length(z)
    rear_nodes = length(z_rear)
    Ts = view(u, 1:nodes)
    Ttube = view(u, (nodes + 1):(nodes + rear_nodes))
    Tcavity = u[nodes + rear_nodes + 1]
    dTs = view(du, 1:nodes)
    dTtube = view(du, (nodes + 1):(nodes + rear_nodes))
    ambient = operating.ambient_temperature(time)

    gas_profile_v12a!(Tg, Qgas, hcell, Qgas_rear, hrear,
                    Ts, Ttube, time, p, operating, z, dx, z_rear, dx_rear)
    fill!(du, 0.0)

    for i in 1:(nodes - 1)
        ki = ks_f(Ts[i])
        kj = ks_f(Ts[i + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        Qcond = kface * A_solid * (Ts[i] - Ts[i + 1]) / dx
        dTs[i] -= Qcond
        dTs[i + 1] += Qcond
        if use_rosseland_radiation
            Tface = 0.5 * (Ts[i] + Ts[i + 1])
            zface = 0.5 * (z[i] + z[i + 1])
            Qrad = rosseland_conductivity_v12a(
                Tface, p; z=zface, profile=rosseland_profile,
            ) * RADIATION_AREA_v12a *
                   (Ts[i] - Ts[i + 1]) / dx
            dTs[i] -= Qrad
            dTs[i + 1] += Qrad
        end
    end

    irradiance = max(0.0, operating.irradiance(time))
    Qsolar = ETA_ABS_FIXED_v12a * irradiance_factor_v12a(irradiance, p) *
             irradiance * A_frt
    radial_conductance = RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v12a
    for i in 1:nodes
        Qside_cavity = radial_conductance * dx * (Ts[i] - Tcavity)
        dTs[i] += Qsolar * solar_weights[i] - Qgas[i] - Qside_cavity
        du[nodes + rear_nodes + 1] += Qside_cavity
    end

    Qfront = H_FRONT_FIXED_v12a * A_frt * (Ts[1] - ambient) +
             EPS_FRONT_FIXED_v12a * sigma_sb * A_frt * (Ts[1]^4 - ambient^4)
    Qreceiver_tube = RECEIVER_TO_TUBE_CONDUCTANCE_v12a * (Ts[end] - Ttube[1])

    dTs[1] -= Qfront
    dTs[end] -= Qreceiver_tube
    dTtube[1] += Qreceiver_tube

    for j in 1:(rear_nodes - 1)
        ki = alumina_conductivity_v12a(Ttube[j])
        kj = alumina_conductivity_v12a(Ttube[j + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        Aface = 0.5 * (rear_tube_solid_area_v12a(z_rear[j]) +
                       rear_tube_solid_area_v12a(z_rear[j + 1]))
        Qcond = kface * Aface * (Ttube[j] - Ttube[j + 1]) / dx_rear
        dTtube[j] -= Qcond
        dTtube[j + 1] += Qcond
    end

    for j in 1:rear_nodes
        G_cavity = rear_tube_cavity_conductance_per_length_v12a(z_rear[j], p)
        G_flange = rear_tube_flange_conductance_per_length_v12a(z_rear[j], Ttube[j])
        Q_cavity = G_cavity * dx_rear * (Ttube[j] - Tcavity)
        Q_flange = G_flange * dx_rear * (Ttube[j] - WATER_FLANGE_TEMPERATURE_v12a)
        dTtube[j] -= Qgas_rear[j] + Q_cavity + Q_flange
        du[nodes + rear_nodes + 1] += Q_cavity
    end

    Qcavity_ambient = cavity_loss_to_ambient_v12a(Tcavity, ambient)
    du[nodes + rear_nodes + 1] -= Qcavity_ambient

    cell_volume = A_solid * dx
    for i in 1:nodes
        capacity = rho_s * Cps_f(Ts[i]) * cell_volume
        dTs[i] /= capacity
    end
    for j in 1:rear_nodes
        dTtube[j] /= max(rear_tube_capacity_v12a(Ttube[j], z_rear[j], dx_rear),
                         eps(Float64))
    end
    du[nodes + rear_nodes + 1] /= max(CAVITY_HEAT_CAPACITY_v12a, eps(Float64))
    return nothing
end

function receiver_ode_v12a!(dTs, Ts, context, time)
    receiver_rhs_v12a!(
        dTs, Ts, time, context.p, context.operating,
        context.z, context.dx, context.solar_weights,
        context.z_rear, context.dx_rear,
        context.Tg, context.Qgas, context.hcell, context.Qgas_rear, context.hrear,
        context.use_rosseland_radiation, context.rosseland_profile,
    )
    return nothing
end

function simulate_v12a(p, operating, save_times;
                     initial_temperature=Tamb, nodes=default_nodes,
                     rear_nodes=REAR_TUBE_DEFAULT_NODES_v12a,
                     initial_rear_temperature=nothing,
                     initial_cavity_temperature=Tamb,
                     use_rosseland_radiation=false,
                     rosseland_profile=:uniform,
                     solver=Rodas5P(autodiff=AutoFiniteDiff()),
                     reltol=1e-6, abstol=1e-7, dtmax=30.0)
    length(p) == 6 || throw(ArgumentError("1D_v12a expects 6 parameters"))
    all(isfinite, p) || throw(ArgumentError("all model parameters must be finite"))
    all(>(0.0), p[1:5]) && 0.0 <= p[6] <= 1.0 ||
        throw(ArgumentError("Nu/entry coefficients and beta_opt must be positive; front_dep must be in [0, 1]"))
    nodes >= 3 || throw(ArgumentError("at least three axial nodes are required"))
    rear_nodes >= 2 || throw(ArgumentError("at least two rear tube nodes are required"))

    time = Float64.(save_times)
    isempty(time) && throw(ArgumentError("save_times cannot be empty"))
    (length(time) == 1 || all(diff(time) .> 0.0)) ||
        throw(ArgumentError("save_times must be strictly increasing"))

    dx = L / nodes
    dx_rear = REAR_TUBE_LENGTH_v12a / rear_nodes
    z = collect(range(dx / 2, step=dx, length=nodes))
    z_rear = collect(range(dx_rear / 2, step=dx_rear, length=rear_nodes))
    zg = vcat(
        collect(range(0.0, L, length=nodes + 1)),
        L .+ collect(range(dx_rear, step=dx_rear, length=rear_nodes)),
    )
    weights = solar_weights_v5(
        p[5], p[6], nodes,
    )
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
    )
    problem = ODEProblem(receiver_ode_v12a!, u_initial, (time[1], time[end]), context)
    solution = solve(
        problem, solver; saveat=time, save_everystep=false, dense=false,
        reltol=reltol, abstol=abstol, dtmax=dtmax, tstops=time,
    )
    successful_retcode(solution) ||
        error("1D_v12a ODE solve failed with retcode $(solution.retcode)")

    state_history = reduce(hcat, solution.u)
    Ts_history = state_history[1:nodes, :]
    rear_temperature_history .= state_history[(nodes + 1):(nodes + rear_nodes), :]
    cavity_temperature_history .= state_history[nodes + rear_nodes + 1, :]
    for output_index in eachindex(time)
        t = time[output_index]
        Ts_now = view(Ts_history, :, output_index)
        Ttube_now = view(rear_temperature_history, :, output_index)
        Tcavity_now = cavity_temperature_history[output_index]
        gas_profile_v12a!(Tg, Qgas, hcell, Qgas_rear, hrear, Ts_now, Ttube_now,
                        t, p, operating, z, dx, z_rear, dx_rear)
        Tg_history[:, output_index] .= Tg
        h_history[:, output_index] .= hcell
        qgas_history[:, output_index] .= Qgas
        flow_now = max(0.0, operating.flow_lpm(t))
        Tin_now = operating.inlet_temperature(t)
        if flow_now <= 1e-12
            fill!(view(nu_history, :, output_index), 0.0)
            fill!(view(effectiveness_history, :, output_index), 0.0)
            fill!(view(ua_over_mcp_history, :, output_index), 0.0)
        else
            mdot_now = m_dot(flow_now, Tin_now)
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
            RECEIVER_TO_TUBE_CONDUCTANCE_v12a * (Ts_now[end] - Ttube_now[1])
        receiver_to_cavity_heat_history[output_index] = sum(
            RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v12a * dx *
            (Ts_now[i] - Tcavity_now)
            for i in 1:nodes
        )
        tube_to_cavity_heat_history[output_index] = sum(
            rear_tube_cavity_conductance_per_length_v12a(z_rear[j], p) * dx_rear *
            (Ttube_now[j] - Tcavity_now)
            for j in 1:rear_nodes
        )
        cavity_ambient_heat_history[output_index] =
            cavity_loss_to_ambient_v12a(Tcavity_now, operating.ambient_temperature(t))
        flange_heat_history[output_index] = sum(
            rear_tube_flange_conductance_per_length_v12a(z_rear[j], Ttube_now[j]) *
            dx_rear * (Ttube_now[j] - WATER_FLANGE_TEMPERATURE_v12a)
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
        rear_tube_heat_transfer_coefficient=h_rear_history,
        receiver_to_tube_heat=receiver_to_tube_heat_history,
        receiver_to_cavity_heat=receiver_to_cavity_heat_history,
        tube_to_cavity_heat=tube_to_cavity_heat_history,
        cavity_ambient_heat_loss=cavity_ambient_heat_history,
        flange_heat_loss=flange_heat_history,
        ode_solution=solution,
    )
end

begin # v12a experimental interface
    paired_observation_v12a(data, simulation_id, obs_a, obs_b) =
        0.5 .* (observation(data, simulation_id, obs_a) .+
                observation(data, simulation_id, obs_b))

    function measured_initial_profile_v12a(data, simulation_id)
        z8 = sensor_positions[:T8]
        z9 = sensor_positions[:T9]
        z10 = sensor_positions[:T10]
        T8 = observation(data, simulation_id, "_T8")[1]
        T9 = paired_observation_v12a(data, simulation_id, "_T9", "_T12")[1]
        T10 = paired_observation_v12a(data, simulation_id, "_T10", "_T11")[1]
        return function (z)
            z <= z8 && return T8
            z >= z10 && return T10
            if z <= z9
                return T8 + (T9 - T8) * (z - z8) / (z9 - z8)
            end
            return T9 + (T10 - T9) * (z - z9) / (z10 - z9)
        end
    end

    function measured_initial_rear_profile_v12a(data, simulation_id)
        receiver_profile = measured_initial_profile_v12a(data, simulation_id)
        T_receiver_exit = receiver_profile(L)
        T_tube_exit = observation(data, simulation_id, "_Tf")[1]
        return function (z_rear)
            fraction = clamp(z_rear / REAR_TUBE_LENGTH_v12a, 0.0, 1.0)
            return T_receiver_exit + fraction * (T_tube_exit - T_receiver_exit)
        end
    end

    function experimental_case_v12a(simulation_id; is_cooling=false)
        data = is_cooling ? measurements_cooling : measurements
        conditions = is_cooling ? simulation_conditions_cooling : simulation_conditions
        time = observation_time(data, simulation_id)
        experiment = hcat(
            observation(data, simulation_id, "_T8"),
            paired_observation_v12a(data, simulation_id, "_T9", "_T12"),
            paired_observation_v12a(data, simulation_id, "_T10", "_T11"),
            observation(data, simulation_id, "_Tf"),
            observation(data, simulation_id, "_T2"),
        )
        return conditions[simulation_id], time, experiment, data
    end

    function gas_at_position_v12a(result, position)
        z = result.z_gas
        values = result.gas_temperature
        position <= z[1] && return vec(values[1, :])
        position >= z[end] && return vec(values[end, :])
        i = searchsortedlast(z, position)
        fraction = (position - z[i]) / (z[i + 1] - z[i])
        return vec((1.0 - fraction) .* values[i, :] .+ fraction .* values[i + 1, :])
    end

    function extract_outputs_v12a(result, p, T3_initial)
        T8 = solid_at_v3(result, sensor_positions[:T8])
        T9_pair = solid_at_v3(result, sensor_positions[:T9])
        T10_pair = solid_at_v3(result, sensor_positions[:T10])
        Tf = gas_at_position_v12a(result, T3_SAMPLE_POSITION_v12a)
        T2 = result.cavity_temperature
        h_mean = vec(mean(result.heat_transfer_coefficient, dims=1))
        return hcat(T8, T9_pair, T10_pair, Tf, T2, h_mean)
    end

    function solve_case_v12a(p, simulation_id; is_cooling=false,
                           nodes=default_nodes,
                           rear_nodes=REAR_TUBE_DEFAULT_NODES_v12a,
                           use_rosseland_radiation=false,
                           rosseland_profile=:uniform,
                           solver=Rodas5P(autodiff=AutoFiniteDiff()),
                           reltol=1e-6, abstol=1e-7, dtmax=30.0)
        conditions, time, experiment, data = experimental_case_v12a(
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
        result = simulate_v12a(
            p, operating, time;
            initial_temperature=measured_initial_profile_v12a(data, simulation_id),
            initial_rear_temperature=measured_initial_rear_profile_v12a(data, simulation_id),
            initial_cavity_temperature=T2[1],
            use_rosseland_radiation=use_rosseland_radiation,
            rosseland_profile=rosseland_profile,
            nodes=nodes, rear_nodes=rear_nodes, solver=solver,
            reltol=reltol, abstol=abstol, dtmax=dtmax,
        )
        outputs = extract_outputs_v12a(result, p, experiment[1, 4])
        return outputs, result, experiment
    end
end

begin # v12a objective helpers
    const TIMING_WEIGHT_v12a = 0.15
    const SLOPE_WEIGHT_v12a = 0.25

    function normalized_slope_mse_v12a(model, experiment)
        length(model) < 2 && return 0.0
        scale = max(maximum(experiment) - minimum(experiment), 1.0)
        return mean(abs2, diff(model) .- diff(experiment)) / scale^2
    end

    function timing_penalty_v12a(time, model, experiment)
        duration = max(time[end] - time[1], 1.0)
        return ((get_t90_v3(time, model) - get_t90_v3(time, experiment)) / duration)^2
    end

    function signal_loss_v12a(time, model, experiment)
        return normalized_mse_v3(model, experiment) +
               SLOPE_WEIGHT_v12a * normalized_slope_mse_v12a(model, experiment) +
               TIMING_WEIGHT_v12a * timing_penalty_v12a(time, model, experiment)
    end

    function ordering_loss_v12a(model, experiment)
        model_offsets = (
            model[end, 2] - model[end, 1],
            model[end, 3] - model[end, 1],
        )
        experiment_offsets = (
            experiment[end, 2] - experiment[end, 1],
            experiment[end, 3] - experiment[end, 1],
        )
        return mean(abs2, collect(model_offsets) .- collect(experiment_offsets)) /
               ORDERING_SCALE_v12a^2
    end

    function loss_cases_v12a(p, keys; is_cooling=false, nodes=15,
                            use_rosseland_radiation=false,
                            rosseland_profile=:uniform)
        total = 0.0
        for simulation_id in keys
            try
                model, result, experiment = solve_case_v12a(
                    p, simulation_id; is_cooling=is_cooling, nodes=nodes,
                    use_rosseland_radiation=use_rosseland_radiation,
                    rosseland_profile=rosseland_profile,
                    reltol=2e-5, abstol=2e-6, dtmax=60.0,
                )
                all(isfinite, model) || return Inf
                total += mean(
                    signal_loss_v12a(result.time, model[:, j], experiment[:, j])
                    for j in 1:5
                )
                if !is_cooling
                    total += ORDERING_WEIGHT_v12a * ordering_loss_v12a(model, experiment)
                end
            catch err
                err isa InterruptException && rethrow()
                return Inf
            end
        end
        return total / length(keys)
    end

    function loss_heating_v12a(p, keys=sim_key_heat; nodes=15)
        regularization = 0.0
        return loss_cases_v12a(p, keys; is_cooling=false, nodes=nodes) + regularization
    end

    function loss_cooling_v12a(p, keys=sim_key_cool; nodes=15)
        regularization = 0.0
        return loss_cases_v12a(p, keys; is_cooling=true, nodes=nodes) + regularization
    end
end

function optimize_parameter_subset_v12a(objective, p_initial, indices, lower, upper;
                                      maximum_iterations=200,
                                      maximum_time_seconds=120.0,
                                      optimizer=NLopt.LN_NELDERMEAD(),
                                      label="subset-v12a")
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

function calibrate_v12a(; heating_keys=sim_key_heat,
                      nodes=15,
                      maximum_iterations=80,
                      maximum_time_seconds=240.0,
                      optimizer=NLopt.LN_NELDERMEAD())
    objective = p -> loss_heating_v12a(p, heating_keys; nodes=nodes)
    return optimize_parameter_subset_v12a(
        objective, pnew_v12a, fit_source_indices_v12a,
        lb_full_v12a, ub_full_v12a;
        maximum_iterations=maximum_iterations,
        maximum_time_seconds=maximum_time_seconds,
        optimizer=optimizer,
        label="v12a-source-stage2",
    )
end

function quick_calibration_v12a()
    return calibrate_v12a(maximum_iterations=20, maximum_time_seconds=60.0)
end

begin # v12a post-processing
    function transient_case_v12a(simulation_id, p=pnew_v12a; is_cooling=false,
                               nodes=default_nodes, use_rosseland_radiation=false,
                               rosseland_profile=:uniform)
        model, result, experiment = solve_case_v12a(
            p, simulation_id; is_cooling=is_cooling, nodes=nodes,
            use_rosseland_radiation=use_rosseland_radiation,
            rosseland_profile=rosseland_profile,
        )
        return (
            time=result.time,
            T8_model=model[:, 1], T8_experiment=experiment[:, 1],
            T9_pair_model=model[:, 2], T9_pair_experiment=experiment[:, 2],
            T10_pair_model=model[:, 3], T10_pair_experiment=experiment[:, 3],
            T3_model=model[:, 4], T3_experiment=experiment[:, 4],
            T2_model=model[:, 5], T2_experiment=experiment[:, 5],
            T2_cavity=result.cavity_temperature,
            receiver_exit_gas=vec(result.gas_temperature[length(result.z_solid) + 1, :]),
            T3_sample_gas=gas_at_position_v12a(result, T3_SAMPLE_POSITION_v12a),
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
            full_result=result,
        )
    end

    function plot_case_v12a(simulation_id, params=pnew_v12a; is_cooling=false,
                          nodes=default_nodes, use_rosseland_radiation=false,
                          rosseland_profile=:uniform)
        @eval using StatsPlots
        data = transient_case_v12a(simulation_id, params;
                                 is_cooling=is_cooling, nodes=nodes,
                                 use_rosseland_radiation=use_rosseland_radiation,
                                 rosseland_profile=rosseland_profile)
        sensors = (:T8, :T9_pair, :T10_pair, :T3, :T2)
        colors = (:blue, :red, :green, :purple, :black)
        plot_object = plot(title="1D_v12a transient: $simulation_id",
                           xlabel="Time (s)", ylabel="Temperature (K)",
                           legend=:outerright,
                           ylims=TEMPERATURE_AXIS_LIMITS_v12a)
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

    function compute_metrics_v12a(p=pnew_v12a; heating_keys=sim_key_heat,
                                cooling_keys=sim_key_cool, nodes=default_nodes,
                                use_rosseland_radiation=false,
                                rosseland_profile=:uniform)
        metrics = NamedTuple[]
        for (keys, cooling) in ((heating_keys, false), (cooling_keys, true))
            for simulation_id in keys
                model, result, experiment = solve_case_v12a(
                    p, simulation_id; is_cooling=cooling, nodes=nodes,
                    use_rosseland_radiation=use_rosseland_radiation,
                    rosseland_profile=rosseland_profile,
                )
                for (sensor, column) in ((:T8, 1), (:T9_pair, 2),
                                         (:T10_pair, 3), (:T3, 4), (:T2, 5))
                    push!(metrics, (
                        simulation_id=simulation_id,
                        phase=cooling ? :cooling : :heating,
                        sensor=sensor,
                        rmse_K=sqrt(mean(abs2, model[:, column] .- experiment[:, column])),
                        steady_error_K=model[end, column] - experiment[end, column],
                        t90_error_s=get_t90_v3(result.time, model[:, column]) -
                                    get_t90_v3(result.time, experiment[:, column]),
                        shape_loss=normalized_slope_mse_v12a(
                            model[:, column], experiment[:, column]
                        ),
                    ))
                end
                for (sensor, model_delta, experiment_delta) in (
                    (:T9_pair_minus_T8, model[:, 2] .- model[:, 1],
                     experiment[:, 2] .- experiment[:, 1]),
                    (:T10_pair_minus_T8, model[:, 3] .- model[:, 1],
                     experiment[:, 3] .- experiment[:, 1]),
                )
                    push!(metrics, (
                        simulation_id=simulation_id,
                        phase=cooling ? :cooling : :heating,
                        sensor=sensor,
                        rmse_K=sqrt(mean(abs2, model_delta .- experiment_delta)),
                        steady_error_K=model_delta[end] - experiment_delta[end],
                        t90_error_s=0.0,
                        shape_loss=normalized_slope_mse_v12a(model_delta, experiment_delta),
                    ))
                end
            end
        end
        return metrics
    end
end

println("""
[1D_v12a] Model loaded.
Use run_1D_v12a.jl for the Stage 2 source-distribution diagnostic.
The default coefficient vector is stored in pnew_v12a.
""")

