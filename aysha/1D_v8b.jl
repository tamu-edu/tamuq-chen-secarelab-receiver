# ============================================================================
# 1D_v8b.jl - 1D receiver model with predicted lumped cavity/T2 state
# ============================================================================
# v8b keeps the v8a rear-tube extension but replaces the measured T2 boundary
# with one lumped cavity thermal state. This is intentionally not a radial mesh:
# the cavity/felt/metal casing is treated as one effective thermal mass whose
# temperature is compared against the measured T2 channel.
#
# Geometry updates supplied for v8b:
#   - Alumina tube length downstream of receiver: 150 mm.
#   - First 46 mm of tube is inside the cavity.
#   - Remaining 104 mm of tube is inside the water-cooled flange.
#   - T3 comparison uses the modeled gas temperature at z = 140 mm.
#   - Alumina adaptor diameter is interpreted as 77.6 mm.
#
# Calibration stages:
#   1. Heating fits gas heat-transfer shape and three irradiance-level factors.
#   2. Cooling fits thermophysical / conduction scales and T3 lag.
#   3. Heating refits gas heat-transfer shape and irradiance-level factors.
#
# The state vector is [receiver solid cells, rear alumina tube cells, T_cavity].
# ============================================================================

include("1D_v5.jl")

begin # v8b fixed constants
    const B_RE_FIXED_v8b = B_RE_FIXED_V5
    const C_PR_FIXED_v8b = C_PR_FIXED_V5
    const EPS_FRONT_FIXED_v8b = EPS_FRONT_FIXED_V5
    const ETA_ABS_FIXED_v8b = ETA_ABS_FIXED_V5
    const BETA_OPT_FIXED_v8b = BETA_OPT_FIXED_V5
    const FRONT_DEPOSITION_FIXED_v8b = FRONT_DEPOSITION_FIXED_V5
    const H_FRONT_FIXED_v8b = 10.0

    # Cavity geometry from the annotated v8b schematic.
    const CAVITY_OUTER_RADIUS_v8b = 75.0e-3
    const METAL_THICKNESS_v8b = 18.0e-3
    const INSULATION_OUTER_RADIUS_v8b = CAVITY_OUTER_RADIUS_v8b - METAL_THICKNESS_v8b
    const METAL_OUTER_RADIUS_v8b = CAVITY_OUTER_RADIUS_v8b
    const CAVITY_LENGTH_v8b = 165.0e-3
    const T3_SAMPLE_POSITION_v8b = 140.0e-3

    const ADAPTOR_DIAMETER_v8b = 77.6e-3
    const ADAPTOR_RADIUS_v8b = ADAPTOR_DIAMETER_v8b / 2.0
    const ADAPTOR_LENGTH_v8b = 57.0e-3
    const ADAPTOR_TUBE_RADIUS_v8b = 6.5e-3
    const ADAPTOR_RECEIVER_OVERLAP_LENGTH_v8b = 28.0e-3
    const ADAPTOR_TUBE_OVERLAP_LENGTH_v8b = ADAPTOR_LENGTH_v8b - ADAPTOR_RECEIVER_OVERLAP_LENGTH_v8b
    const ADAPTOR_OVERLAP_LENGTH_v8b = ADAPTOR_RECEIVER_OVERLAP_LENGTH_v8b
    const ADAPTOR_FREE_LENGTH_v8b = ADAPTOR_TUBE_OVERLAP_LENGTH_v8b
    const METAL_LENGTH_v8b = CAVITY_LENGTH_v8b
    const BACKPLATE_THICKNESS_v8b = METAL_THICKNESS_v8b
    const INSULATION_MONOLITH_LENGTH_v8b = L
    const INSULATION_ADAPTOR_LENGTH_v8b = ADAPTOR_LENGTH_v8b
    const INSULATION_CONDUCTIVITY_v8b = 0.08

    const RECEIVER_EQ_RADIUS_v8b = sqrt(A_frt / pi)
    const RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v8b =
        4.0 * pi * INSULATION_CONDUCTIVITY_v8b /
        log(INSULATION_OUTER_RADIUS_v8b / RECEIVER_EQ_RADIUS_v8b)
    const ADAPTOR_TO_CAVITY_RESISTANCE_v8b =
        log(INSULATION_OUTER_RADIUS_v8b / ADAPTOR_RADIUS_v8b) /
        (4.0 * pi * INSULATION_CONDUCTIVITY_v8b * ADAPTOR_LENGTH_v8b)
    const ADAPTOR_CONTACT_RESISTANCE_v8b =
        1.0 / (100.0 * 4.0 * w_t * ADAPTOR_OVERLAP_LENGTH_v8b)

    # Rear extension geometry supplied for v8b.
    const REAR_TUBE_LENGTH_v8b = 150.0e-3
    const REAR_TUBE_CAVITY_LENGTH_v8b = 46.0e-3
    const REAR_TUBE_FLANGE_LENGTH_v8b = REAR_TUBE_LENGTH_v8b - REAR_TUBE_CAVITY_LENGTH_v8b
    const REAR_TUBE_DEFAULT_NODES_v8b = 12
    const REAR_TUBE_GAS_RADIUS_v8b = ADAPTOR_TUBE_RADIUS_v8b
    const REAR_TUBE_WALL_THICKNESS_v8b = 1.5e-3
    const REAR_TUBE_OUTER_RADIUS_v8b = REAR_TUBE_GAS_RADIUS_v8b + REAR_TUBE_WALL_THICKNESS_v8b
    const REAR_TUBE_FLOW_AREA_v8b = pi * REAR_TUBE_GAS_RADIUS_v8b^2
    const REAR_TUBE_HYDRAULIC_DIAMETER_v8b = 2.0 * REAR_TUBE_GAS_RADIUS_v8b
    const REAR_TUBE_EXCHANGE_PERIMETER_v8b = 2.0 * pi * REAR_TUBE_GAS_RADIUS_v8b

    const RECEIVER_TO_TUBE_CONDUCTANCE_v8b =
        1.0 / ADAPTOR_CONTACT_RESISTANCE_v8b

    # Alumina and aluminum properties follow the older 0D material set. The
    # flange is treated as an isothermal water-cooled sink in v8ba.
    const ALUMINA_DENSITY_v8b = 3900.0
    const ALUMINUM_DENSITY_v8b = 2700.0
    const ALUMINUM_HEAT_CAPACITY_v8b = 900.0
    const ALUMINUM_CONDUCTIVITY_v8b = 205.0
    const METAL_EMISSIVITY_v8b = 0.2
    const H_NAT_CAVITY_v8b = 10.0
    const INSULATION_DENSITY_v8b = 140.0
    const INSULATION_HEAT_CAPACITY_v8b = 1360.0
    const WATER_FLANGE_TEMPERATURE_v8b = 25.0 + 273.15

    const CAVITY_FELT_VOLUME_v8b = max(
        0.0,
        pi * INSULATION_OUTER_RADIUS_v8b^2 * CAVITY_LENGTH_v8b -
        A_frt * L -
        pi * ADAPTOR_RADIUS_v8b^2 * ADAPTOR_LENGTH_v8b,
    )
    const CAVITY_METAL_VOLUME_v8b =
        pi * (METAL_OUTER_RADIUS_v8b^2 - INSULATION_OUTER_RADIUS_v8b^2) *
        CAVITY_LENGTH_v8b +
        pi * METAL_OUTER_RADIUS_v8b^2 * BACKPLATE_THICKNESS_v8b
    const CAVITY_OUTER_AREA_v8b =
        2.0 * pi * METAL_OUTER_RADIUS_v8b * CAVITY_LENGTH_v8b +
        pi * METAL_OUTER_RADIUS_v8b^2
    const CAVITY_HEAT_CAPACITY_v8b =
        INSULATION_DENSITY_v8b * INSULATION_HEAT_CAPACITY_v8b *
        CAVITY_FELT_VOLUME_v8b +
        ALUMINUM_DENSITY_v8b * ALUMINUM_HEAT_CAPACITY_v8b *
        CAVITY_METAL_VOLUME_v8b
end

# Final v8b parameter vector p[1:10]
# Cooling / thermophysical:
#   p[1]  gamma_C       receiver solid heat-capacity multiplier
#   p[2]  k_scale       axial receiver conductivity multiplier
#   p[3]  k_ins_scale   receiver/tube-to-cavity conductance multiplier
#   p[7]  tau_T3        retained for vector compatibility; unused in v8b
# Heating:
#   p[4]  A_Nu          Nusselt multiplier exponent in 10^A_Nu
#   p[5]  h_floor       downstream exchange floor, 0<h_floor<=1
#   p[6]  L_h           axial exchange decay length [m]
#   p[8]  f_I_high      irradiance correction for high-Io group
#   p[9]  f_I_mid       irradiance correction for mid-Io group
#   p[10] f_I_low       irradiance correction for low-Io group
begin # v8b parameter values and bounds
    pnew_v8b = [
        1.0, 1.0, 1.0,
        -1.0, 0.25, 0.050,
        20.0,
        1.0, 1.0, 1.0,
    ]
    lb_full_v8b = [
        0.20, 0.20, 0.25,
        -4.00, 0.02, 0.005,
        0.50,
        0.60, 0.60, 0.60,
    ]
    ub_full_v8b = [
        6.00, 3.00, 4.00,
        2.50, 1.00, 1.000,
        500.0,
        2.00, 2.00, 1.60,
    ]

    fit_heat_transfer_indices_v8b = [4, 5, 6, 8, 9, 10]
    fit_cooling_thermal_indices_v8b = [1, 2, 3]
end

irradiance_factor_v8b(irradiance, p) =
    irradiance >= 3.80e5 ? p[8] :
    irradiance >= 2.80e5 ? p[9] : p[10]

function axial_exchange_shape_v8b(z, p)
    floor_value = clamp(p[5], 0.0, 1.0)
    length_scale = max(p[6], 1e-6)
    return floor_value + (1.0 - floor_value) * exp(-max(0.0, z) / length_scale)
end

alumina_conductivity_v8b(T) =
    5.5 + 34.5 * exp(-0.0033 * (property_temperature(T) - 273.15))

alumina_heat_capacity_v8b(T) =
    max(500.0, (1.00446 + 1.742e-4 * property_temperature(T) -
                2.796e4 * property_temperature(T)^-2) * 1000.0)

function rear_tube_outer_radius_v8b(z_rear)
    return z_rear <= ADAPTOR_TUBE_OVERLAP_LENGTH_v8b ?
           ADAPTOR_RADIUS_v8b : REAR_TUBE_OUTER_RADIUS_v8b
end

function rear_tube_solid_area_v8b(z_rear)
    outer_radius = rear_tube_outer_radius_v8b(z_rear)
    return pi * (outer_radius^2 - REAR_TUBE_GAS_RADIUS_v8b^2)
end

function rear_tube_cavity_conductance_per_length_v8b(z_rear, p)
    z_rear > REAR_TUBE_CAVITY_LENGTH_v8b && return 0.0
    outer_radius = min(rear_tube_outer_radius_v8b(z_rear), 0.98 * INSULATION_OUTER_RADIUS_v8b)
    return p[3] * 2.0 * pi * INSULATION_CONDUCTIVITY_v8b /
           log(INSULATION_OUTER_RADIUS_v8b / outer_radius)
end

function rear_tube_flange_conductance_per_length_v8b(z_rear, Ttube)
    z_rear <= REAR_TUBE_CAVITY_LENGTH_v8b && return 0.0
    tube_wall_resistance =
        log(REAR_TUBE_OUTER_RADIUS_v8b / REAR_TUBE_GAS_RADIUS_v8b) /
        (2.0 * pi * alumina_conductivity_v8b(Ttube))
    aluminum_resistance =
        log(METAL_OUTER_RADIUS_v8b / REAR_TUBE_OUTER_RADIUS_v8b) /
        (2.0 * pi * ALUMINUM_CONDUCTIVITY_v8b)
    return 1.0 / (tube_wall_resistance + aluminum_resistance)
end

function cavity_loss_to_ambient_v8b(Tcavity, ambient)
    return H_NAT_CAVITY_v8b * CAVITY_OUTER_AREA_v8b * (Tcavity - ambient) +
           METAL_EMISSIVITY_v8b * sigma_sb * CAVITY_OUTER_AREA_v8b *
           (Tcavity^4 - ambient^4)
end

function rear_tube_capacity_v8b(Ttube, z_rear, dx_rear)
    return ALUMINA_DENSITY_v8b * alumina_heat_capacity_v8b(Ttube) *
           rear_tube_solid_area_v8b(z_rear) * dx_rear
end

function rear_tube_htc_v8b(Twall, Tgas, flow_lpm, Tin)
    flow_lpm <= 1e-12 && return 0.0
    mdot = m_dot(flow_lpm, Tin)
    Tfilm = 0.5 * (Twall + Tgas)
    cp = cpf_f(Tfilm)
    mu = mu_f_f(Tfilm)
    kf = kf_f(Tfilm)
    reynolds = mdot * REAR_TUBE_HYDRAULIC_DIAMETER_v8b /
               (REAR_TUBE_FLOW_AREA_v8b * mu)
    prandtl = cp * mu / kf
    nusselt = reynolds < 2300.0 ?
              3.66 :
              0.023 * reynolds^0.8 * prandtl^0.4
    return nusselt * kf / REAR_TUBE_HYDRAULIC_DIAMETER_v8b
end

function gas_profile_v8b!(Tg, Qgas, hcell, Qgas_rear, hrear,
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
        nusselt = 10.0^p[4] *
                  reynolds^B_RE_FIXED_v8b *
                  prandtl^(10.0^C_PR_FIXED_v8b)
        hcell[i] = nusselt * kf / Dh * axial_exchange_shape_v8b(z[i], p)
        UA = hcell[i] * P_exchange * dx
        effectiveness = -expm1(-UA / (mdot * cp))
        Tg[i + 1] = Tg[i] + effectiveness * (Ts[i] - Tg[i])
        Qgas[i] = mdot * cp * (Tg[i + 1] - Tg[i])
    end

    for j in eachindex(Ttube)
        inlet_index = receiver_nodes + j
        Tfilm = 0.5 * (Ttube[j] + Tg[inlet_index])
        cp = cpf_f(Tfilm)
        hrear[j] = rear_tube_htc_v8b(Ttube[j], Tg[inlet_index], flow, Tin)
        UA = hrear[j] * REAR_TUBE_EXCHANGE_PERIMETER_v8b * dx_rear
        effectiveness = -expm1(-UA / (mdot * cp))
        Tg[inlet_index + 1] = Tg[inlet_index] + effectiveness *
                              (Ttube[j] - Tg[inlet_index])
        Qgas_rear[j] = mdot * cp * (Tg[inlet_index + 1] - Tg[inlet_index])
    end
    return nothing
end

function receiver_rhs_v8b!(du, u, time, p, operating, z, dx, solar_weights,
                          z_rear, dx_rear,
                          Tg, Qgas, hcell, Qgas_rear, hrear)
    nodes = length(z)
    rear_nodes = length(z_rear)
    Ts = view(u, 1:nodes)
    Ttube = view(u, (nodes + 1):(nodes + rear_nodes))
    Tcavity = u[nodes + rear_nodes + 1]
    dTs = view(du, 1:nodes)
    dTtube = view(du, (nodes + 1):(nodes + rear_nodes))
    ambient = operating.ambient_temperature(time)

    gas_profile_v8b!(Tg, Qgas, hcell, Qgas_rear, hrear,
                    Ts, Ttube, time, p, operating, z, dx, z_rear, dx_rear)
    fill!(du, 0.0)

    for i in 1:(nodes - 1)
        ki = p[2] * ks_f(Ts[i])
        kj = p[2] * ks_f(Ts[i + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        Qcond = kface * A_solid * (Ts[i] - Ts[i + 1]) / dx
        dTs[i] -= Qcond
        dTs[i + 1] += Qcond
    end

    irradiance = max(0.0, operating.irradiance(time))
    Qsolar = ETA_ABS_FIXED_v8b * irradiance_factor_v8b(irradiance, p) *
             irradiance * A_frt
    radial_conductance = p[3] * RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v8b
    for i in 1:nodes
        Qside_cavity = radial_conductance * dx * (Ts[i] - Tcavity)
        dTs[i] += Qsolar * solar_weights[i] - Qgas[i] - Qside_cavity
        du[nodes + rear_nodes + 1] += Qside_cavity
    end

    Qfront = H_FRONT_FIXED_v8b * A_frt * (Ts[1] - ambient) +
             EPS_FRONT_FIXED_v8b * sigma_sb * A_frt * (Ts[1]^4 - ambient^4)
    Qreceiver_tube = RECEIVER_TO_TUBE_CONDUCTANCE_v8b * (Ts[end] - Ttube[1])

    dTs[1] -= Qfront
    dTs[end] -= Qreceiver_tube
    dTtube[1] += Qreceiver_tube

    for j in 1:(rear_nodes - 1)
        ki = alumina_conductivity_v8b(Ttube[j])
        kj = alumina_conductivity_v8b(Ttube[j + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        Aface = 0.5 * (rear_tube_solid_area_v8b(z_rear[j]) +
                       rear_tube_solid_area_v8b(z_rear[j + 1]))
        Qcond = kface * Aface * (Ttube[j] - Ttube[j + 1]) / dx_rear
        dTtube[j] -= Qcond
        dTtube[j + 1] += Qcond
    end

    for j in 1:rear_nodes
        G_cavity = rear_tube_cavity_conductance_per_length_v8b(z_rear[j], p)
        G_flange = rear_tube_flange_conductance_per_length_v8b(z_rear[j], Ttube[j])
        Q_cavity = G_cavity * dx_rear * (Ttube[j] - Tcavity)
        Q_flange = G_flange * dx_rear * (Ttube[j] - WATER_FLANGE_TEMPERATURE_v8b)
        dTtube[j] -= Qgas_rear[j] + Q_cavity + Q_flange
        du[nodes + rear_nodes + 1] += Q_cavity
    end

    Qcavity_ambient = cavity_loss_to_ambient_v8b(Tcavity, ambient)
    du[nodes + rear_nodes + 1] -= Qcavity_ambient

    cell_volume = A_solid * dx
    for i in 1:nodes
        capacity = rho_s * p[1] * Cps_f(Ts[i]) * cell_volume
        dTs[i] /= capacity
    end
    for j in 1:rear_nodes
        dTtube[j] /= max(rear_tube_capacity_v8b(Ttube[j], z_rear[j], dx_rear),
                         eps(Float64))
    end
    du[nodes + rear_nodes + 1] /= max(CAVITY_HEAT_CAPACITY_v8b, eps(Float64))
    return nothing
end

function receiver_ode_v8b!(dTs, Ts, context, time)
    receiver_rhs_v8b!(
        dTs, Ts, time, context.p, context.operating,
        context.z, context.dx, context.solar_weights,
        context.z_rear, context.dx_rear,
        context.Tg, context.Qgas, context.hcell, context.Qgas_rear, context.hrear,
    )
    return nothing
end

function simulate_v8b(p, operating, save_times;
                     initial_temperature=Tamb, nodes=default_nodes,
                     rear_nodes=REAR_TUBE_DEFAULT_NODES_v8b,
                     initial_rear_temperature=nothing,
                     initial_cavity_temperature=Tamb,
                     solver=Rodas5P(autodiff=AutoFiniteDiff()),
                     reltol=1e-6, abstol=1e-7, dtmax=30.0)
    length(p) == 10 || throw(ArgumentError("1D_v8b expects 10 parameters"))
    all(isfinite, p) || throw(ArgumentError("all model parameters must be finite"))
    p[1] > 0.0 && p[2] > 0.0 && p[3] > 0.0 ||
        throw(ArgumentError("capacity, conductivity, and insulation scale must be positive"))
    nodes >= 3 || throw(ArgumentError("at least three axial nodes are required"))
    rear_nodes >= 2 || throw(ArgumentError("at least two rear tube nodes are required"))

    time = Float64.(save_times)
    isempty(time) && throw(ArgumentError("save_times cannot be empty"))
    (length(time) == 1 || all(diff(time) .> 0.0)) ||
        throw(ArgumentError("save_times must be strictly increasing"))

    dx = L / nodes
    dx_rear = REAR_TUBE_LENGTH_v8b / rear_nodes
    z = collect(range(dx / 2, step=dx, length=nodes))
    z_rear = collect(range(dx_rear / 2, step=dx_rear, length=rear_nodes))
    zg = vcat(
        collect(range(0.0, L, length=nodes + 1)),
        L .+ collect(range(dx_rear, step=dx_rear, length=rear_nodes)),
    )
    weights = solar_weights_v5(BETA_OPT_FIXED_v8b, FRONT_DEPOSITION_FIXED_v8b, nodes)
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
    )
    problem = ODEProblem(receiver_ode_v8b!, u_initial, (time[1], time[end]), context)
    solution = solve(
        problem, solver; saveat=time, save_everystep=false, dense=false,
        reltol=reltol, abstol=abstol, dtmax=dtmax, tstops=time,
    )
    successful_retcode(solution) ||
        error("1D_v8b ODE solve failed with retcode $(solution.retcode)")

    state_history = reduce(hcat, solution.u)
    Ts_history = state_history[1:nodes, :]
    rear_temperature_history .= state_history[(nodes + 1):(nodes + rear_nodes), :]
    cavity_temperature_history .= state_history[nodes + rear_nodes + 1, :]
    for output_index in eachindex(time)
        t = time[output_index]
        Ts_now = view(Ts_history, :, output_index)
        Ttube_now = view(rear_temperature_history, :, output_index)
        Tcavity_now = cavity_temperature_history[output_index]
        gas_profile_v8b!(Tg, Qgas, hcell, Qgas_rear, hrear, Ts_now, Ttube_now,
                        t, p, operating, z, dx, z_rear, dx_rear)
        Tg_history[:, output_index] .= Tg
        h_history[:, output_index] .= hcell
        h_rear_history[:, output_index] .= hrear
        receiver_to_tube_heat_history[output_index] =
            RECEIVER_TO_TUBE_CONDUCTANCE_v8b * (Ts_now[end] - Ttube_now[1])
        receiver_to_cavity_heat_history[output_index] = sum(
            p[3] * RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v8b * dx *
            (Ts_now[i] - Tcavity_now)
            for i in 1:nodes
        )
        tube_to_cavity_heat_history[output_index] = sum(
            rear_tube_cavity_conductance_per_length_v8b(z_rear[j], p) * dx_rear *
            (Ttube_now[j] - Tcavity_now)
            for j in 1:rear_nodes
        )
        cavity_ambient_heat_history[output_index] =
            cavity_loss_to_ambient_v8b(Tcavity_now, operating.ambient_temperature(t))
        flange_heat_history[output_index] = sum(
            rear_tube_flange_conductance_per_length_v8b(z_rear[j], Ttube_now[j]) *
            dx_rear * (Ttube_now[j] - WATER_FLANGE_TEMPERATURE_v8b)
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
        rear_tube_heat_transfer_coefficient=h_rear_history,
        receiver_to_tube_heat=receiver_to_tube_heat_history,
        receiver_to_cavity_heat=receiver_to_cavity_heat_history,
        tube_to_cavity_heat=tube_to_cavity_heat_history,
        cavity_ambient_heat_loss=cavity_ambient_heat_history,
        flange_heat_loss=flange_heat_history,
        ode_solution=solution,
    )
end

begin # v8b experimental interface
    paired_observation_v8b(data, simulation_id, obs_a, obs_b) =
        0.5 .* (observation(data, simulation_id, obs_a) .+
                observation(data, simulation_id, obs_b))

    function measured_initial_profile_v8b(data, simulation_id)
        z8 = sensor_positions[:T8]
        z9 = sensor_positions[:T9]
        z10 = sensor_positions[:T10]
        T8 = observation(data, simulation_id, "_T8")[1]
        T9 = paired_observation_v8b(data, simulation_id, "_T9", "_T12")[1]
        T10 = paired_observation_v8b(data, simulation_id, "_T10", "_T11")[1]
        return function (z)
            z <= z8 && return T8
            z >= z10 && return T10
            if z <= z9
                return T8 + (T9 - T8) * (z - z8) / (z9 - z8)
            end
            return T9 + (T10 - T9) * (z - z9) / (z10 - z9)
        end
    end

    function measured_initial_rear_profile_v8b(data, simulation_id)
        receiver_profile = measured_initial_profile_v8b(data, simulation_id)
        T_receiver_exit = receiver_profile(L)
        T_tube_exit = observation(data, simulation_id, "_Tf")[1]
        return function (z_rear)
            fraction = clamp(z_rear / REAR_TUBE_LENGTH_v8b, 0.0, 1.0)
            return T_receiver_exit + fraction * (T_tube_exit - T_receiver_exit)
        end
    end

    function experimental_case_v8b(simulation_id; is_cooling=false)
        data = is_cooling ? measurements_cooling : measurements
        conditions = is_cooling ? simulation_conditions_cooling : simulation_conditions
        time = observation_time(data, simulation_id)
        experiment = hcat(
            observation(data, simulation_id, "_T8"),
            paired_observation_v8b(data, simulation_id, "_T9", "_T12"),
            paired_observation_v8b(data, simulation_id, "_T10", "_T11"),
            observation(data, simulation_id, "_Tf"),
            observation(data, simulation_id, "_T2"),
        )
        return conditions[simulation_id], time, experiment, data
    end

    function gas_at_position_v8b(result, position)
        z = result.z_gas
        values = result.gas_temperature
        position <= z[1] && return vec(values[1, :])
        position >= z[end] && return vec(values[end, :])
        i = searchsortedlast(z, position)
        fraction = (position - z[i]) / (z[i + 1] - z[i])
        return vec((1.0 - fraction) .* values[i, :] .+ fraction .* values[i + 1, :])
    end

    function extract_outputs_v8b(result, p, T3_initial)
        T8 = solid_at_v3(result, sensor_positions[:T8])
        T9_pair = solid_at_v3(result, sensor_positions[:T9])
        T10_pair = solid_at_v3(result, sensor_positions[:T10])
        Tf = gas_at_position_v8b(result, T3_SAMPLE_POSITION_v8b)
        T2 = result.cavity_temperature
        h_mean = vec(mean(result.heat_transfer_coefficient, dims=1))
        return hcat(T8, T9_pair, T10_pair, Tf, T2, h_mean)
    end

    function solve_case_v8b(p, simulation_id; is_cooling=false,
                           nodes=default_nodes,
                           rear_nodes=REAR_TUBE_DEFAULT_NODES_v8b,
                           solver=Rodas5P(autodiff=AutoFiniteDiff()),
                           reltol=1e-6, abstol=1e-7, dtmax=30.0)
        conditions, time, experiment, data = experimental_case_v8b(
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
        result = simulate_v8b(
            p, operating, time;
            initial_temperature=measured_initial_profile_v8b(data, simulation_id),
            initial_rear_temperature=measured_initial_rear_profile_v8b(data, simulation_id),
            initial_cavity_temperature=T2[1],
            nodes=nodes, rear_nodes=rear_nodes, solver=solver,
            reltol=reltol, abstol=abstol, dtmax=dtmax,
        )
        outputs = extract_outputs_v8b(result, p, experiment[1, 4])
        return outputs, result, experiment
    end
end

begin # v8b calibration
    const TIMING_WEIGHT_v8b = 0.15
    const SLOPE_WEIGHT_v8b = 0.25

    heat_shape_regularization_v8b(p) =
        0.001 * (log(p[5] / 0.25)^2 + log(p[6] / 0.050)^2)

    function normalized_slope_mse_v8b(model, experiment)
        length(model) < 2 && return 0.0
        scale = max(maximum(experiment) - minimum(experiment), 1.0)
        return mean(abs2, diff(model) .- diff(experiment)) / scale^2
    end

    function timing_penalty_v8b(time, model, experiment)
        duration = max(time[end] - time[1], 1.0)
        return ((get_t90_v3(time, model) - get_t90_v3(time, experiment)) / duration)^2
    end

    function signal_loss_v8b(time, model, experiment)
        return normalized_mse_v3(model, experiment) +
               SLOPE_WEIGHT_v8b * normalized_slope_mse_v8b(model, experiment) +
               TIMING_WEIGHT_v8b * timing_penalty_v8b(time, model, experiment)
    end

    function loss_cases_v8b(p, keys; is_cooling=false, nodes=15)
        total = 0.0
        for simulation_id in keys
            try
                model, result, experiment = solve_case_v8b(
                    p, simulation_id; is_cooling=is_cooling, nodes=nodes,
                    reltol=2e-5, abstol=2e-6, dtmax=60.0,
                )
                all(isfinite, model) || return Inf
                total += mean(
                    signal_loss_v8b(result.time, model[:, j], experiment[:, j])
                    for j in 1:5
                )
            catch err
                err isa InterruptException && rethrow()
                return Inf
            end
        end
        return total / length(keys)
    end

    function loss_heating_v8b(p, keys=sim_key_heat; nodes=15)
        regularization = heat_shape_regularization_v8b(p) +
                         0.003 * sum(log(p[i])^2 for i in 8:10)
        return loss_cases_v8b(p, keys; is_cooling=false, nodes=nodes) + regularization
    end

    function loss_cooling_v8b(p, keys=sim_key_cool; nodes=15)
        regularization = 0.005 * (log(p[1])^2 + log(p[2])^2 + log(p[3])^2)
        return loss_cases_v8b(p, keys; is_cooling=true, nodes=nodes) + regularization
    end
end

function optimize_parameter_subset_v8b(objective, p_initial, indices, lower, upper;
                                      maximum_iterations=200,
                                      maximum_time_seconds=120.0,
                                      optimizer=NLopt.LN_NELDERMEAD(),
                                      label="subset-v8b")
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

function calibrate_v8b(; cooling_keys=sim_key_cool,
                      heating_keys=sim_key_heat,
                      nodes=15,
                      heating_base_iterations=20,
                      cooling_iterations=20,
                      heating_refit_iterations=20,
                      cooling_time_seconds=120.0,
                      heating_base_time_seconds=120.0,
                      heating_refit_time_seconds=120.0,
                      optimizer=NLopt.LN_NELDERMEAD())
    p0 = copy(pnew_v8b)

    println("\n=== 1D_v8b STAGE 1: heating heat-transfer + irradiance factors ===")
    heating_base_fit = optimize_parameter_subset_v8b(
        p -> loss_heating_v8b(p, heating_keys; nodes=nodes),
        p0, fit_heat_transfer_indices_v8b, lb_full_v8b, ub_full_v8b;
        maximum_iterations=nlopt_evaluation_budget_v3(
            heating_base_iterations, length(fit_heat_transfer_indices_v8b)
        ),
        maximum_time_seconds=heating_base_time_seconds,
        optimizer=optimizer,
        label="heating-base-v8b",
    )
    p_after_heating_base = heating_base_fit.parameters

    println("\n=== 1D_v8b STAGE 2: cooling thermophysical / predicted-cavity calibration ===")
    cooling_fit = optimize_parameter_subset_v8b(
        p -> loss_cooling_v8b(p, cooling_keys; nodes=nodes),
        p_after_heating_base, fit_cooling_thermal_indices_v8b, lb_full_v8b, ub_full_v8b;
        maximum_iterations=nlopt_evaluation_budget_v3(
            cooling_iterations, length(fit_cooling_thermal_indices_v8b)
        ),
        maximum_time_seconds=cooling_time_seconds,
        optimizer=optimizer,
        label="cooling-v8b",
    )
    p_after_cooling = cooling_fit.parameters

    println("\n=== 1D_v8b STAGE 3: heating heat-transfer + irradiance refit ===")
    heating_refit = optimize_parameter_subset_v8b(
        p -> loss_heating_v8b(p, heating_keys; nodes=nodes),
        p_after_cooling, fit_heat_transfer_indices_v8b, lb_full_v8b, ub_full_v8b;
        maximum_iterations=nlopt_evaluation_budget_v3(
            heating_refit_iterations, length(fit_heat_transfer_indices_v8b)
        ),
        maximum_time_seconds=heating_refit_time_seconds,
        optimizer=optimizer,
        label="heating-refit-v8b",
    )
    parameters = heating_refit.parameters

    global pnew_v8b = parameters
    global calibration_v8b = (
        heating_base=heating_base_fit,
        cooling=cooling_fit,
        heating_refit=heating_refit,
        fixed=(
            eta_abs=ETA_ABS_FIXED_v8b,
            beta_opt=BETA_OPT_FIXED_v8b,
            front_deposition_fraction=FRONT_DEPOSITION_FIXED_v8b,
            h_front=H_FRONT_FIXED_v8b,
            eps_front=EPS_FRONT_FIXED_v8b,
            cavity_outer_radius=CAVITY_OUTER_RADIUS_v8b,
            insulation_outer_radius=INSULATION_OUTER_RADIUS_v8b,
            receiver_to_cavity_conductance_per_length=RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v8b,
            receiver_to_tube_conductance=RECEIVER_TO_TUBE_CONDUCTANCE_v8b,
            rear_tube_length=REAR_TUBE_LENGTH_v8b,
            rear_tube_cavity_length=REAR_TUBE_CAVITY_LENGTH_v8b,
            rear_tube_flange_length=REAR_TUBE_FLANGE_LENGTH_v8b,
            t3_sample_position=T3_SAMPLE_POSITION_v8b,
            rear_tube_wall_thickness=REAR_TUBE_WALL_THICKNESS_v8b,
            cavity_heat_capacity=CAVITY_HEAT_CAPACITY_v8b,
            water_flange_temperature=WATER_FLANGE_TEMPERATURE_v8b,
        ),
        parameters=parameters,
    )
    println("\nFinal 1D_v8b parameters:")
    println(parameters)
    return calibration_v8b
end

function quick_calibration_v8b()
    return calibrate_v8b(
        cooling_keys=["C69"], heating_keys=["E74"],
        nodes=11,
        heating_base_iterations=2,
        cooling_iterations=2,
        heating_refit_iterations=2,
    )
end

begin # v8b post-processing
    function transient_case_v8b(simulation_id, p=pnew_v8b; is_cooling=false,
                               nodes=default_nodes)
        model, result, experiment = solve_case_v8b(
            p, simulation_id; is_cooling=is_cooling, nodes=nodes
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
            T3_sample_gas=gas_at_position_v8b(result, T3_SAMPLE_POSITION_v8b),
            rear_tube_inlet=vec(result.rear_tube_temperature[1, :]),
            rear_tube_exit=vec(result.rear_tube_temperature[end, :]),
            receiver_to_tube_heat=result.receiver_to_tube_heat,
            receiver_to_cavity_heat=result.receiver_to_cavity_heat,
            tube_to_cavity_heat=result.tube_to_cavity_heat,
            cavity_ambient_heat_loss=result.cavity_ambient_heat_loss,
            flange_heat_loss=result.flange_heat_loss,
            h_effective=model[:, 6],
            full_result=result,
        )
    end

    function plot_case_v8b(simulation_id, params=pnew_v8b; is_cooling=false,
                          nodes=default_nodes)
        @eval using StatsPlots
        data = transient_case_v8b(simulation_id, params;
                                 is_cooling=is_cooling, nodes=nodes)
        sensors = (:T8, :T9_pair, :T10_pair, :T3, :T2)
        colors = (:blue, :red, :green, :purple, :black)
        plot_object = plot(title="1D_v8b transient: $simulation_id",
                           xlabel="Time (s)", ylabel="Temperature (K)",
                           legend=:outerright)
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

    function compute_metrics_v8b(p=pnew_v8b; heating_keys=sim_key_heat,
                                cooling_keys=sim_key_cool, nodes=default_nodes)
        metrics = NamedTuple[]
        for (keys, cooling) in ((heating_keys, false), (cooling_keys, true))
            for simulation_id in keys
                model, result, experiment = solve_case_v8b(
                    p, simulation_id; is_cooling=cooling, nodes=nodes
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
                        shape_loss=normalized_slope_mse_v8b(
                            model[:, column], experiment[:, column]
                        ),
                    ))
                end
            end
        end
        return metrics
    end
end

println("""
[1D_v8b] Model loaded.
Run quick_calibration_v8b() for a short check.
Run calibration_v8b = calibrate_v8b() for the full staged calibration.
The calibrated vector is stored in pnew_v8b.
""")
