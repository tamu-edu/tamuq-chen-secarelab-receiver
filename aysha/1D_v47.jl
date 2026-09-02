# ============================================================================
# 1D_v47.jl - Entire Converter Model for Monolithic Solar Receivers
# ============================================================================
#
# Key Features & Refinements:
# 1. Geometry & Material Synchronization (1D_v3.jl baseline):
#    - 1.50 mm channel opening, 0.40 mm web, 1.90 mm cell pitch, N_ch = 100
#    - Monolith cross-section 19x19 mm (A_frt = 3.61e-4 m^2), L = 137 mm
#    - Void porosity epsilon = 0.623269, A_solid = 1.361e-4 m^2
#    - Temperature-dependent SiC & Air properties via polynomial functions
#    - Sensor coordinates: T8=11mm, T12=58mm, T11=107mm, T9=58mm, T10=107mm, T3=140mm
# 2. 100% Mass Flow & Enthalpy Conservation:
#    - Entire mass flow m_dot passes through core matrix and rear exit tube (phi_act = 1.0)
# 3. Front Aperture Suction Preheating (z = 0):
#    - Inward flowing air preheats across front aperture face: Q_suction = h_suction * A_frt * (T_perim,1 - T_in)
#    - Bounded h_suction in [10, 150] W/m^2K ensuring volumetric core heating dominates
#    - Radiative loss Q_front_rad radiates to ambient
# 4. Front Perimeter Spillage & Aluminum Casing Conduction:
#    - Spilled power (1 - chi) * Q_delivered is absorbed at front face (z=0, cell 1)
#    - Axial transport governed by aluminum casing conductivity k_perim_ref
# 5. Phase-Invariant Boundary Couplings:
#    - Flange conductance and boundary properties are invariant across heating and cooling
# 6. Physical Zero-Flow Physics:
#    - Proper natural convection & radiation modeling for T3 sensor under zero flow (C80)
# 7. Exact Energy Conservation Audit:
#    - Rigorous first-law audit function compute_energy_balance_v47
# ============================================================================

using DifferentialEquations
using LinearAlgebra
using Statistics
using Optimization
using OptimizationNLopt

include("1D_v5.jl")

begin # physical and geometrical constants for v47
    # Use existing base geometric values from 1D_v3
    const N_channels_v47 = 100           # Number of square channels
    const n_refr_v47 = 1.0               # Refractive index in channels
    const default_nodes_v47 = 15         # Default axial nodes
    const w_ch_v47 = 0.0015              # Channel opening width (1.50 mm)
    const w_cell_v47 = 0.019 / sqrt(N_channels_v47) # 0.0019 m (1.90 mm)
    const t_web_nominal_v47 = w_cell_v47 - w_ch_v47    # 0.0004 m (0.40 mm)
    const Dh_v47 = w_ch_v47              # Hydraulic diameter (1.50 mm)
    const P_exchange_v47 = 4.0 * N_channels_v47 * w_ch_v47 # Internal wetted perimeter (0.60 m)
    const a_v_v47 = P_exchange_v47 / A_frt   # Specific exchange area (1662.05 m^2/m^3)
    const epsilon_v47 = (N_channels_v47 * w_ch_v47^2) / A_frt # Void porosity (0.623269)
    const A_solid_v47 = (1.0 - epsilon_v47) * A_frt       # Solid matrix area (1.361e-4 m^2)
    const rho_s_v47 = 3100.0             # SiC density [kg/m^3]
    const sigma_sb_v47 = 5.670374419e-8  # Stefan-Boltzmann constant [W/m^2 K^4]

    # Packaging & Housing Geometry
    const R_core_v47 = sqrt(A_frt / pi)  # 0.010735 m (10.74 mm)
    const R_felt_v47 = 0.025             # Outer felt radius (25 mm)
    const R_case_v47 = 0.030             # Outer casing radius (30 mm)
    const L_case_v47 = L                 # Casing length (137 mm)
    const A_case_ext_v47 = 2.0 * pi * R_case_v47 * L_case_v47 # External casing area (0.02582 m^2)
    const k_felt_val_v47 = 0.06          # Alumina-silicate felt conductivity [W/m K]

    # Rear Hardware Assembly
    const REAR_TUBE_LENGTH_v47 = 0.063       # Rear exit tube length (63 mm)
    const REAR_TUBE_DEFAULT_NODES_v47 = 15   # Rear tube axial nodes
    const R_tube_inner_v47 = 0.007           # Rear tube gas radius (7 mm)
    const R_tube_outer_v47 = 0.011           # Rear tube outer radius (11 mm)
    const WATER_FLANGE_TEMPERATURE_v47 = 298.15 # Water-cooled flange temperature (25 C)
    const COOLING_ROOM_TEMPERATURE_v47 = 295.15 # Room temperature during cooling (22 C)

    # Capacitances
    const CAVITY_HEAT_CAPACITY_v47 = 4026.0  # Metal cavity / housing shell [J/K]
    const T3_SENSOR_CAPACITY_v47 = 0.05      # Small TC bead capacity [J/K]
    const MEASURED_ASSEMBLY_CAPACITANCE_v47 = 301.0 # Experimental total capacity target [J/K]
    const CAPACITANCE_REGULARIZATION_WEIGHT_v47 = 2.0
    const ORDERING_WEIGHT_v47 = 0.05
    const COOLING_UPTURN_WEIGHT_v47 = 0.05

    # Radiation & Boundary constants
    const EPS_FRONT_FIXED_v47 = 0.85         # Front aperture emissivity
    const EPS_CASE_EXT_v47 = 0.80            # External casing emissivity
    const H_NAT_EXT_v47 = 8.0                # External natural convection HTC [W/m^2 K]
    const NU_FLOOR_v47 = 0.0                 # Nu floor (Option A / clean formulation)

    const T3_SAMPLE_POSITION_v47 = 0.140
    const CAVITY_LENGTH_v47 = 0.0685
    const REAR_TUBE_CAVITY_LENGTH_v47 = 0.030
end

begin # thermophysical material properties (1D_v3.jl exact polynomials)
    # Silicon Carbide (SiC) Monolith
    function ks_f(T)
        Tc = T - 273.15
        return max(1.0, 51.0 * exp(-0.0030 * Tc) + 1.2)
    end

    function Cps_f(T)
        Tc = T - 273.15
        return max(500.0, 900.0 + 0.30 * Tc - 3.0e5 * (T^-2))
    end

    # Fluid (Air) Properties
    const p_atm = 101325.0
    const R_spec_air = 287.05

    rhof_f(T) = p_atm / (R_spec_air * max(T, 100.0))
    cpf_f(T) = 1005.0 + 0.05 * (T - 273.15)
    muf_f(T) = 1.81e-5 * (max(T, 100.0) / 293.15)^0.70
    kf_f(T) = 0.0257 * (max(T, 100.0) / 293.15)^0.80
    Prf_f(T) = (cpf_f(T) * muf_f(T)) / kf_f(T)

    # Standard air mass flow conversion: 1 LPM at standard conditions
    const rho_std_air = 101325.0 / (287.05 * 293.15) # 1.2041 kg/m^3
    m_dot_standard_v47(flow_lpm) = (flow_lpm / 60000.0) * rho_std_air

    # Alumina Rear Tube Properties
    function k_al_f(T)
        Tc = T - 273.15
        return 5.5 + 34.5 * exp(-0.0033 * Tc)
    end

    function cp_al_f(T)
        return max(500.0, (1.00446 + 1.742e-4 * T - 2.796e4 * (T^-2)) * 1000.0)
    end
    const rho_al = 3900.0
end

begin # parameter vector definition & bounds (23 parameters)
    # Calibrated baseline parameter vector for 1D_v47
    const pnew_v47 = [
        1.0,                    # 1: A_Nu (Nusselt baseline constant)
        0.9283026938833806,     # 2: B_Re (Reynolds exponent)
        0.3333,                 # 3: C_Pr (Prandtl exponent, fixed)
        0.4551200986314161,     # 4: front_dep (front-face absorption fraction)
        1.3400,                 # 5: scale_456 (456 kW/m^2 cluster power scale)
        1.5800,                 # 6: scale_304 (304 kW/m^2 cluster power scale)
        1.1100,                 # 7: scale_256 (256 kW/m^2 cluster power scale)
        24.872039140594868,     # 8: G_core_perim (W/m K)
        137.70432244158707,     # 9: C_perim_eff (J/K)
        9.10064389265549,       # 10: k_perim_ref (W/m K)
        193.77838332109653,     # 11: beta_opt (1/m)
        0.7958715381834904,     # 12: chi (core solar fraction)
        0.9728454436569068,     # 13: f_core_rear (rear split core fraction)
        0.19925163375662922,    # 14: flange_scale (water flange conductance scale)
        0.26501771319320255,    # 15: k_core_axial_scale (core axial conduction scale)
        102.20711085813221,     # 16: C_rear_eff (J/K)
        0.8312024110892697,     # 17: G_receiver_rear (W/K)
        6.9838677475793896,     # 18: G_rear_tube (W/K)
        0.01,                   # 19: G_rear_cavity (W/K)
        8.658389709221446,      # 20: G_rear_axial (W/K)
        8.736702133267356e-5,   # 21: delta_web (web thickness correction, m)
        1.529840450127574,      # 22: C_z (entry length scale factor)
        150.0,                  # 23: h_suction (aperture suction HTC, W/m^2 K)
    ]

    const lb_full_v47 = [
        1.0, 0.0, 0.3333, 0.05, 1.0, 1.0, 0.8, 0.1, 50.0, 1.0,
        20.0, 0.5, 0.5, 0.05, 0.05, 20.0, 0.01, 0.01, 0.01, 0.01,
        1.0e-5, 0.1, 10.0
    ]

    const ub_full_v47 = [
        15.0, 1.0, 0.3333, 0.60, 2.0, 2.0, 1.8, 30.0, 400.0, 60.0,
        350.0, 1.0, 1.0, 5.0, 2.0, 300.0, 5.0, 10.0, 5.0, 20.0,
        2.0e-3, 5.0, 150.0
    ]

    const fit_rear_stage_indices_v47 = [14, 17, 18, 19, 20]
    const fit_transport_stage_indices_v47 = [1, 2, 8, 10, 14, 15, 17, 18, 19, 20, 21, 22, 23]
    const fit_full_stage_indices_v47 = [1, 2, 4, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]

    # Parameter Accessors
    nusselt_a_v47(p) = clamp(p[1], lb_full_v47[1], ub_full_v47[1])
    nusselt_b_v47(p) = clamp(p[2], lb_full_v47[2], ub_full_v47[2])
    core_front_deposition_v47(p) = clamp(p[4], lb_full_v47[4], ub_full_v47[4])
    core_perimeter_conductance_v47(p) = clamp(p[8], lb_full_v47[8], ub_full_v47[8])
    perimeter_heat_capacity_total_v47(p) = clamp(p[9], lb_full_v47[9], ub_full_v47[9])
    perimeter_conductivity_ref_v47(p) = clamp(p[10], lb_full_v47[10], ub_full_v47[10])
    core_extinction_coefficient_v47(p) = clamp(p[11], lb_full_v47[11], ub_full_v47[11])
    core_source_fraction_v47(p) = clamp(p[12], lb_full_v47[12], ub_full_v47[12])
    rear_core_fraction_v47(p) = clamp(p[13], lb_full_v47[13], ub_full_v47[13])
    core_axial_conduction_scale_v47(p) = clamp(p[15], lb_full_v47[15], ub_full_v47[15])
    rear_reservoir_heat_capacity_v47(p) = clamp(p[16], lb_full_v47[16], ub_full_v47[16])
    web_thickness_correction_v47(p) = clamp(p[21], lb_full_v47[21], ub_full_v47[21])
    entrance_length_scale_v47(p) = clamp(p[22], lb_full_v47[22], ub_full_v47[22])
    suction_heat_transfer_coefficient_v47(p) = clamp(p[23], lb_full_v47[23], ub_full_v47[23])

    function absorbed_power_scale_v47(irradiance, p)
        I_val = irradiance
        if I_val >= 400000.0
            return clamp(p[5], lb_full_v47[5], ub_full_v47[5])
        elseif I_val >= 280000.0
            return clamp(p[6], lb_full_v47[6], ub_full_v47[6])
        else
            return clamp(p[7], lb_full_v47[7], ub_full_v47[7])
        end
    end
end

begin # constitutive transport closures
    # Developing laminar Nusselt correlation & effective series HTC
    function local_heat_transfer_coefficient_v47(Re, Pr, k_fluid, z_pos, T_solid, p)
        A_nu = nusselt_a_v47(p)
        B_re = nusselt_b_v47(p)
        Cz = entrance_length_scale_v47(p)
        delta_web = web_thickness_correction_v47(p)
        
        # Entrance coordinate
        z_star = max(z_pos, 1.0e-4) / Cz
        Gz = (Dh_v47 / z_star) * Re * Pr
        
        # Developing Nu
        Nu_entry = (0.0668 * Gz) / (1.0 + 0.04 * (Gz^(2.0 / 3.0)))
        Nu_fluid = A_nu + Nu_entry * ((Re * Pr * Dh_v47 / z_star)^B_re)
        Nu_fluid = max(NU_FLOOR_v47, Nu_fluid)
        
        h_fluid = Nu_fluid * k_fluid / Dh_v47
        
        # Intra-strut solid conduction resistance in series
        t_eff = max(1.0e-5, t_web_nominal_v47 + delta_web)
        k_s_val = ks_f(T_solid)
        R_solid = t_eff / (4.0 * k_s_val)
        
        h_eff = 1.0 / (1.0 / h_fluid + R_solid)
        return h_eff, Nu_fluid
    end

    # Core axial effective conductivity (solid conduction + Rosseland radiation)
    function core_effective_axial_conductivity_v47(T, p)
        k_scale = core_axial_conduction_scale_v47(p)
        beta_opt = core_extinction_coefficient_v47(p)
        k_cond = k_scale * (1.0 - epsilon_v47) * ks_f(T)
        k_rad = (16.0 * sigma_sb_v47 * (n_refr_v47^2) * (T^3)) / (3.0 * beta_opt)
        return k_cond + k_rad
    end

    # Perimeter casing axial conductivity
    function perimeter_axial_conductivity_v47(T, p)
        k_ref = perimeter_conductivity_ref_v47(p)
        return k_ref * (1.0 + 0.0005 * (T - 293.15))
    end

    # Perimeter to external cavity conductance
    const G_felt_per_length_v47 = (2.0 * pi * k_felt_val_v47) / log(R_felt_v47 / R_core_v47)

    # Cavity to ambient loss
    function cavity_loss_to_ambient_v47(T_cavity, T_amb)
        Q_conv = H_NAT_EXT_v47 * A_case_ext_v47 * (T_cavity - T_amb)
        Q_rad = EPS_CASE_EXT_v47 * sigma_sb_v47 * A_case_ext_v47 * (T_cavity^4 - T_amb^4)
        return Q_conv + Q_rad
    end

    # Rear tube flange conductance per unit length
    function rear_tube_flange_conductance_per_length_v47(z_rear, T_tube, p)
        s_flange = clamp(p[14], 0.05, 5.0)
        fraction = clamp(z_rear / REAR_TUBE_LENGTH_v47, 0.0, 1.0)
        # Activation towards water flange (z_rear -> L_tube)
        sigmoid = 1.0 / (1.0 + exp(-30.0 * (fraction - 0.70)))
        G_al = (2.0 * pi * k_al_f(T_tube)) / log(R_tube_outer_v47 / R_tube_inner_v47)
        return s_flange * G_al * sigmoid
    end

    # Rear tube sensible heat capacity
    function rear_tube_capacity_v47(T_tube, z_rear, dx_rear)
        V_al = pi * (R_tube_outer_v47^2 - R_tube_inner_v47^2) * dx_rear
        return rho_al * cp_al_f(T_tube) * V_al
    end

    # Distributed rear contact weights
    function rear_contact_weights_v47(z_solid)
        weights = [exp(12.0 * (z / L - 1.0)) for z in z_solid]
        return weights ./ sum(weights)
    end
end

begin # gas profile marching (100% mass flow & zero-flow physics)
    function gas_profile_v47!(Tg, Qgas, hcell, Qgas_rear, hrear,
                             Tcore, Tperim, Ttube, Tin, flow_lpm, p,
                             z_solid, dx, z_rear, dx_rear, zg)
        nodes = length(Tcore)
        rear_nodes = length(Ttube)
        flow = max(0.0, flow_lpm)
        mdot = m_dot_standard_v47(flow)
        h_suction = suction_heat_transfer_coefficient_v47(p)

        if flow > 1e-12
            # 1. Suction Preheating at Aperture Face (z = 0)
            Q_suction = h_suction * A_frt * (Tperim[1] - Tin)
            Tg[1] = Tin + Q_suction / (mdot * cpf_f(Tin))

            # 2. Marching through Core Honeycomb Channels (z in [0, L])
            for i in 1:nodes
                T_solid = Tcore[i]
                T_g_in = Tg[i]
                T_film = 0.5 * (T_solid + T_g_in)
                
                cp_f = cpf_f(T_film)
                k_f = kf_f(T_film)
                mu_f = muf_f(T_film)
                Pr_f = (cp_f * mu_f) / k_f
                
                v_ch = mdot / (N_channels_v47 * rhof_f(T_film) * (w_ch_v47^2))
                Re_ch = (rhof_f(T_film) * v_ch * Dh) / mu_f
                
                h_eff, _ = local_heat_transfer_coefficient_v47(Re_ch, Pr_f, k_f, z_solid[i], T_solid, p)
                hcell[i] = h_eff
                
                NTU = (h_eff * P_exchange_v47 * dx) / (mdot * cp_f)
                T_g_out = T_solid - (T_solid - T_g_in) * exp(-NTU)
                
                Tg[i + 1] = T_g_out
                Qgas[i] = mdot * cp_f * (T_g_out - T_g_in)
            end

            # 3. Marching through Alumina Exit Tube (z in [L, L + L_tube])
            P_tube = 2.0 * pi * R_tube_inner_v47
            for j in 1:rear_nodes
                T_tube_j = Ttube[j]
                T_g_in = Tg[nodes + j]
                T_film = 0.5 * (T_tube_j + T_g_in)
                
                cp_f = cpf_f(T_film)
                k_f = kf_f(T_film)
                mu_f = muf_f(T_film)
                Pr_f = (cp_f * mu_f) / k_f
                
                v_tube = mdot / (pi * (R_tube_inner_v47^2) * rhof_f(T_film))
                Re_tube = (rhof_f(T_film) * v_tube * (2.0 * R_tube_inner_v47)) / mu_f
                
                # Fully developed / developing pipe flow Nu
                Nu_tube = 4.36 + (0.0668 * (2.0 * R_tube_inner_v47 / max(z_rear[j], 1e-4)) * Re_tube * Pr_f) /
                                 (1.0 + 0.04 * ((2.0 * R_tube_inner_v47 / max(z_rear[j], 1e-4)) * Re_tube * Pr_f)^(2.0 / 3.0))
                h_tube = Nu_tube * k_f / (2.0 * R_tube_inner_v47)
                hrear[j] = h_tube
                
                NTU_rear = (h_tube * P_tube * dx_rear) / (mdot * cp_f)
                T_g_out = T_tube_j - (T_tube_j - T_g_in) * exp(-NTU_rear)
                
                Tg[nodes + j + 1] = T_g_out
                Qgas_rear[j] = mdot * cp_f * (T_g_out - T_g_in)
            end
        else
            # Zero-flow physical state (pure stationary conduction/diffusion)
            Tg .= Tin
            Qgas .= 0.0
            hcell .= 0.0
            Qgas_rear .= 0.0
            hrear .= 0.0
        end
        return Tg
    end
end

begin # ODE RHS System
    function receiver_rhs_v47!(du, u, p_tuple, t)
        p, operating, z_solid, dx, core_weights, rear_weights,
        z_rear, dx_rear, Tg, Qgas, hcell, Qgas_rear, hrear, zg = p_tuple

        nodes = length(z_solid)
        rear_nodes = length(z_rear)

        # Unpack state variables
        Tcore = @view u[1:nodes]
        Tperim = @view u[nodes+1:2*nodes]
        Ttube = @view u[2*nodes+1:2*nodes+rear_nodes]
        cavity_index = 2 * nodes + rear_nodes + 1
        Tcavity = u[cavity_index]
        Trear = @view u[cavity_index+1:cavity_index+nodes]
        T3_sensor = u[cavity_index+nodes+1]

        # Derivatives views
        dTcore = @view du[1:nodes]
        dTperim = @view du[nodes+1:2*nodes]
        dTtube = @view du[2*nodes+1:2*nodes+rear_nodes]
        dTrear = @view du[cavity_index+1:cavity_index+nodes]

        # Operating Conditions at time t
        irradiance = operating.irradiance(t)
        flow_lpm = max(0.0, operating.flow_lpm(t))
        Tin = operating.inlet_temperature(t)
        Tamb = operating.ambient_temperature(t)

        # Delivered Power & Source Partition
        Q_aperture = max(0.0, irradiance) * A_frt
        M = absorbed_power_scale_v47(irradiance, p)
        Q_delivered = M * Q_aperture
        chi = core_source_fraction_v47(p)
        Qcore_solar = chi * Q_delivered
        Qperim_solar = (1.0 - chi) * Q_delivered

        # March gas temperatures
        gas_profile_v47!(Tg, Qgas, hcell, Qgas_rear, hrear,
                         Tcore, Tperim, Ttube, Tin, flow_lpm, p,
                         z_solid, dx, z_rear, dx_rear, zg)

        # 1. CORE SOLID MATRIX (Nodes 1:nodes)
        G_cp = core_perimeter_conductance_v47(p)
        G_rec_rear = clamp(p[17], 0.01, 5.0)
        f_core_rear = rear_core_fraction_v47(p)
        cell_vol = A_solid_v47 * dx

        for i in 1:nodes
            C_core_i = rho_s_v47 * Cps_f(Tcore[i]) * cell_vol
            Q_solar_i = Qcore_solar * core_weights[i]
            Q_gas_i = Qgas[i]
            Q_radial_i = G_cp * dx * (Tcore[i] - Tperim[i])
            Q_rear_i = G_rec_rear * rear_weights[i] * f_core_rear * (Tcore[i] - Trear[i])

            # Axial conduction
            Q_axial = 0.0
            if i > 1
                k_prev = 0.5 * (core_effective_axial_conductivity_v47(Tcore[i], p) + core_effective_axial_conductivity_v47(Tcore[i-1], p))
                Q_axial += (k_prev * A_solid_v47 / dx) * (Tcore[i-1] - Tcore[i])
            end
            if i < nodes
                k_next = 0.5 * (core_effective_axial_conductivity_v47(Tcore[i], p) + core_effective_axial_conductivity_v47(Tcore[i+1], p))
                Q_axial += (k_next * A_solid_v47 / dx) * (Tcore[i+1] - Tcore[i])
            end

            dTcore[i] = (Q_solar_i - Q_gas_i - Q_radial_i - Q_rear_i + Q_axial) / C_core_i
        end

        # 2. PERIMETER HOUSING (Nodes 1:nodes)
        C_perim_cell = perimeter_heat_capacity_total_v47(p) / nodes
        A_perim_cs = pi * (R_case_v47^2 - R_felt_v47^2) # Aluminum shell conduction area
        h_suction = suction_heat_transfer_coefficient_v47(p)
        mdot = m_dot_standard_v47(flow_lpm)

        Q_suction = flow_lpm > 1e-12 ? h_suction * A_frt * (Tperim[1] - Tin) : 0.0
        Q_front_rad = EPS_FRONT_FIXED_v47 * sigma_sb_v47 * A_frt * (Tperim[1]^4 - Tamb^4)

        for i in 1:nodes
            Q_solar_perim_i = i == 1 ? Qperim_solar : 0.0
            Q_radial_i = G_cp * dx * (Tcore[i] - Tperim[i])
            Q_cavity_i = G_felt_per_length_v47 * dx * (Tperim[i] - Tcavity)
            Q_rear_perim_i = G_rec_rear * rear_weights[i] * (1.0 - f_core_rear) * (Tperim[i] - Trear[i])

            Q_axial_perim = 0.0
            if i > 1
                k_prev = 0.5 * (perimeter_axial_conductivity_v47(Tperim[i], p) + perimeter_axial_conductivity_v47(Tperim[i-1], p))
                Q_axial_perim += (k_prev * A_perim_cs / dx) * (Tperim[i-1] - Tperim[i])
            end
            if i < nodes
                k_next = 0.5 * (perimeter_axial_conductivity_v47(Tperim[i], p) + perimeter_axial_conductivity_v47(Tperim[i+1], p))
                Q_axial_perim += (k_next * A_perim_cs / dx) * (Tperim[i+1] - Tperim[i])
            end

            Q_boundary = i == 1 ? -(Q_suction + Q_front_rad) : 0.0

            dTperim[i] = (Q_solar_perim_i + Q_radial_i - Q_cavity_i - Q_rear_perim_i + Q_axial_perim + Q_boundary) / C_perim_cell
        end

        # 3. DISTRIBUTED REAR RAIL (Nodes 1:nodes)
        C_rear_total = rear_reservoir_heat_capacity_v47(p)
        G_rear_tube = clamp(p[18], 0.01, 10.0)
        G_rear_cavity = clamp(p[19], 0.01, 5.0)
        G_rear_axial = clamp(p[20], 0.01, 20.0)
        Q_rear_to_tube = G_rear_tube * (Trear[nodes] - Ttube[1])

        for i in 1:nodes
            C_rear_i = C_rear_total * rear_weights[i]
            Q_from_core_i = G_rec_rear * rear_weights[i] * f_core_rear * (Tcore[i] - Trear[i])
            Q_from_perim_i = G_rec_rear * rear_weights[i] * (1.0 - f_core_rear) * (Tperim[i] - Trear[i])
            Q_to_cavity_i = G_rear_cavity * rear_weights[i] * (Trear[i] - Tcavity)

            Q_axial_rear = 0.0
            if i > 1
                Q_axial_rear += (G_rear_axial / dx) * (Trear[i-1] - Trear[i])
            end
            if i < nodes
                Q_axial_rear += (G_rear_axial / dx) * (Trear[i+1] - Trear[i])
            end

            Q_exit_coupling = i == nodes ? -Q_rear_to_tube : 0.0

            dTrear[i] = (Q_from_core_i + Q_from_perim_i - Q_to_cavity_i + Q_axial_rear + Q_exit_coupling) / C_rear_i
        end

        # 4. ALUMINA EXIT TUBE (Nodes 1:rear_nodes)
        A_tube_cs = pi * (R_tube_outer_v47^2 - R_tube_inner_v47^2)
        for j in 1:rear_nodes
            C_tube_j = rear_tube_capacity_v47(Ttube[j], z_rear[j], dx_rear)
            Q_gas_rear_j = Qgas_rear[j]
            Q_flange_j = rear_tube_flange_conductance_per_length_v47(z_rear[j], Ttube[j], p) * dx_rear * (Ttube[j] - WATER_FLANGE_TEMPERATURE_v47)
            Q_tube_cavity_j = G_felt_per_length_v47 * dx_rear * (Ttube[j] - Tcavity)

            Q_axial_tube = 0.0
            if j > 1
                k_prev = 0.5 * (k_al_f(Ttube[j]) + k_al_f(Ttube[j-1]))
                Q_axial_tube += (k_prev * A_tube_cs / dx_rear) * (Ttube[j-1] - Ttube[j])
            end
            if j < rear_nodes
                k_next = 0.5 * (k_al_f(Ttube[j]) + k_al_f(Ttube[j+1]))
                Q_axial_tube += (k_next * A_tube_cs / dx_rear) * (Ttube[j+1] - Ttube[j])
            end

            Q_entry_coupling = j == 1 ? Q_rear_to_tube : 0.0

            dTtube[j] = (Q_entry_coupling - Q_gas_rear_j - Q_flange_j - Q_tube_cavity_j + Q_axial_tube) / C_tube_j
        end

        # 5. EXTERNAL CAVITY SHELL
        Q_perim_to_cavity_total = sum(G_felt_per_length_v47 * dx * (Tperim[i] - Tcavity) for i in 1:nodes)
        Q_rear_to_cavity_total = sum(G_rear_cavity * rear_weights[i] * (Trear[i] - Tcavity) for i in 1:nodes)
        Q_tube_to_cavity_total = sum(G_felt_per_length_v47 * dx_rear * (Ttube[j] - Tcavity) for j in 1:rear_nodes)
        Q_cavity_to_amb = cavity_loss_to_ambient_v47(Tcavity, Tamb)

        du[cavity_index] = (Q_perim_to_cavity_total + Q_rear_to_cavity_total + Q_tube_to_cavity_total - Q_cavity_to_amb) / CAVITY_HEAT_CAPACITY_v47

        # 6. EXIT GAS THERMOCOUPLE SENSOR T3
        T_gas_sensor = Tg[nodes + 1] # Entrance of rear tube
        T_tube_sensor = Ttube[1]
        
        if flow_lpm > 1e-12
            h_tc = 40.0 * (flow_lpm / 15.0)^0.60
            Q_conv_sensor = h_tc * 1.0e-5 * (T_gas_sensor - T3_sensor)
            Q_rad_sensor = 0.85 * sigma_sb_v47 * 1.0e-5 * (T_tube_sensor^4 - T3_sensor^4)
            du[cavity_index + nodes + 1] = (Q_conv_sensor + Q_rad_sensor) / T3_SENSOR_CAPACITY_v47
        else
            # Zero-flow natural cooling
            h_tc_nat = 8.0
            Q_conv_sensor = h_tc_nat * 1.0e-5 * (T_tube_sensor - T3_sensor)
            Q_rad_sensor = 0.85 * sigma_sb_v47 * 1.0e-5 * (T_tube_sensor^4 - T3_sensor^4)
            du[cavity_index + nodes + 1] = (Q_conv_sensor + Q_rad_sensor) / T3_SENSOR_CAPACITY_v47
        end

        return nothing
    end
end

begin # simulation driver & solution extraction
    function simulate_v47(p, operating, time_history;
                          nodes=default_nodes_v47,
                          rear_nodes=REAR_TUBE_DEFAULT_NODES_v47,
                          initial_temperature=300.0,
                          initial_perim_temperature=300.0,
                          initial_rear_temperature=300.0,
                          initial_rear_reservoir_temperature=300.0,
                          initial_cavity_temperature=300.0,
                          solver=Rodas5P(autodiff=AutoFiniteDiff()),
                          reltol=1e-6, abstol=1e-7, dtmax=30.0)
        dx = L / nodes
        dx_rear = REAR_TUBE_LENGTH_v47 / rear_nodes
        z_solid = [(i - 0.5) * dx for i in 1:nodes]
        z_rear = [(j - 0.5) * dx_rear for j in 1:rear_nodes]
        z_gas = vcat([0.0], [i * dx for i in 1:nodes], [L + j * dx_rear for j in 1:rear_nodes])

        # Core solar weights
        core_weights = solar_weights_v5(p[11], p[4], nodes)
        rear_weights = rear_contact_weights_v47(z_solid)

        # Initial conditions vector
        Tcore0 = initial_temperature isa Function ? [initial_temperature(z) for z in z_solid] : fill(Float64(initial_temperature), nodes)
        Tperim0 = initial_perim_temperature isa Function ? [initial_perim_temperature(z) for z in z_solid] : fill(Float64(initial_perim_temperature), nodes)
        Ttube0 = initial_rear_temperature isa Function ? [initial_rear_temperature(zr) for zr in z_rear] : fill(Float64(initial_rear_temperature), rear_nodes)
        Trear0 = initial_rear_reservoir_temperature isa Function ? [initial_rear_reservoir_temperature(z) for z in z_solid] : fill(Float64(initial_rear_reservoir_temperature), nodes)
        Tcavity0 = Float64(initial_cavity_temperature)
        T3_0 = Ttube0[1]

        u0 = vcat(Tcore0, Tperim0, Ttube0, [Tcavity0], Trear0, [T3_0])

        # Context buffers
        context = (
            p=p, operating=operating, z=z_solid, dx=dx,
            core_weights=core_weights, rear_weights=rear_weights,
            z_rear=z_rear, dx_rear=dx_rear,
            Tg=zeros(nodes + rear_nodes + 1),
            Qgas=zeros(nodes),
            hcell=zeros(nodes),
            Qgas_rear=zeros(rear_nodes),
            hrear=zeros(rear_nodes),
            zg=z_gas,
        )

        f_ode = ODEFunction{true}(receiver_rhs_v47!)
        tspan = (time_history[1], time_history[end])
        prob = ODEProblem(f_ode, u0, tspan, context)

        sol = solve(prob, solver;
                    saveat=time_history,
                    reltol=reltol,
                    abstol=abstol,
                    dtmax=dtmax,
                    maxiters=1_000_000)

        # Extract matrices over time
        nt = length(time_history)
        core_temp = zeros(nodes, nt)
        perim_temp = zeros(nodes, nt)
        tube_temp = zeros(rear_nodes, nt)
        cavity_temp = zeros(nt)
        rear_temp = zeros(nodes, nt)
        T3_vec = zeros(nt)

        gas_temp = zeros(nodes + rear_nodes + 1, nt)
        htc_mat = zeros(nodes, nt)
        htc_rear_mat = zeros(rear_nodes, nt)
        Qgas_mat = zeros(nodes, nt)
        Qgas_rear_mat = zeros(rear_nodes, nt)

        cavity_index = 2 * nodes + rear_nodes + 1
        for (k, t_val) in enumerate(time_history)
            uk = sol.u[k]
            core_temp[:, k] .= uk[1:nodes]
            perim_temp[:, k] .= uk[nodes+1:2*nodes]
            tube_temp[:, k] .= uk[2*nodes+1:2*nodes+rear_nodes]
            cavity_temp[k] = uk[cavity_index]
            rear_temp[:, k] .= uk[cavity_index+1:cavity_index+nodes]
            T3_vec[k] = uk[cavity_index+nodes+1]

            flow = max(0.0, operating.flow_lpm(t_val))
            Tin = operating.inlet_temperature(t_val)
            
            gas_profile_v47!(
                context.Tg, context.Qgas, context.hcell, context.Qgas_rear, context.hrear,
                view(core_temp, :, k), view(perim_temp, :, k), view(tube_temp, :, k),
                Tin, flow, p, z_solid, dx, z_rear, dx_rear, z_gas
            )
            
            gas_temp[:, k] .= context.Tg
            htc_mat[:, k] .= context.hcell
            htc_rear_mat[:, k] .= context.hrear
            Qgas_mat[:, k] .= context.Qgas
            Qgas_rear_mat[:, k] .= context.Qgas_rear
        end

        return (
            time=time_history,
            z_solid=z_solid,
            z_gas=z_gas,
            z_rear=z_rear,
            core_temperature=core_temp,
            perim_temperature=perim_temp,
            rear_temperature=tube_temp,
            cavity_temperature=cavity_temp,
            rear_reservoir_temperature=rear_temp,
            T3_sensor=T3_vec,
            gas_temperature=gas_temp,
            heat_transfer_coefficient=htc_mat,
            rear_heat_transfer_coefficient=htc_rear_mat,
            receiver_gas_heat=Qgas_mat,
            rear_tube_gas_heat=Qgas_rear_mat,
            core_source_weights=core_weights,
            rear_contact_weights=rear_weights,
            ode_solution=sol,
        )
    end
end

begin # experimental case interface & loss function
    const WALL_CHAIN_POSITIONS_v47 = (
        T8=sensor_positions[:T8],   # 0.011 m (11 mm)
        T12=sensor_positions[:T9],  # 0.058 m (58 mm)
        T11=sensor_positions[:T10], # 0.107 m (107 mm)
    )

    function measured_initial_perim_profile_v47(data, simulation_id)
        id_str = string(simulation_id)
        z8 = WALL_CHAIN_POSITIONS_v47.T8
        z12 = WALL_CHAIN_POSITIONS_v47.T12
        z11 = WALL_CHAIN_POSITIONS_v47.T11
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

    function measured_initial_core_profile_v47(data, simulation_id)
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

    function measured_initial_rear_profile_v47(data, simulation_id, p)
        id_str = string(simulation_id)
        core_profile = measured_initial_core_profile_v47(data, id_str)
        perim_profile = measured_initial_perim_profile_v47(data, id_str)
        f_core_rear = rear_core_fraction_v47(p)
        T_receiver_exit = f_core_rear * core_profile(L) + (1.0 - f_core_rear) * perim_profile(L)
        T_tube_exit = observation(data, id_str, "_Tf")[1]
        return function (z_rear)
            fraction = clamp(z_rear / REAR_TUBE_LENGTH_v47, 0.0, 1.0)
            return T_receiver_exit + fraction * (T_tube_exit - T_receiver_exit)
        end
    end

    function measured_initial_rear_reservoir_profile_v47(data, simulation_id, p)
        core_profile = measured_initial_core_profile_v47(data, simulation_id)
        perim_profile = measured_initial_perim_profile_v47(data, simulation_id)
        f_core = rear_core_fraction_v47(p)
        return function (z)
            return f_core * core_profile(z) + (1.0 - f_core) * perim_profile(z)
        end
    end

    function experimental_case_v47(simulation_id; is_cooling=false)
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

    function core_at_v47(result, position)
        z = result.z_solid
        values = result.core_temperature
        position <= z[1] && return vec(values[1, :])
        position >= z[end] && return vec(values[end, :])
        i = searchsortedlast(z, position)
        fraction = (position - z[i]) / (z[i + 1] - z[i])
        return vec((1.0 - fraction) .* values[i, :] .+ fraction .* values[i + 1, :])
    end

    function perim_at_v47(result, position)
        z = result.z_solid
        values = result.perim_temperature
        position <= z[1] && return vec(values[1, :])
        position >= z[end] && return vec(values[end, :])
        i = searchsortedlast(z, position)
        fraction = (position - z[i]) / (z[i + 1] - z[i])
        return vec((1.0 - fraction) .* values[i, :] .+ fraction .* values[i + 1, :])
    end

    function extract_outputs_v47(result, p)
        T8 = perim_at_v47(result, sensor_positions[:T8])
        T12_perim = perim_at_v47(result, WALL_CHAIN_POSITIONS_v47.T12)
        T11_perim = perim_at_v47(result, WALL_CHAIN_POSITIONS_v47.T11)
        T9_core = core_at_v47(result, sensor_positions[:T9])
        T10_core = core_at_v47(result, sensor_positions[:T10])
        Tf = result.T3_sensor
        T2 = result.cavity_temperature
        h_mean = vec(mean(result.heat_transfer_coefficient, dims=1))
        return hcat(T8, T12_perim, T11_perim, T9_core, T10_core, Tf, T2, h_mean)
    end

    function align_cooling_initial_outputs_v47!(outputs, experiment)
        size(outputs, 1) >= 1 || return outputs
        outputs[1, 1:7] .= experiment[1, 1:7]
        return outputs
    end

    function solve_case_v47(p, simulation_id; is_cooling=false,
                           nodes=default_nodes,
                           rear_nodes=REAR_TUBE_DEFAULT_NODES_v47,
                           solver=Rodas5P(autodiff=AutoFiniteDiff()),
                           reltol=1e-6, abstol=1e-7, dtmax=30.0)
        id_str = string(simulation_id)
        conditions, time, experiment, data = experimental_case_v47(
            id_str; is_cooling=is_cooling
        )
        flow = observation(data, id_str, "_flow")
        Tin = observation(data, id_str, "_Tin")
        ambient = observation(data, id_str, "_Tamb")
        if is_cooling
            Tin = fill(COOLING_ROOM_TEMPERATURE_v47, length(time))
            ambient = fill(COOLING_ROOM_TEMPERATURE_v47, length(time))
        end
        T2 = observation(data, id_str, "_T2")
        irradiance = is_cooling ? zeros(length(time)) : fill(conditions[Io], length(time))
        operating = operating_condition_v3(
            irradiance=linear_history_v3(time, irradiance),
            flow_lpm=linear_history_v3(time, flow),
            inlet_temperature=linear_history_v3(time, Tin),
            ambient_temperature=linear_history_v3(time, ambient),
        )
        result = simulate_v47(
            p, operating, time;
            nodes=nodes,
            rear_nodes=rear_nodes,
            initial_temperature=measured_initial_core_profile_v47(data, id_str),
            initial_perim_temperature=measured_initial_perim_profile_v47(data, id_str),
            initial_rear_temperature=measured_initial_rear_profile_v47(data, id_str, p),
            initial_rear_reservoir_temperature=measured_initial_rear_reservoir_profile_v47(data, id_str, p),
            initial_cavity_temperature=T2[1],
            solver=solver,
            reltol=reltol,
            abstol=abstol,
            dtmax=dtmax,
        )
        outputs = extract_outputs_v47(result, p)
        is_cooling && align_cooling_initial_outputs_v47!(outputs, experiment)
        return outputs, result, experiment
    end
end

begin # objective loss & calibration functions
    function signal_loss_v47(time, model, exp_v)
        scale = max(maximum(exp_v) - minimum(exp_v), 20.0)
        norm_sq = mean(abs2, (model .- exp_v) ./ scale)
        dt = diff(time)
        grad_model = diff(model) ./ dt
        grad_exp = diff(exp_v) ./ dt
        grad_scale = max(maximum(abs, grad_exp), 0.05)
        grad_loss = mean(abs2, (grad_model .- grad_exp) ./ grad_scale)
        return norm_sq + 0.15 * grad_loss
    end

    function ordering_loss_v47(outputs, experiment)
        loss = 0.0
        T9_m = outputs[:, 4]
        T10_m = outputs[:, 5]
        loss += mean(max.(0.0, T10_m .- T9_m) .^ 2)
        
        T8_m = outputs[:, 1]
        T11_m = outputs[:, 3]
        loss += mean(max.(0.0, T11_m .- T8_m) .^ 2)
        
        return loss
    end

    function cooling_upturn_loss_v47(outputs)
        loss = 0.0
        for j in 1:7
            diffs = diff(outputs[:, j])
            loss += mean(max.(0.0, diffs) .^ 2)
        end
        return loss
    end

    function participating_total_heat_capacity_v47(p, nodes=15)
        C_core = rho_s_v47 * Cps_f(900.0) * A_solid_v47 * L
        C_perim = perimeter_heat_capacity_total_v47(p)
        C_rear = rear_reservoir_heat_capacity_v47(p)
        return C_core + C_perim + C_rear
    end

    function capacitance_regularization_v47(p; nodes=15)
        capacitance = participating_total_heat_capacity_v47(p, nodes)
        return CAPACITANCE_REGULARIZATION_WEIGHT_v47 *
               ((capacitance - MEASURED_ASSEMBLY_CAPACITANCE_v47) / MEASURED_ASSEMBLY_CAPACITANCE_v47)^2
    end
    
    function power_scale_regularization_v47(p)
        min_upper = min(p[5], p[6])
        if p[7] < min_upper - 0.1
            return 100.0 * (min_upper - 0.1 - p[7])^2
        end
        return 0.0
    end

    const keys_heating = sim_key_heat
    const keys_cooling = sim_key_cool

    function loss_cases_v47(p, keys; is_cooling=false, nodes=15)
        total_loss = 0.0
        data = is_cooling ? measurements_cooling : measurements
        for key in keys
            try
                outputs, _, experiment = solve_case_v47(
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
                    case_loss += signal_loss_v47(time, outputs[:, j], experiment[:, j])
                end
                case_loss /= 7.0
                if is_cooling
                    case_loss += COOLING_UPTURN_WEIGHT_v47 * cooling_upturn_loss_v47(outputs)
                else
                    case_loss += ORDERING_WEIGHT_v47 * ordering_loss_v47(outputs, experiment)
                end
                total_loss += case_loss
            catch e
                return 1e6
            end
        end
        return total_loss / length(keys)
    end

    function calibration_loss_v47(p; stage=:full, nodes=15)
        for i in eachindex(p)
            (p[i] < lb_full_v47[i] || p[i] > ub_full_v47[i]) && return 1e6 + 1e3 * (abs(p[i] - clamp(p[i], lb_full_v47[i], ub_full_v47[i])))
        end
        
        heating_loss = loss_cases_v47(p, keys_heating; is_cooling=false, nodes=nodes)
        reg = capacitance_regularization_v47(p; nodes=nodes) + power_scale_regularization_v47(p)
        return heating_loss + reg
    end
end

begin # exact energy conservation audit function (v47)
    function compute_energy_balance_v47(result, operating, time_val, p; nodes=default_nodes_v47, rear_nodes=REAR_TUBE_DEFAULT_NODES_v47)
        dx = L / nodes
        dx_rear = REAR_TUBE_LENGTH_v47 / rear_nodes
        irradiance = operating.irradiance(time_val)
        ambient = operating.ambient_temperature(time_val)
        flow = max(0.0, operating.flow_lpm(time_val))
        Tin = operating.inlet_temperature(time_val)

        Q_aperture = max(0.0, irradiance) * A_frt
        M = absorbed_power_scale_v47(irradiance, p)
        Q_delivered = M * Q_aperture
        chi = core_source_fraction_v47(p)
        Qcore_solar = chi * Q_delivered
        Qperim_solar = (1.0 - chi) * Q_delivered

        # Temperatures at this time
        idx = argmin(abs.(result.time .- time_val))
        Tcore = result.core_temperature[:, idx]
        Tperim = result.perim_temperature[:, idx]
        Ttube = result.rear_temperature[:, idx]
        Tcavity = result.cavity_temperature[idx]
        Tg = result.gas_temperature[:, idx]

        mdot = m_dot_standard_v47(flow)
        
        # Enthalpy transferred into gas
        h_suction = suction_heat_transfer_coefficient_v47(p)
        Q_suction = flow > 1e-12 ? h_suction * A_frt * (Tperim[1] - Tin) : 0.0
        Q_gas_core = sum(result.receiver_gas_heat[:, idx])
        Q_gas_rear = sum(result.rear_tube_gas_heat[:, idx])
        Q_gas_total = Q_suction + Q_gas_core + Q_gas_rear
        
        Q_front_rad = EPS_FRONT_FIXED_v47 * sigma_sb_v47 * A_frt * (Tperim[1]^4 - ambient^4)
        Q_cavity_amb = cavity_loss_to_ambient_v47(Tcavity, ambient)
        
        Q_flange = 0.0
        for j in 1:rear_nodes
            G_fl = rear_tube_flange_conductance_per_length_v47(result.z_rear[j], Ttube[j], p)
            Q_flange += G_fl * dx_rear * (Ttube[j] - WATER_FLANGE_TEMPERATURE_v47)
        end

        # Calculate instantaneous sensible storage rate dE_stored/dt for the solid network
        u_k = result.ode_solution.u[idx]
        du_k = similar(u_k)
        core_weights = solar_weights_v5(p[11], p[4], nodes)
        rear_weights = rear_contact_weights_v47(result.z_solid)
        context = (
            p=p, operating=operating, z=result.z_solid, dx=dx,
            core_weights=core_weights, rear_weights=rear_weights,
            z_rear=result.z_rear, dx_rear=dx_rear,
            Tg=zeros(nodes + rear_nodes + 1), Qgas=zeros(nodes), hcell=zeros(nodes),
            Qgas_rear=zeros(rear_nodes), hrear=zeros(rear_nodes), zg=result.z_gas,
        )
        receiver_rhs_v47!(
            du_k, u_k, context, time_val
        )
        
        dE_stored = 0.0
        cell_volume = A_solid_v47 * dx
        perim_cap_cell = perimeter_heat_capacity_total_v47(p) / nodes
        for i in 1:nodes
            core_cap = rho_s_v47 * Cps_f(Tcore[i]) * cell_volume
            dE_stored += core_cap * du_k[i]
            dE_stored += perim_cap_cell * du_k[nodes + i]
        end
        for j in 1:rear_nodes
            tube_cap = rear_tube_capacity_v47(Ttube[j], result.z_rear[j], dx_rear)
            dE_stored += tube_cap * du_k[2 * nodes + j]
        end
        cavity_index = 2 * nodes + rear_nodes + 1
        dE_stored += CAVITY_HEAT_CAPACITY_v47 * du_k[cavity_index]
        C_rear_tot = rear_reservoir_heat_capacity_v47(p)
        for i in 1:nodes
            dE_stored += (C_rear_tot * rear_weights[i]) * du_k[cavity_index + i]
        end

        instantaneous_residual = Q_delivered - (Q_gas_total + Q_front_rad + Q_cavity_amb + Q_flange + dE_stored)
        steady_flux_gap = Q_delivered - (Q_gas_total + Q_front_rad + Q_cavity_amb + Q_flange)

        return (
            delivered_W=Q_delivered,
            core_solar_W=Qcore_solar,
            perim_solar_W=Qperim_solar,
            suction_gas_heat_W=Q_suction,
            core_gas_heat_W=Q_gas_core,
            rear_tube_gas_heat_W=Q_gas_rear,
            gas_heat_W=Q_gas_total,
            front_rad_loss_W=Q_front_rad,
            cavity_amb_loss_W=Q_cavity_amb,
            flange_loss_W=Q_flange,
            dE_stored_dt_W=dE_stored,
            instantaneous_residual_W=instantaneous_residual,
            residual_W=instantaneous_residual,
            steady_flux_gap_W=steady_flux_gap,
        )
    end

    function fit_indices_for_stage_v47(stage)
        stage_symbol = Symbol(stage)
        stage_symbol == :rear && return fit_rear_stage_indices_v47
        stage_symbol == :transport && return fit_transport_stage_indices_v47
        stage_symbol == :full && return fit_full_stage_indices_v47
        stage_symbol == :all && return collect(1:23)
        throw(ArgumentError("unknown v47 calibration stage: $stage"))
    end

    function calibrate_v47(; nodes=15, maximum_iterations=500,
                           maximum_time_seconds=600.0, stage=:full, dataset=:heating)
        p0 = copy(pnew_v47)
        idx_fit = fit_indices_for_stage_v47(stage)
        function obj_free(theta)
            p_full = copy(p0)
            p_full[idx_fit] .= theta
            if dataset == :heating
                return loss_cases_v47(p_full, keys_heating; is_cooling=false, nodes=nodes) +
                       capacitance_regularization_v47(p_full; nodes=nodes) + power_scale_regularization_v47(p_full)
            elseif dataset == :cooling
                return loss_cases_v47(p_full, keys_cooling; is_cooling=true, nodes=nodes) +
                       capacitance_regularization_v47(p_full; nodes=nodes)
            else
                return loss_cases_v47(p_full, keys_heating; is_cooling=false, nodes=nodes) +
                       loss_cases_v47(p_full, keys_cooling; is_cooling=true, nodes=nodes) +
                       capacitance_regularization_v47(p_full; nodes=nodes) + power_scale_regularization_v47(p_full)
            end
        end
        res = optimize_with_nlopt_v3(
            obj_free,
            p0[idx_fit],
            lb_full_v47[idx_fit],
            ub_full_v47[idx_fit];
            maximum_iterations=maximum_iterations,
            maximum_time_seconds=maximum_time_seconds,
            algorithm=NLopt.LN_BOBYQA(),
            label="fit-v47-$(Symbol(stage))",
        )
        p_opt = copy(p0)
        p_opt[idx_fit] .= res.minimizer
        return (objective=res.minimum, parameters=p_opt, minimizer=res.minimizer, retcode=res.retcode)
    end
end
