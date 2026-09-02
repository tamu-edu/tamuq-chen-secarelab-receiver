# 1D_v50.jl - Entire Converter Model (ECM) for Monolithic Solar Receivers
# Version 50: audited repair of the v49 observation submodels, flow-slope objective,
# parameter bounds, and probe initialization. No new plant-physics mechanism is added.

begin # dependencies & base workspace inclusion
    using OrdinaryDiffEq
    using Statistics
    using LinearAlgebra
    using Optimization
    using OptimizationNLopt
    using DelimitedFiles

    if !@isdefined(observations) || isempty(observations)
        include("1D_v5.jl")
    end
end

begin # geometric & physical constants (authoritative v50 specifications)
    const N_channels_v50 = 100           # Number of square honeycomb channels
    const w_ch_v50 = 0.0015              # Channel opening width (1.50 mm)
    const w_cell_v50 = 0.019 / sqrt(N_channels_v50) # 0.0019 m (1.90 mm)
    const t_web_nominal_v50 = w_cell_v50 - w_ch_v50    # 0.0004 m (0.40 mm)
    const Dh_v50 = w_ch_v50              # Hydraulic diameter (1.50 mm)
    const L_v50 = 0.137                  # Monolith length (137 mm)
    const P_exchange_v50 = 4.0 * N_channels_v50 * w_ch_v50 # Internal wetted perimeter (0.60 m)
    
    # Frontal Aperture & Solid Porosity
    const A_frt_v50 = 3.6305e-4          # Front circular aperture area (pi * 0.010735^2 = 3.6305e-4 m^2)
    const R_core_v50 = sqrt(A_frt_v50 / pi) # 0.010735 m (10.74 mm)
    const a_v_v50 = P_exchange_v50 / A_frt_v50 # Specific volumetric exchange area (1652.66 m^2/m^3)
    const epsilon_v50 = (N_channels_v50 * w_ch_v50^2) / A_frt_v50 # Void porosity (0.6196)
    const A_solid_v50 = (1.0 - epsilon_v50) * A_frt_v50       # Solid matrix area (1.3809e-4 m^2)
    
    # Fundamental Physical Constants
    const sigma_sb_v50 = 5.670374419e-8  # Stefan-Boltzmann constant [W/m^2 K^4]
    const n_refr_v50 = 1.0               # Refractive index of air
    const default_nodes_v50 = 15         # Default axial nodes

    # Packaging & Housing Geometry
    const R_felt_v50 = 0.025             # Outer felt radius (25 mm)
    const R_case_v50 = 0.030             # Outer casing radius (30 mm)
    const L_case_v50 = L_v50             # Casing length (137 mm)
    const A_case_ext_v50 = 2.0 * pi * R_case_v50 * L_case_v50 # External casing area (0.02582 m^2)
    const k_felt_val_v50 = 0.06          # Alumina-silicate felt conductivity [W/m K]

    # Rear Hardware Assembly
    const REAR_TUBE_LENGTH_v50 = 0.063       # Rear exit tube length (63 mm)
    const REAR_TUBE_DEFAULT_NODES_v50 = 15   # Rear tube axial nodes
    const R_tube_inner_v50 = 0.007           # Rear tube gas radius (7 mm)
    const R_tube_outer_v50 = 0.011           # Rear tube outer radius (11 mm)

    # Thermal Capacitances & Sinks
    const CAVITY_HEAT_CAPACITY_v50 = 4000.0        # Outer cavity shell capacity [J/K]
    const T3_SENSOR_CAPACITY_v50 = 0.05            # T3 thermocouple junction capacity [J/K]
    const T3_PROBE_AREA_v50 = 5.0e-5               # T3 thermocouple effective surface area [m^2]
    const SLOPE_LOSS_WEIGHT_v50 = 1.0              # Flow slope derivative loss weight
    const WATER_FLANGE_TEMPERATURE_v50 = 293.15    # Chilled water flange boundary [K]
    const MEASURED_ASSEMBLY_CAPACITANCE_v50 = 301.0# Reference measured thermal capacitance [J/K]
    const CAPACITANCE_REGULARIZATION_WEIGHT_v50 = 5.0
    const COOLING_UPTURN_WEIGHT_v50 = 10.0
    const ORDERING_WEIGHT_v50 = 5.0
    const COOLING_ROOM_TEMPERATURE_v50 = 296.15

    # Radiation & Loss Emissivities
    const EPS_FRONT_FIXED_v50 = 0.85     # Front SiC / rim aperture emissivity
    const EPS_CASE_EXT_v50 = 0.80        # Outer shell casing emissivity
    const H_NAT_EXT_v50 = 5.0            # External natural convection to ambient [W/m^2 K]
    const NU_INFINITY_SQUARE_v50 = 3.61   # Fully developed laminar square duct Nu baseline
end

begin # material properties (unified polynomial functions)
    # SiC Solid Matrix Properties
    function ks_v50(T)
        Tc = T - 273.15
        return max(1.0, 51.0 * exp(-0.0030 * Tc) + 1.2)
    end

    function Cps_v50(T)
        Tc = T - 273.15
        return max(500.0, 900.0 + 0.30 * Tc - 3.0e5 * (T^-2))
    end
    const rho_s_v50 = 3100.0             # SiC density [kg/m^3]

    # Fluid (Air) Properties
    const p_atm_v50 = 101325.0
    const R_spec_air_v50 = 287.05
    rhof_v50(T) = p_atm_v50 / (R_spec_air_v50 * max(T, 100.0))
    cpf_v50(T) = 1005.0 + 0.05 * (T - 273.15)
    muf_v50(T) = 1.81e-5 * (max(T, 100.0) / 293.15)^0.70
    kf_v50(T) = 0.0257 * (max(T, 100.0) / 293.15)^0.80
    Prf_v50(T) = (cpf_v50(T) * muf_v50(T)) / kf_v50(T)

    # Standard air mass flow conversion: 1 LPM at standard conditions
    const rho_std_air_v50 = 101325.0 / (287.05 * 293.15) # 1.2041 kg/m^3
    m_dot_standard_v50(flow_lpm) = (flow_lpm / 60000.0) * rho_std_air_v50

    # Alumina Rear Tube Properties
    function k_al_v50(T)
        Tc = T - 273.15
        return 5.5 + 34.5 * exp(-0.0033 * Tc)
    end

    function cp_al_v50(T)
        return max(500.0, (1.00446 + 1.742e-4 * T - 2.796e4 * (T^-2)) * 1000.0)
    end
    const rho_al_v50 = 3900.0
end

begin # parameter vector definition & bounds (28 parameters)
    # V50 starts from the saved v49 vector. These values are seeds, not a v50 calibration.
    const pnew_v50 = [
        0.35508897745473206,     # 1: C1_Nu (Nusselt entry scale factor)
        0.013182183431474639,    # 2: C2_Nu (Nusselt entry denominator scale factor)
        3.6100,                  # 3: Nu_inf (fixed asymptotic fully developed square duct Nu)
        0.45289250852056845,     # 4: front_dep (front-face direct absorption fraction)
        1.3400,                  # 5: scale_456 (456 kW/m^2 cluster power scale)
        1.5800,                  # 6: scale_304 (304 kW/m^2 cluster power scale)
        1.1100,                  # 7: scale_256 (256 kW/m^2 cluster power scale)
        53.971081967474575,      # 8: G_core_perim (W/m K, core-to-perimeter radial conductance)
        85.27824863197442,       # 9: C_perim_eff (J/K, perimeter housing effective capacity)
        10.108991005714476,      # 10: k_perim_ref (W/m K, perimeter casing axial conductivity)
        274.49912776906086,      # 11: beta_opt (1/m, solar band optical extinction)
        0.7542695414237681,      # 12: chi (core solar absorption fraction)
        0.9743146639117137,      # 13: f_core_rear (rear split core fraction)
        5.124666398428759,       # 14: flange_scale (water flange conductance scale)
        0.9997602289810131,      # 15: k_core_axial_scale (core axial solid conduction scale)
        152.06272947224176,      # 16: C_rear_eff (J/K, rear hardware rail heat capacity)
        1.8532169665844016,      # 17: G_receiver_rear (W/K, receiver-to-rear rail contact conductance)
        6.828083484701864,       # 18: G_rear_tube (W/K, rear rail to alumina tube coupling)
        778.2455084200903,       # 19: beta_rad (1/m, thermal IR Rosseland extinction)
        0.11591364722106841,     # 20: kA_rear_eff (W m / K, rear rail axial conduction product)
        3.5463420889289515e-5,   # 21: delta_web (m, web thickness correction)
        0.003940910909642221,    # 22: F_LoS (front-to-rear line-of-sight radiation view factor)
        9.712355465108256,       # 23: h_suction (W/m^2 K, aperture suction preheating HTC)
        46.340318660513205,      # 24: h_probe_ref (W/m^2 K, convective HTC to T3 probe at 15 LPM)
        0.8242722722705073,      # 25: w_probe_rad (radiation enclosure factor of T3 probe to tube wall)
        0.022400324356015293,    # 26: G_probe_stem (W/K, probe stem conduction to mounting flange)
        0.047073136973866524,    # 27: w10_stem (sheath conduction weight for core sensor T10)
        0.0,                     # 28: w11_stem (sheath conduction weight for perimeter sensor T11)
    ]

    const lb_full_v50 = [
        0.01, 0.005, 3.61, 0.05, 1.0, 1.0, 0.8, 0.1, 50.0, 1.0,
        20.0, 0.5, 0.5, 0.05, 0.05, 20.0, 0.01, 0.01, 50.0, 0.001,
        1.0e-5, 0.0001, 0.0, 5.0, 0.05, 0.0001, 0.0, 0.0
    ]

    const ub_full_v50 = [
        0.50, 0.20, 3.61, 0.95, 2.0, 2.0, 2.0, 100.0, 300.0, 40.0,
        600.0, 1.0, 1.0, 20.0, 1.0, 250.0, 5.0, 15.0, 1000.0, 0.5,
        3.0e-4, 0.05, 150.0, 300.0, 1.0, 0.50, 0.50, 0.50
    ]

    const fit_full_stage_indices_v50 = vcat([1, 2], collect(4:28))
    const fit_plant_stage_indices_v50 = vcat([1, 2], collect(4:23))
    const fit_observation_stage_indices_v50 = collect(24:28)

    function fit_stage_indices_v50(stage)
        stage_name = lowercase(string(stage))
        stage_name == "observation" && return fit_observation_stage_indices_v50
        stage_name == "plant" && return fit_plant_stage_indices_v50
        stage_name == "full" && return fit_full_stage_indices_v50
        throw(ArgumentError("unknown v50 fit stage: $stage"))
    end
end

begin # parameter accessor helpers
    nu_c1_v50(p) = p[1]
    nu_c2_v50(p) = p[2]
    nu_inf_v50(p) = NU_INFINITY_SQUARE_v50
    front_deposition_fraction_v50(p) = clamp(p[4], 0.05, 0.95)
    core_perimeter_conductance_v50(p) = clamp(p[8], 0.1, 100.0)
    perimeter_heat_capacity_total_v50(p) = clamp(p[9], 50.0, 300.0)
    perimeter_conductivity_ref_v50(p) = clamp(p[10], 1.0, 40.0)
    core_optical_extinction_v50(p) = clamp(p[11], 20.0, 600.0)
    core_source_fraction_v50(p) = clamp(p[12], 0.5, 1.0)
    rear_core_fraction_v50(p) = clamp(p[13], 0.5, 1.0)
    flange_scale_v50(p) = clamp(p[14], 0.05, 20.0)
    core_axial_conduction_scale_v50(p) = clamp(p[15], 0.05, 1.0)
    rear_reservoir_heat_capacity_v50(p) = clamp(p[16], 20.0, 250.0)
    receiver_rear_conductance_v50(p) = clamp(p[17], 0.01, 5.0)
    rear_tube_conductance_v50(p) = clamp(p[18], 0.01, 15.0)
    core_thermal_rad_extinction_v50(p) = clamp(p[19], 50.0, 1000.0)
    rear_axial_kA_v50(p) = clamp(p[20], 0.001, 0.5)
    web_thickness_correction_v50(p) = clamp(p[21], 1.0e-5, 3.0e-4)
    line_of_sight_view_factor_v50(p) = clamp(p[22], 0.0001, 0.05)
    suction_heat_transfer_coefficient_v50(p) = clamp(p[23], 0.0, 150.0)
    probe_heat_transfer_coefficient_ref_v50(p) = clamp(p[24], 5.0, 300.0)
    probe_radiation_weight_v50(p) = clamp(p[25], 0.05, 1.0)
    probe_stem_conductance_v50(p) = clamp(p[26], 0.0001, 0.50)
    T10_rear_observation_weight_v50(p) = clamp(p[27], 0.0, 0.50)
    T11_stem_observation_weight_v50(p) = clamp(p[28], 0.0, 0.50)

    function absorbed_power_scale_v50(irradiance, p)
        irradiance <= 0.0 && return 0.0
        irradiance >= 400_000.0 && return p[5] # 456 kW/m^2 cluster
        irradiance >= 280_000.0 && return p[6] # 304 kW/m^2 cluster
        return p[7]                            # 256 kW/m^2 cluster
    end
end

begin # physical closures & correlations
    # Developing laminar Nusselt correlation (Shah-London single-group form) & series HTC
    function local_heat_transfer_coefficient_v50(Re, Pr, k_fluid, z_pos, T_solid, p)
        C1 = nu_c1_v50(p)
        C2 = nu_c2_v50(p)
        Nu_inf = NU_INFINITY_SQUARE_v50
        delta_web = web_thickness_correction_v50(p)
        
        # Entrance coordinate
        z_star = max(z_pos, 1.0e-4)
        Gz = (Dh_v50 / z_star) * Re * Pr
        
        # Developing Nu (Shah-London entry without compound singularities)
        Nu_entry = (C1 * Gz) / (1.0 + C2 * (Gz^(2.0 / 3.0)))
        Nu_fluid = Nu_inf + Nu_entry
        
        h_fluid = Nu_fluid * k_fluid / Dh_v50
        
        # Intra-strut solid conduction resistance in series
        t_eff = max(1.0e-5, t_web_nominal_v50 + delta_web)
        k_s_val = ks_v50(T_solid)
        R_solid = t_eff / (4.0 * k_s_val)
        
        h_eff = 1.0 / (1.0 / h_fluid + R_solid)
        return h_eff, Nu_fluid
    end

    # Core axial effective conductivity (solid conduction + decoupled Rosseland radiation)
    function core_effective_axial_conductivity_v50(T, p)
        k_scale = core_axial_conduction_scale_v50(p)
        beta_rad = core_thermal_rad_extinction_v50(p)
        k_cond = k_scale * (1.0 - epsilon_v50) * ks_v50(T)
        k_rad = (16.0 * sigma_sb_v50 * (n_refr_v50^2) * (T^3)) / (3.0 * beta_rad)
        return k_cond + k_rad
    end

    # Perimeter casing axial conductivity
    function perimeter_axial_conductivity_v50(T, p)
        k_ref = perimeter_conductivity_ref_v50(p)
        return k_ref * (1.0 + 0.0005 * (T - 293.15))
    end

    # Perimeter to external cavity conductance
    const G_felt_per_length_v50 = (2.0 * pi * k_felt_val_v50) / log(R_felt_v50 / R_core_v50)

    # Cavity to ambient loss
    function cavity_loss_to_ambient_v50(T_cavity, T_amb)
        Q_conv = H_NAT_EXT_v50 * A_case_ext_v50 * (T_cavity - T_amb)
        Q_rad = EPS_CASE_EXT_v50 * sigma_sb_v50 * A_case_ext_v50 * (T_cavity^4 - T_amb^4)
        return Q_conv + Q_rad
    end

    # Rear tube flange conductance per unit length
    function rear_tube_flange_conductance_per_length_v50(z_rear, T_tube, p)
        s_flange = flange_scale_v50(p)
        fraction = clamp(z_rear / REAR_TUBE_LENGTH_v50, 0.0, 1.0)
        # Activation towards water flange (z_rear -> L_tube)
        sigmoid = 1.0 / (1.0 + exp(-30.0 * (fraction - 0.70)))
        G_al = (2.0 * pi * k_al_v50(T_tube)) / log(R_tube_outer_v50 / R_tube_inner_v50)
        return s_flange * G_al * sigmoid
    end

    # Rear tube sensible heat capacity
    function rear_tube_capacity_v50(T_tube, z_rear, dx_rear)
        V_al = pi * (R_tube_outer_v50^2 - R_tube_inner_v50^2) * dx_rear
        return rho_al_v50 * cp_al_v50(T_tube) * V_al
    end

    # Distributed rear contact weights
    function rear_contact_weights_v50(z_solid)
        weights = [exp(12.0 * (z / L_v50 - 1.0)) for z in z_solid]
        return weights ./ sum(weights)
    end

    # Solar absorption spatial distribution weights
    function solar_weights_v50(beta_opt, f_front, nodes)
        dx = L_v50 / nodes
        weights = zeros(nodes)
        for i in 1:nodes
            z_left = (i - 1) * dx
            z_right = i * dx
            # Beer-Lambert volume integral
            w_beer = exp(-beta_opt * z_left) - exp(-beta_opt * z_right)
            weights[i] = (1.0 - f_front) * w_beer
        end
        weights[1] += f_front
        return weights ./ sum(weights)
    end
end

begin # gas profile marching (100% mass flow & zero-flow physics)
    function gas_profile_v50!(Tg, Qgas, hcell, Qgas_rear, hrear,
                             Tcore, Tperim, Ttube, Tin, flow_lpm, p,
                             z_solid, dx, z_rear, dx_rear, zg)
        nodes = length(Tcore)
        rear_nodes = length(Ttube)
        flow = max(0.0, flow_lpm)
        mdot = m_dot_standard_v50(flow)
        h_suction = suction_heat_transfer_coefficient_v50(p)

        if flow > 1e-12
            # 1. Suction Preheating at Aperture Face (z = 0)
            Q_suction = h_suction * A_frt_v50 * (Tperim[1] - Tin)
            Tg[1] = Tin + Q_suction / (mdot * cpf_v50(Tin))

            # 2. Marching through Core Honeycomb Channels (z in [0, L])
            for i in 1:nodes
                T_solid = Tcore[i]
                T_g_in = Tg[i]
                T_film = 0.5 * (T_solid + T_g_in)
                
                cp_f = cpf_v50(T_film)
                k_f = kf_v50(T_film)
                mu_f = muf_v50(T_film)
                Pr_f = (cp_f * mu_f) / k_f
                
                v_ch = mdot / (N_channels_v50 * rhof_v50(T_film) * (w_ch_v50^2))
                Re_ch = (rhof_v50(T_film) * v_ch * Dh_v50) / mu_f
                
                h_eff, _ = local_heat_transfer_coefficient_v50(Re_ch, Pr_f, k_f, z_solid[i], T_solid, p)
                hcell[i] = h_eff
                
                NTU = (h_eff * P_exchange_v50 * dx) / (mdot * cp_f)
                T_g_out = T_solid - (T_solid - T_g_in) * exp(-NTU)
                
                Tg[i + 1] = T_g_out
                Qgas[i] = mdot * cp_f * (T_g_out - T_g_in)
            end

            # 3. Marching through Alumina Exit Tube (z in [L, L + L_tube])
            P_tube = 2.0 * pi * R_tube_inner_v50
            for j in 1:rear_nodes
                T_tube_j = Ttube[j]
                T_g_in = Tg[nodes + j]
                T_film = 0.5 * (T_tube_j + T_g_in)
                
                cp_f = cpf_v50(T_film)
                k_f = kf_v50(T_film)
                mu_f = muf_v50(T_film)
                Pr_f = (cp_f * mu_f) / k_f
                
                v_tube = mdot / (pi * (R_tube_inner_v50^2) * rhof_v50(T_film))
                Re_tube = (rhof_v50(T_film) * v_tube * (2.0 * R_tube_inner_v50)) / mu_f
                
                # Fully developed / developing pipe flow Nu
                Nu_tube = 4.36 + (0.0668 * (2.0 * R_tube_inner_v50 / max(z_rear[j], 1e-4)) * Re_tube * Pr_f) /
                                 (1.0 + 0.04 * ((2.0 * R_tube_inner_v50 / max(z_rear[j], 1e-4)) * Re_tube * Pr_f)^(2.0 / 3.0))
                h_tube = Nu_tube * k_f / (2.0 * R_tube_inner_v50)
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
    function receiver_rhs_v50!(du, u, p_tuple, t)
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
        Q_aperture = max(0.0, irradiance) * A_frt_v50
        M = absorbed_power_scale_v50(irradiance, p)
        Q_delivered = M * Q_aperture
        chi = core_source_fraction_v50(p)
        Qcore_solar = chi * Q_delivered
        Qperim_solar = (1.0 - chi) * Q_delivered

        # March gas temperatures
        gas_profile_v50!(Tg, Qgas, hcell, Qgas_rear, hrear,
                         Tcore, Tperim, Ttube, Tin, flow_lpm, p,
                         z_solid, dx, z_rear, dx_rear, zg)

        # Line-of-sight cavity direct radiation exchange (front node 1 to rear node N)
        F_LoS = line_of_sight_view_factor_v50(p)
        Q_rad_LoS = F_LoS * A_frt_v50 * sigma_sb_v50 * (Tcore[1]^4 - Tcore[nodes]^4)

        # 1. CORE SOLID MATRIX (Nodes 1:nodes)
        G_cp = core_perimeter_conductance_v50(p)
        G_rec_rear = receiver_rear_conductance_v50(p)
        f_core_rear = rear_core_fraction_v50(p)
        cell_vol = A_solid_v50 * dx

        for i in 1:nodes
            C_core_i = rho_s_v50 * Cps_v50(Tcore[i]) * cell_vol
            Q_solar_i = Qcore_solar * core_weights[i]
            Q_gas_i = Qgas[i]
            Q_radial_i = G_cp * dx * (Tcore[i] - Tperim[i])
            Q_rear_i = G_rec_rear * rear_weights[i] * f_core_rear * (Tcore[i] - Trear[i])

            # Axial conduction
            Q_axial = 0.0
            if i > 1
                k_prev = 0.5 * (core_effective_axial_conductivity_v50(Tcore[i], p) + core_effective_axial_conductivity_v50(Tcore[i-1], p))
                Q_axial += (k_prev * A_solid_v50 / dx) * (Tcore[i-1] - Tcore[i])
            end
            if i < nodes
                k_next = 0.5 * (core_effective_axial_conductivity_v50(Tcore[i], p) + core_effective_axial_conductivity_v50(Tcore[i+1], p))
                Q_axial += (k_next * A_solid_v50 / dx) * (Tcore[i+1] - Tcore[i])
            end

            # Line-of-sight radiation boundary term
            Q_LoS_i = i == 1 ? -Q_rad_LoS : (i == nodes ? Q_rad_LoS : 0.0)

            dTcore[i] = (Q_solar_i - Q_gas_i - Q_radial_i - Q_rear_i + Q_axial + Q_LoS_i) / C_core_i
        end

        # 2. PERIMETER HOUSING (Nodes 1:nodes)
        C_perim_cell = perimeter_heat_capacity_total_v50(p) / nodes
        A_perim_cs = pi * (R_case_v50^2 - R_felt_v50^2) # Aluminum shell conduction area
        h_suction = suction_heat_transfer_coefficient_v50(p)
        mdot = m_dot_standard_v50(flow_lpm)

        Q_suction = flow_lpm > 1e-12 ? h_suction * A_frt_v50 * (Tperim[1] - Tin) : 0.0
        Q_front_rad = EPS_FRONT_FIXED_v50 * sigma_sb_v50 * A_frt_v50 * (Tperim[1]^4 - Tamb^4)

        for i in 1:nodes
            Q_solar_perim_i = i == 1 ? Qperim_solar : 0.0
            Q_radial_i = G_cp * dx * (Tcore[i] - Tperim[i])
            Q_cavity_i = G_felt_per_length_v50 * dx * (Tperim[i] - Tcavity)
            Q_rear_perim_i = G_rec_rear * rear_weights[i] * (1.0 - f_core_rear) * (Tperim[i] - Trear[i])

            Q_axial_perim = 0.0
            if i > 1
                k_prev = 0.5 * (perimeter_axial_conductivity_v50(Tperim[i], p) + perimeter_axial_conductivity_v50(Tperim[i-1], p))
                Q_axial_perim += (k_prev * A_perim_cs / dx) * (Tperim[i-1] - Tperim[i])
            end
            if i < nodes
                k_next = 0.5 * (perimeter_axial_conductivity_v50(Tperim[i], p) + perimeter_axial_conductivity_v50(Tperim[i+1], p))
                Q_axial_perim += (k_next * A_perim_cs / dx) * (Tperim[i+1] - Tperim[i])
            end

            Q_boundary = i == 1 ? -(Q_suction + Q_front_rad) : 0.0

            dTperim[i] = (Q_solar_perim_i + Q_radial_i - Q_cavity_i - Q_rear_perim_i + Q_axial_perim + Q_boundary) / C_perim_cell
        end

        # 3. DISTRIBUTED REAR RAIL (Nodes 1:nodes)
        C_rear_total = rear_reservoir_heat_capacity_v50(p)
        G_rear_tube = rear_tube_conductance_v50(p)
        kA_rear = rear_axial_kA_v50(p)
        Q_rear_to_tube = G_rear_tube * (Trear[nodes] - Ttube[1])

        for i in 1:nodes
            C_rear_i = C_rear_total * rear_weights[i]
            Q_from_core_i = G_rec_rear * rear_weights[i] * f_core_rear * (Tcore[i] - Trear[i])
            Q_from_perim_i = G_rec_rear * rear_weights[i] * (1.0 - f_core_rear) * (Tperim[i] - Trear[i])

            # Dimensionally consistent rear axial conduction: (kA)/dx * Delta T
            Q_axial_rear = 0.0
            if i > 1
                Q_axial_rear += (kA_rear / dx) * (Trear[i-1] - Trear[i])
            end
            if i < nodes
                Q_axial_rear += (kA_rear / dx) * (Trear[i+1] - Trear[i])
            end

            Q_exit_coupling = i == nodes ? -Q_rear_to_tube : 0.0

            dTrear[i] = (Q_from_core_i + Q_from_perim_i + Q_axial_rear + Q_exit_coupling) / C_rear_i
        end

        # 4. ALUMINA EXIT TUBE (Nodes 1:rear_nodes)
        A_tube_cs = pi * (R_tube_outer_v50^2 - R_tube_inner_v50^2)
        for j in 1:rear_nodes
            C_tube_j = rear_tube_capacity_v50(Ttube[j], z_rear[j], dx_rear)
            Q_gas_rear_j = Qgas_rear[j]
            Q_flange_j = rear_tube_flange_conductance_per_length_v50(z_rear[j], Ttube[j], p) * dx_rear * (Ttube[j] - WATER_FLANGE_TEMPERATURE_v50)
            Q_tube_cavity_j = G_felt_per_length_v50 * dx_rear * (Ttube[j] - Tcavity)

            Q_axial_tube = 0.0
            if j > 1
                k_prev = 0.5 * (k_al_v50(Ttube[j]) + k_al_v50(Ttube[j-1]))
                Q_axial_tube += (k_prev * A_tube_cs / dx_rear) * (Ttube[j-1] - Ttube[j])
            end
            if j < rear_nodes
                k_next = 0.5 * (k_al_v50(Ttube[j]) + k_al_v50(Ttube[j+1]))
                Q_axial_tube += (k_next * A_tube_cs / dx_rear) * (Ttube[j+1] - Ttube[j])
            end

            Q_entry_coupling = j == 1 ? Q_rear_to_tube : 0.0

            dTtube[j] = (Q_entry_coupling - Q_gas_rear_j - Q_flange_j - Q_tube_cavity_j + Q_axial_tube) / C_tube_j
        end

        # 5. EXTERNAL CAVITY SHELL
        Q_perim_to_cavity_total = sum(G_felt_per_length_v50 * dx * (Tperim[i] - Tcavity) for i in 1:nodes)
        Q_tube_to_cavity_total = sum(G_felt_per_length_v50 * dx_rear * (Ttube[j] - Tcavity) for j in 1:rear_nodes)
        Q_cavity_to_amb = cavity_loss_to_ambient_v50(Tcavity, Tamb)

        du[cavity_index] = (Q_perim_to_cavity_total + Q_tube_to_cavity_total - Q_cavity_to_amb) / CAVITY_HEAT_CAPACITY_v50

        # 6. EXIT GAS THERMOCOUPLE SENSOR T3.
        # This is a passive one-way observation model: it does not feed energy
        # back into the plant gas/tube states.
        T_gas_sensor = Tg[nodes + 1] # Entrance of rear tube
        T_tube_sensor = Ttube[1]

        h_probe_ref = probe_heat_transfer_coefficient_ref_v50(p)
        radiation_weight = probe_radiation_weight_v50(p)
        G_probe_stem = probe_stem_conductance_v50(p)

        if flow_lpm > 1e-12
            h_tc = h_probe_ref * (flow_lpm / 15.0)^0.60
            Q_conv_sensor = h_tc * T3_PROBE_AREA_v50 * (T_gas_sensor - T3_sensor)
        else
            # Zero-flow natural cooling
            h_tc_nat = 4.0
            Q_conv_sensor = h_tc_nat * T3_PROBE_AREA_v50 * (T_tube_sensor - T3_sensor)
        end
        Q_rad_sensor = radiation_weight * sigma_sb_v50 * T3_PROBE_AREA_v50 *
                       (T_tube_sensor^4 - T3_sensor^4)
        Q_stem_sensor = G_probe_stem * (WATER_FLANGE_TEMPERATURE_v50 - T3_sensor)
        du[cavity_index + nodes + 1] =
            (Q_conv_sensor + Q_rad_sensor + Q_stem_sensor) / T3_SENSOR_CAPACITY_v50

        return nothing
    end
end

begin # simulation driver & solution extraction
    function simulate_v50(p, operating, time_history;
                          nodes=default_nodes_v50,
                          rear_nodes=REAR_TUBE_DEFAULT_NODES_v50,
                          initial_temperature=300.0,
                          initial_perim_temperature=300.0,
                          initial_rear_temperature=300.0,
                          initial_rear_reservoir_temperature=300.0,
                          initial_cavity_temperature=300.0,
                          initial_probe_temperature=nothing,
                          solver=Rodas5P(autodiff=AutoFiniteDiff()),
                          reltol=1e-6, abstol=1e-7, dtmax=30.0)
        dx = L_v50 / nodes
        dx_rear = REAR_TUBE_LENGTH_v50 / rear_nodes
        z_solid = [(i - 0.5) * dx for i in 1:nodes]
        z_rear = [(j - 0.5) * dx_rear for j in 1:rear_nodes]
        z_gas = vcat([0.0], [i * dx for i in 1:nodes], [L_v50 + j * dx_rear for j in 1:rear_nodes])

        # Core solar weights
        core_weights = solar_weights_v50(core_optical_extinction_v50(p), front_deposition_fraction_v50(p), nodes)
        rear_weights = rear_contact_weights_v50(z_solid)

        # Initial conditions vector
        Tcore0 = initial_temperature isa Function ? [initial_temperature(z) for z in z_solid] : fill(Float64(initial_temperature), nodes)
        Tperim0 = initial_perim_temperature isa Function ? [initial_perim_temperature(z) for z in z_solid] : fill(Float64(initial_perim_temperature), nodes)
        Ttube0 = initial_rear_temperature isa Function ? [initial_rear_temperature(zr) for zr in z_rear] : fill(Float64(initial_rear_temperature), rear_nodes)
        Trear0 = initial_rear_reservoir_temperature isa Function ? [initial_rear_reservoir_temperature(z) for z in z_solid] : fill(Float64(initial_rear_reservoir_temperature), nodes)
        Tcavity0 = Float64(initial_cavity_temperature)
        T3_0 = initial_probe_temperature === nothing ? Ttube0[1] : Float64(initial_probe_temperature)

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

        f_ode = ODEFunction{true}(receiver_rhs_v50!)
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
            
            gas_profile_v50!(
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
    const WALL_CHAIN_POSITIONS_v50 = (
        T8=sensor_positions[:T8],   # 0.011 m (11 mm)
        T12=sensor_positions[:T9],  # 0.058 m (58 mm)
        T11=sensor_positions[:T10], # 0.107 m (107 mm)
    )

    function measured_initial_perim_profile_v50(data, simulation_id)
        id_str = string(simulation_id)
        z8 = WALL_CHAIN_POSITIONS_v50.T8
        z12 = WALL_CHAIN_POSITIONS_v50.T12
        z11 = WALL_CHAIN_POSITIONS_v50.T11
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

    function measured_initial_core_profile_v50(data, simulation_id)
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

    function measured_initial_rear_profile_v50(data, simulation_id, p)
        id_str = string(simulation_id)
        core_profile = measured_initial_core_profile_v50(data, id_str)
        perim_profile = measured_initial_perim_profile_v50(data, id_str)
        f_core_rear = rear_core_fraction_v50(p)
        T_receiver_exit = f_core_rear * core_profile(L_v50) + (1.0 - f_core_rear) * perim_profile(L_v50)
        T_tube_exit = observation(data, id_str, "_Tf")[1]
        return function (z_rear)
            fraction = clamp(z_rear / REAR_TUBE_LENGTH_v50, 0.0, 1.0)
            return T_receiver_exit + fraction * (T_tube_exit - T_receiver_exit)
        end
    end

    function measured_initial_rear_reservoir_profile_v50(data, simulation_id, p)
        id_str = string(simulation_id)
        perim_profile = measured_initial_perim_profile_v50(data, id_str)
        core_profile = measured_initial_core_profile_v50(data, id_str)
        f_core_rear = rear_core_fraction_v50(p)
        return function (z)
            return f_core_rear * core_profile(z) + (1.0 - f_core_rear) * perim_profile(z)
        end
    end

    function experimental_case_v50(simulation_id; is_cooling=false)
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

    function perim_at_v50(result, z_target)
        idx = argmin(abs.(result.z_solid .- z_target))
        return vec(result.perim_temperature[idx, :])
    end

    function core_at_v50(result, z_target)
        idx = argmin(abs.(result.z_solid .- z_target))
        return vec(result.core_temperature[idx, :])
    end

    function rear_reservoir_at_v50(result, z_target)
        idx = argmin(abs.(result.z_solid .- z_target))
        return vec(result.rear_reservoir_temperature[idx, :])
    end

    function extract_outputs_v50(result, p)
        T8 = perim_at_v50(result, sensor_positions[:T8])
        T12_perim = perim_at_v50(result, WALL_CHAIN_POSITIONS_v50.T12)
        T11_local = perim_at_v50(result, WALL_CHAIN_POSITIONS_v50.T11)
        w11 = T11_stem_observation_weight_v50(p)
        T11_perim = (1.0 - w11) .* T11_local .+ w11 .* WATER_FLANGE_TEMPERATURE_v50
        T9_core = core_at_v50(result, sensor_positions[:T9])
        T10_local = core_at_v50(result, sensor_positions[:T10])
        T10_rear = rear_reservoir_at_v50(result, sensor_positions[:T10])
        w10 = T10_rear_observation_weight_v50(p)
        T10_core = (1.0 - w10) .* T10_local .+ w10 .* T10_rear
        Tf = result.T3_sensor
        T2 = result.cavity_temperature
        h_mean = vec(mean(result.heat_transfer_coefficient, dims=1))
        return hcat(T8, T12_perim, T11_perim, T9_core, T10_core, Tf, T2, h_mean)
    end

    function solve_case_v50(p, simulation_id; is_cooling=false,
                           nodes=default_nodes_v50,
                           rear_nodes=REAR_TUBE_DEFAULT_NODES_v50,
                           solver=Rodas5P(autodiff=AutoFiniteDiff()),
                           reltol=1e-6, abstol=1e-7, dtmax=30.0)
        id_str = string(simulation_id)
        conditions, time, experiment, data = experimental_case_v50(
            id_str; is_cooling=is_cooling
        )
        flow = observation(data, id_str, "_flow")
        Tin = observation(data, id_str, "_Tin")
        ambient = observation(data, id_str, "_Tamb")
        if is_cooling
            Tin = fill(COOLING_ROOM_TEMPERATURE_v50, length(time))
            ambient = fill(COOLING_ROOM_TEMPERATURE_v50, length(time))
        end
        T2 = observation(data, id_str, "_T2")
        irradiance = is_cooling ? zeros(length(time)) : fill(conditions[Io], length(time))
        operating = operating_condition_v3(
            irradiance=linear_history_v3(time, irradiance),
            flow_lpm=linear_history_v3(time, flow),
            inlet_temperature=linear_history_v3(time, Tin),
            ambient_temperature=linear_history_v3(time, ambient),
        )
        result = simulate_v50(
            p, operating, time;
            nodes=nodes,
            rear_nodes=rear_nodes,
            initial_temperature=measured_initial_core_profile_v50(data, id_str),
            initial_perim_temperature=measured_initial_perim_profile_v50(data, id_str),
            initial_rear_temperature=measured_initial_rear_profile_v50(data, id_str, p),
            initial_rear_reservoir_temperature=measured_initial_rear_reservoir_profile_v50(data, id_str, p),
            initial_cavity_temperature=T2[1],
            initial_probe_temperature=observation(data, id_str, "_Tf")[1],
            solver=solver,
            reltol=reltol,
            abstol=abstol,
            dtmax=dtmax,
        )
        outputs = extract_outputs_v50(result, p)
        return outputs, result, experiment
    end
end

begin # objective loss & calibration functions
    function signal_loss_v50(time, model, exp_v)
        scale = max(maximum(exp_v) - minimum(exp_v), 20.0)
        norm_sq = mean(abs2, (model .- exp_v) ./ scale)
        dt = diff(time)
        grad_model = diff(model) ./ dt
        grad_exp = diff(exp_v) ./ dt
        grad_scale = max(maximum(abs, grad_exp), 0.05)
        grad_loss = mean(abs2, (grad_model .- grad_exp) ./ grad_scale)
        return norm_sq + 0.15 * grad_loss
    end

    function ordering_loss_v50(outputs, experiment)
        loss = 0.0
        T9_m = outputs[:, 4]
        T10_m = outputs[:, 5]
        loss += mean(max.(0.0, (T10_m .- T9_m) ./ 50.0) .^ 2)
        
        T8_m = outputs[:, 1]
        T11_m = outputs[:, 3]
        loss += mean(max.(0.0, (T11_m .- T8_m) ./ 50.0) .^ 2)
        
        return loss
    end

    function cooling_upturn_loss_v50(outputs)
        loss = 0.0
        for j in 1:7
            diffs = diff(outputs[:, j])
            loss += mean(max.(0.0, diffs) .^ 2)
        end
        return loss
    end

    function participating_total_heat_capacity_v50(p, nodes=15; kwargs...)
        C_core = rho_s_v50 * Cps_v50(900.0) * A_solid_v50 * L_v50
        C_perim = perimeter_heat_capacity_total_v50(p)
        C_rear = rear_reservoir_heat_capacity_v50(p)
        return C_core + C_perim + C_rear
    end
    function participating_total_heat_capacity_v50(p; nodes=15, kwargs...)
        return participating_total_heat_capacity_v50(p, nodes)
    end
    const total_participating_heat_capacity_v50 = participating_total_heat_capacity_v50

    function capacitance_regularization_v50(p; nodes=15)
        capacitance = participating_total_heat_capacity_v50(p, nodes)
        return CAPACITANCE_REGULARIZATION_WEIGHT_v50 *
               ((capacitance - MEASURED_ASSEMBLY_CAPACITANCE_v50) / MEASURED_ASSEMBLY_CAPACITANCE_v50)^2
    end
    
    function power_scale_regularization_v50(p)
        min_upper = min(p[5], p[6])
        if p[7] < min_upper - 0.1
            return 100.0 * (min_upper - 0.1 - p[7])^2
        end
        return 0.0
    end

    const keys_heating_v50 = sim_key_heat
    const keys_cooling_v50 = sim_key_cool

    # Explicit flow-derivative penalty across irradiance clusters.
    const CLUSTER_IRR_MAP_v50 = Dict(
        456000.0 => ["E67", "E68", "E69", "E70", "E71"],
        304000.0 => ["E72", "E73", "E74", "E75", "E76"],
        256000.0 => ["E77", "E78", "E79", "E80", "E81"]
    )
    const SLOPE_NORM_SIGMAS_v50 = [5.0, 3.0, 3.0, 3.0, 3.0, 3.0, 1.0] # K/LPM scale per sensor

    function loss_flow_slopes_v50(case_outputs_dict)
        total_slope_loss = 0.0
        for (irr, case_keys) in CLUSTER_IRR_MAP_v50
            flows = Float64[]
            mod_temps = [Float64[] for _ in 1:7]
            exp_temps = [Float64[] for _ in 1:7]

            for key in case_keys
                haskey(case_outputs_dict, key) || continue
                outputs, experiment, _ = case_outputs_dict[key]
                cond_dict = simulation_conditions[key]
                flow_val = cond_dict[qlpm]
                push!(flows, flow_val)
                for s in 1:7
                    push!(mod_temps[s], outputs[end, s])
                    push!(exp_temps[s], experiment[end, s])
                end
            end

            length(flows) >= 3 || continue
            f_mean = mean(flows)
            f_var = sum((flows .- f_mean) .^ 2)
            f_var > 1e-6 || continue

            for s in [1, 3, 4, 5, 6, 7] # T8, T11, T9, T10, T3, T2
                slope_mod = sum((flows .- f_mean) .* (mod_temps[s] .- mean(mod_temps[s]))) / f_var
                slope_exp = sum((flows .- f_mean) .* (exp_temps[s] .- mean(exp_temps[s]))) / f_var
                sigma_s = SLOPE_NORM_SIGMAS_v50[s]
                slope_diff = (slope_mod - slope_exp) / sigma_s
                total_slope_loss += slope_diff^2
            end
        end
        return SLOPE_LOSS_WEIGHT_v50 * (total_slope_loss / 18.0)
    end

    function dataset_loss_components_v50(p, keys; is_cooling=false, nodes=15)
        total_loss = 0.0
        data = is_cooling ? measurements_cooling : measurements
        case_dict = Dict{String, Tuple{Matrix{Float64}, Matrix{Float64}, Vector{Float64}}}()
        for key in keys
            try
                key_str = string(key)
                outputs, _, experiment = solve_case_v50(
                    p, key_str;
                    is_cooling=is_cooling,
                    nodes=nodes,
                    solver=Rodas5P(autodiff=AutoFiniteDiff()),
                    reltol=1e-4,
                    abstol=1e-5,
                    dtmax=60.0,
                )
                case_loss = 0.0
                time = Float64.(observation_time(data, key_str))
                for j in 1:7
                    case_loss += signal_loss_v50(time, outputs[:, j], experiment[:, j])
                end
                case_loss /= 7.0
                if is_cooling
                    case_loss += COOLING_UPTURN_WEIGHT_v50 * cooling_upturn_loss_v50(outputs)
                else
                    case_loss += ORDERING_WEIGHT_v50 * ordering_loss_v50(outputs, experiment)
                end
                total_loss += case_loss
                case_dict[key_str] = (outputs, experiment, time)
            catch e
                return (signal_loss=1.0e6, slope_loss=0.0, case_outputs=case_dict, error=e)
            end
        end
        mean_signal_loss = total_loss / length(keys)
        slope_loss = (!is_cooling && length(keys) == length(keys_heating_v50)) ?
                     loss_flow_slopes_v50(case_dict) : 0.0
        return (signal_loss=mean_signal_loss, slope_loss=slope_loss,
                case_outputs=case_dict, error=nothing)
    end

    function loss_cases_v50(p, keys; is_cooling=false, nodes=15,
                            include_flow_slope=true)
        components = dataset_loss_components_v50(
            p, keys; is_cooling=is_cooling, nodes=nodes
        )
        components.error === nothing || return 1.0e6
        return components.signal_loss +
               (include_flow_slope ? components.slope_loss : 0.0)
    end

    function calibration_loss_components_v50(p; stage=:full, nodes=15)
        for i in eachindex(p)
            if p[i] < lb_full_v50[i] || p[i] > ub_full_v50[i]
                bound_penalty = 1.0e6 + 1.0e3 *
                    abs(p[i] - clamp(p[i], lb_full_v50[i], ub_full_v50[i]))
                return (signal_loss=bound_penalty, slope_loss=0.0,
                        capacitance_regularization=0.0, power_regularization=0.0,
                        total=bound_penalty, error=:bounds)
            end
        end

        dataset = dataset_loss_components_v50(
            p, keys_heating_v50; is_cooling=false, nodes=nodes
        )
        if dataset.error !== nothing
            return (signal_loss=1.0e6, slope_loss=0.0,
                    capacitance_regularization=0.0, power_regularization=0.0,
                    total=1.0e6, error=dataset.error)
        end
        cap_reg = capacitance_regularization_v50(p; nodes=nodes)
        power_reg = power_scale_regularization_v50(p)
        total = dataset.signal_loss + dataset.slope_loss + cap_reg + power_reg
        return (signal_loss=dataset.signal_loss, slope_loss=dataset.slope_loss,
                capacitance_regularization=cap_reg, power_regularization=power_reg,
                total=total, error=nothing)
    end

    function calibration_loss_v50(p; stage=:full, nodes=15)
        return calibration_loss_components_v50(p; stage=stage, nodes=nodes).total
    end
end

begin # exact energy conservation audit & macroscopic LMTD HTC
    function macroscopic_htc_v50(result, idx)
        Q_core = sum(result.receiver_gas_heat[:, idx])
        T_solid_in = result.core_temperature[1, idx]
        T_solid_out = result.core_temperature[end, idx]
        T_gas_in = result.gas_temperature[1, idx]
        T_gas_out = result.gas_temperature[length(result.z_solid) + 1, idx]
        
        dT1 = max(0.1, T_solid_in - T_gas_in)
        dT2 = max(0.1, T_solid_out - T_gas_out)
        
        LMTD = (abs(dT1 - dT2) < 1e-4) ? dT1 : (dT1 - dT2) / log(dT1 / dT2)
        A_total = P_exchange_v50 * L_v50
        
        return Q_core / (A_total * max(LMTD, 0.1))
    end

    function compute_energy_balance_v50(result, operating, time_val, p; nodes=default_nodes_v50, rear_nodes=REAR_TUBE_DEFAULT_NODES_v50)
        dx = L_v50 / nodes
        dx_rear = REAR_TUBE_LENGTH_v50 / rear_nodes
        irradiance = operating.irradiance(time_val)
        ambient = operating.ambient_temperature(time_val)
        flow = max(0.0, operating.flow_lpm(time_val))
        Tin = operating.inlet_temperature(time_val)

        Q_aperture = max(0.0, irradiance) * A_frt_v50
        M = absorbed_power_scale_v50(irradiance, p)
        Q_delivered = M * Q_aperture
        chi = core_source_fraction_v50(p)
        Qcore_solar = chi * Q_delivered
        Qperim_solar = (1.0 - chi) * Q_delivered

        # Temperatures at this time
        idx = argmin(abs.(result.time .- time_val))
        Tcore = result.core_temperature[:, idx]
        Tperim = result.perim_temperature[:, idx]
        Ttube = result.rear_temperature[:, idx]
        Tcavity = result.cavity_temperature[idx]
        Tg = result.gas_temperature[:, idx]

        mdot = m_dot_standard_v50(flow)
        
        # Enthalpy transferred into gas
        h_suction = suction_heat_transfer_coefficient_v50(p)
        Q_suction = flow > 1e-12 ? h_suction * A_frt_v50 * (Tperim[1] - Tin) : 0.0
        Q_gas_core = sum(result.receiver_gas_heat[:, idx])
        Q_gas_rear = sum(result.rear_tube_gas_heat[:, idx])
        Q_gas_total = Q_suction + Q_gas_core + Q_gas_rear
        
        Q_front_rad = EPS_FRONT_FIXED_v50 * sigma_sb_v50 * A_frt_v50 * (Tperim[1]^4 - ambient^4)
        Q_cavity_amb = cavity_loss_to_ambient_v50(Tcavity, ambient)
        
        Q_flange = 0.0
        for j in 1:rear_nodes
            G_fl = rear_tube_flange_conductance_per_length_v50(result.z_rear[j], Ttube[j], p)
            Q_flange += G_fl * dx_rear * (Ttube[j] - WATER_FLANGE_TEMPERATURE_v50)
        end

        # Calculate instantaneous sensible storage rate dE_stored/dt for the solid network
        u_k = result.ode_solution.u[idx]
        du_k = similar(u_k)
        core_weights = solar_weights_v50(core_optical_extinction_v50(p), front_deposition_fraction_v50(p), nodes)
        rear_weights = rear_contact_weights_v50(result.z_solid)
        context = (
            p=p, operating=operating, z=result.z_solid, dx=dx,
            core_weights=core_weights, rear_weights=rear_weights,
            z_rear=result.z_rear, dx_rear=dx_rear,
            Tg=zeros(nodes + rear_nodes + 1), Qgas=zeros(nodes), hcell=zeros(nodes),
            Qgas_rear=zeros(rear_nodes), hrear=zeros(rear_nodes), zg=result.z_gas,
        )
        receiver_rhs_v50!(
            du_k, u_k, context, time_val
        )
        
        dE_stored = 0.0
        cell_volume = A_solid_v50 * dx
        perim_cap_cell = perimeter_heat_capacity_total_v50(p) / nodes
        for i in 1:nodes
            core_cap = rho_s_v50 * Cps_v50(Tcore[i]) * cell_volume
            dE_stored += core_cap * du_k[i]
            dE_stored += perim_cap_cell * du_k[nodes + i]
        end
        for j in 1:rear_nodes
            tube_cap = rear_tube_capacity_v50(Ttube[j], result.z_rear[j], dx_rear)
            dE_stored += tube_cap * du_k[2 * nodes + j]
        end
        cavity_index = 2 * nodes + rear_nodes + 1
        dE_stored += CAVITY_HEAT_CAPACITY_v50 * du_k[cavity_index]
        C_rear_tot = rear_reservoir_heat_capacity_v50(p)
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
            macro_htc=macroscopic_htc_v50(result, idx),
        )
    end

    function calibrate_v50(; nodes=15, maximum_iterations=500,
                           maximum_time_seconds=600.0, stage=:full,
                           dataset=:heating, initial_parameters=pnew_v50)
        lowercase(string(dataset)) == "heating" ||
            throw(ArgumentError("v50 calibration currently supports dataset=:heating only"))
        p0 = copy(initial_parameters)
        idx_fit = fit_stage_indices_v50(stage)
        function obj_free(theta)
            p_full = copy(p0)
            p_full[idx_fit] .= theta
            return calibration_loss_v50(p_full; stage=stage, nodes=nodes)
        end

        res = optimize_with_nlopt_v3(
            obj_free,
            p0[idx_fit],
            lb_full_v50[idx_fit],
            ub_full_v50[idx_fit];
            maximum_iterations=maximum_iterations,
            maximum_time_seconds=maximum_time_seconds,
            algorithm=OptimizationNLopt.NLopt.LN_BOBYQA(),
            label="fit-v50-$(lowercase(string(stage)))",
        )
        p_opt = copy(p0)
        p_opt[idx_fit] .= res.minimizer
        return (parameters=p_opt, objective=res.minimum, retcode=res.retcode,
                evaluations=res.evaluations, stage=stage, fitted_indices=idx_fit)
    end
end
