# ============================================================================
# 2D_v10.jl - Front-Exchange Conservative 2D Axisymmetric Macro-ECM
# ============================================================================
# Correction relative to v9:
#   1. Flow-sensitive heat transfer from the exposed SiC front web to the
#      entering gas is represented with a bounded effectiveness closure.
#      The same heat is removed from the front solid and added to the gas.
#   2. Natural convection remains on the front felt/casing but is replaced by
#      this internal exchange on the receiver web. Front radiation remains an
#      external loss, and the centered irradiance field is unchanged.
#
# Corrections inherited from v9:
#   1. The Aalborg air-equivalent flow is treated as standard volumetric flow
#      at 101.4 kPa and 294.25 K, then converted once to conserved mass flow.
#   2. Local gas density and channel velocity vary with the computed gas
#      temperature. The mass-flux Reynolds expression remains unchanged.
#   3. A square-channel laminar pressure model predicts the atmospheric-
#      referenced flush-tap DP1 signal. Cold t0 data determine one sensor
#      offset and one shared hydraulic-resistance scale; no bypass is fitted.
#   4. The verified T8 depth is 5 mm rather than 11 mm.
#
# Corrections inherited from v8:
#   1. The 19 x 19 mm receiver area is 3.61e-4 m2 and its area-equivalent
#      radius is 10.7196 mm. This restores porosity near 0.623 and mass near 40 g.
#   2. The full measured gas mass flow is conserved through the channel rings.
#   3. Optical transmission is retained; Beer-Lambert weights are not normalized
#      to force complete absorption. Rim spillage is deposited only at the front.
#   4. The unused cavity state/heat-capacity parameter is removed.
#   5. Rear-tube geometry, terminal flange boundary, T3 interpolation, cooling
#      t90, SiC conductivity, and physical receiver/felt contact resistance are corrected.
#   6. The fitted heat-transfer closure is explicitly an apparent Macro-ECM
#      Nusselt law; the non-identifiable single-channel Nu=3.61 floor is removed.
# ============================================================================

module Receiver2D_v10

using DifferentialEquations
using OrdinaryDiffEq
using SciMLBase
using LinearAlgebra
using Statistics

export Geometry2D, SolidProperties2D, HeatTransferParameters2D, LossParameters2D
export OpticalParameters2D, HydraulicParameters2D, ModelParameters2D
export OperatingCondition2D, SimulationResult2D
export default_parameters2D, simulate2D, sensor_predictions2D, get_t90_2D, normalized_slope_mse_2D
export linear_history, sic_conductivity, sic_heat_capacity
export air_conductivity, air_heat_capacity, air_viscosity, air_density
export standard_air_density2D, standard_mass_flow2D, ideal_square_channel_dp_slope2D
export alumina_conductivity, alumina_heat_capacity
export felt_conductivity_temp, front_nusselt_correlation, build_initial_state_2D
export pack_parameters2D, unpack_parameters2D, calibrate2D, LB_2D, UB_2D, FIT_INDICES_2D
export geometry_invariants2D, solar_power_budget2D, energy_rate_ledger2D
export hydraulic_summary2D

const SIGMA = 5.670374419e-8
const WATER_FLANGE_TEMP = 298.15
const G_ACCEL = 9.81
const AIR_SPECIFIC_GAS_CONSTANT = 287.05
const SQUARE_DUCT_DARCY_POISEUILLE = 56.91

# ============================================================================
# Struct Definitions
# ============================================================================

Base.@kwdef struct Geometry2D
    length::Float64 = 137e-3
    receiver_width::Float64 = 19.0e-3           # Physical square monolith width
    receiver_radius::Float64 = sqrt((19.0e-3)^2 / pi) # Area-equivalent disc radius
    insulation_outer_radius::Float64 = 57.0e-3  # Alumina felt outer radius
    casing_outer_radius::Float64 = 75.0e-3      # Aluminum casing outer radius
    channel_width::Float64 = 1.5e-3             # Hydraulic diameter Dh
    wall_thickness::Float64 = 0.4e-3            # Monolith web wall thickness
    channel_count::Int = 100                    # Total channel count
    nodes_r_rec::Int = 14                       # Mesh-verified radial nodes in SiC disc
    nodes_r_felt::Int = 7                       # Radial nodes in alumina felt insulation
    nodes_r_case::Int = 3                       # Radial nodes in aluminum casing
    nodes_z::Int = 45                           # Refined non-uniform axial mesh
    rear_tube_length::Float64 = 150.0e-3        # Rear alumina exit tube length
    rear_tube_nodes::Int = 30                   # Rear tube axial nodes
    rear_tube_inner_radius::Float64 = 6.5e-3
    rear_tube_outer_radius::Float64 = 8.0e-3
    t3_distance_from_receiver::Float64 = 3.0e-3
end

Base.@kwdef struct SolidProperties2D
    density::Float64 = 2150.0                   # Measured effective SiC solid density [kg/m3]
    radial_conductivity_scale::Float64 = 0.05   # Radial web/tortuosity scale
    axial_conductivity_scale::Float64 = 1.00    # Axial wall-connectivity scale
    rad_extinction_coeff::Float64 = 50.0        # Thermal infrared extinction beta_rad [1/m]
    felt_conductivity_ref::Float64 = 0.06       # Baseline alumina felt thermal conductivity [W/m/K]
    felt_density::Float64 = 140.0               # Alumina felt density [kg/m3]
    felt_heat_capacity::Float64 = 1360.0        # Alumina felt heat capacity [J/kg/K]
    casing_conductivity::Float64 = 205.0        # Aluminum casing thermal conductivity [W/m/K]
    casing_density::Float64 = 2700.0            # Aluminum casing density [kg/m3]
    casing_heat_capacity::Float64 = 900.0       # Aluminum casing heat capacity [J/kg/K]
    receiver_felt_contact_resistance::Float64 = 1.0e-3 # Areal resistance [m2 K/W]
end

Base.@kwdef struct HeatTransferParameters2D
    coefficient::Float64 = 1.0e-3               # Apparent Macro-ECM Nu prefactor A_Nu
    reynolds_exponent::Float64 = 1.440           # Fixed independent apparent exponent B_Re
    prandtl_exponent::Float64 = 0.333           # Prandtl exponent C_Pr
    minimum_nusselt::Float64 = 0.0              # No single-channel floor in apparent closure
    c_radial_flow::Float64 = 0.10               # Preferential central core channel flow fraction
    front_coefficient::Float64 = 0.0             # Front-web Nu prefactor C_f (sensitivity parameter)
    front_reynolds_exponent::Float64 = 0.5       # Fixed front-web Reynolds exponent m_f
end

Base.@kwdef struct LossParameters2D
    front_emissivity::Float64 = 0.85            # Front face emissivity
    casing_emissivity::Float64 = 0.20           # Casing outer emissivity
    front_loss_scale::Float64 = 1.0
    casing_loss_scale::Float64 = 1.0
    rear_adaptor_conductance::Float64 = 0.448   # 100 W/m2/K over channel-wall overlap
    flange_conductance_scale::Float64 = 1.00    # Water-cooled flange sink scale
end

Base.@kwdef struct OpticalParameters2D
    absorbed_fraction::Float64 = 1.00            # Power scales represent net aperture delivery
    extinction_coefficient::Float64 = 110.0      # Deeper optical extinction beta_opt [1/m]
    beam_radius_sigma::Float64 = 14.0e-3         # Gaussian beam radial width sigma [m]
    spillage_fraction::Float64 = 0.10            # Fraction of power spilling to perimeter rim
    front_deposition_fraction::Float64 = 0.20    # Front surface absorption fraction
    scale_456::Float64 = 1.35                    # Initial net delivered-power factor
    scale_304::Float64 = 1.35
    scale_256::Float64 = 0.85
end

Base.@kwdef struct HydraulicParameters2D
    standard_pressure::Float64 = 101.4e3       # Aalborg default reference [Pa absolute]
    standard_temperature::Float64 = 294.25     # Aalborg default reference [K, 70 F]
    atmospheric_pressure::Float64 = 101.325e3  # Receiver open-face pressure [Pa absolute]
    mass_flow_scale::Float64 = 1.0             # Global DP1-validated correction, not per-run
    dp1_zero_offset_mbar::Float64 = -0.614226030202630
    hydraulic_resistance_scale::Float64 = 1.95171196637134
    minor_loss_coefficient::Float64 = 0.0      # Shared hot-excess path-loss coefficient
end

Base.@kwdef struct ModelParameters2D
    geometry::Geometry2D = Geometry2D()
    solid::SolidProperties2D = SolidProperties2D()
    heat_transfer::HeatTransferParameters2D = HeatTransferParameters2D()
    losses::LossParameters2D = LossParameters2D()
    optics::OpticalParameters2D = OpticalParameters2D()
    hydraulics::HydraulicParameters2D = HydraulicParameters2D()
end

struct OperatingCondition2D
    irradiance::Function
    flow_lpm::Function
    inlet_temperature::Function
    ambient_temperature::Function
end

function OperatingCondition2D(;
    irradiance=t -> 0.0,
    flow_lpm=t -> 10.0,
    inlet_temperature=t -> 295.0,
    ambient_temperature=t -> 295.0
)
    irr_fn = irradiance isa Function ? (t -> Float64(irradiance(t))) : (t -> Float64(irradiance))
    flow_fn = flow_lpm isa Function ? (t -> Float64(flow_lpm(t))) : (t -> Float64(flow_lpm))
    tin_fn = inlet_temperature isa Function ? (t -> Float64(inlet_temperature(t))) : (t -> Float64(inlet_temperature))
    tamb_fn = ambient_temperature isa Function ? (t -> Float64(ambient_temperature(t))) : (t -> Float64(ambient_temperature))
    return OperatingCondition2D(irr_fn, flow_fn, tin_fn, tamb_fn)
end

struct SimulationResult2D
    time::Vector{Float64}
    nr_rec::Int
    receiver_radius::Float64
    r_solid::Vector{Float64}
    z_solid::Vector{Float64}
    z_rear::Vector{Float64}
    z_rear_gas::Vector{Float64}
    r_gas::Vector{Float64}
    z_gas::Vector{Float64}
    solid_temperature::Array{Float64, 3}       # (Nr_total, Nz, Nt)
    gas_temperature::Array{Float64, 3}         # (Nr_rec, Nz+1, Nt)
    rear_tube_temperature::Array{Float64, 2}   # (Nz_rear, Nt)
    rear_gas_temperature::Array{Float64, 2}    # (Nz_rear+1, Nt)
    heat_transfer_coefficient::Array{Float64, 3} # (Nr_rec, Nz, Nt)
    front_heat_transfer_coefficient::Array{Float64, 2} # (Nr_rec, Nt)
    front_gas_heat_transfer_W::Array{Float64, 2} # (Nr_rec, Nt), solid -> gas
    gas_density::Array{Float64, 3}             # (Nr_rec, Nz, Nt)
    gas_velocity::Array{Float64, 3}            # (Nr_rec, Nz, Nt)
    gas_reynolds::Array{Float64, 3}            # (Nr_rec, Nz, Nt)
    cell_pressure_drop::Array{Float64, 3}      # (Nr_rec, Nz, Nt) [Pa]
    receiver_pressure_drop_mbar::Vector{Float64}
    dp1_prediction_mbar::Vector{Float64}
    mass_flow_kg_s::Vector{Float64}
    ode_solution::Any
end

# ============================================================================
# Material & Transport Property Functions
# ============================================================================

function linear_history(times::AbstractVector, values::AbstractVector)
    length(times) == length(values) || error("times and values must match in length.")
    t_min, t_max = times[1], times[end]
    return function (t)
        t_clamped = clamp(t, t_min, t_max)
        idx = searchsortedlast(times, t_clamped)
        if idx == 0
            return Float64(values[1])
        elseif idx >= length(times)
            return Float64(values[end])
        else
            dt = times[idx + 1] - times[idx]
            frac = dt <= 1e-12 ? 0.0 : (t_clamped - times[idx]) / dt
            return Float64((1.0 - frac) * values[idx] + frac * values[idx + 1])
        end
    end
end

clamp_temperature(T) = max(200.0, Float64(T))

function sic_conductivity(T)
    Tk = clamp_temperature(T)
    # COMSOL alpha-polycrystalline SiC correlation used by the established
    # 0D/1D receiver lineage. T is in kelvin and k is in W/m/K.
    val = 191.9216 - 0.3261784 * Tk + 2.739462e-4 * Tk^2 - 7.70926e-8 * Tk^3
    return max(2.0, val)
end

function sic_heat_capacity(T)
    Tk = clamp_temperature(T)
    val = 1110.0 + 0.15 * Tk - 4.2e5 / (Tk^2)
    return max(400.0, val)
end

function air_conductivity(T)
    Tk = clamp_temperature(T)
    return 2.414e-3 * sqrt(Tk) / (1.0 + 245.4 * (10.0^(-12.0 / Tk)) / Tk)
end

function air_heat_capacity(T)
    Tk = clamp_temperature(T)
    return (1.0 + 1.983e-4 * Tk - 4.14e-8 * Tk^2) * 1004.0
end

function air_viscosity(T)
    Tk = clamp_temperature(T)
    return 1.458e-6 * (Tk^1.5) / (Tk + 110.4)
end

function air_density(T, pressure=101.325e3)
    Tk = clamp_temperature(T)
    return max(1.0, Float64(pressure)) / (AIR_SPECIFIC_GAS_CONSTANT * Tk)
end

standard_air_density2D(h::HydraulicParameters2D=HydraulicParameters2D()) =
    air_density(h.standard_temperature, h.standard_pressure)

function standard_mass_flow2D(flow_lpm::Real,
                              h::HydraulicParameters2D=HydraulicParameters2D())
    return h.mass_flow_scale * standard_air_density2D(h) *
           max(0.0, Float64(flow_lpm)) / 60000.0
end

function ideal_square_channel_dp_slope2D(
    p::ModelParameters2D=default_parameters2D(),
)
    g = p.geometry
    flow_area = g.channel_count * g.channel_width^2
    velocity_per_lpm = 1.0e-3 / 60.0 / flow_area
    mu_ref = air_viscosity(p.hydraulics.standard_temperature)
    dp_pa_per_lpm = 0.5 * SQUARE_DUCT_DARCY_POISEUILLE * mu_ref *
                    g.length * velocity_per_lpm / g.channel_width^2
    return dp_pa_per_lpm / 100.0
end

function alumina_conductivity(T)
    Tc = clamp_temperature(T) - 273.15
    return max(1.0, 5.5 + 34.5 * exp(-0.0033 * Tc))
end

function alumina_heat_capacity(T)
    Tk = clamp_temperature(T)
    return max(500.0, (1.00446 + 1.742e-4 * Tk - 2.796e4 * (Tk^-2)) * 1000.0)
end

function felt_conductivity_temp(T, ref_k=0.06)
    Tk = clamp_temperature(T)
    return max(0.05, ref_k + 1.2e-10 * Tk^3)
end

function front_nusselt_correlation(T_surface, T_ambient, L_char=0.0678)
    Ts = clamp_temperature(T_surface)
    Tamb = clamp_temperature(T_ambient)
    dT = max(0.1, Ts - Tamb)
    Tfilm = 0.5 * (Ts + Tamb)

    kf = air_conductivity(Tfilm)
    cp = air_heat_capacity(Tfilm)
    mu = air_viscosity(Tfilm)
    rho = 101325.0 / (287.05 * Tfilm)
    nu = mu / rho
    alpha = kf / (rho * cp)
    beta_exp = 1.0 / Tfilm

    Ra = G_ACCEL * beta_exp * dT * (L_char^3) / (nu * alpha)
    Pr = nu / alpha

    Nu_front = (0.825 + (0.387 * Ra^(1/6)) / ((1.0 + (0.492 / Pr)^(9/16))^(8/27)))^2
    h_conv = max(10.0, Nu_front * kf / L_char)
    return h_conv
end

default_parameters2D() = ModelParameters2D()

# ============================================================================
# Multi-Domain 2D Grid Discretization Struct
# ============================================================================

struct DiscretizationGrid2D
    nr_rec::Int
    nr_felt::Int
    nr_case::Int
    nr_total::Int
    nz::Int
    nz_rear::Int
    dz_cells::Vector{Float64}
    dz_rear::Float64
    r_faces::Vector{Float64}
    r_centers::Vector{Float64}
    z_faces::Vector{Float64}
    z_centers::Vector{Float64}
    z_rear_faces::Vector{Float64}
    z_rear_centers::Vector{Float64}
    area_frt::Vector{Float64}
    area_solid::Vector{Float64}
    area_flow::Vector{Float64}
    volume_cell::Matrix{Float64}
    perimeter_ex::Vector{Float64}
    porosity::Vector{Float64}
    dh::Float64
    total_frt::Float64
end

function build_grid2D(p::ModelParameters2D)
    g = p.geometry
    nr_rec, nr_felt, nr_case = g.nodes_r_rec, g.nodes_r_felt, g.nodes_r_case
    nr_total = nr_rec + nr_felt + nr_case
    nz, nz_rear = g.nodes_z, g.rear_tube_nodes
    
    R_rec = g.receiver_radius
    R_felt = g.insulation_outer_radius
    R_case = g.casing_outer_radius
    L_rec = g.length
    L_rear = g.rear_tube_length
    isapprox(pi * R_rec^2, g.receiver_width^2; rtol=1e-10) ||
        error("receiver_radius must be area-equivalent to receiver_width")

    # Refined Non-Uniform Axial Mesh
    z_faces = zeros(nz + 1)
    z_faces[1] = 0.0
    
    n_fine = 8
    z_fine_end = 15.0e-3
    for j in 1:n_fine
        z_faces[j + 1] = z_fine_end * (j / n_fine)
    end
    
    n_coarse = nz - n_fine
    for j in 1:n_coarse
        frac = j / n_coarse
        z_faces[n_fine + j + 1] = z_fine_end + (L_rec - z_fine_end) * frac
    end

    dz_cells = [z_faces[j + 1] - z_faces[j] for j in 1:nz]
    z_centers = [0.5 * (z_faces[j] + z_faces[j + 1]) for j in 1:nz]

    dz_rear = L_rear / nz_rear
    z_rear_faces = collect(range(0.0, L_rear, length=nz_rear + 1))
    z_rear_centers = [0.5 * (z_rear_faces[j] + z_rear_faces[j + 1]) for j in 1:nz_rear]

    r_faces_rec = collect(range(0.0, R_rec, length=nr_rec + 1))
    r_faces_felt = collect(range(R_rec, R_felt, length=nr_felt + 1))[2:end]
    r_faces_case = collect(range(R_felt, R_case, length=nr_case + 1))[2:end]
    r_faces = vcat(r_faces_rec, r_faces_felt, r_faces_case)
    r_centers = [0.5 * (r_faces[i] + r_faces[i + 1]) for i in 1:nr_total]

    total_frt = pi * R_rec^2
    area_frt = [pi * (r_faces[i + 1]^2 - r_faces[i]^2) for i in 1:nr_total]

    ch_area = g.channel_width^2
    total_ch_area = g.channel_count * ch_area
    total_ch_area < total_frt || error("channel flow area must be smaller than receiver frontal area")
    porosity_rec = total_ch_area / total_frt

    porosity = zeros(nr_total)
    area_flow = zeros(nr_total)
    area_solid = zeros(nr_total)
    perimeter_ex = zeros(nr_total)

    for i in 1:nr_total
        if i <= nr_rec
            porosity[i] = porosity_rec
            area_flow[i] = area_frt[i] * porosity_rec
            area_solid[i] = area_frt[i] * (1.0 - porosity_rec)
            perimeter_ex[i] = area_frt[i] * (4.0 * g.channel_count * g.channel_width / total_frt)
        else
            porosity[i] = 0.0
            area_flow[i] = 0.0
            area_solid[i] = area_frt[i]
            perimeter_ex[i] = 0.0
        end
    end

    volume_cell = zeros(nr_total, nz)
    for i in 1:nr_total, j in 1:nz
        volume_cell[i, j] = area_frt[i] * dz_cells[j]
    end

    return DiscretizationGrid2D(
        nr_rec, nr_felt, nr_case, nr_total, nz, nz_rear, dz_cells, dz_rear,
        r_faces, r_centers, z_faces, z_centers,
        z_rear_faces, z_rear_centers,
        area_frt, area_solid, area_flow, volume_cell, perimeter_ex,
        porosity, g.channel_width, total_frt
    )
end

function geometry_invariants2D(p::ModelParameters2D=default_parameters2D())
    grid = build_grid2D(p)
    receiver_area = grid.total_frt
    flow_area = sum(grid.area_flow[1:grid.nr_rec])
    solid_area = sum(grid.area_solid[1:grid.nr_rec])
    receiver_mass = p.solid.density * solid_area * p.geometry.length
    return (
        receiver_area=receiver_area,
        physical_square_area=p.geometry.receiver_width^2,
        equivalent_radius=p.geometry.receiver_radius,
        porosity=flow_area / receiver_area,
        solid_fraction=solid_area / receiver_area,
        receiver_mass=receiver_mass,
    )
end

function solar_weights2D(grid::DiscretizationGrid2D, p::ModelParameters2D)
    nr_rec, nz = grid.nr_rec, grid.nz
    sigma = p.optics.beam_radius_sigma
    beta_opt = p.optics.extinction_coefficient
    f_front = p.optics.front_deposition_fraction

    I_rad = [exp(-0.5 * (grid.r_centers[i] / sigma)^2) for i in 1:nr_rec]
    I_rad ./= sum(I_rad[i] * grid.area_frt[i] for i in 1:nr_rec) / grid.total_frt

    # Do not renormalize Beer-Lambert absorption over the finite receiver.
    # Its sum is 1-exp(-beta*L), and the remainder is transmitted.
    w_z = [exp(-beta_opt * grid.z_faces[j]) - exp(-beta_opt * grid.z_faces[j + 1]) for j in 1:nz]
    w_z_depth = (1.0 - f_front) .* w_z
    w_z_depth[1] += f_front

    weights = zeros(grid.nr_total, nz)
    for i in 1:nr_rec, j in 1:nz
        weights[i, j] = I_rad[i] * (grid.area_frt[i] / grid.total_frt) * w_z_depth[j]
    end
    return weights
end

function build_initial_state_2D(grid::DiscretizationGrid2D, p::ModelParameters2D, t0_sensors::Dict{Symbol, Float64}, T_ambient::Float64)
    nr_total = grid.nr_total
    nz, nz_rear = grid.nz, grid.nz_rear
    u_length = nr_total * nz + nz_rear
    u0 = fill(T_ambient, u_length)

    T8 = get(t0_sensors, :T8, T_ambient)
    T12 = get(t0_sensors, :T12, T_ambient)
    T11 = get(t0_sensors, :T11, T_ambient)
    T9 = get(t0_sensors, :T9, T_ambient)
    T10 = get(t0_sensors, :T10, T_ambient)
    T3 = get(t0_sensors, :T3, T_ambient)
    T2 = get(t0_sensors, :T2, T_ambient)

    T_s0 = reshape(view(u0, 1:(nr_total * nz)), nr_total, nz)
    T_rear0 = view(u0, (nr_total * nz + 1):(nr_total * nz + nz_rear))
    
    R_rec = grid.r_faces[grid.nr_rec + 1]
    R_case = grid.r_faces[end]
    r_t2 = min(R_case, R_rec + 40.0e-3)
    r_core_sample = grid.r_centers[1]
    r_perim_sample = grid.r_centers[grid.nr_rec]

    function axial_sensor_profile(z_val, T_mid, T_deep)
        if z_val <= 5.0e-3
            return T8
        elseif z_val <= 58.0e-3
            frac = (z_val - 5.0e-3) / (58.0e-3 - 5.0e-3)
            return (1.0 - frac) * T8 + frac * T_mid
        elseif z_val <= 107.0e-3
            frac = (z_val - 58.0e-3) / (107.0e-3 - 58.0e-3)
            return (1.0 - frac) * T_mid + frac * T_deep
        else
            return T_deep
        end
    end

    # The sensor depths generally fall between cell centres.  Sampling a
    # piecewise profile with a corner at a sensor and then interpolating the
    # cell-centre values does not reproduce that sensor exactly.  Enforce each
    # measured value on the discrete profile used to initialize the state.
    function enforce_sample!(profile, target_z, target_value)
        j0, j1, w = _axis_interp_indices(grid.z_centers, target_z)
        if j0 == j1
            profile[j0] = target_value
        elseif w <= 0.5
            profile[j0] = (target_value - w * profile[j1]) / (1.0 - w)
        else
            profile[j1] = (target_value - (1.0 - w) * profile[j0]) / w
        end
        return profile
    end

    core_profile = [axial_sensor_profile(z, T9, T10) for z in grid.z_centers]
    perim_profile = [axial_sensor_profile(z, T12, T11) for z in grid.z_centers]
    enforce_sample!(core_profile, 58.0e-3, T9)
    enforce_sample!(core_profile, 107.0e-3, T10)
    enforce_sample!(perim_profile, 5.0e-3, T8)
    enforce_sample!(perim_profile, 58.0e-3, T12)
    enforce_sample!(perim_profile, 107.0e-3, T11)

    for j in 1:nz
        T_core_z = core_profile[j]
        T_perim_z = perim_profile[j]

        for i in 1:nr_total
            r_val = grid.r_centers[i]
            if r_val <= R_rec
                frac_r = clamp((r_val - r_core_sample) /
                               max(r_perim_sample - r_core_sample, eps(Float64)), 0.0, 1.0)^2
                T_s0[i, j] = (1.0 - frac_r) * T_core_z + frac_r * T_perim_z
            elseif r_val <= r_t2
                frac_out = clamp((r_val - R_rec) / max(r_t2 - R_rec, eps(Float64)), 0.0, 1.0)
                T_s0[i, j] = (1.0 - frac_out) * T_perim_z + frac_out * T2
            else
                frac_out = clamp((r_val - r_t2) / max(R_case - r_t2, eps(Float64)), 0.0, 1.0)
                T_s0[i, j] = (1.0 - frac_out) * T2 + frac_out * T_ambient
            end
        end
    end

    for j in 1:nz_rear
        z_val = grid.z_rear_centers[j]
        frac_z = clamp(z_val / max(p.geometry.t3_distance_from_receiver, eps(Float64)), 0.0, 1.0)
        T_rear0[j] = (1.0 - frac_z) * (0.5 * (T10 + T11)) + frac_z * T3
    end

    return u0
end

# ============================================================================
# Multi-Domain Thermal Properties Helper
# ============================================================================

function cell_thermal_properties2D(i, T, p::ModelParameters2D, grid::DiscretizationGrid2D)
    if i <= grid.nr_rec
        k_val = p.solid.radial_conductivity_scale * sic_conductivity(T)
        rho_val = p.solid.density
        cp_val = sic_heat_capacity(T)
    elseif i <= (grid.nr_rec + grid.nr_felt)
        k_val = felt_conductivity_temp(T, p.solid.felt_conductivity_ref)
        rho_val = p.solid.felt_density
        cp_val = p.solid.felt_heat_capacity
    else
        k_val = p.solid.casing_conductivity
        rho_val = p.solid.casing_density
        cp_val = p.solid.casing_heat_capacity
    end
    return (k=k_val, rho=rho_val, cp=cp_val)
end

# ============================================================================
# RHS & Gas Profile Marching Functions
# ============================================================================

function _gas_profile2D!(gas, qgas, hcell, front_qgas, front_h, density, velocity, reynolds_flow,
                         dp_cell, gas_rear, qgas_rear, hrear,
                         solid_temp, T_rear, time, p, op, grid)
    nr_rec, nz = grid.nr_rec, grid.nz
    nz_rear = grid.nz_rear
    flow_lpm = max(0.0, Float64(op.flow_lpm(time)))
    inlet = Float64(op.inlet_temperature(time))
    hp = p.heat_transfer

    mdot_total = standard_mass_flow2D(flow_lpm, p.hydraulics)
    
    dh = grid.dh

    if flow_lpm <= 1e-12
        fill!(qgas, 0.0)
        fill!(hcell, 0.0)
        fill!(front_qgas, 0.0)
        fill!(front_h, 0.0)
        fill!(velocity, 0.0)
        fill!(reynolds_flow, 0.0)
        fill!(dp_cell, 0.0)
        fill!(qgas_rear, 0.0)
        fill!(hrear, 0.0)
        gas .= inlet
        density .= air_density(inlet, p.hydraulics.atmospheric_pressure)
        gas_rear .= inlet
        return nothing
    end

    R_rec = grid.r_faces[nr_rec + 1]
    psi_r = [1.0 - hp.c_radial_flow * (grid.r_centers[i] / R_rec)^2 for i in 1:nr_rec]
    psi_sum = sum(psi_r[i] * grid.area_flow[i] for i in 1:nr_rec)

    for i in 1:nr_rec
        mdot_ring = mdot_total * (psi_r[i] * grid.area_flow[i] / psi_sum)
        area_ring = max(grid.area_flow[i], eps(Float64))

        # Flow-sensitive exchange at the exposed SiC web occurs before the
        # internal channel-wall march. The effectiveness form bounds the gas
        # inlet temperature between the supplied inlet and front-solid values.
        cp_in = air_heat_capacity(inlet)
        k_in = air_conductivity(inlet)
        mu_in = air_viscosity(inlet)
        pr_in = cp_in * mu_in / k_in
        re_front = mdot_ring * dh / (area_ring * mu_in)
        nu_front = max(0.0, hp.front_coefficient) *
                   max(re_front, eps(Float64))^hp.front_reynolds_exponent *
                   max(pr_in, eps(Float64))^(1.0 / 3.0)
        front_h[i] = nu_front * k_in / dh
        front_ntu = front_h[i] * grid.area_solid[i] /
                    max(mdot_ring * cp_in, eps(Float64))
        front_effectiveness = -expm1(-front_ntu)
        gas[i, 1] = inlet + front_effectiveness * (solid_temp[i, 1] - inlet)
        front_qgas[i] = mdot_ring * cp_in * (gas[i, 1] - inlet)
        
        for j in 1:nz
            dz_j = grid.dz_cells[j]
            Tfilm = clamp_temperature(0.5 * (solid_temp[i, j] + gas[i, j]))
            cp = air_heat_capacity(Tfilm)
            kf = air_conductivity(Tfilm)
            mu = air_viscosity(Tfilm)
            
            reynolds = mdot_ring * dh / (max(grid.area_flow[i], eps(Float64)) * mu)
            prandtl = cp * mu / kf
            
            entry_factor = (dh / max(grid.z_centers[j], dh / 2.0))^(1.0 / 3.0)
            Nu_dev = hp.coefficient * max(reynolds, eps(Float64))^hp.reynolds_exponent *
                     max(prandtl, eps(Float64))^hp.prandtl_exponent * entry_factor
            nusselt = max(hp.minimum_nusselt, Nu_dev)
            
            hcell[i, j] = nusselt * kf / dh
            ua = hcell[i, j] * grid.perimeter_ex[i] * dz_j
            effectiveness = -expm1(-ua / max(mdot_ring * cp, eps(Float64)))
            
            gas[i, j + 1] = gas[i, j] + effectiveness * (solid_temp[i, j] - gas[i, j])
            qgas[i, j] = mdot_ring * cp * (gas[i, j + 1] - gas[i, j])

            Tbulk = clamp_temperature(0.5 * (gas[i, j] + gas[i, j + 1]))
            rho_bulk = air_density(Tbulk, p.hydraulics.atmospheric_pressure)
            mu_bulk = air_viscosity(Tbulk)
            u_bulk = mdot_ring / (rho_bulk * area_ring)
            re_bulk = rho_bulk * u_bulk * dh / mu_bulk
            density[i, j] = rho_bulk
            velocity[i, j] = u_bulk
            reynolds_flow[i, j] = re_bulk
            dp_cell[i, j] = p.hydraulics.hydraulic_resistance_scale *
                            0.5 * SQUARE_DUCT_DARCY_POISEUILLE *
                            mu_bulk * dz_j * u_bulk / dh^2
        end
        # The empirical cold multiplier already contains the reference-state
        # inlet/outlet/path loss. Add only the temperature-dependent excess
        # dynamic pressure so a nonzero K does not double count the cold fit.
        rho_ref_path = air_density(
            p.hydraulics.standard_temperature,
            p.hydraulics.atmospheric_pressure,
        )
        u_ref_path = mdot_ring / (rho_ref_path * area_ring)
        dynamic_excess = 0.5 * (
            density[i, end] * velocity[i, end]^2 -
            rho_ref_path * u_ref_path^2
        )
        dp_cell[i, end] += p.hydraulics.minor_loss_coefficient *
                           dynamic_excess
    end

    total_enthalpy = sum(mdot_total * (psi_r[i] * grid.area_flow[i] / psi_sum) * air_heat_capacity(gas[i, end]) * gas[i, end] for i in 1:nr_rec)
    total_mcp = sum(mdot_total * (psi_r[i] * grid.area_flow[i] / psi_sum) * air_heat_capacity(gas[i, end]) for i in 1:nr_rec)
    gas_rear[1] = total_mcp > 0 ? total_enthalpy / total_mcp : inlet

    r_tube_dh = 2.0 * p.geometry.rear_tube_inner_radius
    r_tube_area = pi * p.geometry.rear_tube_inner_radius^2
    r_tube_perim = 2.0 * pi * p.geometry.rear_tube_inner_radius

    for j in 1:nz_rear
        Tfilm = clamp_temperature(0.5 * (T_rear[j] + gas_rear[j]))
        cp = air_heat_capacity(Tfilm)
        kf = air_conductivity(Tfilm)
        mu = air_viscosity(Tfilm)
        
        reynolds = mdot_total * r_tube_dh / (r_tube_area * mu)
        prandtl = cp * mu / kf
        nusselt = reynolds < 2300.0 ? 3.66 : 0.023 * reynolds^0.8 * prandtl^0.4
        
        hrear[j] = nusselt * kf / r_tube_dh
        ua = hrear[j] * r_tube_perim * grid.dz_rear
        effectiveness = -expm1(-ua / max(mdot_total * cp, eps(Float64)))
        
        gas_rear[j + 1] = gas_rear[j] + effectiveness * (T_rear[j] - gas_rear[j])
        qgas_rear[j] = mdot_total * cp * (gas_rear[j + 1] - gas_rear[j])
    end

    return nothing
end

function receiver_ode2D!(du, u, context, time)
    grid = context.grid
    nr_rec, nr_total = grid.nr_rec, grid.nr_total
    nz, nz_rear = grid.nz, grid.nz_rear
    p = context.p
    op = context.op
    weights = context.weights
    work = context.work

    T_s = reshape(view(u, 1:(nr_total * nz)), nr_total, nz)
    T_rear = view(u, (nr_total * nz + 1):(nr_total * nz + nz_rear))

    du_s = reshape(view(du, 1:(nr_total * nz)), nr_total, nz)
    du_rear = view(du, (nr_total * nz + 1):(nr_total * nz + nz_rear))

    g = p.geometry
    solid = p.solid
    loss = p.losses
    ambient = Float64(op.ambient_temperature(time))

    _gas_profile2D!(
        work.gas, work.qgas, work.h, work.front_qgas, work.front_h,
        work.density, work.velocity,
        work.reynolds, work.dp_cell, work.gas_rear, work.qgas_rear,
        work.hrear, T_s[1:nr_rec, :], T_rear, time, p, op, grid,
    )
    fill!(du_s, 0.0)
    fill!(du_rear, 0.0)

    for i in 1:(nr_total - 1), j in 1:nz
        dz_j = grid.dz_cells[j]
        prop_i = cell_thermal_properties2D(i, T_s[i, j], p, grid)
        prop_ip1 = cell_thermal_properties2D(i + 1, T_s[i + 1, j], p, grid)
        dr_i = grid.r_faces[i + 1] - grid.r_faces[i]
        dr_ip1 = grid.r_faces[i + 2] - grid.r_faces[i + 1]
        areal_resistance = 0.5 * dr_i / prop_i.k + 0.5 * dr_ip1 / prop_ip1.k
        if i == nr_rec
            areal_resistance += max(0.0, solid.receiver_felt_contact_resistance)
        end
        r_face = grid.r_faces[i + 1]
        area_r_face = 2.0 * pi * r_face * dz_j
        Qcond_r = area_r_face * (T_s[i, j] - T_s[i + 1, j]) / max(areal_resistance, eps(Float64))
        du_s[i, j] -= Qcond_r
        du_s[i + 1, j] += Qcond_r
    end

    for i in 1:nr_total, j in 1:(nz - 1)
        dz_avg = 0.5 * (grid.dz_cells[j] + grid.dz_cells[j + 1])
        prop_j = cell_thermal_properties2D(i, T_s[i, j], p, grid)
        prop_jp1 = cell_thermal_properties2D(i, T_s[i, j + 1], p, grid)
        
        # Combined solid conduction + High-Temperature Radiative Diffusion in channels
        k_rad_j = i <= nr_rec ? 16.0 * SIGMA * clamp_temperature(T_s[i, j])^3 / (3.0 * max(1.0, p.solid.rad_extinction_coeff)) : 0.0
        k_rad_jp1 = i <= nr_rec ? 16.0 * SIGMA * clamp_temperature(T_s[i, j + 1])^3 / (3.0 * max(1.0, p.solid.rad_extinction_coeff)) : 0.0

        kz_j = i <= nr_rec ? p.solid.axial_conductivity_scale * sic_conductivity(T_s[i, j]) + k_rad_j : prop_j.k
        kz_jp1 = i <= nr_rec ? p.solid.axial_conductivity_scale * sic_conductivity(T_s[i, j + 1]) + k_rad_jp1 : prop_jp1.k
        kz_face = 2.0 * kz_j * kz_jp1 / (kz_j + kz_jp1)
        
        area_ax = i <= nr_rec ? grid.area_solid[i] : grid.area_frt[i]
        Qcond_z = kz_face * area_ax * (T_s[i, j] - T_s[i, j + 1]) / dz_avg
        du_s[i, j] -= Qcond_z
        du_s[i, j + 1] += Qcond_z
    end

    delivered_power = p.optics.absorbed_fraction * max(0.0, Float64(op.irradiance(time))) * grid.total_frt
    spillage_power = p.optics.spillage_fraction * delivered_power
    core_power = delivered_power - spillage_power

    for i in 1:nr_rec, j in 1:nz
        du_s[i, j] += core_power * weights[i, j] - work.qgas[i, j]
    end
    for i in 1:nr_rec
        du_s[i, 1] -= work.front_qgas[i]
    end
    felt_indices = (nr_rec + 1):(nr_rec + grid.nr_felt)
    felt_front_area = sum(grid.area_frt[i] for i in felt_indices)
    for i in felt_indices
        du_s[i, 1] += spillage_power * grid.area_frt[i] / felt_front_area
    end

    for i in 1:nr_total
        area_f = i <= nr_rec ? grid.area_solid[i] : grid.area_frt[i]
        T_surf_clamped = clamp_temperature(T_s[i, 1])
        Qfront_rad = loss.front_emissivity * SIGMA * area_f * (T_surf_clamped^4 - ambient^4)
        # Receiver-web convection is the conservative front-to-gas exchange
        # above. Only the non-receiver front retains natural ambient convection.
        Qfront_conv = i <= nr_rec ? 0.0 :
            front_nusselt_correlation(T_s[i, 1], ambient) *
            area_f * (T_s[i, 1] - ambient)
        du_s[i, 1] -= loss.front_loss_scale * (Qfront_rad + Qfront_conv)
    end

    Qrear_exit = 0.0
    for i in 1:nr_rec
        Qrear_i = loss.rear_adaptor_conductance * (grid.area_solid[i] / sum(grid.area_solid[1:nr_rec])) * (T_s[i, nz] - T_rear[1])
        du_s[i, nz] -= Qrear_i
        Qrear_exit += Qrear_i
    end
    du_rear[1] += Qrear_exit

    A_tube_wall = pi * (g.rear_tube_outer_radius^2 - g.rear_tube_inner_radius^2)
    for j in 1:(nz_rear - 1)
        k_j = alumina_conductivity(T_rear[j])
        k_jp1 = alumina_conductivity(T_rear[j + 1])
        k_face = 2.0 * k_j * k_jp1 / (k_j + k_jp1)
        Qcond_rear = k_face * A_tube_wall * (T_rear[j] - T_rear[j + 1]) / grid.dz_rear
        du_rear[j] -= Qcond_rear
        du_rear[j + 1] += Qcond_rear
    end

    for j in 1:nz_rear
        Q_gas_exchange = work.qgas_rear[j]
        du_rear[j] -= Q_gas_exchange
    end
    Q_flange = loss.flange_conductance_scale *
               (alumina_conductivity(T_rear[end]) * A_tube_wall / (0.5 * grid.dz_rear)) *
               (T_rear[end] - WATER_FLANGE_TEMP)
    du_rear[end] -= Q_flange

    for j in 1:nz
        dz_j = grid.dz_cells[j]
        T_case_clamped = clamp_temperature(T_s[nr_total, j])
        Qcasing_loss = 10.0 * (2.0 * pi * g.casing_outer_radius * dz_j) * (T_s[nr_total, j] - ambient) +
                       loss.casing_emissivity * SIGMA * (2.0 * pi * g.casing_outer_radius * dz_j) * (T_case_clamped^4 - ambient^4)
        du_s[nr_total, j] -= loss.casing_loss_scale * Qcasing_loss
    end

    for i in 1:nr_total, j in 1:nz
        dz_j = grid.dz_cells[j]
        prop = cell_thermal_properties2D(i, T_s[i, j], p, grid)
        cell_vol = i <= nr_rec ? (grid.area_solid[i] * dz_j) : (grid.area_frt[i] * dz_j)
        cap = prop.rho * prop.cp * cell_vol
        du_s[i, j] /= max(cap, eps(Float64))
    end

    for j in 1:nz_rear
        cap_rear = 3900.0 * alumina_heat_capacity(T_rear[j]) * A_tube_wall * grid.dz_rear
        du_rear[j] /= max(cap_rear, eps(Float64))
    end
    return nothing
end

mutable struct Workspace2D
    gas::Array{Float64, 2}
    qgas::Array{Float64, 2}
    h::Array{Float64, 2}
    front_qgas::Vector{Float64}
    front_h::Vector{Float64}
    density::Array{Float64, 2}
    velocity::Array{Float64, 2}
    reynolds::Array{Float64, 2}
    dp_cell::Array{Float64, 2}
    gas_rear::Vector{Float64}
    qgas_rear::Vector{Float64}
    hrear::Vector{Float64}
end

function Workspace2D(grid::DiscretizationGrid2D)
    return Workspace2D(
        zeros(grid.nr_rec, grid.nz + 1),
        zeros(grid.nr_rec, grid.nz),
        zeros(grid.nr_rec, grid.nz),
        zeros(grid.nr_rec),
        zeros(grid.nr_rec),
        zeros(grid.nr_rec, grid.nz),
        zeros(grid.nr_rec, grid.nz),
        zeros(grid.nr_rec, grid.nz),
        zeros(grid.nr_rec, grid.nz),
        zeros(grid.nz_rear + 1),
        zeros(grid.nz_rear),
        zeros(grid.nz_rear)
    )
end

function hydraulic_summary2D(work::Workspace2D, p::ModelParameters2D,
                             op::OperatingCondition2D, grid::DiscretizationGrid2D,
                             time::Real)
    flow_lpm = max(0.0, Float64(op.flow_lpm(time)))
    mdot_total = standard_mass_flow2D(flow_lpm, p.hydraulics)
    if mdot_total <= eps(Float64)
        receiver_dp_mbar = 0.0
    else
        R_rec = grid.r_faces[grid.nr_rec + 1]
        psi_r = [1.0 - p.heat_transfer.c_radial_flow *
                 (grid.r_centers[i] / R_rec)^2 for i in 1:grid.nr_rec]
        psi_sum = sum(psi_r[i] * grid.area_flow[i] for i in 1:grid.nr_rec)
        fractions = [psi_r[i] * grid.area_flow[i] / psi_sum for i in 1:grid.nr_rec]
        ring_dp = vec(sum(work.dp_cell, dims=2))
        receiver_dp_mbar = sum(fractions .* ring_dp) / 100.0
    end
    return (
        standard_flow_lpm=flow_lpm,
        mass_flow_kg_s=mdot_total,
        receiver_pressure_drop_mbar=receiver_dp_mbar,
        dp1_prediction_mbar=p.hydraulics.dp1_zero_offset_mbar + receiver_dp_mbar,
    )
end

function solar_power_budget2D(p::ModelParameters2D, irradiance::Real)
    grid = build_grid2D(p)
    weights = solar_weights2D(grid, p)
    delivered = p.optics.absorbed_fraction * max(0.0, Float64(irradiance)) * grid.total_frt
    spillage = p.optics.spillage_fraction * delivered
    core_incident = delivered - spillage
    core_deposited = core_incident * sum(weights)
    transmitted = max(0.0, core_incident - core_deposited)
    deposited = core_deposited + spillage
    return (
        delivered=delivered,
        deposited=deposited,
        core_deposited=core_deposited,
        spillage_deposited=spillage,
        transmitted=transmitted,
        closure=delivered - deposited - transmitted,
    )
end

function energy_rate_ledger2D(u::AbstractVector, p::ModelParameters2D,
                              op::OperatingCondition2D, time::Real=0.0)
    grid = build_grid2D(p)
    nr_rec, nr_total = grid.nr_rec, grid.nr_total
    nz, nz_rear = grid.nz, grid.nz_rear
    expected_length = nr_total * nz + nz_rear
    length(u) == expected_length || error("state vector length mismatch in energy ledger")

    weights = solar_weights2D(grid, p)
    work = Workspace2D(grid)
    context = (grid=grid, nz_rear=nz_rear, p=p, op=op, weights=weights, work=work)
    du = zeros(Float64, expected_length)
    receiver_ode2D!(du, u, context, Float64(time))

    T_s = reshape(view(u, 1:(nr_total * nz)), nr_total, nz)
    T_rear = view(u, (nr_total * nz + 1):expected_length)
    du_s = reshape(view(du, 1:(nr_total * nz)), nr_total, nz)
    du_rear = view(du, (nr_total * nz + 1):expected_length)

    denergy_dt = 0.0
    for i in 1:nr_total, j in 1:nz
        prop = cell_thermal_properties2D(i, T_s[i, j], p, grid)
        volume = (i <= nr_rec ? grid.area_solid[i] : grid.area_frt[i]) * grid.dz_cells[j]
        denergy_dt += prop.rho * prop.cp * volume * du_s[i, j]
    end
    tube_area = pi * (p.geometry.rear_tube_outer_radius^2 - p.geometry.rear_tube_inner_radius^2)
    for j in 1:nz_rear
        capacity = 3900.0 * alumina_heat_capacity(T_rear[j]) * tube_area * grid.dz_rear
        denergy_dt += capacity * du_rear[j]
    end

    ambient = Float64(op.ambient_temperature(time))
    front_loss = 0.0
    for i in 1:nr_total
        area_f = i <= nr_rec ? grid.area_solid[i] : grid.area_frt[i]
        Ts = T_s[i, 1]
        front_loss += p.losses.front_loss_scale * area_f * (
            p.losses.front_emissivity * SIGMA * (clamp_temperature(Ts)^4 - ambient^4) +
            (i <= nr_rec ? 0.0 :
             front_nusselt_correlation(Ts, ambient) * (Ts - ambient))
        )
    end

    casing_loss = 0.0
    for j in 1:nz
        area = 2.0 * pi * p.geometry.casing_outer_radius * grid.dz_cells[j]
        Ts = T_s[nr_total, j]
        casing_loss += p.losses.casing_loss_scale * area * (
            10.0 * (Ts - ambient) +
            p.losses.casing_emissivity * SIGMA * (clamp_temperature(Ts)^4 - ambient^4)
        )
    end

    flange_loss = p.losses.flange_conductance_scale *
                  (alumina_conductivity(T_rear[end]) * tube_area / (0.5 * grid.dz_rear)) *
                  (T_rear[end] - WATER_FLANGE_TEMP)
    gas_receiver = sum(work.qgas)
    gas_front = sum(work.front_qgas)
    gas_rear = sum(work.qgas_rear)
    solar = solar_power_budget2D(p, op.irradiance(time))
    expected_denergy_dt = solar.deposited - gas_front - gas_receiver - gas_rear -
                          front_loss - casing_loss - flange_loss
    return (
        denergy_dt=denergy_dt,
        expected_denergy_dt=expected_denergy_dt,
        residual=denergy_dt - expected_denergy_dt,
        solar_deposited=solar.deposited,
        solar_transmitted=solar.transmitted,
        gas_front=gas_front,
        gas_receiver=gas_receiver,
        gas_rear=gas_rear,
        front_loss=front_loss,
        casing_loss=casing_loss,
        flange_loss=flange_loss,
    )
end

# ============================================================================
# Main Simulation API
# ============================================================================

function simulate2D(
    p::ModelParameters2D,
    op::OperatingCondition2D,
    save_times::AbstractVector;
    initial_temperature=Float64(op.ambient_temperature(save_times[1])),
    solver=FBDF(autodiff=AutoFiniteDiff()),
    reltol=1e-6,
    abstol=1e-7,
    dtmax=30.0
)
    grid = build_grid2D(p)
    nr_rec, nr_total = grid.nr_rec, grid.nr_total
    nz, nz_rear = grid.nz, grid.nz_rear
    weights = solar_weights2D(grid, p)
    work = Workspace2D(grid)

    context = (grid=grid, nz_rear=nz_rear, p=p, op=op, weights=weights, work=work)
    u_length = nr_total * nz + nz_rear

    if initial_temperature isa AbstractVector
        length(initial_temperature) == u_length || error("initial_temperature vector length mismatch.")
        u_initial = copy(Float64.(initial_temperature))
    else
        u_initial = fill(Float64(initial_temperature), u_length)
    end

    t_span = (Float64(save_times[1]), Float64(save_times[end]))
    prob = ODEProblem(receiver_ode2D!, u_initial, t_span, context)

    sol = solve(prob, solver; saveat=save_times, reltol=reltol, abstol=abstol, dtmax=dtmax,
                isoutofdomain=(u, p, t) -> any(x -> isnan(x) || x < 100.0 || x > 3500.0, u))

    nt = length(sol.t)
    solid_T = zeros(nr_total, nz, nt)
    gas_T = zeros(nr_rec, nz + 1, nt)
    rear_tube_T = zeros(nz_rear, nt)
    rear_gas_T = zeros(nz_rear + 1, nt)
    htc = zeros(nr_rec, nz, nt)
    front_htc = zeros(nr_rec, nt)
    front_qgas = zeros(nr_rec, nt)
    gas_rho = zeros(nr_rec, nz, nt)
    gas_velocity = zeros(nr_rec, nz, nt)
    gas_reynolds = zeros(nr_rec, nz, nt)
    cell_dp = zeros(nr_rec, nz, nt)
    receiver_dp_mbar = zeros(nt)
    dp1_prediction_mbar = zeros(nt)
    mass_flow_kg_s = zeros(nt)

    for k in 1:nt
        u_k = sol.u[k]
        t_k = sol.t[k]
        
        solid_T[:, :, k] .= reshape(view(u_k, 1:(nr_total * nz)), nr_total, nz)
        rear_tube_T[:, k] .= view(u_k, (nr_total * nz + 1):(nr_total * nz + nz_rear))

        _gas_profile2D!(
            work.gas, work.qgas, work.h, work.front_qgas, work.front_h,
            work.density, work.velocity,
            work.reynolds, work.dp_cell, work.gas_rear, work.qgas_rear,
            work.hrear, solid_T[1:nr_rec, :, k], rear_tube_T[:, k],
            t_k, p, op, grid,
        )
        gas_T[:, :, k] .= work.gas
        rear_gas_T[:, k] .= work.gas_rear
        htc[:, :, k] .= work.h
        front_htc[:, k] .= work.front_h
        front_qgas[:, k] .= work.front_qgas
        gas_rho[:, :, k] .= work.density
        gas_velocity[:, :, k] .= work.velocity
        gas_reynolds[:, :, k] .= work.reynolds
        cell_dp[:, :, k] .= work.dp_cell
        hyd = hydraulic_summary2D(work, p, op, grid, t_k)
        receiver_dp_mbar[k] = hyd.receiver_pressure_drop_mbar
        dp1_prediction_mbar[k] = hyd.dp1_prediction_mbar
        mass_flow_kg_s[k] = hyd.mass_flow_kg_s
    end

    return SimulationResult2D(
        Vector{Float64}(sol.t),
        nr_rec,
        p.geometry.receiver_radius,
        grid.r_centers,
        grid.z_centers,
        grid.z_rear_centers,
        grid.z_rear_faces,
        grid.r_centers[1:nr_rec],
        grid.z_faces,
        solid_T,
        gas_T,
        rear_tube_T,
        rear_gas_T,
        htc,
        front_htc,
        front_qgas,
        gas_rho,
        gas_velocity,
        gas_reynolds,
        cell_dp,
        receiver_dp_mbar,
        dp1_prediction_mbar,
        mass_flow_kg_s,
        sol
    )
end

function _axis_interp_indices(coords::AbstractVector, target)
    target <= coords[1] && return (1, 1, 0.0)
    target >= coords[end] && return (length(coords), length(coords), 0.0)
    hi = searchsortedfirst(coords, target)
    lo = hi - 1
    weight = (target - coords[lo]) / (coords[hi] - coords[lo])
    return (lo, hi, weight)
end

function _sample_solid(result::SimulationResult2D, r_target, z_target; receiver_only=false)
    r_coords = receiver_only ? view(result.r_solid, 1:result.nr_rec) : result.r_solid
    i0, i1, wr = _axis_interp_indices(r_coords, r_target)
    j0, j1, wz = _axis_interp_indices(result.z_solid, z_target)
    nt = length(result.time)
    values = zeros(nt)
    for k in 1:nt
        t00 = result.solid_temperature[i0, j0, k]
        t10 = result.solid_temperature[i1, j0, k]
        t01 = result.solid_temperature[i0, j1, k]
        t11 = result.solid_temperature[i1, j1, k]
        values[k] = (1.0 - wz) * ((1.0 - wr) * t00 + wr * t10) +
                    wz * ((1.0 - wr) * t01 + wr * t11)
    end
    return values
end

function _sample_rear_gas(result::SimulationResult2D, z_target)
    j0, j1, wz = _axis_interp_indices(result.z_rear_gas, z_target)
    return vec((1.0 - wz) .* result.rear_gas_temperature[j0, :] .+
               wz .* result.rear_gas_temperature[j1, :])
end

function sensor_predictions2D(result::SimulationResult2D)
    r_core = 0.0
    r_perim = result.receiver_radius
    r_t2 = result.receiver_radius + 40.0e-3

    T8 = _sample_solid(result, r_perim, 5.0e-3; receiver_only=true)
    T12 = _sample_solid(result, r_perim, 58.0e-3; receiver_only=true)
    T11 = _sample_solid(result, r_perim, 107.0e-3; receiver_only=true)
    T9 = _sample_solid(result, r_core, 58.0e-3; receiver_only=true)
    T10 = _sample_solid(result, r_core, 107.0e-3; receiver_only=true)
    T3 = _sample_rear_gas(result, 3.0e-3)
    T2 = _sample_solid(result, r_t2, 58.0e-3)

    return (T8=T8, T12=T12, T11=T11, T9=T9, T10=T10, T3=T3, T2=T2)
end

function get_t90_2D(times::AbstractVector, signal::AbstractVector)
    length(times) == length(signal) || return NaN
    delta = signal[end] - signal[1]
    abs(delta) <= 1e-6 && return 0.0
    target = signal[1] + 0.90 * delta
    idx = delta > 0.0 ? findfirst(y -> y >= target, signal) :
                        findfirst(y -> y <= target, signal)
    idx === nothing && return times[end] - times[1]
    return times[idx] - times[1]
end

function normalized_slope_mse_2D(model::AbstractVector, experiment::AbstractVector)
    length(model) == length(experiment) || return NaN
    length(model) <= 2 && return 0.0
    span = max(maximum(experiment) - minimum(experiment), 1.0)
    d_model = diff(model)
    d_exp = diff(experiment)
    return mean(abs2, (d_model .- d_exp)) / span^2
end

# ============================================================================
# Parameter Optimization & Calibration Interface
# ============================================================================
using Optimization
using OptimizationNLopt

export pack_parameters2D, unpack_parameters2D, calibrate2D, LB_2D, UB_2D

function pack_parameters2D(p::ModelParameters2D)
    return [
        p.heat_transfer.coefficient,
        p.heat_transfer.reynolds_exponent,
        p.solid.radial_conductivity_scale,
        p.solid.axial_conductivity_scale,
        p.optics.beam_radius_sigma * 1000.0,
        p.optics.spillage_fraction,
        p.optics.scale_456,
        p.optics.scale_304,
        p.optics.scale_256,
        p.solid.rad_extinction_coeff,
        p.solid.receiver_felt_contact_resistance,
        p.optics.extinction_coefficient,
    ]
end

function unpack_parameters2D(theta::AbstractVector, p_base::ModelParameters2D=default_parameters2D())
    g = p_base.geometry
    s = SolidProperties2D(
        density = p_base.solid.density,
        radial_conductivity_scale = theta[3],
        axial_conductivity_scale = theta[4],
        rad_extinction_coeff = theta[10],
        receiver_felt_contact_resistance = theta[11],
        felt_conductivity_ref = p_base.solid.felt_conductivity_ref,
        felt_density = p_base.solid.felt_density,
        felt_heat_capacity = p_base.solid.felt_heat_capacity,
        casing_conductivity = p_base.solid.casing_conductivity,
        casing_density = p_base.solid.casing_density,
        casing_heat_capacity = p_base.solid.casing_heat_capacity,
    )
    ht = HeatTransferParameters2D(
        coefficient = theta[1],
        reynolds_exponent = theta[2],
        prandtl_exponent = p_base.heat_transfer.prandtl_exponent,
        minimum_nusselt = p_base.heat_transfer.minimum_nusselt,
        c_radial_flow = p_base.heat_transfer.c_radial_flow,
        front_coefficient = p_base.heat_transfer.front_coefficient,
        front_reynolds_exponent = p_base.heat_transfer.front_reynolds_exponent,
    )
    l = LossParameters2D(
        front_emissivity = p_base.losses.front_emissivity,
        casing_emissivity = p_base.losses.casing_emissivity,
        front_loss_scale = p_base.losses.front_loss_scale,
        casing_loss_scale = p_base.losses.casing_loss_scale,
        rear_adaptor_conductance = p_base.losses.rear_adaptor_conductance,
        flange_conductance_scale = p_base.losses.flange_conductance_scale,
    )
    opt = OpticalParameters2D(
        absorbed_fraction = p_base.optics.absorbed_fraction,
        extinction_coefficient = theta[12],
        beam_radius_sigma = theta[5] / 1000.0,
        spillage_fraction = theta[6],
        front_deposition_fraction = p_base.optics.front_deposition_fraction,
        scale_456 = theta[7],
        scale_304 = theta[8],
        scale_256 = theta[9],
    )
    return ModelParameters2D(g, s, ht, l, opt, p_base.hydraulics)
end

# HEATING-ONLY optimization bounds
const LB_2D = [1.0e-5, 0.20, 0.005, 0.20, 5.0, 0.00, 0.50, 0.50, 0.30, 20.0, 0.0, 10.0]
const UB_2D = [0.10, 1.80, 0.300, 2.00, 35.0, 0.25, 2.20, 2.20, 1.50, 1000.0, 0.020, 300.0]
# Fit only transport and group-level delivered-power terms. Beam width,
# spillage and receiver/felt contact resistance are fixed independently to
# prevent the strongest structural confounding seen in v1-v7.
const FIT_INDICES_2D = [1, 3, 4, 7, 8, 9, 10, 12]

function loss_function_2D(theta, (heating_cases, solve_fn))
    p = unpack_parameters2D(theta)
    total_loss = 0.0
    count = 0
    
    # Fit ONLY heating cases (Cooling cases used as benchmark)
    for sim_id in heating_cases
        try
            case_data = solve_fn(sim_id, p; is_cooling=false)
            model = case_data.model
            exp = case_data.experiment
            times = case_data.times
            
            all(isfinite, model) || return Inf
            
            for j in 1:7
                scale = max(maximum(exp[:, j]) - minimum(exp[:, j]), 1.0)
                squared_error = mean(abs2, model[:, j] .- exp[:, j])
                nmse_span = squared_error / scale^2
                nmse_100K = squared_error / 100.0^2
                steady_100K = ((model[end, j] - exp[end, j]) / 100.0)^2
                shape = normalized_slope_mse_2D(model[:, j], exp[:, j])
                t90 = ((get_t90_2D(times, model[:, j]) - get_t90_2D(times, exp[:, j])) / max(times[end] - times[1], 1.0))^2
                
                # Priority on deep rear (T10, T11, T3) and T12
                w_sensor = (j in [3, 5, 6]) ? 2.0 : 1.0
                total_loss += w_sensor * (
                    0.70 * nmse_100K +
                    0.30 * nmse_span +
                    0.20 * steady_100K +
                    0.10 * shape +
                    0.05 * t90
                )
            end
            
            # Preserve the two measured core/perimeter offsets without imposing a sign.
            t12_t9_model = model[end, 2] - model[end, 4]
            t12_t9_exp = exp[end, 2] - exp[end, 4]
            t10_t11_model = model[end, 5] - model[end, 3] # T10 - T11 model
            t10_t11_exp = exp[end, 5] - exp[end, 3]     # T10 - T11 exp
            total_loss += ((t12_t9_model - t12_t9_exp)^2 +
                           (t10_t11_model - t10_t11_exp)^2) / 30.0^2
            
            count += 1
        catch err
            err isa InterruptException && rethrow()
            return Inf
        end
    end
    return count > 0 ? total_loss / count : Inf
end

function calibrate2D(solve_fn; heating_cases=("E67", "E76", "E80"),
                     max_iters=150, max_time=1800.0,
                     p_init=default_parameters2D())
    theta_base = pack_parameters2D(p_init)
    theta0 = theta_base[FIT_INDICES_2D]
    evals = Ref(0)
    
    function counted_objective(theta_fit, p_data)
        evals[] += 1
        theta = copy(theta_base)
        theta[FIT_INDICES_2D] .= theta_fit
        val = loss_function_2D(theta, p_data)
        if isfinite(val)
            println("[calibrate2D v10] Heating-Only Eval $(evals[]): loss = $(round(val, digits=4))")
            flush(stdout)
        end
        return val
    end

    opt_prob = OptimizationFunction(counted_objective, SciMLBase.NoAD())
    prob = OptimizationProblem(
        opt_prob, theta0, (heating_cases, solve_fn),
        lb=LB_2D[FIT_INDICES_2D], ub=UB_2D[FIT_INDICES_2D],
    )
    
    println("[calibrate2D v10] Starting Heating-Only Parameter Optimization...")
    sol = solve(prob, OptimizationNLopt.NLopt.LN_BOBYQA(); maxiters=max_iters, maxtime=max_time)
    
    theta_opt = copy(theta_base)
    theta_opt[FIT_INDICES_2D] .= sol.u
    p_opt = unpack_parameters2D(theta_opt)
    return (
        objective=sol.objective,
        parameters=p_opt,
        minimizer=theta_opt,
        fitted_indices=copy(FIT_INDICES_2D),
        evaluations=evals[],
        retcode=sol.retcode,
    )
end

end # module Receiver2D_v10
