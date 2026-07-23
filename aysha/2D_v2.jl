# ============================================================================
# 2D_v2.jl - Multi-Domain 2D Continuum Model with Physical Buoyant Plume & LTNE Physics
# ============================================================================
# Key Physics Additions in v2:
#   1. Churchill-Chu Natural Convection Plume Correlation h_front(T_s, T_amb)
#      dynamically increasing front aperture heat dissipation at high T.
#   2. Temperature-Dependent Alumina Felt Insulation Conductivity k_felt(T)
#      accounting for radiative pore diffusion: k_felt(T) = 0.06 + 1.2e-10 * T^3.
#   3. Radial Flow Recruitment Partition u_z(r) = u_0 * (1 - c_flow * (r/R_rec)^2)
#      enabling preferential core channel cooling.
#   4. High spatial resolution disc discretization (17 x 25 grid, 20 rear tube cells).
# ============================================================================

module Receiver2D_v2

using DifferentialEquations
using OrdinaryDiffEq
using SciMLBase
using LinearAlgebra
using Statistics

export Geometry2D, SolidProperties2D, HeatTransferParameters2D, LossParameters2D
export OpticalParameters2D, ModelParameters2D, OperatingCondition2D, SimulationResult2D
export default_parameters2D, simulate2D, sensor_predictions2D, get_t90_2D, normalized_slope_mse_2D
export linear_history, sic_conductivity, sic_heat_capacity
export air_conductivity, air_heat_capacity, air_viscosity, alumina_conductivity, alumina_heat_capacity
export felt_conductivity_temp, front_nusselt_correlation
export pack_parameters2D, unpack_parameters2D, calibrate2D, LB_2D, UB_2D

const SIGMA = 5.670374419e-8
const WATER_FLANGE_TEMP = 298.15
const G_ACCEL = 9.81

# ============================================================================
# Struct Definitions
# ============================================================================

Base.@kwdef struct Geometry2D
    length::Float64 = 137e-3
    receiver_radius::Float64 = 33.9e-3          # Equivalent disc radius sqrt(A_frt / pi)
    insulation_outer_radius::Float64 = 57.0e-3  # Alumina felt outer radius
    casing_outer_radius::Float64 = 75.0e-3      # Aluminum casing outer radius
    channel_width::Float64 = 1.5e-3             # Hydraulic diameter Dh
    wall_thickness::Float64 = 0.4e-3            # Monolith web wall thickness
    channel_count::Int = 100                    # Total channel count
    nodes_r_rec::Int = 10                       # Radial nodes inside SiC receiver disc
    nodes_r_felt::Int = 5                       # Radial nodes in alumina felt insulation
    nodes_r_case::Int = 2                       # Radial nodes in aluminum casing
    nodes_z::Int = 25                           # High-resolution axial nodes along receiver depth
    rear_tube_length::Float64 = 150.0e-3        # Rear alumina exit tube length
    rear_tube_nodes::Int = 20                   # High-resolution rear tube axial nodes
end

Base.@kwdef struct SolidProperties2D
    density::Float64 = 2150.0                   # Dense SiC solid density [kg/m3] (40g porous mass via area_solid)
    radial_conductivity_scale::Float64 = 0.03   # Effective radial conductivity scale (k_r_eff ~ 3 W/m/K)
    axial_conductivity_scale::Float64 = 1.00    # Axial effective conductivity scale
    heat_capacity_scale::Float64 = 1.00         # Heat capacity scale
    felt_conductivity_ref::Float64 = 0.06       # Baseline alumina felt thermal conductivity [W/m/K]
    felt_density::Float64 = 140.0               # Alumina felt density [kg/m3]
    felt_heat_capacity::Float64 = 1360.0        # Alumina felt heat capacity [J/kg/K]
    casing_conductivity::Float64 = 205.0        # Aluminum casing thermal conductivity [W/m/K]
    casing_density::Float64 = 2700.0            # Aluminum casing density [kg/m3]
    casing_heat_capacity::Float64 = 900.0       # Aluminum casing heat capacity [J/kg/K]
end

Base.@kwdef struct HeatTransferParameters2D
    coefficient::Float64 = 2.50                 # Developing flow Nu prefactor A_Nu
    reynolds_exponent::Float64 = 0.333          # Reynolds exponent B_Re
    prandtl_exponent::Float64 = 0.333           # Prandtl exponent C_Pr
    minimum_nusselt::Float64 = 3.61             # Fully developed laminar Nu floor
    phi_0::Float64 = 1.00                       # Active flow fraction at Re = 50
    m_rec::Float64 = 0.00                       # Flow recruitment exponent
    c_radial_flow::Float64 = 0.10               # Preferential central core channel flow fraction
end

Base.@kwdef struct LossParameters2D
    front_emissivity::Float64 = 0.85            # Front face emissivity
    cavity_emissivity::Float64 = 0.20           # Casing outer emissivity
    cavity_heat_capacity::Float64 = 301.0       # Participating thermal mass [J/K]
    rear_adaptor_conductance::Float64 = 15.0    # Adaptor sleeve contact conductance [W/K]
    flange_conductance_scale::Float64 = 1.00    # Water-cooled flange sink scale
end

Base.@kwdef struct OpticalParameters2D
    absorbed_fraction::Float64 = 0.85            # Eta_abs
    extinction_coefficient::Float64 = 184.67     # Beta_opt [1/m]
    beam_radius_sigma::Float64 = 14.0e-3         # Gaussian beam radial width sigma [m]
    spillage_fraction::Float64 = 0.20            # Fraction of power spilling to perimeter rim
    front_deposition_fraction::Float64 = 0.20    # Fraction absorbed specularly at front surface
    scale_456::Float64 = 1.336                   # Total power scale factor at 456 kW/m2
    scale_304::Float64 = 1.374                   # Total power scale factor at 304 kW/m2
    scale_256::Float64 = 0.786                   # Total power scale factor at 256 kW/m2
end

Base.@kwdef struct ModelParameters2D
    geometry::Geometry2D = Geometry2D()
    solid::SolidProperties2D = SolidProperties2D()
    heat_transfer::HeatTransferParameters2D = HeatTransferParameters2D()
    losses::LossParameters2D = LossParameters2D()
    optics::OpticalParameters2D = OpticalParameters2D()
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
    r_solid::Vector{Float64}
    z_solid::Vector{Float64}
    z_rear::Vector{Float64}
    r_gas::Vector{Float64}
    z_gas::Vector{Float64}
    solid_temperature::Array{Float64, 3}       # (Nr_total, Nz, Nt)
    gas_temperature::Array{Float64, 3}         # (Nr_rec, Nz+1, Nt)
    rear_tube_temperature::Array{Float64, 2}   # (Nz_rear, Nt)
    rear_gas_temperature::Array{Float64, 2}    # (Nz_rear+1, Nt)
    cavity_temperature::Vector{Float64}        # (Nt,)
    heat_transfer_coefficient::Array{Float64, 3} # (Nr_rec, Nz, Nt)
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
    Tc = clamp_temperature(T) - 273.15
    val = 1.0 / (1.10e-4 + 3.42e-6 * Tc)
    return max(5.0, val)
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

function alumina_conductivity(T)
    Tc = clamp_temperature(T) - 273.15
    return max(1.0, 5.5 + 34.5 * exp(-0.0033 * Tc))
end

function alumina_heat_capacity(T)
    Tk = clamp_temperature(T)
    return max(500.0, (1.00446 + 1.742e-4 * Tk - 2.796e4 * (Tk^-2)) * 1000.0)
end

# Temperature-dependent alumina felt thermal conductivity (solid + pore radiation)
function felt_conductivity_temp(T, ref_k=0.06)
    Tk = clamp_temperature(T)
    return max(0.05, ref_k + 1.2e-10 * Tk^3)
end

# Churchill-Chu Natural Convection Plume Correlation at Front Aperture
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
    dz::Float64
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
    volume_cell::Vector{Float64}
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

    dz = L_rec / nz
    dz_rear = L_rear / nz_rear

    # Multi-domain radial faces
    r_faces_rec = collect(range(0.0, R_rec, length=nr_rec + 1))
    r_faces_felt = collect(range(R_rec, R_felt, length=nr_felt + 1))[2:end]
    r_faces_case = collect(range(R_felt, R_case, length=nr_case + 1))[2:end]
    r_faces = vcat(r_faces_rec, r_faces_felt, r_faces_case)

    r_centers = [0.5 * (r_faces[i] + r_faces[i + 1]) for i in 1:nr_total]

    z_faces = collect(range(0.0, L_rec, length=nz + 1))
    z_centers = [0.5 * (z_faces[j] + z_faces[j + 1]) for j in 1:nz]

    z_rear_faces = collect(range(0.0, L_rear, length=nz_rear + 1))
    z_rear_centers = [0.5 * (z_rear_faces[j] + z_rear_faces[j + 1]) for j in 1:nz_rear]

    total_frt = pi * R_rec^2
    area_frt = [pi * (r_faces[i + 1]^2 - r_faces[i]^2) for i in 1:nr_total]

    ch_area = g.channel_width^2
    total_ch_area = g.channel_count * ch_area
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

    volume_cell = area_frt .* dz

    return DiscretizationGrid2D(
        nr_rec, nr_felt, nr_case, nr_total, nz, nz_rear, dz, dz_rear,
        r_faces, r_centers, z_faces, z_centers,
        z_rear_faces, z_rear_centers,
        area_frt, area_solid, area_flow, volume_cell, perimeter_ex,
        porosity, g.channel_width, total_frt
    )
end

function solar_weights2D(grid::DiscretizationGrid2D, p::ModelParameters2D)
    nr_rec, nz = grid.nr_rec, grid.nz
    sigma = p.optics.beam_radius_sigma
    beta_opt = p.optics.extinction_coefficient
    dz = grid.dz
    f_front = p.optics.front_deposition_fraction

    I_rad = [exp(-0.5 * (grid.r_centers[i] / sigma)^2) for i in 1:nr_rec]
    I_rad ./= sum(I_rad[i] * grid.area_frt[i] for i in 1:nr_rec) / grid.total_frt

    w_z = [exp(-beta_opt * (j - 1) * dz) - exp(-beta_opt * j * dz) for j in 1:nz]
    w_z_depth = (1.0 - f_front) .* (w_z ./ sum(w_z))
    w_z_depth[1] += f_front

    weights = zeros(grid.nr_total, nz)
    for i in 1:nr_rec, j in 1:nz
        weights[i, j] = I_rad[i] * (grid.area_frt[i] / grid.total_frt) * w_z_depth[j]
    end
    weights ./= sum(weights)
    return weights
end

# ============================================================================
# Multi-Domain Thermal Properties Helper
# ============================================================================

function cell_thermal_properties2D(i, T, p::ModelParameters2D, grid::DiscretizationGrid2D)
    if i <= grid.nr_rec
        k_val = p.solid.radial_conductivity_scale * sic_conductivity(T)
        rho_val = p.solid.density
        cp_val = p.solid.heat_capacity_scale * sic_heat_capacity(T)
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

function _gas_profile2D!(gas, qgas, hcell, gas_rear, qgas_rear, hrear, solid_temp, T_rear, time, p, op, grid)
    nr_rec, nz = grid.nr_rec, grid.nz
    nz_rear = grid.nz_rear
    flow_lpm = max(0.0, Float64(op.flow_lpm(time)))
    inlet = Float64(op.inlet_temperature(time))
    hp = p.heat_transfer

    mdot_total = 1.199 * flow_lpm / 60000.0
    
    # Active Flow Recruitment phi_act(Re)
    reynolds_ref = 50.0
    re_est = mdot_total * grid.dh / (sum(grid.area_flow[1:nr_rec]) * air_viscosity(inlet))
    phi_act = clamp(hp.phi_0 * (max(re_est, 1.0) / reynolds_ref)^hp.m_rec, 0.20, 1.00)
    mdot_active = phi_act * mdot_total
    dh = grid.dh

    if flow_lpm <= 1e-12
        fill!(qgas, 0.0)
        fill!(hcell, 0.0)
        fill!(qgas_rear, 0.0)
        fill!(hrear, 0.0)
        gas .= inlet
        gas_rear .= inlet
        return nothing
    end

    # Preferential radial channel flow partition: psi(r) = 1 - c_flow * (r/R_rec)^2
    R_rec = grid.r_faces[nr_rec + 1]
    psi_r = [1.0 - hp.c_radial_flow * (grid.r_centers[i] / R_rec)^2 for i in 1:nr_rec]
    psi_sum = sum(psi_r[i] * grid.area_flow[i] for i in 1:nr_rec)

    for i in 1:nr_rec
        mdot_ring = mdot_active * (psi_r[i] * grid.area_flow[i] / psi_sum)
        gas[i, 1] = inlet
        
        for j in 1:nz
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
            ua = hcell[i, j] * grid.perimeter_ex[i] * grid.dz
            effectiveness = -expm1(-ua / max(mdot_ring * cp, eps(Float64)))
            
            gas[i, j + 1] = gas[i, j] + effectiveness * (solid_temp[i, j] - gas[i, j])
            qgas[i, j] = mdot_ring * cp * (gas[i, j + 1] - gas[i, j])
        end
    end

    # Mixed gas profile entering rear exit tube
    total_enthalpy = sum(mdot_active * (psi_r[i] * grid.area_flow[i] / psi_sum) * air_heat_capacity(gas[i, end]) * gas[i, end] for i in 1:nr_rec)
    total_mcp = sum(mdot_active * (psi_r[i] * grid.area_flow[i] / psi_sum) * air_heat_capacity(gas[i, end]) for i in 1:nr_rec)
    gas_rear[1] = total_mcp > 0 ? total_enthalpy / total_mcp : inlet

    # March gas through rear alumina tube
    r_tube_dh = 13.0e-3
    r_tube_area = pi * (13.0e-3 / 2.0)^2
    r_tube_perim = pi * 13.0e-3

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
    T_cavity = view(u, (nr_total * nz + nz_rear + 1):(nr_total * nz + nz_rear + 1))

    du_s = reshape(view(du, 1:(nr_total * nz)), nr_total, nz)
    du_rear = view(du, (nr_total * nz + 1):(nr_total * nz + nz_rear))
    du_cavity = view(du, (nr_total * nz + nz_rear + 1):(nr_total * nz + nz_rear + 1))

    g = p.geometry
    solid = p.solid
    loss = p.losses
    ambient = Float64(op.ambient_temperature(time))

    _gas_profile2D!(work.gas, work.qgas, work.h, work.gas_rear, work.qgas_rear, work.hrear, T_s[1:nr_rec, :], T_rear, time, p, op, grid)
    fill!(du_s, 0.0)
    fill!(du_rear, 0.0)
    du_cavity[1] = 0.0

    # 1. Radial Solid/Felt/Casing Conduction across all radial nodes
    for i in 1:(nr_total - 1), j in 1:nz
        prop_i = cell_thermal_properties2D(i, T_s[i, j], p, grid)
        prop_ip1 = cell_thermal_properties2D(i + 1, T_s[i + 1, j], p, grid)
        kr_face = 2.0 * prop_i.k * prop_ip1.k / (prop_i.k + prop_ip1.k)
        
        r_face = grid.r_faces[i + 1]
        area_r_face = 2.0 * pi * r_face * grid.dz
        dr = grid.r_centers[i + 1] - grid.r_centers[i]
        
        Qcond_r = kr_face * area_r_face * (T_s[i, j] - T_s[i + 1, j]) / dr
        du_s[i, j] -= Qcond_r
        du_s[i + 1, j] += Qcond_r
    end

    # 2. Axial Conduction + Rosseland Radiation
    for i in 1:nr_total, j in 1:(nz - 1)
        prop_j = cell_thermal_properties2D(i, T_s[i, j], p, grid)
        prop_jp1 = cell_thermal_properties2D(i, T_s[i, j + 1], p, grid)
        
        kz_j = prop_j.k + (i <= nr_rec ? 16.0 * SIGMA * T_s[i, j]^3 / (3.0 * p.optics.extinction_coefficient) : 0.0)
        kz_jp1 = prop_jp1.k + (i <= nr_rec ? 16.0 * SIGMA * T_s[i, j + 1]^3 / (3.0 * p.optics.extinction_coefficient) : 0.0)
        kz_face = 2.0 * kz_j * kz_jp1 / (kz_j + kz_jp1)
        
        area_ax = i <= nr_rec ? grid.area_solid[i] : grid.area_frt[i]
        Qcond_z = kz_face * area_ax * (T_s[i, j] - T_s[i, j + 1]) / grid.dz
        du_s[i, j] -= Qcond_z
        du_s[i, j + 1] += Qcond_z
    end

    # 3. Solar Absorption & Heat Exchanges
    absorbed_power = p.optics.absorbed_fraction * max(0.0, Float64(op.irradiance(time))) * grid.total_frt
    spillage_power = p.optics.spillage_fraction * absorbed_power
    core_power = absorbed_power - spillage_power

    for i in 1:nr_total, j in 1:nz
        Qsolar = 0.0
        if i <= nr_rec
            Qsolar = core_power * weights[i, j]
            du_s[i, j] += Qsolar - work.qgas[i, j]
        elseif i <= (nr_rec + grid.nr_felt)
            Qsolar = (spillage_power / (grid.nr_felt * nz))
            du_s[i, j] += Qsolar
        end
    end

    # 4. Front Face Losses with Churchill-Chu Buoyant Plume Correlation
    for i in 1:nr_total
        area_f = i <= nr_rec ? grid.area_solid[i] : grid.area_frt[i]
        h_front_i = front_nusselt_correlation(T_s[i, 1], ambient)
        Qfront_rad = loss.front_emissivity * SIGMA * area_f * (T_s[i, 1]^4 - ambient^4)
        Qfront_conv = h_front_i * area_f * (T_s[i, 1] - ambient)
        du_s[i, 1] -= (Qfront_rad + Qfront_conv)
    end

    # 5. Rear Contact Adaptor & Alumina Tube Heat Sink
    Qrear_exit = 0.0
    for i in 1:nr_rec
        Qrear_i = loss.rear_adaptor_conductance * (grid.area_solid[i] / sum(grid.area_solid[1:nr_rec])) * (T_s[i, nz] - T_rear[1])
        du_s[i, nz] -= Qrear_i
        Qrear_exit += Qrear_i
    end
    du_rear[1] += Qrear_exit

    A_tube_wall = pi * (25.0e-3^2 - 19.0e-3^2) / 4.0
    for j in 1:(nz_rear - 1)
        k_j = alumina_conductivity(T_rear[j])
        k_jp1 = alumina_conductivity(T_rear[j + 1])
        k_face = 2.0 * k_j * k_jp1 / (k_j + k_jp1)
        Qcond_rear = k_face * A_tube_wall * (T_rear[j] - T_rear[j + 1]) / grid.dz_rear
        du_rear[j] -= Qcond_rear
        du_rear[j + 1] += Qcond_rear
    end

    # Rear alumina tube loss to gas & water-cooled flange sink (25°C)
    for j in 1:nz_rear
        Q_gas_exchange = work.qgas_rear[j]
        Q_flange = loss.flange_conductance_scale * (alumina_conductivity(T_rear[j]) * A_tube_wall / grid.dz_rear) * (T_rear[j] - WATER_FLANGE_TEMP)
        du_rear[j] -= (Q_gas_exchange + Q_flange)
    end

    # Casing outer convective & radiative loss to ambient
    for j in 1:nz
        Qcasing_loss = 10.0 * (2.0 * pi * g.casing_outer_radius * grid.dz) * (T_s[nr_total, j] - ambient) +
                       loss.cavity_emissivity * SIGMA * (2.0 * pi * g.casing_outer_radius * grid.dz) * (T_s[nr_total, j]^4 - ambient^4)
        du_s[nr_total, j] -= Qcasing_loss
    end

    # 6. Convert Power Rates to Temperature Rates dT/dt
    # Physical Thermal Mass Correction:
    # Use porous solid volume (area_solid * dz) for SiC matrix cells (i <= nr_rec)
    # so total monolith solid mass matches exact 40.0g measured mass.
    for i in 1:nr_total, j in 1:nz
        prop = cell_thermal_properties2D(i, T_s[i, j], p, grid)
        cell_vol = i <= nr_rec ? (grid.area_solid[i] * grid.dz) : grid.volume_cell[i]
        cap = prop.rho * prop.cp * cell_vol
        du_s[i, j] /= max(cap, eps(Float64))
    end

    for j in 1:nz_rear
        cap_rear = 3950.0 * alumina_heat_capacity(T_rear[j]) * A_tube_wall * grid.dz_rear
        du_rear[j] /= max(cap_rear, eps(Float64))
    end
    
    du_cavity[1] = 0.0
    return nothing
end

mutable struct Workspace2D
    gas::Array{Float64, 2}
    qgas::Array{Float64, 2}
    h::Array{Float64, 2}
    gas_rear::Vector{Float64}
    qgas_rear::Vector{Float64}
    hrear::Vector{Float64}
end

function Workspace2D(grid::DiscretizationGrid2D)
    return Workspace2D(
        zeros(grid.nr_rec, grid.nz + 1),
        zeros(grid.nr_rec, grid.nz),
        zeros(grid.nr_rec, grid.nz),
        zeros(grid.nz_rear + 1),
        zeros(grid.nz_rear),
        zeros(grid.nz_rear)
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
    solver=Rodas5P(autodiff=AutoFiniteDiff()),
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

    u_length = nr_total * nz + nz_rear + 1
    
    if initial_temperature isa AbstractVector
        length(initial_temperature) == u_length || error("initial_temperature vector length mismatch.")
        u_initial = copy(Float64.(initial_temperature))
    else
        u_initial = fill(Float64(initial_temperature), u_length)
    end

    t_span = (Float64(save_times[1]), Float64(save_times[end]))
    prob = ODEProblem(receiver_ode2D!, u_initial, t_span, context)

    sol = solve(prob, solver; saveat=save_times, reltol=reltol, abstol=abstol, dtmax=dtmax)

    nt = length(sol.t)
    solid_T = zeros(nr_total, nz, nt)
    gas_T = zeros(nr_rec, nz + 1, nt)
    rear_tube_T = zeros(nz_rear, nt)
    rear_gas_T = zeros(nz_rear + 1, nt)
    cavity_T = zeros(nt)
    htc = zeros(nr_rec, nz, nt)

    for k in 1:nt
        u_k = sol.u[k]
        t_k = sol.t[k]
        
        solid_T[:, :, k] .= reshape(view(u_k, 1:(nr_total * nz)), nr_total, nz)
        rear_tube_T[:, k] .= view(u_k, (nr_total * nz + 1):(nr_total * nz + nz_rear))
        cavity_T[k] = u_k[nr_total * nz + nz_rear + 1]

        _gas_profile2D!(work.gas, work.qgas, work.h, work.gas_rear, work.qgas_rear, work.hrear, solid_T[1:nr_rec, :, k], rear_tube_T[:, k], t_k, p, op, grid)
        gas_T[:, :, k] .= work.gas
        rear_gas_T[:, k] .= work.gas_rear
        htc[:, :, k] .= work.h
    end

    return SimulationResult2D(
        Vector{Float64}(sol.t),
        grid.r_centers,
        grid.z_centers,
        grid.z_rear_centers,
        grid.r_centers[1:nr_rec],
        grid.z_faces,
        solid_T,
        gas_T,
        rear_tube_T,
        rear_gas_T,
        cavity_T,
        htc,
        sol
    )
end

function sensor_predictions2D(result::SimulationResult2D)
    nr_total = length(result.r_solid)
    nz = length(result.z_solid)
    
    r_core_idx = 1
    r_perim_idx = min(10, nr_total)
    
    # T2 thermocouple is embedded 40 mm outside receiver wall (r = R_rec + 40.0 mm = 73.9 mm)
    r_t2_target = 33.9e-3 + 40.0e-3
    r_t2_idx = argmin(abs.(result.r_solid .- r_t2_target))

    z8_idx = argmin(abs.(result.z_solid .- 11.0e-3))
    z12_idx = argmin(abs.(result.z_solid .- 58.0e-3))
    z11_idx = argmin(abs.(result.z_solid .- 107.0e-3))

    T8 = vec(result.solid_temperature[r_perim_idx, z8_idx, :])
    T12 = vec(result.solid_temperature[r_perim_idx, z12_idx, :])
    T11 = vec(result.solid_temperature[r_perim_idx, z11_idx, :])
    T9 = vec(result.solid_temperature[r_core_idx, z12_idx, :])
    T10 = vec(result.solid_temperature[r_core_idx, z11_idx, :])
    
    # Gas exit probe T3 is located inside rear alumina tube at ~140 mm depth
    z_t3_idx = argmin(abs.(result.z_rear .- (140.0e-3 - 137.0e-3)))
    T3 = vec(result.rear_gas_temperature[z_t3_idx, :])
    
    # T2 thermocouple embedded 40 mm deep inside housing insulation
    T2 = vec(result.solid_temperature[r_t2_idx, z12_idx, :])

    return (T8=T8, T12=T12, T11=T11, T9=T9, T10=T10, T3=T3, T2=T2)
end

function get_t90_2D(times::AbstractVector, signal::AbstractVector)
    length(times) == length(signal) || return NaN
    s_min, s_max = minimum(signal), maximum(signal)
    span = s_max - s_min
    span <= 1e-6 && return 0.0
    target = s_min + 0.90 * span
    idx = findfirst(y -> y >= target, signal)
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
        p.losses.cavity_heat_capacity,
        p.optics.scale_456,
        p.optics.scale_304,
        p.optics.scale_256,
        p.heat_transfer.phi_0,
        p.heat_transfer.m_rec,
        p.optics.front_deposition_fraction,
        p.losses.flange_conductance_scale,
        p.heat_transfer.c_radial_flow,
    ]
end

function unpack_parameters2D(theta::AbstractVector, p_base::ModelParameters2D=default_parameters2D())
    g = p_base.geometry
    s = SolidProperties2D(
        density = p_base.solid.density,
        radial_conductivity_scale = theta[3],
        axial_conductivity_scale = theta[4],
        heat_capacity_scale = p_base.solid.heat_capacity_scale,
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
        phi_0 = theta[11],
        m_rec = theta[12],
        c_radial_flow = theta[15],
    )
    l = LossParameters2D(
        front_emissivity = p_base.losses.front_emissivity,
        cavity_emissivity = p_base.losses.cavity_emissivity,
        cavity_heat_capacity = theta[7],
        rear_adaptor_conductance = p_base.losses.rear_adaptor_conductance,
        flange_conductance_scale = theta[14],
    )
    opt = OpticalParameters2D(
        absorbed_fraction = p_base.optics.absorbed_fraction,
        extinction_coefficient = p_base.optics.extinction_coefficient,
        beam_radius_sigma = theta[5] / 1000.0,
        spillage_fraction = theta[6],
        front_deposition_fraction = theta[13],
        scale_456 = theta[8],
        scale_304 = theta[9],
        scale_256 = theta[10],
    )
    return ModelParameters2D(g, s, ht, l, opt)
end

# Unmodified power scale bounds per user directive [0.50, 2.80]
const LB_2D = [0.5, 0.10, 0.001, 0.10, 5.0, 0.01, 100.0, 0.50, 0.50, 0.40, 0.40, 0.00, 0.00, 0.10, 0.00]
const UB_2D = [6.0, 0.60, 0.100, 3.00, 35.0, 0.40, 500.0, 2.80, 2.80, 2.00, 1.00, 0.60, 0.50, 5.00, 0.50]

function loss_function_2D(theta, (heating_cases, cooling_cases, solve_fn))
    p = unpack_parameters2D(theta)
    total_loss = 0.0
    count = 0
    
    for (sim_id, is_cool) in vcat([(id, false) for id in heating_cases], [(id, true) for id in cooling_cases])
        try
            case_data = solve_fn(sim_id, p; is_cooling=is_cool)
            model = case_data.model
            exp = case_data.experiment
            times = case_data.times
            
            all(isfinite, model) || return Inf
            
            for j in 1:7
                scale = max(maximum(exp[:, j]) - minimum(exp[:, j]), 1.0)
                nmse = mean(abs2, model[:, j] .- exp[:, j]) / scale^2
                shape = normalized_slope_mse_2D(model[:, j], exp[:, j])
                t90 = ((get_t90_2D(times, model[:, j]) - get_t90_2D(times, exp[:, j])) / max(times[end] - times[1], 1.0))^2
                total_loss += nmse + 0.50 * shape + 0.25 * t90
            end
            
            if !is_cool
                t12_t9_diff = (model[end, 2] - model[end, 4]) - (exp[end, 2] - exp[end, 4])
                t11_t10_diff = (model[end, 3] - model[end, 5]) - (exp[end, 3] - exp[end, 5])
                total_loss += 0.25 * (t12_t9_diff^2 + t11_t10_diff^2) / 100.0^2
            end
            
            count += 1
        catch err
            err isa InterruptException && rethrow()
            return Inf
        end
    end
    return count > 0 ? total_loss / count : Inf
end

function calibrate2D(solve_fn; heating_cases=("E67", "E76", "E80"),
                     cooling_cases=("C69", "C80", "C81"),
                     max_iters=150, max_time=600.0,
                     p_init=default_parameters2D())
    theta0 = pack_parameters2D(p_init)
    evals = Ref(0)
    
    function counted_objective(theta, p_data)
        evals[] += 1
        val = loss_function_2D(theta, p_data)
        if isfinite(val)
            println("[calibrate2D v2] Eval $(evals[]): loss = $(round(val, digits=4))")
            flush(stdout)
        end
        return val
    end
    
    opt_prob = OptimizationFunction(counted_objective, SciMLBase.NoAD())
    prob = OptimizationProblem(opt_prob, theta0, (heating_cases, cooling_cases, solve_fn), lb=LB_2D, ub=UB_2D)
    
    println("[calibrate2D v2] Starting NLopt.LN_BOBYQA parameter optimization...")
    sol = solve(prob, OptimizationNLopt.NLopt.LN_BOBYQA(); maxiters=max_iters, maxtime=max_time)
    
    p_opt = unpack_parameters2D(sol.u)
    return (objective=sol.objective, parameters=p_opt, minimizer=sol.u, retcode=sol.retcode)
end

end # module Receiver2D_v2
