# ============================================================================
# 2D_v14.jl - symmetry-reduced square-channel thermal-hydraulic network
# ============================================================================

module Receiver2D_v14

include("2D_v12.jl")
const V12 = Receiver2D_v12
const V11 = V12.Receiver2D_v11

using DifferentialEquations
using OrdinaryDiffEq
using SciMLBase
using Statistics

export NetworkParameters2D, ModelParameters2D, NetworkGrid2D
export SimulationResult2D, default_parameters2D, build_network_grid2D
export network_inventory2D, simulate2D, sensor_predictions2D
export energy_rate_ledger2D, hydraulic_feedback_diagnostic2D
export Geometry2D, SolidProperties2D, HeatTransferParameters2D
export LossParameters2D, OpticalParameters2D, HydraulicParameters2D
export OperatingCondition2D, AssemblyParameters2D, ObservationParameters2D
export build_initial_state_2D, linear_history, get_t90_2D

const Geometry2D = V12.Geometry2D
const SolidProperties2D = V12.SolidProperties2D
const HeatTransferParameters2D = V12.HeatTransferParameters2D
const LossParameters2D = V12.LossParameters2D
const OpticalParameters2D = V12.OpticalParameters2D
const HydraulicParameters2D = V12.HydraulicParameters2D
const OperatingCondition2D = V12.OperatingCondition2D
const AssemblyParameters2D = V12.AssemblyParameters2D
const ObservationParameters2D = V12.ObservationParameters2D
const build_initial_state_2D = V12.build_initial_state_2D
const linear_history = V12.linear_history
const get_t90_2D = V12.get_t90_2D

Base.@kwdef struct NetworkParameters2D
    lateral_conductivity_scale::Float64 = 1.0
    edge_felt_contact_scale::Float64 = 1.0
    groove_quadrature_points::Int = 9
end

Base.@kwdef struct ModelParameters2D
    base::V12.ModelParameters2D = V12.default_parameters2D()
    network::NetworkParameters2D = NetworkParameters2D()
end

function Base.getproperty(p::ModelParameters2D, name::Symbol)
    name === :base && return getfield(p, :base)
    name === :network && return getfield(p, :network)
    return getproperty(getfield(p, :base), name)
end

default_parameters2D() = ModelParameters2D()

struct NetworkGrid2D
    base_grid::V11.DiscretizationGrid2D
    channel_side_count::Int
    group_count::Int
    group_keys::Vector{NTuple{2,Int}}
    group_map::Matrix{Int}
    multiplicity::Vector{Int}
    boundary_faces::Vector{Int}
    centroid_x::Vector{Float64}
    centroid_y::Vector{Float64}
    centroid_radius::Vector{Float64}
    groove_overlap::Vector{Float64}
    neighbor_a::Vector{Int}
    neighbor_b::Vector{Int}
    neighbor_count::Vector{Int}
    pitch::Float64
    solid_area_per_channel::Float64
    frontal_area_per_channel::Float64
    side_group::Int
    core_group::Int
    corner_group::Int
end

function _channel_key(i, j, n)
    x2 = abs(2 * i - n - 1)
    y2 = abs(2 * j - n - 1)
    return (max(x2, y2), min(x2, y2))
end

function _channel_groove_overlap(xc, yc, width, free_radius, nq)
    nq = max(nq, 3)
    count_out = 0
    count_all = nq * nq
    for ix in 1:nq, iy in 1:nq
        x = xc - 0.5 * width + (ix - 0.5) * width / nq
        y = yc - 0.5 * width + (iy - 0.5) * width / nq
        x^2 + y^2 > free_radius^2 && (count_out += 1)
    end
    return count_out / count_all
end

function build_network_grid2D(
    p::ModelParameters2D=default_parameters2D(),
)
    b = p.base.base
    g = b.geometry
    n = round(Int, sqrt(g.channel_count))
    n^2 == g.channel_count ||
        error("v14 requires a square channel count")
    iseven(n) ||
        error("v14 midpoint-side mapping currently requires an even grid")
    pitch = g.receiver_width / n
    isapprox(
        pitch, g.channel_width + g.wall_thickness; rtol=1e-10,
    ) || error("receiver width is inconsistent with channel pitch")

    raw_keys = [
        _channel_key(i, j, n) for i in 1:n for j in 1:n
    ]
    keys = sort(unique(raw_keys); by=key -> (
        key[1]^2 + key[2]^2, key[1], key[2],
    ))
    key_to_group = Dict(key => index
                        for (index, key) in enumerate(keys))
    group_map = zeros(Int, n, n)
    multiplicity = zeros(Int, length(keys))
    boundary_faces = zeros(Int, length(keys))
    for i in 1:n, j in 1:n
        group = key_to_group[_channel_key(i, j, n)]
        group_map[i, j] = group
        multiplicity[group] += 1
        boundary_faces[group] +=
            (i == 1) + (i == n) + (j == 1) + (j == n)
    end

    edge_counts = Dict{Tuple{Int,Int},Int}()
    function add_edge!(ga, gb)
        ga == gb && return
        key = ga < gb ? (ga, gb) : (gb, ga)
        edge_counts[key] = get(edge_counts, key, 0) + 1
    end
    for i in 1:n, j in 1:n
        i < n && add_edge!(group_map[i, j], group_map[i + 1, j])
        j < n && add_edge!(group_map[i, j], group_map[i, j + 1])
    end
    edge_keys = sort(collect(Base.keys(edge_counts)))
    neighbor_a = [key[1] for key in edge_keys]
    neighbor_b = [key[2] for key in edge_keys]
    neighbor_count = [edge_counts[key] for key in edge_keys]

    centroid_x = [0.5 * key[1] * pitch for key in keys]
    centroid_y = [0.5 * key[2] * pitch for key in keys]
    centroid_radius = hypot.(centroid_x, centroid_y)
    free_radius = 0.5 * b.hydraulics.groove_free_diameter
    groove_overlap = [
        _channel_groove_overlap(
            centroid_x[group], centroid_y[group],
            g.channel_width, free_radius,
            p.network.groove_quadrature_points,
        ) for group in eachindex(keys)
    ]
    total_area = g.receiver_width^2
    flow_area = g.channel_count * g.channel_width^2
    solid_area_per_channel =
        (total_area - flow_area) / g.channel_count
    frontal_area_per_channel = total_area / g.channel_count
    side_key = (n - 1, 1)
    core_key = (1, 1)
    corner_key = (n - 1, n - 1)
    return NetworkGrid2D(
        V11.build_grid2D(b), n, length(keys), keys, group_map,
        multiplicity, boundary_faces, centroid_x, centroid_y,
        centroid_radius, groove_overlap, neighbor_a, neighbor_b,
        neighbor_count, pitch, solid_area_per_channel,
        frontal_area_per_channel, key_to_group[side_key],
        key_to_group[core_key], key_to_group[corner_key],
    )
end

function network_inventory2D(
    p::ModelParameters2D=default_parameters2D();
    temperature=800.0,
)
    grid = build_network_grid2D(p)
    b = p.base.base
    total_solid_area = grid.solid_area_per_channel *
                       sum(grid.multiplicity)
    mass = b.solid.density * total_solid_area * b.geometry.length
    capacity = mass * V11.sic_heat_capacity(temperature)
    return (
        group_count=grid.group_count,
        multiplicity_sum=sum(grid.multiplicity),
        boundary_face_sum=sum(grid.boundary_faces),
        receiver_solid_area_m2=total_solid_area,
        receiver_mass_kg=mass,
        receiver_capacity_J_K=capacity,
        flow_area_m2=sum(grid.multiplicity) *
                     b.geometry.channel_width^2,
    )
end

mutable struct NetworkWorkspace2D
    gas::Matrix{Float64}
    qgas::Matrix{Float64}
    h::Matrix{Float64}
    front_qgas::Vector{Float64}
    front_h::Vector{Float64}
    density::Matrix{Float64}
    velocity::Matrix{Float64}
    reynolds::Matrix{Float64}
    prandtl::Matrix{Float64}
    graetz::Matrix{Float64}
    nusselt::Matrix{Float64}
    dp_cell::Matrix{Float64}
    channel_mdot::Vector{Float64}
    channel_dp::Vector{Float64}
    groove_dp::Vector{Float64}
    common_dp::Float64
    equal_pressure_error::Float64
    gas_rear::Vector{Float64}
    qgas_rear::Vector{Float64}
    hrear::Vector{Float64}
end

function NetworkWorkspace2D(grid::NetworkGrid2D)
    ng = grid.group_count
    nz = grid.base_grid.nz
    nr = grid.base_grid.nz_rear
    return NetworkWorkspace2D(
        zeros(ng, nz + 1), zeros(ng, nz), zeros(ng, nz),
        zeros(ng), zeros(ng), zeros(ng, nz), zeros(ng, nz),
        zeros(ng, nz), zeros(ng, nz), zeros(ng, nz),
        zeros(ng, nz), zeros(ng, nz), zeros(ng), zeros(ng),
        zeros(ng), 0.0, 0.0, zeros(nr + 1), zeros(nr),
        zeros(nr),
    )
end

function _network_solar_weights(grid::NetworkGrid2D,
                                p::ModelParameters2D)
    b = p.base.base
    sigma = b.optics.beam_radius_sigma
    radial = exp.(
        -0.5 .* (grid.centroid_radius ./ sigma).^2,
    )
    normalization = sum(
        grid.multiplicity[group] * radial[group]
        for group in 1:grid.group_count
    ) / sum(grid.multiplicity)
    radial ./= normalization
    beta = b.optics.extinction_coefficient
    ffront = b.optics.front_deposition_fraction
    bg = grid.base_grid
    wz = [
        exp(-beta * bg.z_faces[j]) -
        exp(-beta * bg.z_faces[j + 1])
        for j in 1:bg.nz
    ]
    wz .*= (1.0 - ffront)
    wz[1] += ffront
    weights = zeros(grid.group_count, bg.nz)
    for group in 1:grid.group_count, j in 1:bg.nz
        weights[group, j] = radial[group] *
            grid.multiplicity[group] /
            sum(grid.multiplicity) * wz[j]
    end
    return weights
end

function _state_layout(grid::NetworkGrid2D)
    bg = grid.base_grid
    nouter = bg.nr_felt + bg.nr_case
    nch = grid.group_count * bg.nz
    nout = nouter * bg.nz
    ntube = bg.nz_rear
    rch = 1:nch
    rout = (last(rch) + 1):(last(rch) + nout)
    rtube = (last(rout) + 1):(last(rout) + ntube)
    iadaptor = last(rtube) + 1
    ihousing = iadaptor + 1
    return (
        channel=rch, outer=rout, tube=rtube,
        adaptor=iadaptor, housing=ihousing,
        total=ihousing, nouter=nouter,
    )
end

function _march_channel_group2D!(
    work::NetworkWorkspace2D,
    group::Int,
    mdot_channel::Float64,
    inlet::Float64,
    channel_temperature,
    p::ModelParameters2D,
    grid::NetworkGrid2D,
)
    b = p.base.base
    bg = grid.base_grid
    hp = b.heat_transfer
    hyd = b.hydraulics
    width = b.geometry.channel_width
    area = width^2
    perimeter = 4.0 * width
    multiplicity = grid.multiplicity[group]

    cp_in = V11.air_heat_capacity(inlet)
    k_in = V11.air_conductivity(inlet)
    mu_in = V11.air_viscosity(inlet)
    pr_in = cp_in * mu_in / k_in
    re_front = mdot_channel * width / (area * mu_in)
    nu_front = max(0.0, hp.front_coefficient) *
        max(re_front, eps(Float64))^hp.front_reynolds_exponent *
        max(pr_in, eps(Float64))^(1.0 / 3.0)
    work.front_h[group] = nu_front * k_in / width
    front_ntu = work.front_h[group] *
        grid.solid_area_per_channel /
        max(mdot_channel * cp_in, eps(Float64))
    front_effectiveness = -expm1(-front_ntu)
    work.gas[group, 1] = inlet + front_effectiveness *
        (channel_temperature[group, 1] - inlet)
    work.front_qgas[group] = multiplicity * mdot_channel * cp_in *
        (work.gas[group, 1] - inlet)

    for j in 1:bg.nz
        Tgas = V11.clamp_temperature(work.gas[group, j])
        Twall = V11.clamp_temperature(
            channel_temperature[group, j],
        )
        cp = V11.air_heat_capacity(Tgas)
        k_pr = V11.air_conductivity(Tgas)
        k_h = V11.air_conductivity(Twall)
        mu = V11.air_viscosity(Tgas)
        reynolds = mdot_channel * width / (area * mu)
        prandtl = cp * mu / k_pr
        correlation = V11.graetz_nusselt2D(
            reynolds, prandtl, bg.z_centers[j], width, hp;
            gas_temperature=Tgas, wall_temperature=Twall,
        )
        nusselt = correlation.nusselt
        work.prandtl[group, j] = prandtl
        work.graetz[group, j] = correlation.graetz
        work.nusselt[group, j] = nusselt
        work.h[group, j] = nusselt * k_h / width
        ua = work.h[group, j] * perimeter * bg.dz_cells[j]
        effectiveness = -expm1(
            -ua / max(mdot_channel * cp, eps(Float64)),
        )
        work.gas[group, j + 1] = Tgas +
            effectiveness * (Twall - Tgas)
        work.qgas[group, j] =
            multiplicity * mdot_channel * cp *
            (work.gas[group, j + 1] - Tgas)

        Tbulk = V11.clamp_temperature(
            0.5 * (
                work.gas[group, j] +
                work.gas[group, j + 1]
            ),
        )
        rho = V11.air_density(
            Tbulk, hyd.atmospheric_pressure,
        )
        mu_bulk = V11.air_viscosity(Tbulk)
        velocity = mdot_channel / (rho * area)
        work.density[group, j] = rho
        work.velocity[group, j] = velocity
        work.reynolds[group, j] =
            rho * velocity * width / mu_bulk
        work.dp_cell[group, j] =
            hyd.hydraulic_resistance_scale *
            0.5 * V11.SQUARE_DUCT_DARCY_POISEUILLE *
            mu_bulk * bg.dz_cells[j] * velocity / width^2
    end

    rho_out = work.density[group, end]
    velocity_out = work.velocity[group, end]
    work.groove_dp[group] = grid.groove_overlap[group] *
        max(0.0, hyd.groove_loss_coefficient) *
        0.5 * rho_out * velocity_out^2
    work.dp_cell[group, end] += work.groove_dp[group]
    work.channel_dp[group] = sum(view(
        work.dp_cell, group, :,
    ))
    return nothing
end

function _flows_at_pressure2D(
    pressure,
    linear_coeff,
    quadratic_coeff,
)
    flows = similar(linear_coeff)
    for group in eachindex(flows)
        a = max(linear_coeff[group], 0.0)
        b = max(quadratic_coeff[group], 0.0)
        flows[group] = if b > eps(Float64)
            (-a + sqrt(a^2 + 4.0 * b * pressure)) /
            (2.0 * b)
        else
            pressure / max(a, eps(Float64))
        end
    end
    return flows
end

function _equal_pressure_candidate2D(
    linear_coeff,
    quadratic_coeff,
    multiplicity,
    mdot_total,
    pressure_hint,
)
    lo = 0.0
    hi = max(pressure_hint, 1.0)
    for _ in 1:30
        flows = _flows_at_pressure2D(
            hi, linear_coeff, quadratic_coeff,
        )
        sum(multiplicity .* flows) >= mdot_total && break
        hi *= 2.0
    end
    for _ in 1:50
        mid = 0.5 * (lo + hi)
        flows = _flows_at_pressure2D(
            mid, linear_coeff, quadratic_coeff,
        )
        if sum(multiplicity .* flows) < mdot_total
            lo = mid
        else
            hi = mid
        end
    end
    return _flows_at_pressure2D(
        0.5 * (lo + hi), linear_coeff, quadratic_coeff,
    )
end

function _solve_channel_flows2D!(
    work::NetworkWorkspace2D,
    mdot_total::Float64,
    inlet::Float64,
    channel_temperature,
    p::ModelParameters2D,
    grid::NetworkGrid2D,
)
    ng = grid.group_count
    hyd = p.base.base.hydraulics
    work.channel_mdot .= mdot_total / sum(grid.multiplicity)
    relaxation = clamp(
        hyd.equal_pressure_relaxation, 0.05, 1.0,
    )
    tolerance = max(hyd.equal_pressure_tolerance, 1e-10)
    max_iterations = max(hyd.equal_pressure_max_iterations, 24)

    for _ in 1:max_iterations
        for group in 1:ng
            _march_channel_group2D!(
                work, group, work.channel_mdot[group], inlet,
                channel_temperature, p, grid,
            )
        end
        common = sum(
            grid.multiplicity[group] *
            work.channel_mdot[group] *
            work.channel_dp[group]
            for group in 1:ng
        ) / mdot_total
        error = maximum(
            abs(work.channel_dp[group] - common)
            for group in 1:ng
        ) / max(common, eps(Float64))
        work.common_dp = common
        work.equal_pressure_error = error
        error <= tolerance && return nothing

        linear_coeff = zeros(ng)
        quadratic_coeff = zeros(ng)
        for group in 1:ng
            mdot = max(work.channel_mdot[group], eps(Float64))
            quadratic_coeff[group] =
                work.groove_dp[group] / mdot^2
            linear_coeff[group] =
                max(0.0, work.channel_dp[group] -
                         work.groove_dp[group]) / mdot
        end
        candidate = _equal_pressure_candidate2D(
            linear_coeff, quadratic_coeff, grid.multiplicity,
            mdot_total, common,
        )
        work.channel_mdot .=
            (1.0 - relaxation) .* work.channel_mdot .+
            relaxation .* candidate
        work.channel_mdot .*= mdot_total /
            sum(grid.multiplicity .* work.channel_mdot)
    end

    for group in 1:ng
        _march_channel_group2D!(
            work, group, work.channel_mdot[group], inlet,
            channel_temperature, p, grid,
        )
    end
    work.common_dp = sum(
        grid.multiplicity[group] *
        work.channel_mdot[group] *
        work.channel_dp[group]
        for group in 1:ng
    ) / mdot_total
    work.equal_pressure_error = maximum(
        abs(work.channel_dp[group] - work.common_dp)
        for group in 1:ng
    ) / max(work.common_dp, eps(Float64))
    return nothing
end

function _gas_profile_network2D!(
    work::NetworkWorkspace2D,
    channel_temperature,
    tube_temperature,
    time,
    p::ModelParameters2D,
    op::OperatingCondition2D,
    grid::NetworkGrid2D,
)
    b = p.base.base
    bg = grid.base_grid
    flow_lpm = max(0.0, Float64(op.flow_lpm(time)))
    inlet = Float64(op.inlet_temperature(time))
    mdot_total = V11.standard_mass_flow2D(
        flow_lpm, b.hydraulics,
    )
    if mdot_total <= eps(Float64)
        fill!(work.qgas, 0.0)
        fill!(work.h, 0.0)
        fill!(work.front_qgas, 0.0)
        fill!(work.front_h, 0.0)
        fill!(work.velocity, 0.0)
        fill!(work.reynolds, 0.0)
        fill!(work.prandtl, 0.0)
        fill!(work.graetz, 0.0)
        fill!(work.nusselt, 0.0)
        fill!(work.dp_cell, 0.0)
        fill!(work.channel_mdot, 0.0)
        fill!(work.channel_dp, 0.0)
        fill!(work.groove_dp, 0.0)
        fill!(work.qgas_rear, 0.0)
        fill!(work.hrear, 0.0)
        work.common_dp = 0.0
        work.equal_pressure_error = 0.0
        work.gas .= inlet
        work.density .= V11.air_density(
            inlet, b.hydraulics.atmospheric_pressure,
        )
        work.gas_rear .= inlet
        return nothing
    end

    _solve_channel_flows2D!(
        work, mdot_total, inlet, channel_temperature,
        p, grid,
    )
    total_enthalpy = sum(
        grid.multiplicity[group] *
        work.channel_mdot[group] *
        V11.air_heat_capacity(work.gas[group, end]) *
        work.gas[group, end]
        for group in 1:grid.group_count
    )
    total_mcp = sum(
        grid.multiplicity[group] *
        work.channel_mdot[group] *
        V11.air_heat_capacity(work.gas[group, end])
        for group in 1:grid.group_count
    )
    work.gas_rear[1] = total_mcp > 0.0 ?
        total_enthalpy / total_mcp : inlet

    radius = b.geometry.rear_tube_inner_radius
    dh = 2.0 * radius
    area = pi * radius^2
    perimeter = 2.0 * pi * radius
    for j in 1:bg.nz_rear
        Tfilm = V11.clamp_temperature(
            0.5 * (
                tube_temperature[j] + work.gas_rear[j]
            ),
        )
        cp = V11.air_heat_capacity(Tfilm)
        kf = V11.air_conductivity(Tfilm)
        mu = V11.air_viscosity(Tfilm)
        reynolds = mdot_total * dh / (area * mu)
        prandtl = cp * mu / kf
        nusselt = reynolds < 2300.0 ?
            3.66 : 0.023 * reynolds^0.8 * prandtl^0.4
        work.hrear[j] = nusselt * kf / dh
        ua = work.hrear[j] * perimeter * bg.dz_rear
        effectiveness = -expm1(
            -ua / max(mdot_total * cp, eps(Float64)),
        )
        work.gas_rear[j + 1] = work.gas_rear[j] +
            effectiveness * (
                tube_temperature[j] - work.gas_rear[j]
            )
        work.qgas_rear[j] = mdot_total * cp *
            (work.gas_rear[j + 1] - work.gas_rear[j])
    end
    return nothing
end

function _outer_conductivity2D(
    outer_index, temperature, p::ModelParameters2D,
    grid::NetworkGrid2D,
)
    bg = grid.base_grid
    if outer_index <= bg.nr_felt
        return V12.felt_conductivity_v12(
            temperature,
            p.base.assembly.felt_conductivity_scale,
        )
    end
    return p.base.base.solid.casing_conductivity
end

function _outer_capacity2D(
    outer_index, j, temperature, p::ModelParameters2D,
    grid::NetworkGrid2D,
)
    bg = grid.base_grid
    original_index = bg.nr_rec + outer_index
    volume = bg.area_frt[original_index] * bg.dz_cells[j]
    if outer_index <= bg.nr_felt
        return p.base.base.solid.felt_density *
            p.base.base.solid.felt_heat_capacity *
            p.base.assembly.felt_heat_capacity_scale * volume
    end
    return p.base.base.solid.casing_density *
        p.base.base.solid.casing_heat_capacity * volume
end

function _housing_capacity2D(p::ModelParameters2D)
    b = p.base.base
    a = p.base.assembly
    extension = max(
        0.0, a.enclosure_internal_length - b.geometry.length,
    )
    annulus_area = pi * (
        a.enclosure_outer_radius^2 -
        a.enclosure_inner_radius^2
    )
    backplate_area = pi * (
        a.enclosure_outer_radius^2 -
        b.geometry.rear_tube_inner_radius^2
    )
    volume = annulus_area * extension +
             backplate_area * a.backplate_thickness
    return a.aluminum_density * a.aluminum_heat_capacity *
           volume
end

function _network_ode2D!(du, u, context, time)
    p = context.p
    grid = context.grid
    bg = grid.base_grid
    work = context.work
    weights = context.weights
    layout = context.layout
    b = p.base.base
    a = p.base.assembly
    ng = grid.group_count
    nz = bg.nz

    Tch = reshape(
        view(u, layout.channel), ng, nz,
    )
    Tout = reshape(
        view(u, layout.outer), layout.nouter, nz,
    )
    Ttube = view(u, layout.tube)
    Tad = u[layout.adaptor]
    Thousing = u[layout.housing]
    duch = reshape(
        view(du, layout.channel), ng, nz,
    )
    duout = reshape(
        view(du, layout.outer), layout.nouter, nz,
    )
    dutube = view(du, layout.tube)
    fill!(du, 0.0)

    _gas_profile_network2D!(
        work, Tch, Ttube, time, p, context.op, grid,
    )
    ambient = Float64(context.op.ambient_temperature(time))

    # Axial conduction within every channel orbit.
    for group in 1:ng, j in 1:(nz - 1)
        dz = 0.5 * (bg.dz_cells[j] + bg.dz_cells[j + 1])
        k0 = b.solid.axial_conductivity_scale *
             V11.sic_conductivity(Tch[group, j]) +
             16.0 * V11.SIGMA *
             V11.clamp_temperature(Tch[group, j])^3 /
             (3.0 * max(1.0, b.solid.rad_extinction_coeff))
        k1 = b.solid.axial_conductivity_scale *
             V11.sic_conductivity(Tch[group, j + 1]) +
             16.0 * V11.SIGMA *
             V11.clamp_temperature(Tch[group, j + 1])^3 /
             (3.0 * max(1.0, b.solid.rad_extinction_coeff))
        kface = 2.0 * k0 * k1 / (k0 + k1)
        area = grid.multiplicity[group] *
               grid.solid_area_per_channel
        Q = kface * area *
            (Tch[group, j] - Tch[group, j + 1]) / dz
        duch[group, j] -= Q
        duch[group, j + 1] += Q
    end

    # Exact square-grid neighbor topology, aggregated by symmetry orbit.
    for edge in eachindex(grid.neighbor_a), j in 1:nz
        ga = grid.neighbor_a[edge]
        gb = grid.neighbor_b[edge]
        k0 = b.solid.radial_conductivity_scale *
             V11.sic_conductivity(Tch[ga, j])
        k1 = b.solid.radial_conductivity_scale *
             V11.sic_conductivity(Tch[gb, j])
        kface = p.network.lateral_conductivity_scale *
                2.0 * k0 * k1 / (k0 + k1)
        Q = grid.neighbor_count[edge] * kface *
            bg.dz_cells[j] *
            (Tch[ga, j] - Tch[gb, j])
        duch[ga, j] -= Q
        duch[gb, j] += Q
    end

    delivered = b.optics.absorbed_fraction *
                max(0.0, Float64(context.op.irradiance(time))) *
                b.geometry.receiver_width^2
    spillage = b.optics.spillage_fraction * delivered
    core_power = delivered - spillage
    for group in 1:ng, j in 1:nz
        duch[group, j] +=
            core_power * weights[group, j] -
            work.qgas[group, j]
    end
    for group in 1:ng
        duch[group, 1] -= work.front_qgas[group]
        area = grid.multiplicity[group] *
               grid.solid_area_per_channel
        Qrad = b.losses.front_loss_scale *
            b.losses.front_emissivity * V11.SIGMA * area *
            (
                V11.clamp_temperature(Tch[group, 1])^4 -
                ambient^4
            )
        duch[group, 1] -= Qrad
    end

    # Exterior square faces contact the felt.
    total_boundary_faces = sum(grid.boundary_faces)
    inner_felt = 1
    dr_felt = bg.r_faces[bg.nr_rec + 2] -
              bg.r_faces[bg.nr_rec + 1]
    for group in 1:ng
        grid.boundary_faces[group] == 0 && continue
        face_fraction =
            grid.boundary_faces[group] / total_boundary_faces
        for j in 1:nz
            krec = b.solid.radial_conductivity_scale *
                   V11.sic_conductivity(Tch[group, j])
            kfelt = _outer_conductivity2D(
                inner_felt, Tout[inner_felt, j], p, grid,
            )
            resistance =
                0.5 * grid.pitch / krec +
                0.5 * dr_felt / kfelt +
                max(0.0,
                    b.solid.receiver_felt_contact_resistance)
            area = face_fraction *
                   4.0 * b.geometry.receiver_width *
                   bg.dz_cells[j]
            Q = p.network.edge_felt_contact_scale * area *
                (Tch[group, j] - Tout[inner_felt, j]) /
                max(resistance, eps(Float64))
            duch[group, j] -= Q
            duout[inner_felt, j] += Q
        end
    end

    # Radial and axial transport through felt/casing.
    for outer in 1:(layout.nouter - 1), j in 1:nz
        i = bg.nr_rec + outer
        dr0 = bg.r_faces[i + 1] - bg.r_faces[i]
        dr1 = bg.r_faces[i + 2] - bg.r_faces[i + 1]
        k0 = _outer_conductivity2D(
            outer, Tout[outer, j], p, grid,
        )
        k1 = _outer_conductivity2D(
            outer + 1, Tout[outer + 1, j], p, grid,
        )
        resistance = 0.5 * dr0 / k0 + 0.5 * dr1 / k1
        area = 2.0 * pi * bg.r_faces[i + 1] *
               bg.dz_cells[j]
        Q = area * (Tout[outer, j] - Tout[outer + 1, j]) /
            max(resistance, eps(Float64))
        duout[outer, j] -= Q
        duout[outer + 1, j] += Q
    end
    for outer in 1:layout.nouter, j in 1:(nz - 1)
        dz = 0.5 * (bg.dz_cells[j] + bg.dz_cells[j + 1])
        k0 = _outer_conductivity2D(
            outer, Tout[outer, j], p, grid,
        )
        k1 = _outer_conductivity2D(
            outer, Tout[outer, j + 1], p, grid,
        )
        kface = 2.0 * k0 * k1 / (k0 + k1)
        original = bg.nr_rec + outer
        Q = kface * bg.area_frt[original] *
            (Tout[outer, j] - Tout[outer, j + 1]) / dz
        duout[outer, j] -= Q
        duout[outer, j + 1] += Q
    end

    felt_front_area = sum(
        bg.area_frt[bg.nr_rec + outer]
        for outer in 1:bg.nr_felt
    )
    for outer in 1:layout.nouter
        original = bg.nr_rec + outer
        if outer <= bg.nr_felt
            duout[outer, 1] +=
                spillage * bg.area_frt[original] /
                felt_front_area
        end
        area = bg.area_frt[original]
        Qrad = b.losses.front_emissivity * V11.SIGMA * area *
            (
                V11.clamp_temperature(Tout[outer, 1])^4 -
                ambient^4
            )
        Qconv = V11.front_nusselt_correlation(
            Tout[outer, 1], ambient,
        ) * area * (Tout[outer, 1] - ambient)
        duout[outer, 1] -=
            b.losses.front_loss_scale * (Qrad + Qconv)
    end

    case_outer = layout.nouter
    for j in 1:nz
        area = 2.0 * pi * b.geometry.casing_outer_radius *
               bg.dz_cells[j]
        Qloss = b.losses.casing_loss_scale * area * (
            10.0 * (Tout[case_outer, j] - ambient) +
            b.losses.casing_emissivity * V11.SIGMA *
            (
                V11.clamp_temperature(
                    Tout[case_outer, j],
                )^4 - ambient^4
            )
        )
        duout[case_outer, j] -= Qloss
    end

    # Loose receiver/adaptor support on the exterior square faces.
    adaptor_power = 0.0
    overlap_start = b.geometry.length -
                    a.adaptor_receiver_overlap
    for group in 1:ng
        grid.boundary_faces[group] == 0 && continue
        face_fraction =
            grid.boundary_faces[group] / total_boundary_faces
        for j in 1:nz
            overlap = V12._overlap_length(
                bg.z_faces[j], bg.z_faces[j + 1],
                overlap_start, b.geometry.length,
            )
            overlap <= 0.0 && continue
            area = face_fraction * 4.0 *
                   b.geometry.receiver_width * overlap
            Q = a.receiver_adaptor_contact_h * area *
                (Tch[group, j] - Tad)
            duch[group, j] -= Q
            adaptor_power += Q
        end
    end

    # Rear tube conduction, gas exchange, adaptor contact and flange.
    tube_area = pi * (
        b.geometry.rear_tube_outer_radius^2 -
        b.geometry.rear_tube_inner_radius^2
    )
    for j in 1:(bg.nz_rear - 1)
        k0 = V11.alumina_conductivity(Ttube[j])
        k1 = V11.alumina_conductivity(Ttube[j + 1])
        kface = 2.0 * k0 * k1 / (k0 + k1)
        Q = kface * tube_area *
            (Ttube[j] - Ttube[j + 1]) / bg.dz_rear
        dutube[j] -= Q
        dutube[j + 1] += Q
    end
    for j in 1:bg.nz_rear
        dutube[j] -= work.qgas_rear[j]
        overlap = V12._overlap_length(
            bg.z_rear_faces[j], bg.z_rear_faces[j + 1],
            0.0, a.adaptor_tube_overlap,
        )
        if overlap > 0.0
            area = 2.0 * pi *
                   b.geometry.rear_tube_outer_radius * overlap
            Q = a.adaptor_tube_contact_h * area *
                (Ttube[j] - Tad)
            dutube[j] -= Q
            adaptor_power += Q
        end
    end
    Qflange = b.losses.flange_conductance_scale *
        V11.alumina_conductivity(Ttube[end]) * tube_area /
        (0.5 * bg.dz_rear) *
        (Ttube[end] - V11.WATER_FLANGE_TEMP)
    dutube[end] -= Qflange

    # Adaptor/felt path.
    kfelt_ad = V12.felt_conductivity_v12(
        Tad, a.felt_conductivity_scale,
    )
    Gadfelt = a.adaptor_felt_conductance_scale *
        2.0 * pi * kfelt_ad * a.adaptor_length /
        log(a.enclosure_inner_radius / a.adaptor_outer_radius)
    rear_weights = [
        V12._overlap_length(
            bg.z_faces[j], bg.z_faces[j + 1],
            overlap_start, b.geometry.length,
        ) for j in 1:nz
    ]
    weight_sum = sum(rear_weights)
    if weight_sum > 0.0
        Tfelt = sum(
            rear_weights[j] * Tout[inner_felt, j]
            for j in 1:nz
        ) / weight_sum
        Q = Gadfelt * (Tad - Tfelt)
        adaptor_power -= Q
        for j in 1:nz
            duout[inner_felt, j] +=
                Q * rear_weights[j] / weight_sum
        end
    end

    # Rear aluminum sleeve/backplate.
    extension = max(
        0.0, a.enclosure_internal_length - b.geometry.length,
    )
    annulus_area = pi * (
        a.enclosure_outer_radius^2 -
        a.enclosure_inner_radius^2
    )
    Gcase = a.aluminum_conductivity * annulus_area /
            max(0.5 * extension, 1e-3)
    Qcase = Gcase * (
        Tout[case_outer, end] - Thousing
    )
    duout[case_outer, end] -= Qcase
    area_loss =
        2.0 * pi * a.enclosure_outer_radius * extension +
        pi * a.enclosure_outer_radius^2
    Qhousing_loss = a.housing_extension_loss_scale *
        area_loss * (
            10.0 * (Thousing - ambient) +
            b.losses.casing_emissivity * V11.SIGMA *
            (
                V11.clamp_temperature(Thousing)^4 -
                ambient^4
            )
        )

    # Convert accumulated powers to temperature derivatives.
    for group in 1:ng, j in 1:nz
        capacity = grid.multiplicity[group] *
            b.solid.density *
            V11.sic_heat_capacity(Tch[group, j]) *
            grid.solid_area_per_channel * bg.dz_cells[j]
        duch[group, j] /= max(capacity, eps(Float64))
    end
    for outer in 1:layout.nouter, j in 1:nz
        duout[outer, j] /= max(
            _outer_capacity2D(
                outer, j, Tout[outer, j], p, grid,
            ), eps(Float64),
        )
    end
    for j in 1:bg.nz_rear
        capacity = 3900.0 *
            V11.alumina_heat_capacity(Ttube[j]) *
            tube_area * bg.dz_rear
        dutube[j] /= max(capacity, eps(Float64))
    end
    Cad = a.adaptor_density * V12.adaptor_volume2D(p.base) *
          V11.alumina_heat_capacity(Tad)
    du[layout.adaptor] =
        adaptor_power / max(Cad, eps(Float64))
    du[layout.housing] =
        (Qcase - Qhousing_loss) /
        max(_housing_capacity2D(p), eps(Float64))
    return nothing
end

function _initial_network_state2D(
    initial_temperature,
    p::ModelParameters2D,
    grid::NetworkGrid2D,
)
    layout = _state_layout(grid)
    initial_temperature isa Number &&
        return fill(Float64(initial_temperature), layout.total)
    initial = Float64.(initial_temperature)
    length(initial) == layout.total && return copy(initial)

    bg = grid.base_grid
    nbase = bg.nr_total * bg.nz + bg.nz_rear
    length(initial) in (nbase, nbase + 2) ||
        error("initial_temperature vector length mismatch in v14")
    Tbase = reshape(
        view(initial, 1:(bg.nr_total * bg.nz)),
        bg.nr_total, bg.nz,
    )
    state = fill(mean(initial), layout.total)
    Tch = reshape(
        view(state, layout.channel), grid.group_count, bg.nz,
    )
    for group in 1:grid.group_count
        radius = min(
            grid.centroid_radius[group],
            bg.r_centers[bg.nr_rec],
        )
        i0, i1, weight = V11._axis_interp_indices(
            bg.r_centers[1:bg.nr_rec], radius,
        )
        for j in 1:bg.nz
            Tch[group, j] =
                (1.0 - weight) * Tbase[i0, j] +
                weight * Tbase[i1, j]
        end
    end
    Tout = reshape(
        view(state, layout.outer), layout.nouter, bg.nz,
    )
    for outer in 1:layout.nouter
        Tout[outer, :] .= Tbase[bg.nr_rec + outer, :]
    end
    state[layout.tube] .= view(
        initial,
        (bg.nr_total * bg.nz + 1):nbase,
    )
    if length(initial) == nbase + 2
        state[layout.adaptor] = initial[nbase + 1]
        state[layout.housing] = initial[nbase + 2]
    else
        state[layout.adaptor] = mean((
            Tch[grid.side_group, end],
            state[first(layout.tube)],
        ))
        state[layout.housing] = Tout[end, end]
    end
    return state
end

struct SimulationResult2D
    time::Vector{Float64}
    group_keys::Vector{NTuple{2,Int}}
    multiplicity::Vector{Int}
    boundary_faces::Vector{Int}
    centroid_radius::Vector{Float64}
    groove_overlap_fraction::Vector{Float64}
    z_solid::Vector{Float64}
    r_outer::Vector{Float64}
    z_rear::Vector{Float64}
    channel_temperature::Array{Float64,3}
    outer_temperature::Array{Float64,3}
    gas_temperature::Array{Float64,3}
    rear_tube_temperature::Matrix{Float64}
    rear_gas_temperature::Matrix{Float64}
    heat_transfer_coefficient::Array{Float64,3}
    gas_density::Array{Float64,3}
    gas_velocity::Array{Float64,3}
    gas_reynolds::Array{Float64,3}
    gas_prandtl::Array{Float64,3}
    gas_graetz::Array{Float64,3}
    gas_nusselt::Array{Float64,3}
    cell_pressure_drop::Array{Float64,3}
    channel_mass_flow_kg_s::Matrix{Float64}
    group_mass_flow_kg_s::Matrix{Float64}
    channel_pressure_drop_Pa::Matrix{Float64}
    groove_pressure_drop_Pa::Matrix{Float64}
    equal_pressure_relative_error::Vector{Float64}
    receiver_pressure_drop_mbar::Vector{Float64}
    dp1_prediction_mbar::Vector{Float64}
    mass_flow_kg_s::Vector{Float64}
    adaptor_temperature::Vector{Float64}
    housing_extension_temperature::Vector{Float64}
    side_group::Int
    core_group::Int
    corner_group::Int
    parameters::ModelParameters2D
    operating_condition::OperatingCondition2D
    ode_solution::Any
end

function simulate2D(
    p::ModelParameters2D,
    op::OperatingCondition2D,
    save_times::AbstractVector;
    initial_temperature=Float64(
        op.ambient_temperature(save_times[1]),
    ),
    solver=FBDF(autodiff=AutoFiniteDiff()),
    reltol=1e-6,
    abstol=1e-7,
    dtmax=30.0,
)
    grid = build_network_grid2D(p)
    layout = _state_layout(grid)
    work = NetworkWorkspace2D(grid)
    weights = _network_solar_weights(grid, p)
    context = (
        p=p, op=op, grid=grid, layout=layout,
        work=work, weights=weights,
    )
    initial = _initial_network_state2D(
        initial_temperature, p, grid,
    )
    problem = ODEProblem(
        _network_ode2D!, initial,
        (Float64(first(save_times)), Float64(last(save_times))),
        context,
    )
    solution = solve(
        problem, solver; saveat=save_times,
        reltol=reltol, abstol=abstol, dtmax=dtmax,
        isoutofdomain=(u, p, t) -> any(
            x -> isnan(x) || x < 100.0 || x > 3500.0, u,
        ),
    )

    bg = grid.base_grid
    nt = length(solution.t)
    ng = grid.group_count
    nchT = zeros(ng, bg.nz, nt)
    noutT = zeros(layout.nouter, bg.nz, nt)
    gasT = zeros(ng, bg.nz + 1, nt)
    tubeT = zeros(bg.nz_rear, nt)
    rearGasT = zeros(bg.nz_rear + 1, nt)
    h = zeros(ng, bg.nz, nt)
    density = zeros(ng, bg.nz, nt)
    velocity = zeros(ng, bg.nz, nt)
    reynolds = zeros(ng, bg.nz, nt)
    prandtl = zeros(ng, bg.nz, nt)
    graetz = zeros(ng, bg.nz, nt)
    nusselt = zeros(ng, bg.nz, nt)
    dp_cell = zeros(ng, bg.nz, nt)
    channel_mdot = zeros(ng, nt)
    group_mdot = zeros(ng, nt)
    channel_dp = zeros(ng, nt)
    groove_dp = zeros(ng, nt)
    equal_error = zeros(nt)
    receiver_dp = zeros(nt)
    dp1 = zeros(nt)
    mass_flow = zeros(nt)
    adaptor = zeros(nt)
    housing = zeros(nt)
    for index in 1:nt
        state = solution.u[index]
        nchT[:, :, index] .= reshape(
            view(state, layout.channel), ng, bg.nz,
        )
        noutT[:, :, index] .= reshape(
            view(state, layout.outer), layout.nouter, bg.nz,
        )
        tubeT[:, index] .= view(state, layout.tube)
        adaptor[index] = state[layout.adaptor]
        housing[index] = state[layout.housing]
        _gas_profile_network2D!(
            work, view(nchT, :, :, index),
            view(tubeT, :, index), solution.t[index],
            p, op, grid,
        )
        gasT[:, :, index] .= work.gas
        rearGasT[:, index] .= work.gas_rear
        h[:, :, index] .= work.h
        density[:, :, index] .= work.density
        velocity[:, :, index] .= work.velocity
        reynolds[:, :, index] .= work.reynolds
        prandtl[:, :, index] .= work.prandtl
        graetz[:, :, index] .= work.graetz
        nusselt[:, :, index] .= work.nusselt
        dp_cell[:, :, index] .= work.dp_cell
        channel_mdot[:, index] .= work.channel_mdot
        group_mdot[:, index] .=
            grid.multiplicity .* work.channel_mdot
        channel_dp[:, index] .= work.channel_dp
        groove_dp[:, index] .= work.groove_dp
        equal_error[index] = work.equal_pressure_error
        receiver_dp[index] = work.common_dp / 100.0
        dp1[index] = p.hydraulics.dp1_zero_offset_mbar +
                     receiver_dp[index]
        flow_lpm = max(0.0, Float64(op.flow_lpm(
            solution.t[index],
        )))
        mass_flow[index] = V11.standard_mass_flow2D(
            flow_lpm, p.hydraulics,
        )
    end
    outer_indices = (
        bg.nr_rec + 1
    ):bg.nr_total
    return SimulationResult2D(
        Float64.(solution.t), copy(grid.group_keys),
        copy(grid.multiplicity), copy(grid.boundary_faces),
        copy(grid.centroid_radius), copy(grid.groove_overlap),
        copy(bg.z_centers), copy(bg.r_centers[outer_indices]),
        copy(bg.z_rear_centers), nchT, noutT, gasT, tubeT,
        rearGasT, h, density, velocity, reynolds, prandtl,
        graetz, nusselt, dp_cell, channel_mdot, group_mdot,
        channel_dp, groove_dp, equal_error, receiver_dp, dp1,
        mass_flow, adaptor, housing, grid.side_group,
        grid.core_group, grid.corner_group, p, op, solution,
    )
end

function _sample_channel2D(result, group, z)
    j0, j1, weight = V11._axis_interp_indices(
        result.z_solid, z,
    )
    return vec(
        (1.0 - weight) .*
        result.channel_temperature[group, j0, :] .+
        weight .* result.channel_temperature[group, j1, :],
    )
end

function _sample_channel_gas2D(result, group, z)
    bulk = 0.5 .* (
        result.gas_temperature[group, 1:end-1, :] .+
        result.gas_temperature[group, 2:end, :]
    )
    j0, j1, weight = V11._axis_interp_indices(
        result.z_solid, z,
    )
    return vec(
        (1.0 - weight) .* bulk[j0, :] .+
        weight .* bulk[j1, :],
    )
end

function _sample_outer2D(result, r, z)
    i0, i1, wr = V11._axis_interp_indices(
        result.r_outer, r,
    )
    j0, j1, wz = V11._axis_interp_indices(
        result.z_solid, z,
    )
    values = zeros(length(result.time))
    for index in eachindex(values)
        t00 = result.outer_temperature[i0, j0, index]
        t10 = result.outer_temperature[i1, j0, index]
        t01 = result.outer_temperature[i0, j1, index]
        t11 = result.outer_temperature[i1, j1, index]
        values[index] =
            (1.0 - wz) * ((1.0 - wr) * t00 + wr * t10) +
            wz * ((1.0 - wr) * t01 + wr * t11)
    end
    return values
end

function _sample_rear_gas2D(result, z)
    zfaces = collect(range(
        0.0,
        result.parameters.geometry.rear_tube_length,
        length=size(result.rear_gas_temperature, 1),
    ))
    j0, j1, weight = V11._axis_interp_indices(zfaces, z)
    return vec(
        (1.0 - weight) .*
        result.rear_gas_temperature[j0, :] .+
        weight .* result.rear_gas_temperature[j1, :],
    )
end

function sensor_predictions2D(result::SimulationResult2D)
    obs = result.parameters.observation
    time = result.time
    side11 = _sample_channel2D(
        result, result.side_group, V11.SIDE_TC_FRONT_Z_2D,
    )
    side58 = _sample_channel2D(
        result, result.side_group, 58e-3,
    )
    side107 = _sample_channel2D(
        result, result.side_group, 107e-3,
    )
    wall58 = _sample_channel2D(
        result, result.core_group, 58e-3,
    )
    wall107 = _sample_channel2D(
        result, result.core_group, 107e-3,
    )
    gas58 = _sample_channel_gas2D(
        result, result.core_group, 58e-3,
    )
    gas107 = _sample_channel_gas2D(
        result, result.core_group, 107e-3,
    )
    wall_fraction = clamp(
        obs.interior_wall_fraction, 0.0, 1.0,
    )
    target9 = wall_fraction .* wall58 +
              (1.0 - wall_fraction) .* gas58
    target10 = wall_fraction .* wall107 +
               (1.0 - wall_fraction) .* gas107
    T3raw = _sample_rear_gas2D(result, 3e-3)
    rT2 = result.parameters.geometry.receiver_radius + 40e-3
    T2 = _sample_outer2D(result, rT2, 58e-3)
    return (
        T8=V12._filter_observation(
            time, side11, obs.side_time_constant_s,
        ),
        T12=V12._filter_observation(
            time, side58, obs.side_time_constant_s,
        ),
        T11=V12._filter_observation(
            time, side107, obs.side_time_constant_s,
        ),
        T9=V12._filter_observation(
            time, target9, obs.interior_time_constant_s,
        ),
        T10=V12._filter_observation(
            time, target10, obs.interior_time_constant_s,
        ),
        T3=V12._filter_observation(
            time, T3raw, obs.outlet_time_constant_s,
        ),
        T2=T2, T9_wall=wall58, T10_wall=wall107,
        T9_gas=gas58, T10_gas=gas107,
    )
end

function energy_rate_ledger2D(
    u::AbstractVector,
    p::ModelParameters2D,
    op::OperatingCondition2D,
    time::Real=0.0,
)
    grid = build_network_grid2D(p)
    layout = _state_layout(grid)
    length(u) == layout.total ||
        error("state vector length mismatch in v14 energy ledger")
    work = NetworkWorkspace2D(grid)
    weights = _network_solar_weights(grid, p)
    context = (
        p=p, op=op, grid=grid, layout=layout,
        work=work, weights=weights,
    )
    du = zeros(Float64, length(u))
    _network_ode2D!(du, u, context, Float64(time))
    bg = grid.base_grid
    b = p.base.base
    a = p.base.assembly
    Tch = reshape(
        view(u, layout.channel), grid.group_count, bg.nz,
    )
    duch = reshape(
        view(du, layout.channel), grid.group_count, bg.nz,
    )
    Tout = reshape(
        view(u, layout.outer), layout.nouter, bg.nz,
    )
    duout = reshape(
        view(du, layout.outer), layout.nouter, bg.nz,
    )
    Ttube = view(u, layout.tube)
    dutube = view(du, layout.tube)
    denergy = 0.0
    for group in 1:grid.group_count, j in 1:bg.nz
        capacity = grid.multiplicity[group] *
            b.solid.density *
            V11.sic_heat_capacity(Tch[group, j]) *
            grid.solid_area_per_channel * bg.dz_cells[j]
        denergy += capacity * duch[group, j]
    end
    for outer in 1:layout.nouter, j in 1:bg.nz
        denergy += _outer_capacity2D(
            outer, j, Tout[outer, j], p, grid,
        ) * duout[outer, j]
    end
    tube_area = pi * (
        b.geometry.rear_tube_outer_radius^2 -
        b.geometry.rear_tube_inner_radius^2
    )
    for j in 1:bg.nz_rear
        capacity = 3900.0 *
            V11.alumina_heat_capacity(Ttube[j]) *
            tube_area * bg.dz_rear
        denergy += capacity * dutube[j]
    end
    Tad = u[layout.adaptor]
    Cad = a.adaptor_density * V12.adaptor_volume2D(p.base) *
          V11.alumina_heat_capacity(Tad)
    denergy += Cad * du[layout.adaptor]
    denergy += _housing_capacity2D(p) * du[layout.housing]

    ambient = Float64(op.ambient_temperature(time))
    receiver_front = sum(
        b.losses.front_loss_scale *
        b.losses.front_emissivity * V11.SIGMA *
        grid.multiplicity[group] *
        grid.solid_area_per_channel *
        (
            V11.clamp_temperature(Tch[group, 1])^4 -
            ambient^4
        ) for group in 1:grid.group_count
    )
    outer_front = 0.0
    for outer in 1:layout.nouter
        original = bg.nr_rec + outer
        area = bg.area_frt[original]
        outer_front += b.losses.front_loss_scale * area * (
            b.losses.front_emissivity * V11.SIGMA *
            (
                V11.clamp_temperature(Tout[outer, 1])^4 -
                ambient^4
            ) +
            V11.front_nusselt_correlation(
                Tout[outer, 1], ambient,
            ) * (Tout[outer, 1] - ambient)
        )
    end
    casing_loss = 0.0
    for j in 1:bg.nz
        area = 2.0 * pi * b.geometry.casing_outer_radius *
               bg.dz_cells[j]
        Tcase = Tout[end, j]
        casing_loss += b.losses.casing_loss_scale * area * (
            10.0 * (Tcase - ambient) +
            b.losses.casing_emissivity * V11.SIGMA *
            (
                V11.clamp_temperature(Tcase)^4 -
                ambient^4
            )
        )
    end
    flange_loss = b.losses.flange_conductance_scale *
        V11.alumina_conductivity(Ttube[end]) * tube_area /
        (0.5 * bg.dz_rear) *
        (Ttube[end] - V11.WATER_FLANGE_TEMP)
    extension = max(
        0.0, a.enclosure_internal_length - b.geometry.length,
    )
    housing_area =
        2.0 * pi * a.enclosure_outer_radius * extension +
        pi * a.enclosure_outer_radius^2
    Thousing = u[layout.housing]
    housing_loss = a.housing_extension_loss_scale *
        housing_area * (
            10.0 * (Thousing - ambient) +
            b.losses.casing_emissivity * V11.SIGMA *
            (
                V11.clamp_temperature(Thousing)^4 -
                ambient^4
            )
        )
    delivered = b.optics.absorbed_fraction *
                max(0.0, Float64(op.irradiance(time))) *
                b.geometry.receiver_width^2
    spillage = b.optics.spillage_fraction * delivered
    core = delivered - spillage
    solar_deposited = spillage + core * sum(weights)
    gas_receiver = sum(work.qgas)
    gas_front = sum(work.front_qgas)
    gas_rear = sum(work.qgas_rear)
    expected = solar_deposited - gas_receiver - gas_front -
               gas_rear - receiver_front - outer_front -
               casing_loss - flange_loss - housing_loss
    return (
        denergy_dt=denergy,
        expected_denergy_dt=expected,
        residual=denergy - expected,
        solar_deposited=solar_deposited,
        gas_receiver=gas_receiver,
        gas_front=gas_front,
        gas_rear=gas_rear,
        receiver_front_loss=receiver_front,
        outer_front_loss=outer_front,
        casing_loss=casing_loss,
        flange_loss=flange_loss,
        housing_loss=housing_loss,
    )
end

function hydraulic_feedback_diagnostic2D(
    p::ModelParameters2D=default_parameters2D();
    flow_lpm=9.0,
    base_temperature=700.0,
    hot_increment=300.0,
)
    grid = build_network_grid2D(p)
    bg = grid.base_grid
    op = OperatingCondition2D(
        irradiance=0.0,
        flow_lpm=flow_lpm,
        inlet_temperature=300.0,
        ambient_temperature=300.0,
    )
    Ttube = fill(300.0, bg.nz_rear)
    uniform = fill(
        Float64(base_temperature), grid.group_count, bg.nz,
    )
    work0 = NetworkWorkspace2D(grid)
    _gas_profile_network2D!(
        work0, uniform, Ttube, 0.0, p, op, grid,
    )
    base_flows = copy(work0.channel_mdot)
    hot = copy(uniform)
    hot[grid.side_group, :] .+= hot_increment
    work1 = NetworkWorkspace2D(grid)
    _gas_profile_network2D!(
        work1, hot, Ttube, 0.0, p, op, grid,
    )
    return (
        side_group=grid.side_group,
        core_group=grid.core_group,
        corner_group=grid.corner_group,
        baseline_side_flow_kg_s=
            base_flows[grid.side_group],
        heated_side_flow_kg_s=
            work1.channel_mdot[grid.side_group],
        heated_to_baseline_side_flow_ratio=
            work1.channel_mdot[grid.side_group] /
            base_flows[grid.side_group],
        baseline_core_to_side_flow_ratio=
            base_flows[grid.core_group] /
            base_flows[grid.side_group],
        heated_core_to_side_flow_ratio=
            work1.channel_mdot[grid.core_group] /
            work1.channel_mdot[grid.side_group],
        baseline_equal_pressure_error=
            work0.equal_pressure_error,
        heated_equal_pressure_error=
            work1.equal_pressure_error,
        baseline_total_mass_flow=sum(
            grid.multiplicity .* base_flows,
        ),
        heated_total_mass_flow=sum(
            grid.multiplicity .* work1.channel_mdot,
        ),
        core_groove_overlap=
            grid.groove_overlap[grid.core_group],
        side_groove_overlap=
            grid.groove_overlap[grid.side_group],
        corner_groove_overlap=
            grid.groove_overlap[grid.corner_group],
    )
end

end
