# ============================================================================
# 2D_v15.jl - v14 square-channel network with an explicit exterior SiC skin
# ============================================================================

module Receiver2D_v15

include("2D_v14.jl")
const V14 = Receiver2D_v14
const V12 = V14.V12
const V11 = V14.V11

using DifferentialEquations
using OrdinaryDiffEq
using SciMLBase

export SkinParameters2D, ModelParameters2D, SimulationResult2D
export default_parameters2D, build_network_grid2D, network_inventory2D
export simulate2D, sensor_predictions2D, energy_rate_ledger2D
export hydraulic_feedback_diagnostic2D
export Geometry2D, SolidProperties2D, HeatTransferParameters2D
export LossParameters2D, OpticalParameters2D, HydraulicParameters2D
export OperatingCondition2D, AssemblyParameters2D, ObservationParameters2D
export NetworkParameters2D, build_initial_state_2D, linear_history, get_t90_2D

const Geometry2D = V14.Geometry2D
const SolidProperties2D = V14.SolidProperties2D
const HeatTransferParameters2D = V14.HeatTransferParameters2D
const LossParameters2D = V14.LossParameters2D
const OpticalParameters2D = V14.OpticalParameters2D
const HydraulicParameters2D = V14.HydraulicParameters2D
const OperatingCondition2D = V14.OperatingCondition2D
const AssemblyParameters2D = V14.AssemblyParameters2D
const ObservationParameters2D = V14.ObservationParameters2D
const NetworkParameters2D = V14.NetworkParameters2D
const build_initial_state_2D = V14.build_initial_state_2D
const linear_history = V14.linear_history
const get_t90_2D = V14.get_t90_2D

Base.@kwdef struct SkinParameters2D
    thickness::Float64 = 0.20e-3
    channel_conductance_scale::Float64 = 0.30
    felt_contact_scale::Float64 = 3.00
end

Base.@kwdef struct ModelParameters2D
    base::V14.ModelParameters2D = V14.default_parameters2D()
    skin::SkinParameters2D = SkinParameters2D()
end

function Base.getproperty(p::ModelParameters2D, name::Symbol)
    name === :base && return getfield(p, :base)
    name === :skin && return getfield(p, :skin)
    return getproperty(getfield(p, :base), name)
end

default_parameters2D() = ModelParameters2D()

build_network_grid2D(p::ModelParameters2D=default_parameters2D()) =
    V14.build_network_grid2D(p.base)

function _skin_area2D(p::ModelParameters2D)
    width = p.geometry.receiver_width
    thickness = p.skin.thickness
    0.0 < thickness < 0.5 * width ||
        error("v15 skin thickness must lie between zero and W/2")
    return width^2 - (width - 2.0 * thickness)^2
end

function _skin_group_areas2D(p::ModelParameters2D, grid)
    total_faces = sum(grid.boundary_faces)
    area = _skin_area2D(p)
    areas = area .* grid.boundary_faces ./ total_faces
    for group in eachindex(areas)
        full = grid.multiplicity[group] *
               grid.solid_area_per_channel
        areas[group] < full ||
            error("v15 skin partition exhausts boundary orbit $group")
    end
    return areas
end

function network_inventory2D(
    p::ModelParameters2D=default_parameters2D();
    temperature=800.0,
)
    grid = build_network_grid2D(p)
    inherited = V14.network_inventory2D(
        p.base; temperature=temperature,
    )
    skin_areas = _skin_group_areas2D(p, grid)
    residual = sum(
        grid.multiplicity[group] *
        grid.solid_area_per_channel - skin_areas[group]
        for group in 1:grid.group_count
    )
    total = residual + sum(skin_areas)
    mass = p.solid.density * total * p.geometry.length
    return (
        group_count=inherited.group_count,
        multiplicity_sum=inherited.multiplicity_sum,
        boundary_face_sum=inherited.boundary_face_sum,
        receiver_solid_area_m2=total,
        residual_channel_area_m2=residual,
        skin_area_m2=sum(skin_areas),
        receiver_mass_kg=mass,
        receiver_capacity_J_K=
            mass * V11.sic_heat_capacity(temperature),
        inherited_receiver_mass_kg=inherited.receiver_mass_kg,
        inherited_receiver_capacity_J_K=
            inherited.receiver_capacity_J_K,
        flow_area_m2=inherited.flow_area_m2,
    )
end

function _state_layout(grid)
    base = V14._state_layout(grid)
    skin = (base.total + 1):(base.total + grid.base_grid.nz)
    return (
        base=base, base_state=1:base.total, skin=skin,
        total=last(skin),
    )
end

function _initial_state2D(initial_temperature, p, grid)
    layout = _state_layout(grid)
    if !(initial_temperature isa Number)
        initial = Float64.(initial_temperature)
        length(initial) == layout.total && return copy(initial)
    end
    base_state = V14._initial_network_state2D(
        initial_temperature, p.base, grid,
    )
    state = zeros(layout.total)
    state[layout.base_state] .= base_state
    Tch = reshape(
        view(base_state, layout.base.channel),
        grid.group_count, grid.base_grid.nz,
    )
    state[layout.skin] .= view(Tch, grid.side_group, :)
    return state
end

function _channel_capacity2D(T, area, dz, p)
    return p.solid.density * V11.sic_heat_capacity(T) * area * dz
end

function _axial_sic_conductivity2D(T, p)
    return p.solid.axial_conductivity_scale *
           V11.sic_conductivity(T) +
           16.0 * V11.SIGMA * V11.clamp_temperature(T)^3 /
           (3.0 * max(1.0, p.solid.rad_extinction_coeff))
end

function _old_edge_exchange2D!(
    Pch, du_base, Tch, Tout, p, grid, layout,
)
    bg = grid.base_grid
    inner_felt = 1
    total_faces = sum(grid.boundary_faces)
    dr_felt = bg.r_faces[bg.nr_rec + 2] -
              bg.r_faces[bg.nr_rec + 1]
    for group in 1:grid.group_count
        grid.boundary_faces[group] == 0 && continue
        fraction = grid.boundary_faces[group] / total_faces
        for j in 1:bg.nz
            krec = p.solid.radial_conductivity_scale *
                   V11.sic_conductivity(Tch[group, j])
            kfelt = V14._outer_conductivity2D(
                inner_felt, Tout[inner_felt, j], p.base, grid,
            )
            resistance = 0.5 * grid.pitch / krec +
                         0.5 * dr_felt / kfelt +
                         max(
                             0.0,
                             p.solid.receiver_felt_contact_resistance,
                         )
            area = fraction * 4.0 * p.geometry.receiver_width *
                   bg.dz_cells[j]
            Qold = p.network.edge_felt_contact_scale * area *
                   (Tch[group, j] - Tout[inner_felt, j]) /
                   max(resistance, eps(Float64))
            Pch[group, j] += Qold
            Cout = V14._outer_capacity2D(
                inner_felt, j, Tout[inner_felt, j], p.base, grid,
            )
            du_base[
                first(layout.outer) +
                (j - 1) * layout.nouter + inner_felt - 1
            ] -= Qold / max(Cout, eps(Float64))
        end
    end
    return nothing
end

function _old_adaptor_exchange2D!(
    Pch, du_base, Tch, Tad, p, grid, layout,
)
    bg = grid.base_grid
    a = p.base.base.assembly
    total_faces = sum(grid.boundary_faces)
    overlap_start = p.geometry.length - a.adaptor_receiver_overlap
    Qtotal = 0.0
    for group in 1:grid.group_count
        grid.boundary_faces[group] == 0 && continue
        fraction = grid.boundary_faces[group] / total_faces
        for j in 1:bg.nz
            overlap = V12._overlap_length(
                bg.z_faces[j], bg.z_faces[j + 1],
                overlap_start, p.geometry.length,
            )
            overlap <= 0.0 && continue
            area = fraction * 4.0 * p.geometry.receiver_width *
                   overlap
            Qold = a.receiver_adaptor_contact_h * area *
                   (Tch[group, j] - Tad)
            Pch[group, j] += Qold
            Qtotal += Qold
        end
    end
    Cad = a.adaptor_density * V12.adaptor_volume2D(p.base.base) *
          V11.alumina_heat_capacity(Tad)
    du_base[layout.adaptor] -=
        Qtotal / max(Cad, eps(Float64))
    return nothing
end

function _skin_ode2D!(du, u, context, time)
    p = context.p
    grid = context.grid
    layout = context.layout
    bg = grid.base_grid
    base_layout = layout.base
    base_u = view(u, layout.base_state)
    base_du = view(du, layout.base_state)
    fill!(du, 0.0)
    V14._network_ode2D!(
        base_du, base_u, context.base_context, time,
    )

    ng = grid.group_count
    nz = bg.nz
    Tch = reshape(
        view(base_u, base_layout.channel), ng, nz,
    )
    duch = reshape(
        view(base_du, base_layout.channel), ng, nz,
    )
    Tout = reshape(
        view(base_u, base_layout.outer),
        base_layout.nouter, nz,
    )
    Tskin = view(u, layout.skin)
    duskin = view(du, layout.skin)
    Tad = base_u[base_layout.adaptor]
    ambient = Float64(context.op.ambient_temperature(time))
    skin_group_area = context.skin_group_area
    skin_area = sum(skin_group_area)

    # Recover the v14 orbit powers before applying the conserved partition.
    Pch = zeros(ng, nz)
    for group in 1:ng, j in 1:nz
        full_area = grid.multiplicity[group] *
                    grid.solid_area_per_channel
        Cfull = _channel_capacity2D(
            Tch[group, j], full_area, bg.dz_cells[j], p,
        )
        Pch[group, j] = Cfull * duch[group, j]
    end
    Pskin = zeros(nz)

    # Remove the v14 direct channel/felt and channel/adaptor paths.
    _old_edge_exchange2D!(
        Pch, base_du, Tch, Tout, p, grid, base_layout,
    )
    _old_adaptor_exchange2D!(
        Pch, base_du, Tch, Tad, p, grid, base_layout,
    )

    # Partition centered solar deposition and the front radiating area.
    delivered = p.optics.absorbed_fraction *
                max(0.0, Float64(context.op.irradiance(time))) *
                p.geometry.receiver_width^2
    core_power = delivered *
                 (1.0 - p.optics.spillage_fraction)
    for group in 1:ng
        skin_group_area[group] == 0.0 && continue
        full_area = grid.multiplicity[group] *
                    grid.solid_area_per_channel
        fraction = skin_group_area[group] / full_area
        for j in 1:nz
            Qsolar = core_power *
                     context.weights[group, j] * fraction
            Pch[group, j] -= Qsolar
            Pskin[j] += Qsolar
        end
        Qrad_full = p.losses.front_loss_scale *
            p.losses.front_emissivity * V11.SIGMA * full_area *
            (
                V11.clamp_temperature(Tch[group, 1])^4 -
                ambient^4
            )
        Qrad_shift = fraction * Qrad_full
        Pch[group, 1] += Qrad_shift
        Pskin[1] -= p.losses.front_loss_scale *
            p.losses.front_emissivity * V11.SIGMA *
            skin_group_area[group] *
            (
                V11.clamp_temperature(Tskin[1])^4 -
                ambient^4
            )
    end

    # Replace the full-area v14 axial conductance by separate residual-solid
    # and skin conductances.
    for group in 1:ng
        skin_group_area[group] == 0.0 && continue
        full_area = grid.multiplicity[group] *
                    grid.solid_area_per_channel
        core_area = full_area - skin_group_area[group]
        for j in 1:(nz - 1)
            dz = 0.5 * (bg.dz_cells[j] + bg.dz_cells[j + 1])
            k0 = _axial_sic_conductivity2D(Tch[group, j], p)
            k1 = _axial_sic_conductivity2D(Tch[group, j + 1], p)
            kface = 2.0 * k0 * k1 / (k0 + k1)
            Qfull = kface * full_area *
                    (Tch[group, j] - Tch[group, j + 1]) / dz
            Qcore = kface * core_area *
                    (Tch[group, j] - Tch[group, j + 1]) / dz
            Pch[group, j] += Qfull - Qcore
            Pch[group, j + 1] -= Qfull - Qcore
        end
    end
    for j in 1:(nz - 1)
        dz = 0.5 * (bg.dz_cells[j] + bg.dz_cells[j + 1])
        k0 = _axial_sic_conductivity2D(Tskin[j], p)
        k1 = _axial_sic_conductivity2D(Tskin[j + 1], p)
        kface = 2.0 * k0 * k1 / (k0 + k1)
        Q = kface * skin_area *
            (Tskin[j] - Tskin[j + 1]) / dz
        Pskin[j] -= Q
        Pskin[j + 1] += Q
    end

    # Continuous-SiC spreading between exterior orbits and the skin.
    total_faces = sum(grid.boundary_faces)
    spreading_distance = 0.5 * grid.pitch
    for group in 1:ng
        grid.boundary_faces[group] == 0 && continue
        fraction = grid.boundary_faces[group] / total_faces
        for j in 1:nz
            k0 = p.solid.radial_conductivity_scale *
                 V11.sic_conductivity(Tch[group, j])
            k1 = p.solid.radial_conductivity_scale *
                 V11.sic_conductivity(Tskin[j])
            kface = 2.0 * k0 * k1 / (k0 + k1)
            area = fraction * 4.0 * p.geometry.receiver_width *
                   bg.dz_cells[j]
            Q = p.skin.channel_conductance_scale * kface * area *
                (Tch[group, j] - Tskin[j]) /
                spreading_distance
            Pch[group, j] -= Q
            Pskin[j] += Q
        end
    end

    # Exterior skin to felt, replacing the removed v14 edge path.
    inner_felt = 1
    dr_felt = bg.r_faces[bg.nr_rec + 2] -
              bg.r_faces[bg.nr_rec + 1]
    for j in 1:nz
        kskin = p.solid.radial_conductivity_scale *
                V11.sic_conductivity(Tskin[j])
        kfelt = V14._outer_conductivity2D(
            inner_felt, Tout[inner_felt, j], p.base, grid,
        )
        resistance = 0.5 * grid.pitch / kskin +
                     0.5 * dr_felt / kfelt +
                     max(
                         0.0,
                         p.solid.receiver_felt_contact_resistance,
                     )
        area = 4.0 * p.geometry.receiver_width * bg.dz_cells[j]
        Q = p.skin.felt_contact_scale * area *
            (Tskin[j] - Tout[inner_felt, j]) /
            max(resistance, eps(Float64))
        Pskin[j] -= Q
        Cout = V14._outer_capacity2D(
            inner_felt, j, Tout[inner_felt, j], p.base, grid,
        )
        outer_index = first(base_layout.outer) +
                      (j - 1) * base_layout.nouter
        base_du[outer_index] += Q / max(Cout, eps(Float64))
    end

    # Exterior skin to the solid adaptor over the measured rear overlap.
    a = p.base.base.assembly
    overlap_start = p.geometry.length - a.adaptor_receiver_overlap
    Qadaptor = 0.0
    for j in 1:nz
        overlap = V12._overlap_length(
            bg.z_faces[j], bg.z_faces[j + 1],
            overlap_start, p.geometry.length,
        )
        overlap <= 0.0 && continue
        area = 4.0 * p.geometry.receiver_width * overlap
        Q = a.receiver_adaptor_contact_h * area *
            (Tskin[j] - Tad)
        Pskin[j] -= Q
        Qadaptor += Q
    end
    Cad = a.adaptor_density * V12.adaptor_volume2D(p.base.base) *
          V11.alumina_heat_capacity(Tad)
    base_du[base_layout.adaptor] +=
        Qadaptor / max(Cad, eps(Float64))

    # Convert the partitioned SiC powers to temperature derivatives.
    for group in 1:ng, j in 1:nz
        full_area = grid.multiplicity[group] *
                    grid.solid_area_per_channel
        core_area = full_area - skin_group_area[group]
        Ccore = _channel_capacity2D(
            Tch[group, j], core_area, bg.dz_cells[j], p,
        )
        duch[group, j] = Pch[group, j] /
                         max(Ccore, eps(Float64))
    end
    for j in 1:nz
        Cskin = _channel_capacity2D(
            Tskin[j], skin_area, bg.dz_cells[j], p,
        )
        duskin[j] = Pskin[j] / max(Cskin, eps(Float64))
    end
    return nothing
end

struct SimulationResult2D
    base_result::V14.SimulationResult2D
    skin_temperature::Matrix{Float64}
    parameters::ModelParameters2D
    ode_solution::Any
end

function Base.getproperty(result::SimulationResult2D, name::Symbol)
    name === :base_result && return getfield(result, :base_result)
    name === :skin_temperature &&
        return getfield(result, :skin_temperature)
    name === :parameters && return getfield(result, :parameters)
    name === :ode_solution && return getfield(result, :ode_solution)
    return getproperty(getfield(result, :base_result), name)
end

function _extract_base_result2D(solution, p, op, grid, layout)
    bg = grid.base_grid
    base_layout = layout.base
    work = V14.NetworkWorkspace2D(grid)
    nt = length(solution.t)
    ng = grid.group_count
    nchT = zeros(ng, bg.nz, nt)
    noutT = zeros(base_layout.nouter, bg.nz, nt)
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
        state = view(solution.u[index], layout.base_state)
        nchT[:, :, index] .= reshape(
            view(state, base_layout.channel), ng, bg.nz,
        )
        noutT[:, :, index] .= reshape(
            view(state, base_layout.outer),
            base_layout.nouter, bg.nz,
        )
        tubeT[:, index] .= view(state, base_layout.tube)
        adaptor[index] = state[base_layout.adaptor]
        housing[index] = state[base_layout.housing]
        V14._gas_profile_network2D!(
            work, view(nchT, :, :, index),
            view(tubeT, :, index), solution.t[index],
            p.base, op, grid,
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
        flow_lpm = max(
            0.0, Float64(op.flow_lpm(solution.t[index])),
        )
        mass_flow[index] = V11.standard_mass_flow2D(
            flow_lpm, p.hydraulics,
        )
    end
    outer_indices = (bg.nr_rec + 1):bg.nr_total
    return V14.SimulationResult2D(
        Float64.(solution.t), copy(grid.group_keys),
        copy(grid.multiplicity), copy(grid.boundary_faces),
        copy(grid.centroid_radius), copy(grid.groove_overlap),
        copy(bg.z_centers), copy(bg.r_centers[outer_indices]),
        copy(bg.z_rear_centers), nchT, noutT, gasT, tubeT,
        rearGasT, h, density, velocity, reynolds, prandtl,
        graetz, nusselt, dp_cell, channel_mdot, group_mdot,
        channel_dp, groove_dp, equal_error, receiver_dp, dp1,
        mass_flow, adaptor, housing, grid.side_group,
        grid.core_group, grid.corner_group, p.base, op, solution,
    )
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
    work = V14.NetworkWorkspace2D(grid)
    weights = V14._network_solar_weights(grid, p.base)
    base_context = (
        p=p.base, op=op, grid=grid, layout=layout.base,
        work=work, weights=weights,
    )
    context = (
        p=p, op=op, grid=grid, layout=layout,
        base_context=base_context, weights=weights,
        skin_group_area=_skin_group_areas2D(p, grid),
    )
    initial = _initial_state2D(initial_temperature, p, grid)
    problem = ODEProblem(
        _skin_ode2D!, initial,
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
    base_result = _extract_base_result2D(
        solution, p, op, grid, layout,
    )
    skin = hcat((
        collect(view(state, layout.skin))
        for state in solution.u
    )...)
    return SimulationResult2D(base_result, skin, p, solution)
end

function _sample_skin2D(result, z)
    j0, j1, weight = V11._axis_interp_indices(
        result.z_solid, z,
    )
    return vec(
        (1.0 - weight) .* result.skin_temperature[j0, :] .+
        weight .* result.skin_temperature[j1, :],
    )
end

function sensor_predictions2D(result::SimulationResult2D)
    obs = result.parameters.observation
    time = result.time
    side5 = _sample_skin2D(result, 5e-3)
    side58 = _sample_skin2D(result, 58e-3)
    side107 = _sample_skin2D(result, 107e-3)
    wall58 = V14._sample_channel2D(
        result, result.core_group, 58e-3,
    )
    wall107 = V14._sample_channel2D(
        result, result.core_group, 107e-3,
    )
    gas58 = V14._sample_channel_gas2D(
        result, result.core_group, 58e-3,
    )
    gas107 = V14._sample_channel_gas2D(
        result, result.core_group, 107e-3,
    )
    wall_fraction = clamp(
        obs.interior_wall_fraction, 0.0, 1.0,
    )
    target9 = wall_fraction .* wall58 +
              (1.0 - wall_fraction) .* gas58
    target10 = wall_fraction .* wall107 +
               (1.0 - wall_fraction) .* gas107
    T3raw = V14._sample_rear_gas2D(result, 3e-3)
    rT2 = result.parameters.geometry.receiver_radius + 40e-3
    T2 = V14._sample_outer2D(result, rT2, 58e-3)
    return (
        T8=V12._filter_observation(
            time, side5, obs.side_time_constant_s,
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
        T2=T2, T8_skin=side5, T12_skin=side58,
        T11_skin=side107, T9_wall=wall58, T10_wall=wall107,
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
        error("state vector length mismatch in v15 energy ledger")
    weights = V14._network_solar_weights(grid, p.base)
    work = V14.NetworkWorkspace2D(grid)
    base_context = (
        p=p.base, op=op, grid=grid, layout=layout.base,
        work=work, weights=weights,
    )
    context = (
        p=p, op=op, grid=grid, layout=layout,
        base_context=base_context, weights=weights,
        skin_group_area=_skin_group_areas2D(p, grid),
    )
    du = zeros(length(u))
    _skin_ode2D!(du, u, context, Float64(time))
    bg = grid.base_grid
    base_u = view(u, layout.base_state)
    base_du = view(du, layout.base_state)
    Tch = reshape(
        view(base_u, layout.base.channel),
        grid.group_count, bg.nz,
    )
    duch = reshape(
        view(base_du, layout.base.channel),
        grid.group_count, bg.nz,
    )
    Tskin = view(u, layout.skin)
    duskin = view(du, layout.skin)
    skin_group_area = context.skin_group_area
    denergy = 0.0
    for group in 1:grid.group_count, j in 1:bg.nz
        full_area = grid.multiplicity[group] *
                    grid.solid_area_per_channel
        core_area = full_area - skin_group_area[group]
        denergy += _channel_capacity2D(
            Tch[group, j], core_area, bg.dz_cells[j], p,
        ) * duch[group, j]
    end
    skin_area = sum(skin_group_area)
    for j in 1:bg.nz
        denergy += _channel_capacity2D(
            Tskin[j], skin_area, bg.dz_cells[j], p,
        ) * duskin[j]
    end
    Tout = reshape(
        view(base_u, layout.base.outer),
        layout.base.nouter, bg.nz,
    )
    duout = reshape(
        view(base_du, layout.base.outer),
        layout.base.nouter, bg.nz,
    )
    for outer in 1:layout.base.nouter, j in 1:bg.nz
        denergy += V14._outer_capacity2D(
            outer, j, Tout[outer, j], p.base, grid,
        ) * duout[outer, j]
    end
    tube_area = pi * (
        p.geometry.rear_tube_outer_radius^2 -
        p.geometry.rear_tube_inner_radius^2
    )
    Ttube = view(base_u, layout.base.tube)
    dutube = view(base_du, layout.base.tube)
    for j in 1:bg.nz_rear
        C = 3900.0 * V11.alumina_heat_capacity(Ttube[j]) *
            tube_area * bg.dz_rear
        denergy += C * dutube[j]
    end
    a = p.base.base.assembly
    Tad = base_u[layout.base.adaptor]
    Cad = a.adaptor_density * V12.adaptor_volume2D(p.base.base) *
          V11.alumina_heat_capacity(Tad)
    denergy += Cad * base_du[layout.base.adaptor]
    denergy += V14._housing_capacity2D(p.base) *
               base_du[layout.base.housing]

    inherited = V14.energy_rate_ledger2D(
        base_u, p.base, op, time,
    )
    ambient = Float64(op.ambient_temperature(time))
    receiver_front = 0.0
    for group in 1:grid.group_count
        full_area = grid.multiplicity[group] *
                    grid.solid_area_per_channel
        core_area = full_area - skin_group_area[group]
        receiver_front += p.losses.front_loss_scale *
            p.losses.front_emissivity * V11.SIGMA * core_area *
            (
                V11.clamp_temperature(Tch[group, 1])^4 -
                ambient^4
            )
    end
    receiver_front += p.losses.front_loss_scale *
        p.losses.front_emissivity * V11.SIGMA * skin_area *
        (V11.clamp_temperature(Tskin[1])^4 - ambient^4)
    expected = inherited.expected_denergy_dt +
               inherited.receiver_front_loss - receiver_front
    return merge(inherited, (
        denergy_dt=denergy,
        expected_denergy_dt=expected,
        residual=denergy - expected,
        receiver_front_loss=receiver_front,
    ))
end

hydraulic_feedback_diagnostic2D(
    p::ModelParameters2D=default_parameters2D(); kwargs...,
) = V14.hydraulic_feedback_diagnostic2D(p.base; kwargs...)

end
