# ============================================================================
# 2D_v19.jl
#
# V18 geometry/source/assembly with:
#   * an orbit-integrated apparent Nu/UA law;
#   * a normalized finite-asymptote axial exchange kernel; and
#   * a gas/tube downstream T3 observation operator.
# ============================================================================

module Receiver2D_v19

include("2D_v18.jl")
const V18 = Receiver2D_v18
const V17 = V18.V17
const V16 = V18.V16
const V15 = V18.V15
const V14 = V18.V14
const V12 = V18.V12
const V11 = V18.V11

using DifferentialEquations
using OrdinaryDiffEq
using SciMLBase
using Statistics

export V18, V17, V16, V15, V14, V12, V11
export IntegratedExchangeParameters2D, OutletObservationParameters2D
export RearTubeFlangeParameters2D
export ModelParameters2D, SimulationResult2D
export default_parameters2D, base_grid2D, build_network_grid2D
export network_inventory2D, simulate2D, sensor_predictions2D
export integrated_ua_diagnostic2D, energy_rate_ledger2D
export exchange_energy_audit2D
export source_weights2D, source_power_invariant2D
export Geometry2D, SolidProperties2D, HeatTransferParameters2D
export LossParameters2D, OpticalParameters2D, HydraulicParameters2D
export OperatingCondition2D, AssemblyParameters2D, ObservationParameters2D
export NetworkParameters2D, SkinParameters2D
export CasingFlangeParameters2D, SourceParameters2D
export build_initial_state_2D, linear_history, get_t90_2D

const Geometry2D = V18.Geometry2D
const SolidProperties2D = V18.SolidProperties2D
const HeatTransferParameters2D = V18.HeatTransferParameters2D
const LossParameters2D = V18.LossParameters2D
const OpticalParameters2D = V18.OpticalParameters2D
const HydraulicParameters2D = V18.HydraulicParameters2D
const OperatingCondition2D = V18.OperatingCondition2D
const AssemblyParameters2D = V18.AssemblyParameters2D
const ObservationParameters2D = V18.ObservationParameters2D
const NetworkParameters2D = V18.NetworkParameters2D
const SkinParameters2D = V18.SkinParameters2D
const CasingFlangeParameters2D = V18.CasingFlangeParameters2D
const SourceParameters2D = V18.SourceParameters2D
const build_initial_state_2D = V18.build_initial_state_2D
const linear_history = V18.linear_history
const get_t90_2D = V18.get_t90_2D

Base.@kwdef struct IntegratedExchangeParameters2D
    # Average Nu = A Re^n.  The defaults are the measured manuscript fit.
    nu_prefactor::Float64 = 3.046900474047653e-4
    reynolds_exponent::Float64 = 1.444917063715293
    prandtl_exponent::Float64 = 0.0
    reference_prandtl::Float64 = 0.70
    kernel_model::Symbol = :graetz
    kernel_floor::Float64 = 1.0e-6
    kernel_iterations::Int = 2
end

Base.@kwdef struct OutletObservationParameters2D
    # Negligible-energy probe, per unit exposed sheath area.
    capacity_areal_J_m2_K::Float64 = 1000.0
    stem_conductance_areal_W_m2_K::Float64 = 20.0
    emissivity::Float64 = 0.30
    probe_diameter_m::Float64 = 1.5e-3
    axial_position_m::Float64 = 3.0e-3
end

Base.@kwdef struct RearTubeFlangeParameters2D
    contact_h_W_m2_K::Float64 = 0.0
    contact_start_m::Float64 = 28.0e-3
end

Base.@kwdef struct ModelParameters2D
    base::V18.ModelParameters2D = V18.default_parameters2D()
    integrated_exchange::IntegratedExchangeParameters2D =
        IntegratedExchangeParameters2D()
    outlet_observation::OutletObservationParameters2D =
        OutletObservationParameters2D()
    rear_tube_flange::RearTubeFlangeParameters2D =
        RearTubeFlangeParameters2D()
end

function Base.getproperty(p::ModelParameters2D, name::Symbol)
    name === :base && return getfield(p, :base)
    name === :integrated_exchange &&
        return getfield(p, :integrated_exchange)
    name === :outlet_observation &&
        return getfield(p, :outlet_observation)
    name === :rear_tube_flange &&
        return getfield(p, :rear_tube_flange)
    return getproperty(getfield(p, :base), name)
end

default_parameters2D() = ModelParameters2D()
base_grid2D(p::ModelParameters2D=default_parameters2D()) =
    V18.base_grid2D(p.base)
build_network_grid2D(p::ModelParameters2D=default_parameters2D()) =
    V18.build_network_grid2D(p.base)
source_weights2D(p::ModelParameters2D, grid) =
    V18.source_weights2D(p.base, grid)
source_power_invariant2D(
    p::ModelParameters2D=default_parameters2D(),
) = V18.source_power_invariant2D(p.base)

function network_inventory2D(
    p::ModelParameters2D=default_parameters2D(); kwargs...,
)
    return merge(V18.network_inventory2D(p.base; kwargs...), (
        integrated_nu_prefactor=
            p.integrated_exchange.nu_prefactor,
        integrated_reynolds_exponent=
            p.integrated_exchange.reynolds_exponent,
        integrated_kernel=p.integrated_exchange.kernel_model,
        outlet_probe_capacity_areal_J_m2_K=
            p.outlet_observation.capacity_areal_J_m2_K,
        outlet_probe_stem_conductance_areal_W_m2_K=
            p.outlet_observation.stem_conductance_areal_W_m2_K,
        rear_tube_flange_contact_h_W_m2_K=
            p.rear_tube_flange.contact_h_W_m2_K,
    ))
end

function _average_nusselt2D(
    mdot_channel, gas_reference, channel_temperature,
    p14, exchange::IntegratedExchangeParameters2D, grid,
)
    width = p14.geometry.channel_width
    area = width^2
    mu = V11.air_viscosity(gas_reference)
    reynolds = mdot_channel * width /
               max(area * mu, eps(Float64))
    cp = V11.air_heat_capacity(gas_reference)
    conductivity = V11.air_conductivity(gas_reference)
    prandtl = cp * mu / conductivity
    nu = max(0.0, exchange.nu_prefactor) *
         max(reynolds, eps(Float64))^
            exchange.reynolds_exponent *
         (
            max(prandtl, eps(Float64)) /
            max(exchange.reference_prandtl, eps(Float64))
         )^exchange.prandtl_exponent
    wall_average = sum(
        channel_temperature[j] * grid.base_grid.dz_cells[j]
        for j in 1:grid.base_grid.nz
    ) / p14.geometry.length
    film = 0.5 * (gas_reference + wall_average)
    kfilm = V11.air_conductivity(film)
    perimeter = 4.0 * width
    ua = nu * kfilm / width *
         perimeter * p14.geometry.length
    return (
        nusselt=nu, reynolds=reynolds, prandtl=prandtl,
        conductivity=kfilm, ua_W_K=ua,
    )
end

function _kernel_integrals2D!(
    kernel, gas, channel_temperature, mdot_channel,
    p14, exchange, grid,
)
    bg = grid.base_grid
    width = p14.geometry.channel_width
    area = width^2
    hp = p14.heat_transfer
    for j in 1:bg.nz
        Tgas = V11.clamp_temperature(gas[j])
        Twall = V11.clamp_temperature(channel_temperature[j])
        cp = V11.air_heat_capacity(Tgas)
        kg = V11.air_conductivity(Tgas)
        mu = V11.air_viscosity(Tgas)
        reynolds = mdot_channel * width /
                   max(area * mu, eps(Float64))
        prandtl = cp * mu / kg
        value = if exchange.kernel_model === :uniform
            bg.dz_cells[j]
        elseif exchange.kernel_model === :graetz
            midpoint = 0.5 * (
                bg.z_faces[j] + bg.z_faces[j + 1]
            )
            halfwidth = 0.5 * bg.dz_cells[j]
            nodes = (-sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0))
            weights = (5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0)
            halfwidth * sum(
                weights[q] * V11.graetz_nusselt2D(
                    reynolds, prandtl,
                    midpoint + halfwidth * nodes[q], width, hp;
                    gas_temperature=Tgas, wall_temperature=Twall,
                ).nusselt for q in 1:3
            )
        else
            error(
                "unknown v19 exchange kernel: " *
                "$(exchange.kernel_model)",
            )
        end
        kernel[j] = max(
            exchange.kernel_floor * bg.dz_cells[j], value,
        )
    end
    return kernel
end

function _march_channel_group2D!(
    work, group, mdot_channel, inlet, channel_temperature,
    p14, exchange::IntegratedExchangeParameters2D, grid,
    ua_per_channel,
)
    bg = grid.base_grid
    hyd = p14.hydraulics
    width = p14.geometry.channel_width
    area = width^2
    perimeter = 4.0 * width
    multiplicity = grid.multiplicity[group]
    work.front_h[group] = 0.0
    work.front_qgas[group] = 0.0
    work.gas[group, :] .= inlet
    kernel = zeros(bg.nz)
    ua_cell = zeros(bg.nz)
    iterations = max(1, exchange.kernel_iterations)
    for _ in 1:iterations
        _kernel_integrals2D!(
            kernel, view(work.gas, group, :),
            view(channel_temperature, group, :),
            mdot_channel, p14, exchange, grid,
        )
        normalization = sum(kernel)
        ua_cell .= ua_per_channel .* kernel ./
                   max(normalization, eps(Float64))
        work.gas[group, 1] = inlet
        for j in 1:bg.nz
            Tgas = V11.clamp_temperature(work.gas[group, j])
            Twall = V11.clamp_temperature(
                channel_temperature[group, j],
            )
            cp = V11.air_heat_capacity(Tgas)
            effectiveness = -expm1(
                -ua_cell[j] /
                max(mdot_channel * cp, eps(Float64)),
            )
            work.gas[group, j + 1] = Tgas +
                effectiveness * (Twall - Tgas)
        end
    end

    for j in 1:bg.nz
        Tgas = V11.clamp_temperature(work.gas[group, j])
        cp = V11.air_heat_capacity(Tgas)
        kg = V11.air_conductivity(Tgas)
        mu = V11.air_viscosity(Tgas)
        work.qgas[group, j] =
            multiplicity * mdot_channel * cp *
            (work.gas[group, j + 1] - Tgas)
        work.h[group, j] = ua_cell[j] /
            max(perimeter * bg.dz_cells[j], eps(Float64))
        work.nusselt[group, j] =
            work.h[group, j] * width / max(kg, eps(Float64))
        work.prandtl[group, j] = cp * mu / kg

        Tbulk = V11.clamp_temperature(0.5 * (
            work.gas[group, j] + work.gas[group, j + 1]
        ))
        rho = V11.air_density(
            Tbulk, hyd.atmospheric_pressure,
        )
        mu_bulk = V11.air_viscosity(Tbulk)
        velocity = mdot_channel / max(rho * area, eps(Float64))
        work.density[group, j] = rho
        work.velocity[group, j] = velocity
        work.reynolds[group, j] =
            rho * velocity * width / mu_bulk
        work.graetz[group, j] =
            work.prandtl[group, j] *
            work.reynolds[group, j] * width /
            max(bg.z_centers[j], eps(Float64))
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
    work.channel_dp[group] = sum(view(work.dp_cell, group, :))
    return nothing
end

function _mean_wall_profile2D(channel_temperature, grid)
    denominator = sum(grid.multiplicity)
    return [
        sum(
            grid.multiplicity[group] *
            channel_temperature[group, j]
            for group in 1:grid.group_count
        ) / denominator
        for j in 1:grid.base_grid.nz
    ]
end

function _mixed_outlet2D(work, grid)
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
    return total_mcp > 0.0 ?
        total_enthalpy / total_mcp : NaN
end

function _solve_channel_flows2D!(
    work, mdot_total, inlet, channel_temperature,
    p14, exchange, grid,
)
    ng = grid.group_count
    hyd = p14.hydraulics
    work.channel_mdot .= mdot_total / sum(grid.multiplicity)
    relaxation = clamp(
        hyd.equal_pressure_relaxation, 0.05, 1.0,
    )
    tolerance = max(hyd.equal_pressure_tolerance, 1e-10)
    max_iterations = max(hyd.equal_pressure_max_iterations, 24)
    gas_reference = inlet
    wall_profile = _mean_wall_profile2D(
        channel_temperature, grid,
    )
    for iteration in 1:max_iterations
        average = _average_nusselt2D(
            mdot_total / sum(grid.multiplicity),
            gas_reference, wall_profile, p14, exchange, grid,
        )
        for group in 1:ng
            _march_channel_group2D!(
                work, group, work.channel_mdot[group], inlet,
                channel_temperature, p14, exchange, grid,
                average.ua_W_K,
            )
        end
        mixed_outlet = _mixed_outlet2D(work, grid)
        new_gas_reference = 0.5 * (inlet + mixed_outlet)
        gas_error = abs(new_gas_reference - gas_reference) /
            max(abs(gas_reference), 1.0)
        common = sum(
            grid.multiplicity[group] *
            work.channel_mdot[group] * work.channel_dp[group]
            for group in 1:ng
        ) / mdot_total
        error = maximum(
            abs(work.channel_dp[group] - common)
            for group in 1:ng
        ) / max(common, eps(Float64))
        work.common_dp = common
        work.equal_pressure_error = error
        if error <= tolerance && gas_error <= tolerance
            return (
                gas_reference_relative_error=gas_error,
                iterations=iteration,
                converged=true,
            )
        end
        linear_coeff = zeros(ng)
        quadratic_coeff = zeros(ng)
        for group in 1:ng
            mdot = max(work.channel_mdot[group], eps(Float64))
            quadratic_coeff[group] =
                work.groove_dp[group] / mdot^2
            linear_coeff[group] = max(
                0.0,
                work.channel_dp[group] - work.groove_dp[group],
            ) / mdot
        end
        candidate = V14._equal_pressure_candidate2D(
            linear_coeff, quadratic_coeff, grid.multiplicity,
            mdot_total, common,
        )
        work.channel_mdot .=
            (1.0 - relaxation) .* work.channel_mdot .+
            relaxation .* candidate
        work.channel_mdot .*= mdot_total /
            sum(grid.multiplicity .* work.channel_mdot)
        gas_reference = new_gas_reference
    end
    average = _average_nusselt2D(
        mdot_total / sum(grid.multiplicity),
        gas_reference, wall_profile, p14, exchange, grid,
    )
    for group in 1:ng
        _march_channel_group2D!(
            work, group, work.channel_mdot[group], inlet,
            channel_temperature, p14, exchange, grid,
            average.ua_W_K,
        )
    end
    work.common_dp = sum(
        grid.multiplicity[group] *
        work.channel_mdot[group] * work.channel_dp[group]
        for group in 1:ng
    ) / mdot_total
    work.equal_pressure_error = maximum(
        abs(work.channel_dp[group] - work.common_dp)
        for group in 1:ng
    ) / max(work.common_dp, eps(Float64))
    mixed_outlet = _mixed_outlet2D(work, grid)
    new_gas_reference = 0.5 * (inlet + mixed_outlet)
    gas_error = abs(new_gas_reference - gas_reference) /
        max(abs(gas_reference), 1.0)
    return (
        gas_reference_relative_error=gas_error,
        iterations=max_iterations,
        converged=(
            work.equal_pressure_error <= tolerance &&
            gas_error <= tolerance
        ),
    )
end

function _gas_profile_integrated2D!(
    work, channel_temperature, tube_temperature, time,
    p14, exchange, op, grid,
)
    bg = grid.base_grid
    flow_lpm = max(0.0, Float64(op.flow_lpm(time)))
    inlet = Float64(op.inlet_temperature(time))
    mdot_total = V11.standard_mass_flow2D(
        flow_lpm, p14.hydraulics,
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
            inlet, p14.hydraulics.atmospheric_pressure,
        )
        work.gas_rear .= inlet
        return (
            gas_reference_relative_error=0.0,
            iterations=0,
            converged=true,
        )
    end
    convergence = _solve_channel_flows2D!(
        work, mdot_total, inlet, channel_temperature,
        p14, exchange, grid,
    )
    work.gas_rear[1] = _mixed_outlet2D(work, grid)

    radius = p14.geometry.rear_tube_inner_radius
    dh = 2.0 * radius
    area = pi * radius^2
    perimeter = 2.0 * pi * radius
    for j in 1:bg.nz_rear
        Tfilm = V11.clamp_temperature(0.5 * (
            tube_temperature[j] + work.gas_rear[j]
        ))
        cp = V11.air_heat_capacity(Tfilm)
        kf = V11.air_conductivity(Tfilm)
        mu = V11.air_viscosity(Tfilm)
        reynolds = mdot_total * dh /
                   max(area * mu, eps(Float64))
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
        work.qgas_rear[j] = mdot_total * cp * (
            work.gas_rear[j + 1] - work.gas_rear[j]
        )
    end
    return convergence
end

function _gas_power_snapshot2D(work)
    return (
        qgas=copy(work.qgas),
        front=copy(work.front_qgas),
        rear=copy(work.qgas_rear),
    )
end

function _v19_ode2D!(du, u, context, time)
    V17._casing_flange_ode2D!(
        du, u, context.casing_context, time,
    )
    old = _gas_power_snapshot2D(context.work)
    grid = context.grid
    bg = grid.base_grid
    layout = context.layout
    base_layout = layout.base
    base_u = view(u, layout.base_state)
    base_du = view(du, layout.base_state)
    Tch = reshape(
        view(base_u, base_layout.channel),
        grid.group_count, bg.nz,
    )
    Ttube = view(base_u, base_layout.tube)
    _gas_profile_integrated2D!(
        context.work, Tch, Ttube, time, context.p14,
        context.exchange, context.op, grid,
    )
    duch = reshape(
        view(base_du, base_layout.channel),
        grid.group_count, bg.nz,
    )
    skin_areas = context.skin_group_area
    for group in 1:grid.group_count, j in 1:bg.nz
        full_area = grid.multiplicity[group] *
                    grid.solid_area_per_channel
        core_area = full_area - skin_areas[group]
        capacity = V15._channel_capacity2D(
            Tch[group, j], core_area, bg.dz_cells[j],
            context.p15,
        )
        delta = context.work.qgas[group, j] -
                old.qgas[group, j]
        j == 1 && (delta +=
            context.work.front_qgas[group] - old.front[group])
        duch[group, j] -= delta /
            max(capacity, eps(Float64))
    end
    tube_area = pi * (
        context.p14.geometry.rear_tube_outer_radius^2 -
        context.p14.geometry.rear_tube_inner_radius^2
    )
    dutube = view(base_du, base_layout.tube)
    for j in 1:bg.nz_rear
        capacity = 3900.0 *
            V11.alumina_heat_capacity(Ttube[j]) *
            tube_area * bg.dz_rear
        dutube[j] -= (
            context.work.qgas_rear[j] - old.rear[j]
        ) / max(capacity, eps(Float64))
        overlap = V12._overlap_length(
            bg.z_rear_faces[j], bg.z_rear_faces[j + 1],
            context.tube_flange.contact_start_m,
            context.p14.geometry.rear_tube_length,
        )
        if overlap > 0.0 &&
           context.tube_flange.contact_h_W_m2_K > 0.0
            area = 2.0 * pi *
                context.p14.geometry.rear_tube_outer_radius *
                overlap
            Q = context.tube_flange.contact_h_W_m2_K *
                area * (Ttube[j] - V11.WATER_FLANGE_TEMP)
            dutube[j] -= Q / max(capacity, eps(Float64))
        end
    end
    return nothing
end

function _model_context2D(p::ModelParameters2D, op)
    grid = build_network_grid2D(p)
    layout = V15._state_layout(grid)
    work = V14.NetworkWorkspace2D(grid)
    weights = source_weights2D(p, grid)
    p17 = p.base.base
    p15 = p17.base
    p14 = p15.base
    base_context = (
        p=p14, op=op, grid=grid, layout=layout.base,
        work=work, weights=weights,
    )
    skin_context = (
        p=p15, op=op, grid=grid, layout=layout,
        base_context=base_context, weights=weights,
        skin_group_area=V15._skin_group_areas2D(p15, grid),
    )
    casing_context = (
        p=p17, op=op, grid=grid, layout=layout,
        skin_context=skin_context,
    )
    return (
        p=p, p14=p14, p15=p15, op=op, grid=grid,
        layout=layout, work=work, weights=weights,
        skin_group_area=skin_context.skin_group_area,
        casing_context=casing_context,
        exchange=p.integrated_exchange,
        tube_flange=p.rear_tube_flange,
    )
end

"""
Independently compare the capacity-weighted derivative introduced by v19
with the gas-exchange and distributed tube/flange powers that should produce
it. This is not derived from the inherited energy-ledger residual and works
on C, M, and F meshes.
"""
function exchange_energy_audit2D(
    u::AbstractVector, p::ModelParameters2D,
    op::OperatingCondition2D, time::Real=0.0,
)
    baseline_context = _model_context2D(p, op)
    v19_context = _model_context2D(p, op)
    du_baseline = zeros(length(u))
    du_v19 = zeros(length(u))
    V17._casing_flange_ode2D!(
        du_baseline, u, baseline_context.casing_context, time,
    )
    _v19_ode2D!(du_v19, u, v19_context, time)

    grid = v19_context.grid
    bg = grid.base_grid
    layout = v19_context.layout
    base_layout = layout.base
    base_u = view(u, layout.base_state)
    Tch = reshape(
        view(base_u, base_layout.channel),
        grid.group_count, bg.nz,
    )
    Ttube = view(base_u, base_layout.tube)
    delta_du = view(du_v19, layout.base_state) .-
               view(du_baseline, layout.base_state)
    delta_channel = reshape(
        view(delta_du, base_layout.channel),
        grid.group_count, bg.nz,
    )
    derivative_power = 0.0
    for group in 1:grid.group_count, j in 1:bg.nz
        full_area = grid.multiplicity[group] *
                    grid.solid_area_per_channel
        core_area = full_area -
                    v19_context.skin_group_area[group]
        capacity = V15._channel_capacity2D(
            Tch[group, j], core_area, bg.dz_cells[j],
            v19_context.p15,
        )
        derivative_power +=
            capacity * delta_channel[group, j]
    end
    tube_area = pi * (
        v19_context.p14.geometry.rear_tube_outer_radius^2 -
        v19_context.p14.geometry.rear_tube_inner_radius^2
    )
    delta_tube = view(delta_du, base_layout.tube)
    for j in 1:bg.nz_rear
        capacity = 3900.0 *
            V11.alumina_heat_capacity(Ttube[j]) *
            tube_area * bg.dz_rear
        derivative_power += capacity * delta_tube[j]
    end

    old = baseline_context.work
    new = v19_context.work
    gas_correction = (
        sum(new.qgas) + sum(new.front_qgas) +
        sum(new.qgas_rear)
    ) - (
        sum(old.qgas) + sum(old.front_qgas) +
        sum(old.qgas_rear)
    )
    tube_flange_loss = 0.0
    for j in 1:bg.nz_rear
        overlap = V12._overlap_length(
            bg.z_rear_faces[j], bg.z_rear_faces[j + 1],
            p.rear_tube_flange.contact_start_m,
            p.geometry.rear_tube_length,
        )
        overlap <= 0.0 && continue
        area = 2.0 * pi *
            p.geometry.rear_tube_outer_radius * overlap
        tube_flange_loss +=
            p.rear_tube_flange.contact_h_W_m2_K *
            area * (Ttube[j] - V11.WATER_FLANGE_TEMP)
    end
    expected_power = -gas_correction - tube_flange_loss
    residual = derivative_power - expected_power
    return (
        capacity_derivative_correction_W=derivative_power,
        expected_correction_W=expected_power,
        gas_correction_W=gas_correction,
        distributed_tube_flange_loss_W=tube_flange_loss,
        residual_W=residual,
        relative_residual=abs(residual) /
            max(abs(expected_power), 1.0),
    )
end

struct SimulationResult2D
    base_result::V17.SimulationResult2D
    parameters::ModelParameters2D
    diagnostics::NamedTuple
    ode_solution::Any
end

function Base.getproperty(result::SimulationResult2D, name::Symbol)
    name === :base_result && return getfield(result, :base_result)
    name === :parameters && return getfield(result, :parameters)
    name === :diagnostics && return getfield(result, :diagnostics)
    name === :ode_solution && return getfield(result, :ode_solution)
    diagnostics = getfield(result, :diagnostics)
    name in propertynames(diagnostics) &&
        return getproperty(diagnostics, name)
    return getproperty(getfield(result, :base_result), name)
end

function _extract_diagnostics2D(
    solution, p, op, grid, layout, p14,
)
    bg = grid.base_grid
    nt = length(solution.t)
    ng = grid.group_count
    work = V14.NetworkWorkspace2D(grid)
    gas = zeros(ng, bg.nz + 1, nt)
    rear_gas = zeros(bg.nz_rear + 1, nt)
    h = zeros(ng, bg.nz, nt)
    nusselt = zeros(ng, bg.nz, nt)
    density = zeros(ng, bg.nz, nt)
    velocity = zeros(ng, bg.nz, nt)
    reynolds = zeros(ng, bg.nz, nt)
    prandtl = zeros(ng, bg.nz, nt)
    graetz = zeros(ng, bg.nz, nt)
    dp_cell = zeros(ng, bg.nz, nt)
    channel_mdot = zeros(ng, nt)
    group_mdot = zeros(ng, nt)
    channel_dp = zeros(ng, nt)
    groove_dp = zeros(ng, nt)
    equal_error = zeros(nt)
    gas_reference_error = zeros(nt)
    flow_iterations = zeros(Int, nt)
    flow_converged = falses(nt)
    receiver_dp = zeros(nt)
    dp1 = zeros(nt)
    mass_flow = zeros(nt)
    for index in 1:nt
        base_u = view(solution.u[index], layout.base_state)
        Tch = reshape(
            view(base_u, layout.base.channel), ng, bg.nz,
        )
        Ttube = view(base_u, layout.base.tube)
        convergence = _gas_profile_integrated2D!(
            work, Tch, Ttube, solution.t[index], p14,
            p.integrated_exchange, op, grid,
        )
        gas[:, :, index] .= work.gas
        rear_gas[:, index] .= work.gas_rear
        h[:, :, index] .= work.h
        nusselt[:, :, index] .= work.nusselt
        density[:, :, index] .= work.density
        velocity[:, :, index] .= work.velocity
        reynolds[:, :, index] .= work.reynolds
        prandtl[:, :, index] .= work.prandtl
        graetz[:, :, index] .= work.graetz
        dp_cell[:, :, index] .= work.dp_cell
        channel_mdot[:, index] .= work.channel_mdot
        group_mdot[:, index] .=
            grid.multiplicity .* work.channel_mdot
        channel_dp[:, index] .= work.channel_dp
        groove_dp[:, index] .= work.groove_dp
        equal_error[index] = work.equal_pressure_error
        gas_reference_error[index] =
            convergence.gas_reference_relative_error
        flow_iterations[index] = convergence.iterations
        flow_converged[index] = convergence.converged
        receiver_dp[index] = work.common_dp / 100.0
        dp1[index] = p14.hydraulics.dp1_zero_offset_mbar +
                     receiver_dp[index]
        mass_flow[index] = V11.standard_mass_flow2D(
            max(0.0, Float64(op.flow_lpm(solution.t[index]))),
            p14.hydraulics,
        )
    end
    return (
        gas_temperature=gas,
        rear_gas_temperature=rear_gas,
        heat_transfer_coefficient=h,
        gas_nusselt=nusselt,
        gas_density=density,
        gas_velocity=velocity,
        gas_reynolds=reynolds,
        gas_prandtl=prandtl,
        gas_graetz=graetz,
        cell_pressure_drop=dp_cell,
        channel_mass_flow_kg_s=channel_mdot,
        group_mass_flow_kg_s=group_mdot,
        channel_pressure_drop_Pa=channel_dp,
        groove_pressure_drop_Pa=groove_dp,
        equal_pressure_relative_error=equal_error,
        gas_reference_relative_error=gas_reference_error,
        flow_solver_iterations=flow_iterations,
        flow_solver_converged=flow_converged,
        receiver_pressure_drop_mbar=receiver_dp,
        dp1_prediction_mbar=dp1,
        mass_flow_kg_s=mass_flow,
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
    context = _model_context2D(p, op)
    grid = context.grid
    layout = context.layout
    p17 = p.base.base
    p15 = context.p15
    p14 = context.p14
    initial = V15._initial_state2D(
        initial_temperature, p15, grid,
    )
    problem = ODEProblem(
        _v19_ode2D!, initial,
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
    base14_result = V15._extract_base_result2D(
        solution, p15, op, grid, layout,
    )
    skin = hcat((
        collect(view(state, layout.skin))
        for state in solution.u
    )...)
    base15_result = V15.SimulationResult2D(
        base14_result, skin, p15, solution,
    )
    base17_result = V17.SimulationResult2D(
        base15_result, p17, solution,
    )
    diagnostics = _extract_diagnostics2D(
        solution, p, op, grid, layout, p14,
    )
    return SimulationResult2D(
        base17_result, p, diagnostics, solution,
    )
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

function _sample_rear_gas2D(result, z)
    zfaces = collect(range(
        0.0, result.parameters.geometry.rear_tube_length,
        length=size(result.rear_gas_temperature, 1),
    ))
    j0, j1, weight = V11._axis_interp_indices(zfaces, z)
    return vec(
        (1.0 - weight) .* result.rear_gas_temperature[j0, :] .+
        weight .* result.rear_gas_temperature[j1, :],
    )
end

function _sample_rear_tube2D(result, z)
    j0, j1, weight = V11._axis_interp_indices(
        result.z_rear, z,
    )
    return vec(
        (1.0 - weight) .* result.rear_tube_temperature[j0, :] .+
        weight .* result.rear_tube_temperature[j1, :],
    )
end

function _probe_gas_coefficient2D(
    temperature, mdot, p::ModelParameters2D,
)
    diameter = p.outlet_observation.probe_diameter_m
    diameter > 0.0 || error(
        "v19 probe diameter must be positive",
    )
    tube_area = pi * p.geometry.rear_tube_inner_radius^2
    rho = V11.air_density(
        temperature, p.hydraulics.atmospheric_pressure,
    )
    velocity = mdot / max(rho * tube_area, eps(Float64))
    mu = V11.air_viscosity(temperature)
    conductivity = V11.air_conductivity(temperature)
    cp = V11.air_heat_capacity(temperature)
    reynolds = rho * velocity * diameter / mu
    prandtl = cp * mu / conductivity
    numerator = 0.62 * sqrt(max(reynolds, 0.0)) *
        prandtl^(1.0 / 3.0)
    denominator = (
        1.0 + (0.4 / max(prandtl, eps(Float64)))^(2.0 / 3.0)
    )^(1.0 / 4.0)
    nusselt = 0.3 + numerator / denominator *
        (
            1.0 + (reynolds / 282000.0)^(5.0 / 8.0)
        )^(4.0 / 5.0)
    return nusselt * conductivity / diameter
end

function _probe_history2D(
    result, Tgas, Ttube; initial_temperature=nothing,
)
    parameters = result.parameters.outlet_observation
    capacity = parameters.capacity_areal_J_m2_K
    stem = parameters.stem_conductance_areal_W_m2_K
    emissivity = parameters.emissivity
    capacity > 0.0 || error(
        "v19 probe areal capacity must be positive",
    )
    stem >= 0.0 || error(
        "v19 probe stem conductance must be nonnegative",
    )
    0.0 <= emissivity <= 1.0 || error(
        "v19 probe emissivity must lie in [0,1]",
    )
    time = result.time
    output = similar(Tgas)
    output[1] = initial_temperature === nothing ?
        Tgas[1] : Float64(initial_temperature)
    sink = V11.WATER_FLANGE_TEMP
    for index in 2:length(time)
        dt = time[index] - time[index - 1]
        dt >= 0.0 || error("v19 probe time must be monotone")
        Tprevious = output[index - 1]
        Tgas_mid = 0.5 * (Tgas[index - 1] + Tgas[index])
        Ttube_mid = 0.5 * (Ttube[index - 1] + Ttube[index])
        mdot_mid = 0.5 * (
            result.mass_flow_kg_s[index - 1] +
            result.mass_flow_kg_s[index]
        )
        hgas = _probe_gas_coefficient2D(
            Tgas_mid, mdot_mid, result.parameters,
        )
        hrad = emissivity * V11.SIGMA *
            (Ttube_mid^2 + Tprevious^2) *
            (Ttube_mid + Tprevious)
        total = hgas + hrad + stem
        target = (
            hgas * Tgas_mid + hrad * Ttube_mid + stem * sink
        ) / max(total, eps(Float64))
        output[index] = target + (Tprevious - target) *
            exp(-total * dt / capacity)
    end
    return output
end

function sensor_predictions2D(
    result::SimulationResult2D; initial_T3=nothing,
)
    obs = result.parameters.observation
    outlet = result.parameters.outlet_observation
    time = result.time
    side11 = _sample_skin2D(
        result, V11.SIDE_TC_FRONT_Z_2D,
    )
    side58 = _sample_skin2D(result, 58e-3)
    side107 = _sample_skin2D(result, 107e-3)
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
    z3 = outlet.axial_position_m
    T3gas = _sample_rear_gas2D(result, z3)
    T3tube = _sample_rear_tube2D(result, z3)
    T3probe = _probe_history2D(
        result, T3gas, T3tube; initial_temperature=initial_T3,
    )
    rT2 = result.parameters.geometry.receiver_radius + 40e-3
    T2 = V14._sample_outer2D(result.base_result, rT2, 58e-3)
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
        T3=T3probe,
        T2=T2,
        T8_skin=side11, T12_skin=side58, T11_skin=side107,
        T9_wall=wall58, T10_wall=wall107,
        T9_gas=gas58, T10_gas=gas107,
        T3_gas_raw=T3gas, T3_tube=T3tube,
        T3_probe=T3probe,
    )
end

function integrated_ua_diagnostic2D(
    p::ModelParameters2D=default_parameters2D();
    flow_lpm=10.0, inlet_temperature=300.0,
    wall_temperature=800.0,
)
    grid = build_network_grid2D(p)
    p14 = p.base.base.base.base
    mdot_total = V11.standard_mass_flow2D(
        flow_lpm, p14.hydraulics,
    )
    mdot_channel = mdot_total / sum(grid.multiplicity)
    wall = fill(wall_temperature, grid.base_grid.nz)
    average = _average_nusselt2D(
        mdot_channel, inlet_temperature, wall,
        p14, p.integrated_exchange, grid,
    )
    work = V14.NetworkWorkspace2D(grid)
    Tch = fill(
        wall_temperature,
        grid.group_count, grid.base_grid.nz,
    )
    _march_channel_group2D!(
        work, 1, mdot_channel, inlet_temperature, Tch,
        p14, p.integrated_exchange, grid, average.ua_W_K,
    )
    perimeter = 4.0 * p14.geometry.channel_width
    recovered_ua = sum(
        work.h[1, j] * perimeter *
        grid.base_grid.dz_cells[j]
        for j in 1:grid.base_grid.nz
    )
    return merge(average, (
        recovered_ua_W_K=recovered_ua,
        relative_ua_error=abs(recovered_ua - average.ua_W_K) /
            max(average.ua_W_K, eps(Float64)),
        outlet_temperature_K=work.gas[1, end],
    ))
end

function energy_rate_ledger2D(
    u::AbstractVector, p::ModelParameters2D,
    op::OperatingCondition2D, time::Real=0.0,
)
    p.geometry.nodes_z == 48 || error(
        "v19 detailed energy ledger is exact only on M",
    )
    inherited = V18.energy_rate_ledger2D(
        u, p.base, op, time,
    )
    grid = build_network_grid2D(p)
    layout = V15._state_layout(grid)
    p14 = p.base.base.base.base
    base_u = view(u, layout.base_state)
    Tch = reshape(
        view(base_u, layout.base.channel),
        grid.group_count, grid.base_grid.nz,
    )
    Ttube = view(base_u, layout.base.tube)
    old_work = V14.NetworkWorkspace2D(grid)
    V14._gas_profile_network2D!(
        old_work, Tch, Ttube, time, p14, op, grid,
    )
    new_work = V14.NetworkWorkspace2D(grid)
    _gas_profile_integrated2D!(
        new_work, Tch, Ttube, time, p14,
        p.integrated_exchange, op, grid,
    )
    old_gas = sum(old_work.qgas) +
              sum(old_work.front_qgas) +
              sum(old_work.qgas_rear)
    new_gas = sum(new_work.qgas) +
              sum(new_work.front_qgas) +
              sum(new_work.qgas_rear)
    correction = new_gas - old_gas
    tube_flange_loss = 0.0
    if p.rear_tube_flange.contact_h_W_m2_K > 0.0
        bg = grid.base_grid
        for j in 1:bg.nz_rear
            overlap = V12._overlap_length(
                bg.z_rear_faces[j], bg.z_rear_faces[j + 1],
                p.rear_tube_flange.contact_start_m,
                p.geometry.rear_tube_length,
            )
            overlap <= 0.0 && continue
            area = 2.0 * pi *
                p.geometry.rear_tube_outer_radius * overlap
            tube_flange_loss +=
                p.rear_tube_flange.contact_h_W_m2_K *
                area * (
                    Ttube[j] - V11.WATER_FLANGE_TEMP
                )
        end
    end
    return merge(inherited, (
        denergy_dt=inherited.denergy_dt - correction -
            tube_flange_loss,
        expected_denergy_dt=
            inherited.expected_denergy_dt - correction -
            tube_flange_loss,
        residual=inherited.residual,
        inherited_gas_removal_W=old_gas,
        integrated_gas_removal_W=new_gas,
        gas_removal_correction_W=correction,
        distributed_tube_flange_loss_W=tube_flange_loss,
    ))
end

end
