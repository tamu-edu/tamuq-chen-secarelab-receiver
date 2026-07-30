# ============================================================================
# 2D_v20.jl
#
# Post-v19 identifiability model:
#   * exact variable-cp enthalpy transport and mixing;
#   * unchanged v19 solid/source/hydraulic architecture and integrated UA law;
#   * discrete, auditable T3 location hypotheses.
# ============================================================================

module Receiver2D_v20

include("2D_v19.jl")
const V19 = Receiver2D_v19
const V18 = V19.V18
const V17 = V19.V17
const V16 = V19.V16
const V15 = V19.V15
const V14 = V19.V14
const V12 = V19.V12
const V11 = V19.V11

using DifferentialEquations
using OrdinaryDiffEq
using SciMLBase

export V19, V18, V17, V16, V15, V14, V12, V11
export T3LocationParameters2D, ModelParameters2D, SimulationResult2D
export default_parameters2D, base_grid2D, build_network_grid2D
export network_inventory2D, simulate2D, sensor_predictions2D
export air_specific_enthalpy2D, air_temperature_from_enthalpy2D
export enthalpy_transport_audit2D, exchange_energy_audit2D
export source_weights2D, source_power_invariant2D
export Geometry2D, SolidProperties2D, HeatTransferParameters2D
export LossParameters2D, OpticalParameters2D, HydraulicParameters2D
export OperatingCondition2D, AssemblyParameters2D, ObservationParameters2D
export NetworkParameters2D, SkinParameters2D
export CasingFlangeParameters2D, SourceParameters2D
export IntegratedExchangeParameters2D, OutletObservationParameters2D
export RearTubeFlangeParameters2D
export build_initial_state_2D, linear_history, get_t90_2D

const Geometry2D = V19.Geometry2D
const SolidProperties2D = V19.SolidProperties2D
const HeatTransferParameters2D = V19.HeatTransferParameters2D
const LossParameters2D = V19.LossParameters2D
const OpticalParameters2D = V19.OpticalParameters2D
const HydraulicParameters2D = V19.HydraulicParameters2D
const OperatingCondition2D = V19.OperatingCondition2D
const AssemblyParameters2D = V19.AssemblyParameters2D
const ObservationParameters2D = V19.ObservationParameters2D
const NetworkParameters2D = V19.NetworkParameters2D
const SkinParameters2D = V19.SkinParameters2D
const CasingFlangeParameters2D = V19.CasingFlangeParameters2D
const SourceParameters2D = V19.SourceParameters2D
const IntegratedExchangeParameters2D = V19.IntegratedExchangeParameters2D
const OutletObservationParameters2D = V19.OutletObservationParameters2D
const RearTubeFlangeParameters2D = V19.RearTubeFlangeParameters2D
const build_initial_state_2D = V19.build_initial_state_2D
const linear_history = V19.linear_history
const get_t90_2D = V19.get_t90_2D

"""
Discrete T3 interpretations. `location` must be one of:

* `:receiver_136` - mass-mixed channel gas at global z=136 mm;
* `:receiver_exit` - mass-mixed receiver outlet at global z=137 mm;
* `:rear_003` - rear alumina-tube gas 3 mm downstream (global z=140 mm).

The location is deliberately discrete because the installation record does
not justify fitting a continuous axial coordinate.
"""
Base.@kwdef struct T3LocationParameters2D
    location::Symbol = :receiver_136
end

Base.@kwdef struct ModelParameters2D
    base::V19.ModelParameters2D = V19.default_parameters2D()
    t3_location::T3LocationParameters2D = T3LocationParameters2D()
end

function Base.getproperty(p::ModelParameters2D, name::Symbol)
    name === :base && return getfield(p, :base)
    name === :t3_location && return getfield(p, :t3_location)
    return getproperty(getfield(p, :base), name)
end

default_parameters2D() = ModelParameters2D()
base_grid2D(p::ModelParameters2D=default_parameters2D()) =
    V19.base_grid2D(p.base)
build_network_grid2D(p::ModelParameters2D=default_parameters2D()) =
    V19.build_network_grid2D(p.base)
source_weights2D(p::ModelParameters2D, grid) =
    V19.source_weights2D(p.base, grid)
source_power_invariant2D(
    p::ModelParameters2D=default_parameters2D(),
) = V19.source_power_invariant2D(p.base)

function network_inventory2D(
    p::ModelParameters2D=default_parameters2D(); kwargs...,
)
    return merge(V19.network_inventory2D(p.base; kwargs...), (
        gas_thermodynamics=:integral_enthalpy,
        t3_location=p.t3_location.location,
    ))
end

# V11 uses cp(T)=1004*(1+aT-bT^2), with T clamped below 200 K.
const _AIR_CP_SCALE = 1004.0
const _AIR_CP_LINEAR = 1.983e-4
const _AIR_CP_QUADRATIC = -4.14e-8
const _AIR_TMIN = 200.0

function _air_enthalpy_polynomial(T)
    return _AIR_CP_SCALE * (
        T + 0.5 * _AIR_CP_LINEAR * T^2 +
        (_AIR_CP_QUADRATIC / 3.0) * T^3
    )
end

"""
Specific air enthalpy relative to 200 K [J/kg], analytically consistent with
`V11.air_heat_capacity`, including its constant-cp continuation below 200 K.
"""
function air_specific_enthalpy2D(T)
    value = Float64(T)
    if value < _AIR_TMIN
        return V11.air_heat_capacity(_AIR_TMIN) *
               (value - _AIR_TMIN)
    end
    return _air_enthalpy_polynomial(value) -
           _air_enthalpy_polynomial(_AIR_TMIN)
end

"""
Invert the monotone air enthalpy relation with safeguarded Newton iteration.
"""
function air_temperature_from_enthalpy2D(
    target; lower=100.0, upper=3500.0,
)
    lo = Float64(lower)
    hi = Float64(upper)
    lo <= hi || ((lo, hi) = (hi, lo))
    hlo = air_specific_enthalpy2D(lo)
    hhi = air_specific_enthalpy2D(hi)
    value = Float64(target)
    value <= hlo && return lo
    value >= hhi && return hi
    T = lo + (hi - lo) * (value - hlo) / max(hhi - hlo, eps())
    for _ in 1:12
        residual = air_specific_enthalpy2D(T) - value
        abs(residual) <= 1e-10 * max(abs(value), 1.0) && return T
        residual > 0.0 ? (hi = T) : (lo = T)
        candidate = T - residual / V11.air_heat_capacity(T)
        T = lo < candidate < hi ? candidate : 0.5 * (lo + hi)
    end
    return T
end

# Antiderivative of cp(T)/(Twall-T), valid in the polynomial range.
function _air_exchange_primitive2D(T, Twall)
    c0 = _AIR_CP_SCALE
    c1 = _AIR_CP_SCALE * _AIR_CP_LINEAR
    c2 = _AIR_CP_SCALE * _AIR_CP_QUADRATIC
    cpwall = c0 + c1 * Twall + c2 * Twall^2
    distance = max(abs(Twall - T), 1e-14)
    return -0.5 * c2 * T^2 -
           (c1 + c2 * Twall) * T -
           cpwall * log(distance)
end

"""
Exact constant-wall, variable-cp plug-flow update:

    UA/mdot = integral(Tin,Tout) cp(T)/(Twall-T) dT

and q=mdot*(h(Tout)-h(Tin)). The returned temperature is bounded between the
inlet gas and wall temperatures.
"""
function _enthalpy_cell_exchange2D(Tin, Twall, ua, mdot)
    gas = V11.clamp_temperature(Tin)
    wall = V11.clamp_temperature(Twall)
    conductance = max(0.0, Float64(ua))
    flow = max(0.0, Float64(mdot))
    if flow <= eps() || conductance <= eps() ||
       abs(wall - gas) <= 1e-12
        return (temperature=gas, heat_W=0.0)
    end

    target = conductance / flow
    cp0 = V11.air_heat_capacity(0.5 * (gas + wall))
    fraction = clamp(-expm1(-target / cp0), 0.0, 1.0 - 1e-12)
    lower_fraction = 0.0
    upper_fraction = 1.0 - 1e-12
    primitive_in = _air_exchange_primitive2D(gas, wall)
    for _ in 1:8
        Tout = gas + fraction * (wall - gas)
        integral = _air_exchange_primitive2D(Tout, wall) -
                   primitive_in
        residual = integral - target
        abs(residual) <= 1e-11 * max(target, 1.0) && break
        residual > 0.0 ?
            (upper_fraction = fraction) :
            (lower_fraction = fraction)
        derivative = V11.air_heat_capacity(Tout) /
                     max(1.0 - fraction, 1e-14)
        candidate = fraction - residual / derivative
        fraction = lower_fraction < candidate < upper_fraction ?
            candidate : 0.5 * (lower_fraction + upper_fraction)
    end
    Tout = gas + fraction * (wall - gas)
    heat = flow * (
        air_specific_enthalpy2D(Tout) -
        air_specific_enthalpy2D(gas)
    )
    return (temperature=Tout, heat_W=heat)
end

function _march_channel_group2D!(
    work, group, mdot_channel, inlet, channel_temperature,
    p14, exchange, grid, ua_per_channel,
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
        V19._kernel_integrals2D!(
            kernel, view(work.gas, group, :),
            view(channel_temperature, group, :),
            mdot_channel, p14, exchange, grid,
        )
        ua_cell .= ua_per_channel .* kernel ./
                   max(sum(kernel), eps())
        work.gas[group, 1] = inlet
        for j in 1:bg.nz
            exchange_cell = _enthalpy_cell_exchange2D(
                work.gas[group, j],
                channel_temperature[group, j],
                ua_cell[j], mdot_channel,
            )
            work.gas[group, j + 1] = exchange_cell.temperature
        end
    end

    for j in 1:bg.nz
        Tgas = V11.clamp_temperature(work.gas[group, j])
        Tout = V11.clamp_temperature(work.gas[group, j + 1])
        cp = V11.air_heat_capacity(Tgas)
        kg = V11.air_conductivity(Tgas)
        mu = V11.air_viscosity(Tgas)
        work.qgas[group, j] = multiplicity * mdot_channel * (
            air_specific_enthalpy2D(Tout) -
            air_specific_enthalpy2D(Tgas)
        )
        work.h[group, j] = ua_cell[j] /
            max(perimeter * bg.dz_cells[j], eps())
        work.nusselt[group, j] =
            work.h[group, j] * width / max(kg, eps())
        work.prandtl[group, j] = cp * mu / kg

        Tbulk = V11.clamp_temperature(0.5 * (Tgas + Tout))
        rho = V11.air_density(Tbulk, hyd.atmospheric_pressure)
        mu_bulk = V11.air_viscosity(Tbulk)
        velocity = mdot_channel / max(rho * area, eps())
        work.density[group, j] = rho
        work.velocity[group, j] = velocity
        work.reynolds[group, j] =
            rho * velocity * width / mu_bulk
        work.graetz[group, j] =
            work.prandtl[group, j] *
            work.reynolds[group, j] * width /
            max(bg.z_centers[j], eps())
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

function _mixed_outlet2D(work, grid)
    total_mdot = sum(
        grid.multiplicity[group] * work.channel_mdot[group]
        for group in 1:grid.group_count
    )
    total_mdot <= eps() && return NaN
    outlet_min = minimum(work.gas[:, end])
    outlet_max = maximum(work.gas[:, end])
    mixed_h = sum(
        grid.multiplicity[group] * work.channel_mdot[group] *
        air_specific_enthalpy2D(work.gas[group, end])
        for group in 1:grid.group_count
    ) / total_mdot
    return air_temperature_from_enthalpy2D(
        mixed_h; lower=outlet_min, upper=outlet_max,
    )
end

function _solve_channel_flows2D!(
    work, mdot_total, inlet, channel_temperature,
    p14, exchange, grid,
)
    ng = grid.group_count
    hyd = p14.hydraulics
    work.channel_mdot .= mdot_total / sum(grid.multiplicity)
    relaxation = clamp(hyd.equal_pressure_relaxation, 0.05, 1.0)
    tolerance = max(hyd.equal_pressure_tolerance, 1e-10)
    max_iterations = max(hyd.equal_pressure_max_iterations, 24)
    gas_reference = inlet
    wall_profile = V19._mean_wall_profile2D(channel_temperature, grid)
    for iteration in 1:max_iterations
        average = V19._average_nusselt2D(
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
        ) / max(common, eps())
        work.common_dp = common
        work.equal_pressure_error = error
        if error <= tolerance && gas_error <= tolerance
            return (
                gas_reference_relative_error=gas_error,
                iterations=iteration, converged=true,
            )
        end
        linear_coeff = zeros(ng)
        quadratic_coeff = zeros(ng)
        for group in 1:ng
            mdot = max(work.channel_mdot[group], eps())
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
    average = V19._average_nusselt2D(
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
    ) / max(work.common_dp, eps())
    new_gas_reference = 0.5 * (inlet + _mixed_outlet2D(work, grid))
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

function _gas_profile_enthalpy2D!(
    work, channel_temperature, tube_temperature, time,
    p14, exchange, op, grid,
)
    bg = grid.base_grid
    flow_lpm = max(0.0, Float64(op.flow_lpm(time)))
    inlet = Float64(op.inlet_temperature(time))
    mdot_total = V11.standard_mass_flow2D(
        flow_lpm, p14.hydraulics,
    )
    if mdot_total <= eps()
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
            iterations=0, converged=true,
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
        reynolds = mdot_total * dh / max(area * mu, eps())
        prandtl = cp * mu / kf
        nusselt = reynolds < 2300.0 ?
            3.66 : 0.023 * reynolds^0.8 * prandtl^0.4
        work.hrear[j] = nusselt * kf / dh
        ua = work.hrear[j] * perimeter * bg.dz_rear
        exchange_cell = _enthalpy_cell_exchange2D(
            work.gas_rear[j], tube_temperature[j],
            ua, mdot_total,
        )
        work.gas_rear[j + 1] = exchange_cell.temperature
        work.qgas_rear[j] = exchange_cell.heat_W
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

function _v20_ode2D!(du, u, context, time)
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
    _gas_profile_enthalpy2D!(
        context.work, Tch, Ttube, time, context.p14,
        context.exchange, context.op, grid,
    )
    duch = reshape(
        view(base_du, base_layout.channel),
        grid.group_count, bg.nz,
    )
    for group in 1:grid.group_count, j in 1:bg.nz
        full_area = grid.multiplicity[group] *
                    grid.solid_area_per_channel
        core_area = full_area -
                    context.skin_group_area[group]
        capacity = V15._channel_capacity2D(
            Tch[group, j], core_area, bg.dz_cells[j],
            context.p15,
        )
        delta = context.work.qgas[group, j] -
                old.qgas[group, j]
        j == 1 && (delta +=
            context.work.front_qgas[group] - old.front[group])
        duch[group, j] -= delta / max(capacity, eps())
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
        ) / max(capacity, eps())
        overlap = V12._overlap_length(
            bg.z_rear_faces[j], bg.z_rear_faces[j + 1],
            context.tube_flange.contact_start_m,
            context.p14.geometry.rear_tube_length,
        )
        if overlap > 0.0 &&
           context.tube_flange.contact_h_W_m2_K > 0.0
            area_contact = 2.0 * pi *
                context.p14.geometry.rear_tube_outer_radius *
                overlap
            Q = context.tube_flange.contact_h_W_m2_K *
                area_contact *
                (Ttube[j] - V11.WATER_FLANGE_TEMP)
            dutube[j] -= Q / max(capacity, eps())
        end
    end
    return nothing
end

function _model_context2D(p::ModelParameters2D, op)
    v19 = V19._model_context2D(p.base, op)
    return merge(v19, (
        p=p,
        exchange=p.integrated_exchange,
        tube_flange=p.rear_tube_flange,
    ))
end

function _enthalpy_residuals2D(work, inlet, mdot_total, grid)
    receiver_gain = sum(work.qgas) + sum(work.front_qgas)
    orbit_gain = sum(
        grid.multiplicity[group] * work.channel_mdot[group] * (
            air_specific_enthalpy2D(work.gas[group, end]) -
            air_specific_enthalpy2D(inlet)
        ) for group in 1:grid.group_count
    )
    mixed = work.gas_rear[1]
    mixing_residual = mdot_total *
        air_specific_enthalpy2D(mixed) - sum(
            grid.multiplicity[group] * work.channel_mdot[group] *
            air_specific_enthalpy2D(work.gas[group, end])
            for group in 1:grid.group_count
        )
    rear_gain = sum(work.qgas_rear)
    rear_stream_gain = mdot_total * (
        air_specific_enthalpy2D(work.gas_rear[end]) -
        air_specific_enthalpy2D(mixed)
    )
    total_gain = receiver_gain + rear_gain
    total_stream_gain = mdot_total * (
        air_specific_enthalpy2D(work.gas_rear[end]) -
        air_specific_enthalpy2D(inlet)
    )
    return (
        receiver_enthalpy_residual_W=receiver_gain - orbit_gain,
        mixing_enthalpy_residual_W=mixing_residual,
        rear_enthalpy_residual_W=rear_gain - rear_stream_gain,
        total_gas_enthalpy_residual_W=
            total_gain - total_stream_gain,
        total_gas_power_W=total_gain,
    )
end

function _extract_diagnostics2D(
    solution, p, op, grid, layout, p14,
)
    # Retain the v19 field names, then replace every gas/hydraulic field with
    # values recomputed by the v20 enthalpy profile.
    inherited = V19._extract_diagnostics2D(
        solution, p.base, op, grid, layout, p14,
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
    receiver_residual = zeros(nt)
    mixing_residual = zeros(nt)
    rear_residual = zeros(nt)
    total_residual = zeros(nt)
    total_gas_power = zeros(nt)
    for index in 1:nt
        base_u = view(solution.u[index], layout.base_state)
        Tch = reshape(
            view(base_u, layout.base.channel), ng, bg.nz,
        )
        Ttube = view(base_u, layout.base.tube)
        convergence = _gas_profile_enthalpy2D!(
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
        mdot = V11.standard_mass_flow2D(
            max(0.0, Float64(op.flow_lpm(solution.t[index]))),
            p14.hydraulics,
        )
        mass_flow[index] = mdot
        residuals = _enthalpy_residuals2D(
            work, Float64(op.inlet_temperature(solution.t[index])),
            mdot, grid,
        )
        receiver_residual[index] =
            residuals.receiver_enthalpy_residual_W
        mixing_residual[index] =
            residuals.mixing_enthalpy_residual_W
        rear_residual[index] =
            residuals.rear_enthalpy_residual_W
        total_residual[index] =
            residuals.total_gas_enthalpy_residual_W
        total_gas_power[index] = residuals.total_gas_power_W
    end
    return merge(inherited, (
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
        receiver_enthalpy_residual_W=receiver_residual,
        mixing_enthalpy_residual_W=mixing_residual,
        rear_enthalpy_residual_W=rear_residual,
        total_gas_enthalpy_residual_W=total_residual,
        total_gas_power_W=total_gas_power,
    ))
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

function simulate2D(
    p::ModelParameters2D,
    op::OperatingCondition2D,
    save_times::AbstractVector;
    initial_temperature=Float64(
        op.ambient_temperature(save_times[1]),
    ),
    solver=FBDF(autodiff=AutoFiniteDiff()),
    reltol=1e-6, abstol=1e-7, dtmax=30.0,
)
    context = _model_context2D(p, op)
    grid = context.grid
    layout = context.layout
    p17 = p.base.base.base
    p15 = context.p15
    p14 = context.p14
    initial = V15._initial_state2D(
        initial_temperature, p15, grid,
    )
    problem = ODEProblem(
        _v20_ode2D!, initial,
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

function _sample_receiver_mixed_gas2D(result, z)
    grid = build_network_grid2D(result.parameters)
    zfaces = grid.base_grid.z_faces
    j0, j1, weight = V11._axis_interp_indices(zfaces, z)
    nt = length(result.time)
    output = zeros(nt)
    for index in 1:nt
        total = sum(result.group_mass_flow_kg_s[:, index])
        if total <= eps()
            output[index] = sum(
                grid.multiplicity[group] *
                ((1.0 - weight) *
                 result.gas_temperature[group, j0, index] +
                 weight *
                 result.gas_temperature[group, j1, index])
                for group in 1:grid.group_count
            ) / sum(grid.multiplicity)
            continue
        end
        mixed_h = sum(
            result.group_mass_flow_kg_s[group, index] *
            air_specific_enthalpy2D(
                (1.0 - weight) *
                result.gas_temperature[group, j0, index] +
                weight *
                result.gas_temperature[group, j1, index]
            ) for group in 1:grid.group_count
        ) / total
        local_temperatures = [
            (1.0 - weight) *
            result.gas_temperature[group, j0, index] +
            weight * result.gas_temperature[group, j1, index]
            for group in 1:grid.group_count
        ]
        output[index] = air_temperature_from_enthalpy2D(
            mixed_h; lower=minimum(local_temperatures),
            upper=maximum(local_temperatures),
        )
    end
    return output
end

function _sample_receiver_mean_wall2D(result, z)
    grid = build_network_grid2D(result.parameters)
    j0, j1, weight = V11._axis_interp_indices(
        result.z_solid, z,
    )
    denominator = sum(grid.multiplicity)
    return vec(sum(
        grid.multiplicity[group] .* (
            (1.0 - weight) .*
            result.channel_temperature[group, j0, :] .+
            weight .*
            result.channel_temperature[group, j1, :]
        ) for group in 1:grid.group_count
    ) ./ denominator)
end

function _t3_targets2D(result)
    location = result.parameters.t3_location.location
    length_receiver = result.parameters.geometry.length
    if location === :receiver_136
        z = min(136e-3, length_receiver)
        return (
            gas=_sample_receiver_mixed_gas2D(result, z),
            surrounding=_sample_receiver_mean_wall2D(result, z),
            global_z_m=z,
        )
    elseif location === :receiver_exit
        z = length_receiver
        return (
            gas=_sample_receiver_mixed_gas2D(result, z),
            surrounding=_sample_receiver_mean_wall2D(result, z),
            global_z_m=z,
        )
    elseif location === :rear_003
        z = 3e-3
        return (
            gas=V19._sample_rear_gas2D(result, z),
            surrounding=V19._sample_rear_tube2D(result, z),
            global_z_m=length_receiver + z,
        )
    end
    error("unknown v20 T3 location: $location")
end

function _probe_gas_coefficient2D(
    temperature, mdot, p::ModelParameters2D,
)
    diameter = p.outlet_observation.probe_diameter_m
    diameter > 0.0 || error(
        "v20 probe diameter must be positive",
    )
    tube_area = pi * p.geometry.rear_tube_inner_radius^2
    rho = V11.air_density(
        temperature, p.hydraulics.atmospheric_pressure,
    )
    velocity = mdot / max(rho * tube_area, eps())
    mu = V11.air_viscosity(temperature)
    conductivity = V11.air_conductivity(temperature)
    cp = V11.air_heat_capacity(temperature)
    reynolds = rho * velocity * diameter / mu
    prandtl = cp * mu / conductivity
    numerator = 0.62 * sqrt(max(reynolds, 0.0)) *
        prandtl^(1.0 / 3.0)
    denominator = (
        1.0 + (0.4 / max(prandtl, eps()))^(2.0 / 3.0)
    )^(1.0 / 4.0)
    nusselt = 0.3 + numerator / denominator *
        (
            1.0 + (reynolds / 282000.0)^(5.0 / 8.0)
        )^(4.0 / 5.0)
    return nusselt * conductivity / diameter
end

function _probe_history2D(
    result, Tgas, Tsurrounding; initial_temperature=nothing,
)
    parameters = result.parameters.outlet_observation
    capacity = parameters.capacity_areal_J_m2_K
    stem = parameters.stem_conductance_areal_W_m2_K
    emissivity = parameters.emissivity
    capacity > 0.0 || error(
        "v20 probe areal capacity must be positive",
    )
    stem >= 0.0 || error(
        "v20 probe stem conductance must be nonnegative",
    )
    0.0 <= emissivity <= 1.0 || error(
        "v20 probe emissivity must lie in [0,1]",
    )
    time = result.time
    output = similar(Tgas)
    output[1] = initial_temperature === nothing ?
        Tgas[1] : Float64(initial_temperature)
    sink = V11.WATER_FLANGE_TEMP
    for index in 2:length(time)
        dt = time[index] - time[index - 1]
        dt >= 0.0 || error("v20 probe time must be monotone")
        Tprevious = output[index - 1]
        Tgas_mid = 0.5 * (Tgas[index - 1] + Tgas[index])
        Tsurrounding_mid = 0.5 * (
            Tsurrounding[index - 1] + Tsurrounding[index]
        )
        mdot_mid = 0.5 * (
            result.mass_flow_kg_s[index - 1] +
            result.mass_flow_kg_s[index]
        )
        hgas = _probe_gas_coefficient2D(
            Tgas_mid, mdot_mid, result.parameters,
        )
        hrad = emissivity * V11.SIGMA *
            (Tsurrounding_mid^2 + Tprevious^2) *
            (Tsurrounding_mid + Tprevious)
        total = hgas + hrad + stem
        target = (
            hgas * Tgas_mid +
            hrad * Tsurrounding_mid + stem * sink
        ) / max(total, eps())
        output[index] = target + (Tprevious - target) *
            exp(-total * dt / capacity)
    end
    return output
end

function sensor_predictions2D(
    result::SimulationResult2D; initial_T3=nothing,
)
    obs = result.parameters.observation
    time = result.time
    side11 = V19._sample_skin2D(
        result, V11.SIDE_TC_FRONT_Z_2D,
    )
    side58 = V19._sample_skin2D(result, 58e-3)
    side107 = V19._sample_skin2D(result, 107e-3)
    wall58 = V19._sample_channel2D(
        result, result.core_group, 58e-3,
    )
    wall107 = V19._sample_channel2D(
        result, result.core_group, 107e-3,
    )
    gas58 = V19._sample_channel_gas2D(
        result, result.core_group, 58e-3,
    )
    gas107 = V19._sample_channel_gas2D(
        result, result.core_group, 107e-3,
    )
    wall_fraction = clamp(obs.interior_wall_fraction, 0.0, 1.0)
    target9 = wall_fraction .* wall58 +
              (1.0 - wall_fraction) .* gas58
    target10 = wall_fraction .* wall107 +
               (1.0 - wall_fraction) .* gas107
    t3 = _t3_targets2D(result)
    # The v19 probe routine is generic: its second thermal target is the
    # surrounding receiver wall for receiver branches and the alumina tube
    # for the downstream branch.
    T3probe = _probe_history2D(
        result, t3.gas, t3.surrounding;
        initial_temperature=initial_T3,
    )
    rT2 = result.parameters.geometry.receiver_radius + 40e-3
    T2 = V14._sample_outer2D(
        result.base_result, rT2, 58e-3,
    )
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
        T3=T3probe, T2=T2,
        T8_skin=side11, T12_skin=side58, T11_skin=side107,
        T9_wall=wall58, T10_wall=wall107,
        T9_gas=gas58, T10_gas=gas107,
        T3_gas_raw=t3.gas,
        T3_surrounding=t3.surrounding,
        T3_probe=T3probe,
        T3_global_z_m=t3.global_z_m,
    )
end

function enthalpy_transport_audit2D(
    result::SimulationResult2D,
)
    absolute_max = maximum(abs, vcat(
        result.receiver_enthalpy_residual_W,
        result.mixing_enthalpy_residual_W,
        result.rear_enthalpy_residual_W,
        result.total_gas_enthalpy_residual_W,
    ))
    scale = max(maximum(abs, result.total_gas_power_W), 1.0)
    return (
        maximum_absolute_residual_W=absolute_max,
        maximum_relative_residual=absolute_max / scale,
        receiver_max_W=maximum(
            abs, result.receiver_enthalpy_residual_W,
        ),
        mixing_max_W=maximum(
            abs, result.mixing_enthalpy_residual_W,
        ),
        rear_max_W=maximum(
            abs, result.rear_enthalpy_residual_W,
        ),
        total_max_W=maximum(
            abs, result.total_gas_enthalpy_residual_W,
        ),
    )
end

function exchange_energy_audit2D(
    u::AbstractVector, p::ModelParameters2D,
    op::OperatingCondition2D, time::Real=0.0,
)
    baseline_context = _model_context2D(p, op)
    v20_context = _model_context2D(p, op)
    du_baseline = zeros(length(u))
    du_v20 = zeros(length(u))
    V17._casing_flange_ode2D!(
        du_baseline, u, baseline_context.casing_context, time,
    )
    _v20_ode2D!(du_v20, u, v20_context, time)
    grid = v20_context.grid
    bg = grid.base_grid
    layout = v20_context.layout
    base_u = view(u, layout.base_state)
    Tch = reshape(
        view(base_u, layout.base.channel),
        grid.group_count, bg.nz,
    )
    Ttube = view(base_u, layout.base.tube)
    delta = view(du_v20, layout.base_state) .-
            view(du_baseline, layout.base_state)
    delta_channel = reshape(
        view(delta, layout.base.channel),
        grid.group_count, bg.nz,
    )
    derivative_power = 0.0
    for group in 1:grid.group_count, j in 1:bg.nz
        full_area = grid.multiplicity[group] *
                    grid.solid_area_per_channel
        core_area = full_area -
                    v20_context.skin_group_area[group]
        capacity = V15._channel_capacity2D(
            Tch[group, j], core_area, bg.dz_cells[j],
            v20_context.p15,
        )
        derivative_power += capacity *
                            delta_channel[group, j]
    end
    tube_area = pi * (
        v20_context.p14.geometry.rear_tube_outer_radius^2 -
        v20_context.p14.geometry.rear_tube_inner_radius^2
    )
    delta_tube = view(delta, layout.base.tube)
    for j in 1:bg.nz_rear
        capacity = 3900.0 *
            V11.alumina_heat_capacity(Ttube[j]) *
            tube_area * bg.dz_rear
        derivative_power += capacity * delta_tube[j]
    end
    old = baseline_context.work
    new = v20_context.work
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
    expected = -gas_correction - tube_flange_loss
    residual = derivative_power - expected
    return (
        capacity_derivative_correction_W=derivative_power,
        expected_correction_W=expected,
        gas_correction_W=gas_correction,
        distributed_tube_flange_loss_W=tube_flange_loss,
        residual_W=residual,
        relative_residual=abs(residual) / max(abs(expected), 1.0),
    )
end

end
