module Receiver1D

"""
Fast conservative 1D model for the small SiC solar receiver.

Coordinate convention: z = 0 is the irradiated air inlet and z = L is the
rear gas outlet. The only dynamic states are axial solid control volumes. The
gas is evaluated as quasi-steady plug flow in every RHS evaluation.
"""

export Geometry, SolidProperties, HeatTransferParameters, LossParameters,
       OpticalParameters, ModelParameters, OperatingCondition, SimulationResult,
       default_parameters, effective_htc_initial, literature_htc,
       constant_history, linear_history, simulate, gas_profile, solid_at,
       heat_transfer_summary, energy_rates, sic_conductivity, sic_heat_capacity,
       air_conductivity, air_heat_capacity, air_viscosity

const SIGMA = 5.670374419e-8

Base.@kwdef struct Geometry
    length::Float64 = 137e-3
    receiver_width::Float64 = 19e-3
    channel_width::Float64 = 1.5e-3
    wall_thickness::Float64 = 0.4e-3
    channel_count::Int = 100
end

Base.@kwdef struct SolidProperties
    density::Float64 = 3200.0
    conductivity_scale::Float64 = 1.0
    heat_capacity_scale::Float64 = 1.0
end

"""
Correlation for the local Nusselt number.

`Nu_bar(z) = C Re^m Pr^n (Tfilm/Tref)^r (1 + a Dh Re Pr/z)^p`

The code differentiates `z*Nu_bar(z)` to obtain a local coefficient. With
`a = 0`, the correlation reduces to a spatially uniform power law.
"""
Base.@kwdef struct HeatTransferParameters
    coefficient::Float64 = 6.1e-4
    reynolds_exponent::Float64 = 1.3835
    prandtl_exponent::Float64 = 0.0
    temperature_exponent::Float64 = -0.4883
    reference_temperature::Float64 = 600.0
    development_coefficient::Float64 = 0.0
    development_exponent::Float64 = 0.45
    minimum_nusselt::Float64 = 1e-5
end

Base.@kwdef struct LossParameters
    front_convection::Float64 = 10.0
    front_emissivity::Float64 = 0.85
    side_conductance_per_length::Float64 = 0.35
    rear_conductance::Float64 = 0.10
    rear_emissivity::Float64 = 0.0
end

Base.@kwdef struct OpticalParameters
    absorbed_fraction::Float64 = 0.85
    extinction_coefficient::Float64 = 50.0
    front_deposition_fraction::Float64 = 0.0
end

Base.@kwdef struct ModelParameters
    geometry::Geometry = Geometry()
    solid::SolidProperties = SolidProperties()
    heat_transfer::HeatTransferParameters = HeatTransferParameters()
    losses::LossParameters = LossParameters()
    optics::OpticalParameters = OpticalParameters()
    nodes::Int = 25
end

struct OperatingCondition{FQ,FF,FT,FA}
    irradiance::FQ
    flow_lpm::FF
    inlet_temperature::FT
    ambient_temperature::FA
end

constant_history(value::Real) = (_ -> Float64(value))
constant_history(f) = f

function OperatingCondition(; irradiance=0.0, flow_lpm=0.0,
                              inlet_temperature=295.6,
                              ambient_temperature=295.6)
    return OperatingCondition(constant_history(irradiance),
                              constant_history(flow_lpm),
                              constant_history(inlet_temperature),
                              constant_history(ambient_temperature))
end

"""Return a piecewise-linear history with constant end extrapolation."""
function linear_history(times::AbstractVector, values::AbstractVector)
    length(times) == length(values) || throw(ArgumentError("times and values must have equal length"))
    length(times) >= 2 || throw(ArgumentError("at least two history points are required"))
    t = Float64.(times)
    y = Float64.(values)
    issorted(t) || throw(ArgumentError("history times must be sorted"))
    return function (x)
        x <= t[1] && return y[1]
        x >= t[end] && return y[end]
        i = searchsortedlast(t, x)
        fraction = (x - t[i]) / (t[i + 1] - t[i])
        return muladd(fraction, y[i + 1] - y[i], y[i])
    end
end

struct SimulationResult
    time::Vector{Float64}
    z_solid::Vector{Float64}
    z_gas::Vector{Float64}
    solid_temperature::Matrix{Float64}
    gas_temperature::Matrix{Float64}
    local_heat_transfer_coefficient::Matrix{Float64}
end

default_parameters() = ModelParameters()

"""
Initial effective correlation inferred from the 15 steady experiments by
inverting the measured outlet gas temperature against the external wall
profile. It is a calibration seed, not a canonical channel-film correlation.
"""
effective_htc_initial() = HeatTransferParameters()

"""Correlation 1 from the manuscript, expressed in local form."""
literature_htc() = HeatTransferParameters(
    coefficient=3.66,
    reynolds_exponent=0.0,
    prandtl_exponent=0.0,
    temperature_exponent=0.0,
    development_coefficient=0.095,
    development_exponent=0.45,
    minimum_nusselt=3.66,
)

frontal_area(g::Geometry) = g.receiver_width^2
channel_area(g::Geometry) = g.channel_count * g.channel_width^2
solid_area(g::Geometry) = frontal_area(g) - channel_area(g)
exchange_perimeter(g::Geometry) = 4.0 * g.channel_width * g.channel_count
hydraulic_diameter(g::Geometry) = g.channel_width

function validate(p::ModelParameters)
    g = p.geometry
    g.length > 0 || throw(ArgumentError("receiver length must be positive"))
    g.channel_count > 0 || throw(ArgumentError("channel_count must be positive"))
    solid_area(g) > 0 || throw(ArgumentError("channel open area exceeds frontal area"))
    p.nodes >= 3 || throw(ArgumentError("at least three axial nodes are required"))
    0.0 <= p.optics.absorbed_fraction <= 1.0 || throw(ArgumentError("absorbed_fraction must be in [0, 1]"))
    0.0 <= p.optics.front_deposition_fraction <= 1.0 || throw(ArgumentError("front_deposition_fraction must be in [0, 1]"))
    return p
end

clamp_temperature(T) = clamp(Float64(T), 250.0, 1600.0)

function sic_conductivity(T)
    x = clamp_temperature(T)
    return max(2.0, 191.9216 - 0.3261784x + 2.739462e-4x^2 - 7.70926e-8x^3)
end

function sic_heat_capacity(T)
    x = clamp_temperature(T)
    return max(400.0, 1110.0 + 0.15x - 425.0exp(-0.003x))
end

function air_conductivity(T)
    x = clamp_temperature(T)
    return max(0.01, -0.00227583562 + 1.15480022e-4x - 7.90252856e-8x^2 + 4.11702505e-11x^3 - 7.43864331e-15x^4)
end

function air_heat_capacity(T)
    x = clamp_temperature(T)
    return max(700.0, 1.93e-10x^4 - 8.0e-7x^3 + 1.14e-3x^2 - 0.449x + 1060.0)
end

function air_viscosity(T)
    x = clamp_temperature(T)
    return max(5e-6, -8.38278e-7 + 8.35717342e-8x - 7.69429583e-11x^2 + 4.6437266e-14x^3 - 1.06585607e-17x^4)
end

air_density_standard(T=295.6) = 352.716 / clamp_temperature(T)

function local_nusselt(reynolds, prandtl, z, dh, film_temperature,
                        hp::HeatTransferParameters)
    base = hp.coefficient * reynolds^hp.reynolds_exponent *
           prandtl^hp.prandtl_exponent *
           (film_temperature / hp.reference_temperature)^hp.temperature_exponent
    b = hp.development_coefficient * dh * reynolds * prandtl
    if b <= 0.0
        return max(hp.minimum_nusselt, base)
    end
    zz = max(z, eps(Float64))
    average_factor = (1.0 + b / zz)^hp.development_exponent
    derivative_factor = 1.0 - hp.development_exponent * b / (zz + b)
    return max(hp.minimum_nusselt, base * average_factor * derivative_factor)
end

function _gas_profile!(gas, qcell, hcell, solid_temperature, time,
                       p::ModelParameters, op::OperatingCondition, z, dx)
    g = p.geometry
    flow = max(0.0, Float64(op.flow_lpm(time)))
    inlet = Float64(op.inlet_temperature(time))
    gas[1] = inlet
    if flow <= 1e-12
        fill!(qcell, 0.0)
        fill!(hcell, 0.0)
        fill!(gas, inlet)
        return nothing
    end

    mass_flow = air_density_standard(inlet) * flow / 60000.0
    area_flow = channel_area(g)
    perimeter = exchange_perimeter(g)
    dh = hydraulic_diameter(g)
    hp = p.heat_transfer

    for i in eachindex(solid_temperature)
        film = clamp_temperature(0.5 * (solid_temperature[i] + gas[i]))
        cp = air_heat_capacity(film)
        kf = air_conductivity(film)
        mu = air_viscosity(film)
        reynolds = mass_flow * dh / (area_flow * mu)
        prandtl = cp * mu / kf
        nusselt = local_nusselt(reynolds, prandtl, z[i], dh, film, hp)
        hcell[i] = nusselt * kf / dh
        ua = hcell[i] * perimeter * dx
        effectiveness = -expm1(-ua / (mass_flow * cp))
        gas[i + 1] = gas[i] + effectiveness * (solid_temperature[i] - gas[i])
        qcell[i] = mass_flow * cp * (gas[i + 1] - gas[i])
    end
    return nothing
end

"""Evaluate the quasi-steady gas profile for a supplied axial solid profile."""
function gas_profile(solid_temperature::AbstractVector, time,
                     p::ModelParameters, op::OperatingCondition)
    validate(p)
    length(solid_temperature) == p.nodes || throw(DimensionMismatch("solid profile length must equal p.nodes"))
    dx = p.geometry.length / p.nodes
    z = collect(range(dx / 2, step=dx, length=p.nodes))
    gas = zeros(p.nodes + 1)
    qcell = zeros(p.nodes)
    hcell = zeros(p.nodes)
    _gas_profile!(gas, qcell, hcell, solid_temperature, time, p, op, z, dx)
    return (temperature=gas, heat_to_gas=qcell, h=hcell)
end

function optical_weights(p::ModelParameters, dx)
    n = p.nodes
    optics = p.optics
    weights = zeros(n)
    beta = optics.extinction_coefficient
    if !isfinite(beta) || beta > 1e8
        weights[1] = 1.0
    elseif beta <= 1e-12
        fill!(weights, 1.0 / n)
    else
        for i in 1:n
            left = (i - 1) * dx
            right = i * dx
            weights[i] = exp(-beta * left) - exp(-beta * right)
        end
        weights ./= sum(weights)
    end
    front = optics.front_deposition_fraction
    weights .*= 1.0 - front
    weights[1] += front
    return weights
end

mutable struct Workspace
    gas::Vector{Float64}
    qgas::Vector{Float64}
    h::Vector{Float64}
    k1::Vector{Float64}
    k2::Vector{Float64}
    k3::Vector{Float64}
    k4::Vector{Float64}
    trial::Vector{Float64}
end

Workspace(n) = Workspace(zeros(n + 1), zeros(n), zeros(n), zeros(n),
                         zeros(n), zeros(n), zeros(n), zeros(n))

function _rhs!(du, temperature, time, p::ModelParameters,
               op::OperatingCondition, z, dx, weights, work::Workspace)
    n = p.nodes
    g = p.geometry
    solid = p.solid
    loss = p.losses
    area_solid = solid_area(g)
    cell_volume = area_solid * dx
    ambient = Float64(op.ambient_temperature(time))

    _gas_profile!(work.gas, work.qgas, work.h, temperature, time, p, op, z, dx)
    fill!(du, 0.0)

    for i in 1:(n - 1)
        ki = solid.conductivity_scale * sic_conductivity(temperature[i])
        kj = solid.conductivity_scale * sic_conductivity(temperature[i + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        qcond = kface * area_solid * (temperature[i] - temperature[i + 1]) / dx
        du[i] -= qcond
        du[i + 1] += qcond
    end

    absorbed_power = p.optics.absorbed_fraction *
                     max(0.0, Float64(op.irradiance(time))) * frontal_area(g)
    for i in 1:n
        side_loss = loss.side_conductance_per_length * dx * (temperature[i] - ambient)
        du[i] += absorbed_power * weights[i] - work.qgas[i] - side_loss
    end

    front_radiation = loss.front_emissivity * SIGMA * frontal_area(g) *
                      (temperature[1]^4 - ambient^4)
    front_convection = loss.front_convection * frontal_area(g) *
                       (temperature[1] - ambient)
    rear_radiation = loss.rear_emissivity * SIGMA * frontal_area(g) *
                     (temperature[end]^4 - ambient^4)
    rear_conduction = loss.rear_conductance * (temperature[end] - ambient)
    du[1] -= front_convection + front_radiation
    du[end] -= rear_conduction + rear_radiation

    for i in 1:n
        capacity = solid.density * solid.heat_capacity_scale *
                   sic_heat_capacity(temperature[i]) * cell_volume
        du[i] /= capacity
    end
    return nothing
end

function stable_step(temperature, p::ModelParameters, dx, maximum_step)
    solid = p.solid
    alpha_max = maximum(solid.conductivity_scale * sic_conductivity(T) /
                        (solid.density * solid.heat_capacity_scale * sic_heat_capacity(T))
                        for T in temperature)
    conduction_limit = 0.35 * dx^2 / max(alpha_max, eps(Float64))
    return min(Float64(maximum_step), conduction_limit)
end

function initial_profile(initial_temperature, z)
    if initial_temperature isa Real
        return fill(Float64(initial_temperature), length(z))
    elseif initial_temperature isa AbstractVector
        length(initial_temperature) == length(z) || throw(DimensionMismatch("initial profile length must equal number of nodes"))
        return Float64.(initial_temperature)
    else
        return Float64[initial_temperature(x) for x in z]
    end
end

"""
Simulate the receiver with an explicit RK4 integrator and a conservative
finite-volume spatial discretization. `maximum_step` is automatically reduced
to satisfy the axial-conduction stability limit.
"""
function simulate(p::ModelParameters, op::OperatingCondition,
                  save_times::AbstractVector;
                  initial_temperature=op.ambient_temperature(first(save_times)),
                  maximum_step=2.0)
    validate(p)
    length(save_times) >= 1 || throw(ArgumentError("save_times cannot be empty"))
    times = Float64.(save_times)
    issorted(times) || throw(ArgumentError("save_times must be sorted"))
    all(diff(times) .> 0.0) || length(times) == 1 || throw(ArgumentError("save_times must be strictly increasing"))

    n = p.nodes
    dx = p.geometry.length / n
    z = collect(range(dx / 2, step=dx, length=n))
    zg = collect(range(0.0, p.geometry.length, length=n + 1))
    weights = optical_weights(p, dx)
    temperature = initial_profile(initial_temperature, z)
    solid_history = Matrix{Float64}(undef, n, length(times))
    gas_history = Matrix{Float64}(undef, n + 1, length(times))
    h_history = Matrix{Float64}(undef, n, length(times))
    work = Workspace(n)

    current_time = times[1]
    solid_history[:, 1] .= temperature
    _gas_profile!(work.gas, work.qgas, work.h, temperature, current_time, p, op, z, dx)
    gas_history[:, 1] .= work.gas
    h_history[:, 1] .= work.h

    for output_index in 2:length(times)
        target_time = times[output_index]
        while current_time < target_time
            dt = min(stable_step(temperature, p, dx, maximum_step), target_time - current_time)

            _rhs!(work.k1, temperature, current_time, p, op, z, dx, weights, work)
            @. work.trial = temperature + 0.5 * dt * work.k1
            _rhs!(work.k2, work.trial, current_time + 0.5dt, p, op, z, dx, weights, work)
            @. work.trial = temperature + 0.5 * dt * work.k2
            _rhs!(work.k3, work.trial, current_time + 0.5dt, p, op, z, dx, weights, work)
            @. work.trial = temperature + dt * work.k3
            _rhs!(work.k4, work.trial, current_time + dt, p, op, z, dx, weights, work)
            @. temperature += dt * (work.k1 + 2.0work.k2 + 2.0work.k3 + work.k4) / 6.0
            current_time += dt
        end
        solid_history[:, output_index] .= temperature
        _gas_profile!(work.gas, work.qgas, work.h, temperature, current_time, p, op, z, dx)
        gas_history[:, output_index] .= work.gas
        h_history[:, output_index] .= work.h
    end
    return SimulationResult(times, z, zg, solid_history, gas_history, h_history)
end

"""Linearly interpolate the predicted solid temperature to a physical z location."""
function solid_at(result::SimulationResult, location::Real)
    z = result.z_solid
    x = clamp(Float64(location), z[1], z[end])
    i = searchsortedlast(z, x)
    i >= length(z) && return copy(result.solid_temperature[end, :])
    i < 1 && return copy(result.solid_temperature[1, :])
    fraction = (x - z[i]) / (z[i + 1] - z[i])
    return @. result.solid_temperature[i, :] + fraction *
              (result.solid_temperature[i + 1, :] - result.solid_temperature[i, :])
end

function heat_transfer_summary(result::SimulationResult, index=length(result.time))
    h = view(result.local_heat_transfer_coefficient, :, index)
    return (mean_h=sum(h) / length(h), minimum_h=minimum(h), maximum_h=maximum(h),
            gas_inlet=result.gas_temperature[1, index],
            gas_outlet=result.gas_temperature[end, index])
end

"""Return instantaneous power terms and the finite-volume energy residual."""
function energy_rates(temperature::AbstractVector, time, p::ModelParameters,
                      op::OperatingCondition)
    validate(p)
    length(temperature) == p.nodes || throw(DimensionMismatch("temperature length must equal p.nodes"))
    n = p.nodes
    g = p.geometry
    dx = g.length / n
    z = collect(range(dx / 2, step=dx, length=n))
    weights = optical_weights(p, dx)
    work = Workspace(n)
    derivative = zeros(n)
    _rhs!(derivative, temperature, time, p, op, z, dx, weights, work)

    ambient = Float64(op.ambient_temperature(time))
    absorbed = p.optics.absorbed_fraction * max(0.0, Float64(op.irradiance(time))) * frontal_area(g)
    gas = sum(work.qgas)
    side = p.losses.side_conductance_per_length * dx * sum(temperature .- ambient)
    front = p.losses.front_convection * frontal_area(g) * (temperature[1] - ambient) +
            p.losses.front_emissivity * SIGMA * frontal_area(g) * (temperature[1]^4 - ambient^4)
    rear = p.losses.rear_conductance * (temperature[end] - ambient) +
           p.losses.rear_emissivity * SIGMA * frontal_area(g) * (temperature[end]^4 - ambient^4)
    cell_volume = solid_area(g) * dx
    storage = sum(p.solid.density * p.solid.heat_capacity_scale *
                  sic_heat_capacity(temperature[i]) * cell_volume * derivative[i]
                  for i in 1:n)
    residual = absorbed - gas - side - front - rear - storage
    return (absorbed=absorbed, gas=gas, side=side, front=front, rear=rear,
            storage=storage, residual=residual)
end

end # module Receiver1D
