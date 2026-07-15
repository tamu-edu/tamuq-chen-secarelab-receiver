# ============================================================================
# 1D_v2.jl - Axial finite-volume solid/gas solar receiver model
# ============================================================================
# Standalone refactor using the research-script flow of 1D_v1.exp.All.jl and
# 0D_v4.jl. This file owns its numerical core and does not include 1D_v1.jl.
# It provides:
#   A. libraries and execution flags
#   B. geometry and parameter ordering
#   C. model construction
#   D. external experimental data import
#   E. solver/remake interface
#   F. cooling-first and heating-second fitting
#   G. post-processing and metrics
#
# Default execution only loads data and defines functions. Set environment
# variables RECEIVER1D_RUN_OPTIMIZATION, RECEIVER1D_RUN_POSTPROCESS,
# RECEIVER1D_RUN_PLOTS, or RECEIVER1D_EXPORT_METRICS to "true" to run those
# expensive/side-effecting sections.
# ============================================================================

begin # libraries
    using Statistics
end

# ============================================================================
# STANDALONE NUMERICAL CORE
# ============================================================================
module Receiver1DV2

export Geometry, SolidProperties, HeatTransferParameters, LossParameters,
       OpticalParameters, ModelParameters, OperatingCondition, SimulationResult,
       linear_history, simulate, solid_at, gas_profile, heat_transfer_summary,
       energy_rates

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

frontal_area(g::Geometry) = g.receiver_width^2
channel_area(g::Geometry) = g.channel_count * g.channel_width^2
solid_area(g::Geometry) = frontal_area(g) - channel_area(g)
exchange_perimeter(g::Geometry) = 4.0 * g.channel_width * g.channel_count
hydraulic_diameter(g::Geometry) = g.channel_width

function validate(parameters::ModelParameters)
    geometry = parameters.geometry
    geometry.length > 0.0 || throw(ArgumentError("receiver length must be positive"))
    geometry.channel_count > 0 || throw(ArgumentError("channel_count must be positive"))
    solid_area(geometry) > 0.0 || throw(ArgumentError("channel area exceeds frontal area"))
    parameters.nodes >= 3 || throw(ArgumentError("at least three axial nodes are required"))
    0.0 <= parameters.optics.absorbed_fraction <= 1.0 ||
        throw(ArgumentError("absorbed_fraction must lie in [0, 1]"))
    return parameters
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
    return max(0.01, -0.00227583562 + 1.15480022e-4x - 7.90252856e-8x^2 +
                     4.11702505e-11x^3 - 7.43864331e-15x^4)
end

function air_heat_capacity(T)
    x = clamp_temperature(T)
    return max(700.0, 1.93e-10x^4 - 8.0e-7x^3 + 1.14e-3x^2 - 0.449x + 1060.0)
end

function air_viscosity(T)
    x = clamp_temperature(T)
    return max(5e-6, -8.38278e-7 + 8.35717342e-8x - 7.69429583e-11x^2 +
                    4.6437266e-14x^3 - 1.06585607e-17x^4)
end

air_density_standard(T=295.6) = 352.716 / clamp_temperature(T)

function local_nusselt(reynolds, prandtl, z, dh, film_temperature,
                        parameters::HeatTransferParameters)
    base = parameters.coefficient * reynolds^parameters.reynolds_exponent *
           prandtl^parameters.prandtl_exponent *
           (film_temperature / parameters.reference_temperature)^parameters.temperature_exponent
    b = parameters.development_coefficient * dh * reynolds * prandtl
    b <= 0.0 && return max(parameters.minimum_nusselt, base)
    zz = max(z, eps(Float64))
    average_factor = (1.0 + b / zz)^parameters.development_exponent
    derivative_factor = 1.0 - parameters.development_exponent * b / (zz + b)
    return max(parameters.minimum_nusselt, base * average_factor * derivative_factor)
end

function gas_profile!(gas, qcell, hcell, solid_temperature, time,
                      parameters::ModelParameters, operating::OperatingCondition,
                      z, dx)
    geometry = parameters.geometry
    flow = max(0.0, Float64(operating.flow_lpm(time)))
    inlet = Float64(operating.inlet_temperature(time))
    gas[1] = inlet
    if flow <= 1e-12
        fill!(gas, inlet)
        fill!(qcell, 0.0)
        fill!(hcell, 0.0)
        return nothing
    end

    mass_flow = air_density_standard(inlet) * flow / 60000.0
    area_flow = channel_area(geometry)
    perimeter = exchange_perimeter(geometry)
    dh = hydraulic_diameter(geometry)
    htc = parameters.heat_transfer

    for i in eachindex(solid_temperature)
        film = clamp_temperature(0.5 * (solid_temperature[i] + gas[i]))
        cp = air_heat_capacity(film)
        kf = air_conductivity(film)
        mu = air_viscosity(film)
        reynolds = mass_flow * dh / (area_flow * mu)
        prandtl = cp * mu / kf
        nu = local_nusselt(reynolds, prandtl, z[i], dh, film, htc)
        hcell[i] = nu * kf / dh
        ua = hcell[i] * perimeter * dx
        effectiveness = -expm1(-ua / (mass_flow * cp))
        gas[i + 1] = gas[i] + effectiveness * (solid_temperature[i] - gas[i])
        qcell[i] = mass_flow * cp * (gas[i + 1] - gas[i])
    end
    return nothing
end

function gas_profile(solid_temperature::AbstractVector, time,
                     parameters::ModelParameters, operating::OperatingCondition)
    validate(parameters)
    length(solid_temperature) == parameters.nodes ||
        throw(DimensionMismatch("solid profile length must equal model nodes"))
    dx = parameters.geometry.length / parameters.nodes
    z = collect(range(dx / 2, step=dx, length=parameters.nodes))
    gas = zeros(parameters.nodes + 1)
    qcell = zeros(parameters.nodes)
    hcell = zeros(parameters.nodes)
    gas_profile!(gas, qcell, hcell, solid_temperature, time,
                 parameters, operating, z, dx)
    return (temperature=gas, heat_to_gas=qcell, h=hcell)
end

function optical_weights(parameters::ModelParameters, dx)
    n = parameters.nodes
    optics = parameters.optics
    weights = zeros(n)
    beta = optics.extinction_coefficient
    if !isfinite(beta) || beta > 1e8
        weights[1] = 1.0
    elseif beta <= 1e-12
        fill!(weights, 1.0 / n)
    else
        for i in 1:n
            weights[i] = exp(-beta * (i - 1) * dx) - exp(-beta * i * dx)
        end
        weights ./= sum(weights)
    end
    weights .*= 1.0 - optics.front_deposition_fraction
    weights[1] += optics.front_deposition_fraction
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

function rhs!(du, temperature, time, parameters::ModelParameters,
              operating::OperatingCondition, z, dx, weights, work::Workspace)
    n = parameters.nodes
    geometry = parameters.geometry
    solid = parameters.solid
    losses = parameters.losses
    area_solid = solid_area(geometry)
    cell_volume = area_solid * dx
    ambient = Float64(operating.ambient_temperature(time))

    gas_profile!(work.gas, work.qgas, work.h, temperature, time,
                 parameters, operating, z, dx)
    fill!(du, 0.0)

    for i in 1:(n - 1)
        ki = solid.conductivity_scale * sic_conductivity(temperature[i])
        kj = solid.conductivity_scale * sic_conductivity(temperature[i + 1])
        kface = 2.0 * ki * kj / (ki + kj)
        qcond = kface * area_solid * (temperature[i] - temperature[i + 1]) / dx
        du[i] -= qcond
        du[i + 1] += qcond
    end

    absorbed = parameters.optics.absorbed_fraction *
               max(0.0, Float64(operating.irradiance(time))) * frontal_area(geometry)
    for i in 1:n
        side_loss = losses.side_conductance_per_length * dx * (temperature[i] - ambient)
        du[i] += absorbed * weights[i] - work.qgas[i] - side_loss
    end

    du[1] -= losses.front_convection * frontal_area(geometry) * (temperature[1] - ambient) +
             losses.front_emissivity * SIGMA * frontal_area(geometry) *
             (temperature[1]^4 - ambient^4)
    du[end] -= losses.rear_conductance * (temperature[end] - ambient) +
               losses.rear_emissivity * SIGMA * frontal_area(geometry) *
               (temperature[end]^4 - ambient^4)

    for i in 1:n
        capacity = solid.density * solid.heat_capacity_scale *
                   sic_heat_capacity(temperature[i]) * cell_volume
        du[i] /= capacity
    end
    return nothing
end

function stable_step(temperature, parameters::ModelParameters, dx, maximum_step)
    solid = parameters.solid
    alpha_max = maximum(solid.conductivity_scale * sic_conductivity(T) /
                        (solid.density * solid.heat_capacity_scale * sic_heat_capacity(T))
                        for T in temperature)
    return min(Float64(maximum_step), 0.35 * dx^2 / max(alpha_max, eps(Float64)))
end

function initial_profile(initial_temperature, z)
    if initial_temperature isa Real
        return fill(Float64(initial_temperature), length(z))
    elseif initial_temperature isa AbstractVector
        length(initial_temperature) == length(z) ||
            throw(DimensionMismatch("initial profile length must equal model nodes"))
        return Float64.(initial_temperature)
    end
    return Float64[initial_temperature(x) for x in z]
end

function simulate(parameters::ModelParameters, operating::OperatingCondition,
                  save_times::AbstractVector; initial_temperature=295.6,
                  maximum_step=2.0)
    validate(parameters)
    times = Float64.(save_times)
    isempty(times) && throw(ArgumentError("save_times cannot be empty"))
    issorted(times) || throw(ArgumentError("save_times must be sorted"))
    (length(times) == 1 || all(diff(times) .> 0.0)) ||
        throw(ArgumentError("save_times must be strictly increasing"))

    n = parameters.nodes
    dx = parameters.geometry.length / n
    z = collect(range(dx / 2, step=dx, length=n))
    zg = collect(range(0.0, parameters.geometry.length, length=n + 1))
    weights = optical_weights(parameters, dx)
    temperature = initial_profile(initial_temperature, z)
    solid_history = Matrix{Float64}(undef, n, length(times))
    gas_history = Matrix{Float64}(undef, n + 1, length(times))
    h_history = Matrix{Float64}(undef, n, length(times))
    work = Workspace(n)

    current_time = times[1]
    solid_history[:, 1] .= temperature
    gas_profile!(work.gas, work.qgas, work.h, temperature, current_time,
                 parameters, operating, z, dx)
    gas_history[:, 1] .= work.gas
    h_history[:, 1] .= work.h

    for output_index in 2:length(times)
        target_time = times[output_index]
        while current_time < target_time
            dt = min(stable_step(temperature, parameters, dx, maximum_step),
                     target_time - current_time)
            rhs!(work.k1, temperature, current_time, parameters, operating,
                 z, dx, weights, work)
            @. work.trial = temperature + 0.5 * dt * work.k1
            rhs!(work.k2, work.trial, current_time + 0.5dt, parameters, operating,
                 z, dx, weights, work)
            @. work.trial = temperature + 0.5 * dt * work.k2
            rhs!(work.k3, work.trial, current_time + 0.5dt, parameters, operating,
                 z, dx, weights, work)
            @. work.trial = temperature + dt * work.k3
            rhs!(work.k4, work.trial, current_time + dt, parameters, operating,
                 z, dx, weights, work)
            @. temperature += dt * (work.k1 + 2.0work.k2 + 2.0work.k3 + work.k4) / 6.0
            current_time += dt
        end
        solid_history[:, output_index] .= temperature
        gas_profile!(work.gas, work.qgas, work.h, temperature, current_time,
                     parameters, operating, z, dx)
        gas_history[:, output_index] .= work.gas
        h_history[:, output_index] .= work.h
    end
    return SimulationResult(times, z, zg, solid_history, gas_history, h_history)
end

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

function energy_rates(temperature::AbstractVector, time, parameters::ModelParameters,
                      operating::OperatingCondition)
    n = parameters.nodes
    length(temperature) == n || throw(DimensionMismatch("temperature length must equal model nodes"))
    geometry = parameters.geometry
    dx = geometry.length / n
    z = collect(range(dx / 2, step=dx, length=n))
    weights = optical_weights(parameters, dx)
    work = Workspace(n)
    derivative = zeros(n)
    rhs!(derivative, temperature, time, parameters, operating, z, dx, weights, work)
    ambient = Float64(operating.ambient_temperature(time))
    absorbed = parameters.optics.absorbed_fraction *
               max(0.0, Float64(operating.irradiance(time))) * frontal_area(geometry)
    gas = sum(work.qgas)
    side = parameters.losses.side_conductance_per_length * dx * sum(temperature .- ambient)
    front = parameters.losses.front_convection * frontal_area(geometry) * (temperature[1] - ambient) +
            parameters.losses.front_emissivity * SIGMA * frontal_area(geometry) *
            (temperature[1]^4 - ambient^4)
    rear = parameters.losses.rear_conductance * (temperature[end] - ambient) +
           parameters.losses.rear_emissivity * SIGMA * frontal_area(geometry) *
           (temperature[end]^4 - ambient^4)
    volume = solid_area(geometry) * dx
    storage = sum(parameters.solid.density * parameters.solid.heat_capacity_scale *
                  sic_heat_capacity(temperature[i]) * volume * derivative[i] for i in 1:n)
    return (absorbed=absorbed, gas=gas, side=side, front=front, rear=rear,
            storage=storage, residual=absorbed - gas - side - front - rear - storage)
end

end # module Receiver1DV2

import .Receiver1DV2

begin # execution controls
    const RUN_OPTIMIZATION = lowercase(get(ENV, "RECEIVER1D_RUN_OPTIMIZATION", "false")) == "true"
    const RUN_POSTPROCESS = lowercase(get(ENV, "RECEIVER1D_RUN_POSTPROCESS", "false")) == "true"
    const RUN_PLOTS = lowercase(get(ENV, "RECEIVER1D_RUN_PLOTS", "false")) == "true"
    const EXPORT_METRICS = lowercase(get(ENV, "RECEIVER1D_EXPORT_METRICS", "false")) == "true"
end

# ============================================================================
# SECTION A: GEOMETRY AND FIXED PARAMETERS
# ============================================================================
begin # fixed parameters
    Tamb = 22.448 + 273.15
    w_t = 19.0e-3
    A_frt = w_t^2
    L = 137e-3
    w_chnl = 1.5e-3
    t_chnl = 0.4e-3
    n_chnl = 10 * 10
    A_chnl_frt_all = n_chnl * w_chnl^2
    A_frt_solid = A_frt - A_chnl_frt_all
    Dh = w_chnl
    Pi_exchange = 4.0 * w_chnl * n_chnl
    Vs = A_frt_solid * L
    sigma_sb = 5.670374419e-8
    rho_s = 3200.0
    N_axial = 25

    geometry_v2 = Receiver1DV2.Geometry(
        length=L,
        receiver_width=w_t,
        channel_width=w_chnl,
        wall_thickness=t_chnl,
        channel_count=n_chnl,
    )

    # Dictionary keys intentionally mirror the legacy import/model interface.
    const Io = :Io
    const qlpm = :qlpm
    const Tinit = :Tinit
    const T_in = :T_in
    const T_amb = :T_amb
end;

# ============================================================================
# SECTION B: FITTED PARAMETER ORDER
# ============================================================================
# Cooling parameters p[1:7]
#   p[1] = k_scale       effective SiC axial-conductivity multiplier
#   p[2] = cp_scale      effective solid heat-capacity multiplier
#   p[3] = U_side_prime  side-loss conductance per length (W/m/K)
#   p[4] = U_rear        rear/adaptor loss conductance (W/K)
#   p[5] = C_Nu          effective Nusselt correlation coefficient
#   p[6] = m_Re          Reynolds-number exponent
#   p[7] = tau_T3        outlet thermocouple/downstream lag (s)
# Heating parameters p[8:11]
#   p[8] = eta_abs       absorbed solar fraction
#   p[9] = beta_opt      optical extinction coefficient (1/m)
#   p[10] = h_front      front convection coefficient (W/m2/K)
#   p[11] = C_h_refine   multiplicative heating refinement of C_Nu
# ============================================================================
begin # parameter definitions
    p_cool_init = [1.0, 1.0, 0.35, 0.10, 6.1e-4, 1.3835, 20.0]
    lb_cool = [0.1, 0.3, 0.02, 0.005, 1e-5, 0.0, 0.5]
    ub_cool = [5.0, 3.0, 5.0, 5.0, 0.1, 3.0, 500.0]

    p_heat_init = [0.85, 50.0, 10.0, 1.0]
    lb_heat = [0.10, 1.0, 0.1, 0.25]
    ub_heat = [0.99, 2000.0, 100.0, 4.0]

    assemble_parameters(p_cool, p_heat) = vcat(Float64.(p_cool), Float64.(p_heat))
end

# ============================================================================
# SECTION C: 1D MODEL CONSTRUCTION
# ============================================================================
function build_model_v2(p; nodes=N_axial)
    length(p) == 11 || throw(ArgumentError("1D_v2 expects 11 parameters; received $(length(p))"))
    all(isfinite, p) || throw(ArgumentError("all model parameters must be finite"))

    solid = Receiver1DV2.SolidProperties(
        density=rho_s,
        conductivity_scale=p[1],
        heat_capacity_scale=p[2],
    )
    htc = Receiver1DV2.HeatTransferParameters(
        coefficient=p[5] * p[11],
        reynolds_exponent=p[6],
        prandtl_exponent=0.0,
        temperature_exponent=-0.4883,
        reference_temperature=600.0,
        development_coefficient=0.0,
        development_exponent=0.45,
        minimum_nusselt=1e-5,
    )
    losses = Receiver1DV2.LossParameters(
        front_convection=p[10],
        front_emissivity=0.85,
        side_conductance_per_length=p[3],
        rear_conductance=p[4],
        rear_emissivity=0.0,
    )
    optics = Receiver1DV2.OpticalParameters(
        absorbed_fraction=p[8],
        extinction_coefficient=p[9],
        front_deposition_fraction=0.0,
    )
    return Receiver1DV2.ModelParameters(
        geometry=geometry_v2,
        solid=solid,
        heat_transfer=htc,
        losses=losses,
        optics=optics,
        nodes=nodes,
    )
end

# ============================================================================
# SECTION D: EXPERIMENTAL DATA IMPORT
# ============================================================================
begin # external import, analogous to import_exp_0D.jl
    include("import_exp_1D_v2.jl")
    sim_key_heat = copy(IDs)
    sim_key_cool = copy(IDs_cooling)
end

# ============================================================================
# SECTION E: MODEL SOLVER AND LEGACY-STYLE REMAKE INTERFACE
# ============================================================================
begin # solver helpers
    const SENSOR_POSITIONS_V2 = Dict("_T8" => 0.011, "_T9" => 0.058, "_T10" => 0.107)
    const MODEL_OUTPUT_COLUMNS = (:T8, :T9, :T10, :Tf, :h_mean)

    function initial_solid_profile(measurement_df, simulation_id)
        zdata = (0.011, 0.058, 0.107)
        tdata = (
            observation(measurement_df, simulation_id, "_T8")[1],
            observation(measurement_df, simulation_id, "_T9")[1],
            observation(measurement_df, simulation_id, "_T10")[1],
        )
        return function (z)
            z <= zdata[1] && return tdata[1]
            z >= zdata[3] && return tdata[3]
            i = z <= zdata[2] ? 1 : 2
            fraction = (z - zdata[i]) / (zdata[i + 1] - zdata[i])
            return tdata[i] + fraction * (tdata[i + 1] - tdata[i])
        end
    end

    function apply_T3_lag(target, measured_initial, time, tau)
        tau <= 1e-12 && return copy(target)
        filtered = similar(target)
        filtered[1] = measured_initial
        for i in 2:length(time)
            gain = -expm1(-(time[i] - time[i - 1]) / tau)
            filtered[i] = filtered[i - 1] + gain * (target[i] - filtered[i - 1])
        end
        return filtered
    end

    function extract_model_outputs(result, pguess; T3_initial=result.gas_temperature[end, 1])
        T8 = Receiver1DV2.solid_at(result, SENSOR_POSITIONS_V2["_T8"])
        T9 = Receiver1DV2.solid_at(result, SENSOR_POSITIONS_V2["_T9"])
        T10 = Receiver1DV2.solid_at(result, SENSOR_POSITIONS_V2["_T10"])
        Tf_true = vec(result.gas_temperature[end, :])
        Tf = apply_T3_lag(Tf_true, T3_initial, result.time, pguess[7])
        h_mean = vec(mean(result.local_heat_transfer_coefficient, dims=1))
        return hcat(T8, T9, T10, Tf, h_mean)
    end

    function solve_model_v2(pguess, operating_condition, tvalues;
                            initial_temperature=operating_condition.ambient_temperature(first(tvalues)),
                            T3_initial=operating_condition.inlet_temperature(first(tvalues)),
                            nodes=N_axial, maximum_step=2.0)
        model = build_model_v2(pguess; nodes=nodes)
        result = Receiver1DV2.simulate(model, operating_condition, tvalues;
                          initial_temperature=initial_temperature,
                          maximum_step=maximum_step)
        outputs = extract_model_outputs(result, pguess; T3_initial=T3_initial)
        return outputs, result
    end

    # Same constant-condition call shape used by 0D_v4.jl and 1D_v1.exp.All.jl.
    function solve_model(pguess, Io_val, qlpm_val, Tinit_val, tvalues;
                         Tin_val=Tamb, Tamb_val=Tamb, initial_profile=Tinit_val,
                         T3_initial=Tinit_val, tolr=1e-7, nodes=N_axial)
        operating_condition = Receiver1DV2.OperatingCondition(
            irradiance=Io_val,
            flow_lpm=qlpm_val,
            inlet_temperature=Tin_val,
            ambient_temperature=Tamb_val,
        )
        outputs, _ = solve_model_v2(pguess, operating_condition, tvalues;
                                    initial_temperature=initial_profile,
                                    T3_initial=T3_initial, nodes=nodes)
        return outputs
    end

    function NLmodeloptim(tvalues, rmp, pguess, tolr=1e-7;
                          initial_profile=rmp[Tinit], T3_initial=rmp[Tinit],
                          nodes=N_axial)
        return solve_model(
            pguess, rmp[Io], rmp[qlpm], rmp[Tinit], tvalues;
            Tin_val=get(rmp, T_in, Tamb), Tamb_val=get(rmp, T_amb, Tamb),
            initial_profile=initial_profile, T3_initial=T3_initial,
            tolr=tolr, nodes=nodes,
        )
    end

    function remakeAysha(pguess, cond_k, time_opt, Tinit_val; tolr=1e-7,
                         initial_profile=Tinit_val, T3_initial=Tinit_val,
                         nodes=N_axial)
        rmp = copy(cond_k)
        rmp[Tinit] = Tinit_val
        return NLmodeloptim(time_opt, rmp, pguess, tolr;
                            initial_profile=initial_profile,
                            T3_initial=T3_initial, nodes=nodes)
    end

    function solve_case_v2(pguess, simulation_id; is_cooling=false,
                           nodes=N_axial, maximum_step=2.0)
        measurement_df = is_cooling ? measurements_cooling : measurements
        conditions = is_cooling ? simulation_conditions_cooling : simulation_conditions
        cond_k = conditions[simulation_id]
        time = observation_time(measurement_df, simulation_id)
        flow = observation(measurement_df, simulation_id, "_flow")
        inlet = observation(measurement_df, simulation_id, "_Tin")
        ambient = observation(measurement_df, simulation_id, "_Tamb")
        irradiance = fill(cond_k[Io], length(time))
        operating_condition = Receiver1DV2.OperatingCondition(
            irradiance=Receiver1DV2.linear_history(time, irradiance),
            flow_lpm=Receiver1DV2.linear_history(time, flow),
            inlet_temperature=Receiver1DV2.linear_history(time, inlet),
            ambient_temperature=Receiver1DV2.linear_history(time, ambient),
        )
        profile = initial_solid_profile(measurement_df, simulation_id)
        T3_initial = observation(measurement_df, simulation_id, "_Tf")[1]
        return solve_model_v2(pguess, operating_condition, time;
                              initial_temperature=profile, T3_initial=T3_initial,
                              nodes=nodes, maximum_step=maximum_step)
    end
end

# ============================================================================
# SECTION F: TWO-STAGE OPTIMIZATION
# ============================================================================
begin # objective functions
    function normalized_mse(model_values, experimental_values)
        scale = max(maximum(experimental_values) - minimum(experimental_values), 20.0)
        return mean(abs2, (model_values .- experimental_values) ./ scale)
    end

    function loss_case_v2(p_full, simulation_id; is_cooling=false, nodes=15)
        measurement_df = is_cooling ? measurements_cooling : measurements
        outputs, _ = solve_case_v2(p_full, simulation_id;
                                   is_cooling=is_cooling, nodes=nodes,
                                   maximum_step=3.0)
        error = 0.0
        for (column, obs_id) in ((1, "_T8"), (2, "_T9"), (3, "_T10"))
            error += normalized_mse(outputs[:, column],
                                    observation(measurement_df, simulation_id, obs_id))
        end
        error += 0.5 * normalized_mse(outputs[:, 4],
                                      observation(measurement_df, simulation_id, "_Tf"))
        return error / 3.5
    end

    function loss_cooling(p_cool, keys=sim_key_cool; nodes=15)
        p_full = assemble_parameters(p_cool, p_heat_init)
        return mean(loss_case_v2(p_full, sm; is_cooling=true, nodes=nodes) for sm in keys)
    end

    function loss_heating(p_heat, keys=sim_key_heat; p_cool=p_cool_init, nodes=15)
        p_full = assemble_parameters(p_cool, p_heat)
        return mean(loss_case_v2(p_full, sm; is_cooling=false, nodes=nodes) for sm in keys)
    end

    function bounded_pattern_search(objective, initial, lower, upper;
                                    initial_step=0.20, maximum_iterations=30,
                                    tolerance=1e-3)
        x = clamp.(Float64.(initial), lower, upper)
        span = Float64.(upper .- lower)
        # Scale each search direction around its initial magnitude. Using only
        # the full bound span is ineffective for C_Nu, which is O(1e-4), while
        # beta_opt is O(10-100).
        step = initial_step .* max.(abs.(x), 0.02 .* span)
        fx = objective(x)
        evaluations = 1
        for iteration in 1:maximum_iterations
            improved = false
            for j in eachindex(x)
                for direction in (-1.0, 1.0)
                    trial = copy(x)
                    trial[j] = clamp(x[j] + direction * step[j], lower[j], upper[j])
                    ftrial = objective(trial)
                    evaluations += 1
                    if ftrial < fx
                        x, fx = trial, ftrial
                        improved = true
                    end
                end
            end
            if !improved
                step .*= 0.5
                maximum(step ./ max.(span, eps())) < tolerance &&
                    return (minimizer=x, minimum=fx, iterations=iteration,
                            evaluations=evaluations, converged=true)
            end
        end
        return (minimizer=x, minimum=fx, iterations=maximum_iterations,
                evaluations=evaluations, converged=false)
    end
end

begin # staged fit driver
    if RUN_OPTIMIZATION
        println("--- STAGE 1: 1D cooling optimization ---")
        cooling_fit_v2 = bounded_pattern_search(
            p -> loss_cooling(p, sim_key_cool),
            p_cool_init, lb_cool, ub_cool;
            maximum_iterations=30,
        )
        p_cool_opt = cooling_fit_v2.minimizer
        println("Cooling objective: ", cooling_fit_v2.minimum)

        println("--- STAGE 2: 1D heating optimization ---")
        heating_fit_v2 = bounded_pattern_search(
            p -> loss_heating(p, sim_key_heat; p_cool=p_cool_opt),
            p_heat_init, lb_heat, ub_heat;
            maximum_iterations=30,
        )
        p_heat_opt = heating_fit_v2.minimizer
        println("Heating objective: ", heating_fit_v2.minimum)
    else
        p_cool_opt = copy(p_cool_init)
        p_heat_opt = copy(p_heat_init)
    end
    pnew = assemble_parameters(p_cool_opt, p_heat_opt)
end

# ============================================================================
# SECTION G: POST-PROCESSING - STEADY STATE, TRANSIENTS, AND METRICS
# ============================================================================
begin # post-processing functions
    function build_steady_results_v2(params=pnew; keys=sim_key_heat, nodes=N_axial)
        results = NamedTuple[]
        for sm in keys
            outputs, _ = solve_case_v2(params, sm; nodes=nodes)
            cond = simulation_conditions[sm]
            push!(results, (
                simulation_id=sm, flow_lpm=cond[qlpm], irradiance=cond[Io],
                T8_model=outputs[end, 1], T8_experiment=observation(measurements, sm, "_T8")[end],
                T9_model=outputs[end, 2], T9_experiment=observation(measurements, sm, "_T9")[end],
                T10_model=outputs[end, 3], T10_experiment=observation(measurements, sm, "_T10")[end],
                T3_model=outputs[end, 4], T3_experiment=observation(measurements, sm, "_Tf")[end],
                h_effective=outputs[end, 5],
            ))
        end
        return results
    end

    function transient_case_v2(simulation_id, params=pnew; is_cooling=false,
                               nodes=N_axial)
        measurement_df = is_cooling ? measurements_cooling : measurements
        outputs, result = solve_case_v2(params, simulation_id;
                                        is_cooling=is_cooling, nodes=nodes)
        return (
            time=copy(result.time),
            T8_model=outputs[:, 1], T8_experiment=observation(measurement_df, simulation_id, "_T8"),
            T9_model=outputs[:, 2], T9_experiment=observation(measurement_df, simulation_id, "_T9"),
            T10_model=outputs[:, 3], T10_experiment=observation(measurement_df, simulation_id, "_T10"),
            T3_model=outputs[:, 4], T3_experiment=observation(measurement_df, simulation_id, "_Tf"),
            h_effective=outputs[:, 5]
        )
    end

    function plot_case_v2(simulation_id, params=pnew; is_cooling=false,
                          nodes=N_axial)
        @eval using StatsPlots
        data = transient_case_v2(simulation_id, params;
                                 is_cooling=is_cooling, nodes=nodes)
        plot_object = plot(title="1D_v2 transient: $simulation_id",
                           xlabel="Time (s)", ylabel="Temperature (K)",
                           legend=:outerright)
        for sensor in (:T8, :T9, :T10, :T3)
            plot!(plot_object, data.time, getproperty(data, Symbol(sensor, "_model")),
                  label="$(sensor) model", lw=2)
            scatter!(plot_object, data.time, getproperty(data, Symbol(sensor, "_experiment")),
                     label="$(sensor) experiment", ms=2, markerstrokewidth=0)
        end
        return plot_object
    end

    function get_t90(time, temperature)
        target = temperature[1] + 0.9 * (temperature[end] - temperature[1])
        index = temperature[end] >= temperature[1] ?
                findfirst(>=(target), temperature) : findfirst(<=(target), temperature)
        return isnothing(index) ? time[end] : time[index]
    end

    function compute_metrics_v2(params=pnew; heating_keys=sim_key_heat,
                                cooling_keys=sim_key_cool, nodes=N_axial)
        metrics = NamedTuple[]
        for (keys, cooling) in ((heating_keys, false), (cooling_keys, true))
            measurement_df = cooling ? measurements_cooling : measurements
            for sm in keys
                outputs, result = solve_case_v2(params, sm;
                                                is_cooling=cooling, nodes=nodes)
                for (sensor, column, obs_id) in (
                    (:T8, 1, "_T8"), (:T9, 2, "_T9"),
                    (:T10, 3, "_T10"), (:T3, 4, "_Tf"),
                )
                    experiment = observation(measurement_df, sm, obs_id)
                    model = outputs[:, column]
                    push!(metrics, (
                        simulation_id=sm,
                        phase=cooling ? :cooling : :heating,
                        sensor=sensor,
                        rmse_K=sqrt(mean(abs2, model .- experiment)),
                        steady_error_K=model[end] - experiment[end],
                        t90_error_s=get_t90(result.time, model) - get_t90(result.time, experiment),
                    ))
                end
            end
        end
        return metrics
    end
end

begin # optional top-level outputs
    if RUN_POSTPROCESS
        steady_results_v2 = build_steady_results_v2(pnew)
        metrics_v2 = compute_metrics_v2(pnew)
        display(steady_results_v2)
        display(metrics_v2)
    end

    if RUN_PLOTS
        display(plot_case_v2(first(sim_key_heat), pnew))
        display(plot_case_v2(first(sim_key_cool), pnew; is_cooling=true))
    end

    if EXPORT_METRICS
        output_path = joinpath(@__DIR__, "summaries", "analysis_results_1D_v2.csv")
        metrics_to_write = compute_metrics_v2(pnew)
        open(output_path, "w") do io
            println(io, "simulation_id,phase,sensor,rmse_K,steady_error_K,t90_error_s")
            for row in metrics_to_write
                println(io, join((row.simulation_id, row.phase, row.sensor,
                                  row.rmse_K, row.steady_error_K, row.t90_error_s), ','))
            end
        end
        println("Saved 1D_v2 metrics to $output_path")
    end
end
