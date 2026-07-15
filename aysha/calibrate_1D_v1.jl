module Receiver1DCalibration

include("1D_v1.jl")
using .Receiver1D
using Statistics

export Experiment, PatternSearchResult, load_heating_experiments,
       load_cooling_experiments, cooling_objective, heating_objective,
       fit_cooling, fit_heating, pattern_search

const SENSOR_POSITIONS = Dict(:T8 => 0.011, :T9 => 0.058, :T10 => 0.107)

const HEATING_FILES = [
    ("E67", "Data_FPT0067_231125_161757.csv", 456e3),
    ("E68", "Data_FPT0068_231126_115725.csv", 456e3),
    ("E69", "Data_FPT0069_231126_140153.csv", 456e3),
    ("E70", "Data_FPT0070_231127_090339.csv", 456e3),
    ("E71", "Data_FPT0071_231128_102707.csv", 456e3),
    ("E72", "Data_FPT0072_231129_104140.csv", 304e3),
    ("E73", "Data_FPT0073_231129_132744.csv", 304e3),
    ("E74", "Data_FPT0074_231130_123228.csv", 304e3),
    ("E75", "Data_FPT0075_231201_162138.csv", 304e3),
    ("E76", "Data_FPT0076_231203_120521.csv", 304e3),
    ("E77", "Data_FPT0077_231203_161315.csv", 256e3),
    ("E78", "Data_FPT0078_231204_132252.csv", 256e3),
    ("E79", "Data_FPT0079_231204_172244.csv", 256e3),
    ("E80", "Data_FPT0080_231205_095122.csv", 256e3),
    ("E81", "Data_FPT0081_231205_135354.csv", 256e3),
]

const COOLING_FILES = [
    ("C69", "Data_FPT0069-Cooling_231126_153148.csv", 0.0),
    ("C80", "Data_FPT0080-cooling_231205_112837.csv", 0.0),
    ("C81", "Data_FPT0081-cooling_231205_153409.csv", 0.0),
]

struct Experiment
    id::String
    phase::Symbol
    time::Vector{Float64}
    flow_lpm::Vector{Float64}
    irradiance::Vector{Float64}
    inlet_temperature::Vector{Float64}
    ambient_temperature::Vector{Float64}
    observations::Dict{Symbol,Vector{Float64}}
end

struct PatternSearchResult
    minimizer::Vector{Float64}
    minimum::Float64
    iterations::Int
    converged::Bool
end

function selected_columns(path)
    rows = Vector{Vector{Float64}}()
    open(path, "r") do io
        first_line = true
        for line in eachline(io)
            if first_line
                first_line = false
                continue
            end
            fields = split(line, ',')
            length(fields) >= 50 || continue
            # time; measured MFCs; T3; T8-T12; two room-temperature channels
            indices = (2, 7, 8, 9, 10, 37, 42, 43, 44, 45, 46, 49, 50)
            try
                push!(rows, [parse(Float64, fields[i]) for i in indices])
            catch error
                error isa ArgumentError || rethrow()
            end
        end
    end
    isempty(rows) && error("no numeric rows found in $path")
    return reduce(vcat, permutedims.(rows))
end

function load_experiment(root, spec, phase; stride=50)
    id, filename, flux = spec
    raw = selected_columns(joinpath(root, "SolarSimulator", "RAW", filename))
    indices = unique(vcat(collect(1:stride:size(raw, 1)), size(raw, 1)))
    data = raw[indices, :]
    time = data[:, 1] .- data[1, 1]
    flow = vec(sum(data[:, 2:5], dims=2))
    ambient = 0.5 .* (data[:, 12] .+ data[:, 13]) .+ 273.15
    observations = Dict(
        :T3 => data[:, 6] .+ 273.15,
        :T8 => data[:, 7] .+ 273.15,
        :T9 => data[:, 8] .+ 273.15,
        :T10 => data[:, 9] .+ 273.15,
        :T11 => data[:, 10] .+ 273.15,
        :T12 => data[:, 11] .+ 273.15,
    )
    return Experiment(id, phase, time, flow, fill(flux, length(time)),
                      copy(ambient), ambient, observations)
end

function load_heating_experiments(root=@__DIR__; stride=50)
    return [load_experiment(root, spec, :heating; stride=stride) for spec in HEATING_FILES]
end

function load_cooling_experiments(root=@__DIR__; stride=50)
    return [load_experiment(root, spec, :cooling; stride=stride) for spec in COOLING_FILES]
end

function experiment_condition(experiment::Experiment)
    return OperatingCondition(
        irradiance=linear_history(experiment.time, experiment.irradiance),
        flow_lpm=linear_history(experiment.time, experiment.flow_lpm),
        inlet_temperature=linear_history(experiment.time, experiment.inlet_temperature),
        ambient_temperature=linear_history(experiment.time, experiment.ambient_temperature),
    )
end

function initial_profile(experiment::Experiment)
    zdata = (SENSOR_POSITIONS[:T8], SENSOR_POSITIONS[:T9], SENSOR_POSITIONS[:T10])
    tdata = (experiment.observations[:T8][1], experiment.observations[:T9][1],
             experiment.observations[:T10][1])
    return function (z)
        z <= zdata[1] && return tdata[1]
        z >= zdata[end] && return tdata[end]
        i = z <= zdata[2] ? 1 : 2
        fraction = (z - zdata[i]) / (zdata[i + 1] - zdata[i])
        return tdata[i] + fraction * (tdata[i + 1] - tdata[i])
    end
end

function gas_sensor_response(target, measured_initial, time, tau)
    tau <= 1e-9 && return copy(target)
    filtered = similar(target)
    filtered[1] = measured_initial
    for i in 2:length(time)
        gain = -expm1(-(time[i] - time[i - 1]) / tau)
        filtered[i] = filtered[i - 1] + gain * (target[i] - filtered[i - 1])
    end
    return filtered
end

function normalized_error(prediction, measurement)
    scale = max(maximum(measurement) - minimum(measurement), 20.0)
    return mean(abs2, (prediction .- measurement) ./ scale)
end

function rebuild(base::ModelParameters; conductivity_scale=base.solid.conductivity_scale,
                 heat_capacity_scale=base.solid.heat_capacity_scale,
                 side_conductance=base.losses.side_conductance_per_length,
                 rear_conductance=base.losses.rear_conductance,
                 htc_coefficient=base.heat_transfer.coefficient,
                 reynolds_exponent=base.heat_transfer.reynolds_exponent,
                 absorbed_fraction=base.optics.absorbed_fraction,
                 extinction_coefficient=base.optics.extinction_coefficient,
                 front_convection=base.losses.front_convection)
    solid = SolidProperties(density=base.solid.density,
                            conductivity_scale=conductivity_scale,
                            heat_capacity_scale=heat_capacity_scale)
    oldh = base.heat_transfer
    htc = HeatTransferParameters(
        coefficient=htc_coefficient,
        reynolds_exponent=reynolds_exponent,
        prandtl_exponent=oldh.prandtl_exponent,
        temperature_exponent=oldh.temperature_exponent,
        reference_temperature=oldh.reference_temperature,
        development_coefficient=oldh.development_coefficient,
        development_exponent=oldh.development_exponent,
        minimum_nusselt=oldh.minimum_nusselt,
    )
    oldloss = base.losses
    losses = LossParameters(
        front_convection=front_convection,
        front_emissivity=oldloss.front_emissivity,
        side_conductance_per_length=side_conductance,
        rear_conductance=rear_conductance,
        rear_emissivity=oldloss.rear_emissivity,
    )
    optics = OpticalParameters(absorbed_fraction=absorbed_fraction,
                               extinction_coefficient=extinction_coefficient,
                               front_deposition_fraction=base.optics.front_deposition_fraction)
    return ModelParameters(geometry=base.geometry, solid=solid, heat_transfer=htc,
                           losses=losses, optics=optics, nodes=base.nodes)
end

function cooling_parameters(theta, base)
    model = rebuild(base,
        conductivity_scale=exp(theta[1]),
        heat_capacity_scale=exp(theta[2]),
        side_conductance=exp(theta[3]),
        rear_conductance=exp(theta[4]),
        htc_coefficient=exp(theta[5]),
        reynolds_exponent=theta[6],
    )
    return model, exp(theta[7])
end

function cooling_objective(theta, experiments, base=ModelParameters(nodes=15))
    model, gas_sensor_tau = cooling_parameters(theta, base)
    total = 0.0
    for experiment in experiments
        result = simulate(model, experiment_condition(experiment), experiment.time;
                          initial_temperature=initial_profile(experiment), maximum_step=3.0)
        case_error = 0.0
        for sensor in (:T8, :T9, :T10)
            case_error += normalized_error(solid_at(result, SENSOR_POSITIONS[sensor]),
                                           experiment.observations[sensor])
        end
        gas = gas_sensor_response(vec(result.gas_temperature[end, :]),
                                  experiment.observations[:T3][1], experiment.time,
                                  gas_sensor_tau)
        case_error += 0.5 * normalized_error(gas, experiment.observations[:T3])
        total += case_error / 3.5
    end
    # Weak physical priors prevent the three cooling runs from fitting non-unique extremes.
    prior = 0.01 * (theta[1]^2 + theta[2]^2) +
            0.005 * (log(exp(theta[3]) / 0.35)^2 + log(exp(theta[4]) / 0.10)^2) +
            0.002 * ((theta[6] - 1.38) / 0.5)^2
    return total / length(experiments) + prior
end

function heating_parameters(theta, cooling_model)
    eta = 1.0 / (1.0 + exp(-theta[1]))
    return rebuild(cooling_model,
        absorbed_fraction=eta,
        extinction_coefficient=exp(theta[2]),
        front_convection=exp(theta[3]),
        htc_coefficient=cooling_model.heat_transfer.coefficient * exp(theta[4]),
    )
end

function heating_objective(theta, experiments, cooling_model; gas_sensor_tau=20.0)
    model = heating_parameters(theta, cooling_model)
    total = 0.0
    for experiment in experiments
        result = simulate(model, experiment_condition(experiment), experiment.time;
                          initial_temperature=initial_profile(experiment), maximum_step=3.0)
        case_error = 0.0
        for sensor in (:T8, :T9, :T10)
            case_error += normalized_error(solid_at(result, SENSOR_POSITIONS[sensor]),
                                           experiment.observations[sensor])
        end
        gas = gas_sensor_response(vec(result.gas_temperature[end, :]),
                                  experiment.observations[:T3][1], experiment.time,
                                  gas_sensor_tau)
        case_error += 0.5 * normalized_error(gas, experiment.observations[:T3])
        total += case_error / 3.5
    end
    prior = 0.005 * theta[4]^2
    return total / length(experiments) + prior
end

function pattern_search(objective, initial; initial_step=fill(0.25, length(initial)),
                        lower=fill(-Inf, length(initial)), upper=fill(Inf, length(initial)),
                        maximum_iterations=60, tolerance=1e-3)
    x = clamp.(Float64.(initial), lower, upper)
    step = Float64.(initial_step)
    fx = objective(x)
    for iteration in 1:maximum_iterations
        improved = false
        for j in eachindex(x)
            for direction in (-1.0, 1.0)
                trial = copy(x)
                trial[j] = clamp(x[j] + direction * step[j], lower[j], upper[j])
                ftrial = objective(trial)
                if ftrial < fx
                    x, fx = trial, ftrial
                    improved = true
                end
            end
        end
        if !improved
            step .*= 0.5
            maximum(step) < tolerance && return PatternSearchResult(x, fx, iteration, true)
        end
    end
    return PatternSearchResult(x, fx, maximum_iterations, false)
end

function fit_cooling(experiments=load_cooling_experiments();
                     base=ModelParameters(nodes=15), maximum_iterations=30)
    initial = [0.0, 0.0, log(0.35), log(0.10),
               log(base.heat_transfer.coefficient),
               base.heat_transfer.reynolds_exponent, log(20.0)]
    lower = [log(0.1), log(0.3), log(0.02), log(0.005), log(1e-5), 0.0, log(0.5)]
    upper = [log(5.0), log(3.0), log(5.0), log(5.0), log(0.1), 3.0, log(200.0)]
    step = [0.25, 0.20, 0.30, 0.30, 0.30, 0.15, 0.30]
    result = pattern_search(x -> cooling_objective(x, experiments, base), initial;
                            initial_step=step, lower=lower, upper=upper,
                            maximum_iterations=maximum_iterations)
    model, gas_sensor_tau = cooling_parameters(result.minimizer, base)
    return (optimization=result, model=model, gas_sensor_tau=gas_sensor_tau)
end

function fit_heating(cooling_fit, experiments=load_heating_experiments();
                     maximum_iterations=30)
    model0 = cooling_fit.model
    eta0 = clamp(model0.optics.absorbed_fraction, 1e-4, 1 - 1e-4)
    initial = [log(eta0 / (1 - eta0)), log(model0.optics.extinction_coefficient),
               log(model0.losses.front_convection), 0.0]
    lower = [log(0.1 / 0.9), log(1.0), log(0.1), log(0.25)]
    upper = [log(0.99 / 0.01), log(2000.0), log(100.0), log(4.0)]
    step = [0.25, 0.35, 0.30, 0.15]
    tau = cooling_fit.gas_sensor_tau
    result = pattern_search(x -> heating_objective(x, experiments, model0;
                                                    gas_sensor_tau=tau), initial;
                            initial_step=step, lower=lower, upper=upper,
                            maximum_iterations=maximum_iterations)
    return (optimization=result, model=heating_parameters(result.minimizer, model0),
            gas_sensor_tau=tau)
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("Loaded staged 1D calibration workflow.")
    println("Call cooling = fit_cooling(); then heating = fit_heating(cooling).")
end

end # module Receiver1DCalibration
