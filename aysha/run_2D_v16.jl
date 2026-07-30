# ============================================================================
# run_2D_v16.jl - experiment loader and cooling-first parameter builders
# ============================================================================

using Statistics

include("run_2D_v14.jl")
include("2D_v16.jl")

const V16 = Receiver2D_v16
const OUTPUT_DIR_2D_v16 = joinpath(@__DIR__, "summaries", "2D_v16")
const HEATING_IDS_2D_v16 = HEATING_IDS_2D_v14
const COOLING_IDS_2D_v16 = COOLING_IDS_2D_v14
const TRAIN_HEATING_2D_v16 = TRAIN_HEATING_2D_v14
const VALID_HEATING_2D_v16 = VALID_HEATING_2D_v14
const SENSOR_NAMES_2D_v16 = SENSOR_NAMES_2D_v14

function _clone_with_2D_v16(destination, source; kwargs...)
    original = NamedTuple(
        name => getfield(source, name)
        for name in fieldnames(typeof(source))
    )
    return destination(; merge(original, (; kwargs...))...)
end

function parameters_2D_v16(;
    nominal_mesh=true,
    screen_mesh=false,
    skin_felt_contact_scale=3.0,
    felt_conductivity_scale=1.20,
    felt_heat_capacity_scale=1.50,
)
    source = inherited_parameters_2D_v14(
        nominal_mesh=nominal_mesh, screen_mesh=screen_mesh,
    )
    source = rebuild_network_2D_v14(
        source; lateral_scale=1.0, edge_felt_scale=3.0,
    )
    source = rebuild_heating_2D_v14(
        source; beta_opt=100.0, beam_sigma=100e-3,
        heat_scale=1.80,
        power_scales=(1.25, 1.2075, 0.8855),
    )
    geometry = _clone_with_2D_v16(V16.Geometry2D, source.geometry)
    solid = _clone_with_2D_v16(V16.SolidProperties2D, source.solid)
    heat = _clone_with_2D_v16(
        V16.HeatTransferParameters2D, source.heat_transfer,
    )
    losses = _clone_with_2D_v16(
        V16.LossParameters2D, source.losses,
    )
    optics = _clone_with_2D_v16(
        V16.OpticalParameters2D, source.optics,
    )
    hydraulics = _clone_with_2D_v16(
        V16.HydraulicParameters2D, source.hydraulics,
    )
    base11 = V16.V11.ModelParameters2D(
        geometry, solid, heat, losses, optics, hydraulics,
    )
    assembly = _clone_with_2D_v16(
        V16.AssemblyParameters2D, source.base.assembly;
        felt_conductivity_scale=felt_conductivity_scale,
        felt_heat_capacity_scale=felt_heat_capacity_scale,
    )
    observation = _clone_with_2D_v16(
        V16.ObservationParameters2D, source.base.observation,
    )
    base12 = V16.V12.ModelParameters2D(
        base=base11, assembly=assembly, observation=observation,
    )
    network = _clone_with_2D_v16(
        V16.NetworkParameters2D, source.network,
    )
    base14 = V16.V14.ModelParameters2D(
        base=base12, network=network,
    )
    skin = V16.SkinParameters2D(
        thickness=0.05e-3,
        channel_conductance_scale=1.0,
        felt_contact_scale=skin_felt_contact_scale,
    )
    return V16.ModelParameters2D(base=base14, skin=skin)
end

case_inputs_2D_v16(args...; kwargs...) =
    case_inputs_2D_v14(args...; kwargs...)

function simulate_case_2D_v16(
    inputs,
    p::V16.ModelParameters2D;
    full_initial_data=true,
    reltol=1e-6,
    abstol=1e-7,
    dtmax=30.0,
)
    scale = inputs.nominal >= 400000.0 ? p.optics.scale_456 :
            (inputs.nominal >= 280000.0 ? p.optics.scale_304 :
             p.optics.scale_256)
    irradiance = inputs.cooling ? zeros(length(inputs.times)) :
                 fill(scale * inputs.nominal, length(inputs.times))
    op = V16.OperatingCondition2D(
        irradiance=V16.linear_history(inputs.times, irradiance),
        flow_lpm=V16.linear_history(inputs.times, inputs.flow),
        inlet_temperature=V16.linear_history(
            inputs.times, inputs.inlet,
        ),
        ambient_temperature=V16.linear_history(
            inputs.times, inputs.ambient,
        ),
    )
    initial = if inputs.cooling && full_initial_data
        t0 = Dict(
            sensor => Float64(observation(
                inputs.data, inputs.id, "_$(sensor)",
            )[1])
            for sensor in SENSOR_NAMES_2D_v16
        )
        grid = V16.V11.build_grid2D(p.base.base.base)
        V16.build_initial_state_2D(
            grid, p.base.base.base, t0, inputs.ambient[1],
        )
    else
        inputs.ambient[1]
    end
    result = V16.simulate2D(
        p, op, inputs.times; initial_temperature=initial,
        reltol=reltol, abstol=abstol, dtmax=dtmax,
    )
    predictions = V16.sensor_predictions2D(result)
    model = hcat((
        getproperty(predictions, sensor)
        for sensor in SENSOR_NAMES_2D_v16
    )...)
    return (
        inputs=inputs, result=result, predictions=predictions,
        model=model, observed=inputs.observed,
    )
end

function _write_namedtuples_csv_2D_v16(path, rows)
    isempty(rows) && return
    names = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(String.(names), ","))
        for row in rows
            println(io, join(
                (getproperty(row, name) for name in names), ",",
            ))
        end
    end
end

