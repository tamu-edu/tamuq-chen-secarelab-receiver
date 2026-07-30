# ============================================================================
# run_2D_v15.jl - experiment loader and builders for the v15 skin model
# ============================================================================

using Statistics

include("run_2D_v14.jl")
include("2D_v15.jl")

const V15 = Receiver2D_v15
const OUTPUT_DIR_2D_v15 = joinpath(@__DIR__, "summaries", "2D_v15")
const HEATING_IDS_2D_v15 = HEATING_IDS_2D_v14
const COOLING_IDS_2D_v15 = COOLING_IDS_2D_v14
const TRAIN_HEATING_2D_v15 = TRAIN_HEATING_2D_v14
const VALID_HEATING_2D_v15 = VALID_HEATING_2D_v14
const SENSOR_NAMES_2D_v15 = SENSOR_NAMES_2D_v14

function _clone_kw_2D_v15(destination, source)
    values = NamedTuple(
        name => getfield(source, name)
        for name in fieldnames(typeof(source))
    )
    return destination(; values...)
end

function selected_v14_base_2D_v15(;
    nominal_mesh=true,
    screen_mesh=false,
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
    geometry = _clone_kw_2D_v15(V15.Geometry2D, source.geometry)
    solid = _clone_kw_2D_v15(V15.SolidProperties2D, source.solid)
    heat = _clone_kw_2D_v15(
        V15.HeatTransferParameters2D, source.heat_transfer,
    )
    losses = _clone_kw_2D_v15(V15.LossParameters2D, source.losses)
    optics = _clone_kw_2D_v15(V15.OpticalParameters2D, source.optics)
    hydraulics = _clone_kw_2D_v15(
        V15.HydraulicParameters2D, source.hydraulics,
    )
    base11 = V15.V11.ModelParameters2D(
        geometry, solid, heat, losses, optics, hydraulics,
    )
    assembly = _clone_kw_2D_v15(
        V15.AssemblyParameters2D, source.base.assembly,
    )
    observation = _clone_kw_2D_v15(
        V15.ObservationParameters2D, source.base.observation,
    )
    base12 = V15.V12.ModelParameters2D(
        base=base11, assembly=assembly, observation=observation,
    )
    network = _clone_kw_2D_v15(
        V15.NetworkParameters2D, source.network,
    )
    return V15.V14.ModelParameters2D(
        base=base12, network=network,
    )
end

function selected_parameters_2D_v15(;
    nominal_mesh=true,
    screen_mesh=false,
    skin_thickness=0.20e-3,
    skin_conductance_scale=0.30,
    skin_felt_contact_scale=3.00,
)
    base = selected_v14_base_2D_v15(
        nominal_mesh=nominal_mesh, screen_mesh=screen_mesh,
    )
    skin = V15.SkinParameters2D(
        thickness=skin_thickness,
        channel_conductance_scale=skin_conductance_scale,
        felt_contact_scale=skin_felt_contact_scale,
    )
    return V15.ModelParameters2D(base=base, skin=skin)
end

case_inputs_2D_v15(args...; kwargs...) =
    case_inputs_2D_v14(args...; kwargs...)

function simulate_case_2D_v15(
    inputs,
    p::V15.ModelParameters2D;
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
    op = V15.OperatingCondition2D(
        irradiance=V15.linear_history(inputs.times, irradiance),
        flow_lpm=V15.linear_history(inputs.times, inputs.flow),
        inlet_temperature=V15.linear_history(
            inputs.times, inputs.inlet,
        ),
        ambient_temperature=V15.linear_history(
            inputs.times, inputs.ambient,
        ),
    )
    initial = if inputs.cooling && full_initial_data
        t0 = Dict(
            sensor => Float64(observation(
                inputs.data, inputs.id, "_$(sensor)",
            )[1])
            for sensor in SENSOR_NAMES_2D_v15
        )
        grid = V15.V11.build_grid2D(p.base.base.base)
        V15.build_initial_state_2D(
            grid, p.base.base.base, t0, inputs.ambient[1],
        )
    else
        inputs.ambient[1]
    end
    result = V15.simulate2D(
        p, op, inputs.times; initial_temperature=initial,
        reltol=reltol, abstol=abstol, dtmax=dtmax,
    )
    predictions = V15.sensor_predictions2D(result)
    model = hcat((
        getproperty(predictions, sensor)
        for sensor in SENSOR_NAMES_2D_v15
    )...)
    return (
        inputs=inputs, result=result, predictions=predictions,
        model=model, observed=inputs.observed,
    )
end

function _write_namedtuples_csv_2D_v15(path, rows)
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
