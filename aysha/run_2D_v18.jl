# ============================================================================
# run_2D_v18.jl - v18 source/exchange screen experiment loader and builders
# ============================================================================

using Statistics

include("run_2D_v14.jl")
include("2D_v18.jl")

const V18 = Receiver2D_v18
const OUTPUT_DIR_2D_v18 = joinpath(@__DIR__, "summaries", "2D_v18")
const HEATING_IDS_2D_v18 = HEATING_IDS_2D_v14
const COOLING_IDS_2D_v18 = COOLING_IDS_2D_v14
const TRAIN_HEATING_2D_v18 = TRAIN_HEATING_2D_v14
const VALID_HEATING_2D_v18 = VALID_HEATING_2D_v14
const SENSOR_NAMES_2D_v18 = SENSOR_NAMES_2D_v14

function _clone_with_2D_v18(destination, source; kwargs...)
    original = NamedTuple(
        name => getfield(source, name)
        for name in fieldnames(typeof(source))
    )
    return destination(; merge(original, (; kwargs...))...)
end

function _geometry_for_mesh_2D_v18(source, mesh)
    mesh in (:screen, :nominal, :refined) ||
        error("unknown v18 mesh: $mesh")
    if mesh === :refined
        return _clone_with_2D_v18(
            V18.Geometry2D, source.geometry;
            nodes_r_rec=14, nodes_r_felt=12,
            nodes_r_case=4, nodes_z=96,
            rear_tube_nodes=60,
        )
    end
    if mesh === :nominal
        return _clone_with_2D_v18(
            V18.Geometry2D, source.geometry;
            nodes_r_rec=14, nodes_r_felt=6,
            nodes_r_case=2, nodes_z=48,
            rear_tube_nodes=30,
        )
    end
    return _clone_with_2D_v18(
        V18.Geometry2D, source.geometry;
        nodes_r_rec=14, nodes_r_felt=3,
        nodes_r_case=1, nodes_z=24,
        rear_tube_nodes=15,
    )
end

function parameters_2D_v18(;
    mesh=:nominal,
    source_model=:beer_lambert,
    deep_fraction=0.0,
    deep_length_m=80.0e-3,
    exchange_multiplier=1.0,
    casing_flange_conductance_W_K=20.0,
    felt_conductivity_scale=1.0,
    felt_heat_capacity_scale=1.0,
    felt_contact_scale=0.30,
    power_scales=(1.20, 1.40, 0.975),
)
    source = inherited_parameters_2D_v14(
        nominal_mesh=mesh !== :screen,
        screen_mesh=mesh === :screen,
    )
    source = rebuild_network_2D_v14(
        source; lateral_scale=1.0, edge_felt_scale=3.0,
    )
    source = rebuild_heating_2D_v14(
        source; beta_opt=100.0, beam_sigma=100e-3,
        heat_scale=1.80 * exchange_multiplier,
        power_scales=power_scales,
    )
    geometry = _geometry_for_mesh_2D_v18(source, mesh)
    solid = _clone_with_2D_v18(
        V18.SolidProperties2D, source.solid,
    )
    heat = _clone_with_2D_v18(
        V18.HeatTransferParameters2D, source.heat_transfer,
    )
    losses = _clone_with_2D_v18(
        V18.LossParameters2D, source.losses,
    )
    optics = _clone_with_2D_v18(
        V18.OpticalParameters2D, source.optics,
    )
    hydraulics = _clone_with_2D_v18(
        V18.HydraulicParameters2D, source.hydraulics,
    )
    base11 = V18.V11.ModelParameters2D(
        geometry, solid, heat, losses, optics, hydraulics,
    )
    assembly = _clone_with_2D_v18(
        V18.AssemblyParameters2D, source.base.assembly;
        felt_conductivity_scale=felt_conductivity_scale,
        felt_heat_capacity_scale=felt_heat_capacity_scale,
    )
    observation = _clone_with_2D_v18(
        V18.ObservationParameters2D, source.base.observation,
    )
    base12 = V18.V12.ModelParameters2D(
        base=base11, assembly=assembly, observation=observation,
    )
    network = _clone_with_2D_v18(
        V18.NetworkParameters2D, source.network,
    )
    base14 = V18.V14.ModelParameters2D(
        base=base12, network=network,
    )
    skin = V18.SkinParameters2D(
        thickness=0.05e-3,
        channel_conductance_scale=1.0,
        felt_contact_scale=felt_contact_scale,
    )
    base16 = V18.V16.ModelParameters2D(
        base=base14, skin=skin,
    )
    flange = V18.CasingFlangeParameters2D(
        conductance_W_K=casing_flange_conductance_W_K,
    )
    base17 = V18.V17.ModelParameters2D(
        base=base16, casing_flange=flange,
    )
    source_parameters = V18.SourceParameters2D(
        model=source_model,
        deep_fraction=deep_fraction,
        deep_length_m=deep_length_m,
    )
    return V18.ModelParameters2D(
        base=base17, source=source_parameters,
    )
end

case_inputs_2D_v18(args...; kwargs...) =
    case_inputs_2D_v14(args...; kwargs...)

function simulate_case_2D_v18(
    inputs,
    p::V18.ModelParameters2D;
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
    op = V18.OperatingCondition2D(
        irradiance=V18.linear_history(inputs.times, irradiance),
        flow_lpm=V18.linear_history(inputs.times, inputs.flow),
        inlet_temperature=V18.linear_history(
            inputs.times, inputs.inlet,
        ),
        ambient_temperature=V18.linear_history(
            inputs.times, inputs.ambient,
        ),
    )
    initial = if inputs.cooling && full_initial_data
        t0 = Dict(
            sensor => Float64(observation(
                inputs.data, inputs.id, "_$(sensor)",
            )[1])
            for sensor in SENSOR_NAMES_2D_v18
        )
        grid = V18.base_grid2D(p)
        V18.build_initial_state_2D(
            grid, p.base.base.base.base.base, t0,
            inputs.ambient[1],
        )
    else
        inputs.ambient[1]
    end
    result = V18.simulate2D(
        p, op, inputs.times; initial_temperature=initial,
        reltol=reltol, abstol=abstol, dtmax=dtmax,
    )
    predictions = V18.sensor_predictions2D(result)
    model = hcat((
        getproperty(predictions, sensor)
        for sensor in SENSOR_NAMES_2D_v18
    )...)
    return (
        inputs=inputs, result=result, predictions=predictions,
        model=model, observed=inputs.observed,
    )
end

function _write_namedtuples_csv_2D_v18(path, rows)
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
