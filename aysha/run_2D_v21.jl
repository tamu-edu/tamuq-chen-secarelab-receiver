# ============================================================================
# run_2D_v21.jl - v21 experiment loader and isolated hot-front builder
# ============================================================================

using Statistics

include("run_2D_v14.jl")
include("2D_v21.jl")

const V21 = Receiver2D_v21
const OUTPUT_DIR_2D_v21 = joinpath(@__DIR__, "summaries", "2D_v21")
const HEATING_IDS_2D_v21 = HEATING_IDS_2D_v14
const COOLING_IDS_2D_v21 = COOLING_IDS_2D_v14
const TRAIN_HEATING_2D_v21 = TRAIN_HEATING_2D_v14
const VALID_HEATING_2D_v21 = VALID_HEATING_2D_v14
const SENSOR_NAMES_2D_v21 = SENSOR_NAMES_2D_v14

function _clone_with_2D_v21(destination, source; kwargs...)
    original = NamedTuple(
        name => getfield(source, name)
        for name in fieldnames(typeof(source))
    )
    return destination(; merge(original, (; kwargs...))...)
end

function parameters_2D_v21(;
    effective_front_area_cm2=3.61,
    power_scales=(1.19, 1.54, 0.94),
    casing_flange_conductance_W_K=20.0,
    felt_conductivity_scale=1.0,
    felt_heat_capacity_scale=1.0,
    nominal_mesh=true,
    screen_mesh=false,
)
    source = inherited_parameters_2D_v14(
        nominal_mesh=nominal_mesh,
        screen_mesh=screen_mesh,
    )
    source = rebuild_network_2D_v14(
        source; lateral_scale=1.0, edge_felt_scale=3.0,
    )
    source = rebuild_heating_2D_v14(
        source; beta_opt=100.0, beam_sigma=100e-3,
        heat_scale=1.80, power_scales=power_scales,
    )
    geometry = _clone_with_2D_v21(
        V21.V20.Geometry2D, source.geometry,
    )
    solid = _clone_with_2D_v21(
        V21.V20.SolidProperties2D, source.solid,
    )
    heat = _clone_with_2D_v21(
        V21.V20.HeatTransferParameters2D, source.heat_transfer,
    )
    losses = _clone_with_2D_v21(
        V21.V20.LossParameters2D, source.losses,
    )
    optics = _clone_with_2D_v21(
        V21.V20.OpticalParameters2D, source.optics,
    )
    hydraulics = _clone_with_2D_v21(
        V21.V20.HydraulicParameters2D, source.hydraulics,
    )
    base11 = V21.V20.V11.ModelParameters2D(
        geometry, solid, heat, losses, optics, hydraulics,
    )
    assembly = _clone_with_2D_v21(
        V21.V20.AssemblyParameters2D, source.base.assembly;
        felt_conductivity_scale=felt_conductivity_scale,
        felt_heat_capacity_scale=felt_heat_capacity_scale,
    )
    observation = _clone_with_2D_v21(
        V21.V20.ObservationParameters2D, source.base.observation,
    )
    base12 = V21.V20.V12.ModelParameters2D(
        base=base11, assembly=assembly,
        observation=observation,
    )
    network = _clone_with_2D_v21(
        V21.V20.NetworkParameters2D, source.network,
    )
    base14 = V21.V20.V14.ModelParameters2D(
        base=base12, network=network,
    )
    skin = V21.V20.SkinParameters2D(
        thickness=0.05e-3,
        channel_conductance_scale=1.0,
        felt_contact_scale=0.30,
    )
    base16 = V21.V20.V16.ModelParameters2D(
        base=base14, skin=skin,
    )
    base = V21.V20.V17.ModelParameters2D(
        base=base16,
        casing_flange=V21.V20.CasingFlangeParameters2D(
            conductance_W_K=
                casing_flange_conductance_W_K,
        ),
    )
    return V21.V20.ModelParameters2D(
        base=base,
        hot_front=V21.V20.HotFrontParameters2D(
            effective_area_m2=1e-4 *
                Float64(effective_front_area_cm2),
        ),
    )
end

case_inputs_2D_v21(args...; kwargs...) =
    case_inputs_2D_v14(args...; kwargs...)

function simulate_case_2D_v21(
    inputs, p::V21.V20.ModelParameters2D;
    full_initial_data=true,
    reltol=1e-6, abstol=1e-7, dtmax=30.0,
)
    scale = inputs.nominal >= 400000.0 ? p.optics.scale_456 :
            (inputs.nominal >= 280000.0 ? p.optics.scale_304 :
             p.optics.scale_256)
    irradiance = inputs.cooling ? zeros(length(inputs.times)) :
                 fill(scale * inputs.nominal, length(inputs.times))
    op = V21.V20.OperatingCondition2D(
        irradiance=V21.V20.linear_history(inputs.times, irradiance),
        flow_lpm=V21.V20.linear_history(inputs.times, inputs.flow),
        inlet_temperature=V21.V20.linear_history(
            inputs.times, inputs.inlet,
        ),
        ambient_temperature=V21.V20.linear_history(
            inputs.times, inputs.ambient,
        ),
    )
    initial = if inputs.cooling && full_initial_data
        t0 = Dict(
            sensor => Float64(observation(
                inputs.data, inputs.id, "_$(sensor)",
            )[1])
            for sensor in SENSOR_NAMES_2D_v21
        )
        grid = V21.V20.V11.build_grid2D(
            p.base.base.base.base.base,
        )
        V21.V20.build_initial_state_2D(
            grid, p.base.base.base.base.base, t0,
            inputs.ambient[1],
        )
    else
        inputs.ambient[1]
    end
    result = V21.V20.simulate2D(
        p, op, inputs.times; initial_temperature=initial,
        reltol=reltol, abstol=abstol, dtmax=dtmax,
    )
    predictions = V21.sensor_predictions2D(result)
    model = hcat((
        getproperty(predictions, sensor)
        for sensor in SENSOR_NAMES_2D_v21
    )...)
    return (
        inputs=inputs, parameters=p, result=result,
        model=V21.V20.sensor_predictions2D(result),
        observed=Matrix(inputs.data[:, SENSOR_NAMES_2D_v21]),
    )
end

function _write_namedtuples_csv_2D_v21(path, rows)
    isempty(rows) && return
    mkpath(dirname(path))
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
