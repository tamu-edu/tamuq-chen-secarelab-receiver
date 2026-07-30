# ============================================================================
# run_2D_v20.jl - v20 experiment loader and parameter builder
# ============================================================================

include("run_2D_v19.jl")
include("2D_v20.jl")

const V20 = Receiver2D_v20
const OUTPUT_DIR_2D_v20 = joinpath(@__DIR__, "summaries", "2D_v20")
const HEATING_IDS_2D_v20 = HEATING_IDS_2D_v19
const COOLING_IDS_2D_v20 = COOLING_IDS_2D_v19
const TRAIN_HEATING_2D_v20 = TRAIN_HEATING_2D_v19
const VALID_HEATING_2D_v20 = VALID_HEATING_2D_v19
const SENSOR_NAMES_2D_v20 = SENSOR_NAMES_2D_v19

function _convert_v19_for_2D_v20(source)
    s18 = source.base
    s17 = s18.base
    s15 = s17.base
    s14 = s15.base
    s12 = s14.base
    s11 = s12.base
    geometry = _clone_for_2D_v19(V20.Geometry2D, s11.geometry)
    solid = _clone_for_2D_v19(V20.SolidProperties2D, s11.solid)
    heat = _clone_for_2D_v19(
        V20.HeatTransferParameters2D, s11.heat_transfer,
    )
    losses = _clone_for_2D_v19(
        V20.LossParameters2D, s11.losses,
    )
    optics = _clone_for_2D_v19(
        V20.OpticalParameters2D, s11.optics,
    )
    hydraulics = _clone_for_2D_v19(
        V20.HydraulicParameters2D, s11.hydraulics,
    )
    d11 = V20.V11.ModelParameters2D(
        geometry, solid, heat, losses, optics, hydraulics,
    )
    assembly = _clone_for_2D_v19(
        V20.AssemblyParameters2D, s12.assembly,
    )
    observation_parameters = _clone_for_2D_v19(
        V20.ObservationParameters2D, s12.observation,
    )
    d12 = V20.V12.ModelParameters2D(
        base=d11, assembly=assembly,
        observation=observation_parameters,
    )
    network = _clone_for_2D_v19(
        V20.NetworkParameters2D, s14.network,
    )
    d14 = V20.V14.ModelParameters2D(
        base=d12, network=network,
    )
    skin = _clone_for_2D_v19(V20.SkinParameters2D, s15.skin)
    d15 = V20.V15.ModelParameters2D(base=d14, skin=skin)
    flange = _clone_for_2D_v19(
        V20.CasingFlangeParameters2D, s17.casing_flange,
    )
    d17 = V20.V17.ModelParameters2D(
        base=d15, casing_flange=flange,
    )
    source_parameters = _clone_for_2D_v19(
        V20.SourceParameters2D, s18.source,
    )
    d18 = V20.V18.ModelParameters2D(
        base=d17, source=source_parameters,
    )
    exchange = _clone_for_2D_v19(
        V20.IntegratedExchangeParameters2D,
        source.integrated_exchange,
    )
    outlet = _clone_for_2D_v19(
        V20.OutletObservationParameters2D,
        source.outlet_observation,
    )
    tube_flange = _clone_for_2D_v19(
        V20.RearTubeFlangeParameters2D,
        source.rear_tube_flange,
    )
    return V20.V19.ModelParameters2D(
        base=d18, integrated_exchange=exchange,
        outlet_observation=outlet,
        rear_tube_flange=tube_flange,
    )
end

function parameters_2D_v20(;
    t3_location=:receiver_136,
    distributed_tube_flange_h_W_m2_K=0.0,
    kwargs...,
)
    base = parameters_2D_v19(
        ; rear_tube_flange_contact_h_W_m2_K=
            distributed_tube_flange_h_W_m2_K,
        kwargs...,
    )
    return V20.ModelParameters2D(
        base=_convert_v19_for_2D_v20(base),
        t3_location=V20.T3LocationParameters2D(
            location=Symbol(t3_location),
        ),
    )
end

case_inputs_2D_v20(args...; kwargs...) =
    case_inputs_2D_v19(args...; kwargs...)

function simulate_case_2D_v20(
    inputs, p::V20.ModelParameters2D;
    full_initial_data=true,
    initialization=(
        full_initial_data ? :full : :ambient
    ),
    reltol=1e-6, abstol=1e-7, dtmax=30.0,
)
    scale = inputs.nominal >= 400000.0 ? p.optics.scale_456 :
            (inputs.nominal >= 280000.0 ? p.optics.scale_304 :
             p.optics.scale_256)
    irradiance = inputs.cooling ? zeros(length(inputs.times)) :
                 fill(scale * inputs.nominal, length(inputs.times))
    op = V20.OperatingCondition2D(
        irradiance=V20.linear_history(inputs.times, irradiance),
        flow_lpm=V20.linear_history(inputs.times, inputs.flow),
        inlet_temperature=V20.linear_history(
            inputs.times, inputs.inlet,
        ),
        ambient_temperature=V20.linear_history(
            inputs.times, inputs.ambient,
        ),
    )
    initial = if inputs.cooling && initialization !== :ambient
        t0_full = Dict(
            sensor => Float64(observation(
                inputs.data, inputs.id, "_$(sensor)",
            )[1])
            for sensor in SENSOR_NAMES_2D_v20
        )
        t0 = if initialization === :full
            t0_full
        elseif initialization === :side_T2_only
            # T3/T9/T10-free plant initialization. The unobserved core is
            # initialized from the measured perimeter profile, and the rear
            # tube from the measured back-side temperature.
            Dict(
                :T8 => t0_full[:T8],
                :T12 => t0_full[:T12],
                :T11 => t0_full[:T11],
                :T9 => t0_full[:T12],
                :T10 => t0_full[:T11],
                :T3 => t0_full[:T11],
                :T2 => t0_full[:T2],
            )
        else
            error(
                "unknown v20 initialization policy: $initialization",
            )
        end
        V20.build_initial_state_2D(
            V20.base_grid2D(p),
            p.base.base.base.base.base.base.base,
            t0, inputs.ambient[1],
        )
    else
        inputs.ambient[1]
    end
    result = V20.simulate2D(
        p, op, inputs.times; initial_temperature=initial,
        reltol=reltol, abstol=abstol, dtmax=dtmax,
    )
    initial_T3 = inputs.cooling && initialization !== :ambient ?
        Float64(observation(
            inputs.data, inputs.id, "_T3",
        )[1]) : nothing
    predictions = V20.sensor_predictions2D(
        result; initial_T3=initial_T3,
    )
    model = hcat((
        getproperty(predictions, sensor)
        for sensor in SENSOR_NAMES_2D_v20
    )...)
    return (
        inputs=inputs, result=result, predictions=predictions,
        model=model, observed=inputs.observed,
    )
end

function _write_namedtuples_csv_2D_v20(path, rows)
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
