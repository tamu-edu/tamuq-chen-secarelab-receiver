# ============================================================================
# run_2D_v19.jl - v19 experiment loader and parameter builder
# ============================================================================

using Statistics

include("run_2D_v18.jl")
include("2D_v19.jl")

const V19 = Receiver2D_v19
const OUTPUT_DIR_2D_v19 = joinpath(@__DIR__, "summaries", "2D_v19")
const HEATING_IDS_2D_v19 = HEATING_IDS_2D_v18
const COOLING_IDS_2D_v19 = COOLING_IDS_2D_v18
const TRAIN_HEATING_2D_v19 = TRAIN_HEATING_2D_v18
const VALID_HEATING_2D_v19 = VALID_HEATING_2D_v18
const SENSOR_NAMES_2D_v19 = SENSOR_NAMES_2D_v18

function _clone_for_2D_v19(destination, source; kwargs...)
    original = NamedTuple(
        name => getfield(source, name)
        for name in fieldnames(typeof(source))
    )
    return destination(; merge(original, (; kwargs...))...)
end

function _convert_base_2D_v19(
    source; terminal_tube_flange_scale=nothing,
)
    s17 = source.base
    s15 = s17.base
    s14 = s15.base
    s12 = s14.base
    s11 = s12.base
    geometry = _clone_for_2D_v19(
        V19.Geometry2D, s11.geometry,
    )
    solid = _clone_for_2D_v19(
        V19.SolidProperties2D, s11.solid,
    )
    heat = _clone_for_2D_v19(
        V19.HeatTransferParameters2D, s11.heat_transfer,
    )
    losses = terminal_tube_flange_scale === nothing ?
        _clone_for_2D_v19(
            V19.LossParameters2D, s11.losses,
        ) :
        _clone_for_2D_v19(
            V19.LossParameters2D, s11.losses;
            flange_conductance_scale=
                Float64(terminal_tube_flange_scale),
        )
    optics = _clone_for_2D_v19(
        V19.OpticalParameters2D, s11.optics,
    )
    hydraulics = _clone_for_2D_v19(
        V19.HydraulicParameters2D, s11.hydraulics,
    )
    d11 = V19.V11.ModelParameters2D(
        geometry, solid, heat, losses, optics, hydraulics,
    )
    assembly = _clone_for_2D_v19(
        V19.AssemblyParameters2D, s12.assembly,
    )
    observation = _clone_for_2D_v19(
        V19.ObservationParameters2D, s12.observation,
    )
    d12 = V19.V12.ModelParameters2D(
        base=d11, assembly=assembly, observation=observation,
    )
    network = _clone_for_2D_v19(
        V19.NetworkParameters2D, s14.network,
    )
    d14 = V19.V14.ModelParameters2D(
        base=d12, network=network,
    )
    skin = _clone_for_2D_v19(
        V19.SkinParameters2D, s15.skin,
    )
    d15 = V19.V15.ModelParameters2D(
        base=d14, skin=skin,
    )
    flange = _clone_for_2D_v19(
        V19.CasingFlangeParameters2D, s17.casing_flange,
    )
    d17 = V19.V17.ModelParameters2D(
        base=d15, casing_flange=flange,
    )
    source_parameters = _clone_for_2D_v19(
        V19.SourceParameters2D, source.source,
    )
    return V19.V18.ModelParameters2D(
        base=d17, source=source_parameters,
    )
end

function parameters_2D_v19(;
    mesh=:screen,
    source_model=:beer_lambert,
    deep_fraction=0.0,
    deep_length_m=80e-3,
    nu_prefactor=3.046900474047653e-4,
    reynolds_exponent=1.444917063715293,
    prandtl_exponent=0.0,
    kernel_model=:graetz,
    probe_capacity_areal_J_m2_K=1000.0,
    probe_stem_conductance_areal_W_m2_K=20.0,
    probe_emissivity=0.30,
    probe_diameter_m=1.5e-3,
    rear_tube_flange_contact_h_W_m2_K=0.0,
    rear_tube_flange_contact_start_m=28e-3,
    terminal_tube_flange_scale=nothing,
    felt_conductivity_scale=1.0,
    felt_heat_capacity_scale=1.0,
    felt_contact_scale=0.30,
    casing_flange_conductance_W_K=20.0,
    power_scales=(1.195, 1.405, 0.975),
)
    source = parameters_2D_v18(
        mesh=mesh,
        source_model=source_model,
        deep_fraction=deep_fraction,
        deep_length_m=deep_length_m,
        exchange_multiplier=1.0,
        casing_flange_conductance_W_K=
            casing_flange_conductance_W_K,
        felt_conductivity_scale=felt_conductivity_scale,
        felt_heat_capacity_scale=felt_heat_capacity_scale,
        felt_contact_scale=felt_contact_scale,
        power_scales=power_scales,
    )
    terminal_scale = terminal_tube_flange_scale === nothing ?
        (rear_tube_flange_contact_h_W_m2_K > 0.0 ? 0.0 : 1.0) :
        terminal_tube_flange_scale
    base = _convert_base_2D_v19(
        source; terminal_tube_flange_scale=terminal_scale,
    )
    exchange = V19.IntegratedExchangeParameters2D(
        nu_prefactor=nu_prefactor,
        reynolds_exponent=reynolds_exponent,
        prandtl_exponent=prandtl_exponent,
        kernel_model=kernel_model,
    )
    outlet = V19.OutletObservationParameters2D(
        capacity_areal_J_m2_K=probe_capacity_areal_J_m2_K,
        stem_conductance_areal_W_m2_K=
            probe_stem_conductance_areal_W_m2_K,
        emissivity=probe_emissivity,
        probe_diameter_m=probe_diameter_m,
        axial_position_m=3e-3,
    )
    tube_flange = V19.RearTubeFlangeParameters2D(
        contact_h_W_m2_K=rear_tube_flange_contact_h_W_m2_K,
        contact_start_m=rear_tube_flange_contact_start_m,
    )
    return V19.ModelParameters2D(
        base=base, integrated_exchange=exchange,
        outlet_observation=outlet,
        rear_tube_flange=tube_flange,
    )
end

case_inputs_2D_v19(args...; kwargs...) =
    case_inputs_2D_v18(args...; kwargs...)

function simulate_case_2D_v19(
    inputs, p::V19.ModelParameters2D;
    full_initial_data=true,
    reltol=1e-6, abstol=1e-7, dtmax=30.0,
)
    scale = inputs.nominal >= 400000.0 ? p.optics.scale_456 :
            (inputs.nominal >= 280000.0 ? p.optics.scale_304 :
             p.optics.scale_256)
    irradiance = inputs.cooling ? zeros(length(inputs.times)) :
                 fill(scale * inputs.nominal, length(inputs.times))
    op = V19.OperatingCondition2D(
        irradiance=V19.linear_history(inputs.times, irradiance),
        flow_lpm=V19.linear_history(inputs.times, inputs.flow),
        inlet_temperature=V19.linear_history(
            inputs.times, inputs.inlet,
        ),
        ambient_temperature=V19.linear_history(
            inputs.times, inputs.ambient,
        ),
    )
    initial = if inputs.cooling && full_initial_data
        t0 = Dict(
            sensor => Float64(observation(
                inputs.data, inputs.id, "_$(sensor)",
            )[1])
            for sensor in SENSOR_NAMES_2D_v19
        )
        V19.build_initial_state_2D(
            V19.base_grid2D(p),
            p.base.base.base.base.base.base,
            t0, inputs.ambient[1],
        )
    else
        inputs.ambient[1]
    end
    result = V19.simulate2D(
        p, op, inputs.times; initial_temperature=initial,
        reltol=reltol, abstol=abstol, dtmax=dtmax,
    )
    initial_T3 = inputs.cooling && full_initial_data ?
        Float64(observation(
            inputs.data, inputs.id, "_T3",
        )[1]) : nothing
    predictions = V19.sensor_predictions2D(
        result; initial_T3=initial_T3,
    )
    model = hcat((
        getproperty(predictions, sensor)
        for sensor in SENSOR_NAMES_2D_v19
    )...)
    return (
        inputs=inputs, result=result, predictions=predictions,
        model=model, observed=inputs.observed,
    )
end

function _write_namedtuples_csv_2D_v19(path, rows)
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
