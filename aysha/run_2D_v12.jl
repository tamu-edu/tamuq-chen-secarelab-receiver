# ============================================================================
# run_2D_v12.jl - shared case loader, baseline matrix and metrics
# ============================================================================

using Statistics
using DelimitedFiles

begin
    const Io = :Io
    const qlpm = :qlpm
    const inlet_temperature = :inlet_temperature
    const ambient_temperature = :ambient_temperature
    const time_data = :time_data
    const Tinit = :Tinit
    const T_in = :T_in
    const T_amb = :T_amb
end

include("2D_v12.jl")
include("import_exp_1D_v2.jl")
using .Receiver2D_v12
const V11_v12run = Receiver2D_v12.Receiver2D_v11

const OUTPUT_DIR_2D_v12 = joinpath(@__DIR__, "summaries", "2D_v12")
const HEATING_IDS_2D_v12 = Tuple(IDs)
const COOLING_IDS_2D_v12 = ("C69", "C80", "C81")
const TRAIN_HEATING_2D_v12 = (
    "E67", "E69", "E71", "E72", "E74", "E76",
    "E77", "E79", "E81",
)
const VALID_HEATING_2D_v12 = (
    "E68", "E70", "E73", "E75", "E78", "E80",
)
const SENSOR_NAMES_2D_v12 = (
    :T8, :T12, :T11, :T9, :T10, :T3, :T2,
)

function base_parameters_2D_v12(; nominal_mesh=true, screen_mesh=false)
    g0 = Geometry2D()
    geometry = if screen_mesh
        Geometry2D(
            length=g0.length,
            receiver_width=g0.receiver_width,
            receiver_radius=g0.receiver_radius,
            insulation_outer_radius=g0.insulation_outer_radius,
            casing_outer_radius=g0.casing_outer_radius,
            channel_width=g0.channel_width,
            wall_thickness=g0.wall_thickness,
            channel_count=g0.channel_count,
            nodes_r_rec=6,
            nodes_r_felt=3,
            nodes_r_case=1,
            nodes_z=18,
            rear_tube_length=g0.rear_tube_length,
            rear_tube_nodes=10,
            rear_tube_inner_radius=g0.rear_tube_inner_radius,
            rear_tube_outer_radius=g0.rear_tube_outer_radius,
            t3_distance_from_receiver=g0.t3_distance_from_receiver,
        )
    elseif nominal_mesh
        g0
    else
        Geometry2D(
            length=g0.length,
            receiver_width=g0.receiver_width,
            receiver_radius=g0.receiver_radius,
            insulation_outer_radius=g0.insulation_outer_radius,
            casing_outer_radius=g0.casing_outer_radius,
            channel_width=g0.channel_width,
            wall_thickness=g0.wall_thickness,
            channel_count=g0.channel_count,
            nodes_r_rec=10,
            nodes_r_felt=5,
            nodes_r_case=2,
            nodes_z=30,
            rear_tube_length=g0.rear_tube_length,
            rear_tube_nodes=20,
            rear_tube_inner_radius=g0.rear_tube_inner_radius,
            rear_tube_outer_radius=g0.rear_tube_outer_radius,
            t3_distance_from_receiver=g0.t3_distance_from_receiver,
        )
    end
    solid = SolidProperties2D(
        density=2150.0,
        radial_conductivity_scale=0.05,
        axial_conductivity_scale=1.0,
        rad_extinction_coeff=50.0,
        felt_conductivity_ref=0.06, # superseded by v12 correction
        felt_density=140.0,
        felt_heat_capacity=1360.0,
        casing_conductivity=205.0,
        casing_density=2700.0,
        casing_heat_capacity=900.0,
        receiver_felt_contact_resistance=1.0e-3,
    )
    heat = HeatTransferParameters2D(
        coefficient=1.0e-3,
        reynolds_exponent=1.44,
        prandtl_exponent=1.0 / 3.0,
        minimum_nusselt=0.0,
        c_radial_flow=0.10,
        front_coefficient=0.0,
        front_reynolds_exponent=0.5,
        nusselt_model=:graetz,
        fully_developed_nusselt=3.61,
        temperature_correction_exponent=0.0,
        heat_transfer_scale=1.142,
    )
    optics = OpticalParameters2D(
        absorbed_fraction=1.0,
        extinction_coefficient=218.8,
        beam_radius_sigma=14.0e-3,
        spillage_fraction=0.10,
        front_deposition_fraction=0.20,
        scale_456=1.25,
        scale_304=1.051748,
        scale_256=0.773518,
    )
    hydraulics = HydraulicParameters2D(
        standard_pressure=101.4e3,
        standard_temperature=294.25,
        atmospheric_pressure=101.325e3,
        mass_flow_scale=1.0,
        dp1_zero_offset_mbar=-0.5428447496201336,
        hydraulic_resistance_scale=0.9701974617588522,
        minor_loss_coefficient=0.0,
        flow_distribution_model=:equal_pressure,
        groove_free_diameter=13.0e-3,
        groove_loss_coefficient=184.15800344228472,
        equal_pressure_max_iterations=16,
        equal_pressure_tolerance=1.0e-5,
        equal_pressure_relaxation=0.55,
    )
    losses = LossParameters2D(
        front_emissivity=0.85,
        casing_emissivity=0.20,
        front_loss_scale=1.0,
        casing_loss_scale=1.0,
        rear_adaptor_conductance=0.0,
        flange_conductance_scale=1.0,
    )
    base = V11_v12run.ModelParameters2D(
        geometry, solid, heat, losses, optics, hydraulics,
    )
    return ModelParameters2D(base=base)
end

function rebuild_parameters_2D_v12(
    p::ModelParameters2D;
    felt_k_scale=p.assembly.felt_conductivity_scale,
    felt_cp_scale=p.assembly.felt_heat_capacity_scale,
    receiver_adaptor_h=p.assembly.receiver_adaptor_contact_h,
    adaptor_tube_h=p.assembly.adaptor_tube_contact_h,
    adaptor_felt_scale=p.assembly.adaptor_felt_conductance_scale,
    side_tau=p.observation.side_time_constant_s,
    interior_tau=p.observation.interior_time_constant_s,
    outlet_tau=p.observation.outlet_time_constant_s,
    interior_wall_fraction=p.observation.interior_wall_fraction,
    beta_opt=p.optics.extinction_coefficient,
    heat_scale=p.heat_transfer.heat_transfer_scale,
    power_scales=(
        p.optics.scale_456, p.optics.scale_304, p.optics.scale_256,
    ),
)
    a0 = p.assembly
    assembly = AssemblyParameters2D(
        enclosure_inner_radius=a0.enclosure_inner_radius,
        enclosure_outer_radius=a0.enclosure_outer_radius,
        enclosure_internal_length=a0.enclosure_internal_length,
        backplate_thickness=a0.backplate_thickness,
        adaptor_outer_radius=a0.adaptor_outer_radius,
        adaptor_length=a0.adaptor_length,
        adaptor_receiver_overlap=a0.adaptor_receiver_overlap,
        adaptor_tube_overlap=a0.adaptor_tube_overlap,
        adaptor_density=a0.adaptor_density,
        receiver_adaptor_contact_h=receiver_adaptor_h,
        adaptor_tube_contact_h=adaptor_tube_h,
        adaptor_felt_conductance_scale=adaptor_felt_scale,
        felt_conductivity_scale=felt_k_scale,
        felt_heat_capacity_scale=felt_cp_scale,
        aluminum_density=a0.aluminum_density,
        aluminum_heat_capacity=a0.aluminum_heat_capacity,
        aluminum_conductivity=a0.aluminum_conductivity,
        housing_extension_loss_scale=a0.housing_extension_loss_scale,
    )
    observation = ObservationParameters2D(
        side_time_constant_s=side_tau,
        interior_time_constant_s=interior_tau,
        outlet_time_constant_s=outlet_tau,
        interior_wall_fraction=interior_wall_fraction,
    )
    h0 = p.heat_transfer
    heat = HeatTransferParameters2D(
        coefficient=h0.coefficient,
        reynolds_exponent=h0.reynolds_exponent,
        prandtl_exponent=h0.prandtl_exponent,
        minimum_nusselt=h0.minimum_nusselt,
        c_radial_flow=h0.c_radial_flow,
        front_coefficient=h0.front_coefficient,
        front_reynolds_exponent=h0.front_reynolds_exponent,
        nusselt_model=h0.nusselt_model,
        fully_developed_nusselt=h0.fully_developed_nusselt,
        temperature_correction_exponent=h0.temperature_correction_exponent,
        heat_transfer_scale=heat_scale,
    )
    o0 = p.optics
    optics = OpticalParameters2D(
        absorbed_fraction=o0.absorbed_fraction,
        extinction_coefficient=beta_opt,
        beam_radius_sigma=o0.beam_radius_sigma,
        spillage_fraction=o0.spillage_fraction,
        front_deposition_fraction=o0.front_deposition_fraction,
        scale_456=power_scales[1],
        scale_304=power_scales[2],
        scale_256=power_scales[3],
    )
    b0 = p.base
    base = V11_v12run.ModelParameters2D(
        b0.geometry, b0.solid, heat, b0.losses, optics, b0.hydraulics,
    )
    return ModelParameters2D(
        base=base, assembly=assembly, observation=observation,
    )
end

function _history_2D_v12(data, id, suffix)
    return Float64.(observation(data, id, suffix))
end

function case_inputs_2D_v12(id; cooling=false, max_points=nothing)
    data = cooling ? measurements_cooling : measurements
    conditions = cooling ? simulation_conditions_cooling :
                 simulation_conditions
    times = Float64.(observation_time(data, id))
    indices = if max_points === nothing || length(times) <= max_points
        collect(eachindex(times))
    else
        unique(round.(Int, range(1, length(times); length=max_points)))
    end
    times = times[indices]
    flow = _history_2D_v12(data, id, "_flow")[indices]
    inlet = _history_2D_v12(data, id, "_Tin")[indices]
    ambient = _history_2D_v12(data, id, "_Tamb")[indices]
    observed = hcat((
        _history_2D_v12(data, id, "_$(sensor)")[indices]
        for sensor in SENSOR_NAMES_2D_v12
    )...)
    nominal = Float64(conditions[id][Io])
    dp1 = _history_2D_v12(data, id, "_DP1")[indices]
    return (
        id=id, cooling=cooling, data=data, times=times, flow=flow,
        inlet=inlet, ambient=ambient, observed=observed,
        nominal=nominal, dp1=dp1,
    )
end

function simulate_case_2D_v12(
    inputs, p::ModelParameters2D;
    full_initial_data=true,
)
    scale = inputs.nominal >= 400000.0 ? p.optics.scale_456 :
            (inputs.nominal >= 280000.0 ? p.optics.scale_304 :
             p.optics.scale_256)
    irradiance = inputs.cooling ? zeros(length(inputs.times)) :
                 fill(scale * inputs.nominal, length(inputs.times))
    op = OperatingCondition2D(
        irradiance=linear_history(inputs.times, irradiance),
        flow_lpm=linear_history(inputs.times, inputs.flow),
        inlet_temperature=linear_history(inputs.times, inputs.inlet),
        ambient_temperature=linear_history(inputs.times, inputs.ambient),
    )
    initial = if inputs.cooling && full_initial_data
        # Use the complete t0 observations even when the fitting trace is
        # downsampled.
        data = inputs.data
        t0 = Dict(
            sensor => Float64(observation(
                data, inputs.id, "_$(sensor)",
            )[1])
            for sensor in SENSOR_NAMES_2D_v12
        )
        grid = V11_v12run.build_grid2D(p.base)
        build_initial_state_2D(
            grid, p.base, t0, inputs.ambient[1],
        )
    else
        inputs.ambient[1]
    end
    result = simulate2D(
        p, op, inputs.times; initial_temperature=initial,
    )
    predictions = sensor_predictions2D(result)
    model = hcat((
        getproperty(predictions, sensor)
        for sensor in SENSOR_NAMES_2D_v12
    )...)
    return (
        inputs=inputs, result=result, predictions=predictions,
        model=model, observed=inputs.observed,
    )
end

function case_metrics_2D_v12(case)
    rows = NamedTuple[]
    for (index, sensor) in enumerate(SENSOR_NAMES_2D_v12)
        model = case.model[:, index]
        observed = case.observed[:, index]
        push!(rows, (
            simulation_id=case.inputs.id,
            phase=case.inputs.cooling ? "cooling" : "heating",
            sensor=String(sensor),
            rmse_K=sqrt(mean(abs2, model .- observed)),
            steady_error_K=model[end] - observed[end],
            t90_error_s=get_t90_2D(case.inputs.times, model) -
                get_t90_2D(case.inputs.times, observed),
        ))
    end
    return rows
end

function _write_namedtuples_csv_2D_v12(path, rows)
    isempty(rows) && return
    names = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(String.(names), ","))
        for row in rows
            println(io, join((getproperty(row, name) for name in names), ","))
        end
    end
end

function run_baseline_2D_v12(; nominal_mesh=true, max_points=121)
    mkpath(OUTPUT_DIR_2D_v12)
    p = base_parameters_2D_v12(; nominal_mesh=nominal_mesh)
    specs = vcat(
        [(id=id, cooling=false) for id in HEATING_IDS_2D_v12],
        [(id=id, cooling=true) for id in COOLING_IDS_2D_v12],
    )
    cases = Vector{Any}(undef, length(specs))
    Threads.@threads for index in eachindex(specs)
        spec = specs[index]
        inputs = case_inputs_2D_v12(
            spec.id; cooling=spec.cooling, max_points=max_points,
        )
        cases[index] = simulate_case_2D_v12(inputs, p)
    end
    sensor_rows = reduce(vcat, case_metrics_2D_v12(case) for case in cases)
    final_rows = NamedTuple[]
    for case in cases
        model = case.model[end, :]
        observed = case.observed[end, :]
        push!(final_rows, (
            simulation_id=case.inputs.id,
            phase=case.inputs.cooling ? "cooling" : "heating",
            flow_lpm=mean(case.inputs.flow),
            T12_minus_T9_model=model[2] - model[4],
            T12_minus_T9_observed=observed[2] - observed[4],
            T11_minus_T10_model=model[3] - model[5],
            T11_minus_T10_observed=observed[3] - observed[5],
            T12_minus_T8_model=model[2] - model[1],
            T12_minus_T8_observed=observed[2] - observed[1],
            dp1_model_mbar=case.result.dp1_prediction_mbar[end],
            dp1_observed_mbar=case.inputs.dp1[end],
            adaptor_T_K=case.result.adaptor_temperature[end],
            housing_extension_T_K=
                case.result.housing_extension_temperature[end],
        ))
    end
    _write_namedtuples_csv_2D_v12(
        joinpath(OUTPUT_DIR_2D_v12, "baseline_sensor_metrics_2D_v12.csv"),
        sensor_rows,
    )
    _write_namedtuples_csv_2D_v12(
        joinpath(OUTPUT_DIR_2D_v12, "baseline_final_profiles_2D_v12.csv"),
        final_rows,
    )
    inventory = hardware_inventory2D(p)
    open(joinpath(OUTPUT_DIR_2D_v12, "hardware_inventory_2D_v12.txt"), "w") do io
        for name in propertynames(inventory)
            println(io, name, "=", getproperty(inventory, name))
        end
    end
    println("baseline mean sensor RMSE = ",
        mean(row.rmse_K for row in sensor_rows))
    println("baseline mean |t90 error| = ",
        mean(abs(row.t90_error_s) for row in sensor_rows))
    println("baseline radial sign mid = ",
        count(row -> row.phase == "heating" &&
            sign(row.T12_minus_T9_model) ==
            sign(row.T12_minus_T9_observed), final_rows), "/15")
    println("baseline radial sign deep = ",
        count(row -> row.phase == "heating" &&
            sign(row.T11_minus_T10_model) ==
            sign(row.T11_minus_T10_observed), final_rows), "/15")
    return (parameters=p, cases=cases, sensor_rows=sensor_rows,
            final_rows=final_rows)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_baseline_2D_v12()
end
