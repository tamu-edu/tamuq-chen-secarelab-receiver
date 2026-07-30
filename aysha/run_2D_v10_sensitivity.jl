# ============================================================================
# run_2D_v10_sensitivity.jl
# Pre-calibration sweep of the conservative front-web-to-inlet-gas coefficient.
# The v9 fitted thermal/hydraulic parameters and centered irradiance are fixed.
# ============================================================================

using Statistics

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

include("2D_v10.jl")
include("import_exp_1D_v2.jl")
using .Receiver2D_v10

const OUTPUT_DIR_2D_v10 = joinpath(@__DIR__, "summaries", "2D_v10")
const FRONT_COEFFICIENT_SWEEP_2D_v10 = (
    0.0, 0.05, 0.10, 0.25, 0.50, 1.00, 2.00, 4.00, 8.00,
)
const HEATING_GROUPS_2D_v10 = (
    ("456", IDs[1:5]),
    ("304", IDs[6:10]),
    ("256", IDs[11:15]),
)

function fitted_v9_theta_2D_v10()
    path = joinpath(
        @__DIR__, "summaries", "2D_v9", "parameters_fitted_2D_v9.csv",
    )
    rows = split.(readlines(path)[2:end], ',')
    values = Dict(parse(Int, row[1]) => parse(Float64, row[3]) for row in rows)
    return [values[index] for index in 1:12]
end

function with_v9_hot_hydraulics_2D_v10(p)
    h0 = p.hydraulics
    hydraulics = HydraulicParameters2D(
        standard_pressure=h0.standard_pressure,
        standard_temperature=h0.standard_temperature,
        atmospheric_pressure=h0.atmospheric_pressure,
        mass_flow_scale=h0.mass_flow_scale,
        dp1_zero_offset_mbar=h0.dp1_zero_offset_mbar,
        hydraulic_resistance_scale=h0.hydraulic_resistance_scale,
        minor_loss_coefficient=118.479,
    )
    return ModelParameters2D(
        p.geometry, p.solid, p.heat_transfer, p.losses, p.optics, hydraulics,
    )
end

function with_front_coefficient_2D_v10(p, coefficient)
    h0 = p.heat_transfer
    heat = HeatTransferParameters2D(
        coefficient=h0.coefficient,
        reynolds_exponent=h0.reynolds_exponent,
        prandtl_exponent=h0.prandtl_exponent,
        minimum_nusselt=h0.minimum_nusselt,
        c_radial_flow=h0.c_radial_flow,
        front_coefficient=Float64(coefficient),
        front_reynolds_exponent=0.5,
    )
    return ModelParameters2D(
        p.geometry, p.solid, heat, p.losses, p.optics, p.hydraulics,
    )
end

function sensitivity_mesh_2D_v10(p)
    g0 = p.geometry
    geometry = Geometry2D(
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
    return ModelParameters2D(
        geometry, p.solid, p.heat_transfer, p.losses, p.optics, p.hydraulics,
    )
end

function operating_case_2D_v10(id, p)
    times = Float64.(observation_time(measurements, id))
    flow = Float64.(observation(measurements, id, "_flow"))
    inlet = Float64.(observation(measurements, id, "_Tin"))
    ambient = Float64.(observation(measurements, id, "_Tamb"))
    nominal = Float64(simulation_conditions[id][Io])
    scale = nominal >= 400000.0 ? p.optics.scale_456 :
            (nominal >= 280000.0 ? p.optics.scale_304 : p.optics.scale_256)
    op = OperatingCondition2D(
        irradiance=linear_history(times, fill(scale * nominal, length(times))),
        flow_lpm=linear_history(times, flow),
        inlet_temperature=linear_history(times, inlet),
        ambient_temperature=linear_history(times, ambient),
    )
    return (
        times=times,
        flow=flow,
        inlet=inlet,
        ambient=ambient,
        nominal=nominal,
        op=op,
    )
end

function front_external_terms_2D_v10(result, p, ambient)
    grid = Receiver2D_v10.build_grid2D(p)
    Ts = view(result.solid_temperature, :, 1, size(result.solid_temperature, 3))
    radiation = 0.0
    convection = 0.0
    for i in 1:grid.nr_total
        area = i <= grid.nr_rec ? grid.area_solid[i] : grid.area_frt[i]
        radiation += p.losses.front_loss_scale * p.losses.front_emissivity *
                     Receiver2D_v10.SIGMA * area *
                     (max(200.0, Ts[i])^4 - ambient^4)
        if i > grid.nr_rec
            convection += p.losses.front_loss_scale *
                          front_nusselt_correlation(Ts[i], ambient) * area *
                          (Ts[i] - ambient)
        end
    end
    return (radiation=radiation, convection=convection)
end

function final_sensitivity_row_2D_v10(id, coefficient, p)
    case = operating_case_2D_v10(id, p)
    result = simulate2D(
        p, case.op, case.times; initial_temperature=case.ambient[1],
    )
    prediction = sensor_predictions2D(result)
    external = front_external_terms_2D_v10(result, p, case.ambient[end])
    observed(name) = Float64(observation(measurements, id, name)[end])
    return (
        front_coefficient=coefficient,
        simulation_id=id,
        irradiance_kW_m2=case.nominal / 1000.0,
        mean_flow_lpm=mean(case.flow),
        model_T12_minus_T8_K=prediction.T12[end] - prediction.T8[end],
        observed_T12_minus_T8_K=observed("_T12") - observed("_T8"),
        model_T12_minus_T9_K=prediction.T12[end] - prediction.T9[end],
        observed_T12_minus_T9_K=observed("_T12") - observed("_T9"),
        model_T11_minus_T10_K=prediction.T11[end] - prediction.T10[end],
        observed_T11_minus_T10_K=observed("_T11") - observed("_T10"),
        model_T3_K=prediction.T3[end],
        observed_T3_K=observed("_T3"),
        front_to_gas_W=sum(result.front_gas_heat_transfer_W[:, end]),
        mean_front_h_W_m2K=mean(result.front_heat_transfer_coefficient[:, end]),
        mean_gas_preheat_K=mean(
            result.gas_temperature[:, 1, end] .- case.inlet[end],
        ),
        front_radiation_W=external.radiation,
        remaining_front_convection_W=external.convection,
        mass_flow_kg_s=result.mass_flow_kg_s[end],
        dp1_model_mbar=result.dp1_prediction_mbar[end],
    )
end

function linear_slope_2D_v10(x, y)
    xm = mean(x)
    ym = mean(y)
    slope = sum((x .- xm) .* (y .- ym)) / sum(abs2, x .- xm)
    return (
        slope=slope,
        correlation=cor(x, y),
        intercept=ym - slope * xm,
    )
end

function write_namedtuple_csv_2D_v10(path, rows)
    isempty(rows) && error("cannot write an empty sensitivity table")
    open(path, "w") do io
        println(io, join(propertynames(rows[1]), ','))
        for row in rows
            println(io, join(values(row), ','))
        end
    end
end

function main_2D_v10_sensitivity()
    mkpath(OUTPUT_DIR_2D_v10)
    base = unpack_parameters2D(
        fitted_v9_theta_2D_v10(), default_parameters2D(),
    )
    base = with_v9_hot_hydraulics_2D_v10(base)

    case_rows = NamedTuple[]
    slope_rows = NamedTuple[]
    for coefficient in FRONT_COEFFICIENT_SWEEP_2D_v10
        p = sensitivity_mesh_2D_v10(
            with_front_coefficient_2D_v10(base, coefficient),
        )
        println("v10 sensitivity C_f=$coefficient")
        for id in IDs
            println("  $id")
            push!(case_rows, final_sensitivity_row_2D_v10(id, coefficient, p))
        end
        selected = filter(row -> row.front_coefficient == coefficient, case_rows)
        for (irradiance, group_ids) in HEATING_GROUPS_2D_v10
            group = [
                only(filter(row -> row.simulation_id == id, selected))
                for id in group_ids
            ]
            flow = [row.mean_flow_lpm for row in group]
            model_axial = [row.model_T12_minus_T8_K for row in group]
            observed_axial = [row.observed_T12_minus_T8_K for row in group]
            model_fit = linear_slope_2D_v10(flow, model_axial)
            observed_fit = linear_slope_2D_v10(flow, observed_axial)
            push!(slope_rows, (
                front_coefficient=coefficient,
                irradiance_kW_m2=parse(Float64, irradiance),
                model_slope_K_per_lpm=model_fit.slope,
                model_correlation=model_fit.correlation,
                observed_slope_K_per_lpm=observed_fit.slope,
                observed_correlation=observed_fit.correlation,
            ))
        end
    end

    case_path = joinpath(
        OUTPUT_DIR_2D_v10, "front_sensitivity_cases_2D_v10.csv",
    )
    slope_path = joinpath(
        OUTPUT_DIR_2D_v10, "front_sensitivity_slopes_2D_v10.csv",
    )
    write_namedtuple_csv_2D_v10(case_path, case_rows)
    write_namedtuple_csv_2D_v10(slope_path, slope_rows)

    summary_path = joinpath(
        OUTPUT_DIR_2D_v10, "front_sensitivity_summary_2D_v10.txt",
    )
    open(summary_path, "w") do io
        println(io, "2D_v10 conservative front-exchange pre-fit sensitivity")
        println(io, "mesh=(nr_rec=10,nr_felt=5,nr_case=2,nz=30,nz_rear=20)")
        println(io, "front_reynolds_exponent=0.5")
        println(io, "front_coefficients=$(collect(FRONT_COEFFICIENT_SWEEP_2D_v10))")
        println(io, "v9 thermal and hydraulic parameters fixed")
        println(io, "irradiance pattern unchanged and centered")
        for row in slope_rows
            println(
                io,
                "Cf=$(row.front_coefficient), I=$(row.irradiance_kW_m2): " *
                "model slope=$(row.model_slope_K_per_lpm), " *
                "observed slope=$(row.observed_slope_K_per_lpm)",
            )
        end
        println(io, "selected aggregate diagnostics:")
        for coefficient in (0.0, 1.0, 2.0, 4.0, 8.0)
            group = filter(
                row -> row.front_coefficient == coefficient, case_rows,
            )
            axial_error = [
                row.model_T12_minus_T8_K - row.observed_T12_minus_T8_K
                for row in group
            ]
            t3_error = [
                row.model_T3_K - row.observed_T3_K for row in group
            ]
            println(
                io,
                "Cf=$coefficient: mean axial=" *
                "$(mean(row.model_T12_minus_T8_K for row in group)) K, " *
                "axial RMSE=$(sqrt(mean(abs2, axial_error))) K, " *
                "mean mid radial=" *
                "$(mean(row.model_T12_minus_T9_K for row in group)) K, " *
                "mean deep radial=" *
                "$(mean(row.model_T11_minus_T10_K for row in group)) K, " *
                "T3 MAE=$(mean(abs, t3_error)) K, " *
                "mean h_front=" *
                "$(mean(row.mean_front_h_W_m2K for row in group)) W/m2/K, " *
                "mean preheat=" *
                "$(mean(row.mean_gas_preheat_K for row in group)) K",
            )
        end
        println(io, "acceptance=FAIL")
        println(
            io,
            "reason=axial mean can be forced only while measured flow slopes " *
            "collapse/reverse; radial signs remain wrong",
        )
        println(io, "calibration_action=do not fit C_f with m_f=0.5")
    end
    println("case_table=$case_path")
    println("slope_table=$slope_path")
    println("summary=$summary_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_2D_v10_sensitivity()
end
