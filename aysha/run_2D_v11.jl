# ============================================================================
# run_2D_v11.jl
# No-refit Graetz / equal-pressure / measured-groove model-form matrix.
# ============================================================================

using LinearAlgebra
using Statistics
using Optim

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

include("2D_v11.jl")
include("import_exp_1D_v2.jl")
using .Receiver2D_v11

const OUTPUT_DIR_2D_v11 = joinpath(@__DIR__, "summaries", "2D_v11")
const REPRESENTATIVE_IDS_2D_v11 = ("E67", "E72", "E77")
const HEATING_GROUPS_2D_v11 = (
    ("456", IDs[1:5]),
    ("304", IDs[6:10]),
    ("256", IDs[11:15]),
)
const V9_HOT_EXCESS_K_2D_v11 = 118.47937747850334
const COLD_TOLERANCE_K_2D_v11 = 5.0

function fitted_v9_theta_2D_v11()
    path = joinpath(
        @__DIR__, "summaries", "2D_v9", "parameters_fitted_2D_v9.csv",
    )
    rows = split.(readlines(path)[2:end], ',')
    values = Dict(
        parse(Int, row[1]) => parse(Float64, row[3]) for row in rows
    )
    return [values[index] for index in 1:12]
end

function sensitivity_mesh_2D_v11(p)
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
        geometry, p.solid, p.heat_transfer, p.losses,
        p.optics, p.hydraulics,
    )
end

function with_heat_model_2D_v11(
    p;
    model=:graetz,
    nu_fd=3.61,
    temperature_exponent=0.0,
    scale=1.0,
)
    h0 = p.heat_transfer
    heat = HeatTransferParameters2D(
        coefficient=h0.coefficient,
        reynolds_exponent=h0.reynolds_exponent,
        prandtl_exponent=h0.prandtl_exponent,
        minimum_nusselt=h0.minimum_nusselt,
        c_radial_flow=h0.c_radial_flow,
        front_coefficient=0.0,
        front_reynolds_exponent=h0.front_reynolds_exponent,
        nusselt_model=model,
        fully_developed_nusselt=Float64(nu_fd),
        temperature_correction_exponent=Float64(temperature_exponent),
        heat_transfer_scale=Float64(scale),
    )
    return ModelParameters2D(
        p.geometry, p.solid, heat, p.losses, p.optics, p.hydraulics,
    )
end

function with_hydraulics_2D_v11(
    p;
    offset=p.hydraulics.dp1_zero_offset_mbar,
    resistance_scale=p.hydraulics.hydraulic_resistance_scale,
    old_hot_K=0.0,
    distribution=:equal_pressure,
    groove_diameter=13.0e-3,
    groove_K=0.0,
    max_iterations=16,
    tolerance=1.0e-5,
)
    h0 = p.hydraulics
    hydraulics = HydraulicParameters2D(
        standard_pressure=h0.standard_pressure,
        standard_temperature=h0.standard_temperature,
        atmospheric_pressure=h0.atmospheric_pressure,
        mass_flow_scale=h0.mass_flow_scale,
        dp1_zero_offset_mbar=Float64(offset),
        hydraulic_resistance_scale=Float64(resistance_scale),
        minor_loss_coefficient=Float64(old_hot_K),
        flow_distribution_model=distribution,
        groove_free_diameter=Float64(groove_diameter),
        groove_loss_coefficient=Float64(groove_K),
        equal_pressure_max_iterations=max_iterations,
        equal_pressure_tolerance=Float64(tolerance),
        equal_pressure_relaxation=0.55,
    )
    return ModelParameters2D(
        p.geometry, p.solid, p.heat_transfer,
        p.losses, p.optics, hydraulics,
    )
end

function cold_rows_2D_v11()
    rows = NamedTuple[]
    for meta in filter(
        row -> row.phase == :heating,
        hydraulic_t0_metadata,
    )
        selected =
            abs(meta.T3_K - meta.Tamb_K) <= COLD_TOLERANCE_K_2D_v11 &&
            abs(meta.Tsolid_mean_K - meta.Tamb_K) <=
                COLD_TOLERANCE_K_2D_v11
        push!(rows, merge(meta, (selected=selected,)))
    end
    return rows
end

function cold_receiver_dp_mbar_2D_v11(
    p, flow_lpm, Tsolid, Tamb,
)
    grid = Receiver2D_v11.build_grid2D(p)
    work = Receiver2D_v11.Workspace2D(grid, p)
    op = OperatingCondition2D(
        irradiance=0.0,
        flow_lpm=flow_lpm,
        inlet_temperature=Tamb,
        ambient_temperature=Tamb,
    )
    solid = fill(Float64(Tsolid), grid.nr_rec, grid.nz)
    rear = fill(Float64(Tsolid), grid.nz_rear)
    Receiver2D_v11._gas_profile2D!(
        work, solid, rear, 0.0, p, op, grid,
    )
    return work.common_dp / 100.0
end

function cold_hydraulic_objective_2D_v11(
    theta, selected, base;
    allow_groove,
    groove_diameter=13.0e-3,
)
    offset = theta[1]
    scale = theta[2]
    groove_K = allow_groove ? theta[3] : 0.0
    p = with_hydraulics_2D_v11(
        base;
        offset=offset,
        resistance_scale=scale,
        old_hot_K=0.0,
        distribution=:equal_pressure,
        groove_diameter=groove_diameter,
        groove_K=groove_K,
        max_iterations=20,
        tolerance=1.0e-7,
    )
    residuals = [
        offset + cold_receiver_dp_mbar_2D_v11(
            p, row.flow_lpm, row.Tsolid_mean_K, row.Tamb_K,
        ) - row.dp1_mbar for row in selected
    ]
    return mean(abs2, residuals)
end

function fit_cold_hydraulics_2D_v11(
    base;
    allow_groove,
    groove_diameter=13.0e-3,
)
    selected = filter(row -> row.selected, cold_rows_2D_v11())
    if allow_groove
        initial = [-0.60, 1.5, 5.0]
        lower = [-1.0, 0.1, 0.0]
        upper = [0.2, 4.0, 300.0]
    else
        initial = [-0.61, 1.95]
        lower = [-1.0, 0.1]
        upper = [0.2, 4.0]
    end
    objective = theta -> cold_hydraulic_objective_2D_v11(
        theta,
        selected,
        base;
        allow_groove=allow_groove,
        groove_diameter=groove_diameter,
    )
    result = Optim.optimize(
        objective,
        lower,
        upper,
        initial,
        Optim.Fminbox(Optim.NelderMead()),
        Optim.Options(
            iterations=500,
            f_reltol=1.0e-12,
            x_abstol=1.0e-9,
            show_trace=false,
        ),
    )
    theta = Optim.minimizer(result)
    fitted = with_hydraulics_2D_v11(
        base;
        offset=theta[1],
        resistance_scale=theta[2],
        old_hot_K=0.0,
        distribution=:equal_pressure,
        groove_diameter=groove_diameter,
        groove_K=allow_groove ? theta[3] : 0.0,
        max_iterations=16,
        tolerance=1.0e-5,
    )
    rows = NamedTuple[]
    for row in selected
        receiver_dp = cold_receiver_dp_mbar_2D_v11(
            fitted,
            row.flow_lpm,
            row.Tsolid_mean_K,
            row.Tamb_K,
        )
        predicted = fitted.hydraulics.dp1_zero_offset_mbar + receiver_dp
        push!(rows, (
            model=allow_groove ? "groove" : "no_groove",
            simulation_id=row.simulation_id,
            flow_lpm=row.flow_lpm,
            observed_dp1_mbar=row.dp1_mbar,
            predicted_dp1_mbar=predicted,
            residual_mbar=predicted - row.dp1_mbar,
        ))
    end
    n = length(rows)
    k = length(theta)
    sse = sum(abs2, row.residual_mbar for row in rows)
    aicc = n * log(max(sse / n, eps(Float64))) + 2k +
           2k * (k + 1) / max(n - k - 1, 1)
    return (
        parameters=fitted,
        theta=theta,
        rmse_mbar=sqrt(sse / n),
        sse=sse,
        aicc=aicc,
        converged=Optim.converged(result),
        iterations=Optim.iterations(result),
        rows=rows,
    )
end

function profile_groove_K_2D_v11(
    base,
    K_values;
    groove_diameter=13.0e-3,
)
    selected = filter(row -> row.selected, cold_rows_2D_v11())
    rows = NamedTuple[]
    for K in K_values
        objective = theta -> begin
            full = [theta[1], theta[2], Float64(K)]
            cold_hydraulic_objective_2D_v11(
                full,
                selected,
                base;
                allow_groove=true,
                groove_diameter=groove_diameter,
            )
        end
        result = Optim.optimize(
            objective,
            [-1.0, 0.1],
            [0.2, 4.0],
            [-0.60, 1.5],
            Optim.Fminbox(Optim.NelderMead()),
            Optim.Options(iterations=250, show_trace=false),
        )
        theta = Optim.minimizer(result)
        push!(rows, (
            groove_K=Float64(K),
            offset_mbar=theta[1],
            resistance_scale=theta[2],
            rmse_mbar=sqrt(Optim.minimum(result)),
        ))
    end
    return rows
end

function operating_case_2D_v11(id, p)
    times = Float64.(observation_time(measurements, id))
    flow = Float64.(observation(measurements, id, "_flow"))
    inlet = Float64.(observation(measurements, id, "_Tin"))
    ambient = Float64.(observation(measurements, id, "_Tamb"))
    nominal = Float64(simulation_conditions[id][Io])
    scale = nominal >= 400000.0 ? p.optics.scale_456 :
            (nominal >= 280000.0 ? p.optics.scale_304 :
             p.optics.scale_256)
    op = OperatingCondition2D(
        irradiance=linear_history(
            times, fill(scale * nominal, length(times)),
        ),
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

function observed_final_2D_v11(id, name)
    return Float64(observation(measurements, id, name)[end])
end

function linear_slope_2D_v11(x, y)
    xmean = mean(x)
    ymean = mean(y)
    sxx = sum((xi - xmean)^2 for xi in x)
    slope = sum(
        (x[i] - xmean) * (y[i] - ymean) for i in eachindex(x)
    ) / sxx
    correlation = cor(x, y)
    return (slope=slope, correlation=correlation)
end

function variants_2D_v11(base, cold_no_groove, cold_groove)
    historical = with_heat_model_2D_v11(base; model=:apparent)
    historical = with_hydraulics_2D_v11(
        historical;
        offset=-0.614226030202630,
        resistance_scale=1.95171196637134,
        old_hot_K=V9_HOT_EXCESS_K_2D_v11,
        distribution=:prescribed,
        groove_K=0.0,
    )

    function graetz_prescribed(nu_fd, exponent=0.0)
        p = with_heat_model_2D_v11(
            base;
            model=:graetz,
            nu_fd=nu_fd,
            temperature_exponent=exponent,
        )
        return with_hydraulics_2D_v11(
            p;
            offset=-0.614226030202630,
            resistance_scale=1.95171196637134,
            old_hot_K=0.0,
            distribution=:prescribed,
            groove_K=0.0,
        )
    end

    equal_open = with_heat_model_2D_v11(
        cold_no_groove.parameters;
        model=:graetz,
        nu_fd=3.61,
    )
    equal_groove = with_heat_model_2D_v11(
        cold_groove.parameters;
        model=:graetz,
        nu_fd=3.61,
    )

    return (
        ("apparent_prescribed_v9", historical),
        ("graetz_fd361_prescribed", graetz_prescribed(3.61)),
        ("graetz_fd298_prescribed", graetz_prescribed(2.98)),
        ("graetz_capuano366_prescribed", graetz_prescribed(3.66, 0.45)),
        ("graetz_fd361_equal_open", equal_open),
        ("graetz_fd361_equal_groove", equal_groove),
    )
end

function final_case_row_2D_v11(label, id, case, result)
    pred = sensor_predictions2D(result)
    grid_nz = size(result.gas_nusselt, 2)
    nearest_z(target) = argmin(abs.(result.z_solid .- target))
    j5 = nearest_z(5.0e-3)
    j58 = nearest_z(58.0e-3)
    j107 = nearest_z(107.0e-3)
    overlap = result.groove_overlap_fraction
    core_indices = findall(x -> x == 0.0, overlap)
    edge_indices = findall(x -> x == 1.0, overlap)
    area_flux = result.ring_mass_flow_kg_s[:, end] ./
                max.(result.r_gas .* 0.0 .+
                     Receiver2D_v11.build_grid2D(
                         # Geometry is reconstructed below only to obtain
                         # ring flow areas; no state is changed.
                         current_parameters_2D_v11[],
                     ).area_flow[1:result.nr_rec], eps(Float64))
    core_flux = isempty(core_indices) ? NaN :
                mean(area_flux[core_indices])
    edge_flux = isempty(edge_indices) ? NaN :
                mean(area_flux[edge_indices])
    return (
        variant=label,
        simulation_id=id,
        irradiance_kW_m2=case.nominal / 1000.0,
        mean_flow_lpm=mean(case.flow),
        model_T8_K=pred.T8[end],
        observed_T8_K=observed_final_2D_v11(id, "_T8"),
        model_T12_K=pred.T12[end],
        observed_T12_K=observed_final_2D_v11(id, "_T12"),
        model_T11_K=pred.T11[end],
        observed_T11_K=observed_final_2D_v11(id, "_T11"),
        model_T9_K=pred.T9[end],
        observed_T9_K=observed_final_2D_v11(id, "_T9"),
        model_T10_K=pred.T10[end],
        observed_T10_K=observed_final_2D_v11(id, "_T10"),
        model_T3_K=pred.T3[end],
        observed_T3_K=observed_final_2D_v11(id, "_T3"),
        model_T2_K=pred.T2[end],
        observed_T2_K=observed_final_2D_v11(id, "_T2"),
        model_T12_minus_T8_K=pred.T12[end] - pred.T8[end],
        observed_T12_minus_T8_K=
            observed_final_2D_v11(id, "_T12") -
            observed_final_2D_v11(id, "_T8"),
        model_T12_minus_T9_K=pred.T12[end] - pred.T9[end],
        observed_T12_minus_T9_K=
            observed_final_2D_v11(id, "_T12") -
            observed_final_2D_v11(id, "_T9"),
        model_T11_minus_T10_K=pred.T11[end] - pred.T10[end],
        observed_T11_minus_T10_K=
            observed_final_2D_v11(id, "_T11") -
            observed_final_2D_v11(id, "_T10"),
        model_dp1_mbar=result.dp1_prediction_mbar[end],
        observed_dp1_mbar=observed_final_2D_v11(id, "_DP1"),
        mean_Nu_5mm=mean(result.gas_nusselt[:, j5, end]),
        mean_Nu_58mm=mean(result.gas_nusselt[:, j58, end]),
        mean_Nu_107mm=mean(result.gas_nusselt[:, j107, end]),
        core_mass_flux_kg_m2s=core_flux,
        edge_mass_flux_kg_m2s=edge_flux,
        core_to_edge_mass_flux_ratio=core_flux / edge_flux,
        max_groove_dp_Pa=maximum(
            result.groove_pressure_drop_Pa[:, end],
        ),
        pressure_equalization_error=
            result.equal_pressure_relative_error[end],
        nz=grid_nz,
    )
end

# The case-row helper needs the exact active mesh areas without passing a
# second large argument through every summary call.
const current_parameters_2D_v11 = Ref{ModelParameters2D}()

function append_transient_rows_2D_v11!(
    rows, label, id, case, result,
)
    pred = sensor_predictions2D(result)
    observed = Dict(
        :T8 => Float64.(observation(measurements, id, "_T8")),
        :T12 => Float64.(observation(measurements, id, "_T12")),
        :T11 => Float64.(observation(measurements, id, "_T11")),
        :T9 => Float64.(observation(measurements, id, "_T9")),
        :T10 => Float64.(observation(measurements, id, "_T10")),
        :T3 => Float64.(observation(measurements, id, "_T3")),
        :T2 => Float64.(observation(measurements, id, "_T2")),
    )
    for k in eachindex(result.time)
        push!(rows, (
            variant=label,
            simulation_id=id,
            time_s=result.time[k],
            flow_lpm=case.flow[k],
            model_T8_K=pred.T8[k],
            observed_T8_K=observed[:T8][k],
            model_T12_K=pred.T12[k],
            observed_T12_K=observed[:T12][k],
            model_T11_K=pred.T11[k],
            observed_T11_K=observed[:T11][k],
            model_T9_K=pred.T9[k],
            observed_T9_K=observed[:T9][k],
            model_T10_K=pred.T10[k],
            observed_T10_K=observed[:T10][k],
            model_T3_K=pred.T3[k],
            observed_T3_K=observed[:T3][k],
            model_T2_K=pred.T2[k],
            observed_T2_K=observed[:T2][k],
            model_dp1_mbar=result.dp1_prediction_mbar[k],
            observed_dp1_mbar=Float64(
                observation(measurements, id, "_DP1")[k],
            ),
        ))
    end
    return nothing
end

function append_ring_rows_2D_v11!(rows, label, id, result, p)
    grid = Receiver2D_v11.build_grid2D(p)
    for i in 1:grid.nr_rec
        push!(rows, (
            variant=label,
            simulation_id=id,
            ring=i,
            radius_mm=1000.0 * grid.r_centers[i],
            groove_overlap=result.groove_overlap_fraction[i],
            ring_mass_flow_kg_s=result.ring_mass_flow_kg_s[i, end],
            ring_mass_flux_kg_m2s=
                result.ring_mass_flow_kg_s[i, end] /
                grid.area_flow[i],
            ring_dp_Pa=result.ring_pressure_drop_Pa[i, end],
            groove_dp_Pa=result.groove_pressure_drop_Pa[i, end],
            outlet_gas_K=result.gas_temperature[i, end, end],
        ))
    end
    return nothing
end

function append_axial_rows_2D_v11!(rows, label, id, result)
    id in REPRESENTATIVE_IDS_2D_v11 || return nothing
    core_i = 1
    edge_i = result.nr_rec
    for j in eachindex(result.z_solid)
        push!(rows, (
            variant=label,
            simulation_id=id,
            z_mm=1000.0 * result.z_solid[j],
            core_solid_K=result.solid_temperature[core_i, j, end],
            edge_solid_K=result.solid_temperature[edge_i, j, end],
            core_gas_K=result.gas_temperature[core_i, j, end],
            edge_gas_K=result.gas_temperature[edge_i, j, end],
            core_Nu=result.gas_nusselt[core_i, j, end],
            edge_Nu=result.gas_nusselt[edge_i, j, end],
        ))
    end
    return nothing
end

function write_namedtuple_csv_2D_v11(path, rows)
    isempty(rows) && error("cannot write empty table: $path")
    open(path, "w") do io
        println(io, join(string.(keys(first(rows))), ','))
        for row in rows
            println(io, join(values(row), ','))
        end
    end
    return path
end

function slope_rows_2D_v11(case_rows)
    rows = NamedTuple[]
    labels = unique(row.variant for row in case_rows)
    for label in labels
        for (irradiance, ids) in HEATING_GROUPS_2D_v11
            group = [
                row for row in case_rows
                if row.variant == label &&
                   row.simulation_id in ids
            ]
            length(group) >= 3 || continue
            flow = [row.mean_flow_lpm for row in group]
            model = [row.model_T12_minus_T8_K for row in group]
            observed = [row.observed_T12_minus_T8_K for row in group]
            mf = linear_slope_2D_v11(flow, model)
            of = linear_slope_2D_v11(flow, observed)
            push!(rows, (
                variant=label,
                irradiance_kW_m2=parse(Float64, irradiance),
                model_slope_K_per_lpm=mf.slope,
                observed_slope_K_per_lpm=of.slope,
                model_correlation=mf.correlation,
                observed_correlation=of.correlation,
            ))
        end
    end
    return rows
end

function aggregate_rows_2D_v11(case_rows, slopes)
    rows = NamedTuple[]
    for label in unique(row.variant for row in case_rows)
        group = filter(row -> row.variant == label, case_rows)
        sg = filter(row -> row.variant == label, slopes)
        axial_error = [
            row.model_T12_minus_T8_K - row.observed_T12_minus_T8_K
            for row in group
        ]
        mid_error = [
            row.model_T12_minus_T9_K - row.observed_T12_minus_T9_K
            for row in group
        ]
        deep_error = [
            row.model_T11_minus_T10_K - row.observed_T11_minus_T10_K
            for row in group
        ]
        t3_error = [
            row.model_T3_K - row.observed_T3_K for row in group
        ]
        dp_error = [
            row.model_dp1_mbar - row.observed_dp1_mbar for row in group
        ]
        slope_error = [
            row.model_slope_K_per_lpm -
            row.observed_slope_K_per_lpm for row in sg
        ]
        push!(rows, (
            variant=label,
            axial_mean_K=mean(
                row.model_T12_minus_T8_K for row in group
            ),
            axial_rmse_K=sqrt(mean(abs2, axial_error)),
            mid_radial_mean_K=mean(
                row.model_T12_minus_T9_K for row in group
            ),
            mid_radial_rmse_K=sqrt(mean(abs2, mid_error)),
            deep_radial_mean_K=mean(
                row.model_T11_minus_T10_K for row in group
            ),
            deep_radial_rmse_K=sqrt(mean(abs2, deep_error)),
            T3_mae_K=mean(abs, t3_error),
            dp1_rmse_mbar=sqrt(mean(abs2, dp_error)),
            dp1_bias_mbar=mean(dp_error),
            slope_rmse_K_per_lpm=isempty(slope_error) ? NaN :
                sqrt(mean(abs2, slope_error)),
            mean_core_edge_flow_ratio=mean(
                row.core_to_edge_mass_flux_ratio for row in group
                if isfinite(row.core_to_edge_mass_flux_ratio)
            ),
            max_pressure_equalization_error=maximum(
                row.pressure_equalization_error for row in group
            ),
        ))
    end
    return rows
end

function main_2D_v11()
    mkpath(OUTPUT_DIR_2D_v11)
    base = unpack_parameters2D(
        fitted_v9_theta_2D_v11(), default_parameters2D(),
    )
    base = sensitivity_mesh_2D_v11(base)

    println("Fitting cold t0 hydraulics without groove...")
    cold_open = fit_cold_hydraulics_2D_v11(
        base; allow_groove=false,
    )
    println("Fitting cold t0 hydraulics with measured groove...")
    cold_groove = fit_cold_hydraulics_2D_v11(
        base; allow_groove=true,
    )
    profile_rows = profile_groove_K_2D_v11(
        base, (0.0, 2.0, 5.0, 10.0, 20.0, 40.0, 80.0, 160.0),
    )
    write_namedtuple_csv_2D_v11(
        joinpath(OUTPUT_DIR_2D_v11, "cold_dp1_cases_2D_v11.csv"),
        vcat(cold_open.rows, cold_groove.rows),
    )
    write_namedtuple_csv_2D_v11(
        joinpath(OUTPUT_DIR_2D_v11, "cold_groove_profile_2D_v11.csv"),
        profile_rows,
    )

    println(
        "cold no-groove theta=$(cold_open.theta), " *
        "RMSE=$(cold_open.rmse_mbar), AICc=$(cold_open.aicc)",
    )
    println(
        "cold groove theta=$(cold_groove.theta), " *
        "RMSE=$(cold_groove.rmse_mbar), AICc=$(cold_groove.aicc)",
    )
    flush(stdout)

    case_rows = NamedTuple[]
    transient_rows = NamedTuple[]
    ring_rows = NamedTuple[]
    axial_rows = NamedTuple[]
    all_variants = variants_2D_v11(base, cold_open, cold_groove)
    requested = get(ENV, "V11_VARIANTS", "")
    if !isempty(requested)
        selected_labels = Set(split(requested, ','))
        all_variants = Tuple(
            item for item in all_variants if item[1] in selected_labels
        )
    end
    requested_ids = get(ENV, "V11_IDS", "")
    run_ids = isempty(requested_ids) ? IDs : split(requested_ids, ',')

    for (label, p) in all_variants
        current_parameters_2D_v11[] = p
        println("Running v11 variant $label")
        flush(stdout)
        for id in run_ids
            println("  $id")
            flush(stdout)
            case = operating_case_2D_v11(id, p)
            result = simulate2D(
                p,
                case.op,
                case.times;
                initial_temperature=case.ambient[1],
            )
            push!(
                case_rows,
                final_case_row_2D_v11(label, id, case, result),
            )
            append_ring_rows_2D_v11!(
                ring_rows, label, id, result, p,
            )
            append_axial_rows_2D_v11!(
                axial_rows, label, id, result,
            )
            if id in REPRESENTATIVE_IDS_2D_v11
                append_transient_rows_2D_v11!(
                    transient_rows, label, id, case, result,
                )
            end
        end
    end

    slopes = slope_rows_2D_v11(case_rows)
    aggregate = aggregate_rows_2D_v11(case_rows, slopes)
    write_namedtuple_csv_2D_v11(
        joinpath(OUTPUT_DIR_2D_v11, "model_form_cases_2D_v11.csv"),
        case_rows,
    )
    if !isempty(slopes)
        write_namedtuple_csv_2D_v11(
            joinpath(
                OUTPUT_DIR_2D_v11, "model_form_slopes_2D_v11.csv",
            ),
            slopes,
        )
    end
    write_namedtuple_csv_2D_v11(
        joinpath(OUTPUT_DIR_2D_v11, "model_form_summary_2D_v11.csv"),
        aggregate,
    )
    write_namedtuple_csv_2D_v11(
        joinpath(OUTPUT_DIR_2D_v11, "ring_profiles_2D_v11.csv"),
        ring_rows,
    )
    write_namedtuple_csv_2D_v11(
        joinpath(OUTPUT_DIR_2D_v11, "axial_profiles_2D_v11.csv"),
        axial_rows,
    )
    write_namedtuple_csv_2D_v11(
        joinpath(OUTPUT_DIR_2D_v11, "representative_transients_2D_v11.csv"),
        transient_rows,
    )

    summary_path = joinpath(
        OUTPUT_DIR_2D_v11, "model_form_summary_2D_v11.txt",
    )
    open(summary_path, "w") do io
        println(io, "cold_open_theta=$(cold_open.theta)")
        println(io, "cold_open_rmse_mbar=$(cold_open.rmse_mbar)")
        println(io, "cold_open_aicc=$(cold_open.aicc)")
        println(io, "cold_groove_theta=$(cold_groove.theta)")
        println(io, "cold_groove_rmse_mbar=$(cold_groove.rmse_mbar)")
        println(io, "cold_groove_aicc=$(cold_groove.aicc)")
        println(io, "groove_profile=$(profile_rows)")
        for row in aggregate
            println(io, row)
        end
    end
    println("v11 complete: $summary_path")
    return (
        cold_open=cold_open,
        cold_groove=cold_groove,
        case_rows=case_rows,
        slopes=slopes,
        aggregate=aggregate,
        summary_path=summary_path,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_2D_v11()
end
