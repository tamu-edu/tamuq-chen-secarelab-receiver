# ============================================================================
# run_2D_v10_exponent_sensitivity.jl
# Reference-normalized test of stronger front HTC Reynolds dependence.
# ============================================================================

include("run_2D_v10_sensitivity.jl")

const FRONT_EQUIVALENT_STRENGTHS_2D_v10 = (2.0, 4.0)
const FRONT_EXPONENTS_2D_v10 = (0.5, 0.8, 1.0, 1.25, 1.5)
const FRONT_REFERENCE_FLOW_LPM_2D_v10 = 10.0

function representative_front_reynolds_2D_v10(p)
    grid = Receiver2D_v10.build_grid2D(p)
    radius = grid.r_faces[grid.nr_rec + 1]
    psi = [
        1.0 - p.heat_transfer.c_radial_flow *
        (grid.r_centers[i] / radius)^2 for i in 1:grid.nr_rec
    ]
    psi_sum = sum(
        psi[i] * grid.area_flow[i] for i in 1:grid.nr_rec
    )
    fractions = [
        psi[i] * grid.area_flow[i] / psi_sum for i in 1:grid.nr_rec
    ]
    mdot_total = standard_mass_flow2D(
        FRONT_REFERENCE_FLOW_LPM_2D_v10, p.hydraulics,
    )
    mu_ref = air_viscosity(p.hydraulics.standard_temperature)
    ring_reynolds = [
        mdot_total * fractions[i] * grid.dh /
        (grid.area_flow[i] * mu_ref) for i in 1:grid.nr_rec
    ]
    return sum(fractions .* ring_reynolds)
end

function with_normalized_front_law_2D_v10(
    p, equivalent_strength, exponent, reference_reynolds,
)
    h0 = p.heat_transfer
    raw_coefficient = Float64(equivalent_strength) *
                      reference_reynolds^(0.5 - Float64(exponent))
    heat = HeatTransferParameters2D(
        coefficient=h0.coefficient,
        reynolds_exponent=h0.reynolds_exponent,
        prandtl_exponent=h0.prandtl_exponent,
        minimum_nusselt=h0.minimum_nusselt,
        c_radial_flow=h0.c_radial_flow,
        front_coefficient=raw_coefficient,
        front_reynolds_exponent=Float64(exponent),
    )
    parameters = ModelParameters2D(
        p.geometry, p.solid, heat, p.losses, p.optics, p.hydraulics,
    )
    return (parameters=parameters, raw_coefficient=raw_coefficient)
end

function exponent_case_row_2D_v10(
    id, equivalent_strength, exponent, raw_coefficient,
    reference_reynolds, p,
)
    case = operating_case_2D_v10(id, p)
    result = simulate2D(
        p, case.op, case.times; initial_temperature=case.ambient[1],
    )
    prediction = sensor_predictions2D(result)
    observed(name) = Float64(observation(measurements, id, name)[end])
    return (
        equivalent_strength=equivalent_strength,
        exponent=exponent,
        raw_coefficient=raw_coefficient,
        reference_reynolds=reference_reynolds,
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
        mass_flow_kg_s=result.mass_flow_kg_s[end],
        dp1_model_mbar=result.dp1_prediction_mbar[end],
        dp1_observed_mbar=observed("_DP1"),
    )
end

function combination_summary_2D_v10(case_rows, slope_rows)
    rows = NamedTuple[]
    for equivalent_strength in FRONT_EQUIVALENT_STRENGTHS_2D_v10
        for exponent in FRONT_EXPONENTS_2D_v10
            group = filter(
                row -> row.equivalent_strength == equivalent_strength &&
                       row.exponent == exponent,
                case_rows,
            )
            slopes = filter(
                row -> row.equivalent_strength == equivalent_strength &&
                       row.exponent == exponent,
                slope_rows,
            )
            slope(irradiance) = only(
                row.model_slope_K_per_lpm for row in slopes
                if row.irradiance_kW_m2 == irradiance
            )
            observed_slope(irradiance) = only(
                row.observed_slope_K_per_lpm for row in slopes
                if row.irradiance_kW_m2 == irradiance
            )
            axial_errors = [
                row.model_T12_minus_T8_K - row.observed_T12_minus_T8_K
                for row in group
            ]
            mid_errors = [
                row.model_T12_minus_T9_K - row.observed_T12_minus_T9_K
                for row in group
            ]
            deep_errors = [
                row.model_T11_minus_T10_K - row.observed_T11_minus_T10_K
                for row in group
            ]
            t3_errors = [
                row.model_T3_K - row.observed_T3_K for row in group
            ]
            dp1_errors = [
                row.dp1_model_mbar - row.dp1_observed_mbar for row in group
            ]
            slope_errors = [
                slope(irradiance) - observed_slope(irradiance)
                for irradiance in (456.0, 304.0, 256.0)
            ]
            push!(rows, (
                equivalent_strength=equivalent_strength,
                exponent=exponent,
                raw_coefficient=first(group).raw_coefficient,
                reference_reynolds=first(group).reference_reynolds,
                axial_mean_K=mean(
                    row.model_T12_minus_T8_K for row in group
                ),
                axial_rmse_K=sqrt(mean(abs2, axial_errors)),
                mid_radial_mean_K=mean(
                    row.model_T12_minus_T9_K for row in group
                ),
                mid_radial_rmse_K=sqrt(mean(abs2, mid_errors)),
                deep_radial_mean_K=mean(
                    row.model_T11_minus_T10_K for row in group
                ),
                deep_radial_rmse_K=sqrt(mean(abs2, deep_errors)),
                T3_mae_K=mean(abs, t3_errors),
                slope_456_K_per_lpm=slope(456.0),
                slope_304_K_per_lpm=slope(304.0),
                slope_256_K_per_lpm=slope(256.0),
                slope_rmse_K_per_lpm=sqrt(mean(abs2, slope_errors)),
                mean_front_h_W_m2K=mean(
                    row.mean_front_h_W_m2K for row in group
                ),
                mean_front_to_gas_W=mean(
                    row.front_to_gas_W for row in group
                ),
                mean_gas_preheat_K=mean(
                    row.mean_gas_preheat_K for row in group
                ),
                mean_dp1_model_mbar=mean(
                    row.dp1_model_mbar for row in group
                ),
                dp1_rmse_mbar=sqrt(mean(abs2, dp1_errors)),
                dp1_bias_mbar=mean(dp1_errors),
            ))
        end
    end
    return rows
end

function main_2D_v10_exponent_sensitivity()
    mkpath(OUTPUT_DIR_2D_v10)
    base = unpack_parameters2D(
        fitted_v9_theta_2D_v10(), default_parameters2D(),
    )
    base = with_v9_hot_hydraulics_2D_v10(base)
    base = sensitivity_mesh_2D_v10(base)
    reference_reynolds = representative_front_reynolds_2D_v10(base)
    println("reference_reynolds=$reference_reynolds")

    case_rows = NamedTuple[]
    slope_rows = NamedTuple[]
    for equivalent_strength in FRONT_EQUIVALENT_STRENGTHS_2D_v10
        for exponent in FRONT_EXPONENTS_2D_v10
            law = with_normalized_front_law_2D_v10(
                base, equivalent_strength, exponent, reference_reynolds,
            )
            println(
                "v10 exponent sensitivity C_eq=$equivalent_strength, " *
                "m=$exponent, raw_C=$(law.raw_coefficient)",
            )
            flush(stdout)
            for id in IDs
                println("  $id")
                push!(
                    case_rows,
                    exponent_case_row_2D_v10(
                        id, equivalent_strength, exponent,
                        law.raw_coefficient, reference_reynolds,
                        law.parameters,
                    ),
                )
            end
            selected = filter(
                row -> row.equivalent_strength == equivalent_strength &&
                       row.exponent == exponent,
                case_rows,
            )
            for (irradiance, group_ids) in HEATING_GROUPS_2D_v10
                group = [
                    only(filter(row -> row.simulation_id == id, selected))
                    for id in group_ids
                ]
                flow = [row.mean_flow_lpm for row in group]
                model_axial = [row.model_T12_minus_T8_K for row in group]
                observed_axial = [
                    row.observed_T12_minus_T8_K for row in group
                ]
                model_fit = linear_slope_2D_v10(flow, model_axial)
                observed_fit = linear_slope_2D_v10(flow, observed_axial)
                push!(slope_rows, (
                    equivalent_strength=equivalent_strength,
                    exponent=exponent,
                    irradiance_kW_m2=parse(Float64, irradiance),
                    model_slope_K_per_lpm=model_fit.slope,
                    model_correlation=model_fit.correlation,
                    observed_slope_K_per_lpm=observed_fit.slope,
                    observed_correlation=observed_fit.correlation,
                ))
            end
        end
    end

    aggregate_rows = combination_summary_2D_v10(case_rows, slope_rows)
    case_path = joinpath(
        OUTPUT_DIR_2D_v10, "front_exponent_cases_2D_v10.csv",
    )
    slope_path = joinpath(
        OUTPUT_DIR_2D_v10, "front_exponent_slopes_2D_v10.csv",
    )
    aggregate_path = joinpath(
        OUTPUT_DIR_2D_v10, "front_exponent_summary_2D_v10.csv",
    )
    write_namedtuple_csv_2D_v10(case_path, case_rows)
    write_namedtuple_csv_2D_v10(slope_path, slope_rows)
    write_namedtuple_csv_2D_v10(aggregate_path, aggregate_rows)

    text_path = joinpath(
        OUTPUT_DIR_2D_v10, "front_exponent_summary_2D_v10.txt",
    )
    open(text_path, "w") do io
        println(io, "2D_v10 reference-normalized front exponent sensitivity")
        println(io, "reference_flow_lpm=$FRONT_REFERENCE_FLOW_LPM_2D_v10")
        println(io, "reference_reynolds=$reference_reynolds")
        println(
            io,
            "equivalent_strengths=" *
            "$(collect(FRONT_EQUIVALENT_STRENGTHS_2D_v10))",
        )
        println(io, "exponents=$(collect(FRONT_EXPONENTS_2D_v10))")
        println(io, "v9 thermal/hydraulic and centered optical fields fixed")
        for row in aggregate_rows
            println(
                io,
                "Ceq=$(row.equivalent_strength), m=$(row.exponent): " *
                "axial mean=$(row.axial_mean_K) K, " *
                "axial RMSE=$(row.axial_rmse_K) K, " *
                "slope RMSE=$(row.slope_rmse_K_per_lpm) K/(L/min), " *
                "slopes=($(row.slope_456_K_per_lpm), " *
                "$(row.slope_304_K_per_lpm), " *
                "$(row.slope_256_K_per_lpm)), " *
                "T3 MAE=$(row.T3_mae_K) K, " *
                "DP1 RMSE=$(row.dp1_rmse_mbar) mbar, " *
                "mid/deep radial means=($(row.mid_radial_mean_K), " *
                "$(row.deep_radial_mean_K)) K",
            )
        end
    end
    println("case_table=$case_path")
    println("slope_table=$slope_path")
    println("aggregate_table=$aggregate_path")
    println("summary=$text_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_2D_v10_exponent_sensitivity()
end
