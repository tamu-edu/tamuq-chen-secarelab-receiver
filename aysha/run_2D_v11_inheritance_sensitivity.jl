# ============================================================================
# Broad one-at-a-time screen of parameters inherited from the v9 thermal fit.
# This is not a calibration. It tests whether v11's profile failures are
# robust to plausible solid/optical reparameterization.
# ============================================================================

include("run_2D_v11.jl")

const INHERITANCE_CASES_2D_v11 = (
    "E67", "E71", "E72", "E76", "E77", "E81",
)
const INHERITANCE_GROUPS_2D_v11 = (
    (456.0, ("E67", "E71")),
    (304.0, ("E72", "E76")),
    (256.0, ("E77", "E81")),
)

function with_solid_sensitivity_2D_v11(
    p;
    radial_scale=p.solid.radial_conductivity_scale,
    axial_scale=p.solid.axial_conductivity_scale,
    rad_extinction=p.solid.rad_extinction_coeff,
)
    s0 = p.solid
    solid = SolidProperties2D(
        density=s0.density,
        radial_conductivity_scale=Float64(radial_scale),
        axial_conductivity_scale=Float64(axial_scale),
        rad_extinction_coeff=Float64(rad_extinction),
        felt_conductivity_ref=s0.felt_conductivity_ref,
        felt_density=s0.felt_density,
        felt_heat_capacity=s0.felt_heat_capacity,
        casing_conductivity=s0.casing_conductivity,
        casing_density=s0.casing_density,
        casing_heat_capacity=s0.casing_heat_capacity,
        receiver_felt_contact_resistance=
            s0.receiver_felt_contact_resistance,
    )
    return ModelParameters2D(
        p.geometry, solid, p.heat_transfer,
        p.losses, p.optics, p.hydraulics,
    )
end

function with_optical_sensitivity_2D_v11(
    p;
    extinction=p.optics.extinction_coefficient,
    power_factor=1.0,
)
    o0 = p.optics
    optics = OpticalParameters2D(
        absorbed_fraction=o0.absorbed_fraction,
        extinction_coefficient=Float64(extinction),
        beam_radius_sigma=o0.beam_radius_sigma,
        spillage_fraction=o0.spillage_fraction,
        front_deposition_fraction=o0.front_deposition_fraction,
        scale_456=Float64(power_factor) * o0.scale_456,
        scale_304=Float64(power_factor) * o0.scale_304,
        scale_256=Float64(power_factor) * o0.scale_256,
    )
    return ModelParameters2D(
        p.geometry, p.solid, p.heat_transfer,
        p.losses, optics, p.hydraulics,
    )
end

function inheritance_variants_2D_v11(base)
    variants = Pair{String, ModelParameters2D}[]
    push!(variants, "baseline" => base)
    for scale in (0.10, 0.25, 0.50)
        push!(
            variants,
            "Sh_$(replace(string(scale), "." => "p"))" =>
                with_heat_model_2D_v11(
                    base;
                    model=:graetz,
                    nu_fd=3.61,
                    scale=scale,
                ),
        )
    end
    for value in (0.05, 0.30)
        push!(
            variants,
            "kr_$(replace(string(value), "." => "p"))" =>
                with_solid_sensitivity_2D_v11(
                    base; radial_scale=value,
                ),
        )
    end
    for value in (0.50, 1.00)
        push!(
            variants,
            "kz_$(replace(string(value), "." => "p"))" =>
                with_solid_sensitivity_2D_v11(
                    base; axial_scale=value,
                ),
        )
    end
    for value in (10.0, 50.0, 110.0)
        push!(
            variants,
            "betaopt_$(Int(value))" =>
                with_optical_sensitivity_2D_v11(
                    base; extinction=value,
                ),
        )
    end
    for value in (20.0, 300.0)
        push!(
            variants,
            "betarad_$(Int(value))" =>
                with_solid_sensitivity_2D_v11(
                    base; rad_extinction=value,
                ),
        )
    end
    for factor in (0.80, 1.20)
        push!(
            variants,
            "power_$(replace(string(factor), "." => "p"))" =>
                with_optical_sensitivity_2D_v11(
                    base; power_factor=factor,
                ),
        )
    end
    return variants
end

function inheritance_summary_2D_v11(case_rows)
    summaries = NamedTuple[]
    for label in unique(row.variant for row in case_rows)
        group = filter(row -> row.variant == label, case_rows)
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
        slopes = Float64[]
        observed_slopes = Float64[]
        for (_, ids) in INHERITANCE_GROUPS_2D_v11
            pair = [
                row for row in group if row.simulation_id in ids
            ]
            push!(
                slopes,
                linear_slope_2D_v11(
                    [row.mean_flow_lpm for row in pair],
                    [row.model_T12_minus_T8_K for row in pair],
                ).slope,
            )
            push!(
                observed_slopes,
                linear_slope_2D_v11(
                    [row.mean_flow_lpm for row in pair],
                    [row.observed_T12_minus_T8_K for row in pair],
                ).slope,
            )
        end
        push!(summaries, (
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
            slope_456_K_per_lpm=slopes[1],
            slope_304_K_per_lpm=slopes[2],
            slope_256_K_per_lpm=slopes[3],
            observed_slope_456_K_per_lpm=observed_slopes[1],
            observed_slope_304_K_per_lpm=observed_slopes[2],
            observed_slope_256_K_per_lpm=observed_slopes[3],
            all_slope_signs_correct=all(
                sign(slopes[i]) == sign(observed_slopes[i])
                for i in eachindex(slopes)
            ),
            T3_mae_K=mean(
                abs(row.model_T3_K - row.observed_T3_K)
                for row in group
            ),
            dp1_rmse_mbar=sqrt(mean(
                (row.model_dp1_mbar - row.observed_dp1_mbar)^2
                for row in group
            )),
        ))
    end
    return summaries
end

function main_inheritance_sensitivity_2D_v11()
    fitted_v9 = sensitivity_mesh_2D_v11(
        unpack_parameters2D(
            fitted_v9_theta_2D_v11(), default_parameters2D(),
        ),
    )
    graetz = with_heat_model_2D_v11(
        fitted_v9;
        model=:graetz,
        nu_fd=3.61,
        temperature_exponent=0.0,
        scale=1.0,
    )
    base = with_hydraulics_2D_v11(
        graetz;
        offset=-0.5428447496201336,
        resistance_scale=0.9701974617588522,
        old_hot_K=0.0,
        distribution=:equal_pressure,
        groove_diameter=13.0e-3,
        groove_K=184.15800344228472,
        max_iterations=16,
        tolerance=1.0e-5,
    )

    rows = NamedTuple[]
    for (label, p) in inheritance_variants_2D_v11(base)
        current_parameters_2D_v11[] = p
        println("Inheritance sensitivity: $label")
        flush(stdout)
        for id in INHERITANCE_CASES_2D_v11
            case = operating_case_2D_v11(id, p)
            result = simulate2D(
                p,
                case.op,
                case.times;
                initial_temperature=case.ambient[1],
            )
            push!(
                rows,
                final_case_row_2D_v11(label, id, case, result),
            )
        end
    end
    summaries = inheritance_summary_2D_v11(rows)
    write_namedtuple_csv_2D_v11(
        joinpath(
            OUTPUT_DIR_2D_v11,
            "inheritance_sensitivity_cases_2D_v11.csv",
        ),
        rows,
    )
    write_namedtuple_csv_2D_v11(
        joinpath(
            OUTPUT_DIR_2D_v11,
            "inheritance_sensitivity_summary_2D_v11.csv",
        ),
        summaries,
    )
    println("Inheritance sensitivity summary:")
    foreach(println, summaries)
    return summaries
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_inheritance_sensitivity_2D_v11()
end
