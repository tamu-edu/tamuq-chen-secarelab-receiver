# ============================================================================
# diagnose_2D_v20_front_energy.jl
#
# No-refit audit requested by the independent v20 energy-closure assessment.
# This reports the front terms that are actually present in the v20 equations:
# receiver-web radiation, outer felt/casing radiation, and outer-front natural
# convection.  The outer radiative term is important because it was omitted
# from the model-free assessment's description of the implemented model.
# ============================================================================

using Statistics

include("run_2D_v20.jl")

const FRONT_AUDIT_POWERS_2D_v20 = (1.05, 1.05, 0.70)

function front_audit_parameters_2D_v20()
    exponent = 0.50
    nu50_reference = 3.07786e-4 * 50.0^1.44346
    nu50 = 2.00 * nu50_reference
    return parameters_2D_v20(
        mesh=:nominal,
        source_model=:near_deep,
        deep_fraction=0.90,
        deep_length_m=0.12,
        nu_prefactor=nu50 / 50.0^exponent,
        reynolds_exponent=exponent,
        distributed_tube_flange_h_W_m2_K=0.0,
        probe_capacity_areal_J_m2_K=30000.0,
        probe_stem_conductance_areal_W_m2_K=0.0,
        felt_conductivity_scale=0.15,
        felt_heat_capacity_scale=0.55,
        felt_contact_scale=0.30,
        power_scales=FRONT_AUDIT_POWERS_2D_v20,
        t3_location=:receiver_exit,
    )
end

function _front_energy_row_2D_v20(case)
    result = case.result
    p = result.parameters
    grid = V20.build_network_grid2D(p)
    bg = grid.base_grid
    p14 = p.base.base.base.base.base
    base = p14.base.base
    ambient = case.inputs.ambient
    final = max(1, floor(Int, 0.9 * length(result.time))):length(result.time)

    receiver_rad = zeros(length(result.time))
    outer_rad = zeros(length(result.time))
    outer_conv = zeros(length(result.time))
    front_fourth = zeros(length(result.time))
    outer_fourth = zeros(length(result.time))
    receiver_area = sum(grid.multiplicity) *
                    grid.solid_area_per_channel
    outer_area = sum(
        bg.area_frt[bg.nr_rec + outer]
        for outer in axes(result.outer_temperature, 1)
    )

    for index in eachindex(result.time)
        Ta = ambient[index]
        receiver_fourth_sum = 0.0
        for group in 1:grid.group_count
            area = grid.multiplicity[group] *
                   grid.solid_area_per_channel
            T = V20.V11.clamp_temperature(
                result.channel_temperature[group, 1, index],
            )
            receiver_fourth_sum += area * T^4
            receiver_rad[index] +=
                base.losses.front_loss_scale *
                base.losses.front_emissivity * V20.V11.SIGMA *
                area * (T^4 - Ta^4)
        end
        front_fourth[index] =
            (receiver_fourth_sum / receiver_area)^(0.25)

        outer_fourth_sum = 0.0
        for outer in axes(result.outer_temperature, 1)
            original = bg.nr_rec + outer
            area = bg.area_frt[original]
            T = V20.V11.clamp_temperature(
                result.outer_temperature[outer, 1, index],
            )
            outer_fourth_sum += area * T^4
            outer_rad[index] +=
                base.losses.front_loss_scale *
                base.losses.front_emissivity * V20.V11.SIGMA *
                area * (T^4 - Ta^4)
            outer_conv[index] +=
                base.losses.front_loss_scale *
                V20.V11.front_nusselt_correlation(T, Ta) *
                area * (T - Ta)
        end
        outer_fourth[index] =
            (outer_fourth_sum / outer_area)^(0.25)
    end

    scale = case.inputs.nominal >= 400000.0 ?
        p.optics.scale_456 :
        (case.inputs.nominal >= 280000.0 ?
         p.optics.scale_304 : p.optics.scale_256)
    absorbed = p.optics.absorbed_fraction * scale *
               case.inputs.nominal * p.geometry.receiver_width^2
    total_rad = receiver_rad .+ outer_rad
    total_front = total_rad .+ outer_conv
    return (
        id=case.inputs.id,
        nominal_kW_m2=case.inputs.nominal / 1000.0,
        flow_L_min=mean(case.inputs.flow[final]),
        receiver_front_T_K=mean(front_fourth[final]),
        outer_front_T_K=mean(outer_fourth[final]),
        absorbed_power_W=absorbed,
        receiver_front_radiation_W=mean(receiver_rad[final]),
        outer_front_radiation_W=mean(outer_rad[final]),
        total_front_radiation_W=mean(total_rad[final]),
        outer_front_convection_W=mean(outer_conv[final]),
        total_front_loss_W=mean(total_front[final]),
        front_radiation_share=mean(total_rad[final]) / absorbed,
        total_front_loss_share=mean(total_front[final]) / absorbed,
        receiver_radiating_area_cm2=1e4 * receiver_area,
        outer_radiating_area_cm2=1e4 * outer_area,
    )
end

function diagnose_2D_v20_front_energy(; max_points=61)
    p = front_audit_parameters_2D_v20()
    rows = Vector{NamedTuple}(undef, length(HEATING_IDS_2D_v20))
    Threads.@threads for index in eachindex(HEATING_IDS_2D_v20)
        id = HEATING_IDS_2D_v20[index]
        println("v20 front-energy audit $index/15: $id")
        flush(stdout)
        inputs = case_inputs_2D_v20(
            id; cooling=false, max_points=max_points,
        )
        case = simulate_case_2D_v20(
            inputs, p; initialization=:ambient,
            reltol=5e-4, abstol=1e-4, dtmax=120.0,
        )
        rows[index] = _front_energy_row_2D_v20(case)
    end
    path = joinpath(
        OUTPUT_DIR_2D_v20, "front_energy_no_refit_2D_v20.csv",
    )
    _write_namedtuples_csv_2D_v20(path, rows)
    println("\nNo-refit v20 front-energy audit")
    println("ID  Tfront  Touter  Qabs  Qrad,rec  Qrad,outer  Qrad,total  Qconv  rad/Qabs  front/Qabs")
    for row in rows
        println(
            row.id, "  ",
            round(row.receiver_front_T_K; digits=1), "  ",
            round(row.outer_front_T_K; digits=1), "  ",
            round(row.absorbed_power_W; digits=1), "  ",
            round(row.receiver_front_radiation_W; digits=1), "  ",
            round(row.outer_front_radiation_W; digits=1), "  ",
            round(row.total_front_radiation_W; digits=1), "  ",
            round(row.outer_front_convection_W; digits=1), "  ",
            round(100 * row.front_radiation_share; digits=1), "%  ",
            round(100 * row.total_front_loss_share; digits=1), "%",
        )
    end
    println("wrote $path")
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    points = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 61
    diagnose_2D_v20_front_energy(; max_points=points)
end
