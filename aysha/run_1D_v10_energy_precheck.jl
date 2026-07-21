# ============================================================================
# run_1D_v10_energy_precheck.jl - v10 energy and axial-transport diagnostics
# ============================================================================
# This runner checks where the fixed absorbed power goes before v10b fitting.
# It intentionally uses eta_opt = 1.0 and front_dep = 1.0 from pnew_v10.
# ============================================================================

begin
    using Statistics
    using StatsPlots
end

include("1D_v10.jl")

begin # output configuration
    const PRECHECK_OUTPUT_DIR_v10 = joinpath(@__DIR__, "summaries", "1D_v10", "precheck")
    const PRECHECK_PLOT_DIR_v10 = joinpath(PRECHECK_OUTPUT_DIR_v10, "plots")
    const PRECHECK_AXIAL_DIR_v10 = joinpath(PRECHECK_OUTPUT_DIR_v10, "axial")
    mkpath(PRECHECK_OUTPUT_DIR_v10)
    mkpath(PRECHECK_PLOT_DIR_v10)
    mkpath(PRECHECK_AXIAL_DIR_v10)
end

begin # variants and cases
    precheck_variant_params_v10(beta_tr) = begin
        p = copy(pnew_v10)
        p[2] = 1.0 / beta_tr
        p
    end

    const PRECHECK_VARIANTS_v10 = (
        (name="no_rad", params=copy(pnew_v10), use_rad=false, profile=:uniform),
        (name="uniform_beta1000", params=precheck_variant_params_v10(1000.0),
         use_rad=true, profile=:uniform),
        (name="front_strong", params=copy(pnew_v10),
         use_rad=true, profile=:front_strong),
        (name="rear_strong", params=copy(pnew_v10),
         use_rad=true, profile=:rear_strong),
    )

    const PRECHECK_CASES_v10 = (
        (simulation_id="E67", is_cooling=false),
        (simulation_id="E76", is_cooling=false),
        (simulation_id="E80", is_cooling=false),
        (simulation_id="C80", is_cooling=true),
    )
end

safe_string_v10(value) = replace(String(value), ":" => "_")

function write_csv_v10(path, header, rows)
    open(path, "w") do io
        println(io, join(header, ','))
        for row in rows
            println(io, join((getproperty(row, Symbol(column)) for column in header), ','))
        end
    end
    return path
end

function gas_heat_gain_cells_v10(Tg, first_index, cell_count, mdot)
    qgas = zeros(cell_count)
    if mdot <= 1e-12
        return qgas
    end
    for i in 1:cell_count
        left = first_index + i - 1
        Tfilm = 0.5 * (Tg[left] + Tg[left + 1])
        cp = cpf_f(Tfilm)
        qgas[i] = mdot * cp * (Tg[left + 1] - Tg[left])
    end
    return qgas
end

function final_energy_diagnostic_v10(params, simulation_id;
                                     is_cooling=false,
                                     use_rosseland_radiation=false,
                                     rosseland_profile=:uniform,
                                     nodes=default_nodes)
    model, result, experiment = solve_case_v10(
        params, simulation_id; is_cooling=is_cooling, nodes=nodes,
        use_rosseland_radiation=use_rosseland_radiation,
        rosseland_profile=rosseland_profile,
    )
    conditions, _time, _experiment_full, data = experimental_case_v10(
        simulation_id; is_cooling=is_cooling,
    )

    final_index = length(result.time)
    z = result.z_solid
    z_rear = result.z_rear_tube .- L
    dx = L / length(z)
    dx_rear = REAR_TUBE_LENGTH_v10 / length(z_rear)
    Ts = result.solid_temperature[:, final_index]
    Ttube = result.rear_tube_temperature[:, final_index]
    Tg = result.gas_temperature[:, final_index]
    h = result.heat_transfer_coefficient[:, final_index]
    Tcavity = result.cavity_temperature[final_index]
    ambient = observation(data, simulation_id, "_Tamb")[final_index]
    flow = observation(data, simulation_id, "_flow")[final_index]
    Tin = observation(data, simulation_id, "_Tin")[final_index]
    mdot = flow <= 1e-12 ? 0.0 : m_dot(flow, Tin)

    weights = solar_weights_v5(params[5], params[4], length(z))
    absorbed_solar = ETA_ABS_FIXED_v10 * irradiance_factor_v10(conditions[Io], params) *
                     conditions[Io] * A_frt
    qsolar = absorbed_solar .* weights
    qgas_receiver = gas_heat_gain_cells_v10(Tg, 1, length(z), mdot)
    qgas_rear = gas_heat_gain_cells_v10(Tg, length(z) + 1, length(z_rear), mdot)

    qside_receiver = [
        RECEIVER_TO_CAVITY_CONDUCTANCE_PER_LENGTH_v10 * dx * (Ts[i] - Tcavity)
        for i in eachindex(Ts)
    ]
    qcond_face = Float64[]
    qrad_face = Float64[]
    for i in 1:(length(z) - 1)
        kface = 2.0 * ks_f(Ts[i]) * ks_f(Ts[i + 1]) / (ks_f(Ts[i]) + ks_f(Ts[i + 1]))
        push!(qcond_face, kface * A_solid * (Ts[i] - Ts[i + 1]) / dx)
        zface = 0.5 * (z[i] + z[i + 1])
        krad = use_rosseland_radiation ?
               rosseland_conductivity_v10(
                   0.5 * (Ts[i] + Ts[i + 1]), params;
                   z=zface, profile=rosseland_profile,
               ) :
               0.0
        push!(qrad_face, krad * RADIATION_AREA_v10 * (Ts[i] - Ts[i + 1]) / dx)
    end

    qrear_cavity = [
        rear_tube_cavity_conductance_per_length_v10(z_rear[j], params) * dx_rear *
        (Ttube[j] - Tcavity)
        for j in eachindex(Ttube)
    ]
    qrear_flange = [
        rear_tube_flange_conductance_per_length_v10(z_rear[j], Ttube[j]) * dx_rear *
        (Ttube[j] - WATER_FLANGE_TEMPERATURE_v10)
        for j in eachindex(Ttube)
    ]

    front_loss = H_FRONT_FIXED_v10 * A_frt * (Ts[1] - ambient) +
                 EPS_FRONT_FIXED_v10 * sigma_sb * A_frt * (Ts[1]^4 - ambient^4)
    receiver_to_tube = RECEIVER_TO_TUBE_CONDUCTANCE_v10 * (Ts[end] - Ttube[1])
    receiver_to_cavity = sum(qside_receiver)
    tube_to_cavity = sum(qrear_cavity)
    flange_loss = sum(qrear_flange)
    cavity_ambient = cavity_loss_to_ambient_v10(Tcavity, ambient)

    beta_cells = [
        rosseland_beta_v10(params; z=zi, profile=rosseland_profile)
        for zi in z
    ]
    krad_cells = use_rosseland_radiation ?
        [rosseland_conductivity_v10(Ts[i], params; z=z[i], profile=rosseland_profile)
         for i in eachindex(Ts)] :
        zeros(length(Ts))
    nu_cells = [
        h[i] * Dh / kf_f(0.5 * (Ts[i] + Tg[i]))
        for i in eachindex(Ts)
    ]

    axial_rows = NamedTuple[]
    for i in eachindex(z)
        push!(axial_rows, (
            z_mm=z[i] * 1000.0,
            solid_K=Ts[i],
            gas_in_K=Tg[i],
            gas_out_K=Tg[i + 1],
            qsolar_W=qsolar[i],
            qgas_W=qgas_receiver[i],
            qside_cavity_W=qside_receiver[i],
            qcond_to_next_W=i < length(z) ? qcond_face[i] : NaN,
            qrad_to_next_W=i < length(z) ? qrad_face[i] : NaN,
            h_W_m2_K=h[i],
            Nu=nu_cells[i],
            beta_tr_1_m=beta_cells[i],
            ell_rad_m=1.0 / beta_cells[i],
            k_rad_W_m_K=krad_cells[i],
            tau_cell=beta_cells[i] * dx,
            tau_channel=beta_cells[i] * Dh,
        ))
    end

    summary = (
        simulation_id=simulation_id,
        phase=is_cooling ? "cooling" : "heating",
        absorbed_solar_W=absorbed_solar,
        front_loss_W=front_loss,
        gas_gain_receiver_W=sum(qgas_receiver),
        gas_gain_rear_W=sum(qgas_rear),
        receiver_to_cavity_W=receiver_to_cavity,
        receiver_to_tube_W=receiver_to_tube,
        tube_to_cavity_W=tube_to_cavity,
        flange_heat_loss_W=flange_loss,
        cavity_ambient_loss_W=cavity_ambient,
        receiver_storage_residual_W=absorbed_solar - front_loss -
                                    sum(qgas_receiver) - receiver_to_cavity -
                                    receiver_to_tube,
        rear_tube_storage_residual_W=receiver_to_tube - sum(qgas_rear) -
                                     tube_to_cavity - flange_loss,
        cavity_storage_residual_W=receiver_to_cavity + tube_to_cavity -
                                  cavity_ambient,
        T8_steady_error_K=model[end, 1] - experiment[end, 1],
        T9_pair_steady_error_K=model[end, 2] - experiment[end, 2],
        T10_pair_steady_error_K=model[end, 3] - experiment[end, 3],
        T3_steady_error_K=model[end, 4] - experiment[end, 4],
        T2_steady_error_K=model[end, 5] - experiment[end, 5],
        mean_Nu=mean(nu_cells),
        min_h_W_m2_K=minimum(h),
        max_h_W_m2_K=maximum(h),
        min_beta_tr_1_m=minimum(beta_cells),
        max_beta_tr_1_m=maximum(beta_cells),
        max_k_rad_W_m_K=maximum(krad_cells),
    )
    return summary, axial_rows
end

function energy_bar_plot_v10(summary, variant_name)
    labels = [
        "solar", "front loss", "gas rx", "gas rear", "rx-cavity",
        "rx-tube", "tube-cavity", "flange", "cavity loss",
    ]
    values = [
        summary.absorbed_solar_W,
        -summary.front_loss_W,
        -summary.gas_gain_receiver_W,
        -summary.gas_gain_rear_W,
        -summary.receiver_to_cavity_W,
        -summary.receiver_to_tube_W,
        -summary.tube_to_cavity_W,
        -summary.flange_heat_loss_W,
        -summary.cavity_ambient_loss_W,
    ]
    return bar(
        labels, values; legend=false, xrotation=35,
        ylabel="Power contribution (W)",
        title="v10 energy partition: $(summary.simulation_id) $variant_name",
    )
end

function axial_heat_plot_v10(axial_rows, simulation_id, variant_name)
    x = [row.z_mm for row in axial_rows]
    plot_object = plot(
        title="v10 axial heat terms: $simulation_id $variant_name",
        xlabel="Receiver length (mm)", ylabel="Cell / face power (W)",
        legend=:outerright,
    )
    plot!(plot_object, x, [row.qsolar_W for row in axial_rows];
          label="solar cell", lw=2, color=:orange)
    plot!(plot_object, x, [row.qgas_W for row in axial_rows];
          label="gas gain cell", lw=2, color=:green)
    plot!(plot_object, x, [row.qside_cavity_W for row in axial_rows];
          label="side to cavity", lw=2, color=:purple)
    plot!(plot_object, x, [row.qcond_to_next_W for row in axial_rows];
          label="solid cond to next", lw=2, color=:blue, linestyle=:dash)
    plot!(plot_object, x, [row.qrad_to_next_W for row in axial_rows];
          label="Rosseland to next", lw=2, color=:red, linestyle=:dot)
    return plot_object
end

begin # run diagnostics
    println("[v10-precheck] Running energy partition diagnostics.")
    summary_rows = NamedTuple[]
    summary_header = (
        "variant", "rosseland_profile", "use_rosseland_radiation",
        "simulation_id", "phase", "absorbed_solar_W", "front_loss_W",
        "gas_gain_receiver_W", "gas_gain_rear_W", "receiver_to_cavity_W",
        "receiver_to_tube_W", "tube_to_cavity_W", "flange_heat_loss_W",
        "cavity_ambient_loss_W", "receiver_storage_residual_W",
        "rear_tube_storage_residual_W", "cavity_storage_residual_W",
        "T8_steady_error_K", "T9_pair_steady_error_K",
        "T10_pair_steady_error_K", "T3_steady_error_K", "T2_steady_error_K",
        "mean_Nu", "min_h_W_m2_K", "max_h_W_m2_K", "min_beta_tr_1_m",
        "max_beta_tr_1_m", "max_k_rad_W_m_K",
    )
    axial_header = (
        "z_mm", "solid_K", "gas_in_K", "gas_out_K", "qsolar_W", "qgas_W",
        "qside_cavity_W", "qcond_to_next_W", "qrad_to_next_W", "h_W_m2_K",
        "Nu", "beta_tr_1_m", "ell_rad_m", "k_rad_W_m_K", "tau_cell",
        "tau_channel",
    )

    for variant in PRECHECK_VARIANTS_v10
        for case in PRECHECK_CASES_v10
            println("[v10-precheck] $(case.simulation_id) $(variant.name)")
            summary, axial_rows = final_energy_diagnostic_v10(
                variant.params, case.simulation_id;
                is_cooling=case.is_cooling,
                use_rosseland_radiation=variant.use_rad,
                rosseland_profile=variant.profile,
            )
            push!(summary_rows, merge(summary, (
                variant=variant.name,
                rosseland_profile=String(variant.profile),
                use_rosseland_radiation=variant.use_rad,
            )))

            stem = "$(case.simulation_id)_$(variant.name)"
            axial_path = joinpath(PRECHECK_AXIAL_DIR_v10, "axial_terms_$(stem)_1D_v10.csv")
            write_csv_v10(axial_path, axial_header, axial_rows)

            savefig(energy_bar_plot_v10(summary, variant.name),
                    joinpath(PRECHECK_PLOT_DIR_v10, "energy_partition_$(stem)_1D_v10.png"))
            savefig(axial_heat_plot_v10(axial_rows, case.simulation_id, variant.name),
                    joinpath(PRECHECK_PLOT_DIR_v10, "axial_heat_terms_$(stem)_1D_v10.png"))
        end
    end

    summary_path = joinpath(PRECHECK_OUTPUT_DIR_v10, "energy_partition_summary_1D_v10.csv")
    write_csv_v10(summary_path, summary_header, summary_rows)
    println("[v10-precheck] Saved summary: $summary_path")
    println("[v10-precheck] Complete.")
end
