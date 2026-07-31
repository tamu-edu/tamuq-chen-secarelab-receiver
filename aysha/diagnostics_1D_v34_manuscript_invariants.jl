include("1D_v34.jl")

const INVARIANT_OUTPUT_DIR_1D = joinpath(@__DIR__, "summaries", "1D_v34")
mkpath(INVARIANT_OUTPUT_DIR_1D)

function linear_regression_1D(xs, ys)
    n = length(xs)
    n == length(ys) || throw(ArgumentError("x/y length mismatch"))
    n < 2 && return (intercept=NaN, slope=NaN, r2=NaN)
    xbar = mean(xs)
    ybar = mean(ys)
    sxx = sum((x - xbar)^2 for x in xs)
    sxy = sum((xs[i] - xbar) * (ys[i] - ybar) for i in eachindex(xs))
    slope = sxx <= eps(Float64) ? NaN : sxy / sxx
    intercept = ybar - slope * xbar
    residual = sum((ys[i] - (intercept + slope * xs[i]))^2 for i in eachindex(xs))
    total = sum((y - ybar)^2 for y in ys)
    r2 = total <= eps(Float64) ? NaN : 1.0 - residual / total
    return (intercept=intercept, slope=slope, r2=r2)
end

function reynolds_case_1D(flow_lpm, Tin)
    mdot = m_dot_standard_v34(flow_lpm)
    mu = mu_f_f(Tin)
    return mdot * Dh / (A_flow * mu)
end

function interp_profile_1D(z, values, position)
    position <= z[1] && return values[1]
    position >= z[end] && return values[end]
    i = searchsortedlast(z, position)
    fraction = (position - z[i]) / (z[i + 1] - z[i])
    return (1.0 - fraction) * values[i] + fraction * values[i + 1]
end

function tube_loss_totals_1D(result, p, time, irradiance)
    z_local = result.z_rear_tube .- L
    dx_rear = REAR_TUBE_LENGTH_v34 / length(z_local)
    Ttube = result.rear_tube_temperature[:, end]
    Tcavity = result.cavity_temperature[end]
    tube_to_cavity = 0.0
    tube_to_flange = 0.0
    for j in eachindex(z_local)
        G_cavity = rear_tube_cavity_conductance_per_length_v34(z_local[j], p)
        G_flange = rear_tube_flange_conductance_per_length_v34(
            z_local[j], Ttube[j], p, time, irradiance,
        )
        tube_to_cavity += G_cavity * dx_rear * (Ttube[j] - Tcavity)
        tube_to_flange += G_flange * dx_rear * (Ttube[j] - WATER_FLANGE_TEMPERATURE_v34)
    end
    return tube_to_cavity, tube_to_flange
end

function rear_rail_heat_totals_1D(result, p)
    Tcore = result.core_temperature[:, end]
    Tperim = result.perim_temperature[:, end]
    Trear = result.rear_reservoir_temperature[:, end]
    Tcavity = result.cavity_temperature[end]
    Ttube_in = result.rear_tube_temperature[1, end]
    f_core = rear_core_fraction_v34(p)
    G_receiver_rear = receiver_rear_conductance_v34(p)
    Qcore_rear = sum(
        G_receiver_rear * result.rear_contact_weights[i] * f_core *
        (Tcore[i] - Trear[i])
        for i in eachindex(Tcore)
    )
    Qperim_rear = sum(
        G_receiver_rear * result.rear_contact_weights[i] * (1.0 - f_core) *
        (Tperim[i] - Trear[i])
        for i in eachindex(Tperim)
    )
    Qrear_cavity = sum(
        rear_cavity_conductance_v34(p) * result.rear_contact_weights[i] *
        (Trear[i] - Tcavity)
        for i in eachindex(Trear)
    )
    Qrear_tube = rear_tube_conductance_v34(p) * (Trear[end] - Ttube_in)
    return Qcore_rear, Qperim_rear, Qrear_tube, Qrear_cavity
end

function invariant_row_1D(simulation_id, p=(@isdefined(fitted_params) ? fitted_params : pnew_v34); nodes=default_nodes)
    outputs, result, experiment = solve_case_v34(p, simulation_id; nodes=nodes)
    conditions = simulation_conditions[simulation_id]
    data = measurements
    time = result.time[end]
    flow = conditions[qlpm]
    irradiance = conditions[Io]
    Tin = observation(data, simulation_id, "_Tin")[end]
    ambient = observation(data, simulation_id, "_Tamb")[end]
    reynolds = reynolds_case_1D(flow, Tin)

    Tcore = result.core_temperature[:, end]
    Tperim = result.perim_temperature[:, end]
    Tgas = result.gas_temperature[:, end]
    q_receiver_gas = sum(result.receiver_gas_heat[:, end])
    q_rear_gas = sum(result.rear_tube_gas_heat[:, end])

    side_wall_basis = mean((
        interp_profile_1D(result.z_solid, Tperim, sensor_positions[:T8]),
        interp_profile_1D(result.z_solid, Tperim, sensor_positions[:T9]),
        interp_profile_1D(result.z_solid, Tperim, sensor_positions[:T10]),
    ))
    gas_bulk_basis = 0.5 * (Tin + Tgas[end])
    deltaT_app = side_wall_basis - gas_bulk_basis
    h_app = q_receiver_gas / (P_exchange * L * max(deltaT_app, 1.0))
    Tfilm = 0.5 * (side_wall_basis + gas_bulk_basis)
    Nu_app = h_app * Dh / max(kf_f(Tfilm), eps(Float64))

    z_peak_core = result.z_solid[argmax(Tcore)]
    z_peak_perim = result.z_solid[argmax(Tperim)]
    epsilon_mean = mean(result.receiver_effectiveness[:, end])
    epsilon_outlet = result.receiver_effectiveness[end, end]

    z107 = sensor_positions[:T10]
    Tcore_107 = interp_profile_1D(result.z_solid, Tcore, z107)
    Tperim_107 = interp_profile_1D(result.z_solid, Tperim, z107)
    lambda_107_raw_K = Tcore_107 - Tperim_107
    lambda_107 = lambda_107_raw_K / max(Tcore_107 - Tin, 1.0)

    Tfront = Tperim[1]
    q_front = H_FRONT_FIXED_v34 * A_frt * (Tfront - ambient) +
              EPS_FRONT_FIXED_v34 * sigma_sb * A_frt * (Tfront^4 - ambient^4)
    q_cavity_ambient = cavity_loss_to_ambient_v34(result.cavity_temperature[end], ambient)
    q_tube_cavity, q_tube_flange = tube_loss_totals_1D(result, p, time, irradiance)
    q_core_rear, q_perim_rear, q_rear_tube, q_rear_cavity = rear_rail_heat_totals_1D(result, p)
    q_loss_external = q_front + q_cavity_ambient + q_tube_flange
    deltaT_loss = side_wall_basis - ambient
    K_loss = q_loss_external / max(deltaT_loss, 1.0)

    return (
        case=simulation_id,
        flow_lpm=flow,
        reynolds=reynolds,
        irradiance=irradiance,
        power_scale=absorbed_power_scale_v34(irradiance, p),
        q_receiver_gas_W=q_receiver_gas,
        q_rear_gas_W=q_rear_gas,
        h_app_W_m2_K=h_app,
        Nu_app=Nu_app,
        deltaT_app_K=deltaT_app,
        epsilon_mean=epsilon_mean,
        epsilon_outlet=epsilon_outlet,
        z_peak_core_mm=1000.0 * z_peak_core,
        z_peak_perim_mm=1000.0 * z_peak_perim,
        lambda_107=lambda_107,
        lambda_107_raw_K=lambda_107_raw_K,
        Tcore_107_K=Tcore_107,
        Tperim_107_K=Tperim_107,
        K_loss_W_K=K_loss,
        q_loss_external_W=q_loss_external,
        q_front_W=q_front,
        q_cavity_ambient_W=q_cavity_ambient,
        q_tube_cavity_W=q_tube_cavity,
        q_tube_flange_W=q_tube_flange,
        q_core_rear_W=q_core_rear,
        q_perim_rear_W=q_perim_rear,
        q_rear_tube_W=q_rear_tube,
        q_rear_cavity_W=q_rear_cavity,
        C_receiver_participating_J_K=participating_assembly_heat_capacity_v34(p, nodes),
        C_total_with_rear_J_K=participating_total_heat_capacity_v34(p, nodes),
        T3_error_K=outputs[end, 6] - experiment[end, 6],
        T2_error_K=outputs[end, 7] - experiment[end, 7],
    )
end

function write_rows_1D(path, rows)
    fields = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(fields, ','))
        for row in rows
            println(io, join((getproperty(row, field) for field in fields), ','))
        end
    end
    return path
end

function write_summary_1D(path, rows)
    re = [row.reynolds for row in rows]
    nu = [row.Nu_app for row in rows]
    valid = [isfinite(re[i]) && isfinite(nu[i]) && re[i] > 0.0 && nu[i] > 0.0 for i in eachindex(re)]
    log_fit = linear_regression_1D(log.(re[valid]), log.(nu[valid]))
    lambda_fit = linear_regression_1D(re, [row.lambda_107 for row in rows])
    K_values = [row.K_loss_W_K for row in rows]
    eps_values = [row.epsilon_mean for row in rows]
    open(path, "w") do io
        println(io, "metric,value,target,status")
        println(io, "Nu_app_prefactor,$(exp(log_fit.intercept)),3.1e-4,diagnostic")
        println(io, "Nu_app_Re_exponent,$(log_fit.slope),1.44,diagnostic")
        println(io, "Nu_app_log_r2,$(log_fit.r2),0.97,diagnostic")
        println(io, "Lambda107_intercept,$(lambda_fit.intercept),0.038,diagnostic")
        println(io, "Lambda107_slope,$(lambda_fit.slope),8.3e-4,diagnostic")
        println(io, "Lambda107_r2,$(lambda_fit.r2),linear trend,diagnostic")
        println(io, "K_loss_min,$(minimum(K_values)),0.10-0.16,diagnostic")
        println(io, "K_loss_max,$(maximum(K_values)),0.10-0.16,diagnostic")
        println(io, "epsilon_mean_min,$(minimum(eps_values)),epsilon*=0.66,diagnostic")
        println(io, "epsilon_mean_max,$(maximum(eps_values)),epsilon*=0.66,diagnostic")
        receiver_capacity = participating_assembly_heat_capacity_v34((@isdefined(fitted_params) ? fitted_params : pnew_v34))
        total_capacity = participating_total_heat_capacity_v34((@isdefined(fitted_params) ? fitted_params : pnew_v34))
        receiver_status = 278.0 <= receiver_capacity <= 324.0 ? "pass" : "below assembly target"
        total_status = 278.0 <= total_capacity <= 324.0 ? "pass" : "failed"
        println(io, "C_receiver_participating,$receiver_capacity,301 +/- 23,$receiver_status")
        println(io, "C_total_with_rear,$total_capacity,301 +/- 23,$total_status")
    end
    return path
end

rows = [invariant_row_1D(simulation_id, (@isdefined(fitted_params) ? fitted_params : pnew_v34)) for simulation_id in sim_key_heat]
write_rows_1D(joinpath(INVARIANT_OUTPUT_DIR_1D, "invariants_1D_v34.csv"), rows)
write_summary_1D(joinpath(INVARIANT_OUTPUT_DIR_1D, "invariant_summary_1D_v34.csv"), rows)

println("[diagnostics_1D_manuscript_invariants] wrote summaries/1D_v34/invariants_1D_v34.csv")
println("[diagnostics_1D_manuscript_invariants] wrote summaries/1D_v34/invariant_summary_1D_v34.csv")
