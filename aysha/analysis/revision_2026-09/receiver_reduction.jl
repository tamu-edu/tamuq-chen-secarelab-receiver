"""
Corrected reduction pipeline for the SiC honeycomb volumetric receiver.

Run with the Aysha project:

    julia --project=. analysis/revision_2026-09/receiver_reduction.jl

The script follows the same linear scientific-workflow structure as
`0D_v1.jl` and `1D_v1.exp.All.jl`: libraries, fixed quantities, property
relations, data reduction functions, uncertainty analysis, and final output.
Running the file directly or with VS Code's "Execute File in REPL" starts the
workflow.  Set `RECEIVER_REDUCTION_INCLUDE_ONLY = true` before including it
only when another script needs the definitions without starting the workflow.
"""

begin # libraries
    using CSV
    using DataFrames
    using ForwardDiff
    using JSON3
    using LinearAlgebra
    using Optim
    using Printf
    using Random
    using Roots
    using Statistics
end

begin # campaign definition
    HEATING = [
        "E67" => ("Data_FPT0067_231125_161757", 456e3),
        "E68" => ("Data_FPT0068_231126_115725", 456e3),
        "E69" => ("Data_FPT0069_231126_140153", 456e3),
        "E70" => ("Data_FPT0070_231127_090339", 456e3),
        "E71" => ("Data_FPT0071_231128_102707", 456e3),
        "E72" => ("Data_FPT0072_231129_104140", 304e3),
        "E73" => ("Data_FPT0073_231129_132744", 304e3),
        "E74" => ("Data_FPT0074_231130_123228", 304e3),
        "E75" => ("Data_FPT0075_231201_162138", 304e3),
        "E76" => ("Data_FPT0076_231203_120521", 304e3),
        "E77" => ("Data_FPT0077_231203_161315", 256e3),
        "E78" => ("Data_FPT0078_231204_132252", 256e3),
        "E79" => ("Data_FPT0079_231204_172244", 256e3),
        "E80" => ("Data_FPT0080_231205_095122", 256e3),
        "E81" => ("Data_FPT0081_231205_135354", 256e3),
    ]

    # Replicates are archived but excluded from fitted quantities.
    REPLICATES = [
        "E82" => ("Data_FPT0082_231210_130825", 256e3),
        "E83" => ("Data_FPT0083_231211_122053", 256e3),
    ]

    COOLING = [
        "C69" => "Data_FPT0069-Cooling_231126_153148",
        "C80" => "Data_FPT0080-cooling_231205_112837",
        "C81" => "Data_FPT0081-cooling_231205_153409",
    ]

    COOL_PROVENANCE = Dict("C69" => "E69", "C80" => "E80", "C81" => "E81")
    FLUXES = [456, 304, 256]
end

begin # fixed geometry and instrument parameters
    # Receiver geometry
    W_CH = 1.5e-3
    T_WEB = 0.4e-3
    N_CH = 100
    L_REC = 0.137
    SIDE = 10 * (W_CH + T_WEB)
    A_FRT = SIDE^2
    A_CH = W_CH^2
    D_H = W_CH
    PER = 4 * W_CH
    POROSITY = N_CH * A_CH / A_FRT
    A_SOLID = A_FRT - N_CH * A_CH
    M_MONO = 0.040
    K_SIC = 40.0

    # Sensor locations
    Z_WALL = Dict("T8" => 0.011, "T12" => 0.058, "T11" => 0.107)
    Z_INT = Dict("T9" => 0.058, "T10" => 0.107)
    wall_sensor_order = ("T8", "T12", "T11")
    wall_positions = [Z_WALL[sensor] for sensor in wall_sensor_order]
    wall_boundaries = [0.0,
                       0.5 * (wall_positions[1] + wall_positions[2]),
                       0.5 * (wall_positions[2] + wall_positions[3]),
                       L_REC]
    WTS = Dict(sensor => (wall_boundaries[index + 1] - wall_boundaries[index]) / L_REC
               for (index, sensor) in enumerate(wall_sensor_order))

    # Flow and pressure instrumentation
    RHO_STD = 101325.0 / (287.05 * 294.25)
    DP_FS = 200.0
    DP_ACC = 0.001
    DT3_BAND = 25.0
    MFC_FS = 5.722
    MFC_A_FS = 0.0025
    MFC_B_REL = 0.025
    TOL_COVERAGE = 2.0
    RHO_MFC_CASES = (1.0, 0.0)

    # The ambient reference can be replaced from the command line with --tamb.
    TAMB_CHANNELS = ("T15", "T16")
    SENSORS = ("T2", "T3", "T8", "T9", "T10", "T11", "T12")
    COOL_SENS = ("T8", "T12", "T11", "T9", "T10", "T3")
    DEEP_SENS = ("T11", "T10", "T3")

    # Python used zero-based column indices. These are the equivalent Julia
    # one-based positions in the raw logger CSV files.
    COLS = (
        t=2, mfc=(7, 8, 9, 10), dp1=17, dp2=18,
        T1=35, T2=36, T3=37, T4=38, T5=39, T6=40, T7=41,
        T8=42, T9=43, T10=44, T11=45, T12=46, T15=49, T16=50,
    )
end

begin # air properties
    # CoolProp 8.0.0 dry-air properties at 1 atm are represented by Chebyshev
    # fits over 200-1600 K. The fits use the same 1 K grid as the Python code.
    # Maximum relative errors are 6.8e-7 (cp), 1.9e-7 (mu), and 6.9e-7 (k).
    P_AIR = 101325.0
    T_AIR_MIN = 200.0
    T_AIR_MAX = 1600.0
    T_AIR_CENTER = 900.0
    T_AIR_HALF_RANGE = 700.0

    CP_COEFFICIENTS = [
        1112.326521833738, 118.55461309246334, -3.0939714544722343,
        -11.51766499051778, 5.007902032742699, -0.4336469410171911,
        -0.5715402549958347, 0.3291775773182228, -0.05376833127680081,
        -0.04277682928194402, 0.042991940328635486, -0.022589009760191594,
        0.008516557171015625, -0.0025377659290189594, 0.0008259448213744829,
        -0.00046141563197419873, 0.00039221815770116344,
        -0.0002517802429630133, 0.00016907740337569307,
    ]
    MU_COEFFICIENTS = [
        3.8358695937796555e-5, 2.2119539563518034e-5, -2.1660765718609892e-6,
        5.461802807570149e-7, -1.4053274184219948e-7, 3.8105626169435703e-8,
        -1.06172172073623e-8, 2.9841508428715123e-9, -8.300841968246233e-10,
        2.2278412693467136e-10, -5.498643577165391e-11, 1.1156897759967311e-11,
        -8.553954229889346e-13, -6.723651882206778e-13, 6.573218038638565e-13,
    ]
    K_COEFFICIENTS = [
        0.060166159950534565, 0.03823890120257938, -0.0025345413742716484,
        0.0006606840752001102, -0.00016874205111542314, 4.4697238515949784e-5,
        -1.1921418643514947e-5, 3.095505654610283e-6, -7.367646494682498e-7,
        1.3574247993189645e-7, -7.715369142426102e-10, -1.7838980556699903e-8,
        1.3841373884409037e-8, -7.085291519288978e-9, 3.587160504666681e-9,
    ]

    function chebyshev_value(coefficients, x)
        b1 = zero(x)
        b2 = zero(x)
        for index in length(coefficients):-1:2
            b0 = 2x * b1 - b2 + coefficients[index]
            b2 = b1
            b1 = b0
        end
        return x * b1 - b2 + coefficients[1]
    end

    function air_property(T, coefficients, name)
        if T < T_AIR_MIN || T > T_AIR_MAX
            @warn "air $name requested outside 200-1600 K; endpoint used" temperature=T
        end
        Tb = clamp(Float64(T), T_AIR_MIN, T_AIR_MAX)
        x = (Tb - T_AIR_CENTER) / T_AIR_HALF_RANGE
        return chebyshev_value(coefficients, x)
    end

    cp_air(T::Real) = air_property(T, CP_COEFFICIENTS, "c_p")
    mu_air(T::Real) = air_property(T, MU_COEFFICIENTS, "mu")
    k_air(T::Real) = air_property(T, K_COEFFICIENTS, "k")
    cp_air(T::AbstractArray) = cp_air.(T)
    mu_air(T::AbstractArray) = mu_air.(T)
    k_air(T::AbstractArray) = k_air.(T)

    function h_gas(T_lo, T_hi; n=64)
        temperature = collect(range(Float64(T_lo), Float64(T_hi), length=n))
        heat_capacity = cp_air(temperature)
        return sum(0.5 .* (heat_capacity[1:end-1] .+ heat_capacity[2:end]) .*
                   diff(temperature))
    end
end

begin # common numerical functions
    function linear_fit(x_values, y_values)
        x = Float64.(x_values)
        y = Float64.(y_values)
        X = hcat(ones(length(x)), x)
        beta = X \ y
        residual = y - X * beta
        total = sum(abs2, y .- mean(y))
        r2 = total > 0 ? 1.0 - sum(abs2, residual) / total : 1.0
        dof = length(x) - 2
        stderr = if dof > 0
            covariance = (sum(abs2, residual) / dof) .* inv(X' * X)
            sqrt(max(covariance[2, 2], 0.0))
        else
            NaN
        end
        return (slope=beta[2], intercept=beta[1], r2=r2,
                stderr=stderr, residual=residual)
    end

    function power_law_fit(x, y)
        fit = linear_fit(log.(Float64.(x)), log.(Float64.(y)))
        return (prefactor=exp(fit.intercept), exponent=fit.slope,
                r2=fit.r2, stderr=fit.stderr)
    end

    function interpolate_clamped(x_grid, y_grid, x)
        x <= x_grid[1] && return y_grid[1]
        x >= x_grid[end] && return y_grid[end]
        index = searchsortedlast(x_grid, x)
        fraction = (x - x_grid[index]) / (x_grid[index + 1] - x_grid[index])
        return muladd(fraction, y_grid[index + 1] - y_grid[index], y_grid[index])
    end

    finite_values(values) = filter(isfinite, Float64.(collect(skipmissing(values))))
    finite_mean(values) = isempty(finite_values(values)) ? NaN : mean(finite_values(values))
    finite_std(values) = length(finite_values(values)) <= 1 ? 0.0 :
                         std(finite_values(values); corrected=false)

    function finite_quantile(values, probability)
        kept = sort(finite_values(values))
        isempty(kept) && return NaN
        position = 1.0 + (length(kept) - 1) * probability
        lower = floor(Int, position)
        upper = ceil(Int, position)
        lower == upper && return kept[lower]
        return kept[lower] + (position - lower) * (kept[upper] - kept[lower])
    end
end

begin # raw logger import
    numeric_column(data, index) = Float64.(coalesce.(data[!, index], NaN))

    function load_data(raw_dir, filename)
        path = joinpath(raw_dir, filename * ".csv")
        isfile(path) || error("raw logger file not found: $path")
        data = DataFrame(CSV.File(path; silencewarnings=true, strict=false))
        raw_time = numeric_column(data, COLS.t)
        keep = findall(isfinite, raw_time)
        time = raw_time[keep] .- raw_time[keep[1]]
        mfc = [numeric_column(data, index)[keep] for index in COLS.mfc]

        output = Dict{String,Any}(
            "t" => time,
            "flow" => reduce(+, mfc),
            "mfc" => mfc,
            "dp1" => numeric_column(data, COLS.dp1)[keep],
            "dp2" => numeric_column(data, COLS.dp2)[keep],
        )
        for sensor in ("T1", "T2", "T3", "T4", "T5", "T6", "T7",
                       "T8", "T9", "T10", "T11", "T12", "T15", "T16")
            output[sensor] = numeric_column(data, getproperty(COLS, Symbol(sensor)))[keep] .+ 273.15
        end
        output["Tamb"] = reduce(+, [output[channel] for channel in TAMB_CHANNELS]) ./
                         length(TAMB_CHANNELS)
        return output
    end

    function tail_mean(values, time; window=120.0)
        return finite_mean(values[time .>= time[end] - window])
    end

    function wall_temperature(data)
        return reduce(+, [WTS[sensor] .* data[sensor] for sensor in wall_sensor_order])
    end
end

begin # steady-state and dimensionless reduction
    function reduce_steady(raw_dir, runs)
        rows = NamedTuple[]
        for (ID, (filename, irradiance)) in runs
            data = load_data(raw_dir, filename)
            time = data["t"]
            flow = tail_mean(data["flow"], time)
            ambient = tail_mean(data["Tamb"], time)
            temperature = Dict(sensor => tail_mean(data[sensor], time) for sensor in SENSORS)
            mass_flow = RHO_STD * flow / 60000.0
            wall = sum(WTS[sensor] * temperature[sensor] for sensor in wall_sensor_order)
            gas_power = mass_flow * h_gas(ambient, temperature["T3"])
            controller_flow = [tail_mean(signal, time) for signal in data["mfc"]]
            shares = sum(controller_flow) > 0 ? controller_flow ./ sum(controller_flow) : fill(0.25, 4)

            push!(rows, (
                ID=ID, Io_kWm2=irradiance / 1e3, q_slpm=flow,
                Tamb=ambient, dur_s=time[end], mfc_rss=sqrt(sum(abs2, shares)),
                mfc_f1=shares[1], mfc_f2=shares[2],
                mfc_f3=shares[3], mfc_f4=shares[4],
                mdot_gs=mass_flow * 1e3, Tw_K=wall, Q_gas_W=gas_power,
                Q_nom_W=irradiance * A_FRT,
                dp1_mbar=tail_mean(data["dp1"], time),
                dp2_mbar=tail_mean(data["dp2"], time),
                T2_ss=temperature["T2"], T3_ss=temperature["T3"],
                T8_ss=temperature["T8"], T9_ss=temperature["T9"],
                T10_ss=temperature["T10"], T11_ss=temperature["T11"],
                T12_ss=temperature["T12"],
            ))
        end
        return DataFrame(rows)
    end

    function dimensionless(ss; dT3=0.0)
        output = copy(ss)
        T3 = output.T3_ss .+ dT3
        output.Tg_bar = 0.5 .* (output.Tamb .+ T3)
        output.mdot_ch = output.mdot_gs .* 1e-3 ./ N_CH
        output.Re = output.mdot_ch .* D_H ./ (A_CH .* mu_air(output.Tg_bar))
        output.Pr = cp_air(output.Tg_bar) .* mu_air(output.Tg_bar) ./ k_air(output.Tg_bar)
        output.Gz_L = D_H .* output.Re .* output.Pr ./ L_REC
        output.Pe_LD = output.Re .* output.Pr .* L_REC ./ D_H
        output.eps = (T3 .- output.Tamb) ./ (output.Tw_K .- output.Tamb)
        output.NTU = -log.(1.0 .- output.eps)
        output.h_app = output.NTU .* output.mdot_ch .* cp_air(output.Tg_bar) ./ (PER * L_REC)
        output.T_film = 0.5 .* (output.Tw_K .+ output.Tg_bar)
        output.Nu = output.h_app .* D_H ./ k_air(output.T_film)
        output.Bi = output.h_app .* T_WEB ./ (2 * K_SIC)
        output.N_rc = 4 * 5.670e-8 .* output.Tw_K.^3 .* D_H ./ k_air(output.Tw_K)
        output.Lam58 = (output.T12_ss .- output.T9_ss) ./ (output.T12_ss .- output.Tamb)
        output.Lam107 = (output.T11_ss .- output.T10_ss) ./ (output.T11_ss .- output.Tamb)
        output.I_vol = output.T12_ss .- output.T8_ss
        output.Q_gas_W = output.mdot_gs .* 1e-3 .*
                         [h_gas(a, b) for (a, b) in zip(output.Tamb, T3)]
        output.eta_nom = output.Q_gas_W ./ output.Q_nom_W
        return output
    end
end

begin # eigenvalue identification
    function eigen_cooling(data; sensors=COOL_SENS, thresh=5.0)
        time = data["t"]
        ambient = tail_mean(data["Tamb"], time)
        second_half = time .> time[end] / 2
        eigenvalues = Float64[]
        for sensor in sensors
            excess = data[sensor] .- ambient
            selected = second_half .& (excess .> thresh)
            count(selected) < 20 && continue
            fit = linear_fit(time[selected], log.(excess[selected]))
            push!(eigenvalues, -fit.slope)
        end
        return finite_mean(eigenvalues), finite_std(eigenvalues), length(eigenvalues)
    end

    function eigen_heating(data, sensors; u_lo=0.07, u_hi=0.45, r2_min=0.95)
        time = data["t"]
        eigenvalues = Float64[]
        for sensor in sensors
            temperature = data[sensor]
            steady = tail_mean(temperature, time)
            initial = mean(temperature[1:min(10, length(temperature))])
            steady - initial < 20 && continue
            deficit = (steady .- temperature) ./ (steady - initial)
            selected = (deficit .> u_lo) .& (deficit .< u_hi) .&
                       (steady .- temperature .> 0)
            count(selected) < 20 && continue
            fit = linear_fit(time[selected], log.((steady .- temperature)[selected]))
            fit.r2 >= r2_min && fit.slope < 0 && push!(eigenvalues, -fit.slope)
        end
        isempty(eigenvalues) && return NaN, NaN, 0
        return mean(eigenvalues), std(eigenvalues; corrected=false), length(eigenvalues)
    end

    function identify(exchange, eigenvalue)
        fit = linear_fit(exchange, eigenvalue)
        C_eff = 1 / fit.slope
        K_loss = fit.intercept / fit.slope
        return C_eff, K_loss, fit.r2, fit.stderr / fit.slope^2
    end
end

begin # inversion crossing, pressure drop, and delivered-power closure
    function crossings(group)
        order = sortperm(group.q_slpm)
        flow = Float64.(group.q_slpm[order])
        inversion = Float64.(group.I_vol[order])
        effectiveness = Float64.(group.eps[order])

        flow_fit = linear_fit(flow, inversion)
        q_global = -flow_fit.intercept / flow_fit.slope
        eps_fit = linear_fit(flow, effectiveness)
        eps_global = eps_fit.intercept + eps_fit.slope * q_global
        crossing_index = findfirst(!=(0.0), diff(sign.(inversion)))
        if isnothing(crossing_index)
            return Dict(
                "q_local" => NaN, "eps_local" => NaN,
                "q_global" => q_global, "eps_global" => eps_global,
                "bracketed" => false, "q_min" => minimum(flow), "q_max" => maximum(flow),
            )
        end

        index = crossing_index
        fraction = -inversion[index] / (inversion[index + 1] - inversion[index])
        q_local = flow[index] + fraction * (flow[index + 1] - flow[index])
        eps_local = effectiveness[index] +
                    fraction * (effectiveness[index + 1] - effectiveness[index])
        return Dict(
            "q_local" => q_local, "eps_local" => eps_local,
            "q_global" => q_global, "eps_global" => eps_global,
            "bracketed" => minimum(flow) <= q_global <= maximum(flow),
            "q_min" => minimum(flow), "q_max" => maximum(flow),
        )
    end

    function dp_laminar(mdot_gs, Tg)
        Po_square = 56.91
        mass_flux = (mdot_gs * 1e-3 / N_CH) / A_CH
        density = 101325.0 / (287.05 * Tg)
        return 0.5 * Po_square * mu_air(Tg) * L_REC * (mass_flux / density) /
               D_H^2 / 100.0
    end

    closure(data, K_loss) =
        (data.Q_gas_W .+ K_loss .* (data.Tw_K .- data.Tamb)) ./ data.Q_nom_W

    function reconciling_dT3(data, K_loss)
        output = Dict{String,Any}()
        for group in groupby(data, :Io_kWm2)
            required = group.Q_nom_W .- K_loss .* (group.Tw_K .- group.Tamb) .-
                       group.Q_gas_W
            shift = required ./ (group.mdot_gs .* 1e-3 .* cp_air(group.Tg_bar))
            output[string(round(Int, group.Io_kWm2[1]))] =
                [mean(shift), minimum(shift), maximum(shift)]
        end
        return output
    end
end

begin # grouped regression and profile-corrected transfer units
    function grouped_powerlaw(data, xcol, ycol; gcol=:Io_kWm2)
        x = log.(Float64.(data[!, xcol]))
        y = log.(Float64.(data[!, ycol]))
        labels = sort(unique(round.(Int, data[!, gcol])); rev=true)
        X = zeros(length(x), 1 + length(labels))
        X[:, 1] .= x
        for (column, label) in enumerate(labels)
            X[:, column + 1] .= round.(Int, data[!, gcol]) .== label
        end
        beta = X \ y
        residual = y - X * beta
        dof = length(y) - size(X, 2)
        stderr = sqrt(max((sum(abs2, residual) / dof) * inv(X' * X)[1, 1], 0.0))
        r2 = 1 - sum(abs2, residual) / sum(abs2, y .- mean(y))
        prefactors = Dict(string(label) => exp(beta[index + 1])
                          for (index, label) in enumerate(labels))
        return Dict("exponent" => beta[1], "stderr" => stderr,
                    "r2" => r2, "prefactors" => prefactors)
    end

    function wall_profile(row, zeta; rear="const", front="const")
        z = zeta * L_REC
        z_nodes = [Z_WALL["T8"], Z_WALL["T12"], Z_WALL["T11"]]
        temperatures = [row.T8_ss, row.T12_ss, row.T11_ss]
        wall = interpolate_clamped(z_nodes, temperatures, z)
        rear_slope = (row.T11_ss - row.T12_ss) / (z_nodes[3] - z_nodes[2])
        front_slope = (row.T12_ss - row.T8_ss) / (z_nodes[2] - z_nodes[1])
        rear == "linear" && z > z_nodes[3] &&
            return row.T11_ss + rear_slope * (z - z_nodes[3])
        rear == "half" && z > z_nodes[3] &&
            return row.T11_ss + 0.5 * rear_slope * (z - z_nodes[3])
        front == "linear" && z < z_nodes[1] &&
            return row.T8_ss - front_slope * (z_nodes[1] - z)
        return wall
    end

    function Tg_exit(N, row; rear="const", front="const", steps=400)
        Tg = Float64(row.Tamb)
        dz = 1 / steps
        for index in 0:(steps - 1)
            left = index * dz
            right = (index + 1) * dz
            k1 = N * (wall_profile(row, left; rear, front) - Tg)
            k2 = N * (wall_profile(row, right; rear, front) - (Tg + dz * k1))
            Tg += dz * (k1 + k2) / 2
        end
        return Tg
    end

    function solve_N(row; dT3=0.0, rear="const", front="const")
        residual(N) = Tg_exit(N, row; rear, front) - (row.T3_ss + dT3)
        return find_zero(residual, (1e-4, 120.0), Roots.Brent())
    end

    function ntu_profile_corrected(data; dT3_band=DT3_BAND)
        solve(offset=0.0) = [solve_N(row; dT3=offset) for row in eachrow(data)]
        corrected = solve()
        augmented = copy(data)
        augmented.NTU_corr = corrected
        fit = power_law_fit(augmented.Re, augmented.NTU_corr)
        band = sort([power_law_fit(data.Re, solve(offset)).exponent
                     for offset in (-dT3_band, dT3_band)])
        return Dict(
            "exponent" => fit.exponent, "stderr" => fit.stderr,
            "r2" => fit.r2,
            "grouped" => grouped_powerlaw(augmented, :Re, :NTU_corr),
            "ratio_to_NTU_app" => [minimum(corrected ./ data.NTU),
                                    maximum(corrected ./ data.NTU)],
            "exponent_T3_band" => band,
            "single_stream_requirement" => -1.0,
            "NTU_corr" => corrected,
        )
    end
end

begin # fixed-conductance falsification
    function Tg_exit_h(profile, nodes, row, scale; rear="const", front="const", steps=400)
        Tg = Float64(row.Tamb)
        dz = 1 / steps
        maximum_density = 1.5 / dz
        for index in 0:(steps - 1)
            left = index * dz
            right = (index + 1) * dz
            density_left = clamp(interpolate_clamped(nodes, profile, left) * scale,
                                 0.0, maximum_density)
            density_right = clamp(interpolate_clamped(nodes, profile, right) * scale,
                                  0.0, maximum_density)
            k1 = density_left * (wall_profile(row, left; rear, front) - Tg)
            k2 = density_right *
                 (wall_profile(row, right; rear, front) - (Tg + dz * k1))
            Tg += dz * (k1 + k2) / 2
        end
        return Tg
    end

    function fit_profile(data, indices, node_count, scale, rng;
                         n_starts=40, warm=nothing)
        nodes = collect(range(0.0, 1.0, length=node_count))
        rows = [data[index, :] for index in indices]
        base = mean(solve_N(row) / scale[index] for (row, index) in zip(rows, indices))

        function residual(log_profile)
            profile = exp.(log_profile)
            return [Tg_exit_h(profile, nodes, row, scale[index]) - row.T3_ss
                    for (row, index) in zip(rows, indices)]
        end
        objective(log_profile) = sum(abs2, residual(log_profile))
        gradient!(storage, log_profile) = ForwardDiff.gradient!(storage, objective, log_profile)

        starts = [log.(fill(base, node_count))]
        if !isnothing(warm)
            warm_nodes, warm_values = warm
            values = [interpolate_clamped(warm_nodes, warm_values, node) for node in nodes]
            push!(starts, log.(max.(values, base * 1e-9)))
        end
        append!(starts, [log.(base .* (0.1 .+ 9.9 .* rand(rng, node_count)))
                         for _ in 1:n_starts])

        best_profile = nothing
        best_value = Inf
        last_error = nothing
        options = Optim.Options(iterations=3000, f_reltol=1e-13,
                                g_tol=1e-10, show_trace=false)
        for initial in starts
            try
                solution = Optim.optimize(objective, gradient!, initial,
                                          Optim.LBFGS(), options)
                if Optim.minimum(solution) < best_value
                    best_value = Optim.minimum(solution)
                    best_profile = exp.(Optim.minimizer(solution))
                end
            catch error
                last_error = error
            end
        end
        isnothing(best_profile) && error("conductance-profile fit failed: " *
                                         sprint(showerror, last_error))
        return best_profile, residual(log.(best_profile))
    end

    function fixed_profile_test(data; n_starts=40, seed=20260904)
        rng = MersenneTwister(seed)
        heat_capacity = cp_air(data.Tg_bar)
        conductivity = k_air(data.Tg_bar)
        h_scale = PER * L_REC ./ (data.mdot_ch .* heat_capacity)
        nu_scale = h_scale .* conductivity ./ D_H
        log_reynolds = log.(data.Re)
        all_indices = collect(1:nrow(data))
        output = Dict(
            "shared_h" => Dict{String,Any}(),
            "per_flux" => Dict{String,Any}(),
            "shared_Nu" => Dict{String,Any}(),
            "run_ids" => String.(data.ID),
        )

        warm = nothing
        warm_five = nothing
        rms_values = Float64[]
        for node_count in (2, 3, 5, 7)
            profile, residual = fit_profile(data, all_indices, node_count, h_scale, rng;
                                            n_starts, warm)
            nodes = collect(range(0.0, 1.0, length=node_count))
            warm = (nodes, profile)
            node_count == 5 && (warm_five = warm)
            fit = linear_fit(log_reynolds, residual)
            rms = sqrt(mean(abs2, residual))
            push!(rms_values, rms)
            output["shared_h"][string(node_count)] = Dict(
                "nodes" => profile, "rms_K" => rms,
                "max_abs_K" => maximum(abs, residual),
                "r_lnRe" => cor(log_reynolds, residual),
                "slope_K_per_lnRe" => fit.slope,
                "residuals_K" => residual,
            )
        end
        any(diff(rms_values) .> 1e-6) &&
            @warn "shared-h residual is not monotone; increase --profile-starts" rms_values

        for group in groupby(data, :Io_kWm2)
            irradiance = round(Int, group.Io_kWm2[1])
            indices = findall(data.Io_kWm2 .== group.Io_kWm2[1])
            profile, residual = fit_profile(data, indices, 5, h_scale, rng;
                                            n_starts, warm=warm_five)
            log_group_reynolds = log.(Float64.(group.Re))
            fit = linear_fit(log_group_reynolds, residual)
            output["per_flux"][string(irradiance)] = Dict(
                "nodes" => profile, "rms_K" => sqrt(mean(abs2, residual)),
                "max_abs_K" => maximum(abs, residual),
                "r_lnRe" => cor(log_group_reynolds, residual),
                "slope_K_per_lnRe" => fit.slope, "n" => nrow(group),
                "run_ids" => String.(group.ID), "residuals_K" => residual,
            )
        end

        warm = nothing
        for node_count in (3, 5)
            profile, residual = fit_profile(data, all_indices, node_count, nu_scale, rng;
                                            n_starts, warm)
            nodes = collect(range(0.0, 1.0, length=node_count))
            warm = (nodes, profile)
            fit = linear_fit(log_reynolds, residual)
            output["shared_Nu"][string(node_count)] = Dict(
                "nodes" => profile, "rms_K" => sqrt(mean(abs2, residual)),
                "max_abs_K" => maximum(abs, residual),
                "r_lnRe" => cor(log_reynolds, residual),
                "slope_K_per_lnRe" => fit.slope,
                "residuals_K" => residual,
            )
        end
        return output
    end

    function wall_extrapolation_sensitivity(data)
        output = Dict{String,Any}()
        for rear in ("const", "half", "linear"), front in ("const", "linear")
            transfer = Float64[]
            reynolds = Float64[]
            infeasible = String[]
            for row in eachrow(data)
                if wall_profile(row, 1.0; rear, front) <= row.T3_ss
                    push!(infeasible, row.ID)
                    continue
                end
                push!(transfer, solve_N(row; rear, front))
                push!(reynolds, row.Re)
            end
            key = "$(rear)_$(front)"
            if length(transfer) == nrow(data)
                fit = power_law_fit(reynolds, transfer)
                output[key] = Dict("exponent" => fit.exponent,
                                   "stderr" => fit.stderr,
                                   "n" => length(transfer),
                                   "infeasible" => String[])
            else
                output[key] = Dict("exponent" => nothing, "stderr" => nothing,
                                   "n" => length(transfer), "infeasible" => infeasible,
                                   "note" => "exit wall below measured T3 for these runs")
            end
        end
        return output
    end
end

begin # Monte Carlo uncertainty propagation
    mfc_sigma(reading) = sqrt.((MFC_A_FS * MFC_FS)^2 .+ (MFC_B_REL .* reading).^2)

    function mfc_rel_perturbation(shares, total_flow, unit_draws)
        reading = shares .* total_flow
        flow_error = vec(sum(mfc_sigma(reading) .* reshape(unit_draws, 1, :), dims=2))
        return flow_error ./ total_flow
    end

    function monte_carlo(ss, eigenvalues; n=40, seed=20260902, rho=1.0)
        rng = MersenneTwister(seed)
        quantity_names = [
            "Nu_a", "Nu_b", "Nu_b_grouped", "NTU_corr_b",
            "eps_star_456", "eps_star_304", "eps_star_256",
            "Lam107_slope", "Lam107_int",
            "C_cool", "K_cool", "C_deep", "K_deep",
            "C_all", "K_all", "C_match", "K_match",
        ]
        samples = Dict(name => Float64[] for name in quantity_names)
        tolerance = Dict(sensor =>
            max.(1.5, 0.004 .* abs.(ss[!, Symbol(sensor, "_ss")] .- 273.15)) ./
            TOL_COVERAGE for sensor in SENSORS)
        steady_shares = Matrix{Float64}(ss[:, [:mfc_f1, :mfc_f2, :mfc_f3, :mfc_f4]])
        eigen_shares = Matrix{Float64}(eigenvalues[:, [:mfc_f1, :mfc_f2, :mfc_f3, :mfc_f4]])

        for realization in 1:n
            perturbed = copy(ss)
            sensor_draw = Dict{String,Float64}()
            for sensor in SENSORS
                sensor_draw[sensor] = randn(rng)
                column = Symbol(sensor, "_ss")
                perturbed[!, column] = perturbed[!, column] .+
                    sensor_draw[sensor] .* tolerance[sensor] .+
                    0.5 .* randn(rng, nrow(perturbed))
            end
            sensor_draw["Tamb"] = randn(rng)

            unit_draws = sqrt(rho) .* randn(rng) .+ sqrt(1 - rho) .* randn(rng, 4)
            perturbed.q_slpm = perturbed.q_slpm .* (1 .+
                mfc_rel_perturbation(steady_shares, ss.q_slpm, unit_draws))
            perturbed.mdot_gs = RHO_STD .* perturbed.q_slpm ./ 60000.0 .* 1e3
            perturbed.Tw_K = sum(WTS[s] .* perturbed[!, Symbol(s, "_ss")]
                                 for s in wall_sensor_order)
            perturbed.Q_gas_W = perturbed.mdot_gs .* 1e-3 .*
                [h_gas(a, b) for (a, b) in zip(perturbed.Tamb, perturbed.T3_ss)]

            groups = dimensionless(perturbed)
            groups.Re ./= 1 + 0.02 * randn(rng)
            groups.Nu ./= 1 + 0.01 * randn(rng)
            try
                fit = power_law_fit(groups.Re, groups.Nu)
                push!(samples["Nu_a"], fit.prefactor)
                push!(samples["Nu_b"], fit.exponent)
                push!(samples["Nu_b_grouped"], grouped_powerlaw(groups, :Re, :Nu)["exponent"])
            catch
                continue
            end
            try
                corrected = [solve_N(row) for row in eachrow(groups)]
                push!(samples["NTU_corr_b"], power_law_fit(groups.Re, corrected).exponent)
            catch
            end
            for group in groupby(groups, :Io_kWm2)
                irradiance = round(Int, group.Io_kWm2[1])
                push!(samples["eps_star_$irradiance"], crossings(DataFrame(group))["eps_local"])
            end

            labels = sort(unique(round.(Int, groups.Io_kWm2)); rev=true)
            X = zeros(nrow(groups), 1 + length(labels))
            X[:, 1] .= groups.Re
            for (column, label) in enumerate(labels)
                X[:, column + 1] .= round.(Int, groups.Io_kWm2) .== label
            end
            beta = X \ groups.Lam107
            push!(samples["Lam107_slope"], beta[1])
            push!(samples["Lam107_int"], mean(beta[2:end]))

            eig = copy(eigenvalues)
            deviation = [isfinite(value) ? value : 0.0 for value in eig.lam_sd]
            eig.lam = eig.lam .+ randn(rng, nrow(eig)) .* deviation
            eps_q = linear_fit(groups.q_slpm, groups.eps)
            flow_factor = 1 .+ mfc_rel_perturbation(eigen_shares, eig.q, unit_draws)
            q_e = eig.q .* flow_factor
            ambient_tolerance = max.(1.5, 0.004 .* abs.(eig.Tamb_e .- 273.15)) ./ TOL_COVERAGE
            outlet_tolerance = max.(1.5, 0.004 .* abs.(eig.T3_e .- 273.15)) ./ TOL_COVERAGE
            temperature_shift = sensor_draw["Tamb"] .* ambient_tolerance .+
                                sensor_draw["T3"] .* outlet_tolerance .+
                                0.5 .* randn(rng, nrow(eig)) .+
                                0.5 .* randn(rng, nrow(eig))
            effectiveness = eps_q.intercept .+ eps_q.slope .* q_e
            eig.x = effectiveness .* (RHO_STD .* q_e ./ 60000.0) .*
                    cp_air(0.5 .* (eig.Tamb_e .+ eig.T3_e) .+ 0.5 .* temperature_shift)

            for (tag, phases) in (("cool", ("cool",)),
                                  ("deep", ("heat",)),
                                  ("all", ("cool", "heat")))
                selected = in.(eig.phase, Ref(phases)) .& isfinite.(eig.lam) .& isfinite.(eig.x)
                count(selected) < 2 && continue
                C_eff, K_loss, _, _ = identify(eig.x[selected], eig.lam[selected])
                push!(samples["C_$tag"], C_eff)
                push!(samples["K_$tag"], K_loss)
            end

            cooling = findall(eig.phase .== "cool")
            matched_exchange = Float64[]
            matched_eigenvalue = Float64[]
            for index in cooling
                source = groups[findfirst(==(COOL_PROVENANCE[eig.ID[index]]), groups.ID), :]
                cp = cp_air(0.5 * (eig.Tamb_e[index] + eig.T3_e[index]) +
                            0.5 * temperature_shift[index])
                push!(matched_exchange, source.eps * RHO_STD * q_e[index] / 60000.0 * cp)
                push!(matched_eigenvalue, eig.lam[index])
            end
            C_match, K_match, _, _ = identify(matched_exchange, matched_eigenvalue)
            push!(samples["C_match"], C_match)
            push!(samples["K_match"], K_match)
        end

        rows = NamedTuple[]
        for quantity in quantity_names
            values = finite_values(samples[quantity])
            push!(rows, (
                quantity=quantity, value=finite_mean(values), sd=finite_std(values),
                ci_lo=finite_quantile(values, 0.025),
                ci_hi=finite_quantile(values, 0.975), n=length(values),
            ))
        end
        return DataFrame(rows)
    end
end

begin # table and JSON output
    function write_markdown(path, lines)
        open(path, "w") do io
            foreach(line -> println(io, line), lines)
        end
    end

    function write_tables(report, data, uncertainty, out_dir)
        order = sortperm(collect(zip(data.Io_kWm2, data.q_slpm)))
        measured = [
            "| Run | \$G_0\$ [kW m\$^{-2}\$] | \$q\$ [sL min\$^{-1}\$] | \$T_{\\rm amb}\$ [K] | \$\\bar T_w\$ [K] | \$T_3\$ [K] | \$T_{12}-T_8\$ [K] |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
        for index in order
            row = data[index, :]
            push!(measured, @sprintf("| %s | %.0f | %.2f | %.1f | %.0f | %.0f | %+.1f |",
                  row.ID, row.Io_kWm2, row.q_slpm, row.Tamb, row.Tw_K,
                  row.T3_ss, row.T12_ss - row.T8_ss))
        end
        write_markdown(joinpath(out_dir, "table_measured_envelope.md"), measured)

        corrected = report["ntu_profile"]["NTU_corr"]
        reduced = [
            "| Run | \$Re_{\\rm nom}\$ | \$Gz_L\$ | \$\\varepsilon\$ | \$NTU_{\\rm app}\$ | \$N_{\\rm prof}\$ | \$Nu_{\\rm app}\$ | \$\\Lambda_{58}\$ | \$\\Lambda_{107}\$ | \$\\eta_{\\rm nom}\$ |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for index in order
            row = data[index, :]
            push!(reduced, @sprintf("| %s | %.1f | %.3f | %.3f | %.3f | %.3f | %.4f | %.4f | %.4f | %.3f |",
                  row.ID, row.Re, row.Gz_L, row.eps, row.NTU, corrected[index],
                  row.Nu, row.Lam58, row.Lam107, row.eta_nom))
        end
        write_markdown(joinpath(out_dir, "table_reduced_envelope.md"), reduced)

        mc = Dict(row.quantity => row for row in eachrow(uncertainty))
        ci(name; scale=1.0, digits=3) = haskey(mc, name) ?
            "[$(round(mc[name].ci_lo * scale; digits)), $(round(mc[name].ci_hi * scale; digits))]" : "—"
        constants = [
            "| Constant | Value | s.d. | 95% interval | Unit | Notes |",
            "|---|---|---|---|---|---|",
            @sprintf("| \$Nu_{\\rm app}\$ prefactor \$a\$ | %.2f×10\$^{-4}\$ | %.2f×10\$^{-4}\$ | %s | ×10\$^{-4}\$ | 15 steady runs |",
                     mc["Nu_a"].value * 1e4, mc["Nu_a"].sd * 1e4,
                     ci("Nu_a"; scale=1e4, digits=2)),
            @sprintf("| \$Nu_{\\rm app}\$ exponent, grouped | %.3f | %.3f | %s | – | primary fit |",
                     report["nusselt"]["grouped"]["exponent"], mc["Nu_b_grouped"].sd,
                     ci("Nu_b_grouped")),
            @sprintf("| \$N_{\\rm prof}\$ exponent | %+.3f | %.3f | %s | – | profile integration |",
                     report["ntu_profile"]["exponent"], mc["NTU_corr_b"].sd,
                     ci("NTU_corr_b")),
        ]
        for (label, key, mc_key) in (
            ("cooling matched-\$\\varepsilon\$", "cooling_matched_eps", "C_match"),
            ("cooling pooled-\$\\varepsilon\$", "cooling", "C_cool"),
            ("heating deep probes", "heating_deep", "C_deep"),
            ("joint eigenvalues", "joint", "C_all"),
        )
            entry = report["identification"][key]
            push!(constants, @sprintf("| \$C_{\\rm eff}\$, %s | %.0f | %.0f | %s | J K\$^{-1}\$ | \$r^2\$=%.3f |",
                  label, entry["C_eff"], mc[mc_key].sd, ci(mc_key; digits=0), entry["r2"]))
        end
        write_markdown(joinpath(out_dir, "table_constants.md"), constants)
    end

    function write_supplementary_tables(report, data, out_dir)
        heating = [
            "| deficit window \$u\$ | \$C_{\\rm eff}\$ [J K\$^{-1}\$] | \$K_{\\rm loss}\$ [W K\$^{-1}\$] | \$r^2\$ |",
            "|---|---:|---:|---:|",
        ]
        for entry in report["heating_window_sensitivity"]
            push!(heating, @sprintf("| (%.2f, %.2f) | %.1f | %.4f | %.4f |",
                  entry["u_lo"], entry["u_hi"], entry["C_eff"], entry["K_loss"], entry["r2"]))
        end
        swing = report["sensor_selection_swing"]
        push!(heating, "")
        push!(heating, @sprintf("Sensor selection: all six probes %.1f J K\$^{-1}\$; deep probes %.1f J K\$^{-1}\$; ratio %.2f.",
              swing["C_all6"], swing["C_deep"], swing["ratio"]))
        write_markdown(joinpath(out_dir, "tableS3_heating_conditionality.md"), heating)

        wall = [
            "| Rear / front extrapolation | exponent | s.e. | runs | infeasible |",
            "|---|---:|---:|---:|---|",
        ]
        for key in ("const_const", "const_linear", "half_const", "half_linear",
                    "linear_const", "linear_linear")
            entry = report["wall_extrapolation"][key]
            exponent = isnothing(entry["exponent"]) ? "—" : @sprintf("%+.4f", entry["exponent"])
            stderr = isnothing(entry["stderr"]) ? "—" : @sprintf("%.4f", entry["stderr"])
            infeasible = isempty(entry["infeasible"]) ? "none" : join(entry["infeasible"], ", ")
            push!(wall, "| $(replace(key, "_" => " / ")) | $exponent | $stderr | $(entry["n"]) | $infeasible |")
        end
        append!(wall, ["", "| Family | nodes | RMS residual [K] | max abs [K] | \$r\$ vs \$\\ln Re\$ | slope [K per e-fold] |",
                       "|---|---:|---:|---:|---:|---:|"])
        fixed = report["fixed_profile_test"]
        for (family, label) in (("shared_h", "shared \$h(z)\$"),
                                ("shared_Nu", "shared \$Nu(z)\$"))
            for node_count in sort(parse.(Int, collect(keys(fixed[family]))))
                entry = fixed[family][string(node_count)]
                push!(wall, @sprintf("| %s | %d | %.1f | %.1f | %+.3f | %+.0f |",
                      label, node_count, entry["rms_K"], entry["max_abs_K"],
                      entry["r_lnRe"], entry["slope_K_per_lnRe"]))
            end
        end
        write_markdown(joinpath(out_dir, "tableS4_wall_and_falsification.md"), wall)

        groups = [
            "| Run | \$G_0\$ [kW m\$^{-2}\$] | \$q\$ [sL min\$^{-1}\$] | \$Re_{\\rm nom}\$ | \$Pr\$ | \$Gz_L\$ | \$x^*=1/Gz_L\$ | \$Re\\,Pr\\,L/D_h\$ | \$Bi\$ | \$N_{rc}\$ |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for row in eachrow(data)
            push!(groups, @sprintf("| %s | %.0f | %.2f | %.1f | %.4f | %.3f | %.2f | %.0f | %.2e | %.2f |",
                  row.ID, row.Io_kWm2, row.q_slpm, row.Re, row.Pr, row.Gz_L,
                  1 / row.Gz_L, row.Pe_LD, row.Bi, row.N_rc))
        end
        write_markdown(joinpath(out_dir, "tableS5_auxiliary_groups.md"), groups)
    end

    function json_safe(value)
        value isa AbstractFloat && !isfinite(value) && return nothing
        value isa AbstractDict && return Dict(string(key) => json_safe(item)
                                              for (key, item) in value)
        value isa NamedTuple && return Dict(string(key) => json_safe(getproperty(value, key))
                                            for key in keys(value))
        value isa AbstractArray && return [json_safe(item) for item in value]
        return value
    end

    function write_json(path, value)
        open(path, "w") do io
            JSON3.pretty(io, json_safe(value))
            println(io)
        end
    end
end

begin # linear reduction workflow and export
    RECEIVER_REDUCTION_RUN = !isdefined(@__MODULE__, :RECEIVER_REDUCTION_INCLUDE_ONLY) ||
                             !getfield(@__MODULE__, :RECEIVER_REDUCTION_INCLUDE_ONLY)

    if RECEIVER_REDUCTION_RUN
        let options = Dict(
                "raw" => normpath(joinpath(@__DIR__, "..", "RAW")),
                "out" => joinpath(@__DIR__, "outputs"),
                "nmc" => "4000",
                "profile-starts" => "40",
                "tamb" => "T15,T16",
            ), position = 1
            while position <= length(ARGS)
                argument = ARGS[position]
                startswith(argument, "--") || error("unknown argument: $argument")
                if occursin('=', argument)
                    key, value = split(argument[3:end], '='; limit=2)
                    options[key] = value
                else
                    position == length(ARGS) && error("missing value for $argument")
                    options[argument[3:end]] = ARGS[position + 1]
                    position += 1
                end
                position += 1
            end
            global raw_dir = abspath(options["raw"])
            global out_dir = abspath(options["out"])
            global n_mc = parse(Int, options["nmc"])
            global profile_starts = parse(Int, options["profile-starts"])
            global TAMB_CHANNELS = Tuple(strip.(split(options["tamb"], ',')))
        end

        println("Running receiver_reduction.jl")
        println("  raw data: ", raw_dir)
        println("  outputs:  ", out_dir)
        println("  Monte Carlo samples: ", n_mc)
        println("  profile random starts: ", profile_starts)
        flush(stdout)

        mkpath(out_dir)
        output(filename) = joinpath(out_dir, filename)
        report = Dict{String,Any}()

        # Steady measurements and dimensionless groups
        println("[1/6] Reducing steady heating and replicate measurements ...")
        flush(stdout)
        steady = reduce_steady(raw_dir, HEATING)
        steady_replicates = reduce_steady(raw_dir, REPLICATES)
        groups = dimensionless(steady)
        groups_replicates = dimensionless(steady_replicates)
        CSV.write(output("groups.csv"), groups)
        CSV.write(output("groups_replicates.csv"), groups_replicates)

        report["air_properties"] = Dict(
            "source" => "CoolProp 'Air' 8.0.0 Chebyshev surrogate",
            "coolprop_version" => "8.0.0",
            "pressure_Pa" => P_AIR,
            "range_K" => [T_AIR_MIN, T_AIR_MAX],
            "maximum_relative_fit_error" => Dict(
                "cp" => 6.8e-7, "viscosity" => 1.9e-7, "conductivity" => 6.9e-7),
        )
        report["geometry"] = Dict(
            "side_mm" => SIDE * 1e3,
            "A_frt_cm2" => A_FRT * 1e4,
            "porosity" => POROSITY,
            "A_solid_m2" => A_SOLID,
            "rho_eff" => M_MONO / (A_SOLID * L_REC),
            "wall_weights" => Dict(sensor => round(WTS[sensor]; digits=4)
                                   for sensor in wall_sensor_order),
        )
        envelope_columns = (:Re, :Pr, :Gz_L, :eps, :NTU, :Nu, :Bi,
                            :N_rc, :Lam58, :Lam107, :Pe_LD, :eta_nom)
        report["envelope"] = Dict(string(column) =>
            [minimum(groups[!, column]), maximum(groups[!, column])]
            for column in envelope_columns)

        # Apparent Nusselt law and transfer-unit structure
        nusselt = power_law_fit(groups.Re, groups.Nu)
        per_flux = Dict{String,Any}()
        for group in groupby(groups, :Io_kWm2)
            fit = power_law_fit(group.Re, group.Nu)
            per_flux[string(round(Int, group.Io_kWm2[1]))] = Dict(
                "exponent" => fit.exponent,
                "prefactor" => fit.prefactor,
                "r2" => fit.r2,
            )
        end
        Nu_fd = Dict("T" => 2.976, "H2" => 3.091, "H1" => 3.608)
        report["nusselt"] = Dict(
            "prefactor" => nusselt.prefactor,
            "exponent" => nusselt.exponent,
            "r2" => nusselt.r2,
            "stderr_exponent" => nusselt.stderr,
            "n" => nrow(groups),
            "fd_reference" => "H2",
            "Nu_fd" => Nu_fd,
            "ratio_to_fd" => [Nu_fd["H2"] / maximum(groups.Nu),
                               Nu_fd["H2"] / minimum(groups.Nu)],
            "ratio_to_fd_all" => Dict(key => [value / maximum(groups.Nu),
                                               value / minimum(groups.Nu)]
                                      for (key, value) in Nu_fd),
            "grouped" => grouped_powerlaw(groups, :Re, :Nu),
            "per_flux" => per_flux,
            "x_star_exit" => [1 / maximum(groups.Gz_L), 1 / minimum(groups.Gz_L)],
        )
        report["ntu_profile"] = ntu_profile_corrected(groups)
        ntu = power_law_fit(groups.Re, groups.NTU)
        fixed_Nu = linear_fit(log.(groups.Re),
            log.(k_air(groups.Tg_bar) ./ (groups.mdot_ch .* cp_air(groups.Tg_bar))))
        fixed_h = linear_fit(log.(groups.Re), -log.(groups.mdot_ch))
        report["ntu_structure"] = Dict(
            "exponent" => ntu.exponent,
            "r2" => ntu.r2,
            "fixed_Nu_requirement" => fixed_Nu.slope,
            "fixed_h_requirement" => fixed_h.slope,
            "stderr" => ntu.stderr,
            "single_stream_requirement" => -1.0,
            "gap_in_Re_power" => ntu.exponent + 1,
        )

        # Crossings and local thermal nonequilibrium
        println("[2/6] Fitting steady correlations and wall-profile transfer units ...")
        flush(stdout)
        report["crossings"] = Dict(string(round(Int, group.Io_kWm2[1])) =>
            crossings(DataFrame(group)) for group in groupby(groups, :Io_kWm2))
        ltne = Dict{String,Any}()
        for column in (:Lam58, :Lam107)
            pooled = linear_fit(groups.Re, groups[!, column])
            per_flux_ltne = Dict{String,Any}()
            for group in groupby(groups, :Io_kWm2)
                fit = linear_fit(group.Re, group[!, column])
                per_flux_ltne[string(round(Int, group.Io_kWm2[1]))] = Dict(
                    "slope" => fit.slope, "intercept" => fit.intercept,
                    "r2" => fit.r2,
                    "range" => [minimum(group[!, column]), maximum(group[!, column])],
                )
            end
            ltne[string(column)] = Dict(
                "pooled" => Dict("slope" => pooled.slope,
                                 "intercept" => pooled.intercept,
                                 "r2" => pooled.r2),
                "per_flux" => per_flux_ltne,
            )
        end
        report["ltne"] = ltne

        # Transient eigenvalues and identified assembly constants
        println("[3/6] Identifying heating and cooling eigenvalues ...")
        flush(stdout)
        eps_q = linear_fit(groups.q_slpm, groups.eps)
        eigenvalue_rows = NamedTuple[]
        for (ID, (filename, irradiance)) in HEATING
            data = load_data(raw_dir, filename)
            time = data["t"]
            flow = tail_mean(data["flow"], time)
            ambient = tail_mean(data["Tamb"], time)
            outlet = tail_mean(data["T3"], time)
            mdot = RHO_STD * flow / 60000.0
            effectiveness = eps_q.intercept + eps_q.slope * flow
            exchange = effectiveness * mdot * cp_air(0.5 * (ambient + outlet))
            controller_flow = [tail_mean(signal, time) for signal in data["mfc"]]
            shares = sum(controller_flow) > 0 ?
                     controller_flow ./ sum(controller_flow) : fill(0.25, 4)

            for (phase, sensors) in (("heat", DEEP_SENS), ("heat6", COOL_SENS))
                eigenvalue, deviation, count_sensors = eigen_heating(data, sensors)
                push!(eigenvalue_rows, (
                    ID=ID, phase=phase, q=flow, x=exchange,
                    lam=eigenvalue, lam_sd=deviation, n=count_sensors,
                    Tamb_e=ambient, T3_e=outlet,
                    mfc_f1=shares[1], mfc_f2=shares[2],
                    mfc_f3=shares[3], mfc_f4=shares[4],
                ))
            end
        end

        for (ID, filename) in COOLING
            data = load_data(raw_dir, filename)
            time = data["t"]
            flow = mean(data["flow"])
            ambient = tail_mean(data["Tamb"], time)
            outlet = tail_mean(data["T3"], time)
            mdot = RHO_STD * flow / 60000.0
            effectiveness = eps_q.intercept + eps_q.slope * flow
            exchange = effectiveness * mdot * cp_air(0.5 * (ambient + outlet))
            controller_flow = [mean(signal) for signal in data["mfc"]]
            shares = sum(controller_flow) > 0 ?
                     controller_flow ./ sum(controller_flow) : fill(0.25, 4)
            eigenvalue, deviation, count_sensors = eigen_cooling(data)
            push!(eigenvalue_rows, (
                ID=ID, phase="cool", q=flow, x=exchange,
                lam=eigenvalue, lam_sd=deviation, n=count_sensors,
                Tamb_e=ambient, T3_e=outlet,
                mfc_f1=shares[1], mfc_f2=shares[2],
                mfc_f3=shares[3], mfc_f4=shares[4],
            ))
        end

        eigenvalues = DataFrame(eigenvalue_rows)
        eigenvalues.x_matched = fill(NaN, nrow(eigenvalues))
        for index in findall(eigenvalues.phase .== "cool")
            source_id = COOL_PROVENANCE[eigenvalues.ID[index]]
            source = groups[findfirst(==(source_id), groups.ID), :]
            eigenvalues.x_matched[index] = source.eps * RHO_STD *
                eigenvalues.q[index] / 60000.0 *
                cp_air(0.5 * (eigenvalues.Tamb_e[index] + eigenvalues.T3_e[index]))
        end
        CSV.write(output("eigenvalues.csv"), eigenvalues)

        identification = Dict{String,Any}()
        identification_selections = [
            "cooling" => (eigenvalues.phase .== "cool"),
            "heating_deep" => (eigenvalues.phase .== "heat"),
            "heating_all6" => (eigenvalues.phase .== "heat6"),
            "joint" => in.(eigenvalues.phase, Ref(("cool", "heat"))),
        ]
        for (name, selected_rows) in identification_selections
            selected_rows .&= isfinite.(eigenvalues.lam) .& isfinite.(eigenvalues.x)
            C_current, K_current, r2_current, _ = identify(
                eigenvalues.x[selected_rows], eigenvalues.lam[selected_rows],
            )
            identification[name] = Dict(
                "C_eff" => C_current, "K_loss" => K_current, "r2" => r2_current,
                "n" => count(selected_rows), "dof" => count(selected_rows) - 2,
            )
        end
        matched_rows = (eigenvalues.phase .== "cool") .&
                       isfinite.(eigenvalues.lam) .& isfinite.(eigenvalues.x_matched)
        C_matched, K_matched, r2_matched, _ = identify(
            eigenvalues.x_matched[matched_rows], eigenvalues.lam[matched_rows],
        )
        identification["cooling_matched_eps"] = Dict(
            "C_eff" => C_matched, "K_loss" => K_matched, "r2" => r2_matched,
            "n" => count(matched_rows), "dof" => count(matched_rows) - 2,
        )
        report["identification"] = identification
        report["sensor_selection_swing"] = Dict(
            "C_deep" => identification["heating_deep"]["C_eff"],
            "C_all6" => identification["heating_all6"]["C_eff"],
            "ratio" => identification["heating_deep"]["C_eff"] /
                       identification["heating_all6"]["C_eff"],
        )
        heating_windows = Dict{String,Any}[]
        heating_eigenvalues = eigenvalues[eigenvalues.phase .== "heat", :]
        for (lower, upper) in ((0.05, 0.35), (0.07, 0.45), (0.10, 0.50),
                               (0.15, 0.60), (0.20, 0.70))
            estimates = Float64[]
            for (ID, (filename, irradiance)) in HEATING
                data = load_data(raw_dir, filename)
                eigenvalue, _, _ = eigen_heating(data, DEEP_SENS;
                                                  u_lo=lower, u_hi=upper)
                push!(estimates, eigenvalue)
            end
            selected_rows = isfinite.(estimates)
            C_window, K_window, r2_window, _ = identify(
                heating_eigenvalues.x[selected_rows], estimates[selected_rows],
            )
            push!(heating_windows, Dict(
                "u_lo" => lower, "u_hi" => upper,
                "C_eff" => C_window, "K_loss" => K_window, "r2" => r2_window,
            ))
        end
        report["heating_window_sensitivity"] = heating_windows
        report["monolith"] = Dict("C_monolith_600K" => M_MONO * 1050.0,
                                  "C_monolith_900K" => M_MONO * 1170.0)

        # Power closure and systematic sensitivities
        println("[4/6] Evaluating closure, probe, pressure, and sensitivity cases ...")
        flush(stdout)
        K_lo = identification["cooling_matched_eps"]["K_loss"]
        K_hi = identification["heating_deep"]["K_loss"]
        delivered_per_flux = Dict{String,Any}()
        for group in groupby(groups, :Io_kWm2)
            irradiance = group.Io_kWm2[1]
            f_lo = mean(closure(group, K_lo))
            f_hi = mean(closure(group, K_hi))
            candidates = (irradiance, irradiance * f_lo, irradiance * f_hi)
            delivered_per_flux[string(round(Int, irradiance))] = Dict(
                "f_Klo" => f_lo, "f_Khi" => f_hi,
                "G_closure_lo" => irradiance * f_lo,
                "G_closure_hi" => irradiance * f_hi,
                "G_nominal" => irradiance,
                "G_interval" => [minimum(candidates), maximum(candidates)],
                "eta_nom" => [minimum(group.eta_nom), maximum(group.eta_nom)],
            )
        end
        report["delivered_power"] = Dict(
            "K_bracket" => [K_lo, K_hi],
            "per_flux" => delivered_per_flux,
            "reconciling_dT3" => reconciling_dT3(groups, K_hi),
        )
        T3_cases = Dict{String,Any}()
        T3_offsets = (-DT3_BAND, 0.0, DT3_BAND)
        for offset in T3_offsets
            offset_groups = dimensionless(steady; dT3=offset)
            offset_nusselt = power_law_fit(offset_groups.Re, offset_groups.Nu)
            eps_star = Dict{String,Any}()
            for group in groupby(offset_groups, :Io_kWm2)
                eps_star[string(round(Int, group.Io_kWm2[1]))] =
                    crossings(DataFrame(group))["eps_local"]
            end
            T3_entry = Dict{String,Any}(
                "eps" => [minimum(offset_groups.eps), maximum(offset_groups.eps)],
                "Nu_prefactor" => offset_nusselt.prefactor,
                "Nu_exponent" => offset_nusselt.exponent,
                "NTU_exponent" => power_law_fit(
                    offset_groups.Re, offset_groups.NTU,
                ).exponent,
                "eps_star" => eps_star,
                "eta_nom" => [minimum(offset_groups.eta_nom),
                               maximum(offset_groups.eta_nom)],
            )

            cooling_rows = eigenvalues[eigenvalues.phase .== "cool", :]
            matched_exchange = Float64[]
            for row in eachrow(cooling_rows)
                source = offset_groups[
                    findfirst(==(COOL_PROVENANCE[row.ID]), offset_groups.ID), :,
                ]
                push!(matched_exchange,
                      source.eps * RHO_STD * row.q / 60000.0 *
                      cp_air(0.5 * (row.Tamb_e + row.T3_e + offset)))
            end
            C_match, K_match, r2_match, _ = identify(
                matched_exchange, cooling_rows.lam,
            )

            heating_rows = eigenvalues[
                (eigenvalues.phase .== "heat") .& isfinite.(eigenvalues.lam), :,
            ]
            offset_eps_q = linear_fit(offset_groups.q_slpm, offset_groups.eps)
            deep_exchange = [
                (offset_eps_q.intercept + offset_eps_q.slope * row.q) *
                RHO_STD * row.q / 60000.0 *
                cp_air(0.5 * (row.Tamb_e + row.T3_e + offset))
                for row in eachrow(heating_rows)
            ]
            C_deep, K_deep, r2_deep, _ = identify(deep_exchange, heating_rows.lam)
            T3_entry["identification"] = Dict(
                "C_match" => C_match, "K_match" => K_match,
                "r2_match" => r2_match, "C_deep" => C_deep,
                "K_deep" => K_deep, "r2_deep" => r2_deep,
            )
            T3_cases[@sprintf("%+.0f", offset)] = T3_entry
        end
        T3_cases["band"] = Dict{String,Any}()
        for quantity in ("C_match", "K_match", "C_deep", "K_deep")
            values = [T3_cases[@sprintf("%+.0f", offset)]["identification"][quantity]
                      for offset in T3_offsets]
            T3_cases["band"][quantity] = [minimum(values), maximum(values)]
        end
        report["T3_sensitivity"] = T3_cases
        reference_definitions = [
            "wall quadrature (T8,T12,T11)" =>
                (data -> sum(WTS[s] .* data[!, Symbol(s, "_ss")]
                             for s in wall_sensor_order)),
            "interior probes (T9,T10)" =>
                (data -> 0.5 .* (data.T9_ss .+ data.T10_ss)),
            "front wall only (T8)" => (data -> copy(data.T8_ss)),
            "mid wall only (T12)" => (data -> copy(data.T12_ss)),
            "rear wall only (T11)" => (data -> copy(data.T11_ss)),
            "wall+interior mean" => (data -> 0.5 .* (
                sum(WTS[s] .* data[!, Symbol(s, "_ss")] for s in wall_sensor_order) .+
                0.5 .* (data.T9_ss .+ data.T10_ss))),
        ]
        reference_rows = NamedTuple[]
        for (name, reference_temperature) in reference_definitions
            alternative = copy(steady)
            alternative.Tw_K = reference_temperature(alternative)
            alternative_groups = dimensionless(alternative)
            fit = power_law_fit(alternative_groups.Re, alternative_groups.Nu)
            alternative_crossings = Dict{Int,Float64}()
            for group in groupby(alternative_groups, :Io_kWm2)
                alternative_crossings[round(Int, group.Io_kWm2[1])] =
                    crossings(DataFrame(group))["eps_local"]
            end
            push!(reference_rows, (
                reference=name,
                eps_min=minimum(alternative_groups.eps),
                eps_max=maximum(alternative_groups.eps),
                NTU_min=minimum(alternative_groups.NTU),
                NTU_max=maximum(alternative_groups.NTU),
                Nu_min=minimum(alternative_groups.Nu),
                Nu_max=maximum(alternative_groups.Nu),
                Nu_prefactor=fit.prefactor,
                Nu_exponent=fit.exponent,
                eps_star_456=get(alternative_crossings, 456, NaN),
                eps_star_304=get(alternative_crossings, 304, NaN),
                eps_star_256=get(alternative_crossings, 256, NaN),
            ))
        end
        references = DataFrame(reference_rows)
        CSV.write(output("reference_sensitivity.csv"), references)
        report["reference_sensitivity"] =
            [Dict(string(name) => row[name] for name in names(references))
             for row in eachrow(references)]

        pressure = select(groups, :ID, :Io_kWm2, :q_slpm, :mdot_gs,
                          :Tg_bar, :dp1_mbar, :dp2_mbar)
        pressure.dp_pred_mbar = [dp_laminar(mass_flow, temperature)
                                 for (mass_flow, temperature) in
                                 zip(pressure.mdot_gs, pressure.Tg_bar)]
        pressure.ratio = pressure.dp1_mbar ./ pressure.dp_pred_mbar
        CSV.write(output("pressure_drop.csv"), pressure)
        report["pressure_drop"] = Dict(
            "resolution_mbar" => DP_FS * DP_ACC,
            "pred_range" => [minimum(pressure.dp_pred_mbar), maximum(pressure.dp_pred_mbar)],
            "meas_range" => [minimum(pressure.dp1_mbar), maximum(pressure.dp1_mbar)],
            "ratio_range" => [minimum(pressure.ratio), maximum(pressure.ratio)],
            "dp2_range" => [minimum(pressure.dp2_mbar), maximum(pressure.dp2_mbar)],
        )

        C_similarity = identification["cooling_matched_eps"]["C_eff"]
        K_similarity = identification["cooling_matched_eps"]["K_loss"]
        half_times = Dict("wall" => Float64[], "gas" => Float64[])
        for (ID, (filename, irradiance)) in HEATING
            data = load_data(raw_dir, filename)
            time = data["t"]
            flow = tail_mean(data["flow"], time)
            mdot = RHO_STD * flow / 60000.0
            ambient = tail_mean(data["Tamb"], time)
            row = groups[findfirst(==(ID), groups.ID), :]
            outlet = tail_mean(data["T3"], time)
            tau = C_similarity /
                  (row.eps * mdot * cp_air(0.5 * (ambient + outlet)) + K_similarity)
            for (name, signal) in (("wall", wall_temperature(data)), ("gas", data["T3"]))
                normalized = (signal .- signal[1]) ./
                             (tail_mean(signal, time) - signal[1])
                half_index = findfirst(>=(0.5), normalized)
                isnothing(half_index) || push!(half_times[name], time[half_index] / tau)
            end
        end
        report["similarity"] = Dict(
            name => Dict("t_half_mean" => mean(values),
                         "cv_percent" => 100 * std(values; corrected=false) / mean(values))
            for (name, values) in half_times
        )

        # Instrument uncertainty, fixed-profile falsification, and tables
        println("[5/6] Propagating instrument uncertainty with $n_mc samples ...")
        flush(stdout)
        uncertainty_cases = Dict{Float64,DataFrame}()
        selected_eigenvalues = eigenvalues[
            in.(eigenvalues.phase, Ref(("cool", "heat"))), :]
        for rho in RHO_MFC_CASES
            uncertainty_case = monte_carlo(steady, selected_eigenvalues; n=n_mc, rho)
            uncertainty_case.rho_mfc = fill(rho, nrow(uncertainty_case))
            uncertainty_cases[rho] = uncertainty_case
        end
        uncertainty = uncertainty_cases[1.0]
        CSV.write(output("uncertainty.csv"),
                  vcat([uncertainty_cases[rho] for rho in RHO_MFC_CASES]...))
        report["controller_correlation_bounds"] = Dict(
            "rho_$rho" => Dict(row.quantity =>
                Dict("value" => row.value, "sd" => row.sd)
                for row in eachrow(uncertainty_cases[rho]))
            for rho in RHO_MFC_CASES)
        report["monte_carlo"] =
            [Dict(string(name) => row[name] for name in names(uncertainty))
             for row in eachrow(uncertainty)]

        println("[6/6] Fitting fixed conductance profiles with $profile_starts random starts ...")
        flush(stdout)
        report["fixed_profile_test"] = fixed_profile_test(groups; n_starts=profile_starts)
        report["wall_extrapolation"] = wall_extrapolation_sensitivity(groups)
        write_tables(report, groups, uncertainty, out_dir)
        write_supplementary_tables(report, groups, out_dir)
        write_json(output("results.json"), report)
        println("Reduction complete. Outputs written to ", out_dir)
        flush(stdout)

        summary_keys = (
            "geometry", "nusselt", "ntu_structure", "identification",
            "sensor_selection_swing", "heating_window_sensitivity",
            "crossings", "ltne", "delivered_power", "T3_sensitivity",
            "pressure_drop", "similarity", "envelope",
        )
        JSON3.pretty(stdout, Dict(key => json_safe(report[key]) for key in summary_keys))
        println()
    end
end
