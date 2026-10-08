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
    using CoolProp_jll
    using DataFrames
    using ForwardDiff
    using GLM
    using Interpolations
    using JSON3
    using Optim
    using Printf
    using Random
    using Roots
    using Statistics
    using Trapz
end

begin # campaign definition 
    # Imposed aperture irradiance levels [W m^-2]. Define them once here so
    # every campaign entry and downstream analysis uses the same values.
    FLUXES = [456e3, 304e3, 256e3]

    # Nominal cooling-flow labels [standard L min^-1]. Reduction always uses
    # the measured sum of the four MFC signals; these values are plot labels.
    COOLING_FLOWS = Dict("C69" => 10.5, "C80" => 6.6, "C81" => 4.5)

    # Heating logger basename and imposed irradiance [W m^-2] for each run.
    HEATING = [
        "E67" => ("Data_FPT0067_231125_161757", FLUXES[1]),
        "E68" => ("Data_FPT0068_231126_115725", FLUXES[1]),
        "E69" => ("Data_FPT0069_231126_140153", FLUXES[1]),
        "E70" => ("Data_FPT0070_231127_090339", FLUXES[1]),
        "E71" => ("Data_FPT0071_231128_102707", FLUXES[1]),
        "E72" => ("Data_FPT0072_231129_104140", FLUXES[2]),
        "E73" => ("Data_FPT0073_231129_132744", FLUXES[2]),
        "E74" => ("Data_FPT0074_231130_123228", FLUXES[2]),
        "E75" => ("Data_FPT0075_231201_162138", FLUXES[2]),
        "E76" => ("Data_FPT0076_231203_120521", FLUXES[2]),
        "E77" => ("Data_FPT0077_231203_161315", FLUXES[3]),
        "E78" => ("Data_FPT0078_231204_132252", FLUXES[3]),
        "E79" => ("Data_FPT0079_231204_172244", FLUXES[3]),
        "E80" => ("Data_FPT0080_231205_095122", FLUXES[3]),
        "E81" => ("Data_FPT0081_231205_135354", FLUXES[3]),
    ]

    # Replicates are archived but excluded from fitted quantities.
    REPLICATES = [
        "E82" => ("Data_FPT0082_231210_130825", FLUXES[3]),
        "E83" => ("Data_FPT0083_231211_122053", FLUXES[3]),
    ]

    # Cooling run IDs and raw logger basenames. Their measured flow is read
    # from the logger; COOLING_FLOWS above only supplies human-readable labels.
    COOLING = [
        "C69" => "Data_FPT0069-Cooling_231126_153148",
        "C80" => "Data_FPT0080-cooling_231205_112837",
        "C81" => "Data_FPT0081-cooling_231205_153409",
    ]

    # Heating run that supplies the steady initial state for each cooling run.
    COOL_PROVENANCE = Dict("C69" => "E69", "C80" => "E80", "C81" => "E81")
end

begin # fixed geometry and instrument parameters
    # Receiver geometry
    W_CH = 1.5e-3                    # square-channel width [m]
    T_WEB = 0.4e-3                   # SiC web thickness [m]
    N_CH = 100                       # number of parallel channels [-]
    L_REC = 0.137                    # receiver axial length [m]
    SIDE = 10 * (W_CH + T_WEB)       # illuminated square side length [m]
    A_FRT = SIDE^2                   # illuminated frontal area [m^2]
    A_CH = W_CH^2                    # open area of one channel [m^2]
    D_H = W_CH                       # square-channel hydraulic diameter [m]
    PER = 4 * W_CH                   # wetted perimeter of one channel [m]
    POROSITY = N_CH * A_CH / A_FRT   # open frontal-area fraction [-]
    A_SOLID = A_FRT - N_CH * A_CH    # solid frontal area [m^2]
    M_MONO = 0.040                   # measured monolith mass [kg]
    K_SIC = 40.0                     # effective SiC conductivity [W m^-1 K^-1]
    D_CAVITY = 0.142                 # aluminium cavity bore [m]
    A_APPROACH = pi / 4 * D_CAVITY^2 # flow approach area ahead of the face [m^2]
    SIGMA_FACE = N_CH * A_CH / A_APPROACH # face contraction ratio [-]

    # Flow-through monolith pressure-drop constants. The Poiseuille number and
    # the incremental pressure-drop number are the square-duct values tabulated
    # by Shah and London; K_INF and C_SHAH are the constants of Shah's apparent
    # friction-factor correlation for the hydrodynamic entrance region. The two
    # face coefficients are the sharp-edged contraction and the full dissipation
    # of the exit velocity head, both at the sigma -> 0 limit of this geometry.
    PO_SQUARE = 56.91                # fully developed square-duct f_D*Re [-]
    K_INF = 1.55                     # square-duct K(infinity) [-]
    C_SHAH = 0.00029                 # square-duct constant of Shah (1978) [-]
    KC_FACE = 0.50                   # face contraction loss coefficient [-]
    KE_FACE = 1.00                   # face expansion loss coefficient [-]

    # Sensor locations
    Z_WALL = Dict("T8" => 0.011, "T12" => 0.058, "T11" => 0.107) # wall TC z [m]
    Z_INT = Dict("T9" => 0.058, "T10" => 0.107)                   # interior TC z [m]
    wall_sensor_order = ("T8", "T12", "T11")                    # front-to-rear order
    wall_positions = [Z_WALL[sensor] for sensor in wall_sensor_order] # z [m]

    # Midpoints between adjacent thermocouples are control-volume boundaries.
    # With the receiver faces as the end boundaries, WTS is the fraction of
    # the receiver length represented by each wall temperature measurement.
    wall_boundaries = [0.0,
                       0.5 * (wall_positions[1] + wall_positions[2]),
                       0.5 * (wall_positions[2] + wall_positions[3]),
                       L_REC] # axial control-volume boundaries [m]
    WTS = Dict(sensor => (wall_boundaries[index + 1] - wall_boundaries[index]) / L_REC
               for (index, sensor) in enumerate(wall_sensor_order)) # length weights [-]

    # Flow and pressure instrumentation
    P_STD = 101325.0                  # reference/property pressure [Pa]
    T_STD = 294.25                    # Aalborg reference temperature [K]
    R_AIR = 287.05                    # dry-air gas constant [J kg^-1 K^-1]
    RHO_STD = P_STD / (R_AIR * T_STD) # standard air density [kg m^-3]
    DP_FS = 200.0                     # pressure-transducer full scale [mbar]
    DP_ACC = 0.001                    # pressure accuracy fraction of full scale [-]
    DT3_BAND = 25.0                   # outlet-TC systematic uncertainty [K]
    MFC_FS = 5.722                    # each MFC full scale [standard L min^-1]
    MFC_A_FS = 0.0025                 # MFC accuracy fraction of full scale [-]
    MFC_B_REL = 0.025                 # MFC accuracy fraction of reading [-]
    TOL_COVERAGE = 2.0                # stated tolerance coverage factor [-]
    RHO_MFC_CASES = (1.0, 0.0)        # correlated/independent MFC error cases [-]

    # The ambient reference can be replaced from the command line with --tamb.
    TAMB_CHANNELS = ("T15", "T16") # inlet-reference thermocouples [-]
    SENSORS = ("T2", "T3", "T8", "T9", "T10", "T11", "T12") # reduced TCs [-]
    COOL_SENS = ("T8", "T12", "T11", "T9", "T10", "T3") # cooling-fit TCs [-]
    DEEP_SENS = ("T11", "T10", "T3") # deep heating-fit TCs [-]

    # columns/positions in the raw logger CSV files.
    COLS = (
        t=2, mfc=(7, 8, 9, 10), dp1=17, dp2=18,
        T1=35, T2=36, T3=37, T4=38, T5=39, T6=40, T7=41,
        T8=42, T9=43, T10=44, T11=45, T12=46, T15=49, T16=50,
    )
end

begin # air properties
    # CoolProp 8.0 supplies dry-air c_p, viscosity and conductivity at 1 atm
    # with the same EOS and transport correlations as the Python reduction.
    # The one-time 1 K table keeps the Monte Carlo fast while retaining less
    # than 1e-6 relative interpolation error over the 200-1600 K range.
    T_AIR_GRID = collect(200.0:1.0:1600.0) # property-table temperature [K]
    COOLPROP_VERSION = string(pkgversion(CoolProp_jll)) # binary package version

    # Evaluate one scalar property through CoolProp's official C interface.
    # `output` is a PropsSI key and the returned value is in the associated SI unit.
    function coolprop_air(output, temperature)
        value = ccall((:PropsSI, CoolProp_jll.libcoolprop), Cdouble,
                      (Cstring, Cstring, Cdouble, Cstring, Cdouble, Cstring),
                      output, "T", Float64(temperature), "P", P_STD, "Air")
        isfinite(value) && abs(value) < 1e100 ||
            error("CoolProp failed for $output at $temperature K")
        return value
    end

    CP_AIR_GRID = coolprop_air.("CPMASS", T_AIR_GRID)       # c_p [J kg^-1 K^-1]
    MU_AIR_GRID = coolprop_air.("VISCOSITY", T_AIR_GRID)   # viscosity [Pa s]
    K_AIR_GRID = coolprop_air.("CONDUCTIVITY", T_AIR_GRID) # conductivity [W m^-1 K^-1]
    CP_AIR = linear_interpolation(T_AIR_GRID, CP_AIR_GRID;
                                  extrapolation_bc=Interpolations.Flat())
    MU_AIR = linear_interpolation(T_AIR_GRID, MU_AIR_GRID;
                                  extrapolation_bc=Interpolations.Flat())
    K_AIR = linear_interpolation(T_AIR_GRID, K_AIR_GRID;
                                 extrapolation_bc=Interpolations.Flat())

    cp_air(T::Real) = CP_AIR(T)              # dry-air c_p(T) [J kg^-1 K^-1]
    mu_air(T::Real) = MU_AIR(T)              # dry-air dynamic viscosity [Pa s]
    k_air(T::Real) = K_AIR(T)                # dry-air conductivity [W m^-1 K^-1]
    cp_air(T::AbstractArray) = CP_AIR.(T)
    mu_air(T::AbstractArray) = MU_AIR.(T)
    k_air(T::AbstractArray) = K_AIR.(T)

    # Integrate c_p dT to obtain the dry-air specific enthalpy change [J kg^-1].
    # Trapz supplies the standard trapezoidal integration used by the Python code.
    function h_gas(T_lo, T_hi; n=64)
        temperature = collect(range(Float64(T_lo), Float64(T_hi), length=n))
        return trapz(temperature, cp_air(temperature))
    end
end

begin # common numerical functions
    # Fit y = intercept + slope*x by ordinary least squares. GLM supplies the
    # coefficients, residuals, R^2 and coefficient standard errors.
    function linear_fit(x_values, y_values)
        x = Float64.(x_values)
        y = Float64.(y_values)
        model = lm(hcat(ones(length(x)), x), y)
        coefficients = coef(model)
        errors = length(x) > 2 ? stderror(model) : [NaN, NaN]
        residual = residuals(model)
        total = sum(abs2, y .- mean(y))
        r_squared = total > 0 ? 1.0 - sum(abs2, residual) / total : 1.0
        return (slope=coefficients[2], intercept=coefficients[1], r2=r_squared,
                stderr=errors[2], residual=residual)
    end

    # Fit y = prefactor*x^exponent by applying linear regression in log space.
    # Returned uncertainty is the standard error of the fitted exponent.
    function power_law_fit(x, y)
        fit = linear_fit(log.(Float64.(x)), log.(Float64.(y)))
        return (prefactor=exp(fit.intercept), exponent=fit.slope,
                r2=fit.r2, stderr=fit.stderr)
    end

    # Linearly interpolate tabulated values and hold the nearest endpoint
    # outside the grid. This helper also supports ForwardDiff dual values.
    function interpolate_clamped(x_grid, y_grid, x)
        x <= x_grid[1] && return y_grid[1]
        x >= x_grid[end] && return y_grid[end]
        index = searchsortedlast(x_grid, x)
        fraction = (x - x_grid[index]) / (x_grid[index + 1] - x_grid[index])
        return muladd(fraction, y_grid[index + 1] - y_grid[index], y_grid[index])
    end

    # Strip missing and non-finite entries before summary statistics.
    finite_values(values) = filter(isfinite, Float64.(collect(skipmissing(values))))
    # Return the finite-only mean, or NaN when no usable value exists.
    finite_mean(values) = isempty(finite_values(values)) ? NaN : mean(finite_values(values))
    # Return the finite-only population standard deviation.
    finite_std(values) = length(finite_values(values)) <= 1 ? 0.0 :
                         std(finite_values(values); corrected=false)
    # Use Statistics.quantile on finite entries, or NaN for an empty sample.
    finite_quantile(values, probability) = isempty(finite_values(values)) ? NaN :
                                           quantile(finite_values(values), probability)
end

begin # raw logger import
    # Convert one DataFrame column to Float64 and map missing entries to NaN.
    numeric_column(data, index) = Float64.(coalesce.(data[!, index], NaN))

    # Read one logger CSV and return time, MFC, pressure and Kelvin-temperature
    # signals in a dictionary with consistent names used by the reduction.
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

    # Average the finite samples in the final `window` seconds of a run. This
    # defines the steady-state value used throughout the reduction.
    function tail_mean(values, time; window=120.0)
        return finite_mean(values[time .>= time[end] - window])
    end

    # Form the axial mean wall temperature at every time using the control-
    # volume length fractions WTS. Each thermocouple represents its local span.
    function wall_temperature(data)
        return reduce(+, [WTS[sensor] .* data[sensor] for sensor in wall_sensor_order])
    end
end

begin # steady-state and dimensionless reduction
    # Reduce each heating run to one steady-state row. Flow is converted from
    # standard L/min to mass flow with the MFC reference density.
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
            total_controller_flow = sum(controller_flow)
            total_controller_flow > 0 ||
                error("$ID has non-positive total MFC flow; controller shares are undefined")
            shares = controller_flow ./ total_controller_flow

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

    # Calculate dimensionless groups and apparent transfer coefficients from
    # steady data. `T3_offset` is zero nominally and nonzero only for the
    # declared systematic outlet-thermocouple sensitivity calculation.
    function dimensionless(ss; T3_offset=0.0)
        output = copy(ss)
        T3 = output.T3_ss .+ T3_offset
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
    # Estimate the cooling decay eigenvalue by fitting log(T-Tamb) versus time
    # for eligible sensors in the latter half of the cooling record.
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

    # Estimate heating eigenvalues from the exponential approach to steady
    # state within the declared normalized-deficit and R^2 acceptance window.
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

    # Identify effective heat capacity and loss conductance from
    # lambda=(exchange+K_loss)/C_eff using a straight-line fit.
    function identify(exchange, eigenvalue)
        fit = linear_fit(exchange, eigenvalue)
        C_eff = 1 / fit.slope
        K_loss = fit.intercept / fit.slope
        return C_eff, K_loss, fit.r2, fit.stderr / fit.slope^2
    end
end

begin # inversion crossing, pressure drop, and delivered-power closure
    # Locate the flow where the volumetric-inversion index crosses zero. Both
    # local bracket interpolation and a global linear estimator are returned.
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

    # Predict laminar square-duct pressure loss [mbar]. For Darcy Poiseuille
    # number Po=f_D*Re=56.91, Darcy-Weisbach reduces to
    # Delta p=(Po/2)*mu*L*u/Dh^2 with u=mass_flux/rho; division by 100 converts Pa.
    function dp_laminar(mdot_gs, Tg)
        mass_flux = (mdot_gs * 1e-3 / N_CH) / A_CH # channel mass flux [kg m^-2 s^-1]
        density = P_STD / (R_AIR * Tg)              # ideal-gas density [kg m^-3]
        return 0.5 * PO_SQUARE * mu_air(Tg) * L_REC * (mass_flux / density) /
               D_H^2 / 100.0
    end

    # Normalize the gas-plus-loss power by nominal incident aperture power.
    closure(data, K_loss) =
        (data.Q_gas_W .+ K_loss .* (data.Tw_K .- data.Tamb)) ./ data.Q_nom_W

    # Compute the uniform T3 shift [K] that would close the steady energy
    # balance at each flux for a specified loss conductance.
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
    # Fit one common log-log exponent with a separate intercept for every
    # irradiance group, then return the group-specific power-law prefactors.
    function grouped_powerlaw(data, xcol, ycol; gcol=:Io_kWm2)
        x = log.(Float64.(data[!, xcol]))
        y = log.(Float64.(data[!, ycol]))
        labels = sort(unique(round.(Int, data[!, gcol])); rev=true)
        X = zeros(length(x), 1 + length(labels))
        X[:, 1] .= x
        for (column, label) in enumerate(labels)
            X[:, column + 1] .= round.(Int, data[!, gcol]) .== label
        end
        model = lm(X, y)
        coefficients = coef(model)
        residual = residuals(model)
        prefactors = Dict(string(label) => exp(coefficients[index + 1])
                          for (index, label) in enumerate(labels))
        r_squared = 1.0 - sum(abs2, residual) / sum(abs2, y .- mean(y))
        return Dict("exponent" => coefficients[1], "stderr" => stderror(model)[1],
                    "r2" => r_squared, "prefactors" => prefactors)
    end

    # Interpolate the three measured wall temperatures along z/L. The front
    # and rear keywords select constant, half-slope or linear extrapolation.
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

    # Integrate dTg/d(z/L)=N[Tw(z/L)-Tg], where N is the total transfer-unit
    # count. Each step is Heun's predictor-corrector (explicit trapezoidal) rule.
    function Tg_exit(N, row; rear="const", front="const", steps=400)
        Tg = Float64(row.Tamb)
        dz = 1 / steps
        for index in 0:(steps - 1)
            left = index * dz
            right = (index + 1) * dz
            k1 = N * (wall_profile(row, left; rear, front) - Tg)
            k2 = N * (wall_profile(row, right; rear, front) - (Tg + dz * k1))
            # k1 predicts the right-end temperature; averaging k1 and k2
            # provides the second-order Heun update over this axial step.
            Tg += dz * (k1 + k2) / 2
        end
        return Tg
    end

    # Invert the axial gas-energy equation for N so its predicted exit
    # temperature matches measured T3 plus an optional systematic offset.
    function solve_N(row; T3_offset=0.0, rear="const", front="const")
        residual(N) = Tg_exit(N, row; rear, front) - (row.T3_ss + T3_offset)
        return find_zero(residual, (1e-4, 120.0), Roots.Brent())
    end

    # Recompute transfer units using the measured nonuniform wall profile and
    # propagate the declared T3 systematic band to the fitted Re exponent.
    function ntu_profile_corrected(data; dT3_band=DT3_BAND)
        # Solve all runs for one common outlet-temperature offset [K].
        solve(offset=0.0) = [solve_N(row; T3_offset=offset) for row in eachrow(data)]
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

begin # composite monolith pressure drop
    # Record the axial gas temperature of the same Heun march that Tg_exit
    # performs. The two are kept separate because Tg_exit is the residual of a
    # root find that also runs inside the Monte Carlo loop, where returning
    # arrays would allocate on every evaluation.
    function Tg_profile(N, row; rear="const", front="const", steps=400)
        z = collect(range(0.0, L_REC; length=steps + 1)) # axial stations [m]
        Tg = zeros(steps + 1)
        Tg[1] = Float64(row.Tamb)
        dz = 1 / steps
        for index in 1:steps
            left = (index - 1) * dz
            right = index * dz
            k1 = N * (wall_profile(row, left; rear, front) - Tg[index])
            k2 = N * (wall_profile(row, right; rear, front) - (Tg[index] + dz * k1))
            Tg[index + 1] = Tg[index] + dz * (k1 + k2) / 2
        end
        return z, Tg
    end

    # Apparent Fanning f*Re over the whole channel length, Shah (1978). The
    # fully developed Fanning value is PO_SQUARE/4 = 14.227 and x_plus is the
    # dimensionless hydrodynamic development length.
    function f_apparent_Re(Re)
        x_plus = L_REC / (D_H * Re)
        developed = PO_SQUARE / 4
        entrance = 3.44 / sqrt(x_plus)
        return entrance +
               (K_INF / (4 * x_plus) + developed - entrance) / (1 + C_SHAH * x_plus^-2)
    end

    # Pressure drop of a flow-through monolith in which every channel is open,
    # as the sum of face contraction, core friction, hydrodynamic entrance,
    # thermal acceleration and face expansion. Core friction is integrated
    # along the channel at constant mass flux G with temperature-dependent
    # viscosity: with u=G*R*Tg/P the Darcy-Weisbach integrand is proportional
    # to mu(Tg)*Tg, so the one-point evaluation of dp_laminar at Tg_bar is
    # replaced by the axial integral of the reconstructed gas temperature.
    function dp_monolith(row; rear="const", front="const", steps=400)
        N = solve_N(row; rear, front)
        z, Tg = Tg_profile(N, row; rear, front, steps)
        mass_flux = (row.mdot_gs * 1e-3 / N_CH) / A_CH # channel mass flux [kg m^-2 s^-1]
        rho_in = P_STD / (R_AIR * row.Tamb)            # inlet density [kg m^-3]
        rho_out = P_STD / (R_AIR * row.T3_ss)          # outlet density [kg m^-3]
        Re_in = mass_flux * D_H / mu_air(row.Tamb)     # inlet channel Reynolds number [-]

        friction = 0.5 * PO_SQUARE * mass_flux * R_AIR / (D_H^2 * P_STD) *
                   trapz(z, mu_air(Tg) .* Tg)
        entrance = K_INF * mass_flux^2 / (2 * rho_in)
        acceleration = mass_flux^2 * (1 / rho_out - 1 / rho_in)
        contraction = (1 - SIGMA_FACE^2 + KC_FACE) * mass_flux^2 / (2 * rho_in)
        expansion = -(1 - SIGMA_FACE^2 - KE_FACE) * mass_flux^2 / (2 * rho_out)

        return (
            NTU_dp=N,
            Re_in=Re_in,
            Tg_len_mean=trapz(z, Tg) / L_REC,
            f_app_ratio=f_apparent_Re(Re_in) / (PO_SQUARE / 4),
            dp_in_mbar=contraction / 100.0,
            dp_fric_mbar=friction / 100.0,
            dp_dev_mbar=entrance / 100.0,
            dp_mom_mbar=acceleration / 100.0,
            dp_out_mbar=expansion / 100.0,
            dp_mono_mbar=(contraction + friction + entrance + acceleration + expansion) /
                         100.0,
        )
    end
end

begin # fixed-conductance falsification
    # Integrate the gas equation for a spatially varying conductance profile.
    # `scale` converts profile values to local transfer-unit density.
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

    # Fit a positive piecewise-linear conductance profile to selected runs by
    # multi-start nonlinear least squares in log-profile coordinates.
    function fit_profile(data, indices, node_count, scale, rng;
                         n_starts=40, warm=nothing)
        nodes = collect(range(0.0, 1.0, length=node_count))
        rows = [data[index, :] for index in indices]
        base = mean(solve_N(row) / scale[index] for (row, index) in zip(rows, indices))

        # Return predicted-minus-measured exit temperatures [K] for one profile.
        function residual(log_profile)
            profile = exp.(log_profile)
            return [Tg_exit_h(profile, nodes, row, scale[index]) - row.T3_ss
                    for (row, index) in zip(rows, indices)]
        end
        # Scalar sum of squared exit-temperature residuals minimized by Optim.
        objective(log_profile) = sum(abs2, residual(log_profile))
        # Fill Optim's gradient buffer by automatic differentiation.
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

    # Test whether one flow-independent h(z) or Nu(z) profile can reproduce all
    # exits, and compare it with separate five-node profiles for each flux.
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

    # Repeat the profile-corrected NTU fit for plausible front/rear wall
    # extrapolations and record infeasible cases where Tw,exit <= Tg,exit.
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
    # Return one-sigma MFC error [standard L min^-1] from full-scale and
    # reading-proportional contributions combined in quadrature.
    mfc_sigma(reading) = sqrt.((MFC_A_FS * MFC_FS)^2 .+ (MFC_B_REL .* reading).^2)

    # Combine the four controller errors into relative total-flow errors for
    # all Monte Carlo draws while preserving their measured flow shares.
    function mfc_rel_perturbation(shares, total_flow, unit_draws)
        reading = shares .* total_flow
        flow_error = vec(sum(mfc_sigma(reading) .* reshape(unit_draws, 1, :), dims=2))
        return flow_error ./ total_flow
    end

    # Propagate thermocouple and MFC errors through steady, transient and
    # fitted quantities. `rho` selects correlated or independent MFC errors.
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
            coefficients = coef(lm(X, groups.Lam107))
            push!(samples["Lam107_slope"], coefficients[1])
            push!(samples["Lam107_int"], mean(coefficients[2:end]))

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
    # Write preformatted Markdown lines to one UTF-8 text file.
    function write_markdown(path, lines)
        open(path, "w") do io
            foreach(line -> println(io, line), lines)
        end
    end

    # Generate the compact main-text tables from the result dictionary and
    # Monte Carlo confidence intervals.
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
        # Monte Carlo 95% intervals, formatted with literal format strings so
        # the emitted table is byte-reproducible; round() left a trailing ".0"
        # on the integer-valued rows. Extended 2026-10-01 to carry every row
        # the manuscript's Table 4 reports, so the table is generated in full
        # rather than part-generated and part-transcribed.
        has(name) = haskey(mc, name)
        ci2(name; scale=1.0) = has(name) ?
            @sprintf("[%.2f, %.2f]", mc[name].ci_lo * scale, mc[name].ci_hi * scale) : "—"
        ci3(name) = has(name) ? @sprintf("[%.3f, %.3f]", mc[name].ci_lo, mc[name].ci_hi) : "—"
        ci0(name) = has(name) ? @sprintf("[%.0f, %.0f]", mc[name].ci_lo, mc[name].ci_hi) : "—"
        sdv(name; scale=1.0) = has(name) ? mc[name].sd * scale : NaN
        nus = report["nusselt"]
        npf = report["ntu_profile"]
        nst = report["ntu_structure"]
        idn = report["identification"]
        constants = [
            "| Constant | Value | s.d. | 95% interval | Unit | Notes |",
            "| ------------------------------------------------------- | -------------- | -------------- | -------------- | ---------- | ------------------------------------------------------------------------- |",
            @sprintf("| \$Nu_{\\rm app}\$ prefactor \$a\$ (pooled) | %.2f×10\$^{-4}\$ | %.2f×10\$^{-4}\$ | %s | ×10\$^{-4}\$ | 15 steady runs |",
                     mc["Nu_a"].value * 1e4, sdv("Nu_a"; scale=1e4), ci2("Nu_a"; scale=1e4)),
            @sprintf("| \$Nu_{\\rm app}\$ exponent, pooled | %.3f | %.3f | %s | – | instrumental MC; regression SE \$\\pm\$%.3f, \$r^2\$=%.3f |",
                     nus["exponent"], sdv("Nu_b"), ci3("Nu_b"),
                     nus["stderr_exponent"], nus["r2"]),
            @sprintf("| \$Nu_{\\rm app}\$ exponent, grouped (primary) | %.3f | %.3f | %s | – | instrumental MC; regression SE \$\\pm\$%.3f, \$r^2\$=%.4f |",
                     nus["grouped"]["exponent"], sdv("Nu_b_grouped"), ci3("Nu_b_grouped"),
                     nus["grouped"]["stderr"], nus["grouped"]["r2"]),
            @sprintf("| \$N_{\\rm prof}\$ exponent (primary) | %+.3f | %.3f | %s | – | instrumental MC; regression SE \$\\pm\$%.3f; fixed-\$Nu\$ requirement %.3f |",
                     npf["exponent"], sdv("NTU_corr_b"), ci3("NTU_corr_b"),
                     npf["stderr"], nst["fixed_Nu_requirement"]),
            @sprintf("| \$NTU_{\\rm app}\$ exponent (identity, superseded) | %+.3f | — | — | – | isothermal-wall identity; retained for comparison only |",
                     nst["exponent"]),
        ]
        # Row labels are abbreviated after the first of each family, matching the
        # author's Table 4 layout (2026-10-07); regenerating must not revert it.
        for (idx, flux) in enumerate(("456", "304", "256"))
            label = idx == 1 ? "Inversion marker \$\\varepsilon^*\$ @ $flux kW m\$^{-2}\$" :
                               "@ $flux kW m\$^{-2}\$"
            note = idx == 1 ? "operational marker under the adopted wall convention; see §5.1" : "same"
            push!(constants,
                  @sprintf("| %s | %.3f | %.3f | %s | – | %s |",
                           label, mc["eps_star_" * flux].value, sdv("eps_star_" * flux),
                           ci3("eps_star_" * flux), note))
        end
        push!(constants,
              @sprintf("| \$\\Lambda_{107}\$ slope [\$Re^{-1}\$] | %.2f×10\$^{-4}\$ | %.2f×10\$^{-4}\$ | %s | ×10\$^{-4}\$ | common slope, per-flux intercepts |",
                       mc["Lam107_slope"].value * 1e4, sdv("Lam107_slope"; scale=1e4),
                       ci2("Lam107_slope"; scale=1e4)))
        for (label, key, mc_key, note) in (
            ("cooling matched-\$\\varepsilon\$ (primary)", "cooling_matched_eps", "C_match",
             @sprintf("\$n\$=3, \$r^2\$=%.3f", idn["cooling_matched_eps"]["r2"])),
            ("cooling pooled-\$\\varepsilon\$", "cooling", "C_cool",
             @sprintf("\$r^2\$=%.3f", idn["cooling"]["r2"])),
            ("joint 18 eigenvalues", "joint", "C_all",
             @sprintf("\$r^2\$=%.3f", idn["joint"]["r2"])),
            ("heating deep probes", "heating_deep", "C_deep",
             @sprintf("%.0f J K\$^{-1}\$ if all six probes used", idn["heating_all6"]["C_eff"])),
        )
            prefix = label == "cooling matched-\$\\varepsilon\$ (primary)" ? "\$C_{\\rm eff}\$, " : ""
            push!(constants, @sprintf("| %s%s | %.0f | %.0f | %s | J K\$^{-1}\$ | %s |",
                  prefix, label, idn[key]["C_eff"], sdv(mc_key), ci0(mc_key), note))
        end
        for (label, key, mc_key, note) in (
            ("cooling matched-\$\\varepsilon\$ (primary)", "cooling_matched_eps", "K_match",
             "secant conductance"),
            ("heating deep probes", "heating_deep", "K_deep", "tangent conductance"),
        )
            prefix = label == "cooling matched-\$\\varepsilon\$ (primary)" ? "\$K_{\\rm loss}\$, " : ""
            push!(constants, @sprintf("| %s%s | %.3f | %.3f | %s | W K\$^{-1}\$ | %s |",
                  prefix, label, idn[key]["K_loss"], sdv(mc_key), ci3(mc_key), note))
        end
        push!(constants,
              @sprintf("| Monolith capacitance (measured mass) | %.1f – %.1f | — | — | J K\$^{-1}\$ | 40 g \$\\times\\,c_p\$(600–900 K) |",
                       report["monolith"]["C_monolith_600K"], report["monolith"]["C_monolith_900K"]))
        write_markdown(joinpath(out_dir, "table_constants.md"), constants)
    end

    # Generate supplementary tables documenting conditionality, wall-reference
    # sensitivity, fixed-profile tests and auxiliary dimensionless groups.
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

    # Recursively replace non-finite floating-point values with `nothing` so
    # the report contains valid JSON null values instead of NaN or infinity.
    function json_safe(value)
        value isa AbstractFloat && !isfinite(value) && return nothing
        value isa AbstractDict && return Dict(string(key) => json_safe(item)
                                              for (key, item) in value)
        value isa NamedTuple && return Dict(string(key) => json_safe(getproperty(value, key))
                                            for key in keys(value))
        value isa AbstractArray && return [json_safe(item) for item in value]
        return value
    end

    # Serialize a Julia object as indented, standards-compliant JSON.
    function write_json(path, value)
        open(path, "w") do io
            JSON3.pretty(io, json_safe(value))
            println(io)
        end
    end
end

begin # MAIN linear reduction workflow and export
    RECEIVER_REDUCTION_RUN = !isdefined(@__MODULE__, :RECEIVER_REDUCTION_INCLUDE_ONLY) ||
                             !getfield(@__MODULE__, :RECEIVER_REDUCTION_INCLUDE_ONLY)

    if RECEIVER_REDUCTION_RUN # do not run if module is simply included and not executed
        # default options and processing of arguments
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
        # Resolve an output filename inside the selected output directory.
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
            "source" => "CoolProp 'Air' through the official CoolProp_jll binary",
            "coolprop_version" => COOLPROP_VERSION,
            "pressure_Pa" => P_STD,
            "range_K" => [first(T_AIR_GRID), last(T_AIR_GRID)],
            "interpolation_step_K" => T_AIR_GRID[2] - T_AIR_GRID[1],
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
        # Standard thermally fully developed laminar square-duct values from
        # Shah & London (1978) for the conventional T, H2 and H1 wall boundary
        # conditions. H2 is the primary comparison used in this study.
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
            total_controller_flow = sum(controller_flow)
            total_controller_flow > 0 ||
                error("$ID has non-positive total MFC flow; controller shares are undefined")
            shares = controller_flow ./ total_controller_flow

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
            total_controller_flow = sum(controller_flow)
            total_controller_flow > 0 ||
                error("$ID has non-positive total MFC flow; controller shares are undefined")
            shares = controller_flow ./ total_controller_flow
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
            offset_groups = dimensionless(steady; T3_offset=offset)
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

        # The composite flow-through monolith prediction, resolved into its
        # five terms, is the reference the measurement is judged against;
        # dp_pred_mbar is retained so the two can be compared run by run.
        monolith = [dp_monolith(row) for row in eachrow(groups)]
        for name in keys(first(monolith))
            pressure[!, name] = [component[name] for component in monolith]
        end
        pressure.uplift = pressure.dp_mono_mbar ./ pressure.dp_pred_mbar
        pressure.ratio_mono = pressure.dp1_mbar ./ pressure.dp_mono_mbar
        pressure.friction_share = pressure.dp_fric_mbar ./ pressure.dp_mono_mbar
        CSV.write(output("pressure_drop.csv"), pressure)
        report["pressure_drop"] = Dict(
            "resolution_mbar" => DP_FS * DP_ACC,
            "pred_range" => [minimum(pressure.dp_pred_mbar), maximum(pressure.dp_pred_mbar)],
            "meas_range" => [minimum(pressure.dp1_mbar), maximum(pressure.dp1_mbar)],
            "ratio_range" => [minimum(pressure.ratio), maximum(pressure.ratio)],
            "dp2_range" => [minimum(pressure.dp2_mbar), maximum(pressure.dp2_mbar)],
            "mono_range" => [minimum(pressure.dp_mono_mbar),
                             maximum(pressure.dp_mono_mbar)],
            "mono_ratio_range" => [minimum(pressure.ratio_mono),
                                   maximum(pressure.ratio_mono)],
            "uplift_range" => [minimum(pressure.uplift), maximum(pressure.uplift)],
            "friction_share_range" => [minimum(pressure.friction_share),
                                       maximum(pressure.friction_share)],
            "f_app_ratio_range" => [minimum(pressure.f_app_ratio),
                                    maximum(pressure.f_app_ratio)],
            "Re_channel_range" => [minimum(pressure.Re_in), maximum(pressure.Re_in)],
            "Tg_len_mean_range" => [minimum(pressure.Tg_len_mean),
                                    maximum(pressure.Tg_len_mean)],
            "permeability_m2" => 2 * POROSITY * D_H^2 / PO_SQUARE,
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
