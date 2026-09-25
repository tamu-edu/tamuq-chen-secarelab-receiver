"""
Figures 3--7 for the SiC volumetric-receiver manuscript.

This script reads the archive written by receiver_reduction.jl, exports the
three plotting trace files, and writes the manuscript figures.  It deliberately
keeps the top-to-bottom, sectioned scientific-script structure used in
0D_v1.jl and 1D_v1.exp.All.jl.

Usage:
    julia --project=. analysis/revision_2026-09/make_figures.jl
    julia --project=. analysis/revision_2026-09/make_figures.jl \
        --raw analysis/RAW \
        --out analysis/revision_2026-09/outputs \
        --fig analysis/revision_2026-09/figures
"""

begin # libraries and shared reduction equations
    using CSV
    using DataFrames
    using JSON3
    using Printf
    using PythonPlot
    using Statistics

    # Import the reduction equations without launching the reduction. This flag
    # is needed because VS Code executes Julia files through `include`.
    let previous = isdefined(@__MODULE__, :RECEIVER_REDUCTION_INCLUDE_ONLY) ?
                   getfield(@__MODULE__, :RECEIVER_REDUCTION_INCLUDE_ONLY) : false
        global RECEIVER_REDUCTION_INCLUDE_ONLY = true
        try
            include(joinpath(@__DIR__, "receiver_reduction.jl"))
        finally
            global RECEIVER_REDUCTION_INCLUDE_ONLY = previous
        end
    end

    plt = PythonPlot
    mpl = PythonPlot.matplotlib
end

begin # plotting constants
    series_color = Dict(
        456 => "#1f4e79",
        304 => "#c1666b",
        256 => "#8f9779",
    )
    flux_levels = (456, 304, 256)

    cooling_color = Dict(
        "C69" => "#1f4e79",
        "C80" => "#c1666b",
        "C81" => "#8f9779",
    )
    cooling_flow = Dict("C69" => 10.5, "C80" => 6.6, "C81" => 4.5)

    reference_color = Dict(
        "front wall only" => "#7a3b2e",
        "mid wall only" => "#c1666b",
        "wall quadrature" => "#1f4e79",
        "interior probes" => "#8f9779",
        "rear wall only" => "#3f5c3a",
    )
    reference_order = (
        "front wall only",
        "mid wall only",
        "wall quadrature",
        "interior probes",
        "rear wall only",
    )
end

begin # manuscript plotting style
    function manuscript_style()
        mpl.rcParams.update(Dict(
            "font.size" => 8,
            "axes.titlesize" => 9,
            "axes.labelsize" => 8,
            "xtick.labelsize" => 7,
            "ytick.labelsize" => 7,
            "legend.fontsize" => 7,
            "axes.spines.top" => false,
            "axes.spines.right" => false,
            "axes.linewidth" => 0.8,
            "xtick.major.width" => 0.8,
            "ytick.major.width" => 0.8,
            "figure.dpi" => 110,
            "savefig.dpi" => 300,
            "lines.markeredgewidth" => 0,
            "font.family" => "sans-serif",
            "mathtext.default" => "regular",
        ))
    end

    function panel_letter(axis, letter)
        axis.text(
            -0.16, 1.06, letter;
            transform=axis.transAxes,
            fontsize=10,
            fontweight="bold",
            va="bottom",
            ha="left",
        )
    end

    function finish_figure(figure, path)
        figure.tight_layout()
        figure.savefig(path; bbox_inches="tight")
        plt.close()
        println("wrote ", basename(path))
    end

    function suppress_minor_labels(axis; x=true, y=false)
        x && axis.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
        y && axis.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    end
end

begin # plotting traces derived directly from the raw logger files
    function export_traces(raw_dir, out_dir)
        mkpath(out_dir)
        steady = reduce_steady(raw_dir, HEATING)
        groups = dimensionless(steady)
        eps_fit = linear_fit(groups.q_slpm, groups.eps)

        cooling_rows = NamedTuple[]
        for (ID, filename) in COOLING
            data = load_data(raw_dir, filename)
            time = data["t"]
            ambient = tail_mean(data["Tamb"], time)
            stride = max(1, length(time) ÷ 400)
            selected = 1:stride:length(time)
            for sensor in COOL_SENS
                theta = (data[sensor] .- ambient) ./ (data[sensor][1] - ambient)
                for index in selected
                    push!(cooling_rows, (
                        ID=ID,
                        sensor=sensor,
                        t=time[index],
                        theta=theta[index],
                    ))
                end
            end
        end
        CSV.write(joinpath(out_dir, "cooling_decays.csv"), DataFrame(cooling_rows))

        report = JSON3.read(
            read(joinpath(out_dir, "results.json"), String),
            Dict{String,Any},
        )
        cooling_fit = report["identification"]["cooling"]
        capacitance = Float64(cooling_fit["C_eff"])
        heat_loss = Float64(cooling_fit["K_loss"])

        master_rows = NamedTuple[]
        for (ID, (filename, irradiance)) in HEATING
            data = load_data(raw_dir, filename)
            time = data["t"]
            flow = tail_mean(data["flow"], time)
            ambient = tail_mean(data["Tamb"], time)
            mass_flow = RHO_STD * flow / 60000.0
            gas_temperature = 0.5 * (ambient + tail_mean(data["T3"], time))
            effectiveness = eps_fit.intercept + eps_fit.slope * flow
            timescale = capacitance /
                        (effectiveness * mass_flow * cp_air(gas_temperature) + heat_loss)
            stride = max(1, length(time) ÷ 400)
            selected = 1:stride:length(time)

            for (signal_name, signal) in
                (("wall", wall_temperature(data)), ("gas", data["T3"]))
                final_value = tail_mean(signal, time)
                theta = (signal .- signal[1]) ./ (final_value - signal[1])
                for index in selected
                    push!(master_rows, (
                        ID=ID,
                        Io=irradiance / 1e3,
                        signal=signal_name,
                        tstar=time[index] / timescale,
                        theta=theta[index],
                    ))
                end
            end
        end
        CSV.write(joinpath(out_dir, "master_curves.csv"), DataFrame(master_rows))

        reference_frames = DataFrame[]
        for reference in reference_order
            alternative = copy(steady)
            alternative.Tw_K = if reference == "wall quadrature"
                sum(WTS[sensor] .* alternative[!, Symbol(sensor, "_ss")]
                    for sensor in wall_sensor_order)
            elseif reference == "interior probes"
                0.5 .* (alternative.T9_ss .+ alternative.T10_ss)
            elseif reference == "front wall only"
                alternative.T8_ss
            elseif reference == "mid wall only"
                alternative.T12_ss
            else
                alternative.T11_ss
            end
            reduced = dimensionless(alternative)
            push!(reference_frames, DataFrame(
                reference=fill(reference, nrow(reduced)),
                Re=reduced.Re,
                Nu=reduced.Nu,
                eps=reduced.eps,
                NTU=reduced.NTU,
                Io=reduced.Io_kWm2,
            ))
        end
        CSV.write(
            joinpath(out_dir, "reference_points.csv"),
            reduce(vcat, reference_frames),
        )
        return nothing
    end
end

begin # author-supplied apparatus montage
    function convert_apparatus(fig_dir)
        source = joinpath(fig_dir, "Figure01_setup.tif")
        destination = joinpath(fig_dir, "fig1_apparatus.png")
        if !isfile(source)
            println("note: no Figure01_setup.tif; leaving ", destination, " as is")
            return nothing
        end

        image_library = PythonPlot.pyimport("PIL.Image")
        image_library.open(source).convert("RGB").save(destination; dpi=(300, 300))
        println("wrote fig1_apparatus.png from Figure01_setup.tif")
        return nothing
    end
end

begin # small data-selection helpers
    function flux_subset(groups, flux)
        selected = groups[isapprox.(groups.Io_kWm2, flux), :]
        sort!(selected, :q_slpm)
        return selected
    end

    function json_number(dictionary, key)
        value = get(dictionary, key, nothing)
        return isnothing(value) ? NaN : Float64(value)
    end

    function trace_files_are_stale(out_dir)
        result_path = joinpath(out_dir, "results.json")
        result_stamp = mtime(result_path)
        derived = ("cooling_decays.csv", "master_curves.csv", "reference_points.csv")
        return any(!isfile(joinpath(out_dir, filename)) ||
                   mtime(joinpath(out_dir, filename)) < result_stamp
                   for filename in derived)
    end
end

begin # Figure 3: steady temperature field
    function figure_steady_field(groups, fig_dir)
        figure, axes = plt.subplots(1, 3; figsize=(7.4, 2.9), sharey=true)
        for (axis, flux) in zip(axes, flux_levels)
            data = flux_subset(groups, flux)
            for (column, marker) in
                ((:T8_ss, "o"), (:T12_ss, "s"), (:T11_ss, "^"), (:T3_ss, "D"))
                fill_color = column == :T3_ss ? "white" : series_color[flux]
                axis.plot(
                    data.q_slpm, data[!, column] .- 273.15, string(marker, "-");
                    color=series_color[flux], ms=3.8, lw=1.1, mfc=fill_color,
                )
            end
            axis.set_title(string(flux, raw" kW m$^{-2}$"); loc="left")
            axis.set_xlabel(raw"$q$ [sL min$^{-1}$]")
        end
        axes[1].set_ylabel(raw"Steady temperature [$^\circ$C]")

        handles = [
            mpl.lines.Line2D(
                [], [];
                color="0.3", marker=marker, ls="-", ms=3.8,
                mfc=fill_color, lw=1.1,
            )
            for (marker, fill_color) in
                (("o", "0.3"), ("s", "0.3"), ("^", "0.3"), ("D", "white"))
        ]
        axes[1].legend(
            handles,
            ["wall 11 mm", "wall 58 mm", "wall 107 mm", "gas outlet"];
            frameon=false, loc="lower left", handlelength=1.5,
        )
        for (axis, letter) in zip(axes, ("a", "b", "c"))
            panel_letter(axis, letter)
        end
        finish_figure(figure, joinpath(fig_dir, "fig3_steady_field.png"))
    end
end

begin # Figure 4: assembly-scale heat-transfer limitation
    function figure_assembly_limitation(report, groups, fig_dir)
        apparent = report["nusselt"]
        prefactor = Float64(apparent["prefactor"])
        exponent = Float64(apparent["exponent"])
        reynolds_line = collect(range(21.0, 102.0; length=60))

        figure, axes = plt.subplots(1, 2; figsize=(7.4, 3.1))
        Nu_axis, NTU_axis = axes

        for flux in flux_levels
            data = flux_subset(groups, flux)
            Nu_axis.plot(
                data.Re, data.Nu, "o";
                color=series_color[flux], ms=4.5,
                label=string(flux, raw" kW m$^{-2}$"),
            )
        end

        grouped = apparent["grouped"]
        grouped_exponent = Float64(grouped["exponent"])
        for flux in flux_levels
            data = flux_subset(groups, flux)
            local_reynolds = collect(range(
                minimum(data.Re) * 0.9,
                maximum(data.Re) * 1.1;
                length=40,
            ))
            grouped_prefactor = Float64(grouped["prefactors"][string(flux)])
            Nu_axis.plot(
                local_reynolds,
                grouped_prefactor .* local_reynolds .^ grouped_exponent,
                "-";
                color=series_color[flux], lw=1.5, zorder=2,
            )
        end
        Nu_axis.plot(
            reynolds_line, prefactor .* reynolds_line .^ exponent, "--";
            color="0.35", lw=1.2, zorder=1,
        )
        Nu_axis.axhspan(2.976, 3.608; color="0.85", zorder=0)
        Nu_axis.axhline(3.091; ls="--", color="0.55", lw=1.4)
        Nu_axis.text(
            23, 3.75,
            "fully developed laminar square duct,\n" *
            raw"$Nu_{H2}=3.09$ (band: $Nu_T=2.98$ to $Nu_{H1}=3.61$)";
            color="0.45", va="bottom", fontsize=7,
        )
        primary_text = @sprintf(
            "grouped (primary, solid):\n\$\\propto Re_{\\rm nom}^{%.3f}\$, \$r^2\$=%.3f\npooled (dashed):\n\$%.2f\\times10^{-4}Re_{\\rm nom}^{%.2f}\$, \$r^2\$=%.3f",
            grouped_exponent,
            Float64(grouped["r2"]),
            prefactor * 1e4,
            exponent,
            Float64(apparent["r2"]),
        )
        Nu_axis.text(21.5, 1.55, primary_text; color="0.15", ha="left", va="top", fontsize=7.5)
        Nu_axis.set_xscale("log")
        Nu_axis.set_yscale("log")
        Nu_axis.set_xlim(20, 110)
        Nu_axis.set_ylim(0.02, 6.5)
        Nu_axis.set_xticks([25, 50, 100])
        Nu_axis.set_xticklabels(["25", "50", "100"])
        suppress_minor_labels(Nu_axis)
        Nu_axis.set_xlabel(raw"Reynolds number, $Re$")
        Nu_axis.set_ylabel(raw"Apparent Nusselt number, $Nu$")
        Nu_axis.set_title("Exchange is limited at the assembly scale"; loc="left")
        Nu_axis.legend(; loc="lower right", frameon=false, handlelength=1.0)

        groups.NTU_corr = Float64.(report["ntu_profile"]["NTU_corr"])
        corrected_exponent = Float64(report["ntu_profile"]["exponent"])
        corrected_fit = linear_fit(log.(groups.Re), log.(groups.NTU_corr))
        for flux in flux_levels
            data = flux_subset(groups, flux)
            NTU_axis.plot(
                data.Re, data.NTU, "o";
                mfc="none", mec=series_color[flux], ms=4.5, mew=1.0,
            )
            NTU_axis.plot(
                data.Re, data.NTU_corr, "o";
                color=series_color[flux], ms=4.5,
            )
        end
        NTU_axis.plot(
            reynolds_line,
            exp(corrected_fit.intercept) .* reynolds_line .^ corrected_exponent,
            "-";
            color="0.25", lw=1.6,
        )
        NTU_axis.plot(
            reynolds_line,
            1.31 .* (reynolds_line ./ 72.5) .^ -1.0,
            "--";
            color="0.55", lw=1.4,
        )
        NTU_axis.set_xscale("log")
        NTU_axis.set_yscale("log")
        NTU_axis.set_xlim(20, 110)
        NTU_axis.set_ylim(0.28, 3.0)
        NTU_axis.set_xticks([25, 50, 100])
        NTU_axis.set_xticklabels(["25", "50", "100"])
        NTU_axis.set_yticks([0.5, 1, 2])
        NTU_axis.set_yticklabels(["0.5", "1", "2"])
        suppress_minor_labels(NTU_axis; y=true)
        NTU_axis.set_xlabel(raw"Reynolds number, $Re$")
        NTU_axis.set_ylabel(raw"Transfer units, $NTU$")
        NTU_axis.set_title("Transfer units rise with flow"; loc="left")
        NTU_axis.text(
            102, 2.35,
            @sprintf("measured, \$\\propto Re^{%+.2f}\$", corrected_exponent);
            color="0.15", ha="right",
        )
        NTU_axis.text(
            21.5, 0.335,
            "filled: wall-profile integration\n" * raw"open: $-\ln(1-\varepsilon)$";
            color="0.35", ha="left", va="bottom", fontsize=7,
        )
        NTU_axis.text(
            102, 0.335,
            "conductance fixed in \$z\$,\n" * raw"$\propto Re^{-1}$";
            color="0.5", ha="right",
        )

        panel_letter(Nu_axis, "a")
        panel_letter(NTU_axis, "b")
        finish_figure(figure, joinpath(fig_dir, "fig4_assembly_limitation.png"))
    end
end

begin # Figure 5: temperature inversion and local nonequilibrium
    function figure_inversion_ltne(report, groups, fig_dir)
        figure, axes = plt.subplots(1, 3; figsize=(7.4, 2.9))
        flow_axis, effectiveness_axis, deficit_axis = axes
        crossing_effectiveness = Float64[]

        for flux in flux_levels
            data = flux_subset(groups, flux)
            crossing = report["crossings"][string(flux)]
            q_local = json_number(crossing, "q_local")
            eps_local = json_number(crossing, "eps_local")

            flow_axis.plot(
                data.q_slpm, data.I_vol, "o-";
                color=series_color[flux], ms=4, lw=1.1,
            )
            effectiveness_axis.plot(
                data.eps, data.I_vol, "o";
                color=series_color[flux], ms=4.5, label=string(flux),
            )
            deficit_axis.plot(
                data.Re, data.Lam107, "o-";
                color=series_color[flux], ms=4, lw=1.1,
            )
            deficit_axis.plot(
                data.Re, data.Lam58, "s--";
                color=series_color[flux], ms=3.4, lw=1.0, mfc="white",
            )

            if isfinite(q_local)
                flow_axis.plot([q_local], [0], "^"; color=series_color[flux], ms=7, zorder=5)
            end
            if isfinite(eps_local)
                effectiveness_axis.plot(
                    [eps_local], [0], "^";
                    color=series_color[flux], ms=7, zorder=5,
                )
                push!(crossing_effectiveness, eps_local)
            end
        end

        flow_axis.axhline(0; color="0.6", lw=0.8, zorder=0)
        effectiveness_axis.axhline(0; color="0.6", lw=0.8, zorder=0)
        flow_axis.set_xlabel(raw"$q$ [sL min$^{-1}$]")
        flow_axis.set_ylabel(raw"$T_{12}-T_8$ [K]")
        flow_axis.set_title("Front-to-mid side-wall crossing"; loc="left")
        flow_axis.annotate(
            "single negative\npoint only";
            xy=(5.08, 0), xytext=(6.4, -88), fontsize=7, color="0.3",
            arrowprops=Dict("arrowstyle" => "-", "lw" => 0.7, "color" => "0.45"),
        )

        if !isempty(crossing_effectiveness)
            lower = minimum(crossing_effectiveness)
            upper = maximum(crossing_effectiveness)
            effectiveness_axis.axvspan(lower, upper; color="0.85", zorder=0)
            effectiveness_axis.text(
                0.5 * (lower + upper), 86, raw"$\varepsilon^*$";
                ha="center", color="0.3",
            )
        end
        effectiveness_axis.set_xlabel(raw"Effectiveness, $\varepsilon$")
        effectiveness_axis.set_title(raw"Crossings collapse on $\varepsilon$"; loc="left")
        effectiveness_axis.legend(
            ; title=raw"kW m$^{-2}$", loc="lower right", frameon=false,
            handlelength=0.8, fontsize=7, title_fontsize=7,
        )

        deficit_axis.set_xlabel(raw"Reynolds number, $Re$")
        deficit_axis.set_ylabel(raw"Apparent wall-to-interior deficit, $\Lambda$")
        deficit_axis.set_title(raw"Wall-to-interior deficit grows with $Re$"; loc="left")
        deficit_axis.text(90, 0.108, raw"$\Lambda_{107}$"; ha="right", color="0.25")
        deficit_axis.text(90, 0.040, raw"$\Lambda_{58}$"; ha="right", color="0.45")

        for (axis, letter) in zip(axes, ("a", "b", "c"))
            panel_letter(axis, letter)
        end
        finish_figure(figure, joinpath(fig_dir, "fig5_inversion_ltne.png"))
    end
end

begin # Figure 6: transient eigenvalue identification
    function figure_transient_identification(report, eigenvalues, cooling_data,
                                             master_curves, fig_dir)
        figure, axes = plt.subplots(1, 3; figsize=(7.4, 2.9))
        cooling_axis, eigenvalue_axis, collapse_axis = axes

        for run_group in groupby(cooling_data, :ID)
            ID = String(first(run_group.ID))
            for sensor_group in groupby(run_group, :sensor)
                cooling_axis.semilogy(
                    sensor_group.t ./ 60,
                    max.(sensor_group.theta, 1e-2);
                    color=cooling_color[ID], lw=0.7, alpha=0.75,
                )
            end
        end
        cooling_axis.set_xlabel("Time [min]")
        cooling_axis.set_ylabel(raw"Normalised excess, $\theta/\theta_0$")
        cooling_axis.set_title("One shared slow mode"; loc="left")
        cooling_axis.set_ylim(1e-1, 1.4)
        cooling_axis.set_xlim(0, 105)
        cooling_axis.set_yticks([0.1, 0.2, 0.5, 1.0])
        cooling_axis.set_yticklabels(["0.1", "0.2", "0.5", "1"])
        suppress_minor_labels(cooling_axis; x=false, y=true)
        for ID in ("C69", "C80", "C81")
            cooling_axis.plot(
                [], [];
                color=cooling_color[ID], lw=1.4,
                label=string(cooling_flow[ID], raw" sL min$^{-1}$"),
            )
        end
        cooling_axis.legend(; loc="lower left", frameon=false, handlelength=1.2, fontsize=7)

        identification_series = (
            ("cool", "o", raw"cooling, matched $\varepsilon$ (primary)",
             "#1f4e79", "cooling_matched_eps"),
            ("heat", "s", "heating, deep 3", "#c1666b", "heating_deep"),
            ("heat6", "^", "heating, all 6", "#8f9779", "heating_all6"),
        )
        for (phase, marker, label, color, key) in identification_series
            selected = eigenvalues[
                (eigenvalues.phase .== phase) .& .!ismissing.(eigenvalues.lam),
                :,
            ]
            lambda = Float64.(selected.lam)
            x_values = phase == "cool" ? Float64.(selected.x_matched) : Float64.(selected.x)
            if phase == "cool"
                eigenvalue_axis.plot(
                    Float64.(selected.x), lambda .* 1e3, marker;
                    mfc="none", mec=color, ms=4.2, mew=1.0,
                    label=raw"cooling, pooled $\varepsilon$",
                )
            end
            eigenvalue_axis.plot(
                x_values, lambda .* 1e3, marker;
                color=color, ms=4.2, label=label,
            )

            exchange = collect(range(0.04, 0.32; length=20))
            fit = report["identification"][key]
            eigenvalue_axis.plot(
                exchange,
                (exchange .+ Float64(fit["K_loss"])) ./ Float64(fit["C_eff"]) .* 1e3,
                "-";
                color=color, lw=1.2,
            )
        end
        eigenvalue_axis.set_xlabel(raw"$\varepsilon\,\dot m c_p$ [W K$^{-1}$]")
        eigenvalue_axis.set_ylabel(raw"$\lambda$ [$10^{-3}$ s$^{-1}$]")
        eigenvalue_axis.set_title(raw"Probe set sets $C_{\rm eff}$"; loc="left")
        eigenvalue_axis.legend(; loc="upper left", frameon=false, handlelength=1.2, fontsize=7)

        for (signal, line_style, color) in
            (("wall", "-", "0.35"), ("gas", "--", "#c1666b"))
            signal_data = master_curves[master_curves.signal .== signal, :]
            for run_group in groupby(signal_data, :ID)
                collapse_axis.plot(
                    run_group.tstar, run_group.theta, line_style;
                    color=color, lw=0.6, alpha=0.8,
                )
            end
        end
        collapse_axis.set_xlim(0, 2.2)
        collapse_axis.set_ylim(0, 1.08)
        collapse_axis.set_xlabel(raw"Rescaled time, $t^*$")
        collapse_axis.set_ylabel(raw"Normalised rise, $\theta^*$")
        collapse_axis.set_title("Transients collapse"; loc="left")
        collapse_axis.text(1.45, 0.45, "wall"; color="0.3")
        collapse_axis.text(1.45, 0.26, "gas outlet"; color="#c1666b")

        for (axis, letter) in zip(axes, ("a", "b", "c"))
            panel_letter(axis, letter)
        end
        finish_figure(figure, joinpath(fig_dir, "fig6_transient_identification.png"))
    end
end

begin # Figure 7: consequences of under-instrumentation
    function figure_under_instrumentation(report, reference_points, fig_dir)
        reynolds_line = collect(range(21.0, 102.0; length=60))
        figure, axes = plt.subplots(1, 2; figsize=(7.4, 3.1))
        reference_axis, capacitance_axis = axes

        for reference in reference_order
            data = reference_points[reference_points.reference .== reference, :]
            fit = power_law_fit(data.Re, data.Nu)
            reference_axis.plot(
                data.Re, data.Nu, "o";
                color=reference_color[reference], ms=3.2, alpha=0.7,
            )
            reference_axis.plot(
                reynolds_line,
                fit.prefactor .* reynolds_line .^ fit.exponent,
                "-";
                color=reference_color[reference], lw=1.5,
                label=@sprintf("%s: \$Re^{%.2f}\$", reference, fit.exponent),
            )
        end
        reference_axis.axhline(3.091; ls="--", color="0.55", lw=1.2)
        reference_axis.text(21, 3.4, "laminar duct"; color="0.45", va="bottom", fontsize=7)
        reference_axis.set_xscale("log")
        reference_axis.set_yscale("log")
        reference_axis.set_xlim(20, 112)
        reference_axis.set_ylim(0.012, 9.0)
        reference_axis.set_xticks([25, 50, 100])
        reference_axis.set_xticklabels(["25", "50", "100"])
        suppress_minor_labels(reference_axis)
        reference_axis.set_xlabel(raw"Reynolds number, $Re$")
        reference_axis.set_ylabel(raw"Apparent Nusselt number, $Nu$")
        reference_axis.set_title("Same data, five reference probes"; loc="left")
        reference_axis.legend(; loc="lower right", frameon=false, handlelength=1.2, fontsize=6.5)

        identification_keys = ("cooling", "heating_deep", "heating_all6", "joint")
        labels = ["cooling\n6 probes", "heating\ndeep 3", "heating\nall 6", "joint\n18 eigenvalues"]
        values = [Float64(report["identification"][key]["C_eff"])
                  for key in identification_keys]
        colors = ["#1f4e79", "#c1666b", "#8f9779", "0.35"]
        positions = collect(0:3)
        capacitance_axis.vlines(positions, zeros(4), values; color=colors, lw=1.2)
        for (position, value, color) in zip(positions, values, colors)
            capacitance_axis.plot([position], [value], "o"; ms=8, color=color, zorder=4)
            capacitance_axis.text(
                position, value + 16, @sprintf("%.0f", value);
                ha="center", fontsize=7, color="0.2",
            )
        end

        window_values = [Float64(item["C_eff"])
                         for item in report["heating_window_sensitivity"]]
        lower, upper = extrema(window_values)
        capacitance_axis.plot(
            [1, 1], [lower, upper];
            lw=6, color="#c1666b", alpha=0.25,
            solid_capstyle="butt", zorder=1,
        )
        capacitance_axis.annotate(
            @sprintf("fit-window\nrange %.0f–%.0f", lower, upper);
            xy=(1.07, 0.5 * (lower + upper)), xytext=(1.38, 150),
            fontsize=6.5, color="#8a4a4e",
            arrowprops=Dict("arrowstyle" => "-", "lw" => 0.7, "color" => "#c1666b"),
        )
        monolith = report["monolith"]
        capacitance_axis.axhspan(
            Float64(monolith["C_monolith_600K"]),
            Float64(monolith["C_monolith_900K"]);
            color="0.8", zorder=0,
        )
        capacitance_axis.text(3.45, 58, "bare monolith"; ha="right", fontsize=7, color="0.35")
        capacitance_axis.set_xticks(positions)
        capacitance_axis.set_xticklabels(labels; fontsize=7)
        capacitance_axis.set_xlim(-0.5, 3.5)
        capacitance_axis.set_ylim(0, 345)
        capacitance_axis.set_ylabel(raw"Effective capacitance, $C_{\rm eff}$ [J K$^{-1}$]")
        capacitance_axis.set_title("Same runs, four defensible answers"; loc="left")

        panel_letter(reference_axis, "a")
        panel_letter(capacitance_axis, "b")
        finish_figure(figure, joinpath(fig_dir, "fig7_under_instrumentation.png"))
    end
end

begin # complete figure workflow
    function make_figures(raw_dir, out_dir, fig_dir)
        manuscript_style()
        mkpath(fig_dir)

        println("Running make_figures.jl")
        println("  reduction archive: ", out_dir)
        println("  figures:           ", fig_dir)
        flush(stdout)

        if trace_files_are_stale(out_dir)
            println("Refreshing plotting traces from the raw logger files ...")
            flush(stdout)
            export_traces(raw_dir, out_dir)
        end

        report = JSON3.read(
            read(joinpath(out_dir, "results.json"), String),
            Dict{String,Any},
        )
        groups = CSV.read(joinpath(out_dir, "groups.csv"), DataFrame)
        eigenvalues = CSV.read(joinpath(out_dir, "eigenvalues.csv"), DataFrame)
        cooling_data = CSV.read(joinpath(out_dir, "cooling_decays.csv"), DataFrame)
        master_curves = CSV.read(joinpath(out_dir, "master_curves.csv"), DataFrame)
        reference_points = CSV.read(joinpath(out_dir, "reference_points.csv"), DataFrame)

        figure_steady_field(groups, fig_dir)
        figure_assembly_limitation(report, groups, fig_dir)
        figure_inversion_ltne(report, groups, fig_dir)
        figure_transient_identification(
            report, eigenvalues, cooling_data, master_curves, fig_dir,
        )
        figure_under_instrumentation(report, reference_points, fig_dir)
        convert_apparatus(fig_dir)

        println("wrote fig3--fig7 to ", fig_dir)
        return nothing
    end
end

begin # command-line execution
    function figure_command_line_options(arguments)
        options = Dict(
            "raw" => normpath(joinpath(@__DIR__, "..", "RAW")),
            "out" => joinpath(@__DIR__, "outputs"),
            "fig" => joinpath(@__DIR__, "figures"),
        )
        index = 1
        while index <= length(arguments)
            argument = arguments[index]
            startswith(argument, "--") || error("unknown argument: $argument")
            if occursin('=', argument)
                key, value = split(argument[3:end], '='; limit=2)
                options[key] = value
            else
                index == length(arguments) && error("missing value for $argument")
                options[argument[3:end]] = arguments[index + 1]
                index += 1
            end
            index += 1
        end
        return options
    end

    function figures_main(arguments=ARGS)
        options = figure_command_line_options(arguments)
        make_figures(
            abspath(options["raw"]),
            abspath(options["out"]),
            abspath(options["fig"]),
        )
    end
end

if !isdefined(@__MODULE__, :MAKE_FIGURES_INCLUDE_ONLY) ||
   !getfield(@__MODULE__, :MAKE_FIGURES_INCLUDE_ONLY)
    figures_main()
end
