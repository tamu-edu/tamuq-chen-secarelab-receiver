# ============================================================================
# Diagnostic plots for the v20 gate-free stress study.
#
# Stress-test parameters are deliberately broader than the independently
# supported physical ranges.  These plots diagnose attainable model-data
# agreement; they do not validate fitted material, power, UA, or probe values.
# Missing optional stress validation/source/power files are skipped.
# ============================================================================

ENV["GKSwstype"] = "100"

using CSV
using DataFrames
using Plots
using Plots.PlotMeasures
using StatsPlots
using Printf
using Statistics

const V20_SUMMARY_DIR = joinpath(@__DIR__, "summaries", "2D_v20")
const V20_PLOT_DIR = joinpath(V20_SUMMARY_DIR, "plots")
mkpath(V20_PLOT_DIR)

const ORIGINAL_THRESHOLDS = (
    heating_side=60.0,
    heating_T2=15.0,
    cooling_side=35.0,
    cooling_T2=12.0,
)

function read_csv_safe(path)
    isfile(path) || return nothing
    try
        return CSV.read(path, DataFrame)
    catch error
        @warn "Skipping unreadable CSV" path exception=(error, catch_backtrace())
        return nothing
    end
end

has_columns(frame, required) =
    frame !== nothing && all(String(column) in names(frame) for column in required)

function matching_csvs(predicate)
    isdir(V20_SUMMARY_DIR) || return String[]
    return sort([
        joinpath(V20_SUMMARY_DIR, name)
        for name in readdir(V20_SUMMARY_DIR)
        if endswith(lowercase(name), ".csv") && predicate(name)
    ])
end

function finite_float(value)
    result = try
        Float64(value)
    catch
        NaN
    end
    return result
end

function plot_ua_profile()
    paths = matching_csvs(name ->
        startswith(name, "stress_ua_profile") &&
        !occursin("summary", lowercase(name))
    )
    best = Dict{Tuple{Float64,Float64},NamedTuple}()
    for path in paths
        frame = read_csv_safe(path)
        required = ("reynolds_exponent", "Nu50_ratio", "objective")
        has_columns(frame, required) || continue
        for row in eachrow(frame)
            n = finite_float(row.reynolds_exponent)
            ratio = finite_float(row.Nu50_ratio)
            objective = finite_float(row.objective)
            all(isfinite, (n, ratio, objective)) || continue
            key = (n, ratio)
            candidate = (
                n=n, ratio=ratio, objective=objective,
                source=basename(path),
            )
            if !haskey(best, key) || objective < best[key].objective
                best[key] = candidate
            end
        end
    end
    if isempty(best)
        @warn "No gate-free UA profile CSVs were available; UA plot omitted."
        return nothing, nothing
    end

    rows = collect(values(best))
    exponents = sort(unique(row.n for row in rows))
    palette = cgrad(:viridis, max(length(exponents), 2), categorical=true)
    figure = plot(
        xlabel="Nu50 / measured-correlation Nu50",
        ylabel="T3-free plant objective",
        title="v20 gate-free UA profile — DIAGNOSTIC STRESS TEST",
        legend=:outerright,
        gridalpha=0.25,
        size=(980, 620),
        dpi=160,
    )
    for (index, exponent) in enumerate(exponents)
        selected = sort(
            filter(row -> row.n == exponent, rows);
            by=row -> row.ratio,
        )
        plot!(
            figure,
            [row.ratio for row in selected],
            [row.objective for row in selected];
            marker=:circle,
            markersize=5,
            linewidth=2,
            color=palette[index],
            label=@sprintf("n = %.2f", exponent),
        )
    end
    annotation = "Broad stress ranges; minima are not validated UA parameters"
    annotate!(figure, minimum(row.ratio for row in rows),
              maximum(row.objective for row in rows),
              text(annotation, 9, :left, :top))
    output = joinpath(V20_PLOT_DIR, "gatefree_ua_objective_profile_2D_v20.png")
    savefig(figure, output)
    global_best = first(sort(rows; by=row -> row.objective))
    return output, global_best
end

function summary_candidates(paths)
    required = (
        "heating_side_rmse_K", "heating_T2_rmse_K",
        "cooling_side_rmse_K", "cooling_T2_rmse_K",
    )
    candidates = NamedTuple[]
    for path in paths
        frame = read_csv_safe(path)
        has_columns(frame, required) || continue
        for (row_index, row) in enumerate(eachrow(frame))
            values = (
                finite_float(row.heating_side_rmse_K),
                finite_float(row.heating_T2_rmse_K),
                finite_float(row.cooling_side_rmse_K),
                finite_float(row.cooling_T2_rmse_K),
            )
            all(isfinite, values) || continue
            normalized = mean(values ./ (
                ORIGINAL_THRESHOLDS.heating_side,
                ORIGINAL_THRESHOLDS.heating_T2,
                ORIGINAL_THRESHOLDS.cooling_side,
                ORIGINAL_THRESHOLDS.cooling_T2,
            ))
            push!(candidates, (
                values=values,
                normalized_score=normalized,
                source=basename(path),
                row_index=row_index,
            ))
        end
    end
    return candidates
end

function plot_rmse_comparison()
    original_paths = matching_csvs(name ->
        startswith(name, "plant_nuisance_selected") &&
        !startswith(name, "stress_")
    )
    # Optional future validation/source/power summaries are harmless here:
    # only files containing all four aggregate plant RMSE columns are eligible.
    stress_paths = matching_csvs(name ->
        startswith(name, "stress_") && (
            occursin("selected", name) ||
            occursin("validation", name) ||
            occursin("summary", name) ||
            occursin("source", name) ||
            occursin("power", name)
        )
    )
    original = summary_candidates(original_paths)
    stress = summary_candidates(stress_paths)
    if isempty(original) || isempty(stress)
        @warn "Bounded or gate-free four-metric summary absent; RMSE comparison omitted." original_count=length(original) stress_count=length(stress)
        return nothing, nothing, nothing
    end
    original_best = first(sort(original; by=row -> row.normalized_score))
    stress_best = first(sort(stress; by=row -> row.normalized_score))

    categories = ["Heating\nside", "Heating\nT2", "Cooling\nside", "Cooling\nT2"]
    thresholds = [
        ORIGINAL_THRESHOLDS.heating_side,
        ORIGINAL_THRESHOLDS.heating_T2,
        ORIGINAL_THRESHOLDS.cooling_side,
        ORIGINAL_THRESHOLDS.cooling_T2,
    ]
    values = hcat(collect(original_best.values), collect(stress_best.values))
    figure = groupedbar(
        categories,
        values;
        bar_position=:dodge,
        label=[
            "Original bounded representative screen" "Gate-free representative screen (3 anchors)"
        ],
        color=[:steelblue :darkorange],
        ylabel="Equal-case mean RMSE (K)",
        title="v20 representative-screen RMSE — not full validation",
        legend=:topright,
        gridalpha=0.25,
        size=(1000, 620),
        dpi=160,
    )
    scatter!(
        figure, 1:4, thresholds;
        marker=:diamond,
        markersize=8,
        markercolor=:black,
        markerstrokecolor=:white,
        label="Original thresholds",
    )
    for (index, threshold) in enumerate(thresholds)
        plot!(
            figure, [index - 0.36, index + 0.36], [threshold, threshold];
            color=:black, linestyle=:dash, linewidth=1.5, label=false,
        )
    end
    output = joinpath(V20_PLOT_DIR, "bounded_vs_gatefree_rmse_2D_v20.png")
    savefig(figure, output)
    return output, original_best, stress_best
end

function bool_value(value)
    value isa Bool && return value
    return lowercase(strip(String(value))) == "true"
end

function load_alltrain_validation()
    phase_path = joinpath(
        V20_SUMMARY_DIR,
        "gatefree_validation_phases_alltrain_power_2D_v20.csv",
    )
    summary_path = joinpath(
        V20_SUMMARY_DIR,
        "gatefree_validation_summary_alltrain_power_2D_v20.csv",
    )
    phases = read_csv_safe(phase_path)
    summary = read_csv_safe(summary_path)
    phase_required = (
        "phase", "side_pooled_rmse_K", "T2_pooled_rmse_K",
        "original_side_limit_K", "original_T2_limit_K",
        "original_plant_threshold_pass",
    )
    has_columns(phases, phase_required) || return nothing
    nrow(phases) > 0 || return nothing
    summary_row = summary !== nothing && nrow(summary) > 0 ?
                  summary[1, :] : nothing
    return (
        phases=phases,
        summary=summary_row,
        phase_path=phase_path,
        summary_path=summary_path,
    )
end

function plot_train_holdout_rmse()
    validation = load_alltrain_validation()
    if validation === nothing
        @warn "Final all-training-power phase CSV absent; train/holdout plot omitted."
        return nothing, nothing
    end
    order = (
        "heating_training", "heating_holdout",
        "cooling_training", "cooling_holdout",
    )
    labels = ["Heating\ntrain", "Heating\nholdout",
              "Cooling\ntrain", "Cooling\nholdout"]
    rows = NamedTuple[]
    for phase in order
        indices = findall(String.(validation.phases.phase) .== phase)
        isempty(indices) && continue
        row = validation.phases[first(indices), :]
        push!(rows, (
            phase=phase,
            side=finite_float(row.side_pooled_rmse_K),
            T2=finite_float(row.T2_pooled_rmse_K),
            side_limit=finite_float(row.original_side_limit_K),
            T2_limit=finite_float(row.original_T2_limit_K),
            pass=bool_value(row.original_plant_threshold_pass),
        ))
    end
    length(rows) == length(order) || begin
        @warn "Final phase CSV did not contain all four expected phases; plot omitted."
        return nothing, validation
    end

    phase_colors = [:steelblue, :lightskyblue, :darkorange, :navajowhite3]
    x = collect(eachindex(rows))
    common = (
        legend=false,
        gridalpha=0.25,
        xticks=(x, labels),
        xrotation=0,
    )
    p_side = bar(
        x, [row.side for row in rows];
        color=phase_colors,
        ylabel="Side pooled RMSE (K)",
        left_margin=12mm,
        top_margin=8mm,
        common...,
    )
    p_T2 = bar(
        x, [row.T2 for row in rows];
        color=phase_colors,
        ylabel="T2 pooled RMSE (K)",
        left_margin=8mm,
        top_margin=8mm,
        common...,
    )
    ylims!(
        p_side, 0.0,
        1.08 * maximum(vcat(
            [row.side for row in rows],
            [row.side_limit for row in rows],
        )),
    )
    ylims!(
        p_T2, 0.0,
        1.10 * maximum(vcat(
            [row.T2 for row in rows],
            [row.T2_limit for row in rows],
        )),
    )
    for (plot_object, metric, limit) in (
        (p_side, :side, :side_limit),
        (p_T2, :T2, :T2_limit),
    )
        limits = [getproperty(row, limit) for row in rows]
        scatter!(
            plot_object, x, limits;
            marker=:diamond, markersize=7,
            markercolor=:black, markerstrokecolor=:white,
            label=false,
        )
        for index in x
            threshold = limits[index]
            plot!(
                plot_object,
                [index - 0.34, index + 0.34],
                [threshold, threshold];
                color=:black, linestyle=:dash,
                linewidth=1.5, label=false,
            )
            value = getproperty(rows[index], metric)
            metric_pass = value <= threshold
            annotate!(
                plot_object, index, value,
                text(metric_pass ? "PASS" : "FAIL", 8,
                     metric_pass ? :darkgreen : :darkred,
                     :center, :bottom),
            )
        end
    end
    figure = plot(
        p_side, p_T2;
        layout=(1, 2),
        size=(1400, 650),
        dpi=160,
        plot_titlefontsize=16,
        plot_title=(
            "v20 all-training-power train/holdout RMSE — " *
            "DIAGNOSTIC STRESS TEST\n" *
            "Diamonds/dashes are original thresholds; labels are metric-level"
        ),
    )
    output = joinpath(
        V20_PLOT_DIR, "gatefree_train_holdout_rmse_2D_v20.png",
    )
    savefig(figure, output)
    return output, validation
end

function choose_t3_branch()
    selected_paths = matching_csvs(name ->
        startswith(name, "stress_t3_selected")
    )
    candidates = NamedTuple[]
    for path in selected_paths
        frame = read_csv_safe(path)
        required = (
            "location", "capacity_areal_J_m2_K",
            "stem_conductance_areal_W_m2_K", "training_objective",
        )
        has_columns(frame, required) || continue
        for row in eachrow(frame)
            objective = finite_float(row.training_objective)
            isfinite(objective) || continue
            push!(candidates, (
                location=String(row.location),
                stem=finite_float(row.stem_conductance_areal_W_m2_K),
                capacity=finite_float(row.capacity_areal_J_m2_K),
                objective=objective,
                source=basename(path),
                diagnostic=true,
            ))
        end
    end
    if !isempty(candidates)
        return first(sort(candidates; by=row -> row.objective))
    end

    original_path = joinpath(
        V20_SUMMARY_DIR, "t3_observer_profile_2D_v20.csv",
    )
    frame = read_csv_safe(original_path)
    required = (
        "location", "capacity_areal_J_m2_K",
        "stem_conductance_areal_W_m2_K", "training_objective",
    )
    has_columns(frame, required) || return nothing
    index = argmin(Float64.(frame.training_objective))
    row = frame[index, :]
    return (
        location=String(row.location),
        stem=finite_float(row.stem_conductance_areal_W_m2_K),
        capacity=finite_float(row.capacity_areal_J_m2_K),
        objective=finite_float(row.training_objective),
        source=basename(original_path),
        diagnostic=false,
    )
end

function plot_t3_profile()
    branch = choose_t3_branch()
    if branch === nothing
        @warn "No T3 observer selection/profile was available; T3 plot omitted."
        return nothing, nothing
    end
    profile_paths = branch.diagnostic ?
        matching_csvs(name -> startswith(name, "stress_t3_profile")) :
        [joinpath(V20_SUMMARY_DIR, "t3_observer_profile_2D_v20.csv")]
    points = Dict{Float64,NamedTuple}()
    for path in profile_paths
        frame = read_csv_safe(path)
        required = (
            "location", "capacity_areal_J_m2_K",
            "stem_conductance_areal_W_m2_K", "training_objective",
            "training_T3_rmse_K", "validation_T3_rmse_K",
        )
        has_columns(frame, required) || continue
        for row in eachrow(frame)
            location = String(row.location)
            stem = finite_float(row.stem_conductance_areal_W_m2_K)
            location == branch.location || continue
            isapprox(stem, branch.stem; atol=1e-9, rtol=1e-9) || continue
            capacity = finite_float(row.capacity_areal_J_m2_K)
            objective = finite_float(row.training_objective)
            training_rmse = finite_float(row.training_T3_rmse_K)
            validation_rmse = finite_float(row.validation_T3_rmse_K)
            all(isfinite, (capacity, objective, training_rmse, validation_rmse)) ||
                continue
            capacity > 0.0 || continue
            candidate = (
                capacity=capacity,
                objective=objective,
                training_rmse=training_rmse,
                validation_rmse=validation_rmse,
                source=basename(path),
            )
            if !haskey(points, capacity) ||
               objective < points[capacity].objective
                points[capacity] = candidate
            end
        end
    end
    if isempty(points)
        @warn "Winning T3 branch had no matching profile rows; T3 plot omitted." branch
        return nothing, branch
    end
    rows = sort(collect(values(points)); by=row -> row.capacity)
    diagnostic_text = branch.diagnostic ? "DIAGNOSTIC STRESS TEST" :
                      "bounded observer profile"
    branch_text = @sprintf(
        "%s, stem = %.1f W m^-2 K^-1", branch.location, branch.stem,
    )
    p_objective = plot(
        [row.capacity for row in rows],
        [row.objective for row in rows];
        xscale=:log10,
        marker=:circle,
        linewidth=2,
        color=:purple,
        label="Training objective",
        xlabel="Probe areal capacity (J m^-2 K^-1)",
        ylabel="Objective",
        title="T3 objective: $branch_text",
        gridalpha=0.25,
    )
    p_rmse = plot(
        [row.capacity for row in rows],
        [row.training_rmse for row in rows];
        xscale=:log10,
        marker=:circle,
        linewidth=2,
        color=:steelblue,
        label="C69/C80 training",
        xlabel="Probe areal capacity (J m^-2 K^-1)",
        ylabel="T3 RMSE (K)",
        title="T3 RMSE: $diagnostic_text",
        gridalpha=0.25,
    )
    plot!(
        p_rmse,
        [row.capacity for row in rows],
        [row.validation_rmse for row in rows];
        marker=:diamond,
        linewidth=2,
        color=:darkorange,
        label="C81 validation",
    )
    figure = plot(
        p_objective, p_rmse;
        layout=(1, 2), size=(1200, 520), dpi=160,
        plot_title="v20 T3 probe profile — $diagnostic_text",
    )
    output = joinpath(V20_PLOT_DIR, "t3_capacity_profile_2D_v20.png")
    savefig(figure, output)
    return output, branch
end

function write_summary(
    ua_output, ua_best, rmse_output, original_best, stress_best,
    train_holdout_output, alltrain_validation,
    t3_output, t3_branch,
)
    output = joinpath(
        V20_SUMMARY_DIR, "gatefree_plot_summary_2D_v20.txt",
    )
    open(output, "w") do io
        println(io, "V20 GATE-FREE DIAGNOSTIC PLOT SUMMARY")
        println(io, "Stress-test parameters are diagnostic and are not validated fits.")
        println(io)
        println(io, "UA plot: ", something(ua_output, "omitted (data absent)"))
        if ua_best !== nothing
            println(io, @sprintf(
                "UA diagnostic minimum: n=%.6g, Nu50 ratio=%.6g, objective=%.6g",
                ua_best.n, ua_best.ratio, ua_best.objective,
            ))
            println(io, "UA source row: ", ua_best.source)
        end
        println(io)
        println(io, "RMSE plot: ", something(rmse_output, "omitted (data absent)"))
        if original_best !== nothing
            println(io, "Original bounded source: ", original_best.source)
            println(io, "Original RMSE [Hside, HT2, Cside, CT2]: ",
                    join(original_best.values, ", "))
        end
        if stress_best !== nothing
            println(io, "Gate-free representative-screen source: ",
                    stress_best.source)
            println(io, "Gate-free representative-screen RMSE ",
                    "[Hside, HT2, Cside, CT2]: ",
                    join(stress_best.values, ", "))
        end
        println(io, "The preceding orange comparison values are ",
                "representative-screen metrics, not full validation.")
        println(io, "Original thresholds [Hside, HT2, Cside, CT2]: ",
                join((
                    ORIGINAL_THRESHOLDS.heating_side,
                    ORIGINAL_THRESHOLDS.heating_T2,
                    ORIGINAL_THRESHOLDS.cooling_side,
                    ORIGINAL_THRESHOLDS.cooling_T2,
                ), ", "))
        println(io)
        println(io, "All-training-power train/holdout plot: ",
                something(train_holdout_output, "omitted (data absent)"))
        if alltrain_validation !== nothing
            summary = alltrain_validation.summary
            if summary !== nothing
                println(io, "FINAL ALL-TRAINING-POWER DIAGNOSTIC PARAMETERS")
                for name in (
                    :mesh, :reynolds_exponent, :Nu50_ratio, :Nu50,
                    :nu_prefactor, :felt_conductivity_scale,
                    :felt_heat_capacity_scale, :felt_contact_scale,
                    :power_scale_456, :power_scale_304, :power_scale_256,
                    :T3_location, :T3_capacity_areal_J_m2_K,
                    :T3_stem_conductance_areal_W_m2_K,
                )
                    name in propertynames(summary) || continue
                    println(io, name, "=", getproperty(summary, name))
                end
                for name in (
                    :all_ode_success, :all_flow_converged,
                    :original_heating_training_gate,
                    :original_heating_holdout_gate,
                    :original_cooling_training_gate,
                    :original_cooling_holdout_gate,
                    :observer_inside_original_bracket,
                    :original_cooling_observer_gate,
                    :original_all_posthoc_thresholds_pass,
                )
                    name in propertynames(summary) || continue
                    println(io, name, "=", getproperty(summary, name))
                end
            end
            println(io, "FINAL PHASE METRICS")
            for row in eachrow(alltrain_validation.phases)
                println(io, @sprintf(
                    "%s: side_RMSE=%.6g K (limit %.6g), T2_RMSE=%.6g K (limit %.6g), plant_pass=%s",
                    String(row.phase),
                    finite_float(row.side_pooled_rmse_K),
                    finite_float(row.original_side_limit_K),
                    finite_float(row.T2_pooled_rmse_K),
                    finite_float(row.original_T2_limit_K),
                    string(bool_value(row.original_plant_threshold_pass)),
                ))
                if "T3_observer_pooled_rmse_K" in names(alltrain_validation.phases)
                    println(io, @sprintf(
                        "  T3_RMSE=%.6g K, T3_final_bias=%.6g K, T3_t90_MAE=%.6g s, observer_pass=%s",
                        finite_float(row.T3_observer_pooled_rmse_K),
                        finite_float(row.T3_observer_final_bias_K),
                        finite_float(row.T3_observer_t90_mae_s),
                        string(bool_value(
                            row.original_T3_performance_threshold_pass,
                        )),
                    ))
                end
            end
            println(io, "Final phase source: ",
                    basename(alltrain_validation.phase_path))
            println(io, "Final parameter source: ",
                    basename(alltrain_validation.summary_path))
        end
        println(io)
        println(io, "T3 plot: ", something(t3_output, "omitted (data absent)"))
        if t3_branch !== nothing
            println(io, "T3 branch source: ", t3_branch.source)
            println(io, @sprintf(
                "T3 branch: location=%s, stem=%.6g, selected capacity=%.6g",
                t3_branch.location, t3_branch.stem, t3_branch.capacity,
            ))
        end
    end
    return output
end

function main()
    ua_output, ua_best = plot_ua_profile()
    rmse_output, original_best, stress_best = plot_rmse_comparison()
    train_holdout_output, alltrain_validation =
        plot_train_holdout_rmse()
    t3_output, t3_branch = plot_t3_profile()
    summary_output = write_summary(
        ua_output, ua_best, rmse_output, original_best, stress_best,
        train_holdout_output, alltrain_validation,
        t3_output, t3_branch,
    )
    println("v20 diagnostic plot directory: ", V20_PLOT_DIR)
    println("v20 diagnostic summary: ", summary_output)
end

main()
