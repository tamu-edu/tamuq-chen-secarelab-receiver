# ============================================================================
# calibrate_2D_v21.jl
#
# Isolated v21 fit.  Only the common effective hot-front area and one power
# multiplier per lamp group move.  All transport, source, felt, contact, and
# observation parameters remain frozen.
# ============================================================================

using Statistics

include("run_2D_v21.jl")

const FRONT_AREA_GRID_CM2_2D_v21 = (3.61, 7.0, 12.0, 20.0, 40.0)
const POWER_GRID_2D_v21 = Dict(
    456 => (1.05, 1.19, 1.34, 1.50, 1.65),
    304 => (1.05, 1.23, 1.40, 1.58, 1.80),
    256 => (0.70, 0.84, 0.94, 1.11, 1.25),
)
const GROUP_IDS_2D_v21 = Dict(
    456 => ("E67", "E69", "E71"),
    304 => ("E72", "E74", "E76"),
    256 => ("E77", "E79", "E81"),
)

function _final_rows_2D_v21(times)
    length(times) >= 2 || error(
        "v21 calibration requires at least two samples",
    )
    threshold = first(times) +
        0.90 * (last(times) - first(times))
    start = something(
        findfirst(time -> time >= threshold, times),
        length(times),
    )
    return start:length(times)
end

function _calibration_metrics_2D_v21(case)
    final = _final_rows_2D_v21(case.inputs.times)
    residual = case.model .- case.observed
    side_final = [
        mean(residual[final, sensor]) for sensor in 1:3
    ]
    T2_final = mean(residual[final, 7])
    T3_final = mean(residual[final, 6])
    # Equal high-level emphasis on side profile, felt/casing, and outlet air.
    final_objective =
        mean((side_final ./ 50.0).^2) / 3.0 +
        (T2_final / 35.0)^2 / 3.0 +
        (T3_final / 60.0)^2 / 3.0
    # A smaller full-history term prevents a steady-only fit from selecting
    # a clearly wrong transient while keeping the requested observables central.
    transient_objective =
        mean((residual[:, 1:3] ./ 75.0).^2) / 3.0 +
        mean((residual[:, 7] ./ 50.0).^2) / 3.0 +
        mean((residual[:, 6] ./ 90.0).^2) / 3.0
    return (
        objective=0.8 * final_objective +
                  0.2 * transient_objective,
        side_final_sse=sum(abs2, side_final),
        T2_final_sse=T2_final^2,
        T3_final_sse=T3_final^2,
        side_final_bias_K=mean(side_final),
        side_final_rmse_K=sqrt(mean(abs2, side_final)),
        T2_final_bias_K=T2_final,
        T3_final_bias_K=T3_final,
        transient_side_sse=sum(abs2, residual[:, 1:3]),
        transient_side_n=length(residual[:, 1:3]),
        transient_T2_sse=sum(abs2, residual[:, 7]),
        transient_T2_n=length(residual[:, 7]),
        transient_T3_sse=sum(abs2, residual[:, 6]),
        transient_T3_n=length(residual[:, 6]),
    )
end

function _group_power_tuple_2D_v21(group, power)
    return group == 456 ? (power, 1.0, 1.0) :
           group == 304 ? (1.0, power, 1.0) :
                          (1.0, 1.0, power)
end

function _profile_candidate_2D_v21(
    area_cm2, group, power; max_points=31,
)
    p = parameters_2D_v21(
        effective_front_area_cm2=area_cm2,
        power_scales=_group_power_tuple_2D_v21(group, power),
    )
    metrics = NamedTuple[]
    for id in GROUP_IDS_2D_v21[group]
        inputs = case_inputs_2D_v21(
            id; cooling=false, max_points=max_points,
        )
        case = simulate_case_2D_v21(
            inputs, p; full_initial_data=false,
            reltol=5e-4, abstol=1e-4, dtmax=120.0,
        )
        push!(metrics, _calibration_metrics_2D_v21(case))
    end
    return (
        effective_front_area_cm2=area_cm2,
        group_kW_m2=group,
        power_scale=power,
        objective=mean(row.objective for row in metrics),
        side_final_rmse_K=sqrt(
            sum(row.side_final_sse for row in metrics) /
            (3 * length(metrics)),
        ),
        T2_final_rmse_K=sqrt(mean(
            row.T2_final_sse for row in metrics
        )),
        T3_final_rmse_K=sqrt(mean(
            row.T3_final_sse for row in metrics
        )),
        side_final_bias_K=mean(
            row.side_final_bias_K for row in metrics
        ),
        T2_final_bias_K=mean(
            row.T2_final_bias_K for row in metrics
        ),
        T3_final_bias_K=mean(
            row.T3_final_bias_K for row in metrics
        ),
        transient_side_rmse_K=sqrt(
            sum(row.transient_side_sse for row in metrics) /
            sum(row.transient_side_n for row in metrics),
        ),
        transient_T2_rmse_K=sqrt(
            sum(row.transient_T2_sse for row in metrics) /
            sum(row.transient_T2_n for row in metrics),
        ),
        transient_T3_rmse_K=sqrt(
            sum(row.transient_T3_sse for row in metrics) /
            sum(row.transient_T3_n for row in metrics),
        ),
    )
end

function calibrate_2D_v21(; max_points=31)
    mkpath(OUTPUT_DIR_2D_v21)
    candidates = [
        (area=area, group=group, power=power)
        for area in FRONT_AREA_GRID_CM2_2D_v21
        for group in (456, 304, 256)
        for power in POWER_GRID_2D_v21[group]
    ]
    rows = Vector{NamedTuple}(undef, length(candidates))
    Threads.@threads for index in eachindex(candidates)
        candidate = candidates[index]
        println(
            "v21 profile $index/$(length(candidates)): A=",
            candidate.area, " group=", candidate.group,
            " power=", candidate.power,
        )
        flush(stdout)
        rows[index] = _profile_candidate_2D_v21(
            candidate.area, candidate.group, candidate.power;
            max_points=max_points,
        )
    end
    _write_namedtuples_csv_2D_v21(
        joinpath(OUTPUT_DIR_2D_v21, "area_power_profile_2D_v21.csv"),
        rows,
    )

    area_rows = NamedTuple[]
    selections = Dict{Tuple{Float64,Int},NamedTuple}()
    for area in FRONT_AREA_GRID_CM2_2D_v21
        selected = NamedTuple[]
        for group in (456, 304, 256)
            subset = filter(
                row -> row.effective_front_area_cm2 == area &&
                       row.group_kW_m2 == group,
                rows,
            )
            best = subset[argmin(row.objective for row in subset)]
            selections[(area, group)] = best
            push!(selected, best)
        end
        push!(area_rows, (
            effective_front_area_cm2=area,
            objective=mean(row.objective for row in selected),
            power_scale_456=selected[1].power_scale,
            power_scale_304=selected[2].power_scale,
            power_scale_256=selected[3].power_scale,
            side_final_rmse_K=sqrt(mean(
                row.side_final_rmse_K^2 for row in selected
            )),
            T2_final_rmse_K=sqrt(mean(
                row.T2_final_rmse_K^2 for row in selected
            )),
            T3_final_rmse_K=sqrt(mean(
                row.T3_final_rmse_K^2 for row in selected
            )),
        ))
    end
    _write_namedtuples_csv_2D_v21(
        joinpath(OUTPUT_DIR_2D_v21, "area_profile_2D_v21.csv"),
        area_rows,
    )
    best_area = area_rows[argmin(row.objective for row in area_rows)]
    _write_namedtuples_csv_2D_v21(
        joinpath(OUTPUT_DIR_2D_v21, "selected_2D_v21.csv"),
        [best_area],
    )
    println("\nv21 selected coarse fit")
    println(best_area)
    return (profiles=rows, areas=area_rows, selected=best_area)
end

if abspath(PROGRAM_FILE) == @__FILE__
    points = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 31
    calibrate_2D_v21(; max_points=points)
end
