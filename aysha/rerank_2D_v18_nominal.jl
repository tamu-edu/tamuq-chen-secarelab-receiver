# ============================================================================
# rerank_2D_v18_nominal.jl - nominal-mesh reranking of profiled structures
# ============================================================================

using CSV
using DataFrames

include("calibrate_2D_v18_staged.jl")

function _parse_candidate_2D_v18(label)
    match_result = match(
        r"^(beer_lambert|near_deep)_f([^_]+)_L([^_]+)_x([^_]+)$",
        String(label),
    )
    match_result === nothing &&
        error("cannot parse v18 candidate label: $label")
    return (
        source_model=Symbol(match_result.captures[1]),
        deep_fraction=parse(Float64, match_result.captures[2]),
        deep_length_m=1e-3 * parse(
            Float64, match_result.captures[3],
        ),
        exchange_multiplier=parse(
            Float64, match_result.captures[4],
        ),
    )
end

function rerank_2D_v18_nominal()
    path = joinpath(
        OUTPUT_DIR_2D_v18, "structural_profiled_2D_v18.csv",
    )
    isfile(path) || error(
        "run calibrate_2D_v18_staged.jl before nominal reranking",
    )
    profiles_frame = CSV.read(path, DataFrame)
    # Structural mesh transfer uses one interleaved case per power group.
    # All nine training cases already entered each nested C-mesh power fit.
    training_inputs = [
        case_inputs_2D_v18(id; max_points=61)
        for id in ("E69", "E74", "E79")
    ]
    results = NamedTuple[]
    rows = NamedTuple[]
    for (index, row) in enumerate(eachrow(profiles_frame))
        candidate = _parse_candidate_2D_v18(row.candidate)
        scales = (
            Float64(row.power_scale_456),
            Float64(row.power_scale_304),
            Float64(row.power_scale_256),
        )
        result = _evaluate_candidate_2D_v18(
            candidate, training_inputs;
            power_scales=scales, mesh=:nominal,
        )
        profile = (
            candidate=candidate, power_scales=scales,
            loss=result.loss.total,
        )
        push!(results, profile)
        push!(rows, (
            candidate=String(row.candidate),
            objective=result.loss.total,
            curve=result.loss.curve,
            level=result.loss.level,
            side_shape=result.loss.side_shape,
            effectiveness=result.loss.effectiveness,
            power_scale_456=scales[1],
            power_scale_304=scales[2],
            power_scale_256=scales[3],
        ))
        println(
            "v18 nominal rerank $index/$(nrow(profiles_frame)) ",
            row.candidate, " J=",
            round(result.loss.total; digits=4),
        )
        flush(stdout)
    end
    _write_namedtuples_csv_2D_v18(
        joinpath(
            OUTPUT_DIR_2D_v18,
            "structural_nominal_rerank_2D_v18.csv",
        ), rows,
    )
    ordered = sort(results; by=x -> x.loss)
    winner = ordered[1]
    runner = ordered[min(2, length(ordered))]
    felt = _profile_felt_2D_v18(winner)
    _write_namedtuples_csv_2D_v18(
        joinpath(OUTPUT_DIR_2D_v18, "felt_profile_nominal_winner_2D_v18.csv"),
        felt.rows,
    )
    open(joinpath(
        OUTPUT_DIR_2D_v18, "parameters_selected_2D_v18.txt",
    ), "w") do io
        println(io, "source_model=", winner.candidate.source_model)
        println(io, "deep_fraction=", winner.candidate.deep_fraction)
        println(io, "deep_length_m=", winner.candidate.deep_length_m)
        println(io, "exchange_multiplier=",
                winner.candidate.exchange_multiplier)
        println(io, "power_scale_456=", winner.power_scales[1])
        println(io, "power_scale_304=", winner.power_scales[2])
        println(io, "power_scale_256=", winner.power_scales[3])
        println(io, "felt_conductivity_scale=", felt.best.conductivity)
        println(io, "felt_heat_capacity_scale=", felt.best.capacity)
        println(io, "felt_contact_scale=0.30")
        println(io, "selection_mesh=nominal")
        println(io, "T9_T10_in_objective=false")
    end
    open(joinpath(
        OUTPUT_DIR_2D_v18, "runner_up_2D_v18.txt",
    ), "w") do io
        println(io, "source_model=", runner.candidate.source_model)
        println(io, "deep_fraction=", runner.candidate.deep_fraction)
        println(io, "deep_length_m=", runner.candidate.deep_length_m)
        println(io, "exchange_multiplier=",
                runner.candidate.exchange_multiplier)
        println(io, "power_scale_456=", runner.power_scales[1])
        println(io, "power_scale_304=", runner.power_scales[2])
        println(io, "power_scale_256=", runner.power_scales[3])
    end
    println("v18 nominal winner=", winner, " felt=", felt.best)
    println("v18 nominal runner=", runner)
    return (winner=winner, runner=runner, felt=felt.best)
end

if abspath(PROGRAM_FILE) == @__FILE__
    rerank_2D_v18_nominal()
end
