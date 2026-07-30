# ============================================================================
# confirm_2D_v18_refined.jl - winner/runner transfer on M and F meshes
# ============================================================================

using Statistics

include("calibrate_2D_v18_staged.jl")

function _read_keyvalues_2D_v18(path)
    values = Dict{String,String}()
    for line in eachline(path)
        occursin("=", line) || continue
        key, value = split(line, "="; limit=2)
        values[key] = value
    end
    return values
end

function _candidate_from_values_2D_v18(values)
    return (
        source_model=Symbol(values["source_model"]),
        deep_fraction=parse(Float64, values["deep_fraction"]),
        deep_length_m=parse(Float64, values["deep_length_m"]),
        exchange_multiplier=parse(
            Float64, values["exchange_multiplier"],
        ),
    )
end

function confirm_2D_v18_refined()
    selected_values = _read_keyvalues_2D_v18(joinpath(
        OUTPUT_DIR_2D_v18, "parameters_selected_2D_v18.txt",
    ))
    runner_values = _read_keyvalues_2D_v18(joinpath(
        OUTPUT_DIR_2D_v18, "runner_up_2D_v18.txt",
    ))
    number(key) = parse(Float64, selected_values[key])
    scales = (
        number("power_scale_456"),
        number("power_scale_304"),
        number("power_scale_256"),
    )
    conductivity = number("felt_conductivity_scale")
    capacity = number("felt_heat_capacity_scale")
    candidates = (
        (name="winner",
         value=_candidate_from_values_2D_v18(selected_values)),
        (name="runner",
         value=_candidate_from_values_2D_v18(runner_values)),
    )
    inputs = [
        case_inputs_2D_v18(id; max_points=61)
        for id in ("E69", "E74", "E79")
    ]
    simulations = Dict{Tuple{String,Symbol},Any}()
    objective_rows = NamedTuple[]
    for candidate in candidates
        for mesh in (:nominal, :refined)
            p = _parameter_tuple_2D_v18(
                candidate.value; mesh=mesh, power_scales=scales,
                felt_conductivity_scale=conductivity,
                felt_heat_capacity_scale=capacity,
            )
            cases = Vector{Any}(undef, length(inputs))
            Threads.@threads for index in eachindex(inputs)
                cases[index] = simulate_case_2D_v18(
                    inputs[index], p; reltol=1e-5,
                    abstol=1e-6, dtmax=45.0,
                )
            end
            simulations[(candidate.name, mesh)] = cases
            loss = set_objective_2D_v18(cases)
            push!(objective_rows, (
                candidate=candidate.name, mesh=String(mesh),
                objective=loss.total, curve=loss.curve,
                level=loss.level, side_shape=loss.side_shape,
                effectiveness=loss.effectiveness,
            ))
            println(
                "v18 refined confirmation ", candidate.name,
                " ", mesh, " J=", round(loss.total; digits=4),
            )
            flush(stdout)
        end
    end
    transfer_rows = NamedTuple[]
    primary = [1, 2, 3, 6, 7]
    for candidate in candidates
        medium = simulations[(candidate.name, :nominal)]
        fine = simulations[(candidate.name, :refined)]
        for index in eachindex(inputs)
            delta = medium[index].model[:, primary] .-
                    fine[index].model[:, primary]
            push!(transfer_rows, (
                candidate=candidate.name,
                simulation_id=inputs[index].id,
                history_rms_K=sqrt(mean(abs2, delta)),
                history_max_K=maximum(abs, delta),
                final_max_K=maximum(abs,
                    medium[index].model[end, primary] .-
                    fine[index].model[end, primary],
                ),
            ))
        end
    end
    _write_namedtuples_csv_2D_v18(joinpath(
        OUTPUT_DIR_2D_v18,
        "winner_runner_refined_objective_2D_v18.csv",
    ), objective_rows)
    _write_namedtuples_csv_2D_v18(joinpath(
        OUTPUT_DIR_2D_v18,
        "winner_runner_mesh_transfer_2D_v18.csv",
    ), transfer_rows)
    foreach(println, objective_rows)
    foreach(println, transfer_rows)
    return (objectives=objective_rows, transfer=transfer_rows)
end

if abspath(PROGRAM_FILE) == @__FILE__
    confirm_2D_v18_refined()
end
