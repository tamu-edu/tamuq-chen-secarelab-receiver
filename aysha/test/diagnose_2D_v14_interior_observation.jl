using Statistics

include(joinpath(@__DIR__, "..", "validate_2D_v14_staged.jl"))

p = selected_parameters_2D_v14(
    nominal_mesh=false, screen_mesh=true,
)
rows = NamedTuple[]
for id in TRAIN_HEATING_2D_v14
    inputs = case_inputs_2D_v14(id; max_points=61)
    case = simulate_case_2D_v14(
        inputs, p; reltol=3e-4, abstol=1e-4, dtmax=120.0,
    )
    pred = case.predictions
    function required_fraction(observed, wall, gas)
        abs(wall - gas) < 1e-9 && return NaN
        return (observed - gas) / (wall - gas)
    end
    push!(rows, (
        simulation_id=id,
        wall_minus_gas_58_K=
            pred.T9_wall[end] - pred.T9_gas[end],
        wall_minus_gas_107_K=
            pred.T10_wall[end] - pred.T10_gas[end],
        required_wall_fraction_58=required_fraction(
            case.observed[end, 4],
            pred.T9_wall[end], pred.T9_gas[end],
        ),
        required_wall_fraction_107=required_fraction(
            case.observed[end, 5],
            pred.T10_wall[end], pred.T10_gas[end],
        ),
        model_side_minus_wall_58_K=
            case.model[end, 2] - pred.T9_wall[end],
        model_side_minus_gas_58_K=
            case.model[end, 2] - pred.T9_gas[end],
        observed_side_minus_interior_58_K=
            case.observed[end, 2] - case.observed[end, 4],
    ))
end
foreach(println, rows)
_write_namedtuples_csv_2D_v14(
    joinpath(
        OUTPUT_DIR_2D_v14,
        "interior_observation_diagnostic_2D_v14.csv",
    ), rows,
)
