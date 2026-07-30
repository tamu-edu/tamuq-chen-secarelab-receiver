include(joinpath(@__DIR__, "..", "run_2D_v12.jl"))

p = base_parameters_2D_v12(; nominal_mesh=false)
for id in ("E67", "E72", "E77")
    inputs = case_inputs_2D_v12(id; max_points=81)
    case = simulate_case_2D_v12(inputs, p)
    pred = case.predictions
    println((
        id=id,
        T12=pred.T12[end],
        T9_observed=case.observed[end, 4],
        T9_reported=pred.T9[end],
        T9_wall=pred.T9_wall[end],
        T9_gas=pred.T9_gas[end],
        T11=pred.T11[end],
        T10_observed=case.observed[end, 5],
        T10_reported=pred.T10[end],
        T10_wall=pred.T10_wall[end],
        T10_gas=pred.T10_gas[end],
    ))
end

