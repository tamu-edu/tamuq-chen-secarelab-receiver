include("1D_v23.jl")

p = [
    4.924815049912292,
    0.5162678835987395,
    0.3333333333333333,
    1.0,
    0.0,
    1.0,
    2.4155667071465876,
    2.497550759733723,
    1.1904357456021397,
    16.044477788186573,
    214.29323236142696,
    7.591183315111345,
    184.67243519237965,
    0.5251021463644217,
    1.8019923585149422,
    0.9993097154917425,
    0.1013817714666468,
    3.336402359922935,
    0.0,
    0.9997676098325551,
    0.0,
]

for simulation_id in sim_key_cool
    model, result, experiment = solve_case_v23(p, simulation_id; is_cooling=true, nodes=11)
    println(simulation_id)
    for (j, name) in enumerate(("T8", "T12", "T11", "T9", "T10", "T3", "T2"))
        count = min(4, size(model, 1))
        model_deltas = diff(model[:, j])
        experiment_deltas = diff(experiment[:, j])
        println(name,
                " model=", collect(model[1:count, j]),
                " experiment=", collect(experiment[1:count, j]),
                " first_delta_model=", model_deltas[1],
                " max_upturn_model=", maximum(max.(model_deltas, 0.0)),
                " first_delta_experiment=", experiment_deltas[1],
                " max_upturn_experiment=", maximum(max.(experiment_deltas, 0.0)))
    end
end


