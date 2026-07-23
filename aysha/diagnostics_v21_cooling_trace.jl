include("1D_v21.jl")

p = [
    4.9207360837689365,
    0.464177148000973,
    0.3333333333333333,
    1.0,
    0.0,
    1.0,
    2.4053684447462973,
    2.5,
    1.4621744961095464,
    28.162239066275774,
    208.18861322274748,
    1.7069790339469213,
    184.67243519237965,
    0.4856581989090208,
    0.37121222436610385,
    0.9992576201178448,
    0.10139749224843751,
    8.858688184209326,
    0.0,
    1.0,
]

for simulation_id in sim_key_cool
    model, result, experiment = solve_case_v21(p, simulation_id; is_cooling=true, nodes=11)
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
