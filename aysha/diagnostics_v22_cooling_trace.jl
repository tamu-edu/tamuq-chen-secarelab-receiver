include("1D_v22.jl")

p = [
    4.921574226427226,
    0.511742950159627,
    0.3333333333333333,
    1.0,
    0.0,
    1.0,
    2.4106297663693166,
    2.4976089017065792,
    1.2277904875626875,
    17.144680365658584,
    220.12006336173621,
    3.660125389206832,
    184.67243519237965,
    0.48091400158208236,
    0.37545813528806443,
    0.999331816726757,
    0.10137111965578832,
    3.5143023555974002,
    0.0,
    0.9997650895622375,
    0.010077580462147243,
]

for simulation_id in sim_key_cool
    model, result, experiment = solve_case_v22(p, simulation_id; is_cooling=true, nodes=11)
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

