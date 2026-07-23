include("1D_v22.jl")

base = [
    4.928347639281704,
    0.5116375454619312,
    0.3333333333333333,
    1.0,
    0.0,
    1.0,
    2.381343080302037,
    2.497597482905488,
    1.324980204605796,
    18.571880646047624,
    210.5771105999091,
    3.04859988090074,
    184.67243519237965,
    0.4977451571213417,
    0.3683475947014179,
    0.9994425303007114,
    0.10140036970180075,
    3.567973794774853,
    0.0,
    1.0,
    0.3037437114125107,
]

println("core_axial_scale,max_T10_upturn_K,C69_T10_first_delta_K,C80_T10_first_delta_K,C81_T10_first_delta_K")
for scale in 0.0:0.025:0.40
    p = copy(base)
    p[21] = scale
    first_deltas = Float64[]
    max_upturns = Float64[]
    for simulation_id in sim_key_cool
        model, _, _ = solve_case_v22(p, simulation_id; is_cooling=true, nodes=11)
        deltas = diff(model[:, 5])
        push!(first_deltas, deltas[1])
        push!(max_upturns, maximum(max.(deltas, 0.0)))
    end
    println(join((scale, maximum(max_upturns), first_deltas...), ','))
end
