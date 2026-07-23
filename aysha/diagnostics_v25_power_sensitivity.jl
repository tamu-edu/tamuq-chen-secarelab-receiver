include("1D_v25.jl")

const BASE_PARAMS_v25_POWER_SENSITIVITY = copy(pnew_v25)

objective_v25_power_sensitivity(p; nodes=11) =
    loss_heating_v25(p; nodes=nodes) + loss_cooling_v25(p; nodes=nodes)

println("case,scale_factor,scale_456,scale_304,scale_256,objective")

for factor in (0.80, 0.90, 1.00, 1.05, 1.10, 1.20)
    p = copy(BASE_PARAMS_v25_POWER_SENSITIVITY)
    p[7:9] .*= factor
    println(join(("all", factor, p[7], p[8], p[9], objective_v25_power_sensitivity(p)), ','))
end

for (index, label) in zip(7:9, ("456", "304", "256"))
    base = BASE_PARAMS_v25_POWER_SENSITIVITY[index]
    for factor in (0.80, 0.90, 1.00, 1.05, 1.10, 1.20, 1.40)
        p = copy(BASE_PARAMS_v25_POWER_SENSITIVITY)
        p[index] = base * factor
        println(join((label, factor, p[7], p[8], p[9], objective_v25_power_sensitivity(p)), ','))
    end
end
