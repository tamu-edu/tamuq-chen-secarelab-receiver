using OrdinaryDiffEq

include(joinpath(@__DIR__, "..", "2D_v8.jl"))
using .Receiver2D_v8

p = default_parameters2D()
op = OperatingCondition2D(
    irradiance=p.optics.scale_304 * 304000.0,
    flow_lpm=10.0,
    inlet_temperature=295.0,
    ambient_temperature=295.0,
)
times = [0.0, 300.0, 600.0]
algorithms = (
    ("Rodas5P-dense-FD", Rodas5P(autodiff=AutoFiniteDiff())),
    ("FBDF-dense-FD", FBDF(autodiff=AutoFiniteDiff())),
    ("Rosenbrock23-dense-FD", Rosenbrock23(autodiff=AutoFiniteDiff())),
)

reference = nothing
for (name, algorithm) in algorithms
    elapsed = @elapsed result = simulate2D(p, op, times; solver=algorithm)
    pred = sensor_predictions2D(result)
    values = [pred.T8[end], pred.T12[end], pred.T11[end], pred.T9[end],
              pred.T10[end], pred.T3[end], pred.T2[end]]
    if isnothing(reference)
        global reference = values
    end
    println(name, " elapsed_s=", round(elapsed, digits=3),
            " max_delta_K=", round(maximum(abs.(values .- reference)), digits=6))
end
