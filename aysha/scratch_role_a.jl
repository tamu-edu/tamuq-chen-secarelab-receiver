include("plot_2D_v20_gatefree_full.jl")
println("Loading complete.")
p = _gatefree_plot_parameters_2D_v20()
id = HEATING_IDS_2D_v20[1]
println("Running case: ", id)
inputs = case_inputs_2D_v20(id; cooling=false, max_points=21)
case = simulate_case_2D_v20(inputs, p; initialization=:ambient, reltol=5e-4, abstol=1e-4, dtmax=120.0)

println(fieldnames(typeof(case.result)))
println(fieldnames(typeof(case.result.base_result)))
println(fieldnames(typeof(case.result.base_result.base_result)))

# find how to access solid temperatures:
if hasproperty(case.result, :diagnostics)
    println("diagnostics keys: ", keys(case.result.diagnostics))
end

if hasproperty(case.result, :solid_temperature)
    println("solid_temperature size: ", size(case.result.solid_temperature))
end
