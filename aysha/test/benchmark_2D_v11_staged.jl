include(joinpath(@__DIR__, "..", "calibrate_2D_v11_staged.jl"))

BLAS.set_num_threads(1)
base = staged_base_2D_v11()
power = [
    base.optics.scale_456,
    base.optics.scale_304,
    base.optics.scale_256,
]
p = staged_parameters_2D_v11(
    base, [110.0, 1.0, 1.0], power,
)
elapsed = @elapsed outputs = compact_training_predictions_2D_v11(p)
println("threads=$(Threads.nthreads())")
println("nine_case_elapsed_s=$elapsed")
println(
    "profile_objective=",
    mean(abs2, profile_residuals_2D_v11(outputs)),
)
