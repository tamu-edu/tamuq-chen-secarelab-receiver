include("calibrate_2D_v12_staged.jl")

result = fit_cooling_stage_2D_v12(
    physical_candidates=PHYSICAL_REFINEMENT_2D_v12,
    output_tag="cooling_refined",
)
println("selected refined cooling index = ", result.selected_index)
println("selected refined cooling loss = ", result.loss)
