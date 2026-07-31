# ============================================================================
# plot_role_a.jl
#
# Generates "Role A" negative proof plots using the 2D_v20 (base) model.
# Demonstrates that a strictly mass-and-energy conserving 2D continuum model
# with a centered beam inherently predicts T_core > T_perim, contradicting 
# the experimental observation of T_perim > T_core deep in the receiver.
# ============================================================================

using Statistics
using Plots

include("run_2D_v20.jl")

function plot_negative_proof()
    outdir = joinpath(@__DIR__, "summaries", "Role_A_Plots")
    mkpath(outdir)
    
    # Run the base 2D model (v20) on a high power case (e.g. E71)
    id = "E71"
    println("Running negative proof simulation for $id...")
    p = parameters_2D_v20(
        power_scales=(1.19, 1.54, 0.94)
    )
    inputs = case_inputs_2D_v20(id; cooling=false, max_points=21)
    
    # We just run one case
    case = simulate_case_2D_v20(inputs, p; full_initial_data=false)
    
    # Extract final time window
    times = case.inputs.times
    threshold = first(times) + 0.90 * (last(times) - first(times))
    final = something(findfirst(time -> time >= threshold, times), length(times)):length(times)
    
    # In V20, SENSOR_NAMES_2D_v20 = ("T8", "T12", "T11", "T9", "T10", "T3", "T2")
    # T11 is deep perimeter (z=107mm)
    # T10 is deep core (z=107mm)
    
    idx_T11 = 3
    idx_T10 = 5
    
    # Model predictions
    T11_model = mean(case.model[final, idx_T11])
    T10_model = mean(case.model[final, idx_T10])
    
    # Experimental observations
    T11_exp = mean(case.observed[final, idx_T11])
    T10_exp = mean(case.observed[final, idx_T10])
    
    # Generate bar chart to highlight the inversion
    categories = ["Perimeter (T11)", "Core (T10)"]
    model_vals = [T11_model, T10_model]
    exp_vals = [T11_exp, T10_exp]
    
    p1 = bar(
        [categories categories],
        [model_vals exp_vals],
        labels=["2D Model (Centered)" "Experiment"],
        color=[:blue :red],
        ylabel="Deep Temperature at z=107mm (K)",
        title="Spatial Temperature Inversion (Negative Proof)",
        ylim=(300, max(maximum(model_vals), maximum(exp_vals)) + 150),
        legend=:topleft,
        bar_width=0.6,
        framestyle=:box
    )
    
    # Annotate contradiction
    annotate!(p1, 1.5, maximum(model_vals) + 120, text("Contradiction:\nExp: T_perim > T_core\nModel: T_core > T_perim", :black, :center, 10))
    
    # Combine and save
    fig = plot(p1, size=(600, 500), left_margin=5Plots.mm, bottom_margin=5Plots.mm)
    savefig(fig, joinpath(outdir, "Role_A_Negative_Proof.png"))
    println("Saved Role A negative proof plot to $outdir/Role_A_Negative_Proof.png")
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_negative_proof()
end
