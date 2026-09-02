include("1D_v42.jl")
using GlobalSensitivity
using QuasiMonteCarlo
using CSV
using DataFrames

function run_sensitivity()
    println("[sensitivity_1D_v42] Starting sensitivity analysis...")
    
    # Define the 19 active parameters and their bounds
    param_names = [
        "A_Nu", "B_Re", "m_rec", "G_core_perim", "C_perim_eff", 
        "k_perim_ref", "beta_opt", "spill_capture", "beta_perim", "f_core_rear", 
        "flange_scale", "flange_cool_gain", "flange_cool_tau_s", "k_core_axial_scale", "C_rear_eff", 
        "G_recv_rear", "G_rear_tube", "G_rear_cavity", "G_rear_axial"
    ]
    
    # Indices in the full p array
    param_indices = [1, 2, 5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
    
    lb = lb_full_v42[param_indices]
    ub = ub_full_v42[param_indices]

    base_p = copy(pnew_v42)
    
    function obj_wrapper(x)
        p = copy(base_p)
        for (i, idx) in enumerate(param_indices)
            p[idx] = x[i]
        end
        try
            loss_h = loss_heating_v42(p, sim_key_heat; nodes=15)
            loss_c = loss_cooling_v42(p, sim_key_cool; nodes=15)
            
            loss = loss_h + loss_c
            if !isfinite(loss)
                return 1e6
            end
            return loss
        catch e
            return 1e6
        end
    end

    println("[sensitivity_1D_v42] Running Morris Sensitivity (may take a while)...")
    # Using Morris method with 100 trajectories, p=4 grid levels
    res = gsa(obj_wrapper, Morris(total_num_trajectory=50, num_trajectory=10, p_steps=fill(5, length(lb))), [[l, u] for (l, u) in zip(lb, ub)])
    
    println("[sensitivity_1D_v42] Analysis complete. Saving results.")
    
    df = DataFrame(
        Parameter = param_names,
        Index = param_indices,
        Mu_Star = res.means_star[1, :],
        Sigma = res.variances[1, :]
    )
    
    # Sort by importance (Mu_Star)
    sort!(df, :Mu_Star, rev=true)
    
    out_dir = joinpath("summaries", "1D_v42")
    mkpath(out_dir)
    out_path = joinpath(out_dir, "sensitivity_morris_1D_v42.csv")
    CSV.write(out_path, df)
    
    println("Top 10 most sensitive parameters:")
    for row in eachrow(first(df, 10))
        println("$(row.Parameter) (Idx $(row.Index)): Mu_Star = $(round(row.Mu_Star, sigdigits=4))")
    end
end

run_sensitivity()







