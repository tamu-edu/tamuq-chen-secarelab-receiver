include("1D_v41.jl")
using GlobalSensitivity
using QuasiMonteCarlo
using CSV
using DataFrames

function run_sensitivity()
    println("[sensitivity_1D_v41] Starting sensitivity analysis...")
    
    # Define the 19 active parameters and their bounds
    param_names = [
        "A_Nu", "B_Re", "m_rec", "G_core_perim", "C_perim_eff", 
        "k_perim_ref", "beta_opt", "spill_capture", "beta_perim", "f_core_rear", 
        "flange_scale", "flange_cool_gain", "flange_cool_tau_s", "k_core_axial_scale", "C_rear_eff", 
        "G_recv_rear", "G_rear_tube", "G_rear_cavity", "G_rear_axial"
    ]
    
    # Indices in the full p array
    param_indices = [1, 2, 5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
    
    # Bounds based on lb_full_v41 and ub_full_v41
    lb = [
        0.5,    # A_Nu
        0.0,    # B_Re
        0.0,    # m_rec
        2.0,    # G_core_perim
        80.0,   # C_perim_eff
        1.0,    # k_perim_ref
        1.0,    # beta_opt
        0.0,    # spill_capture
        0.5,    # beta_perim
        0.95,   # f_core_rear
        0.05,   # flange_scale
        1.0,    # flange_cool_gain
        30.0,   # flange_cool_tau
        0.01,   # k_core_axial_scale
        40.0,   # C_rear_eff
        0.0,    # G_recv_rear
        0.0,    # G_rear_tube
        0.0,    # G_rear_cavity
        0.0     # G_rear_axial
    ]
    ub = [
        10.0,   # A_Nu
        0.6,    # B_Re
        2.0,    # m_rec
        20.0,   # G_core_perim
        200.0,  # C_perim_eff
        20.0,   # k_perim_ref
        2000.0, # beta_opt
        1.0,    # spill_capture
        300.0,  # beta_perim
        1.0,    # f_core_rear
        0.2,    # flange_scale
        20.0,   # flange_cool_gain
        250.0,  # flange_cool_tau
        0.1,    # k_core_axial_scale
        120.0,  # C_rear_eff
        5.0,    # G_recv_rear
        5.0,    # G_rear_tube
        15.0,   # G_rear_cavity
        15.0    # G_rear_axial
    ]

    base_p = copy(pnew_v41)
    
    function obj_wrapper(x)
        p = copy(base_p)
        for (i, idx) in enumerate(param_indices)
            p[idx] = x[i]
        end
        try
            loss_h = loss_heating_v41(p, sim_key_heat; nodes=15)
            loss_c = loss_cooling_v41(p, sim_key_cool; nodes=15)
            
            loss = loss_h + loss_c
            if !isfinite(loss)
                return 1e6
            end
            return loss
        catch e
            return 1e6
        end
    end

    println("[sensitivity_1D_v41] Running Morris Sensitivity (may take a while)...")
    # Using Morris method with 100 trajectories, p=4 grid levels
    res = gsa(obj_wrapper, Morris(total_num_trajectory=50, num_trajectory=10, p_steps=fill(5, length(lb))), [[l, u] for (l, u) in zip(lb, ub)])
    
    println("[sensitivity_1D_v41] Analysis complete. Saving results.")
    
    df = DataFrame(
        Parameter = param_names,
        Index = param_indices,
        Mu_Star = res.means_star[1, :],
        Sigma = res.variances[1, :]
    )
    
    # Sort by importance (Mu_Star)
    sort!(df, :Mu_Star, rev=true)
    
    out_dir = joinpath("summaries", "1D_v41")
    mkpath(out_dir)
    out_path = joinpath(out_dir, "sensitivity_morris_1D_v41.csv")
    CSV.write(out_path, df)
    
    println("Top 10 most sensitive parameters:")
    for row in eachrow(first(df, 10))
        println("$(row.Parameter) (Idx $(row.Index)): Mu_Star = $(round(row.Mu_Star, sigdigits=4))")
    end
end

run_sensitivity()
