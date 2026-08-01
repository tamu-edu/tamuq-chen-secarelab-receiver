include("1D_v41.jl"); println("Testing base pnew_v41 loss..."); try; loss = loss_heating_v41(pnew_v41, sim_key_heat; nodes=15); println("Loss = "); catch e; println(e); end
