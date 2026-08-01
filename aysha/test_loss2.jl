include("1D_v41.jl"); try; loss_heating_v41(pnew_v41, sim_key_heat; nodes=15); catch e; println(e); end
