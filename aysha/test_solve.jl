include("1D_v41.jl"); model, result, experiment = solve_case_v41(pnew_v41, sim_key_heat[1]; is_cooling=false, nodes=15); println("Success? ", successful_retcode(result))
