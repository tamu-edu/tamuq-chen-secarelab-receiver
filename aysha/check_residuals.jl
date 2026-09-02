include("1D_v43.jl")

using Printf

p_opt = [2.087730202135983, 0.49940644166095516, 0.3333333333333333, 0.861956527420877, 0.2691043290536539, 1.0, 1.34, 1.58, 1.11, 11.000547307967606, 180.71299861548547, 6.344221043258331, 184.67243519237965, 0.4537519078039935, 12.607965629872641, 0.9992076236151953, 0.10145026438772616, 4.710520400565495, 121.23490607136135, 0.05077605326790914, 87.2262384009021, 0.14161292617715687, 2.188212572243105, 1.4096784900537245, 4.525763338199959, 0.00010056420362922036, 0.4104457219696818]

println("Evaluating SS for heating runs...")
println("--- STEADY STATE RESIDUALS (Mod - Exp) [K] ---")
println("Run   | Irrad | Flow | err(T_gas_out) | err(T2) | err(T8)")
for k in sim_key_heat
    outputs, result, experiment = solve_case_v43(p_opt, k; is_cooling=false, nodes=15, rear_nodes=5)
    
    # We can get irrad/flow directly from heating_runs!
    run_data = heating_runs[k]
    # heating_runs is loaded in 1D_v3.jl
    
    irrad = maximum(run_data.input.irradiance.data)
    flow = maximum(run_data.input.flow_lpm.data)
    
    # Model steady state
    mod_gas_out = result.gas.T_out
    mod_T2 = result.solid.T2
    mod_T8 = result.solid.T8
    
    # Exp steady state
    exp_gas_out = experiment.gas.T_out
    exp_T2 = experiment.solid.T2
    exp_T8 = experiment.solid.T8
    
    err_gas = mod_gas_out - exp_gas_out
    err_T2 = mod_T2 - exp_T2
    err_T8 = mod_T8 - exp_T8
    
    @printf("%-5s | %5.0f | %4.0f | %+14.2f | %+7.2f | %+7.2f\n", k, irrad/1000, flow, err_gas, err_T2, err_T8)
end
