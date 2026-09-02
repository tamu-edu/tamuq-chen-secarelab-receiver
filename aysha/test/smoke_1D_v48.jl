# test/smoke_1D_v48.jl
# Automated test suite for 1D_v48 Entire Converter Model, mesh sensitivity, and conservation audits.

using Test
include("../1D_v48.jl")

@testset "1D_v48 Basic Model & Parameter Tests" begin
    @test length(pnew_v48) == 23
    @test length(lb_full_v48) == 23
    @test length(ub_full_v48) == 23
    @test all(lb_full_v48 .<= pnew_v48 .<= ub_full_v48)
    
    # Check physical geometric properties
    @test isapprox(A_frt_v48, pi * R_core_v48^2, rtol=1e-5)
    @test isapprox(epsilon_v48, (100 * 0.0015^2) / A_frt_v48, rtol=1e-5)
    @test isapprox(A_solid_v48, (1.0 - epsilon_v48) * A_frt_v48, rtol=1e-5)
    @test isapprox(P_exchange_v48, 4 * 100 * 0.0015, rtol=1e-5)
    
    # Check participating assembly capacitance near 301 J/K
    C_total = participating_total_heat_capacity_v48(pnew_v48, 15)
    @test 250.0 < C_total < 350.0
    @test isapprox(C_total, 301.0, atol=25.0)
end

@testset "1D_v48 Forward Simulation & Conservation Residual (Heating E67)" begin
    outputs, result, experiment = solve_case_v48(pnew_v48, "E67"; nodes=15)
    @test size(outputs, 1) == size(experiment, 1)
    @test size(outputs, 2) == 8
    @test all(isfinite, outputs)
    
    # Check exact First-Law energy conservation at endpoint
    conditions, time, _, _ = experimental_case_v48("E67")
    t_end = time[end]
    operating = operating_condition_v3(
        irradiance=linear_history_v3(time, fill(conditions[Io], length(time))),
        flow_lpm=linear_history_v3(time, observation(measurements, "E67", "_flow")),
        inlet_temperature=linear_history_v3(time, observation(measurements, "E67", "_Tin")),
        ambient_temperature=linear_history_v3(time, observation(measurements, "E67", "_Tamb")),
    )
    
    balance = compute_energy_balance_v48(result, operating, t_end, pnew_v48; nodes=15)
    println("E67 Energy Balance at t = $(t_end) s (nodes=15):")
    println("  Delivered Solar Input : $(round(balance.delivered_W, digits=2)) W")
    println("  Suction Gas Preheat   : $(round(balance.suction_gas_heat_W, digits=2)) W")
    println("  Core Gas Heat         : $(round(balance.core_gas_heat_W, digits=2)) W")
    println("  Rear Tube Gas Heat    : $(round(balance.rear_tube_gas_heat_W, digits=2)) W")
    println("  Total Gas Enthalpy    : $(round(balance.gas_heat_W, digits=2)) W")
    println("  Front Radiative Loss  : $(round(balance.front_rad_loss_W, digits=2)) W")
    println("  Cavity Loss to Amb    : $(round(balance.cavity_amb_loss_W, digits=2)) W")
    println("  Flange Loss to Water  : $(round(balance.flange_loss_W, digits=2)) W")
    println("  Sensible Storage Rate : $(round(balance.dE_stored_dt_W, digits=2)) W")
    println("  First-Law Residual    : $(round(balance.instantaneous_residual_W, digits=6)) W")
    println("  Macroscopic LMTD HTC  : $(round(balance.macro_htc, digits=2)) W/m^2 K")
    
    # Residual must be machine zero (< 1e-10 W)
    @test abs(balance.instantaneous_residual_W) < 1.0e-10
    # Core gas heat must be positive under nominal flow
    @test balance.core_gas_heat_W > 0.0
end

@testset "1D_v48 Mesh Sensitivity & Macroscopic HTC Convergence" begin
    # Forward simulation at N = 15, 25, 50 nodes for E67
    htc_list = Float64[]
    q_gas_list = Float64[]
    t_core_exit_list = Float64[]
    
    node_counts = [15, 25, 50]
    for N in node_counts
        _, res, _ = solve_case_v48(pnew_v48, "E67"; nodes=N)
        idx = length(res.time)
        htc_val = macroscopic_htc_v48(res, idx)
        push!(htc_list, htc_val)
        push!(q_gas_list, sum(res.receiver_gas_heat[:, idx]))
        push!(t_core_exit_list, res.core_temperature[end, idx])
    end
    
    println("\nE67 Mesh Convergence Results across N in $(node_counts):")
    for (k, N) in enumerate(node_counts)
        println("  N = $(N): Macro HTC = $(round(htc_list[k], digits=2)) W/m^2K | Core Gas Heat = $(round(q_gas_list[k], digits=2)) W | Core Exit T = $(round(t_core_exit_list[k], digits=2)) K")
    end
    
    # Check that macroscopic HTC is stable within 15% and Core Gas Heat within 2% across N=15, 25, 50
    htc_rel_diff = abs(htc_list[end] - htc_list[1]) / htc_list[1]
    q_rel_diff = abs(q_gas_list[end] - q_gas_list[1]) / q_gas_list[1]
    println("  Macro HTC max relative difference (N=15 to 50): $(round(htc_rel_diff * 100, digits=2))%")
    println("  Core Gas Heat max relative difference (N=15 to 50): $(round(q_rel_diff * 100, digits=2))%")
    
    @test htc_rel_diff < 0.15
    @test q_rel_diff < 0.02
end

@testset "1D_v48 Forward Simulation (Cooling C69 & C80)" begin
    # C69 (60 LPM forced cooling)
    outputs_c69, _, _ = solve_case_v48(pnew_v48, "C69"; is_cooling=true, nodes=15)
    @test all(isfinite, outputs_c69)
    println("Cooling C69 solved successfully.")
    
    # C80 (0 LPM natural cooling)
    outputs_c80, _, _ = solve_case_v48(pnew_v48, "C80"; is_cooling=true, nodes=15)
    @test all(isfinite, outputs_c80)
    println("Cooling C80 (0 LPM natural cooling) solved successfully.")
    
    # Verify non-trivial cooling temperature decay
    @test outputs_c80[end, 4] < outputs_c80[1, 4] # T9 decreases
    @test outputs_c80[end, 1] < outputs_c80[1, 1] # T8 decreases
end

println("\nAll 1D_v48 smoke tests passed!")
