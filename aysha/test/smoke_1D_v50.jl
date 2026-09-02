# test/smoke_1D_v50.jl
# Automated test suite for 1D_v50 Entire Converter Model, mesh sensitivity, and conservation audits.

using Test
include("../1D_v50.jl")

@testset "1D_v50 Basic Model & Parameter Tests" begin
    @test length(pnew_v50) == 28
    @test length(lb_full_v50) == 28
    @test length(ub_full_v50) == 28
    @test all(lb_full_v50 .<= pnew_v50 .<= ub_full_v50)

    # Optimizer bounds and active accessors must agree exactly at the seed and
    # at both advertised endpoints. This guards against the hidden clipping
    # that made part of the v49 search domain ineffective.
    @test core_perimeter_conductance_v50(pnew_v50) == pnew_v50[8]
    @test flange_scale_v50(pnew_v50) == pnew_v50[14]
    @test rear_tube_conductance_v50(pnew_v50) == pnew_v50[18]
    @test suction_heat_transfer_coefficient_v50(pnew_v50) == pnew_v50[23]
    for (index, accessor) in (
        (8, core_perimeter_conductance_v50),
        (14, flange_scale_v50),
        (18, rear_tube_conductance_v50),
        (23, suction_heat_transfer_coefficient_v50),
    )
        p_endpoint = copy(pnew_v50)
        p_endpoint[index] = lb_full_v50[index]
        @test accessor(p_endpoint) == lb_full_v50[index]
        p_endpoint[index] = ub_full_v50[index]
        @test accessor(p_endpoint) == ub_full_v50[index]
    end
    @test fit_stage_indices_v50(:observation) == collect(24:28)
    @test all(i in fit_stage_indices_v50(:full) for i in (5, 6, 7, 24, 28))
    
    # Check physical geometric properties
    @test isapprox(A_frt_v50, pi * R_core_v50^2, rtol=1e-5)
    @test isapprox(epsilon_v50, (100 * 0.0015^2) / A_frt_v50, rtol=1e-5)
    @test isapprox(A_solid_v50, (1.0 - epsilon_v50) * A_frt_v50, rtol=1e-5)
    @test isapprox(P_exchange_v50, 4 * 100 * 0.0015, rtol=1e-5)
    
    # Check participating assembly capacitance near 301 J/K
    C_total = participating_total_heat_capacity_v50(pnew_v50, 15)
    @test 250.0 < C_total < 350.0
    @test isapprox(C_total, 301.0, atol=25.0)
end

@testset "1D_v50 Observation Parameter Activity" begin
    baseline, _, _ = solve_case_v50(pnew_v50, "E67"; nodes=15)
    perturbations = [
        (24, 200.0, 6), # T3 convection
        (25, 0.10, 6),  # T3 radiation
        (26, 0.30, 6),  # T3 stem conduction
        (27, 0.40, 5),  # T10 observation
        (28, 0.40, 3),  # T11 observation
    ]
    for (parameter_index, perturbed_value, signal_index) in perturbations
        p_test = copy(pnew_v50)
        p_test[parameter_index] = perturbed_value
        changed, _, _ = solve_case_v50(p_test, "E67"; nodes=15)
        @test maximum(abs.(changed[:, signal_index] .- baseline[:, signal_index])) > 1.0e-6
    end
end

@testset "1D_v50 Objective Composition & Active Flow-Slope Loss" begin
    components = calibration_loss_components_v50(pnew_v50; nodes=15)
    @test components.error === nothing
    @test isfinite(components.total)
    @test components.slope_loss > 0.0
    @test isapprox(
        components.total,
        components.signal_loss + components.slope_loss +
        components.capacitance_regularization + components.power_regularization;
        rtol=1.0e-12,
    )
    @test isapprox(calibration_loss_v50(pnew_v50; nodes=15), components.total; rtol=1.0e-12)
end

@testset "1D_v50 Forward Simulation & Conservation Residual (Heating E67)" begin
    outputs, result, experiment = solve_case_v50(pnew_v50, "E67"; nodes=15)
    @test size(outputs, 1) == size(experiment, 1)
    @test size(outputs, 2) == 8
    @test all(isfinite, outputs)
    
    # Check exact First-Law energy conservation at endpoint
    conditions, time, _, _ = experimental_case_v50("E67")
    t_end = time[end]
    operating = operating_condition_v3(
        irradiance=linear_history_v3(time, fill(conditions[Io], length(time))),
        flow_lpm=linear_history_v3(time, observation(measurements, "E67", "_flow")),
        inlet_temperature=linear_history_v3(time, observation(measurements, "E67", "_Tin")),
        ambient_temperature=linear_history_v3(time, observation(measurements, "E67", "_Tamb")),
    )
    
    balance = compute_energy_balance_v50(result, operating, t_end, pnew_v50; nodes=15)
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

@testset "1D_v50 Mesh Sensitivity & Macroscopic HTC Convergence" begin
    # Forward simulation at N = 15, 25, 50 nodes for E67
    htc_list = Float64[]
    q_gas_list = Float64[]
    t_core_exit_list = Float64[]
    
    node_counts = [15, 25, 50]
    for N in node_counts
        _, res, _ = solve_case_v50(pnew_v50, "E67"; nodes=N)
        idx = length(res.time)
        htc_val = macroscopic_htc_v50(res, idx)
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

@testset "1D_v50 Forward Simulation (Cooling C69, C80 & C81)" begin
    # C69 (60 LPM forced cooling)
    outputs_c69, _, experiment_c69 = solve_case_v50(pnew_v50, "C69"; is_cooling=true, nodes=15)
    @test all(isfinite, outputs_c69)
    @test isapprox(outputs_c69[1, 6], experiment_c69[1, 6]; atol=1.0e-8)
    println("Cooling C69 solved successfully.")
    
    # C80 (0 LPM natural cooling)
    outputs_c80, _, _ = solve_case_v50(pnew_v50, "C80"; is_cooling=true, nodes=15)
    @test all(isfinite, outputs_c80)
    println("Cooling C80 (0 LPM natural cooling) solved successfully.")
    
    # Verify non-trivial cooling temperature decay
    @test outputs_c80[end, 4] < outputs_c80[1, 4] # T9 decreases
    @test outputs_c80[end, 1] < outputs_c80[1, 1] # T8 decreases

    # C81 is the third exported cooling case and must also remain executable.
    outputs_c81, _, _ = solve_case_v50(pnew_v50, "C81"; is_cooling=true, nodes=15)
    @test all(isfinite, outputs_c81)
end

println("\nAll 1D_v50 smoke tests passed!")
