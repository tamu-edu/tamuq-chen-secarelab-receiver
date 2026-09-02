# ============================================================================
# smoke_1D_v46.jl - Fast verification test for 1D_v46 Entire Converter Model
# ============================================================================

using Test
include("../1D_v46.jl")

println("Running smoke tests for 1D_v46...")

@testset "1D_v46 Basic Model & Parameter Tests" begin
    @test length(pnew_v46) == 23
    @test length(lb_full_v46) == 23
    @test length(ub_full_v46) == 23
    @test all(lb_full_v46 .<= pnew_v46 .<= ub_full_v46)
    @test NU_FLOOR_v46 == 0.0
end

@testset "1D_v46 Forward Simulation & Conservation Residual (Heating E67)" begin
    p = copy(pnew_v46)
    outputs, result, experiment = solve_case_v46(
        p, :E67;
        is_cooling=false,
        nodes=15,
        reltol=1e-5,
        abstol=1e-6,
        dtmax=30.0,
    )
    @test size(outputs, 1) == length(result.time)
    @test size(outputs, 2) == 8 # 7 signals + mean HTC
    @test all(outputs[:, 1:7] .> 250.0) # all temperatures physically positive (> 250 K)
    @test all(outputs[:, 1:7] .< 2000.0)
    
    # Check energy conservation at steady state (end of run)
    t_end = result.time[end]
    conditions, time_vec, _, data = experimental_case_v46("E67"; is_cooling=false)
    flow = observation(data, "E67", "_flow")
    Tin = observation(data, "E67", "_Tin")
    ambient = observation(data, "E67", "_Tamb")
    irradiance = fill(conditions[Io], length(time_vec))
    operating = operating_condition_v3(
        irradiance=linear_history_v3(time_vec, irradiance),
        flow_lpm=linear_history_v3(time_vec, flow),
        inlet_temperature=linear_history_v3(time_vec, Tin),
        ambient_temperature=linear_history_v3(time_vec, ambient),
    )
    
    balance = compute_energy_balance_v46(result, operating, t_end, p; nodes=15)
    println("E67 Energy Balance at t = $(t_end) s:")
    println("  Delivered Solar Input : $(round(balance.delivered_W, digits=2)) W")
    println("  Gas Enthalpy Gain     : $(round(balance.gas_heat_W, digits=2)) W")
    println("  Front Radiative Loss  : $(round(balance.front_rad_loss_W, digits=2)) W")
    println("  Cavity Loss to Amb    : $(round(balance.cavity_amb_loss_W, digits=2)) W")
    println("  Flange Loss to Water  : $(round(balance.flange_loss_W, digits=2)) W")
    println("  Sensible Storage Rate : $(round(balance.dE_stored_dt_W, digits=2)) W")
    println("  First-Law Residual    : $(round(balance.instantaneous_residual_W, digits=6)) W")
    
    # Exact instantaneous conservation
    @test abs(balance.instantaneous_residual_W) < 1e-4
end

@testset "1D_v46 Forward Simulation (Cooling C69)" begin
    p = copy(pnew_v46)
    outputs, result, experiment = solve_case_v46(
        p, :C69;
        is_cooling=true,
        nodes=15,
        reltol=1e-5,
        abstol=1e-6,
        dtmax=30.0,
    )
    @test size(outputs, 1) == length(result.time)
    @test size(outputs, 2) == 8
    @test all(outputs[:, 1:7] .> 250.0)
    println("Cooling C69 solved successfully.")
end

println("All 1D_v46 smoke tests passed!")
