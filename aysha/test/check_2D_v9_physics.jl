using Test
using Statistics

include(joinpath(@__DIR__, "..", "run_2D_v9.jl"))

@testset "2D_v9 DP1 import and cold-t0 calibration" begin
    for id in vcat(IDs, IDs_cooling)
        data = id in IDs ? measurements : measurements_cooling
        dp1 = observation(data, id, "_DP1")
        flow = observation(data, id, "_flow")
        @test length(dp1) == length(flow)
        @test all(isfinite, dp1)
    end

    calibration = cold_t0_dp1_calibration_2D_v9()
    selected_ids = [row.simulation_id for row in calibration.selected]
    @test selected_ids == [
        "E67", "E68", "E70", "E72", "E74",
        "E75", "E76", "E78", "E80",
    ]
    @test calibration.intercept_mbar ≈ -0.614226030202630 rtol=1e-8
    @test calibration.slope_mbar_per_lpm ≈ 0.0455544997548929 rtol=1e-8
    @test calibration.r2 > 0.98
    @test calibration.rmse_mbar < 0.03
    @test calibration.hydraulic_resistance_scale ≈ 1.95171196637134 rtol=1e-6

    p = with_t0_hydraulics_2D_v9(default_parameters2D(), calibration)
    @test p.hydraulics.mass_flow_scale == 1.0
    @test p.hydraulics.minor_loss_coefficient == 0.0
    @test p.hydraulics.dp1_zero_offset_mbar ≈ calibration.intercept_mbar
    @test p.hydraulics.hydraulic_resistance_scale ≈
          calibration.hydraulic_resistance_scale
end
