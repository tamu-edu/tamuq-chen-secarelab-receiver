using Test
using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "calibrate_2D_v11_staged.jl"))

function read_staged_theta_2D_v11()
    table = CSV.read(
        joinpath(
            OUTPUT_DIR_2D_v11,
            "staged_parameters_2D_v11.csv",
        ),
        DataFrame,
    )
    values = Float64.(table.value)
    return values[1:3], values[4:6]
end

profile_theta, power_theta = read_staged_theta_2D_v11()
base = staged_base_2D_v11()
p = staged_parameters_2D_v11(base, profile_theta, power_theta)
outputs = compact_training_predictions_2D_v11(p)
profile_objective = mean(abs2, profile_residuals_2D_v11(outputs))
level_objective = mean(abs2, level_residuals_2D_v11(outputs))

grid = Receiver2D_v11.build_grid2D(p)
uniform_op = OperatingCondition2D(
    irradiance=0.0,
    flow_lpm=12.0,
    inlet_temperature=300.0,
    ambient_temperature=300.0,
)
ledger = energy_rate_ledger2D(
    fill(300.0, grid.nr_total * grid.nz + grid.nz_rear),
    p,
    uniform_op,
    0.0,
)
micro = simulate2D(
    p,
    uniform_op,
    [0.0, 1.0e-6];
    initial_temperature=300.0,
)
mesh = CSV.read(
    joinpath(
        OUTPUT_DIR_2D_v11,
        "staged_mesh_comparison_2D_v11.csv",
    ),
    DataFrame,
)

@testset "v11 staged calibration verification" begin
    @test isempty(
        intersect(
            Set(TRAIN_HEATING_2D_v11),
            Set(VALID_HEATING_2D_v11),
        ),
    )
    @test all(PROFILE_LOWER_2D_v11 .<= profile_theta) &&
          all(profile_theta .<= PROFILE_UPPER_2D_v11)
    @test all(POWER_LOWER_2D_v11 .<= power_theta) &&
          all(power_theta .<= POWER_UPPER_2D_v11)
    @test all(isfinite, profile_theta)
    @test all(isfinite, power_theta)
    @test isfinite(profile_objective)
    @test isfinite(level_objective)
    @test profile_objective < 5.193912960052663
    @test abs(ledger.residual) < 1.0e-8
    @test sum(micro.ring_mass_flow_kg_s[:, end]) ≈
          micro.mass_flow_kg_s[end] rtol=1.0e-12
    @test micro.equal_pressure_relative_error[end] < 1.0e-5
    @test maximum(mesh.max_abs_sensor_delta_K) < 3.0
end

println("final_training_profile_objective=$profile_objective")
println("final_training_level_objective=$level_objective")
println("energy_residual_W=$(ledger.residual)")
println(
    "max_nominal_mesh_sensor_delta_K=",
    maximum(mesh.max_abs_sensor_delta_K),
)
