using Test

include(joinpath(@__DIR__, "..", "run_2D_v18.jl"))

@testset "2D v18 conservative source" begin
    baseline = parameters_2D_v18(
        source_model=:beer_lambert,
    )
    zero_deep = parameters_2D_v18(
        source_model=:near_deep, deep_fraction=0.0,
        deep_length_m=80e-3,
    )
    deep = parameters_2D_v18(
        source_model=:near_deep, deep_fraction=0.65,
        deep_length_m=100e-3,
    )
    g = V18.build_network_grid2D(baseline)
    wb = V18.source_weights2D(baseline, g)
    wz = V18.source_weights2D(zero_deep, g)
    invariant = V18.source_power_invariant2D(deep)
    @test maximum(abs.(wb .- wz)) < 1e-15
    @test invariant.absolute_error < 1e-14
    @test maximum(abs.(
        invariant.radial_reference .-
        invariant.radial_candidate
    )) < 1e-14
    @test all(V18.source_weights2D(deep, g) .>= 0.0)
end

@testset "2D v18 exact mesh family" begin
    coarse = V18.base_grid2D(parameters_2D_v18(mesh=:screen))
    nominal = V18.base_grid2D(parameters_2D_v18(mesh=:nominal))
    refined = V18.base_grid2D(parameters_2D_v18(mesh=:refined))
    @test (coarse.nz, nominal.nz, refined.nz) == (24, 48, 96)
    @test coarse.z_faces == nominal.z_faces[1:2:end]
    @test nominal.z_faces == refined.z_faces[1:2:end]
    @test coarse.r_faces[15:end] ==
          nominal.r_faces[15:2:end]
    @test nominal.r_faces[15:end] ==
          refined.r_faces[15:2:end]
    @test coarse.z_rear_faces ==
          nominal.z_rear_faces[1:2:end]
    @test nominal.z_rear_faces ==
          refined.z_rear_faces[1:2:end]
end

@testset "2D v18 exchange and energy" begin
    low = parameters_2D_v18(exchange_multiplier=0.40)
    high = parameters_2D_v18(exchange_multiplier=1.00)
    @test low.heat_transfer.heat_transfer_scale ≈ 0.72
    @test high.heat_transfer.heat_transfer_scale ≈ 1.80
    p = parameters_2D_v18(
        source_model=:near_deep, deep_fraction=0.6,
        deep_length_m=100e-3,
    )
    grid = V18.build_network_grid2D(p)
    layout = V18.V15._state_layout(grid)
    state = fill(700.0, layout.total)
    op = V18.OperatingCondition2D(
        irradiance=250000.0, flow_lpm=9.0,
        inlet_temperature=300.0, ambient_temperature=300.0,
    )
    ledger = V18.energy_rate_ledger2D(state, p, op, 0.0)
    @test abs(ledger.residual) < 1e-7
    @test ledger.source_power_error < 1e-14
end
