using Test
using Statistics

include(joinpath(@__DIR__, "..", "2D_v9.jl"))
using .Receiver2D_v9

@testset "2D_v9 standard flow, velocity, pressure, and sensor smoke tests" begin
    p = default_parameters2D()
    grid = Receiver2D_v9.build_grid2D(p)

    rho_ref = standard_air_density2D(p.hydraulics)
    @test rho_ref ≈ 1.20018 rtol=5e-4
    @test standard_mass_flow2D(10.0, p.hydraulics) ≈
          rho_ref * 10.0 / 60000.0 rtol=1e-14
    @test air_density(600.0) < air_density(300.0)
    @test ideal_square_channel_dp_slope2D(p) ≈ 0.0233408 rtol=1e-4
    @test p.hydraulics.hydraulic_resistance_scale *
          ideal_square_channel_dp_slope2D(p) ≈ 0.0455545 rtol=1e-4

    sensors = Dict(
        :T8 => 760.0, :T12 => 750.0, :T11 => 620.0,
        :T9 => 730.0, :T10 => 600.0, :T3 => 550.0, :T2 => 320.0,
    )
    u0 = build_initial_state_2D(grid, p, sensors, 295.0)
    op = OperatingCondition2D(
        irradiance=0.0,
        flow_lpm=10.0,
        inlet_temperature=295.0,
        ambient_temperature=295.0,
    )
    result = simulate2D(p, op, [0.0, 1.0e-6]; initial_temperature=u0)
    pred = sensor_predictions2D(result)
    @test pred.T8[1] ≈ sensors[:T8] atol=1e-8
    @test pred.T12[1] ≈ sensors[:T12] atol=1e-8
    @test pred.T11[1] ≈ sensors[:T11] atol=1e-8
    @test pred.T9[1] ≈ sensors[:T9] atol=1e-8
    @test pred.T10[1] ≈ sensors[:T10] atol=1e-8

    @test size(result.gas_density) == (grid.nr_rec, grid.nz, 2)
    @test size(result.gas_velocity) == (grid.nr_rec, grid.nz, 2)
    @test size(result.gas_reynolds) == (grid.nr_rec, grid.nz, 2)
    @test size(result.cell_pressure_drop) == (grid.nr_rec, grid.nz, 2)
    @test all(result.gas_density .> 0.0)
    @test all(result.gas_velocity .> 0.0)
    @test all(result.gas_reynolds .> 0.0)
    @test all(result.cell_pressure_drop .> 0.0)
    @test all(result.receiver_pressure_drop_mbar .> 0.0)
    @test result.dp1_prediction_mbar ≈
          p.hydraulics.dp1_zero_offset_mbar .+
          result.receiver_pressure_drop_mbar

    mdot_from_cells = sum(
        result.gas_density[i, 1, 1] * result.gas_velocity[i, 1, 1] *
        grid.area_flow[i] for i in 1:grid.nr_rec
    )
    @test mdot_from_cells ≈ result.mass_flow_kg_s[1] rtol=1e-12

    i, j = 1, 1
    re_rhou = result.gas_density[i, j, 1] * result.gas_velocity[i, j, 1] *
              grid.dh / air_viscosity(
                  0.5 * (result.gas_temperature[i, j, 1] +
                         result.gas_temperature[i, j + 1, 1]),
              )
    @test re_rhou ≈ result.gas_reynolds[i, j, 1] rtol=1e-12

    function uniform_case(T, flow)
        operating = OperatingCondition2D(
            irradiance=0.0,
            flow_lpm=flow,
            inlet_temperature=T,
            ambient_temperature=T,
        )
        return simulate2D(p, operating, [0.0, 1.0e-6]; initial_temperature=T)
    end

    cold = uniform_case(300.0, 10.0)
    hot = uniform_case(700.0, 10.0)
    low_flow = uniform_case(300.0, 5.0)
    zero_flow = uniform_case(300.0, 0.0)
    @test mean(hot.gas_velocity[:, :, 1]) > mean(cold.gas_velocity[:, :, 1])
    @test hot.receiver_pressure_drop_mbar[1] > cold.receiver_pressure_drop_mbar[1]
    @test cold.receiver_pressure_drop_mbar[1] >
          low_flow.receiver_pressure_drop_mbar[1] > 0.0
    @test zero_flow.receiver_pressure_drop_mbar[1] == 0.0
    @test all(iszero, zero_flow.gas_velocity)
    @test all(iszero, zero_flow.gas_reynolds)
    @test all(iszero, zero_flow.cell_pressure_drop)

    h0 = p.hydraulics
    h_minor = HydraulicParameters2D(
        standard_pressure=h0.standard_pressure,
        standard_temperature=h0.standard_temperature,
        atmospheric_pressure=h0.atmospheric_pressure,
        mass_flow_scale=h0.mass_flow_scale,
        dp1_zero_offset_mbar=h0.dp1_zero_offset_mbar,
        hydraulic_resistance_scale=h0.hydraulic_resistance_scale,
        minor_loss_coefficient=1.0,
    )
    p_minor = ModelParameters2D(
        p.geometry, p.solid, p.heat_transfer, p.losses, p.optics, h_minor,
    )
    reference_op = OperatingCondition2D(
        irradiance=0.0,
        flow_lpm=10.0,
        inlet_temperature=h0.standard_temperature,
        ambient_temperature=h0.standard_temperature,
    )
    reference_base = simulate2D(
        p, reference_op, [0.0, 1.0e-6];
        initial_temperature=h0.standard_temperature,
    )
    reference_minor = simulate2D(
        p_minor, reference_op, [0.0, 1.0e-6];
        initial_temperature=h0.standard_temperature,
    )
    hot_minor_op = OperatingCondition2D(
        irradiance=0.0,
        flow_lpm=10.0,
        inlet_temperature=700.0,
        ambient_temperature=700.0,
    )
    hot_minor = simulate2D(
        p_minor, hot_minor_op, [0.0, 1.0e-6];
        initial_temperature=700.0,
    )
    @test reference_minor.receiver_pressure_drop_mbar[1] ≈
          reference_base.receiver_pressure_drop_mbar[1] atol=1e-12
    @test hot_minor.receiver_pressure_drop_mbar[1] >
          hot.receiver_pressure_drop_mbar[1]

    ledger = energy_rate_ledger2D(u0, p, op, 0.0)
    @test abs(ledger.residual) < 1e-8
    @test length(pack_parameters2D(p)) == 12
    @test unpack_parameters2D(pack_parameters2D(p), p).hydraulics == p.hydraulics
end
