using Test

include(joinpath(@__DIR__, "..", "run_2D_v19.jl"))

@testset "2D v19 integrated UA invariants" begin
    for mesh in (:screen, :nominal, :refined)
        p = parameters_2D_v19(mesh=mesh)
        inventory = V19.network_inventory2D(p)
        @test inventory.group_count == 15
        @test inventory.multiplicity_sum == 100
        diagnostic = V19.integrated_ua_diagnostic2D(p)
        @test diagnostic.relative_ua_error < 1e-12
        @test diagnostic.ua_W_K > 0.0
    end
    pzero = parameters_2D_v19(nu_prefactor=0.0)
    zero = V19.integrated_ua_diagnostic2D(pzero)
    @test zero.ua_W_K == 0.0
    @test zero.outlet_temperature_K ≈ 300.0 atol=1e-10
    invariant = V19.source_power_invariant2D(
        parameters_2D_v19(
            source_model=:near_deep,
            deep_fraction=0.6, deep_length_m=0.05,
        ),
    )
    @test invariant.absolute_error < 1e-12
end

@testset "2D v19 case and probe" begin
    p = parameters_2D_v19(mesh=:screen)
    input = case_inputs_2D_v19("E72"; max_points=5)
    case = simulate_case_2D_v19(
        input, p; reltol=5e-4, abstol=1e-4, dtmax=120.0,
    )
    @test all(isfinite, case.model)
    @test case.result.mass_flow_kg_s[end] > 0.0
    @test maximum(case.result.equal_pressure_relative_error) < 1e-4
    @test all(case.result.flow_solver_converged)
    @test maximum(case.result.gas_reference_relative_error) < 1e-4
    grid = V19.build_network_grid2D(p)
    final = length(case.result.time)
    perimeter = 4.0 * p.geometry.channel_width
    orbit_ua = [
        sum(
            case.result.heat_transfer_coefficient[group, j, final] *
            perimeter * grid.base_grid.dz_cells[j]
            for j in 1:grid.base_grid.nz
        )
        for group in 1:grid.group_count
    ]
    @test maximum(orbit_ua) - minimum(orbit_ua) < 1e-10
    @test isapprox(
        sum(grid.multiplicity .* orbit_ua),
        100.0 * first(orbit_ua); rtol=1e-12,
    )
    @test isapprox(
        sum(
            grid.multiplicity .*
            case.result.channel_mass_flow_kg_s[:, final]
        ),
        case.result.mass_flow_kg_s[final]; rtol=1e-10,
    )

    cooling_input = case_inputs_2D_v19(
        "C81"; cooling=true, max_points=5,
    )
    cooling = simulate_case_2D_v19(
        cooling_input, p;
        reltol=5e-4, abstol=1e-4, dtmax=120.0,
    )
    @test cooling.model[1, 6] ≈ cooling.observed[1, 6] atol=1e-10
    @test all(isfinite, cooling.model)
end

@testset "2D v19 nominal energy ledger" begin
    p = parameters_2D_v19(mesh=:nominal)
    input = case_inputs_2D_v19("E72"; max_points=3)
    case = simulate_case_2D_v19(
        input, p; reltol=5e-4, abstol=1e-4, dtmax=120.0,
    )
    ledger = V19.energy_rate_ledger2D(
        case.result.ode_solution.u[end], p,
        case.result.operating_condition,
        case.result.time[end],
    )
    @test abs(ledger.residual) < 1e-8
    audit = V19.exchange_energy_audit2D(
        case.result.ode_solution.u[end], p,
        case.result.operating_condition,
        case.result.time[end],
    )
    @test audit.relative_residual < 1e-10
end
