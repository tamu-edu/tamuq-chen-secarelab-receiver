# ============================================================================
# 2D_v21.jl
#
# Phase 2 Reduced-Order Model:
#   * Inherits v20 exact enthalpy and identifiability structures
#   * Introduces Phase2Parameters2D for missing physics:
#       - Explicit Spillage Heating (Orbit 15)
#       - Flow Maldistribution (core_preference parameter)
#   * Drops T3 from the calibration objective by default
# ============================================================================

module Receiver2D_v21

include("2D_v20.jl")
include("2D_v21_Properties.jl")
const V20 = Receiver2D_v20
const V21_Props = Receiver2D_v21_Properties
const V19 = V20.V19
const V18 = V19.V18
const V17 = V19.V17
const V16 = V19.V16
const V15 = V19.V15
const V14 = V19.V14
const V12 = V19.V12
const V11 = V19.V11

using DifferentialEquations
using OrdinaryDiffEq
using SciMLBase

export Phase2Parameters2D, ModelParameters2D_v21
export simulate2D_v21, apply_v21_property_fixes!

function apply_v21_property_fixes!()
    @eval V11 begin
        # 1. Corrected air properties
        air_conductivity(T) = Main.Receiver2D_v21.V21_Props.air_conductivity_v21(T)
        air_heat_capacity(T) = Main.Receiver2D_v21.V21_Props.air_heat_capacity_v21(T)
        air_viscosity(T) = Main.Receiver2D_v21.V21_Props.air_viscosity_v21(T)
        air_density(T, p=101325.0) = Main.Receiver2D_v21.V21_Props.air_density_v21(T, p)
        
        # 2. Corrected SiC properties
        sic_conductivity(T) = Main.Receiver2D_v21.V21_Props.sic_conductivity_v21(T)
        sic_heat_capacity(T) = Main.Receiver2D_v21.V21_Props.sic_heat_capacity_v21(T)
        
        # 3. Corrected felt properties
        felt_conductivity_temp(T, ref_k=0.06) = Main.Receiver2D_v21.V21_Props.felt_conductivity_v21(T)
        
        # Override keyword constructors to fix constants (rad_extinction_coeff, felt_heat_capacity, fully_developed_nusselt)
        function SolidProperties2D(;
            density = 2150.0,
            radial_conductivity_scale = 0.05,
            axial_conductivity_scale = 1.00,
            rad_extinction_coeff = 400.0, # FIXED from 50.0
            felt_conductivity_ref = 0.06,
            felt_density = 140.0,
            felt_heat_capacity = 1000.0, # FIXED from 1360.0
            casing_conductivity = 205.0,
            casing_density = 2700.0,
            casing_heat_capacity = 900.0,
            receiver_felt_contact_resistance = 1.0e-3
        )
            return SolidProperties2D(density, radial_conductivity_scale, axial_conductivity_scale, rad_extinction_coeff, felt_conductivity_ref, felt_density, felt_heat_capacity, casing_conductivity, casing_density, casing_heat_capacity, receiver_felt_contact_resistance)
        end
        
        function HeatTransferParameters2D(;
            coefficient = 1.0e-3,
            reynolds_exponent = 1.440,
            prandtl_exponent = 0.333,
            minimum_nusselt = 0.0,
            c_radial_flow = 0.10,
            front_coefficient = 0.0,
            front_reynolds_exponent = 0.5,
            nusselt_model = :graetz,
            fully_developed_nusselt = 2.98, # FIXED from 3.61
            temperature_correction_exponent = 0.0,
            heat_transfer_scale = 1.0
        )
            return HeatTransferParameters2D(coefficient, reynolds_exponent, prandtl_exponent, minimum_nusselt, c_radial_flow, front_coefficient, front_reynolds_exponent, nusselt_model, fully_developed_nusselt, temperature_correction_exponent, heat_transfer_scale)
        end
    end
end

Base.@kwdef struct Phase2Parameters2D
    spillage_power_W::Float64 = 0.0
    core_preference::Float64 = 0.0
    spillage_axial_spread_m::Float64 = 5.0e-3
end

Base.@kwdef struct ModelParameters2D_v21
    base::V20.ModelParameters2D = V20.default_parameters2D()
    phase2::Phase2Parameters2D = Phase2Parameters2D()
end

function Base.getproperty(p::ModelParameters2D_v21, name::Symbol)
    name === :base && return getfield(p, :base)
    name === :phase2 && return getfield(p, :phase2)
    return getproperty(getfield(p, :base), name)
end

function _solve_channel_flows2D_v21!(
    work, mdot_total, inlet, channel_temperature,
    p14, exchange, grid, phase2
)
    ng = grid.group_count
    
    # If core_preference > 0.0, we manually allocate flow towards the core
    # and bypass the equal pressure solver for this simple reduced order model
    if phase2.core_preference > 1e-6
        pref = phase2.core_preference
        weights = zeros(ng)
        # Orbit 1 is core, Orbit 15 is perimeter. 
        # A simple linear weighting where outer orbits get less flow.
        for group in 1:ng
            r_frac = (group - 1) / (ng - 1) # 0 at core, 1 at perimeter
            weights[group] = 1.0 - pref * r_frac
        end
        weights .= max.(weights, 0.01) # Minimum 1% flow to avoid div by zero
        
        weighted_multiplicity = grid.multiplicity .* weights
        work.channel_mdot .= mdot_total .* weights ./ sum(weighted_multiplicity)
        
        # We still need to march the channels to compute heat transfer
        gas_reference = inlet
        wall_profile = V19._mean_wall_profile2D(channel_temperature, grid)
        average = V19._average_nusselt2D(
            mdot_total / sum(grid.multiplicity),
            gas_reference, wall_profile, p14, exchange, grid,
        )
        for group in 1:ng
            V20._march_channel_group2D!(
                work, group, work.channel_mdot[group], inlet,
                channel_temperature, p14, exchange, grid,
                average.ua_W_K,
            )
        end
        
        # We skip the iterative pressure solver
        work.common_dp = 0.0
        work.equal_pressure_error = 0.0
        return (gas_reference_relative_error=0.0, iterations=1, converged=true)
    else
        # Fall back to v20 common-pressure solver
        return V20._solve_channel_flows2D!(
            work, mdot_total, inlet, channel_temperature,
            p14, exchange, grid
        )
    end
end

function _gas_profile_enthalpy2D_v21!(
    work, channel_temperature, tube_temperature, time,
    p14, exchange, op, grid, phase2
)
    bg = grid.base_grid
    flow_lpm = max(0.0, Float64(op.flow_lpm(time)))
    inlet = Float64(op.inlet_temperature(time))
    mdot_total = V11.standard_mass_flow2D(flow_lpm, p14.hydraulics)
    
    if mdot_total <= eps()
        fill!(work.qgas, 0.0)
        fill!(work.h, 0.0)
        fill!(work.front_qgas, 0.0)
        fill!(work.front_h, 0.0)
        fill!(work.velocity, 0.0)
        fill!(work.reynolds, 0.0)
        fill!(work.prandtl, 0.0)
        fill!(work.graetz, 0.0)
        fill!(work.nusselt, 0.0)
        fill!(work.dp_cell, 0.0)
        fill!(work.channel_mdot, 0.0)
        fill!(work.channel_dp, 0.0)
        fill!(work.groove_dp, 0.0)
        fill!(work.qgas_rear, 0.0)
        fill!(work.hrear, 0.0)
        work.common_dp = 0.0
        work.equal_pressure_error = 0.0
        work.gas .= inlet
        work.density .= V11.air_density(inlet, p14.hydraulics.atmospheric_pressure)
        work.gas_rear .= inlet
        return (gas_reference_relative_error=0.0, iterations=0, converged=true)
    end
    
    convergence = _solve_channel_flows2D_v21!(
        work, mdot_total, inlet, channel_temperature,
        p14, exchange, grid, phase2
    )
    work.gas_rear[1] = V20._mixed_outlet2D(work, grid)

    radius = p14.geometry.rear_tube_inner_radius
    dh = 2.0 * radius
    area = pi * radius^2
    perimeter = 2.0 * pi * radius
    for j in 1:bg.nz_rear
        Tfilm = V11.clamp_temperature(0.5 * (tube_temperature[j] + work.gas_rear[j]))
        cp = V11.air_heat_capacity(Tfilm)
        kf = V11.air_conductivity(Tfilm)
        mu = V11.air_viscosity(Tfilm)
        reynolds = mdot_total * dh / max(area * mu, eps())
        prandtl = cp * mu / kf
        nusselt = reynolds < 2300.0 ? 3.66 : 0.023 * reynolds^0.8 * prandtl^0.4
        work.hrear[j] = nusselt * kf / dh
        ua = work.hrear[j] * perimeter * bg.dz_rear
        exchange_cell = V20._enthalpy_cell_exchange2D(
            work.gas_rear[j], tube_temperature[j], ua, mdot_total,
        )
        work.gas_rear[j + 1] = exchange_cell.temperature
        work.qgas_rear[j] = exchange_cell.heat_W
    end
    return convergence
end

function _v21_ode2D!(du, u, context, time)
    V17._casing_flange_ode2D!(du, u, context.casing_context, time)
    
    old = V20._gas_power_snapshot2D(context.work)
    grid = context.grid
    bg = grid.base_grid
    layout = context.layout
    base_layout = layout.base
    base_u = view(u, layout.base_state)
    base_du = view(du, layout.base_state)
    
    Tch = reshape(view(base_u, base_layout.channel), grid.group_count, bg.nz)
    Ttube = view(base_u, base_layout.tube)
    
    _gas_profile_enthalpy2D_v21!(
        context.work, Tch, Ttube, time, context.p14,
        context.exchange, context.op, grid, context.p.phase2
    )
    
    duch = reshape(view(base_du, base_layout.channel), grid.group_count, bg.nz)
    for group in 1:grid.group_count, j in 1:bg.nz
        full_area = grid.multiplicity[group] * grid.solid_area_per_channel
        core_area = full_area - context.skin_group_area[group]
        capacity = V15._channel_capacity2D(Tch[group, j], core_area, bg.dz_cells[j], context.p15)
        
        delta = context.work.qgas[group, j] - old.qgas[group, j]
        if j == 1
            delta += context.work.front_qgas[group] - old.front[group]
            
            # Phase 2 Physics: Correct front-face radiation area
            # The V15 baseline only radiates from the solid web area. We must radiate
            # from the channel flow area as well (acting as blackbody cavities).
            ambient = Float64(context.op.ambient_temperature(time))
            extra_area = grid.multiplicity[group] * (grid.frontal_area_per_channel - grid.solid_area_per_channel)
            extra_rad = context.p14.losses.front_loss_scale * context.p14.losses.front_emissivity * V11.SIGMA * extra_area * (V11.clamp_temperature(Tch[group, 1])^4 - ambient^4)
            delta += extra_rad
        end
        
        # Phase 2 Physics: Spillage Heating on Perimeter (Orbit 15)
        if group == grid.group_count && context.p.phase2.spillage_power_W > 0.0
            # Apply heating exponentially in first few mm
            L_spill = context.p.phase2.spillage_axial_spread_m
            z = bg.z_centers[j]
            dz = bg.dz_cells[j]
            fraction = exp(-z / L_spill) * (dz / L_spill) # Approximate integral
            spill_q = context.p.phase2.spillage_power_W * fraction
            
            # delta is heat removed by gas. So subtracting spill_q from delta INCREASES solid T
            delta -= spill_q
        end
        
        duch[group, j] -= delta / max(capacity, eps())
    end
    
    tube_area = pi * (context.p14.geometry.rear_tube_outer_radius^2 - context.p14.geometry.rear_tube_inner_radius^2)
    dutube = view(base_du, base_layout.tube)
    for j in 1:bg.nz_rear
        capacity = 3900.0 * V11.alumina_heat_capacity(Ttube[j]) * tube_area * bg.dz_rear
        dutube[j] -= (context.work.qgas_rear[j] - old.rear[j]) / max(capacity, eps())
        
        overlap = V12._overlap_length(
            bg.z_rear_faces[j], bg.z_rear_faces[j + 1],
            context.tube_flange.contact_start_m,
            context.p14.geometry.rear_tube_length,
        )
        if overlap > 0.0 && context.tube_flange.contact_h_W_m2_K > 0.0
            area_contact = 2.0 * pi * context.p14.geometry.rear_tube_outer_radius * overlap
            Q = context.tube_flange.contact_h_W_m2_K * area_contact * (Ttube[j] - V11.WATER_FLANGE_TEMP)
            dutube[j] -= Q / max(capacity, eps())
        end
    end
    return nothing
end

function simulate2D_v21(
    p::ModelParameters2D_v21,
    op::V20.OperatingCondition2D,
    save_times::AbstractVector;
    initial_temperature=Float64(op.ambient_temperature(save_times[1])),
    solver=FBDF(autodiff=AutoFiniteDiff()),
    reltol=1e-6, abstol=1e-7, dtmax=30.0,
)
    apply_v21_property_fixes!()
    
    context = V20._model_context2D(p.base, op)
    context_v21 = merge(context, (p=p,))
    
    grid = context.grid
    layout = context.layout
    p17 = p.base.base.base.base
    p15 = context.p15
    p14 = context.p14
    
    initial = V15._initial_state2D(initial_temperature, p15, grid)
    problem = ODEProblem(
        _v21_ode2D!, initial,
        (Float64(first(save_times)), Float64(last(save_times))),
        context_v21,
    )
    solution = solve(
        problem, solver; saveat=save_times,
        reltol=reltol, abstol=abstol, dtmax=dtmax,
        isoutofdomain=(u, _p, t) -> any(x -> isnan(x) || x < 100.0 || x > 3500.0, u)
    )
    
    base14_result = V15._extract_base_result2D(solution, p15, op, grid, layout)
    skin = hcat((collect(view(state, layout.skin)) for state in solution.u)...)
    base15_result = V15.SimulationResult2D(base14_result, skin, p15, solution)
    base17_result = V17.SimulationResult2D(base15_result, p17, solution)
    
    # Use V20 diagnostics extraction logic
    diagnostics = V20._extract_diagnostics2D(solution, p.base, op, grid, layout, p14)
    
    return V20.SimulationResult2D(base17_result, p.base, diagnostics, solution)
end

end # module
