# ============================================================================
# 2D_v13.jl - Conservative two-solid-population extension of v12
# ============================================================================

module Receiver2D_v13

include("2D_v12.jl")
const V12 = Receiver2D_v12
const V11 = V12.Receiver2D_v11

using DifferentialEquations
using OrdinaryDiffEq
using Statistics

export PopulationParameters2D, ModelParameters2D, SimulationResult2D
export default_parameters2D, simulate2D, sensor_predictions2D
export population_inventory2D
export energy_rate_ledger2D
export Geometry2D, SolidProperties2D, HeatTransferParameters2D
export LossParameters2D, OpticalParameters2D, HydraulicParameters2D
export OperatingCondition2D, AssemblyParameters2D, ObservationParameters2D
export build_initial_state_2D, linear_history, get_t90_2D

const Geometry2D = V12.Geometry2D
const SolidProperties2D = V12.SolidProperties2D
const HeatTransferParameters2D = V12.HeatTransferParameters2D
const LossParameters2D = V12.LossParameters2D
const OpticalParameters2D = V12.OpticalParameters2D
const HydraulicParameters2D = V12.HydraulicParameters2D
const OperatingCondition2D = V12.OperatingCondition2D
const AssemblyParameters2D = V12.AssemblyParameters2D
const ObservationParameters2D = V12.ObservationParameters2D
const build_initial_state_2D = V12.build_initial_state_2D
const linear_history = V12.linear_history
const get_t90_2D = V12.get_t90_2D

Base.@kwdef struct PopulationParameters2D
    side_solid_fraction::Float64 = 0.30
    side_gas_heat_fraction::Float64 = 0.05
    active_side_conductance_per_length::Float64 = 1.0 # W/K/m
    side_felt_contact_h::Float64 = 5.0                # W/m2/K
end

Base.@kwdef struct ModelParameters2D
    base::V12.ModelParameters2D = V12.default_parameters2D()
    population::PopulationParameters2D = PopulationParameters2D()
end

function Base.getproperty(p::ModelParameters2D, name::Symbol)
    name === :base && return getfield(p, :base)
    name === :population && return getfield(p, :population)
    return getproperty(getfield(p, :base), name)
end

default_parameters2D() = ModelParameters2D()

struct SimulationResult2D
    base_result::V12.SimulationResult2D
    side_temperature::Matrix{Float64} # Nz x Nt
    parameters::ModelParameters2D
end

function Base.getproperty(result::SimulationResult2D, name::Symbol)
    if name === :base_result || name === :side_temperature ||
       name === :parameters
        return getfield(result, name)
    end
    return getproperty(getfield(result, :base_result), name)
end

function population_inventory2D(p::ModelParameters2D=default_parameters2D();
                                temperature=800.0)
    inventory = V12.hardware_inventory2D(p.base; temperature=temperature)
    fraction = clamp(p.population.side_solid_fraction, 0.01, 0.80)
    return (
        total_receiver_capacity_J_K=inventory.receiver_capacity_J_K,
        active_receiver_capacity_J_K=
            (1.0 - fraction) * inventory.receiver_capacity_J_K,
        side_receiver_capacity_J_K=
            fraction * inventory.receiver_capacity_J_K,
        side_solid_fraction=fraction,
    )
end

function _population_ode2D!(du, u, context, time)
    p = context.p
    p12 = p.base
    grid = context.grid
    n12 = context.n12
    V12._assembly_ode2D!(
        view(du, 1:n12), view(u, 1:n12),
        context.v12_context, time,
    )
    T_s = reshape(
        view(u, 1:(grid.nr_total * grid.nz)),
        grid.nr_total, grid.nz,
    )
    du_s = reshape(
        view(du, 1:(grid.nr_total * grid.nz)),
        grid.nr_total, grid.nz,
    )
    Tside = view(u, (n12 + 1):(n12 + grid.nz))
    duside = view(du, (n12 + 1):(n12 + grid.nz))
    fill!(duside, 0.0)

    pop = p.population
    fraction = clamp(pop.side_solid_fraction, 0.01, 0.80)
    gas_fraction = clamp(pop.side_gas_heat_fraction, 0.0, fraction)
    active_fraction = 1.0 - fraction
    work = context.v12_context.base_context.work
    base = p12.base

    total_solid_area = sum(grid.area_solid[1:grid.nr_rec])
    capacities = zeros(grid.nz)
    active_mean = zeros(grid.nz)
    for j in 1:grid.nz
        capacities[j] = sum(
            p12.solid.density * V11.sic_heat_capacity(T_s[i, j]) *
            grid.area_solid[i] * grid.dz_cells[j]
            for i in 1:grid.nr_rec
        )
        active_mean[j] = sum(
            grid.area_solid[i] * T_s[i, j]
            for i in 1:grid.nr_rec
        ) / total_solid_area
    end

    irradiance = max(0.0, Float64(
        context.v12_context.base_context.op.irradiance(time),
    ))
    delivered = p12.optics.absorbed_fraction * irradiance *
                grid.total_frt
    core_power = (1.0 - p12.optics.spillage_fraction) * delivered
    weights = context.v12_context.base_context.weights

    side_power = zeros(grid.nz)
    for j in 1:grid.nz
        solar_j = core_power * sum(weights[i, j] for i in 1:grid.nr_rec)
        gas_cell = collect(view(work.qgas, :, j))
        j == 1 && (gas_cell .+= work.front_qgas)
        gas_j = sum(gas_cell)

        # Reallocate gas extraction between the two physical fractions.
        # The original v12 derivative corresponds to proportional allocation.
        for i in 1:grid.nr_rec
            cap_active = active_fraction *
                p12.solid.density * V11.sic_heat_capacity(T_s[i, j]) *
                grid.area_solid[i] * grid.dz_cells[j]
            du_s[i, j] += (gas_fraction - fraction) * gas_cell[i] /
                          max(cap_active, eps(Float64))
        end
        side_power[j] += fraction * solar_j - gas_fraction * gas_j

        # The front face exists in both solid populations.  V12 applies the
        # receiver-web radiation to the active field; its conceptual energy
        # therefore retains (1-f) of that loss.  Add the complementary side
        # loss at the side-population temperature.
        if j == 1
            ambient = Float64(
                context.v12_context.base_context.op.ambient_temperature(time),
            )
            Qsidefront = fraction * base.losses.front_loss_scale *
                base.losses.front_emissivity * V11.SIGMA *
                total_solid_area *
                (V11.clamp_temperature(Tside[j])^4 - ambient^4)
            side_power[j] -= Qsidefront
        end

        Gaw = pop.active_side_conductance_per_length *
              grid.dz_cells[j]
        Qaw = Gaw * (active_mean[j] - Tside[j])
        side_power[j] += Qaw
        for i in 1:grid.nr_rec
            share = grid.area_solid[i] / total_solid_area
            cap_active = active_fraction *
                p12.solid.density * V11.sic_heat_capacity(T_s[i, j]) *
                grid.area_solid[i] * grid.dz_cells[j]
            du_s[i, j] -= share * Qaw /
                          max(cap_active, eps(Float64))
        end

        # Correct receiver/felt exchange: the active field owns only its
        # mass fraction; the side population has its own felt contact.
        irec = grid.nr_rec
        ifelt = irec + 1
        dr_rec = grid.r_faces[irec + 1] - grid.r_faces[irec]
        dr_felt = grid.r_faces[ifelt + 1] - grid.r_faces[ifelt]
        krec = p12.solid.radial_conductivity_scale *
               V11.sic_conductivity(T_s[irec, j])
        kfelt = V12.felt_conductivity_v12(
            T_s[ifelt, j],
            p12.assembly.felt_conductivity_scale,
        )
        resistance = 0.5 * dr_rec / krec + 0.5 * dr_felt / kfelt +
                     p12.solid.receiver_felt_contact_resistance
        area = 2.0 * pi * grid.r_faces[irec + 1] * grid.dz_cells[j]
        Qoriginal = area * (T_s[irec, j] - T_s[ifelt, j]) /
                    resistance
        Qsidefelt = pop.side_felt_contact_h *
                    (4.0 * base.geometry.receiver_width *
                     grid.dz_cells[j]) *
                    (Tside[j] - T_s[ifelt, j])
        side_power[j] -= Qsidefelt
        delta_felt = -fraction * Qoriginal + Qsidefelt
        V12._add_power_v12!(
            du_s, ifelt, j, delta_felt, T_s, p12, grid,
        )

        # The loose adaptor support acts on the side/support population, not
        # on the gas-exposed active-channel population.  Remove the inherited
        # V12 outer-cell/adaptor path and replace it conservatively below.
        overlap_start = base.geometry.length -
                        p12.assembly.adaptor_receiver_overlap
        overlap = V12._overlap_length(
            grid.z_faces[j], grid.z_faces[j + 1],
            overlap_start, base.geometry.length,
        )
        if overlap > 0.0
            contact_area = 4.0 * base.geometry.receiver_width * overlap
            Qoriginal_adaptor =
                p12.assembly.receiver_adaptor_contact_h * contact_area *
                (T_s[irec, j] - u[context.v12_context.nbase + 1])
            Qside_adaptor =
                p12.assembly.receiver_adaptor_contact_h * contact_area *
                (Tside[j] - u[context.v12_context.nbase + 1])
            cap_active = active_fraction * p12.solid.density *
                V11.sic_heat_capacity(T_s[irec, j]) *
                grid.area_solid[irec] * grid.dz_cells[j]
            du_s[irec, j] += (1.0 - fraction) * Qoriginal_adaptor /
                             max(cap_active, eps(Float64))
            side_power[j] -= Qside_adaptor

            Tad = u[context.v12_context.nbase + 1]
            Cad = p12.assembly.adaptor_density *
                  V12.adaptor_volume2D(p12) *
                  V11.alumina_heat_capacity(Tad)
            du[context.v12_context.nbase + 1] +=
                (-Qoriginal_adaptor + Qside_adaptor) /
                max(Cad, eps(Float64))
        end
    end

    # Axial conduction carried by the side fraction.
    for j in 1:(grid.nz - 1)
        dz = 0.5 * (grid.dz_cells[j] + grid.dz_cells[j + 1])
        k0 = p12.solid.axial_conductivity_scale *
             V11.sic_conductivity(Tside[j])
        k1 = p12.solid.axial_conductivity_scale *
             V11.sic_conductivity(Tside[j + 1])
        kface = 2.0 * k0 * k1 / (k0 + k1)
        Q = fraction * kface * total_solid_area *
            (Tside[j] - Tside[j + 1]) / dz
        side_power[j] -= Q
        side_power[j + 1] += Q
    end
    for j in 1:grid.nz
        duside[j] = side_power[j] /
            max(fraction * capacities[j], eps(Float64))
    end
    return nothing
end

function _contexts_2D(p::ModelParameters2D,
                      op::OperatingCondition2D)
    p12 = p.base
    base = p12.base
    grid = V11.build_grid2D(base)
    weights = V11.solar_weights2D(grid, base)
    work = V11.Workspace2D(grid, base)
    nbase = grid.nr_total * grid.nz + grid.nz_rear
    base_context = (
        grid=grid, nz_rear=grid.nz_rear, p=base, op=op,
        weights=weights, work=work,
    )
    v12_context = (
        grid=grid, p=p12, nbase=nbase, base_context=base_context,
    )
    n12 = nbase + 2
    context = (
        p=p, grid=grid, n12=n12, v12_context=v12_context,
    )
    return context
end

function simulate2D(
    p::ModelParameters2D,
    op::OperatingCondition2D,
    save_times::AbstractVector;
    initial_temperature=Float64(op.ambient_temperature(save_times[1])),
    solver=FBDF(autodiff=AutoFiniteDiff()),
    reltol=1e-6, abstol=1e-7, dtmax=30.0,
)
    p12 = p.base
    context = _contexts_2D(p, op)
    grid = context.grid
    work = context.v12_context.base_context.work
    nbase = context.v12_context.nbase
    n12 = context.n12
    initial12 = V12._initial_augmented_state(
        initial_temperature, p12, grid,
    )
    T_s0 = reshape(
        view(initial12, 1:(grid.nr_total * grid.nz)),
        grid.nr_total, grid.nz,
    )
    side0 = [T_s0[grid.nr_rec, j] for j in 1:grid.nz]
    initial = vcat(initial12, side0)
    problem = ODEProblem(
        _population_ode2D!, initial,
        (Float64(first(save_times)), Float64(last(save_times))),
        context,
    )
    solution = solve(
        problem, solver; saveat=save_times, reltol=reltol,
        abstol=abstol, dtmax=dtmax,
        isoutofdomain=(u, p, t) -> any(
            x -> isnan(x) || x < 100.0 || x > 3500.0, u,
        ),
    )
    base_result = V12._base_result_from_solution(
        solution, p12, op, grid, work,
    )
    adaptor = [state[nbase + 1] for state in solution.u]
    housing = [state[nbase + 2] for state in solution.u]
    result12 = V12.SimulationResult2D(
        base_result, adaptor, housing, p12,
    )
    side = hcat((
        collect(view(state, (n12 + 1):(n12 + grid.nz)))
        for state in solution.u
    )...)
    return SimulationResult2D(result12, side, p)
end

function energy_rate_ledger2D(
    u::AbstractVector,
    p::ModelParameters2D,
    op::OperatingCondition2D,
    time::Real=0.0,
)
    context = _contexts_2D(p, op)
    grid = context.grid
    nbase = context.v12_context.nbase
    n12 = context.n12
    length(u) == n12 + grid.nz ||
        error("state vector length mismatch in v13 energy ledger")

    du = zeros(Float64, length(u))
    _population_ode2D!(du, u, context, Float64(time))
    p12 = p.base
    fraction = clamp(p.population.side_solid_fraction, 0.01, 0.80)
    active_fraction = 1.0 - fraction
    solid_count = grid.nr_total * grid.nz
    T_s = reshape(view(u, 1:solid_count), grid.nr_total, grid.nz)
    du_s = reshape(view(du, 1:solid_count), grid.nr_total, grid.nz)

    denergy_dt = 0.0
    for i in 1:grid.nr_total, j in 1:grid.nz
        capacity = V12._cell_capacity_v12(
            i, j, T_s[i, j], p12, grid,
        )
        i <= grid.nr_rec && (capacity *= active_fraction)
        denergy_dt += capacity * du_s[i, j]
    end

    tube_area = pi * (
        p12.geometry.rear_tube_outer_radius^2 -
        p12.geometry.rear_tube_inner_radius^2
    )
    Ttube = view(u, (solid_count + 1):nbase)
    dutube = view(du, (solid_count + 1):nbase)
    for j in 1:grid.nz_rear
        capacity = 3900.0 * V11.alumina_heat_capacity(Ttube[j]) *
                   tube_area * grid.dz_rear
        denergy_dt += capacity * dutube[j]
    end

    Tad = u[nbase + 1]
    Cad = p12.assembly.adaptor_density * V12.adaptor_volume2D(p12) *
          V11.alumina_heat_capacity(Tad)
    denergy_dt += Cad * du[nbase + 1]

    a = p12.assembly
    extension_length = max(
        0.0, a.enclosure_internal_length - p12.geometry.length,
    )
    annulus_area = pi * (
        a.enclosure_outer_radius^2 - a.enclosure_inner_radius^2
    )
    backplate_area = pi * (
        a.enclosure_outer_radius^2 -
        p12.geometry.rear_tube_inner_radius^2
    )
    Vhousing = annulus_area * extension_length +
               backplate_area * a.backplate_thickness
    Chousing = a.aluminum_density * a.aluminum_heat_capacity * Vhousing
    denergy_dt += Chousing * du[nbase + 2]

    Tside = view(u, (n12 + 1):(n12 + grid.nz))
    duside = view(du, (n12 + 1):(n12 + grid.nz))
    for j in 1:grid.nz
        Cside = fraction * sum(
            p12.solid.density * V11.sic_heat_capacity(T_s[i, j]) *
            grid.area_solid[i] * grid.dz_cells[j]
            for i in 1:grid.nr_rec
        )
        denergy_dt += Cside * duside[j]
    end

    base_ledger = V12.energy_rate_ledger2D(
        view(u, 1:n12), p12, op, time,
    )
    ambient = Float64(op.ambient_temperature(time))
    total_solid_area = sum(grid.area_solid[1:grid.nr_rec])
    base_receiver_front = p12.losses.front_loss_scale *
        p12.losses.front_emissivity * V11.SIGMA *
        sum(
            grid.area_solid[i] *
            (V11.clamp_temperature(T_s[i, 1])^4 - ambient^4)
            for i in 1:grid.nr_rec
        )
    side_receiver_front = fraction * p12.losses.front_loss_scale *
        p12.losses.front_emissivity * V11.SIGMA * total_solid_area *
        (V11.clamp_temperature(Tside[1])^4 - ambient^4)
    expected = base_ledger.expected_denergy_dt +
               fraction * base_receiver_front -
               side_receiver_front
    return (
        denergy_dt=denergy_dt,
        expected_denergy_dt=expected,
        residual=denergy_dt - expected,
        base_external_rate=base_ledger.expected_denergy_dt,
        active_receiver_front_loss=
            active_fraction * base_receiver_front,
        side_receiver_front_loss=side_receiver_front,
    )
end

function _sample_side(result::SimulationResult2D, z)
    base = result.base_result.base_result
    j0, j1, w = V11._axis_interp_indices(base.z_solid, z)
    return (1.0 - w) .* result.side_temperature[j0, :] .+
           w .* result.side_temperature[j1, :]
end

function sensor_predictions2D(result::SimulationResult2D)
    p12raw = V12.sensor_predictions2D(result.base_result)
    obs = result.parameters.observation
    time = result.time
    return (
        T8=V12._filter_observation(
            time, _sample_side(result, 5e-3),
            obs.side_time_constant_s,
        ),
        T12=V12._filter_observation(
            time, _sample_side(result, 58e-3),
            obs.side_time_constant_s,
        ),
        T11=V12._filter_observation(
            time, _sample_side(result, 107e-3),
            obs.side_time_constant_s,
        ),
        T9=p12raw.T9, T10=p12raw.T10, T3=p12raw.T3, T2=p12raw.T2,
    )
end

end
