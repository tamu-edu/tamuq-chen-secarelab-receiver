# ============================================================================
# 2D_v12.jl - Assembly-capacitance and thermocouple-observation extension
# ============================================================================
#
# V12 deliberately builds on the conservative v11 gas/receiver solver.  It
# adds the hardware and observation physics identified after the v11 rejection:
#
#   * the 77.6 mm OD, 57 mm long dense-alumina adaptor (29 mm loose receiver
#     overlap and 28 mm tube overlap);
#   * the missing 28 mm aluminum sleeve extension and 18 mm rear backplate;
#   * finite receiver/adaptor, adaptor/tube and adaptor/felt coupling;
#   * datasheet-shaped, bounded effective felt k(T) and Cp scaling;
#   * gas/solid-aware interior-thermocouple observations and shared sensor lag.
#
# The full hardware heat capacities are retained.  Their participation is
# determined by finite conductances; they are not added to the SiC Cp.
# ============================================================================

module Receiver2D_v12

include("2D_v11.jl")
const V11 = Receiver2D_v11

using DifferentialEquations
using OrdinaryDiffEq
using SciMLBase
using Statistics

export Geometry2D, SolidProperties2D, HeatTransferParameters2D
export LossParameters2D, OpticalParameters2D, HydraulicParameters2D
export OperatingCondition2D, AssemblyParameters2D, ObservationParameters2D
export ModelParameters2D, SimulationResult2D
export default_parameters2D, simulate2D, sensor_predictions2D
export build_initial_state_2D, linear_history, get_t90_2D
export geometry_invariants2D, hardware_inventory2D
export felt_conductivity_v12, adaptor_volume2D
export energy_rate_ledger2D

const Geometry2D = V11.Geometry2D
const SolidProperties2D = V11.SolidProperties2D
const HeatTransferParameters2D = V11.HeatTransferParameters2D
const LossParameters2D = V11.LossParameters2D
const OpticalParameters2D = V11.OpticalParameters2D
const HydraulicParameters2D = V11.HydraulicParameters2D
const OperatingCondition2D = V11.OperatingCondition2D
const linear_history = V11.linear_history
const get_t90_2D = V11.get_t90_2D
const build_initial_state_2D = V11.build_initial_state_2D
const geometry_invariants2D = V11.geometry_invariants2D

Base.@kwdef struct AssemblyParameters2D
    enclosure_inner_radius::Float64 = 57.0e-3
    enclosure_outer_radius::Float64 = 75.0e-3
    enclosure_internal_length::Float64 = 165.0e-3
    backplate_thickness::Float64 = 18.0e-3

    adaptor_outer_radius::Float64 = 38.8e-3
    adaptor_length::Float64 = 57.0e-3
    adaptor_receiver_overlap::Float64 = 29.0e-3
    adaptor_tube_overlap::Float64 = 28.0e-3
    adaptor_density::Float64 = 3900.0

    # The receiver is supported but the contact is explicitly not firm.
    receiver_adaptor_contact_h::Float64 = 25.0       # W/m2/K
    adaptor_tube_contact_h::Float64 = 80.0           # W/m2/K
    adaptor_felt_conductance_scale::Float64 = 1.0

    # The felt grade is unknown.  The nominal curve follows the supplied
    # alumina-felt sheet only as a shape prior; both k and Cp may move modestly.
    felt_conductivity_scale::Float64 = 1.0
    felt_heat_capacity_scale::Float64 = 1.0

    aluminum_density::Float64 = 2700.0
    aluminum_heat_capacity::Float64 = 900.0
    aluminum_conductivity::Float64 = 205.0
    housing_extension_loss_scale::Float64 = 1.0
end

Base.@kwdef struct ObservationParameters2D
    side_time_constant_s::Float64 = 45.0
    interior_time_constant_s::Float64 = 25.0
    outlet_time_constant_s::Float64 = 20.0
    # Interior probes are flow exposed.  This is the wall fraction in their
    # equilibrium target; the balance is the local channel-gas temperature.
    interior_wall_fraction::Float64 = 0.55
end

Base.@kwdef struct ModelParameters2D
    base::V11.ModelParameters2D = V11.default_parameters2D()
    assembly::AssemblyParameters2D = AssemblyParameters2D()
    observation::ObservationParameters2D = ObservationParameters2D()
end

function ModelParameters2D(
    geometry::Geometry2D,
    solid::SolidProperties2D,
    heat_transfer::HeatTransferParameters2D,
    losses::LossParameters2D,
    optics::OpticalParameters2D,
    hydraulics::HydraulicParameters2D;
    assembly=AssemblyParameters2D(),
    observation=ObservationParameters2D(),
)
    return ModelParameters2D(
        V11.ModelParameters2D(
            geometry, solid, heat_transfer, losses, optics, hydraulics,
        ),
        assembly,
        observation,
    )
end

function default_parameters2D()
    b0 = V11.default_parameters2D()
    l0 = b0.losses
    losses = LossParameters2D(
        front_emissivity=l0.front_emissivity,
        casing_emissivity=l0.casing_emissivity,
        front_loss_scale=l0.front_loss_scale,
        casing_loss_scale=l0.casing_loss_scale,
        # V12 resolves the physical adaptor explicitly.  The v11 direct
        # receiver-to-tube surrogate must be disabled to avoid a bypass.
        rear_adaptor_conductance=0.0,
        flange_conductance_scale=l0.flange_conductance_scale,
    )
    base = V11.ModelParameters2D(
        b0.geometry, b0.solid, b0.heat_transfer, losses,
        b0.optics, b0.hydraulics,
    )
    return ModelParameters2D(base=base)
end

function Base.getproperty(p::ModelParameters2D, name::Symbol)
    if name === :base || name === :assembly || name === :observation
        return getfield(p, name)
    end
    return getproperty(getfield(p, :base), name)
end

"""
Nominal felt conductivity curve used in v12.

The supplied RS-3000 sheet reports 0.08 W/m/K at 500 C, 0.11 at 800 C,
0.17 at 1100 C, 0.26 at 1400 C and 0.32 at 1600 C.  The installed grade is
unknown, so this curve is a temperature-shape prior and a bounded scale is
identified separately.
"""
function felt_conductivity_v12(T::Real, scale::Real=1.0)
    Tc = clamp(Float64(T) - 273.15, 20.0, 1600.0)
    knots_T = [20.0, 500.0, 800.0, 1100.0, 1400.0, 1600.0]
    knots_k = [0.050, 0.080, 0.110, 0.170, 0.260, 0.320]
    index = searchsortedlast(knots_T, Tc)
    nominal = if index <= 0
        knots_k[1]
    elseif index >= length(knots_T)
        knots_k[end]
    else
        fraction = (Tc - knots_T[index]) /
                   (knots_T[index + 1] - knots_T[index])
        knots_k[index] + fraction *
            (knots_k[index + 1] - knots_k[index])
    end
    return max(0.025, Float64(scale) * nominal)
end

function adaptor_volume2D(p::ModelParameters2D=default_parameters2D())
    a = p.assembly
    g = p.geometry
    receiver_area = g.receiver_width^2
    tube_area = pi * g.rear_tube_inner_radius^2
    overlap_volume = (
        pi * a.adaptor_outer_radius^2 - receiver_area
    ) * a.adaptor_receiver_overlap
    tube_volume = (
        pi * a.adaptor_outer_radius^2 - tube_area
    ) * a.adaptor_tube_overlap
    return overlap_volume + tube_volume
end

function hardware_inventory2D(p::ModelParameters2D=default_parameters2D();
                              temperature=800.0)
    base = p.base
    a = p.assembly
    grid = V11.build_grid2D(base)

    receiver_mass = V11.geometry_invariants2D(base).receiver_mass
    receiver_capacity =
        receiver_mass * V11.sic_heat_capacity(temperature)

    # Actual felt-filled union: bare receiver for 108 mm, adaptor OD for
    # 57 mm.  The tube is inside the adaptor or outside the cavity.
    bare_receiver_length = max(
        0.0,
        base.geometry.length - a.adaptor_receiver_overlap,
    )
    felt_volume = (
        pi * a.enclosure_inner_radius^2 *
            a.enclosure_internal_length -
        base.geometry.receiver_width^2 * bare_receiver_length -
        pi * a.adaptor_outer_radius^2 * a.adaptor_length
    )
    felt_mass = base.solid.felt_density * felt_volume
    felt_capacity = felt_mass * base.solid.felt_heat_capacity *
                    a.felt_heat_capacity_scale

    adaptor_volume = adaptor_volume2D(p)
    adaptor_mass = a.adaptor_density * adaptor_volume
    adaptor_capacity =
        adaptor_mass * V11.alumina_heat_capacity(temperature)

    sleeve_volume = pi * (
        a.enclosure_outer_radius^2 -
        a.enclosure_inner_radius^2
    ) * a.enclosure_internal_length
    backplate_volume = pi * (
        a.enclosure_outer_radius^2 -
        base.geometry.rear_tube_inner_radius^2
    ) * a.backplate_thickness
    sleeve_mass = a.aluminum_density * sleeve_volume
    backplate_mass = a.aluminum_density * backplate_volume

    tube_area = pi * (
        base.geometry.rear_tube_outer_radius^2 -
        base.geometry.rear_tube_inner_radius^2
    )
    tube_mass = 3900.0 * tube_area * base.geometry.rear_tube_length
    tube_capacity =
        tube_mass * V11.alumina_heat_capacity(temperature)

    return (
        receiver_mass_kg=receiver_mass,
        receiver_capacity_J_K=receiver_capacity,
        felt_volume_m3=felt_volume,
        felt_mass_kg=felt_mass,
        felt_capacity_J_K=felt_capacity,
        adaptor_volume_m3=adaptor_volume,
        adaptor_mass_kg=adaptor_mass,
        adaptor_capacity_J_K=adaptor_capacity,
        aluminum_sleeve_mass_kg=sleeve_mass,
        aluminum_sleeve_capacity_J_K=
            sleeve_mass * a.aluminum_heat_capacity,
        aluminum_backplate_mass_kg=backplate_mass,
        aluminum_backplate_capacity_J_K=
            backplate_mass * a.aluminum_heat_capacity,
        rear_tube_mass_kg=tube_mass,
        rear_tube_capacity_J_K=tube_capacity,
    )
end

struct SimulationResult2D
    base_result::V11.SimulationResult2D
    adaptor_temperature::Vector{Float64}
    housing_extension_temperature::Vector{Float64}
    parameters::ModelParameters2D
end

function Base.getproperty(result::SimulationResult2D, name::Symbol)
    if name === :base_result || name === :adaptor_temperature ||
       name === :housing_extension_temperature || name === :parameters
        return getfield(result, name)
    end
    return getproperty(getfield(result, :base_result), name)
end

function _cell_capacity_v12(i, j, T, p, grid)
    base = p.base
    prop = V11.cell_thermal_properties2D(i, T, base, grid)
    volume = (i <= grid.nr_rec ? grid.area_solid[i] :
              grid.area_frt[i]) * grid.dz_cells[j]
    scale = if i > grid.nr_rec &&
               i <= grid.nr_rec + grid.nr_felt
        p.assembly.felt_heat_capacity_scale
    else
        1.0
    end
    return prop.rho * prop.cp * volume * scale
end

function _desired_k_v12(i, T, p, grid)
    if i <= grid.nr_rec
        return p.solid.radial_conductivity_scale *
               V11.sic_conductivity(T)
    elseif i <= grid.nr_rec + grid.nr_felt
        return felt_conductivity_v12(
            T, p.assembly.felt_conductivity_scale,
        )
    end
    return p.solid.casing_conductivity
end

function _base_k_v12(i, T, p, grid)
    return V11.cell_thermal_properties2D(
        i, T, p.base, grid,
    ).k
end

function _add_power_v12!(du_s, i, j, power, T_s, p, grid)
    du_s[i, j] += power /
        max(_cell_capacity_v12(i, j, T_s[i, j], p, grid), eps(Float64))
end

function _correct_felt_transport_v12!(du_s, T_s, p, grid)
    nr_rec = grid.nr_rec
    first_felt = nr_rec + 1
    last_felt = nr_rec + grid.nr_felt

    # Base v11 has already divided the felt equations by the nominal Cp.
    # Rescale every felt energy balance to the v12 Cp before adding flux
    # corrections.
    cp_scale = max(p.assembly.felt_heat_capacity_scale, eps(Float64))
    for i in first_felt:last_felt, j in 1:grid.nz
        du_s[i, j] /= cp_scale
    end

    for i in nr_rec:(last_felt), j in 1:grid.nz
        # Only interfaces touching felt need correction.
        (i >= first_felt || i + 1 >= first_felt) || continue
        (i <= last_felt || i + 1 <= last_felt) || continue
        dr_i = grid.r_faces[i + 1] - grid.r_faces[i]
        dr_ip1 = grid.r_faces[i + 2] - grid.r_faces[i + 1]

        kb_i = _base_k_v12(i, T_s[i, j], p, grid)
        kb_ip1 = _base_k_v12(i + 1, T_s[i + 1, j], p, grid)
        kd_i = _desired_k_v12(i, T_s[i, j], p, grid)
        kd_ip1 = _desired_k_v12(i + 1, T_s[i + 1, j], p, grid)

        rb = 0.5 * dr_i / kb_i + 0.5 * dr_ip1 / kb_ip1
        rd = 0.5 * dr_i / kd_i + 0.5 * dr_ip1 / kd_ip1
        if i == nr_rec
            contact = max(
                0.0,
                p.solid.receiver_felt_contact_resistance,
            )
            rb += contact
            rd += contact
        end
        area = 2.0 * pi * grid.r_faces[i + 1] *
               grid.dz_cells[j]
        delta_power = area * (
            T_s[i, j] - T_s[i + 1, j]
        ) * (1.0 / rd - 1.0 / rb)
        _add_power_v12!(
            du_s, i, j, -delta_power, T_s, p, grid,
        )
        _add_power_v12!(
            du_s, i + 1, j, delta_power, T_s, p, grid,
        )
    end

    for i in first_felt:last_felt, j in 1:(grid.nz - 1)
        dz = 0.5 * (grid.dz_cells[j] + grid.dz_cells[j + 1])
        kb0 = _base_k_v12(i, T_s[i, j], p, grid)
        kb1 = _base_k_v12(i, T_s[i, j + 1], p, grid)
        kd0 = _desired_k_v12(i, T_s[i, j], p, grid)
        kd1 = _desired_k_v12(i, T_s[i, j + 1], p, grid)
        kb = 2.0 * kb0 * kb1 / (kb0 + kb1)
        kd = 2.0 * kd0 * kd1 / (kd0 + kd1)
        delta_power = (kd - kb) * grid.area_frt[i] *
                      (T_s[i, j] - T_s[i, j + 1]) / dz
        _add_power_v12!(
            du_s, i, j, -delta_power, T_s, p, grid,
        )
        _add_power_v12!(
            du_s, i, j + 1, delta_power, T_s, p, grid,
        )
    end
    return nothing
end

function _overlap_length(z0, z1, a0, a1)
    return max(0.0, min(z1, a1) - max(z0, a0))
end

function _assembly_ode2D!(du, u, context, time)
    p = context.p
    base = p.base
    grid = context.grid
    nbase = context.nbase

    V11.receiver_ode2D!(
        view(du, 1:nbase),
        view(u, 1:nbase),
        context.base_context,
        time,
    )

    nr_total = grid.nr_total
    nz = grid.nz
    nz_rear = grid.nz_rear
    solid_count = nr_total * nz
    T_s = reshape(view(u, 1:solid_count), nr_total, nz)
    T_tube = view(u, (solid_count + 1):(solid_count + nz_rear))
    du_s = reshape(view(du, 1:solid_count), nr_total, nz)
    du_tube = view(du, (solid_count + 1):(solid_count + nz_rear))

    _correct_felt_transport_v12!(du_s, T_s, p, grid)

    Tad = u[nbase + 1]
    Thousing = u[nbase + 2]
    du[nbase + 1] = 0.0
    du[nbase + 2] = 0.0
    a = p.assembly

    adaptor_power = 0.0
    overlap_start = base.geometry.length -
                    a.adaptor_receiver_overlap
    for j in 1:nz
        overlap = _overlap_length(
            grid.z_faces[j], grid.z_faces[j + 1],
            overlap_start, base.geometry.length,
        )
        overlap <= 0.0 && continue
        area = 4.0 * base.geometry.receiver_width * overlap
        Q = a.receiver_adaptor_contact_h * area *
            (T_s[grid.nr_rec, j] - Tad)
        _add_power_v12!(
            du_s, grid.nr_rec, j, -Q, T_s, p, grid,
        )
        adaptor_power += Q
    end

    # The first 28 mm of the tube is inside the alumina adaptor.
    tube_wall_area = pi * (
        base.geometry.rear_tube_outer_radius^2 -
        base.geometry.rear_tube_inner_radius^2
    )
    for j in 1:nz_rear
        overlap = _overlap_length(
            grid.z_rear_faces[j], grid.z_rear_faces[j + 1],
            0.0, a.adaptor_tube_overlap,
        )
        overlap <= 0.0 && continue
        area = 2.0 * pi * base.geometry.rear_tube_outer_radius *
               overlap
        Q = a.adaptor_tube_contact_h * area * (T_tube[j] - Tad)
        tube_capacity = 3900.0 *
            V11.alumina_heat_capacity(T_tube[j]) *
            tube_wall_area * grid.dz_rear
        du_tube[j] -= Q / max(tube_capacity, eps(Float64))
        adaptor_power += Q
    end

    # Felt fills the remaining cavity around the adaptor.  The represented
    # rear felt cells carry this exchange; the 57/29 area factor accounts for
    # the free rear adaptor section without inventing tube/felt contact.
    felt_indices = (grid.nr_rec + 1):(
        grid.nr_rec + grid.nr_felt
    )
    interface_radius = a.adaptor_outer_radius
    kfelt = felt_conductivity_v12(
        Tad, a.felt_conductivity_scale,
    )
    G_adaptor_felt = a.adaptor_felt_conductance_scale *
        2.0 * pi * kfelt * a.adaptor_length /
        log(a.enclosure_inner_radius / interface_radius)
    rear_weights = zeros(nz)
    for j in 1:nz
        rear_weights[j] = _overlap_length(
            grid.z_faces[j], grid.z_faces[j + 1],
            overlap_start, base.geometry.length,
        )
    end
    weight_sum = sum(rear_weights)
    if weight_sum > 0.0
        inner_felt = first(felt_indices)
        Tfelt = sum(
            rear_weights[j] * T_s[inner_felt, j]
            for j in 1:nz
        ) / weight_sum
        Q = G_adaptor_felt * (Tad - Tfelt)
        adaptor_power -= Q
        for j in 1:nz
            rear_weights[j] <= 0.0 && continue
            _add_power_v12!(
                du_s, inner_felt, j,
                Q * rear_weights[j] / weight_sum,
                T_s, p, grid,
            )
        end
    end

    Cad = a.adaptor_density * adaptor_volume2D(p) *
          V11.alumina_heat_capacity(Tad)
    du[nbase + 1] = adaptor_power / max(Cad, eps(Float64))

    # Missing 28 mm sleeve plus 18 mm backplate, represented as one aluminum
    # rear-housing state coupled to the end of the existing casing field.
    extension_length = max(
        0.0,
        a.enclosure_internal_length - base.geometry.length,
    )
    annulus_area = pi * (
        a.enclosure_outer_radius^2 -
        a.enclosure_inner_radius^2
    )
    backplate_area = pi * (
        a.enclosure_outer_radius^2 -
        base.geometry.rear_tube_inner_radius^2
    )
    Vhousing = annulus_area * extension_length +
               backplate_area * a.backplate_thickness
    Chousing = a.aluminum_density * a.aluminum_heat_capacity *
               Vhousing
    case_index = grid.nr_total
    Tcase_end = T_s[case_index, end]
    Gcase = a.aluminum_conductivity * annulus_area /
            max(0.5 * extension_length, 1.0e-3)
    Qcase = Gcase * (Tcase_end - Thousing)
    _add_power_v12!(
        du_s, case_index, nz, -Qcase, T_s, p, grid,
    )

    ambient = Float64(context.base_context.op.ambient_temperature(time))
    area_loss = (
        2.0 * pi * a.enclosure_outer_radius * extension_length +
        pi * a.enclosure_outer_radius^2
    )
    Qloss = a.housing_extension_loss_scale * area_loss * (
        10.0 * (Thousing - ambient) +
        base.losses.casing_emissivity * V11.SIGMA *
        (V11.clamp_temperature(Thousing)^4 - ambient^4)
    )
    du[nbase + 2] = (Qcase - Qloss) /
                    max(Chousing, eps(Float64))
    return nothing
end

function _initial_augmented_state(initial_temperature, p, grid)
    nbase = grid.nr_total * grid.nz + grid.nz_rear
    if initial_temperature isa AbstractVector
        if length(initial_temperature) == nbase + 2
            return copy(Float64.(initial_temperature))
        end
        length(initial_temperature) == nbase ||
            error("initial_temperature vector length mismatch")
        base_initial = copy(Float64.(initial_temperature))
        T_s = reshape(
            view(base_initial, 1:(grid.nr_total * grid.nz)),
            grid.nr_total, grid.nz,
        )
        T_tube = view(
            base_initial,
            (grid.nr_total * grid.nz + 1):nbase,
        )
        Tad = mean((
            T_s[grid.nr_rec, end],
            T_tube[1],
        ))
        Thousing = T_s[grid.nr_total, end]
        return vcat(base_initial, Tad, Thousing)
    end
    return fill(Float64(initial_temperature), nbase + 2)
end

function _base_result_from_solution(sol, p, op, grid, work)
    base = p.base
    nr_rec, nr_total = grid.nr_rec, grid.nr_total
    nz, nz_rear = grid.nz, grid.nz_rear
    nt = length(sol.t)
    solid_count = nr_total * nz

    solid_T = zeros(nr_total, nz, nt)
    gas_T = zeros(nr_rec, nz + 1, nt)
    rear_tube_T = zeros(nz_rear, nt)
    rear_gas_T = zeros(nz_rear + 1, nt)
    htc = zeros(nr_rec, nz, nt)
    front_htc = zeros(nr_rec, nt)
    front_qgas = zeros(nr_rec, nt)
    gas_rho = zeros(nr_rec, nz, nt)
    gas_velocity = zeros(nr_rec, nz, nt)
    gas_reynolds = zeros(nr_rec, nz, nt)
    gas_prandtl = zeros(nr_rec, nz, nt)
    gas_graetz = zeros(nr_rec, nz, nt)
    gas_nusselt = zeros(nr_rec, nz, nt)
    cell_dp = zeros(nr_rec, nz, nt)
    ring_mdot = zeros(nr_rec, nt)
    ring_dp = zeros(nr_rec, nt)
    groove_dp = zeros(nr_rec, nt)
    equal_pressure_error = zeros(nt)
    receiver_dp_mbar = zeros(nt)
    dp1_prediction_mbar = zeros(nt)
    mass_flow_kg_s = zeros(nt)

    for k in 1:nt
        state = sol.u[k]
        solid_T[:, :, k] .= reshape(
            view(state, 1:solid_count), nr_total, nz,
        )
        rear_tube_T[:, k] .= view(
            state, (solid_count + 1):(solid_count + nz_rear),
        )
        V11._gas_profile2D!(
            work, solid_T[1:nr_rec, :, k], rear_tube_T[:, k],
            sol.t[k], base, op, grid,
        )
        gas_T[:, :, k] .= work.gas
        rear_gas_T[:, k] .= work.gas_rear
        htc[:, :, k] .= work.h
        front_htc[:, k] .= work.front_h
        front_qgas[:, k] .= work.front_qgas
        gas_rho[:, :, k] .= work.density
        gas_velocity[:, :, k] .= work.velocity
        gas_reynolds[:, :, k] .= work.reynolds
        gas_prandtl[:, :, k] .= work.prandtl
        gas_graetz[:, :, k] .= work.graetz
        gas_nusselt[:, :, k] .= work.nusselt
        cell_dp[:, :, k] .= work.dp_cell
        ring_mdot[:, k] .= work.ring_mdot
        ring_dp[:, k] .= work.ring_dp
        groove_dp[:, k] .= work.groove_dp
        equal_pressure_error[k] = work.equal_pressure_error
        hyd = V11.hydraulic_summary2D(
            work, base, op, grid, sol.t[k],
        )
        receiver_dp_mbar[k] = hyd.receiver_pressure_drop_mbar
        dp1_prediction_mbar[k] = hyd.dp1_prediction_mbar
        mass_flow_kg_s[k] = hyd.mass_flow_kg_s
    end

    return V11.SimulationResult2D(
        Vector{Float64}(sol.t), nr_rec, base.geometry.receiver_radius,
        grid.r_centers, grid.z_centers, grid.z_rear_centers,
        grid.z_rear_faces, grid.r_centers[1:nr_rec], grid.z_faces,
        solid_T, gas_T, rear_tube_T, rear_gas_T, htc, front_htc,
        front_qgas, gas_rho, gas_velocity, gas_reynolds, gas_prandtl,
        gas_graetz, gas_nusselt, cell_dp, ring_mdot, ring_dp,
        groove_dp, copy(work.groove_overlap), equal_pressure_error,
        receiver_dp_mbar, dp1_prediction_mbar, mass_flow_kg_s, sol,
    )
end

function simulate2D(
    p::ModelParameters2D,
    op::OperatingCondition2D,
    save_times::AbstractVector;
    initial_temperature=Float64(op.ambient_temperature(save_times[1])),
    solver=FBDF(autodiff=AutoFiniteDiff()),
    reltol=1e-6,
    abstol=1e-7,
    dtmax=30.0,
)
    base = p.base
    grid = V11.build_grid2D(base)
    weights = V11.solar_weights2D(grid, base)
    work = V11.Workspace2D(grid, base)
    nbase = grid.nr_total * grid.nz + grid.nz_rear
    base_context = (
        grid=grid, nz_rear=grid.nz_rear, p=base, op=op,
        weights=weights, work=work,
    )
    context = (
        grid=grid, p=p, nbase=nbase, base_context=base_context,
    )
    initial = _initial_augmented_state(
        initial_temperature, p, grid,
    )
    tspan = (Float64(first(save_times)), Float64(last(save_times)))
    problem = ODEProblem(_assembly_ode2D!, initial, tspan, context)
    solution = solve(
        problem, solver; saveat=save_times, reltol=reltol,
        abstol=abstol, dtmax=dtmax,
        isoutofdomain=(u, p, t) -> any(
            x -> isnan(x) || x < 100.0 || x > 3500.0, u,
        ),
    )
    base_result = _base_result_from_solution(
        solution, p, op, grid, work,
    )
    adaptor = [state[nbase + 1] for state in solution.u]
    housing = [state[nbase + 2] for state in solution.u]
    return SimulationResult2D(
        base_result, adaptor, housing, p,
    )
end

function energy_rate_ledger2D(
    u::AbstractVector,
    p::ModelParameters2D,
    op::OperatingCondition2D,
    time::Real=0.0,
)
    base = p.base
    grid = V11.build_grid2D(base)
    nbase = grid.nr_total * grid.nz + grid.nz_rear
    length(u) == nbase + 2 ||
        error("state vector length mismatch in v12 energy ledger")
    weights = V11.solar_weights2D(grid, base)
    work = V11.Workspace2D(grid, base)
    base_context = (
        grid=grid, nz_rear=grid.nz_rear, p=base, op=op,
        weights=weights, work=work,
    )
    context = (
        grid=grid, p=p, nbase=nbase, base_context=base_context,
    )
    du = zeros(Float64, length(u))
    _assembly_ode2D!(du, u, context, Float64(time))

    solid_count = grid.nr_total * grid.nz
    T_s = reshape(view(u, 1:solid_count), grid.nr_total, grid.nz)
    du_s = reshape(view(du, 1:solid_count), grid.nr_total, grid.nz)
    T_tube = view(u, (solid_count + 1):nbase)
    du_tube = view(du, (solid_count + 1):nbase)

    denergy_dt = 0.0
    for i in 1:grid.nr_total, j in 1:grid.nz
        denergy_dt += _cell_capacity_v12(
            i, j, T_s[i, j], p, grid,
        ) * du_s[i, j]
    end
    tube_area = pi * (
        base.geometry.rear_tube_outer_radius^2 -
        base.geometry.rear_tube_inner_radius^2
    )
    for j in 1:grid.nz_rear
        capacity = 3900.0 * V11.alumina_heat_capacity(T_tube[j]) *
                   tube_area * grid.dz_rear
        denergy_dt += capacity * du_tube[j]
    end
    Tad = u[nbase + 1]
    Cad = p.assembly.adaptor_density * adaptor_volume2D(p) *
          V11.alumina_heat_capacity(Tad)
    denergy_dt += Cad * du[nbase + 1]

    a = p.assembly
    extension_length = max(
        0.0,
        a.enclosure_internal_length - base.geometry.length,
    )
    annulus_area = pi * (
        a.enclosure_outer_radius^2 -
        a.enclosure_inner_radius^2
    )
    backplate_area = pi * (
        a.enclosure_outer_radius^2 -
        base.geometry.rear_tube_inner_radius^2
    )
    Vhousing = annulus_area * extension_length +
               backplate_area * a.backplate_thickness
    Chousing = a.aluminum_density * a.aluminum_heat_capacity *
               Vhousing
    denergy_dt += Chousing * du[nbase + 2]

    base_ledger = V11.energy_rate_ledger2D(
        view(u, 1:nbase), base, op, time,
    )
    ambient = Float64(op.ambient_temperature(time))
    Thousing = u[nbase + 2]
    extra_area = (
        2.0 * pi * a.enclosure_outer_radius * extension_length +
        pi * a.enclosure_outer_radius^2
    )
    housing_loss = a.housing_extension_loss_scale * extra_area * (
        10.0 * (Thousing - ambient) +
        base.losses.casing_emissivity * V11.SIGMA *
        (V11.clamp_temperature(Thousing)^4 - ambient^4)
    )
    expected = base_ledger.expected_denergy_dt - housing_loss
    return (
        denergy_dt=denergy_dt,
        expected_denergy_dt=expected,
        residual=denergy_dt - expected,
        housing_extension_loss=housing_loss,
        base_external_rate=base_ledger.expected_denergy_dt,
    )
end

function _filter_observation(time, target, tau)
    tau <= 0.0 && return copy(Float64.(target))
    result = zeros(Float64, length(target))
    result[1] = target[1]
    for i in 2:length(target)
        dt = max(0.0, time[i] - time[i - 1])
        gain = 1.0 - exp(-dt / tau)
        result[i] = result[i - 1] +
                    gain * (target[i] - result[i - 1])
    end
    return result
end

function _sample_local_gas(result::SimulationResult2D, z_target)
    base = result.base_result
    j0, j1, weight = V11._axis_interp_indices(
        base.z_solid, z_target,
    )
    # The innermost annulus represents the measured central channel group.
    values = zeros(length(base.time))
    for k in eachindex(base.time)
        g0 = 0.5 * (
            base.gas_temperature[1, j0, k] +
            base.gas_temperature[1, j0 + 1, k]
        )
        g1 = 0.5 * (
            base.gas_temperature[1, j1, k] +
            base.gas_temperature[1, j1 + 1, k]
        )
        values[k] = (1.0 - weight) * g0 + weight * g1
    end
    return values
end

function sensor_predictions2D(result::SimulationResult2D)
    raw = V11.sensor_predictions2D(result.base_result)
    obs = result.parameters.observation
    time = result.time

    Tgas58 = _sample_local_gas(result, 58.0e-3)
    Tgas107 = _sample_local_gas(result, 107.0e-3)
    wall_fraction = clamp(obs.interior_wall_fraction, 0.0, 1.0)
    T9target = wall_fraction .* raw.T9 .+
               (1.0 - wall_fraction) .* Tgas58
    T10target = wall_fraction .* raw.T10 .+
                (1.0 - wall_fraction) .* Tgas107

    return (
        T8=_filter_observation(
            time, raw.T8, obs.side_time_constant_s,
        ),
        T12=_filter_observation(
            time, raw.T12, obs.side_time_constant_s,
        ),
        T11=_filter_observation(
            time, raw.T11, obs.side_time_constant_s,
        ),
        T9=_filter_observation(
            time, T9target, obs.interior_time_constant_s,
        ),
        T10=_filter_observation(
            time, T10target, obs.interior_time_constant_s,
        ),
        T3=_filter_observation(
            time, raw.T3, obs.outlet_time_constant_s,
        ),
        T2=raw.T2,
        T9_wall=raw.T9,
        T10_wall=raw.T10,
        T9_gas=Tgas58,
        T10_gas=Tgas107,
    )
end

end # module Receiver2D_v12
