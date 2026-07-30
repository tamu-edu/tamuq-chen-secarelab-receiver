# ============================================================================
# 2D_v18.jl - conservative axial source and exchange-magnitude screen
# ============================================================================

module Receiver2D_v18

include("2D_v17.jl")
const V17 = Receiver2D_v17
const V16 = V17.V16
const V15 = V17.V15
const V14 = V17.V14
const V12 = V17.V12
const V11 = V17.V11

using DifferentialEquations
using OrdinaryDiffEq
using SciMLBase

export V17, V16, V15, V14, V12, V11
export SourceParameters2D, ModelParameters2D, SimulationResult2D
export default_parameters2D, base_grid2D, build_network_grid2D
export network_inventory2D
export simulate2D, sensor_predictions2D, energy_rate_ledger2D
export source_weights2D, source_power_invariant2D
export Geometry2D, SolidProperties2D, HeatTransferParameters2D
export LossParameters2D, OpticalParameters2D, HydraulicParameters2D
export OperatingCondition2D, AssemblyParameters2D, ObservationParameters2D
export NetworkParameters2D, SkinParameters2D
export CasingFlangeParameters2D
export build_initial_state_2D, linear_history, get_t90_2D

const Geometry2D = V17.Geometry2D
const SolidProperties2D = V17.SolidProperties2D
const HeatTransferParameters2D = V17.HeatTransferParameters2D
const LossParameters2D = V17.LossParameters2D
const OpticalParameters2D = V17.OpticalParameters2D
const HydraulicParameters2D = V17.HydraulicParameters2D
const OperatingCondition2D = V17.OperatingCondition2D
const AssemblyParameters2D = V17.AssemblyParameters2D
const ObservationParameters2D = V17.ObservationParameters2D
const NetworkParameters2D = V17.NetworkParameters2D
const SkinParameters2D = V17.SkinParameters2D
const CasingFlangeParameters2D = V17.CasingFlangeParameters2D
const build_initial_state_2D = V17.build_initial_state_2D
const linear_history = V17.linear_history
const get_t90_2D = V17.get_t90_2D

Base.@kwdef struct SourceParameters2D
    model::Symbol = :beer_lambert
    deep_fraction::Float64 = 0.0
    deep_length_m::Float64 = 80.0e-3
end

Base.@kwdef struct ModelParameters2D
    base::V17.ModelParameters2D = V17.default_parameters2D()
    source::SourceParameters2D = SourceParameters2D()
end

function Base.getproperty(p::ModelParameters2D, name::Symbol)
    name === :base && return getfield(p, :base)
    name === :source && return getfield(p, :source)
    return getproperty(getfield(p, :base), name)
end

default_parameters2D() = ModelParameters2D()

function _front_node_count2D(nz)
    nz == 24 && return 4
    nz == 48 && return 8
    nz == 96 && return 16
    return 8
end

function base_grid2D(
    p::ModelParameters2D=default_parameters2D(),
)
    b = p.base.base.base.base.base
    original = V11.build_grid2D(b)
    nz = b.geometry.nodes_z
    nfront = _front_node_count2D(nz)
    nrear = nz - nfront
    nrear > 0 || error("v18 axial mesh requires rear cells")
    front_end = 15.0e-3
    z_faces = vcat(
        collect(range(
            0.0, front_end; length=nfront + 1,
        )),
        collect(range(
            front_end, b.geometry.length; length=nrear + 1,
        ))[2:end],
    )
    dz = diff(z_faces)
    zcenters = 0.5 .* (
        z_faces[1:end-1] .+ z_faces[2:end]
    )
    volume = [
        original.area_frt[i] * dz[j]
        for i in 1:original.nr_total, j in 1:nz
    ]
    return V11.DiscretizationGrid2D(
        original.nr_rec, original.nr_felt,
        original.nr_case, original.nr_total,
        nz, original.nz_rear, dz, original.dz_rear,
        original.r_faces, original.r_centers,
        z_faces, zcenters, original.z_rear_faces,
        original.z_rear_centers, original.area_frt,
        original.area_solid, original.area_flow,
        volume, original.perimeter_ex, original.porosity,
        original.dh, original.total_frt,
    )
end

function build_network_grid2D(
    p::ModelParameters2D=default_parameters2D(),
)
    inherited = V17.build_network_grid2D(p.base)
    values = (
        getfield(inherited, name)
        for name in fieldnames(typeof(inherited))[2:end]
    )
    return V14.NetworkGrid2D(base_grid2D(p), values...)
end

function network_inventory2D(
    p::ModelParameters2D=default_parameters2D(); kwargs...,
)
    return merge(
        V17.network_inventory2D(p.base; kwargs...),
        (
            source_model=p.source.model,
            source_deep_fraction=p.source.deep_fraction,
            source_deep_length_m=p.source.deep_length_m,
        ),
    )
end

function source_weights2D(p::ModelParameters2D, grid)
    reference = V14._network_solar_weights(
        grid, p.base.base.base,
    )
    p.source.model === :beer_lambert && return reference
    p.source.model === :near_deep ||
        error("unknown v18 source model: $(p.source.model)")
    fraction = p.source.deep_fraction
    0.0 <= fraction <= 1.0 ||
        error("v18 deep source fraction must lie in [0,1]")
    length_scale = p.source.deep_length_m
    length_scale > 0.0 ||
        error("v18 deep source length must be positive")
    bg = grid.base_grid
    deep = [
        exp(-bg.z_faces[j] / length_scale) -
        exp(-bg.z_faces[j + 1] / length_scale)
        for j in 1:bg.nz
    ]
    deep ./= sum(deep)
    weights = similar(reference)
    for group in 1:grid.group_count
        group_total = sum(view(reference, group, :))
        near = view(reference, group, :) ./ group_total
        weights[group, :] .= group_total .* (
            (1.0 - fraction) .* near .+ fraction .* deep
        )
    end
    return weights
end

function source_power_invariant2D(
    p::ModelParameters2D=default_parameters2D(),
)
    grid = build_network_grid2D(p)
    reference = V14._network_solar_weights(
        grid, p.base.base.base,
    )
    candidate = source_weights2D(p, grid)
    return (
        reference_sum=sum(reference),
        candidate_sum=sum(candidate),
        absolute_error=abs(sum(reference) - sum(candidate)),
        radial_reference=vec(sum(reference; dims=2)),
        radial_candidate=vec(sum(candidate; dims=2)),
    )
end

struct SimulationResult2D
    base_result::V17.SimulationResult2D
    parameters::ModelParameters2D
    ode_solution::Any
end

function Base.getproperty(result::SimulationResult2D, name::Symbol)
    name === :base_result && return getfield(result, :base_result)
    name === :parameters && return getfield(result, :parameters)
    name === :ode_solution && return getfield(result, :ode_solution)
    return getproperty(getfield(result, :base_result), name)
end

function simulate2D(
    p::ModelParameters2D,
    op::OperatingCondition2D,
    save_times::AbstractVector;
    initial_temperature=Float64(
        op.ambient_temperature(save_times[1]),
    ),
    solver=FBDF(autodiff=AutoFiniteDiff()),
    reltol=1e-6,
    abstol=1e-7,
    dtmax=30.0,
)
    grid = build_network_grid2D(p)
    layout = V15._state_layout(grid)
    work = V14.NetworkWorkspace2D(grid)
    weights = source_weights2D(p, grid)
    base_context = (
        p=p.base.base.base, op=op, grid=grid,
        layout=layout.base, work=work, weights=weights,
    )
    skin_context = (
        p=p.base.base, op=op, grid=grid, layout=layout,
        base_context=base_context, weights=weights,
        skin_group_area=V15._skin_group_areas2D(
            p.base.base, grid,
        ),
    )
    casing_context = (
        p=p.base, op=op, grid=grid, layout=layout,
        skin_context=skin_context,
    )
    initial = V15._initial_state2D(
        initial_temperature, p.base.base, grid,
    )
    problem = ODEProblem(
        V17._casing_flange_ode2D!, initial,
        (Float64(first(save_times)), Float64(last(save_times))),
        casing_context,
    )
    solution = solve(
        problem, solver; saveat=save_times,
        reltol=reltol, abstol=abstol, dtmax=dtmax,
        isoutofdomain=(u, p, t) -> any(
            x -> isnan(x) || x < 100.0 || x > 3500.0, u,
        ),
    )
    base14_result = V15._extract_base_result2D(
        solution, p.base.base, op, grid, layout,
    )
    skin = hcat((
        collect(view(state, layout.skin))
        for state in solution.u
    )...)
    base15_result = V15.SimulationResult2D(
        base14_result, skin, p.base.base, solution,
    )
    base17_result = V17.SimulationResult2D(
        base15_result, p.base, solution,
    )
    return SimulationResult2D(base17_result, p, solution)
end

sensor_predictions2D(result::SimulationResult2D) =
    V17.sensor_predictions2D(result.base_result)

function energy_rate_ledger2D(
    u::AbstractVector,
    p::ModelParameters2D,
    op::OperatingCondition2D,
    time::Real=0.0,
)
    # The inherited v17 ledger reconstructs the v17 axial faces.  Those faces
    # are identical to v18 only on M (48 cells, 8/40 front/rear).  Refuse a
    # numerically false closure report on the new C/F face distributions.
    p.geometry.nodes_z == 48 || error(
        "v18 energy ledger is exact only on the nominal M mesh; " *
        "C/F require the custom-grid ledger",
    )
    inherited = V17.energy_rate_ledger2D(
        u, p.base, op, time,
    )
    invariant = source_power_invariant2D(p)
    return merge(inherited, (
        source_reference_sum=invariant.reference_sum,
        source_candidate_sum=invariant.candidate_sum,
        source_power_error=invariant.absolute_error,
    ))
end

end
