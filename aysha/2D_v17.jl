# ============================================================================
# 2D_v17.jl - verified rear aluminum casing-to-cooling-flange path
# ============================================================================

module Receiver2D_v17

include("2D_v16.jl")
const V16 = Receiver2D_v16
const V15 = V16.V15
const V14 = V16.V14
const V12 = V16.V12
const V11 = V16.V11

using DifferentialEquations
using OrdinaryDiffEq
using SciMLBase

export V16, V15, V14, V12, V11
export CasingFlangeParameters2D, ModelParameters2D, SimulationResult2D
export default_parameters2D, build_network_grid2D, network_inventory2D
export simulate2D, sensor_predictions2D, energy_rate_ledger2D
export hydraulic_feedback_diagnostic2D
export Geometry2D, SolidProperties2D, HeatTransferParameters2D
export LossParameters2D, OpticalParameters2D, HydraulicParameters2D
export OperatingCondition2D, AssemblyParameters2D, ObservationParameters2D
export NetworkParameters2D, SkinParameters2D
export build_initial_state_2D, linear_history, get_t90_2D

const Geometry2D = V16.Geometry2D
const SolidProperties2D = V16.SolidProperties2D
const HeatTransferParameters2D = V16.HeatTransferParameters2D
const LossParameters2D = V16.LossParameters2D
const OpticalParameters2D = V16.OpticalParameters2D
const HydraulicParameters2D = V16.HydraulicParameters2D
const OperatingCondition2D = V16.OperatingCondition2D
const AssemblyParameters2D = V16.AssemblyParameters2D
const ObservationParameters2D = V16.ObservationParameters2D
const NetworkParameters2D = V16.NetworkParameters2D
const SkinParameters2D = V16.SkinParameters2D
const build_initial_state_2D = V16.build_initial_state_2D
const linear_history = V16.linear_history
const get_t90_2D = V16.get_t90_2D

Base.@kwdef struct CasingFlangeParameters2D
    # Effective contact + spreading + water-side conductance.  The aluminum
    # bulk path is already represented by the rear housing state.
    conductance_W_K::Float64 = 0.0
end

Base.@kwdef struct ModelParameters2D
    base::V16.ModelParameters2D = V16.default_parameters2D()
    casing_flange::CasingFlangeParameters2D =
        CasingFlangeParameters2D()
end

function Base.getproperty(p::ModelParameters2D, name::Symbol)
    name === :base && return getfield(p, :base)
    name === :casing_flange && return getfield(p, :casing_flange)
    return getproperty(getfield(p, :base), name)
end

default_parameters2D() = ModelParameters2D()

build_network_grid2D(p::ModelParameters2D=default_parameters2D()) =
    V16.build_network_grid2D(p.base)

function network_inventory2D(
    p::ModelParameters2D=default_parameters2D(); kwargs...,
)
    return merge(
        V16.network_inventory2D(p.base; kwargs...),
        (
            casing_flange_conductance_W_K=
                p.casing_flange.conductance_W_K,
        ),
    )
end

function _casing_flange_power2D(u, p, layout)
    Thousing = u[layout.base.housing]
    return p.casing_flange.conductance_W_K *
           (Thousing - V11.WATER_FLANGE_TEMP)
end

function _casing_flange_ode2D!(du, u, context, time)
    V15._skin_ode2D!(
        du, u, context.skin_context, time,
    )
    Qflange = _casing_flange_power2D(
        u, context.p, context.layout,
    )
    Chousing = V14._housing_capacity2D(
        context.p.base.base,
    )
    du[context.layout.base.housing] -=
        Qflange / max(Chousing, eps(Float64))
    return nothing
end

struct SimulationResult2D
    base_result::V15.SimulationResult2D
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
    weights = V14._network_solar_weights(grid, p.base.base)
    base_context = (
        p=p.base.base, op=op, grid=grid, layout=layout.base,
        work=work, weights=weights,
    )
    skin_context = (
        p=p.base, op=op, grid=grid, layout=layout,
        base_context=base_context, weights=weights,
        skin_group_area=V15._skin_group_areas2D(p.base, grid),
    )
    context = (
        p=p, op=op, grid=grid, layout=layout,
        skin_context=skin_context,
    )
    initial = V15._initial_state2D(
        initial_temperature, p.base, grid,
    )
    problem = ODEProblem(
        _casing_flange_ode2D!, initial,
        (Float64(first(save_times)), Float64(last(save_times))),
        context,
    )
    solution = solve(
        problem, solver; saveat=save_times,
        reltol=reltol, abstol=abstol, dtmax=dtmax,
        isoutofdomain=(u, p, t) -> any(
            x -> isnan(x) || x < 100.0 || x > 3500.0, u,
        ),
    )
    base_result = V15._extract_base_result2D(
        solution, p.base, op, grid, layout,
    )
    skin = hcat((
        collect(view(state, layout.skin))
        for state in solution.u
    )...)
    inherited = V15.SimulationResult2D(
        base_result, skin, p.base, solution,
    )
    return SimulationResult2D(inherited, p, solution)
end

sensor_predictions2D(result::SimulationResult2D) =
    V15.sensor_predictions2D(result.base_result)

function energy_rate_ledger2D(
    u::AbstractVector,
    p::ModelParameters2D,
    op::OperatingCondition2D,
    time::Real=0.0,
)
    inherited = V16.energy_rate_ledger2D(
        u, p.base, op, time,
    )
    grid = build_network_grid2D(p)
    layout = V15._state_layout(grid)
    Qflange = _casing_flange_power2D(u, p, layout)
    return merge(inherited, (
        denergy_dt=inherited.denergy_dt - Qflange,
        expected_denergy_dt=
            inherited.expected_denergy_dt - Qflange,
        residual=(
            inherited.denergy_dt - Qflange -
            (inherited.expected_denergy_dt - Qflange)
        ),
        tube_flange_loss=inherited.flange_loss,
        casing_flange_loss=Qflange,
        flange_loss=inherited.flange_loss + Qflange,
    ))
end

hydraulic_feedback_diagnostic2D(
    p::ModelParameters2D=default_parameters2D(); kwargs...,
) = V16.hydraulic_feedback_diagnostic2D(p.base; kwargs...)

end
