# ============================================================================
# 2D_v16.jl - cooling-identified felt path for the v15 exterior SiC skin
# ============================================================================

module Receiver2D_v16

include("2D_v15.jl")
const V15 = Receiver2D_v15
const V14 = V15.V14
const V12 = V15.V12
const V11 = V15.V11

export V15, V14, V12, V11
export SkinParameters2D, ModelParameters2D, SimulationResult2D
export default_parameters2D, build_network_grid2D, network_inventory2D
export simulate2D, sensor_predictions2D, energy_rate_ledger2D
export hydraulic_feedback_diagnostic2D
export Geometry2D, SolidProperties2D, HeatTransferParameters2D
export LossParameters2D, OpticalParameters2D, HydraulicParameters2D
export OperatingCondition2D, AssemblyParameters2D, ObservationParameters2D
export NetworkParameters2D, build_initial_state_2D, linear_history, get_t90_2D

const SkinParameters2D = V15.SkinParameters2D
const ModelParameters2D = V15.ModelParameters2D
const SimulationResult2D = V15.SimulationResult2D
const Geometry2D = V15.Geometry2D
const SolidProperties2D = V15.SolidProperties2D
const HeatTransferParameters2D = V15.HeatTransferParameters2D
const LossParameters2D = V15.LossParameters2D
const OpticalParameters2D = V15.OpticalParameters2D
const HydraulicParameters2D = V15.HydraulicParameters2D
const OperatingCondition2D = V15.OperatingCondition2D
const AssemblyParameters2D = V15.AssemblyParameters2D
const ObservationParameters2D = V15.ObservationParameters2D
const NetworkParameters2D = V15.NetworkParameters2D
const build_initial_state_2D = V15.build_initial_state_2D
const linear_history = V15.linear_history
const get_t90_2D = V15.get_t90_2D

default_parameters2D() = V15.default_parameters2D()
build_network_grid2D(args...; kwargs...) =
    V15.build_network_grid2D(args...; kwargs...)
network_inventory2D(args...; kwargs...) =
    V15.network_inventory2D(args...; kwargs...)
simulate2D(args...; kwargs...) = V15.simulate2D(args...; kwargs...)
sensor_predictions2D(args...; kwargs...) =
    V15.sensor_predictions2D(args...; kwargs...)
energy_rate_ledger2D(args...; kwargs...) =
    V15.energy_rate_ledger2D(args...; kwargs...)
hydraulic_feedback_diagnostic2D(args...; kwargs...) =
    V15.hydraulic_feedback_diagnostic2D(args...; kwargs...)

end
