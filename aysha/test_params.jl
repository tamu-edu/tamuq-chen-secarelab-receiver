include("run_2D_v22.jl")

# 1. Test parameter construction logic in run_2D_v22
p_base = parameters_2D_v22()
println("Default Nu: ", p_base.base.base.base.base.base.heat_transfer.fully_developed_nusselt)

# Update global Nu
Main.Receiver2D_v22.Receiver2D_v22_Properties.GLOBAL_NU_MULTIPLIER[] = 5.0

# Re-create params
p_new = parameters_2D_v22()
println("Nu after multiplier change (should be 5 * 3.66 = 18.3): ", p_new.base.base.base.base.base.base.heat_transfer.fully_developed_nusselt)

# 2. Look at calibrate_2D_v22_phase4.jl
include("calibrate_2D_v22_phase4.jl")
p_cal = build_base_parameters(5.0, 5.0)
println("Nu in calibrate phase 4 (should be 5 * 3.66 = 18.3): ", p_cal.base.base.base.base.base.base.heat_transfer.fully_developed_nusselt)
