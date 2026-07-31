module Receiver2D_v21_Properties

export air_conductivity_v21, air_heat_capacity_v21, air_viscosity_v21, air_density_v21
export felt_conductivity_v21, felt_heat_capacity_v21
export sic_conductivity_v21, sic_heat_capacity_v21
export fully_developed_nusselt_v21

const _T_TAB = [300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0]
const _CP_TAB = [1007.0, 1014.0, 1030.0, 1051.0, 1075.0, 1099.0, 1121.0, 1141.0, 1159.0, 1175.0]
const _SIC_CP_TAB = [690.0, 830.0, 960.0, 1040.0, 1090.0, 1120.0, 1140.0, 1160.0, 1170.0, 1180.0]

function _interp(x, xs, ys)
    x <= xs[1] && return ys[1]
    x >= xs[end] && return ys[end]
    idx = searchsortedlast(xs, x)
    idx == length(xs) && return ys[end]
    f = (x - xs[idx]) / (xs[idx+1] - xs[idx])
    return ys[idx] + f * (ys[idx+1] - ys[idx])
end

# 1. Corrected air conductivity (Sutherland, fixes 8% bias)
function air_conductivity_v21(T)
    T_clamp = clamp(T, 100.0, 3500.0)
    return 2.646e-3 * T_clamp^1.5 / (T_clamp + 245.4 * 10^(-12.0/T_clamp))
end

# 2. Corrected air heat capacity (interpolated from reference table)
function air_heat_capacity_v21(T)
    return _interp(T, _T_TAB, _CP_TAB)
end

function air_viscosity_v21(T)
    T_clamp = clamp(T, 100.0, 3500.0)
    return 1.458e-6 * T_clamp^1.5 / (T_clamp + 110.4)
end

function air_density_v21(T, p=101325.0)
    return p / (287.05 * clamp(T, 100.0, 3500.0))
end

# 3. Corrected felt conductivity (v11 functional form, correctly tracking local T)
function felt_conductivity_v21(T)
    return 0.06 + 1.2e-10 * T^3
end

# 4. Corrected felt heat capacity (fixed from 1360 to 1000)
function felt_heat_capacity_v21(T)
    return 1000.0
end

# 5. Corrected SiC capacity and conductivity scale
function sic_conductivity_v21(T, scale=0.55)
    # dense k = 143.0 - 0.09*T + 0.0001*T^2 roughly?
    # Actually, we can use the original dense conductivity but multiplied by porosity scale
    return scale * max(10.0, 206.0 - 0.354 * T + 1.83e-4 * T^2)
end

function sic_heat_capacity_v21(T)
    return _interp(T, _T_TAB, _SIC_CP_TAB)
end

# 6. Corrected Nusselt asymptote for constant wall temperature
function fully_developed_nusselt_v21()
    return 2.98
end

end
