# 1D_v10 Revision Notes and Strategy

This document outlines the fundamental physical and architectural revisions planned for the `1D_v10` receiver model. Following a review of recent experimental discrepancies and relevant literature on multi-scale monolith modeling, v10 will shift away from empirical error-correction scalars towards a strict continuum physical model (Entire Converter Model, ECM). 

The ultimate goal of this study is to **obtain and validate effective macroscopic heat transfer coefficients (convective, radiative, and conductive)** for the structured monolithic receiver. The model serves to bridge the gap between detailed channel-scale physics and full-reactor behavior.

## 1. Strict Physical Properties (Removal of Empirical Scalars)
Previous versions introduced arbitrary scaling factors to correct for transient cooling or steady-state heat loss misalignments. In v10, we will strictly enforce physically grounded material definitions:
- **`k_scale` is removed:** The axial effective solid conductivity will no longer be multiplied by an arbitrary scalar. It will rely purely on the temperature-dependent dense material property $k_s(T)$ and the solid volume fraction $(1-\phi)$.
- **`k_ins_scale` is removed:** The radial conductance to the cavity will be fixed by the geometric dimensions and literature thermal conductivity of the insulation.
- **`f_exchange` is removed:** The gas-solid heat exchange area will be strictly bounded by the physical geometric perimeter ($P_{exchange}$).

## 2. Radiative Heat Transfer: Rosseland Diffusion Approximation
Internal thermal radiation significantly enhances effective axial heat transfer at high temperatures, which smooths the axial temperature gradients (as highlighted in Hayes, 2021). Rather than artificially inflating `k_scale`, v10 will explicitly model internal radiation using the **Rosseland diffusion approximation** (Filho, 2020; Mendes, 2014).

The effective solid axial conductivity will be modeled as the sum of pure solid conduction and a radiative contribution:
$$k_{eff,s}(T_s) = (1 - \phi) k_s(T_s) + k_{rad}(T_s)$$
where the radiative conductivity is given by:
$$k_{rad}(T_s) = \frac{16 \sigma T_s^3}{3 \beta_R}$$
Here, $\sigma$ is the Stefan-Boltzmann constant, and $\beta_R$ is the **Rosseland extinction coefficient**. 
- $\beta_R$ will be introduced as a primary physically-meaningful fitting parameter in the optimization loops. Typical values for such structures range between $300$ and $2700 \text{ m}^{-1}$ depending on porosity and pore/channel scale.

## 3. Optimization Strategy Refactoring
Without the artificial parameters (`gamma_C`, `k_scale`, etc.), the calibration stages must be redesigned around fundamental transport physics:
1.  **Cooling Stage (Thermophysical Calibration):** 
    During cooling, there is no forced flow and no solar source. The transient thermal decay is governed exclusively by sensible heat capacity, axial conduction, natural convection losses, and internal radiation. 
    - **Fitted Parameter:** $\beta_R$ (Rosseland extinction coefficient). 
    - This stage isolates the effective macroscopic thermal transport of the solid matrix.
2.  **Heating Stage (Convection and Source Calibration):** 
    During heating, forced internal convection and external solar irradiation dominate. 
    - **Fitted Parameters:** $A_{Nu}$, $B_{Re}$, $C_{Pr}$ (developing-flow Nusselt correlation), and solar source distribution variables (e.g., `front_dep`, `beta_opt`, and any irradiance correction factors).
    - This stage extracts the convective heat transfer coefficients now that the baseline solid thermal backbone is correctly parameterized.

## 4. Code Cleanup and Legacy Component Removal
To cleanly implement the ECM approach, the `1D_v10.jl` code will be purged of all unused and legacy artifacts inherited from v3 through v5. The following components will be entirely removed from the parameter vectors, structs, and code logic:
- `gamma_C` (capacity multiplier)
- `U_side` and `U_rear` (legacy arbitrary boundary losses)
- `h_floor` and `L_h` (legacy empirical axial heat-exchange shape)
- `tau_T3` (legacy temporal lag for gas measurement)
- Hardcoded legacy constants like `B_RE_FIXED_V5` and `C_PR_FIXED_V5`

By enforcing these fundamental constraints, `1D_v10` will ensure that the optimized parameters are not merely mathematical curve fits, but valid, interpretable transport coefficients for square-channel monolith reactors.
