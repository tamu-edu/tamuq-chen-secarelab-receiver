# 2D Continuum Solar Receiver Model Development Journal (Macro-ECM)

This journal records the model iterations, mathematical formulations, calibration outcomes, validation metrics, and spatial physical insights for 2D Axisymmetric Continuum Macro-ECM models of the monolithic SiC solar receiver with square channels.

---

## Overall Study Objective

The overarching goal of this study and 2D model development is to **obtain and validate effective macroscopic heat transfer coefficients (convective, radiative, and conductive)** for a structured monolithic solar receiver with square channels. The 2D model serves as a continuum representation (Entire Converter Model / Macro-ECM) where fundamental transport parameters (such as anisotropic effective thermal conductivities $k_{s,r}^{\text{eff}}$ and $k_{s,z}^{\text{eff}}$, Nusselt number correlations, and Rosseland radiation extinction coefficients) are extracted from experimental data, bridging the gap between detailed single-channel physics and full-reactor behavior.

---

## 2026-07-22 — 2D_v1 Axisymmetric Continuum Macro-ECM Model Initial Formulation

Files:
- `2D_v1.jl`
- `run_2D_v1.jl`
- `test/smoke_2D_v1.jl`
- `summaries/2D_v1/`

Mathematical Formulation:
- **2D Solid Continuum $T_s(r, z, t)$**: Discretized on an $N_r \times N_z$ finite-volume radial-axial grid ($8 \times 15$ default).
- **Governing Solid Energy Equation**:
  $$\rho_s (1 - \varepsilon) C_{p,s}(T_s) \frac{\partial T_s}{\partial t} = \frac{1}{r} \frac{\partial}{\partial r} \left( r k_{s,r}^{\text{eff}} \frac{\partial T_s}{\partial r} \right) + \frac{\partial}{\partial z} \left( k_{s,z}^{\text{eff}} \frac{\partial T_s}{\partial z} \right) + q_{\text{solar}}^{\prime\prime\prime}(r, z) - a_v h_c(r, z) (T_s - T_g)$$
- **Gaussian Beam Radial Distribution**:
  $$I_0(r) = I_{\text{peak}} \exp\left(-\frac{r^2}{2\sigma_{\text{beam}}^2}\right)$$
- **Quasi-Steady 2D Fluid Profile $T_g(r, z, t)$**:
  $$\rho_g u_z(r) C_{p,g} \frac{\partial T_g}{\partial z} = a_v h_c(r, z) (T_s - T_g)$$
- **ODE Solver Framework**: Uses `DifferentialEquations.jl` (`OrdinaryDiffEq.jl` / `Rodas5P` / `Tsit5` / `RK4`) for robust state integration.

Validation:
- `test/smoke_2D_v1.jl` test suite passing (25/25 checks).
- Executed `run_2D_v1.jl` across all 15 heating and 3 cooling runs.

Physical Heatmap & Radial Gradient Analysis:
1. **Cause of Homogeneous Heatmaps in Initial Draft**: In the first draft of `2D_v1.jl`, the radial conductivity scale `radial_conductivity_scale` was set to 0.35 multiplying dense SiC conductivity ($k_s \sim 100\text{ W/m/K}$), giving $k_{s,r}^{\text{eff}} \approx 35\text{ W/m/K}$. In a porous square-channel matrix ($\varepsilon \approx 0.65$), radial heat transfer must cross web walls and gas pores, yielding a much lower effective medium radial conductivity ($k_{s,r}^{\text{eff}} \approx 1-4\text{ W/m/K}$). Setting $k_{s,r}^{\text{eff}} \approx 35\text{ W/m/K}$ caused unphysically fast radial thermal smearing, flattening the radial gradient.
2. **Correction ($k_{s,r}^{\text{eff}} \sim 3\text{ W/m/K}$)**: Adjusting `radial_conductivity_scale` to 0.03 ($k_{s,r}^{\text{eff}} \approx 3\text{ W/m/K}$) and applying the Gaussian beam width ($\sigma_{\text{beam}} = 14\text{ mm}$) and delivered power factors ($f_{456}=1.336, f_{304}=1.374, f_{256}=0.786$) reveals a **sharp, realistic 2D thermal gradient** ($T_{\text{core}} > T_{\text{perim}}$ near front, with perimeter heating via rim spillage).

Quantitative Sensor RMSE & Bias Summary (Updated Baseline 2D_v1):

```text
Sensor        Heating RMSE (K)    Heating Steady Error (K)    Cooling RMSE (K)    Cooling Steady Error (K)
T8                 385.6               +486.2                      119.2                -9.4
T12_perim          432.2               +547.4                      162.0               -17.8
T11_perim          607.3               +762.6                      164.3               -27.3
T9_core            487.6               +617.9                      164.4               -20.4
T10_core           625.5               +781.7                      163.9               -31.3
T3                 693.3               +850.5                      162.9               -34.3
T2                  16.9                -25.7                       25.5               -15.0
```

---

## 2026-07-23 — 2D_v2 Multi-Domain Model Formulation & Calibration Results

Physical & Mathematical Formulations Introduced in 2D_v2:
1. **Churchill-Chu Natural Convection Plume Correlation at Front Aperture**:
   - Dynamically calculates front face buoyant plume convection:
     $$\text{Ra}_L = \frac{g \beta (T_{\text{front}} - T_{\text{amb}}) L^3}{\nu \alpha}$$
     $$\text{Nu}_{\text{front}} = \left\{0.825 + \frac{0.387 \text{Ra}_L^{1/6}}{\left[1 + (0.492 / \text{Pr})^{9/16}\right]^{8/27}}\right\}^2$$
     $$h_{\text{front}}(T) = \max\left(10.0, \frac{\text{Nu}_{\text{front}} k_{\text{air}}}{L}\right)$$
   - Dynamically increases front aperture heat dissipation from $10\text{ W/m}^2\text{K}$ at room temperature up to **$35 - 50\text{ W/m}^2\text{K}$ at $1000 - 1300\text{ K}$**.
2. **Temperature-Dependent Alumina Felt Insulation Thermal Conductivity**:
   - Accounts for matrix conduction plus radiative pore transport through alumina fibers:
     $$k_{\text{felt}}(T) = 0.06 + 1.2 \times 10^{-10} T^3\text{ W/m/K}$$
   - Increases radial dissipation into the outer housing from $0.06\text{ W/m/K}$ at $300\text{ K}$ up to **$0.32\text{ W/m/K}$ at $1300\text{ K}$**.
3. **Radial Preferential Flow Distribution**:
   - Preferential central core channel cooling:
     $$\psi(r_i) = 1 - c_{\text{radial\_flow}} \left(\frac{r_i}{R_{\text{rec}}}\right)^2$$
4. **Cooling Transients Initial Condition ($t_0$) & 40 g Solid Mass**:
   - Initializes cooling from solved hot steady state $u_{\text{hot}} = u(t_{\text{end}})$.
   - Corrected porous solid cell volume to $A_{\text{solid}} \cdot dz$, matching exact $40.0\text{ g}$ measured mass.

Quantitative Sensor Metric Summary (2D_v2 Fitted Results):

```text
Phase/Case        Sensor      2D_v1 Fitted RMSE (K)    2D_v2 Fitted RMSE (K)    2D_v2 Steady Error (K)
Heating (E81)     T8                 171.5                    171.5                    +217.2
Heating (E81)     T12_perim           57.0                     57.0                     +70.8
Heating (E81)     T11_perim           27.5                     25.5                     +25.0
Heating (E81)     T9_core             69.1                     64.6                     +84.1
Heating (E81)     T10_core            46.6                     44.6                     +45.6
Heating (E81)     T3                  30.8                     30.0                     +15.9
Heating (E81)     T2                  53.9                     53.9                     +76.7

Cooling (C81)     T8                 101.9                     98.8                      +1.0
Cooling (C81)     T12_perim           61.3                     58.4                      -1.6
Cooling (C81)     T11_perim           26.7                     25.0                      -5.2
Cooling (C81)     T9_core             60.8                     57.8                      -2.5
Cooling (C81)     T10_core            29.5                     27.8                      -6.2
Cooling (C81)     T3                  21.5                     21.4                      -7.1
Cooling (C81)     T2                  41.3                     42.6                     +10.4
```

Generated 2D_v2 Artifacts (`fitted` output only):
- Fitted Metrics CSV: [analysis_results_2D_v2.csv](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v2/analysis_results_2D_v2.csv)
- Fitted Steady State Results CSV: [steady_results_2D_v2.csv](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v2/steady_results_2D_v2.csv)
- Fitted Flow Slopes CSV: [flow_slopes_2D_v2.csv](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v2/flow_slopes_2D_v2.csv)
- Fitted Parameters CSV: [parameters_fitted_2D_v2.csv](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v2/parameters_fitted_2D_v2.csv)
- Fitted Steady Parity Plot: [steady_comparison_2D_v2.png](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v2/plots/steady_comparison_2D_v2.png)
- Fitted Representative Axial Profiles (E67): [axial_profile_E67_2D_v2.png](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v2/plots/axial_profiles/axial_profile_E67_2D_v2.png)
- Fitted Representative 2D Heatmap (E67): [heatmap_2D_E67_2D_v2.png](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v2/plots/2d_profiles/heatmap_2D_E67_2D_v2.png)
- Fitted Cooling Transient Plot (C81): [transient_C81_cooling_2D_v2.png](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v2/plots/transients/transient_C81_cooling_2D_v2.png)
- Complete directory of figures: [summaries/2D_v2/plots/](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v2/plots)








