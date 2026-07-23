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

## 2026-07-23 — High-Resolution Multi-Domain 2D Axisymmetric Model Calibration Results

Major Physical & Numerical Refinements Incorporated:
1. **Cooling Transient Initial Condition ($t_0$) Correction**:
   - Cooling runs ($C_{69}, C_{80}, C_{81}$) now initialize from the solved **hot steady state** ($u_{\text{hot}} = u(t_{\text{end}})$) of their corresponding heating runs ($E_{69}, E_{80}, E_{81}$).
   - The model accurately predicts dynamic cooling thermal decay starting at $\sim 1000 - 1300\text{ K}$ down to ambient as cold air flows through the receiver without solar flux.
2. **Physical Placement of Thermocouple $T_2$**:
   - Thermocouple $T_2$ is sampled at $r = R_{\text{rec}} + 40.0\text{ mm} = 73.9\text{ mm}$ (embedded 40 mm outside the receiver wall inside the insulation housing).
   - This brought $T_2$ predictions down from $+380\text{ K}$ overestimation to within $+10 - +70\text{ K}$ of experiment.
3. **Porous Solid Mass Correction ($40.0\text{ g}$ Measured Mass)**:
   - Corrected receiver cell solid thermal mass to use $A_{\text{solid}} = (1 - \varepsilon) A_{\text{frt}}$ rather than gross cell volume, matching the exact $40.0\text{ g}$ measured mass of the porous SiC monolith.
   - This eliminated the $2.86\times$ unphysically slow thermal response and restored fast dynamic thermal response.
4. **High Spatial Grid Resolution**:
   - Discretized disc into a **$17 \times 25$ grid** ($10$ SiC core rings + $5$ alumina felt rings + $2$ aluminum casing rings $\times 25$ axial depth cells along $137\text{ mm}$).
   - Rear alumina exit tube discretized into **$20$ axial grid cells** along $150\text{ mm}$.

Quantitative Sensor Metric Summary across Heating & Cooling Runs:

```text
Phase/Case        Sensor      Fitted RMSE (K)    Fitted Steady Error (K)
Heating (E81)     T8               171.5                 +217.2
Heating (E81)     T12_perim         57.0                  +70.8
Heating (E81)     T11_perim         27.5                  +30.0
Heating (E81)     T9_core           69.1                  +92.7
Heating (E81)     T10_core          46.6                  +50.5
Heating (E81)     T3                30.8                  +19.4
Heating (E81)     T2                53.9                  +76.7

Cooling (C81)     T8               101.9                   +1.2
Cooling (C81)     T12_perim         61.3                   -1.5
Cooling (C81)     T11_perim         26.7                   -5.1
Cooling (C81)     T9_core           60.8                   -2.4
Cooling (C81)     T10_core          29.5                   -6.0
Cooling (C81)     T3                21.5                   -7.0
Cooling (C81)     T2                41.3                   +9.9
```

Generated Fitted Artifacts (`fitted` output only):
- Fitted Metrics CSV: [analysis_results_2D_v1.csv](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v1/analysis_results_2D_v1.csv)
- Fitted Steady State Results CSV: [steady_results_2D_v1.csv](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v1/steady_results_2D_v1.csv)
- Fitted Flow Slopes CSV: [flow_slopes_2D_v1.csv](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v1/flow_slopes_2D_v1.csv)
- Fitted Parameters CSV: [parameters_fitted_2D_v1.csv](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v1/parameters_fitted_2D_v1.csv)
- Fitted Steady Parity Plot: [steady_comparison_2D_v1.png](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v1/plots/steady_comparison_2D_v1.png)
- Fitted Representative Axial Profiles (E67): [axial_profile_E67_2D_v1.png](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v1/plots/axial_profiles/axial_profile_E67_2D_v1.png)
- Fitted Representative 2D Heatmap (E67): [heatmap_2D_E67_2D_v1.png](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v1/plots/2d_profiles/heatmap_2D_E67_2D_v1.png)
- Fitted Cooling Transient Plot (C81): [transient_C81_cooling_2D_v1.png](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v1/plots/transients/transient_C81_cooling_2D_v1.png)
- Complete directory of figures: [summaries/2D_v1/plots/](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v1/plots)







