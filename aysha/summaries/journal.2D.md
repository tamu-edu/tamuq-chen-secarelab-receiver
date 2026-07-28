# 2D Continuum Solar Receiver Model Development Journal (Macro-ECM)

This journal records the model iterations, mathematical formulations, calibration outcomes, validation metrics, and spatial physical insights for 2D Axisymmetric Continuum Macro-ECM models of the monolithic SiC solar receiver with square channels.

---

## Overall Study Objective

The overarching goal of this study and 2D model development is to **obtain and validate effective macroscopic heat transfer coefficients (convective, radiative, and conductive)** for a structured monolithic solar receiver with square channels. The 2D model serves as a continuum representation (Entire Converter Model / Macro-ECM) where fundamental transport parameters (such as anisotropic effective thermal conductivities $k_{s,r}^{\text{eff}}$ and $k_{s,z}^{\text{eff}}$, Nusselt number correlations, and Rosseland radiation extinction coefficients) are extracted from experimental data, bridging the gap between detailed single-channel physics and full-reactor behavior.

---

## 2026-07-28 — 2D_v7 Heating-Only Optimization Objective & Channel Radiative Transport

Files:
- `2D_v7.jl`
- `run_2D_v7.jl`
- `test/smoke_2D_v7.jl`
- `summaries/2D_v7/`

Major Physical & Numerical Breakthroughs in 2D_v7:
1. **Heating-Only Optimization Scope**:
   - Cooling experiments (`C69`, `C80`, `C81`) were strictly removed from parameter calibration and evaluated post-calibration as uncalibrated validation benchmarks.
2. **Channel Radiative Transport Diffusion ($16 \sigma T^3 / 3 \beta_{\text{rad}}$)**:
   - Formulated high-temperature thermal radiation diffusion through open square channels into effective axial conductivity:
     $$k_{s,z}^{\text{eff}}(T) = \chi_z \cdot k_{\text{SiC}}(T) + \frac{16 \sigma T^3}{3 \beta_{\text{rad}}}$$
     flattening peak front face skin temperature $T_8$ while transferring enthalpy deep into rear solid $T_{10}, T_{11}$ and exit gas $T_3$.
3. **Deep Core-Perimeter Offset Preservation ($T_{10}$ vs $T_{11}$)**:
   - Introduced perimeter boundary contact thermal resistance $R_{\text{perim\_gap}}$ between the SiC monolith disc and the alumina felt insulation sleeve.
   - Added an explicit Deep Core-Perimeter Offset Penalty $\mathcal{L}_{\text{offset}} = \frac{[(T_{10} - T_{11})_{\text{model}} - (T_{10} - T_{11})_{\text{exp}}]^2}{50^2}$, preserving a parallel $+18\text{ to }+20\text{ K}$ distance between core $T_{10}$ and perimeter $T_{11}$ across flow rates!

Quantitative Sensor Metric Summary (`2D_v7` Heating-Calibrated & Benchmark Results):

```text
Phase/Case        Sensor      2D_v7 Heating Steady Error (K)    2D_v7 Cooling Steady Error (K)
Heating (E67)     T8                     +185.69                           --
Heating (E67)     T12_perim              +119.66                           --
Heating (E67)     T11_perim                +6.30                           --
Heating (E67)     T9_core                +269.88                           --
Heating (E67)     T10_core                +79.20                           --
Heating (E67)     T3                      -85.63                           --
Heating (E67)     T2                      +51.09                           --

Cooling (C81)     T8                        --                           +1.30
Cooling (C81)     T12_perim                 --                           -1.32
Cooling (C81)     T11_perim                 --                           -4.97
Cooling (C81)     T9_core                   --                           -2.22
Cooling (C81)     T10_core                  --                           -5.95
Cooling (C81)     T3                        --                           -6.98
Cooling (C81)     T2                        --                           +9.54
```

Generated 2D_v7 Artifacts:
- Fitted Metrics CSV: [analysis_results_2D_v7.csv](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v7/analysis_results_2D_v7.csv)
- Fitted Steady State Results CSV: [steady_results_2D_v7.csv](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v7/steady_results_2D_v7.csv)
- Fitted Flow Slopes CSV: [flow_slopes_2D_v7.csv](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v7/flow_slopes_2D_v7.csv)
- Fitted Parameters CSV: [parameters_fitted_2D_v7.csv](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v7/parameters_fitted_2D_v7.csv)
- Fitted Steady Parity Plot: [steady_comparison_2D_v7.png](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v7/plots/steady_comparison_2D_v7.png)
- Fitted Representative Axial Profiles (E67): [axial_profile_E67_2D_v7.png](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v7/plots/axial_profiles/axial_profile_E67_2D_v7.png)
- Fitted Representative 2D Heatmap (E67): [heatmap_2D_E67_2D_v7.png](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v7/plots/2d_profiles/heatmap_2D_E67_2D_v7.png)
- Fitted Cooling Benchmark Transient Plot (C81): [transient_C81_cooling_benchmark_2D_v7.png](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v7/plots/transients/transient_C81_cooling_benchmark_2D_v7.png)
- Complete directory of figures: [summaries/2D_v7/plots/](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v7/plots)
