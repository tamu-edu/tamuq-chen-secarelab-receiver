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

---

## 2026-07-28 - Retrospective error audit of 2D_v1 through 2D_v7

> **Superseding notice:** The historical v7 claims above are not accepted
> validation evidence. The v7 run ended with `MaxTime`, most parameters did
> not move from their seeds, and all prior versions inherited a receiver
> frontal-area error. The old outputs are retained only as development history.

### Audit conclusion

None of the fitted coefficients or manuscript-readiness claims from
`2D_v1`-`2D_v7` should be treated as validated. The numerical LTNE scaffold
was useful, but a fundamental receiver-area error propagated into geometry,
mass, porosity, incident power and fitted transport parameters. Additional
conservation, observation-map and calibration defects prevented the old fits
from being physically interpretable.

`2D_v8` is the first version in this line that enforces geometry, solar-power
and mass-flow closure with executable tests.

### Numbered error ledger

| ID | Affected model/manual statement | Error and quantitative consequence | Correction in `2D_v8` | Verification |
|---|---|---|---|---|
| E01 | `2D_v1`-`2D_v7`; v2 theory manual Sec. 2 | The area of a 19 mm by 19 mm face was coded/documented as `3.61e-3 m2`. The correct value is `3.61e-4 m2`. The area-equivalent radius is 10.7196 mm, not 33.9 mm. | `receiver_width=19e-3`; `receiver_radius=sqrt(receiver_width^2/pi)`; grid construction rejects inconsistent values. | Exact area/radius tests pass. |
| E02 | v1-v7 geometry and power source | E01 forced porosity to about 0.0623 instead of 0.6233, receiver mass to about 0.997 kg instead of 0.04006 kg, and aperture power to ten times its physical value. Every old time constant and coefficient was conditional on the wrong plant. | Porosity derives from 100 channels of 1.5 mm square area; solid mass derives from corrected solid area, 137 mm length and 2150 kg/m3 density; power uses corrected aperture area. | Tests give porosity 0.6232687 and mass 0.0400588 kg. |
| E03 | v1-v7 axisymmetric domain | Felt began outside the erroneous 33.9 mm receiver, so receiver/felt volumes and radial paths did not represent the stated monolith. | Receiver ends at 10.7196 mm; felt and casing retain measured 57 and 75 mm outer radii. | Grid and mesh tests pass. |
| E04 | v2-v7 active-flow model | Only a Reynolds-dependent fraction of measured mass flow was heated in the receiver, but the rear tube used full measured flow starting from the active-flow outlet. This violates gas enthalpy conservation at the junction. | Full measured mass flow is conservatively partitioned among rings and the same total continues through the rear tube. Radial preference changes only partition. | Global energy-rate residual is below `1e-8 W`; zero-source equilibrium passes. |
| E05 | v2-v7 cavity state and `C_cavity` | An extra cavity state and fitted heat capacity were present, but its derivative was zero and it did not close an energy balance. The parameter was inert. | Remove the cavity state and parameter; retain distributed felt, casing and tube capacities. | State-length test confirms no inert state. |
| E06 | v1-v7/manual Beer-Lambert source | Axial absorption weights were normalized over 137 mm, forcing all non-front/non-spill power to be absorbed. This erased transmission and made `beta_opt` mainly a shape parameter. | Do not normalize finite-length Beer-Lambert fractions; retain `exp(-beta_opt L)` as transmitted power. | Beta monotonicity and exact solar-closure tests pass. |
| E07 | v1-v7 spillage source | Rim spillage was deposited through every felt axial cell, creating an artificial volumetric heater. | Deposit spillage only in the illuminated front felt/rim cell, area-weighted. | Solar ledger separates core, rim and transmission terms. |
| E08 | v1-v7 SiC conductivity | The old correlation produced unrealistically large values and differed from the established 0D/1D property basis. | Use the established polynomial: about 116.6 W/m/K at 300 K and 61.8 W/m/K at 1200 K. | Range and monotonicity tests pass. |
| E09 | v7 axial conductivity | v7 used `kz_scale*prop_j.k` after `prop_j.k` already included the radial scale. Axial conductivity depended on both `chi_z` and `chi_r`, contrary to its equation. | Use `chi_z*k_SiC(T) + 16*sigma*T^3/(3*beta_rad)` axially and `chi_r*k_SiC(T)` radially. | Code path, energy and mesh tests pass. |
| E10 | v7 perimeter gap | The gap term was an unqualified additive number rather than dimensional area-normalized resistance, so its magnitude was not an interpretable coefficient. | Use finite-volume half-cell resistances plus `R''_contact` in `m2 K/W`; fix contact resistance in the first v8 fit to reduce confounding. | Units are explicit; fit-mask test confirms it is fixed. |
| E11 | v1-v7 rear tube | Rear gas/tube areas used 25/19 mm diameters, inconsistent with documented 6.5/8.0 mm inner/outer radii. | Derive all rear areas from 6.5 and 8.0 mm radii. | Smoke and energy tests exercise the corrected geometry. |
| E12 | v1-v7 flange boundary | The water-cooled flange sink acted on every rear-tube cell, a distributed sink rather than a terminal boundary. | Apply it only to the final rear cell with a half-cell conduction distance. | Adiabatic uniform-state test passes when disabled. |
| E13 | v1-v7 observation map | Nearest-cell sampling, an off-by-one T3 location, and geometry-inconsistent T2/receiver positions biased comparisons and made them mesh-dependent. | Bilinear solid interpolation; rear-gas-face interpolation for T3 at 3 mm; T2 sampled 40 mm outside the corrected receiver radius. | Sensor-map and T3 identity tests pass. |
| E14 | v1-v7 cooling initialization | The constructed cooling field did not reproduce measured solid sensors through the observation map, adding a fictitious initial residual. | Constrain discrete axial profiles so the initial map exactly returns T8, T12, T11, T9 and T10. | Five tests pass to `1e-8 K`. T3 remains an algebraic gas output, not an independently assignable state. |
| E15 | v1-v7 response-time metric | `t90` searched only for rising crossings, so cooling response times were wrong or defaulted to record duration. | Select the crossing direction from the sign of `T_final-T_initial`. | Rising and falling tests pass. |
| E16 | v1-v7 calibration | High-dimensional single-start fits used only three heating cases, routinely ended at `MaxTime`, and often returned seed/bound-active values. v7 ended at objective 1.94295 with most entries unchanged. | Production v8 fits eight transport/power terms; fixes the independently reported `B_Re=1.44`, beam width, spillage and contact resistance; uses all 15 heating cases; records evaluations and return code; validates on nominal mesh. | Fit-mask tests pass; pilot and rejection evidence are recorded below. |
| E17 | v6/v7 interpretation | Aggregate v6 evidence contradicted readiness claims: heating RMSE mean 92.8 K, cooling RMSE mean 71.0 K, heating steady MAE 80.5 K and cooling steady MAE 20.0 K. All 15 heating cases had wrong mid/deep core-perimeter ordering, and 8/21 flow-slope signs were wrong. v7 still had errors up to about +270 K. | Report aggregate metrics and both spatial offsets/flow-slope signs without equating execution with validation. | Baseline and fitted metrics are separated. |
| E18 | v1-v7 tests | Smoke tests checked execution/shapes but not geometry, conservation, transmission, sensor coordinates, cooling timing or grid sensitivity. | Add geometry, property, solar, energy-rate, equilibrium, sensor-map, cooling-t90 and three-mesh tests. | 53 checks pass: 34 smoke, 14 physics, 5 mesh. |
| E19 | v1-v7 solver cost | Dense finite-difference Jacobians made each stiff solve expensive and contributed to time-limited fits. | Use FBDF finite differences. On a 600 s nominal run it matched Rodas5P within 0.00021 K and reduced 8.897 s to 1.557 s. | All suites pass with FBDF. |
| E20 | v1-v7 and first v8 attempt, convective closure | The Macro-ECM imposed the single-channel fully developed floor `Nu=3.61` even though the independently reported whole-receiver apparent correlation gives Nu of order 0.03-0.2 in this Reynolds range. The first v8 attempt also added the developing term to the floor. Gas/solid NTU saturated, so fitted `A_Nu` and `B_Re` had exactly zero local sensitivity across the complete transient response. | Define the fitted closure explicitly as an apparent Macro-ECM Nusselt law, remove the single-channel floor, use `max(Nu_floor,Nu_dev)` with `Nu_floor=0`, seed `A_Nu=1e-3`, and bound it in `[1e-5,0.1]`. Add a regression test that doubling `A_Nu` changes receiver gas heating. | Apparent-Nu sensitivity test passes; the pilot fit moves both convective parameters to interior values. |

### Corrections required in `2D_v2_theory_manual.md`

The v2 manual is historical and must not be reused as the current model:

1. Replace the 33.9 mm receiver radius and `3.61e-3 m2` area by 10.7196 mm
   and `3.61e-4 m2`.
2. Remove active-flow fraction unless a bypass stream and mixing balance are
   modeled explicitly.
3. Remove the normalized Beer-Lambert denominator and report transmission.
4. Separate optical `beta_opt` from thermal-radiation `beta_rad`.
5. Remove the inert cavity state and claimed 301 J/K inference.
6. Replace nearest-cell sensor assignments and old T3 position by the v8 map.
7. Replace old rear-tube dimensions and distributed flange sink.
8. Remove the local `Nu=3.61` floor if the fitted coefficient is interpreted
   as the observed whole-receiver apparent closure.
9. Treat v2/v3 coefficient and validation tables as invalidated by E01-E20.

### Approximation status after correction

V8 remains an area-equivalent axisymmetric representation of a square
monolith. Radial flow preference, beam width, spillage, front deposition,
contact resistance and external-loss closures remain effective assumptions,
not independently measured fields. Conservation and mesh tests establish
numerical/bookkeeping integrity; they do not alone validate coefficients.

---

## 2026-07-28 - 2D_v8 corrected-geometry conservative model

Files:

- `2D_v8.jl`
- `run_2D_v8.jl`
- `test/smoke_2D_v8.jl`
- `test/check_2D_v8_physics.jl`
- `test/check_2D_v8_mesh.jl`
- `test/check_2D_v8_identifiability.jl`
- `test/diagnose_2D_v8_fitted_mesh.jl`
- `test/diagnose_2D_v8_radial_ordering.jl`
- `test/pilot_calibrate_2D_v8.jl`
- `test/benchmark_2D_v8_solvers.jl`
- `test/benchmark_2D_v8.jl`
- `summaries/2D_v8/`

Pre-calibration verification:

- 34/34 smoke/regression checks passed.
- 14/14 physical-invariant and observation-map checks passed.
- 5/5 mesh-sensitivity checks passed.
- At 600 s with the corrected apparent-Nu closure, maximum
  coarse-to-nominal sensor change was 15.33 K and nominal-to-fine change was
  11.12 K. The nominal mesh is `(14,7,3,45,30)` and the refinement gate passes.
- Solar power closure and instantaneous global energy-rate closure pass.
- The current apparent-Nu defaults are deliberately uncalibrated: heating
  sensor-RMSE mean is 104.98 K and steady MAE is 99.54 K; cooling
  sensor-RMSE mean is 37.29 K and steady MAE is 14.82 K.

### Rejected v8 fit A: inherited `Nu=3.61` saturation

This run is preserved in the primary `summaries/2D_v8` CSV/plot artifacts so
that the failed hypothesis remains auditable:

- Objective: 2.00928; return code `MaxIters`; 120 evaluations.
- Heating RMSE mean: 118.70 K; heating steady MAE: 123.68 K.
- Cooling RMSE mean: 68.50 K; cooling steady MAE: 19.82 K.
- Mid-depth and deep core/perimeter signs: 0/15 correct for each.
- Flow-slope signs: 20/21 correct, but several magnitudes were far too steep.
- `A_Nu=5.976` and `B_Re=1.790` were near their upper bounds;
  `beta_opt=300 1/m` was at its upper bound.
- Full transient local sensitivity rank was only 7/9; condition number was
  `3.32e17`. `A_Nu` and `B_Re` sensitivity columns were exactly zero because
  exchange effectiveness was saturated.
- A separate fitted-point synthetic mesh gate passed (3.27 K
  coarse-to-nominal and 3.21 K nominal-to-fine), proving that the bad absolute
  RMSE was caused by the closure/loss formulation rather than mesh transfer.

This vector is rejected and must not be quoted as extracted transport
coefficients.

### Corrected v8 apparent-Nu pilot

The corrected closure and absolute-error-aware loss were exercised on six
heating cases spanning the low/high flow endpoints at all three irradiances.
The pilot used 60 evaluations on the former nominal mesh and was validated on
the refined nominal mesh across all 15 heating and three cooling cases.

Pilot fitted vector:

```text
A_Nu       = 0.00200748
B_Re       = 1.11359
chi_r      = 0.0536842
chi_z      = 0.908613
scale_456  = 0.858066
scale_304  = 0.776186
scale_256  = 0.516722
beta_rad   = 54.0814 1/m
beta_opt   = 96.1096 1/m
```

All nine pilot-fitted parameters were interior. Full-transient local
sensitivity improved to rank 8/9 with condition number 1055.7, but
`corr(A_Nu,B_Re)=0.9934`. The production v8 fit mask therefore fixes the
independently reported apparent exponent at `B_Re=1.44` and fits only `A_Nu`;
beam width, spillage and contact resistance also remain fixed. Quantitative
pilot results:

```text
objective                  12.6711 (MaxIters, 60 evaluations)
heating RMSE mean          76.84 K
heating steady MAE         85.36 K
heating t90 MAE           677.53 s
cooling RMSE mean          36.93 K
cooling steady MAE         11.50 K
cooling t90 MAE           440.48 s
mid ordering correct        0/15
deep ordering correct       0/15
```

The pilot improves heating RMSE relative to v6 (92.8 K) and substantially
improves cooling RMSE, while restoring identifiable convective leverage.
It is nevertheless rejected as a final coefficient set because spatial
ordering fails every case and steady heating error remains large.

### Remaining structural failure and expert status

For E81, increasing central-flow preference from 0.1 to 0.9 changed the
predicted mid-depth perimeter-minus-core offset only from -2.34 K to -1.58 K,
while the measured offset is +21.40 K. The deep measured offset is +21.21 K
and also remains negative in the model. Therefore the current centered
Gaussian source plus radial conduction/flow topology cannot generate the
observed inversion.

Do not add an unconstrained rim heater merely to force the sign. Before a
publishable next fit, obtain or justify at least one of:

1. an independently measured radial flux map or defensible annular/rim
   absorption model;
2. a conservative manifold/bypass model with measured flow partition;
3. confirmation of T9/T12 and T10/T11 physical placement and labels.

After that structural choice, run the all-15-case nominal-grid fit, repeat
mesh and full-transient sensitivity diagnostics, and require interior
parameters plus correct spatial ordering before extracting coefficients.

Current verdict: **v8 implementation PASS; v8 coefficient validation FAIL**.
See `summaries/2D_v8/acceptance_status_2D_v8.txt`.
