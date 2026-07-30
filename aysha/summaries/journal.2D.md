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

---

## 2026-07-29 - Implemented 2D_v9 temperature-dependent velocity and DP1 hydraulic extension

> **Acceptance notice:** `2D_v9` now implements the formulation below and
> inherits the conservative thermal and optical equations of `2D_v8`.
> Reference-flow conversion, T8 mapping, mass conservation, local velocity,
> pressure scaling and axial-mesh invariance pass their tests. The cold DP1
> closure is accepted as a hydraulic constraint. The hot-path DP1 closure is
> not yet accepted, and no new thermal transport coefficient has been fitted
> or accepted in this step.

### Verified thermocouple-coordinate correction

The receiver thermocouple axial planes have now been independently verified
as 5, 58 and 107 mm from the illuminated front face. Consequently, the T8
perimeter-skin coordinate inherited from v8 must be corrected from 11 mm to
5 mm in both the observation map and cooling-field initialization. The
mid-depth and deep planes remain:

```text
front/perimeter T8                         z =   5 mm
mid-depth core T9 and perimeter T12        z =  58 mm
deep core T10 and perimeter T11            z = 107 mm
```

This is a coordinate correction, not a fitted parameter. It may materially
affect inference of the approximately 10 mm optical absorption length, while
the independently verified 58 and 107 mm locations preserve the evidence for
the unresolved perimeter-hotter-than-core spatial inversion.

### MFC quantity and reference conditions

The Aalborg GFC manual reports calibration reference conditions of
14.7 psia (approximately 101.4 kPa absolute) and 70 degF (21.1 degC or
294.25 K). The experimental MFC values have already been adjusted using the
manufacturer's gas correction factors and were checked separately with the
bubble flow meter. They are therefore treated as air-equivalent standard
volumetric flow and must not receive a second gas-factor correction.

Using the ideal-gas relation for air at the manual reference condition gives

$$
\rho_{\mathrm{ref}} =
\frac{p_{\mathrm{ref}}}{R_{\mathrm{air}}T_{\mathrm{ref}}}
\simeq 1.200\ \mathrm{kg\,m^{-3}},
$$

and the prescribed receiver mass flow is

$$
\dot m_{\mathrm{MFC}} =
\rho_{\mathrm{ref}}\frac{Q_{\mathrm{std}}}{60000}.
$$

`2D_v9` retains this prescribed mass flow as its primary flow input. The
hydraulic extension does not infer an unmeasured per-run bypass branch.

### Temperature-dependent density and local velocity

The v8 Reynolds-number expression was already consistent with prescribed mass
flow because density cancels from
$\mathrm{Re}=\rho uD_h/\mu=\dot mD_h/(A_{\mathrm{flow}}\mu)$. It did not,
however, calculate the actual local channel velocity or pressure loss.
`2D_v9` explicitly evaluates

$$
\rho_{i,j}(t)=
\frac{p_{i,j}(t)}{R_{\mathrm{air}}T_{g,i,j}(t)}
\simeq
\frac{p_{\mathrm{atm}}}{R_{\mathrm{air}}T_{g,i,j}(t)}
$$

and

$$
u_{i,j}(t)=
\frac{\dot m_i}{\rho_{i,j}(t)A_{\mathrm{flow},i}}.
$$

The atmospheric-pressure approximation is adequate for the first hydraulic
implementation because the measured DP1 magnitude is only a few millibar.
Mass flow remains fixed, while hot-gas density decreases and local actual
volumetric flow and velocity increase.

### DP1 geometry and cold-$t_0$ empirical constraint

DP1 is a flush wall static-pressure tap in the flow path just inside the
water jacket. Its reference port and the receiver front are both open to
atmosphere. DP1 therefore measures local static gauge pressure relative to
the atmospheric receiver face; it is not a Pitot or stagnation-pressure
measurement.

Nine near-ambient heating starts were selected as the cold-$t_0$ hydraulic
set:

```text
E67, E68, E70, E72, E74, E75, E76, E78, E80
```

Flow, DP1, T3, the five receiver temperatures and ambient temperature are
averaged over the raw inclusive interval $0\le t\le10$ s before data
decimation. A case is selected only when both T3 and the mean receiver-solid
temperature are within 5 K of ambient.

The regression against corrected standard MFC flow is

$$
DP1_{\mathrm{raw}}[\mathrm{mbar}]
= -0.614226
+ 0.0455545\,Q_{\mathrm{std}}[\mathrm{L\,min^{-1}}],
$$

with

```text
R2    = 0.98144
RMSE  = 0.02768 mbar
```

The intercept is an empirical sensor/path zero for this dataset, not a
separately verified transducer-zero measurement. At the MFC reference
temperature of 294.25 K, the ideal fully developed square-channel slope is
0.0233408 mbar/(standard L/min). The observed cold hydraulic slope is therefore
1.95171, approximately 1.95 times, the ideal receiver-channel-only estimate.
This multiplier plausibly aggregates
developing-channel, inlet, outlet and short water-jacket path losses contained
between the atmospheric face and the DP1 tap.

### Pressure-drop closure and use of DP1

For laminar flow in a square channel,

$$
f_D\mathrm{Re}=56.91,
$$

so the axial-cell channel pressure loss may be written as

$$
\Delta p_{i,j} =
f_{D,i,j}\frac{\Delta z_j}{D_h}
\frac{\rho_{i,j}u_{i,j}^2}{2}
=
28.455\,
\frac{\mu_{i,j}u_{i,j}\Delta z_j}{D_h^2}.
$$

The implemented v9 comparison uses the geometry/property prediction multiplied
by the cold empirical hydraulic factor:

$$
DP1_{\mathrm{model}}
= p_{0,\mathrm{DP1}}
+ C_{\mathrm{hyd}}\Delta p_{\mathrm{square}},
\qquad
p_{0,\mathrm{DP1}}\simeq-0.614226\ \mathrm{mbar},
\qquad
C_{\mathrm{hyd}}\simeq1.95,
$$

with the recorded DP1 sign convention preserved. Additional quadratic minor
losses are to be introduced only if the heating and cooling residuals require
them; they are not part of the initial cold calibration.

DP1 is an independent hydraulic validation output and, if supported by the
full histories, may constrain one common global flow scale multiplying the
corrected MFC mass flow. It must not be used to fit an independent flow or
bypass fraction for every run. This ordering prevents flow uncertainty from
being absorbed into the apparent Nusselt prefactor and exponent.

### Bubble-meter evidence and exclusion

The bubble flow-meter characterization was performed **without the receiver**.
Its restricted-flow response therefore characterizes the MFC/bubble-meter test
arrangement and cannot calibrate flow through the installed receiver. It is
excluded from the v9 receiver-flow calibration. The corrected MFC standard
flow remains the nominal input, with the normal-configuration DP1 data used to
test its consistency.

### Implemented files and verification results

The implementation is in `2D_v9.jl`; `run_2D_v9.jl` imports DP1, performs the
cold-$t_0$ regression and runs the 15 heating plus 3 cooling cases. The
importer now exposes the original spreadsheet DP1 column without changing the
legacy observation indices.

The following checks passed on 2026-07-29:

```text
test/smoke_2D_v9.jl          34/34
test/check_2D_v9_physics.jl  46/46
test/check_2D_v9_mesh.jl       2/2
test/smoke_2D_v8.jl          34/34
```

The tests cover the MFC reference density and mass-flow conversion, corrected
T8 mapping and cooling initialization, exact ring mass-flow closure,
$\rho uD_h/\mu$ consistency, zero-flow behavior, velocity and pressure
temperature response, thermal energy conservation, preservation of the
12-parameter v8 thermal pack, DP1 import/calibration and hydraulic axial-mesh
invariance.

The full-run artifacts are:

```text
summaries/2D_v9/dp1_cold_t0_calibration_2D_v9.csv
summaries/2D_v9/dp1_summary_2D_v9.csv
```

### Full-run DP1 assessment

With the corrected MFC mass flow fixed (`mass_flow_scale = 1`), the cold
calibration predicts all 15 heating starts with bias -0.0159 mbar and RMSE
0.0392 mbar. The model mean local channel-velocity range is
0.3404--1.4701 m/s at heating/cooling start and 0.3395--2.3111 m/s at the end
of the histories. The increase at hot conditions follows directly from the
decrease in $\rho(T)$ at conserved mass flow.

The common cold linear closure does **not** reproduce every hot condition:

```text
comparison                  bias (model-data)   RMSE
15 heating final points       -0.4507 mbar      0.5233 mbar
3 cooling initial points      -0.4218 mbar      0.4811 mbar
3 cooling final points        -0.0336 mbar      0.0397 mbar
```

The matching cold endpoints and underpredicted hot endpoints show that a
temperature/velocity-dependent path contribution remains unresolved. Likely
candidates are developing-flow and contraction/expansion losses, the short
water-jacket path represented only by the cold multiplier, and uncertainty in
the modeled gas temperature used by the hydraulic diagnostic. The residual
does not justify an experiment-specific bypass fraction. Nor does it yet
justify changing the independently corrected MFC flow: in laminar flow, a
global flow scale and the hydraulic multiplier are strongly confounded.

The reserved single global quadratic minor-loss coefficient should be tested
next against complete time histories, using separate calibration and
validation cases and a thermally calibrated v9 state. It must not be fitted
together freely with DP1 offset, hydraulic scale and mass-flow scale.

Initial pre-full-test verdict: **implementation PASS; cold hydraulic
constraint PASS; hot-path hydraulic validation FAIL; thermal-coefficient
validation remains at the v8 FAIL status**. This is superseded by the
train/validation result below.

---

## 2026-07-29 - Full 2D_v9 train/validation test

### Leakage-free test design

The full test fixes the corrected MFC mass-flow multiplier at 1.0 and retains
the cold DP1 offset and resistance obtained independently from the raw
$0\le t\le10$ s starts. Thermal/optical parameters are fitted on nine heating
runs:

```text
E67, E69, E71, E72, E74, E76, E77, E79, E81
```

The following data are not used by the thermal objective:

```text
held-out heating: E68, E70, E73, E75, E78, E80
cooling validation: C69, C80, C81
```

Each irradiance group therefore contributes three training and two held-out
heating runs. Cooling remains validation-only and starts from its measured
temperature field.

The fitted thermal/optical vector is:

```text
A_Nu       = 0.00123739
B_Re       = 1.44       fixed
chi_r      = 0.180685
chi_z      = 0.200000   active lower bound
scale_456  = 1.050672
scale_304  = 1.051748
scale_256  = 0.621298
beta_rad   = 102.894 1/m
beta_opt   = 21.2914 1/m
```

The seed objective decreased from 14.4361 to 10.3625, but the derivative-free
fit ended at `MaxIters` after 120 evaluations.

### Hot-excess DP1 closure

The cold hydraulic multiplier already contains the reference-temperature
channel and path loss. Adding an ordinary $K\rho u^2/2$ term would count that
cold contribution twice. V9 therefore tests one shared incremental term:

$$
DP1_{\mathrm{model}} =
p_{0,\mathrm{DP1}}+
C_{\mathrm{hyd}}\Delta p_{\mathrm{square}}+
K_{\mathrm{hot}}\left[
\overline{\frac{\rho u^2}{2}}-
\overline{\frac{\rho_{\mathrm{ref}}u_{\mathrm{ref}}^2}{2}}
\right],
$$

where the overbar is the existing mass-flow-weighted ring aggregation and
$u_{\mathrm{ref}}$ uses the same conserved ring mass flow at 294.25 K. Only
$K_{\mathrm{hot}}$ is fitted, analytically, to complete DP1 histories from the
nine heating-training cases. Offset, cold hydraulic scale and MFC mass-flow
scale remain fixed.

The result is

```text
K_hot = 118.479
```

This number is an effective coefficient referenced to the small channel-end
dynamic pressure. It aggregates contractions, expansions, developing flow,
the water-jacket path and thermal-state error; it must not be quoted as the
minor-loss coefficient of one identified fitting.

Transient DP1 results transfer to the held-out cases:

| Phase | Mean base RMSE (mbar) | Mean augmented RMSE (mbar) | Augmented mean bias (mbar) |
| :--- | ---: | ---: | ---: |
| Nine heating training runs | 0.4398 | 0.1708 | +0.0044 |
| Six held-out heating runs | 0.4106 | 0.1932 | +0.0113 |
| Three cooling validation runs | 0.1648 | 0.0958 | -0.0840 |

The DP1 improvement on unseen data supports a common hot-path contribution
and does not require a per-run bypass. It is a provisional predictive closure
because the sensor accuracy/zero specification and the detailed tap-to-face
path geometry are still unavailable.

### Thermal and spatial validation

| Phase | Mean sensor RMSE (K) | Steady MAE (K) | t90 MAE (s) |
| :--- | ---: | ---: | ---: |
| Heating training | 76.21 | 78.72 | 700.65 |
| Held-out heating | 63.92 | 68.46 | 746.43 |
| Cooling validation | 34.71 | 15.75 | 535.71 |

The apparent improvement over the v8 pilot is not sufficient for coefficient
acceptance:

- mid-depth perimeter-minus-core sign: 0/15 correct;
- deep perimeter-minus-core sign: 0/15 correct;
- flow-slope signs: 20/21 correct, with the 456 kW/m2 T3 slope incorrect;
- `chi_z` is exactly at its lower bound;
- full-transient local sensitivity rank is 7/8;
- sensitivity condition number is 21461 and the weakest normalized column is
  radial conductivity.

The spatial failure is decisive. Measured perimeter-minus-core offsets are
roughly +20 to +53 K, whereas v9 predicts approximately -0.2 to -2.7 K. The
centered source/radial-flow architecture cannot reproduce the observation by
adjusting the declared coefficients.

### Numerical gates and artifacts

The fitted E73 mesh comparison gives:

```text
coarse -> nominal maximum sensor change = 5.8760 K
nominal -> fine maximum sensor change   = 4.9088 K
nominal -> fine DP1 change              = 0.00943 mbar
```

The fitted-point energy-rate residual is $1.42\times10^{-14}$ W. All 87 v9
execution/invariant/mesh/diagnostic assertions pass. The scientific
acceptance gates above are reported separately and are not converted into
software-test failures.

Primary artifacts are:

```text
run_2D_v9_full.jl
plot_2D_v9_full.jl
summaries/2D_v9/full_test_summary_2D_v9.txt
summaries/2D_v9/thermal_metrics_2D_v9.csv
summaries/2D_v9/dp1_transient_metrics_2D_v9.csv
summaries/2D_v9/transient_predictions_2D_v9.csv
summaries/2D_v9/steady_ordering_2D_v9.csv
summaries/2D_v9/flow_slopes_2D_v9.csv
summaries/2D_v9/identifiability_summary_2D_v9.txt
summaries/2D_v9/fitted_mesh_summary_2D_v9.txt
summaries/2D_v9/acceptance_status_2D_v9.txt
```

Post-processing produces 51 verified PNG figures under
`summaries/2D_v9/plots/`: two parity plots, 18 four-panel temperature/DP1
transient plots, 15 final axial profiles, 15 final 2D temperature fields and
one identifiability-correlation plot. The complete list and counts are recorded
in `summaries/2D_v9/plots/plot_manifest_2D_v9.txt`.

Final full-test verdict: **v9 implementation and MFC/$\rho(T)$ velocity
calculation PASS; cold hydraulics PASS; common hot DP1 prediction PROVISIONAL
PASS; thermal coefficient extraction FAIL**.

---

## 2026-07-29 - Planned 2D_v10 flow-sensitive front solid-to-gas exchange

> **Pre-implementation record:** This section defines the physical hypothesis,
> equations and acceptance checks before `2D_v10` is coded. The centered radial
> irradiation pattern is retained. No side-weighted or fitted rim irradiation
> term is introduced.

### Confirmed side-thermocouple installation

The perimeter thermocouples were installed at the middle of a square receiver
side wall after making a shallow dip to retain each junction. They are
therefore side-surface/near-surface SiC measurements rather than corner
measurements or freely exposed beads. In the area-equivalent axisymmetric
model they remain boundary/perimeter observations. The physical square-side
distance from the center is 9.5 mm, whereas the area-equivalent circular
boundary is 10.7196 mm; this mapping uncertainty will be checked separately
but cannot explain the current result because v9 predicts almost no radial
temperature difference.

### Experimental evidence for a front-flow mechanism

At fixed irradiance, the measured mid-depth perimeter temperature relative to
the 5 mm front perimeter temperature is strongly flow-sensitive:

| Irradiance | Slope of $T_{12}-T_8$ | Correlation with flow |
| :--- | ---: | ---: |
| 456 kW/m2 | +17.34 K/(L/min) | 0.960 |
| 304 kW/m2 | +10.84 K/(L/min) | 0.896 |
| 256 kW/m2 | +6.90 K/(L/min) | 0.889 |

The deep perimeter-minus-core offset also increases with flow within each
irradiance group:

```text
456 kW/m2: +2.18 K/(L/min)
304 kW/m2: +1.80 K/(L/min)
256 kW/m2: +0.97 K/(L/min)
```

By contrast, the 58 mm perimeter-minus-core offset is nearly independent of
flow within each group. This supports testing a flow-mediated front cooling
and gas-preheating mechanism before altering the incident radial flux.

### Defect in the v9 front boundary

V9 applies natural convection from every front solid/felt/casing cell directly
to ambient:

$$
\dot Q_{\mathrm{front,conv}}=
h_{\mathrm{natural}}A_{\mathrm{front}}
(T_s-T_{\mathrm{amb}}),
$$

where `front_nusselt_correlation` has no MFC-flow dependence and enforces
$h\ge10$ W/m2/K. This energy is treated as an external loss. The channel gas,
however, starts at the measured inlet temperature and receives heat only from
the internal channel-wall perimeter. Heat transfer from the exposed SiC front
web into the air being drawn into the channels is absent.

### Proposed v10 front correlation

For receiver ring $i$, define a local inlet Reynolds number using the same
conserved ring mass flow and local actual velocity as v9:

$$
\mathrm{Re}_{f,i}
=\frac{\rho_{\mathrm{in}}u_{\mathrm{in},i}D_h}
{\mu_{\mathrm{in}}}.
$$

The front forced-convection closure is

$$
\mathrm{Nu}_{f,i}
=C_f\,\mathrm{Re}_{f,i}^{m_f}\mathrm{Pr}_{\mathrm{in}}^{1/3},
\qquad
h_{f,i}=\frac{\mathrm{Nu}_{f,i}k_{\mathrm{in}}}{D_h}.
$$

The initial sensitivity study fixes $m_f=0.5$ and varies only the common
prefactor $C_f$. The exponent is not fitted simultaneously with $C_f$.

### Conservative gas preheating

The exposed frontal SiC web area in each ring is
$A_{\mathrm{front,solid},i}$. A bounded effectiveness formulation is used:

$$
\mathrm{NTU}_{f,i}
=\frac{h_{f,i}A_{\mathrm{front,solid},i}}
{\dot m_i c_{p,\mathrm{in}}},
\qquad
\epsilon_{f,i}=1-\exp(-\mathrm{NTU}_{f,i}),
$$

$$
T_{g,i}(0^+)=
T_{\mathrm{in}}+
\epsilon_{f,i}\left[T_{s,i}(z_1)-T_{\mathrm{in}}\right],
$$

$$
\dot Q_{f\rightarrow g,i}
=\dot m_i c_{p,\mathrm{in}}
\left[T_{g,i}(0^+)-T_{\mathrm{in}}\right].
$$

The same $\dot Q_{f\rightarrow g,i}$ is subtracted from the front receiver
solid cell and added to the channel gas through its updated inlet
temperature. It is a loss from the solid but not from the complete receiver;
the global energy ledger must therefore remain closed.

### External front losses and double-counting rule

For the SiC receiver opening:

- front thermal radiation to the surroundings remains an external loss;
- natural-convection loss is replaced by the forced front-to-inlet-gas
  transfer;
- no independent flow-dependent ambient-loss coefficient is fitted.

Natural convection remains on the exposed front felt/rim/casing regions that
do not feed the receiver channels. A captured/escaping-air fraction is not
introduced in the first implementation. This prevents the same forced
convective heat from being counted both as ambient loss and gas gain.

The internal channel-wall exchange begins downstream of the frontal web and
retains the apparent Macro-ECM Nusselt law. Frontal area and internal wetted
perimeter are distinct surfaces, so this is not geometrical double counting.

### Required pre-calibration sensitivity test

Before fitting any v10 coefficient, sweep one common $C_f$ over a range that
includes zero and report, for all heating cases:

```text
T12 - T8       front-to-mid axial response
T12 - T9       58 mm perimeter-core response
T11 - T10      107 mm perimeter-core response
T3             outlet-gas response
Q_front_to_gas useful front exchange
front radiation and remaining external convection
```

The hypothesis advances to calibration only if increasing $C_f$:

1. makes $T_{12}-T_8$ more positive, with stronger response at higher flow;
2. moves both radial offsets toward the measured positive sign;
3. increases gas enthalpy by exactly the heat removed from the front solid;
4. leaves the MFC mass-flow and DP1 hydraulic closures unchanged;
5. preserves zero-flow behavior, uniform equilibrium, global energy
   conservation and mesh stability.

If the front mechanism corrects the axial profile but does not materially
change the radial signs, the remaining radial discrepancy will be retained as
a separate structural problem rather than forcing it through irradiation or
an unmeasured bypass.

Current v10 status at the time this section was written: **hypothesis
documented before implementation; code and sensitivity results pending**.

## 2026-07-29 - 2D_v10 implementation and pre-fit sensitivity result

### Implementation

The documented front-exchange hypothesis was implemented in `2D_v10.jl`
without changing the centered irradiance distribution. Two fields were added
to `HeatTransferParameters2D`:

```text
front_coefficient          C_f, default 0
front_reynolds_exponent    m_f, default 0.5
```

For every receiver ring and every ODE evaluation, the code now:

1. uses the conserved v9 ring mass flow to calculate the inlet Reynolds
   number;
2. evaluates the documented front Nusselt and heat-transfer coefficient;
3. calculates the bounded effectiveness and post-front gas temperature;
4. subtracts `front_gas_heat_transfer_W` from the first receiver solid cell;
5. begins the internal channel-wall gas march at that post-front gas
   temperature.

The heat is therefore transferred from one modeled subsystem to another. It
is not classified as an external receiver loss. Receiver-front radiation is
still external. Natural convection remains only on the front felt/casing
cells that do not feed the channels. The parameter pack remains the same 12
v9 thermal parameters; $C_f$ was deliberately excluded from calibration.

The simulation result now exposes:

```text
front_heat_transfer_coefficient[nr_rec, time]
front_gas_heat_transfer_W[nr_rec, time]
gas_temperature[:, 1, :] = gas temperature immediately after the front web
```

### Conservation and limiting-case tests

`test/smoke_2D_v10.jl` passes **21/21** tests. The checks include:

- $C_f=0$ gives zero front exchange;
- finite flow and a hot front give positive heat transfer and bounded gas
  preheating;
- reconstructed $\dot m_i c_p(T_{g,i}(0^+)-T_\mathrm{in})$ equals the reported
  ring heat rate;
- the energy ledger includes the equal front solid loss and gas gain and has
  residual below $10^{-8}$ W;
- uniform solid/gas/ambient temperature remains an equilibrium;
- zero flow gives zero front $h$, zero front heat, zero velocity, and zero
  channel pressure drop;
- the standard mass flow is unchanged;
- pack/unpack preserves the non-fitted front parameters.

### Pre-fit sweep definition

The sensitivity runner is `run_2D_v10_sensitivity.jl`. It holds fixed:

- all fitted v9 thermal parameters;
- the v9 cold hydraulic calibration and hot-excess coefficient;
- the centered v9 irradiation distribution;
- $m_f=0.5$.

All 15 heating experiments were simulated on the established calibration
mesh for:

```text
C_f = 0, 0.05, 0.10, 0.25, 0.50, 1, 2, 4, 8
```

The upper values are intentionally extreme and probe saturation rather than
representing plausible fitted values.

### Main sensitivity result

The front exchange moves the mean axial offset in the required direction but
does not reproduce its flow dependence or either radial profile:

| $C_f$ | Mean model $T_{12}-T_8$ | Axial-offset RMSE | Mean $h_f$ | Mean gas preheat | Mean front-to-gas heat |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | -128.1 K | 143.1 K | 0 W/m2/K | 0 K | 0 W |
| 1 | -77.7 K | 94.2 K | 124.8 W/m2/K | 50.8 K | 9.45 W |
| 2 | -38.1 K | 60.4 K | 249.7 W/m2/K | 89.3 K | 16.62 W |
| 4 | +19.0 K | 48.3 K | 499.3 W/m2/K | 142.7 K | 26.55 W |
| 8 | +83.9 K | 99.8 K | 998.7 W/m2/K | 200.4 K | 37.39 W |

The experimental mean $T_{12}-T_8$ is **+10.88 K**. Thus an average sign
change can be forced near $C_f=4$, but this requires approximately
500 W/m2/K and 143 K mean inlet-gas preheat. It does not give a good
case-by-case profile match.

More decisively, increasing $C_f$ weakens rather than strengthens the modeled
flow slope:

| Irradiance | Measured slope | Model slope, $C_f=0$ | Model slope, $C_f=4$ | Model slope, $C_f=8$ |
| --- | ---: | ---: | ---: | ---: |
| 456 kW/m2 | +17.34 | +5.39 | +1.83 | -1.45 |
| 304 kW/m2 | +10.48 | +7.10 | +1.11 | -1.97 |
| 256 kW/m2 | +6.90 | +8.78 | +0.92 | -1.78 |

Units are K/(L/min). At the strongest exchange the slope reverses, contrary
to the observed positive correlations.

The mean modeled radial offsets remain almost invariant over the entire
sweep:

```text
T12 - T9:  -1.23 K at C_f=0; -1.25 K at C_f=8
            measured mean = +24.48 K

T11 - T10: -0.70 K at C_f=0; -0.80 K at C_f=8
            measured mean = +35.97 K
```

The mean T3 absolute error improves from 78.4 K at $C_f=0$ to a minimum near
65.8 K at $C_f=2$, then worsens to 69.1 K at $C_f=4$ and 76.0 K at
$C_f=8$. The exchange increases the modeled hot DP1 through the predicted gas
temperature/density change, but the conserved mass flow remains exactly
unchanged. This is a thermal consequence, not a refit of the hydraulic
closure.

### Acceptance decision

The fixed-$m_f=0.5$ front-exchange hypothesis **fails the pre-calibration
gate**:

- conservative energy transfer: **pass**;
- zero-flow and equilibrium limits: **pass**;
- axial mean moves in the correct direction: **pass**;
- observed axial flow sensitivity: **fail**;
- positive mid and deep radial offsets: **fail**;
- plausible simultaneous T3/profile improvement: **fail**.

Accordingly, $C_f$ must **not** be added to the calibration vector on the
evidence above. Doing so would trade a forced average axial sign for worse
flow ordering while leaving the radial structural error untouched.

Artifacts:

```text
summaries/2D_v10/front_sensitivity_cases_2D_v10.csv
summaries/2D_v10/front_sensitivity_slopes_2D_v10.csv
summaries/2D_v10/front_sensitivity_summary_2D_v10.txt
summaries/2D_v10/plots/front_sensitivity_profile_offsets_2D_v10.png
summaries/2D_v10/plots/front_sensitivity_flow_slopes_2D_v10.png
summaries/2D_v10/plots/front_sensitivity_exchange_2D_v10.png
```

## 2026-07-29 - Planned v10 reference-normalized flow-exponent test

The fixed $m_f=0.5$ test showed that increasing the mean front-exchange
strength cools the front but weakens the modeled $T_{12}-T_8$ flow slope. The
next diagnostic therefore tests whether a more flow-sensitive front
coefficient can correct that failure. It is a sensitivity test, not a fit.

Directly comparing the raw coefficient $C_f$ at different exponents is
invalid because changing $m_f$ also changes the coefficient magnitude at
every Reynolds number. Define a representative Reynolds number
$\mathrm{Re}_{f,\mathrm{ref}}$ from the v9 conserved ring-flow distribution at
10 standard L/min and the MFC reference temperature. For an equivalent
strength $C_{f,\mathrm{eq}}$, use

$$
C_f(m_f)
=
C_{f,\mathrm{eq}}\,
\mathrm{Re}_{f,\mathrm{ref}}^{0.5-m_f}.
$$

Then

$$
C_f(m_f)\,\mathrm{Re}_{f,\mathrm{ref}}^{m_f}
=
C_{f,\mathrm{eq}}\,
\mathrm{Re}_{f,\mathrm{ref}}^{0.5},
$$

so all exponents have the same front Nusselt number and $h_f$ at the
reference condition. Differences across experiments are caused by the flow
dependence, not by an arbitrary change in reference magnitude.

Planned grid:

```text
equivalent strengths C_f,eq = 2 and 4
exponents m_f              = 0.5, 0.8, 1.0, 1.25, 1.5
heating cases              = all E67-E81
mesh                       = established v9 calibration mesh
```

The $m_f=0.5$ rows reproduce the preceding sensitivity. Values above one are
included as deliberate stress tests; they are not assumed physically
plausible merely because they are simulated.

Every combination will report:

```text
case-level T12-T8, T12-T9, T11-T10, T3
flow slope of T12-T8 within each irradiance group
front h, heat rate and gas preheat
mass flow and DP1 prediction
```

The more-sensitive law advances only if it improves the three irradiance
group flow slopes without requiring an extreme exponent or sacrificing the
case-level axial error and T3 response. The radial offsets remain an
independent acceptance check: exponent changes must not be credited with
capturing them unless their signs and magnitudes actually move.

Current status at the time this section was written: **test design documented
before implementation; exponent sweep pending**.

## 2026-07-29 - V10 reference-normalized exponent sensitivity result

### Execution

`run_2D_v10_exponent_sensitivity.jl` executed the documented matrix:

```text
equivalent strengths = 2, 4
front Re exponents   = 0.5, 0.8, 1.0, 1.25, 1.5
heating cases        = E67-E81
reference flow       = 10 standard L/min
reference Re         = 73.4110
```

The raw $C_f$ was rescaled for every exponent so the front Nusselt number was
identical at $\mathrm{Re}_{f,\mathrm{ref}}$. All v9 thermal, hydraulic and
centered optical parameters remained fixed. The output tables include final
observed DP1 as well as the spatial temperatures.

### Effect of exponent at the profile-matching strength

For $C_{f,\mathrm{eq}}=4$:

| $m_f$ | Mean $T_{12}-T_8$ | Axial RMSE | Flow-slope RMSE | T3 MAE | Mean $h_f$ | DP1 RMSE |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.50 | +18.96 K | 48.31 K | 11.01 K/(L/min) | 69.10 K | 499 W/m2/K | 0.289 mbar |
| 0.80 | +15.76 K | 40.91 K | 8.40 K/(L/min) | 66.30 K | 504 W/m2/K | 0.279 mbar |
| 1.00 | +13.82 K | 37.17 K | 6.72 K/(L/min) | 64.48 K | 510 W/m2/K | 0.274 mbar |
| 1.25 | +11.62 K | 34.18 K | 4.71 K/(L/min) | 62.29 K | 523 W/m2/K | 0.268 mbar |
| 1.50 | +9.70 K | 33.20 K | 2.93 K/(L/min) | 60.23 K | 540 W/m2/K | 0.263 mbar |

The experimental mean axial difference is +10.88 K. Increasing $m_f$
therefore improves the mean axial offset, case-level axial RMSE, T3 MAE and
the flow-slope error simultaneously within the tested matrix.

At $C_{f,\mathrm{eq}}=4$, $m_f=1.5$, the modeled slopes are:

| Irradiance | Model | Experiment | Model correlation |
| --- | ---: | ---: | ---: |
| 456 kW/m2 | +12.67 | +17.34 | 0.974 |
| 304 kW/m2 | +8.70 | +10.48 | 0.916 |
| 256 kW/m2 | +7.84 | +6.90 | 0.927 |

Units are K/(L/min). This confirms that a sufficiently superlinear front law
can restore much of the observed flow response.

### Why this is not yet an acceptable coefficient

The numerical improvement is real, but it is not yet a defensible extracted
front heat-transfer correlation:

1. The best point is the upper stress-test exponent, $m_f=1.5$. The range
   through $m_f=1.0$ remains well below the measured slopes, especially at
   456 and 304 kW/m2.
2. Although the reference strength is normalized, the $m_f=1.5$ law produces
   case-level mean $h_f$ from 156 to 1275 W/m2/K across the experiments.
   Mean gas preheat spans 51.6 to 208.8 K and front-to-gas heat spans 9.38 to
   44.94 W.
3. The remaining axial errors are irradiance-structured rather than random.
   At the best point, mean $T_{12}-T_8$ bias is approximately:

   ```text
   456 kW/m2: +34.7 K
   304 kW/m2:  +3.4 K
   256 kW/m2: -41.6 K
   ```

   A common front-flow law therefore does not by itself explain the
   irradiance dependence.
4. The radial discrepancy is unchanged:

   ```text
   model mean T12-T9  = -1.27 K; measured = +24.48 K
   model mean T11-T10 = -0.77 K; measured = +35.97 K
   ```

5. The best axial point has final DP1 RMSE 0.263 mbar and positive bias
   0.128 mbar. The lower-strength $C_{f,\mathrm{eq}}=2$, $m_f=1.5$ point has
   better DP1 RMSE, 0.183 mbar, but its mean axial difference remains
   -43.47 K with 58.31 K axial RMSE. No tested point simultaneously gives
   the stronger axial correction and the better DP1 behavior.

The conserved mass flow is unchanged in every run; the DP1 changes arise
only from the calculated gas temperature, density and velocity.

### Decision

The exponent sensitivity gives a **directional success but parameter-
extraction failure**:

- reference-normalized experiment design: **pass**;
- stronger flow dependence improves the measured axial trend: **pass**;
- adequate result with $m_f\le1$: **fail**;
- plausible front $h_f$ and gas-preheat range: **not established**;
- common behavior across irradiance groups: **fail**;
- radial profile recovery: **fail**;
- simultaneous axial and DP1 improvement: **fail**.

Consequently, neither $C_f$ nor $m_f$ should be added to the calibration
vector yet. The result should instead be interpreted as evidence that the
missing flow-linked mechanism is stronger than a classical
$\mathrm{Re}^{0.5}$ front closure and may represent more than ordinary
front-web convection. Treating the fitted $m_f=1.5$ value as a fundamental
heat-transfer exponent would over-interpret the data.

Artifacts:

```text
summaries/2D_v10/front_exponent_cases_2D_v10.csv
summaries/2D_v10/front_exponent_slopes_2D_v10.csv
summaries/2D_v10/front_exponent_summary_2D_v10.csv
summaries/2D_v10/front_exponent_summary_2D_v10.txt
summaries/2D_v10/plots/front_exponent_flow_slopes_2D_v10.png
summaries/2D_v10/plots/front_exponent_profile_tradeoff_2D_v10.png
summaries/2D_v10/plots/front_exponent_radial_exchange_2D_v10.png
```

Post-run regression checks:

```text
test/smoke_2D_v10.jl  21/21 pass
test/smoke_2D_v9.jl   34/34 pass
```

## 2026-07-29 - Ceramic receiver literature audit and v11 strategy

Full assessment:

```text
summaries/2D_v11_literature_strategy.md
```

### Literature-supported findings

The local literature set was reviewed with emphasis on ordered ceramic
honeycombs and monolith channels rather than transferring foam correlations
to the present geometry.

1. Capuano and Fend's review of the Becker ceramic square-channel receiver
   uses a local Kays/Graetz entry law with a finite downstream Nusselt
   asymptote. It treats exposed-web exchange separately and uses the web
   thickness as that problem's characteristic length.
2. Fend's detailed SiC channel model places most wall-to-air transfer within
   approximately the first 20--30 mm. This is a distributed channel-entry
   layer, not an instantaneous front-face gas-temperature jump.
3. Cornejo and Hayes show that local properties, inlet contraction and
   heating rate materially affect monolith entry transfer. They require a
   correlation to converge to a fully developed value and warn that one
   universal Graetz curve can fail under strong variable-property heating.
4. Fend and Hoffschmidt experimentally show that a parallel-channel carrier
   can retain temperature nonuniformity that a foam rapidly compensates.
   Avila-Marin and Kribus describe thermally induced flow maldistribution:
   hotter, more viscous paths can receive less flow, reinforcing hot zones.
5. For honeycombs the channels prevent radial mass exchange. The appropriate
   macroscopic radial-flow question is therefore how the conserved inlet flow
   divides among parallel channels, not how much flow bypasses the receiver
   or crosses channel walls.
6. The Avila-Marin review reports porous-inlet convective losses as normally
   small relative to thermal losses. This weighs against interpreting the
   v10 high-exponent exposed-front sensitivity as a physical loss or
   heat-transfer coefficient.

### Quantitative conflict with the current internal closure

V9/v10 uses

$$
\mathrm{Nu}
=
A_{\mathrm{Nu}}\mathrm{Re}^{1.44}\mathrm{Pr}^{1/3}
\left(D_h/z\right)^{1/3}
=
A_{\mathrm{Nu}}\mathrm{Re}^{1.107}\mathrm{Gz}^{1/3}.
$$

With `minimum_nusselt=0`, this law approaches zero downstream and has a
strong independent Reynolds dependence in addition to Graetz scaling. Using
the fitted `A_Nu=0.00123739`, the experimental flow range and the verified
5, 58 and 107 mm depths gives approximately:

| Condition | Current Nu | Receiver Kays-form Nu |
| --- | ---: | ---: |
| 4.53 L/min, 5 mm | 0.114 | 4.335 |
| 4.53 L/min, 58 mm | 0.050 | 3.722 |
| 4.53 L/min, 107 mm | 0.041 | 3.694 |
| 18.34 L/min, 5 mm | 0.856 | 6.047 |
| 18.34 L/min, 58 mm | 0.378 | 3.905 |
| 18.34 L/min, 107 mm | 0.308 | 3.795 |

The Kays values are not assumed correct merely because they are published.
The discrepancy demonstrates that the fitted v9/v10 number is an apparent
whole-model compensator and cannot yet be called a channel heat-transfer
coefficient.

The Capuano front-web correlation is not adopted. Its stated lower range is
approximately `Re=100`, whereas the present web-thickness Reynolds range is
only about 9--36. The v10 front sensitivity also used the channel hydraulic
diameter rather than the web thickness.

### Approved v11 sequence

The next model-form test is:

1. Freeze the corrected MFC total mass flow, geometry, optics, solid
   properties, external losses and the t0-calibrated hydraulic basis.
2. Disable the v10 exposed-front gas jump.
3. Replace the apparent internal law by a local, temperature-dependent
   square-channel Graetz law. Test fixed downstream alternatives
   `Nu_fd=3.61`, `2.98`, and the Capuano `3.66` form; do not fit a new
   exponent or virtual origin.
4. After the axial kernel is tested, replace the prescribed radial factor
   `1-c_r(r/R)^2` by a common-pressure parallel-ring solve:

   $$
   \Delta p_i(\dot m_i,T_i)=\Delta p_{\mathrm{common}},
   \qquad
   \sum_i\dot m_i=\dot m_{\mathrm{MFC}}.
   $$

   Every ring remains axially isolated; no bypass and no radial mass transfer
   are introduced.
5. DP1 becomes the common predicted path pressure plus its established
   sensor offset, not an average of inconsistent ring pressure drops. Cold
   t0 data test the hydraulic base; hot DP1 tests the calculated
   temperature-dependent density, viscosity and velocity.
6. Permit at most one order-unity heat-transfer scale only if a fixed
   literature form captures profile shapes with a common magnitude error.
   Do not refit heat-transfer exponent, entry length, optics and
   conductivity together.

### Predeclared decision gates

V11 must conserve mass and energy, converge with first-cell refinement,
approach a finite downstream Nusselt value, improve the v10 axial-profile
benchmark (`T12-T8` RMSE below 33.2 K), preserve the observed flow-slope
directions, and not degrade held-out DP1 beyond the v9 reference RMSE of
about 0.193 mbar without an explicit pressure-model improvement.

The equal-pressure stage must also improve the positive measured radial
offsets `T12-T9` and `T11-T10`; otherwise thermal maldistribution is rejected
as their explanation. Side thermocouples remain mapped to the middle of the
side wall at the verified 5, 58 and 107 mm depths. Their shallow mounting dips
do not justify a side-weighted irradiance pattern.

Decision: implement v11 as a no-refit model-form matrix first. A failure of
the literature closures triggers a single-channel surrogate or source/sensor
audit, not another free front-Reynolds exponent.

### Geometry clarification: open outer channels with rear groove obstruction

The observed outer flow restriction is not a set of blocked monolith
channels. All outer channels remain open. A subsequent gauge of the rear
support geometry found an approximately 13 mm diameter region that is
completely free; discharge paths outside that diameter may face partial
obstruction from the support groove.

This supplies a physical radial-resistance seed for the common-pressure v11
model:

$$
\Delta p_{\mathrm{groove},i}
=
I_{\mathrm{groove},i}K_{\mathrm{groove}}
\frac{\rho_{i,\mathrm{out}}u_{\mathrm{groove},i}^2}{2}.
$$

Consequences for implementation:

1. Retain all 100 channels and their full internal heat-transfer perimeter.
2. Apply the additional quadratic loss only to the groove-overlapped outer
   rings.
3. Set the nominal free radius to 6.5 mm. For the corrected 361 mm2 receiver
   face, the free circular area is 132.73 mm2 or 36.77% of the face. The
   remaining groove-exposed region is 63.23%. This supersedes the earlier
   outermost-row/2.14 mm nominal mapping.
4. The 63.23% figure is an overlap area, not a closed-channel fraction.
   Retain flow and heat-transfer area throughout it and apply only a partial
   quadratic restriction.
5. Calculate each axisymmetric ring's groove weight from its exact area
   overlap outside 6.5 mm. Use 12 and 14 mm free diameters as geometry
   sensitivities around the gauged value.
6. Keep the corrected MFC total mass flow conserved. The groove redistributes
   flow from the perimeter toward the core; it does not reduce receiver flow
   unless an independently demonstrated external leak or MFC control limit
   exists.
7. Use cold t0 DP1-versus-flow data to identify or bound one combined groove
   strength. Because the loss is quadratic, the core/edge flow contrast
   should weaken as total flow decreases, although temperature-dependent
   resistance may partially oppose it.
8. Disable the old empirical heating-only hot-excess loss when testing the
   groove law to avoid double counting. The cold-identified groove model must
   predict the hot DP1 response through local $\rho(T)$ and velocity.

This observed groove makes the perimeter-hot mechanism substantially more
specific:

$$
\text{rear groove restriction}
\rightarrow
\dot m_{\mathrm{edge}}\downarrow
\rightarrow
T_{\mathrm{edge}}\uparrow
\rightarrow
(\mu/\rho)_{\mathrm{edge}}\uparrow
\rightarrow
\dot m_{\mathrm{edge}}\downarrow.
$$

It is therefore incorporated as measured geometry plus a pressure-loss law,
not as a fitted radial irradiance pattern, blocked-channel fraction or
unmeasured bypass.

## 2026-07-29: v11 implementation and no-refit model-form results

### What was implemented

`2D_v11.jl` implements the predeclared literature test without refitting the
v9 optical, solid or external-loss parameters:

1. The v10 exposed-front gas-temperature jump is disabled.
2. The internal apparent closure can be replaced by the local Kays/Graetz
   form

   $$
   \mathrm{Gz}=\mathrm{Re}\,\mathrm{Pr}\frac{D_h}{z},\qquad
   \mathrm{Nu}
   =
   S_h\left[
   \mathrm{Nu}_{fd}
   \frac{0.104\,\mathrm{Gz}}
        {1+0.016\,\mathrm{Gz}^{0.8}}
   \right]
   \left(\frac{T_g}{T_w}\right)^n.
   $$

   Local gas temperature supplies $\mu$, $c_p$ and $\mathrm{Pr}$; local wall
   temperature supplies the conductivity used in
   $h=\mathrm{Nu}k_w/D_h$. The production branch uses
   `Nu_fd=3.61`, `n=0` and `S_h=1`. Alternatives `Nu_fd=2.98` and
   `Nu_fd=3.66, n=0.45` were also tested.
3. A parallel-ring iteration conserves the corrected MFC mass flow exactly
   while solving for a common path pressure drop.
4. The rear support is represented by an exact annular overlap outside the
   measured 13 mm free diameter:

   $$
   \Delta p_{\mathrm{groove},i}
   =
   I_iK_{\mathrm{groove}}\frac{\rho_{i,\mathrm{out}}u_{i,\mathrm{out}}^2}{2}.
   $$

   All channels remain open and retain their complete internal flow and
   heat-transfer areas. The historical heating-only excess-loss coefficient
   is zero in this branch.
5. Results now retain ring flow, ring pressure, groove pressure, overlap,
   local `Pr`, `Gz`, `Nu`, velocity and the common-pressure convergence
   residual.

Reproducible entry points are `run_2D_v11.jl`,
`plot_2D_v11.jl`, `test/smoke_2D_v11.jl` and
`test/mesh_2D_v11.jl`. Numerical tables and figures are under
`summaries/2D_v11/`.

### Cold t0 identification: the groove is not uniquely identified by DP1

Nine heating records satisfying the predeclared cold-state tolerance were
used. The two cold fits were:

| Cold hydraulic form | Offset (mbar) | Shared resistance scale | Groove K | RMSE (mbar) | AICc |
| --- | ---: | ---: | ---: | ---: | ---: |
| Equal pressure, no groove | -0.61383 | 1.94402 | 0 | 0.02541 | -60.11 |
| Equal pressure, 13 mm free diameter | -0.54284 | 0.97020 | 184.16 | 0.02452 | -55.95 |

The groove fit has a physically interesting decomposition: it changes the
shared channel-resistance scale from about 1.94 to about 0.97 and assigns the
remaining pressure loss to the measured rear restriction. It is not,
however, a statistical improvement. RMSE improves by only 0.00088 mbar and
the extra parameter worsens AICc by 4.16.

Profiling fixed groove coefficients demonstrates the tradeoff:

| Groove K | 0 | 10 | 40 | 80 | 160 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| Fitted resistance scale | 1.944 | 1.751 | 1.384 | 1.163 | 0.996 |
| Cold RMSE (mbar) | 0.02541 | 0.02492 | 0.02476 | 0.02465 | 0.02453 |

Therefore cold DP1 supports a combined path resistance but cannot separately
determine how much belongs to channel friction and how much belongs to the
rear groove. `K=184.16` is a useful model-form test value, not a validated
groove coefficient.

### Full 15-case no-refit matrix

Six variants were run over E67--E81 on the same sensitivity mesh:

- historical v9 apparent heat transfer with prescribed radial flow;
- fixed-flow Graetz with `Nu_fd=3.61`;
- fixed-flow Graetz with `Nu_fd=2.98`;
- fixed-flow Capuano alternative (`Nu_fd=3.66`, `n=0.45`);
- Graetz with equal-pressure open discharge; and
- Graetz with equal pressure plus the cold-fitted rear groove.

This is intentionally a **frozen-v9 compatibility test**, not a fair
recalibrated competition between v9 and v11. V10 also retained the v9 fit and
varied only its new front-exchange terms. In the nominal Graetz branch,
v9's fitted `A_Nu` and its fixed `B_Re` are inactive, and the old hot-excess
hydraulic coefficient is disabled. However, v11 still inherits the
v9-fitted radial and axial conductivity scales, three irradiance-group power
scales, radiative extinction and optical extinction. In particular,
`k_scale_z=0.2` was bound-active in the v9 fit. Those parameters can be
compensators for the v9 apparent heat-transfer law and cannot be presumed
optimal for Graetz heat transfer.

Key aggregate results are:

| Variant | Axial `T12-T8` RMSE (K) | Radial `T12-T9` RMSE (K) | Radial `T11-T10` RMSE (K) | T3 MAE (K) | hot DP1 RMSE (mbar) | flow-slope RMSE (K/LPM) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| v9 apparent, prescribed | 143.08 | 25.85 | 38.01 | 78.38 | 0.161 | 7.25 |
| Graetz 3.61, prescribed | 178.64 | 25.45 | 38.08 | 99.95 | 0.523 | 22.72 |
| Graetz 2.98, prescribed | 176.94 | 25.45 | 38.08 | 99.84 | 0.523 | 22.77 |
| Capuano alternative, prescribed | 177.88 | 25.45 | 38.08 | 99.66 | 0.522 | 22.87 |
| Graetz 3.61, equal-pressure open | 178.78 | 25.50 | 38.09 | 99.95 | 0.523 | 22.73 |
| Graetz 3.61, equal-pressure groove | 176.14 | 24.83 | 37.97 | 99.91 | 0.570 | 23.01 |

The v11 Graetz forms fail the axial-profile gate by a large margin. More
importantly, they reverse every observed flow trend. The measured slopes of
`T12-T8` are `+17.34`, `+10.48` and `+6.90 K/(L/min)` at 456, 304 and
256 kW/m2. The nominal Graetz/groove predictions are `-14.59`, `-9.40`
and `-6.29 K/(L/min)`.

Changing the downstream asymptote or applying the published temperature
correction has almost no practical effect. The literature-scale interfacial
exchange is already strong enough that the gas/solid coupling is near its
effective NTU saturation in this continuum model.

### The rear-groove hypothesis produces flow maldistribution but not the measured solid radial profile

The nominal groove branch predicts a mean core-to-edge mass-flux ratio of
2.18, reaching approximately 3.1 at the highest flow. The ratio decreases
toward roughly 1.5 at the lowest flow, so the quadratic restriction produces
the expected convergence of channel flow distributions as total flow falls.
Mass is conserved and ring pressure drops agree to better than
`8.7e-6` relative error.

Despite that strong hydraulic effect, the predicted solid radial offsets
remain nearly zero:

| Quantity | Measured mean scale | Groove-model mean | Model RMSE |
| --- | ---: | ---: | ---: |
| `T12-T9` at 58 mm | about +24.5 K | -0.25 K | 24.83 K |
| `T11-T10` at 107 mm | about +36 K | -0.70 K | 37.97 K |

For example, E67 has a predicted core/edge flow ratio of 2.67, but predicts
`T12-T9=+0.68 K` and `T11-T10=-0.36 K`, compared with measured
`+24.44 K` and `+53.06 K`. The current continuum solid therefore smooths
the imposed hydraulic contrast before it appears in the side-versus-core
thermocouple observables.

This rejects **flow redistribution alone within the current axisymmetric
solid model** as the explanation of the radial thermocouple profile. It does
not show that the observed support groove is absent.

Repeating the complete 15-case groove branch for 12, 13 and 14 mm free
diameters does not change this conclusion:

| Free diameter | Cold-fitted `K` | Mean core/edge mass flux | `T12-T9` mean / RMSE | `T11-T10` mean / RMSE |
| ---: | ---: | ---: | ---: | ---: |
| 12 mm | 291.68 | 3.33 | +0.19 / 24.38 K | -0.63 / 37.88 K |
| 13 mm | 184.16 | 2.18 | -0.25 / 24.83 K | -0.70 / 37.97 K |
| 14 mm | 207.18 | 2.26 | -0.25 / 24.82 K | -0.70 / 37.97 K |

Even the 12 mm case, with the largest hydraulic contrast, leaves the radial
solid signal near zero. The fitted `K` values remain non-unique decompositions
of the cold path resistance.

### Pressure result

Replacing the empirical hot-excess v9 pressure term with the cold-fitted
groove does not transfer to heating: hot DP1 RMSE rises from 0.161 mbar for
the historical comparison to 0.570 mbar, with a -0.373 mbar bias. Because
the hot pressure calculation uses the predicted gas-temperature field, its
failure is coupled to the rejected thermal model. No heating-only groove
refit is justified.

### Verification

- v9 regression: 34/34 tests passed.
- v10 regression: 21/21 tests passed.
- v11 unit/smoke suite: 27/27 tests passed.
- v11 nominal-mesh confirmation: 4/4 tests passed.
- Mass conservation is at floating-point precision.
- The groove branch pressure-equality residual is below `8.7e-6`.
- The uniform-state energy ledger residual is below `1e-8 W`.
- The implemented local law approaches its prescribed finite
  `Nu_fd` downstream.
- Re-running E67, E72 and E77 on the nominal 14/7/3 radial,
  45-node receiver and 30-node rear mesh changes any final sensor by at most
  1.26 K. Axial offsets change by at most 0.36 K and radial offsets by at
  most 0.04 K. The rejection is not a sensitivity-mesh artifact.

### Decision

V11 is retained as a tested model-form branch and **is not accepted as the
new calibrated physical model**. More precisely, the frozen-v9 parameter
point is rejected for v11; the run does not yet reject every independently
recalibrated Graetz model. No heat-transfer scale was fitted in this first
stage because the predeclared condition for a simple one-parameter scale was
correct profile shape with a common magnitude error, whereas the inherited
point has the wrong axial slope signs and essentially no radial signal.

The next stage must therefore be a v11-specific sensitivity and
identifiability audit followed, only if supported, by a constrained
training-only refit. It must use the existing nine heating training cases and
six held-out heating cases. At minimum it should examine `S_h`,
`k_scale_r`, `k_scale_z`, `beta_opt` and the three irradiance-group delivery
scales, while avoiding a simultaneous unconstrained refit of every v9
parameter. Cooling remains validation. This is needed to distinguish a
rejected transport form from rejection caused by v9 compensation parameters.

### Inheritance sensitivity: the v9 optical fit materially confounds the axial verdict

A broad one-at-a-time screen was subsequently run using the high- and
low-flow records in every irradiance group. It tested:

- `S_h=0.10, 0.25, 0.50`;
- `k_scale_r=0.05, 0.30`;
- `k_scale_z=0.50, 1.00`;
- `beta_opt=10, 50, 110 1/m`;
- `beta_rad=20, 300 1/m`; and
- a common 0.8/1.2 factor on delivered power.

The screen confirms that the user's concern is substantive. At the inherited
v9 value `beta_opt=21.2914 1/m`, the six-endpoint axial RMSE is 203.73 K
and all slopes are negative. Restoring `beta_opt=110 1/m` changes all three
slopes to positive and reduces endpoint axial RMSE to 67.57 K. None of the
other one-at-a-time inherited-parameter changes restores all three signs.

The complete 15-case confirmation at `beta_opt=110 1/m` gives:

| Metric | Frozen v9 `beta_opt=21.29` | `beta_opt=110` confirmation |
| --- | ---: | ---: |
| Axial `T12-T8` RMSE | 176.14 K | 58.17 K |
| Mean `T12-T8` | 149.87 K | 35.06 K |
| 456-kW/m2 slope | -14.59 | +0.19 K/(L/min) |
| 304-kW/m2 slope | -9.40 | +1.05 K/(L/min) |
| 256-kW/m2 slope | -6.29 | +0.52 K/(L/min) |
| T3 MAE | 99.91 K | 91.19 K |
| hot DP1 RMSE | 0.570 mbar | 0.437 mbar |
| mid radial mean / RMSE | -0.25 / 24.83 K | -1.01 / 25.61 K |
| deep radial mean / RMSE | -0.70 / 37.97 K | -0.88 / 38.19 K |

The positive full-group correlations are weak (`0.05`, `0.23`, `0.24`) and
the slopes remain far below the measured `17.34`, `10.48` and
`6.90 K/(L/min)`. The result does not yet validate v11. It does show that
the earlier statement “Graetz reverses the flow trends” was conditional on
the v9 optical extinction and cannot be used to reject Graetz before a
v11-specific refit.

The radial verdict is more robust. Every broad sensitivity, including
`k_scale_r=0.05`, leaves the radial means near zero or negative, and
`beta_opt=110` does not improve them. Thus:

- **axial/Graetz verdict:** reopened; controlled v11 recalibration is
  justified;
- **groove-only radial verdict:** still rejected within the current
  axisymmetric solid/observation model.

The next calibration must profile optical extinction jointly with the
smallest identifiable subset of `S_h`, axial conductivity and group power
scales, use only the nine established heating training cases, and report the
six held-out cases and cooling unchanged. Radial conductivity should not be
freed merely to chase the side TCs because its broad sensitivity has the
wrong direction and the current annular observation model remains suspect.

The next model should not add a side-weighted irradiance profile or another
free Reynolds exponent. The evidence instead directs the next test toward
the mapping between the three-dimensional square receiver/support/side-wall
geometry and the thermocouple observations. In particular, the present
axisymmetric perimeter state is an annular average, while the side
thermocouples were fixed in shallow wall dips and the rear groove is a local
mechanical boundary. A bounded square-sector/contact or two-solid-population
test should separate:

1. internal SiC radial conduction;
2. local receiver-to-holder/groove thermal contact and rear obstruction; and
3. the side-thermocouple observation/contact model.

The conserved MFC mass flow, local $\rho(T)$ velocity, common-pressure
hydraulics and measured 13 mm free-diameter geometry should remain. The
groove coefficient must remain profiled or externally constrained rather
than presented as independently identified.

## 2026-07-29: completed staged v11-specific calibration

### Fitting separation and method

The agreed staged calibration was completed without validation leakage:

```text
training heating:
E67 E69 E71 E72 E74 E76 E77 E79 E81

held-out heating:
E68 E70 E73 E75 E78 E80

cooling validation:
C69 C80 C81
```

Stage 1 fitted only optical extinction `beta_opt`, the Graetz multiplier
`S_h`, and axial SiC conductivity scale `k_scale_z`. Its objective used
axial differences at the verified 5, 58 and 107 mm locations, the core
58-to-107 mm difference, the receiver/T3 relation and within-irradiance flow
ordering. It did not include the structurally unresolved side-minus-core
radial offsets.

Twenty-two broad true-simulation design points were evaluated. A quadratic
response surface proposed a candidate, after which seven true local
evaluations selected the result. The response surface chose only where to
evaluate; every selection loss came from actual v11 transient integrations.

Stage 2 held the transport parameters fixed and profiled each irradiance
group's power scale separately. It used `T8`, radially averaged temperatures
at 58 and 107 mm, T3 and T2. This prevents the known side/core mismatch from
controlling delivered power.

### Selected parameters

| Parameter | v9 inherited | Staged v11 | Status |
| --- | ---: | ---: | --- |
| `beta_opt` | 21.2914 1/m | 218.8 1/m | interior |
| `S_h` | not applicable | 1.142 | interior, order unity |
| `k_scale_z` | 0.2 | 0.2 | lower-bound active |
| power scale, 456 | 1.05067 | 1.25 | interior |
| power scale, 304 | 1.05175 | 1.05175 | unchanged |
| power scale, 256 | 0.62130 | 0.77352 | interior |

The final normalized training profile and level objectives are 3.5724 and
1.4331. The six-parameter residual Jacobian has singular values

```text
31.70, 27.68, 17.77, 11.27, 5.06, 3.38
```

giving local rank 6/6, condition number 9.37 and maximum absolute column
correlation 0.454. This is a cleaner local numerical rank than v9, but does
not validate the parameters: `k_scale_z` is still bound-active and the
independent validations fail.

### What improved

The complete E67--E81 axial `T12-T8` RMSE falls from 176.14 K at the
frozen-v9 v11 point to 44.14 K. Training and held-out axial RMSE are similar,
45.47 and 42.06 K, so this particular improvement transfers.

All three full-group flow-slope signs are now correct:

| Irradiance | Model slope | Observed slope | Model correlation |
| ---: | ---: | ---: | ---: |
| 456 kW/m2 | +8.02 | +17.34 | 0.848 |
| 304 kW/m2 | +4.62 | +10.48 | 0.678 |
| 256 kW/m2 | +4.24 | +6.90 | 0.762 |

Thus the frozen-v9 negative-slope result was partly a wrong-parameter
rejection. The staged model still predicts only about half the observed
slope at 456 and 304 kW/m2 and remains above the predeclared 33.2 K axial
RMSE target.

Cooling DP1 RMSE is also good at 0.049 mbar. Mass conservation, energy
closure and equal-pressure convergence remain satisfied.

### What failed

The axial-profile improvement does not transfer to the complete thermal
response:

| Phase | v9 mean sensor RMSE | staged v11 | v9 steady MAE | staged v11 |
| --- | ---: | ---: | ---: | ---: |
| Heating training | 76.21 K | 114.39 K | 78.72 K | 98.68 K |
| Held-out heating | 63.92 K | 98.63 K | 68.46 K | 81.63 K |
| Cooling validation | 34.71 K | 70.99 K | 15.75 K | 19.94 K |

The time histories expose the primary failure. Staged v11 heats most
receiver states to near plateau within a few hundred seconds, while the
experiments continue evolving for roughly 1,000--3,000 seconds. It also
cools too quickly. Mean absolute t90 errors are approximately 1,186 s in
training, 1,280 s in held-out heating and 1,167 s in cooling. This is a
missing thermal-timescale or observation-response issue, not a small optical
power adjustment.

Radial ordering remains wrong in every heating case:

```text
mean modeled T12-T9  = -1.22 K; sign correct 0/15
mean modeled T11-T10 = -1.03 K; sign correct 0/15
```

Training/held-out radial RMSE remain about 25--26 K and 38 K. The staged
transport/optical fit therefore does not rescue the groove-only radial
mechanism.

Heating DP1 also fails to transfer:

```text
training DP1 RMSE = 0.482 mbar
held-out DP1 RMSE = 0.366 mbar
v9 held-out reference = 0.193 mbar
```

The cooling DP1 improvement cannot compensate for this heating failure.

### Numerical verification and decision

`test/check_2D_v11_staged.jl` passes 12/12 checks. Parameters are finite and
within bounds; the training profile improves from the `beta_opt=110` anchor;
the uniform-state energy residual is zero; total mass is conserved; pressure
equality remains below `1e-5`; and the largest sensitivity-to-nominal mesh
final-sensor change over E67/E72/E77 is 2.58 K.

The staged fit is **rejected for coefficient extraction**.

Two conclusions can now be separated:

1. A v11-specific optical/transport fit can restore the observed axial
   flow-slope direction. The frozen-v9 v11 result was partly a
   wrong-parameter rejection, exactly as questioned.
2. The current v11 state and observation equations still cannot
   simultaneously reproduce transient timescales, absolute held-out
   temperatures, positive radial offsets and heating DP1. Therefore
   `S_h=1.142`, `beta_opt=218.8 1/m` and the power scales must not be
   reported as validated physical coefficients.

Further fitting of `S_h`, `beta_opt`, power and conductivity inside the same
state/observation model is not justified. The next test must independently
constrain the missing time-scale/observation physics while addressing the
square side-wall/support geometry:

- a thermocouple/contact response model or independently measured sensor
  time constant;
- an audit of effective thermal mass and holder/contact storage using
  cooling data;
- a square-sector or two-solid-population representation for side/core
  temperatures and local support contact; and
- retention of conserved MFC mass flow, local `rho(T)` velocity,
  equal-pressure channel allocation and the measured groove geometry.

The formal result is in
`summaries/2D_v11/staged_acceptance_status_2D_v11.txt`.

---

## 2026-07-29: v12 assembly-capacitance and observation-model formulation

### Motivation and evidence carried forward

Staged v11 improved the axial flow-trend direction but heated and cooled about
1,200 s too quickly, predicted the wrong core/perimeter ordering in all 15
heating cases and failed held-out absolute temperatures.  Independently, the
experimental slow mode gives

```text
C_eff = 301 +/- 23 J/K
```

compared with approximately 42--47 J/K for the measured 40 g receiver.  The
large fitted SiC `Cp` and `k` values used in earlier three-dimensional trials
are therefore treated as model-form evidence for missing participating
hardware, redistribution and observation physics, not as material-property
measurements.

### Geometry and hardware audit

Annotated assembly images and subsequent clarification establish:

```text
SiC receiver                 19 x 19 mm, length 137 mm
aluminum enclosure OD        150 mm
aluminum radial thickness     18 mm
enclosure inner radius        57 mm
enclosure internal length    165 mm
rear aluminum backplate       18 mm
alumina adaptor OD            77.6 mm
alumina adaptor length        57 mm
receiver overlap              29 mm
tube overlap                  28 mm
alumina tube bore diameter    13 mm
```

The adaptor is solid apart from the square receiver opening in the overlap
section and the tube bore in the rear section.  Its calculated volume and mass
are

```text
overlap volume = 126.686 mL
rear volume    = 128.709 mL
total volume   = 255.395 mL
mass at 3900 kg/m3 = 0.99604 kg
```

The felt fills all remaining cavity space.  Receiver/adaptor contact exists
but is not firm.  The tube does not touch the felt: its internal portion is
covered by the adaptor and the remaining portion is outside the cavity inside
the continuously water-cooled flange.

The v11 hardware representation was incomplete:

1. the 0.996 kg alumina adaptor was absent and replaced by a direct
   receiver-to-tube conductance;
2. the aluminum sleeve ended at 137 rather than 165 mm;
3. the 18 mm rear backplate was absent;
4. the rear felt/adaptor topology and finite loose contact were absent; and
5. interior flow-exposed thermocouples were compared directly with solid
   temperature.

Nominal v12 inventory at representative temperature is:

| Component | Mass | Nominal capacity |
| --- | ---: | ---: |
| SiC receiver | 0.0401 kg | about 45 J/K |
| cavity felt | 0.1926 kg | about 262 J/K at the nominal Cp |
| alumina adaptor | 0.9960 kg | about 0.9--1.1 kJ/K |
| aluminum sleeve | 3.325 kg | 2.993 kJ/K |
| aluminum rear backplate | 0.852 kg | 0.767 kJ/K |
| rear alumina tube | about 0.040 kg | about 36--45 J/K |

These nominal capacities greatly exceed the observed 301 J/K slow-mode
capacity.  V12 must therefore use finite contacts so only a fraction
participates on the experimental time scale.  The capacities must not be
lumped into the SiC state.

### Felt-property correction and uncertainty

The installed felt grade is unknown because the supplier did not provide a
specification sheet.  The earlier v11 law,

```text
k_felt = 0.06 + 1.2e-10 T^3
```

gave about 0.21 W/m/K at 800 degC.  The supplied RS-3000 sheet, used only as
a shape prior, reports 0.11 W/m/K at the same temperature.  V12 replaces the
old cubic rise with interpolation through

```text
20, 500, 800, 1100, 1400, 1600 degC
0.050, 0.080, 0.110, 0.170, 0.260, 0.320 W/m/K
```

and exposes one bounded global conductivity scale.  Density remains
140 kg/m3.  A separate modest global felt-Cp scale is allowed because the
grade is unknown.  Neither may vary by experiment.

### V12 state and contact topology

`2D_v12.jl` retains the conservative v11 receiver, gas, DP1 and
equal-pressure ring equations and adds:

```text
T_adaptor(t)       physical 0.996 kg dense-alumina adaptor
T_housing,rear(t)  missing 28 mm sleeve plus 18 mm backplate
```

The adaptor couples to:

- the perimeter receiver cells over the measured rear 29 mm through a loose
  receiver/adaptor contact coefficient;
- the first 28 mm of rear tube through an adaptor/tube contact coefficient;
  and
- the rear felt through the cylindrical adaptor-to-felt resistance.

The inherited v11 direct receiver-to-tube adaptor surrogate is disabled; it
must not operate in parallel with the explicit adaptor state.

There is explicitly no tube-to-felt path.  The extra aluminum housing state
couples conductively to the existing sleeve and loses heat from its missing
side/rear external area.  The continuously cooled flange remains the terminal
tube boundary.

### Thermocouple observations

Side thermocouples remain solid-contact observations at 5, 58 and 107 mm with
one shared response time.  T9 and T10 are flow-exposed channel probes, so their
equilibrium target is

```text
T_TC,target = w_wall T_wall + (1-w_wall) T_gas
```

with one shared wall fraction and one shared response time.  T3 receives a
shared outlet-probe response time.  These are observation equations and do
not add thermal mass to the receiver energy balance.

### Pre-fit verification and staged fit rules

`test/smoke_2D_v12.jl` passes 31/31 checks covering the hardware inventory,
felt-property knots, uniform equilibrium, augmented-state response, interior
observation leverage, conserved MFC mass flow, equal-pressure convergence and
the augmented instantaneous energy balance.

The calibration order is fixed before fitting:

1. cooling only: felt `k`/`Cp`, loose receiver/adaptor contact and shared
   sensor time constants;
2. heating training only: interior wall fraction and, only after stage 1,
   optical/heat-transfer levels;
3. unchanged held-out heating and cooling validation;
4. rejection if the full response still needs experiment-specific material
   properties, violates the 301 +/- 23 J/K participation constraint or fails
   radial ordering.

The physical SiC density, `Cp(T)` and `k(T)` are not fitted per experiment.

### Completed v12 staged calibration and nominal validation

An implementation audit before the definitive run found that the first v12
draft still inherited v11's direct receiver-to-tube adaptor conductance.  That
path bypassed the newly explicit adaptor while operating in parallel with it.
It was removed, the smoke/energy tests were repeated, and every cooling,
heating and validation stage below was rerun.  Earlier provisional v12 outputs
are superseded.

Cooling used only C69/C80/C81 and selected:

| Quantity | Selected value | Qualification |
| --- | ---: | --- |
| felt conductivity scale | 1.20 | global effective grade/packing scale |
| felt heat-capacity scale | 1.50 | refinement upper bound active |
| receiver/adaptor contact `h` | 15 W/m2/K | weak/loose contact |
| shared side-probe time | 180 s | tested upper bound active |
| shared interior-probe time | 60 s | tested upper bound active |
| shared outlet-probe time | 180 s | tested upper bound active |
| interior wall fraction | 0.80 | tested upper bound active |

The contact corresponds to only about

```text
G_rec-adaptor = h A_overlap
               = 15 * (4 * 0.019 * 0.029)
               = 0.0331 W/K.
```

Thus the confirmed 0.996 kg adaptor is present but only weakly participates
on the experimental time scale.  The cooling fit assigns most of the
additional accessible inventory to the felt.  The active bounds on felt Cp
and all probe times show that cooling does not uniquely separate material
storage from observation delay.

Heating then froze all cooling quantities and the physical SiC properties.
The nine established heating-training cases selected:

```text
beta_opt = 150 1/m
Graetz heat-transfer scale S_h = 1.20
power scales = 1.25, 1.05, 0.77 for 456, 304, 256 kW/m2
```

No experiment-specific `k` or `Cp` was used.  The six held-out heating cases
and three cooling cases were evaluated on the nominal mesh without refitting.

### Definitive nominal-mesh results

| Phase | Mean sensor RMSE | Steady MAE | t90 MAE | Axial RMSE | DP1 RMSE |
| --- | ---: | ---: | ---: | ---: | ---: |
| heating training | 91.70 K | 89.69 K | 756 s | 53.52 K | 0.459 mbar |
| held-out heating | 90.21 K | 86.15 K | 848 s | 51.80 K | 0.349 mbar |
| cooling validation | 31.17 K | 30.06 K | 700 s | 9.88 K | 0.058 mbar |

Compared with prior models:

- cooling RMSE improves from 70.99 K in staged v11 and 34.71 K in v9 to
  31.17 K;
- held-out heating improves from 98.63 K in staged v11 but remains worse than
  63.92 K in v9;
- transient timing improves markedly relative to staged v11's roughly
  1,200--1,300 s errors, but remains unacceptable; and
- heating DP1 remains worse than v9's 0.193 mbar held-out reference.

The radial verdict does not improve:

```text
T12-T9 sign correct  = 0/15
T11-T10 sign correct = 0/15
mid radial RMSE      = 29.2--29.9 K
deep radial RMSE     = 46.7--46.8 K
```

The flow-exposed observation model cannot create the missing profile because
the calculated local gas and wall temperatures differ by only a few kelvin.
The heating profiles remain almost radially uniform and much flatter/cooler
than the measured square-receiver profiles.

The very coarse screening mesh changes representative final sensors by
5.94--16.04 K relative to the nominal mesh.  This is within the explicit
20 K diagnostic gate but too large to support precise coefficient extraction.
The nominal validation results themselves remain the decision basis.

### V12 decision

V12 is **rejected for coefficient extraction** but retained as a useful
mechanism result.

What is supported:

1. the installed adaptor is approximately 0.996 kg under the confirmed solid
   77.6 mm geometry;
2. it is only weakly coupled through the loose receiver contact;
3. correcting the felt/inventory topology substantially improves cooling; and
4. the large artificial SiC Cp used previously was partly substituting for
   accessible felt/housing storage and probe response.

What is not supported:

1. the selected felt `Cp` scale is validated, because it and the probe delays
   are bound-active and confounded;
2. the v12 optical and Graetz values are validated transport coefficients;
3. an axisymmetric one-solid receiver can explain the side/core profile; or
4. further fitting within this topology is justified.

The next defensible model-form test is a square-sector or explicit two-solid
population receiver.  It must distinguish gas-exposed channel-wall material
from side/corner/support-connected material while retaining the v12 hardware
inventory and loose adaptor contact.  No side-weighted irradiance source or
experiment-specific material property should be introduced.

Definitive artifacts are under `summaries/2D_v12/`, including all 18 temporal
plots, 15 axial plots, parity, full transient/final CSVs, calibration traces,
mesh transfer and the formal rejection file.

## V13: conservative active-channel and side/support populations

### Motivation and pre-fit decision

V12 showed that observation lag and missing hardware inventory improve
cooling but cannot produce the measured side-to-interior temperature
difference.  The calculated channel gas and wall temperatures are too close,
and an axisymmetric single receiver temperature field forces the shallow-dip
side thermocouples to observe the same local population that transfers heat
to the flowing gas.

V13 therefore tests a model-form change, not another irradiance fit.  The
measured SiC inventory is partitioned into two co-located macroscopic
populations:

1. an active channel-wall population that receives most of the gas heat
   transfer; and
2. a side/corner/support population observed by T8, T12 and T11, with reduced
   gas participation, its own felt contact and the loose rear-adaptor support.

The split is conservative:

```text
C_active = (1-f_side) C_SiC
C_side   = f_side C_SiC
C_active + C_side = C_SiC
```

No SiC mass or heat capacity is added.  The same radial/axial optical source
is divided in direct proportion to the two solid fractions; there is no
side-weighted illumination.  The populations exchange

```text
dQ_active-side = G'_as dz (T_active,mean - T_side).
```

The side population has a perimeter felt path

```text
dQ_side-felt = h_side-felt (4 w_rec dz)
               (T_side - T_felt),
```

and it owns the measured loose receiver/adaptor contact over the final 29 mm.
The active population retains the channel flow, equal-pressure ring
allocation and Graetz heat-transfer equations.  A single global
`f_gas,side <= f_side` allocates the fraction of receiver-to-gas heat transfer
borne by the side population; the remainder is removed from the active
population.  This explicitly tests reduced flow-accessible area without
inventing or fitting an unmeasured bypass mass flow.

Front radiation, receiver/felt exchange, active/side exchange and adaptor
exchange are all split conservatively.  `test/smoke_2D_v13.jl` verifies the
capacity inventory, uniform equilibrium, heated response, conserved MFC mass
flow, equal-pressure convergence and instantaneous whole-assembly energy
closure.

### Frozen and fitted quantities

The definitive v12 hardware, felt, observation, optical, Graetz, power and
hydraulic values are frozen for the first v13 test.  Only four shared
population quantities may change:

```text
f_side       side/support fraction of measured SiC inventory
f_gas,side   absolute fraction of gas heat transfer assigned to that population
G'_as        active-to-side conductance per receiver length
h_side-felt  side/support-to-felt contact coefficient
```

They are fitted on the established nine heating-training cases with a
cooling-retention penalty from C69/C80/C81.  The six held-out heating cases
remain untouched until nominal-mesh validation.  Experiment-specific
properties, side-weighted irradiance and a fitted bypass branch remain
prohibited.

### Completed v13 staged identification

The initial 18-point population screen and 29-point refinement selected

```text
f_side       = 0.075
f_gas,side   = 0.0225
G'_as        = 30 W/K/m
h_side-felt  = 5 W/m2/K
```

with the side fraction at the refinement lower bound and active/side
conductance at its upper bound.  A formal bound expansion to
`f_side = 0.025--0.10` and `G'_as = 30--100 W/K/m` moved the optimum farther
toward

```text
f_side       = 0.025
f_gas,side   = 0.0125
G'_as        = 100 W/K/m
h_side-felt  = 5 W/m2/K.
```

At 800 K the full receiver capacity is 49.25 J/K, so the selected side
capacity is only about 1.23 J/K.  Its total active/side conductance is

```text
G_as = G'_as L = 100 * 0.137 = 13.7 W/K,
```

giving a characteristic equilibration time of only about 0.09 s.  The
optimization therefore collapses the new field back into the v12
single-solid limit; it does not identify a physically meaningful second slow
population.

Candidate diagnostics nevertheless show why the mechanism initially looked
promising:

| Population test | Mid radial sign | Deep radial sign | Heating RMSE | Cooling RMSE |
| --- | ---: | ---: | ---: | ---: |
| collapsed selected limit | 0/9 | 0/9 | 111.5 K | 32.0 K |
| moderate split | 7/9 | 0/9 | 121.1 K | 41.6 K |
| strongly separated | 9/9 | 3/9 | 163.4 K | 86.3 K |

A separated side population can make the middle side thermocouple hotter than
the interior probe, but it simultaneously produces the wrong deep profile
and destroys overall heating/cooling levels.

To avoid rejecting the topology merely because it inherited coefficients
fitted to v12, a second heating stage re-screened every defensible population
candidate over

```text
beta_opt = 75, 150, 250, 400 1/m
S_h      = 0.75, 1.20, 1.80
```

and then refitted the three shared irradiance-group scales.  This independent
stage again selected the collapsed population and

```text
beta_opt = 150 1/m
S_h = 1.80
power scales = 1.25, 1.2075, 0.77.
```

The Graetz scale and 304-kW/m2 power scale are at the tested upper bounds.
Meaningful population splits were not rescued by the optical-depth or
heat-transfer refit.

### Definitive v13 nominal-mesh validation

The six held-out heating cases and three cooling cases were run without
further adjustment:

| Phase | Mean sensor RMSE | Steady MAE | t90 MAE | Axial RMSE | Mid radial RMSE | Deep radial RMSE | DP1 RMSE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| heating training | 88.00 K | 85.74 K | 735 s | 45.35 K | 27.50 K | 43.07 K | 0.402 mbar |
| held-out heating | 82.60 K | 78.15 K | 820 s | 45.59 K | 28.27 K | 43.17 K | 0.287 mbar |
| cooling validation | 31.04 K | 31.21 K | 688 s | 10.59 K | 4.96 K | 9.77 K | 0.062 mbar |

Relative to v12, held-out heating RMSE improves from 90.21 to 82.60 K,
held-out axial RMSE from 51.80 to 45.59 K and held-out DP1 RMSE from 0.349 to
0.287 mbar.  Cooling is essentially unchanged.  These small gains come from
the refitted heating coefficients, not an identified side population.

Both heating profile-sign tests remain

```text
T12-T9 sign correct  = 0/15
T11-T10 sign correct = 0/15.
```

The representative screen-to-nominal maximum final-sensor changes are
7.86 K for E67, 3.07 K for E72 and 2.39 K for E77, passing the 20 K
diagnostic gate and improving on v12 mesh transfer.  The rejection is
therefore not a coarse-mesh artifact.

`test/smoke_2D_v13.jl` passes 22 checks and
`test/mesh_2D_v13_staged.jl` passes 2 checks.  All 34 requested plots were
generated: one parity plot, 18 temporal plots and 15 axial plots.

### V13 decision and physical interpretation

V13 is **rejected for coefficient extraction**.

The test does not support a globally co-located, weakly gas-coupled side
thermal mass as the missing physics.  Forcing a significant split can create
the mid-depth sign, which supports the broader idea that side and core
channels experience different thermal-fluid conditions.  But one global
side fraction and one global conductance cannot reproduce both the 58 and
107 mm profiles while retaining cooling and temperature level.

The next defensible step is an explicit square-channel network or square
sector, not another lumped side population.  It should:

1. represent center, edge and corner channel groups explicitly;
2. solve their flow rates from a common pressure drop with local
   temperature-dependent density and viscosity;
3. apply the measured radial groove obstruction only at the rear exit;
4. retain lateral SiC conduction between channel groups and the confirmed
   felt/adaptor contacts; and
5. map side thermocouples to the exterior wall and interior probes to their
   actual channel groups.

That model can test the proposed positive feedback—hotter channels have
higher resistance, receive less flow and become hotter—without inventing an
unmeasured total-flow bypass.  The v13 results show that this feedback must be
spatially and axially resolved; a global gas-participation fraction is too
coarse.

Definitive v13 artifacts are under `summaries/2D_v13/`.

## V14: square-channel orbit network

### Pre-fit formulation and acceptance rules

V13 established that a global co-located side population is too coarse.  V14
therefore resolves the actual 10 by 10 square-channel topology while
retaining symmetry.  The 100 channels are reduced to 15 exact square-symmetry
orbits: channels related by horizontal, vertical and diagonal reflection
share one solid and gas solution.  Orbit multiplicities sum to 100, so the
reduction changes neither flow area nor receiver inventory.

Each orbit has:

- one SiC channel-cell temperature per axial cell;
- one isolated axial gas stream;
- its exact channel multiplicity;
- its square-grid radius and rear-groove overlap;
- the number of exterior receiver faces it owns; and
- lateral conductive connections to neighboring channel orbits.

The measured SiC solid area

```text
A_SiC = w_receiver^2 - 100 b_channel^2
```

is divided equally among the 100 channel cells and then multiplied by each
orbit multiplicity.  Thus the total SiC mass and heat capacity remain exactly
the photographed receiver values.

#### Common-pressure channel flow

Every physical channel in orbit `g` carries the same per-channel mass flow.
At every thermal residual evaluation, v14 solves

```text
Delta p_g(mdot_g, Tgas_g(z)) = Delta p_common
sum_g multiplicity_g mdot_g = mdot_MFC.
```

The local gas density, viscosity, heat capacity and conductivity use the
calculated gas temperature.  The channel friction and Graetz heat-transfer
laws are inherited from v11/v12.  The rear groove adds its cold-DP-calibrated
quadratic loss according to the actual overlap of each square channel with
the region outside the measured 6.5 mm free radius.  It does not close or
remove outer channels.

This is the first 2D version that can represent the proposed feedback at the
individual channel-group level:

```text
hot channel -> lower density/higher resistance
            -> less flow at common pressure drop
            -> less convective cooling -> hotter channel.
```

The corrected MFC total flow remains conserved and no bypass is fitted.

#### Solid network and observations

Adjacent square-channel cells exchange heat through the inherited effective
transverse SiC conductivity.  Each orbit also conducts axially.  Exterior
channel groups exchange with the felt and the rear alumina support according
to their exact number of exterior square faces.  The v12 felt, adaptor, tube,
rear housing and water-cooled flange states remain explicit.

The centered Gaussian/Beer-Lambert source already used by v12 is evaluated at
the channel-orbit centroids and normalized over the physical square face.
There is no side-weighted source.

The thermocouple mapping is fixed before fitting:

```text
T8, T12, T11 -> midpoint exterior-side orbit at z=5, 58, 107 mm
T9, T10      -> central-channel orbit at z=58, 107 mm,
                with the established wall/gas observation blend
T3           -> mixed rear gas at 3 mm
T2           -> felt field at z=58 mm
```

The midpoint side orbit, central orbit and corner orbit are separately
reported.  This mapping uses the verified axial positions and the confirmed
shallow-dip placement in the middle of the side wall.

#### Staged test

V14 first runs with the definitive v12 hardware, observation, optics, Graetz,
power and cold-DP parameters.  Before calibration it must pass:

1. exact orbit multiplicity, flow area and SiC inventory closure;
2. whole-assembly instantaneous energy closure;
3. uniform equilibrium;
4. conserved total MFC mass and common channel pressure drop;
5. decreasing flow in a deliberately hotter channel group; and
6. stronger groove loss outside the measured free diameter.

Only after those checks will training fit a bounded transverse-conduction
scale and, if the profile mechanism is present, independently re-evaluate
optical depth, the single Graetz multiplier and shared power-group scales.
The six held-out heating cases remain untouched until the nominal-mesh
confirmation.  Experiment-specific properties and a bypass branch remain
prohibited.

### V14 numerical verification

`test/smoke_2D_v14.jl` passes 35 checks.  The implemented network has:

```text
15 symmetry orbits
sum of orbit multiplicities = 100 channels
sum of exterior square faces = 40
receiver mass = 0.0400588 kg
```

The hot-channel diagnostic deliberately raises the midpoint-side orbit by
300 K.  Its common-pressure solution reduces that orbit's flow while
preserving total MFC mass.  Uniform equilibrium, equal pressure and
whole-assembly energy-rate closure also pass.

### Required square-network cold hydraulic recalibration

The first provisional v14 runs transferred v11's annular groove coefficient
directly.  That is not valid because the square channel overlaps differ from
annular overlap fractions.  Those provisional thermal results are
superseded.

Using only the 15 heating-experiment `t0` DP1/flow points, with the established
physical sensor offset frozen, the square network selected:

```text
hydraulic resistance scale = 0.80
groove loss coefficient K  = 335
DP1 zero offset             = -0.542845 mbar
t0 DP1 RMSE                 = 0.02347 mbar
```

The best nested no-groove model gives 0.03359 mbar RMSE with resistance scale
1.661.  The groove term therefore improves cold DP1, although the two
resistance contributions remain partially correlated.  All v14 thermal
screens and validation were repeated after freezing the square-network
values.

### Completed staged v14 calibration

The first thermal topology screen varied only:

```text
lateral conductivity multiplier = 0.05, 0.20, 1.00
edge/felt contact multiplier     = 0.30, 1.00, 3.00
centered beam sigma              = 14, 30, 100 mm
```

Increasing beam sigma broadens the centered Gaussian toward uniform
illumination; it never biases power toward the side.

Weak lateral coupling can produce the observed direction at 58 mm, but the
tradeoff is severe:

| Network candidate | Mid sign | Deep sign | Heating RMSE | Cooling RMSE |
| --- | ---: | ---: | ---: | ---: |
| inherited/homogenized | 0/9 | 0/9 | 114.6 K | 33.4 K |
| moderate isolation | 7/9 | 0/9 | 136.9 K | 42.0 K |
| strong isolation | 7/9 | 4/9 | 142.3 K | 47.1 K |

The independent heating screen then tested every error-optimal and
profile-producing network over `beta_opt=50--400 1/m` and
`S_h=0.75--1.80`, followed by the three shared power-group scales.  It
selected:

```text
lateral conductivity multiplier = 1.00
edge/felt contact multiplier     = 3.00
centered beam sigma              = 100 mm
beta_opt                         = 100 1/m
Graetz scale S_h                 = 1.80
power scales                     = 1.25, 1.2075, 0.8855
```

The lateral, edge/felt, beam-width and Graetz quantities sit at tested bounds,
as do the 304 and 256 kW/m2 power refinements.  The fit therefore prefers a
nearly uniform beam and a strongly homogenized receiver, not the isolated
channel configuration that generated some radial signs.

### Definitive v14 nominal-mesh validation

| Phase | Sensor RMSE | Steady MAE | t90 MAE | Axial RMSE | Mid radial RMSE | Deep radial RMSE | DP1 RMSE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| heating training | 83.10 K | 81.40 K | 630 s | 46.38 K | 28.85 K | 46.56 K | 0.413 mbar |
| held-out heating | 81.57 K | 77.80 K | 717 s | 44.28 K | 29.70 K | 46.73 K | 0.280 mbar |
| cooling validation | 32.73 K | 36.87 K | 607 s | 12.63 K | 5.47 K | 11.50 K | 0.064 mbar |

The square network modestly improves held-out overall RMSE relative to v12
and v13, but cooling and radial errors degrade relative to v13.  Both
heating radial signs remain wrong in all 15 cases:

```text
T12-T9 sign correct  = 0/15
T11-T10 sign correct = 0/15.
```

The hydraulic effect itself is large and has the expected flow dependence.
The final per-channel core/side flow ratio decreases from approximately
4.26 at 18.7 L/min to 1.81--1.95 at 4.53 L/min.  Thus the quadratic groove
does make the profiles hydraulically converge as total flow decreases.
Nevertheless, the modeled `T12-T9` remains negative, approximately
-1.8 to -8.6 K, while every experiment is positive by 20--27 K.  The final
channel maps retain a core that is about 2--7 K hotter than the side/corner
channels.

The explicit interior observation diagnostic also fails to rescue the
profile.  At 58 and 107 mm the local central-channel gas is generally
1--14 K hotter than its local wall after being heated upstream.  Reducing the
wall fraction would therefore make the modeled interior thermocouple hotter,
not colder.  The wall fractions required to match the experiments lie
outside the physical interval `[0,1]` in nearly every case.

Representative screen-to-nominal final-sensor changes are 6.94 K for E67,
3.85 K for E72 and 2.98 K for E77.  The mesh diagnostic passes its 20 K gate.
All 49 nonzero plots were generated: parity, 18 transients, 15 axial profiles
and 15 square channel/flow maps.

### V14 decision

V14 is **rejected for coefficient extraction**, but it gives a strong
mechanism conclusion.

The proposed common-pressure feedback is present and quantitatively large:
outer channels receive much less flow, and the contrast weakens at lower
total flow exactly as expected for the quadratic groove.  It is not,
however, sufficient to reverse the predicted solid temperature ordering
while also matching cooling and absolute temperatures.  Thermohydraulic
maldistribution is therefore real but is not the missing radial-temperature
mechanism by itself.

The next model-form change should resolve the **exterior receiver wall
actually touched by T8/T12/T11** separately from the orbit-averaged edge
channel wall.  A perimeter-skin state must use physical outer-web mass,
axial SiC conduction, conduction to the adjacent edge-channel cells,
felt/adaptor contact and the same centered optical source.  It must not add
mass or side-weighted irradiance.  Before fitting such a model, the exact
transverse channel locations of T9/T10 and the approximate depth of the
side-wall dips would materially constrain the observation mapping.

Definitive artifacts are under `summaries/2D_v14/`.

## V15 pre-fit plan: explicit exterior SiC perimeter wall

### Motivation

V14 demonstrated that the measured 13 mm free rear opening produces a real,
large channel-flow redistribution and the expected reduction of that
redistribution at lower total flow.  It nevertheless predicted the
orbit-averaged outer channel wall colder than the interior measurements in
all 15 heating experiments.  Increasing transverse isolation could reverse
some profiles only by degrading the absolute temperatures and cooling
transients.

The side thermocouples were installed in shallow dips at the middle of the
physical exterior side wall.  They therefore need not measure the
orbit-averaged temperature of the complete outer channel solid assigned to
the `(9,1)` symmetry orbit.  V15 tests the narrower hypothesis that these
thermocouples follow a distinct exterior SiC perimeter-wall state.

### Frozen physics

V15 inherits the definitive v14 model without refitting:

- the 10 by 10 square-channel network and its 15 exact symmetry orbits;
- the common-pressure-drop channel-flow solution and temperature-dependent
  gas properties;
- the cold-flow hydraulic calibration
  (`hydraulic_resistance_scale = 0.80`,
  `groove_loss_coefficient = 335`, and
  `DP1_offset = -0.54284475 mbar`);
- the centered, nearly uniform irradiance distribution
  (`beam_sigma = 100 mm`) and Beer-Lambert axial deposition;
- the Graetz heat-transfer law, optical extinction, group power scales,
  felt, solid adaptor, alumina tube, aluminum enclosure, and cooling flange;
- the verified thermocouple axial locations of 5, 58, and 107 mm.

No side-weighted optical source, bypass branch, artificial experiment-wise
heat capacity, or experiment-wise conductivity is introduced.

### New state and strict mass partition

One exterior-perimeter SiC temperature, `Tskin(z,t)`, is added per axial
cell.  Its cross-sectional area is the area of a square shell of effective
thickness `delta_skin`,

`A_skin = W^2 - (W - 2 delta_skin)^2`.

That area is allocated among the exterior channel orbits in proportion to
their physical boundary-face counts.  The same area is subtracted from the
corresponding orbit solid areas.  Consequently,

`A_channel,residual + A_skin = A_SiC,v14`

at every axial location.  V15 must reproduce the v14 SiC mass and heat
capacity to numerical precision; the new state is a partition, not added
thermal mass.

The v14 direct exterior-channel/felt and exterior-channel/adaptor exchanges
are removed.  V15 instead uses

`channel solid <-> exterior SiC skin <-> felt`

and, over the measured rear overlap,

`exterior SiC skin <-> solid alumina adaptor`.

The skin has its own axial SiC conduction.  The removed share of the v14
axial conduction, local solar deposition, and front radiative loss is
transferred consistently to the skin state.  Internal channel/gas exchange
remains on the residual channel solid because the exterior wall is not an
additional gas passage.

The only new fitted structural quantities are:

1. the effective exterior-skin thickness; and
2. a shared channel-to-skin spreading-conductance multiplier.

The multiplier represents the unresolved two-dimensional spreading distance
between the orbit-average channel solid and the shallow-dip exterior-wall
location.  It is global, not experiment-specific.

### Observation model

`T8`, `T12`, and `T11` observe the filtered skin temperature at 5, 58, and
107 mm.  `T9` and `T10` retain the v14 central-channel wall/gas observation
model.  `T3`, `T2`, and all hardware states are unchanged.

### Pre-declared acceptance tests

Before inspecting a fitted result, v15 is accepted for coefficient
extraction only if all of the following hold:

1. receiver mass and heat capacity equal v14 to numerical precision;
2. uniform-temperature equilibrium and whole-assembly energy closure pass;
3. hydraulic flow splits and DP1 predictions are unchanged from v14 for the
   same channel-temperature field;
4. at least 12 of 15 heating cases have the correct sign for both
   `T12 - T9` and `T11 - T10`;
5. held-out heating mid- and deep-radial RMSE are each below 20 K;
6. held-out mean sensor RMSE does not exceed the v14 value of 81.57 K;
7. cooling mean sensor RMSE does not exceed 37.73 K, five kelvin above v14;
8. screen-to-nominal maximum final-sensor change remains below 20 K.

If a profile improvement requires the skin thickness or spreading
conductance to remain at an implausible search boundary, v15 will be treated
as a diagnostic model rather than evidence for a transferable heat-transfer
coefficient.

### V15 implementation and numerical verification

`2D_v15.jl` implements the documented mass-conserving exterior-wall
partition as a wrapper around the v14 channel-network equations.  At every
ODE evaluation it:

1. recovers the v14 energy rates;
2. removes the former direct outer-channel/felt and
   outer-channel/adaptor transfers;
3. partitions the local SiC source, axial conductance, front radiating area,
   mass, and heat capacity between residual channel solid and skin;
4. applies channel/skin, skin/felt, and skin/adaptor exchanges; and
5. advances the exterior skin as one additional state per axial cell.

This construction retains the v14 lateral channel topology and gas exchange.
The hydraulic calculation sees the same channel temperatures and equations
as v14; adding the skin does not alter the cold-flow calibration.

`test/smoke_2D_v15.jl` passes 28 checks:

- the 15 orbits, 100 channels, and 40 exterior faces are retained;
- the residual-channel and skin areas sum exactly to the v14 SiC area;
- receiver mass and common-temperature heat capacity equal v14 to
  `1e-14` relative tolerance;
- hydraulic flow splits are bitwise identical for the same channel field;
- uniform equilibrium includes the new skin and all seven observations;
- heated states are finite and conserve total mass flow; and
- whole-assembly energy closure is below `1e-7 W`.

The selected screen-to-nominal mesh test also passes.  Maximum final-sensor
changes were 6.929 K for E67, 3.850 K for E72, and 2.979 K for E77, all below
the pre-declared 20 K limit.

### V15 structural screen

The screen froze every definitive v14 parameter and evaluated:

- `delta_skin = 0.05, 0.15, 0.30 mm`; and
- channel/skin spreading scale `= 0.01, 0.03, 0.10, 0.30, 1.00`.

The first stage used E67, E72, and E77, one experiment from each irradiance
group.  None of the 15 configurations gave the correct sign at either 58 or
107 mm.  The four best objective values were then evaluated on all nine
heating-training cases and the three cooling experiments.

The numerical optimum was:

| Quantity | Selected value | Search position |
|---|---:|---|
| effective skin thickness | 0.05 mm | lower bound |
| channel/skin spreading scale | 1.00 | upper bound |
| training mid sign | 0/9 | failed |
| training deep sign | 0/9 | failed |

Both selections collapse the new state toward v14: minimum mass reduces its
dynamic independence and maximum coupling forces it toward the outer-channel
temperature.  The fit therefore does not support a separately evolving
exterior SiC wall under the inherited felt contact.

### Definitive nominal-mesh validation

| Metric | Heating training | Heating held out | Cooling |
|---|---:|---:|---:|
| mean sensor RMSE (K) | 83.03 | 81.50 | 32.93 |
| steady MAE (K) | 81.30 | 77.90 | 36.99 |
| t90 MAE (s) | 630.8 | 713.1 | 602.4 |
| axial RMSE (K) | 50.16 | 47.59 | 12.50 |
| mid radial RMSE (K) | 30.19 | 31.00 | 5.57 |
| deep radial RMSE (K) | 49.34 | 49.49 | 12.85 |
| DP1 RMSE (mbar) | 0.413 | 0.280 | 0.063 |

The heating sign tests remain:

- `T12 - T9`: 0/15 correct;
- `T11 - T10`: 0/15 correct.

Compared with v14, held-out mean sensor RMSE changes only from 81.57 to
81.50 K, while mid-radial RMSE worsens from 29.70 to 31.00 K and deep-radial
RMSE worsens from 46.73 to 49.49 K.  Cooling RMSE changes from 32.73 to
32.93 K.  These differences are small and consistent with the selected skin
state collapsing back toward v14.

For all heating experiments the modeled radial differences remain negative:
`T12 - T9` ranges approximately from -2.2 to -11.1 K and
`T11 - T10` from -5.6 to -19.7 K, whereas the measurements are positive.
The E72 axial plot makes the cause visible: the skin and outer-channel curves
nearly coincide at about 525--580 K while the measured side thermocouples are
about 740--800 K.

### V15 decision and physical interpretation

V15 is **rejected for coefficient extraction**.  The result is stronger than
simply finding a poor parameter set:

1. With the centered source, the skin receives only its mass-proportional
   share of local SiC heating.
2. The skin has no independent internal gas-heated surface.
3. It can lose heat to the felt and adaptor.
4. Weakening channel/skin coupling therefore makes the skin cooler, not
   hotter; strengthening it returns the v14 orbit temperature.

Consequently, a mass-conserving exterior wall alone cannot explain the hot
side measurements while the inherited skin/felt exchange is retained.  This
also rules out adding an arbitrary exterior thermal mass as the solution.

The next physically identifiable test should not alter the centered optical
field.  It should recalibrate the thermal path that directly changed when
the skin was introduced:

1. use C69/C80/C81 first to constrain a global skin/felt contact resistance,
   felt conductivity scale, and felt heat-capacity scale from T8/T12/T11 and
   T2 during zero-solar cooling;
2. freeze those cooling-derived values;
3. rerun the heating sign test with the v15 mass partition; and
4. only then consider a shared heating-coefficient refit.

This is physically motivated by the observed full but non-firm felt contact.
It also prevents heating errors from forcing the edge/felt multiplier to the
v14 upper bound of 3.0.  If a cooling-constrained weak contact still cannot
produce the signs, the remaining explanation is more likely an observation
or local-deposition issue than a missing bulk SiC thermal state.

Definitive artifacts are under `summaries/2D_v15/`; 49 nonempty plots include
one parity plot, 18 temporal plots, 15 axial profiles, and 15 channel maps.

## V16 pre-fit plan: cooling-first skin/felt identification

### Why this stage is necessary

V15 inherited the v14 exterior/felt contact multiplier of 3.0 even though
v14 fitted that multiplier while the side thermocouples still observed an
orbit-averaged channel state.  Once v15 introduced an explicit exterior
wall, that inherited coefficient became the direct loss path from the newly
observed state and was no longer guaranteed to be transferable.

The physical observation is that the alumina-silicate felt fills the cavity
but its contact with the receiver is not firm.  The cooling data provide the
only solar-independent way to identify this path.  V16 therefore fits it
before looking again at heating profiles.

The v15 cooling residuals also prevent assuming in advance that contact must
be weaker.  At the ends of C69/C80/C81, the model is generally 12--66 K
hotter than the measured SiC and gas temperatures, whereas the final T2 felt
errors are only about 0.1--4.3 K.  Reducing contact alone would isolate the
receiver and make its late cooling still slower.  A joint fit is required to
separate contact, felt conduction, and felt storage.

### Frozen model

V16 retains:

- the v14 square-channel hydraulics and cold DP1 calibration;
- the v15 mass-conserving exterior SiC skin;
- `delta_skin = 0.05 mm` and channel/skin spreading scale `= 1.0`;
- the centered optical field and every heating, Graetz, loss, power, and
  hardware parameter;
- the measured sensor locations and observation time constants.

No heating data enter the cooling fit.  No experiment-specific coefficient
is allowed.

### Cooling parameters and bounds

Only three global quantities may vary:

| Parameter | Screen values |
|---|---|
| skin/felt contact multiplier | 0.10, 0.30, 1.0, 3.0, 10.0 |
| felt conductivity multiplier | 0.60, 1.20, 2.40 |
| felt heat-capacity multiplier | 0.75, 1.50, 3.00 |

The contact multiplier acts only on the new skin/felt interface.  It does not
alter the felt/adaptor path or retrospectively alter the v14 channel
equations.  Conductivity and capacity scales remain shared by the entire
felt inventory.

C69 and C80 are the calibration cases.  C81 is held out during selection.
The primary cooling sensors are T8, T12, T11, and T2.  T9, T10, and T3 enter
only as a lower-weight guard against transferring error elsewhere in the
assembly.

### Pre-declared acceptance

The cooling identification is considered useful only if:

1. all v15 mass, equilibrium, hydraulic, and energy tests remain valid;
2. the selected screen configuration is finite and screen-to-nominal changes
   remain below 20 K;
3. mean primary-sensor cooling RMSE over all three cases improves by at least
   10% from the v15 value of 26.87 K;
4. held-out C81 primary-sensor RMSE does not exceed its v15 value of
   approximately 28.14 K;
5. mean T2 cooling RMSE remains below 20 K;
6. one-step neighboring parameters do not cause a qualitative collapse,
   and any boundary optimum is reported as non-identifiable.

After selecting from C69/C80, the cooling-derived parameters will be frozen
and all 15 heating experiments will be simulated.  The overall model is
accepted for coefficient extraction only if the earlier spatial criteria
also hold: at least 12/15 correct signs at both depths, held-out radial RMSE
below 20 K, held-out mean sensor RMSE no worse than 81.57 K, and cooling mean
sensor RMSE no worse than 37.73 K.

### V16 implementation and cooling screen

V16 retains the v15 equations but separates
`skin.felt_contact_scale` from the historical v14
`network.edge_felt_contact_scale`.  The old v14 transfer is still removed
exactly before the new interface is applied, so the new multiplier changes
only the skin/felt path.

`test/smoke_2D_v16.jl` passes nine checks:

- the v15 receiver mass and heat capacity remain conserved;
- uniform equilibrium is unchanged;
- changing the skin/felt multiplier changes skin cooling without changing
  the inherited network multiplier; and
- whole-assembly energy closure remains below `1e-7 W`.

The initial 45-point screen fit C69 and C80 and held C81 out.  Its best
screen-mesh point was:

| Quantity | Initial selected value |
|---|---:|
| skin/felt contact multiplier | 1.0 |
| felt conductivity multiplier | 2.4 |
| felt heat-capacity multiplier | 1.5 |
| C69/C80 primary RMSE | 23.45 K |
| held-out C81 primary RMSE | 23.06 K |

This was a meaningful transfer from calibration to C81, but the
conductivity multiplier was at the original upper bound.

### Boundary extension and nominal-mesh correction

The conductivity range was extended to 3.6, 4.8, and 7.2 while testing
neighboring contact and heat-capacity values.  The objective continued to
improve toward weaker contact and larger conductivity.  The extended screen
selected contact 0.3, conductivity scale 7.2, and capacity scale 1.5, again
at the conductivity ceiling.

The original screen-selected configuration also failed the pre-declared
mesh-transfer rule: C69 changed by 24.61 K.  Therefore the ten best and
boundary-diagnostic configurations were re-ranked directly on the nominal
mesh.  The nominal mesh again selected:

| Quantity | Definitive nominal selection |
|---|---:|
| skin/felt contact multiplier | 0.30 |
| felt conductivity multiplier | 7.20 |
| felt heat-capacity multiplier | 1.50 |
| C69/C80 primary RMSE at 61 points | 25.76 K |
| held-out C81 primary RMSE at 61 points | 25.50 K |

This is not an identifiable felt property.  A multiplier of 7.2 corresponds
to approximately 0.36 W/m/K near room temperature and 0.58 W/m/K at 500 C,
far above the temperature-shaped felt prior.  The optimum remains at the
extended ceiling, indicating that it represents a missing parallel cooling
path or compensating model error.

The definitive screen-to-nominal comparison also fails:

| Case | Maximum final-sensor change |
|---|---:|
| C69 | 28.68 K |
| C81 | 11.60 K |
| E72 | 36.61 K after the power refit |

Thus the v16 cooling coefficients cannot be reported as mesh-stable
macroscopic properties.

### Frozen-parameter nominal validation

With the nominally selected cooling parameters frozen, the 121-point
all-case results are:

| Metric | Heating training | Heating held out | Cooling |
|---|---:|---:|---:|
| mean sensor RMSE (K) | 113.19 | 124.91 | 33.59 |
| steady MAE (K) | 123.78 | 138.96 | 27.92 |
| t90 MAE (s) | 754.6 | 826.2 | 892.9 |
| axial RMSE (K) | 81.29 | 83.91 | 9.25 |
| mid-radial RMSE (K) | 36.66 | 37.63 | 5.06 |
| deep-radial RMSE (K) | 49.61 | 49.87 | 12.60 |
| DP1 RMSE (mbar) | 0.523 | 0.424 | 0.052 |

The primary cooling sensors T8/T12/T11/T2 have an aggregate RMSE of
23.12 K, compared with 26.87 K in v15, so the targeted cooling objective
improves by about 14%.  T2 improves to 9.34 K.  However, T9/T10/T3 cooling
errors increase enough that overall cooling RMSE becomes 33.59 K, slightly
worse than v15's 32.93 K.

Heating deteriorates decisively.  Held-out mean sensor RMSE rises from
81.50 K in v15 to 124.91 K.  Both radial signs remain 0/15 correct.
The high-conductivity, weak-contact boundary point tested separately on
E67/E72/E77 also retains negative radial gaps and produces sensor RMSEs of
153.1, 157.0, and 69.8 K, respectively.

### V16 decision and new spatial diagnosis

V16 is **rejected for coefficient extraction** for three independent
reasons:

1. the fitted felt conductivity runs to the extended upper boundary;
2. the selected quantities fail mesh transfer; and
3. cooling-derived parameters strongly degrade heating and do not change the
   radial signs.

The full/non-firm felt contact is real, and a lower interface multiplier is
consistent with it.  Nevertheless, the data demand an unrealistically
conductive parallel path to cool the assembly, so felt properties alone are
not the missing physics behind the heating profiles.

The axial data now give a sharper next direction.  In 10 of 15 heating
experiments the measured side temperature is highest at 58 mm; the model is
highest at 5 mm in all 15 cases.  For E72, for example, the measured
5/58/107 mm side temperatures are approximately 738/798/734 K, while the
v16 skin prediction is approximately 562/534/502 K.  Increasing radial or
felt losses cannot move an axial heat maximum from the entrance toward the
middle.

The next model should therefore test the **axial optical/volumetric source
law**, not a side-weighted radial beam:

- retain the centered radial irradiance and v14 hydraulics;
- replace the single Beer-Lambert solid source with a two-component axial
  source: a near-entrance wall-absorption component and a deeper
  channel-radiation component;
- conserve total absorbed power exactly;
- fit the shared source split and deep extinction length using the
  5/58/107 mm axial profile only;
- freeze those optical parameters before evaluating radial profiles and
  cooling.

This tests the structured-receiver volumetric effect without imposing extra
power on the measured side.  It is also directly falsifiable: the new source
must reproduce the middle-depth maximum without destroying total-temperature
parity or cooling.

Definitive artifacts are under `summaries/2D_v16/`; the 49 nonempty plots
include parity, 18 temporal comparisons, 15 axial profiles, and 15 channel
maps.

### V16 correction: phase-resolved parity and required power refit

The initial v16 parity figure can be misread because it pools three different
phases and all seven sensors on one set of axes.  Its x-axis is observed
temperature and its y-axis is modeled temperature.  The apparent
overprediction near the lower-left is primarily cooling:

- final cooling errors are positive for every sensor on average;
- heating receiver and gas sensors are strongly negative;
- heating T2 is the main positive-error exception.

For held-out heating, mean final errors range from approximately -117 K for
T10 to -238 K for T12, while T2 is about +15 K.  For cooling, mean final
errors are about +4 to +52 K.  Thus the transient heating underprediction and
the parity plot are consistent; the pooled plot contains both regimes.

More importantly, v16 retained the v14 irradiance multipliers
`(1.25, 1.2075, 0.8855)` while changing the felt conductivity and skin/felt
loss path.  That freeze was appropriate for testing whether a cooling-only
fit transferred without help, but it does not constitute a fully refitted
heating model.  The earlier v16 heating metrics must therefore be labeled
**frozen-power diagnostics**, not the final optimized v16 heating result.

Before changing the axial source shape, v16 will receive one corrective
calibration stage:

1. freeze the nominally selected cooling parameters
   `(contact, felt-k, felt-Cp) = (0.3, 7.2, 1.5)`;
2. use the nominal mesh because the screen mesh failed transfer;
3. fit one shared irradiance multiplier for each experimental irradiance
   group using only its heating-training experiments;
4. fit absolute temperature level using T8/T12/T11/T9/T10/T3, with T2 as a
   low-weight guard;
5. evaluate E68/E70, E73/E75, and E78/E80 without refitting; and
6. rerun all cooling cases unchanged.

Only the three group power multipliers may move.  The centered radial beam,
single Beer-Lambert axial law, Graetz coefficient, hydraulic network, and
all cooling-derived quantities remain frozen.  Consequently, this stage can
correct temperature level but cannot be credited with fixing an axial or
radial profile unless the corresponding held-out differences also improve.

### Completed v16 power refit

The nominal-mesh group fits and local refinements selected:

| Irradiance group | Previous scale | Refitted scale | Held-out receiver RMSE |
|---|---:|---:|---:|
| 456 kW/m2 | 1.2500 | 1.65 | 87.1 K |
| 304 kW/m2 | 1.2075 | 1.80 | 82.5 K |
| 256 kW/m2 | 0.8855 | 1.25 | 46.8 K |

The 256-group minimum was checked at 1.20, 1.25, 1.30, and 1.35; the complete
objective has an internal minimum at 1.25.  The fit therefore confirms that
the heating underprediction was partly a total-power calibration error after
the v16 loss path changed.

The corrected 121-point validation is:

| Metric | Heating training | Heating held out | Cooling |
|---|---:|---:|---:|
| mean sensor RMSE (K) | 65.62 | 61.46 | 33.59 |
| steady MAE (K) | 54.94 | 50.04 | 27.92 |
| t90 MAE (s) | 810.2 | 875.0 | 892.9 |
| axial RMSE (K) | 106.70 | 107.75 | 9.25 |
| mid-radial RMSE (K) | 42.58 | 43.72 | 5.06 |
| deep-radial RMSE (K) | 55.30 | 55.70 | 12.60 |
| DP1 RMSE (mbar) | 0.250 | 0.181 | 0.052 |

Refitting power reduces held-out mean sensor RMSE from the frozen-power
124.91 K to 61.46 K.  The updated parity plot therefore has heating points
on both sides of the 1:1 line instead of systematic heating underprediction.
Distinct markers are now used for training heating, held-out heating, and
cooling so the regimes cannot be confused.

This improvement strengthens, rather than removes, the spatial diagnosis.
The extra power magnifies the wrong axial shape:

- held-out axial RMSE increases from 83.91 to 107.75 K;
- held-out mid/deep radial RMSE become 43.72/55.70 K;
- both radial signs remain 0/15 correct; and
- the modeled side maximum remains at 5 mm in 15/15 cases.

For E72, the power-refitted model predicts side temperatures of approximately
688/648/599 K at 5/58/107 mm, compared with measured
738/798/734 K.  Total temperature level is much closer, but the model still
cools monotonically downstream instead of reaching the measured middle
maximum.  The power-refitted E72 screen-to-nominal change also rises to
36.61 K, so v16 remains mesh-unstable.

The corrected conclusion is:

- power scaling **did require refitting** and explains much of the absolute
  heating bias;
- felt-path parameters remain non-identifiable and mesh-sensitive; and
- the remaining dominant error is the axial and radial distribution of
  absorbed energy and heat removal, not insufficient total power.

V16 remains rejected for coefficient extraction, while the evidence for the
next power-conserving two-component axial source test is stronger.
