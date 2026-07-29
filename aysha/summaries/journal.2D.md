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
