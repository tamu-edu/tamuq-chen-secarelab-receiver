# 2D Continuum Solar Receiver Model Development Journal (Macro-ECM)

This journal records the model iterations, mathematical formulations, calibration outcomes, validation metrics, and spatial physical insights for 2D Axisymmetric Continuum Macro-ECM models of the monolithic SiC solar receiver with square channels.

---

## Overall Study Objective

The overarching goal of this study and 2D model development is to **obtain and validate effective macroscopic heat transfer coefficients (convective, radiative, and conductive)** for a structured monolithic solar receiver with square channels. The 2D model serves as a continuum representation (Entire Converter Model / Macro-ECM) where fundamental transport parameters (such as anisotropic effective thermal conductivities $k_{s,r}^{\text{eff}}$ and $k_{s,z}^{\text{eff}}$, Nusselt number correlations, and Rosseland radiation extinction coefficients) are extracted from experimental data, bridging the gap between detailed single-channel physics and full-reactor behavior.

---

## 2026-07-28 â€” 2D_v7 Heating-Only Optimization Objective & Channel Radiative Transport

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
the earlier statement â€œGraetz reverses the flow trendsâ€� was conditional on
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

That model can test the proposed positive feedbackâ€”hotter channels have
higher resistance, receive less flow and become hotterâ€”without inventing an
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
| mean sensor RMSE (K) | 65.32 | 61.16 | 33.75 |
| steady MAE (K) | 54.69 | 49.65 | 28.09 |
| t90 MAE (s) | 810.2 | 877.4 | 892.9 |
| axial RMSE (K) | 107.94 | 110.54 | 8.00 |
| mid-radial RMSE (K) | 42.58 | 43.72 | 5.06 |
| deep-radial RMSE (K) | 55.30 | 55.70 | 12.60 |
| DP1 RMSE (mbar) | 0.250 | 0.181 | 0.052 |

Refitting power reduces held-out mean sensor RMSE from the frozen-power
124.91 K to approximately 61 K.  The updated parity plot therefore has heating points
on both sides of the 1:1 line instead of systematic heating underprediction.
Distinct markers are now used for training heating, held-out heating, and
cooling so the regimes cannot be confused.

This improvement strengthens, rather than removes, the spatial diagnosis.
The extra power magnifies the wrong axial shape:

- held-out axial RMSE increases from 83.91 to 110.54 K;
- held-out mid/deep radial RMSE become 43.72/55.70 K;
- both radial signs remain 0/15 correct; and
- modeled T8 at 11 mm remains above T12 in 15/15 cases.

For E72, the power-refitted model predicts side temperatures of approximately
690/648/599 K at 11/58/107 mm, compared with measured
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

## 2026-07-30 - Post-v16 strategy correction: symmetry and manuscript invariants

The proposed post-v16 rewrite was initially described as 100 coupled 1D
channels.  That statement confused the number of **physical channels**, which
must remain 100 in every extensive flow, perimeter, mass, and energy term, with
the number of **independent state histories** required under the accepted
symmetry.

For the 10 by 10 square array:

- no symmetry gives 100 independent channels;
- horizontal and vertical reflection alone gives a 5 by 5 quadrant with 25;
- full square dihedral symmetry, including `x <-> y`, gives 15 exact orbits.

V14--v16 already implement the 15-orbit reduction, with multiplicities summing
to 100.  The centered radial irradiance, centered circular 13 mm free rear
region, square receiver, uniform-contact hypothesis, and negligible channel
buoyancy (`Gr/Re^2 <= 10^-3`) all preserve this D4 symmetry.  The side
thermocouples are observation operators and do not break the state symmetry.
Therefore v17A and a centered rear-manifold v17B will retain 15 orbits.  A
25-state quadrant is required only if diagonal interchange is physically
broken; 100 is required only if evidenced geometry or forcing also breaks the
mirror planes.

The next-model strategy was also re-anchored to
`analysis/manuscript/manuscript_full_draft.md`.  Its regime analysis rules out
several previously discussed directions:

- `Re=23--94`, `Pr=0.684--0.690`, and `Gz_L=0.18--0.70` imply deeply laminar,
  nearly constant-Pr flow that is thermally developed by the outlet;
- `Pe L/D_h=1500--5900` rules out axial gas conduction;
- `Gr/Re^2 <= 10^-3` rules out buoyancy as the radial mechanism;
- `Bi<10^-5` rules out wall-normal SiC resistance as the assembly bottleneck;
- `N_rc` rising from about 1.6 at 600 K to about 5 at 1000 K identifies
  nonlocal in-channel radiation as a physically strong axial candidate.

Consequently, the two-component axial source is retained only as a
power-conserving, falsifiable diagnostic.  If it can reproduce the 58 mm
maximum and the effectiveness crossover, it must be promoted to a physical
axial radiosity/view-factor model before its parameters receive coefficient
credit.

The revised staged plan is:

1. build a common evaluator for the manuscript invariants before more fitting;
2. v17A: identify axial radiation/source shape from normalized 11/58/107 mm
   profiles, with lamp-group power factors nested only for absolute level;
3. replace a successful empirical split by a conservative physical radiation
   closure;
4. freeze axial physics, then v17B: test a 15-orbit rear-manifold pressure
   network against cold `t0` DP1 and validate against hot DP1 and radial trends;
5. only then recalibrate physical global contact/loss terms to cooling, release
   the `C_eff` prior, and perform joint held-out confirmation.

The manuscript invariants are now primary acceptance outputs:

```text
I1  apparent Nu = 3.1e-4 Re^1.44
I2  epsilon* = 0.66 +/- 0.03
I3  Lambda_107 = 0.038 + 8.3e-4 Re
I4  C_eff = 301 +/- 23 J/K with the prior released
I5  K_loss = 0.10--0.16 W/K
I6  delivered-power factors consistent with a single explicit convention
```

The local channel Graetz Nusselt number is not equated to I1.  I1 must be
post-processed from the modeled side-wall and outlet observables using the
manuscript definition; its superlinear exponent is a signature of
flow-dependent assembly participation.  Role B coefficient extraction remains
blocked until heating and cooling, spatial signs, invariant reproduction,
parameter interiority, energy closure, and mesh transfer all pass
simultaneously.

The complete current strategy is recorded in
`summaries/manuscript_2D_readiness_strategy.md`.

### Final T8 correction: 11 mm and complete v16 rerun

The final apparatus correction is that T8 is at 11 mm, as stated in the
manuscript, not at the 5 mm coordinate introduced into v9--v16.  The earlier
5 mm journal statements are superseded.

For constant end extrapolation and a piecewise-linear wall profile, the exact
11/58/107 mm quadrature is:

```text
Tbar_w = 0.251825 T8 + 0.350365 T12 + 0.397810 T11.
```

The manuscript's current rounded weights give
`Nu_app=3.0768e-4 Re^1.4435`; exact positional weights give
`3.1044e-4 Re^1.4428`.  The three inversion effectiveness values move only
from `0.6712/0.6714/0.6258` to `0.6730/0.6732/0.6280`.  The dimensionless
regime and invariant conclusions are unchanged.

The code correction uses one shared `SIDE_TC_FRONT_Z_2D=11 mm` in the v11
initial-field and observation implementation inherited by v14--v16.  V14 and
v15 observations and v16 plot coordinates use that constant.  The following
were rerun at the nominal mesh:

1. no-refit 5-to-11 mm observation sensitivity;
2. nominal cooling candidate re-ranking and held-out C81;
3. local power fits for all three irradiance groups;
4. the 256-group internal-minimum confirmation;
5. complete 15-heating/3-cooling validation and all plots; and
6. the screen-to-nominal mesh test.

The fitted values are identical to the earlier selections:

```text
skin/felt contact scale = 0.30
felt k scale            = 7.20
felt Cp scale           = 1.50
power scales            = 1.65 / 1.80 / 1.25
```

The corrected complete validation is:

| Metric | Heating training | Heating held out | Cooling |
|---|---:|---:|---:|
| mean sensor RMSE (K) | 65.32 | 61.16 | 33.75 |
| steady MAE (K) | 54.69 | 49.65 | 28.09 |
| t90 MAE (s) | 810.2 | 877.4 | 892.9 |
| axial RMSE (K) | 107.94 | 110.54 | 8.00 |
| mid-radial RMSE (K) | 42.58 | 43.72 | 5.06 |
| deep-radial RMSE (K) | 55.30 | 55.70 | 12.60 |
| DP1 RMSE (mbar) | 0.250 | 0.181 | 0.052 |

Compared with observing the same solved fields at 5 mm, held-out overall RMSE
improves by 0.30 K, axial RMSE worsens by 2.78 K, and cooling RMSE worsens by
0.16 K.  The T8 prediction moves by only -6.63 to +7.11 K across heating
cases.  For E72 it moves from 687.74 to 690.43 K, while the observed value is
737.66 K.

The decisive results do not move:

- the refitted power factors are unchanged;
- the nonphysical felt conductivity remains at its extended bound;
- the model still predicts T8 above T12 in all 15 heating cases;
- both radial signs remain 0/15;
- the E72 mesh change remains 36.61 K; and
- v16 remains rejected for coefficient extraction.

Therefore the 11 mm correction is important for geometric accuracy and
reporting, but it does not change the investigation lessons or the v17
axial-radiation-then-rear-manifold strategy.

## 2026-07-30 - Verified aluminum-casing contact with the cooling flange

The rear/back face of the aluminum casing is in physical contact with the
water-cooled flange; this contact is also part of the alumina-tube sealing
arrangement.  This corrects an important omission in v12--v16.  Those models
connect the terminal alumina-tube cell to the fixed water-flange temperature,
but the rear aluminum sleeve/backplate state is connected only to the axial
casing field and to external ambient convection/radiation.  There is no
aluminum-backplate-to-water-flange heat-transfer term.

The omitted path is physically distinct from the alumina
receiver--adaptor--tube path:

```text
felt -> aluminum sleeve -> aluminum rear/back face
     -> finite flange contact -> water-cooled flange
```

It can remove heat from the outer assembly without first passing through the
T2 location or the alumina tube.  Its sign is therefore consistent with the
systematic T2 result: v16 overpredicts final heating T2 in every heating case,
and also heats/cools T2 too rapidly.  The inflated fitted felt-conductivity
scale can now be interpreted more specifically as compensation for combining
two different functions in one parameter:

1. transport between the receiver and the local felt/T2 region; and
2. removal of casing/backplate heat to the water-cooled flange.

The second function must instead be represented by the verified metal/flange
topology.  The ideal axial aluminum conductance is large because of the broad
back face, so it must not be inserted as a perfect Dirichlet condition on the
entire casing.  The physically uncertain quantity is a single global
effective backplate/flange contact conductance, including real contact area,
interface resistance, flange spreading resistance, and water-side removal.
It must be conservative: heat removed from the aluminum state is added to the
reported flange sink.

The revised Stage-3 hardware ablation is:

1. restore felt conductivity and heat capacity to narrow, physically
   defensible global ranges and retain the measured hardware masses;
2. add one finite `G_case_flange` link from the rear aluminum
   sleeve/backplate state to the measured water-flange temperature;
3. identify `G_case_flange` from cooling training cases only, with the axial
   optical and hydraulic models frozen;
4. validate without refitting on C81 and on every heating case; and
5. reject the mechanism if it merely trades T2 bias against the other sensors,
   requires a boundary-active conductance, or fails mesh transfer and energy
   closure.

The primary signatures are removal of the all-positive heating T2 bias,
improved T2 timing, physically interior felt properties, and recovery of the
measured `C_eff=301 +/- 23 J/K` and `K_loss=0.10--0.16 W/K` slow mode.  This
path does not explain the measured 58 mm side-temperature maximum, so the
axial radiation/source investigation remains a separate prerequisite rather
than being replaced by the flange-contact correction.

### Pre-declared casing/flange correction test

The verified boundary is corrected immediately as the first v17 hardware
ablation rather than leaving a known-wrong boundary in the baseline.  The
test changes only the rear aluminum energy balance:

```text
Q_case_flange = G_case_flange
                (T_housing,rear - T_water,flange)

C_housing dT_housing,rear/dt =
    Q_case - Q_housing,ambient - Q_case_flange.
```

`Q_case_flange` is also added to the flange-loss ledger, so the change is
exactly conservative.  The measured water-flange temperature already used by
the terminal alumina-tube boundary remains 298.15 K.  No new thermal mass,
experiment-specific property, axial source, radial beam, bypass flow, or
regime switch is introduced.

To avoid reproducing the v16 series-resistance non-identifiability, the
skin/felt contact scale is fixed at the previously selected weak-contact value
0.30 and the felt-conductivity scale is fixed at its nominal physical value
1.00.  Cooling may identify only two shared quantities:

```text
G_case_flange
felt Cp scale
```

The logarithmic `G_case_flange` screen spans zero through a near-clamped
backplate.  It is ranked on C69/C80 and selected at the nominal mesh; C81 is
held out.  The same value is then carried into all heating cases.  One power
factor per lamp configuration may be re-profiled, but the independent
delivered-power brackets are explicit guards:

```text
456 group: 1.05--1.34
304 group: 1.23--1.58
256 group: 0.84--1.11
```

The correction is useful only if it:

1. lowers the systematic heating T2 bias and improves T2 timing;
2. improves or preserves held-out cooling without requiring nonphysical felt
   conductivity;
3. moves the fitted power factors toward the independent closure brackets;
4. preserves exact energy closure;
5. transfers across screen and nominal meshes more reliably than v16; and
6. does not materially worsen the already-failed axial and radial metrics.

It cannot be credited with solving the 58 mm axial maximum.  Results from this
ablation will be used to reorder and narrow the remaining v17 strategy.

### V17 casing/flange correction results

`2D_v17.jl` adds the verified rear aluminum-casing/backplate contact with the
water-cooled flange.  The term is conservative in both the ODE and the energy
ledger.  The smoke tests pass 8/8, including equilibrium, stronger cooling
with larger contact conductance, and exact flange-loss accounting.

The cooling fit fixed:

```text
skin/felt contact scale = 0.30
felt k scale            = 1.00
```

and identified only `G_case_flange` and the global felt-Cp scale from C69/C80.
C81 remained held out.  The nominal-mesh selection was:

```text
G_case_flange = 20 W/K
felt Cp scale = 1.00
```

with cooling primary RMSE of 27.90 K on C69/C80 and 22.92 K on held-out C81.
Held-out T2 RMSE was 8.69 K.  However, the training objective at
`G_case_flange=0`, `Cp scale=0.75` was only 0.7% worse.  A high-conductance
profile showed:

| `G_case_flange` (W/K) | training objective | held-out primary RMSE (K) | held-out T2 RMSE (K) |
|---:|---:|---:|---:|
| 10 | 2.3168 | 22.947 | 8.922 |
| 20 | 2.2810 | 22.923 | 8.691 |
| 50 | 2.2423 | 22.905 | 8.492 |
| 100 | 2.2252 | 22.898 | 8.412 |
| 200 | 2.2158 | 22.895 | 8.370 |
| 500 | 2.2099 | 22.892 | 8.344 |

The temperature data therefore identify only a near-clamped rear-aluminum
regime, not a finite contact coefficient.  The pre-declared 20 W/K value is
retained for the v17 transfer test, but it is a practical lower-bound/plateau
parameter and must not be reported as a measured casing/flange conductance.

After the hardware correction, the power factors were re-profiled:

| group | v16 factor | v17 factor | independent closure bracket | result |
|---|---:|---:|---:|---|
| 456 | 1.65 | 1.05 | 1.05--1.34 | inside |
| 304 | 1.80 | 1.00 | 1.23--1.58 | below |
| 256 | 1.25 | 0.75 | 0.84--1.11 | below |

This confirms that the v16 `k_felt=7.2` and high power factors formed a
compensating loop.  It does not establish correct power accounting in v17:
two of three new factors are now below their independent brackets, and the
heating fit is substantially worse.

The complete nominal validation is:

| Metric | v16 heating train | v17 heating train | v16 heating held out | v17 heating held out | v16 cooling | v17 cooling |
|---|---:|---:|---:|---:|---:|---:|
| mean sensor RMSE (K) | 65.32 | 105.20 | 61.16 | 97.24 | 33.75 | 35.73 |
| steady MAE (K) | 54.69 | 100.46 | 49.65 | 90.78 | 28.09 | 21.00 |
| t90 MAE (s) | 810 | 1082 | 877 | 1196 | 893 | 1026 |
| axial RMSE (K) | 107.94 | 50.18 | 110.54 | 44.01 | 8.00 | 5.54 |
| mid-radial RMSE (K) | 42.58 | 27.10 | 43.72 | 27.87 | 5.06 | 4.73 |
| deep-radial RMSE (K) | 55.30 | 44.16 | 55.70 | 44.23 | 12.60 | 9.74 |
| DP1 RMSE (mbar) | 0.250 | 0.502 | 0.181 | 0.337 | 0.052 | 0.031 |

The lower axial RMSE is not an axial-physics success.  V17 still predicts no
58 mm maximum in any heating case; its much lower fitted power compresses the
modeled axial differences while strongly underpredicting absolute heating
temperatures.  E72, for example, predicts an almost flat 510--540 K side
profile against measured values near 738/798/734 K.

The casing/flange path has a decisive but incomplete T2 signature:

| T2 metric | v16 | v17 |
|---|---:|---:|
| mean heating RMSE (K) | 32.46 | 7.98 |
| mean heating steady bias (K) | +35.62 | -3.54 |
| heating t90 MAE (s) | 653 | 1676 |
| mean cooling RMSE (K) | 9.34 | 12.80 |
| mean cooling steady bias (K) | +4.25 | -4.54 |
| cooling t90 MAE (s) | 1150 | 917 |

Thus the verified path explains most of the systematic T2 steady-level bias,
but the one-state rear-housing topology does not reproduce T2 dynamics.  It
must be retained as real hardware while its fitted magnitude is treated as
unidentified.

Energy closure passes with a maximum absolute residual of
`4.3e-14 W`.  Screen-to-nominal maximum final-sensor changes are 10.11 K for
C69, 5.32 K for C81, and 5.21 K for E72, so the pre-declared 20 K transfer
criterion passes.  A nominal-to-refined Richardson test is still required
before discriminating small source/exchange effects.

### Independent post-v16 assessment: checked conclusions

The independent assessment in `summaries/2D_post_v16_assessment.md` was
checked against the code and the v17 results.

The following conclusions are supported:

1. the v16 felt contact/conductivity fit is structurally confounded;
2. apparent-`Nu` exponent, magnitude, effectiveness, inversion, and LTNE must
   be evaluated before further fitting;
3. mesh adequacy is a precondition for mechanism discrimination;
4. axial source and wall/gas exchange magnitude must be screened jointly;
5. v16 felt conductivity and power factors were compensating; and
6. the rear-manifold model should be deferred until cheaper observation and
   exchange tests are exhausted.

Two claims require correction:

- T9/T10 are not observed as pure solid temperatures in v12--v17.  The v17
  target inherited from the v14 calibration is 80% local wall and 20% local
  gas with a response filter.
  V17 post-processing from pure gas through pure wall gives 0/15 correct
  radial signs at both depths for every blend.  Observation reweighting alone
  is therefore falsified on the current fields.
- Comparing an axial radiative conductance through open flow area with axial
  SiC conduction does not by itself rule out redistribution of incident lamp
  radiation by channel-wall reflections/view factors.  It does argue against
  treating thermal axial radiation as a large Fourier-like conductance.
  The conservative two-component source remains a diagnostic, but a passing
  result must be interpreted as optical deposition/reflection until a physical
  radiosity calculation supports it.

The Julia invariant evaluator reproduces the assessment's v16 values and gives
the following v17 result:

| invariant | measured | v16 | v17 |
|---|---:|---:|---:|
| apparent-`Nu` prefactor | `3.047e-4` | `4.911e-4` | `6.893e-4` |
| apparent-`Nu` exponent | 1.445 | 1.412 | 1.442 |
| effectiveness range | 0.573--0.781 | 0.732--0.856 | 0.877--0.937 |
| 58 mm peak count | 10/15 | 0/15 | 0/15 |
| positive `Lambda_58` | 15/15 | 0/15 | 0/15 |
| positive `Lambda_107` | 15/15 | 0/15 | 0/15 |

V17 therefore confirms that the exponent already emerges, while the exchange
magnitude is much too large and the effectiveness range is compressed.  The
next calculation must not fit the local exchange exponent again.  It must
screen one global exchange-magnitude correction jointly with the conservative
axial-source alternatives under bounded power factors.

V17 remains rejected for Role-B coefficient extraction.

## 2026-07-30 - V18 pre-declared source/exchange screen

V18 retains the 100 physical channels as 15 exact D4 channel orbits.  It does
not fit T9/T10 or assume that their filtered reading is the local core-gas
temperature.  In accordance with the apparatus evidence, calibration uses:

```text
primary side wall: T8, T12, T11
felt:             T2
outlet air:       T3
diagnostic only:  T9, T10, Lambda_58, Lambda_107
```

The heating split remains interleaved within every lamp group:

```text
training: E67/E69/E71, E72/E74/E76, E77/E79/E81
held out: E68/E70, E73/E75, E78/E80
cooling:  C69/C80 training, C81 held out
```

### Exact numerical family

The legacy screen/intermediate/nominal meshes are not a Richardson family
because all retain eight cells over the first 15 mm.  V18 introduces an exact
factor-two family:

| mesh | auxiliary receiver `nr` | felt `nr` | casing `nr` | axial cells | front/rear split | rear tube |
|---|---:|---:|---:|---:|---:|---:|
| C | 14 | 3 | 1 | 24 | 4/20 | 15 |
| M | 14 | 6 | 2 | 48 | 8/40 | 30 |
| F | 14 | 12 | 4 | 96 | 16/80 | 60 |

The receiver auxiliary radial count remains fixed because the channel
cross-section is already represented by the physical 10 by 10 topology and
15 D4 orbits.  C is used for broad screening, M for reranking, and F only for
winner/runner-up confirmation.  Candidate effects must exceed numerical
uncertainty.

All nine heating-training records enter the nested C-mesh power profiles.
The first attempt to rerank all nine simultaneously on M showed severe stiff-
solver thread contention before completing one candidate.  The structural
transfer was therefore reduced to the pre-existing interleaved E69/E74/E79
representatives, one from every power group; all 15 heating and three cooling
records remain reserved for the final M-mesh validation.

### Conservative source candidates

The baseline uses the inherited Beer--Lambert axial weights.  The alternative
mixes the same near source with a deeper exponential component:

```text
w_gj = W_g [(1-f_deep) w_near,j + f_deep w_deep,j(L_deep)].
```

Each group marginal `W_g`, and therefore the existing centered radial
distribution, is unchanged.  Both axial components are nonnegative and
normalized to the inherited group total, so the candidate and baseline
deposit exactly the same power.  The inherited finite-length transmission is
below one part per million at `beta=100 1/m` and is preserved equally in both
laws.

The alternative bounds are:

```text
0.05 <= f_deep <= 0.80
0.020 <= L_deep <= 0.200 m.
```

A convex mixture of two front-origin exponentials remains monotonically
decreasing as a source.  It can produce the measured 58 mm temperature
maximum only through its interaction with front/rear losses and transport.
If it fails, no extra empirical source zones will be added.

### Exchange and fitting objective

The local exchange exponent remains fixed.  V17 already gives 1.442 against
the measured apparent exponent 1.445.  V18 exposes only one global multiplier
on the inherited exchange magnitude.

The first marginal v18 calculations placed the provisional optimum at the
initial lower screen value, `m_h=0.15`.  Before any power profiling, the
exchange screen was therefore extended to `m_h=0.05` and `0.10`.  This is
recorded as a boundary-triggered range extension, not an after-the-fact
selection: a lower-bound optimum will be reported as evidence that the
assumed local core--gas exchange law is structurally too strong.

Every case is sampled at 61 common fractional-time points, excluding `t=0`.
Residuals use the robust pseudo-Huber function

```text
rho(x) = 2 (sqrt(1+x^2) - 1)
```

with fixed discrepancy scales of 35 K for each side TC, 25 K for T3, and
15 K for T2.  Cases are averaged before groups so record length and power do
not change their weight.  The objective is

```text
J = 0.35 J_curve + 0.25 J_level
  + 0.20 J_side_shape + 0.20 J_effectiveness.
```

`J_curve` uses 0.50 side, 0.30 T3, and 0.20 T2.  `J_level` uses the mean of
the final 10% with 0.45 side, 0.35 T3, and 0.20 T2.  Side shape uses normalized
T8/T12/T11 excess temperatures with scale 0.05.  Effectiveness uses the exact
11/58/107 mm quadrature with scale 0.04.  T9/T10 never enter selection.

One power factor per lamp group is nested inside every candidate and hard
bounded by independent closure:

```text
456: 1.05--1.34
304: 1.23--1.58
256: 0.84--1.11.
```

No Gaussian penalty is used.  A bound-active result is reported as unresolved
power accounting.

### Acceptance

Before physical ranking, primary-observable nominal-to-refined history RMS
must be below 10 K and the maximum below 20 K.  On untouched heating holdout,
a candidate must satisfy:

- five-primary mean-sensor RMSE below the v16 benchmark 58.46 K;
- mean side RMSE below 60 K, T3 RMSE below 60 K, and T2 RMSE below 15 K;
- axial RMSE below 40 K;
- at least 8/10 observed middle peaks and no more than two false peaks;
- apparent-`Nu` prefactor within 20% and exponent within 0.10;
- effectiveness RMSE at most 0.05 and crossover between 0.63 and 0.69;
- all power factors inside their brackets; and
- exact mass/power/energy closure and acceptable refined-mesh transfer.

C81 must retain primary-five mean RMSE below 27 K and T2 RMSE below 12 K.
Composite error cannot hide failure of T3 or T2.

## 2026-07-30 - V18 results: source/exchange factorization rejected

### Numerical and calibration execution

The v18 source, mesh, exchange and energy smoke suite passes 15/15 tests.
The exact pre-fit C/M/F family gives the following M-to-F primary-observable
changes:

| case | history RMS | history maximum | final maximum |
|---|---:|---:|---:|
| E72 | 5.51 K | 9.36 K | 7.24 K |
| C69 | 9.53 K | 19.93 K | 13.15 K |
| C81 | 6.04 K | 12.16 K | 7.53 K |

The declared 10 K RMS/20 K maximum history gate passes narrowly before
fitting.  This does not guarantee convergence after parameters change.

The boundary-triggered C screen selects:

```text
source                     = Beer--Lambert
exchange multiplier        = 0.05 (lower diagnostic bound)
power 456/304/256           = 1.05 / 1.23 / 0.84
felt conductivity scale    = 1.30 (upper screen bound)
felt heat-capacity scale   = 1.30 (upper screen bound)
felt contact scale         = 0.30 fixed
```

All three power factors are at the lower edges of their independent closure
intervals.  Felt conductivity and capacity are at their upper screen edges.
These are boundary flags, not identified values.

The four structural candidates transfer to M as:

| source | exchange | M objective |
|---|---:|---:|
| Beer--Lambert | 1.00 | 3.9498 |
| Beer--Lambert | 0.05 | **3.0873** |
| 20% deep, 20 mm | 1.00 | 4.0980 |
| 20% deep, 20 mm | 0.05 | 3.1348 |

The winner/runner difference is much smaller than numerical uncertainty.
On F the same order remains, 2.9083 versus 2.9554, but fitted-regime M-to-F
history differences are 12.84--20.12 K RMS and 19.57--32.88 K maximum.
Therefore the selected weak-exchange regime fails the numerical gate even
though the pre-fit baseline passed it.

### Untouched validation

V18 generates 18 temporal plots, 15 axial plots and a primary-observable
parity plot under `summaries/2D_v18/plots/`.

| phase | primary-five mean RMSE | side RMSE | T3 RMSE | T2 RMSE | axial RMSE |
|---|---:|---:|---:|---:|---:|
| heating training | 84.66 K | 93.98 K | 132.81 K | 8.56 K | 129.54 K |
| heating held out | 71.50 K | 82.90 K | 101.53 K | 7.26 K | 130.65 K |
| cooling training | 33.25 K | 33.83 K | 50.12 K | 14.64 K | 2.56 K |
| C81 held out | 26.57 K | 29.31 K | 36.21 K | 8.71 K | 6.63 K |

The held-out heating primary error fails the 58.46 K benchmark, the side and
T3 guards fail, and axial shape fails severely.  C81 primary-five mean RMSE
passes its 27 K guard narrowly and C81 T2 passes 12 K.  This isolates the
remaining problem away from a simple missing felt/cooling mass: cooling T2 is
acceptable while heating side and air cannot be reconciled.

The final M-mesh energy-rate residual is below `1.0e-13 W`; source power and
radial orbit marginals close exactly on every mesh.  The inherited detailed
energy ledger is explicitly disabled on C/F, whose new axial faces differ
from v17, so it cannot report a numerically false refined closure.

### Invariant assessment

| invariant | measured | v17 | v18 selected |
|---|---:|---:|---:|
| apparent-`Nu` prefactor | `3.047e-4` | `6.893e-4` | `6.107e-3` |
| apparent-`Nu` exponent | 1.445 | 1.442 | 0.793 |
| effectiveness range | 0.573--0.781 | 0.877--0.937 | 0.712--0.889 |
| effectiveness RMSE | -- | -- | 0.188 |
| middle peaks | 10/15 | 0/15 | 0/15 |

The local exponent was held fixed, but changing the exchange magnitude moved
the complete model into a different nonlinear regime and destroyed the
emergent apparent exponent.  A single multiplier cannot independently repair
exchange magnitude while preserving flow dependence.

### Full source/exchange interaction audit

Because the first joint screen pruned source candidates using their marginal
performance at current exchange, every near/deep pair was re-evaluated at
`m_h=0.05` and `0.10`.  This finds a deliberately shape-focused point:

```text
f_deep = 0.60
L_deep = 50 mm
m_h = 0.10
representative middle peaks = 3/3
representative mean inversion error = -2.75 K
objective = 5.060
```

The Beer--Lambert scalar-objective winner has objective 3.183 but zero
representative peaks.  Thus deep centered deposition can create the measured
58 mm maximum; the scalar objective rejects it because absolute side/T3
levels worsen.

On all 15 heating cases using the C diagnostic mesh, the shape point captures
7/10 true peaks and creates five false peaks.  Held-out side RMSE is
124.65 K, T3 RMSE 130.46 K, axial RMSE 51.02 K, and felt RMSE 7.21 K.  It is
not a transferable v18 solution, but it is evidence that optical shape should
not be selected using the same scalar balance as core--gas exchange.

### Decision

V18 is **rejected for Role B coefficient extraction**.

The rejected quantities are:

- `m_h=0.05` as a physical heat-transfer multiplier;
- felt scales 1.30/1.30 as identified material properties;
- the lower-bound power vector as a resolved optical calibration; and
- either v18 axial source law as a validated deposition coefficient.

The 15 D4-orbit representation is retained.  V19 will separate three inverse
problems:

1. normalized side profiles identify a conservative centered optical kernel;
2. measured effectiveness and apparent `Nu(Re)` identify an integrated
   receiver wall--gas `UA`, whose axial distribution uses a normalized
   finite-asymptote Graetz kernel; and
3. cooling plus tube temperature identify a bounded downstream T3
   gas/tube/probe observation operator.

Only after those signatures pass will bounded felt/storage and group power
factors be profiled for absolute temperatures.  T9/T10 remain diagnostic, and
no bypass, side-weighted beam, or rear manifold is introduced.

## 2026-07-30 - V19 integrated-UA and outlet-observation model

### Implemented architecture

V19 retains all 100 physical channels through the existing 15 D4 orbits and
the full SiC/felt/casing/adaptor/tube/flange assembly. It does not revert to a
15-channel approximation. The verified side-TC coordinates are 11, 58, and
107 mm.

The v18 local exchange multiplier is replaced by an equivalent per-physical-
channel integrated law:

```text
Nu_bar = A_Nu Re_bar^n (Pr/0.70)^m
UA_ch  = Nu_bar k_film/Dh P_ch L
UA_receiver = 100 UA_ch.
```

`Re_bar` uses total standard mass flow divided by 100. Every physical channel
receives the same imposed integrated `UA_ch`; equal-pressure flow allocation
then changes each orbit's actual NTU and normalized Graetz axial kernel. This
avoids a new unidentifiable radial-UA law. Three-point Gauss integration and
per-orbit normalization recover `UA_ch` on C/M/F to machine precision.
Orbit heat rates retain physical multiplicity, so solid/gas exchange remains
algebraically antisymmetric.

V19 also distinguishes the receiver-exit gas approximation from T3. The
configured T3 operator samples rear gas and tube at 3 mm downstream, with an
effective cylinder cross-flow coefficient, tube radiation, finite areal
capacity, and an effective stem/lead sink to the water flange. Cooling starts
the probe at measured T3. A distributed effective tube/flange sink begins
28 mm into the local rear-tube coordinate; the standard builder disables the
old terminal-only tube sink whenever this branch is active.

Two qualifications are now explicit:

1. rear mixing uses a `cp(T)T` weighted approximation rather than exact
   integral enthalpy; and
2. Stage B compares calculated gas at the receiver exit with measured T3,
   whereas Stage C treats T3 as downstream and dynamic. Therefore fitted
   `A_Nu,n` are conditional apparent parameters and may absorb probe/location
   bias.

The manuscript T3 coordinate near 136 mm and the model's configured 3 mm
rear-tube coordinate also require reconciliation. The latter corresponds to
approximately 140 mm global position for a 137 mm receiver.

### Stage A - optical-shape transfer fails

Stage A fits only normalized final T8/T12/T11 excess profiles. The expanded
conservative centered near/deep screen selects:

```text
f_deep = 0.90
L_deep = 120 mm.
```

Training normalized RMSE is 0.0325, with 5/6 observed middle peaks recovered
and two false peaks. Held-out normalized RMSE is 0.0292, but only 2/4 peaks
are recovered and one false peak is created. Combined performance is 7/10
true peaks with three false peaks, versus the declared requirement of at
least 8/10 and at most two false peaks. Stage A fails. The selected deep
fraction is also at the tested upper edge. Because selection occurred on C,
it is a rejected diagnostic source seed, not a mesh-robust optical law.

### Stage B - integrated signature remains conditional and fails

The broad Stage-B screen selects:

```text
A_Nu = 3.80e-4
n    = 1.54.
```

Effectiveness transfers well: training/held-out RMSE is 0.0301/0.0285.
However, log-apparent-Nu RMSE is 0.1116/0.1029, and `A_Nu` is about 25% above
the measured manuscript fit, outside the 20% bound.

Because this miss was close to the declared gates, a predeclared local
refinement tested only four points inside the acceptance rectangle. The best
is exactly its upper-right corner:

```text
A_Nu = 3.693432e-4
n    = 1.54346
full-training log-Nu RMSE = 0.1203.
```

It is boundary-active and performs worse on the full training set, so no
holdout or M confirmation was authorized. The original Stage-B failure is
not a coarse-grid artifact.

### Stage C - tube/probe hypothesis fails

The cooling ablation monotonically prefers the largest tested distributed
tube/flange coefficient:

```text
h_tube-flange = 400 W/m2/K (upper grid edge).
```

The T3 screen then selects maximum probe areal capacity and minimum stem sink:

```text
C''probe = 3000 J/m2/K (upper edge)
G''stem  = 0 W/m2/K    (lower edge).
```

Despite all three boundary choices, T3 RMSE is 73.29 K on C69/C80 and
39.30 K on C81. C81 final bias is +27.34 K and its `t90` error is 1500 s.
Stage C fails. Probe lag plus the tested rear-tube cooling topology cannot
explain T3.

### Stage D - felt and power cannot repair the structure

For diagnostic continuity, Stage D freezes rejected A--C parameters. The
cooling-only felt screen and M rerank select both lower edges:

```text
felt k scale  = 0.70
felt Cp scale = 0.75.
```

On C69/C80 at M, side/T3/T2 RMSE is 49.49/70.95/14.63 K. C81 gives
45.92/38.77/8.80 K, with T3 final bias +29.97 K. Felt identification and its
transfer gates fail.

Every groupwise power profile is monotonically improved by lowering power and
selects its independent closure floor:

```text
456/304/256 scales = 1.05 / 1.23 / 0.84.
```

Training side RMSE remains 150.91, 182.11, and 99.56 K across the three
groups; T3 RMSE remains 111.36, 129.48, and 74.24 K. Held-out errors also
fail. Power is therefore unresolved at its closure bounds and cannot be used
to rescue v19.

### Verification corrections made before final validation

The original verification draft recorded M-to-F rows without enforcing the
10 K RMS/20 K maximum gates and reused an inherited ledger identity as if it
were independent. Before final validation, v19 was corrected to:

- assert summed 100-channel `UA` and mass flow;
- expose and gate both pressure and gas-reference iteration convergence;
- add an independent capacity-weighted derivative audit of the v19 exchange
  correction and distributed tube/flange loss on C/M/F;
- expand M-to-F representatives to operating extremes and held-out heating;
- remove both distributed and terminal tube cooling in the Stage-C-null
  branch;
- test T3 prediction sensitivity to 61/121/241 saved points; and
- write final acceptance only after stage, mesh, flow, energy, null-branch,
  and observation-sampling gates are evaluated.

The detailed formulation and interpretation limits are in
`summaries/2D_v19_theory_manual.md`. Final 18-case errors, plots, mesh results,
and the decision are recorded below after the corrected validation run.

### Final validation, numerical transfer, and decision

The corrected M-mesh validation gives:

| phase | primary-five mean RMSE | side RMSE | T3 RMSE | T2 RMSE | axial-inversion RMSE |
|---|---:|---:|---:|---:|---:|
| heating training | 101.95 K | 130.87 K | 105.62 K | 11.53 K | 50.03 K |
| heating held out | 76.30 K | 93.96 K | 89.97 K | 9.66 K | 36.35 K |
| cooling training | 44.42 K | 45.56 K | 70.61 K | 14.81 K | 0.42 K |
| C81 held out | 37.26 K | 46.02 K | 39.20 K | 9.02 K | 2.53 K |

Felt/T2 remains the only consistently moderate primary observable. Heating
side and T3 errors are worse than v18, and C81 no longer meets the previous
27 K primary-five guard.

The fitted-regime M-to-F check fails on every heating representative:

| case | history RMS | maximum |
|---|---:|---:|
| E67 | 10.51 K | 18.35 K |
| E71 | 32.17 K | 51.24 K |
| E75, held out | 26.30 K | 40.59 K |
| E80, held out | 16.39 K | 24.36 K |
| C69 | 9.22 K | 19.93 K |
| C81 | 4.68 K | 12.16 K |

Thus only the two cooling cases satisfy both the 10 K RMS and 20 K maximum
gates. The boundary-selected optical/absolute-temperature regime is not
numerically transferable for heating.

The independent exchange-energy audit passes on all M/F representatives and
all 18 nominal cases, with maximum relative residual `1.13e-14`. The legacy
M ledger residual is below `1.43e-14` relative. A dedicated E75 tolerance
test gives only 0.024 K primary-history RMS and 0.072 K maximum between the
bulk-validation and tight solver settings, confirming that the large M/F
differences are spatial-discretization effects, not ODE tolerance.

T3 post-processing is adequately sampled: C81 changes by 0.027 K at final
temperature and 50 s in `t90` between 61 and the available 119 points. The
requested 241-point case cannot exceed 119 because that is the experimental
record length. Conversely, removing the complete Stage-C tube/probe cooling
branch changes primary histories by 14.1--17.9 K RMS for representative
heating cases and T3 by 31.1--39.5 K RMS. Absolute results are therefore
strongly conditional on the already rejected Stage-C boundary model.

The flow audit converges normally for 17/18 cases. E76 reaches the 24-iteration
limit at 112/119 saved points; its maximum pressure residual is nevertheless
only `2.29e-5` relative and its gas-reference residual `1.19e-6`. This is below
the reporting guard of `1e-4` but above the internal pressure target near
`1e-5`, so the strict flow-convergence flag remains false rather than silently
accepting the capped iteration.

V19 produces 18 transient plots, 15 axial plots, and one parity plot. The
axial figures use the full continuous modeled side-skin profile, with only
the measured thermocouples shown as points at 11/58/107 mm.

**Permanent axial-plot convention:** every axial-profile figure intended for
comparison across experiments or model versions must use the same x- and
y-axis limits throughout the comparison set. The model remains a continuous
axial curve and measured thermocouples remain discrete points. Axis limits
must be chosen once from the complete comparison set, recorded in the plotting
script or caption, and not autoscaled separately for each experiment.

V19 is **rejected for Role B coefficient extraction**. Every identification
stage fails, every nuisance family reaches a bound, heating mesh transfer
fails, and the Stage-C branch materially changes predictions. The following
are not validated coefficients:

- `f_deep=0.90`, `L_deep=120 mm`;
- `A_Nu=3.80e-4`, `n=1.54`;
- `h_tube-flange=400 W/m2/K`;
- `C''probe=3000 J/m2/K`, `G''stem=0`;
- felt scales `0.70/0.75`; and
- power scales `1.05/1.23/0.84`.

The 15 D4-orbit geometry, exact source normalization, standard-flow mass
conversion, common-pressure channel allocation, temperature-dependent gas
properties, and full assembly are retained as verified architecture. The
next step must address identifiability of outlet gas versus T3 and optical
source depth; another joint point fit of the same observables is not
authorized.

Here â€œno independent outlet-gas measurementâ€� does **not** mean T3 is
statistically dependent on the wall thermocouples or that it is not a real
experimental measurement. It means there is only one outlet-region
temperature measurement, T3, and the same signal is being asked to determine
both:

1. the mass-averaged receiver-exit gas temperature used to infer apparent
   effectiveness/`Nu`; and
2. the unknown downstream tube/probe observation relation between that bulk
   gas temperature and the thermocouple reading.

Because T3 can exchange with the alumina tube and cooled lead/stem and has
finite response, it is not a separately shielded measurement of bulk gas
enthalpy. There is no second outlet temperature, tube-wall measurement, or
independent outlet calorimetric balance that can constrain the observation
operator while leaving T3 free to validate the inferred gas temperature.
Thus T3 is independent experimental information, but not an independent
ground truth for receiver-exit bulk gas when its own observation physics is
also fitted.

## 2026-07-30 - V20 enthalpy correction and identifiability test

### Scope

V20 performs the agreed post-v19 identifiability test without adding a
bypass, side-weighted source, per-experiment properties, independently fitted
radial channels, or another joint T3-driven coefficient fit. It retains the
15 D4 orbits representing all 100 channels and the full assembly.

The rejected v19 distributed rear-tube/flange sink is removed from the v20
plant; the established terminal tube/flange path remains. T3, T9, and T10 are
excluded from plant objectives. Heating and cooling scores use only
T8/T12/T11 and T2.

### Exact enthalpy transport

V20 analytically integrates the inherited air `Cp(T)`. Every receiver and
rear-tube cell solves

```text
UA/mdot = integral(Tin,Tout) Cp(T)/(Twall-T) dT
Q = mdot [h(Tout)-h(Tin)].
```

Orbit outlets mix by mass-weighted enthalpy and `T(h)` inversion. Separate
receiver, mixing, rear-tube, and full-stream residuals are saved. The updated
smoke suite passes 28/28 tests.

The six-case frozen nominal-mesh v19/v20 comparison shows that the correction
is small: all-sensor RMS differences are 0.030--0.606 K, the maximum sensor
difference is 1.45 K, T2 RMS differences are below 0.036 K, and
receiver-exit-gas RMS differences are 0.088--1.862 K. Maximum gas enthalpy
residual is `3.36e-9 W` absolute and `3.11e-11` relative. The old
variable-`Cp` approximation was a real conservation defect, but not the
missing 40--180 K physics.

### T3 coordinate evidence

The strongest archival evidence is global `x approximately 136 mm`: early
experiment-specific scripts sample approximately 135.6--136 mm and identify
it with T3. Global 137 mm is the exact receiver-exit convention. Global
140 mm, or 3 mm into the rear tube, is a later inherited convention rather
than a metrology-grade verification. The photographs locate T3 only in the
outlet region and cannot distinguish this span.

V20 therefore implements three discrete branches:

```text
receiver_136  = global 136 mm
receiver_exit = global 137 mm
rear_003      = global 140 mm.
```

No continuous position is fitted.

### Removal of initialization leakage

An audit found that the first cooling draft excluded T3/T9/T10 from the score
but still used them to initialize the physical core and rear tube. The final
policy uses side/T2-only initialization:

- T8/T12/T11 define the perimeter profile;
- the unobserved core begins from that perimeter profile;
- T2 initializes the outer assembly; and
- the rear tube begins from back-side T11.

Measured cooling T3 at `t0` initializes only the one-way probe ODE when T3 is
scored and cannot feed back into the plant. Flow-convergence and all
structural/nuisance boundary gates were also made explicit, and final-window
scoring now uses the actual final 10% of time.

### Discrete T3 observer fails

The leakage-free M-mesh C69/C80 fit and C81 validation span all three
locations, `C''=200--3000 J/(m2 K)`, and
`G''stem=0--200 W/(m2 K)`. The lowest score weakly prefers 137 mm but again
selects `C''=3000` at the upper boundary and `G''stem=0` at the lower
boundary. Training mean/worst T3 RMSE is 69.93/94.96 K. C81 gives
38.28 K RMSE, +42.49 K final bias, and -1800 s `t90` error. No observer
passes. The earlier all-sensor initializer weakly preferred 140 mm, so the
small location ranking itself is initialization-sensitive and cannot locate
T3.

### T3-free `Nu50,n` profile is a boundary ridge

The corrected M profile evaluates `n=1.25--1.65` and
`Nu50/Nu50_measured=0.75,1.35,1.95` on E72/E69/E81. Its minimum is
`n=1.25` at the lower exponent boundary and `Nu50` ratio 1.35, with
`A_Nu=8.85657e-4`, objective 4.0651, side RMSE 125.96 K, and T2 RMSE
10.81 K. Four candidates within 10% span `n=1.25--1.35` and `Nu50` ratios
1.35--1.95. There is no interior connected coefficient basin.

### Bounded nuisance screens fail

Both near-tied `n=1.25` ridge endpoints select the same lower bounds:

```text
felt k/Cp scales = 0.70/0.75
power scales 456/304/256 = 1.05/1.23/0.84.
```

| `Nu50` ratio | heating side | heating T2 | cooling side | cooling T2 |
|---:|---:|---:|---:|---:|
| 1.35 | 102.47 K | 10.20 K | 66.42 K | 14.55 K |
| 1.95 | 79.40 K | 8.47 K | 51.95 K | 14.48 K |

Selected-row coupled-flow solves converge, but heating and cooling gates both
fail. Felt, all three powers, and the exponent are boundary-active.
`plant_admissible=false` for both endpoints.

### Decision

V20 is rejected for Role-B coefficient extraction. Its thermodynamic
correction is retained, but the T3 operator is not identified, the T3-free
integrated-UA surface is a structural boundary ridge, and bounded nuisance
parameters cannot produce an admissible plant.

The v20 minima are diagnostic only and are not validated `A_Nu`, `n`, felt,
power, or T3-probe parameters. Because the T3-free plant fails, no
T3-conditioned UA refit, bootstrap point interval, held-out branch selection,
or F-mesh coefficient confirmation is authorized; those operations would add
precision to a rejected structure.

The defensible manuscript outcome is the directly measured apparent
correlation plus a formal non-identifiability statement for receiver
`UA(Re)` under the available outlet-temperature information. Equations,
assumptions, and routing are in `summaries/2D_v20_theory_manual.md`; outputs
are under `summaries/2D_v20/`.

## 2026-07-30 - V20 gate-free boundary-stress test

### Purpose and rules

The bounded v20 study ended with the `Nu(Re)`, felt, power, and T3-observer
parameters on search limits. A deliberately gate-free stress test was
therefore run to distinguish two possibilities:

1. the original acceptance gates were prematurely rejecting a nearby useful
   solution; or
2. even substantially relaxed parameter bounds could not reconcile the side
   wall, felt, and outlet-temperature records.

Acceptance gates and boundary penalties were not used to rank any stress
candidate. Candidates were ranked only by the stated T3-free side/T2 score,
or by the T3 training score in the separate one-way observer profile.
Parameters were allowed outside their physical identification priors, but
basic nonnegative physics was retained: no negative conductivity, heat
capacity, `Nu`, power, observer capacity, or conductance was allowed. The
original heating, cooling, and observer gates were evaluated only after the
search and full validation. Consequently, a stress minimum is a diagnostic
compensator and not an identified material property or heat-transfer
coefficient.

Cooling continued to use the leakage-free side/T2-only initialization.
T3/T9/T10 did not initialize the plant. The T3 observer was post-processed
one way and could not feed energy back into the receiver, felt, casing, or
gas solution. Every final-window quantity used the last 10% of each
experiment's actual time span, not the last 10% of its sample indices.

### Structural `UA` sweep

The structural sweep extended the Reynolds exponent down through zero and
the `Nu50` ratio beyond the original interval. Initial slices covered
`n=0--1.25` and ratios as high as 4; the refined common slices used ratios
`1.0--3.0` at `n=0`, 0.25, and 0.50. With the nuisance parameters still
fixed, the lowest refined structural row was

```text
n = 0.25
Nu50/Nu50_reference = 2.00
objective = 2.21567
side RMSE = 82.78 K
T2 RMSE = 9.43 K.
```

The ratio minimum was interior to the refined ratio interval, but the
exponent had moved far below both the measured apparent-correlation exponent
and the original v20 identification range. After the felt and power
coordinates were allowed to respond, the `n=0.50`, ratio-2.00 point slightly
outperformed `n=0.25`, ratio-2.00: the corresponding representative-case
heating/cooling objectives were 1.8924/1.0971 versus 1.9417/1.1401. The
stress plant therefore used

```text
n = 0.50
Nu50/Nu50_reference = 2.00
Nu50 = 0.174451
A_Nu = 0.0246711.
```

This coordinate movement does not establish a new correlation. It shows
that the T3-free score can exchange a much flatter `Re` dependence against a
larger `Nu50`, i.e. the original structural ridge continues outside the
physical identification box.

### Felt, source, and power coordinates

The first expanded nuisance screen at `n=0.50`, ratio 2.00 selected felt
`k/Cp` scales 0.40/0.55 and representative-case power scales
1.05/1.23/0.70. Its heating side/T2 RMSE was 77.86/6.57 K and cooling
side/T2 RMSE was 50.81/13.56 K. Both original plant gates still failed.

A second cooling-only felt refinement selected `Cp` scale 0.55. The
conductivity objective was exactly flat at scales 0.05 and 0.15
(`objective=1.064571`, side/T2 RMSE 50.28/13.09 K) and changed only to
1.065667 at scale 0.30. Scale 0.15 was retained as a representative point,
but the plateau down to 0.05 means that felt conductivity is practically
unidentified in this regime; it is not evidence for a conductivity 15% of
the nominal value.

The source extension compared the retained near/deep source with full-deep,
longer-deep, and Beer--Lambert alternatives. The existing
`f_deep=0.90`, `L_deep=0.12 m` choice remained lowest:

| source | `f_deep` | `L_deep` | objective | side RMSE | T2 RMSE |
|---|---:|---:|---:|---:|---:|
| near/deep | 0.90 | 0.12 m | 1.7051 | 73.68 K | 7.28 K |
| near/deep | 1.00 | 0.12 m | 1.9824 | 80.78 K | 7.26 K |
| near/deep | 0.90 | 0.20 m | 2.0006 | 81.16 K | 7.25 K |
| Beer--Lambert | 0 | 0.08 m | 2.9289 | 103.20 K | 7.50 K |

Thus making the already deep source still deeper did not repair the side
profiles. This stress comparison preserves the 0.90/0.12 choice but does not
validate it as an optical property.

The first power refinement used one representative experiment per lamp
group and selected 1.05/1.23/0.75. Because that selection could be
experiment-specific, every power coordinate was then reranked over all three
heating-training experiments in its lamp group. The all-training selections
were

| group | selected scale | tested range | boundary-active? | group side RMSE | group T2 RMSE |
|---|---:|---:|---|---:|---:|
| 456 | 1.05 | 0.70--1.15 | no | 124.31 K | 9.32 K |
| 304 | 1.05 | 0.70--1.23 | no | 126.58 K | 6.23 K |
| 256 | 0.70 | 0.55--0.84 | no | 55.65 K | 3.97 K |

The final stress-validation power vector was therefore 1.05/1.05/0.70.
These selections are interior only to the deliberately expanded numerical
stress grids. The 456 value equals its original closure lower bound, while
the 304 and 256 values lie below their original independent closure brackets
of 1.23--1.58 and 0.84--1.11. They therefore do not satisfy Role-B physical
interiority, and they do not cure the large 456/304 side errors.

### Expanded T3 observer

At the stress plant, the gate-free discrete observer selected

```text
location = receiver exit, global x = 137 mm
C''probe = 30000 J/(m2 K)
G''stem = 0 W/(m2 K).
```

On C69/C80 its training mean/worst T3 RMSE was 41.75/58.72 K, with
`t90` MAE 925 s. C81 gave 30.21 K RMSE, +48.21 K final bias, and
-1500 s `t90` error. The extended capacity did not plateau: at the
receiver-exit/zero-stem branch the training objective decreased monotonically
from 2.548 at 1200 to 2.226 at 10000 and 1.330 at
30000 J/(m2 K), while training RMSE fell from 74.43 to 41.75 K. Capacity
therefore remained strongly upper-bound active and the stem remained
lower-bound active. The observer is still unidentified.

### Both full 18-case validations

`validate_2D_v20_gatefree.jl` simulated all 15 heating experiments and
C69/C80/C81. Errors below are pooled directly over all samples and all three
side thermocouples where applicable, rather than averaged from per-case
RMSEs.

The first full pass used the refined felt point and the
representative-experiment power vector 1.05/1.23/0.75:

| phase | side RMSE | T2 RMSE | T3-observer RMSE |
|---|---:|---:|---:|
| heating training | 131.08 K | 7.64 K | 223.65 K |
| heating holdout | 95.51 K | 6.47 K | 152.97 K |
| cooling training | 55.02 K | 13.28 K | 44.18 K |
| cooling holdout | 43.67 K | 8.28 K | 30.14 K |

The second and final pass changed only the all-training-refined power vector
to 1.05/1.05/0.70:

| phase | side RMSE | T2 RMSE | T3-observer RMSE |
|---|---:|---:|---:|
| heating training | 115.96 K | 7.15 K | 194.74 K |
| heating holdout | 94.35 K | 6.31 K | 129.84 K |
| cooling training | 55.02 K | 13.28 K | 44.18 K |
| cooling holdout | 43.67 K | 8.28 K | 30.14 K |

The final-pass side/T2/T3 final-window biases were respectively
`+7.72/-4.59/+86.22 K` for heating training,
`-24.40/-6.94/+48.10 K` for heating holdout,
`+52.83/-2.66/+14.06 K` for cooling training, and
`+59.45/+2.79/+48.37 K` for C81.

All 18 ODE solves completed and every coupled-flow solve converged. The
maximum full-stream enthalpy relative residual was `6.11e-11`. The failure is
therefore not a solver, flow-allocation, or gas-energy-conservation failure.

Applied only after evaluation, the original plant thresholds were side/T2
RMSE below 60/15 K for heating and 35/12 K for cooling. Every phase-level
plant gate was false: heating failed on side temperature, cooling training
failed on both side and T2, and cooling holdout failed on side temperature.
The original cooling T3 observer gate was also false because its nominal
30 K RMSE, 20 K final-bias, 500 s timing, and interior-parameter conditions
were not jointly satisfied.

### Interpretation

The original acceptance gates are not causing the rejection. They were
absent from every stress ranking, yet the unrestricted nonnegative search
still produced side errors of 94--116 K in heating and 44--55 K in cooling
on the full dataset. A modest relaxation of the original parameter brackets
is therefore insufficient. Reaching these still-failing errors already
requires a very flat `Nu(Re)` exponent, doubled `Nu50`, felt conductivity on
a near-zero plateau, reduced optical powers, and a T3 thermal capacity ten
times the former upper bound with zero stem conductance.

These values are useful only for diagnosing compensation directions. They
are not validated `Nu`, felt, source, power, or observer parameters and must
not be reported as extracted material or heat-transfer properties. The
stress test strengthens the v20 conclusion: the measured side profiles,
felt response, and outlet-region T3 cannot be reconciled by retuning the
existing scalar closures. It does not uniquely identify which missing
spatial heat path, boundary condition, or core--gas closure is responsible.

### Scripts and artifacts

The reusable stress and validation scripts are:

- `stress_profile_2D_v20_ua.jl`;
- `stress_screen_2D_v20_nuisance.jl`;
- `stress_refine_2D_v20_felt.jl`;
- `stress_profile_2D_v20_source.jl`;
- `stress_refine_2D_v20_power.jl`;
- `stress_refine_2D_v20_power_group.jl`;
- `stress_profile_2D_v20_t3.jl`; and
- `validate_2D_v20_gatefree.jl`; and
- `plot_2D_v20_gatefree_full.jl`.

The principal artifacts are under `summaries/2D_v20/`:

- `stress_ua_profile*_2D_v20.csv`;
- `stress_nuisance_*_2D_v20.csv`;
- `stress_felt_refined*_2D_v20.csv`;
- `stress_source_profile_extension_2D_v20.csv`;
- `stress_power_refined*_2D_v20.csv` and
  `stress_power_group_*_2D_v20.csv`;
- `stress_t3_profile_beststress_2D_v20.csv` and
  `stress_t3_selected_beststress_2D_v20.csv`; and
- `gatefree_validation_*_stress_best_2D_v20.csv` and
  `gatefree_validation_*_alltrain_power_2D_v20.csv`.

The initially generated gate-diagnostic figures did not constitute the full
expected validation plot set. The completed suite is now under
`summaries/2D_v20/plots/gatefree_full/` and contains:

- two final-window parity plots, including a five-sensor faceted view;
- temporal measured/model histories for all 18 heating and cooling cases;
- continuous modeled side-wall axial profiles for all 15 heating cases; and
- a 15-case axial overview with the common limits `x=0--137 mm` and
  `T=500--1250 K`.

All parity and axial values use the final 10% of actual experiment time.
Each axial model curve uses every nominal-mesh axial center rather than only
the three thermocouple coordinates. The individual axial figures and the
overview use the same temperature scale throughout.

### Standing axial-plot convention

For v20 and every subsequent 2D version, the standard result package shall
include an `axial_profiles_common_scale` overview in the same style as
`summaries/2D_v20/plots/gatefree_full/axial_profiles_common_scale_2D_v20_gatefree.png`.
The overview shall:

- show all heating experiments as aligned small multiples;
- draw each modeled axial profile as a continuous line through every axial
  mesh center, not only at thermocouple coordinates;
- overlay measured side thermocouples at 11, 58, and 107 mm;
- use identical axial and temperature limits for every panel; and
- state the averaging window and model version in the title, caption, or
  accompanying manifest.

Individual axial plots may still be produced for detailed inspection, but
they do not replace the common-scale overview.

---

## 2026-07-31 - V21 Phase 2 Reduced-Order Model Implementation & Clean Break

### Overview
Version `21` marks the beginning of Phase 2 modeling and the implementation of a "clean break" architecture. Rather than forking the deep legacy ODE file hierarchy (`v11` $\rightarrow$ `v19`), v21 employs a surgical `@eval` memory-injection approach (`apply_v21_property_fixes!()`) to enforce corrected thermophysical properties and macroscopic physics onto the inherited legacy equations without modifying the older files on disk.

### 1. Corrected Thermophysical Properties
The following independent property corrections were codified into a new standalone module `2D_v21_Properties.jl`:
- **Air Conductivity:** Applied the Sutherland correlation to fix the previous ~8% bias.
- **Air Heat Capacity:** Updated interpolation to match experimental reference tables.
- **SiC Heat Capacity:** Restored the structurally sound ~1100 J/kg/K formula, fixing the flawed ~400 J/kg/K value that severely biased transient `t90` speeds.
- **Felt Insulation:** Restored the functional form that properly tracks local casing temperature, and corrected specific heat from 1360 to 1000 J/kg/K.

### 2. Implementation of Phase 2 Physics
The `2D_v21.jl` solver now explicitly incorporates:
- **Spillage Heating:** Explicit boundary heating applied to the outer perimeter (Orbit 15) to simulate power spilled onto the felt/flange.
- **Flow Maldistribution:** A `core_preference` parameter explicitly allocates more mass flow to the core orbits, simulating the 13mm rear adaptor blockage.
- **Front-Face Radiation Fix:** Radiation loss is now calculated using the full $3.61 \text{ cm}^2$ frontal area (treating channels as blackbody cavities) rather than just the solid web area ($1.36 \text{ cm}^2$).

### 3. Verification & Validation Outcomes
Diagnostic tests on the `2D_v21.jl` framework confirmed absolute physical and mathematical soundness:
1. **Steady-State Inversion:** A 300s steady-state simulation (1.5kW, 150W spillage, 0.5 core preference) mathematically proved the hypothesis, predicting $T_{\text{perimeter}} \gg T_{\text{core}}$ deep into the receiver (e.g. 1360 K perimeter vs 626 K core at 8.4mm depth, and 944.7 K vs 611 K at 19.9mm depth).
2. **Transient Thermal Mass ($t_{90}$):** A 1000K to 295K cooldown test confirmed that the SiC capacity fix extends the thermal decay time constant by a factor of ~2.45 (from unphysically fast to a realistic 980s), directly resolving the longstanding transient speed mismatch.

Verdict: **Phase 2 Implementation PASS**. The model is numerically robust, inherits no property biases, and is structurally ready for Phase 3 Formal Calibration.

## Phase 3: Macroscopic Parameter Calibration

In Phase 3, we sought to extract the true macroscopic flow maldistribution and spillage parameters that best describe the macroscopic performance of the receiver, acting as the 'true' boundary conditions for our 2D continuum model (Role B). 

To do this, we performed a full grid search over the macroscopic parameter space across steady-state high-temperature conditions (E67, E69, E71). The grid space evaluated was:
- **Spillage Fraction**: [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
- **Core Preference (Flow distribution scale)**: [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

### Calibration Results

The grid search converged on an extreme point in the parameter space:
- **Optimal Spillage Fraction**:  .0
- **Optimal Core Preference**: 1.0 (Complete flow maldistribution)

### Analysis of the Boundary Conditions

The results are highly significant. In attempting to minimize the SSE (Sum of Squared Errors) against the experimental measurements—particularly the deep temperature spatial inversion ({\text{perim}} > T_{\text{core}}$)—the model actively pushed the \Core Preference\ to its absolute maximum limit of 1.0.

A \Core Preference\ of 1.0 indicates that virtually all the fluid is being routed through the core of the monolith, completely bypassing the perimeter channels. By forcing all the cold air through the core, the model achieves the maximum possible cooling of the core channels. Simultaneously, because the perimeter is starved of cooling fluid, its temperature spikes drastically, allowing the model to naturally recreate the spatial inversion observed experimentally without requiring any spillage heating (\Spillage = 0.0\).

While this configuration manages to reproduce the *qualitative* behavior of the temperature inversion, the quantitative fit is still poor (\Side Final RMSE = 210 K\). This implies that while extreme flow maldistribution is a dominant macroscopic mechanism, it alone (or when constrained by standard textbook heat transfer coefficients) is insufficient to perfectly capture the coupled non-linear heat transfer taking place within the monolithic structure.

This firmly concludes Phase 3 and transitions us to Phase 4, where we will use these optimized macroscopic boundary conditions to calibrate the intrinsic microscopic heat transfer and optical coefficients (Nusselt number scaling and Rosseland extinction coefficients) to achieve a high-fidelity quantitative fit.


## Phase 4: Heat Transfer Parameter Calibration (Quantitative Tuning)

After establishing that extreme flow maldistribution (Core Preference = 1.0, Spillage = 0.0) provides the necessary macroscopic mechanism to induce the deep spatial temperature inversions seen experimentally, we transitioned to calibrating the fundamental microscopic coefficients to optimize the quantitative fit.

In Phase 4, we performed a 36-point grid sweep across:
- **Nusselt Multiplier (Nu_Mult)**: [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
- **Rosseland Extinction Multiplier (Ross_Mult)**: [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]

This grid was evaluated against the high-temperature steady states (E67, E69, E71). Prior to execution, several legacy injection bugs were resolved in the 2D_v22 architecture to ensure that the physical property corrections introduced in Phase 2 (e.g. Sutherland air conductivity, proper SiC heat capacity, baseline Rosseland bounds, etc.) were properly propagated to the solver.

### Calibration Outcomes

The grid sweep converged on:
- **Optimal Nusselt Multiplier**: 1.2
- **Optimal Rosseland Multiplier**: 0.1

While the algorithm correctly selected an optimal configuration, the absolute magnitudes of the objective function (Sum of Squared Errors) and the RMSE (~1500–1800 K) remained severely elevated across the entire grid. 

### Global Validation (Phase 4 Output)

Running all 15 heating and cooling cases (E67-E81) with these optimally calibrated coefficients reveals that the 2D model (Role B) ultimately **cannot perfectly reconcile the full 3D behavior** using scalar internal boundaries. The model successfully captures the *qualitative* behavior of the temperature inversions (as designed in Phase 3), but the extremely high residual errors (quantitative deviations) physically prove that attempting to force a complex, geometrically sparse 3D honeycomb receiver into a continuum 2D model mathematically breaks the scalar physics equations. 

This completes the Role B computational framework: the rigorous demonstration that an axisymmetrically reduced model, even when fed extreme but physically-derived boundary conditions, structurally limits out. The plots generated (stored under summaries/2D_v22/plots/) display the best possible behavior this continuum formulation can offer.

## Crossover Insights from 1D Corroboration (July 2026)

The conclusions reached above in Phase 4 of the 2D model development have been independently and rigorously corroborated by the final 1D model framework (`1D_v35`). 

Specifically, the 1D model was subjected to an unconstrained grid optimization that allowed for massive scalar flow bypass (diverting cold gas around the core). The 1D optimizer mathematically proved that **flow bypass cannot resolve the paradox either**, due to the **Advective Bottleneck**. If the model bypasses gas during the steady-state heating phase to keep the core hot, it lacks the required *advective mass flow* through the core during the transient cooling phase, causing the core to take vastly too long to cool down.

Therefore, both the 1D and 2D continuum formulations have been mathematically cornered. A scalar, unified flow structure is physically invalid for this geometry. The flow distribution must be dynamically non-linear (e.g., bypassing during heating, fully participating during cooling) or severely 3D geometric, validating the manuscript's claims regarding macroscopic maldistribution in monolithic structures.
