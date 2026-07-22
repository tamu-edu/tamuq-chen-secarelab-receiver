# 1D Model Journal

Brief development notes for the receiver 1D model family, from the early v1 scripts through the staged rear-mass v6 model.

## v1: Baseline 1D Receiver Formulation

Files:

- `1D_v1.jl`
- `1D_v1.exp.67.jl`
- `1D_v1.exp.71.jl`
- `1D_v1.exp.All.jl`
- `1D_v1.exp.All.PEtab.jl`

Purpose:

- Establish a first axial receiver model.
- Test 1D transient solid/gas behavior against selected experiments.
- Explore fitting workflows and PEtab-style calibration attempts.

Main characteristics:

- More exploratory structure.
- Mixed model, experiment, fitting, and plotting logic.
- Useful for proving that axial resolution mattered, but not yet a clean reusable workflow.

Lessons:

- A lumped model misses axial thermocouple behavior.
- The code needed clearer separation between model physics, data import, calibration, and post-processing.

## v2: Standalone Refactor

Files:

- `1D_v2.jl`
- `import_exp_1D_v2.jl`
- `test/test_1D_v2.jl`

Purpose:

- Refactor the 1D model into a cleaner research-script flow.
- Remove dependency on `1D_v1.jl`.
- Add a reusable external data import contract.

Main characteristics:

- Own numerical core.
- Structured sections for geometry, model construction, data import, solver interface, optimization, and post-processing.
- Introduced a more explicit finite-volume solid/gas model.
- Added tests for coexistence with earlier model versions.

Lessons:

- A clean data import layer is essential.
- The model needed a more conservative heat balance and simpler top-level execution behavior.

## v3: Conservative Finite-Volume Model

Files:

- `1D_v3.jl`
- `run_1D_v3.jl`
- `test/smoke_1D_v3.jl`
- `test/test_1D_v3_static.jl`

Purpose:

- Create a conservative finite-volume receiver model with a clean model/function file and a separate runner.
- Use Julia `Optimization` / `OptimizationNLopt` for calibration.

Main characteristics:

- Dynamic states are axial solid temperatures.
- Gas is quasi-steady plug flow, marched cell-by-cell with exact NTU updates.
- Conservative axial conduction between neighboring cells.
- Cooling-first and heating-second calibration.
- Separate `run_1D_v3.jl` for full workflow: calibration, plots, metrics.
- Later updated from direct empirical `h_ref * flow^n` to the A/B/C Nusselt form:

```text
Nu = 10^A Re^B Pr^(10^C)
h = Nu k_f / Dh
```

Lessons:

- v3 gave a much cleaner model/workflow separation.
- The model captured front thermocouples reasonably well.
- Deeper thermocouples and gas temperature showed stronger simulated flow dependence than the experiments.
- The discrepancy likely required changing heat-exchange structure, not only retuning scalar parameters.

## v4: Axial Heat-Exchange Shape and Paired Thermocouples

Files:

- `1D_v4.jl`
- `run_1D_v4.jl`
- `test/smoke_1D_v4.jl`

Purpose:

- Test whether axial gas-solid exchange weakening could reduce the excessive downstream flow dependence.
- Include paired thermocouple averages because the 1D model has no radial coordinate.

Main characteristics:

- Kept the A/B/C `Nu(Re,Pr)` heat-transfer law.
- Added a Graetz-style axial exchange multiplier inspired by earlier `0D_v3.jl` experiments.
- Used paired experimental channels:

```text
T9_pair  = (T9 + T12) / 2
T10_pair = (T10 + T11) / 2
```

Main outcome:

- The fitted axial shape strongly suppressed downstream exchange, often down to the imposed lower clamp.
- Several parameters ran to bounds:
  - very high front convection/emissivity
  - high optical attenuation
  - very low side/rear loss
- Heating fits remained structurally poor for many cases.

Lessons:

- The model does appear to need weaker downstream exchange.
- The v4 Graetz-style parameterization was too flexible/awkward and encouraged bound-hitting fits.
- Paired thermocouple averaging is useful and should remain.

## v5: Reduced-Parameter Staged Model

Files:

- `1D_v5.jl`
- `run_1D_v5.jl`
- `test/smoke_1D_v5.jl`

Purpose:

- Reduce fitted degrees of freedom.
- Use cooling data to identify thermophysical and gas-exchange behavior before fixing those parameters in heating simulations.
- Replace the v4 axial exchange function with a simpler interpretable shape.

Main characteristics:

- Reduced to 12 total parameters.
- Evolved through two calibration designs:
  - first: cooling-first thermophysical calibration followed by heating optical/front calibration
  - later: heating base, cooling thermophysical, heating refit
- Cooling/thermal parameters:

```text
gamma_C, k_scale, U_side, U_rear,
A_Nu, h_floor, L_h, tau_T3
```

- Fixed optical/input-energy parameters after the gamma-source trial:

```text
eta_abs = 0.80
beta_opt = 50.0
front_dep = 0.50
```

- Other fixed values:

```text
B_Re = 1.0
C_Pr = 0.5
eps_front = 0.95
```

- Axial exchange shape:

```text
s(z) = h_floor + (1 - h_floor) exp(-z / L_h)
```

- Optical deposition was reset to the `1D_v1.jl` parameterization:

```text
Beer-Lambert cell-integrated weights
+ front-deposition fraction added to the first cell
```

- Keeps paired thermocouple channels:

```text
T9_pair  = (T9 + T12) / 2
T10_pair = (T10 + T11) / 2
```

Intended lesson:

- Heating should identify gas heat-transfer behavior in the irradiated receiver.
- Cooling should identify thermophysical properties and hidden thermal storage/losses.
- Incoming energy should be adjusted manually rather than fitted through optical parameters during these structural tests.

Observed issue:

- Cooling profiles improved relative to earlier attempts but still showed excessive cooling speed, especially in `T10_pair` and `T3`.
- Some cooling cases, notably C69, retained nonphysical shape features.
- This suggested that a single global heat-capacity multiplier was insufficient; a missing rear/deep thermal reservoir was likely.

## v6: Rear Thermal Mass and Fixed Front Convection

Files:

- `1D_v6.jl`
- `run_1D_v6.jl`
- `test/smoke_1D_v6.jl`

Purpose:

- Add a rear/adaptor thermal mass to slow deep cooling/heating dynamics without making the whole receiver uniformly heavier.
- Keep `h_front` fixed to a literature-style natural-convection value rather than fitting it.
- Keep optical/input-energy parameters fixed so incoming energy can be adjusted manually.

Main characteristics:

- Adds one dynamic state:

```text
state = [T_s,1 ... T_s,N, T_rear]
```

- Rear mass coupling:

```text
Q_rear_coupling = K_rear (T_s,N - T_rear)
Q_rear_mass_loss = U_rear_mass (T_rear - T_amb)
```

- Fixed:

```text
h_front = 10 W/m2/K
eps_front = 0.95
eta_abs = 0.80
beta_opt = 50.0
front_dep = 0.50
B_Re = 1.0
C_Pr = 0.5
```

- Heating-stage fitted parameters:

```text
A_Nu, h_floor, L_h
```

- Cooling-stage fitted parameters:

```text
gamma_C, k_scale, U_side, U_rear,
tau_T3, C_rear_scale, K_rear, U_rear_mass
```

- Calibration sequence:

```text
1. heating heat-transfer calibration
2. cooling thermophysical/rear-mass calibration
3. heating heat-transfer refit
```

Intended lesson:

- If deep cooling is too fast, the model needs extra thermal storage near the rear/deeper region, not just a larger uniform `gamma_C`.
- If the rear mass improves cooling tails without ruining front behavior, the missing physics is likely adaptor/holder/insulation/outlet thermal inertia.

## Current Working Hypotheses

1. The front/depth temperature reversal, such as `T9 > T8` at high irradiance, likely needs correct optical deposition and front-region losses, but optical fitting is currently held fixed to avoid hiding other structural errors.
2. The weak experimental flow dependence deeper in the receiver likely requires reduced downstream effective gas-solid exchange or inactive/bypassed flow.
3. The paired thermocouple averages are more appropriate for 1D fitting than using only one radial thermocouple at each axial station.
4. Excessively fast cooling in `T10_pair` and `T3` points toward missing rear/deep thermal mass.
5. Overly flexible calibration can hide missing physics by pushing loss and optical parameters to bounds.

## Recommended Next Checks

1. Run quick then full `1D_v7` calibration.
2. Compare cooling-only profiles for C69/C80/C81 before judging heating.
3. Inspect whether `k_ins_scale`, `gamma_C`, and `k_scale` stay inside credible ranges.
4. Check whether measured `T2` as a boundary removes the unphysical v6 rear-mass warm-up during cooling.
5. If v7 improves cooling but heating still misses `T8/T9_pair`, assess whether the fitted irradiance factors are physically credible before changing heat-transfer physics.

## v7 - T2 Boundary and Irradiance-Level Factors

Files:

- `1D_v7.jl`
- `run_1D_v7.jl`
- `test/smoke_1D_v7.jl`

Purpose:

- Replace the free v6 rear thermal mass with measured insulation temperature `T2(t)` as a boundary condition.
- Bring the old 0D_v3 cavity geometry back in as a resistance from the receiver/adaptor region to the measured `T2` location.
- Move the three irradiance-level correction factors out of the data importer and into the calibration vector.
- Add timing/shape terms to the objective so the optimizer cannot hide poor transient response behind acceptable steady-state errors.

Main characteristics:

- No rear/adaptor dynamic state is solved.
- `T2` is treated as a measured boundary about 40 mm radially outside the equivalent receiver wall.
- Receiver cells lose radially through a geometry-based felt conductance:

```text
Q_side,T2 = k_ins_scale G_receiver_to_T2 dx (T_s,i - T2(t))
```

- The receiver rear/adaptor path is represented by a fixed contact-plus-felt resistance to `T2(t)`:

```text
Q_rear,T2 = k_ins_scale G_rear_to_T2 (T_s,N - T2(t))
```

- Heating-stage fitted parameters:

```text
A_Nu, h_floor, L_h, f_I_high, f_I_mid, f_I_low
```

- Cooling-stage fitted parameters:

```text
gamma_C, k_scale, k_ins_scale, tau_T3
```

- Removed from calibration:

```text
U_side, U_rear, C_rear_scale, K_rear, U_rear_mass
```

Quick-run notes:

- The quick v7 runner completed successfully.
- The smoke test passed with Julia 1.12.6.
- The quick-fit cooling behavior no longer shows the v6 free rear-mass warm-up problem.
- Heating remains imperfect in the quick run, especially front/depth steady-state bias, so the full v7 calibration is needed before judging the irradiance factors.

Full-run assessment:

- Latest saved v7 parameters:

```text
gamma_C = 2.72
k_scale = 1.17
k_ins_scale = 0.97
A_Nu = 1.73
h_floor = 0.257
L_h = 0.050 m
tau_T3 = 20.4 s
f_I_high = 1.29
f_I_mid = 1.37
f_I_low = 0.945
```

- The insulation/T2 conductance scale stayed credible, close to unity. This supports using measured `T2(t)` as a boundary rather than a free hidden rear mass.
- The fitted irradiance factors are large for the high and middle irradiance groups. Because the effective absorbed factor is `eta_abs * f_I`, the high and middle groups correspond to roughly `1.03` and `1.09` of nominal irradiance after `eta_abs = 0.80`. That is plausible only if the nominal aperture irradiance or projected area is low; it should not be interpreted as pure absorptivity.
- Compared with v6, cooling no longer suffers from the artificial rear-mass warm-up. Cooling `T8` and `T9_pair` are reasonable, with mean RMSE near 22-25 K.
- Cooling `T10_pair` and `T3` still show early model overshoot and then too-fast decay/tail underprediction. This suggests either excessive axial redistribution from the hot front during early cooling, an incompatible cooling initial profile, or missing outlet/adaptor thermal inertia.
- Heating remains the main failure. Mean steady errors are about:

```text
T8       model - experiment = -114 K
T9_pair  model - experiment = -110 K
T10_pair model - experiment =   -1 K
T3       model - experiment =  +54 K
```

- This sign pattern is important: the front/mid solid is too cold while the outlet gas/T3 can be too hot. A bulk irradiance correction alone cannot resolve that contradiction.
- Worst heating deviations are concentrated in E67-E76. Later cases E77-E81 are substantially better, so the model error has a case/flow/irradiance dependence rather than a uniform offset.

Interpretation:

- v7 is a better structural model than v6 for cooling because it removes the unphysical free rear energy reservoir.
- v7 does not yet solve the key heating physics. The remaining mismatch points toward the coupling between gas heat removal, axial exchange distribution, and what the T3 thermocouple actually measures.
- The combination "solid too cold, gas too hot" is consistent with too much modeled gas-solid energy extraction from the hot solid path and/or too little bypass/mixing before the T3 measurement.

## Revalidation of Gemini Comments Across v3-v7

Accepted and implemented in v6/v7:

- Heat-shape regularization now belongs in the heating objective when `h_floor` and `L_h` are fitted there. This was a real issue for v6. v7 also uses the corrected placement.
- Calibration loss functions in v6/v7 now catch ODE solve failures and return `Inf` rather than crashing a full optimizer run.
- Zero-flow gas diagnostics in v6/v7 now set the internal gas outlet profile to local solid temperatures while keeping `Qgas = 0`, avoiding an artificial instant drop to inlet temperature.

Still present in older reference versions:

- v3-v5 still hard-error on failed solves in the loss path. If those versions are used again for calibration, port the v6/v7 loss `try/catch` pattern back.
- v3-v5 still set stagnant gas to inlet temperature at zero flow. If old cooling comparisons are revisited, remember that their T3 behavior is biased by this assumption.
- v5 regularizes `h_floor` and `L_h` during cooling because v5 actually fits those parameters in the cooling set. That is not the same logical bug as v6, but the modeling premise is now superseded by the v6/v7 staged strategy.

Deferred:

- Boundary-face temperatures remain approximated by first/last cell-center temperatures in all finite-volume versions. This is a real mesh-sensitivity caveat, but a proper fix requires a boundary face energy balance and should be handled as a dedicated model revision.
- Calibration functions still mutate global `pnew_v*` values for interactive REPL convenience. For batch reproducibility, add a no-mutation wrapper later rather than changing the current workflow abruptly.

## Proposed v8 Strategy

1. Regenerate v7 results after the latest robustness fixes.

The saved full-run v7 plots were generated before the most recent regularization/zero-flow/loss-handling edits. The expected qualitative behavior should be similar, but the official comparison should be rerun before making publication-grade conclusions.

2. Add a gas bypass / active-flow fraction model.

The current v7 failure pattern suggests that all measured flow is being forced to exchange too strongly with the hot receiver path. Introduce a heating-stage parameter such as `f_active`:

```text
mdot_active = f_active mdot_total
Qgas = mdot_active cp (Tg_active,out - Tin)
T3_true = Tin + f_active (Tg_active,out - Tin)
```

This can raise solid temperatures by reducing total heat removal while also lowering the mixed outlet temperature seen by T3. That is the right direction for the observed `T8/T9_pair` underprediction and `T3` overprediction.

3. Add a cooling overshoot/monotonicity penalty.

The early cooling bump in `T10_pair` and `T3` is not acceptable if the experiment decays monotonically. Add a targeted cooling penalty:

```text
penalty += mean(max(model(t) - model(t0), 0)^2)
```

for cooling sensors that should not warm after shutdown. This should help the optimizer reject excessive early axial heat redistribution.

4. Revisit axial conductivity after adding the cooling penalty.

If the bump remains, reduce or restructure the axial conduction path. The fitted `k_scale` near 1.17 is not extreme, but the finite-volume model may still over-connect hot front material to deeper thermocouple locations because the real porous/channeled structure has weaker effective axial conduction than dense SiC.

5. Replace the simple T3 lag with a physical outlet/sensor node.

The fitted `tau_T3` is now small, yet T3 timing errors remain large. A single first-order lag is probably not representing the outlet tube, adaptor, and thermocouple environment. A small outlet gas/metal measurement node connected to `T2(t)` or ambient would be more physical than continuing to tune `tau_T3`.

6. Keep optical shape as a controlled sensitivity, not a broad fit.

Do not reopen unrestricted optical fitting yet. Instead, run a small manual/sensitivity grid over `front_dep` and `beta_opt` after the gas-bypass test. The question is whether front-weighted deposition can improve `T8/T9_pair` without worsening `T10_pair` and T3.

7. Judge versions by signed residual structure, not objective alone.

For v8 comparisons, always report:

```text
mean RMSE
mean signed steady error
t90 error
cooling overshoot
T8/T9_pair intersection behavior versus flow
```

The present objective can improve while preserving the wrong physical sign pattern.

## v8/v8b Result Assessment and Next Revision Plan

Files:

- `1D_v8.jl`
- `run_1D_v8.jl`
- `test/smoke_1D_v8.jl`
- `1D_v8b.jl`
- `run_1D_v8b.jl`
- `test/smoke_1D_v8b.jl`
- `summaries/1D_v8_theory_manual.md`

Purpose:

- v8a added an explicit rear alumina tube and water-cooled flange while keeping measured T2 as a boundary.
- v8b kept the rear tube/flange domain but replaced measured T2 with one predicted lumped cavity state.
- v8b compares T3 to the gas temperature at 140 mm and compares T2 directly to the predicted cavity temperature.

Saved v8b parameters:

```text
gamma_C = 1.50
k_scale = 1.28
k_ins_scale = 0.856
A_Nu = 1.956
h_floor = 0.251
L_h = 0.050 m
tau_T3 = 20.0  # retained but unused in v8b
f_I_high = 1.600
f_I_mid = 1.600
f_I_low = 1.185
```

Key result:

- v8b is a substantial improvement over v7/v8a in heating solid temperatures and T2 prediction.
- T2 is now predicted well with a single lumped cavity state:

```text
heating T2 mean RMSE = 4.3 K, mean steady error = -5.5 K
cooling T2 mean RMSE = 2.3 K, mean steady error =  3.3 K
```

- Cooling remains reasonable:

```text
cooling T8       RMSE 14.9 K
cooling T9_pair  RMSE 29.1 K
cooling T10_pair RMSE 19.0 K
cooling T3       RMSE 26.5 K
```

- Heating is improved but still structurally biased:

```text
heating T8       mean steady error = -62 K
heating T9_pair  mean steady error = -83 K
heating T10_pair mean steady error = -13 K
heating T3       mean steady error = +11 K
```

Irradiance-band residuals show the remaining problem more clearly:

```text
high-Io: T8 -86 K, T9_pair -83 K, T10_pair +13 K, T3 +46 K
mid-Io:  T8 -100 K, T9_pair -116 K, T10_pair -39 K, T3 -11 K
low-Io:  T8  ~0 K, T9_pair -51 K, T10_pair -14 K, T3  -1 K
```

Interpretation:

- The cavity/T2 problem is largely resolved by v8b. The single lumped cavity is good enough to proceed without a radial mesh.
- The remaining heating error is not primarily a missing-cavity-mass problem.
- High and mid irradiance factors hit the upper bound (`f_I_high ~= f_I_mid ~= 1.6`), while the front/mid receiver solid is still too cold. This points to an input-energy/source-distribution or gas-removal structure issue rather than insufficient optimizer effort.
- The low-irradiance cases are much closer, so the next revision should target irradiance-dependent behavior, not a uniform offset.
- Heat-flow diagnostics show the cooled flange is a strong sink and the rear tube exit wall stays near the water temperature, but T3 is sampled at 140 mm, so the main T3 comparison is upstream of most flange cooling.

Recommended next revision:

1. Freeze the v8b cavity structure.

Do not add a radial mesh next. Keep the lumped cavity/T2 model and use it as the baseline because it predicts T2 well.

2. Address the irradiance/source limit.

The fitted `f_I_high` and `f_I_mid` are at the upper bound. Before adding more thermal states, test whether the fixed optical assumptions are too restrictive:

```text
eta_abs = 0.80
beta_opt = 50 1/m
front_dep = 0.50
```

The next revision should allow a controlled optical/source sensitivity, not a broad unconstrained fit. Candidate tests:

```text
front_dep sweep
beta_opt sweep
eta_abs or aperture-power scale by irradiance group
```

3. Add an active-flow or gas-bypass model after the source check.

If increasing/reshaping energy input fixes solid temperatures but breaks T3, add a mixed outlet model:

```text
mdot_active = f_active mdot_total
T3_mix = Tin + f_active (Tg_active(140 mm) - Tin)
```

This can decouple solid heat removal from the measured gas temperature.

4. Keep T3 at 140 mm for now.

The v8b T3 definition is more physical than a fitted lag. Do not reintroduce `tau_T3` until the source and active-flow questions have been tested.

5. Add a source/flow diagnostic runner.

Before creating a large new version, make a small diagnostic script that runs v8b over a grid:

```text
front_dep
beta_opt
f_active
```

and reports signed steady residuals by irradiance group. This is faster and more informative than immediately adding more fitted states.

6. Judge the next revision by signed residuals.

The target should be:

```text
T2 remains within ~5-10 K
high/mid T8 and T9_pair underprediction is reduced
T10_pair does not flip into large overprediction
T3 at 140 mm remains within a plausible signed residual
```

### 2026-07-20 - v8b axial profile plots after widened irradiance bounds

Added steady axial `T vs Length` plots to `run_1D_v8b.jl` for every heating and cooling experiment. Each profile overlays:

- continuous receiver/tube wall temperature,
- continuous gas temperature,
- experimental solid points at T8/T9/T10 axial locations,
- experimental T3 at 140 mm,
- T2 experiment/model markers at the lumped cavity location for context.

The plotting pass was run with `RECEIVER1D_V8B_RUNNER_REUSE_PARAMETERS=true`, so it reused the saved full-calibration parameters instead of refitting. The generated files are:

```text
summaries/1D_v8b/plots/axial_profile_*_1D_v8b.png
```

Saved parameters with the updated bounds:

```text
gamma_C     = 1.5021
k_scale     = 1.2804
k_ins_scale = 0.8562
A_Nu        = 2.3628
h_floor     = 0.2506
L_h         = 0.0500 m
f_I_high    = 1.6131
f_I_mid     = 1.7809
f_I_low     = 1.1849
```

Aggregate residuals from the regenerated metrics:

```text
heating: T8 RMSE 63 K, steady -44 K
heating: T9_pair RMSE 72 K, steady -65 K
heating: T10_pair RMSE 59 K, steady +2 K
heating: T3 RMSE 65 K, steady +25 K
heating: T2 RMSE 4 K, steady -5 K

cooling: T8 RMSE 15 K, steady +9 K
cooling: T9_pair RMSE 29 K, steady +8 K
cooling: T10_pair RMSE 19 K, steady +13 K
cooling: T3 RMSE 27 K, steady +21 K
cooling: T2 RMSE 2 K, steady +3 K
```

Heating steady residuals by irradiance group:

```text
high-Io: T8 -81 K, T9_pair -78 K, T10_pair +18 K, T3 +50 K, T2 -10 K
mid-Io:  T8 -50 K, T9_pair -67 K, T10_pair  +3 K, T3 +26 K, T2  -4 K
low-Io:  T8  ~0 K, T9_pair -51 K, T10_pair -14 K, T3  -1 K, T2  ~0 K
```

Interpretation:

- The widened high/mid irradiance factors improved the front/mid solid underprediction relative to the earlier 1.6-limited run.
- The price is a warmer gas prediction, especially T3 in high and mid irradiance cases.
- T2 remains well captured, so the lumped cavity thermal mass should stay frozen for the next revision.
- The axial profiles make the remaining structural issue clearer: the model can lift the receiver profile, but it then tends to overheat the sampled gas in low-flow/mid-irradiance cases such as E76 and E71.

Next revision recommendation:

Keep the v8b cavity/rear-tube structure and test a source-shape or active-flow/mixing diagnostic before adding more thermal masses. The next useful diagnostic is a small grid over front deposition / beta_opt / active gas fraction, judged by signed residuals and the new axial profile plots.

# 2026-07-20 - v9a heat-transfer-law revision

Implemented `1D_v9a.jl`, `run_1D_v9a.jl`, and `test/smoke_1D_v9a.jl`.

v9a keeps the v8b rear tube, water-cooled flange, and predicted lumped cavity/T2 state, but changes the receiver gas-solid heat-transfer model:

```text
Nu(z) = max(Nu_fd, 10^A_Nu Re^B_Re Pr^C_Pr (Dh / z)^1/3)
Nu_fd = 3.61
UA = h f_exchange P_exchange dx
```

Removed from the fitted vector:

```text
gamma_C
tau_T3
h_floor
L_h
```

Fitted vector:

```text
k_scale
k_ins_scale
A_Nu
B_Re
C_Pr
f_exchange
f_I_high
f_I_mid
f_I_low
```

Full v9a runner completed and generated:

```text
summaries/1D_v9a/parameters_1D_v9a.csv
summaries/1D_v9a/analysis_results_1D_v9a.csv
summaries/1D_v9a/plots/
```

Full-run fitted parameters:

```text
k_scale     = 1.3716
k_ins_scale = 0.7982
A_Nu        = 0.2372
B_Re        = 0.7109
C_Pr        = 0.3920
f_exchange  = 0.9908
f_I_high    = 1.5735
f_I_mid     = 1.7319
f_I_low     = 1.1516
```

Aggregate residuals:

```text
heating: T8       RMSE 72.9 K, steady -54.7 K
heating: T9_pair  RMSE 83.3 K, steady -74.8 K
heating: T10_pair RMSE 68.7 K, steady  -2.2 K
heating: T3       RMSE 74.9 K, steady +25.2 K
heating: T2       RMSE  4.5 K, steady  -6.4 K

cooling: T8       RMSE 28.1 K, steady  +8.4 K
cooling: T9_pair  RMSE 49.0 K, steady  +7.3 K
cooling: T10_pair RMSE 33.6 K, steady +11.2 K
cooling: T3       RMSE 27.5 K, steady +17.7 K
cooling: T2       RMSE  2.3 K, steady  +2.6 K
```

Heating steady residuals by irradiance group:

```text
high-Io: T8 -93 K, T9_pair -90 K, T10_pair +11 K, T3 +50 K, T2 -12 K
mid-Io:  T8 -62 K, T9_pair -77 K, T10_pair  -1 K, T3 +27 K, T2  -6 K
low-Io:  T8  -9 K, T9_pair -58 K, T10_pair -16 K, T3  -1 K, T2  -1 K
```

Interpretation:

- v9a is a useful negative result. It did not improve the central signed residual pattern relative to the latest v8b run.
- `f_exchange` calibrated near unity, so the optimizer did not use reduced exchange area to improve the fit in this formulation.
- `B_Re = 0.71` and `C_Pr = 0.39` are physically more reasonable than the old fixed `B_Re = 1.0`, but the axial profiles remain similar to v8b.
- Removing `gamma_C` did not break T2, but cooling T8/T9/T10 worsened relative to v8b, so the missing transient storage/coupling question may still matter for cooling even if it is not the main heating error.
- The remaining heating problem is increasingly likely to involve source distribution and/or active flow/T3 mixing rather than only the receiver heat-transfer exponent.

Recommended next move:

Do not continue tuning v9a heat-transfer bounds alone. The next informative test should be one of:

```text
v9b-source: fit/sweep front_dep and beta_opt with v9a or v8b baseline
v9b-flow: replace f_exchange with an active-flow/mixed-T3 model
```

Because v9a did not use `f_exchange < 1`, the source-distribution test should probably come first unless independent flow evidence supports bypass/mixing.

# 2026-07-21 — Gemini assessment of v9a results

#### What v9a tested and what it decided

v9a implemented the heat-transfer-law revision from the combined Gemini+GPT plan in `1D_v8b_geminicomments.md`:

- **Replaced** the empirical `Nu = 10^A_Nu × Re^1.0 × Pr^3.16 × s(z)` with a developing laminar flow form: `Nu(z) = max(3.61, 10^A_Nu × Re^B_Re × Pr^C_Pr × (Dh/z)^(1/3))`
- **Freed** `B_Re` and `C_Pr` as fitted parameters (previously fixed at 1.0 and 0.5 via the `10^C_Pr` encoding)
- **Added** `f_exchange` (exchange-area fraction) multiplying `P_exchange` in the UA term
- **Removed** `gamma_C`, `tau_T3`, `h_floor`, `L_h`

The code exactly followed the GPT plan rather than the Gemini counter-proposal on `gamma_C` (Gemini recommended keeping it; v9a removed it).

#### Validation of predictions from the Gemini review

**1. Effectiveness saturation — confirmed by `f_exchange ≈ 1.0`.**

The Gemini review predicted that with v8b parameters (NTU ≈ 380 per cell), the gas-solid exchange was fully saturated and that changing the Re exponent alone would have zero effect. It also predicted that with classical laminar Nu ≈ 3.6, the per-cell NTU would drop to ~2.1 and ε to ~0.88.

v9a results:
- `A_Nu = 0.237`, so the Nusselt pre-factor is `10^0.237 ≈ 1.72`, and the developing-flow term at the inlet (where `(Dh/z)^(1/3)` is large) gives Nu ≈ 10–20 in the first few cells, decaying toward `Nu_fd = 3.61` downstream. This is consistent with the predicted ε ≈ 0.88 regime for most of the receiver.
- `f_exchange = 0.991` — the optimizer did not reduce the exchange area. This means even at the lower NTU ≈ 2 regime, the gas-solid coupling is still too strong to explain the residual pattern. The remaining error is not about how much exchange happens per cell — it is about what the gas encounters (source distribution) or how much gas actually enters the channels (active flow).
- `B_Re = 0.71` — physically plausible for a developing laminar/transitional regime, but the parameter had little leverage because the NTU is above the threshold where ε sensitivity saturates (for NTU > 2, ε changes slowly).

**2. Removing `gamma_C` worsened cooling — as predicted.**

| Sensor | v8b cooling RMSE | v9a cooling RMSE | Change |
|---|---|---|---|
| T8 | 15 K | 28 K | +87% |
| T9_pair | 29 K | 49 K | +69% |
| T10_pair | 19 K | 34 K | +77% |
| T3 | 27 K | 28 K | ~same |
| T2 | 2.3 K | 2.3 K | ~same |

The Gemini review warned:

> Removing gamma_C simultaneously with changing the Nusselt law makes it impossible to separate the effects. If the new Nusselt law changes the steady-state balance, gamma_C will need to re-adjust the transient response.

This is exactly what happened. The receiver's thermal response became too fast without the capacity multiplier, producing larger cooling errors. The v8b `gamma_C = 1.50` was not just a fudge — it was absorbing real uncertainty in effective thermal mass (adaptor, channel geometry, contact interfaces).

**3. Heating residuals did not improve — as predicted for a Nusselt-only change.**

v9a heating signed residuals are essentially indistinguishable from v8b:

| Group | v8b T8 | v9a T8 | v8b T9_pair | v9a T9_pair | v8b T3 | v9a T3 |
|---|---|---|---|---|---|---|
| high | -81 K | -93 K | -78 K | -90 K | +50 K | +50 K |
| mid | -50 K | -62 K | -67 K | -77 K | +26 K | +27 K |
| low | ~0 K | -9 K | -51 K | -58 K | -1 K | -1 K |

v9a is actually slightly **worse** in every irradiance group for solid temperatures, while T3 is unchanged. The irradiance factors remain near bounds: `f_I_high = 1.57`, `f_I_mid = 1.73`.

This confirms that the Nusselt-law revision was a necessary structural correction (the exponents are now physical), but it was not the cause of the steady-state bias.

#### What this means for the model

The v9a result narrows the diagnosis significantly:

1. **The gas-solid exchange law is no longer the primary suspect.** Whether the Re exponent is 1.0 or 0.71, whether the Nusselt is 2566 or 5, the gas outlet temperature and solid temperatures barely change. The model is in a regime where gas-solid coupling is strong enough to equilibrate within a few cells, and the residual pattern is set by the overall energy balance, not the exchange rate.

2. **The irradiance factors at bounds point to insufficient energy input.** The optimizer wants `f_I × eta_abs > 1.0` for high/mid irradiance, which means either:
   - The nominal irradiance is underestimated
   - `eta_abs = 0.80` is too low
   - The actual absorbed power is higher than `0.80 × Io × A_frt`
   - Or: the energy is correctly absorbed but distributed too far forward (front_dep = 0.50), causing excessive front radiative loss that reduces the *net* energy available to heat the deeper solid

3. **The source distribution is the remaining structural lever.** With `front_dep = 0.50`, roughly half the absorbed solar power goes into the first 5.5 mm cell. At high irradiance, this front cell reaches ~1100–1200 K, producing a front radiative loss of `0.95 × 5.67e-8 × 361e-6 × (1200⁴ - 295⁴) ≈ 30 W`. With a total absorbed power of ~60–80 W, this is a 40–50% front loss fraction, which is enormous. Reducing `front_dep` to 0.15 and/or lowering `beta_opt` to 20 1/m would spread the absorption deeper, reduce the front peak temperature, reduce front radiative loss, and make more net energy available for heating the mid/rear receiver — without changing the gas exchange at all.

#### Concrete recommendation for v9b

```text
v9b should be: v9a + gamma_C restored + front_dep/beta_opt freed

Fitted vector (11 parameters):
  p[1]  gamma_C        [0.50, 3.00]
  p[2]  k_scale        [0.20, 3.00]
  p[3]  k_ins_scale    [0.25, 4.00]
  p[4]  A_Nu           [-2.00, 1.00]
  p[5]  B_Re           [0.00, 1.00]
  p[6]  C_Pr           [0.20, 0.50]
  p[7]  f_exchange     [0.05, 1.00]
  p[8]  front_dep      [0.05, 0.50]
  p[9]  beta_opt       [10.0, 100.0]
  p[10] f_I_high       [0.60, 2.00]
  p[11] f_I_mid        [0.60, 2.00]
  p[12] f_I_low        [0.60, 1.60]
```

Calibration stages:

```text
Stage 1 — Heating: fit A_Nu, B_Re, C_Pr, f_exchange, front_dep,
          beta_opt, f_I_high, f_I_mid, f_I_low (9 parameters)
Stage 2 — Cooling: fit gamma_C, k_scale, k_ins_scale (3 parameters)
Stage 3 — Heating refit (same 9 parameters)
```

Success criteria:

```text
1. front_dep calibrates below 0.30 (evidence for deeper absorption)
2. f_I factors move closer to 1.0 (less need for energy compensation)
3. T8/T9_pair steady error decreases in high/mid groups
4. Cooling RMSE returns to v8b levels with gamma_C restored
5. f_exchange either stays near 1.0 or drops — either result is
   informative for deciding whether v9c needs an active-flow model
```

Alternative if v9b-source doesn't close the gap: the next change should be an active-flow/mixed-T3 model (v9c-flow). But the source-distribution test should come first because it has independent physical motivation and doesn't add new conceptual complexity.

### axial profile thermocouple markers

Updated the axial profile plot helpers in `run_1D_v8b.jl` and `run_1D_v9a.jl` so the experimental solid markers show the raw thermocouple channels:

```text
T8 at z = 11 mm
T9 and T12 at z = 58 mm
T10 and T11 at z = 107 mm
```

The model curves are unchanged. The previous profile plots used the 1D comparison channels `T8`, `T9_pair`, and `T10_pair`; the new plots make the radial spread at the T9/T12 and T10/T11 axial stations explicit.

Regenerated both saved-parameter plot sets:

```text
summaries/1D_v8b/plots/axial_profile_*_1D_v8b.png
summaries/1D_v9a/plots/axial_profile_*_1D_v9a.png
```

### v10a comparison-first ECM baseline

Implemented `1D_v10.jl`, `run_1D_v10.jl`, and `test/smoke_1D_v10.jl` as the
first fundamental v10 revision. The model keeps the v8b rear tube/flange/cavity
domain, but removes the broad fitted scalars:

```text
removed/replaced: gamma_C, k_scale, k_ins_scale, h_floor, L_h, f_exchange,
                  tau_T3
kept/fitted-ready at first: C_Nu_model, ell_rad_m, eta_opt, front_dep, beta_opt
```

The v10 runner is intentionally comparison-first rather than calibration-first.
It reports:

```text
M1_no_rad:        fixed material/geometry properties, monolith Nu, no Rosseland
M2_rad_beta1000:  M1 + Rosseland axial radiation with beta_tr = 1000 1/m
M3_rad_beta300:   Rosseland sensitivity with beta_tr = 300 1/m
M3_rad_beta2700:  Rosseland sensitivity with beta_tr = 2700 1/m
```

Smoke test:

```text
test/smoke_1D_v10.jl: 44/44 passed
```

The first full matrix shows that the fixed-property v10a baseline is much too
cold during heating, especially in the receiver solid:

| variant | phase | T8 steady | T9_pair steady | T10_pair steady | T3 steady | T2 steady |
|---|---:|---:|---:|---:|---:|---:|
| M1_no_rad | heating | -219 K | -241 K | -156 K | -118 K | -11 K |
| M2_rad_beta1000 | heating | -219 K | -242 K | -156 K | -118 K | -11 K |
| M3_rad_beta300 | heating | -219 K | -242 K | -157 K | -118 K | -11 K |
| M3_rad_beta2700 | heating | -219 K | -241 K | -156 K | -118 K | -11 K |

Cooling remains moderate, because those runs start near the measured hot state:

| variant | phase | T8 steady | T9_pair steady | T10_pair steady | T3 steady | T2 steady |
|---|---:|---:|---:|---:|---:|---:|
| M1_no_rad | cooling | +3 K | +2 K | +7 K | +15 K | +3 K |
| M2/M3 | cooling | nearly unchanged | nearly unchanged | nearly unchanged | nearly unchanged | nearly unchanged |

Interpretation:

1. Rosseland axial diffusion across `beta_tr = 300-2700 1/m` barely changes the
   result. At this magnitude it is not the missing dominant mechanism.
2. Removing `k_scale`, `k_ins_scale`, and `f_exchange` exposes a large heating
   energy/source deficit. The model is not getting enough useful absorbed heat
   into the receiver interior.
3. Because T2 stays close while the receiver is cold, the next v10 step should
   not be a hidden cavity/insulation scalar. The problem is upstream of cavity
   heat loss: absorbed-power level, source distribution, or inlet/front thermal
   boundary physics.

Recommended v10b:

```text
keep eta_opt = 1.0 fixed
keep front_dep = 1.0 fixed
add Nu A/B/C coefficients
add low-dimensional axial Rosseland parameters
defer per-irradiance absorbed-power calibration to a later stage
```

Follow-up change after discussion:

```text
eta_opt   fixed to 1.0
front_dep fixed to 1.0
```

This removes global absorbed-power compensation and forces all shortwave source
into the front cell under the current `solar_weights_v5` convention. Since
`front_dep = 1.0` makes `beta_opt` inactive, the current v10 fit-ready
heat-transfer index list was reduced to `C_Nu_model` only; Rosseland remains a
separate radiation index.

Regenerated `summaries/1D_v10` after this change. The smoke test now passes
`49/49` checks. Mean heating signed steady errors for `M1_no_rad` became:

| sensor | steady error |
|---|---:|
| T8 | -192 K |
| T9_pair | -248 K |
| T10_pair | -164 K |
| T3 | -125 K |
| T2 | -11 K |

Relative to the previous `front_dep = 0.50` v10a matrix, T8 becomes less cold
but T9/T10/T3 become slightly colder. This is a useful source-location signal:
front-only deposition preferentially helps the front thermocouple and does not
solve the deeper receiver deficit.

### v10 energy-partition precheck

Implemented `run_1D_v10_energy_precheck.jl` and added optional named Rosseland
axial profiles to `1D_v10.jl`:

```text
:uniform
:front_strong
:rear_strong
:weak_gradient
```

The default model behavior remains uniform Rosseland unless a profile is
explicitly requested. The smoke test now passes `52/52` checks.

Precheck outputs:

```text
summaries/1D_v10/precheck/energy_partition_summary_1D_v10.csv
summaries/1D_v10/precheck/axial/axial_terms_*_1D_v10.csv
summaries/1D_v10/precheck/plots/energy_partition_*_1D_v10.png
summaries/1D_v10/precheck/plots/axial_heat_terms_*_1D_v10.png
```

Representative cases checked:

```text
E67, E76, E80, C80
```

Mean heating energy partition over E67/E76/E80:

| variant | solar | front loss | gas gain receiver | receiver-cavity | receiver-tube | flange | receiver residual |
|---|---:|---:|---:|---:|---:|---:|---:|
| no_rad | 97.8 W | 7.8 W | 44.3 W | 26.6 W | 18.6 W | 20.6 W | 0.5 W |
| uniform_beta1000 | 97.8 W | 7.8 W | 44.4 W | 26.5 W | 18.6 W | 20.6 W | 0.5 W |
| front_strong | 97.8 W | 7.8 W | 44.4 W | 26.5 W | 18.6 W | 20.6 W | 0.5 W |
| rear_strong | 97.8 W | 7.8 W | 44.4 W | 26.5 W | 18.6 W | 20.6 W | 0.5 W |

Mean heating steady errors remain essentially unchanged by spatial Rosseland:

| variant | T8 | T9_pair | T10_pair | T3 | T2 |
|---|---:|---:|---:|---:|---:|
| no_rad | -181 K | -235 K | -150 K | -113 K | -10 K |
| uniform_beta1000 | -181 K | -235 K | -150 K | -113 K | -10 K |
| front_strong | -182 K | -235 K | -150 K | -112 K | -10 K |
| rear_strong | -181 K | -235 K | -150 K | -112 K | -10 K |

For E76, axial solid conduction is roughly 20 W through most of the receiver,
while Rosseland face flux remains very small:

| profile | max \|Qrad\| | mean \|Qrad\| | max \|Qcond\| | mean \|Qcond\| |
|---|---:|---:|---:|---:|
| uniform_beta1000 | 0.17 W | 0.08 W | 23.6 W | 20.3 W |
| front_strong | 0.56 W | 0.19 W | 23.3 W | 20.2 W |
| rear_strong | 0.15 W | 0.12 W | 23.7 W | 20.3 W |

Interpretation:

1. The energy balance is numerically consistent: the receiver residual is only
   about 0.5 W at the final sample.
2. Spatial Rosseland profiles in the current physically bounded range do not
   materially change the energy partition or the deeper receiver/T3 residuals.
3. With `front_dep = 1.0`, all shortwave source is assigned to the first cell.
   The first cell passes a large amount to the gas; downstream cells then
   mostly see hot-gas exchange, side loss to the cavity, and solid conduction.
4. The next fit should not rely on Rosseland as the main missing mechanism.

Recommended next revision:

```text
v10b-Nu:
  keep eta_opt = 1.0
  keep front_dep = 1.0
  keep Rosseland off or fixed as a diagnostic
  replace C_Nu_model with A_Nu, B_Re, C_Pr
  report the same energy partition after fitting
```

If Nu A/B/C cannot reduce the deeper receiver and T3 residuals without
unphysical coefficients, then the following stage should revisit absorbed-power
level by irradiance group, including the fixed `ETA_ABS_FIXED_v10 = 0.8`
assumption.

### experimental T8 versus T9/T12 flow crossover

Added `run_1D_sensor_pattern.jl` to summarize the raw steady thermocouple
ordering and produce flow-trend plots:

```text
summaries/1D_sensor_pattern/steady_thermocouple_pattern.csv
summaries/1D_sensor_pattern/plots/temperature_vs_flow_I*_1D.png
summaries/1D_sensor_pattern/plots/internal_minus_front_I*_1D.png
```

The diagnostic confirms the user's observation and rules out a simple "T9 direct
light" explanation, because side thermocouple T12 follows the same pattern.
The diagnostic was then expanded to include raw T10/T11 and their pair average.

At high flow, the internal station is hotter than T8:

| irradiance | high-flow case | flow | T9_pair - T8 | T12 - T8 |
|---:|---|---:|---:|---:|
| 456 kW/m2 | E67 | 15.28 L/min | +47 K | +60 K |
| 304 kW/m2 | E72 | 18.71 L/min | +47 K | +60 K |
| 256 kW/m2 | E77 | 13.85 L/min | +47 K | +58 K |

As flow is reduced, the sign reverses:

| irradiance | low-flow case | flow | T9_pair - T8 | T12 - T8 |
|---:|---|---:|---:|---:|
| 456 kW/m2 | E71 | 7.13 L/min | -96 K | -83 K |
| 304 kW/m2 | E76 | 4.53 L/min | -109 K | -97 K |
| 256 kW/m2 | E81 | 4.53 L/min | -21 K | -10 K |

The deeper T10/T11 station behaves differently: `T10_pair - T8` is negative in
all heating cases and becomes more negative as flow is reduced:

| irradiance | high-flow case | low-flow case |
|---:|---:|---:|
| 456 kW/m2 | -80 K | -338 K |
| 304 kW/m2 | -29 K | -329 K |
| 256 kW/m2 | -6 K | -152 K |

Interpretation:

1. The sign reversal is flow-coupled. A purely flow-independent ETC or
   Rosseland correction would not naturally create the observed crossover by
   itself.
2. A more plausible primary mechanism is **front-region inlet cooling plus gas
   preheating**. At high flow, the front section sees the coldest gas and a
   high entry/developing-flow heat-transfer coefficient, so T8 can be depressed
   while the downstream gas is already warmer and removes heat less aggressively
   from the T9/T12 station. At low flow, the gas quickly heats and the total heat
   removal capacity drops, so the front source dominates and T8 becomes the
   hottest solid station.
3. This does not eliminate ETC/Rosseland. It means ETC/Rosseland should be
   secondary: it can smooth or redistribute the axial field, but the crossover
   itself is a diagnostic for flow-dependent gas/solid exchange and inlet
   thermal development.
4. The current v10a Nu closure is too weakly structured to test this properly
   because it is effectively near `Nu_fd = 3.61` over much of the receiver.

Recommended next check:

```text
v10b-flow-shape:
  keep eta_opt = 1.0
  keep front_dep = 1.0
  keep Rosseland off/fixed
  replace C_Nu_model with A_Nu, B_Re, C_Pr
  retain an axial entry/developing-flow dependence
  add diagnostics for predicted T9_pair - T8, T12 - T8, and T10_pair - T8
```

If v10b-flow-shape cannot reproduce the high-flow internal-hot / low-flow
front-hot crossover with credible Nu coefficients, then test an explicit axial
ETC term as a diagnostic, but require it to be reported as an equivalent
conductivity and compared against the physical solid/Rosseland scales.

### front-only deposition reruns and flow-slope interpretation

The project constant was changed to:

```text
FRONT_DEPOSITION_FIXED_V5 = 1.0
```

Saved parameter files show that v8b, v9a, and v10 are now using
`front_dep = 1.0`. The older plain `summaries/1D_v8/parameters_1D_v8.csv`
still records `front_dep = 0.5`, so the directly comparable front-only reruns
are v8b, v9a, and v10.

Mean signed steady errors after the front-only reruns:

| model | T8 | T9_pair | T10_pair | T3 | T2 |
|---|---:|---:|---:|---:|---:|
| v8b | -7 K | -76 K | -10 K | +14 K | -5 K |
| v9a | -18 K | -83 K | -12 K | +17 K | -6 K |
| v10 M1 | -192 K | -248 K | -164 K | -125 K | -11 K |

Observations:

1. Front-only deposition strongly helps the T8 level in v8b/v9a, but it leaves a
   large T9/T12 cold bias. This is consistent with the experimental pattern:
   the front source can correct the front thermocouple without reproducing the
   internal hot ridge.
2. v8b remains better than v9a after the change. The more physical v9a Nu law
   still worsens T8/T9/T3 relative to the empirical v8b heat-transfer shape.
3. v10 remains much too cold because it removed the irradiance-level
   compensation and broad thermal scalars. This is not a Rosseland failure alone;
   it is also a source/input/closure deficit under the fixed `eta_abs = 0.8`,
   `eta_opt = 1.0` assumptions.

The sensor-pattern slopes versus flow provide an important constraint. Linear
slopes within each irradiance group are:

| irradiance | dT8/dflow | dT9_pair/dflow | dT10_pair/dflow |
|---:|---:|---:|---:|
| 456 kW/m2 | -34 K/(L/min) | -17 K/(L/min) | -3 K/(L/min) |
| 304 kW/m2 | -24 K/(L/min) | -13 K/(L/min) | -3 K/(L/min) |
| 256 kW/m2 | -21 K/(L/min) | -14 K/(L/min) | -6 K/(L/min) |

This agrees with the observation that the downstream temperatures become
flatter with flow rate, especially at T10/T11. However, it does **not** support
the statement that front heat transfer has less temperature impact. T8 is the
most flow-sensitive channel. The more consistent interpretation is:

```text
T8:
  cold inlet gas + entry/developing-flow exchange makes the front strongly
  flow-sensitive.

T9/T12:
  gas has already been heated, so the local driving temperature difference and
  flow sensitivity are reduced; this station forms the internal hot ridge at
  high flow.

T10/T11:
  source/available heat and gas-solid driving force are both weaker, so the
  temperature response to flow is much flatter and always below T8 in the
  current data.
```

Recommended correlation form for v10b:

```text
Nu(z) = Nu_fd + A * Re^B * Pr^C * F_entry(z, Re, Pr)

F_entry(z, Re, Pr) = (Dh / max(z, z0))^m
                     * exp(-z / L_entry)

or equivalently

Nu(z) = max(Nu_fd, A * Re^B * Pr^C * (1 + E * Gz(z)^m))
Gz(z) = Re * Pr * Dh / max(z, z0)
```

The fit should constrain `B > 0` so heat transfer increases with flow and should
report signed diagnostics for:

```text
T9_pair - T8
T12 - T8
T10_pair - T8
```

If this entry/developing-flow Nu shape cannot reproduce the crossover and the
flattening toward T10/T11, then add an explicit axial ETC diagnostic. But that
ETC should be required to report an equivalent conductivity and should not be
allowed to silently replace source/input calibration.

Actionable v10b implementation suggestions:

```text
1. Keep fixed during the first v10b-flow-shape test:
   eta_opt = 1.0
   front_dep = 1.0
   ETA_ABS_FIXED_v10 = 0.8

2. Keep Rosseland off or fixed as a secondary diagnostic:
   do not fit ell_rad/front-rear Rosseland in the first Nu test
   still report Qrad and equivalent k_rad if enabled

3. Replace the single C_Nu_model parameter with:
   A_Nu
   B_Re
   C_Pr

4. Add one explicit entry/development shape:
   either F_entry(z) = (Dh / max(z, z0))^m * exp(-z / L_entry)
   or     F_entry(z) = 1 + E * Gz(z)^m

5. Use physically bounded coefficients:
   B_Re > 0
   0.2 <= C_Pr <= 0.5
   Nu(z) should remain near credible square-channel/monolith ranges

6. Optimize against the normal sensor residuals plus signed ordering metrics:
   T9_pair - T8
   T12 - T8
   T10_pair - T8

7. Require the fitted model to reproduce the qualitative flow pattern:
   high flow: T9/T12 hotter than T8
   low flow:  T8 hotter than T9/T12
   all flows: T10/T11 generally below T8, with flatter flow response

8. Only if this fails, add an axial ETC diagnostic:
   report k_ETC, equivalent Q_ETC, and compare them with solid conduction and
   Rosseland flux before accepting it as physics.
```

## 2026-07-21 - v11 Stage 1 Flow-Shape Isolation

Implemented `1D_v11.jl` and `run_1D_v11.jl` from the final
`1D_v10_geminicomments.md` plan. This revision intentionally keeps the optical
and source terms locked:

```text
eta_abs = 0.8
eta_opt = 1.0
beta_opt = 50 1/m
front_dep = 1.0
Rosseland off during the Stage 1 fit
```

The fitted vector is now:

```text
p[1] A_Nu
p[2] B_Re
p[3] C_Pr
p[4] L_entry_m
p[5] ell_rad_m fixed diagnostic only
```

The receiver Nusselt model is:

```text
Nu(z) = Nu_fd + A_Nu * Re^B_Re * Pr^C_Pr
        * (Dh / max(z, z0))^(1/3)
        * exp(-z / L_entry)
```

New diagnostics are written for each representative axial profile:

```text
Nu(z)
h(z)
UA/(m cp)
cell effectiveness
Qgas(z)
solid and gas axial temperatures
```

Generated artifacts:

```text
summaries/1D_v11/analysis_results_all_variants_1D_v11.csv
summaries/1D_v11/steady_results_baseline_front_source_1D_v11.csv
summaries/1D_v11/steady_results_stage1_nu_fit_1D_v11.csv
summaries/1D_v11/diagnostics/cell_diagnostics_*_1D_v11.csv
summaries/1D_v11/plots/*.png
```

Validation:

```text
test/smoke_1D_v11.jl passed: 62/62 tests
run_1D_v11.jl completed
```

Fitted Stage 1 parameters:

| parameter | value |
|---|---:|
| A_Nu | 1.0996 |
| B_Re | 1.0000 |
| C_Pr | 0.3183 |
| L_entry_m | 0.0379 |
| ell_rad_m | 0.0010 |

The optimizer reached the upper bound for `B_Re`, but the objective barely
changed:

```text
baseline objective = 0.20282
fitted objective   = 0.20222
```

Mean signed steady errors for Stage 1 fitted v11:

| phase | sensor | mean steady error |
|---|---|---:|
| heating | T8 | -189 K |
| heating | T9_pair | -245 K |
| heating | T10_pair | -162 K |
| heating | T3 | -123 K |
| heating | T2 | -11 K |
| heating | T9_pair - T8 | -56 K |
| heating | T10_pair - T8 | +28 K |
| cooling | T8 | +3 K |
| cooling | T9_pair | +2 K |
| cooling | T10_pair | +7 K |
| cooling | T3 | +15 K |
| cooling | T2 | +3 K |

Representative fitted heating cases:

| case | observation | model | experiment |
|---|---|---:|---:|
| E67 high flow | T8 | 670 K | 932 K |
| E67 high flow | T9_pair - T8 | -33 K | +47 K |
| E67 high flow | T10_pair - T8 | -78 K | -80 K |
| E76 low flow | T8 | 825 K | 1068 K |
| E76 low flow | T9_pair - T8 | -107 K | -109 K |
| E80 mid/low flow | T8 | 676 K | 710 K |
| E80 mid/low flow | T9_pair - T8 | -64 K | +14 K |

Interpretation:

1. v11 Stage 1 fails to reproduce the high-flow internal hot ridge. Under the
   fixed front-source assumption, the fitted model still keeps T9 below T8 in
   high-flow cases where the experiment has T9/T12 hotter than T8.
2. The cell effectiveness diagnostic explains why the Nu shape has limited
   leverage. For the representative fitted cases, the front cell effectiveness
   is essentially 1.0 and the receiver mean effectiveness is already high:

```text
E67 mean effectiveness = 0.875
E76 mean effectiveness = 0.993
E80 mean effectiveness = 0.969
```

3. The sensitivity check confirms the failure is not simply an optimizer corner.
   Deliberately lowering the entry enhancement increased T8 by only about 2 K
   in E67 and made the loss worse, while T9 remained below T8.
4. Cooling remains fairly reasonable, especially T2. The major mismatch is in
   heating energy/source distribution, not in the rear/cavity thermal mass alone.

Recommended next revision:

```text
Proceed to Stage 2.

Release one source/input family, not several at once:
  preferred v12a: keep eta_abs = 0.8 and eta_opt = 1.0 fixed,
                  fit beta_opt and front_dep only.

Keep the v11 Nu diagnostic fields active so any source-distribution gain can be
checked against Nu(z), effectiveness, Qgas(z), and the signed ordering metrics.

Do not fit Rosseland/ETC yet. The v11 result first proves that the isolated
entry-Nu mechanism cannot create the hot T9/T12 ridge with front-only source.
```

## 2026-07-21 - v12a Stage 2 Source-Distribution Diagnostic

Motivation:

The user correctly questioned whether releasing source distribution can address
the observed pattern:

```text
front: flow dependent
back/rear: nearly flow independent
```

The answer going into v12a was: no, not directly. Stage 2 can only test whether
wrong axial placement of absorbed power is responsible for the internal levels
and hot ridge. It should not be expected to create the front/back flow-slope
structure by itself.

Implementation:

```text
1D_v12a.jl
run_1D_v12a.jl
test/smoke_1D_v12a.jl
```

v12a freezes the v11 fitted Nu shape:

| parameter | frozen value |
|---|---:|
| A_Nu | 1.0996 |
| B_Re | 1.0000 |
| C_Pr | 0.3183 |
| L_entry_m | 0.0379 |

and fits only:

```text
beta_opt
front_dep
```

with:

```text
eta_abs = 0.8 fixed
eta_opt = 1.0 fixed
Rosseland off during the fit
```

Validation:

```text
test/smoke_1D_v12a.jl passed: 65/65 tests
run_1D_v12a.jl completed
```

Generated artifacts:

```text
summaries/1D_v12a/analysis_results_all_variants_1D_v12a.csv
summaries/1D_v12a/steady_results_baseline_v11_nu_front_source_1D_v12a.csv
summaries/1D_v12a/steady_results_stage2_source_fit_1D_v12a.csv
summaries/1D_v12a/diagnostics/cell_diagnostics_*_1D_v12a.csv
summaries/1D_v12a/plots/*.png
```

Fitted source parameters:

| parameter | value |
|---|---:|
| beta_opt | 184.67 1/m |
| front_dep | 0.0 |
| first-cell source fraction | 0.637 |
| downstream source fraction | 0.363 |

This means the optimizer moved away from pure front deposition, but only to an
effective distribution where about 64% of the power is still in the first cell.

Objective change:

```text
baseline v11 Nu + front source = 0.20216
v12a source fit                = 0.20035
```

The improvement is small.

Mean signed steady errors:

| phase | sensor | baseline | v12a source fit |
|---|---|---:|---:|
| heating | T8 | -189 K | -192 K |
| heating | T9_pair | -245 K | -241 K |
| heating | T10_pair | -162 K | -159 K |
| heating | T3 | -123 K | -121 K |
| heating | T2 | -11 K | -11 K |
| heating | T9_pair - T8 | -56 K | -50 K |
| heating | T10_pair - T8 | +28 K | +33 K |

Representative source-fit cases:

| case | observation | model | experiment |
|---|---|---:|---:|
| E67 high flow | T8 | 664 K | 932 K |
| E67 high flow | T9_pair - T8 | -26 K | +47 K |
| E67 high flow | T10_pair - T8 | -71 K | -80 K |
| E76 low flow | T8 | 827 K | 1068 K |
| E76 low flow | T9_pair - T8 | -102 K | -109 K |
| E80 mid/low flow | T8 | 674 K | 710 K |
| E80 mid/low flow | T9_pair - T8 | -60 K | +14 K |

Flow-slope diagnostic for the fitted v12a source:

| irradiance | sensor | model slope | experimental slope |
|---:|---|---:|---:|
| 456 kW/m2 | T8 | -31 K/(L/min) | -34 K/(L/min) |
| 456 kW/m2 | T9_pair | -24 K/(L/min) | -17 K/(L/min) |
| 456 kW/m2 | T10_pair | -16 K/(L/min) | -2 K/(L/min) |
| 456 kW/m2 | T3 | -13 K/(L/min) | +1 K/(L/min) |
| 304 kW/m2 | T8 | -22 K/(L/min) | -24 K/(L/min) |
| 304 kW/m2 | T9_pair | -16 K/(L/min) | -13 K/(L/min) |
| 304 kW/m2 | T10_pair | -11 K/(L/min) | -3 K/(L/min) |
| 304 kW/m2 | T3 | -9 K/(L/min) | 0 K/(L/min) |
| 256 kW/m2 | T8 | -24 K/(L/min) | -21 K/(L/min) |
| 256 kW/m2 | T9_pair | -17 K/(L/min) | -14 K/(L/min) |
| 256 kW/m2 | T10_pair | -11 K/(L/min) | -6 K/(L/min) |
| 256 kW/m2 | T3 | -9 K/(L/min) | -3 K/(L/min) |

Interpretation:

1. v12a confirms the user's concern. Source redistribution acts mostly as an
   axial level correction and does not fix the front/back flow-dependence
   mismatch.
2. The optimizer does prefer moving some energy downstream, improving
   `T9_pair - T8` by about 6 K on average, but this is far too small to create
   the high-flow internal hot ridge.
3. The model rear remains much too flow-sensitive. In the experiments T10/T11
   and T3 are nearly flat with flow at 304 and 456 kW/m2, while the model still
   gives large negative slopes.
4. The result suggests the next missing mechanism is not only source placement.
   The gas/rear response needs an additional structure that decouples downstream
   measured temperatures from the through-flow more strongly than the current
   1D gas-solid model allows.

Next recommended check:

```text
Do not proceed to per-irradiance power scaling yet if the immediate question is
the flow-slope pattern.

Instead, build a diagnostic v12b/v13 that explicitly targets downstream
flow-slope decoupling, for example:

1. Add an axial bypass / active-flow fraction model:
   only a fraction phi(z, Re) of the gas exchanges with the solid downstream,
   while the remainder bypasses thermally. This can make T3 and rear sensors
   less flow-sensitive without forcing all temperatures horizontally.

2. Or add a two-temperature gas/solid exchange diagnostic:
   a near-wall/interstitial gas branch exchanges strongly with the solid,
   while a core/bypass gas branch carries most of the mass flow.

3. Keep source distribution fixed at the v12a fitted values during that test
   so flow-decoupling is not hidden by another optical refit.
```

## 2026-07-21 - v13 Axial Radiation and Thermocouple Diagnostics

Motivation:

After rejecting a two-gas-branch model as likely overfitting, v13 tested two
more defensible possibilities without adding fitted parameters:

```text
1. nonlocal axial surface-to-surface radiative exchange inside the receiver
2. a fixed front-biased thermocouple-wire measurement model
```

Implementation:

```text
1D_v13.jl
run_1D_v13.jl
test/smoke_1D_v13.jl
```

v13 freezes the v11/v12a closure:

| parameter | fixed value |
|---|---:|
| A_Nu | 1.0996 |
| B_Re | 1.0000 |
| C_Pr | 0.3183 |
| L_entry_m | 0.0379 |
| beta_opt | 184.67 1/m |
| front_dep | 0.0 |

The effective source distribution is the v12a result:

```text
first-cell source fraction = 0.637
downstream source fraction = 0.363
```

The axial view-radiation diagnostic uses a fixed geometry-tied kernel:

```text
emissivity = 0.85
maximum row-sum view fraction = 0.06
axial decay length = 50 mm
```

The thermocouple diagnostic is only a measurement-layer correction; it does not
change the energy balance:

```text
T_tc(z) = T_solid(z) - f_wire(z) * (T_solid(z) - 298.15 K)
f_wire(z) = 0.25 * exp(-z / 25 mm)
```

Validation:

```text
test/smoke_1D_v13.jl passed: 83/83 tests
run_1D_v13.jl completed
```

Generated artifacts:

```text
summaries/1D_v13/analysis_results_all_variants_1D_v13.csv
summaries/1D_v13/steady_results_*_1D_v13.csv
summaries/1D_v13/diagnostics/cell_diagnostics_*_1D_v13.csv
summaries/1D_v13/plots/*.png
```

Variants:

```text
baseline_v12a_source
axial_view_rad_fixed
axial_view_rad_upper
tc_wire_fixed
combined_fixed
```

Mean signed heating steady errors:

| variant | T8 | T9_pair | T10_pair | T3 | T9_pair - T8 |
|---|---:|---:|---:|---:|---:|
| baseline_v12a_source | -192 K | -241 K | -159 K | -121 K | -50 K |
| axial_view_rad_fixed | -193 K | -242 K | -159 K | -120 K | -49 K |
| axial_view_rad_upper | -194 K | -243 K | -159 K | -120 K | -49 K |
| tc_wire_fixed | -256 K | -250 K | -160 K | -121 K | +6 K |
| combined_fixed | -256 K | -250 K | -160 K | -120 K | +6 K |

Compact representative ordering:

| variant | E67 T9-T8 | E80 T9-T8 | E76 T9-T8 | mean |Q_view| | mean T8 TC shift |
|---|---:|---:|---:|---:|---:|
| baseline_v12a_source | -25.7 K | -59.6 K | -102.0 K | 0.000 W | 0.0 K |
| axial_view_rad_fixed | -25.8 K | -59.4 K | -100.8 K | 0.586 W | 0.0 K |
| axial_view_rad_upper | -26.0 K | -58.8 K | -98.6 K | 1.711 W | 0.0 K |
| tc_wire_fixed | +24.9 K | -6.9 K | -27.4 K | 0.000 W | -63.9 K |
| combined_fixed | +24.8 K | -6.7 K | -26.5 K | 0.586 W | -63.8 K |

Flow-slope diagnostic for the combined fixed variant:

| irradiance | sensor | model slope | experimental slope |
|---:|---|---:|---:|
| 456 kW/m2 | T8 | -26 K/(L/min) | -34 K/(L/min) |
| 456 kW/m2 | T9_pair | -23 K/(L/min) | -17 K/(L/min) |
| 456 kW/m2 | T10_pair | -16 K/(L/min) | -2 K/(L/min) |
| 456 kW/m2 | T3 | -13 K/(L/min) | +1 K/(L/min) |
| 304 kW/m2 | T8 | -18 K/(L/min) | -24 K/(L/min) |
| 304 kW/m2 | T9_pair | -15 K/(L/min) | -13 K/(L/min) |
| 304 kW/m2 | T10_pair | -11 K/(L/min) | -3 K/(L/min) |
| 304 kW/m2 | T3 | -9 K/(L/min) | 0 K/(L/min) |

Interpretation:

1. The axial view-radiation model is too weak to explain the mismatch. Even the
   upper-bound variant produces only about `1-2 W` of total redistributed
   axial radiative exchange in representative heating cases. It barely changes
   T9-T8 or the steady-state errors.
2. The thermocouple-wire diagnostic can flip the high-flow E67 ordering toward
   `T9 > T8`, but only by imposing a large T8 downward measurement shift
   (`~64 K` mean, and up to `~85 K` in low-flow hotter cases). This improves
   the ordering metric while making the absolute T8 error much worse.
3. Combining fixed axial radiation with the thermocouple model does not solve
   the rear flow-slope issue. T10/T11 and T3 remain far too flow-sensitive.
4. Therefore the internal hot-ridge mismatch is unlikely to be solved by
   conservative axial radiation alone or by a defensible simple thermocouple
   measurement correction.

Recommendation:

```text
The next defendable direction should not be another fitted transport branch.

Before adding more 1D closure flexibility, perform an energy/input and
measurement-consistency audit:

1. Compare measured gas enthalpy rise + estimated losses with the modeled
   absorbed power for each irradiance.
2. Revisit whether the experimental "irradiance" corresponds to aperture
   average, peak, illuminated area average, or receiver-absorbed equivalent.
3. Add a diagnostic per-irradiance power-scale solve, not as final physics, but
   to estimate the missing watts required after v12a/v13 fixed mechanisms.
4. If the required power scale is reasonable and consistent by irradiance,
   then the model is energy-starved.
5. If the required power scale varies mainly with flow, then the remaining
   issue is still flow/measurement coupling rather than optics.
```

## 2026-07-21 - v14 Absorbed-Power Scale Audit

Motivation:

v13 showed that fixed axial radiation and a simple thermocouple correction do
not explain the mismatch. The next question was whether the model is simply
energy-starved under the nominal irradiance/absorptivity inputs, and whether
the required correction is grouped by irradiance level or by flow rate.

Implementation:

```text
1D_v14.jl
run_1D_v14.jl
test/smoke_1D_v14.jl
```

v14 keeps the v11/v12a closure fixed:

| parameter | fixed value |
|---|---:|
| A_Nu | 1.0996 |
| B_Re | 1.0000 |
| C_Pr | 0.3183 |
| L_entry_m | 0.0379 |
| beta_opt | 184.67 1/m |
| front_dep | 0.0 |

The only fitted parameters are absorbed-power scale factors:

```text
scale_456
scale_304
scale_256
```

The model still uses:

```text
eta_abs = 0.8 fixed
eta_opt = 1.0 fixed before scale
Rosseland off
axial view radiation off
thermocouple measurement model off
```

Validation:

```text
test/smoke_1D_v14.jl passed: 92/92 tests
run_1D_v14.jl completed
```

Generated artifacts:

```text
summaries/1D_v14/analysis_results_all_variants_1D_v14.csv
summaries/1D_v14/steady_results_baseline_scale1_1D_v14.csv
summaries/1D_v14/steady_results_per_irradiance_power_fit_1D_v14.csv
summaries/1D_v14/per_case_power_scale_audit_1D_v14.csv
summaries/1D_v14/flow_slopes_*_1D_v14.csv
summaries/1D_v14/diagnostics/cell_diagnostics_*_1D_v14.csv
summaries/1D_v14/plots/*.png
```

Per-irradiance fitted scales:

| irradiance | scale | nominal absorbed power | scaled absorbed power | added absorbed power |
|---:|---:|---:|---:|---:|
| 456 kW/m2 | 1.670 | 131.7 W | 219.9 W | +88.2 W |
| 304 kW/m2 | 1.717 | 87.8 W | 150.8 W | +63.0 W |
| 256 kW/m2 | 0.983 | 73.9 W | 72.7 W | -1.3 W |

Objective change:

```text
baseline v12a/v13 scale=1 = 0.20035
per-irradiance power fit  = 0.09531
```

Mean signed heating steady errors:

| sensor | baseline | per-irradiance power fit |
|---|---:|---:|
| T8 | -192 K | -28 K |
| T9_pair | -241 K | -92 K |
| T10_pair | -159 K | -34 K |
| T3 | -121 K | -15 K |
| T2 | -11 K | -2 K |
| T9_pair - T8 | -50 K | -65 K |
| T10_pair - T8 | +33 K | -6 K |

Interpretation of these errors:

1. A large input-energy deficit is strongly indicated for the 456 and
   304 kW/m2 datasets. Raising absorbed power by about `+60 to +90 W` greatly
   improves absolute temperatures.
2. The 256 kW/m2 group does not want added power. Its fitted scale is almost
   exactly 1.0.
3. Power scaling improves the global level but worsens the internal hot-ridge
   metric `T9_pair - T8`: the model becomes too front-hot after the added power.
4. This confirms that the energy deficit and the axial/flow-shape mismatch are
   separate problems.

Per-case power-scale audit:

| irradiance | mean scale | min scale | max scale | mean added power |
|---:|---:|---:|---:|---:|
| 456 kW/m2 | 1.677 | 1.550 | 1.776 | +89 W |
| 304 kW/m2 | 1.755 | 1.617 | 1.972 | +66 W |
| 256 kW/m2 | 1.008 | 0.898 | 1.152 | +1 W |

Per-case scale slope versus flow:

| irradiance | d(scale)/d(flow) |
|---:|---:|
| 456 kW/m2 | -0.0126 per L/min |
| 304 kW/m2 | +0.0243 per L/min |
| 256 kW/m2 | +0.0274 per L/min |

The scales cluster more by irradiance than by flow, especially for the two
higher irradiance levels. This supports an irradiance/input calibration issue
for 304 and 456 kW/m2.

Flow-slope consequence after per-irradiance power scaling:

| irradiance | sensor | model slope | experimental slope |
|---:|---|---:|---:|
| 456 kW/m2 | T8 | -41 K/(L/min) | -34 K/(L/min) |
| 456 kW/m2 | T9_pair | -31 K/(L/min) | -17 K/(L/min) |
| 456 kW/m2 | T10_pair | -21 K/(L/min) | -2 K/(L/min) |
| 456 kW/m2 | T3 | -17 K/(L/min) | +1 K/(L/min) |
| 304 kW/m2 | T8 | -32 K/(L/min) | -24 K/(L/min) |
| 304 kW/m2 | T9_pair | -23 K/(L/min) | -13 K/(L/min) |
| 304 kW/m2 | T10_pair | -15 K/(L/min) | -3 K/(L/min) |
| 304 kW/m2 | T3 | -12 K/(L/min) | 0 K/(L/min) |
| 256 kW/m2 | T8 | -23 K/(L/min) | -21 K/(L/min) |
| 256 kW/m2 | T9_pair | -16 K/(L/min) | -14 K/(L/min) |
| 256 kW/m2 | T10_pair | -11 K/(L/min) | -6 K/(L/min) |
| 256 kW/m2 | T3 | -8 K/(L/min) | -3 K/(L/min) |

This is the key diagnostic result:

```text
Power scaling fixes much of the absolute heating level,
but it makes the downstream flow-slope problem more obvious.
```

Recommended next direction:

```text
1. Treat per-irradiance absorbed-power calibration as necessary for the
   304 and 456 kW/m2 datasets, but not sufficient.

2. Do not use a single global power scale:
   the 256 kW/m2 group rejects it.

3. Do not interpret the remaining mismatch as purely optical:
   after adding the missing watts, T10/T11 and T3 are still much too
   flow-sensitive.

4. Next model revision should combine:
   a) per-irradiance power-scale inputs fixed from v14, and
   b) a conservative downstream flow-sensitivity reduction that is not a
      two-branch overfit.

Candidate conservative form:

   replace the gas-solid exchange effectiveness downstream with a bounded
   thermal-equilibration ceiling:

   epsilon_eff(z, flow) = min(epsilon_Nu(z, flow), epsilon_cap(z))

   where epsilon_cap(z) is monotonic and physically interpreted as unresolved
   nonuniform exchange / finite thermal participation, not a new heat source.

This should be tested as a diagnostic with one or two global parameters and
reported directly through epsilon_eff(z), Qgas(z), and the flow-slope tables.
```
