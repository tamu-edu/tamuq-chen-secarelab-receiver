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
