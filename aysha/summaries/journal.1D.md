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
