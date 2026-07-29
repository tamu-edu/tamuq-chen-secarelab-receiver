# 1D Manuscript-Gap Strategy

Date: 2026-07-28

## Current Status

The 1D model remains a **Role A** manuscript asset: it can corroborate
mechanisms and stress-test interpretations, but it should not yet be used as a
**Role B** source of validated macroscopic coefficients.

The v31 refactor moved the model much closer to the manuscript strategy by
removing the worst artificial closures and enforcing measured thermal
inventory:

- no direct distributed rear-core sink,
- no hard rear/flange cooling switch,
- T3 is a gas-temperature observable,
- total participating heat capacity is constrained near 301 J/K,
- staged calibration separates rear topology, transport, and full source-power
  release,
- invariant diagnostics are now written to CSV.

The remaining gap is not a coding cleanup problem; it is a model-identification
problem. v31 still improves fit by driving some physically meaningful paths or
power terms to implausible edges.

## Distance From Manuscript Gate

| Gate | Target | v31 Status |
|---|---:|---|
| Fit quality | near empirical baseline without surrogate sink/switch/T3 blend | improved but still ~3.3x v27 objective |
| Total heat capacity | 301 +/- 23 J/K | pass: 302.5 J/K |
| Apparent Nu law | `3.1e-4 Re^1.44` | partial: prefactor 2.23e-4, exponent 1.334 |
| Nu trend quality | high log-linearity | partial: R2 0.961 vs 0.97 |
| Volumetric inversion threshold | epsilon* = 0.66 +/- 0.03 | fail: epsilon_mean 0.876-0.983 |
| Lambda_107 | `0.038 + 8.3e-4 Re` | fail: wrong sign/trend |
| External loss | 0.10-0.16 W/K basis | fail: 0.250-0.341 |
| Power convention | 456/304 elevated, 256 nominal/lower but not collapsed | fail: 256 scale = 0.636 |
| Parameter interiority | no critical edge collapse | fail: C bounds, flange scale, G_rear_tube = 0 |
| Heating/cooling simultaneous | no cooling-only patch | partial: smooth gate retained, but edge parameters remain |

## Specific v32 Strategy

1. **Preserve v31 as the scaffold.**
   Do not return to v27/v25 sink closures. Keep the capacity gate, staged
   calibration, pure-gas T3, and invariant diagnostics.

2. **Repair the rear-tube/adaptor heat path.**
   v31 sets `G_rear_tube = 0`, which leaves T3 too cold and routes too much
   heat through the cavity/loss path. v32 should replace the free zero-capable
   scalar with either:
   - a geometry-derived lower-bounded conductance prior, or
   - a short manifold/tube gas-control-volume exchange with a physically
     bounded wall-gas conductance.

3. **Constrain the source-power convention explicitly.**
   Add a weak regularization or staged bound policy so `scale_256` cannot
   collapse while `scale_456` and `scale_304` remain elevated. The target is not
   to force equal scales, but to prevent the optimizer from using low-flux power
   suppression as a substitute for missing rear/outlet physics.

4. **Promote invariants into the objective only after diagnostic agreement.**
   Add soft penalties for:
   - `Nu_app_Re_exponent`,
   - `Lambda_107` slope/sign,
   - `K_loss` range,
   - total heat capacity.

   Keep these penalties reported separately so the scientific cost is visible.

5. **Treat T3 as a downstream gas-path validation, not as a wall blend.**
   The next T3 fix should not blend wall and gas temperatures. It should make
   the downstream gas path physically warm enough through tube/manifold heat
   exchange and source/loss balance.

6. **Only then extract coefficients.**
   Convective/radiative/conductive coefficients should be reported as validated
   only when the same parameter set passes fit, inventory, Lambda, Nu, K-loss,
   and power-convention gates without critical edge collapse.

## Immediate Implementation Tasks

For v32:

- fork from v31,
- set a physically positive rear-tube coupling floor or prior,
- add a source-scale regularization term with explicit output,
- add optional invariant penalties behind named weights,
- keep runner no-plot mode,
- recalibrate in stages: `rear`, `transport`, `full`,
- compare v31 and v32 with the invariant CSVs before considering another
  topology change.
