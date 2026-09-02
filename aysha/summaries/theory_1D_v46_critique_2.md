Rechecking the complete September 1 artifacts changes part of my previous assessment. The current full results are not the optimized v46 result described in the theory document—they are a later forward run using the default/manual parameter vector.

No files were modified.

## Corrected overall verdict

V46 successfully repairs several structural conservation problems and substantially improves cooling predictions over v45. However, the currently saved “full results” do not constitute a calibrated or validated v46 model. They are best classified as a conservative forward-model baseline.

The most immediate problem is artifact provenance:

| Item | Theory document | Current full artifacts |
|---|---:|---:|
| Run status | `MaxTime` optimization | `Manual`, calibration skipped |
| Objective | 1.973576 | 2.430907 |
| \(A_{Nu}\) | 3.388 | 2.500 |
| \(G_{core,perim}\) | 12.347 | 6.000 |
| \(\chi\) | 1.000 | 0.450 |
| \(h_{suction}\) | 331.6 W/m²K | 250 W/m²K |
| Mean HTC range | 44.7–84.3 W/m²K | 30.7–55.2 W/m²K |

The current run explicitly says `return_code=Manual` in [optimization_summary_1D_v46.txt](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/1D_v46/optimization_summary_1D_v46.txt). The runner defaults to skipping calibration in [run_1D_v46.jl](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/run_1D_v46.jl:487).

Consequently, filenames containing `fitted` are currently misleading.

## What v46 has genuinely addressed

- Total measured mass flow is used throughout the receiver and rear tube.
- Front suction heat is transferred consistently into gas enthalpy.
- Core plus perimeter source power exactly equals delivered power.
- Lamp-off-specific conductance multipliers were removed.
- The instantaneous first-law ledger closes to machine precision.
- Cooling is excluded from the stated heating objective.
- Perimeter deposition is now active in the current result: \(\chi=0.45\), so 45% of power enters the core and 55% the perimeter.
- Selected-grid tests show reasonably stable predicted sensor temperatures. Changing from 15 to 50 axial nodes changed final temperatures by generally only a few kelvin in E67, E71, and E81.

These are real improvements in model architecture.

## Full-result performance

### Heating

The current v46 forward result is materially worse than v45 for every sensor group:

| Sensor | V45 mean RMSE | V46 mean RMSE |
|---|---:|---:|
| T8 | 51.4 K | 53.1 K |
| T12 | 69.0 K | 130.1 K |
| T11 | 52.9 K | 70.2 K |
| T9 | 59.9 K | 80.7 K |
| T10 | 35.8 K | 41.3 K |
| T3 | 22.3 K | 43.6 K |
| T2 | 4.9 K | 8.3 K |

Across all 105 heating signal/case combinations:

- Mean RMSE: **61.1 K**
- Median RMSE: **47.4 K**
- Maximum RMSE: **237.3 K**
- Mean absolute endpoint error: **59.0 K**

The strongest failure is T12, which is systematically underpredicted:

- Mean endpoint bias: **−120.8 K**
- Worst endpoint error: **−217.5 K**

The high-power cases E67–E72 are especially poor. Their mean per-case RMSEs range from approximately 74 to 115 K, except E71 at 60 K.

### Cooling

Cooling is a genuine improvement over v45:

| Sensor | V45 mean RMSE | V46 mean RMSE |
|---|---:|---:|
| T8 | 24.3 K | 20.9 K |
| T12 | 49.4 K | 21.6 K |
| T11 | 54.8 K | 25.5 K |
| T9 | 48.1 K | 10.9 K |
| T10 | 57.4 K | 27.3 K |
| T3 | 57.6 K | 35.8 K |
| T2 | 9.2 K | 10.5 K |

Case-level results:

- C80: mean RMSE **12.9 K**
- C81: mean RMSE **15.2 K**
- C69: mean RMSE **37.3 K**

Thus moderate- and zero-flow cooling are encouraging, but C69 remains weak. Cooling is also initialized from measured core, perimeter, tube, and cavity temperatures, so this is not a blind forecast.

## Important new finding: gas heating is dominated by the front boundary

The cell diagnostics show that the current parameterization does not behave primarily as a core-to-gas volumetric receiver.

For each heating case, 20–24 of the 25 core cells have negative local gas heat transfer: after rapid front heating, the gas becomes hotter than the downstream solid and returns heat to it. That reversal can be physically possible, but it changes the interpretation of the inferred coefficient.

Examples of the final gas-energy decomposition:

| Case | Suction heat | Net core-to-gas | Rear tube-to-gas |
|---|---:|---:|---:|
| E67 | +70.2 W | +64.3 W | −13.5 W |
| E71 | +60.9 W | +24.2 W | −12.8 W |
| E76 | +44.4 W | +2.1 W | −9.7 W |
| E81 | +29.7 W | +0.23 W | −6.6 W |

For E81, only about 1% of reported total gas heating originates as net receiver-core convection. Almost all comes through the fitted suction boundary.

This does not violate energy conservation, but it means:

- Total gas heat cannot be presented as monolith convective extraction.
- The fitted Nusselt correlation is weakly constrained by net core heat transfer at low flow.
- The suction coefficient, source partition, and Nusselt parameters remain strongly confounded.
- Rear-tube heat transfer is negative in every heating case.

## Axial and flow behavior remains unresolved

### Wall inversion

The model reproduces the sign of \(T_{12}-T_8\) in only **5 of 15 cases**. It incorrectly predicts \(T_{12}<T_8\) in ten cases where the experiment often shows a warmer mid-perimeter region.

### Core axial temperature drop

For \(T_9-T_{10}\):

- Mean absolute error: **54.1 K**
- Worst error: **98.4 K**

The model systematically produces an axial core drop that is too small.

### Flow slopes

Front and external measurements are relatively reasonable:

- T8 slope errors: approximately 2.7–4.4 K/LPM.
- T2 slope errors: approximately 0.9–1.1 K/LPM.

Deep measurements remain much too flow-sensitive. At 456 kW/m²:

- T11: model −21.84 versus experiment −1.37 K/LPM.
- T10: model −20.86 versus experiment −3.54 K/LPM.
- T3: model −14.95 versus experiment +0.54 K/LPM.

This suggests that the external/front energy pathways are becoming plausible, while the internal gas–solid and rear transport network remains structurally incorrect.

## Conservation is correct, but “steady state” is not

The instantaneous residual is below \(10^{-13}\) W, confirming internal bookkeeping. However, the endpoint energy balances retain substantial storage:

- Minimum storage fraction: **12.3%**
- Maximum storage fraction: **24.0%**
- Mean storage fraction: **18.8%**

For E67:

- Delivered: 220.59 W
- Accounted external/gas losses: 174.77 W
- Continuing storage: 45.82 W

Therefore:

- The files named `steady_results` contain experimental end points, not verified steady states.
- `max_steady_residual_W` is actually the instantaneous conservation residual.
- The theory’s statement that steady-state residuals are below \(10^{-13}\) W is false for these runs.

## Coefficient extraction is not yet defensible

The current mean HTC range is 30.7–55.2 W/m²K, not the range reported in the theory. In addition:

- The saved HTC is the raw fluid coefficient \(h\), not the effective series coefficient \(U_{eff}\).
- Its arithmetic axial mean is grid-dependent because \(h\propto z^{-C_z}\) near the inlet.
- E67 mean HTC rises from 51.9 to 59.1 W/m²K when changing from 15 to 50 nodes—about 14%.
- Sensor temperatures are much less grid-sensitive, meaning output convergence does not imply coefficient convergence.
- `NU_FLOOR=0` still provides no fully developed asymptote.
- Most core optical power is deposited extremely near the front: 80.5% in the first 5.48 mm and 94.2% within approximately 11 mm.

The current result therefore does not validate either a representative Nusselt correlation or a volumetric optical extinction coefficient.

## Pending work, in priority order

1. **Resolve result provenance.** Decide whether the manual parameter vector or the earlier MaxTime optimum is the intended v46 result. Do not mix them.

2. **Regenerate documentation from the authoritative run.** The theory, journal, parameter table, energy table, cooling table, and HTC discussion currently describe a different parameter vector.

3. **Run an actual calibration.** The present full artifacts are forward baseline results, despite their `fitted` filenames.

4. **Fix the heating model before interpreting cooling success.** Prioritize T12, T9, wall inversion, axial core drop, and deep flow slopes.

5. **Investigate suction dominance.** Separate:
   - gas heat acquired at the aperture,
   - gas heat acquired in the monolith,
   - gas heat lost in the rear tube.

6. **Constrain source partition physically.** The internal flux calculation reports a direct receiver fraction of only 0.069, while the manual model uses \(\chi=0.45\). This may be defensible only if reflection or redistribution is explicitly modeled and documented.

7. **Define steady-state criteria.** Report \(dE/dt\), the true steady-flux residual, and a normalized threshold before labelling results steady.

8. **Make HTC/Nu grid-invariant.** Report an area- or heat-transfer-weighted effective coefficient and demonstrate mesh convergence of the extracted coefficient, not just temperatures.

9. **Repair zero-flow T3 physics.** At zero flow, the gas array is assigned local solid temperatures despite zero gas heat transfer. C80 T3 is therefore not a solved gas-energy prediction.

10. **Perform identifiability analysis.** Multistart fitting, profiles, sensitivity matrices, uncertainty intervals, and parameter-correlation analysis remain absent.

11. **Reconcile total capacity.** The current 270.53 J/K remains below the stated 278–324 J/K interval, while capacity regularization is disabled.

12. **Create condition-level holdouts.** Cooling validation should be supplemented with withheld heating flux/flow conditions.

The complete artifacts therefore strengthen a narrower conclusion: v46 has a much cleaner conservative architecture and encouraging cooling behavior, but its current full result is an uncalibrated baseline with unresolved internal heat-transfer physics and inconsistent documentation.