# 1D Model Lessons Learned and Strategy Forward to v29

## Purpose

This note pauses implementation before `1D_v29.jl` and consolidates the lessons from the previous 1D model series. The goal is to prevent v29 from becoming another parameter patch and instead make it a targeted test of the remaining physical hypothesis: the receiver cannot be represented only by axial core/perimeter temperatures plus an external tube/cavity path; it needs a physically bounded rear/adaptor thermal reservoir between the ceramic exit region and the measured/external rear hardware.

The overarching scientific objective remains the extraction and validation of effective macroscopic heat-transfer coefficients for the structured monolithic receiver. Therefore, v29 should improve physical identifiability, not merely reduce the residual.

## Executive Strategy

Do not proceed by restoring the v27 direct rear-core sink unchanged. v27 remains the best numerical fit, but its split direct rear sink is an effective surrogate. Do not proceed from v28 as if its failed fit were a calibrated model; v28 is a negative control showing that removing the distributed rear sink without adding a physical rear reservoir makes the topology under-complete.

The next version should be a constrained topology test:

1. Start from the v27/v25 energy-accounting framework, not the collapsed v28 fitted parameters.
2. Remove the direct distributed receiver rear sink.
3. Add one explicit rear/adaptor thermal reservoir state with bounded heat capacity and physically interpretable conductances.
4. Keep `T3` as gas temperature in the first v29 implementation.
5. Keep the smooth operating-state gate, but apply it only to boundary hardware loss paths, not to an arbitrary receiver-volume sink.
6. Treat irradiance-dependent power scales as experimental input-calibration nuisance parameters, not as absorptivity or intrinsic material properties.
7. Judge v29 by residual structure, heat-path plausibility, parameter interiority, and cooling behavior, not by objective value alone.

## Lessons by Model Phase

### v1-v3: Numerical Backbone

The finite-volume solid receiver and exact quasi-steady gas march are the durable core of the 1D model. These versions established the most important numerical principle: the gas energy exchange should be computed conservatively cell by cell, and the solid residual should close the energy balance without hidden source/sink terms.

Retain:

- Axial finite-volume energy accounting.
- Exact NTU-style gas update within each cell.
- Clear separation between model equations, runner/calibration scripts, and diagnostics.
- Conservative bookkeeping of absorbed input, gas uptake, storage, and losses.

Avoid:

- Replacing the gas march with an approximate local linearization unless it is proven energy-equivalent.
- Adding residual correction terms that do not appear in the reported energy accounting.

### v4-v5: Shape Flexibility and Overfitting

The early flexible axial exchange/source-shape experiments showed that arbitrary floors, free redistribution factors, and broad empirical profiles can reduce residuals while damaging physical credit. They are useful for diagnostics but weak as final scientific parameters.

Lesson:

- Shape freedom should be used only as a hypothesis test. Once a shape is indicated, v29 should convert it into a physical state or geometric coupling.

Risk:

- A model can appear improved because a profile compensates for missing rear/topological inventory.

### v6-v8: Rear and Cavity Inventory

The rear, cavity, tube, and adaptor region cannot be ignored. Versions in this range showed that adding rear/cavity thermal inventory improves the timing and downstream thermocouple behavior, especially when T2/cavity behavior is predicted instead of prescribed.

Lesson:

- Rear hardware is not just a boundary condition; it stores heat and mediates heat flow.
- A rear thermal mass can improve physical realism, but if its capacity and conductances are unconstrained it becomes another empirical reservoir.

Retain:

- Prediction of rear/cavity-like states where possible.
- Comparison against measured T2/cavity signals when available.

Avoid:

- Free rear masses with no geometry or heat-capacity plausibility check.

### v9-v15: Convective Law, Radiation, and Source Diagnostics

Changing the heat-transfer law alone did not solve the signed residual pattern. Rosseland axial radiation within tested bounds was too small to explain the missing behavior. Source redistribution and apparent-Nusselt variants helped diagnostically but converged toward effective compensation rather than a validated physical coefficient.

Lessons:

- High fitted Nusselt coefficients are not automatically physical. They may be compensating for missing source distribution, side/perimeter exchange, rear losses, or sensor mapping.
- Radiation must be reported with optical-thickness and heat-flow diagnostics. A fitted radiation coefficient without a magnitude check has weak scientific credit.
- Source distribution should be constrained by flux/optical accounting and not allowed to absorb every residual.

Retain:

- The ability to report actual Nu ranges, Reynolds ranges, and gas-side heat-flow fractions.
- Irradiance-specific source/power scale diagnostics.

Avoid:

- Claiming fitted `A_Nu` or `B_Re` as a validated channel-scale Nusselt correlation unless other compensating pathways are controlled.

### v16-v20: Heterogeneous Topology Tests

Heterogeneous active/wall, side-wall, and measurement-bias hypotheses were either falsified or only partially successful. The major constructive result was the emergence of a two-zone macro-topology: core and perimeter behave differently and must exchange heat.

Lessons:

- A two-zone core/perimeter model is necessary for macroscopic interpretation.
- Perimeter coupling, perimeter capacity, and perimeter conduction are real candidates for effective transport parameters.
- T3/outlet behavior cannot be fixed by side heating alone; rear/outlet topology matters.

Retain:

- Core/perimeter states and conductive coupling.
- Separate accounting for core-to-gas, core-to-perimeter, perimeter-to-hardware, and axial conduction/radiation terms.

Avoid:

- One-zone receiver models as final candidates.
- Pure measurement-bias explanations unless they are tested after the physical topology is adequate.

### v21-v25: Cooling Artifacts and Energy Accounting

Cooling exposed artificial behavior that heating fits could hide. Nonphysical cooling upturns in T3/T10/T11 were decisive diagnostics. The energy-accounting family also corrected a hidden absorptivity/power-scale issue: the old effective absorption factor was not defensible as a material property.

Lessons:

- Cooling trajectories are strong physical constraints because there is no lamp source to mask topology errors.
- Initial-condition alignment and cooling start time affect apparent losses and must be explicitly disclosed.
- Irradiance power scales around 2 for 456 W and 304 W cases are input/participating-flux calibration factors, not intrinsic absorptivity.
- Low-irradiance behavior needs deeper energy/topology interpretation, not merely uniform power rescaling.

Retain:

- Cooling monotonicity and signed-residual checks.
- Participating heat-capacity reporting.
- Energy-accounting output tables.

Avoid:

- Hidden eta/absorptivity corrections.
- Cooling-only fixes that create heating overcooling.

### v26-v28: Rear Sink and T3 Lessons

v26 improved cooling-side handling but overcooled heating because the same rear-loss representation acted in incompatible operating regimes. v27 split the rear sink into heating and cooling terms, achieving the best fit but with weak scientific credit because the direct rear sink is a surrogate. v28 removed the direct sink, made T3 pure gas, and introduced a smooth operating-state gate on the explicit rear hardware path; the fit collapsed badly.

Lessons:

- v27 proves that some rear/outlet heat-removal pathway is needed.
- v28 proves that the existing explicit tube/cavity/flange path is insufficient by itself.
- The missing element is likely a rear/adaptor/holder thermal reservoir and contact network, not a distributed volumetric sink inside the receiver.
- T3 as pure gas is scientifically cleaner, but it may fail unless the rear/outlet thermal topology is physically complete.
- The smooth operating-state gate is preferable to a hard lamp-on/off switch, but it should modulate real hardware conductances only.

Retain:

- Smooth operating-state gates.
- Higher allowed `A_Nu` upper bound for diagnostic freedom.
- T3 as gas temperature for v29a.

Avoid:

- Direct receiver-volume rear sinks.
- Hard lamp-on/off rear-sink switches.
- Immediate restoration of a T3 gas-wall blend before testing a better rear topology.

## Non-Negotiable Constraints for v29

v29 should satisfy these modeling rules from the start:

- No arbitrary direct distributed sink from receiver core cells to ambient/flange.
- No hard operating-state switch.
- No hidden absorptivity multiplier.
- No T3 wall blending in the baseline v29a model.
- No unreported heat path.
- No broad, unconstrained artificial source-shape correction.
- All fitted parameters must have physical units, bounds, and interpretation.
- Report whether each fitted parameter is interior or boundary-active.
- Preserve exact gas energy accounting.
- Preserve core/perimeter macro-topology.
- Preserve the 17 C cooling inlet/ambient correction and explicit cooling time alignment, with documentation.

## Proposed v29 Physical Model

### State Topology

Use the v27/v25 energy-accounting receiver model as the stable base, but replace the direct rear sink with a rear/adaptor reservoir:

- `T_core[i]`: active receiver/core solid temperature.
- `T_perim[i]`: perimeter or less-active solid temperature.
- `T_tube[j]`: explicit rear tube or downstream hardware state if already present.
- `T_cavity`: rear/cavity effective temperature if already present.
- `T_rear`: new rear/adaptor/holder reservoir state.

The new `T_rear` state should represent solid material and contact inventory near the receiver exit: ceramic rear face, holder/adaptor contact region, and thermally participating rear support mass. It should not represent an unlimited environmental sink.

### Heat Paths

The minimum v29a path should be:

```text
core[N] / perimeter[N] -> rear/adaptor reservoir -> tube/cavity/flange hardware -> ambient/water boundary
```

Candidate equations:

```text
Q_core_rear  = G_core_rear  * (T_core[N]  - T_rear)
Q_perim_rear = G_perim_rear * (T_perim[N] - T_rear)
Q_rear_tube  = G_rear_tube  * (T_rear - T_tube[1])
Q_rear_cav   = G_rear_cav   * (T_rear - T_cavity)
```

The rear reservoir ODE:

```text
C_rear * dT_rear/dt =
    Q_core_rear
  + Q_perim_rear
  - Q_rear_tube
  - Q_rear_cav
  - Q_rear_flange_optional
```

For v29a, omit `Q_rear_flange_optional` unless the minimal reservoir cannot reproduce cooling without restoring an artificial sink. If a direct reservoir-to-flange path is later needed, add it as v29b with a documented physical contact interpretation and smooth operating-state gate.

### Smooth Operating-State Gate

Use a continuous lamp-off/cooling gate only on hardware boundary conductance:

```text
s_off(t, I) = [1 / (1 + (I / I_gate)^n)] * [1 - exp(-t / tau_gate)]
G_eff = G_base * (1 + gain_off * s_off)
```

Initial choices:

- `I_gate = 1000 W/m2`, fixed unless diagnostics demand otherwise.
- `n = 4`, fixed to keep the transition smooth but sharp.
- `tau_gate` fitted or bounded by cooling transient onset.
- `gain_off` fitted with conservative upper bound.

Interpretation:

- The gate represents operating-state-dependent rear hardware heat transfer, such as enhanced effective contact/convection during no-flux cooling.
- It must not act as a hidden receiver-volume loss term.

### T3 Mapping

Baseline v29a:

```text
T3_model = gas temperature at the measured downstream axial position
```

Do not blend T3 with wall/rear solid temperature in v29a. If T3 remains systematically biased after the rear reservoir is tested, create a later explicitly labeled sensor-submodel branch. That later branch must distinguish between:

- true gas thermocouple position,
- radiation/conduction contamination of the bead,
- rear hardware recirculation or mixing,
- outlet manifold thermal exchange.

T3 wall blending should not be reintroduced as a quiet calibration parameter.

## Parameterization Strategy

### Seed Values

Use v27 fitted values as seeds for stable non-rear parameters. Do not seed from the collapsed v28 optimum except for structural changes such as the smooth gate form.

Carry forward:

- `A_Nu`, `B_Re`, `C_Pr`, with upper `A_Nu` bound relaxed.
- irradiance power scales as nuisance input factors.
- `G_core_perim`, `C_perim_eff`, `k_perim_ref`, and `beta_perim`.
- `f_core_tube` or equivalent source-to-hardware partition, if still structurally needed.
- `flange_scale`, but reinterpret it only as explicit hardware boundary scaling.
- `k_core_axial_scale`, with reporting if it remains near zero.

Remove from v29:

- `G_rear_core_heat`.
- `G_rear_core_cool`.
- `G_rear_perim` if it was only part of the distributed sink.
- `rear_sink_shape`.
- `f_T3_wall`.

Add for v29:

- `C_rear_eff`: rear/adaptor participating heat capacity.
- `G_core_rear`: conductance from final core cell or rear-face core region to rear reservoir.
- `G_perim_rear` or `f_rear_core`: perimeter participation in rear contact.
- `G_rear_tube`: rear reservoir to explicit tube/downstream hardware.
- `G_rear_cav`: rear reservoir to cavity/housing state, if physically represented.
- `gate_gain_off`: smooth cooling/off-state conductance gain.
- `gate_tau_off`: smooth gate onset time.

### Bounds and Priors

Bounds should be narrow enough to prevent a new empirical reservoir:

- `C_rear_eff`: estimate from ceramic/adaptor/holder geometry before fitting; use a bounded range around that estimate. If geometry is uncertain, start with an order-of-magnitude range such as 50-500 J/K and report the implied mass.
- Conductances: bound by plausible contact and conduction areas; if broad exploratory bounds are used, label v29a as diagnostic.
- `A_Nu`: allow values above v27's old ceiling, but report if the fitted value becomes implausibly high or boundary-active.
- Power scales: keep bounded and report as participating-flux calibration factors.
- Axial source depth: hold fixed in v29a unless the rear-reservoir test succeeds and residuals still demand source redistribution.

## Calibration Sequence

### Stage 0: Forward Consistency Test

Implement the new state and run with v27 non-rear parameters plus geometry-seeded rear parameters. This stage checks numerical stability and heat-accounting closure.

Required outputs:

- Energy closure table.
- Heat-path time integrals.
- T3, T10, T11 cooling traces.
- Bound/infeasible parameter report.

### Stage 1: Rear-Only Calibration

Fit only the new rear reservoir parameters and smooth gate parameters while holding convective/source/perimeter parameters close to v27 values.

Purpose:

- Determine whether the missing v28 heat path is genuinely a rear topology issue.

Pass condition:

- Major improvement over v28 without collapsing `C_perim_eff`, `G_core_perim`, `k_perim_ref`, or `k_core_axial_scale` to artificial bounds.

### Stage 2: Coupled Physical Calibration

Release the main physical parameters:

- `A_Nu`, `B_Re`.
- power scales.
- core/perimeter exchange.
- rear reservoir conductances.
- explicit hardware boundary scale.

Still hold source-depth freedom fixed unless Stage 1 and Stage 2 reveal a stable residual requiring it.

### Stage 3: Optional Source-Depth Test

Only after the rear reservoir has passed should source-depth or flux-penetration freedom be reopened. If reopened, it should be tested as an ablation:

- fixed front-local source,
- smooth exponential/Beer-like source depth,
- irradiance-dependent source depth only if optically justified.

## Acceptance Criteria

v29 should not be accepted simply because it improves the objective. It should meet most of the following:

- Objective is dramatically better than v28 and approaches v27 without restoring a direct receiver-volume sink.
- T3 heating residual improves while T3 remains pure gas.
- Cooling T3/T10/T11 do not show artificial late upturns.
- Fitted `C_rear_eff` implies a plausible participating rear mass.
- Rear heat-flow integral is plausible relative to absorbed input and gas uptake.
- `A_Nu` is not merely pinned to a new high upper bound.
- `G_core_perim`, `C_perim_eff`, and `k_perim_ref` remain physically interpretable and are not compensating for missing rear physics.
- Irradiance power scales remain interpretable as flux/input calibration factors.
- Residuals by irradiance and flow show improved structure, not only reduced aggregate error.
- Energy closure remains explicit and conservative.

## Diagnostics Required for v29

Add or extend diagnostics to include:

- Parameter table with value, units, bounds, and boundary-active flag.
- Heat-path integral table:
  - absorbed input,
  - gas uptake,
  - core-to-perimeter transfer,
  - core/perimeter-to-rear transfer,
  - rear-to-tube/cavity/flange transfer,
  - ambient/water losses,
  - storage terms.
- Residual table by thermocouple, irradiance, and flow rate.
- Heating/cooling split objective.
- Flow-slope diagnostic, especially for T3 and downstream thermocouples.
- Cooling monotonicity/upturn diagnostic.
- Participating heat-capacity table:
  - receiver measured mass basis,
  - perimeter effective capacity,
  - rear reservoir implied capacity and mass.
- Ablation comparison:
  - v27 baseline,
  - v28 negative control,
  - v29a rear reservoir without optional direct flange path,
  - v29b rear reservoir with gated reservoir-to-flange path if required.

## Scientific Credit Rules

The model should be described conservatively:

- v27: best empirical/effective baseline with split rear sink; useful but not final coefficient-validation model.
- v28: negative control showing that explicit tube/flange hardware alone is insufficient.
- v29: topology-identification model testing whether a bounded rear/adaptor reservoir can replace the empirical sink.

Only after v29 passes should the fitted convective/radiative/conductive coefficients be presented as candidate macroscopic receiver coefficients. If v29 still requires strong empirical gates, wall-blended T3, or boundary-active conductances, it should be reported as another diagnostic stage rather than a validated coefficient model.

## Recommended v29 Implementation Plan

1. Copy v28 into `1D_v29.jl` and runner/test files, but seed from v27 non-rear parameters.
2. Add `T_rear` as a new state with heat capacity and energy accounting.
3. Remove any direct receiver-volume rear sink from the RHS.
4. Connect final core/perimeter cells to `T_rear`.
5. Connect `T_rear` to explicit rear tube/cavity hardware.
6. Keep the smooth operating-state gate only on external hardware loss conductance.
7. Keep T3 as pure gas at the measured position.
8. Add rear-reservoir heat-path diagnostics before optimization.
9. Run smoke tests and a forward run before fitting.
10. Fit in staged sequence: rear-only, then coupled physical parameters, then optional source-depth ablation.
11. Update `summaries/journal.1D.md` with the v29 hypothesis, calibration outcome, and whether the rear reservoir replaced the v27 surrogate.

## Decision Before Coding v29

The recommended next version is:

```text
1D_v29a = v27 energy-accounting baseline
         - direct distributed rear sink
         - T3 wall blend
         + bounded rear/adaptor reservoir
         + smooth gated hardware boundary loss
         + pure-gas T3 mapping
```

The key question for v29a is not "can it fit better than v28?" It almost certainly should. The key question is:

```text
Can a physically bounded rear/adaptor reservoir recover most of v27's performance
without direct receiver-volume sinks or T3 wall blending?
```

If yes, v29 becomes a scientifically credible bridge toward coefficient validation. If no, the next hypothesis should be a measured-outlet/sensor submodel or a more explicit 2D-informed rear/manifold topology, not another hidden sink.
