# Strategy: Bringing the 1D and 2D Models to Manuscript-Incorporation Readiness

## 0. Audience and Intent

This note is written to direct the modeling agent(s) working on `1D_vNN.jl` and
`2D_vNN.jl`. It sits above `1D_v29_strategy_forward.md`: that document governs
the next 1D coding step, this document defines what "done" means relative to the
manuscript under preparation and what each model must deliver before any of its
fitted coefficients can appear in the paper.

The manuscript is currently a **model-free experimental characterization** paper.
Its strength is that its results are model-independent. Modeling must not
dilute that. The models have exactly two admissible roles in the paper:

- **Role A (available now): external corroboration.** Support the delivered-power
  accounting and the qualitative structural claims, without reporting fitted
  coefficients as physical.
- **Role B (the target): validated coefficient identification.** Report the
  effective macroscopic heat-transfer coefficients as candidate physical values.

No model output may be presented in Role B until it clears the acceptance gate in
Section 3. Both theory manuals already say this; this strategy makes the gate
explicit and testable.

---

## 1. Current Status (as of 2026-07-23)

**1D — mid rear-topology identification.**
- `v27` (split heating/cooling rear sink): best empirical fit, objective ≈ 0.175,
  but the direct rear sink is an acknowledged surrogate → weak scientific credit.
- `v28` (direct sink removed, explicit tube/flange only): negative control,
  objective ≈ 117.7, parameters collapsed to bounds. Falsifies the
  hardware-path-only hypothesis.
- `v29a` (bounded rear/adaptor reservoir): objective ≈ 53.1;
  `C_total_with_rear` ≈ 336 J/K vs 301 J/K reference. Better than v28, far from
  v27, T3 still biased low, several conductances bound-active.
- Durable result: **two-zone core/perimeter macro-topology with radial coupling
  is necessary.** One-zone models are closed out. Coefficients not yet validated.

**2D — genuine axisymmetric continuum LTNE Macro-ECM.**
- `2D_v2`: 446-state finite-volume LTNE model with multi-domain geometry,
  anisotropic `k_s,r^eff`/`k_s,z^eff` (+Rosseland axial term), developing-flow Nu,
  Churchill-Chu front plume, T-dependent felt conductivity, radial preferential
  flow.
- Performance: **cooling fits well** (single-digit to ~tens of K steady error),
  **heating overpredicts substantially** (steady error T8 +217 K, T9 +84 K,
  T2 +77 K at E81). Several parameters bound-active (`chi_z`=1.0, `m_rec`=0,
  `s_flange`=1.0). `C_cavity` fitted essentially to the 301 J/K target.
- This is the correct model *class* for the paper's Section 5 claims, but it is a
  formulation milestone, not a converged result.

**Theory manuals — current and honest.** `1D_v28_theory_manual.md`,
`2D_v2_theory_manual.md`, `1D_v29_strategy_forward.md` are up to date and
correctly refuse to present bound-active fitted values as physical. Keep this
discipline.

---

## 2. The Validation Target: Reproduce the Manuscript's Measured Invariants

The manuscript already owns a small set of model-independent, quantitative
results. These are the "portable constraints on any future model" named in the
Discussion. **A model earns Role B only by reproducing these as emergent
outputs, with the corresponding prior/regularization released.** This is the
single most important instruction in this document.

| # | Measured invariant | Value (from manuscript) | How the model must reproduce it |
|---|---|---|---|
| I1 | Apparent global Nusselt | `Nu = 3.1e-4 Re^1.44`, 23 ≤ Re ≤ 94, r²=0.97, 15–100× below duct theory | Post-process the model's steady solid/gas fields into the *same* apparent global `h` definition used in the paper (instrumented side-wall temperature, full nominal exchange area), regress `Nu` vs `Re`, compare prefactor and exponent. The exponent > 1 is the key signature (flow-dependent solid recruitment). |
| I2 | Volumetric-inversion threshold | `epsilon* = 0.66 ± 0.03`, flux-independent | Model must show the axial wall-temperature peak migrating front→interior when gas-side effectiveness crosses ε*. Report the model's crossover ε* and its flux-independence. |
| I3 | Deep nonequilibrium slope | `Lambda_107 = 0.038 + 8.3e-4 Re` (interior-to-wall deficit at 107 mm, linear in Re, ~zero intercept) | Compute the model's `T_core(107) − T_perim(107)` deficit vs Re; match slope and near-zero intercept. Present state: 1D two-zone matches sign/trend; 2D underpredicts in heating. This is the sharpest discriminator. |
| I4 | Effective capacitance | `C_eff = 301 ± 23 J/K`, ≈ 7× monolith | Release the `C_ref`/`C_cavity` prior (set `w_C = 0` / widen bounds) and check that the *fitted* participating capacity still lands near 301 J/K. Only then is C_eff a model-independent confirmation. Today both models ingest 301 as a prior → not independent. |
| I5 | Parasitic loss conductance | `K_loss = 0.10–0.16 W/K` | Report the model's total ambient/flange loss conductance in the same excess-temperature basis and confirm it falls in this band. |
| I6 | Delivered-power factors | +30–60% over nominal at 456/304 kW·m⁻²; ≈ nominal or less at 256 | Already used as corroboration. Reconcile 1D and 2D (Section 4, fix F2). |

**Rule:** invariants that a model currently *ingests* (I4 via the 301 J/K prior;
I6 via fitted power scales) cannot be claimed as independent confirmations until
the prior is released and the value re-emerges. State this explicitly in every
results table.

---

## 3. Acceptance Gate for Role B (validated coefficients)

A model may report fitted coefficients as physical only if **all** of the
following hold. Report each as an explicit pass/fail line.

1. Fit quality: objective within a stated factor of the best empirical baseline
   **without** any surrogate sink, hard operating switch, hidden absorptivity,
   or T3 wall-blend.
2. Parameter interiority: no coefficient intended as physical (`A_Nu`, `B_Re`,
   `G_core_perim`, `C_perim_eff`, `k_perim_ref`, `k_s,r^eff`, `k_s,z^eff`,
   `beta_opt`) is bound-active. Flag every bound-active parameter.
3. Invariant reproduction: I1, I2, I3 reproduced as emergent outputs (not
   ingested), within stated tolerances; I4 re-emerges near 301 J/K with the prior
   released; I5 in band.
4. Heating **and** cooling both satisfied simultaneously — no regime-specific
   patch that helps one and damages the other.
5. Energy closure explicit and conservative (input = gas uptake + storage +
   losses), reported per case.

Until the gate is cleared, the model stays in Role A and its coefficients are
described as "candidate/effective, not yet validated."

---

## 4. Three Consistency Fixes (required before either model is elevated)

**F1 — C_eff must not be a circular confirmation.** Both models currently take
the manuscript's 301 J/K as a soft prior (`1D C_ref=301, w_C=0.10`;
`2D C_cavity` fitted to 301). The paper cannot claim the models "independently
confirm" C_eff while feeding them C_eff. Action: run each model with the C prior
released (I4) and report the fitted value. If it re-emerges near 301, that is a
real confirmation; if not, drop the confirmation claim.

**F2 — 1D and 2D disagree on absolute power scale.** 1D fits `scale_456 ≈ 2.0`
with `eta_abs = 1`; 2D fits `f_456 = 1.336` with `eta_abs = 0.85` (effective
≈ 1.14). The manuscript's model-calibration column uses the 2D numbers
(1.336 / 1.374 / 0.786). Action: adopt one absorptivity/spill convention across
both model families, restate both sets of scales in that convention, and confirm
they agree within the ±8% calibration scatter. Document the convention in both
theory manuals.

**F3 — T3 is an unresolved outlet measurement in both models.** 1D v29 leaves T3
biased low; the 2D exit-gas T3 rests on the unconverged heating field. Action:
do not build any manuscript claim on T3 until the rear/outlet topology (1D) and
heating convergence (2D) are resolved. Keep T3 as pure gas until a physically
justified sensor submodel is warranted (per `1D_v29_strategy_forward.md`).

---

## 5. 1D Workstream

Goal: convert the two-zone macro-topology into a validated coefficient model by
closing the rear/outlet topology.

1. Execute `1D_v29_strategy_forward.md` staged plan (rear/adaptor reservoir,
   staged calibration, diagnostics). Do not restore the v27 direct sink.
2. If v29a's single terminal-cell reservoir remains insufficient (current
   evidence), proceed to **v29b: a distributed rear-contact / manifold topology**
   or an explicitly labeled outlet-sensor submodel — not another hidden sink.
3. Target: recover v27-level fit (~0.175) with all rear/perimeter parameters
   interior, then release the C prior (I4) and derive I1/I2/I3 as emergent
   outputs.
4. Deliverable: a coefficient table (`A_Nu`, `B_Re`, `G_core_perim`,
   `C_perim_eff`, `k_perim_ref`, rear conductances) each tagged interior/bound,
   with the invariant-reproduction report.

---

## 6. 2D Workstream

Goal: converge the heating calibration and demonstrate the continuum LTNE model
reproduces the paper's invariants.

1. Diagnose the heating overprediction (T8 +217 K). Likely suspects, in order:
   front/rim spillage partition and front-plume loss magnitude; radial effective
   conductivity `chi_r` splitting heat outward too weakly/strongly; power-scale ×
   `eta_abs` convention (see F2). The cooling fit being good while heating is far
   off points to a source/loss-magnitude problem, not a thermal-mass problem.
2. Free the bound-active parameters (`chi_z`=1.0, `m_rec`=0, `s_flange`=1.0) with
   physically justified bounds and confirm they leave the bounds; if they don't,
   report and interpret.
3. Reproduce I2 (ε* peak migration) and I3 (Λ₁₀₇ slope) directly from the 2D
   field — these are the model's natural strengths and the strongest validation.
   The 2D model failing to grow the center/wall split downstream has been the
   persistent structural symptom; treat I3 as the primary convergence metric.
4. Release the `C_cavity`=301 prior (I4) and report the emergent assembly capacity.
5. Deliverable: validated anisotropic `k_s,r^eff`, `k_s,z^eff`, `Nu(Re,Pr,z)`
   closure, `beta_opt`, felt `k(T)`, and front `h(T)` — each with interior/bound
   status and the invariant-reproduction report.

---

## 7. Theory-Manual and Scientific-Credit Requirements

- Every theory manual must carry a **Validation Status** section stating which of
  I1–I6 are reproduced-as-output vs ingested-as-prior vs not-yet, and which
  parameters are bound-active.
- Never present a bound-active fitted value as a physical coefficient.
- Maintain the ablation record (v27 baseline / v28 negative control / v29 rear
  reservoir; 2D_v1 / 2D_v2) so reviewers can see the falsification chain.
- Keep power scales described as participating-flux calibration factors, not
  absorptivity.

---

## 8. What to Hand Back for the Manuscript

Two possible deliverables, depending on which gate is cleared:

- **If Role A only (near-term):** a short "Model corroboration" supplementary
  note citing both theory manuals, stating precisely (a) that the two-zone and
  LTNE structures independently support the assembly-scale / nonequilibrium
  claims of Sections 5.1 and 5.3, (b) the reconciled delivered-power factors
  (F2), and (c) an explicit statement that fitted coefficients are not yet
  validated and C_eff is a model input, not a model confirmation (F1).

- **If Role B cleared (target):** validated coefficient tables for 1D (two-zone
  effective transport) and/or 2D (anisotropic continuum coefficients), plus the
  invariant-reproduction report (I1–I6), suitable for a modeling Results section
  or a companion modeling paper. Do not merge into the experimental paper in a
  way that makes any experimental result depend on a fitted coefficient.

---

## 9. One-Line Summary for the Modeling Agent

Do not optimize the objective further in isolation. The next milestone is to make
each model **reproduce the manuscript's measured invariants (Nu–Re exponent, ε*,
Λ₁₀₇ slope, C_eff with the prior released) as emergent outputs, with all physical
parameters interior and heating+cooling satisfied together.** That, not a lower
loss, is what promotes a coefficient from "effective/fitted" to "validated" and
into the paper.
