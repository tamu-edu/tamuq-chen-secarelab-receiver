# Manuscript Readiness Assessment — 1D_v36, 1D_v37, 1D_v39

Assessed 2026-08-01 against the gate table in `1D_manuscript_gap_strategy.md`,
using the saved artifacts in `summaries/1D_v3{5,6,7,8,9}/` and the branch
narrative in `journal.1D.md` (§ "Multiphysics Branch Exploration").

**Note on sources:** a *Theory Manual v39* is not present in either connected
folder (`aysha/analysis`, `aysha/summaries`). The latest 1D theory manual on
disk is `1D_v28_theory_manual.md`. This assessment therefore reconstructs the
v36/v37/v39 formulations from the journal entry, the parameter CSVs, and the
saved diagnostics. If the v39 manual exists elsewhere in the repo, the
formulation-level conclusions below should be re-checked against it.

---

## 1. Headline verdict

| Version | Readiness | One-line |
|---|---|---|
| **1D_v36** (dynamic bypass) | **3 / 5 — citable with framing** | Only branch with invariant diagnostics; reproduces the target Nu exponent to 0.04%. Blocked by bound-locked capacities and a non-monotonic power convention. |
| **1D_v37** (optical redistribution) | **2 / 5 — supporting evidence only** | Better power convention, materially worse residuals than v36, no invariants computed. Its own mechanism undercuts the volumetric framing. |
| **1D_v39** (combined) | **1 / 5 — not citable as saved** | Saved diagnostics are ~200 K worse than the v35 baseline on every heating channel, which contradicts the reported objective of 0.629. Must be re-run before any v39 number appears in the manuscript. |

---

## 2. Gate-by-gate scoring

Invariant diagnostics (`invariants_*.csv`, `invariant_summary_*.csv`) exist
**only for v36**. v37, v38 and v39 wrote parameters, steady results, flow slopes
and residuals but no invariants. That alone prevents v37/v39 from being scored
on five of the ten gates.

| Gate | Target | v31 (ref.) | **v36** | **v37** | **v39** |
|---|---|---|---|---|---|
| Total heat capacity | 301 ± 23 J/K | 302.5 pass | **302.5 pass** | 302.5 pass | 302.7 pass |
| Apparent Nu exponent | Re^1.44 | 1.334 partial | **1.4394 pass** | not computed | not computed |
| Apparent Nu prefactor | 3.1e-4 | 2.23e-4 | **1.74e-4 fail** (1.8× low) | not computed | not computed |
| Nu trend quality | R² ≥ 0.97 | 0.961 | **0.935 partial** | not computed | not computed |
| Inversion threshold | ε* = 0.66 ± 0.03 | 0.876–0.983 fail | **0.885–0.983 fail** | 0.900–0.992 fail | 0.966–0.998 fail |
| Λ₁₀₇ | 0.038 + 8.3e-4·Re | wrong sign fail | **slope +5.7e-5, R² 0.083 fail** | not computed | not computed |
| External loss K_loss | 0.10–0.16 W/K | 0.250–0.341 fail | **0.150–0.227 partial** | not computed | not computed |
| Power convention | 456/304 elevated, 256 not collapsed | 0.636 fail | **1.369 / 1.484 / 0.629 fail** (non-monotonic) | 1.325 / 0.996 / 0.611 partial | 1.270 / 0.978 / 0.551 fail |
| Parameter interiority | no edge collapse | fail | **fail** (C_perim, C_rear, flange_scale, front_dep all at bounds) | fail (C_perim, C_rear, G_rear_tube = 0) | fail (m_rec = 0, G_rear_tube = 0, C's at bounds) |
| Heating/cooling simultaneous | no cooling-only patch | partial | **pass** (balanced, no switch) | partial | fail |

### Residual summary (mean signed steady error, K)

| Sensor | v35 base | **v36** | **v37** | v38 | **v39** |
|---|---:|---:|---:|---:|---:|
| T8 (11 mm) | −27.9 | **+13.8** | −49.8 | +481.8 | **−223.8** |
| T12_perim (58 mm) | −72.0 | **−40.1** | −97.1 | +527.9 | **−294.9** |
| T9_core (58 mm) | −65.6 | **−21.4** | −86.2 | +339.1 | **−265.9** |
| T11_perim (107 mm) | −31.5 | **+16.3** | −35.6 | +389.5 | **−236.8** |
| T10_core (107 mm) | −30.2 | **+28.7** | −22.7 | +233.9 | **−267.4** |
| T3 (gas, 140 mm) | −56.6 | **−0.1** | −24.9 | +102.9 | **−245.7** |
| T2 (cavity) | +31.1 | **+19.8** | +18.6 | +193.4 | **+37.5** |
| Cooling, worst channel RMSE | 30.8 | **49.9** | 50.9 | 55.9 | **62.7** |

v36 is the best-balanced heating fit ever produced in this family: every solid
channel is inside ±40 K and T3 is essentially exact, without a cooling-only
patch. That is the single strongest argument for using v36 as the manuscript
model.

---

## 3. Version-specific findings

### 1D_v36 — Dynamic bypass

**What earns readiness.** `Nu_app ∝ Re^1.4394` against a target of 1.44 is a
0.04% match, and it is obtained with a *physically motivated* mechanism
(temperature-dependent flow resistance, `m_rec = 0.0606`) rather than by
inflating a fit scalar. The total participating capacity closes at 302.5 J/K
against the 301 ± 23 J/K gravimetric bound, and K_loss has dropped from v31's
0.250–0.341 to 0.150–0.227 — the lower end is now inside the measured band.
Heating and cooling are satisfied simultaneously with no switch.

**What blocks it.**

1. **Bound-locked thermal inventory.** `C_perim_eff = 150.0` and
   `C_rear_eff = 80.0` sit *exactly* on their lower bounds, and `front_dep = 1.0`
   and `flange_scale = 0.1014` are also at edges. The capacity gate passes only
   because the bounds were placed to make it pass. A reviewer will read the
   inventory closure as imposed, not identified.
2. **Non-monotonic power convention.** `scale_304 = 1.484 > scale_456 = 1.369`,
   with `scale_256 = 0.629`. There is no optical or calorimetric story in which
   the 304 kW/m² case needs *more* power correction than the 456 kW/m² case.
   This is the clearest indication that a source term is still absorbing a
   structural error, and it is the item most likely to be caught in review.
3. **Λ₁₀₇ has no trend.** The slope sign is now correct (+5.7e-5 vs +8.3e-4)
   but it is 14× too small and R² = 0.083 — i.e. the model produces no
   Re-dependence in the deep-station lag at all. Any manuscript claim built on
   Λ₁₀₇ cannot be supported by v36.
4. **Depth flow sensitivity is 6–16× too steep.** At 456 kW/m², model
   dT₁₁/dṁ = −22.9 K/(L/min) against −1.4 experimental; dT₁₀/dṁ = −22.0 against
   −3.5. The experimental signature that deep stations go *flat* with flow is
   not reproduced.
5. **The crossover is still missing** (see §4).

**Framing that works:** present v36 as an *inverse-identification consistency
check* that independently recovers the experimentally-extracted Nu exponent,
with the bound-locked capacities disclosed. Do not present it as a validated
predictive coefficient-extraction model.

### 1D_v37 — Optical redistribution

The journal records v37 as a "massive success" on objective (0.2388), but the
saved residuals are worse than v36 on **every** heating channel — T9_core −86 K
vs −21 K, T12_perim −97 K vs −40 K, T8 −50 K vs +14 K. The objective ranking and
the physical residual ranking disagree, which means the v37 objective is being
carried by shape/transient terms rather than by steady accuracy.

Its one real advantage is the power convention: `1.325 / 0.996 / 0.611` is
monotonic in irradiance, unlike v36. That is worth preserving.

**The deeper problem is scientific, not numerical.** The fitted optics —
`front_dep = 0.571` with `beta_opt = 190 m⁻¹` — deposit 57% of the absorbed
power on the front face and the remainder within roughly 5 mm. That is a
*surface* receiver. Using it to fix the axial profile directly undercuts the
manuscript's volumetric-absorption framing, and 190 m⁻¹ has no independent
optical characterization behind it. v37 should be reported as a
sensitivity/hypothesis test that establishes how much of the residual is
attributable to source placement, not as a preferred formulation.

### 1D_v39 — Combined formulation

v39 is **not usable in its saved state**, and there is an internal
inconsistency that must be resolved first.

1. **Reported objective vs saved diagnostics disagree.** The journal records
   objective = 0.629, i.e. ~3× better than v35's 1.866. The saved
   `analysis_results_fitted_energy_accounting_1D_v39.csv` shows heating steady
   errors of −224 to −295 K on every solid channel and −246 K on T3 — roughly
   200 K *worse* than v35 on every channel. Either the diagnostics were written
   from a non-optimal parameter vector, or the v39 loss normalization differs
   from v35/v37. Until this is closed, no v39 number should be cited.
2. **v39 is not an independent combined fit.** `beta_opt = 190.23` is exactly
   v37's fitted value truncated to two decimals, indicating it was fixed/seeded
   rather than re-identified. The "degeneracy proof" in the journal is therefore
   a statement about a partly-frozen parameter space, not about the full
   combined space.
3. **The degeneracy claim is directionally right but over-stated.** `m_rec` did
   go to exactly 0.0, disabling dynamic bypass. But the journal's inference —
   that front-loaded optics "adequately covers for both the spatial and temporal
   errors" — is contradicted by the same run's residuals, which are the worst
   heating errors in the series. The defensible version of the claim is: *the
   optimizer sacrificed the transient mechanism and failed to recover the steady
   field*, i.e. the two mechanisms are not additive, not that optics subsumes
   flow.
4. **Regressions on gates v31 was supposed to fix.** `G_rear_tube` returns to
   0.0 (floored to 0.5), and `scale_256 = 0.551` is the deepest low-flux power
   collapse in the family. `ε_mean = 0.966–0.998` and `Nu` up to 17.4 are the
   furthest any version has been from the ε* = 0.66 inversion threshold.
5. One genuine positive: v39's flow slopes are the closest to experiment at
   depth (T₁₀ −1.85 vs −3.54; T₁₂ −16.7 vs −16.8 at 456 kW/m²). The *shape* of
   the flow response improves even as the absolute level collapses — worth
   diagnosing rather than discarding.

---

## 4. The result that is actually manuscript-grade

Two failures are shared by v35, v36, v37 and v39 and are more valuable than any
of the individual fits:

**(a) The high-flow thermal crossover is never reproduced.**

| Case | I (kW/m²) | ṁ (L/min) | T9−T8 model | T9−T8 exp |
|---|---:|---:|---:|---:|
| E67 | 456 | 15.28 | −9 (v36) / −26 (v39) | **+35** |
| E72 | 304 | 18.71 | −13 (v36) / −19 (v39) | **+33** |
| E77 | 256 | 13.85 | −12 (v36) / −16 (v39) | **+37** |
| E71 | 456 | 7.13 | −74 (v36) | −108 |
| E76 | 304 | 4.53 | −176 (v36) | −120 |

Every formulation gets the low-flow ordering roughly right and the high-flow
sign wrong. No scalar closure — empirical Nu, developing-flow Nu, exchange-area
fraction, static bypass, dynamic bypass, two-stream gas, or front-loaded optics
— produces an internal station hotter than the front face at high flow.

**(b) Cell effectiveness never approaches the inversion threshold.**

Fitted ε_mean is 0.885–0.983 (v36), 0.900–0.992 (v37), 0.966–0.998 (v39),
against the experimentally-inferred ε* = 0.66 ± 0.03. With ε = 1 − exp(−NTU),
ε = 0.98 implies NTU ≈ 3.9 while ε* = 0.66 implies NTU ≈ 1.08 — the model needs
roughly a **3.6× reduction in UA/(ṁ·c_p)** to reach volumetric inversion. Neither
optical redistribution nor dynamic bypass moves this by more than a few percent,
because both leave the *product* of exchange coefficient and swept capacity
essentially intact. This is the quantitative statement of the "Advective
Bottleneck" and it is the strongest physics result in the branch.

Together (a) and (b) constitute a clean inverted proof: the LTNE 1D continuum,
with any closure tested here, is structurally incapable of volumetric inversion
at the measured effectiveness. That is a defensible manuscript conclusion and
does not depend on which of v36/v37/v39 fits best.

---

## 5. Recommended path to manuscript

1. **Adopt v36 as the quantitative 1D model**, with v37 and v38 reported as the
   two falsified alternative hypotheses and v39 as the non-additivity test.
2. **Re-run v39 diagnostics before citing anything from it.** Confirm whether
   the saved CSV corresponds to the optimum, and whether the v39 loss uses the
   same normalization as v35/v37. This is the highest-priority action item.
3. **Fix v36's power convention.** Re-fit with a monotonicity constraint
   `scale_456 ≥ scale_304 ≥ scale_256` and report the objective penalty. If v36
   survives with the Nu exponent intact, the result becomes substantially harder
   to attack.
4. **Re-open the capacity bounds** (`C_perim_eff`, `C_rear_eff`) by ±30% and
   re-fit. If 302.5 J/K is recovered from the interior, the inventory gate
   becomes a genuine validation rather than an imposed one.
5. **Compute invariants for v37 and v39.** Five gates are currently unscorable
   for two of the three versions under assessment.
6. **Lead the modelling section with the inverted proof (§4)**, not with a fit
   quality table. The crossover failure and the ε gap are reproducible,
   mechanism-independent, and quantitative; the fit rankings are neither.

## 6. Open items to verify

- Locate or regenerate the *Theory Manual v39*; the formulation-level statements
  above are inferred from parameter CSVs and the journal, not from the manual.
- Confirm whether `beta_opt` was fitted or fixed in v39, and whether `front_dep`
  in v36 (= 1.0) is a bound or a project constant
  (`FRONT_DEPOSITION_FIXED_V5 = 1.0` suggests the latter).
- Reconcile the v37 objective ranking (0.2388, better than v36's 0.2173 claim
  ordering in the journal narrative) against v37's uniformly worse steady
  residuals — the loss weighting should be documented in the manuscript.
