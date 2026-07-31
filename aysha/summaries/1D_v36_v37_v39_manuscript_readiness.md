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

**Provenance of the gate targets — read this first.** None of the targets used
below are literature values or general porous-media correlations. All are
identified from *this* campaign and are stated in `claude_data_analysis_sections_v2.md`:

| Target | Source | What it is |
|---|---|---|
| Nu = (3.1 ± 0.1)e-4 · Re^1.44, r² = 0.97 | R3 | Apparent global exchange, 15 runs, 23 ≤ Re ≤ 94 |
| ε* = 0.66 ± 0.03 | R4 | Flux-independent effectiveness at which the wall peak moves inside |
| Λ₁₀₇ = 0.038 + 8.3e-4·Re, r² = 0.90 | R5 | Deep-station gas–solid depression (a stated *lower bound*) |
| C_eff = 301 ± 23 J/K, K_loss = 0.10–0.16 W/K | R7 / M4 | Slow-mode eigenvalue identification from cooling decays |
| Envelope: ε = 0.57–0.78, NTU = 0.85–1.51, Nu_app = 0.028–0.212 | R1 | Measured operating range |

They are therefore legitimate *validation* targets — the model should reproduce
the receiver's own measured behaviour — but they are **not** physical bounds
constraining the formulation, and nothing here should be read as importing an
external correlation onto a receiver geometry that has none.

**They must also be compared like-for-like.** ε and Nu_app are *reduced*
quantities with specific definitions (M2):

$$\varepsilon = \frac{T_3 - T_{amb}}{\bar T_w - T_{amb}},\qquad
\bar T_w = 0.248\,T_8 + 0.365\,T_{12} + 0.387\,T_{11},\qquad NTU = -\ln(1-\varepsilon)$$

with T9/T10 explicitly excluded (radiation-biased interior probes), and the
inversion indicator defined on the wall chain as `I_vol = T12 − T8`. The model's
saved `mean_effectiveness` column is a *per-cell* effectiveness and is **not**
the same quantity; the same applies to `mean_nu` versus Nu_app. Section 4 below
recomputes both from the model output using the manuscript definitions.

---

## 1. Headline verdict

| Version | Readiness | One-line |
|---|---|---|
| **1D_v36** (dynamic bypass) | **3.5 / 5 — citable with framing** | Only branch with invariant diagnostics; reproduces the measured ε envelope (0.583–0.777 vs 0.573–0.781) and the Nu exponent to 0.04%. Blocked by bound-locked capacities, a non-monotonic power convention, and an inversion threshold shifted to ε ≈ 0.75. |
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
| ε envelope (0.57–0.78) | manuscript defn. | — | **0.583–0.777 pass** | 0.656–0.770 partial (range compressed) | 0.351–0.488 fail |
| Inversion threshold, ε at I_vol sign change | 0.66 ± 0.03, flux-independent | — | **0.72–0.78, flux-dependent, no crossing at 256 — fail** | 0.73 (456 only), no crossing at 304/256 — fail | no crossing — fail |
| Λ₁₀₇ | 0.038 + 8.3e-4·Re | wrong sign fail | **slope +5.7e-5, R² 0.083 fail** | not computed | not computed |
| External loss K_loss | 0.10–0.16 W/K | 0.250–0.341 fail | **0.150–0.227 partial** | not computed | not computed |
| Power convention (vs R6 closure, §7) | f = 1.34 / 1.58 / 1.11 | 0.636 | **1.369 / 1.484 / 0.629 — 2 of 3 within 6%, 256 group 43% low** | 1.325 / 0.996 / 0.611 — only 456 agrees | 1.270 / 0.978 / 0.551 — none agree |
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
2. **Power convention — criticism retracted, one group still off.** An earlier
   draft of this file called `scale_304 = 1.484 > scale_456 = 1.369`
   "physically incoherent." That was wrong: R6's model-free delivered-power
   closure produces the *same* non-monotonic pattern from an independent route
   (see §7). v36 agrees with it to 2.2% at 456 and 6.1% at 304. Only the 256
   group remains a genuine outlier, at 43% below the R6 estimate.
3. **Λ₁₀₇ has no trend.** The slope sign is now correct (+5.7e-5 vs +8.3e-4)
   but it is 14× too small and R² = 0.083 — i.e. the model produces no
   Re-dependence in the deep-station lag at all. Any manuscript claim built on
   Λ₁₀₇ cannot be supported by v36.
4. **Depth flow sensitivity is 6–16× too steep.** At 456 kW/m², model
   dT₁₁/dṁ = −22.9 K/(L/min) against −1.4 experimental; dT₁₀/dṁ = −22.0 against
   −3.5. The experimental signature that deep stations go *flat* with flow is
   not reproduced.
5. **The inversion threshold is displaced and flux-dependent** (see §4b) — v36
   needs ε ≈ 0.72–0.78 to move the wall peak inside, against a measured
   flux-independent 0.66, and never inverts at all in the 256 kW/m² group.
   Since R4 is one of the paper's headline claims, this is the most
   scientifically consequential of v36's gaps, and it is *not* an effectiveness
   magnitude problem — v36 gets ε right (§4a).

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

## 4. Recomputed with the manuscript definitions

An earlier draft of this assessment compared the model's per-cell
`mean_effectiveness` (0.885–0.998) against the measured global ε* = 0.66 and
concluded the models were nowhere near the inversion regime. **That was a
category error** — two different quantities. Recomputing ε from the model output
using M2's definition, ε = (T3 − T_amb)/(T̄_w − T_amb) with the trapezoid wall
weights, gives a very different and much more interesting picture.

### (a) v36 reproduces the measured ε envelope almost exactly

| | ε range | NTU range |
|---|---|---|
| **Experiment (R1)** | 0.57–0.78 | 0.85–1.51 |
| Experiment (recomputed here) | 0.573–0.781 | 0.85–1.52 |
| **v36** | **0.583–0.777** | **0.88–1.50** |
| v37 | 0.656–0.770 | 1.07–1.47 |
| v35 | 0.520–0.732 | 0.73–1.32 |
| v39 | 0.351–0.488 | 0.43–0.67 |

v36 lands inside the measured envelope on both ends and tracks it run-by-run
(E67 0.743/0.735, E72 0.777/0.773, E77 0.775/0.781, E76 0.583/0.573). This is a
strong, previously unreported validation result and it materially raises v36's
standing. v37's range is compressed at the low-ε end — it cannot get down to the
0.57 low-flow points. v39 is uniformly ~0.25 too low, consistent with its
collapsed temperature field.

### (b) The real failure is the *location* of the inversion threshold

R4's claim is that the wall peak moves inside the receiver when ε exceeds
≈0.66, **independent of flux**. Recomputing the ε at which `I_vol = T12 − T8`
changes sign:

| Flux (kW/m²) | Experiment | v36 | v37 |
|---:|---|---|---|
| 456 | 0.644–0.673 | 0.720–0.743 | 0.728–0.732 |
| 304 | 0.625–0.668 | 0.747–0.777 | no crossing |
| 256 | 0.628–0.679 | **no crossing** | no crossing |

The experimental crossings bracket 0.66 in all three groups — the flux
independence in R4 is clean. v36 inverts, but only at ε ≈ 0.72–0.78, and the
threshold drifts with flux; at 256 kW/m² it never inverts at all (max I_vol
= −9.5 K where the experiment reaches +57.8 K).

So the correct statement is not "the model can't invert." It is:

> The model reaches the measured effectiveness envelope but requires roughly
> 0.08–0.10 more effectiveness than the real receiver to move the wall peak
> inside, and its threshold is flux-dependent where the measurement's is not.

That is a much sharper and more diagnosable claim. It says the deficiency is in
how axial *source* and axial *exchange* are distributed relative to each other —
at fixed ε the model puts too much of its wall temperature at the front face —
rather than in the overall magnitude of gas–solid coupling, which v36 gets
right.

### (c) What this does to the "Advective Bottleneck" argument

The journal's §"Final Proof" concludes that no scalar flow structure can satisfy
heating and cooling simultaneously. v36's ε agreement weakens the strong form of
that claim: a *dynamic* (temperature-dependent) flow structure with `m_rec` =
0.0606 does satisfy both, and lands the effectiveness envelope. What it does not
do is relocate the inversion threshold. The defensible version of the argument
is therefore narrower than the journal states, and should be rewritten before it
goes into the manuscript.

### (d) Caveat on the remaining reduced quantities

Given that ε required correction, the `Nu_app_prefactor` (1.74e-4 vs 3.1e-4) and
`Nu_app_Re_exponent` (1.4394 vs 1.44) in `invariant_summary_1D_v36.csv` should
be verified to use M2's reduction — h_app from the same ε·ṁc_p route and the
same wetted-area convention — before either is cited. The exponent match to
0.04% is striking enough that it deserves that check. Note that `mean_nu`
(3.07–8.50) is a local channel Nusselt number and is three orders of magnitude
above Nu_app (0.028–0.212); they are not comparable and neither should be
substituted for the other.

---

## 5. Recommended path to manuscript

1. **Adopt v36 as the quantitative 1D model**, with v37 and v38 reported as the
   two falsified alternative hypotheses and v39 as the non-additivity test.
2. **Re-run v39 diagnostics before citing anything from it.** Confirm whether
   the saved CSV corresponds to the optimum, and whether the v39 loss uses the
   same normalization as v35/v37. This is the highest-priority action item.
3. ~~Fix v36's power convention with a monotonicity constraint.~~ **Superseded
   by §7** — do *not* impose monotonicity on the scale factors; R6 shows the
   real delivered-power factors are non-monotonic. Fix the aperture fluxes to
   the independent estimate instead and re-fit the remainder.
4. **Re-open the capacity bounds** (`C_perim_eff`, `C_rear_eff`) by ±30% and
   re-fit. If 302.5 J/K is recovered from the interior, the inventory gate
   becomes a genuine validation rather than an imposed one.
5. **Compute invariants for v37 and v39.** Five gates are currently unscorable
   for two of the three versions under assessment.
6. **Lead the modelling section with §4, not with a fit-quality table.** The
   pairing is the story: v36 matches the measured ε envelope run-by-run yet
   needs ~0.1 more ε to invert, and loses flux independence doing it. That is a
   quantitative, well-posed statement about axial source/exchange placement.
7. **Rewrite the "Advective Bottleneck" conclusion** in the journal before it
   reaches the manuscript — v36 falsifies its strong form (§4c).
8. **Do not import external Nusselt correlations as validation targets.** The
   receiver's apparent Nu is 15–100× below fully developed duct theory (R3); a
   monolith-channel correlation from the literature is a *contrast*, not a
   bound. The literature in `analysis/literature/` (Cornejo & Hayes 2020 on
   monolith Nu, Fend 2004/2013, Avila-Marin 2011/2019, Kribus 2014) has not been
   read into this assessment and should be used only to position the result, not
   to constrain the model.

## 7. Fixing the irradiance basis instead of fitting it

**Planned change:** replace the three fitted `scale_*` parameters with
independently estimated aperture irradiances, and restrict the calibration to
the remaining parameters. The power-scale fitting is not presented explicitly,
on the grounds of a known experimental flux discrepancy.

**This is the single best available move, and it is already supported.** R6 of
`claude_data_analysis_sections_v2.md` derives a per-flux delivered-power factor
with no thermal model in the loop —
f = [Q_gas + K·(T̄_w − T_amb)]/(G₀·A_frt), with the cooling-identified K_loss
bracket (0.096–0.164 W/K) carried as a systematic band:

| G₀ [kW/m²] | f (K_cool) | f (K_heat) | **v36 fitted** | v36 vs f(K_heat) |
|---:|---:|---:|---:|---:|
| 456 | 1.05 | 1.34 | **1.369** | **+2.2%** |
| 304 | 1.23 | 1.58 | **1.484** | **−6.1%** |
| 256 | 0.84 | 1.11 | **0.629** | **−43.4%** |

Two consequences:

1. **The non-monotonic pattern is corroborated, not an artifact.** R6's own
   f values peak at 304 (1.58 > 1.34 > 1.11). The optimizer independently found
   the same ordering. Reviewer risk on this point is far lower than it appeared.
   Checklist item 3 in the draft already commits to reporting Io·A_frt as a
   *lower bound* — the modelling section becomes a second, independent
   confirmation of exactly that statement.
2. **The 256 group is the real discrepancy.** v36 wants 161 kW/m² effective
   against R6's 215–284. Fixing the flux there adds 26–76% more power to that
   group.

**Why the validation targets are safe under this change.** ε, NTU, Nu_app,
Λ₁₀₇, C_eff and K_loss are all constructed from temperatures, metered flow and
gas properties — none contains an irradiance term. Changing the flux basis moves
the *model*, not the targets. The comparison therefore stays valid, and becomes
stricter, since three degrees of freedom are removed.

**Predicted effect on the open gaps.** The 256 group is currently the weakest in
v36: it is under-predicted (E77 T8 496 K vs 566 K; E81 714 K vs 760 K) and it is
the one group that never inverts (max I_vol = −9.5 K against +57.8 K measured).
Both defects point the same way, and raising that group's power is the correct
direction for both. There is a real chance that fixing the flux basis closes
§4b — the most consequential remaining gap — rather than costing anything. That
prediction should be checked explicitly, not assumed.

**What must still be re-verified after the refit**, since every headline number
is conditional on the fitted source:

```text
1. Nu_app exponent (currently 1.4394 vs 1.44) — does it survive?
2. ε envelope (0.583-0.777 vs measured 0.573-0.781) — does it survive?
3. Inversion threshold per flux group, especially 256 — does it appear?
4. C_total_with_rear — still 301 +/- 23 J/K without bound-locking?
5. Whether C_perim_eff / C_rear_eff leave their lower bounds once three
   source parameters no longer compete with them.
```

**Disclosure requirement.** Not presenting the *fitting* is fine. Not presenting
the *values* is not. The modelling section needs one line stating the aperture
fluxes actually used, their provenance (R6 closure, or whatever independent
estimate replaces it) and their uncertainty band. Note that the draft's header
currently reads "All irradiance values used in the calculations are the nominal
measured aperture fluxes (456, 304, 256 kW/m²); no correction factors of any
kind are applied." That statement is true and important for the *experimental*
reduction — and it must be scoped to it, or it will read as false once the
modelling section uses corrected fluxes. A reviewer who finds an undisclosed
flux correction will discount everything else in the paper; a reviewer who finds
a disclosed one, corroborated by R6, will read it as a strength.

## 6. Open items to verify

- Locate or regenerate the *Theory Manual v39*; the formulation-level statements
  above are inferred from parameter CSVs and the journal, not from the manual.
- Confirm whether `beta_opt` was fitted or fixed in v39, and whether `front_dep`
  in v36 (= 1.0) is a bound or a project constant
  (`FRONT_DEPOSITION_FIXED_V5 = 1.0` suggests the latter).
- Confirm that `Nu_app_*` in the invariant diagnostics uses M2's reduction
  (§4d). One reduced-quantity mismatch has already been found in this file.
- Check whether the model's `mean_effectiveness` diagnostic should simply be
  replaced by the manuscript-definition ε in future runners, to prevent the
  per-cell/global confusion recurring.
- Reconcile the v37 objective ranking (0.2388, better than v36's 0.2173 claim
  ordering in the journal narrative) against v37's uniformly worse steady
  residuals — the loss weighting should be documented in the manuscript.
