# Project status and forward plan — 28 August 2026

**Revision 3** (rev. 1 preceded the 1D_v41 full run; rev. 2 documented its
source-scale failure; rev. 3 incorporates the documentation and implementation
audit). Restart briefing after a ~4 week pause. Reconstructed from
`journal.1D.md`, `journal.2D.md`, `1D_v36_v37_v39_manuscript_readiness.md`
(1 Aug), `manuscript_2D_readiness_strategy.md`, `1D_manuscript_gap_strategy.md`,
`analysis/manuscript/manuscript_full_draft.md`, and the saved CSVs under
`summaries/1D_v3{5,6,7,8,9}/`, `1D_v4{0,1}/` and `2D_v2{0,1,2}/`.

**What changed in rev. 3:** the manuscript addendum was removed from after the
table captions, conservative Role-A language was moved into the Discussion,
both journals received superseding audit entries, and the v22 Phase-4 wiring
was checked against its saved grid. The experimental power discussion now
treats the T3-based closure and model-derived factors as provisional rather
than independent validation. Temperature-only results are unchanged.

---

## 1. Where the three tracks stand

### 1D — a preferred model, and a failed test that located the real problem

| Version | Mechanism | Objective | Verdict |
|---|---|---:|---|
| v35 | baseline, scalar bypass | 1.866 | superseded; "Advective Bottleneck" argument comes from here |
| **v36** | **dynamic bypass, phi(T_core)** | **0.2173** | **best balanced saved fit; Role-A manuscript candidate** |
| v37 | optical redistribution | 0.2388 | source-placement sensitivity; no identified optical coefficient |
| v38 | two-stream gas | 5.701 | rejected within the tested parameterization |
| v39 | v36+v37 combined | 0.629 reported | **not citable** — saved residuals contradict the objective |
| v40 | v39 variant, m_rec active | 1.476 | undocumented; worse than v36 on nearly every gate |
| **v41** | **R6 fluxes frozen, no dynamic bypass** | **5780.6** | **did not converge; see §2 — the most useful failure so far** |

**v36 remains the strongest result the campaign has produced.** Mean steady
heating error inside ±40 K on every solid channel (T8 +13.8, T9 −21.4,
T12 −40.1, T11 +16.3, T10 +28.7 K), T3 essentially exact (−0.1 K), heating and
cooling satisfied simultaneously with no switch. Nu exponent 1.4394 against the
measured 1.44. The measured effectiveness envelope is reproduced run-by-run
(0.583–0.777 vs 0.573–0.781). C_total 302.5 J/K against 301 ± 23.

**What still blocks v36:**
- `C_perim_eff = 150.0` and `C_rear_eff = 80.0` sit *exactly* on their lower
  bounds; `flange_scale` and `front_dep` also at edges. The capacity gate
  passes because the bounds were placed to make it pass.
- Nu prefactor 1.74e-4 vs 3.1e-4 (1.8x low).
- Lambda_107 slope 5.7e-5 vs 8.3e-4, r2 = 0.083 — no Re-trend at all.
- The inversion threshold sits at eps ~ 0.72–0.78 and drifts with flux,
  against a measured flux-independent 0.66; at 256 kW/m2 it never inverts.
- The power-convention gap is **not** what it was thought to be — see §2.

**v39/v40 remain unresolved loose ends.** v39's saved CSVs show mean steady
errors of −180 to −243 K on every channel, contradicting its reported objective
of 0.629. v40 flips the sign (+236 K at T8) with no journal entry and no theory
manual. Neither is documented.

### 2D — v20 closed for Role B; v22 result invalid pending wiring repair

The v20 gate-free stress test is the trustworthy endpoint. It confirms that
its rejection is not caused by solver failure or acceptance penalties (maximum
enthalpy residual about 6e-11). Receiver `UA(Re)` remains non-identifiable
without an independent bulk outlet-enthalpy measurement or verified T3
transfer function. The standing v20 decision is to stop Role-B coefficient
extraction from that inverse formulation and retain it as a Role-A diagnostic.

The later v21/v22 record does not strengthen that conclusion as previously
claimed. Phase 3 selected `core_preference = 1.0` and zero spillage with side
RMSE near 211 K. Phase 4 then used `core_preference = 0.6` and 525 W spillage,
so it did not inherit the selected macroscopic point. More importantly, every
Rosseland multiplier at a fixed Nusselt multiplier produces exactly identical
output. The patch changes `V11.felt_conductivity_temp`, while inherited v12/v14
paths call `V12.felt_conductivity_v12`. The 1500--1800 K v22 result therefore
diagnoses a configuration/wiring failure and cannot be called a formal
falsification or best-possible continuum result.

### Manuscript — model contradictions removed; bibliography still open

`manuscript_full_draft.md` runs from abstract through Section 6 conclusions
with figure and table captions. The orphaned model addendum and its capacity,
v34/v36, run-count, and encoding contradictions have been removed. Section 5.6
now states that the model coefficients are not validated. The standalone
`Role_A_Corroboration.md` has been rewritten to match the v20/v36 evidence.

The power-dependent discussion is now explicitly conditional: T3 is a local
outlet-region probe with an unidentified transfer relation, the R6 closure and
1D fit use the same campaign, and neither provides an independent optical power
calibration. This does not affect Nu, epsilon, Lambda, C_eff, K_loss, or the
master curves.

Outstanding manuscript work:

- finalise the references against the v5 bibliography;
- write the acknowledgments;
- incorporate the literature set for positioning only. Cui & Kaer (2018) and
  Schlereth & Hinrichsen (2014) provide precedent for 2D continuum monolith
  reduction and its limitations; Cornejo & Hayes (2020) is a channel-Nu
  contrast, not a validation bound; and
- decide whether provisional delivered-basis efficiency belongs in the main
  text or supplementary material pending direct flux/spillage and outlet
  enthalpy measurements.

---

## 2. The v41 result — what it actually shows

`1D_v41` froze the three source scales at the R6 K_heat closure values
(1.34 / 1.58 / 1.11), disabled the dynamic bypass (`phi_0 = 1.0`, `m_rec = 0`),
kept v40's optics (`beta_opt = 184.67`), and refitted the rest. Outcome:
objective 5780.6, `return_code = MaxTime`, mean steady heating errors of
**+909 K (T8), +914 K (T12), +773 K (T11), +567 K (T10), +346 K (T3)**,
K_loss 0.49–0.86 W/K against a measured 0.10–0.16, and eps_mean ≈ 1.000.

**This is not a verdict on the fixed-flux strategy. It is an energy-bookkeeping
error, and finding it is worth more than the run cost.**

The model's absorbed power is not `scale × G0 × A_frt`. Only 6.92% of the beam
enters the aperture (`flux_receiver_fraction`); 93.08% spills, and
`spill_capture` returns a fraction of that to the perimeter. So the true
delivered-to-nominal-aperture ratio is

    M = scale × (1 + 13.45 × spill_capture)

| Version | spill_capture | spill multiplier | M at 456 / 304 / 256 |
|---|---:|---:|---|
| v36 | 0.0280 | 1.376 | **1.88 / 2.04 / 0.87** |
| v37 | 0.0362 | 1.487 | 1.97 / 1.48 / 0.91 |
| v40 | 0.0921 | 2.238 | 2.98 / 2.30 / 1.41 |
| **v41** | **0.500** | **7.725** | **10.35 / 12.21 / 8.57** |
| R6 closure (K_heat) | — | — | 1.34 / 1.58 / 1.11 |
| R6 closure (K_cool) | — | — | 1.05 / 1.23 / 0.84 |

v41 therefore drove **5.5 times more total power** into the assembly than v36 at
456 kW/m2, and 10.4 times the nominal aperture power. A 900 K over-prediction is
the arithmetically expected consequence. The optimizer then did what it could to
shed it — `k_perim_ref` 7.5 → 48.4 W/m/K (aluminium-like), `G_core_perim`
10.6 → 27.1, K_loss to 3–5x the measured band — and still could not.

Two consequences, both correcting rev. 1 of this note and §6 of the Aug-1
readiness assessment:

1. **RETRACTION — `scale_456` was never comparable to R6's f.** The readiness
   note's headline agreement ("v36 fitted 1.369 vs R6 1.34, +2.2%") compares a
   partial source scale against a total delivered-power ratio. The comparable
   quantity is M. On that basis v36 sits at 1.88 / 2.04 / 0.87 against R6's
   1.34 / 1.58 / 1.11 — roughly +40% / +29% / −22%, not +2% / −6% / −43%. The
   qualitative story survives (non-monotonic, 256 low) but every number in that
   table needs restating before it goes near the manuscript. R6 is itself a
   T3-based provisional closure, not an independent optical measurement.
2. **RETRACTION — the "five dead parameters" finding is a screening artifact.**
   Rev. 1 reported that the Morris screen gives exactly zero elementary effects
   for `spill_capture`, `beta_perim`, `C_rear_eff` and `G_rear_axial`. v41 shows
   `spill_capture` moving from 0.028 to 0.5 changes absorbed power by a factor
   of 5.6 and the wall temperature by ~700 K. A parameter with that leverage
   cannot have mu* = 0. Identically zero mu* *and* sigma for four parameters is
   the signature of parameters that were never perturbed — the screen is buggy
   for those four. Do not prune on it; fix it or discard it.

There is also a real degeneracy here worth naming: `scale` and `spill_capture`
enter the absorbed power only through their product M. They are not separately
identifiable from temperature data alone. That degeneracy is the likely origin
of both the `scale_256` collapse and the erratic `spill_capture` values across
versions, and it is a cleaner explanation than any of the physical stories so
far attached to the 256 group.

---

## 3. Recommended sequence

### A. Freeze the story and document the audits — completed 2026-08-28

The journals now record v40, v41, the source-product degeneracy, the v22 wiring
failure, and the narrowed Role-A conclusions. Adopt v36 as the best balanced 1D
diagnostic, v37/v38 as tested alternatives, v39/v40 as unresolved or rejected
combinations, v41 as the power-bookkeeping diagnostic, and v20 as the
trustworthy 2D non-identifiability endpoint.

### B. Reparametrise the source, then redo the fixed-flux test (v42)

Replace `(scale, spill_capture)` with `(M, chi)`, where M is the total
delivered-to-aperture ratio and chi is the core/perimeter partition of it. Then:

- initially fix M per flux group to the R6 closure as a provisional sensitivity
  band, carrying the K_cool–K_heat spread and the unresolved T3 transfer
  relation explicitly, and let chi and the transport parameters fit;
- state explicitly what M excludes (front-face reradiation, flange cooling)
  when comparing to R6, since R6's f counts only gas enthalpy plus assembly
  loss;
- re-verify afterwards: Nu exponent, eps envelope, per-group inversion
  threshold, C_total, and whether the capacities leave their bounds once the
  source degeneracy is gone.

This is still the highest-payoff run. v41 does not argue against it — v41 is
what happens when it is done on the wrong variable, and it removes a degeneracy
that has been contaminating the source parameters since at least v31.

Note also that v41's optimisation timed out (`MaxTime`) at an objective four
orders of magnitude worse than v36. Whatever time budget it used is not enough
for a cold start; seed v42 from v36 and give it a longer budget.

### C. Bound release, and fix the sensitivity screen (v43)

Open `C_perim_eff` and `C_rear_eff` by ±30% and refit. If 302.5 J/K is
recovered from the interior, the inventory closure becomes a genuine validation
and §5.5 of the manuscript can be rewritten honestly. If not, report it as
imposed — still publishable, just a different claim. Repair the Morris screen
before using it to justify dropping anything.

### D. The inversion threshold — the one open piece of physics

The most scientifically consequential gap, and it is well-posed: at the correct
effectiveness the model needs ~0.08–0.10 more eps than the real receiver to move
the wall peak inside, and its threshold drifts with flux where the measurement's
does not. The deficiency is in *how axial source and axial exchange are placed
relative to each other*, not in the magnitude of gas–solid coupling, which v36
gets right.

Specific hypothesis, unchanged by v41 and now better supported: **v36's flux
dependence is imported by the mechanism itself.** Its active-flow fraction is a
function of `T_core`, and `T_core` scales with flux, so a temperature-driven
bypass cannot produce a flux-independent threshold. The measurement says the
threshold is set by something scale-free in flux — geometric or Re-driven.
Replacing phi(T_core) with an axially varying, temperature-independent
participation phi(z) (or phi(Re)) is the structurally correct move: flow
spreading past the front-face obstruction and contracting toward the 13 mm rear
adaptor is documented hardware, it redistributes exchange axially without a
surface-absorber optical model, and it matches the manuscript's own reading of
the exponent above unity as flow-dependent recruitment of solid participation
rather than a boundary-layer property. A Graetz/entry-region argument does *not*
work here — at Re 23–94 and Pr 0.69 the thermal entry length is a few
millimetres against a 137 mm receiver — which is itself a useful negative result
for the discussion.

### E. Manuscript work (parallel, independent of B–D)

1. Keep the corrected limitations language and do not restore the removed
   addendum claims.
2. Restate any future capacitance claim according to whatever step C returns.
3. If a delivered-power set is retained, report one provisional M basis with
   provenance, structural uncertainty, and the T3 caveat.
4. Lead any modelling subsection with the eps-envelope / inversion-threshold
   pairing, not a fit-quality table.
5. Preserve the revised "Advective Bottleneck" wording: v36 weakens its strong
   form, and the data do not uniquely identify bypass as the mechanism.
6. Read the literature set for positioning only. The receiver's apparent Nu is
   15–100x below duct theory; a monolith correlation is a contrast, not a bound.
7. Finalise references against the v5 bibliography; write acknowledgments.

---

## 4. Suggested immediate next step

A is complete. Next is B: implement the conserved `(M, chi)` source basis and
verify it with forward energy-ledger tests before optimization. The v41
post-mortem converts a failed run into a clean structural finding: source
magnitude and spill capture have not been separately identifiable, and the
power-convention anomalies from v31 onward are consistent with that degeneracy.
