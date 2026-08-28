# Project status and forward plan — 28 August 2026

**Revision 2** (rev. 1 written earlier the same day, before the 1D_v41 full run
appeared). Restart briefing after a ~4 week pause. Reconstructed from
`journal.1D.md`, `journal.2D.md`, `1D_v36_v37_v39_manuscript_readiness.md`
(1 Aug), `manuscript_2D_readiness_strategy.md`, `1D_manuscript_gap_strategy.md`,
`analysis/manuscript/manuscript_full_draft.md`, and the saved CSVs under
`summaries/1D_v3{5,6,7,8,9}/`, `1D_v4{0,1}/` and `2D_v2{0,1,2}/`.

**What changed since rev. 1:** a complete `1D_v41` run has appeared — the
fixed-flux refit. It is a failed optimisation but a very informative
diagnostic, and it overturns two conclusions in rev. 1. Both retractions are
marked below. `journal.2D.md` shows a timestamp change with no content change.
No other file has moved since 5 Aug.

---

## 1. Where the three tracks stand

### 1D — a preferred model, and a failed test that located the real problem

| Version | Mechanism | Objective | Verdict |
|---|---|---:|---|
| v35 | baseline, scalar bypass | 1.866 | superseded; "Advective Bottleneck" argument comes from here |
| **v36** | **dynamic bypass, phi(T_core)** | **0.2173** | **best model in the family — manuscript candidate** |
| v37 | optical redistribution | 0.2388 | falsified as a formulation (surface receiver); keep as sensitivity test |
| v38 | two-stream gas | 5.701 | falsified |
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

### 2D — closed out, deliberately

Phases 2–4 (v21/v22) finished 31 Jul. The grid searches drove core preference
to 1.0 and spillage to 0.0, reproducing the temperature inversion
*qualitatively* while leaving side RMSE ~210 K; the Nu/Rosseland sweep left
RMSE at 1500–1800 K. The v20 gate-free stress test independently confirmed the
failure is structural, not a gate or solver problem (max enthalpy residual
6e-11). Standing decision (`manuscript_2D_readiness_strategy.md` §12): **stop
coefficient extraction from the 2D inverse model**, report it as a
structural-limit (Role A) demonstration, and state that receiver UA(Re) is not
uniquely identifiable without an independent bulk outlet-enthalpy measurement
or a verified T3 transfer function. Nothing to re-open — it needs writing up.

### Manuscript — complete draft, unmerged addendum, three contradictions

`manuscript_full_draft.md` runs from abstract through Section 6 conclusions with
figure and table captions. Outstanding:

- **References not finalised** (bracketed placeholders against the v5
  bibliography); **acknowledgments empty**.
- **An addendum block sits after the table captions, unmerged** into the body.
- **Three internal contradictions in that block, in order of reviewer risk:**
  1. §5.5 claims the model "organically identifies 302.5 J/K" under "a
     completely uninformative prior", calling it independent verification of
     C_eff. §3 of the same block states C_eff is ingested as a prior, and the
     capacities sit on their bounds. A reviewer can falsify this from the
     parameter file.
  2. §5.5 quotes v34 numbers (prefactor 3.6e-4, exponent 1.11) as the model
     corroboration. v36 gives 1.74e-4 and 1.4394 — the exponent match is far
     *better* than what is written, the prefactor worse.
  3. Three different delivered-power factor sets are in circulation: §2 of the
     block (1.34 / 1.37 / 0.79), the R6 closure (1.34 / 1.58 / 1.11 on K_heat;
     1.05 / 1.23 / 0.84 on K_cool), and the models' fitted `scale_*`. §2 below
     shows the last of these was never the right comparison quantity.
- Encoding damage in the appended text ("3.6A-10^-4", "301 A 23 J/K").
- The 21 literature PDFs are cited nowhere. The four added 5 Aug matter for the
  2D write-up: Cui & Kaer (2018) and Schlereth & Hinrichsen (2014) are
  precedent for 2D continuum monolith reduction and its limits; Cornejo & Hayes
  (2020) is the monolith-Nu *contrast* for §5.1.

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

Two consequences, both correcting rev. 1 of this note and §7 of the Aug-1
readiness assessment:

1. **RETRACTION — `scale_456` was never comparable to R6's f.** The readiness
   note's headline agreement ("v36 fitted 1.369 vs R6 1.34, +2.2%") compares a
   partial source scale against a total delivered-power ratio. The comparable
   quantity is M. On that basis v36 sits at 1.88 / 2.04 / 0.87 against R6's
   1.34 / 1.58 / 1.11 — roughly +40% / +29% / −22%, not +2% / −6% / −43%. The
   qualitative story survives (non-monotonic, 256 low) but every number in that
   table needs restating before it goes near the manuscript, and the third
   delivered-power set in circulation should be dropped rather than reconciled.
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

### A. Freeze the story, then write the three missing entries (half a day)

Adopt: **v36 = the manuscript 1D model**; v37 and v38 = the two falsified
alternative hypotheses; v39/v40 = the non-additivity test; v41 = the
power-bookkeeping diagnostic; 2D = the structural-limit demonstration. Write
journal entries for v40, v41 and the M-reparametrisation finding now. This is
the piece that decays fastest and three of the last four runs are undocumented.

### B. Reparametrise the source, then redo the fixed-flux test (v42)

Replace `(scale, spill_capture)` with `(M, chi)`, where M is the total
delivered-to-aperture ratio and chi is the core/perimeter partition of it. Then:

- fix M per flux group to the R6 closure, carrying the K_cool–K_heat spread as
  a declared systematic band, and let chi and the transport parameters fit;
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

1. Merge the addendum into §5.5/§5.6 and fix the three contradictions.
2. Restate the capacitance claim according to whatever C returns in C.
3. Report **one** delivered-power set, on the M basis, with provenance and band.
4. Lead the modelling subsection with the eps-envelope / inversion-threshold
   pairing, not a fit-quality table.
5. Rewrite the "Advective Bottleneck" conclusion in `journal.1D.md` before it
   migrates into the manuscript: v36 falsifies its strong form, since a
   *dynamic* flow structure does satisfy heating and cooling simultaneously.
   The defensible claim is narrower — no *scalar* structure works, and nothing
   tried so far relocates the inversion threshold.
6. Read the literature set for positioning only. The receiver's apparent Nu is
   15–100x below duct theory; a monolith correlation is a contrast, not a bound.
7. Finalise references against the v5 bibliography; write acknowledgments.

---

## 4. Suggested immediate next step

A, then B. The v41 post-mortem and the M reparametrisation are the same piece of
work, and writing it down is what converts a failed run into the campaign's
cleanest structural finding: the source magnitude and the spillage capture have
never been separately identifiable, and every power-convention anomaly from v31
onward is consistent with that single degeneracy.
