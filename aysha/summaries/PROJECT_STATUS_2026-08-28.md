# Project status and forward plan — 28 August 2026

Restart briefing after a ~4 week pause. Reconstructed from `journal.1D.md`,
`journal.2D.md`, `1D_v36_v37_v39_manuscript_readiness.md` (1 Aug),
`manuscript_2D_readiness_strategy.md`, `1D_manuscript_gap_strategy.md`,
`analysis/manuscript/manuscript_full_draft.md`, and the saved CSVs under
`summaries/1D_v3{5,6,7,8,9}/`, `1D_v4{0,1}/` and `2D_v2{0,1,2}/`.

**Last substantive activity:** 1 Aug 2026 (1D_v39/v40 runs, v39 theory manual,
readiness assessment). Only change since: four literature PDFs added on 5 Aug
(Zaversky 2018, Cui & Kaer 2018, Schlereth 2014, Wang 2021) — not yet read into
any document.

---

## 1. Where the three tracks stand

### 1D — converged on a preferred model, with open blockers

| Version | Mechanism | Objective | Verdict |
|---|---|---:|---|
| v35 | baseline, scalar bypass | 1.866 | superseded; "Advective Bottleneck" argument comes from here |
| **v36** | **dynamic bypass, phi(T_core)** | **0.2173** | **best model in the family — manuscript candidate** |
| v37 | optical redistribution | 0.2388 | falsified as a formulation (surface receiver); keep as sensitivity test |
| v38 | two-stream gas | 5.701 | falsified |
| v39 | v36+v37 combined | 0.629 reported | **not citable** — saved residuals contradict the objective |
| v40 | v39 variant, m_rec active | 1.476 | undocumented; worse than v36 on nearly every gate |
| v41 | — | — | only a Morris sensitivity CSV; no run, no notes |

**v36 is the strongest result the campaign has produced.** Mean steady heating
error inside ±40 K on every solid channel (T8 +13.8, T9 −21.4, T12 −40.1,
T11 +16.3, T10 +28.7 K), T3 essentially exact (−0.1 K), heating and cooling
satisfied simultaneously with no switch. Nu exponent 1.4394 against the
measured 1.44. The measured effectiveness envelope is reproduced run-by-run
(0.583–0.777 vs 0.573–0.781). C_total 302.5 J/K against 301 ± 23.

**What still blocks v36:**
- `C_perim_eff = 150.0` and `C_rear_eff = 80.0` sit *exactly* on their lower
  bounds; `flange_scale` and `front_dep` also at edges. The capacity gate
  passes because the bounds were placed to make it pass.
- Nu prefactor 1.74e-4 vs 3.1e-4 (1.8x low).
- Lambda_107 slope 5.7e-5 vs 8.3e-4, r2 = 0.083 — no Re-trend at all.
- `scale_256 = 0.629` vs the R6 closure estimate 1.11 — 43% low.
- The inversion threshold sits at eps ~ 0.72–0.78 and drifts with flux,
  against a measured flux-independent 0.66; at 256 kW/m2 it never inverts.

**v39/v40 are unresolved loose ends.** v39's saved CSVs show mean steady errors
of −180 to −243 K on every channel, roughly 200 K worse than v35, which
contradicts the reported objective of 0.629 — either the diagnostics were
written from a non-optimal vector or the loss normalisation differs. v40 flips
the sign of the error (T8 +236, T12 +156, T3 −100 K) with no journal entry and
no theory manual explaining what changed. Two of the three most recent runs are
undocumented.

**v41 contains an unread and useful result.** The Morris screen shows
`spill_capture`, `beta_perim`, `C_rear_eff` and `G_rear_axial` have *exactly
zero* elementary effects, and `G_rear_tube` is ~1e-2 — five of twenty-five
parameters are non-identifiable in the present formulation. `C_perim_eff`
dominates (mu* = 1.7e4) with an enormous sigma, i.e. strongly interacting.
This is direct evidence for pruning the parameter vector before any further
refit, and it has not been written into any note.

### 2D — closed out, deliberately

Phases 2–4 (v21/v22) finished 31 Jul. The grid searches drove core preference
to 1.0 and spillage to 0.0, reproducing the temperature inversion
*qualitatively* while leaving side RMSE ~210 K; the Nu/Rosseland sweep left
RMSE at 1500–1800 K. The v20 gate-free stress test independently confirmed that
the failure is structural, not a gate or a solver problem (max enthalpy
residual 6e-11). The standing decision in
`manuscript_2D_readiness_strategy.md` §12 is: **stop coefficient extraction
from the 2D inverse model**, report it as a structural-limit demonstration
(Role A) and state that receiver UA(Re) is not uniquely identifiable without an
independent bulk outlet-enthalpy measurement or a verified T3 transfer function.
Nothing here needs re-opening. It needs writing up.

### Manuscript — complete draft, unmerged addendum, three contradictions

`manuscript_full_draft.md` runs from abstract through Section 6 conclusions with
figure and table captions. Outstanding:

- **References not finalised** (bracketed placeholders against the v5
  bibliography); **acknowledgments empty**.
- **An addendum block sits after the table captions, unmerged** into the body:
  "Energy Conservation and Structural Support", "Reconciled Delivered-Power
  Factors", "Explicit Disclaimer on Validation", and a §5.5 "Quantitative Model
  Corroboration".
- **Three internal contradictions in that block, in order of reviewer risk:**
  1. §5.5 claims the model "organically identifies 302.5 J/K" under "a
     completely uninformative prior", calling it independent verification of
     C_eff. §3 of the same block states C_eff is ingested as a prior, and the
     readiness assessment shows the capacities sit on their bounds. As written
     this is a claim a reviewer can falsify from the parameter file.
  2. §5.5 quotes v34 numbers (prefactor 3.6e-4, exponent 1.11) as the model
     corroboration. The current model, v36, gives 1.74e-4 and 1.4394 — the
     exponent match is far *better* than what is written, the prefactor worse.
  3. Three different delivered-power factor sets are in circulation:
     §2 of the block (1.34 / 1.37 / 0.79), the R6 closure
     (1.34 / 1.58 / 1.11 on K_heat; 1.05 / 1.23 / 0.84 on K_cool), and v36's
     fitted values (1.369 / 1.484 / 0.629). Only one set can appear.
- Encoding damage in the appended text ("3.6A-10^-4", "301 A 23 J/K").
- The 21 literature PDFs have not been cited anywhere. The four added on 5 Aug
  are directly relevant to the 2D write-up: Cui & Kaer (2018) and Schlereth &
  Hinrichsen (2014) are precedent for 2D continuum monolith reduction and its
  limits; Cornejo & Hayes (2020) is the monolith-Nu *contrast* for §5.1.

---

## 2. Recommended sequence

The Aug-1 readiness note listed eight actions; none were executed except the
v39 theory manual. The list below supersedes it, reordered by payoff and with
the v40/v41 findings folded in.

### A. Freeze the story before running anything else (half a day)

Adopt: **v36 = the manuscript 1D model**; v37 and v38 = the two falsified
alternative hypotheses; v39/v40 = the non-additivity test; 2D = the structural
-limit demonstration. Write the missing journal entries for v40 and v41 now,
including the Morris result, and mark v41's dead parameters. Everything after
this is either a refit of v36 or manuscript prose.

### B. Fixed-flux refit — the single highest-payoff run (v42)

Replace the three fitted `scale_*` with the R6 delivered-power factors and
refit the remainder. Rationale (readiness §7): eps, NTU, Nu_app, Lambda_107,
C_eff and K_loss contain no irradiance term, so fixing the flux basis moves the
model and not the targets — three degrees of freedom are removed and the
comparison becomes *stricter*, not looser. The 256 group is simultaneously the
most under-predicted and the only group that never inverts; both defects point
the same way, and raising that group's power by 26–76% is the correct direction
for both. Verify afterwards: Nu exponent, eps envelope, per-group inversion
threshold, C_total, and whether the capacities leave their bounds once three
source parameters stop competing with them.

### C. Bound-release and parameter pruning (v43)

Open `C_perim_eff` and `C_rear_eff` by ±30% and drop the five parameters the
Morris screen shows are dead. If 302.5 J/K is recovered from the interior, the
inventory closure becomes a genuine validation and §5.5 can be rewritten
honestly. If it is not, the capacity claim must be reported as imposed — which
is still a publishable statement, just a different one.

### D. The inversion threshold — the one open piece of physics

This is the most scientifically consequential gap and it is well-posed: at the
correct effectiveness, the model needs ~0.08–0.10 more eps than the real
receiver to move the wall peak inside, and its threshold drifts with flux where
the measurement's does not. The deficiency is in *how axial source and axial
exchange are placed relative to each other*, not in the magnitude of gas–solid
coupling, which v36 gets right.

A specific hypothesis worth testing before more optics: **the flux dependence
of v36's threshold is imported by the mechanism itself.** v36's active-flow
fraction is a function of `T_core`, and `T_core` scales with flux, so a
temperature-driven bypass cannot produce a flux-independent threshold. The
measurement says the threshold is set by something scale-free in flux — i.e.
geometric or Re-driven. Replacing phi(T_core) with an axially varying,
temperature-independent participation phi(z) (or phi(Re)) is the structurally
correct move: flow spreading past the front-face obstruction and contracting
toward the 13 mm rear adaptor is a documented feature of this hardware, it
redistributes exchange axially without invoking a surface-absorber optical
model, and it is consistent with the manuscript's own reading of the exponent
above unity as flow-dependent recruitment of solid participation rather than a
boundary-layer property. Note that a Graetz/entry-region argument does *not*
work here — with Re = 23–94 and Pr = 0.69 the thermal entry length is a few
millimetres against a 137 mm receiver — which is itself a useful negative
result for the discussion.

### E. Manuscript work (parallel, independent of B–D)

1. Merge the addendum into §5.5/§5.6 and fix the three contradictions above.
2. Restate the capacitance claim according to whatever C returns.
3. Choose one delivered-power factor set and state its provenance and band.
4. Lead the modelling subsection with the eps-envelope / inversion-threshold
   pairing, not a fit-quality table — that pairing is a quantitative,
   well-posed statement, and it is the most interesting thing the models found.
5. Rewrite the "Advective Bottleneck" conclusion in `journal.1D.md` before it
   migrates into the manuscript: v36 falsifies its strong form, since a
   *dynamic* flow structure does satisfy heating and cooling simultaneously.
   The defensible claim is narrower — no *scalar* structure works, and no
   structure tried so far relocates the inversion threshold.
6. Read the literature set for positioning only. The receiver's apparent Nu is
   15–100x below duct theory; a monolith correlation is a contrast, not a bound.
7. Finalise references against the v5 bibliography; write acknowledgments.

---

## 3. Suggested immediate next step

A then B: write the two missing journal entries and the pruning note (a
morning's work, and it is the piece that decays fastest), then launch the
fixed-flux refit, since it is the one run that could close the inversion-
threshold gap for free and it removes three fitted degrees of freedom from the
manuscript's most exposed claim.
