# Options for closing the first manuscript

Date: 28 August 2026. Constraint: no further experiments are possible, and new
computation should be minimised. Companion to `PROJECT_STATUS_2026-08-28.md`.

---

## 1. The decisive fact

**Every headline claim in the current draft is model-free.** R1, R3, R4, R5, R6
and R7/M4 are algebraic reductions of measured temperatures, metered flow and
gas properties. None of them contains a fitted coefficient, and none of them
depends on the 1D or 2D model being right:

| Result | What it is | Depends on a model? |
|---|---|---|
| R1 | envelope: eps 0.57–0.78, NTU 0.85–1.51, Re 23–94 | no |
| R3 | Nu_app = 3.1e-4 Re^1.44, r2 = 0.97, 15–100x below duct theory | no |
| R4 | inversion criterion eps* = 0.671/0.671/0.626, flux-independent | no |
| R5 | Lambda_107 = 0.038 + 8.3e-4 Re (lower bound); LTE fails everywhere | no |
| R6 | delivered-power closure; gas output exceeds nominal by 23 ± 4% | no |
| R7/M4 | C_eff = 301 ± 23 J/K, K_loss = 0.10–0.16 W/K, 15 confirming eigenvalues | no |
| M5 | 15 transients collapse on C_eff/(mdot cp + K_loss), CV 9–13% | no |

The manuscript is therefore **already complete as a scientific object**. What is
missing is editorial (references, acknowledgments) and the removal of an
addendum that makes model-dependent claims the data do not need.

This reframes the whole question. It is not "how much modelling work is needed
to publish"; it is "how much modelling should be allowed into a paper that does
not require any."

## 2. What must not be claimed, at any option

These are contaminated and no reasonable amount of new work fixes them for a
first paper:

- **Any fitted transport coefficient presented as validated** (Role B). Neither
  1D nor 2D passed that gate, and the 2D identifiability analysis says UA(Re) is
  not reachable from this instrumentation at all.
- **"The model independently confirms C_eff = 301 J/K."** The capacities sit on
  their lower bounds and C_eff is also ingested as a prior. This is the single
  most falsifiable sentence in the draft.
- **Delivered-power factors read off `scale_*`.** `scale` and `spill_capture`
  enter absorbed power only through M = scale x (1 + 13.45 x spill_capture);
  they are not separately identifiable. The three competing factor sets in the
  draft are an artifact of this.
- **Any v39 or v40 number.** Internally inconsistent, undocumented.
- **The model's Nu prefactor.** The reduction convention has not been verified
  against the manuscript definition.

## 3. What is robust from the modelling, at zero additional computation

Three classes of model result survive everything above, because none of them
requires a fitted coefficient to be correct:

1. **Ratio-level agreement.** v36 reproduces the measured effectiveness
   envelope run-by-run (0.583–0.777 vs 0.573–0.781, NTU 0.88–1.50 vs 0.85–1.52).
   eps is a ratio of temperature differences and is therefore insensitive to the
   absolute source magnitude — which is precisely why this result survives the
   source degeneracy that invalidates the power claims.
2. **A sharp, quantitative discrepancy.** At the correct effectiveness, the
   model needs roughly 0.08–0.10 *more* eps than the real receiver to move the
   wall peak inside, and its threshold drifts with flux where the measurement's
   does not (v36: 0.720–0.743 at 456, 0.747–0.777 at 304, no crossing at 256).
   This is a statement about the relative axial placement of source and
   exchange, not about coupling magnitude — and it is a better result than
   agreement would have been.
3. **Falsifications and identifiability limits.** Each is a negative result that
   is *strengthened*, not weakened, by the fits failing:
   - a two-stream core/perimeter gas split is rejected by the optimizer (v38);
   - no *scalar* flow structure satisfies heating and cooling simultaneously
     (the advective bottleneck), while a temperature-dependent one does — so the
     correct claim is narrower than the journal currently states;
   - a 2D axisymmetric continuum fed extreme, physically-derived maldistribution
     boundary conditions reproduces the inversion qualitatively and still leaves
     ~210 K side RMSE, with all conservation checks passing to 1e-11;
   - source magnitude and spillage capture are not separately identifiable from
     temperature data when 93% of the beam spills — a general, quantitative
     backing for the draft's existing recommendation that future campaigns
     integrate the flux map over the receiver face.

---

## 4. The options

### Option A — Experiment-only paper

Cut all modelling. Delete the addendum block. Scope the paper as a
measurement-and-reduction study: the title already says so.

- **New computation:** none.
- **Work:** finalise references against the v5 bibliography, acknowledgments,
  scope the "no correction factors applied" statement to the experimental
  reduction, tidy §5.6.
- **Strength:** every claim is model-free and defensible. Nothing in the paper
  can be attacked through a model.
- **Weakness:** a year of modelling stays unpublished, and a referee may ask why
  no model was attempted — answerable in one sentence, but it is a question.
- **Verdict:** the safe floor. Ship-ready in the shortest time of any option.

### Option B — Experiment + a compact modelling section framed as falsification and identifiability *(recommended)*

Keep Option A's backbone and add one Results/Discussion section built entirely
from §3 above. It reports what the continuum models *cannot* do and what the
instrumentation *cannot* identify, and reports no fitted coefficient as physical.

- **New computation:** one plotting script over existing CSVs — no new
  simulations. Two panels: (i) modelled vs measured eps parity, run by run;
  (ii) I_vol = T12 − T8 against eps for measurement and model, showing the
  threshold displacement and its flux dependence.
- **Work:** ~1500–2500 words plus a half-page formulation summary (the v36
  theory manual already has the equations), one figure, one table.
- **Strength:** converts the modelling campaign into publishable content without
  a single new fit; inoculates against the "did you model this?" referee; and
  the identifiability results are genuinely useful to anyone building receivers
  of this class. The eps-envelope/threshold pairing is the most interesting
  single thing either model produced.
- **Weakness:** a negative result about a model is only as credible as the model
  is. The section must include enough formulation detail to make the falsifi-
  cation land, which is real writing and gives a referee a surface to attack.
  Mitigate by keeping the model description tight and the claims conditional:
  "within this class of formulations", never "no model can".
- **Verdict:** best value per unit of work. Recommended.

### Option C — Full Role B coefficient-extraction paper

Complete the (M, chi) reparametrisation, the bound release, the invariant
recomputation, and the v39/v40 reconciliation; then report validated
coefficients.

- **New computation:** several full calibration campaigns, uncertain outcome.
- **Verdict:** not a first paper. The 2D identifiability analysis already
  concludes UA(Re) is not uniquely identifiable from this instrumentation, and
  no further experiments are possible — so the gate may be unreachable in
  principle with this dataset, not merely unreached. Pursuing it now risks
  spending months to arrive back at Option B.

### Option D — Option B now, identifiability paper later

Publish B, then develop the modelling material into a second, methods-oriented
paper: what can and cannot be identified from wall-thermocouple data in a
small-aperture volumetric receiver, using the 1D and 2D campaigns as the
worked case. The degeneracy, the T3 observer non-identifiability, the ridge
structure in (Nu50, n), and the failed hypothesis branches are the content.

- **Verdict:** the natural home for everything Option B cannot fit. Requires no
  new experiments either. Decide after B is submitted, not before.

---

## 5. Recommendation

**Option B, then D.** Concretely, and in order:

1. Delete the addendum block wholesale. Every claim in it is either already
   made better in the body or is one of the contaminated claims in §2.
2. Write the modelling section from §3 above. Lead with the eps-envelope
   agreement, then the threshold displacement as the diagnostic, then the
   falsifications, then the identifiability limits. Close with the flux-map
   recommendation, now quantitatively backed.
3. Add one paragraph to §5.6 Limitations stating plainly that the transport
   coefficients are effective parameters, that C_eff enters the models as a
   prior, and that the source magnitude and spillage capture are degenerate.
   Said once, clearly, in Limitations, this is a strength; said as a disclaimer
   attached to a corroboration claim, it reads as a retraction.
4. Rewrite the advective-bottleneck conclusion in `journal.1D.md` to its
   defensible form before it migrates into the text.
5. Finalise references and acknowledgments.
6. Use the literature set for positioning only — Cornejo & Hayes (2020) as the
   monolith-Nu contrast; Cui & Kaer (2018) and Schlereth & Hinrichsen (2014) as
   precedent for 2D continuum monolith reduction and its limits; Kribus (2014),
   Avila-Marin (2011, 2019) and Fend (2004, 2013) for the volumetric framing.

Nothing in this list requires a new simulation, and the only new computation is
one plotting script.

## 6. What the paper's scientific contribution then is

Stated plainly, so the abstract can be written against it:

1. A flux-independent operating-point criterion for the volumetric inversion,
   eps* ~ 2/3, resolved to ±0.005 — a design condition on the heat-exchange
   operating point rather than on the optics.
2. An apparent assembly-scale exchange law, Nu = 3.1e-4 Re^1.44, lying 15–100x
   below duct theory with an exponent above unity — evidence that the limiting
   mechanism is flow-dependent recruitment of participating solid, not the
   channel boundary layer.
3. Quantified local thermal nonequilibrium throughout the envelope, as a
   measured lower bound.
4. A one-parameter transient identification, C_eff and K_loss, cross-validated
   by fifteen independent eigenvalues, showing the housing carries about six
   sevenths of the participating thermal mass.
5. A model-free delivered-power closure that localises a delivery error by lamp
   configuration and reports every efficiency as a bound.
6. The identifiability boundary: which of these quantities a continuum model of
   this class can reproduce (the effectiveness envelope), which it cannot (the
   location and flux independence of the inversion threshold), and which the
   instrumentation cannot determine at all (receiver UA(Re); source magnitude
   versus spillage capture).

Item 6 is what the modelling year bought. It is worth stating, and it costs one
figure.
