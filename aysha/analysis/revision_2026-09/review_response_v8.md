# Response to review v8 of manuscript_revised_v2

All four must-fix items and all seven wording items are accepted and implemented. Nothing in this
round required a scientific judgement against the reviewer; every item was work owed. The
manuscript is now internally consistent with a single archived pipeline run, and we consider the
revision closed.

Line numbers refer to `manuscript_revised_v2.md` as submitted with this response. All values come
from one 4000-realization run of `receiver_reduction.py`, archived in `outputs/`.

## 1. Controller uncertainty — implemented as specified

Accepted in full, and the reviewer's diagnosis was exact: the code drew a fresh aggregate error for
each steady run and a separate one for the transient records, which is a run-to-run repeatability
the manufacturer does not specify, and it contradicted both the manuscript and our own previous
letter.

Implemented as directed. Four per-unit errors are drawn once per realization, each as a correlated
part of weight √ρ plus an independent part of weight √(1−ρ) so the per-unit total stays at the 2.5%
specification for any ρ. They are combined with each record's own four measured shares, and the
same four errors reach the corresponding steady and transient records. The raw cooling files do
contain the four controller signals, as the reviewer noted, so cooling shares are now read from
them rather than borrowed from the steady mean (§3.4, line 123).

This changed a conclusion, and in the reviewer's favour. With errors correctly persistent across
runs, the two correlation cases become nearly indistinguishable:

| Quantity | ρ = 1 | ρ = 0 |
|---|---|---|
| Nu exponent, pooled | ±0.0029 | ±0.0031 |
| Nu exponent, grouped | ±0.0031 | ±0.0031 |
| N_prof exponent | ±0.0029 | ±0.0030 |
| C_eff matched-ε | ±37.0 J K⁻¹ | ±36.6 J K⁻¹ |

The reason is that the four shares vary little between runs, so a persistent set of unit errors
produces a per-run flow error that is almost purely systematic whatever ρ is. It therefore displaces
a correlation's prefactor and barely tilts its exponent. The consequence for the paper is that the
two-limit reporting we introduced last round is no longer needed: single instrumental values are
quoted, the agreement between cases is stated, and both remain archived (§3.4, line 121).
Our previous letter's claim that the independent component "does reach an exponent" was an artefact
of the defective implementation and is withdrawn.

Both cases were rerun at 4000 realizations and `uncertainty.csv`, `results.json` and Table 2 were
regenerated. The instrumental term on the Nusselt exponent moved from ±0.005 to ±0.003 and the
matched-ε Monte Carlo mean from 291 to 280 J K⁻¹ at ρ = 0; every affected number in the text was
updated.

## 2. Matched estimator — synchronized, and the provenance is now fixed

All three sub-points accepted.

*Methods.* §3.2 no longer says the pooled ε(q) correlation is in use. It states that each cooling
transient is referred to the heating run it decays from, and that the pooled correlation is reported
for comparison only (line 109).

*Provenance.* The reviewer is right that a nearest-flow match is unstable: E81 and E76 differ by
0.001 sL min⁻¹, so a perturbed flow could switch a cooling run between them inside the Monte Carlo,
and that instability was silently inflating the matched-ε spread. Provenance is now a fixed map
`COOL_PROVENANCE = {"C69": "E69", "C80": "E80", "C81": "E81"}` taken from the experiment log, used
identically in the point estimate and in every realization.

*Similarity.* §3.3 and §4.5 now rescale on the primary identification and on each run's own measured
effectiveness rather than on the pooled constants (line 115). The values become t*½ = 0.205
with CV 16.8% for the wall and 0.634 with CV 7.8% for the outlet — matching the reviewer's
independent check to three digits and two, which we take as confirmation the rescaling is now the
intended one (line 201).

The stale claim that the wall half-rise "landed at 0.20 rather than near ln 2" is removed. The
half-rise under this rescaling has no expected value; ln 2 is what a single-node system would give,
and the departure from it is itself the measure of how far the receiver is from one state. That is
now what the text says.

Figure 5b was also inconsistent with the primary estimator — it plotted the cooling points at the
pooled abscissa. It now plots them at the matched abscissa as filled symbols with the pooled ones
open behind, and the fitted line is the primary identification. The matched abscissa is archived
per run in `eigenvalues.csv` as `x_matched` so the figure is not recomputing it.

## 3. Exact Re⁻¹ — removed from both remaining places

Accepted. §4.2 (line 169) now gives the computed requirement, −0.9996 for a fixed Nusselt
number and −1.017 for a conductance independent of temperature as well as flow, and adopts −1.00 as
the working value. The conclusion (line 340) writes Re⁻¹·⁰⁰, gives −0.9996 rather than
"exact", and now carries the reviewer's own high-level formulation naming all three tested families.

## 4. Reproducibility claim — archive now matches it

Accepted. `results.json → fixed_profile_test` archives per-run residuals and run ids for all three
families, not only the shared-conductance one, and `wall_extrapolation` records which runs make a
full linear rear continuation inadmissible. A `README_reproduction.md` accompanies the package,
stating the two commands that regenerate everything, the runtime, the seeds that make both scripts
deterministic, what each archived file supports, and the nominal/apparent conventions a reader needs
before using any number.

## Wording items

All seven accepted: the unquantified-bias sentence now points at §4.2 where the bias is quantified
(line 91); the web Biot statement is now about through-web conduction not being the limiting
resistance (line 171); the "everywhere, one to two orders" phrasing is replaced by the tested
scope and the 36 to 68 K range (line 292); the §4.1 crossing is described as the side-wall
signature with the caveat that it is a two-probe difference and not the interior peak; the
conclusion uses the same conventions throughout; the sentence-case error in the data-availability
paragraph is fixed and Figure 1 is now identified as drawn rather than computed (line 352); and
a semicolon splice introduced by the §3.2 edit was corrected.

## Verification

Sixty pipeline scalars were checked against `results.json` and `uncertainty.csv`, all present and
correct in the text. Reference cross-referencing is consistent in both directions at 46 references.
The reduction and figure scripts run clean from the raw files with logs archived.

## On the reviewer's framing

We agree with the assessment that this is a measurement, identifiability and model-class diagnostic
paper, and that no further conceptual expansion or literature work is needed. We also agree with the
recommendation to state the falsification within its tested scope and to leave the physical cause
unidentified, which is what §4.2, §4.6 and §5.2 now do.
