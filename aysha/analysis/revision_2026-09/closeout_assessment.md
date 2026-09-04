# Close-out assessment: what to cut for submission

The revision cycle is closed on correctness. Nothing in the manuscript is now unfixable, so nothing
needs removing on those grounds. What does need a decision is length, which every review round has
excluded from scope and which is now the only barrier to submission.

## Where the manuscript stands

The body is 18,338 words with a 725-word abstract, 46 references, 6 figures and 2 tables. Typical
limits are 8,000 to 10,000 words for *Solar Energy*, *Applied Thermal Engineering* or *Renewable
Energy*, and 250 to 300 for the abstract. The manuscript is roughly twice a submittable length, and
the abstract between two and three times.

Section shares of the body:

| Section | Words | Share |
|---|---:|---:|
| 5.2 Against uniform-temperature or uniform-flux coefficients | 3137 | 17.1% |
| 5.1 Against validating on a gas outlet temperature | 1461 | 8.0% |
| 1. Introduction | 1514 | 8.3% |
| 4.6 Delivered aperture irradiance | 1377 | 7.5% |
| 4.2 Assembly-scale exchange and the structural constraint | 1284 | 7.0% |
| 4.4 Transient identification and the energy budget | 1165 | 6.4% |
| 6. Conclusions | 1009 | 5.5% |
| 3.4 Uncertainty | 836 | 4.6% |
| 4.1 Steady field and the inversion criterion | 831 | 4.5% |
| Abstract | 725 | 4.0% |
| 4.7 Pressure drop as installed | 660 | 3.6% |
| 3.1 Dimensionless groups | 649 | 3.5% |
| 4.3 Local nonequilibrium | 590 | 3.2% |
| 3.2 Transient identification | 449 | 2.4% |
| 5.3 What would resolve it | 426 | 2.3% |
| 6 others (2.1–2.3, 3.3, 4.5, 5.4, back matter) | 2225 | 12.1% |

## Recommendation: one paper, roughly 9,000 words

The paper has one argument — that the standard receiver validation dataset does not constrain the
model class it is used to calibrate, and that this receiver's flow response falsifies a specific,
widely used closure. Everything that serves that argument stays. Everything that is a second paper
in embryo goes to supplementary material, where the pipeline already archives it.

**Move to supplementary, keeping a two-to-three sentence summary and a pointer in the main text.**

*§5.2, cut from 3137 to about 1200 words.* This is the single largest saving available and the
section least damaged by it. Its core is short: the measured flow response requires a
flow-dependent participating contact; a flow-dependent split requires inertial resistance, a
viscosity difference or buoyancy; the two candidate mechanisms are viscosity-driven redistribution
and peripheral bypass; the remedy is a multi-stream formulation. The literature scaffolding around
that core — the three-agent generalisation across packed beds, molten-salt panels and this
honeycomb, the parallel-channel boiling analogies cited for structure only, the Cornejo-group
flow-invariance derivation, the criterion lineage from 1996 through 2006 to 2022 — is a review-grade
discussion that would carry a paper of its own. Keep the flow-invariance result, which is the one
that does load-bearing work, and move the rest.

*§4.6, cut from 1377 to about 700 words.* The result the paper needs is that the radiometric and
enthalpy determinations disagree by a flow-dependent amount that cannot presently be assigned. The
gauge spatial-averaging calculation, the per-configuration factor table and the reconciling-δT3
analysis are supporting detail; `flux_geometry.json` and `results.json → delivered_power` already
hold them.

*§3.4, cut from 836 to about 350 words.* State the three kinds of uncertainty, the coverage-factor
convention and the controller model in a paragraph. The four-part justification of what was wrong in
earlier versions is correspondence, not method, and belongs in the response letters — which are
part of the record anyway.

*§4.7, cut from 660 to about 250 words.* The finding is one sentence: the installed transducer
cannot resolve the monolith pressure drop, for reasons of both range and tap placement, so the
absorber cannot be placed on either side of the stability criterion. The dual-cell-density and
misaligned-monolith comparisons and the zero-offset diagnosis go to supplementary.

*Conclusions, cut from 1009 to about 500 words.* It currently restates §4 and §5 in sequence. Three
paragraphs suffice: what was measured, what it rules out, what would resolve it.

**Nothing to delete outright.** I considered four candidates and reject all four. The pressure-drop
null result motivates the paper's central recommendation and would leave §5.3 unsupported. The
energy-closure discrepancy is the strongest independent evidence for a flow-dependent term. The
disequilibrium indices are the only direct check available on the two-equation closures the field
uses. The similarity collapse is what a start-up model needs and is one paragraph.

**The abstract needs rewriting to 300 words, not trimming.** At 725 words it currently carries the
falsification protocol, three families of fitted profiles, four capacitances and the closure spans.
A 300-word abstract holds the geometry and conditions, the two contradicting results with their
exponents, the identifiability finding, and the recommendation. Everything else is in the paper.

## Estimated result

| Element | Now | Target |
|---|---:|---:|
| Abstract | 725 | 300 |
| Body | 18,338 | ~9,000 |
| Figures | 6 | 6 |
| Tables | 2 | 2 |
| References | 46 | 46 |

The supplementary document would run to roughly 6,000 words of material already written and already
reproducible from the archived pipeline, so the cut is a redistribution rather than a loss.

## What is genuinely finished

The numerical results, the uncertainty model, the falsification and its scope statement, the figure
set, the reproduction package and the reference list all appear closed. Sixty pipeline scalars match
the text, cross-references are consistent in both directions, and both scripts run clean from the
raw files. No open technical item remains from any of the four review rounds.
