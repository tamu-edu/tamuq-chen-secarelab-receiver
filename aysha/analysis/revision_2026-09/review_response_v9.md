# Response to review v9 (manuscript_revised_v3)

Every item is accepted. Two required reruns of the pipeline; one required a correction to the reviewer's own diagnosis, and that correction came from information not in the manuscript.

## 1. Flow-controller uncertainty — accepted, with a correction to the proposed remedy

The reviewer is right that the manuscript attributed 2.5% of reading to a "specification" and that the Aalborg GFC17 specification says something different: ±1.0% of full scale over the whole 0–100% range, with ±0.5% full-scale repeatability. That is an additive bound in absolute flow. The manuscript was wrong to call 2.5% a specification, and §2.2 and §3.4 now state the specification correctly.

The proposed remedy — replace the proportional model with a purely additive one — is not right either, for a reason the manuscript failed to record. The factory calibration is gas-specific, so running air rather than the calibration gas requires a gas-conversion factor whose error is multiplicative. Each unit was therefore calibrated in place against a bubble flowmeter. That substitution removes the conversion factor and replaces it with the residual uncertainty of a directly measured characteristic — also proportional to reading. A purely additive model would discard the term the in-house calibration actually constrains.

The model is now two-term:

$$\sigma_u = \left[(a_{\rm FS} F)^2 + (b_{\rm rel}\, q_u)^2\right]^{1/2},\qquad a_{\rm FS} = 0.0025,\ b_{\rm rel} = 0.025,\ F = 10\ {\rm sL\,min^{-1}},$$

the additive coefficient being the full-scale repeatability read as a 95% bound and the proportional coefficient the calibration residual. One standard-normal draw per unit per realization is combined with each record's own measured shares, as before.

Two points on the reviewer's arithmetic. The units are not 0–5 sL min⁻¹: the largest per-unit reading in the campaign is 5.72 sL min⁻¹, which a 5 sL min⁻¹ unit cannot deliver, and the GFC17 air/N₂ maximum is 10 sL min⁻¹. Per-unit readings therefore span 0.69 to 5.72 sL min⁻¹, or 7 to 57% of full scale, and the additive term at 0.5% FS is 0.025 sL min⁻¹ against a proportional term of 0.017 sL min⁻¹ at the lowest reading — so the additive floor exceeds the proportional term for 48% of the per-unit readings in the campaign, and neither term can be dropped.

The consequence is larger than "regenerate the intervals". The two terms have opposite flow dependence, so the relative uncertainty on the summed flow now falls from 3.4% at the lowest campaign flow to 2.6% at the highest (1.7% to 1.3% for independent units), where the previous model gave a flow-independent 2.5% at every run. A flow-independent relative error cannot tilt a log-log slope by construction, which is why the previous draft could assert that the flow error "displaces a prefactor rather than tilting an exponent". Under the corrected model the error retains a weak gradient that does reach an exponent, and the manuscript now says so.

Regenerated instrumental terms, all point estimates unchanged:

| quantity | old sd | new sd (ρ=1) | new sd (ρ=0) |
|---|---:|---:|---:|
| $Nu$ exponent, pooled | 0.003 | 0.004 | 0.003 |
| $Nu$ exponent, grouped | 0.003 | 0.004 | 0.003 |
| $N_{\rm prof}$ exponent | 0.003 | 0.004 | 0.003 |
| $Nu$ prefactor [10⁻⁴] | 0.11 | 0.12 | 0.11 |
| $\Lambda_{107}$ slope [10⁻⁴] | 0.36 | 0.35 | 0.30 |
| $C_{\rm eff}$, matched [J K⁻¹] | 37 | 37 | 37 |

Two dependent claims were corrected rather than left standing. The two correlation limits no longer "agree to within 0.0002 on every exponent" — they agree to within 0.0007, which is what the manuscript now states. And the regression term no longer exceeds the instrumental term "more than twentyfold" for the pooled correlation: the ratio is 16.8, stated as roughly seventeenfold. For the profile-corrected transfer-unit count it is 11.2, but for the grouped Nusselt correlation it is only 2.6, and §3.4 and S2 now say that instead of claiming an order of magnitude throughout.

**Is the bubble-flowmeter calibration sufficient?** It is sufficient to justify keeping a proportional term, and it is the reason the reviewer's purely additive proposal was declined. It is not sufficient on its own, for three reasons now addressed: it cannot remove the unit's own zero, linearity and repeatability floor, which is additive and dominates at the bottom of this campaign's range; the number was undocumented, so §2.2 now records the calibration and §3.4 the model; and 2.5% was being used as the whole error rather than as one of two terms. If the calibration record gives a residual scatter materially different from 2.5%, `MFC_B_REL` is a single constant in `receiver_reduction.py` and the intervals regenerate in one run.

## 2. T3 band on the identified constants — accepted and implemented

The claim was ahead of the archive: the pipeline propagated the ±25 K band through ε, Nu, the transfer-unit exponent, ε* and η, but not through the identification. It now does, recomputing the abscissa at each offset with the same estimator the archived value uses, so the zero-offset case reproduces the archived point estimate exactly for both branches.

| quantity | −25 K | 0 | +25 K | band |
|---|---:|---:|---:|---|
| $C_{\rm eff}$, matched-ε (primary) [J K⁻¹] | 271.7 | 275.9 | 280.5 | 271.7 – 280.5 |
| $K_{\rm loss}$, matched-ε (primary) [W K⁻¹] | 0.0844 | 0.0798 | 0.0754 | 0.0754 – 0.0844 |
| $C_{\rm eff}$, heating deep probes [J K⁻¹] | 262.2 | 280.7 | 299.4 | 262.2 – 299.4 |
| $K_{\rm loss}$, heating deep probes [W K⁻¹] | 0.1073 | 0.1137 | 0.1200 | 0.1073 – 0.1200 |

The primary capacitance moves 1.6% across the whole band, against an 11% spread among estimators, so it is the quantity least sensitive to the outlet-probe systematic — a result worth having rather than a formality. The loss bracket widens from 0.080–0.114 to 0.075–0.120 W K⁻¹, and §4.4 now reports that.

## 3. Claims contradicted by v2 — resolved by retirement

The v2 items are real and are not fixed in v2. v2 is retired as the working master, as the reviewer recommends; v3 plus the supplement is the submission package, and no claim in v3 rests on v2.

The "scripts run clean" claim was not accurate and is now made true rather than softened. The overflow warnings came from the RK2 integrator inside the multi-start optimizer of §5.2, where the search probes log-h values for which the explicit step is unstable. The local transfer-unit density is now clipped at 1.5/Δζ, i.e. N ≈ 600, where the gas is already at wall temperature to machine precision, so the clip changes no attainable outlet temperature. Every falsification statistic is bit-identical before and after — RMS residuals 64.2/55.4/54.2/54.2, 55.8/54.6 and 36.0/67.7/52.1 K, correlations −0.943 to −0.997 — and the run log is now empty of warnings and archived as `run_reduction.log`.

## 4. Supplementary items — all accepted

The wall-extrapolation table header was reversed: the code keys are `{rear}_{front}` and the table said front/rear. Rather than edit the header, the table is now written by the pipeline from those keys, so the ordering cannot drift again.

"Every table is generated" and "no value is hand-entered" were stronger than the workflow supported — the supplementary tables had been assembled in an analysis session, not by the script. A `write_supplementary_tables()` step now emits `tableS3_heating_conditionality.md`, `tableS4_wall_and_falsification.md` and `tableS5_auxiliary_groups.md`, which the supplement reproduces unaltered, and the claim is restated as "no measured value in this document is hand-entered".

"Bit for bit" is gone. The data-availability section now names the stack the archive was produced with — Python 3.11.16, NumPy 2.4.6, SciPy 1.17.1, pandas 2.3.3, Matplotlib 3.11.1 — and identifies the multi-start optimizer as the one step whose last digits may differ across BLAS builds.

## 5. Wording pass — accepted

"Proves to be" in the abstract is now "is". "Immune to the identifiability problems" is now "does not inherit" them. "The gap measures how little of the structure participates" is now "is a measure of", with the governing assumptions named. "The gas falls further behind the solid" is now "the wall-to-interior probe difference widens", which is what is measured. The five-sixths attribution now states that it is obtained by subtraction and that the identification does not resolve how the remainder divides between holder, insulation and plenum hardware; the second occurrence is now "the greater part".

On the pressure measurement the reviewer is right and the recommendation is strengthened, not merely hedged: the two mechanisms differ in the flow dependence of the characteristic, not in its value at one flow, so §5.3 and §6 now specify a pressure-drop series across the flow range with a sealed-clearance control run, and state explicitly that a single point would not discriminate them.

The web Biot number claim in the conclusion now excludes "a through-web conduction resistance as the limiting one, though not every solid-side or assembly-scale resistance".

## 6. The 1D/2D boundary

We concur, and nothing changed. The full project models remain contextual discussion and companion work; the only thermal integrations in the paper are the single-stream diagnostics of §4.2 and §5.2, used to test the model class the paper argues against.

## Status

Sections 1–6 are 10,671 words, up 467 from v3 because §2.2 and §3.4 gained the controller specification, the calibration provenance and the two-term model, and §4.4 the identification band. All 46 references remain cited in the main text, the pipeline runs without warnings, and every number quoted in the manuscript reproduces from the regenerated archive.
