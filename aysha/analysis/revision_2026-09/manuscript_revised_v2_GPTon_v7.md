Excluding the abstract and overall length as before, **the response letter and v7 are not yet ready to submit**. The underlying revision is stronger, and I independently reproduced the headline shared-\(h(z)\) failure approximately, but most response items are only partially closed because the manuscript, code and archived outputs do not yet agree.

## Point-by-point audit

| Response item | Verdict | Main reason |
|---|---|---|
| 1. Structural proof | Partially resolved | New test is numerically credible, but absent from the reproducibility pipeline and still stated too universally. |
| 2. Uncertainty | Partially resolved | Regression/MC separation is fixed, but the controller error model double-counts the 2.5% specification and transient propagation remains inconsistent. |
| 3. Transient estimator | Partially resolved | Matched effectiveness is correctly primary, but Methods, loss bracket and conclusions retain old values. |
| 4. Causal language | Partially resolved | Section 5.2 is improved, but earlier sections still identify flow distribution too strongly. |
| 5. Numerical synchronization | Not fully resolved | Tables are now generated, but several manuscript values remain stale and archived MC contains 1200—not 4000—realizations. |
| 6. Observational qualifications | Partially resolved | Figures are corrected, but “peak migration” and directly measured nonequilibrium remain in the text and conclusion. |
| 7. Energy closure | Partially resolved | The conceptual qualification is fixed, but the conclusion reports superseded spans and conductances. |

## 1. The new shared-\(h(z)\) test is valuable, but not yet a general proof

The response correctly accepts the former logical gap and adds a simultaneous fit of a nonnegative, piecewise-linear \(h(z)\) ([response:9](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/review_response_v6.md:9); [manuscript:290](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/manuscript_revised_v2.md:290)).

I independently reconstructed this test. I obtained approximately:

- 2 nodes: 64 K RMS
- 3 nodes: 56 K RMS
- 5 nodes: 55 K RMS
- 7 nodes: 54 K RMS
- residual slope: approximately −124 to −126 K per \(\ln Re\)

These agree closely with the response, so the result appears numerically genuine.

However:

- No shared-profile fitting code exists in [receiver_reduction.py](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/receiver_reduction.py).
- No corresponding result exists in `results.json`.
- Rear-extrapolation sensitivity is likewise absent from the pipeline.
- Therefore the data-availability claim that every scalar is generated and archived is currently false ([manuscript:346](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/manuscript_revised_v2.md:346)).

The test should also be strengthened in two ways:

1. Fit a separate fixed profile within each irradiance configuration. This avoids rejecting a profile merely because conductance differs with temperature or lamp configuration. My independent five-node fits still left approximately 36, 68 and 53 K RMS at 256, 304 and 456 kW m\(^{-2}\), respectively, suggesting the conclusion will survive.

2. Test a shared \(Nu(z)\), with \(h_i(z)=Nu(z)k_i/D_h\), not only a shared dimensional \(h(z)\). That is the more direct test of the constant-Nusselt model class.

Until those tests are archived, change “no such profile” to “none of the tested two- to seven-node profiles.”

The correction from exactly −1 is properly introduced in Section 5.2, but old “exactly \(Re^{-1}\)” statements remain in Section 4.2 and the conclusion ([line 169](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/manuscript_revised_v2.md:169); [line 334](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/manuscript_revised_v2.md:334)).

## 2. The revised uncertainty model still needs correction

The separation of instrumental Monte Carlo uncertainty from regression SE is now clear and appropriate.

The controller model, however, uses the same 2.5% accuracy specification twice:

- a common 2.5% calibration multiplier; and
- another independent 2.5% per-controller contribution, reduced to approximately 1.3% after summation.

That gives an individual-controller uncertainty equivalent to approximately \(\sqrt{2}\times2.5\%\). An accuracy specification cannot be decomposed into two full-variance components without independent evidence.

A defensible approach would be either:

- draw four persistent controller errors with a stated covariance matrix; or
- report uncertainty bounds for fully correlated and independent-controller cases.

The independent component should not be treated as run-to-run repeatability unless such repeatability is specified.

The transient branch also remains inconsistent with the description:

- cooling flow receives the common factor but not the proposed independent run component;
- T3/ambient errors are redrawn rather than sharing the steady-data sensor calibration draws;
- the matched-effectiveness calculation uses the mean temperature perturbation across eigenvalue rows for every cooling case ([code:408–435](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/receiver_reduction.py:408)).

Finally, the archived [uncertainty.csv](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/outputs/uncertainty.csv:1) contains **1200 realizations**, whereas the manuscript and response state 4000.

## 3. The matched transient estimator is correct, but old material remains

Adopting the matched-effectiveness result as primary is well justified:

\[
C_{\rm eff}=276\pm38\ {\rm J\,K^{-1}},\qquad
K_{\rm loss}=0.080\pm0.024\ {\rm W\,K^{-1}}.
\]

The one-state contradiction is also correctly removed from the conclusion.

Remaining inconsistencies:

- Methods still state that \(\varepsilon\) is obtained from the pooled campaign correlation ([line 109](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/manuscript_revised_v2.md:109)).
- Section 4.4 still gives the loss bracket as 0.098–0.118 W K\(^{-1}\), rather than the current 0.080–0.114 ([line 197](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/manuscript_revised_v2.md:197)).
- The conclusion retains 0.098–0.118 ([line 332](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/manuscript_revised_v2.md:332)).

The response’s claim that these values are synchronized everywhere is therefore incorrect.

## 4. Mechanistic language remains stronger than current project evidence

Section 5.2 is substantially improved by distinguishing observation from hypothesis ([line 296](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/manuscript_revised_v2.md:296)).

But two passages retain the old causal conclusion:

- Section 4.3 says both observations “require” participating contact to increase ([line 185](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/manuscript_revised_v2.md:185)).
- Section 4.6 says what remains is bypass or maldistribution and calls this independently required ([line 227](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/manuscript_revised_v2.md:227)).

These should be conditional, matching Section 5.2.

This matters because current project progress still says:

- 1D v50 is repaired and testable, but not calibrated or physically explanatory;
- the power-dependent T3 response remains unresolved;
- the trustworthy 2D result establishes \(UA(Re)\) non-identifiability, not unique maldistribution.

Thus the manuscript may argue that bypass and redistribution are leading hypotheses, but not that the project has selected them.

## 5. Numerical synchronization is incomplete

The new generated [Table 2](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/outputs/table2_constants.md:1) is much improved. Nevertheless:

- Conclusion loss bracket: 0.098–0.118 instead of 0.080–0.114.
- Conclusion irradiance span: 456–517 instead of 451–517 at 456 kW m\(^{-2}\).
- Conclusion irradiance span: 216–256 instead of 201–256 at 256 kW m\(^{-2}\) ([line 338](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/manuscript_revised_v2.md:338)).
- Table 2 does not include the stated 95% intervals, despite the table description claiming it does.
- Shared-profile and extrapolation results are absent from `results.json`.
- Figure 3 labels the grouped model as primary, but the plotted black fit remains the pooled single-prefactor line. Three parallel grouped lines should be shown.

The response’s claimed “36-value closing audit” should be removed or repeated after these corrections.

## 6. Observational wording is still inconsistent

Figure 4 is now appropriately labelled. The following manuscript statements remain:

- The definition still says the peak is on the front and then migrates inward ([line 101](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/manuscript_revised_v2.md:101)).
- Section 4.1 still converts the observed side-wall crossing into “the volumetric effect” ([line 149](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/manuscript_revised_v2.md:149)).
- The conclusion again calls the quantity measured “local nonequilibrium” and says it is measured rather than assumed ([line 330](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/manuscript_revised_v2.md:330)).
- The conclusion says volumetric inversion is “governed by a single dimensionless criterion” ([line 328](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/manuscript_revised_v2.md:328)).

Use “front-to-mid side-wall crossing,” “apparent wall-to-interior disequilibrium,” and “operational marker under the adopted convention” consistently.

## 7. Energy-closure response is conceptually right but numerically stale

The closure table itself now uses the appropriate conditional-estimate language and updated values. The conclusion correctly says the two estimates do not bound the irradiance.

But the conclusion still reports the superseded spans, and Section 4.6 still makes an exhaustive causal inference. It also says power-normalized quantities are reported as “intervals” at [line 229](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/manuscript_revised_v2.md:229). Use “reported spans” throughout.

## Additional response-letter issues

- The pressure correction was not actually completed: the Methods still call ±0.2 mbar a “resolution” ([line 51](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/manuscript_revised_v2.md:51)).
- The limitations section is stale. It says the non-isothermal-wall bias is unquantified and that no result uses a thermal model, although the manuscript now quantifies profile correction and solves two single-stream thermal models ([line 320](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision_2026-09/manuscript_revised_v2.md:320)).
- The opening tally in the response letter is internally confusing: nine fully accepted plus two differently resolved plus one contested totals twelve, while it describes seven principal and four additional items.
- “The one item contested” is unnecessarily adversarial because the original criticism was accepted and the suggested test was performed. Better: “Accepted; broader claim retained following an additional test.”
- For an external response letter, include the exact revised passage and manuscript page/line location beneath each response.

## Overall decision

The new fixed-profile calculation materially strengthens the paper and appears numerically plausible. However, its code and outputs must be archived, the controller uncertainty model needs reconsideration, and the remaining stale statements require a global consistency pass.

The paper’s defensible high-level conclusion remains:

> Within the tested single-stream formulations and measured side-wall reconstruction, neither a uniform conductance nor the tested two- to seven-node flow-independent axial conductance profiles reproduce the outlet-temperature flow response. The physical cause remains unidentified.

That formulation is consistent with both the new empirical test and the current 1D/2D project status.