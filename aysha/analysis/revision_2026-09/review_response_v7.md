# Response to review v7 of manuscript_revised_v2

We accept every finding in this review. Nine were matters of work not yet done rather than of
disagreement, and three changed a result. We note with thanks that the reviewer independently
reconstructed the shared-conductance test and obtained 64, 56, 55 and 54 K at two, three, five and
seven nodes with a residual slope of −124 to −126 K per ln Re; our archived values are 64.2, 55.4,
54.2 and 54.2 K with slopes of −148 to −126 K per ln Re. That independent reproduction is the
strongest evidence we have that the result is real, and we have acted on the reviewer's judgement
that it was not yet properly archived or properly bounded.

All values below are read from `outputs/results.json` and `outputs/uncertainty.csv` produced by a
single 4000-realization run of `receiver_reduction.py`. Line numbers refer to
`manuscript_revised_v2.md` as submitted with this response.

## 1. The fixed-conductance test — accepted; broader claim retained following the additional tests

The reviewer is right that the test existed only in scratch code, that the data-availability claim
was therefore false, and that two stronger variants were owed. All three are now in the pipeline as
`fixed_profile_test()`, archived under `results.json → fixed_profile_test` with fitted nodes and
per-run residuals, and the wall-extrapolation sweep is archived as `wall_extrapolation`.

Both requested variants were run, and both fail in the same way.

*Per-configuration profiles.* Fitting a separate five-node profile within each irradiance
configuration gives 36.0, 67.7 and 52.1 K rms at 256, 304 and 456 kW m⁻², closely matching the
reviewer's independent 36, 68 and 53 K. Allowing each configuration its own profile does not rescue
the hypothesis; it sharpens the diagnosis, because the residual ordering becomes nearly perfect,
r = −0.996, −0.982 and −0.997. Within one lamp configuration, at fixed flux distribution and
comparable temperature level, a flow-independent profile fails monotonically in flow.

*Shared Nu(z).* The direct test of the constant-Nusselt class, h_i(z) = Nu(z) k_i/D_h with each run's
own conductivity, gives 55.8 and 54.6 K rms at three and five nodes with slopes of −131.6 and
−129.1 K per ln Re. It is indistinguishable in outcome from the dimensional family. We note that our
first attempt at this fit did not converge, because the multi-start range had been scaled for
dimensional h and Nu is some two hundred times smaller; the fit is now started from the
uniform-conductance solution for the runs being fitted, and the correction is documented in the
code.

The manuscript now states the claim at the scope the evidence supports (§4.2, lines 292 ff.):
none of the tested families — shared h(z) at two to seven nodes, shared Nu(z) at three and five,
and per-configuration h(z) at five — reproduces the outlet data, with the two structural limits
stated explicitly, namely that the single-stream form itself is assumed and that the wall
reconstruction is inherited.

On "exactly Re⁻¹": the §4.2 and conclusion instances the reviewer flagged are gone; the computed
requirement, −0.9996 for a fixed Nusselt number and −1.017 for a temperature-independent
conductance, is now the only form used.

## 2. Controller uncertainty — accepted; the specification was being used twice

This is the most substantive correction in the round and the reviewer is entirely right. Applying a
full-variance common component and a full-variance independent component gave each unit
√2 × 2.5%, which no accuracy specification supports.

Each controller's error is now written as a correlated part of weight √ρ and an independent part of
weight √(1−ρ), so the per-unit total is 2.5% for any ρ. Since the manufacturer states no
correlation between units, the pipeline runs both limits and archives both
(`results.json → controller_correlation_bounds`; `uncertainty.csv` carries a `rho_mfc` column).
The instrumental terms are consequently reported as bounds (§3.4, line 121; §4.2, line 157):

| Quantity | ρ = 1 (correlated) | ρ = 0 (independent) |
|---|---|---|
| Nu exponent, pooled | ±0.003 | ±0.005 |
| Nu exponent, grouped | ±0.003 | ±0.005 |
| N_prof exponent | ±0.003 | ±0.005 |
| C_eff matched-ε | ±37 J K⁻¹ | ±40 J K⁻¹ |

We have also removed the claim that the independent part is run-to-run repeatability. It is the
between-unit part of the accuracy specification, and it varies between runs only because the shares
do (§3.4, line 123).

The three transient-branch inconsistencies are fixed. The cooling flow now carries the independent
component as well as the common one, using the mean steady share factor since per-cooling-run
shares are not logged. The T3 and ambient calibration offsets now reuse the steady frame's
per-sensor draws instead of being redrawn. The matched-effectiveness abscissa now takes each
cooling run's own temperature perturbation rather than the mean across eigenvalue rows.

The archived Monte Carlo is now 4000 realizations, as the text states. The 1200-realization file the
reviewer found was a working run left in place; that was our error.

## 3. Matched estimator — accepted; the stale material is removed

Methods no longer present the pooled effectiveness as the one in use (§4.4, line 191). The loss
bracket reads 0.080 to 0.114 W K⁻¹ in both places the reviewer identified. The primary values are
C_eff = 276 ± 37 J K⁻¹ and K_loss = 0.080 ± 0.022 W K⁻¹ at ρ = 1 (§4.4, line 191).

Our statement that these were synchronized everywhere was wrong, and the reason is worth recording:
the audit we ran checked the values we had changed rather than sweeping the document, so it could
not find the places we had missed. The audit has been redone in the other direction, over sixty
pipeline-sourced scalars matched against `results.json` and `uncertainty.csv`, and it now also
covers the fixed-profile and extrapolation results. One further discrepancy surfaced that the
review did not list: the inverse Graetz number at the lowest flow was written 5.72 against a
computed 5.71.

## 4. Mechanistic language — accepted

§4.3 (line 185) now says the two observations point to increasing participating contact taken at
face value, and defers the conditional to §5.2. §4.6 (line 227) no longer treats bypass or
maldistribution as what remains: it states that the discrepancy requires a flow-dependent term,
that this does not identify which, that a stationarity departure or a flow-dependent outlet
response would supply one equally well, and that the convergence with §5.2 is a consistency rather
than an independent confirmation. The manuscript nowhere claims the project has selected a
mechanism.

## 5. Synchronization — accepted

Conclusion loss bracket and irradiance spans corrected to 0.080–0.114 W K⁻¹ and 451–517, 304–402,
201–256 kW m⁻² (line 344). Table 2 now carries the 95% percentile intervals its description
promised, alongside the standard deviation and, for fitted slopes, the regression standard error.
Figure 3a now plots the three parallel grouped lines as the primary fit, colour-matched to the
configurations, with the pooled single-prefactor line dashed behind them.

## 6. Observational wording — accepted

The §3.1 definition (line 101) now describes the indicator as negative while the front-face probe is
the hotter and positive once the mid-depth probe has overtaken it, with no reference to a migrating
peak. §4.1 no longer converts the crossing into the volumetric effect: it reports the side-wall
signature of volumetric behaviour and states that the measured quantity is a difference between two
side-wall probes rather than the location of the internal maximum. The conclusion (line 334) uses
"front-to-mid side-wall crossing", "apparent wall-to-interior disequilibrium" and "operational
marker under the adopted convention", and no longer says the inversion is governed by a single
dimensionless criterion.

## 7. Energy closure — accepted

"Reported as intervals" is now "reported as spans between the two estimates rather than as intervals
containing the true value" (line 229), and the exhaustive causal inference in §4.6 is replaced as
described under item 4.

## Additional items

*Transducer.* Accepted, and the datasheet figure was already in the sentence: the Keller PD33X is
±200 mbar range at ±0.1% full scale, which is where the ±0.2 mbar comes from. Only the word
"resolution" was wrong, and it now reads "an accuracy floor of ±0.2 mbar" (line 51). The caveat in
our previous letter, that the figure needed confirmation against the datasheet, is withdrawn as
unnecessary.

*Limitations.* Accepted; the section was stale in both respects the reviewer names (§5.4, line 326).
It now states that the isothermal-wall bias is quantified at 1.5% to 8.7% in the transfer units and
0.048 in the exponent, that the correction inherits an extrapolated wall reconstruction, and that
two single-stream thermal models are in fact solved — deliberately of the class the paper argues
against, and used to test that class rather than to correct any measurement. Two limits not
previously listed are added: thermal stationarity is assumed rather than demonstrated, and a
flow-dependent outlet-probe bias is unbounded and would act directly on the exponents.

*Letter tally and tone.* Both accepted. The previous letter's arithmetic was inconsistent, and this
letter states no tally. We have adopted the reviewer's formulation for item 1 — accepted, broader
claim retained following an additional test — which is a more accurate description of what happened
than "contested", since the criticism was correct and the remedy was the reviewer's own suggestion.
Passage locations are given beneath each item as requested.

## On the reviewer's proposed high-level conclusion

We adopt it, with one amendment. The reviewer's formulation reads:

> Within the tested single-stream formulations and measured side-wall reconstruction, neither a
> uniform conductance nor the tested two- to seven-node flow-independent axial conductance profiles
> reproduce the outlet-temperature flow response. The physical cause remains unidentified.

Both clauses are right and the manuscript now says both. The amendment is that the tested set is
larger than the reviewer had seen: it also includes a shared Nusselt profile and per-configuration
profiles, the latter being the variant most likely to rescue the hypothesis and the one that fails
most sharply. We would therefore write "neither a uniform conductance, nor a shared conductance or
Nusselt profile of two to seven nodes, nor a separate profile fitted within each irradiance
configuration".
