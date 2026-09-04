# Response to review of manuscript_revised_v2 (v6)

Every checkable claim in this review was tested against the data and the reduction code before
being acted on. Of the seven highest-priority findings and four additional corrections, nine are
accepted in full, two are accepted with a different resolution than the one proposed, and one is
contested on substance. Two of the accepted items changed the paper's central argument rather than
its wording, and one of those changes strengthens it.

## 1. Scope of the structural proof — ACCEPTED, and resolved by the test the review proposes

The review is right that the profile inversion assumes a spatially uniform h when solving for
N = hPL/(m-dot c_p), and therefore on its own establishes incompatibility with a uniform
constant-Nusselt closure rather than with every fixed h(z). This was a genuine gap.

Rather than narrow the claim, we ran the second resolution the review offers, which is decisive.
Taking h(z) as piecewise linear and non-negative on evenly spaced nodes, a single shared profile
was fitted to all fifteen runs simultaneously -- one h(z) common to every flow, which is what
"fixed" means -- and tested against the fifteen measured outlet temperatures. No such profile
reproduces them. With two nodes the best fit leaves 64 K rms outlet residual; with three, five and
seven nodes it saturates at 55, 54 and 54 K, so the failure is not a want of freedom. The residuals
are monotone in flow, correlating with ln Re at r = -0.95 and falling 126 K per e-fold, from about
+80 K at the lowest flow of each configuration to about -70 K at the highest. A uniform outlet-bias
offset cannot produce a trend in flow, and the thermocouple uncertainty is one to two orders of
magnitude smaller than the residual. The general fixed-conductance hypothesis is therefore excluded
empirically, and section 4.2 now rests the broad claim on this test instead of on the uniform-h
inversion. Multi-start optimisation was used throughout; the best-fitting profiles push all
conductance to the channel rear, which is where the family best imitates a rising participation,
and they still fail.

The rear-extrapolation concern is also accepted and quantified. Replacing the constant rear
extrapolation with a continuation at half the measured T12-to-T11 gradient moves the exponent from
+0.341 to +0.331; a linear front back-extrapolation moves it by under 0.003. A full linear
continuation of the rear gradient is not admissible: for five of the fifteen runs it places the exit
wall below the measured outlet gas temperature, which no single-stream model can produce, so the
outlet measurement itself bounds how steeply the rear wall can fall. Sensitivity is thus at the 0.01
level against a discrepancy of 1.34.

On "exactly Re^-1": accepted, and resolved by computing the requirement rather than asserting it or
substituting (Re Pr)^-1. From the measured k, c_p and m-dot, a fixed Nusselt number requires
d ln N / d ln Re = -0.9996 with r^2 = 0.9999; a conductance independent of temperature as well as of
flow requires -1.017. Pr varies only over 0.6834-0.6889, which is why the value is so close to -1.

## 2. Uncertainty statement versus implementation — ACCEPTED in full

All four sub-points were verified in the code and all four were wrong as implemented or as described.

Regression standard errors and instrumental Monte Carlo uncertainties are now reported separately
and labelled, in the text and in Table 2. They differ by an order of magnitude and answer different
questions: for the profile-corrected exponent the instrumental term is +/-0.004 and the regression
term +/-0.040. The "without exception" sentence in section 3.4 has been replaced by an explicit
statement of which quantity carries which, and of the fact that for the section 4.2 correlations the
regression term is the larger and is the one to read as the uncertainty of the correlation.

The four-controller error is corrected. The metered flow is the sum of four Aalborg units carrying
roughly 25, 33, 28 and 14% of the total; the text described them as one. The 2.5% specification is
now applied in two components: a common multiplier shared by all four units and all runs
(calibration accuracy, treated as fully correlated -- the conservative choice -- and displacing a
prefactor rather than an exponent), plus an independent per-run component of
2.5% x sqrt(sum f_i^2) = 1.3% at the measured shares, which does act on an exponent. The previous
single common multiplier suppressed the component that reaches the exponent.

The thermocouple class error is now drawn once per sensor per realization and shared across runs,
scaled by each run's temperature-dependent tolerance, with the 0.5 K drift independent per run. It
was previously drawn independently for every run, which contradicted calling it systematic.

One further defect was found while implementing this and is not in the review: the common controller
factor was being drawn independently for the steady frame and for the transient abscissa, though it
is one physical error. Sharing it moved the Monte Carlo mean of the matched capacitance from 302 to
286 J K^-1.

## 3. Transient estimator — ACCEPTED in full

"For continuity" was not a defensible reason and the matched-effectiveness identification is now
primary: C_eff = 276 +/- 38 J K^-1 and K_loss = 0.080 +/- 0.024 W K^-1 at r^2 = 0.986, against 0.964
for the pooled form. Its uncertainty is propagated in the Monte Carlo, and the energy-closure
bracket has been recomputed from it, which changes the closure factors to 0.989-1.134, 1.149-1.324
and 0.784-0.916 at 456, 304 and 256 kW m^-2 and the closure estimates to 451-517, 349-402 and
201-234 kW m^-2. The conclusion that assembly capacitance substantially exceeds monolith
capacitance survives, as the review anticipated: 276 J K^-1 against 42-47 J K^-1.

One point of presentation: the Monte Carlo mean of the matched capacitance (286 J K^-1) sits above
its point estimate (276) because C_eff is the reciprocal of a slope fitted to three points and is
therefore a nonlinear function of it. We quote the point estimate with the Monte Carlo standard
deviation and state the bias rather than averaging it away.

The one-parameter contradiction is fixed: the conclusion now says a single slow time scale organizes
late behaviour and explicitly denies a one-state response, matching section 4.5.

## 4. Causal language on the mechanism — ACCEPTED

"What remains is that the participating contact must grow" is now explicitly conditional on the
measured flow dependence being physical, and the two alternatives that would produce the same
signature without any change in participation are named in the same sentence: a flow-dependent
outlet-probe bias, and a departure from thermal stationarity varying with run duration and hence
with flow. Neither is bounded by this campaign, which is why the two measurements ranked first and
second in section 5.3 target them. Bypass and viscosity-driven redistribution are now described as
the leading hypotheses consistent with the observations rather than as identified causes, and
"a flow-distribution degree of freedom supplies one" has become an argument for preferring a
hypothesis with an explicit statement that preferring is not confirming.

We have not added optical deposition to that list, since nothing in this campaign speaks to it and
naming an unevidenced mechanism alongside two evidenced ones would misrepresent the balance.

## 5. Numerical synchronization — ACCEPTED, and fixed structurally rather than by hand

The review is right on every instance, and right that a single-source audit was needed. Rather than
reconcile by hand, table generation has been moved into the reduction script, so Tables 1 and 2 are
now written by the same pass that computes the numbers and cannot drift again. Table 2 previously
carried the superseded +0.389 exponent and the pooled Lambda_107 slope because it was
hand-maintained; it now reports point estimates with instrumental Monte Carlo standard deviations
and regression standard errors labelled separately, on one consistent basis.

The conclusion and Figure 3 now lead with the grouped exponent 1.470 and label the pooled 1.443 as
such. The joint capacitance is 269 +/- 31 everywhere. The deep-probe estimate is 281 J K^-1 as a
point estimate with 37 J K^-1 as its Monte Carlo standard deviation, and the 284 and 286 values -- a
Monte Carlo mean and a stale figure respectively -- are gone. A closing audit checked 36 quoted
values against results.json and uncertainty.csv, all of which now agree, and a reverse sweep for
superseded strings returns nothing.

## 6. Observational qualifications — ACCEPTED in full

The Lambda definition in section 3.1 and the conclusion both now describe an apparent
wall-to-interior disequilibrium and a bounded surrogate for gas-solid nonequilibrium, matching
section 4.3. "Peak migration" is replaced by "front-to-mid side-wall crossing" throughout, including
in section 5.2 and in the Figure 4a title, and Figure 4c's axis label now reads "apparent
wall-to-interior deficit". The threshold is described as an operational marker under the stated
wall-averaging convention wherever it appears. Figure 4a's "bounded, not resolved" annotation, which
contradicted the locally bracketed crossing reported in the text, now reads "single negative point
only", which is what the 256 kW m^-2 group actually shows.

## 7. Energy closure in the conclusion — ACCEPTED

The conclusion no longer says the irradiance is "determined to an interval". It states that the
irradiance is not determined to a value and that the two estimates do not bound it, matching section
4.6, and the closure table's column headings now read "closure estimate" and "span reported" rather
than "interval adopted".

## Additional corrections

The section 3.1 / section 4.2 tension is accepted and resolved by reframing rather than removal. The
"category error" wording is replaced by the specific statement that comparing this exponent against
channel-correlation values near 0.3-0.6 compares unlike objects, and the packed-bed comparison in
section 4.2 is now explicitly justified on the grounds that those coefficients were themselves
identified from measurement on a porous assembly rather than prescribed from channel geometry. The
comparison no longer claims our value should agree with theirs, only that a superlinear exponent is
a recurring outcome of measuring rather than assuming this coefficient at low Reynolds number.

The pressure transducer wording is accepted with a caveat. "Resolution" has been replaced by an
accuracy floor of order +/-0.2 mbar, noted as 0.1% of the +/-200 mbar span. We have not seen the
datasheet, so this figure should be confirmed against it before submission; if the manufacturer
states a different accuracy class the sentence needs revising, though the conclusion that the
transducer cannot resolve a 0.2-1.0 mbar signal is unlikely to change.

Figure 3 and Figure 4 labels are corrected as described under items 5 and 6.

## The one item contested

The review's framing of item 4 -- that the 2D endpoint establishes structural non-identifiability
"not unique maldistribution, bypass, or 3D physics" -- is accepted and implemented. But we do not
accept the implication that the positive transfer-unit slope should be presented only as evidence
against an axially uniform constant-Nu closure. That was the correct reading of the previous
version, in which the only test performed assumed uniform h. It is no longer the state of the
evidence: the shared-profile fit described under item 1 tests the general fixed-conductance
hypothesis directly, with up to seven free non-negative nodes, and excludes it with a residual
structure that is systematic in flow. The broader claim is now an empirical result rather than an
inference from a special case, and it is stated as such.

Two limitations of that test are stated in the manuscript and are worth restating here. It assumes
the single-stream form itself -- one gas stream exchanging with a wall whose temperature is the
measured profile -- so it excludes fixed conductances within that form, not every conceivable
reduced-order model. And it inherits the wall reconstruction, whose rear extrapolation is bounded
but not measured. Neither affects the sign of the result.
