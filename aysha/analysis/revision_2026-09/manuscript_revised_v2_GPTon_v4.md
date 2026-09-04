\## Overall verdict



\*\*Major revision; not yet submission-ready.\*\* The manuscript contains a strong experimental paper, but its most emphatic conclusion—that zero- and one-dimensional single-stream models are structurally incapable of reproducing the receiver and that flow redistribution is required—is broader than the analysis and the project’s current modelling evidence support.



The most defensible contribution is sharper:



> This receiver exhibits a reproducible, flow-dependent apparent exchange response and a front-to-mid wall-temperature inversion, while sparse thermocouple measurements do not uniquely identify physical convective, radiative, or conductive coefficients.



That finding is important, consistent with the project’s progress, and publishable.



\## Priority findings



1\. \*\*The central “structural impossibility” proof is not valid as currently written.\*\*



The manuscript acknowledges that \\(NTU=-\\ln(1-\\varepsilon)\\) is not exact for the measured, non-isothermal wall, but later treats the resulting \\(NTU\\) as the exact integral \\(\\int hP\\,dz/\\dot mc\_p\\) (\[caveat](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision\_2026-09/manuscript\_revised\_v2.md:89), \[claimed proof](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision\_2026-09/manuscript\_revised\_v2.md:274)). Because the axial wall profile changes strongly with flow, the two quantities cannot be equated without calculation.



I performed that calculation using the measured piecewise T8–T12–T11 wall profile. The profile-corrected transfer-unit exponent remains positive, approximately \*\*+0.34 rather than +0.39\*\*, and individual transfer units change by about 1–9%. This is encouraging: the sign result appears robust. That calculation should be incorporated.



Even then, the result rules out only a \*\*flow-independent conductance\*\* or conventional laminar duct closure. It does not rule out every local correlation: an empirical \\(h\\propto Re^{1.34}\\), for example, can reproduce the exponent by construction. Replace “structurally incapable at any calibration” with “inconsistent with conventional single-channel laminar closures under an equal-flow, single-stream assumption.”



2\. \*\*The Graetz-number interpretation is reversed.\*\*



The manuscript states that \\(Gz\_L\\le0.704\\) means the channels remain thermally developing over their entire length (\[line 165](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision\_2026-09/manuscript\_revised\_v2.md:165)). Using the usual estimate,



\\\[

L\_t/L \\simeq 0.05\\,Re\\,Pr\\,D\_h/L=0.05Gz\_L,

\\]



the thermal entrance length occupies at most about \*\*3.5% of the monolith\*\*, not the whole length. The project’s earlier dimensional analysis reached the same conclusion. This also means the benchmark should distinguish the fully developed square-duct limits—approximately 2.98 for uniform wall temperature and 3.61 for uniform heat flux—rather than use an unexplained value of 3.09.



3\. \*\*The paper does not identify a physical macroscopic heat-transfer coefficient.\*\*



The reported Reynolds number assumes all metered flow divides equally among 100 channels, while the interpretation invokes bypass and maldistribution. It is therefore a \*\*nominal Reynolds number\*\*. Likewise, \\(Nu\_{\\rm app}\\) is algebraically constructed from T3 and an axial wall-temperature average; it is not an independently measured interfacial Nusselt number.



This distinction is essential given the project objective of extracting effective convective, radiative, and conductive coefficients. The latest project evidence says:



\- 2D v20 found \\(UA(Re)\\) non-identifiable under the available observations (\[2D audit](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/journal.2D.md:4497)).

\- The apparent HTC in 1D v50 still changes about 11.2% from 15 to 50 cells.

\- v50 is repaired and testable, but explicitly \*\*not calibrated or validated\*\* (\[v50 status](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/theory\_1D\_v50.md:123)).



Use \\(Re\_{\\rm nom}\\), \\(NTU\_{\\rm app}\\), \\(h\_{\\rm app}\\), and \\(Nu\_{\\rm app}\\) consistently. Present them as installed-receiver response indices, not validated transport coefficients.



4\. \*\*T3 uncertainty is flow-dependent, so a uniform ±25 K shift does not protect the exponent.\*\*



The manuscript correctly recognizes that T3 is governed by convection, radiation and stem conduction. Such a bias necessarily varies with flow and enclosure temperature. Applying the same ±25 K shift to every run tests an offset, not a slope error.



This matters because T3 controls effectiveness, apparent NTU/Nu, gas power, inversion thresholds, and the transient-regression abscissa. The v50 work reinforces this problem: its T3 observation model is explicitly flow-dependent, and its observation-only fits drove parameters to bounds rather than identifying a trustworthy correction.



The structural result should be described as robust to a \*\*uniform calibration offset\*\*, not to outlet-probe bias generally. A shielded or aspirated outlet measurement remains necessary to establish the physical gas response.



5\. \*\*The “measured LTNE” claim is too strong.\*\*



T9 and T10 are exposed thermocouples, not direct bulk-gas measurements. The manuscript itself notes radiative perturbation of an internal monolith thermocouple (\[instrumentation](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision\_2026-09/manuscript\_revised\_v2.md:49)). Their difference from the side-wall probes is valuable, but it is an \*\*apparent wall-to-interior temperature disequilibrium\*\*, containing gas convection, radiative loading, sheath conduction and spatial sampling effects.



Rename \\(\\Lambda\\) accordingly and state that it is an observation-weighted LTNE surrogate or bound. “Measured rather than assumed gas–solid nonequilibrium” will invite a decisive reviewer objection.



6\. \*\*The delivered-irradiance closure is neither model-free nor demonstrably steady.\*\*



Equation \\(G\_0A=Q\_{\\rm gas}+K\_{\\rm loss}(\\bar T\_w-T\_{\\rm amb})\\) uses a loss conductance obtained from a lumped transient model and omits storage. Yet the project’s own experimental analysis reports that T2—the cavity/insulation probe—is still charging when the heating runs end (\[existing diagnosis](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/exp\_analysis/exp\_analysis\_report.md:29)). The manuscript currently omits T2 from its instrumentation description.



Consequently:



\- The two irradiance values are not guaranteed upper and lower bounds.

\- Their disagreement cannot presently be assigned to bypass or flow maldistribution.

\- \\(K\_{\\rm loss}\\) from a cooling secant need not equal the steady heating loss, especially with radiative losses and an evolving insulation temperature.



Report these as “radiometric and apparent energy-closure estimates,” add a measured stationarity/storage assessment, and avoid calling their span a delivered-irradiance interval until storage and T3 bias are independently bounded.



7\. \*\*The statistical reduction needs correction before the printed uncertainties are trusted.\*\*



The reduction script does not fully implement the stated uncertainty model:



\- The declared 2.5% flow uncertainty is implemented as a fixed 0.1 sL min\\(^{-1}\\) perturbation, not 2.5% of each run (\[implementation](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision\_2026-09/receiver\_reduction.py:341)).

\- The IEC thermocouple tolerance is divided by two without explaining whether it is being interpreted as a 95% bound (\[implementation](/D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/analysis/revision\_2026-09/receiver\_reduction.py:335)).

\- Eigenvalue-abscissa uncertainty is represented by a common 1% multiplier rather than recomputing it from flow and temperature realizations.

\- The quoted \\(\\Lambda\_{107}\\) slope, \\(8.27\\times10^{-4}\\), is the pooled slope despite the claim of separate flux-dependent intercepts. A fixed-intercept-by-flux regression gives approximately \*\*\\(8.50\\times10^{-4}\\,Re^{-1}\\)\*\*.

\- The global Nu regression hides irradiance-level offsets. The per-level exponents are approximately 1.522, 1.473 and 1.440, each with \\(r^2\\approx0.999\\); a common-slope, flux-intercept model gives about 1.470.



This grouped analysis is scientifically better than the current pooled fit. Also, ±0.069 is a regression standard error, not “run-to-run exponent scatter.”



8\. \*\*Several headline descriptions exceed what the sensors establish.\*\*



The quantity \\(T\_{12}-T\_8\\) detects a \*\*front-to-mid side-wall crossing\*\*. It does not locate the three-dimensional solid-temperature maximum, particularly under a Gaussian beam with expected radial gradients. Replace “the axial peak migrates” with the directly observed statement unless a spatial model or surface map supports the former.



Similarly, the reference-sensitivity table shows that \\(\\varepsilon^\*\\) varies from approximately 0.57 to 0.79 depending on the temperature reference. That is inconsistent with presenting \\(\\varepsilon^\*\\approx0.65\\) as a receiver property. It is a useful operational criterion under the explicitly adopted wall-average convention.



9\. \*\*The transient claims need narrowing.\*\*



The abstract says all eighteen transients collapse within 14%, whereas the Results report a collapse only for fifteen heating runs. Moreover, the wall half-rise is \\(t^\*=0.203\\) while the outlet is \\(0.630\\), and the manuscript recognizes a separate fast wall component. This supports a shared flow-dependent time scale, but not a universal one-state response.



Also test the cooling identification against the matched preceding heating effectiveness rather than the pooled \\(\\varepsilon(q)\\). Using the matched E69/E80/E81 values changes the point estimate from roughly \\(C\_{\\rm eff}=299\\) to \*\*276 J K\\(^{-1}\\)\*\* and \\(K\_{\\rm loss}=0.095\\) to \*\*0.080 W K\\(^{-1}\\)\*\*. That sensitivity belongs beside the sensor/window sensitivity already reported.



\## Publishing and presentation



The abstract is approximately \*\*609 words\*\*, and the main text before the references is approximately \*\*13,600 words\*\*. Both need substantial reduction. Section 5.2 and the pressure-drop literature discussion are especially long and speculative relative to fifteen primary runs.



Other necessary corrections:



\- “Energy-weighted wall temperature” is actually a length- or surface-area-weighted wall temperature.

\- The instrumentation count is inconsistent, and T2 is omitted.

\- Figure 4 says “bounded, not resolved,” while the text reports a locally bracketed 256 kW m\\(^{-2}\\) crossing.

\- The data statement says one script produces all figures, but the README requires both the reduction and figure scripts.

\- Figure 1 is not present in the revision figure directory.

\- The reference-completion note cannot remain in a submitted manuscript.

\- “The two most recent reviews” is obsolete: a dedicated 2023 review of porous volumetric receivers and transient research now exists and should be incorporated (\[Energy Reviews, 2023](https://www.sciencedirect.com/science/article/pii/S2772970223000226)).



\## Recommended disposition



Keep this as an experimental and identifiability paper. Retain the measured temperature slopes, apparent grouped transfer scaling, reference-sensor sensitivity, late-time eigenvalues, and pressure-instrument diagnosis. Recast bypass and viscosity-driven maldistribution as competing hypotheses.



Do not add v49/v50 fitted coefficients to this paper. Reserve physical coefficient extraction and a mechanistic flow-distribution conclusion for the companion modelling study after v50 multi-start calibration, mesh convergence, parameter profiling, and a corrected 2D implementation.



With those changes, the manuscript can become a strong paper. In its current form, the empirical evidence is valuable, but the claims outrun both the measurement model and the project’s audited modelling status. No files were changed during this review.

