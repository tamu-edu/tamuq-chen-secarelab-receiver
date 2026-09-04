# Supplementary Material

## Transient Thermal Characterization and Dimensionless Correlations for a Structured Silicon Carbide Volumetric Solar Receiver, with Consequences for Reduced-Order Model Validation

A. Melhim, A. G. Konstandopoulos, K. J. Kakosimos

This document holds the campaign-by-campaign literature detail, the full uncertainty propagation, and the sensitivity and per-run tables supporting the main text. Section numbers S1 to S8 are those cited in the manuscript. The tables in sections S3, S4 and S5 are written by `receiver_reduction.py` as `tableS3_heating_conditionality.md`, `tableS4_wall_and_falsification.md` and `tableS5_auxiliary_groups.md` and are reproduced here unaltered; no measured value in this document is hand-entered.

---

## S1 Experimental campaigns on porous and structured volumetric absorbers

The response has largely been to build better absorbers. Fend and co-workers measured thermophysical and heat transfer properties of cordierite, SSiC and CBSiC foams, SiC fibre meshes and catalyst carriers, obtaining efficiencies of 75–92% at 300–550 °C for combined material systems and showing that 45 ppi foams transfer far more heat than 10–20 ppi foams by virtue of specific surface areas of up to 2500 m² m⁻³ [4]; they also reported double-layer SiC foams reaching 98% efficiency at 672 °C outlet [5]. Ávila-Marín and colleagues built and tested twenty-six graded-porosity wire-mesh absorbers against the TSA and SOLAIR baselines [6], and later five periodic and non-periodic SiC morphologies including Voronoi structures that outperformed a conventional SiC honeycomb [7]. Zaversky and co-workers optimized single- against multi-layer ceramic foams and concluded that only the first layer's properties matter thermally, with a broad efficiency plateau between 30 and 50 PPI [8]. Wang and co-workers tested four porous SiC absorbers and found that reducing pore diameter from 5.83 to 2.22 mm at fixed porosity raised efficiency by nearly three points [9]. The pattern continues in the most recent work. Patil and co-workers built and tested a 5 kW reticulated-porous-ceramic cavity receiver for process heat above 1000 °C, comparing SiSiC, alumina and ceria structures across mean pore sizes [31], and Broeske and co-workers tested screen-printed ceramic and metallic-foil three-dimensionally shaped absorbers in a sixteen-module matrix at DLR's Synlight facility, reaching air exit temperatures to 735 °C and thermal efficiencies up to 91.8% at 650 °C against 84.6 to 85.3% for the 80 cpsi HiTRec SiSiC honeycomb benchmark [29]. D'Souza and co-workers took the complementary route of a calorimetric test bench, measuring the heat transfer to pressurised gas in absorbers of varying rectangular channel geometry under a high-flux solar simulator [32]. The recurring output of this family of work is a design recommendation — a porosity, a cell density, a layering, a channel aspect ratio — obtained from steady-state efficiency measurements.

The recurring output of these campaigns is a design recommendation — porosity, cell density, layering, channel aspect ratio — inferred from steady-state efficiency and outlet-temperature measurements. Section 5.1 of the manuscript quantifies what such a dataset constrains.

---

## S2 Uncertainty propagation

Three kinds of uncertainty are propagated and reported separately, because they behave differently.

**Random and calibration error** is propagated by Monte Carlo over 4000 realizations. Reported ± values are Monte Carlo standard deviations and bracketed ranges are 95% percentile intervals. The type N tolerance of IEC 60584 class 1 — the greater of 1.5 K and 0.004|t| — is a maximum permissible error rather than a standard deviation, so it is read as a 95% bound and halved to a 1σ term. That term is a per-sensor calibration offset, so it is drawn once per sensor per realization and shared across all runs, scaled by each run's temperature-dependent tolerance, with an independent 0.5 K steady-window drift term per run. Viscosity and conductivity carry 2% and 1% perturbations. The abscissa of the transient identification, x = ε ṁ c_p, is rebuilt in every realization from that realization's perturbed flow and temperatures, with the controller factor shared between the steady frame and the abscissa because it is one physical error.

**The controller model is two-term.** The four Aalborg GFC17 units carry roughly 25, 33, 28 and 14% of the total flow, and two independent pieces of information constrain their error. The manufacturer specification is ±1.0% of full scale over the whole 0 to 100% range with ±0.5% full-scale repeatability — an additive bound in absolute flow, not a percent of reading — and the factory calibration is gas-specific. Because the receiver runs air rather than the calibration gas, each unit was calibrated in place against a bubble flowmeter; that substitution removes the gas-conversion factor, which is a percent-of-reading systematic, and replaces it with the residual uncertainty of the measured characteristic, which is also proportional to reading. What it cannot remove is the unit's own zero, linearity and repeatability floor, which remains additive. Both terms are therefore carried,

$$\sigma_u = \left[(a_{\rm FS} F)^2 + (b_{\rm rel}\, q_u)^2\right]^{1/2},\qquad a_{\rm FS} = 0.0025,\ b_{\rm rel} = 0.025,\ F = 10\ {\rm sL\,min^{-1}},$$

with the additive coefficient the full-scale repeatability read as a 95% bound and the proportional coefficient the calibration residual. One standard-normal draw per unit per realization is combined with each record's own measured shares, so a single physical unit error reaches every steady and transient record at the magnitude its own reading implies.

The two terms have opposite flow dependence, and in this campaign neither is negligible: per-unit readings run 0.69 to 5.72 sL min⁻¹, that is 7 to 57% of full scale, and the additive term exceeds the proportional one for 48% of the per-unit readings. On the summed flow the relative 1σ falls from 3.4% at the lowest campaign flow to 2.6% at the highest for fully correlated units, and from 1.7% to 1.3% for independent ones. A purely proportional model would give a flow-independent 2.5% and would therefore, by construction, be unable to reach a fitted exponent at all; the correct model is predominantly systematic, so it displaces a prefactor, but retains a weak gradient that does reach an exponent. Because the manufacturer states no correlation between units, every instrumental term is computed in both limits: they agree to within 0.0007 on every exponent and within 0.5 J K⁻¹ on every capacitance, so the larger is quoted and both are archived in `uncertainty.csv` under `rho_mfc`.

Where a quantity is a fitted slope, its regression standard error is reported alongside and labelled as such, because the two measure different things: the Monte Carlo term is what the instruments contribute and the regression term is the scatter of the runs about the fitted law. For the pooled Nusselt correlation the regression term is larger by a factor of 16.8 and for the profile-corrected transfer-unit count by 11.2, so in both cases it is the one to read as the uncertainty of the correlation; for the grouped Nusselt correlation, which is much better determined, the ratio is only 2.6 and the instrumental term is not negligible.

**Systematic bias in the outlet gas probe** is propagated as a declared band, because it is not a calibration error and is much larger than one. T3 sits in the exit plenum at a convective–radiative equilibrium between the gas and the surrounding hardware and reads neither bulk gas nor wall temperature; the bias is of order tens of kelvin and is not characterized, a difficulty identified as a principal source of the model–experiment gap in this field [10] and demonstrated for thermocouples in monolith channels specifically [19]. Since ε, NTU, h_app, Nu, ε*, η and — through ε — C_eff and K_loss all derive from T3, each is reported with a δT3 = ±25 K band alongside its Monte Carlo interval. The band is carried through the identification as well as through the steady groups, by recomputing the abscissa at each offset with the same estimator the archived value uses, so that the zero-offset case reproduces the archived point estimate exactly:

| quantity | δT3 = −25 K | δT3 = 0 | δT3 = +25 K | band |
|---|---:|---:|---:|---|
| $C_{\rm eff}$, cooling matched-ε (primary) [J K⁻¹] | 271.7 | 275.9 | 280.5 | 271.7 – 280.5 |
| $K_{\rm loss}$, cooling matched-ε (primary) [W K⁻¹] | 0.0844 | 0.0798 | 0.0754 | 0.0754 – 0.0844 |
| $C_{\rm eff}$, heating deep probes [J K⁻¹] | 262.2 | 280.7 | 299.4 | 262.2 – 299.4 |
| $K_{\rm loss}$, heating deep probes [W K⁻¹] | 0.1073 | 0.1137 | 0.1200 | 0.1073 – 0.1200 |

The identified capacitance is thus the quantity least sensitive to the outlet-probe systematic — the primary determination moves by 1.6% over the whole band, against the 11% spread among estimators and the ±37 J K⁻¹ Monte Carlo standard deviation — while the loss bracket widens from 0.080–0.114 to 0.075–0.120 W K⁻¹. The structural result of section 5.2 is deliberately constructed to survive the band. Source: `results.json`, key `T3_sensitivity`.

**Estimator choice** is propagated for the inversion threshold, where it dominates: section 4.1 reports ε* located both by interpolation between the bracketing runs and by a global linear regression, and carries the difference as a systematic.

---

## S3 Conditionality of the heating identification

The heating approach to steady state decays with the same eigenvalue in principle, and log(T_ss − T) is regressed over a deficit window in u = (T_ss − T)/(T_ss − T₀). Two choices are unresolvable within this dataset.

Sensor selection. With the lamps on, the front probes are held by a local balance between absorbed flux and aperture reradiation with its own short time constant, so their approach to steady state is bi-exponential and does not carry the assembly-scale mode. Using all six receiver probes gives C_eff = 121.6 J K⁻¹; using the three deep and outlet probes gives 280.7 J K⁻¹, a ratio of 2.31.

Fit window. Sweeping the deficit window moves the identified capacitance monotonically while the coefficient of determination improves, so goodness of fit does not select the window.

| deficit window $u$ | $C_{\rm eff}$ [J K$^{-1}$] | $K_{\rm loss}$ [W K$^{-1}$] | $r^2$ |
|---|---:|---:|---:|
| (0.05, 0.35) | 294.9 | 0.1274 | 0.8405 |
| (0.07, 0.45) | 280.7 | 0.1137 | 0.9033 |
| (0.10, 0.50) | 263.9 | 0.0999 | 0.9358 |
| (0.15, 0.60) | 228.1 | 0.0740 | 0.9621 |
| (0.20, 0.70) | 199.6 | 0.0564 | 0.9666 |

Sensor selection: all six receiver probes give 121.6 J K$^{-1}$, the three deep and outlet probes 280.7 J K$^{-1}$, a ratio of 2.31.

Source: `results.json`, keys `sensor_selection_swing` and `heating_window_sensitivity`.

---

## S4 Wall-profile reconstruction and the fixed-conductance falsification

The length-averaged wall temperature and the profile-corrected transfer-unit count both use a piecewise-linear axial wall profile through the probes at 11, 58 and 107 mm over a 137 mm receiver. The front 11 mm and rear 30 mm are unmeasured. Four admissible extrapolation combinations give exponents within 0.013 of one another; a full linear continuation of the rear gradient is inadmissible, because for five of the fifteen runs it places the exit wall below the measured outlet gas temperature, which no single-stream model can produce.

| Rear / front extrapolation | exponent | s.e. | runs | infeasible |
|---|---:|---:|---:|---|
| constant / constant | +0.3407 | 0.0405 | 15 | none |
| constant / linear | +0.3429 | 0.0409 | 15 | none |
| half / constant | +0.3306 | 0.0410 | 15 | none |
| half / linear | +0.3328 | 0.0414 | 15 | none |
| linear / constant | — | — | 10 | E71, E75, E76, E80, E81 |
| linear / linear | — | — | 10 | E71, E75, E76, E80, E81 |

The primary value is the constant/constant case, +0.3407 ± 0.0405, against the single-stream requirement of -1.00. Under the common-slope grouping used for the Nusselt correlation the exponent is +0.3560 ± 0.0099 with r² = 0.9922.

The fixed-conductance falsification of section 5.2 fits a piecewise-linear non-negative h(z) by multi-start least squares in log h and asks whether one flow-independent profile reproduces the fifteen measured outlet temperatures. It does not, in any family or at any node count, and the discriminating statistic is the ordering of the residual in flow rather than its magnitude.

| Family | nodes | RMS residual [K] | max abs [K] | $r$ vs $\ln Re$ | slope [K per e-fold] |
|---|---:|---:|---:|---:|---:|
| shared dimensional $h(z)$ | 2 | 64.2 | 111.3 | -0.943 | -148 |
| shared dimensional $h(z)$ | 3 | 55.4 | 96.4 | -0.946 | -128 |
| shared dimensional $h(z)$ | 5 | 54.2 | 95.0 | -0.947 | -126 |
| shared dimensional $h(z)$ | 7 | 54.2 | 95.0 | -0.947 | -126 |
| shared $Nu(z)$ | 3 | 55.8 | 96.0 | -0.945 | -132 |
| shared $Nu(z)$ | 5 | 54.6 | 94.4 | -0.947 | -129 |
| per-configuration, 256 kW m$^{-2}$ | 5 | 36.0 | 56.9 | -0.996 | -90 |
| per-configuration, 304 kW m$^{-2}$ | 5 | 67.7 | 97.3 | -0.982 | -128 |
| per-configuration, 456 kW m$^{-2}$ | 5 | 52.1 | 79.1 | -0.997 | -197 |

Source: `results.json`, keys `wall_extrapolation`, `ntu_profile` and `fixed_profile_test`; the latter also holds the fitted nodes and per-run residuals with run identifiers.

---

## S5 Auxiliary dimensionless groups, per run

The four auxiliary groups close off the local explanations for the assembly-scale exchange deficit of section 4.2: the web Biot number, the Graetz number with its reciprocal x* = 1/Gz_L (the dimensionless axial coordinate of the thermal-entrance problem at the channel exit), the axial Péclet group, and the radiation–conduction number.

| Run | $G_0$ [kW m$^{-2}$] | $q$ [sL min$^{-1}$] | $Re_{\rm nom}$ | $Pr$ | $Gz_L$ | $x^*=1/Gz_L$ | $Re\,Pr\,L/D_h$ | $Bi$ | $N_{rc}$ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| E67 | 456 | 15.28 | 72.6 | 0.6837 | 0.543 | 1.84 | 4531 | 2.53e-05 | 4.37 |
| E68 | 456 | 12.50 | 58.9 | 0.6838 | 0.441 | 2.27 | 3680 | 1.88e-05 | 4.86 |
| E69 | 456 | 10.50 | 49.6 | 0.6838 | 0.372 | 2.69 | 3100 | 1.48e-05 | 5.09 |
| E70 | 456 | 9.11 | 42.7 | 0.6839 | 0.320 | 3.12 | 2670 | 1.19e-05 | 5.64 |
| E71 | 456 | 7.13 | 34.0 | 0.6837 | 0.254 | 3.93 | 2122 | 8.38e-06 | 5.66 |
| E72 | 304 | 18.32 | 94.0 | 0.6845 | 0.704 | 1.42 | 5876 | 3.39e-05 | 2.68 |
| E73 | 304 | 13.17 | 65.5 | 0.6835 | 0.490 | 2.04 | 4086 | 2.17e-05 | 3.41 |
| E74 | 304 | 9.03 | 44.8 | 0.6834 | 0.335 | 2.99 | 2794 | 1.25e-05 | 3.96 |
| E75 | 304 | 6.95 | 34.9 | 0.6838 | 0.261 | 3.83 | 2179 | 8.53e-06 | 4.11 |
| E76 | 304 | 4.53 | 23.3 | 0.6846 | 0.175 | 5.71 | 1460 | 4.81e-06 | 4.08 |
| E77 | 256 | 13.85 | 78.8 | 0.6889 | 0.594 | 1.68 | 4956 | 2.60e-05 | 1.50 |
| E78 | 256 | 10.01 | 55.5 | 0.6876 | 0.418 | 2.39 | 3486 | 1.68e-05 | 1.84 |
| E79 | 256 | 8.04 | 44.4 | 0.6875 | 0.334 | 2.99 | 2788 | 1.22e-05 | 2.01 |
| E80 | 256 | 6.61 | 36.5 | 0.6874 | 0.274 | 3.64 | 2289 | 9.31e-06 | 2.12 |
| E81 | 256 | 4.53 | 25.2 | 0.6878 | 0.190 | 5.27 | 1582 | 5.54e-06 | 2.23 |

Source: `groups.csv`. Properties are evaluated as stated in section 3.1 of the manuscript.

---

## S6 Delivered aperture irradiance: the two determinations in full

The aperture irradiance is determined two ways, and they disagree.

The first is radiometric flux mapping, using the three-step characterization procedure established for this simulator [22,23]: an alumina-coated Lambertian target on a precision three-axis slider, a Gardon-type circular foil radiometer traversed across the focal region on a graded grid with NIST-traceable calibration and calorimetric correction for spectral mismatch, and a parallel high-dynamic-range CCD reconstruction whose grayscale-to-flux transfer fits a power law with a = 1117.5 W m⁻² and exponent 0.903. Gaussian fits to the resulting profiles validate the focal distribution with normalized RMSE of about 0.07 and 0.12 on the two axes [22]. Integrating the calibrated map over the receiver face gives nominal values of 456, 304 and 256 kW m⁻² for the three configurations.

The second is an apparent steady energy closure, not model-free: it takes K_loss from the lumped identification of section 4.4 and omits storage by assuming the assembly is stationary at run end. Neither assumption is safe: insulation probe T2 rises monotonically as flow falls and run duration lengthens, from 308 to 356 K across the campaign, with decay far slower than any run window, so insulation is still charging when runs end, leaving a residual storage term of unknown sign. Nor need a K_loss from a cooling secant equal the steady loss under irradiation, since part of the path is radiative and insulation temperature differs between states. At steady state:

$$f = \frac{Q_{\rm gas} + K_{\rm loss}(\bar T_w - T_{\rm amb})}{G_0 A_{\rm frt}}$$

which should equal unity if the radiometric flux is correct. Evaluated at both ends of the identified K_loss bracket and averaged per configuration, it does not: f runs 0.989–1.134 at 456 kW m⁻², 1.149–1.324 at 304 and 0.784–0.916 at 256, giving closure estimates of 451–517, 349–402 and 201–234 kW m⁻², with reported spans 451–517, 304–402 and 201–256 kW m⁻².

We report both estimates, with every power-normalized quantity spanning them. This is not a bounding interval: closure omits storage, inherits a lumped loss model, and carries T3 bias directly into Q_gas, so neither determination is an established bound and the disagreement cannot be assigned to one cause. Apparent efficiency on the nominal basis ranges 0.252–1.224, exceeding unity in the 304 kW m⁻² group.

One contributor is geometric. The Gardon gauge reports flux over its 12.6 mm-radius exposed area [22], 499 mm² against the receiver's 361 mm² face, a ratio of 1.38. In a focused Gaussian spot the gauge averages a larger, more peripheral region than the receiver intercepts, under-reading the face flux: the face-to-gauge ratio for a centred Gaussian is 1.28, 1.17, 1.11 and 1.04 for 1/e² radii of 10, 14, 18 and 30 mm, bracketing the 1.13 required at 456 kW m⁻² and covering roughly half the 1.32 required at 304. This is structural to characterizing a 19 mm receiver with a 25 mm gauge, not execution error, but it does not explain 256 kW m⁻², where closure sits 8–16% below nominal and a broader spot pushes the factor toward unity, so aperture spillage dominates oppositely. Comparably, optical modelling of a directly irradiated reactor on the same simulator found only 1–6% of aperture-incident power reaching the catalytic bed [23]. Supplementary section S6 gives the gauge-averaging calculation and remaining instrument contributions.

Outlet-probe bias alone cannot explain the discrepancy: the uniform T3 shift bringing closure to f = 1 is −75 K at 456 kW m⁻², −126 K at 304 and +75 K at 256 — sign-varying and varying 250 K within a configuration across flow, so no uniform correction reconciles the two determinations. A flow-dependent energy-budget term is required, without identifying which: bypass around the monolith periphery and flow maldistribution among channel groups are candidates consistent with section 4.2, but non-stationarity or flow-dependent probe response are not excluded. Convergence with section 5.2 is a consistency, not independent confirmation.

No temperature-based result depends on this choice: effectiveness, transfer units, the Nusselt correlation, the inversion criterion, disequilibrium indices, eigenvalue identification and similarity collapse are all ratios of measured temperatures and flows. Only efficiency and delivered-power-per-unit-mass figures shift, reported as spans rather than intervals containing the true value.

Source: `results.json`, key `delivered_power`; `flux_geometry.json` for the gauge-averaging calculation.

---

## S7 What a flow-dependent split requires: three destabilising agents

That prescription has been carried out elsewhere, which makes the recommendation concrete rather than aspirational. Nagarajan, Ström and Sjöblom build a reduced-order pseudochannel model of a catalytic converter precisely because single-channel models cannot predict accurately under flow maldistribution in realistic geometries [35] — the remedy being a few representative channels carrying different flows rather than one channel with a better closure. The same group whose channel-scale Nusselt work anchors section 4.2 has developed flow-distribution models for dual-cell-density monoliths [36] and dual-zone packed beds [43]. Within the solar field the step has been taken at panel and module scale: Schwager and co-workers resolve the mass-flow distribution among the parallel tubes of a molten-salt receiver panel under flux gradients [37], and the multi-zone unsteady model of Kumar and co-workers resolves exchange between absorbers, casing and return air rather than treating the receiver as one stream [34].

Those two flow-distribution studies do more than set a precedent: they identify what a flow-dependent split requires. Reinao and Cornejo model a dual-cell-density monolith as two concentric regions of differing apparent permeability with a Darcy momentum sink, isothermal, and find the flow fraction through each region depends only on relative permeability and core size, not on flow rate, across several flow rates [36]. That invariance is exact: when every parallel path has a resistance linear in velocity and shares one viscosity, flow rate cancels from the equal-pressure-drop condition. Díaz and co-workers repeat the exercise for a dual-zone packed bed, where the Ergun equation adds an inertial term, and there the split does depend on Reynolds number — strongly sensitive to the resistance contrast at low Re where the viscous term dominates, converging as inertia takes over [43]. A split that changes with flow rate therefore cannot arise in parallel linear-viscous paths of common viscosity; it requires an inertial contribution or a path-to-path viscosity difference. This receiver's participation does change with flow, as Re_nom^+0.341, and both admissible causes are present: the heated absorber carries a non-uniform viscosity field [9], and the viscous-to-inertial transition is the same axis as the linear-versus-quadratic stability criterion of section 4.7 [27,28]. The direction matches too, the packed-bed split growing more uniform as flow rises.

A third agent gives the same structure. In Schwager's molten-salt panels the destabilizing force is buoyancy rather than inertia or viscosity — hotter, less dense salt gains flow in an upward-flowing panel but loses it in a downward-flowing one, where the effect is self-reinforcing and can drive tubes to zero flow — and what arrests it is friction, so reducing overall flow intensifies the maldistribution: below 20% of nominal flow mild flux spreads suffice to trigger reversal, and at 10% no spread is needed [37]. That mechanism is not this receiver's, but in all three cases a driving force for maldistribution competes with a frictional resistance whose relative weight grows with flow, so the split is most differentiated at low flow and most uniform at high flow. The flow-rate dependence of participation is therefore a general property of heated parallel-path arrays rather than an artefact of one geometry.

---

## S8 Pressure drop as installed

Steady-window means of the differential pressure channel compare poorly with the laminar square-duct prediction, computed with a Poiseuille number of 56.91, all metered flow through 100 channels, and properties at T̄_g. The measured signal exceeds the monolith-only prediction by 1.4 to 2.8 times at high flow (2.74 against 0.99 mbar at 15.3 sL min⁻¹), falls below it at low flow, and goes negative at the lowest flow (−0.02 against 0.20 mbar at 4.5 sL min⁻¹) — the signature of a zero offset on a transducer whose ±0.2 mbar accuracy floor, 0.1% of its ±200 mbar span, is the same order as the entire expected monolith drop of 0.20 to 1.0 mbar. The taps evidently span more than the monolith faces, and the second differential channel reads a flat −0.26 to −0.30 mbar in every run and appears unconnected.

This matters because pressure drop is how the volumetric-receiver literature diagnoses honeycomb flow behaviour: absorbers with linear velocity dependence tend inherently towards flow instability [27]; high-thermal-conductivity materials with a quadratic characteristic are stable [28]; instability becomes possible when one pressure-drop level maps to more than one outlet temperature, letting a region drift and fail if temperature exceeds material limits [26]. Ávila-Marín's review extends this to honeycombs, identifying a critical flux above which highly porous honeycombs become unstable while wire meshes and low-porosity foams stay stable via their quadratic characteristic [1]; permeability and Forchheimer coefficients from pressure-drop data feed homogeneous equivalent models [13]. Recent unheated monolith measurements show a linear Re–pressure-drop relationship throughout, with progressive underestimation as Re rises through entrance/exit effects [45], large enough at channel Re 50–320 to need a spacing-dependent correction factor increased by misalignment [44]. Our isothermal honeycomb sits on the instability-prone, linear side, but our transducer, spanning monolith plus plumbing, cannot separate the absorber's characteristic from comparable entrance, exit and fitting losses. We report this because the fix is inexpensive, and section 5.3 shows this measurement would discriminate between the two candidate mechanisms above.

Source: `pressure_drop.csv`, which carries the per-run measured differential pressure on both channels alongside the laminar square-duct prediction.

---

## Reproduction

Two commands regenerate every table, scalar and figure in the manuscript and in this document from the raw logger files: `python receiver_reduction.py` writes the archived outputs, and `python make_figures.py` draws Figures 2 to 6 from them. Both are seeded and deterministic. Figure 1 is a schematic of the apparatus and sensor positions, drawn rather than computed, and is the only figure not produced by the pipeline. The accompanying `README_reproduction.md` documents the runtime, the seeds, what each archived file supports, and the nominal and apparent conventions a reader needs before using any number.

## References

References are numbered as in the manuscript.
