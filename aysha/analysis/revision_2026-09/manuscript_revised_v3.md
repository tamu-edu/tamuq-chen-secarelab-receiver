# Transient Thermal Characterization and Dimensionless Correlations for a Structured Silicon Carbide Volumetric Solar Receiver, with Consequences for Reduced-Order Model Validation

A. Melhim, A. G. Konstandopoulos, K. J. Kakosimos

---





## Abstract

A monolithic silicon carbide honeycomb receiver — 100 square channels of 1.5 mm side, 0.40 mm webs, 137 mm long, 0.623 open frontal fraction — was characterized in a high-flux solar simulator at three lamp configurations and air flows of 4.5 to 18.3 sL min⁻¹, giving nominal channel Reynolds numbers of 23 to 94. Fifteen heating and three cooling transients were reduced to dimensionless form without recourse to a thermal model, and the receiver's constants identified from the transient eigenstructure rather than by inverse simulation. Gas–solid exchange is limited at the assembly scale rather than in the channel: the apparent Nusselt number, referenced to the length-averaged wall temperature, scales as Re_nom^1.470 ± 0.011, one to two orders of magnitude below the fully developed laminar square-duct values, and the slow thermal eigenvalue identifies an effective capacitance of 276 ± 37 J K⁻¹, six times that of the bare monolith.

Two results contradict standard practice in receiver modelling. First, the transfer-unit count increases with flow: obtained by integrating the measured axial wall profile, it scales as Re_nom^+0.341 ± 0.041, where a single-stream formulation whose volumetric conductance is any fixed function of axial position requires Re_nom^−1.00. Flow-independent axial conductance profiles fitted instead — shared across runs, shared as Nusselt profiles, and fitted within each irradiance configuration, with up to seven free non-negative nodes — leave outlet residuals of 36 to 68 K root-mean-square that are monotone in flow, so no tested family reproduces the measured response. Second, the same fifteen runs yield a correlation prefactor spanning a factor of 29, an inversion threshold spanning 0.57 to 0.79 and a capacitance spanning 122 to 281 J K⁻¹, according only to which solid thermocouple is nominated as the reference and which probe subset and fit window are used. A model validated against a gas outlet temperature and one or two solid temperatures is therefore matched to a quantity that is not a property of the receiver, which offers a concrete explanation for the field's long-standing observation that measured absorber performance falls short of prediction while calibrated models appear to fit. We identify the instrumentation that would convert these installed-receiver response indices into physical coefficients.

**Keywords:** volumetric solar receiver; silicon carbide honeycomb; local thermal non-equilibrium; transient identification; model validation; dimensionless correlation

---





## 1. Introduction

Volumetric receivers distribute concentrated radiation into the depth of a porous or structured absorber so the solid temperature peak sits inside the material rather than on the irradiated face, suppressing reradiation from the hot front surface and raising attainable outlet temperature. Reviews anticipate thermal efficiencies above 75% with air outlet temperatures beyond 1000 °C, silicon carbide absorbers usable to 1500 °C where metallic meshes are confined to 800–1000 °C [1]. Three decades of prototypes have delivered unevenly: the wire-mesh TSA absorber reached 85% efficiency at 700 °C outlet and the recrystallized SiC SOLAIR absorber 70–75% at 750 °C [3,6], but Sulzer 1 delivered 68% against 80% predicted by its own design code, and Catrec 2 was destroyed attempting to raise its 460 °C outlet temperature [1,6]. Reviews, most recently in 2023 [46], state that the theoretical volumetric effect has not been validated by experiment, and most tested absorbers perform worse than predicted [3].

The response has largely been to build better absorbers: foams and fibre meshes characterized for thermophysical and heat transfer properties [4,5], twenty-six graded-porosity wire meshes and five SiC morphologies including Voronoi structures [6,7], single- against multi-layer ceramic foams [8], and pore-diameter variation at fixed porosity [9]. Current work continues the pattern — a 5 kW reticulated-porous-ceramic cavity receiver compared across SiSiC, alumina and ceria [31]; screen-printed ceramic and metallic-foil three-dimensionally shaped absorbers in a sixteen-module matrix at DLR's Synlight facility, reaching 735 °C exit and 91.8% efficiency at 650 °C against 84.6–85.3% for the 80 cpsi HiTRec SiSiC benchmark [29]; and a calorimetric bench measuring heat transfer to pressurised gas across rectangular channel geometries [32]. The recurring output is a design recommendation — porosity, cell density, layering, channel aspect ratio — from steady-state efficiency measurements. Supplementary section S1 gives campaign-by-campaign detail.

In parallel, a modelling literature treats the absorber as a porous continuum with local thermal non-equilibrium and a radiation model. An idealized absorber in local thermal equilibrium predicts above 95% efficiency at 800 kW m⁻² and 1000 °C, unrealistically optimistic against the roughly 70% real absorbers achieve [2]; an extruded SiC honeycomb was matched to within 1–2% in outlet temperature and efficiency once its extinction coefficient of 2300 m⁻¹ was measured independently [11]; and homogeneous equivalent models with a P1 radiation approximation close the gas–solid coupling with volumetric correlations of the form Nu_lv = a Re^b Pr^c [13]. Superposing separately computed conduction and radiation contributions fails to recover the coupled effective conductivity [14], and a Rosseland treatment can differ from discrete ordinates by up to 76% [15]. Our own group showed specular rather than diffuse channel reflectivity flattens the source term and raises modelled efficiency to 96% against 85% for the commercial HiTRec-II geometry [16]. Local thermal non-equilibrium is throughout a premise of the equations rather than a measured quantity, and the volumetric coupling coefficient is taken from a correlation obtained elsewhere or fitted.

The most rigorous heat transfer analysis of honeycomb channels comes from outside the solar field. Cornejo, Nikrityuk and Hayes computed laminar convection in monolith channels under constant wall temperature and constant wall heat flux across Re of 50–600, showing that forming the Péclet number on an inlet rather than local temperature introduces errors of 50% or more, and that at high heating rates a minimum appears in the Nusselt–Graetz curve some 16% below the asymptotic value [17]. Effective axial and radial conductivities are established for honeycomb substrates [18], and the multi-scale review notes radiation is required to reproduce measured axial temperature profiles and that a thermocouple inside a monolith channel can be significantly affected by radiation even without blocking the channel [19]. A 0.05 mm wall gap substantially reduces effective radial conductivity [20], and a pseudo-continuous two-dimensional reduction reproduces full three-dimensional CFD to better than 0.5 K [21]. These results assume a heated wall and uniform inflow — flow equally distributed over all channels [21], gas treated as static [20] — conditions an installed receiver under a concentrated, spilling beam does not satisfy.

Neither family measures an installed receiver through its transients and reads the governing quantities directly off the measurement, which matters because calibration data for a reduced-order model are almost always a gas outlet temperature with one or two embedded solid temperatures. Capuano and co-workers illustrate the consequence: reviewing four established DLR models against one reference experiment, they report predicted outlet air temperatures of 1050, 1048, 1057 and 1073 K against 972 K measured, and efficiencies of 0.78, 0.77, 0.79 and 0.81 against 0.69, attributing the gap partly to the difficulty of experimentally evaluating a mean air flow temperature [10]. Four independently developed models agreeing with each other to within 25 K yet disagreeing with the experiment by 80–100 K signals a validation target that fails to constrain what the models disagree about.

This paper instead reduces the measurements to dimensionless groups, identifies the receiver's constants from the transient eigenstructure, and tests the structure of reduced-order formulations against a measured exponent. This structural test does not inherit the identifiability problems afflicting inverse calibration, since it constrains the sign of a derivative rather than a coefficient's value; it also quantifies how much identified coefficients depend on the analyst's choice of reference sensor, and thus how much information a conventional validation exercise carries. The contributions are a dimensionless characterization of an installed SiC honeycomb receiver including a directly instrumented rather than assumed wall-to-interior disequilibrium; identification of effective capacitance and loss conductance from transient eigenvalues with sensor-selection sensitivity quantified; a similarity rescaling on which fifteen heating transients collapse; and two methodological results — that validation against a gas outlet temperature plus one or two solid temperatures is largely uninformative for this class of receiver, and that a uniform-temperature or uniform-flux coefficient evaluated on the fixed channel geometry cannot reproduce its flow response at any calibration.

## 2. Experiment





### 2.1 Receiver and apparatus

The absorber is a monolithic silicon carbide honeycomb with a ten-by-ten array of square channels of 1.5 mm side and 0.40 mm web thickness, giving a 1.9 mm pitch, a 19.0 mm square frontal face of area 3.61 cm², an open frontal fraction of 0.6233 and an axial length of 137 mm. The measured mass is 40 g, which with the solid cross-section of 1.360×10⁻⁴ m² corresponds to an effective bulk density of 2147 kg m⁻³, consistent with a recrystallized SiC of roughly one third internal porosity; this is the value used for every capacitance estimate below. The geometry is of the same family as the extruded SiC honeycombs and catalyst-carrier substrates characterized in the solar literature [4,11] and modelled in our earlier work [16], with a channel size between the two samples of Fend and co-workers (2.18 and 1.43 mm) and a comparable web thickness.

The monolith is held in an insulated cavity and irradiated through a circular aperture by the seven-lamp high-flux solar simulator at Texas A&M University at Qatar, each lamp a 6 kW xenon short-arc unit; the apparatus and its optical characterization are described in detail elsewhere [22,23]. Air is metered by four Aalborg GFC17 thermal mass-flow controllers, all carrying air, referenced to 21.1 °C and 1 atm, and the total flow is their sum. Their manufacturer specification is ±1.0% of full scale over the whole 0 to 100% range with ±0.5% full-scale repeatability, and the factory calibration is gas-specific; because the receiver runs air rather than the calibration gas, each unit was calibrated in place against a bubble flowmeter, which replaces the gas-conversion factor with a directly measured characteristic and leaves a residual proportional to reading. Both terms are carried in the uncertainty model of section 3.4. Per-unit readings span 0.69 to 5.72 sL min⁻¹, that is 7 to 57% of full scale, so neither term is negligible against the other. Because these are mass-flow devices reporting standard volumetric units, the mass flow follows from the reference density ρ_std = 1.1996 kg m⁻³ rather than from the density at the prevailing ambient temperature.





### 2.2 Instrumentation

Nine thermocouples define the measurement. Three sheathed probes contact the monolith side wall at axial positions of 11, 58 and 107 mm, denoted T8, T12 and T11. Two probes sit in the interior, flow-exposed region at 58 and 107 mm, denoted T9 and T10. One probe, T3, sits in the exit plenum and reports the outlet gas temperature. One further probe, T2, sits in the surrounding insulation and is used in section 4.6 to assess thermal stationarity; the wall-to-insulation drop T12 − T2 is 310 to 770 K, so the insulation resistance dominates radial transfer. The remaining two probes give the cavity ambient, whose mean is T_amb. The distinction between the wall chain and the interior pair is not cosmetic: sections 4.3 and 5.1 show that the two report materially different temperatures at the same depth, and that which is taken to be the solid temperature changes every identified coefficient in this paper. That the interior probe reading is itself perturbed by the radiation field is established for this geometry [19], which is why we treat the reference-temperature choice as a systematic rather than a detail.

Differential pressure across the monolith is monitored with a Keller PD33X of ±200 mbar range and ±0.1% full-scale accuracy, hence an accuracy floor of ±0.2 mbar. Section 4.7 shows this to be the same order as the entire expected monolith pressure drop, so the channel is unusable as installed; we report it because the corrective action that follows is concrete, and because pressure drop is the quantity that discriminates between the two mechanisms left open by section 5.2.

Thermocouple tolerance is treated per IEC 60584 for type N as the greater of 1.5 K and 0.004|t|, which at the front probe near 1200 K is ±4.4 K. Installation error for the three sheathed wall probes, comprising contact resistance against the monolith and radiative exchange with the cavity, is not quantified, and we declare it as such rather than absorbing it into the stated budget.





### 2.3 Campaign and reduction

Fifteen heating runs span three lamp configurations at nominal aperture irradiances of 456, 304 and 256 kW m⁻²: runs E67–E71 at 7.1 to 15.3 sL min⁻¹, E72–E76 at 4.5 to 18.3 sL min⁻¹, and E77–E81 at 4.5 to 13.8 sL min⁻¹ respectively. Three cooling transients, C69, C80 and C81 at 10.5, 6.6 and 4.5 sL min⁻¹, follow lamp shutdown with flow maintained. Two additional runs, E82 and E83, repeat E77 and E81; they are retained in the archive for transparency and excluded from every fit reported here.

Steady-state values are 120 s means at the end of each heating run. The length-averaged wall temperature is the length average of the piecewise-linear profile through the three wall probes with constant extrapolation to the end faces,

$$\bar T_w = 0.2518\,T_8 + 0.3504\,T_{12} + 0.3978\,T_{11},$$

these being the exact trapezoid coefficients for probes at 11, 58 and 107 mm over 137 mm. Gas enthalpy differences are computed by trapezoidal integration of c_p(T) between T_amb and the relevant gas temperature rather than as c_p(T̄)ΔT.

---

## 3. Analysis framework





### 3.1 Dimensionless groups

Per-channel mass flow is ṁ_ch = ṁ/100 and the hydraulic diameter is the channel side, D_h = 1.5 mm. Gas properties are evaluated at T̄_g = (T_amb + T_3)/2 except the thermal conductivity in the Nusselt number, which uses the film temperature (T̄_w + T̄_g)/2 — not innocuous here, since forming the Péclet number on an inlet rather than a local temperature shifts the apparent Nusselt–Graetz curve by 50% or more [17], so it is stated explicitly and held fixed across every variant in section 5.1.

The receiver is characterized as a single-stream heat exchanger following Kays and London [25], with effectiveness and transfer-unit count

$$\varepsilon = \frac{T_3 - T_{\rm amb}}{\bar T_w - T_{\rm amb}}, \qquad NTU = -\ln(1-\varepsilon),$$

giving an apparent assembly-level coefficient and Nusselt number:

$$h_{\rm app} = \frac{NTU\,\dot m_{\rm ch} c_p}{P L}, \qquad Nu = \frac{h_{\rm app} D_h}{k(T_{\rm film})}.$$

Re is formed on ṁ/100, assuming equal division of metered flow among the hundred channels — an assumption sections 5.2 and 5.3 argue is violated — so it is a nominal Reynolds number, Re_nom, where the distinction matters. Similarly, NTU_app, h_app and Nu_app are built algebraically from the outlet-plenum probe T3 and an axial wall average, not independently measured interfacial coefficients. All four are reproducible installed-receiver response indices, appropriate for the structural arguments below but not validated transport coefficients and not transferable to another geometry; extracting physical coefficients is left to a companion modelling study.

Two qualifications matter. First, the ε–NTU relation is exact only for a wall at fixed temperature, whereas here the wall spans up to 321 K between T8 and T11 and is non-monotonic in the inverted regime, making the inversion approximate, with bias quantified in section 4.2. Second, since h_app is defined through NTU, Nu(Re) is algebraically tied to ε(Re):

$$\frac{d\ln Nu}{d\ln Re} = 1 + \frac{d\ln NTU}{d\ln Re} + \frac{d\ln[c_p/k]}{d\ln Re},$$

so a fitted Nu(Re) exponent is not independent corroboration of the measured effectiveness range, enabling the structural test of section 5.2: comparing that exponent against channel-correlation values near 0.3–0.6 compares unlike objects, since the relevant null here is −1 for the transfer-unit count, not 0.5 for Nu.

The apparent wall-to-interior disequilibrium — a bounded surrogate for gas–solid nonequilibrium (section 4.3) — is formed at each instrumented depth:

$$\Lambda_{58} = \frac{T_{12}-T_9}{T_{12}-T_{\rm amb}}, \qquad \Lambda_{107} = \frac{T_{11}-T_{10}}{T_{11}-T_{\rm amb}},$$

with volumetric inversion indicator I = T₁₂ − T₈, negative while the front-face probe is hotter, positive once mid-depth overtakes it. Four auxiliary groups, used in section 5.2, bound the candidate physics: web Biot number Bi = h_app t_w/2k_SiC, Graetz number Gz_L = D_h Re Pr/L, axial Péclet group Re·Pr·L/D_h, and radiation–conduction number N_rc = 4σT̄_w³D_h/k(T̄_w).



### 3.2 Transient identification

After lamp shutdown the assembly decays as a single slow mode: for each sensor the late-window excess temperature is regressed as log(T − T_amb) against time, and λ is the mean of the per-sensor slopes. Lumping the assembly into one capacitance C_eff losing heat to ambient through a conductance K_loss while delivering enthalpy to the gas gives

$$\lambda = \frac{\varepsilon \dot m c_p + K_{\rm loss}}{C_{\rm eff}},$$

so regressing λ against the gas advective conductance x = ε ṁ c_p identifies both constants from slope and intercept with no simulation in the loop. For each cooling transient, ε is fixed by the experiment log rather than by flow proximity, since two heating runs differ by 0.001 sL min⁻¹. The approach is a transient identification in the spirit of the temperature-wave method used to obtain volumetric heat transfer coefficients of porous absorbers from phase shift and attenuation [4], here applied to the installed assembly rather than a material sample.

The heating approach to steady state decays with the same eigenvalue in principle, with log(T_ss − T) regressed over the deficit window 0.07 < u < 0.45, where u = (T_ss − T)/(T_ss − T₀). That identification is conditional in two ways, on sensor selection and on the fit window: all six receiver probes give C_eff = 122 J K⁻¹ where the three deep and outlet probes give 281 J K⁻¹, and sweeping the deficit window from (0.05, 0.35) to (0.20, 0.70) moves C_eff monotonically from 295 to 200 J K⁻¹ while r² improves from 0.84 to 0.97, so fit quality cannot select the window (supplementary section S3). We therefore report cooling as the primary identification, needing no sensor rule once the source is removed; heating as a conditional consistency check; and a joint fit to all eighteen eigenvalues as the pooled estimate. Their disagreement is addressed in section 5.1.


### 3.3 Similarity rescaling

If the lumped description holds, every transient should collapse when time is rescaled on the identified eigenvalue as t* = t/τ with τ = C_eff/(ε ṁ c_p + K_loss). The ε factor belongs in τ because it is the effectiveness that converts the metered flow into the conductance the eigenvalue actually sees; both the reference constants and the effectiveness are taken from the primary identification and from each run's own measured effectiveness. The half-rise value that results is a property of the receiver's response shape and carries no expected value: a single-node lumped system would give ln 2, and the departure from it is itself the measure of how far the receiver is from one state.





### 3.4 Uncertainty

Three kinds of uncertainty are propagated and reported separately, since conflating them misleads.

Random and calibration error is propagated by Monte Carlo over 4000 realizations. The IEC 60584 class 1 type N tolerance — the greater of 1.5 K and 0.004|t| — is a maximum permissible error rather than a standard deviation, so it is read as a 95% bound and halved to a 1σ systematic term, drawn once per sensor per realization and shared across runs, with an independent 0.5 K drift term per run; viscosity and conductivity carry 2% and 1% perturbations. The metered flow sums four controllers carrying roughly 25, 33, 28 and 14% of the total, and their error is two-term because the two available pieces of information constrain different things: the manufacturer specification is ±1.0% of full scale with ±0.5% full-scale repeatability, an additive bound in absolute flow rather than a percent of reading, while the in-house bubble-flowmeter calibration constrains the slope of each unit's characteristic and so contributes a term proportional to reading. Each unit therefore carries σ = [(0.0025 F)² + (0.025 q_u)²]^{1/2} with F = 10 sL min⁻¹ the per-unit full scale, and one draw per unit per realization is combined with each record's own measured shares. The two terms have opposite flow dependence, so the relative uncertainty on the summed flow falls from 3.4% at the lowest campaign flow to 2.6% at the highest for fully correlated units: the flow error is predominantly systematic, displacing a correlation's prefactor, with a weak residual gradient that does reach an exponent. Both correlation limits are computed since none is specified; they agree to within 0.0007 on every exponent and 0.5 J K⁻¹ on every capacitance, so the larger is quoted and both archived. Fitted-slope quantities report the regression standard error alongside, labelled as such: for the pooled correlation and the transfer-unit count of section 4.2 it is larger by an order of magnitude and is the one to read as the uncertainty, while for the much better determined grouped correlation it exceeds the instrumental term only 2.6-fold. Supplementary section S2 details the propagation.

Systematic bias in the outlet gas probe is propagated as a declared band, much larger than calibration error. T3, in the exit plenum, sits at a convective–radiative equilibrium between gas and hardware and reads neither bulk gas nor wall temperature; the bias is of order tens of kelvin and uncharacterized — a principal source of the model–experiment gap in this field [10], demonstrated for thermocouples in monolith channels specifically [19]. Since ε, NTU, h_app, Nu, ε*, η and, through ε, C_eff and K_loss all derive from T3, each carries a δT3 = ±25 K band alongside its Monte Carlo interval; the band is propagated through the identification as well as through the steady groups, giving 271.7 to 280.5 J K⁻¹ on the primary capacitance and 0.0754 to 0.0844 W K⁻¹ on its loss conductance, so the identified capacitance is the quantity least sensitive to the outlet-probe systematic. The structural result of section 5.2 is deliberately constructed to survive the band.

Estimator choice is propagated for the inversion threshold, where it dominates: section 4.1 reports ε* located both by interpolation between bracketing runs and by global linear regression, carrying the difference as a systematic.


## 4. Results





### 4.1 Steady field, flow response, and the inversion criterion

The steady field against flow at each configuration is shown in Figure 2, and the response is strongly graded in depth. Front-face wall temperature falls steeply with flow, at −33.2, −24.3 and −20.8 K per sL min⁻¹ at 456, 304 and 256 kW m⁻² with r² of at least 0.97; the mid-depth wall falls at half to two-thirds that rate; the rear wall barely falls, at −1.1, −1.9 and −5.2 K per sL min⁻¹ with r² of 0.07, 0.12 and 0.76. Outlet gas temperature is nearly flow-invariant at the two higher fluxes (+0.6 and −0.1 K per sL min⁻¹, r² ≤ 0.03) and falls at the lowest, at −3.1 K per sL min⁻¹. At 456 kW m⁻² the outlet gas spans 26 K across the flow range against 264 K for the front-face wall, a ratio of ten; at 304 and 256 kW m⁻² the spans are 51 K against 331 K and 36 K against 194 K.

The side-wall profile inverts as flow increases, crossing zero within the tested range at all three configurations (Figure 4a), and the crossings collapse when expressed in effectiveness rather than flow (Figure 4b). Two caveats apply: the indicator, a difference between two side-wall probes at 11 and 58 mm, detects a front-to-mid crossing along the instrumented wall chain rather than the three-dimensional solid temperature maximum, which under a Gaussian beam with radial gradients need not lie on the wall — the radial traverse of section 5.3 would settle that. The threshold is also convention-dependent: ε* is formed on the length-averaged wall temperature, and the sensitivity table of section 5.1 shows it moving between 0.57 and 0.79 as the reference probe changes. What follows is an operational criterion under the stated convention, not a receiver property.

| Nominal $G_0$ [kW m⁻²] | $q^*$ local [sL min⁻¹] | $\varepsilon^*$ local | $q^*$ global | $\varepsilon^*$ global |
|---:|---:|---:|---:|---:|
| 456 | 10.12 | 0.666 ± 0.002 | 11.08 | 0.673 |
| 304 | 8.36 | 0.655 ± 0.002 | 10.32 | 0.673 |
| 256 | 5.08 | 0.642 ± 0.003 | 3.69 | 0.628 |

The two estimators disagree materially. At 256 kW m⁻² the global regression places the crossing at 3.69 sL min⁻¹, 18% below the lowest run at that flux; only one point in that group is negative and the indicator is strongly concave, so a global fit under-slopes near the bottom and pushes the crossing outside the data. Local bracketing keeps all three crossings inside the measured range and halves the ε* spread from 0.045 to 0.024. The defensible statement is that inversion occurs at ε* ≈ 0.65, with a weak monotonic increase from 0.642 to 0.666 over a 1.8-fold flux range; neither flux independence nor a resolved flux dependence is claimed, since the 0.024 spread exceeds the 0.002–0.003 Monte Carlo intervals yet is comparable to the estimator systematic. The criterion is robust in form: inversion is set by how much of the available wall-to-ambient temperature difference the gas has recovered, not by flow rate or flux separately. Under the δT3 = ±25 K band, ε* at 256 kW m⁻² moves from 0.579 to 0.704, so the threshold is conditional on outlet-probe calibration.

This is the experimental counterpart of a result the modelling literature has approached from the other side. Our single-channel study found a receiver could be moved into genuinely volumetric behaviour by making the surface reflectivity specular, flattening the source term without changing porosity or cell density [16], and the field has generally sought the volumetric effect through material and geometric design [6,7,8,9]. On one plain, ungraded SiC honeycomb the side-wall profile inverts simply by raising flow past a threshold expressible in a single dimensionless number: the side-wall signature of volumetric behaviour is at least as much an operating condition as a design property — a statement about the signature, since the measured quantity is a two-probe wall difference, not the interior maximum's location.

### 4.2 Assembly-scale exchange, and a structural constraint

The apparent Nusselt number follows a clean power law over the full campaign (Figure 3a),

$$Nu_{\rm app} = (3.10 \pm 0.12)\times10^{-4}\, Re_{\rm nom}^{1.443}, \qquad r^2 = 0.971,\ n = 15,$$

with pooled regression standard error ±0.069 and instrumental Monte Carlo term ±0.004; the regression term exceeds the instrumental one roughly seventeenfold, so residual dispersion is physical and structured, not metrological or random. Fitted separately per irradiance configuration the exponent is 1.440, 1.473 and 1.522 at 256, 304 and 456 kW m⁻² (r² ≥ 0.9996 each), so pooled scatter reflects a flux-ordered offset plus a weak monotonic trend in the exponent. A common-slope model with a separate prefactor per configuration gives exponent 1.470, regression SE ±0.011, instrumental term ±0.004, r² = 0.9994, and prefactors 3.20, 2.70 and 2.55 ×10⁻⁴ in the same order; we report this grouped model as the primary correlation alongside the pooled fit, being better determined and the honest representation of data acquired at three discrete irradiance levels.

Absolute values run from 0.028 to 0.213. The fully developed laminar square-duct reference depends on thermal boundary condition; Shah and London tabulate all three for unity aspect ratio: Nu_T = 2.976, Nu_H2 = 3.091 and Nu_H1 = 3.608, with fRe = 14.227 throughout [24]. We adopt Nu_H2 as primary, since interior channels see similar radiative and convective conditions on all four walls and thin webs make a peripherally uniform flux the better idealization than a uniform temperature. Measured Nu_app is a factor of 14.5 to 111 below Nu_H2, moving only to 14.0–107 or 17.0–130 under the other two limits — a one-to-two-order-of-magnitude deficit regardless of choice, comparably far below channel-level Nusselt numbers of order three to four computed for monolith channels under either boundary condition [17]. Exchange is limited at the assembly scale, not by the channel film.

A Reynolds exponent above unity might seem to disqualify this as a heat-transfer correlation, but such exponents are typical of interfacial coefficients measured in porous media at these Re. Song and co-workers compile eight correlations for Re ≤ O(10²) with exponents from 0.2 to 1.97, adding their own, Nu = 0.001 Re^1.66 Pr^1/3, r² = 0.93 over 13 ≤ Re ≤ 133 from fifty-four tests [33]; our 1.443 over Re = 23–94 sits inside that range. The comparison is one of exponent and magnitude only — their coefficient is volumetric and particle-diameter-referenced in a packed bed, ours apparent and wall-referenced in a honeycomb — showing that a superlinear exponent recurs when this coefficient is measured rather than assumed at low Re, not that the values should agree. This is not a competing channel correlation: the channel-scale work correctly answers its own question, while ours is an assembly-scale quantity in the same non-dimensional form, and the gap between the two is a measure of how little of the structure participates under the assumptions set out in section 3.1.

The exponent itself is not the interesting quantity, since the identity of section 3.1 shows it to be ε(Re) rewritten. The interesting quantity is the transfer-unit count (Figure 3b),

$$N_{\rm prof} \propto Re_{\rm nom}^{+0.341},\qquad \pm0.004\ \text{(instrumental)},\ \pm0.040\ \text{(regression SE)},\qquad r^2 = 0.845,$$

obtained by integrating the gas energy equation through the measured piecewise-linear wall profile with uniform h, solving for N = hPL/ṁc_p reproducing each run's measured outlet temperature. Profile-corrected counts exceed isothermal-wall identity values by 1.5% to 8.7%, largest at low flow where the axial wall gradient is steepest; the identity gives NTU_app ∝ Re_nom^(+0.389 ± 0.045), retained only for comparison. The front 11 mm and rear 30 mm are unmeasured; the exponent moves only to +0.331 continuing the rear gradient at half slope, and by less than 0.003 back-extrapolating the front gradient. A full linear continuation of the rear gradient is inadmissible: for five of the fifteen runs it places the exit wall below the measured outlet gas temperature, which no single-stream model can produce, so the outlet measurement bounds the rear wall's permissible steepness. The exponent is insensitive to extrapolation at the 0.01 level, against a discrepancy of 1.34 (supplementary section S4).

Transfer units increase with flow. For a single stream exchanging with a solid whose volumetric conductance is any fixed function of axial position — h(z) independent of mass flow — the transfer-unit count is

$$N = \frac{1}{\dot m c_p}\int_0^L h(z)\,P\,dz \;\propto\; \dot m^{-1} \;\propto\; Re^{-1}$$

with no adjustable constant. On the measured properties this evaluates to −0.9996 for fixed Nusselt number and −1.017 for conductance independent of temperature and flow; we write it as −1.00 below. The profile-corrected exponent differs from it by 1.34 in the power of Re and in sign. Under a uniform δT3 = ±25 K offset it moves only from +0.320 to +0.377, remaining positive and over eight standard errors from zero, let alone −1.00 — not an artefact of outlet-probe calibration. This is the paper's central structural result, developed in section 5.2.

The auxiliary groups close off local explanations. The web Biot number never exceeds 3.4×10⁻⁵, so the solid is isothermal across web thickness and through-web conduction is not limiting, consistent with high effective conductivities reported for honeycomb substrates [18] and with SiC at the conductive end of substrates compared by Cui and Kær [20]. The Graetz number reaches only 0.704; its reciprocal, x* = 1/Gz_L, runs from 1.42 at highest flow to 5.71 at lowest, against a fully developed asymptote reached at x* of order 10⁻¹ [24]. The channel is thus thermally fully developed over all but a few percent of its length at every condition, so the constant fully developed Nusselt number is the appropriate closure; a developing-flow closure would raise h while scaling no faster than about Re^0.5, driving N toward Re^−0.5 and worsening the sign discrepancy. The axial Péclet group is at least 1460, so axial gas conduction is negligible, and the radiation–conduction number is 1.50 to 5.66, so the solid redistributes absorbed energy axially by radiation and an imposed uniform flux is not the correct boundary condition — consistent with radiation being required to reproduce measured axial temperature profiles in monolith channels [19]. Supplementary section S5 tabulates all four groups per run.

### 4.3 Local nonequilibrium

Both disequilibrium indices rise linearly with Reynolds number within every flux group (Figure 4c) and both are small in magnitude, Λ₅₈ spanning 0.031 to 0.064 and Λ₁₀₇ spanning 0.055 to 0.114. At 107 mm the three per-flux slopes agree to within ±2% — 8.72, 8.51 and 8.35 ×10⁻⁴ per unit Re at 456, 304 and 256 kW m⁻², with r² of 0.997, 0.999 and 1.000 — while offsets are ordered by flux: Λ₁₀₇ = c(G₀) + (8.51 ± 0.35)×10⁻⁴ Re, a flux-invariant slope with a flux-dependent offset falling as flux rises, intercepts 0.0440, 0.0342 and 0.0313 at 256, 304 and 456 kW m⁻². A pooled single-intercept fit averages the offsets away, gives r² = 0.897 against the 0.997–1.000 achieved within groups, and inflates the slope's standard error ninefold. At 58 mm the same structure holds at roughly a third of the slope — 1.81, 2.65 and 3.11 ×10⁻⁴ per unit Re, r² of 0.984, 0.999 and 0.992 — with Λ₅₈ rising by 22%, 55% and 35% across the flow range at the three fluxes: disequilibrium is a smooth monotonic function of depth and flow rather than a switch.

T9 and T10 are flow-exposed sheathed thermocouples, not bulk-gas measurements: each sits at a convective–radiative–conductive equilibrium among gas, surrounding solid and its own stem, and for this geometry a probe inside a monolith channel is known to be significantly perturbed by radiation even unblocked [19]. Λ is thus an apparent wall-to-interior disequilibrium rather than a true gas–solid nonequilibrium — an observation-weighted surrogate, most usefully a bound, since radiative loading biases the interior probe toward the solid and understates the true gas–solid gap. Local thermal non-equilibrium is a premise of the two-equation porous-absorber formulations [2,11,13,15] and of heterogeneous monolith models [18,21], its magnitude following from an assumed or fitted volumetric coefficient; a directly instrumented disequilibrium at two depths gives that literature something to check against, provided the comparison uses the same observation model rather than bulk-gas temperature. That the gap grows with Reynolds number is the microscopic counterpart of section 4.2: the wall-to-interior probe difference widens as flow rises, yet the integrated transfer-unit count also rises. Taken at face value both point to a participating gas–solid contact that increases with flow; the conditional is discharged in section 5.2.

### 4.4 Transient identification and the energy budget

The three cooling transients decay as a single exponential in all six receiver probes, with a per-sensor eigenvalue spread of 3.7% (Figure 5a). The abscissa x = ε ṁ c_p requires an effectiveness for each cooling run; referring each run to the heating state it decays from, with c_p on that run's own tail temperatures, gives the primary identification (Figure 5b),

$$C_{\rm eff} = 276 \pm 37\ {\rm J\,K^{-1}}, \qquad K_{\rm loss} = 0.080 \pm 0.023\ {\rm W\,K^{-1}}, \qquad r^2 = 0.986.$$

A pooled ε(q) correlation across all fifteen heating runs, averaging over the flux ordering of ε, gives 299 ± 41 J K⁻¹ and 0.095 ± 0.025 W K⁻¹ at r² = 0.964; the matched form is better-conditioned and 8% lower in capacitance. Both rest on three eigenvalues and hence one degree of freedom, the main caveat. The joint fit to all eighteen eigenvalues gives 269 ± 30 J K⁻¹ and 0.101 ± 0.022 W K⁻¹ with sixteen degrees of freedom and r² = 0.899; heating-only identification on the deep and outlet probes gives 281 ± 36 J K⁻¹ and 0.114 ± 0.028 W K⁻¹. The spread, 269 to 299 J K⁻¹ (11%), is the honest measure of how well the dataset determines the capacitance (Figure 6b). The Monte Carlo mean of the matched capacitance, 280 J K⁻¹, sits above the point estimate because C_eff is the reciprocal of a slope fitted to three points; we quote the point estimate with the Monte Carlo standard deviation and note the bias rather than average it away.

Whichever estimate is taken, C_eff exceeds the bare monolith capacitance of 42.0 to 46.8 J K⁻¹ (from the measured 40 g mass and c_p at 600–900 K) by a factor of six to seven, so on this lumped description roughly five sixths of the thermal mass governing the receiver's transient is attributable to the holder and insulation rather than to the absorber. The attribution is by subtraction — the absorber's own capacitance is measured and the rest is the remainder — and the identification does not resolve how that remainder divides between holder, insulation and the plenum hardware. This matters for start-up and cycling estimates, is invisible to steady-state characterization — one reason it is absent from a literature built on steady efficiency measurements [4,5,6,7,8,9] — and cannot be recovered from absorber material properties alone. It also qualifies He, Du and Shen's claim that porous volumetric receivers have quick transient response and small thermal inertia owing to the short flow path of the heat transfer fluid [46]: true for the absorber alone (42–47 J K⁻¹ is indeed small), but not for the installed receiver, where the greater part of the responding mass belongs to hardware the absorber's properties say nothing about.

Transient behaviour has typically been characterized via a receiver time constant from outlet air temperature rather than by identifying the underlying capacitance. The same review reports time constants of about 90 s (Sulzer), 70 s (HiTRec II), 365 s (Sandia foam), 660 s (CeramTec) and 600–840 s (SOLAIR), noting no standard exists [46]. A time constant, C_eff/(x + K_loss), measured at one flow leaves numerator and denominator unseparated — the 70–840 s spread reflects that conflation. Here, separation comes from varying advective conductance, attributing the six- to sevenfold excess specifically to holder and insulation. Detailed resolution elsewhere has been computational: a coupled transient model for start-up, shut-down, clear-sky and cloud passage [41], and ray tracing with pore-scale simulation and PID control giving 7.0 s (air) versus 45.9 s (molten salt), rising 129% as mass flow falls threefold [42]. A predicted time constant implies a predicted capacitance—the factor found here cannot be captured unless holder and insulation are modelled explicitly.

The loss conductance is bracketed between 0.080 and 0.114 W K⁻¹, with the heating tangent exceeding the matched cooling secant by 1.41× — consistent with a partly radiative loss path, though only weak evidence for T³ dominance, and we claim no more than the sign. Carried through the δT3 = ±25 K band, the primary determination moves over 271.7 to 280.5 J K⁻¹ and 0.0754 to 0.0844 W K⁻¹, and the heating determination over 262.2 to 299.4 J K⁻¹ and 0.1073 to 0.1200 W K⁻¹, so the band widens the loss bracket to 0.075 to 0.120 W K⁻¹ while leaving the capacitance spread dominated by estimator choice rather than by outlet-probe bias.


### 4.5 Similarity collapse

Rescaling time on τ = C_eff/(ε ṁ c_p + K_loss) collapses all fifteen heating transients (Figure 5c). The length-averaged wall temperature reaches half its total rise at t* = 0.205 with a coefficient of variation of 16.8% across runs spanning a fourfold flow range and a 1.8-fold flux range, and the outlet gas reaches half rise at t* = 0.634 with a coefficient of variation of 7.8%. A single slow mode therefore organizes the late heating transient, and the offset between the two half-rise values quantifies the phase lag between solid and gas. This does not establish a one-state response: the wall half-rise sits well below ln 2 while the outlet sits close to it, so the wall retains a faster component the lumped model does not represent, consistent with section 3.2's finding that the front probes carry their own short time constant. The defensible statement is that the heating transients share a common flow-dependent time scale, set by the slow eigenvalue, on which their late behaviour collapses — which is what a start-up or control model needs — and not that the receiver is a one-parameter dynamic system.





### 4.6 Delivered aperture irradiance

The aperture irradiance is determined two ways, and they disagree. The first is radiometric flux mapping by the three-step characterization procedure established for this simulator [22,23]: an alumina-coated Lambertian target on a three-axis slider, a Gardon-type circular foil radiometer traversed on a graded grid with NIST-traceable calibration and a calorimetric correction for spectral mismatch, and a parallel high-dynamic-range CCD reconstruction, with Gaussian fits validating the focal distribution to normalized RMSE of about 0.07 and 0.12 on the two axes [22]. Integrating the calibrated map over the receiver face gives the nominal values 456, 304 and 256 kW m⁻² that label the three configurations. The second is an apparent steady energy closure, which is not model-free: it takes K_loss from section 4.4 and omits storage by assuming the assembly is stationary at run end, and neither assumption is safe, since the insulation probe rises from 308 to 356 K across the campaign with a decay far slower than any run window. At steady state

$$f = \frac{Q_{\rm gas} + K_{\rm loss}(\bar T_w - T_{\rm amb})}{G_0 A_{\rm frt}}$$

should equal unity if the radiometric flux is correct. Evaluated at both ends of the K_loss bracket and averaged per configuration it does not: f runs 0.989–1.134 at 456 kW m⁻², 1.149–1.324 at 304 and 0.784–0.916 at 256, and apparent efficiency on the nominal basis ranges 0.252–1.224, exceeding unity in the 304 kW m⁻² group. We therefore report spans of 451–517, 304–402 and 201–256 kW m⁻², and every power-normalized quantity as a band across them — not as a bounding interval, since closure omits storage, inherits a lumped loss model and carries T3 bias directly into Q_gas, so neither determination is an established bound.

Supplementary section S6 shows that no uniform outlet-probe shift and no flux-only correction reconciles the two, because the required correction is flow-dependent; that part of the gap at the two higher fluxes is structural, a Gardon gauge sensing over 499 mm² being unable to resolve the flux a 361 mm² face intercepts in a focused spot, with a face-to-gauge ratio of 1.04 to 1.28 depending on spot width; and that the same difficulty of converting incident onto delivered power has been reported for this simulator elsewhere [23]. What the residual requires is a flow-dependent term in the energy budget, and peripheral bypass and channel-group maldistribution are two candidates consistent with section 4.2 — though non-stationarity and a flow-dependent probe response are not excluded, so the convergence with section 5.2 is a consistency and not an independent confirmation.

No temperature-based result depends on the choice: effectiveness, transfer units, the Nusselt correlation, the inversion criterion, the disequilibrium indices, the eigenvalue identification and the similarity collapse are all ratios of measured temperatures and flows. Only the efficiency and delivered-power-per-unit-mass figures shift, and those are reported as spans rather than as intervals containing the true value.


### 4.7 Pressure drop as installed

Steady-window means of the differential pressure channel compare poorly with the laminar square-duct prediction, computed with a Poiseuille number of 56.91, all metered flow through 100 channels, and properties at T̄_g. The measured signal exceeds the monolith-only prediction by 1.4 to 2.8 times at high flow (2.74 against 0.99 mbar at 15.3 sL min⁻¹), falls below it at low flow, and goes negative at the lowest (−0.02 against 0.20 mbar at 4.5 sL min⁻¹) — the signature of a zero offset on a transducer whose ±0.2 mbar accuracy floor, 0.1% of its ±200 mbar span, is the same order as the entire expected monolith drop of 0.20 to 1.0 mbar. The taps evidently span more than the monolith faces, and a second differential channel reads a flat −0.26 to −0.30 mbar in every run and appears unconnected.

This matters because pressure drop is how the volumetric-receiver literature diagnoses honeycomb flow behaviour. Absorbers whose pressure loss is linear in velocity tend inherently towards flow instability [27], high-conductivity materials with a predominantly quadratic characteristic are the stable ones [28], and instability becomes possible when one pressure-drop level maps to more than one outlet temperature, letting a region drift away from the operating point and fail locally if it exceeds the material limit [26]; Ávila-Marín's review extends the criterion to highly porous honeycombs specifically [1], and permeability and Forchheimer coefficients from pressure-drop data are the standard inputs to homogeneous equivalent models [13]. Recent unheated monolith measurements place a honeycomb on the linear, instability-prone side of that axis: the first experimental validation of a dual-cell-density pressure-drop model finds a linear relationship between Reynolds number and pressure drop in every configuration tested while progressively underestimating the drop as Re rises through entrance and exit effects [45], and those effects are large enough at channel Re of 50 to 320 to require an explicit correction factor, strongly dependent on spacing and increased by misalignment [44]. Our transducer, spanning monolith plus plumbing, cannot separate the absorber's characteristic from those comparable losses. We report it because the fix is inexpensive and because section 5.3 shows this to be the measurement that would discriminate between the two candidate mechanisms above; supplementary section S8 gives the per-run comparison.

## 5. Two consequences for reduced-order receiver modelling

### 5.1 Against validating receiver models on a gas outlet temperature and one or two solid temperatures

The standard validation dataset for a volumetric receiver model is an outlet gas temperature together with one or two embedded solid temperatures. Almost every campaign surveyed above reports precisely this — outlet air temperature and efficiency, with front-face temperature by pyrometer or infrared camera where available [4,5,6,7,8,9] — and the convention is explicit at review level: He and co-workers name outlet air temperature and thermal efficiency as the two key indicators for solar receivers, noting that a non-uniform radiation distribution produces a non-uniform air temperature distribution requiring more thermocouples than campaigns typically deploy [46]. Modelling studies calibrate and validate against those same quantities [10,11,12,13], and recent practice has not moved past it [29,31,34]. This campaign contains the same measurements several times over, permitting a direct test of how much such a dataset constrains. The answer is very little, far less than the apparent quality of the fit suggests.

The gas outlet temperature is nearly blind to the internal field. At 456 kW m⁻² T3 varies by 26 K across the flow range while the front-face wall varies by 264 K; at 304 kW m⁻² the ratio is 51 K against 331 K. Worse, T3's flow slope changes sign between configurations, so the qualitative behaviour a model must match is inconsistent across the campaign.

Which solid probe is called the solid temperature changes every identified coefficient. Recomputing the entire analysis — same fifteen runs, reduction, and fitting procedure, varying only which thermocouple defines T_s — gives Figure 6a:

| Solid reference | $\varepsilon$ range | $Nu$ prefactor | $Nu$ exponent | $\varepsilon^*$ (456/304/256) |
|---|---|---:|---:|---|
| front wall only (T8) | 0.454 – 0.851 | 5.95×10⁻⁵ | 1.841 | 0.605 / 0.591 / 0.572 |
| mid wall only (T12) | 0.520 – 0.710 | 2.67×10⁻⁴ | 1.426 | 0.604 / 0.590 / 0.570 |
| wall quadrature (adopted) | 0.573 – 0.781 | 3.10×10⁻⁴ | 1.443 | 0.666 / 0.655 / 0.642 |
| wall + interior mean | 0.608 – 0.803 | 3.70×10⁻⁴ | 1.419 | 0.691 / 0.682 / 0.671 |
| interior probes (T9, T10) | 0.648 – 0.826 | 4.49×10⁻⁴ | 1.393 | 0.719 / 0.711 / 0.702 |
| rear wall only (T11) | 0.771 – 0.823 | 1.72×10⁻³ | 1.092 | 0.787 / 0.787 / 0.791 |

The prefactor spans a factor of 28.9, the exponent 1.09 to 1.84, the inversion threshold 0.57 to 0.79. Every row is a defensible unremarked choice, and every row is a different receiver. Interior and wall probes at the same depth differ by 21–55 K at 107 mm and 21–27 K at 58 mm, gaps growing with flow, so no single solid temperature represents the section. The same sensitivity is documented channel-side [17,19], but at assembly scale the effect is a factor of twenty-nine, not tens of percent. The transient identification is equally under-determined: the same eighteen transients yield C_eff of 122, 269, 276 or 281 J K⁻¹ depending on probe subset and fit window, the window alone spanning 200–295 J K⁻¹, with the coefficient of determination increasing as the estimate degrades (Figure 6b). This factor of 2.3 is the identifiability basin's width, not measurement scatter.

Agreement between a receiver model and a gas outlet temperature plus one or two solid temperatures is therefore not evidence of correct internal physics; the identified coefficient reflects sensor placement rather than the receiver, and will not transfer to a different geometry, scale or operating envelope. This offers a deflationary explanation for a two-decade puzzle: most tested volumetric absorbers underperformed predictions [1,3], four established models cluster within 25 K of each other and 80–100 K from experiment with mean-air-temperature measurement named a principal suspect [10], and the idealized equilibrium model is optimistic by over twenty efficiency points against practice [2]. If the calibration target does not constrain the internal field, models can be simultaneously well fitted and structurally wrong.

Two external results sharpen this. Iding and co-workers reproduce the front temperature of the 1080-cup Jülich receiver to 6.78 K RMS and behind-cup air temperature to 7.51 K, on exactly the two-quantity dataset above, yet report this per-cup performance does not carry to the array [30]: the field resolves the irradiated surface to a few kelvin but not depth, where identified coefficients live. Conversely, Song and co-workers resolve the interfacial coefficient in porous media at low Reynolds number by measuring solid and gas temperatures simultaneously at each of a series of axial stations and solving the discretized gas energy equation directly, rather than inverting against a temperature field, which they note is ill-posed [33]; their term for co-located temperatures is data integrity, and the present results quantify the error from choosing the solid reference instead.

We recommend that reported receiver-model validations state explicitly which solid temperature defines the reference and, where the instrumentation permits, report the identified coefficient's sensitivity to that choice; where it does not permit, the coefficient should be reported as conditional rather than as a receiver property. This is not an argument against reduced-order models but that the standard validation dataset cannot discriminate among them — the discriminating measurements are different ones.

### 5.2 Against uniform-temperature or uniform-flux heat transfer coefficients in zero- and one-dimensional receiver models

Zero- and one-dimensional receiver models close the gas–solid coupling with a coefficient from a correlation derived for a duct with uniform wall temperature or uniform wall flux — a Churchill cross-flow correlation in the inlet region [10], volumetric forms Nu_lv = a Re^b Pr^c [13], the Fu and Kuwahara forms [15], Graetz-problem Nusselt and Sherwood numbers at constant wall temperature [21]. That closure is inadmissible, and no model is needed to show it. Every such model reduces to the single-stream form

$$\dot m c_p \frac{dT_g}{dz} = h(z)\,P\,[T_s(z) - T_g(z)],$$

and if h is fixed in position — as a Nusselt number on unchanging local geometry implies — integration gives N = (1/ṁc_p)∫h P dz, whose flow dependence is that of 1/ṁ alone. Section 4.2 evaluates this on measured properties as −0.9996 but measures +0.341 ± 0.040 instead: the sign is wrong, and each local repair is excluded — a developing-flow correlation [17,24] scales h no faster than Re^0.5 and worsens the disagreement; solid-side resistance is excluded by Bi ≤ 3.4×10⁻⁵ and axial gas conduction by Re·Pr·L/D_h ≥ 1460; a uniform imposed flux is also wrong, since N_rc of 1.50 to 5.66 means the solid redistributes absorbed energy axially by radiation [14,15,19].

That inversion assumes spatially uniform h, excluding only a constant-Nusselt closure, not every fixed h(z). The broader claim is tested across three families in which h(z) is piecewise linear and non-negative on evenly spaced nodes, fitted by multi-start least squares in log h, asking whether one flow-independent profile reproduces the fifteen measured outlet temperatures: a single dimensional h(z) shared by all runs; a separate five-node profile per irradiance configuration (so a conductance differing with lamp setting or temperature level is not held against the hypothesis); and a shared Nusselt profile, h_i(z) = Nu(z) k_i/D_h with each run's own conductivity — the direct test of the constant-Nusselt class.

None succeeds; the discriminating statistic is the ordering of the residual in flow, not its magnitude. For the shared conductance profile, the RMS outlet residual is 64 K with two nodes and saturates at 55, 54 and 54 K with three, five and seven — freedom is not the issue — while residuals correlate with ln Re at r = −0.95, falling 126–148 K per e-fold, from about +80 K at the lowest flow of each configuration to −70 K at the highest. The shared Nusselt profile behaves identically, at 56 and 55 K with three and five nodes. The per-configuration profile lowers the residual in one group and raises it in another (36, 68, 52 K at 256, 304, 456 kW m⁻²) but sharpens the ordering to r = −0.996, −0.982, −0.997: a flow-independent profile fails monotonically even within a single lamp configuration. A fixed conductance over-predicts outlet temperature at low flow and under-predicts it at high flow in every family, a trend the ±25 K outlet-bias band cannot absorb. Two limits apply: the test presumes the single-stream form, excluding only flow-independent conductance within it; and it inherits the wall reconstruction's bounded but unmeasured extrapolations. Neither affects the sign or ordering. All three families, fitted nodes and per-run residuals are archived under `fixed_profile_test`.

Taken at face value, the measurement says participating gas–solid contact grows with Reynolds number. That conditional is doing work: an outlet-probe bias varying systematically with flow, or a departure from thermal stationarity varying with run duration and hence with flow, would produce the same signature without any change in participation, and neither is bounded by this campaign — which is why section 5.3's first two measurements target them. Read physically, what a single-stream model fixes but this receiver need not is the distribution of flow among channels. Two mechanisms, advanced as leading hypotheses rather than identified causes, produce the required sign. Viscosity-driven redistribution at common pressure drop: hot core channels choke relative to cooler perimeter channels; as total flow rises the field cools, the viscosity contrast weakens, and the participating fraction grows — predicting transfer units and the wall-to-interior deficit both rising with Re and a front-to-mid wall crossing, all observed in sections 4.1–4.3, matching the coupling reported for porous SiC absorbers where non-uniform hydrodynamics are significantly influenced by increased viscosity at the hot spot [9]. Peripheral bypass: a leakage path through the housing clearance whose share falls as total flow rises, being cooler with resistance a weaker function of temperature, gives the same sign while inflating a nominal Reynolds number. The two are indistinguishable from temperature data alone, both require a flow-distribution degree of freedom absent from a single-stream model, and differ in that redistribution conserves the flow path while bypass short-circuits it — the distinction section 5.3 measures.

The first mechanism is not a postulate invented for this measurement: it is the flow-instability criterion carried by the volumetric-receiver literature for three decades (section 4.7 [26,27,28]), restated as continuous rather than catastrophic — the linear regime is viscosity-dominated, using the viscosity of the heated gas, so the contrast driving runaway overheating in an unstable absorber must, in a stable one, still redistribute flow reversibly, weakening as total flow rises. That is the sign measured, and the criterion attaches specifically to highly porous honeycombs [1]. The parallel-channel literature bears on it structurally rather than quantitatively, both results being for flow boiling and therefore cited for structure and not for transfer of values: lateral channel-to-channel conductance alleviates maldistribution with a single eigenvalue setting its onset [39], consistent with radial heat transfer reducing instability probability in porous absorbers [28], while a maldistributed state can be thermodynamically favoured over an even one [40]. The radiation–conduction number measured here is such a lateral conductance, so transfer units rising with Reynolds number occur in the presence of, not the absence of, a stabilizing mechanism.

The remedy is a multi-stream or distributed formulation with a pressure-coupled flow split, not a better Nusselt correlation — precedented in pseudochannel models built because single-channel models cannot predict accurately under flow maldistribution [35], flow-distribution models for dual-cell-density monoliths [36] and dual-zone packed beds [43], and solar-field models at panel and module scale [34,37]. A dual-cell-density monolith with a Darcy momentum sink, isothermal, gives a flow fraction depending only on relative permeability and core size, not flow rate [36], since flow rate cancels from the equal-pressure-drop condition when every path is linear in velocity with one shared viscosity. Adding an inertial term via the Ergun equation makes the split Reynolds-dependent [43]. A flow-dependent split thus needs an inertial contribution or a path-to-path viscosity difference, and both are present: the heated absorber carries a non-uniform viscosity field [9], and the viscous-to-inertial transition is the axis of the stability criterion itself [27,28]. Buoyancy supplies a third such agent in molten-salt panels [37], and in all three the split is most differentiated at low flow — flow-dependent participation is a general property of heated parallel-path arrays. Supplementary section S7 develops the comparison.

The scope claimed is precise: not that no continuum model can reproduce these data, nor an impossibility over all parameters, but that a single-stream formulation whose volumetric conductance is a fixed function of axial position cannot produce a transfer-unit count increasing with flow, and this receiver's does. Any such model, however calibrated, will misrepresent the flow response, and since such models are typically calibrated against the dataset section 5.1 shows uninformative, the misrepresentation can survive an apparently successful validation. This does not contradict the catalysis literature's continuum reductions but delimits them: a pseudo-continuous model reproduces three-dimensional CFD to better than 0.5 K assuming flow equally distributed over all channels [21], and effective radial conductivities are obtained with the gas static [20] — reasonable for a reactor with a heated wall and metered uniform inflow. What fails here is not the reduction but its use under a concentrated, spilling beam with an insulated housing and a flow split the boundary conditions do not control.

### 5.3 What would resolve it

Three measurements would settle the questions left open, ranked by value per unit cost. Most valuable is a differential pressure measurement across the monolith faces alone with a ±2 to ±5 mbar transducer. The two candidate mechanisms of section 5.2 have different pressure-drop signatures at the same total flow, since redistribution conserves the flow path while bypass short-circuits it, so a pressure-drop series across the flow range — preferably with a sealed-clearance control run to remove the bypass path — would discriminate them; a single point would not, since the two differ in the flow dependence of the characteristic rather than in its value at one flow. The expected signal is 0.20 to 1.0 mbar, resolvable to a few percent with an appropriately ranged transducer. It would also connect this receiver to the permeability and Forchheimer characterization the homogeneous-equivalent modelling literature depends on [13], and place the absorber on one side or the other of the linear-versus-quadratic flow-stability criterion [26,27,28]. Its practicability at this scale is shown by a calorimetric test bench for pressurised gas receivers carrying absolute pressure sensors before and after the absorber alongside the thermocouples [32]. No measurement of channel-to-channel flow split in an irradiated honeycomb absorber appears in the recent literature: Broeske and co-workers state that individual mass-flow measurements for each of their sixteen absorber modules were impossible for reasons of space, so efficiencies were evaluated through radiative energy balance rather than air enthalpy flow [29]; the nearest published measurements isolate the effect in surrogates, varying air velocities in a two-tube fluidized-particle receiver at ambient temperature to impose controlled thermal heterogeneity and observe the resulting particle mass-flow split [38]. The measurement proposed here is not a replication.

Second is an independent determination of outlet enthalpy, by mixing calorimeter or shielded aspirated probe, which would break the T3 dependence propagating into effectiveness, transfer units, the Nusselt correlation, the inversion threshold, efficiency, and the identified constants, addressing the mean-air-temperature difficulty named as a principal source of the field's model–experiment gap [10]. Third is a radial thermocouple traverse at one depth to resolve the core-to-perimeter solid profile, since section 5.1 identifies a single solid temperature per depth as the binding limitation on identifiability; two radial positions at one depth would test the redistribution mechanism directly and replace the reference-choice sensitivity of Figure 6a with a measured profile.


### 5.4 Limitations

Several limitations bound what has been established. The isothermal-wall ε–NTU relation approximates a wall spanning up to 321 K; the bias is quantified — profile-corrected transfer units exceed identity values by 1.5% to 8.7%, shifting the exponent from +0.389 to +0.341 — but that correction inherits the wall reconstruction, whose front 11 mm and rear 30 mm are extrapolated rather than measured, the rear extrapolation bounded by the outlet temperature but not determined by it. The wall probes' installation error, comprising contact resistance and cavity radiation, is unquantified. The outlet-probe systematic bias is bounded by declaration rather than measurement, and only against a uniform offset: a bias varying with flow is unbounded and would act directly on the exponents carrying the paper's argument. Thermal stationarity at the end of each heating window is assumed, not demonstrated. The primary transient identification rests on one degree of freedom. The inversion threshold at 256 kW m⁻² is sensitive to the crossing estimator at a level comparable to its flux dependence, and that group contains a single negative point. Flux-map uncertainty is characterized by its discrepancy with the energy closure rather than propagated from the instrument chain. Pressure drop is unusable as installed.

Two thermal models are solved: the profile-corrected transfer-unit count and the fixed-conductance falsification of section 4.2 both integrate a single-stream gas energy equation through the measured wall profile, deliberately of the class the paper argues against, used to test that class rather than correct any measurement. Nothing here resolves the internal flow field, and no reported measurement is corrected using a model of it.


## 6. Conclusions

Gas–solid exchange in this structured SiC receiver is limited at the assembly scale rather than in the channel. The apparent Nusselt number referenced to the length-averaged wall temperature follows Nu_app ∝ Re_nom^1.470 ± 0.011 with r² = 0.9994 under a common-slope, per-configuration-prefactor fit over Reynolds numbers of 23 to 94, prefactors 2.55 to 3.20×10⁻⁴; the pooled form gives 1.443 ± 0.069 at r² = 0.971. Either way this is 14.5 to 111 below the fully developed laminar square-duct value Nu_H2 = 3.091 and far below channel-level monolith-literature values, with a web Biot number of at most 3.4×10⁻⁵ excluding a through-web conduction resistance as the limiting one, though not every solid-side or assembly-scale resistance. The front-to-mid side-wall crossing collapses onto a single dimensionless marker under the adopted wall-averaging convention, the mid-depth wall overtaking the front face once the gas has recovered about two thirds of the available wall-to-ambient temperature difference, at ε* = 0.666, 0.655 and 0.642; the criterion is robust in form, but its weak flux dependence, comparable to the estimator systematic, is unresolved. On this plain, ungraded honeycomb the side-wall signature of volumetric behaviour is thus at least as much an operating condition as a design property. The apparent wall-to-interior disequilibrium grows linearly with flow at both depths, with a flux-invariant slope of (8.51 ± 0.35)×10⁻⁴ per unit Re at 107 mm and roughly a third of that at 58 mm — a smooth function of depth and flow, not a threshold phenomenon — giving the two-equation modelling literature a measured number to check against.

The receiver's transient is set by six to seven times the absorber's own thermal mass. The slow eigenvalue is linear in the gas advective conductance, giving C_eff = 276 ± 37 J K⁻¹ from cooling and 269 ± 30 J K⁻¹ from a joint fit to all eighteen eigenvalues against a bare monolith capacitance of 42 to 47 J K⁻¹, with loss conductance bracketed at 0.080 to 0.114 W K⁻¹. All fifteen heating transients rescaled on this eigenvalue collapse to a wall half-rise of t* = 0.205 within 16.8%, though not a one-state response since the wall retains a faster component the lumped model omits.

A transfer-unit count increasing with flow rules out a widely used reduced-order model class: integrated from the measured axial wall profile, it scales as Re_nom^+0.341 ± 0.041, opposite in sign to the Re_nom^−1.00 required by any single-stream formulation with position-fixed volumetric conductance, since channels are fully developed over at least 96% of their length, giving −0.9996 with developing-flow closure only widening the gap. Neither uniform, shared, nor separately fitted two-to-seven-node conductance/Nusselt profiles reproduce the outlet-temperature flow response. An empirical h ∝ Re^1.34 reproduces the exponent by construction but lacks mechanism or predictive content beyond the fitted range; a flow-distribution degree of freedom would supply one, but this campaign cannot discriminate it from flow-dependent probe bias or non-stationarity.

Validation against a gas outlet temperature and one or two solid temperatures does not constrain receiver physics: from the same fifteen runs, changing the reference solid thermocouple moves the identified Nusselt prefactor by a factor of 28.9, the exponent from 1.09 to 1.84, and the inversion threshold from 0.57 to 0.79; changing probe subset and fit window moves identified capacitance from 122 to 281 J K⁻¹. At the highest flux the outlet gas spans 26 K across the flow range while the front-face solid spans 264 K.

Delivered aperture irradiance is not determined to a value. Radiometric flux mapping and steady energy closure disagree by −1 to +13%, +15 to +32%, and −22 to −8% across the three configurations, with flow-dependent discrepancy precluding a uniform correction; a Gardon gauge over 499 mm² cannot resolve flux from a 361 mm² face in a focused spot, giving a geometric bias of 1.04 to 1.28. Reported spans are 451–517, 304–402, and 201–256 kW m⁻². No temperature-based result depends on this choice.

The discriminating measurement is a properly ranged differential pressure across the monolith faces, taken as a series across the flow range: with an accuracy floor of ±0.2 mbar against an expected 0.20 to 1.0 mbar signal, the installed transducer cannot constrain flow distribution, whereas a ±2 to ±5 mbar device — with a sealed-clearance control run — would separate viscosity-driven redistribution from peripheral bypass.

## Data availability and reproducibility

Every derived table, uncertainty estimate and figure input in this paper and in the supplementary material is produced by one script, `receiver_reduction.py`, from the raw logger files with no hand-entered measured values; the figures are then drawn from those outputs by a second script, `make_figures.py`, which writes Figures 2 to 6. Figure 1 is a schematic of the apparatus and sensor positions, drawn rather than computed, and is the only figure not produced by the pipeline. Both scripts are seeded, so a rerun on the same software stack reproduces the archived outputs; the archive was produced with Python 3.11.16, NumPy 2.4.6, SciPy 1.17.1, pandas 2.3.3 and Matplotlib 3.11.1, and the multi-start optimizer of section 5.2 is the one step whose last digits may differ across BLAS builds.

Archived outputs comprise `groups.csv` with the per-run dimensionless groups, `groups_replicates.csv` for E82 and E83, `eigenvalues.csv` with the eighteen eigenvalues under both heating sensor rules together with the per-run controller shares and the matched identification abscissa, `reference_sensitivity.csv` for section 5.1, `pressure_drop.csv` for section 4.7, `flux_geometry.json` for the gauge-averaging calculation of section 4.6, `uncertainty.csv` for the 4000-realization Monte Carlo in both controller-correlation limits, and `results.json` containing every reported scalar — including `ntu_profile` for the profile-corrected transfer units with their pooled, grouped and T3-band exponents, `nusselt.grouped` and `nusselt.per_flux` for the correlations, `fixed_profile_test` for the three falsification families with their fitted nodes and per-run residuals, and `wall_extrapolation` for the reconstruction sensitivity; the supplementary tables are written alongside as `tableS3_heating_conditionality.md`, `tableS4_wall_and_falsification.md` and `tableS5_auxiliary_groups.md`. The supplementary material and its accompanying README document the reproduction commands and the archived contents in full.





## Tables

Table 1 gives the campaign envelope: per-run flow, dimensionless groups, temperatures, effectiveness, transfer units, apparent Nusselt number, profile-corrected transfer units, wall-to-interior disequilibrium indices, inversion indicator and nominal-basis efficiency for all fifteen heating runs, generated as `table1_envelope.md`. Table 2 gives the identified constants with Monte Carlo standard deviations, 95% percentile intervals and, for fitted slopes, the regression standard error alongside, including all four capacitance determinations and the structural transfer-unit exponent, generated as `table2_constants.md`.





## Figures

Figure 1 shows the apparatus, receiver geometry and thermocouple positions. Figure 2 shows steady wall-chain and outlet gas temperatures against flow at the three configurations, displaying the depth-graded flow response and the near-invariance of the outlet gas. Figure 3 shows the apparent Nusselt number against Reynolds number with the fully developed laminar reference, and the transfer-unit count against Reynolds number with the Re⁻¹ requirement of a position-fixed single-stream conductance. Figure 4 shows the inversion indicator against flow with the located crossings, the same data against effectiveness showing the collapse at ε* ≈ 0.65, and the wall-to-interior disequilibrium indices at both depths against Reynolds number. Figure 5 shows the cooling decays for all six probes at three flows, the slow eigenvalue against gas advective conductance for the three identifications, and the similarity collapse of the fifteen heating transients. Figure 6 shows the apparent Nusselt correlation recomputed against five different solid reference probes from the same fifteen runs, and the effective capacitance from the same eighteen transients under four defensible analysis choices with the fit-window range and the bare monolith capacitance marked.

---

---

## References

[1] A. L. Ávila-Marín, Volumetric receivers in Solar Thermal Power Plants with Central Receiver System technology: A review, *Solar Energy* 85 (2011) 891–910.

[2] A. Kribus, Y. Gray, M. Grijnevich, G. Mittelman, The promise and challenge of solar volumetric absorbers, *Solar Energy* 110 (2014) 463–481.

[3] A. L. Ávila-Marín, C. Caliot, et al., Modelling strategies for porous structures as solar receivers in central receiver systems: A review, *Renewable and Sustainable Energy Reviews* 111 (2019) 15–33.

[4] T. Fend, B. Hoffschmidt, R. Pitz-Paal, O. Reutter, P. Rietbrock, Porous materials as open volumetric solar receivers: Experimental determination of thermophysical and heat transfer properties, *Energy* 29 (2004) 823–833.

[5] T. Fend, R. Pitz-Paal, O. Reutter, J. Bauer, B. Hoffschmidt, Two novel high-porosity materials as volumetric receivers for concentrated solar radiation, *Solar Energy Materials & Solar Cells* 84 (2004) 291–304.

[6] A. L. Ávila-Marín, M. Alvarez de Lara, J. Fernández-Reche, Experimental results of gradual porosity volumetric air receivers with wire meshes, *Renewable Energy* 122 (2018) 339–353.

[7] A. L. Ávila-Marín, J. Fernández-Reche, S. Gianella, et al., Experimental study of innovative periodic cellular structures as air volumetric absorbers, *Renewable Energy* 184 (2022).

[8] F. Zaversky, L. Aldaz, et al., Numerical and experimental evaluation and optimization of ceramic foam as solar absorber – Single-layer vs multi-layer configurations, *Applied Energy* 210 (2018) 351–375.

[9] P. Wang, J. B. Li, R. N. Xu, et al., Non-uniform and volumetric effect on the hydrodynamic and thermal characteristic in a unit solar absorber, *Energy* 225 (2021).

[10] R. Capuano, T. Fend, P. Schwarzbözl, et al., Numerical models of advanced ceramic absorbers for volumetric solar receivers, *Renewable and Sustainable Energy Reviews* 58 (2016) 656–665.

[11] T. Fend, P. Schwarzbözl, O. Smirnova, D. Schöllgen, C. Jakob, Numerical investigation of flow and heat transfer in a volumetric solar receiver, *Renewable Energy* 60 (2013) 655–661.

[12] T. Fend, O. Smirnova, D. Schöllgen, Numeric modelling of a compact high temperature heat exchanger, German Aerospace Center, Institute of Solar Research (COMSOL Conference contribution), 2011.

[13] A. L. Ávila-Marín, C. Caliot, M. Alvarez de Lara, et al., Homogeneous equivalent model coupled with P1-approximation for dense wire meshes volumetric air receivers, *Renewable Energy* 135 (2019) 908–919.

[14] M. A. A. Mendes, P. Talukdar, S. Ray, D. Trimis, Detailed and simplified models for evaluation of effective thermal conductivity of open-cell porous foams at high temperatures in presence of thermal radiation, *International Journal of Heat and Mass Transfer* 68 (2014) 612–624.

[15] R. C. Moro Filho, W. Malalasekera, An Analysis of Thermal Radiation in Porous Media Under Local Thermal Non-equilibrium, *Transport in Porous Media* 132 (2020) 683–705.

[16] F. Salih, K. E. Kakosimos, Exploring and enhancing the volumetric behaviour in solar receivers through specular reflectivity and simplified design, *Energy Conversion and Management* 291 (2023) 117306.

[17] I. Cornejo, P. Nikrityuk, R. E. Hayes, Improved Nu number correlations for gas flow in monolith reactors using temperature-dependent fluid properties, *International Journal of Thermal Sciences* 155 (2020) 106419.

[18] R. E. Hayes, A. Rojas, J. Mmbaga, The effective thermal conductivity of monolith honeycomb structures, *Catalysis Today* 147 (2009) S113–S119.

[19] R. E. Hayes, I. Cornejo, Multi-scale modelling of monolith reactors: A 30-year perspective from 1990 to 2020, *The Canadian Journal of Chemical Engineering* 99 (2021) 2589–2606.

[20] X. Cui, S. K. Kær, Two-dimensional thermal analysis of radial heat transfer of monoliths in small-scale steam methane reforming, *International Journal of Hydrogen Energy* 43 (2018) 11952–11968.

[21] D. Schlereth, O. Hinrichsen, Comparison of a Pseudocontinuous, Heterogeneous 2D Conductive Monolith Reactor Model to a 3D Computational Fluid Dynamics Model, *Industrial & Engineering Chemistry Research* 53 (2014).

[22] D. Haseler, A. M. Ali, K. E. Kakosimos, Positioning and focusing of light sources in light concentrating systems using convolutional neural network modelling, *Solar Energy* 218 (2021) 445–454.

[23] T. M. Abdellateif, J. Sarwar, et al., K. E. Kakosimos, Optical and experimental evaluation of a directly irradiated solar reactor for the catalytic dry reforming of methane, *Chemical Engineering Journal* 452 (2023) 139190.

[24] R. K. Shah, A. L. London, *Laminar Flow Forced Convection in Ducts*, Academic Press, New York, 1978.

[25] W. M. Kays, A. L. London, *Compact Heat Exchangers*, 3rd ed., McGraw-Hill, New York, 1984 (reprint ed., Krieger, 2018).

[26] R. T. Broeske, P. Schwarzbözl, B. Hoffschmidt, Multi-dimensional numerical analysis of flow instabilities in 3D-shaped honeycomb absorbers, *Solar Energy* 247 (2022) 86-95. doi:10.1016/j.solener.2022.10.007

[27] A. Kribus, H. Ries, W. Spirkl, Inherent Limitations of Volumetric Solar Receivers, *Journal of Solar Energy Engineering* 118 (1996) 151-155. doi:10.1115/1.2870891

[28] M. Becker, T. Fend, B. Hoffschmidt, R. Pitz-Paal, O. Reutter, V. Stamatov, M. Steven, D. Trimis, Theoretical and numerical investigation of flow stability in porous materials applied as volumetric solar receivers, *Solar Energy* 80 (2006) 1241-1248. doi:10.1016/j.solener.2005.11.006

[29] R. T. Broeske, P. Schwarzbözl, L. Birkigt, S. Vasic, S. Dung, T. Doerbeck, B. Hoffschmidt, Experimentally assessed efficiency improvement of innovative 3D-shaped structures as volumetric absorbers, *Renewable Energy* 218 (2023) 119220. doi:10.1016/j.renene.2023.119220

[30] K. Iding, D. Zanger, D. Maldonado Quinto, R. Pitz-Paal, A Real-Time Capable Simulation of Open Volumetric Receiver Surface Temperatures With Spatially High Resolution, *SolarPACES Conference Proceedings* 1 (2024) art. 885. doi:10.52825/solarpaces.v1i.885

[31] V. R. Patil, F. Kiener, A. Grylka, A. Steinfeld, Experimental testing of a solar air cavity-receiver with reticulated porous ceramic absorbers for thermal processing at above 1000 °C, *Solar Energy* 214 (2021) 72-85. doi:10.1016/j.solener.2020.11.045

[32] D. D'Souza, M. Montes, M. Romero, J. González-Aguilar, Experimental assessment of different compact flow channel geometries on pressurised gas solar receivers, *Applied Thermal Engineering* 266 (2025) 125634. doi:10.1016/j.applthermaleng.2025.125634

[33] Z. Song, J. Wang, S. Hui, R. Dai, Data integrity resolving gas-solid interfacial heat transfer coefficient in porous media under low Reynolds number, *International Journal of Heat and Mass Transfer* 255 (2026) 127803. doi:10.1016/j.ijheatmasstransfer.2025.127803

[34] V. D. Kumar, G. Singh, L. Chandra, S. Mukhopadhyay, R. Shekhar, A Multi-Zone Unsteady Heat Transfer Model for an Open Volumetric Air Receiver: A Step Towards Scale-Up and Design Optimization, *International Journal of Heat and Mass Transfer* 191 (2022) 122747. doi:10.1016/j.ijheatmasstransfer.2022.122747

[35] P. C. Nagarajan, H. Ström, J. Sjöblom, A reduced order pseudochannel model accounting for flow maldistribution in automotive catalysis, *Scientific Reports* 15 (2025) 5082. doi:10.1038/s41598-025-89756-w

[36] C. Reinao, I. Cornejo, A Model for the Flow Distribution in Dual Cell Density Monoliths, *Processes* 11 (2023) 827. doi:10.3390/pr11030827

[37] C. Schwager, J. Schulte, M. Binder, C. J. T. Boura, U. Herrmann, Flow Distribution in Molten Salt Receiver Panels at High Flux Gradients, *SolarPACES Conference Proceedings* 2 (2024) art. 821. doi:10.52825/solarpaces.v2i.821

[38] G. Sahuquet, R. Gueguen, L. Fontalvo, S. Mer, A. Toutant, F. Bataille, G. Flamant, Particle Flow Distribution in a Fluidized-Particles Multitube Solar Receiver, *SolarPACES Conference Proceedings* 1 (2024) art. 717. doi:10.52825/solarpaces.v1i.717

[39] M. E. Rahman, J. A. Weibel, The mechanism by which thermal coupling leads to an asymptotic maldistribution stability boundary for flow boiling in parallel channels, *International Journal of Heat and Mass Transfer* 228 (2024) 125671. doi:10.1016/j.ijheatmasstransfer.2024.125671

[40] T. Aka, S. Narayan, An entropic understanding of flow maldistribution in thermally isolated parallel channels, *International Journal of Heat and Mass Transfer* 227 (2024) 125564. doi:10.1016/j.ijheatmasstransfer.2024.125564

[41] S. Sharma, P. Talukdar, Dynamic Performance Characteristics of a Porous Volumetric Solar Receiver Under Transient Flux Conditions, *Journal of Solar Energy Engineering* 145 (2023) 041009. doi:10.1115/1.4056622

[42] S. Du, Y. He, M. Li, Z. Liu, Transient response performances and control strategy of porous volumetric solar receivers with different heat transfer fluids, *Solar Energy* 283 (2024) 112984. doi:10.1016/j.solener.2024.112984

[43] P. Díaz, C. Reinao, I. Cornejo, A predictive model for flow distribution in dual-zone packed bed reactors, *Chemical Engineering Journal* 515 (2025) 162962. doi:10.1016/j.cej.2025.162962

[44] I. Cornejo, A new model for pressure drop correction for series-arranged misaligned monoliths, *Chemical Engineering Science* 299 (2024) 120515. doi:10.1016/j.ces.2024.120515

[45] C. Reinao, P. Diaz, I. Cornejo, A new model for the pressure drop in dual-cell density monoliths: Experimental validation and CFD, *Chemical Engineering Journal* 516 (2025) 163599. doi:10.1016/j.cej.2025.163599

[46] Y.-L. He, S. Du, S. Shen, Advances in porous volumetric solar receivers and enhancement of volumetric absorption, *Energy Reviews* 2 (2023) 100035. doi:10.1016/j.enrev.2023.100035
