# Transient Thermal Characterization and Dimensionless Correlations for a Structured Silicon Carbide Volumetric Solar Receiver

Aysha Melhim¹, Athanasios G. Konstandopoulos²,³, Konstantinos Kakosimos¹,*

¹ Chemical Engineering Program, Texas A&M University at Qatar, Doha, Qatar
² Chemical Process and Energy Resources Institute, CERTH, Thessaloniki, P.O. Box 60361, 57001, Greece
³ Aerosol & Particle Technology Laboratory, Department of Chemical Engineering, Aristotle University, University Campus, Thessaloniki, 54124, Greece
\* Corresponding author: k.kakosimos@qatar.tamu.edu

---

## Abstract

Transient and steady-state temperature measurements from a laboratory-scale structured silicon carbide (SiC) volumetric receiver are analyzed to obtain dimensionless heat-transfer characteristics without recourse to detailed simulation. The receiver, a square-channel honeycomb monolith (137 mm length, 19 mm width, 100 channels of 1.5 mm), was operated in a high-flux solar simulator at three aperture irradiances (456, 304, and 256 kW m⁻²) and five air flow rates per irradiance, with three additional cooling experiments. The data reduce to a small set of quantitative results. The apparent global Nusselt number follows Nu = 3.1×10⁻⁴ Re^1.44 for 23 ≤ Re ≤ 94, one to two orders of magnitude below fully developed laminar duct theory, which shows that heat exchange is limited at the scale of the monolith assembly rather than the channel boundary layer. The axial wall-temperature peak moves from the front face into the interior when the gas-side heat-exchange effectiveness exceeds a flux-independent threshold, ε* = 0.66 ± 0.03, providing an operating-point criterion for volumetric behavior. The interior-to-wall temperature deficit deep in the receiver grows linearly with Reynolds number, Λ₁₀₇ = 0.038 + 8.3×10⁻⁴ Re, quantifying local thermal nonequilibrium. Analysis of the slow cooling eigenmode identifies an effective thermal capacitance C_eff = 301 ± 23 J K⁻¹ and a parasitic loss conductance K_loss = 0.10–0.16 W K⁻¹; the capacitance is confirmed within 1% by fifteen independent estimates from the heating transients and is seven times that of the 40 g monolith, showing that the housing dominates the thermal inertia. All fifteen heating transients collapse onto master curves under a single time scale constructed from these constants. An energy closure demonstrates that the nominal aperture-power accounting underestimates delivered power by 30–60% at the two higher irradiances, and per-configuration calibration factors consistent with an independent model calibration are applied where power enters the results. A Monte Carlo propagation of instrument uncertainties accompanies all reported quantities.

**Keywords:** concentrated solar power; volumetric receiver; silicon carbide monolith; transient analysis; Nusselt correlation; local thermal nonequilibrium

---

## 1. Introduction

Concentrated solar technology can supply industrial process heat at temperatures above 1000 °C and is therefore a candidate replacement for fossil-fired heating in high-temperature processes [1]. Among receiver concepts for air as the heat-transfer fluid, volumetric receivers absorb concentrated radiation within a porous or channeled structure through which the air flows, with the aim of transferring heat to the gas over a distributed volume rather than at a single irradiated surface [2,3]. The intended operating condition, often called the volumetric effect, is one in which the irradiated front face remains cooler than the interior of the absorber, so that thermal emission losses from the aperture are reduced [4].

Experimental studies of structured absorbers have concentrated on steady-state performance, and the heat-transfer correlations proposed for monolithic and foam absorbers scatter widely: published Nusselt correlations of the form Nu = a Re^b report prefactors spanning two orders of magnitude and exponents from about 0.3 to 0.6 [5–8]. Transient measurements on structured monoliths are scarce, and the few available datasets have been used mainly to validate numerical models rather than to extract transport characteristics directly. As a result, several basic questions about this receiver class lack direct experimental answers: which resistance limits heat extraction; under what operating condition the volumetric effect appears; how far the gas departs from local thermal equilibrium with the solid; and which thermal mass governs the transient response of a receiver installed in its housing.

This work addresses these questions with a data-driven analysis of a transient experimental campaign on a laboratory-scale SiC honeycomb receiver. The approach deliberately avoids detailed simulation. Instead, the measurements are reduced to dimensionless groups, and the quantities of interest are obtained from algebraic energy balances, regression on derived groups, and the eigenstructure of the measured transients. The products are: (i) an apparent global Nusselt correlation and its interpretation through a regime analysis; (ii) a flux-independent effectiveness criterion for the volumetric inversion; (iii) a measured local-thermal-nonequilibrium relation; (iv) an effective capacitance and loss conductance identified from the cooling decay and verified against the heating transients; (v) a similarity collapse of all heating transients under a single time scale; and (vi) a delivered-power assessment based on energy closure. Each quantity is accompanied by a Monte Carlo uncertainty estimate. The complete dataset and the reduction scripts are publicly available.

The paper is organized as follows. Section 2 describes the experimental system and campaign. Section 3 presents the data reduction, the dimensionless framework, the identification methods, and the uncertainty analysis. Section 4 reports the results, Section 5 discusses their implications and limitations, and Section 6 summarizes the conclusions.

---

## 2. Experimental system

### 2.1 High-flux solar simulator

The high-flux solar simulator (HFSS) comprises seven 6 kWₑ xenon short-arc lamps with ellipsoidal aluminum-coated reflectors arranged in a compressed hexagonal pattern with the seventh lamp at the center. Each reflector has focal points at 90 mm and 2790 mm from the vertex, a 735 mm opening, and a 450 mm length; the lamp sits at the first focal point and the receiver aperture at the second. Details of the design and its characterization are given in [9].

Irradiance on the target plane was characterized by a three-step flux-mapping method [9]: a water-cooled Lambertian target (250 × 250 mm² stainless steel slab, plasma-coated with Al₂O₃) mounted on a three-axis precision slider provides the spatial distribution, and a Gardon-type circular-foil radiometer (Vatell TG1000-0, colloidal graphite coating, transducer calibration accuracy ±3%) placed at the target center provides the absolute scale. Three lamp configurations were used, giving nominal aperture irradiances of 456, 304, and 256 kW m⁻².

### 2.2 Receiver and instrumentation

The receiver is a commercial extruded honeycomb monolith of recrystallized silicon carbide with silica and secondary SiC-polytype additives. It is 137 mm long with a 19 × 19 mm² frontal cross-section and contains 100 square channels of 1.5 mm width separated by 0.4 mm walls. The measured monolith mass is 40 g, which corresponds to an effective solid density of approximately 2150 kg m⁻³, about one third below that of dense SiC; the extrudate is internally porous. The monolith is mounted horizontally inside an aluminum cylindrical cavity (142 mm diameter) and is surrounded by alumina felt insulation (thermal conductivity approximately 0.078 W m⁻¹ K⁻¹). Ambient air is drawn through the channels by a rotary-vane vacuum pump downstream of the flow controllers.

Temperatures were measured with type-N thermocouples. Three sheathed probes (1.5 mm OD) contact the external side wall of the monolith at 11, 58, and 107 mm from the irradiated face; these are denoted T8, T12, and T11 respectively. Two miniature probes (0.5 mm OD) are inserted into interior channels at 58 mm (T9) and 107 mm (T10) and are exposed to the channel flow. The gas outlet temperature T3 is measured at approximately 136 mm, at the receiver exit and off the optical axis, so that the probe has no direct view of the radiation source. A further probe (T2) is embedded 40 mm deep in the insulation, and two probes (T15, T16) record ambient temperature. It is noted that an earlier description of this apparatus interchanged the wall and interior probe designations; the assignment stated here is the correct one and is used throughout.

The air flow is set by four mass-flow controllers (Aalborg GFC17, 0–5 standard L min⁻¹ each, accuracy ±1% of full scale), operated in parallel to improve resolution; all four channels carry air, and the total standard volumetric flow is their sum. The manufacturer reference conditions are 21.1 °C and 101.325 kPa. Pressure drop across the monolith is monitored with a differential transducer (Keller PD33X, ±200 mbar, ±0.1% FS). All signals are logged at 1 Hz.

### 2.3 Campaign

Fifteen heating experiments cover the three irradiance levels at five flow rates each, between 4.5 and 18.3 standard L min⁻¹ (Table 1). Each experiment starts from thermal equilibrium with the laboratory and continues until the temperature drift is below 0.5 K min⁻¹ over the final ten minutes; run durations are between 3000 and 7150 s. Three cooling experiments follow selected heating runs: the lamps are switched off at time zero and the air flow is maintained (10.5, 6.6, and 4.5 L min⁻¹). One replicate heating run at the low-flux condition (13.9 L min⁻¹) agrees with its counterpart within 2 K on every sensor, which is within the thermocouple class accuracy; a second replicate was excluded after a 25–30 K discrepancy attributed to a change in lamp condition. Per-channel Reynolds numbers span 23 to 94.

**Table 1.** Campaign summary. G₀ is the nominal aperture irradiance; q the total standard volumetric flow; Re the per-channel Reynolds number at mean gas temperature (Section 3.2).

| G₀ [kW m⁻²] | q [sL min⁻¹] | Re | runs |
|---:|---|---|---|
| 456 | 7.1–15.3 | 34–73 | 5 heating, 1 cooling (10.5) |
| 304 | 4.5–18.3 | 23–94 | 5 heating |
| 256 | 4.5–13.9 | 25–79 | 5 heating + 1 replicate, 2 cooling (6.6, 4.5) |

---

## 3. Data reduction and analysis methods

### 3.1 Reduction

Temperatures are converted to kelvin and the ambient reference is T_amb = (T15 + T16)/2. The mass flow is ṁ = ρ_std q/60000 with ρ_std = 1.199 kg m⁻³ evaluated at the controller reference conditions. Steady-state values are means over the final 120 s of each run; the residual drift over this window is below 0.5 K and is carried as an uncertainty component. Air properties μ(T), k(T), and c_p(T) are interpolated from standard tabulations between 300 and 1200 K, and gas enthalpy differences are computed as ∫c_p dT rather than c_pΔT.

The wall-temperature profile along the receiver is represented by the three side-wall probes. An energy-representative wall temperature is formed by trapezoidal quadrature of the piecewise-linear profile over the receiver length,
$$\bar T_w = 0.248\,T_8 + 0.365\,T_{12} + 0.387\,T_{11},$$
where the weights are the exact trapezoid coefficients for the probe positions with constant extrapolation to the end faces. The interior probes T9 and T10 are used only for the nonequilibrium metrics of Section 3.2 and enter no other quantity.

### 3.2 Dimensionless framework

Per-channel quantities use ṁ_ch = ṁ/100, channel width w = 1.5 mm, hydraulic diameter D_h = w (square section), wetted perimeter P = 4w, flow area A_ch = w², and heated length L = 137 mm. The evaluation temperature is stated for each group because gas properties vary along the receiver.

The Reynolds number is Re = ṁ_ch D_h/(A_ch μ(T̄_g)) with T̄_g = (T_amb + T3)/2. Since ṁ_ch is constant along the channel while μ increases with temperature, Re decreases from inlet to outlet by roughly a factor of two; the reported value at T̄_g lies between the two.

The Prandtl number Pr = c_p μ/k is computed at T̄_g for every run. It spans 0.684 to 0.690 over the campaign, a variation of ±0.4%, for two reasons. First, the viscosity and thermal conductivity of a gas share nearly the same temperature dependence, so their ratio is almost independent of temperature; the residual dependence enters through c_p and is weak (by the Eucken relation, Pr ≈ c_p/(c_p + 1.25R/M) changes by about 2% while c_p itself changes by 17% between 300 and 1200 K). Second, the mean-gas reference temperatures of all runs fall within 410–530 K. The Prandtl number is therefore constant in fact although the properties vary by design, and the Nusselt correlation of Section 4.3 is reported as Nu(Re) at fixed Pr.

The Graetz number at the outlet, Gz_L = D_h Re Pr/L, measures the residual thermal development.

The receiver is treated as a single-stream heat exchanger between the solid matrix and the gas. The heat-exchange effectiveness and number of transfer units are

$$\varepsilon = \frac{T_3 - T_{\rm amb}}{\bar T_w - T_{\rm amb}},\qquad NTU = -\ln(1-\varepsilon),$$

the second relation being exact for a single stream exchanging with a locally fixed wall temperature. The effectiveness is a purely measured quantity; the logarithm amplifies its uncertainty as ε approaches unity, which the uncertainty propagation of Section 3.6 accounts for.

The apparent global heat-transfer coefficient and Nusselt number are

$$h = \frac{NTU\,\dot m_{\rm ch} c_p(\bar T_g)}{P L},\qquad Nu = \frac{h D_h}{k(T_{\rm film})},\qquad T_{\rm film} = \tfrac12(\bar T_w + \bar T_g).$$

The qualifier *apparent* is used deliberately: h is defined with respect to the instrumented side-wall temperature and the full nominal exchange area, and therefore lumps the channel-film resistance together with cross-sectional solid nonuniformity, flow maldistribution among channels, and any peripheral bypass. It is the coefficient accessible to external measurement, and its relation to channel-film theory is examined in Sections 4.3 and 5.1.

The wall Biot number is Bi = h t_w/(2k_SiC) with wall thickness t_w = 0.4 mm.

Local thermal nonequilibrium is quantified by the normalized interior-probe deficit at the two instrumented stations,

$$\Lambda_{58} = \frac{T_{12} - T_9}{T_{12} - T_{\rm amb}},\qquad \Lambda_{107} = \frac{T_{11} - T_{10}}{T_{11} - T_{\rm amb}}.$$

Because the interior probes exchange radiation with the surrounding channel walls, any radiative bias raises their readings toward the solid temperature; the measured Λ values are therefore lower bounds on the gas-side deficit.

The operating coordinate of the campaign is the power ratio PR = P_del/ṁ in kJ kg⁻¹, where P_del is the delivered power defined in Section 3.5.

### 3.3 Regime classification

Before any correlation is interpreted, the dimensionless envelope is used to establish which transport mechanisms can be rate-limiting: laminarity and entrance effects (Re, Gz), axial gas conduction (Péclet number Pe = Re Pr), buoyancy (Gr/Re²), wall-normal solid resistance (Bi), gas rarefaction (Knudsen number), and radiative exchange within the channels (radiation-conduction number N_rc = 4σT³D_h/k). The resulting classification is reported in Section 4.1.

### 3.4 Transient identification

After the fast internal modes of the assembly decay, the energy balance reduces to a single node,

$$C_{\rm eff}\,\frac{d\theta}{dt} = -\left(\varepsilon\,\dot m c_p + K_{\rm loss}\right)\theta,\qquad \theta = \bar T - T_{\rm amb},$$

whose observable is the slow eigenvalue

$$\lambda = \frac{\varepsilon\,\dot m c_p + K_{\rm loss}}{C_{\rm eff}}.$$

Two independent estimates of λ are available at each operating point. In the cooling experiments, λ is the log-slope of (T − T_amb) over the late half of the record (excess temperature above 5 K), evaluated for each of the six receiver sensors; the spread among sensors, below 5% in every run, verifies that a single mode governs the decay before the model is used. In the heating experiments, linearization about the steady state gives T_ss − T(t) ∝ e^(−λt) with the same eigenvalue; the deficit fraction u = (T_ss − T)/(T_ss − T₀) is fitted on the interval 0.07 < u < 0.45, excluding early times so that the fast modes have decayed and excluding the tail where the logarithm amplifies noise, and sensors with r² < 0.95 are discarded.

The eigenvalues are regressed against the gas advective conductance ε ṁc_p, with ε evaluated at each flow from the steady correlation ε(q); the slope and intercept identify C_eff and K_loss without any discretized model or optimizer. The cooling and heating datasets provide two independent identifications whose agreement serves as internal validation.

Early-window single-exponential fits to the wall probes give position-dependent time constants τ(z), from which pairwise apparent axial diffusivities α_app = Δz²/Δτ are formed.

Finally, the heating transients are tested for similarity. Each record is normalized as θ*(t) = (T − T₀)/(T_ss − T₀) and time is rescaled as t* = t(ṁc_p + K_loss)/C_eff using only constants identified from the cooling data. The quality of the collapse is quantified by the coefficient of variation of t* at θ* = 0.5 across the fifteen runs.

### 3.5 Delivered power

The nominal power accounting, P_nom = G₀A_frt with A_frt = 3.61 cm², is shown in Section 4.6 to underestimate the power reaching the receiver system: the measured gas enthalpy output alone exceeds P_nom at the highest flow. Where power enters a reported quantity explicitly, namely in η and PR, delivered-power calibration factors are applied per lamp configuration:

$$P_{\rm del} = f\,G_0 A_{\rm frt},\qquad f_{456} = 1.336,\quad f_{304} = 1.374,\quad f_{256} = 0.786.$$

These factors originate from an independent calibration of a one-dimensional transport model of the same receiver in which per-irradiance power scales were the only free parameters. They are adopted here only because a model-free steady energy closure on the present data, f = [Q_gas + K_loss(T̄_w − T_amb)]/(G₀A_frt), brackets the same values for every group (Section 4.6). Their uncertainty is taken as ±8%, the per-case scatter of the calibration, and treated as systematic per configuration. Three provisions protect the analysis from circularity: the factors affect no temperature-based result (ε, Nu, Λ, ε*, C_eff, K_loss, and the master curves are factor-free); only their per-configuration grouping, which the closure confirms independently, is relied upon; and all nominal-basis values are retained in the published data products so that a future flux-map measurement can replace the factors without repeating the reduction.

### 3.6 Uncertainty analysis

A Monte Carlo propagation with 4000 realizations treats instrument errors as systematic within a run and independent between instruments. Each thermocouple carries a bias drawn from a normal distribution with standard deviation 1.1 K (class accuracy ±2.2 K) plus a steady-window drift of standard deviation 0.5 K; each of the four flow controllers carries a bias of standard deviation 0.025 sL min⁻¹; property tabulations are perturbed by ±2% (μ, k) and ±1% (c_p), correlated across runs; the radiometer calibration by ±3%; and the delivered-power factors by ±8%. For every realization the complete analysis chain is recomputed, including all per-run groups, the Nusselt-law fit, the three inversion crossings, and the eigenvalue regressions with each λ additionally perturbed by its fit standard error. Nonlinear amplifications, in particular dNTU = dε/(1 − ε), are thereby captured without linearization. Reported uncertainties are Monte Carlo standard deviations, and parameter intervals are 95% percentile ranges.

---

### 3.7 Continuum model used for comparison

The model used in Section 4.9 is a two-zone quasi-one-dimensional formulation of the installed assembly. A finite-volume energy balance is solved along the receiver axis for a participating core solid and for a perimeter zone that lumps the monolith periphery, the alumina-silicate felt and the aluminium casing; the two are coupled by a linear radial conductance and each carries its own axial conduction. The gas is advanced quasi-steadily through the core, T_{g,i+1} = T_{g,i} + (1 − e^{−NTU_i})(T_{s,i} − T_{g,i}), with NTU_i built from a Reynolds-dependent effective Nusselt number; gas storage is neglected because the channel residence time is several orders of magnitude below the ceramic response time. The aperture flux is deposited in the core with Beer–Lambert attenuation, and the fraction of the spillage that is captured is deposited on the perimeter front face, since neither the dense felt nor the casing wall admits axial penetration. A lumped cavity node and a rear adaptor/tube/water-cooled-flange path close the assembly. The active flow fraction is allowed to depend on the core temperature, which is the only mechanism in the formulation capable of satisfying the heating and cooling phases simultaneously (Section 5.6).

The parameters are identified in a single simultaneous calibration against all fifteen heating and three cooling runs; heating and cooling are never fitted separately and no phase switch is used. Two properties of this calibration must be stated because they bound what the comparison can support. First, the participating heat capacity enters through bounds derived from the value identified in Section 4.7, so the model's capacitance is imposed and its agreement with that value is not an independent confirmation. Second, the absorbed power is not identifiable from the temperature data alone, for the reason given in Section 5.6. Every model quantity compared in Section 4.9 is reduced from the model output using the identical definitions of Section 3.2, including the same wall quadrature and the same exclusion of the interior probes.

## 4. Results

### 4.1 Operating envelope and transport regimes

The campaign spans Re = 23–94, Pr = 0.684–0.690, Gz_L = 0.18–0.70, ε = 0.57–0.78, NTU = 0.85–1.51, apparent Nu = 0.028–0.212, and PR = 262–1665 kJ kg⁻¹ on the delivered-power basis. The per-run values and their uncertainties are tabulated in the supplementary data.

These ranges fix the transport regime unambiguously. The flow is deeply laminar throughout, two orders of magnitude below transition. The hydrodynamic and thermal entrance lengths, L_h ≈ 0.05 Re D_h ≤ 7 mm and L_t ≈ 0.05 Re Pr D_h ≤ 5 mm, occupy at most the first few percent of the channel length, and Gz_L < 0.7 confirms that the flow is thermally developed at the exit; fully developed laminar duct theory, with the constant-Nusselt limit Nu ≈ 3 for a square duct, is therefore the appropriate single-channel reference. Axial conduction in the gas is negligible, since Pe(L/D_h) = 1500–5900. Buoyancy is irrelevant inside the channels despite the horizontal orientation: Gr ≈ 0.4–0.7 at film conditions, so Gr/Re² ≤ 10⁻³. The channel walls are radially lumped, with Bi < 10⁻⁵, so any solid-side nonuniformity must occur at the scale of the monolith cross-section rather than the wall thickness. The flow is continuum, with Kn ≈ 10⁻⁴. One mechanism is not negligible: the radiation-conduction number N_rc = 4σT³D_h/k reaches 1.6 at 600 K and about 5 at 1000 K, so radiative exchange between axial wall locations within a channel is comparable to or larger than gas conduction, and radiation is available to redistribute heat along the solid without gas participation.

### 4.2 Steady temperature maps

Figure 2 shows the steady wall and gas temperatures against flow for the three irradiances. At fixed irradiance, the front wall temperature T8 decreases steeply with flow (−33.2, −24.3, and −20.8 K per sL min⁻¹ at 456, 304, and 256 kW m⁻², r² ≥ 0.97), the mid-length wall T12 at about half that rate, and the rear wall T11 hardly at all (−1.1 to −5.2 K per sL min⁻¹, r² ≤ 0.76); the outlet gas temperature is statistically independent of flow at the two higher irradiances. The flow sensitivity of the wall temperature therefore decays monotonically with depth and is essentially extinguished by 107 mm. The insulation probe rises only 10–60 K above ambient while the adjacent wall runs 300–770 K above ambient, indicating that the radial heat path is limited by the insulation, and the negligible decay of the insulation temperature during the cooling runs shows that the felt is still charging when the heating runs end.

The replicate pair at the low-flux condition agrees within 2 K on every sensor, and within 2% in ε and Nu. Run-to-run repeatability therefore does not limit any conclusion that follows.

### 4.3 Apparent global exchange correlation

Across all fifteen heating runs the apparent Nusselt number follows a single power law (Figure 4, left),

$$Nu = (3.1 \pm 0.1)\times 10^{-4}\,Re^{1.44},\qquad r^2 = 0.97,\qquad 23 \le Re \le 94,\ Pr = 0.69.$$

The uncertainty of the exponent separates into an instrument contribution of ±0.004 (Monte Carlo) and a run-to-run scatter contribution of ±0.069 (regression); the correlation is limited by real operating-point variability rather than by measurement noise.

Two features of this result carry the physical content. The first is the magnitude. The measured Nu of 0.03–0.21 lies a factor of 15–100 below the fully developed laminar value. If the channel film were the controlling resistance, the corresponding NTU would be approximately 45 and the gas would reach the wall temperature within a few millimeters of entry; the measured NTU is close to 1 and the outlet gas leaves 150–300 K below the wall. Heat extraction in this receiver is therefore not limited by the channel boundary layer. Since the regime analysis excludes wall-normal conduction, entrance effects, buoyancy, and rarefaction, the limitation must arise at the scale of the assembly (Section 5.1). The second feature is the exponent. A value of 1.44 is well above the range 0.3–0.6 characteristic of film-controlled laminar correlations. An exponent above unity means that the apparent conductance grows faster than the flow rate itself, which cannot be produced by a boundary-layer coefficient and instead indicates that the fraction of the solid participating effectively in exchange increases with flow.

### 4.4 Effectiveness criterion for the volumetric inversion

At low flow the wall-temperature maximum sits at the front face (T12 − T8 < 0); at high flow it moves into the receiver (T12 − T8 > 0). The crossing flows differ strongly between irradiance groups (11.1, 10.3, and 3.7 sL min⁻¹) and do not collapse on the Reynolds number (Re* = 52, 52, and 20). They collapse on the effectiveness (Figure 3):

$$\varepsilon^* = 0.671 \pm 0.003,\quad 0.671 \pm 0.003,\quad 0.626 \pm 0.005\qquad (456,\ 304,\ 256\ \mathrm{kW\,m^{-2}}),$$

with 95% intervals of ±0.005–0.010. The volumetric inversion therefore occurs when the gas stream recovers approximately two thirds of the wall excess temperature, independent of irradiance over the 1.8-fold flux range tested.

The mechanism follows from the front-face energy balance. The face loses heat by thermal emission and by convective quenching from the entering gas; both the absorbed flux and the quench term scale essentially linearly with irradiance at fixed geometry, so the balance that determines the location of the temperature maximum reduces to a condition on the gas-side recovery, that is, on ε. The small but statistically significant reduction of ε* at the lowest flux (0.626 against 0.671, nine times its uncertainty interval) is consistent with the growing relative weight of the flux-independent parasitic loss identified in Section 4.7: a fixed loss removes proportionally more from the face balance when the input is small, allowing the inversion at a slightly lower recovery. On the delivered-power coordinate the crossings lie at PR* = 993, 731, and 985 kJ kg⁻¹, a 1.4-fold spread compared with 2.4-fold on the nominal basis; the delivered-power correction improves the power-ratio coordinate but no power ratio collapses the inversion as tightly as the effectiveness.

As a design statement: a receiver of this class operates volumetrically when its operating point satisfies ε > ε* ≈ 2/3, a condition on the heat-exchange operating point rather than on the optics.

### 4.5 Local thermal nonequilibrium

The interior probes read below their wall counterparts at the same axial station in every run. At 58 mm the normalized deficit is small and independent of flow, Λ₅₈ = 0.03–0.06 (21–27 K). At 107 mm it grows linearly with the Reynolds number (Figure 4, center),

$$\Lambda_{107} = 0.038 + 8.3\times10^{-4}\,Re,\qquad r^2 = 0.90,$$

reaching 0.11 (55 K) at Re = 94; within individual irradiance groups the linearity is nearly exact (r² up to 0.998), and the per-point uncertainty of ±0.003–0.006 is an order of magnitude smaller than the flow trend. Because radiative exchange with the surrounding channel walls biases the interior probes toward the solid temperature, the true gas deficit is at least as large as measured, and the conclusion is conservative: deep in the receiver the gas has not equilibrated with the solid, and the shortfall grows in proportion to velocity. Equivalently, the thermal equilibration length exceeds the receiver length at every operating point, and increasingly so at high flow. This is the local counterpart, measured by an independent set of sensors, of the global NTU ≈ 1 found in Section 4.3.

### 4.6 Energy accounting and delivered power

The gas enthalpy output Q_gas = ṁ∫c_p dT (from T_amb to T3) is a bias-free measurement: the outlet probe sits at the receiver exit, off the optical axis, without a direct view of the source. Compared with the nominal input P_nom = G₀A_frt, the apparent efficiency rises nearly linearly with flow within each irradiance group and reaches 1.23 ± 0.04 at the highest flow (304 kW m⁻², 18.3 sL min⁻¹). The Monte Carlo interval includes the radiometer calibration, all thermocouple biases, and the metering uncertainty, and alternative defensible metering conventions raise rather than lower the value. An output exceeding the nominal input by five standard deviations admits one conclusion: the power reaching the receiver system exceeds G₀A_frt. The route is geometric, since the concentrated beam overfills the 3.61 cm² frontal face and spillage absorbed by the cavity front and insulation reaches the gas through the cavity air path and the heated rear hardware. The nominal accounting must therefore be reported as a lower bound on delivered power, with a shortfall of at least 23 ± 4% at 304 kW m⁻².

Combining the two quantities identified in this work, Q_gas and K_loss, gives a per-run estimate of the delivered power itself, f = [Q_gas + K_loss(T̄_w − T_amb)]/(G₀A_frt). No thermal model enters this closure. Group means, with the K_loss bracket carried as a systematic band, are:

| G₀ [kW m⁻²] | f, closure (K lower) | f, closure (K upper) |
|---:|---:|---:|
| 456 | 1.05 | 1.34 |
| 304 | 1.23 | 1.58 |
| 256 | 0.84 | 1.11 |

Two statements survive the width of the band: the 456 and 304 kW m⁻² configurations deliver of order 30–60% more power than the nominal accounting, and the 256 kW m⁻² configuration delivers approximately the nominal value or less. The residual flow dependence of the estimate (slopes of +0.03 to +0.06 per sL min⁻¹) shows that it is not a purely optical quantity; part of the factor compensates flow-dependent physics outside the accounting. The per-configuration grouping, not the third digit, is the reliable content. This closure is the only route to the delivered power available from these measurements: as shown in Section 5.6, the source magnitude cannot be separated from the captured spillage by inverse modelling, so a model-calibrated delivery factor would restate the present result rather than corroborate it. Different lamp configurations at the three flux levels make a per-configuration delivery error physically plausible, and the near-nominal value at 256 kW m⁻² is consistent with the reduced ε* of that group (Section 4.4).

With the calibration factors of Section 3.5 applied, every efficiency falls below unity (Figure 3, right):

$$\eta_{\rm del} = 0.30\text{–}0.67,\quad 0.21\text{–}0.90,\quad 0.32\text{–}0.87\qquad (456,\ 304,\ 256\ \mathrm{kW\,m^{-2}}),$$

with uncertainties of ±0.06–0.08 dominated by the factor uncertainty. The efficiency rises linearly with flow in every group (slopes 0.045–0.059 per sL min⁻¹, r² ≥ 0.98) and does not saturate within the envelope; the maximum demonstrated recovery is about 90% of delivered power at 18.3 sL min⁻¹. The receiver is flow-starved rather than exchange-saturated at every operating point tested.

### 4.7 Transient identification and cross-validation

In each cooling run the late-time log-slopes of the six receiver sensors agree within 5% (λ = 7.99, 6.38, and 4.84 ×10⁻⁴ s⁻¹ at 10.5, 6.6, and 4.5 sL min⁻¹): after internal redistribution the assembly decays as a single thermal mode, which verifies the premise of the identification before it is applied. Regression of λ against ε ṁc_p yields, from the cooling data alone,

$$C_{\rm eff} = 301 \pm 23\ \mathrm{J\,K^{-1}},\qquad K_{\rm loss} = 0.096 \pm 0.011\ \mathrm{W\,K^{-1}}\qquad (r^2 = 0.96).$$

The heating experiments provide fifteen further eigenvalues through the late approach to steady state, spanning a wider flow range. Regression of these gives C_eff = 305 ± 32 J K⁻¹ and K_loss = 0.165 ± 0.028 W K⁻¹ (r² = 0.79). The two datasets share no measurements, one being recorded with the lamps on and the other with the lamps off, yet they return the same capacitance within 1%. The identification is therefore confirmed by eighteen eigenvalue measurements in total, and no additional experiments are required to establish C_eff.

The loss conductance differs between the two datasets, and the difference is informative rather than contradictory: K_loss linearizes a loss path that is partly radiative and hence proportional to T⁴, and the heating approaches sample the neighborhood of the hot steady states while the cooling windows lie 150–250 K lower. The ratio of 1.7 over that interval is consistent with the T³ sensitivity of a radiative surface loss. The two values bracket the loss law over the operating range: K(T) ≈ 0.10–0.16 W K⁻¹ between roughly 550 and 1050 K.

The measured monolith mass of 40 g, with c_p of SiC between 1050 and 1170 J kg⁻¹ K⁻¹ at 600–900 K, gives a monolith capacitance of 42–47 J K⁻¹. The identified C_eff is 6.4–7.2 times larger: the thermal mass participating in the slow mode is dominated by the insulation and mounting hardware coupled to the monolith, not by the ceramic itself. This ratio is obtainable only from transient data, and it quantifies how far the installed receiver departs from the bare-absorber idealization. The measured mass also shows that property sets based on the density of dense SiC would overestimate the monolith heat capacity by about 50%.

Early-window time constants along the wall increase front to back (τ = 226, 307, and 545 s at 10.5 sL min⁻¹; 593, 845, and 1440 s at 4.5 sL min⁻¹). The pairwise apparent axial diffusivities, α_app = Δz²/Δτ = 0.4–2.7 ×10⁻⁵ m² s⁻¹, increase with flow, indicating that axial heat redistribution in the assembly is partly gas-mediated rather than purely conductive. The front decays fastest at every flow although it carries no preferential flow path; its early transient is dominated by radiative loss through the aperture, the same T⁴ path that produces the temperature dependence of K_loss.

### 4.8 Similarity of the heating transients

Rescaling each heating run by the single time scale identified from the cooling data, t* = t(ṁc_p + K_loss)/C_eff, with no parameter fitted to any heating transient, collapses the fifteen normalized histories (Figure 6). The half-rise points are t* = 0.20 with a coefficient of variation of 13% for the energy-weighted wall temperature, and t* = 0.62 with 9% for the outlet gas, across a 3.4-fold range of flow and a 1.8-fold range of flux.

Two conclusions follow. First, for system-level purposes the receiver exhibits one-parameter dynamics: the ratio C_eff/(ṁc_p + K_loss), measurable from a few decay records, predicts the transient shape of every heating run to within about 10% regardless of power level. Second, the factor of three between the wall and gas half-rise times measures the wall-to-gas thermal lag. The outlet gas is not a proxy for the state of the solid during transients, and any reduced-order description that couples them rigidly will misplace one or the other by this factor. Both numbers are directly applicable to control design and ramp-rate specification.

---

### 4.9 Comparison with a continuum model

The calibrated model of Section 3.7 reproduces the measured effectiveness envelope run by run (Figure 7a). Over the fifteen heating runs the mean absolute deviation in ε is 0.018 and the maximum is 0.038, with a bias of +0.005; the modelled range 0.583–0.776 covers the measured 0.572–0.779, and NTU agrees to a mean absolute deviation of 0.055. Mean steady residuals are +14 K at 11 mm, −21 and −40 K at 58 mm, +16 and +29 K at 107 mm, +20 K on the insulation and −0.1 K on the outlet gas. At the level of the operating point, therefore, the assembly-scale exchange embodied in ε and NTU is captured.

The volumetric inversion is not (Figure 7b). Reducing the model output by the same convention used in Section 4.4 — regression of the indicator I_vol = T12 − T8 on flow to locate the crossing, then evaluation of the effectiveness regression there — the model inverts at ε* = 0.786, 0.764 and 0.726 for the 256, 304 and 456 kW m⁻² groups, against the measured 0.626, 0.671 and 0.671. It therefore requires 0.06–0.16 more gas-side recovery to move the wall peak inside than the receiver does, equivalently 1.2 to 3.6 times more flow (q* = 13.7, 16.3 and 13.5 against 3.8, 10.4 and 11.0 sL min⁻¹). The residual flux dependence is not merely retained but reversed: the measured threshold rises weakly with irradiance while the modelled one falls. The amplitude fails as well. Where the measured indicator reaches +58, +60 and +59 K at the top of each group, the model reaches +0.7, +5.5 and −9.5 K; it never develops an interior wall maximum of more than a few kelvin, and in the lowest-flux group it does not develop one at all.

The two observations are consistent with a single defect, and the residuals localize it. The energy-weighted wall mean is reproduced — which is why ε is right — while the front-to-mid difference is not: the +14 K residual at 11 mm against the −40 K residual at 58 mm is a 54 K error in the axial gradient, and the mean offset of the indicator across the fifteen runs is −54 K. The discrepancy is therefore a redistribution error rather than an error in the magnitude of gas–solid coupling: at the correct effectiveness the model places too much of its wall temperature at the front face. This is a statement about the axial placement of the volumetric source relative to the axial distribution of exchange, and it is the sharpest available diagnosis of what a continuum representation of this receiver is missing.

## 5. Discussion

### 5.1 Interpretation of the apparent Nusselt correlation

The regime analysis of Section 4.1 eliminates the channel-scale explanations for the gap between the apparent Nusselt number and duct theory one by one: wall-normal conduction (Bi < 10⁻⁵), entrance limitation (Gz_L < 0.7), buoyancy (Gr/Re² ≤ 10⁻³), axial gas conduction (Pe L/D_h > 10³), and rarefaction (Kn ≈ 10⁻⁴). What remains is the assembly scale, with two candidate mechanisms that the data support directly. The first is cross-sectional solid nonuniformity: the instrumented side wall is backed by insulation and exposed to cavity radiation, and therefore runs hotter than the average channel wall that the gas contacts; the interior probes reading 21–55 K below the side wall provide independent evidence. The second is maldistribution: flow non-uniformity among the hundred channels, and possibly a peripheral bypass, lowers the mixed outlet temperature at a given wall temperature. Both mechanisms reduce the apparent Nu without any violation of channel physics.

On this reading the super-linear exponent acquires a definite meaning. Increasing the flow deepens and extends the cold interior region, so the fraction of the solid that participates effectively in exchange grows with Re; this multiplies the film scaling and drives the apparent conductance faster than linearly. The wide scatter of published correlations for volumetric absorbers is consistent with this interpretation, since an assembly-level coefficient necessarily reflects assembly-specific nonuniformity, and such correlations should not be expected to transfer between geometries.

The implication for modeling is structural rather than parametric. A single-channel conjugate model equipped with a textbook Nusselt law cannot reproduce these data at any parameter value, because the controlling resistance is not in the channel. Model closures for this receiver class require either a resolved solid-temperature distribution across the section or an explicit flow-dependent active fraction.

### 5.2 Scope of the effectiveness criterion

The strength of the ε* criterion is its economy: one measured operating-point quantity, requiring one gas and three wall thermocouples, predicts whether the receiver operates volumetrically, across all tested fluxes, with a threshold resolved to ±0.005. Its scope is bounded by the dataset: one geometry, laminar flow, and an effectiveness referenced to the side-wall probes. The reduction of ε* at the lowest flux is resolved and is consistent with the fixed parasitic loss claiming a larger share of the front-face balance at low power; a two-term criterion could capture this, but three flux levels are insufficient to fix the additional coefficient. Testing the threshold on other structured absorbers, where inversion data exist but are usually reported against flow rather than effectiveness, is the natural external validation.

### 5.3 The nonequilibrium relation

Λ₁₀₇(Re) is a directly measured counterpart of the interphase nonequilibrium that two-temperature porous-medium models represent, and the radiative bias of the interior probes acts in the conservative direction. Its linearity in Re with a near-zero intercept indicates that the deep-receiver gas deficit is advection-controlled and vanishes only in the static limit; there is no flow regime within the envelope in which local thermal equilibrium holds at 107 mm. Model formulations that assume local equilibrium are outside their domain of validity for this receiver class at Re above roughly 20, and calibrated two-temperature closures should reproduce both the slope and the flux-invariance of the measured relation.

### 5.4 The transient identification as a measurement

The eigenvalue method is deliberately sparing in assumptions: it requires only that a slowest mode exist and be observable, which the 5% agreement among sensors verifies empirically. Its products are portable constraints on any future model of this receiver, of whatever fidelity. The agreement of the capacitance between disjoint datasets, together with the temperature dependence residing in the loss coefficient, could not have been obtained from either dataset alone. The practical consequence of C_eff ≈ 7 C_monolith is that start-up and cloud-transient response of the installed receiver are governed by the housing and insulation rather than by the absorber, so improvements in dynamic performance should address the mounting and insulation architecture rather than the ceramic.

### 5.5 Delivered power

The result that the measured gas output exceeds the nominal input is an uncomfortable one to report, and reporting it as a bound rather than absorbing it silently into a calibration factor keeps the dataset useful: every efficiency herein is stated with its basis, and every absorbed-power figure is a minimum. The convergence of the algebraic closure and the independent model calibration on the same per-configuration grouping justifies the adopted factors, while their residual flow dependence warns against interpreting them as purely optical constants. For future campaigns at this scale, the flux map should be integrated over the actual receiver face and the spillage fraction quantified, since spillage comparable to the face area is unavoidable in small-aperture cavity configurations.

### 5.6 What the continuum comparison establishes

Four statements follow from Section 4.9 and from the wider set of formulations examined, none of which requires any fitted transport coefficient to be correct.

First, the agreement in effectiveness and the failure in inversion are not in tension. The effectiveness is a ratio of temperature differences and is therefore insensitive both to the absolute source magnitude and to the axial distribution of wall temperature; the inversion indicator depends only on that distribution. A formulation can reproduce the assembly-scale operating point and still misplace the source–exchange balance along the axis, and this one does.

Second, the class of internal flow structures compatible with the data is narrow. A two-stream formulation, in which core and perimeter channels carry separately tracked gas streams that mix downstream, is rejected by the calibration: the objective worsens by a factor of three and the optimizer routes the entire flow through the core rather than use the perimeter path. A scalar, operating-point-independent bypass fraction cannot satisfy the heating and cooling phases simultaneously, because a bypass large enough to keep the core hot in steady heating removes the advective capacity required to cool it at the observed rate. A bypass that depends on the core temperature does satisfy both. But a temperature-driven mechanism necessarily inherits a flux dependence, since the core temperature scales with irradiance, and this is precisely the property the measured threshold does not have. The flux independence reported in Section 4.4 is therefore evidence that the partition governing the inversion is geometric or Reynolds-driven rather than thermally driven — a stronger statement about the mechanism than the fit quality of any single formulation.

Third, the limitation is not removed by adding a spatial dimension. An axisymmetric continuum treatment of the same assembly, supplied with the extreme flow maldistribution and peripheral heating that the one-dimensional analysis implies, reproduces the sign of the radial and axial inversions but leaves side-wall errors of 94–116 K in heating and 44–55 K in cooling in its least unfavourable configuration, and about 210 K when the macroscopic maldistribution parameters are selected by grid search, while satisfying every conservation check to a relative enthalpy residual below 10⁻¹⁰. The failure is structural rather than numerical, and it is consistent with the assembly-scale reading of the apparent Nusselt correlation in Section 5.1: what is missing from both formulations is a description of which solid participates, not a better coefficient for exchange with the solid that does.

Fourth, two quantities are not identifiable from this instrumentation, and this bounds what any model of this class can be asked to deliver. The absorbed power is one of them. Because the aperture subtends only about 7% of the beam, the power reaching the assembly is the product of a source scale and the captured fraction of the remaining 93%, and the temperature field responds only to that product; the two factors are individually unconstrained. The delivered-power factors of Section 4.6, which are obtained algebraically from gas enthalpy and the identified loss conductance with no thermal model in the loop, are therefore the only defensible route to that quantity, and any model-derived value of it should be read as a restatement of them. The receiver conductance UA(Re) is the second: no unique identification exists without an independent bulk outlet enthalpy or a verified transfer function for the outlet probe, and admissible parameter sets form a ridge rather than a basin. Both limits are properties of the measurement, not of the optimizer, and they carry a direct design consequence for campaigns of this type: the flux map should be integrated over the receiver face with the spillage fraction quantified separately, and a bulk outlet enthalpy should be measured independently of any single probe.

### 5.7 Limitations

The transient identification rests on three cooling runs, but is confirmed by fifteen independent eigenvalues from the heating data; the remaining open quantity is the shape of K(T) between its two bracketing values. The effectiveness and Nusselt number are referenced to the side-wall probes; referencing the interior probes instead shifts ε by +0.10–0.15 and Nu by about a factor of two without changing any qualitative conclusion, and the choice is stated and reproducible. The radiative bias of the interior probes is signed but not quantified, and only the sign is used. All conclusions are confined to Re = 23–94, Pr ≈ 0.69, one monolith geometry, and one cavity configuration; the ε* threshold and the Nusselt prefactor in particular should be transferred to other geometries as hypotheses to be tested rather than as constants.

The model comparison of Sections 4.9 and 5.6 carries its own limitations, which are stated here rather than qualified in place. The transport parameters of that model are effective quantities identified by inverse fitting against this dataset; they are used to test structural hypotheses and are not reported, and should not be read, as measured material or heat-transfer properties. The participating heat capacity enters the calibration through bounds taken from Section 4.7, so the model's agreement with that value is imposed and is not an independent confirmation of it. The absorbed power is degenerate in the sense described in Section 5.6, and the model is therefore silent on the delivered-power question, which is settled algebraically in Section 4.6. Finally, the negative results in Section 5.6 are conditional on the class of formulations examined — two-zone and axisymmetric continuum representations with scalar internal closures — and establish that this class cannot reproduce the measured inversion, not that no continuum description can.

---

## 6. Conclusions

1. Over 23 ≤ Re ≤ 94 the apparent global exchange of the receiver follows Nu = 3.1×10⁻⁴ Re^1.44 (r² = 0.97), a factor of 15–100 below fully developed duct theory, with an exponent above unity. Heat extraction is limited at the assembly scale, and the exponent measures the flow-dependent recruitment of solid participation rather than a boundary-layer property.
2. The volumetric inversion obeys a flux-independent operating-point criterion: the wall-temperature peak moves into the receiver when the gas-side effectiveness exceeds ε* ≈ 0.66 ± 0.03, resolved to ±0.005 at each flux level.
3. Gas–solid nonequilibrium deep in the receiver grows linearly with Reynolds number, Λ₁₀₇ = 0.038 + 8.3×10⁻⁴ Re, measured as a conservative lower bound. Local thermal equilibrium does not hold anywhere in the tested envelope.
4. The slow-mode eigenvalue analysis identifies C_eff = 301 ± 23 J K⁻¹ and K_loss = 0.10–0.16 W K⁻¹ (temperature-bracketed) from the cooling data, confirmed within 1% by fifteen independent estimates from the heating data. With the measured 40 g monolith mass, the housing accounts for about six sevenths of the participating thermal mass.
5. All fifteen heating transients collapse onto master curves under the single time scale C_eff/(ṁc_p + K_loss), with coefficients of variation of 9–13%; the wall-to-gas lag ratio is a factor of three. The installed receiver exhibits one-parameter transient dynamics.
6. The measured gas output exceeds the nominal aperture-power accounting by 23 ± 4% at the highest flow. An energy closure, corroborated by an independent model calibration, localizes the delivery error by lamp configuration (approximately +34–37% at 456 and 304 kW m⁻², −21% at 256 kW m⁻²). With per-configuration calibration factors applied, all efficiencies are sub-unity and rise linearly with flow to about 90% recovery of delivered power; the receiver is flow-limited, not exchange-limited, throughout the campaign.
7. A two-zone continuum model calibrated simultaneously against all fifteen heating and three cooling runs reproduces the measured effectiveness envelope to a mean absolute deviation of 0.018 in ε, yet fails to reproduce the volumetric inversion: it inverts only at 0.06–0.16 more gas-side recovery than the receiver does, reverses the weak flux ordering of the measured threshold, and reaches an indicator amplitude an order of magnitude too small. The residuals localize the failure to the axial distribution of wall temperature rather than to the magnitude of gas–solid coupling. Separately, the aperture subtends about 7% of the beam, so the source magnitude and the captured spillage fraction enter the temperature field only as a product and are individually unidentifiable; the delivered power and the receiver conductance are therefore not recoverable by inverse modelling of this instrumentation, and are reported here only where they follow algebraically from the measurements.

---

## Data availability

The complete raw dataset, the reduction and analysis scripts, and all derived data products are available at the project repository (github.com/tamu-edu/tamuq-chen-secarelab-receiver). The published files include the corrected thermocouple map and flow-channel designations.

## Acknowledgments

[To be completed.]

## References

[1] International Energy Agency, CO₂ Emissions Report. *(complete citation from v5 library)*
[2] T. Fend et al., volumetric receiver materials and performance. *(v5 refs)*
[3] A.L. Ávila-Marín, volumetric receiver review. *(v5 refs)*
[4] A. Kribus et al., on the volumetric effect. *(v5 refs)*
[5–8] Nu/Sh correlations for structured and porous absorbers, as compiled in Table 1 of the experimental companion description. *(v5 refs 17–24)*
[9] HFSS design, operation, and flux characterization. *(v5 ref 25)*
[10] Lambertian target method. *(v5 ref 26)*
[11] SiC property data. *(v5 ref 31)*

*Reference list to be finalized against the v5 bibliography; bracketed notes indicate the corresponding v5 entries.*

---

## Figure captions

**Figure 1.** Receiver assembly and instrumentation: honeycomb monolith in the insulated cavity, with side-wall probes T8/T12/T11 (11/58/107 mm), interior probes T9/T10 (58/107 mm), gas outlet probe T3 (≈136 mm, off-axis), insulation probe T2, and flow train. *(Revised from the earlier apparatus figure; probe designations corrected.)*

**Figure 2.** Steady side-wall and outlet gas temperatures versus air flow for the three irradiance levels. Flow sensitivity decreases monotonically with depth and vanishes at the rear. *(fig2_ss_vs_flow.png)*

**Figure 3.** Left: wall inversion metric T12 − T8 versus flow with zero crossings. Center: interior-minus-wall temperature differences at 58 and 107 mm. Right: gas thermal efficiency on delivered (solid) and nominal (hollow) power bases. *(fig3_inversion_gaps_eta.png)*

**Figure 4.** Left: apparent Nusselt number versus Reynolds number with the fitted power law and the fully developed laminar reference. Center: nonequilibrium indices Λ₅₈ and Λ₁₀₇ versus Re. Right: effectiveness versus flow. *(fig7_Nu_LTNE_eps.png)*

**Figure 5.** Cooling decays of the normalized excess temperatures (logarithmic scale) for the three cooling runs, and regression of the slow eigenvalue against the gas advective conductance including the fifteen heating-approach eigenvalues. *(fig4_cooling_lin.png; eigenvalue_verification.csv)*

**Figure 6.** Collapse of the fifteen normalized heating transients in rescaled time t* = t(ṁc_p + K_loss)/C_eff for the energy-weighted wall temperature and the outlet gas temperature. *(fig6_master_curve.png)*

**Figure 7.** Continuum model comparison. Left: modelled against measured effectiveness for the fifteen heating runs, with the ±0.03 band; the model spans the measured envelope with a mean absolute deviation of 0.018. Right: inversion indicator I_vol = T12 − T8 against effectiveness, measured (filled symbols, solid) and modelled (open symbols, dashed), by irradiance group; the shaded band marks the measured threshold. The modelled inversion is displaced to higher effectiveness, becomes flux dependent, and is an order of magnitude weaker. *(fig8_model_comparison.png)*

## Table captions

**Table 1.** Experimental campaign and dimensionless envelope.

**Table 2.** Identified constants with 95% intervals: Nusselt law (prefactor and exponent), ε* per flux level, Λ₁₀₇(Re) coefficients, ε(q) regression, C_eff and K_loss from both datasets, monolith capacitance from measured mass, and delivered-power factors with their closure brackets.
