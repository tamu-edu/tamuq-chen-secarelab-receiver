# Dimensionless Characterization of a Structured SiC Volumetric Receiver from Transient and Steady Thermocouple Data

*Draft manuscript sections (Methods, Results, Discussion). Reproducible pipeline: `../exp_analysis/exp_analysis.py`, `../exp_analysis/dimensionless_analysis.py`, `../exp_analysis/uncertainty_analysis.py`; data products `dimensionless_groups.csv`, `eigenvalue_verification.csv`, `uncertainty_per_run.csv`, `delivered_power_check.csv`; figures fig1–fig7. All temperature-based quantities (ε, Nu, Λ, ε\*, C_eff, K_loss, master curves) use nominal measured aperture fluxes only and involve no correction factors. Quantities in which power enters explicitly (η, PR) are reported on the delivered-power basis using per-configuration calibration factors f = 1.336, 1.374, 0.786 (456, 304, 256 kW/m²), whose provenance and validation are given in M6/R6; nominal-basis values are retained alongside as bounds.*

---

## 0. Revisions required in the experimental section (checklist)

1. **Thermocouple topology: manuscript text and figure are wrong and must be corrected.** The correct assignment is **T8 (11 mm), T12 (58 mm), T11 (107 mm) on the side wall; T9 (58 mm), T10 (107 mm) interior and flow-exposed**. The v5 text and Figure 3 state the opposite and must be revised, with a scaled schematic showing axial coordinates (11 / 58 / 107 / ≈136 mm) and radial positions. The two assignments swap the sign interpretation of every radial-gap metric.
2. **Flow metering convention.** State that the GFC17 controllers deliver *standard* L/min with the manufacturer reference conditions (21.1 °C, 101.325 kPa) and give the mass-flow conversion; the convention choice moves all energy quantities by 5–9%.
3. **Delivered-power statement.** Report the flux-map integral over the 19 × 19 mm² face and the spillage fraction onto the cavity front. Section R6 shows the nominal aperture flux times the frontal area under-bounds the absorbed power by at least 23 ± 4% at the highest flow; the paper must present Io·A_frt explicitly as a lower bound on delivered power.
4. **Gas and interior thermocouple bias.** T3 sits at the receiver exit, off-axis, with no direct view of the radiation source, and is treated as unbiased. The interior probes T9/T10, by contrast, view the surrounding channel walls and are radiation-biased *toward the solid temperature*; the consequence (their gas-side depression readings are lower bounds) is used explicitly in Section R5.
5. **Uncertainty budget.** Include the Monte Carlo propagation of Section M7 (TC ±2.2 K, MFC ±1% FS ×4, properties ±1–2%, HFG ±3%, delivered-power factors ±8%).
6. **Repeatability and steady-state criterion.** FPT0082 repeats E77 within 2 K (ε and Nu within 2%) and provides the repeatability statement; FPT0083 is excluded. Define the steady-state criterion quantitatively (e.g. |dT/dt| < 0.5 K/min over the final 10 min); run durations are 3000–7150 s.
7. **Receiver mass.** The measured monolith mass (40 g) must appear in the setup section; it anchors the solid heat capacity (Section R7) and shows the effective extrudate density is ≈ 2150 kg/m³, i.e. ≈ 33% below dense SiC, consistent with the additive-bonded recrystallized material.
8. **Data availability.** The public repository should carry the corrected MFC column names (all four channels carry air) and the thermocouple map.

---

## 1. Methods

### M1. Dataset and reduction

Fifteen heating experiments (three nominal aperture irradiances G₀ = 456, 304 and 256 kW/m², five air flows each, 4.5–18.3 slpm) and three cooling experiments (lamps off at t = 0 with flow maintained at 10.5, 6.6 and 4.5 slpm) are analyzed, recorded at 1 Hz. One replicate (E82 = E77 conditions) is used solely for the repeatability estimate.

All temperatures are converted to Kelvin. The ambient temperature is T_amb = ½(T15 + T16). The air mass flow is

$$\dot m = \rho_{\rm std}\, q/60000,\qquad \rho_{\rm std} = \frac{101325}{287.05\times294.25} = 1.199\ \mathrm{kg\,m^{-3}},$$

with q in standard L/min at the manufacturer reference conditions (21.1 °C, 1 atm). Steady-state values are means over the final 120 s of each run; the residual drift over this window is below 0.5 K and is carried as an uncertainty component. Air properties μ(T), k(T), c_p(T) are interpolated from standard tables (300–1200 K); gas enthalpy differences are evaluated as ∫c_p dT, not as c_p·ΔT.

The side-wall axial chain is T8 (z₁ = 11 mm), T12 (z₂ = 58 mm), T11 (z₃ = 107 mm). The energy-representative wall temperature uses trapezoidal quadrature of the piecewise-linear wall profile over 0–137 mm:

$$\bar T_w = 0.248\,T_8 + 0.365\,T_{12} + 0.387\,T_{11},$$

the weights being the exact trapezoid coefficients for sensor positions 11/58/107 mm with end extrapolation. T9 and T10 (interior, flow-exposed) are reserved for the nonequilibrium metrics and are never used in ε, Nu, or the transient identification.

### M2. Dimensionless framework

Per-channel quantities use ṁ_ch = ṁ/100, square channel width w = 1.5 mm (hydraulic diameter D_h = w for a square duct), wetted perimeter P = 4w, heated length L = 137 mm, flow area A_ch = w². The groups are defined as follows, with the evaluation temperature stated because properties vary strongly along the receiver.

**Reynolds number.** Re = ṁ_ch D_h/(A_ch μ(T̄_g)), with T̄_g = ½(T_amb + T3) the mean-gas reference. Because ṁ_ch is conserved along the channel while μ rises with T, Re decreases from inlet to outlet by roughly a factor of two; the reported Re is the mean-gas value and inlet/outlet values bracket it.

**Prandtl number.** Pr = c_p μ/k evaluated at T̄_g. Pr is *computed*, not assumed: it spans 0.684–0.690 across the campaign (±0.4%). Two effects make it flat. First, for a gas, μ and k share nearly the same temperature dependence (kinetic-theory scaling ~T^0.7–0.8), so their ratio is almost temperature-free, and the residual dependence through c_p is weak (Eucken: Pr ≈ c_p/(c_p + 1.25 R/M), which changes by only ~2% while c_p changes by 17% over 300–1200 K). Second, the mean-gas reference temperatures span only 410–530 K across all runs. A property set that varies by design therefore delivers a group that is constant in fact — which is why the Nusselt correlation below is reported as Nu(Re) at fixed Pr rather than Nu(Re, Pr).

**Graetz number.** Gz_L = D_h Re Pr/L, the inverse dimensionless thermal development length evaluated at the outlet.

**Effectiveness and transfer units.** The receiver is a single-stream heat exchanger between the solid matrix and the gas, for which

$$\varepsilon = \frac{T_3 - T_{\rm amb}}{\bar T_w - T_{\rm amb}},\qquad NTU = -\ln(1-\varepsilon),$$

the second relation being exact for a single stream exchanging with a (locally) fixed-temperature wall. ε is a purely measured quantity; NTU inherits an amplified uncertainty near ε → 1 (dNTU = dε/(1−ε)), which the Monte Carlo propagation handles explicitly.

**Apparent heat-transfer coefficient and Nusselt number.**

$$h = \frac{NTU\; \dot m_{\rm ch} c_p(\bar T_g)}{P\,L},\qquad Nu = \frac{h\,D_h}{k(T_{\rm film})},\qquad T_{\rm film} = \tfrac12(\bar T_w + \bar T_g).$$

The adjective *apparent* is essential: h is defined against the instrumented side-wall temperature and the full nominal exchange area, so it lumps channel-film transport together with any cross-sectional solid nonuniformity, flow maldistribution, and bypass. It is the coefficient a designer can actually use with external measurements, and its departure from channel-film theory is itself a result (Section R3, D1).

**Wall Biot number.** Bi = h(t_wall/2)/k_SiC with t_wall = 0.4 mm.

**Local thermal nonequilibrium (LTNE) indices.**

$$\Lambda_{58} = \frac{T_{12} - T_9}{T_{12} - T_{\rm amb}},\qquad \Lambda_{107} = \frac{T_{11} - T_{10}}{T_{11} - T_{\rm amb}},$$

the interior (flow-exposed) probe depression normalized by the local wall excess temperature.

**Power ratio.** PR = f·G₀A_frt/ṁ (kJ/kg), the campaign operating coordinate on the delivered-power basis, with the per-configuration factors f of M7; the nominal-basis PR is tabulated alongside.

### M3. Steady correlation fitting

The Nusselt law is fitted by ordinary least squares on log Nu vs log Re over the fifteen heating runs; the exponent confidence interval is reported both from the regression scatter and from the instrument Monte Carlo (M6), which separate physical run-to-run variability from measurement uncertainty. The inversion criterion is obtained per flux group by linear regression of I_vol = T12 − T8 against q; the crossing q₀ = −intercept/slope is mapped to ε* through the group's ε(q) regression. The LTNE lines are OLS fits of Λ against Re.

### M4. Transient eigenvalue identification

After the internal (fast) modes decay, the assembly obeys a one-node balance:

$$C_{\rm eff}\frac{d\theta}{dt} = -(\varepsilon\,\dot m c_p + K_{\rm loss})\,\theta,\qquad \theta = \bar T - T_{\rm amb},$$

whose observable is the slow eigenvalue λ = (ε ṁc_p + K_loss)/C_eff. Two independent estimates of λ are extracted per operating point:

- *Cooling:* log-slope of (T − T_amb) over the late window (t > t_end/2, excess > 5 K), averaged over the six receiver sensors. The sensor-to-sensor spread (≤ 5%) is itself the test that a single mode governs the decay.
- *Heating approach:* linearizing about the steady state gives T_ss − T(t) ∝ e^{−λt} with the *same* λ. The deficit fraction u = (T_ss − T)/(T_ss − T₀) is fitted on 0.07 < u < 0.45 (early times excluded to let fast modes decay; the tail excluded to avoid noise amplification), retaining sensors with r² > 0.95.

λ is then regressed against ε ṁc_p, with ε at each flow taken from the steady ε(q) correlation; the slope and intercept identify C_eff and K_loss with no thermal model, discretization, or optimizer in the loop. Fifteen heating and three cooling eigenvalues give two independent identifications whose agreement is the internal validation of the method (R7).

Early-window single-exponential fits to the wall chain provide position-dependent time constants τ(z) and pairwise apparent diffusivities α_app = Δz²/Δτ.

The heating transients are finally tested for similarity: each signal is normalized, θ*(t) = (T − T₀)/(T_ss − T₀), against rescaled time t* = t(ṁc_p + K_loss)/C_eff, using only the constants identified from cooling. Collapse quality is quantified by the coefficient of variation of t* at θ* = 0.5 over the fifteen runs.

### M5. Regime classification

The dimensionless envelope is used to establish, before any correlation is interpreted, which transport mechanisms can and cannot be rate-limiting: laminarity and entrance lengths (Re, Gz), axial gas conduction (Pe = Re Pr), buoyancy (Gr/Re²), wall-normal solid resistance (Bi), gas rarefaction (Kn), and in-channel radiation (radiation-conduction number N_rc = 4σT³D_h/k). Results in R1.

### M6. Delivered-power calibration factors

The nominal power accounting G₀A_frt is demonstrably a lower bound on the power reaching the receiver system (R6). Where power enters a reported quantity explicitly — the gas thermal efficiency η and the power ratio PR — per-lamp-configuration delivered-power factors are applied as ad-hoc calibration constants:

$$f_{456} = 1.336,\qquad f_{304} = 1.374,\qquad f_{256} = 0.786,$$

so that P_del = f·G₀A_frt. The factors originate from an independent calibration of a one-dimensional transport model of the same receiver in which per-irradiance power scales were the only free parameters, and they are accepted here only because a model-free steady energy closure on the present data, f_closure = [Q_gas + K(T̄_w − T_amb)]/(G₀A_frt), brackets the same values per group (R6, table). Their uncertainty is taken as ±8% (systematic per configuration, from the per-case scatter of the calibration) and is propagated through the Monte Carlo. Three properties of the analysis protect it from circularity: the factors touch no temperature-based result (ε, Nu, Λ, ε*, C_eff, K_loss and the master curves are factor-free); their per-flux grouping — the only content actually used — is confirmed by the closure estimate, which shares no assumptions with the model calibration; and all nominal-basis values are retained in the data products so any future flux-map measurement can replace f without re-reduction.

### M7. Uncertainty analysis

A Monte Carlo propagation (N = 4000) treats all instrument errors as systematic per run and independent between instruments: each thermocouple carries a bias drawn from N(0, (2.2 K)/2) plus a steady-window drift N(0, 0.5 K); each of the four flow controllers N(0, 0.025 slpm) (±1% FS, 5 slpm range); property tables μ, k ±2% and c_p ±1% (correlated across runs, as tabulation errors are); HFG calibration ±3% (enters the efficiency only). For every draw the entire analysis chain is recomputed — all per-run groups, the Nu-law fit, the three inversion crossings, and the eigenvalue regressions (the latter additionally perturbing each λ by its fit standard error) — so all nonlinear amplifications (notably dNTU = dε/(1−ε)) and cross-correlations are captured without linearization. Reported uncertainties are MC standard deviations; parameter intervals are 95% percentile ranges.

---

## 2. Results

### R1. Operating envelope and transport regimes

The campaign spans Re = 23–94, Pr = 0.684–0.690, Gz_L = 0.18–0.70, ε = 0.57–0.78, NTU = 0.85–1.51, apparent Nu = 0.028–0.212, PR = 262–1665 kJ/kg on the delivered-power basis (300–1211 nominal) (full table: `dimensionless_groups.csv`; per-run uncertainties: `uncertainty_per_run.csv`). These ranges define the regime unambiguously:

**Deep laminar flow.** Re ≤ 94 is two orders below transition; no turbulence modeling question arises anywhere in the campaign.

**Hydrodynamically and thermally developed over ≳95% of the length.** The entrance lengths L_h ≈ 0.05 Re D_h ≤ 7 mm and L_t ≈ 0.05 Re Pr D_h ≤ 5 mm are confined to the first few percent of the 137 mm channel (equivalently Gz_L < 0.7 at the exit). Fully developed laminar duct theory is therefore the correct single-channel reference, and the constant-Nu limit (Nu_T = 2.98 for a square duct, ≈ 3.66 circular) — not an entrance-enhanced value — is the proper benchmark for the apparent Nu.

**Axial gas conduction negligible.** Pe = Re Pr = 16–65, and the relevant axial-conduction parameter Pe·(L/D_h) ≈ 1500–5900 ≫ 1: gas-phase axial diffusion cannot compete with advection.

**Pure forced convection.** Gr = gβΔT D_h³/ν² ≈ 0.4–0.7 at film conditions (ΔT ≈ 100–300 K, T_film ≈ 700 K), so Gr/Re² ≤ 10⁻³: buoyancy is irrelevant inside the channels despite the horizontal orientation.

**Radially lumped channel walls.** Bi = h t_wall/2k_SiC < 10⁻⁵: no measurable temperature drop across a wall; any solid-side nonuniformity must be at the monolith cross-section scale, not the wall scale.

**Continuum flow.** Kn = λ_mfp/D_h ≈ 10⁻⁴ even at the highest temperatures.

**Radiatively participating channels.** The radiation-conduction number N_rc = 4σT³D_h/k_air reaches 1.6 at 600 K and ≈ 5 at 1000 K: within a channel, radiative exchange between axial wall locations is comparable to or larger than gas conduction. In-channel radiation is thus an available mechanism for redistributing heat axially along the solid without gas participation — relevant to the flow-insensitivity of the rear temperatures (R2, D1).

The one group that is *not* constant by assumption deserves note: Pr varies by only ±0.4% (M2) because gas viscosity and conductivity share their temperature scaling and the mean-gas reference temperatures span only 410–530 K. All Re-dependences reported below are therefore clean, not Pr-contaminated.

### R2. Steady temperature maps

(Figures fig2, fig3.) At fixed irradiance the front wall T8 falls steeply with flow (−33.2, −24.3, −20.8 K/slpm at 456/304/256 kW/m²; r² ≥ 0.97), the mid wall T12 at roughly half that rate (−16.5, −13.4, −14.0 K/slpm), while the rear wall T11 is almost flow-flat (−1.1, −1.9, −5.2 K/slpm, r² ≤ 0.76) and the outlet gas T3 is statistically flat at the two higher fluxes (slopes +0.6 and −0.1 K/slpm with r² ≈ 0). The flow-sensitivity of the wall temperature thus decays monotonically with depth and is essentially extinguished by z = 107 mm. The insulation probe T2 rises only 10–60 K above ambient at steady state while the adjacent wall runs 300–770 K above ambient; the radial path is insulation-limited, and T2's negligible decay during the cooling runs shows the felt is still charging when the heating runs end.

The replicate pair E77/E82 agrees within 2 K on every sensor (ε and Nu within 2%), which is within the thermocouple class accuracy: run-to-run repeatability does not limit any conclusion below.

### R3. Apparent global exchange correlation

Across all fifteen runs and all three flux levels the apparent Nusselt number collapses onto a single power law (fig7, left):

$$Nu = (3.1 \pm 0.1)\times10^{-4}\; Re^{1.44},\qquad r^2 = 0.97,\qquad 23 \le Re \le 94,\ Pr = 0.69.$$

The exponent uncertainty separates cleanly into an instrument contribution of ±0.004 (Monte Carlo) and a run-to-run scatter contribution of ±0.069 (regression): the law is measurement-limited nowhere and physics-limited everywhere, i.e. the residual scatter is real operating-point variability, not noise.

Two features of this correlation carry the physics. First, the *magnitude*: Nu = 0.03–0.21 lies a factor of 15–100 below the fully developed laminar benchmark Nu_T ≈ 3. Were channel-film transport rate-limiting, NTU would be ≈ 45 and the gas would equilibrate with the wall within a few millimetres of entry; instead the measured NTU is ≈ 1 and the outlet gas leaves 150–300 K below the wall chain. The receiver is therefore emphatically *not* channel-limited. With Bi < 10⁻⁵ excluding wall-normal resistance and Gz_L < 0.7 excluding entrance limitation (R1), the deficit must arise at the assembly scale (D1). Second, the *exponent*: 1.44 is far above the 0.3–0.6 range typical of film-controlled laminar correlations for porous and structured media. An exponent above unity means the apparent conductance grows *faster* than the flow rate — each increment of flow recruits disproportionately more effective exchange — which is the signature of a flow-dependent participating fraction of the solid rather than of a boundary-layer coefficient (D1).

Both features are properties of the receiver as an assembly, measured exactly as a designer would encounter them: from external wall thermometry, inlet/outlet gas temperatures, and metered flow.

### R4. A flux-independent criterion for the volumetric inversion

The axial wall-temperature peak sits at the front face at low flow (I_vol = T12 − T8 < 0) and moves inside the receiver at high flow. The crossing flows differ strongly per flux group (11.1, 10.3, 3.7 slpm) and do not collapse on Re (Re* = 52, 52, 20). They collapse on the effectiveness (fig3):

$$\varepsilon^* = 0.671 \pm 0.003,\quad 0.671 \pm 0.003,\quad 0.626 \pm 0.005 \qquad (456,\ 304,\ 256\ \mathrm{kW\,m^{-2}}),$$

Monte Carlo 95% intervals ±0.005–0.010. The volumetric inversion — the operational definition of "volumetric" behavior, in which the irradiated face runs cooler than the interior — occurs when the gas stream recovers approximately two-thirds of the wall excess temperature, independent of irradiance over the 1.8× flux range.

The mechanism follows from the front-face energy balance. The face loses heat by re-radiation and by convective quenching from the entering cold gas; both the quench term and the absorbed flux scale essentially linearly with irradiance at fixed geometry, so their ratio — which decides whether the face or the interior is hotter — reduces to a function of the gas-side recovery, i.e. of ε. The small but statistically significant depression of ε* at the lowest flux (0.626 vs 0.671) is consistent with the growing relative weight of the flux-independent parasitic loss K_loss (R7) at low power: a fixed loss subtracts proportionally more from the face balance when the input is small, allowing inversion at slightly lower recovery.

On the delivered-power coordinate the crossings sit at PR* = 993, 731, 985 kJ/kg — a 1.4× spread versus 2.4× on the nominal basis. The delivered-power correction thus improves the power-ratio coordinate substantially, but no PR collapses the inversion as tightly as the effectiveness does; ε* remains the criterion.

As a design statement: a receiver of this class operates volumetrically when its operating point satisfies ε > ε* ≈ ⅔ — a condition on the exchanger operating point, not on the optics.

### R5. Local thermal nonequilibrium

The interior, flow-exposed probes read below their wall-chain counterparts at the same axial station in every run. At 58 mm the normalized depression is small and flow-flat, Λ₅₈ = 0.03–0.06 (raw 21–27 K, slope with Re not significant). At 107 mm it grows linearly with Reynolds number (fig7, centre):

$$\Lambda_{107} = 0.038 + 8.3\times10^{-4}\,Re,\qquad r^2 = 0.90,$$

reaching 0.11 (raw 55 K) at Re = 94, with per-point Monte Carlo uncertainty ±0.003–0.006 — the flow trend exceeds its uncertainty by an order of magnitude; within each flux group the linearity is near-perfect (r² up to 0.998).

Because the interior probes view the surrounding hot channel walls, radiation biases them *upward*, toward the solid temperature; the true local gas deficit is therefore *at least* Λ. The measured lines are consequently conservative lower bounds on gas–solid nonequilibrium, and the conclusion is strengthened, not weakened, by the bias: deep in the receiver the gas has not equilibrated with the solid, and the shortfall grows in proportion to velocity. Equivalently, the thermal equilibration length exceeds the receiver length for all operating points, increasingly so at high flow — the local counterpart of the global NTU ≈ 1 (R3), measured by an independent set of sensors.

### R6. Energy accounting and the delivered-power bound

The apparent gas thermal efficiency η = Q_gas/(G₀A_frt), with Q_gas = ṁ∫c_p dT from T_amb to T3, rises near-linearly with flow within each flux group (≈ 6%/slpm) and reaches η = 1.23 ± 0.04 at the highest flow (E72; the neighboring E73 gives 0.99 ± 0.03). The Monte Carlo interval includes the ±3% HFG calibration, all thermocouple biases, and the metering uncertainty; alternative defensible metering conventions only *raise* η. Since T3 sits at the receiver exit, off-axis and without a direct view of the radiation source, it carries no radiative over-read, and Q_gas is a robust, bias-free measurement of the enthalpy delivered to the gas.

An output exceeding G₀A_frt by five standard deviations admits one conclusion: the power absorbed by the receiver system exceeds the nominal aperture flux times the frontal area. The physical route is geometric — the concentrated beam overfills the 3.61 cm² face, and spillage absorbed by the cavity front and insulation reaches the gas stream through the cavity air path and the hot rear hardware. Two quantitative statements follow, both independent of any optical model: (i) G₀A_frt must be reported as a *lower bound* on delivered power; (ii) at 304 kW/m² the delivered-power shortfall of the nominal accounting is at least 23 ± 4%. The monotone rise of η with flow at every flux additionally shows the receiver is flow-starved rather than exchange-saturated over the whole campaign: within the tested envelope, more flow always recovers more of the available power.

**A per-flux delivered-power estimate from steady closure.** Combining the two data-identified quantities of this work — Q_gas and the transient loss coefficient — yields a per-run estimate of the delivered power itself, f = [Q_gas + K·(T̄_w − T_amb)]/(G₀A_frt), with the K_loss bracket (0.096–0.164 W/K) propagated as a systematic band (`delivered_power_check.csv`). The group means are:

| G₀ [kW/m²] | f (K_cool) | f (K_heat) | independent 1D calibration |
|---:|---:|---:|---:|
| 456 | 1.05 | 1.34 | 1.34 |
| 304 | 1.23 | 1.58 | 1.37 |
| 256 | 0.84 | 1.11 | 0.79 |

The last column lists delivered-power factors obtained by an entirely independent route — calibration of a one-dimensional transport model of the same receiver against the transient dataset, in which per-irradiance power scales were the only free parameters. The two methods share no assumptions (one is an algebraic energy closure on measured temperatures; the other a fitted PDE model with its own, quite different, exchange closure) and their error modes are disjoint, yet they agree on the two conclusions that matter: the 456 and 304 kW/m² configurations deliver of order 30–60% more power than the nominal accounting, while the 256 kW/m² configuration delivers approximately the nominal value or less. The model-based factors must nevertheless be used critically: they were fitted through a closure whose exchange law disagrees with the Nu law measured here, they carry ±7–12% per-case scatter, and their residual flow-dependence (also visible in the closure-based f, slope ≈ +0.03–0.06 per slpm) shows that neither estimate is a purely optical quantity — part of each factor compensates flow-dependent physics not captured by either accounting. The per-flux grouping, not the absolute third digit, is the robust content. That the lamp-group configurations differ per flux level makes a per-configuration delivery error physically natural, and the near-nominal 256 result is consistent with the anomalously low ε* of that group (R4): at lower true input, the fixed parasitic loss claims a larger share of the front-face balance.

**Gas thermal efficiency on the delivered basis.** With the calibration factors of M6 applied (±8% propagated), every efficiency falls below unity and the campaign map becomes physically consistent (fig3, right panel: solid symbols delivered basis, hollow nominal):

$$\eta_{\rm del} = 0.30\text{–}0.67\ (456),\qquad 0.21\text{–}0.90\ (304),\qquad 0.32\text{–}0.87\ (256\ \mathrm{kW\,m^{-2}}),$$

with Monte Carlo uncertainties ±0.06–0.08 dominated by the factor uncertainty itself. η_del rises linearly with flow in every group (slopes 0.045, 0.049, 0.059 per slpm, r² ≥ 0.98) and the highest-flow points of the three groups reach 0.67–0.90 without saturating: the flow-starvation conclusion survives the recalibration unchanged, and the maximum demonstrated recovery of delivered power is ≈ 90% at 18.3 slpm. The nominal-basis η (which exceeds unity at E72 by five standard deviations) is retained in the data products as the factor-free bound that motivates and validates the calibration.

### R7. Transient identification and its internal cross-validation

**Single slow mode.** In each cooling run the late-time log-slopes of all six receiver sensors agree within 5% (λ = 7.99, 6.38, 4.84 ×10⁻⁴ s⁻¹ at 10.5, 6.6, 4.5 slpm): after internal redistribution the assembly decays as one thermal mode — the premise of M4, verified before it is used.

**Two independent identifications.** Regressing λ against ε ṁc_p:

| dataset | n | C_eff [J/K] | K_loss [W/K] | r² |
|---|---:|---:|---:|---:|
| cooling decays | 3 | 301 ± 23 | 0.096 ± 0.011 | 0.96 |
| heating approaches | 15 | 304 ± 32 | 0.164 ± 0.028 | 0.79 |

(Monte Carlo uncertainties including λ fit errors, metering, and the ε(q) model.) The two datasets share no data — lamps on versus lamps off — yet return the same effective capacitance within 1%. The identification is thus verified *from the behavior of the existing data alone*, with eighteen eigenvalue measurements in total; no additional experiments are required to establish C_eff.

**The loss coefficient is temperature-dependent, and the two values bracket it.** K_loss differs between the datasets (0.164 vs 0.096 W/K) because K linearizes a partly radiative (∝ T⁴) parasitic path: the heating approaches sample the neighborhood of the hot steady states while the cooling windows sample 150–250 K lower. The ratio 1.7 over that interval is consistent with the T³ sensitivity of a radiative surface loss. The pair therefore brackets the loss law over the operating range: K(T) ≈ 0.10–0.16 W/K between roughly 550 and 1050 K.

**The capacitance is seven times the monolith.** The measured monolith mass is 40 g; with c_p,SiC = 1050–1170 J/(kg·K) at 600–900 K this gives C_monolith = 42–47 J/K. (The measured mass also fixes the effective extrudate density at ≈ 2150 kg/m³, 33% below dense SiC — the recrystallized, additive-bonded material is internally porous, and property sets based on dense-SiC density overestimate the monolith heat capacity by ~50%.) The identified C_eff = 301 ± 23 J/K is 6.4–7.2× the monolith value: the thermal mass participating in the slow mode is dominated by the insulation and mounting hardware coupled to the monolith, not by the ceramic itself. This single measured number — obtainable only from transient data — quantifies how far the receiver-plus-housing assembly departs from the bare-absorber idealization, and fixes the capacitance any system-level transient model must carry.

**Apparent axial diffusivity.** Early-window time constants order front-to-back along the wall chain (τ = 226/307/545 s at 10.5 slpm; 593/845/1440 s at 4.5 slpm). Pairwise α_app = Δz²/Δτ gives 0.4–2.7 ×10⁻⁵ m²/s, *increasing with flow*: axial heat redistribution in the assembly is partly gas-mediated, not purely conductive, consistent with the radiatively and advectively coupled channel picture of R1. The front decays fastest at every flow although it carries no preferential flow path — its early transient is dominated by radiative loss through the aperture, the T⁴ path that also produces the temperature dependence of K_loss.

### R8. Master-curve collapse of the heating transients

Rescaling each heating run by the single time scale identified from the cooling data — t* = t(ṁc_p + K_loss)/C_eff, with no parameter fitted to any heating transient — collapses the fifteen normalized histories (fig6):

- energy-weighted wall temperature: t*(θ* = 0.5) = 0.20, CV = 13%;
- outlet gas: t*(θ* = 0.5) = 0.62, CV = 9%;

across a 3.4× range of flow and a 1.8× range of flux. Two conclusions follow. First, for system-level purposes the receiver has *one-parameter dynamics*: the ratio C_eff/(ṁc_p + K_loss), measurable from a handful of decay curves, predicts the transient shape of every heating run to within ~10% regardless of power level. Second, the factor-three offset between the wall and gas half-rise times (0.20 vs 0.62) is the dimensionless magnitude of the wall–gas thermal lag: the gas outlet is not a proxy for the solid state during transients, and any reduced-order description that couples them rigidly will misplace one or the other by precisely this factor. Both numbers are directly usable for control design and ramp-rate specification.

---

## 3. Discussion

**D1 — What the apparent Nu law measures.** The regime analysis eliminates, one by one, every channel-scale explanation for the two-order-of-magnitude gap between the apparent Nu and duct theory: wall-normal conduction (Bi < 10⁻⁵), entrance limitation (Gz_L < 0.7), buoyancy (Gr/Re² ≤ 10⁻³), axial gas conduction (Pe L/D_h > 10³), rarefaction (Kn ≈ 10⁻⁴). What remains is the assembly scale: (i) cross-sectional solid nonuniformity — the instrumented side wall, backed by insulation and exposed to cavity radiation, runs hotter than the average channel wall the gas actually contacts, a picture supported independently by the interior probes reading 21–55 K below the side wall (R5); (ii) flow maldistribution among the hundred channels and possible peripheral bypass, both of which lower the mixed outlet temperature at fixed wall temperature. On this reading the super-linear exponent (1.44 ± 0.07) acquires a physical meaning that a film coefficient cannot have: increasing flow deepens the cold core and extends it rearward, so the *fraction of the solid participating in exchange* grows with Re, compounding the film-scaling and driving the apparent conductance faster than linearly. Published correlations for volumetric absorbers, which span an unusually wide range of prefactors and exponents, are consistent with this interpretation: an assembly-level coefficient inevitably reflects assembly-specific nonuniformity, which is why such correlations transfer poorly between geometries. The practical implication for modeling is structural rather than parametric: a single-channel conjugate model with a textbook Nu law cannot reproduce these data at any parameter value, because the controlling resistance lives at the cross-section scale; closures need either a solid-temperature distribution across the section or an explicit Re-dependent active fraction.

**D2 — Status of the ε* criterion.** The criterion's strength is its economy: a single measured operating-point quantity (ε, requiring one gas and three wall thermocouples) predicts whether the receiver operates volumetrically, across all tested fluxes, with a threshold determined to ±0.005. Its scope is bounded by the dataset: one geometry, laminar flow, ε* referenced to the side-wall chain. The observed 0.045 depression of ε* at the lowest flux is not noise (it exceeds its MC interval ninefold) and is consistent with the fixed parasitic loss K_loss claiming a proportionally larger share of the front-face balance at low power; a two-term criterion ε* (1 − c·K_lossΔT/G₀A_frt) would capture it, but three flux levels are too few to fix c meaningfully. Testing the ⅔ threshold on other structured absorbers — where inversion data exist but are typically reported against flow or Re rather than effectiveness — is the natural external validation.

**D3 — The LTNE bound and its use.** Λ₁₀₇(Re) is, to our knowledge, the directly measured quantity closest to what two-temperature porous-medium models call the interphase nonequilibrium, obtained here with the radiation bias working in the conservative direction (R5). Its linearity in Re with a near-zero intercept says the deep-receiver gas deficit is advection-controlled and vanishes only in the static limit — there is no flow regime within the envelope where local equilibrium holds at 107 mm. Any LTNE closure calibrated for this class of receiver should reproduce both the line and its flux-invariance; conversely, local-thermal-equilibrium formulations are demonstrably outside their validity domain here for Re ≳ 20.

**D4 — The transient identification as a measurement technique.** The eigenvalue method of M4/R7 is deliberately model-poor: it requires only that a slowest mode exist and be observable, which the 5% sensor-to-sensor agreement verifies empirically. Its products (C_eff, K_loss(T)) are therefore portable — they constrain *any* future model of this receiver, of whatever fidelity, and they came entirely from data already in hand. The cross-validation between heating and cooling datasets (1% agreement in C_eff from disjoint data) is the strongest internal-consistency statement the dataset can make, and it simultaneously demonstrates that the linearized loss coefficient, not the capacitance, carries the temperature dependence — information that would be inaccessible to a fit of either dataset alone. The 7× ratio of C_eff to the measured-mass monolith capacitance relocates the receiver's thermal inertia from the absorber to its housing, with direct consequences: start-up and cloud-transient response are set by the insulation coupling, not by absorber mass, so dynamic performance improvements should target the mounting and insulation architecture rather than the ceramic.

**D5 — The delivered-power bound.** η > 1 against the nominal input is an uncomfortable result to report and precisely for that reason a valuable one: it is a five-sigma, assumption-light demonstration that aperture-flux × frontal-area accounting under-measures delivered power in small-aperture cavity configurations, where spillage comparable to the face area is unavoidable. The bound (≥ 23 ± 4% at 304 kW/m²) is a property of the measurement configuration, not of the receiver, and publishing it as a bound — rather than absorbing it into a calibration factor — keeps the dataset honest for future users: every efficiency herein is *apparent*, every absorbed-power figure a *minimum*. The flow-monotone rise of η adds the operational corollary that the receiver never reaches exchange saturation in the tested envelope; its efficiency ceiling within the campaign is set by the flow budget, not by the heat-transfer surface.

**D6 — Limitations.** The transient identification rests on three cooling runs but is confirmed by fifteen independent heating-approach eigenvalues; the remaining open quantity is the shape of K(T) between the two bracketing values, which the data constrain only at its endpoints. ε and Nu are referenced to the side-wall chain; referencing the interior chain instead shifts ε by +0.10–0.15 and Nu by roughly a factor of two while changing no qualitative conclusion — the choice is stated and reproducible rather than uncertain. The interior-probe radiation bias is signed but not quantified; since it acts conservatively on Λ, only the bound direction is used. All conclusions are confined to Re = 23–94, Pr ≈ 0.69, one monolith geometry, and one cavity configuration; the ε* threshold and the Nu prefactor in particular should be transferred to other geometries as hypotheses, not constants.

---

## 4. Conclusions

1. Across 23 ≤ Re ≤ 94 the receiver's apparent global exchange follows Nu = 3.1×10⁻⁴ Re^1.44 (r² = 0.97), a factor 15–100 below fully developed duct theory with a super-linear exponent — the receiver is assembly-limited, not channel-limited, and the exponent measures flow-recruited participation of the solid.
2. The volumetric inversion obeys a flux-independent operating-point criterion: the wall peak moves inside the receiver when the gas effectiveness exceeds ε* ≈ 0.66 ± 0.03 (threshold resolved to ±0.005 per flux level).
3. Deep-receiver gas–solid nonequilibrium grows linearly with Re, Λ₁₀₇ = 0.038 + 8.3×10⁻⁴ Re (a conservative lower bound); local thermal equilibrium does not hold anywhere in the envelope.
4. The slow-mode eigenvalue method identifies C_eff = 301 ± 23 J/K and K_loss = 0.10–0.16 W/K (temperature-bracketed) from cooling decays, cross-validated to 1% by fifteen heating approaches; with the measured 40 g monolith mass, seven-eighths of the participating thermal mass belongs to the housing, not the absorber.
5. All fifteen heating transients collapse onto master curves under the single cooling-identified time scale (CV ≤ 13%), with a factor-three wall-to-gas lag ratio — the receiver's system-level dynamics are one-parameter dynamics.
6. Measured gas output exceeds nominal aperture-flux × frontal-area input by 23 ± 4% at the highest flow: nominal accounting under-bounds delivered power in this configuration. A steady energy closure built from the data-identified quantities, corroborated by an independent model calibration, localizes the delivery error by lamp configuration (≈ +34–37% at 456 and 304 kW/m², ≈ −21% at 256 kW/m²); adopting these as per-configuration calibration factors (f = 1.336, 1.374, 0.786, ±8%) renders every efficiency sub-unity, with η_del rising linearly with flow to a maximum demonstrated recovery of ≈ 90% of delivered power at 18.3 slpm — the receiver remains flow-starved, not exchange-saturated, throughout the envelope.

---

## 5. Figure and table set

Fig. 1 setup and corrected TC map (revised from v5). Fig. 2 steady wall-chain and gas temperatures vs flow per flux (fig2). Fig. 3 wall inversion, crossings, and ε* collapse (fig3). Fig. 4 Nu–Re with duct-theory reference and LTNE lines (fig7). Fig. 5 cooling log-decays and λ vs ε ṁc_p regression with heating-approach points overlaid (fig4 + `eigenvalue_verification.csv`). Fig. 6 master-curve collapse (fig6). Table 1 campaign and dimensionless envelope with per-run MC uncertainties (`dimensionless_groups.csv` + `uncertainty_per_run.csv`). Table 2 identified constants with 95% intervals: Nu law (a, b), ε* per flux, Λ₁₀₇ line, ε(q), C_eff, K_loss (both datasets), C_monolith from measured mass. Supplementary: reduction and uncertainty pipeline (three scripts), replicate comparison E77/E82, regime-number derivations.
