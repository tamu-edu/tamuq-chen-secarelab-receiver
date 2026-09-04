# Review of `manuscript_full_draft.md`

**Manuscript:** *Transient Thermal Characterization and Dimensionless Correlations for a Structured Silicon Carbide Volumetric Solar Receiver* (Melhim, Konstandopoulos, Kakosimos)

**Basis of review.** Every reported number was recomputed from the repository's own derived data products in `aysha/analysis/exp_analysis/` (`dimensionless_groups.csv`, `flow_slopes.csv`, `eigenvalue_verification.csv`, `identified_constants_mc.csv`, `delivered_power_check.csv`, `table2_identified_constants.csv`), from the raw logger files in `aysha/analysis/RAW/`, and cross-checked against the 1D/2D Julia sources (`1D_v2 … 1D_v50`, `2D_v1 … 2D_v22`), the COMSOL summary (`background/comsol_model_analysis_summary.md`), and the modeling progress report (`receiver_analysis_v2.html`).

---

## 0. Summary verdict

The physical argument is sound and the analysis strategy — reduce to dimensionless groups, identify constants from the transient eigenstructure, avoid inverse modelling — is the right one for this dataset. Most reported quantities reproduce exactly.

**Reproduced exactly from the pipeline:**

| Quantity | Manuscript | Recomputed |
|---|---|---|
| Nusselt law (n = 15) | 3.1×10⁻⁴ Re^1.44, r² = 0.97 | 3.104×10⁻⁴ Re^1.443, r² = 0.971 |
| ε* per flux level | 0.671 / 0.671 / 0.626 | 0.673 / 0.673 / 0.628 |
| Crossing flows q* | 11.1 / 10.3 / 3.7 sL min⁻¹ | 11.08 / 10.32 / 3.69 |
| Cooling eigenvalues λ | 7.99 / 6.38 / 4.84 ×10⁻⁴ s⁻¹ | identical |
| C_eff, K_loss (cooling) | 301 ± 23 J K⁻¹, 0.096 ± 0.011 W K⁻¹, r² = 0.96 | 300.9, 0.0963, r² = 0.964 |
| T8 flow slopes | −33.2 / −24.3 / −20.8 K per sL min⁻¹ | −33.20 / −24.31 / −20.83 |
| T11 flow slopes | −1.1 to −5.2, r² ≤ 0.76 | −1.10 / −1.86 / −5.16, r² = 0.07 / 0.12 / 0.76 |
| Envelope (Re, Pr, Gz, ε, NTU, Nu) | as stated | agrees to the printed digits |
| Geometry, porosity, ρ_eff, C_monolith | 19 mm, 3.61 cm², 40 g → 2150 kg m⁻³, 42–47 J K⁻¹ | ε = 0.6233, A_solid = 1.360×10⁻⁴ m², ρ_eff = 2147, 42.0–46.8 |

**Not reproduced — these are blocking.** Six reported numbers in §4.7 and §3.5/§4.6 disagree with the pipeline that produced the paper, and two of them carry headline claims in the abstract and conclusions. Details in §1 below.

There is also one structural framing issue (§3.1) that a referee in this field will raise: the Nusselt exponent is not independent information, and the reference exponent it is compared against is the wrong one.

---

## 1. Blocking: reported values that disagree with the analysis pipeline

### 1.1 The heating-transient identification (§4.7) is wrong

The manuscript reports, from the fifteen heating eigenvalues:

> C_eff = 305 ± 32 J K⁻¹ and K_loss = 0.165 ± 0.028 W K⁻¹ (r² = 0.79)

The pipeline gives, from the same 15 rows of `eigenvalue_verification.csv`:

| | C_eff [J K⁻¹] | K_loss [W K⁻¹] | r² |
|---|---|---|---|
| Manuscript §4.7 | 305 ± 32 | 0.165 ± 0.028 | 0.79 |
| My regression | **287.0** | **0.1189** | **0.898** |
| `identified_constants_mc.csv` | **287.5 ± 26.7** | **0.1194 ± 0.0179** | — |
| `table2_identified_constants.csv` | **287** | **0.119** | — |

Three consequences, all of which propagate to the abstract and conclusions:

1. **The "within 1%" cross-validation claim fails.** 301 vs 287 J K⁻¹ is **4.7%**, not 1%. (Even the manuscript's own 305 vs 301 is 1.3%, not 1%.) The claim appears in the abstract, in §4.7, and in Conclusion 4. The two estimates remain statistically consistent — the 95% intervals [254, 348] and [235, 340] overlap heavily — so the *conclusion* survives, but the number does not. Recommended wording: "agree to within 5%, well inside their respective 95% intervals."
2. **The K_loss bracket 0.10–0.16 W K⁻¹ is not a pipeline result.** The two point estimates are 0.097 and 0.119 W K⁻¹. The value 0.16 is close to the *upper 95% bound* of the heating estimate (0.155), so the stated bracket mixes a point estimate at one end with a confidence bound at the other. `delivered_power_check.csv` confirms the pipeline's own bracket: its columns are named `f_K0.097` and `f_K0.119`. The defensible bracket is **K_loss = 0.097–0.119 W K⁻¹**.
3. **The radiative-loss argument in §4.7 loses its quantitative support.** The text reads "The ratio of 1.7 over that interval is consistent with the T³ sensitivity of a radiative surface loss." The actual ratio is 0.1194/0.0967 = **1.23**. A factor 1.23 across a 150–250 K interval is much weaker evidence for a T³-dominated path; either recompute the expected ratio for the actual temperature window and show that 1.23 is what T³ predicts, or downgrade the sentence to "the heating-window value is the larger of the two, in the direction expected for a partly radiative loss path."

`identified_constants_mc.csv` also contains a **combined 18-eigenvalue fit that the manuscript never reports**: C_eff = 275.6 ± 23.7 J K⁻¹, K_loss = 0.1037 ± 0.0145 W K⁻¹ (my regression: 274.5, 0.1028, r² = 0.894). This is the statistically preferred estimate — one fit, 18 observations, 16 degrees of freedom — and it is 8.4% below the cooling-only value. The three-point cooling fit has **one degree of freedom**, so the ±23 J K⁻¹ attached to the headline value is optimistic. Recommendation: report the combined fit as the primary identification and the two subsets as the consistency check. The C_eff/C_monolith ratio becomes ≈ 6 rather than ≈ 7, which changes Conclusion 4's "six sevenths" to "about five sixths".

### 1.2 The delivered-power factors are not traceable, and three mutually inconsistent sets are in circulation

§3.5 declares f₄₅₆ = 1.336, f₃₀₄ = 1.374, f₂₅₆ = 0.786, attributed to "an earlier one-dimensional inverse calibration". I could not confirm that provenance and the repository contradicts it:

- **Set A — {1.336, 1.374, 0.786}.** These exact values appear in the repository only as hard-coded struct defaults in `2D_v1.jl` and `2D_v2.jl` (`scale_456`, `scale_304`, `scale_256`, commented "Total power scale factor"). No 1D version's `scale_456` ever takes the value 1.336; the 1D history is 1.6695, 2.3548, 2.4106, 2.4156, 1.9325, 2.0025, 0.7022, 1.9928, 1.9812, 1.2704, 1.3400. So these are **2D** seeds, not the output of a 1D inverse calibration.
- **Set B — {1.34, 1.58, 1.11}.** This is the manuscript's own "f, closure (K upper)" column in §4.6 **and, identically, the seed values in `1D_v49.jl` and `1D_v50.jl` lines 116–118.** See §4.2 below for why this matters.
- **Set C — {1.146, 1.343, 0.931}.** This is what the archived pipeline actually applied: `delivered_power_check.csv` column `f_applied`. These are exactly the ratios 523/456, 408/304, 238/256 — i.e. the ratio of the flux values used to label the three levels in `table2_identified_constants.csv` (523, 408, 238 kW m⁻²) to the nominal fluxes in the manuscript (456, 304, 256 kW m⁻²).

**Reproducibility consequence.** The manuscript's η_del ranges reproduce under Set A (I get 0.297–0.668, 0.214–0.881, 0.316–0.874 against the stated 0.30–0.67, 0.21–0.90, 0.32–0.87), so the *text* is internally consistent. But the archived CSV, which the data-availability statement points a reader to, was generated with Set C and yields 0.347–0.779, 0.219–0.902, 0.267–0.738. A referee or reader who downloads the repository will not reproduce §4.6 or Figure 3(right). Similarly §4.1's PR_del range 262–1665 kJ kg⁻¹ does not match the archived column (311–1627).

**Also check the §4.6 table.** The lower-bracket column reproduces exactly (`f_K0.097` group means = 1.052, 1.230, 0.845 vs printed 1.05, 1.23, 0.84 ✓). The upper-bracket column does **not**: `f_K0.119` group means are 1.146, 1.343, 0.931, against the printed 1.34, 1.58, 1.11. The printed upper column corresponds to K_loss = 0.16, which per §1.1 is not a pipeline value.

**Recommendation.** The cleanest fix — and the one most consistent with the abstract's claim of working "without recourse to detailed simulation" — is to **delete Set A from the paper entirely**. §3.5 already concedes the factors "affect no temperature-based result", so nothing that matters is lost. Report η and PR on the nominal basis, present the T3-based closure as a bracket computed with the paper's own K_loss = 0.097–0.119, state the flux-accounting discrepancy as the open problem it is, and re-run the pipeline so the archived products match the printed numbers.

### 1.3 Two irradiance labellings are in use

`table2_identified_constants.csv` labels the three levels **523, 408, 238 kW m⁻²**; the manuscript uses **456, 304, 256 kW m⁻²**. These are presumably peak-versus-aperture-average flux from the Lambertian-target mapping, but the paper never defines a second flux measure, and the difference is what generates Set C above. Resolve, define once in §2.1, and use one set throughout.

---

## 2. Claims contradicted by the data

### 2.1 "Λ₅₈ is small and independent of flow" (§4.5) — not supported

Within each irradiance group Λ₅₈ rises systematically and almost perfectly linearly with Re:

| G₀ [kW m⁻²] | Λ₅₈ range | change | slope [Re⁻¹] | r² |
|---:|---|---:|---:|---:|
| 456 | 0.0306 → 0.0374 | +22% | 1.80×10⁻⁴ | 0.984 |
| 304 | 0.0343 → 0.0531 | +55% | 2.65×10⁻⁴ | 0.999 |
| 256 | 0.0472 → 0.0635 | +35% | 3.11×10⁻⁴ | 0.992 |

An r² of 0.98–0.999 is not scatter. The correct statement is that Λ₅₈ is *small* and rises with Re at roughly a third of the rate seen at 107 mm — which strengthens rather than weakens the paper, because it makes the nonequilibrium a smooth function of depth rather than a switch. This also touches Conclusion 3.

### 2.2 The pooled Λ₁₀₇ fit hides a flux-ordered intercept

Per-group fits show a common slope and a systematically ordered intercept:

| G₀ [kW m⁻²] | intercept | slope [Re⁻¹] |
|---:|---:|---:|
| 456 | 0.0301 | 8.72×10⁻⁴ |
| 304 | 0.0341 | 8.52×10⁻⁴ |
| 256 | 0.0447 | 8.36×10⁻⁴ |

The slopes agree to ±2% — a genuinely flux-invariant result worth stating as such. The intercepts differ by 48% and are monotonic in flux. The pooled fit (0.0395 + 7.94×10⁻⁴ Re, r² = 0.904; the paper prints 0.038 + 8.3×10⁻⁴, r² = 0.90) averages this away and its r² of 0.90 understates the within-group quality. Report **Λ₁₀₇ = c(G₀) + 8.5×10⁻⁴ Re** with the three offsets tabulated. Separately, §4.5's "r² up to 0.998" understates the actual per-group values (0.9964, 0.9994, 0.9997) and "up to" invites the criticism that the best group was quoted; give all three.

### 2.3 ε* is claimed to be both flux-independent and significantly flux-dependent

The abstract and Conclusion 2 assert "a flux-independent threshold, ε* = 0.66 ± 0.03", while §4.4 argues the 256 kW m⁻² value is lower by "nine times its uncertainty interval". Both cannot stand as written. The arithmetic behind the 9σ is correct (0.045 / 0.005), which makes the flux-independence claim the one that has to go: a threshold resolved to ±0.005 at each level that then moves by 0.045 between levels is flux-*dependent* and weakly so. Suggested framing: "ε* = 0.673 at both 456 and 304 kW m⁻², falling to 0.628 at 256 kW m⁻² — invariant across a 1.5-fold flux change and weakly dependent at the lowest flux."

But see §3.2 — the size of that flux dependence is a method artifact.

### 2.4 "The two datasets share no measurements" (§4.7) — overstated

§3.4 states that in the cooling regression "ε [is] evaluated at each flow from the steady correlation ε(q)". That correlation comes from the steady endpoints of the heating runs. The cooling identification therefore *does* use information from the heating dataset, and its independence is partial: the eigenvalues λ are independent, the abscissa ε ṁc_p is not. Say "the two sets of eigenvalues are measured with the lamps on and off respectively and share no temperature record; both use the same steady ε(q) regression to form the advective conductance."

### 2.5 "cannot reproduce these data at any parameter value" (§5.1) — not demonstrated

> A single-channel conjugate model equipped with a textbook Nusselt law cannot reproduce these data at any parameter value, because the controlling resistance is not in the channel.

This is an impossibility claim over a parameter space, and the experimental data alone cannot establish it. It is also precisely the claim your own modeling report warns against ("Do not claim … the impossibility of all continuum models"). Two honest options: soften to "is unlikely to reproduce … without an explicit flow-dependent participating fraction", or support it with the model evidence — which is the argument for the companion paper in §6. Note that §3.1 below gives you a *provable* weaker version of this statement, which I would use instead.

---

## 3. Methods and statistics

### 3.1 The Nusselt exponent is not independent information, and 0.3–0.6 is the wrong reference

Because h is *defined* as h = NTU ṁ_ch c_p /(PL), the correlation Nu(Re) is algebraically locked to ε(Re):

Nu ∝ NTU × ṁ × 1/k(T_film) ⇒ d ln Nu / d ln Re = 1 + d ln NTU / d ln Re

Fitting NTU against Re over the 15 runs gives an exponent of **0.389**, and 1 + 0.389 = 1.39, with the remainder to 1.443 supplied by the property variation in k(T_film). So the exponent 1.44 contains exactly the information already in the measured ε range 0.573–0.783 — it is not a second, corroborating observation.

This matters for two sentences in §4.3 and §5.1. "A value of 1.44 is well above the range 0.3–0.6 characteristic of film-controlled laminar correlations" compares against literature correlations built on a different definition of h. Under *this* definition the correct null is not 0.5 but **1.0**: any conductance that is a fixed function of position, with NTU ∝ 1/ṁ, gives Nu ∝ Re⁰; any *flow-independent* volumetric conductance gives Nu ∝ Re¹.

The strong, defensible, and far more interesting statement is:

> **The measured number of transfer units increases with flow (NTU ∝ Re^0.39), whereas a single-stream exchanger with any conductance that is a fixed function of axial position requires NTU ∝ Re⁻¹.** The receiver's gas-side coupling is therefore not merely weak — it is recruited by flow.

That is a claim about the *sign* of a derivative, it is immune to the T3 calibration question below (a flow-independent T3 bias cannot reverse it), and it is the falsifiable structural constraint that a model must satisfy. It also supersedes §2.5: rather than asserting that no parameterization can work, you can prove that no *flow-independent* one can. Panel (a) of the accompanying figure shows this directly.

### 3.2 The ε* crossings are located by a method that extrapolates outside the data

At 256 kW m⁻² the crossing is placed at q* = 3.69 sL min⁻¹, **18% below the lowest run at that flux** (4.53 sL min⁻¹). Only one point in that group is negative (I_vol = −9.4 K), the next is +25.8 K, and the indicator is strongly concave, so a global linear regression through five points under-slopes near the bottom and pushes the crossing left. Locating each crossing by interpolation between the bracketing runs instead:

| G₀ [kW m⁻²] | global linear q* | local q* | global ε* | local ε* |
|---:|---:|---:|---:|---:|
| 456 | 11.08 | 10.12 | 0.673 | 0.666 |
| 304 | 10.32 | 8.36 | 0.673 | 0.655 |
| 256 | 3.69 | 5.08 | 0.628 | 0.642 |

The spread in ε* halves, from 0.045 to 0.024, and becomes monotonic in flux. So the 9σ flux dependence of §4.4 — and equally the flux-independence headline — is **method-dependent at the current data density**, and neither survives as a resolved result. The Monte Carlo interval of ±0.005 propagates instrument error but not the choice of crossing estimator, which here is the larger term.

Two fixes, either acceptable: (a) restrict the crossing regression to the points bracketing the zero and propagate the estimator choice into the reported interval; or (b) state that the 256 kW m⁻² threshold is bounded rather than resolved (ε* < 0.64), because the receiver was already inverted at every tested flow at that flux. Option (b) is honest and costs little — the central claim, ε* ≈ 2/3 independent of flux over the tested range, holds under both methods. Panel (b) of the figure shows the two estimators.

### 3.3 T3 is load-bearing for almost every result, and its systematic bias is not in the budget

ε, NTU, h, Nu, ε*, η_gas, the closure factor f, and — through ε in the eigenvalue regression — C_eff and K_loss all derive from T3. §3.6 propagates only thermocouple class accuracy (1.1 K standard deviation) and steady-window drift (0.5 K). The probe-physics bias that §4.6 and §5.5 correctly flag as uncharacterized is a *different and much larger* term. With T̄_w − T_amb ≈ 600 K:

| systematic T3 bias | shift in ε | vs quoted ε* uncertainty |
|---:|---:|---|
| 10 K | 0.017 | 3.4× the per-flux ±0.005 |
| 30 K | 0.050 | 10× per-flux; 1.7× the abstract's ±0.03 |
| 50 K | 0.083 | 2.8× the abstract's ±0.03 |

Your own modeling notes state that T3 "reads neither bulk gas temperature nor tube wall temperature" but sits at a convective–radiative equilibrium between them — a bias of tens of kelvin, not 1.1 K. **Recommendation:** add a one-parameter sensitivity band. Report every T3-derived quantity with a stated δT3 = ±25 K (or whatever the probe geometry supports) systematic band alongside the Monte Carlo interval. This costs one paragraph and one extra column in Table 2, and it converts the paper's most obvious point of attack into a demonstration of rigour. Note the argument of §3.1 is deliberately constructed to survive this: a flow-independent T3 bias shifts ε but cannot make dNTU/dRe negative.

### 3.4 The thermocouple uncertainty may be optimistic by ~2×

§3.6 uses "class accuracy ±2.2 K". IEC 60584 class 1 for type N is ±1.5 K or ±0.004|t|, whichever is greater — ±4 K at 1000 °C, and T8 runs near 1200 K. Class 2 is ±2.5 K or ±0.0075|t|. Separately, the dominant error for the three sheathed *wall* probes is not the sensor class but the installation: contact resistance against the monolith side wall and radiative exchange with the cavity. Neither is in the budget. State the class explicitly, use the temperature-dependent tolerance, and either bound the installation term or declare it unquantified.

### 3.5 §3.6 understates the delivered-power uncertainty

The factors are perturbed by ±8%. The candidate values for the 256 kW m⁻² level span 0.786 to 1.11 — a 41% range (§1.2). If the factors are retained at all, the perturbation must cover their actual dispersion.

### 3.6 The master-curve time scale is inconsistent with the identified eigenvalue

§3.4 identifies λ = (**ε** ṁc_p + K_loss)/C_eff. §4.8 rescales time as t* = t(ṁc_p + K_loss)/C_eff — the ε factor is absent. Since ε ≈ 0.57–0.78, the two differ by up to a factor 1.75, which is a large part of why the wall half-rise lands at t* = 0.20 rather than near ln 2. Either it is a typographical omission, or the collapse deliberately uses a different (and undeclared) group. State which, and if the ε factor is intentionally dropped, say why the collapse is better without it.

### 3.7 The ε–NTU relation is applied to a wall that is not isothermal

NTU = −ln(1 − ε) is exact for a single stream exchanging with a wall at *fixed* temperature. Here ε is referenced to the energy-weighted average T̄_w of a wall whose T8→T11 spread reaches several hundred kelvin and whose axial profile is non-monotonic in the inverted regime. The relation is then an approximation of unquantified bias. §3.2's "the second relation being exact" should be qualified, and a one-line estimate of the bias (e.g. from an assumed linear wall profile) would settle it.

### 3.8 The trapezoid weights are not the exact coefficients for the stated positions

§3.1 gives T̄_w = 0.248 T8 + 0.365 T12 + 0.387 T11 and states the weights "are the exact trapezoid coefficients for the probe positions with constant extrapolation to the end faces". For probes at 11 / 58 / 107 mm over L = 137 mm the exact coefficients are **0.2518 / 0.3504 / 0.3978**. The printed weights are the exact coefficients for probes at 11 / 57 / **111** mm. The numerical effect on T̄_w is under 1 K, so no result changes, but the claim of exactness is falsifiable as printed — and given the probe-mapping history noted in §2.2 of the manuscript, an apparent 111 mm is worth double-checking against the drawing.

### 3.9 Nu prefactor uncertainty — checked, no defect

Recorded here because it is the kind of thing a referee will query. The paper reports Nu = (3.1 ± 0.1)×10⁻⁴. The Monte Carlo standard deviation in `identified_constants_mc.csv` is 0.12×10⁻⁴, and §3.6 declares the convention explicitly ("Reported uncertainties are Monte Carlo standard deviations, and parameter intervals are 95% percentile ranges"). The printed ± is therefore the correct quantity under the paper's own convention, rounded from 0.12 to 0.1 — consistent with the exponent's ±0.004 (sd 0.0036) and with the separately bracketed 95% ranges. **No change needed**, beyond optionally printing 0.12 rather than 0.1 for consistency with the two significant figures used elsewhere. The exponent's decomposition into ±0.004 instrumental and ±0.069 run-to-run reproduces and is one of the better parts of the paper.

---

## 4. Cross-document consistency

### 4.1 Sensor positions: the manuscript is right, the progress report is wrong

The modeling progress report (`receiver_analysis_v2.html`, §1 Sensor Network) places T10 and T11 at **z = 91 mm**. The manuscript places them at **107 mm**. The code settles it: `1D_v2.jl`, `1D_v3.jl` (`SENSOR_POSITIONS`), `1D_v46`–`1D_v50` (`T11=sensor_positions[:T10], # 0.107 m`), `2D_v10`–`2D_v19` (`enforce_sample!(…, 107.0e-3, …)`), and the COMSOL export map (`cpt4, y = 107 mm`) all use 107 mm. **No manuscript change needed; correct the progress report** before any of it reaches a draft. The report is also internally inconsistent — its §14 refers to "z = 107 mm".

### 4.2 The progress report's "power scale conflict" is circular

`receiver_analysis_v2.html` §12 tabulates model values scale_456 = 1.34, scale_304 = 1.58, scale_256 = 1.11 against the manuscript's f = 1.336, 1.374, 0.786, labels the last two "CONFLICT (15% gap)" and "CONFLICT (41% gap)", and builds §13 Argument III on them ("The Power Scale Discrepancy at 256 kW/m² Is a Research Finding").

Those three model numbers are the **seed values** in `1D_v49.jl` and `1D_v50.jl` lines 116–118, not fitted outputs — and the report itself states that v50 is "seeded, not calibrated" and that v49's fitted parameters p[24]–p[28] are numerically inactive. They are also numerically identical to the manuscript's own §4.6 "closure (K upper)" column. So the claimed model-versus-experiment conflict is a comparison between two of the manuscript's own power estimates. **Argument III cannot be used as written.**

### 4.3 The progress report's run-to-flux mapping is wrong

The report's §1 table assigns 456 kW m⁻² → E67, E68, E69; 304 → E70–E77; 256 → E78–E81, with flows to 18.7 LPM and Re to 130. The data say 456 → E67–E71, 304 → E72–E76, 256 → E77–E81, with Re ≤ 94 (`dimensionless_groups.csv`), and the COMSOL summary agrees with the data. The manuscript's Table 1 is correct.

### 4.4 E72's flow rate differs between data products

`dimensionless_groups.csv` and `delivered_power_check.csv` give q = 18.320 sL min⁻¹; `eigenvalue_verification.csv` gives 18.708 — a 2.1% difference for the same run, which feeds directly into the abscissa of the eigenvalue regression. Reconcile.

### 4.5 A quantitative model-comparison section was already drafted and removed

`analysis/manuscript/_session_backups_2026-08/manuscript_full_draft_pre_audit_backup.md` contains a model-versus-experiment section reporting modelled ε* = 0.786 / 0.764 / 0.726 against measured 0.626 / 0.671 / 0.671, modelled crossing flows q* = 13.7 / 16.3 / 13.5 against 3.8 / 10.4 / 11.0, and a modelled inversion amplitude of +0.7 / +5.5 / −9.5 K against measured +58 / +60 / +59 K. `analysis/exp_analysis/fig8_model_comparison.png` and `fig8_model_comparison_data.csv` still exist; the current draft has no Figure 8. Those numbers are the strongest model-side result in the whole campaign and they belong in the companion paper (§6), not in a backup directory.

### 4.6 Pressure drop is logged but unusable as instrumented — and this is fixable cheaply

§2.2 states that pressure drop across the monolith is monitored (Keller PD33X, ±200 mbar, ±0.1% FS), but no Δp is reported anywhere. It is in the raw files. Steady-window means of `DP1 (mbar)` against a laminar square-duct prediction (Po = 56.91, all metered flow through 100 channels, μ and ρ at T̄_g):

| run | q [sL min⁻¹] | Δp predicted [mbar] | DP1 measured [mbar] |
|---|---:|---:|---:|
| E67 | 15.28 | 0.99 | 2.74 |
| E72 | 18.32 | 0.99 | 2.49 |
| E77 | 13.85 | 0.58 | 0.98 |
| E76 | 4.53 | 0.24 | 0.16 |
| E81 | 4.53 | 0.20 | **−0.02** |

Three things follow. The measured signal exceeds the monolith-only prediction by 1.7–2.8× at high flow, so DP1 evidently spans more than the monolith faces. It falls *below* prediction at low flow and goes negative, which is diagnostic of a zero offset. And the transducer's ±0.1% FS on a ±200 mbar range is **±0.2 mbar** — the same order as the entire expected monolith Δp (0.20–0.99 mbar) across the whole campaign. `DP2` is a flat −0.28 to −0.32 mbar in every run and appears unconnected.

So Δp as configured cannot constrain flow distribution, and this is the right thing to say in §5.5 rather than omitting the measurement silently. **A ±2 to ±5 mbar differential transducer tapped across the monolith faces alone would resolve the signal to a few percent.** That is a far cheaper intervention than the outlet calorimetry the progress report ranks highest, and — see §6 — it is the measurement that discriminates between the two candidate mechanisms for the whole story. I would put it first on the instrumentation list.

### 4.7 A computed bypass surrogate is never reported

`delivered_power_check.csv` carries `R_leak` = 0.520–0.710, rising at +0.014 per sL min⁻¹ with r² = 0.97–0.99 (`flow_slopes.csv`). Whatever its precise definition, a quantity of that magnitude with that flow trend is directly relevant to §5.1's "possibly a peripheral bypass" and should either be defined and reported or removed from the archive. For context, the COMSOL model requires `qlpm_f_all` = 0.15–0.2, described in its own notes as a "leakage surrogate" — i.e. it reproduces the measurements with 15–20% of the metered flow through the monolith.

---

## 5. Presentation and editorial

- **Abstract contradiction.** "without recourse to detailed simulation" versus §3.5's import of factors from a model calibration. Resolved by §1.2's recommendation.
- **Table 2 is never cited in the text.** Cite it where the constants are introduced (§4.3, §4.4, §4.7).
- **References are placeholders.** [1]–[11] carry "*(v5 refs)*" notes. [1] (an IEA CO₂ report) is also weak support for the specific claim about >1000 °C process heat; cite a process-heat assessment.
- **Figure filenames do not match figure numbers.** Figure 4 → `fig7_Nu_LTNE_eps.png`, Figure 5 → `fig4_cooling_lin.png`. Harmless in the repository, but rename before submission so the archive is navigable.
- **Units.** "L min⁻¹" in §2.3 (cooling flows) versus "sL min⁻¹" elsewhere. Standardize on sL min⁻¹ and state the reference condition once.
- **Run bookkeeping.** §2.3 describes "fifteen heating experiments", plus one replicate that agreed and one that was excluded; Table 1 lists "5 heating + 1 replicate" at 256 kW m⁻². `dimensionless_groups.csv` contains 17 rows (E67–E83), with both E82 and E83 present and no irradiance label. The Nu fit correctly uses n = 15. State explicitly that both replicates are excluded from all fits and that they appear in the archive for transparency.
- **"T12 at about half that rate" (§4.2)** holds at 456 kW m⁻² (16.5/33.2 = 0.50) but not at the lower fluxes (13.4/24.3 = 0.55; 14.0/20.8 = 0.67). Say "half to two-thirds".
- **NTU maximum** printed as 1.51; the campaign maximum is 1.528 (E82, a replicate) and 1.518 among the fifteen. Use 1.52.
- **The T3 flow slope changes sign** (+0.58 at 456, −0.07 at 304, −3.05 at 256, r² = 0.03 / 0.0003 / 0.55). §4.2 mentions only that T3 is flow-independent at the two higher fluxes. The sign change is a real feature that the modelling cannot yet reproduce; worth one sentence.
- **"a factor of 15–100 below"** — the actual range is 13.9–105.9. Use "one to two orders of magnitude", as the abstract already does.

---

## 6. Should the 1D modelling go into this paper, stand alone, or anchor a new 1D + 2D/3D paper?

**Short answer: do not put it in this paper; build the new multi-fidelity modelling paper. There is a strong scientific reason for it, but not the one the progress report gives.**

### 6.1 Why not this paper

Three of the five arguments in `receiver_analysis_v2.html` §13 do not currently support publication:

- **Argument II (C_eff cross-validation)** is circular by construction. The 1D objective explicitly regularizes C_total toward the 301 J K⁻¹ that came from these experiments. The report concedes this and then claims the agreement anyway. A geometry-only estimate computed *without* the regularization term would be a real test; it has not been run.
- **Argument III (power-scale conflict)** is invalid — see §4.2 above. The numbers compared are the manuscript's own closure bracket against the manuscript's own declared factors.
- **Arguments IV and V (assembly-level Nu, LTNE index Λ)** are proposals, not results. Neither has been computed; the report lists Λ post-processing as a gap.

Only **Argument I (rear-temperature flow sensitivity)** stands, and it is genuinely strong.

Beyond that, the model has nothing calibrated to report: v49's flow-slope penalty and all five observation parameters are numerically inactive because of a duplicate loss-function definition, and v50 is repaired but seeded. Attaching an uncalibrated model with a documented non-identifiability result to a clean, model-free experimental paper transfers the model's liabilities onto the paper's strongest asset and invites a referee to question the correlations by association.

The reverse dependency is the one to act on now: this manuscript currently *leans* on unpublished modelling in three places, and all three should be cut regardless of what you decide about the companion paper — the f factors (§1.2), the impossibility claim in §5.1 (§2.5, replaced by the provable statement in §3.1), and the dangling §5.6 reference to "the accompanying 1D or 2D inverse models", which should become a forward citation to a named companion or be deleted.

### 6.2 Why the 1D work cannot stand alone as it is

Two framings are available and both fail at present. As a **coefficient-extraction** paper it is ruled out by your own v20 result — receiver UA(Re) is not uniquely identifiable from this measurement set, and the fitted parameters are boundary-active. As a **negative result** ("our 1D model overpredicts rear flow sensitivity by 4–5×") it is thin: a single model class failing on one dataset, with a known implementation defect in the version that produced the numbers, is a technical note at best.

What converts it into a real paper is a *mechanism plus a demonstrated remedy* — which requires the 2D and 3D models. That is exactly the paper you want to write.

### 6.3 The thesis I would give the multi-fidelity paper

The audit above supplies a sharper organizing claim than the report's five-argument list, and it unifies every discrepancy in the campaign into one:

> **The receiver's transfer-unit count increases with flow (NTU ∝ Re^0.39), and no single-stream formulation can produce that. The missing physics is flow-participation that grows with Reynolds number.**

Why this is one story rather than five findings:

1. The experimental Nu ∝ Re^1.44 *is* the statement NTU ∝ Re^0.39 (§3.1) — the correlation and the anomaly are the same measurement.
2. NTU rising with flow *is* ε rising with flow, which *is* T̄_w falling steeply while T3 stays flat — which *is* the rear-temperature flow-invariance that the 1D model overpredicts by 4–5×. Argument I and the Nusselt exponent are one phenomenon seen from two ends.
3. The ε* ≈ 2/3 inversion criterion is the same statement written at the front face, which is why the 1D model both misses the inversion and gets the rear slopes wrong.
4. The structural theorem is provable, not fitted: in a single-stream LTNE formulation with a volumetric conductance that is any fixed function of position, NTU ∝ Re⁻¹. With a developing-flow closure of the Shah–London type, h grows no faster than about Re^0.5 in the entrance region, giving NTU ∝ Re^−0.5 at best. The gap to the measured Re^+0.39 is roughly **Re^0.9** — nearly a full power. That is a quantitative structural deficit, not a calibration residual, and it is immune to the identifiability problem because it constrains an exponent rather than a coefficient.
5. It names the remedy and makes it falsifiable. An active-participation fraction growing as Re^~0.9 closes the gap. Two mechanisms in your campaign predict exactly that sign: viscosity-driven redistribution between core and perimeter channel groups at common Δp (hot core chokes at low flow, redistribution weakens as temperatures fall at high flow → participation rises with Re), and peripheral bypass (`R_leak` 0.52–0.71 rising with flow; COMSOL's `qlpm_f_all` = 0.15–0.2). The paper's result is which one reproduces the exponent.

### 6.4 Assets you already have for it

- **1D LTNE** (`1D_v50.jl`, repaired) — supplies the structural theorem and the 4–5× rear-slope deficit.
- **2D axisymmetric r–z** (`2D_v22.jl`) plus the pressure-coupled multi-stream architecture — the candidate remedy, and the v20 non-identifiability result as an honest bound on what the data can constrain.
- **3D COMSOL Cav_Hex v7.18** — and its *calibration* is itself evidence. The two knobs it needs are `k_air_mult` = 15 (an explicit "surface area surrogate") and `qlpm_f_all` = 0.15–0.2 (an explicit "leakage surrogate"): a third, independently built model class reaching for exactly the two corrections a missing flow-participation mechanism would force. That is a much better cross-model argument than the C_eff or power-scale ones.
- **A partially independent power-factor corroboration.** COMSOL's separately tuned optical factors (`I_f_high` = 1.15–1.25, `I_f_low` = 0.85) reproduce the same super-nominal/sub-nominal ordering as the experimental T3 closure (1.05–1.15 at 456, 0.85–0.93 at 256, per §1.2 Set C). Different code, different physics, same sign and comparable magnitude. Unlike Arguments II and III this is not circular, and it is the argument I would actually make.
- **The removed model-comparison numbers** in the pre-audit backup (§4.5), already reduced by the same convention as the experimental analysis.

### 6.5 What has to be done before that paper is submittable

1. Repair and **calibrate** v50 from multiple starts to convergence, with the activity test the report specifies (each fitted parameter must move at least one output). Until then no model coefficient can be printed.
2. Implement the multi-stream 2D model and test the single prediction that matters: does it reproduce **NTU ∝ Re^+0.39**, the rear-slope magnitudes, and inversion at ε* ≈ 0.66? If it does, the paper is strong. If it does not, the bypass mechanism is the alternative and the paper becomes about discriminating them.
3. Compute the effective assembly-level Nu from the calibrated model and compare against 3.1×10⁻⁴ Re^1.44 (report Argument IV) and Λ₁₀₇ against 0.038 + 8.3×10⁻⁴ Re (Argument V). Both are cheap post-processing on an existing model and both are new results.
4. **Measure Δp across the monolith faces with a ±2–5 mbar transducer** (§4.6). Redistribution and bypass have different Δp signatures at the same total flow, so this single measurement discriminates the two mechanisms — and it is much cheaper than the outlet calorimetry the report ranks first. Outlet enthalpy is still needed to close the T3 question, but for *this* paper's thesis, pressure drop is the higher-value measurement.
5. Drop Arguments II and III from the framing, or reconstruct them non-circularly (unregularized C_total; genuinely fitted power scales compared against a closure that does not share the same K_loss bracket).

### 6.6 Recommended plan

- **Paper I — this draft, revised.** Experimental and model-free. Fix §1 and §2, restate the Nusselt result per §3.1, add the T3 sensitivity band per §3.3, delete the f factors, re-run the pipeline so the archive reproduces the text. Submit as it stands otherwise; the correlations, the ε* criterion, the Λ relation, the eigenvalue identification, and the similarity collapse are five distinct contributions and that is a complete paper.
- **Paper II — new.** *Structural limits of reduced-order models for a structured SiC volumetric receiver: 1D LTNE, 2D multi-stream, and 3D CFD.* Thesis per §6.3. Cites Paper I for the data. Under the discipline of §6.5 this is a stronger paper than the combined version would have been, because its central claim is a provable structural constraint rather than a set of fitted coefficients that your own analysis has already shown to be non-identifiable.
- **Journals.** Paper I fits *Solar Energy* or *Applied Thermal Engineering*. Paper II is a better fit for *International Journal of Heat and Mass Transfer*, where a structural non-identifiability result with a multi-fidelity remedy is a recognized contribution type.

One caution on Paper II's framing: keep the Role-A language for every fitted coefficient, but do **not** apply it to the exponent argument. "NTU rises with Re and single-stream formulations cannot produce that" is a structural claim about a derivative, and hedging it as a model-dependent effective parameter would discard the paper's strongest result.
