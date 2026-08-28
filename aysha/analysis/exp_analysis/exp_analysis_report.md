# Experimental Data Analysis: SiC Honeycomb Solar Receiver (rev. 2)

> Historical note (2026-08-28): this report records an early reduction. Its
> approximate `K_loss = 0.06 W/K` statement is superseded by the later
> heating/cooling eigenvalue analysis in
> `analysis/manuscript/data_analysis_sections.md`: `0.096 +/- 0.011 W/K` for
> cooling decays and `0.164 +/- 0.028 W/K` for heating approaches. The two later
> values are interpreted as a temperature-dependent bracket, not competing
> estimates of one constant.

Analysis of the 18 heating and 3 cooling runs in `RAW/`, using the corrected thermocouple topology: the side-wall axial chain is **T8 (11 mm) – T12 (58 mm) – T11 (107 mm)**; **T9 (58 mm) and T10 (107 mm) are interior probes exposed to the channel flow**, so they are biased toward the local gas temperature. Radial transfer is read from the pairs T9–T12 and T10–T11 combined with T2 (40 mm into the insulation). Companion files: `steady_state_metrics.csv`, `flow_slopes.csv`, `cooling_time_constants.csv`, `fig1`–`fig5`, and `exp_analysis.py`.

## 1. System comprehension (from background/)

SiC honeycomb monolith, L = 137 mm, frontal width ≈ 19 mm (A_frt ≈ 3.61 cm²), 10×10 square channels of 1.5 mm (porosity ≈ 0.62, Dh = 1.5 mm), end-irradiated at three nominal aperture fluxes (456 / 304 / 256 kW/m²) with air drawn axially through the channels. T3 reads the gas at ≈ 136 mm; T15/T16 give ambient; flow is the sum of the four MFC air streams. The modelling history (0D v1–v4, 1D v3–v12, COMSOL) leaves three open problems: reproducing the volumetric inversion, the over-predicted flow-sensitivity of rear temperatures, and the gas/solid cooling hysteresis.

## 2. Analysis strategy

Steady-state values (last 120 s) for all sensors per run; flow-slope regressions of each quantity within each irradiance group; t90 rise times; two-window exponential fits of the cooling runs to separate fast local decay from the slow system eigenmode; and a model-free identification of effective capacitance and parasitic loss from the slow-mode rate versus gas advective conductance. FPT0082/0083 (absent from both import scripts) serve as replication checks. Pipeline validation: the extracted flow slopes at 456 kW/m² (T8: −33.2, T10: −3.5, T3: +0.6 K/(L/min)) match the values quoted in `journal.1D.md`.

## 3. Findings with the corrected sensor topology

**Wall-chain volumetric inversion is flow-driven.** I_vol,wall = T12 − T8 rises with flow at ≈ +17 / +11 / +7 K/(L/min) (456/304/256 kW/m², r² ≈ 0.8–0.9) and crosses zero near 11 / 10 / 4–5 L/min. The wall peak moves from the front face to ≥ 58 mm only when the through-flow is strong enough to quench the front — the mechanism is convective (entrance-region heat transfer at z < 11 mm, where developing-flow Nu is 2–3× the fully developed value), not purely optical. This is consistent with the COMSOL finding that shifting `source_y` alone cannot produce the inversion.

**Rear wall is almost flow-invariant.** Along the wall chain the flow slopes collapse monotonically with depth: T8 −21 to −33, T12 −13 to −17, T11 only −1 to −5 K/(L/min), and gas T3 ≈ 0. Whatever sets the rear wall temperature is nearly decoupled from the through-flow — the strongest experimental constraint for the bypass / two-temperature closure proposed in `journal.1D.md`. A candidate model must produce (dT11/dq ≈ dT3/dq ≈ 0) while keeping dT8/dq strongly negative.

**Interior-minus-wall depression quantifies gas-solid nonequilibrium.** T9 − T12 ≈ −25 K, essentially independent of flow and flux (slope ≲ 0.2 K/(L/min)). T10 − T11 is also negative but deepens linearly with flow: −2.38, −1.82, −0.99 K/(L/min) at 456/304/256 kW/m² with r² = 0.998, 0.969, 0.938 — the cleanest correlation in the whole dataset. Reading the interior probes as partially gas-coupled, this says the gas is still ≈ 20–55 K below the wall at 107 mm, with the deficit proportional to flow: direct evidence that thermal equilibration lengthens with velocity, i.e. the effective NTU per unit length falls as Re rises. The flux-normalized deepening rate (−2.38/−1.82/−0.99 scaling roughly with wall temperature level) can calibrate the exchange law's flow exponent.

**Radial loss path.** The wall-to-insulation drop T12 − T2 is 310–770 K, scaling with flux and only weakly with flow — the insulation resistance dominates radial transfer, and T2's cooling trace (below) shows the insulation is still charging when heating runs end.

**R_leak scaling law.** (T3 − Tamb)/(T12 − Tamb) = 0.45–0.50 + 0.014·q_lpm across all flux groups (r² ≥ 0.98): one intercept and one slope summarize gas-side energy recovery campaign-wide.

**Cooling: fast local decay + one slow system mode.** Early time constants order along the flow path (C69: T8 226 s, T12 307 s, T11 545 s, T3 1050 s; interior T9/T10 slightly slower than their wall partners), while late-window constants collapse to a common value per run (≈ 1050–1200, 1500–1620, 1950–2430 s at 10.5, 6.6, 4.5 L/min). Regressing the slow-mode rate against ṁcp gives **C_eff ≈ 300 J/K and K_loss ≈ 0.06 W/K**. The bare monolith holds ≈ 60 J/K, so γ_C ≈ 5: the cooling data measure monolith + insulation + adaptor thermal mass acting together — the physical reason every model version fits γ_C > 1, and an argument for an explicit insulation capacitance node. T2 itself barely decays over the cooling windows (its fits diverge), confirming the insulation time scale is much longer than the runs.

**Replication.** FPT0082 duplicates E77 (13.8 L/min): steady T8 agrees within 2 K. FPT0083 duplicates E81 (4.5 L/min) but runs ≈ 25–30 K hotter — check flux/alignment before using it. Both are usable as held-out validation cases.

## 4. Energy-balance audit: why η_gas > 1 is not an arithmetic error

The check for E72 (18.3 L/min, 304 kW/m² nominal): ΔT = T3 − Tamb = 356.3 K. The gas enthalpy rise was recomputed with the exact enthalpy integral ∫cp(T)dT (not cp at mean temperature) under three flow-metering conventions:

| convention | Q_gas [W] | η vs 19 mm face | η vs 20 mm face |
|---|---:|---:|---:|
| volumetric at ambient density | 133.3 | 1.22 | 1.10 |
| SLPM referenced to 0 °C | 144.7 | 1.32 | 1.19 |
| SLPM referenced to 20 °C | 134.8 | 1.23 | 1.11 |

No defensible choice of density, cp treatment, or frontal area brings E72 (or E73, η ≈ 1.0) below unity. The numerator is solid; the discrepancy sits in the denominator and, partly, in T3 itself:

1. **The nominal flux underestimates delivered power.** The project's own artifacts admit this: `import_exp_0D.jl` multiplies the 456 and 304 kW/m² levels by 1.1 and 1.15, and the COMSOL calibration uses `I_f_high` = 1.15–1.25. Applying 1.15 to E72 still leaves η ≈ 0.95–1.06 — necessary but not sufficient.
2. **Power enters the gas without crossing the monolith face.** The simulator beam overfills the 3.6 cm² face; spillage heats the cavity/insulation front, and that energy reaches the gas via the cavity air path and the hot rear adaptor upstream of T3. Io·A_frt is therefore a lower bound on the system's absorbed power, not the total.
3. **T3 is radiation-biased high.** The bead at 136 mm views the rear face at 850–1100 K. A bead balance h(T_tc − Tg) = εσ(T_surr⁴ − T_tc⁴) with h ≈ 50–100 W/m²K, ε ≈ 0.8 gives a plausible +20–60 K over-read at the hot, low-flow end; removing ~40 K from E72 alone lowers η by ~0.11. (The import code already flags the inlet TC as "strongly radiation-biased"; the same physics applies at the outlet, attenuated by higher gas velocity.)

Practical consequence: η_gas computed against nominal Io·A_frt should be reported as an *apparent* efficiency. For model calibration, treat measured Q_gas as a hard lower bound on total absorbed power per run, fit an effective delivered-power factor per flux level (consistent with the existing 1.1–1.25 factors), and either shield-correct T3 or carry a radiation-bias term ∝ (T_rear⁴ − T3⁴) in the observation model.

## 5. Implications for the next model iteration

The corrected topology strengthens the v10 physics-constrained direction: an entrance-region (developing-flow) Nusselt law concentrating exchange in the first few cm to reproduce the flow-driven wall inversion; a rear bypass fraction φ(Re) or two-temperature gas branch to make T11/T3 flow-flat while T10−T11 deepens with flow at the measured −1 to −2.4 K/(L/min); an explicit insulation capacitance node (≈ 240 of the 300 J/K total) behind a contact resistance, sized by the T12−T2 drop; and per-flux delivered-power calibration bounded below by Q_gas with a radiation-bias observation model on T3. The model-free pair (C_eff, K_loss), the R_leak law, and the gapC107(q) line should be imposed as constraints rather than left to the optimizer.
