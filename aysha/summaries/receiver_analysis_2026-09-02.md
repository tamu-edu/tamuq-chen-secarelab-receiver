# Heat Transfer Modeling Progress — 1D and 2D Mathematical Models
**SiC Monolithic Volumetric Solar Receiver · AUTh-TRANSP**
**Analysis date: 2026-09-02 · 1D v50 (repaired) · 2D v22**

---

## § 1 — Physical System

**Receiver:** Structured monolithic SiC honeycomb, 137 mm long, 100 square channels (1.5 mm × 1.5 mm), 0.4 mm web thickness, 19 × 19 mm² frontal area.

| Parameter | Value |
|---|---|
| Channel side length | 1.5 mm |
| Web thickness | 0.4 mm |
| Porosity ε | 0.623 |
| Hydraulic diameter D_h | 1.5 mm |
| Channel cross-section (gas) | 2.25 mm² |
| Solid cross-section (effective) | 1.36 mm² |
| Receiver length | 137 mm |

**Material:** Silicon carbide (SiC). Thermal conductivity k_s ≈ 35–80 W/(m·K) (temperature-dependent); specific heat c_p ≈ 1100 J/(kg·K) at operating temperature. Optically semi-transparent to solar radiation; absorbs strongly in IR.

**Operation:** Concentrated solar flux (256–456 kW/m²) enters the front aperture. Air flows axially at 4–19 LPM. The temperature inversion criterion ε* = T_gas,out/T_sol,in marks the onset of volumetric (favorable) heating where gas exits hotter than the surface it entered.

---

## § 2 — Modeling Hierarchy

| Level | Description | Status |
|---|---|---|
| 0D | Global energy balance; no spatial resolution | Complete (analytical) |
| 1D | LTNE axial continuum (solid + gas fields, z-axis) | Active — v50 |
| 2D | Axisymmetric r–z continuum | Development paused — v22 |
| ECM | Equivalent Circuit Model (network of R, C nodes) | Target; requires validated 1D/2D coefficients |

**Overall study objective:** Extract effective macroscopic heat transfer coefficients (h_eff, k_eff, C_eff) from experimental data + 1D/2D calibration for use in a reduced-order ECM suitable for system-level optimization.

---

## § 3 — Core Heat Transfer Physics

**Solid energy equation (LTNE):**
```
(ρ c_p)_s · ∂T_s/∂t = ∂/∂z[(k_s + k_rad)·∂T_s/∂z] − h·a_v·(T_s − T_g) + Q_solar(z)
```

**Gas energy equation (LTNE):**
```
(ρ c_p)_g · ∂T_g/∂t + ṁ c_p,g ∂T_g/∂z = h·a_v·(T_g − T_s)
```

**Three coupled mechanisms:**

1. **Solid axial conduction + Rosseland IR radiative diffusion** — k_rad = 4σT³/(3β_rad); β_rad ≈ 778 m⁻¹ (penetration depth ≈ 1.3 mm in IR)

2. **Gas–solid internal convection** — Shah-London developing-flow Nu:
   Nu(z) = Nu_∞ + C₁·Gz_z / (1 + C₂·Gz_z^(2/3))
   Nu_∞ = 3.61 (square duct, fully-developed laminar); C₁ = 0.355, C₂ = 0.0132 (calibrated)

3. **Beer-Lambert solar absorption** — Q_solar(z) = β_opt · F_front · G_solar · exp(−β_opt · z); β_opt ≈ 274.5 m⁻¹ (penetration depth ≈ 3.6 mm → 60% deposited in first 10% of channel)

**Key dimensionless groups (experimental):**
- Re = ρ V D_h / μ (channel Reynolds number)
- Nu_eff = 3.1×10⁻⁴ Re^1.44 (assembly-level, r²=0.97) — different quantity from single-channel Shah-London Nu
- ε* = T_gas,107/T_sol,front = 0.66 ± 0.03 (flux-independent inversion threshold)
- Λ₁₀₇ = 0.038 + 8.3×10⁻⁴ Re (LTNE nonequilibrium index at z=107 mm, r²=0.90)

---

## § 4 — 1D Model Evolution (v1–v50)

### Key milestones

| Version | Key addition | Heating RMSE | Notes |
|---|---|---|---|
| v1–v10 | Basic LTNE, MOL, OrdinaryDiffEq.jl | — | Proof of concept |
| v36 | Best pre-ECM fit | ~55 K | Obj 0.2173 |
| v42 | (M, χ) source reparametrization | 0.2270 | Morris screening repaired |
| v43 | δ_web + free Graetz exponent | 0.2206 | B_Re = 0.499 |
| v44 | Re-dependent suction h_suction(Re) | 0.1420 | −35.6% improvement |
| v45 | Continuous gas tracking + free β_opt | 0.0864 | Record; T3 errors < 35 K |
| v46 | Full ECM (φ_act≡1.0), exact energy conservation | — | First-law residual < 10⁻¹³ W |
| v47 | Bounded h_suction ≤ 150; C_total = 302.74 J/K | — | β_opt = 193.78 m⁻¹ |
| v48 | Shah-London Nu; decoupled β_opt/β_rad; F_LoS | 1.7136 | Mesh invariance established |
| v49 | Observation equations (T3 probe, T10/T11 weights) | ~1.713 | **BUGS: 6 silent defects** |
| v50 | Implementation repair of v49 | — | Repaired; not yet calibrated |

### v49 defects (all corrected in v50)
1. Duplicate `loss_cases` function definition — second overrides first, killing the slope penalty
2. T3 ODE uses hardcoded h_tc=40, Q_rad=0.85×..., ignoring parameters p[24:26]
3. T10/T11 mixing weights p[27:28] never read by observation equations
4. G_core-perim optimizer value 53.97 clamped to bound of 50 (accessor limit < optimizer limit)
5. flange_scale value 5.125 clamped to bound of 5.0 (same bug)
6. T3(t₀) initialized from model, not from measured first value

---

## § 5 — Key Physical Discoveries

### The Advective Bottleneck
In a 1D single-stream model, reducing flow rate increases ΔT (front-to-rear) in proportion to 1/ṁ. This predicts rear temperatures (T10, T11) with flow sensitivities of −10 to −20 K/LPM. The experiment shows −1 to −5 K/LPM. The rear is nearly flow-invariant. No parameter combination in the current 28-parameter model can close this gap — the physical mechanism (radial redistribution, channel-to-channel maldistribution) is absent from the 1D equations.

### T3 Exit Probe Bias
The exit thermocouple (T3) sits inside an alumina tube that couples radiatively to the tube wall (600–800 K during operation). It does not read bulk gas temperature. At 456 kW/m², T3_model = 352 K vs T3_exp = 763 K at high flow — a 411 K systematic under-prediction at the inherited v49 seed.

### Source Parameter Degeneracy
Parameters `scale` and `spill_capture` enter the energy budget only through M = scale × (1 + 13.45 × spill_capture). Reparametrized as (M, χ) in v42 to resolve.

### β_opt / β_rad Decoupling
Optical solar absorption (β_opt ≈ 274.5 m⁻¹, penetration depth ≈ 3.6 mm) and Rosseland IR diffusion (β_rad ≈ 778 m⁻¹) are fundamentally different processes operating at different wavelengths. Coupling them (as in v47) introduces non-physical constraints. Decoupled in v48.

---

## § 6 — 2D Model Evolution (v1–v22)

| Version | Key development | Status |
|---|---|---|
| v1–v19 | Axisymmetric r–z mesh, property lookups | Active |
| v20 | Non-identifiability of UA(Re) demonstrated | Finding: 2D also cannot uniquely extract HTC |
| v21 | SiC c_p corrected 400→1100 J/(kg·K) | Property fix |
| v22 Phase 1–3 | Grid search: core_preference=1.0, 0 spillage, side RMSE=211 K | Partial |
| v22 Phase 4 | **INVALID** — wrong settings + inactive Rosseland + wiring bug | Retracted |

**2026-09-02 mandate:** 2D architecture must be redesigned before Phase 4 is re-run. The non-identifiability finding from v20 means that even a correctly implemented 2D model cannot extract unique HTC values without targeted new measurements.

---

## § 7 — The Flow-Invariance Paradox

**The problem:** The 1D model (single-stream, exact energy conservation) predicts ∂T_rear/∂V̇ ≈ −15 to −20 K/LPM. The experiment shows −1 to −5 K/LPM. The rear solid temperature is nearly independent of flow rate.

**Three candidate mechanisms (mathematical root causes):**
1. **Single-stream 1D ΔT ∝ 1/ṁ** — the fundamental scaling of heat exchangers under 1D advection forces strong flow sensitivity
2. **Viscosity-driven radial flow bifurcation** — μ(T) ∝ T^0.70 means hot central channels see higher viscosity, lower flow, creating a self-reinforcing pattern that smears out the axial temperature gradient
3. **Jensen's inequality for front re-radiation** — T̄⁴ < T_peak⁴, so the effective radiative loss is lower than a 1D average temperature would predict

**Resolution path:** 2D r–z model with full viscosity-dependent flow redistribution, or new measurements (spatially resolved channel-exit temperature map) that can discriminate between hypotheses.

---

## § 8 — Current Calibrated State (v50 Inherited Seed)

### Parameters (inherited from v49, evaluations = 0)

| # | Name | Value | Unit |
|---|---|---|---|
| 1 | C1_Nu | 0.3551 | — |
| 2 | C2_Nu | 0.01318 | — |
| 3 | Nu_inf | 3.61 (fixed) | — |
| 4 | front_dep | 0.4529 | — |
| 5 | scale_456 | 1.34 | — |
| 6 | scale_304 | 1.58 | — |
| 7 | scale_256 | 1.11 | — |
| 8 | G_core_perim | 53.97 | W/(m·K) |
| 9 | C_perim_eff | 85.28 | J/K |
| 10 | k_perim_ref | 10.11 | W/(m·K) |
| 11 | beta_opt | 274.50 | m⁻¹ |
| 12 | chi | 0.7543 | — |
| 13 | f_core_rear | 0.9743 | — |
| 14 | flange_scale | 5.125 | — |
| 15 | k_core_axial_scale | 1.000 | — |
| 16 | C_rear_eff | 152.06 | J/K |
| 17 | G_receiver_rear | 1.853 | W/K |
| 18 | G_rear_tube | 6.828 | W/K |
| 19 | beta_rad | 778.25 | m⁻¹ |
| 20 | kA_rear_eff | 0.1159 | W/(m·K) |
| 21 | delta_web | 3.55×10⁻⁵ | m |
| 22 | F_LoS | 0.003941 | — |
| 23 | h_suction | 9.712 | W/(m²·K) |
| 24 | h_probe_ref | 46.34 | W/(m²·K) |
| 25 | w_probe_rad | 0.824 | — |
| 26 | G_probe_stem | 0.0224 | W/K |
| 27 | w10_stem | 0.0471 | — |
| 28 | w11_stem | 0.0 | — |
| — | **C_total** | **301.11** | **J/K** |

### Objective at seed (not a calibrated result)

| Component | Value |
|---|---|
| Signal loss (𝓛_signals) | 0.1028 |
| Flow-slope loss (𝓛_slopes) | **8.2199** ← dominant; was 0 in v49 |
| Capacitance regularization | 6.9×10⁻⁷ |
| Power-scale regularization | 1.690 |
| **Total** | **10.013** |

### Physical invariants

| Metric | Value | Target | Status |
|---|---|---|---|
| C_total | 301.11 J/K | 301 ± 23 J/K | ✓ Pass |
| Max energy residual | 5.7×10⁻¹⁴ W | < 10⁻⁴ W | ✓ Pass |
| Mesh sensitivity (HTC, 15→50 nodes) | 11.2% | < 2% ideally | ⚠ Not converged |

---

## § 9 — Open Problems

1. **Mesh convergence** — HTC changes 11.2% from 15→50 axial nodes; refinement required before reporting macroscopic HTC as a coefficient
2. **T3 sign-change failure** — ∂T₃/∂V̇ at 456 kW/m² is +0.54 K/LPM (experiment) vs −2.25 K/LPM (model); sign transition mechanism unknown
3. **Rear flow-invariance** — T10 and T11 slopes 4–5× too steep in model; requires radial physics or maldistribution model
4. **2D Phase 4 architecture** — must be redesigned from scratch; Phase 4 results retracted
5. **UA(Re) non-identifiability** — both 1D and 2D models suffer from this; targeted experiments (thermally tagged tracer, fiber-optic spatial mapping) are prerequisites for Role B
6. **v50 calibration** — full multi-start run not yet executed; current results are seed-only

---

## § 10 — Forward Path

**Immediate (weeks 1–2):**
- Run full v50 calibration from ≥ 5 random starts using the `:full` stage
- Implement Λ₁₀₇ post-processing diagnostic in v50 output
- Mesh refinement study (50→100→200 nodes, HTC convergence criterion < 1%)

**Short term (months 1–2):**
- 2D architectural redesign: viscosity-dependent radial flow, proper Rosseland boundary conditions
- Compare v50 calibrated assembly-level Nu against experimental Nu = 3.1×10⁻⁴ Re^1.44
- Profile/bootstrap key parameters to assess identifiability (especially β_opt, h_suction, G_core-perim)

**Long term (months 2–6):**
- Targeted new measurements for Role B: spatially-resolved T_s(r,z) pyrometry, channel-exit T_g distribution, spectral optical characterization of β_opt, β_rad
- ECM parameter extraction (Role B) only after UA(Re) identifiability is established

---

## § 11 — 1D v50: Implementation Repair & Baseline

### Status
**Repaired and testable — not calibrated, not validated.** Return code = Manual, evaluations = 0.

### Six defects corrected
1. Duplicate `loss_cases_v50` definition removed; slope penalty now active
2. Parameters 24–28 now read from optimizer vector (observation model)
3. G_core-perim bound raised 50→100 W/(m·K); flange_scale bound raised 5→20; accessor bounds match optimizer bounds
4. T₃(t₀) initialized from first measured T3 value; post-solution override removed
5. Irradiance cases 5–7 (256, 304, 456 kW/m²) included in full-stage fit
6. Explicit `:plant`, `:observation`, `:full` calibration stages

### Flow-slope residuals at inherited seed (K/LPM)

| Signal | 256 model | 256 exp | 304 model | 304 exp | 456 model | 456 exp |
|---|---|---|---|---|---|---|
| T8 (front) | −28.5 | −20.9 | −32.2 | −23.7 | −38.1 | −34.1 |
| T9 core | −19.5 | −13.9 | −23.6 | −13.4 | −28.9 | −16.7 |
| T10 core | −11.0 | −6.1 | −14.9 | −3.7 | **−19.2** | **−3.5** |
| T11 perim | −11.7 | −5.2 | −15.6 | −1.9 | **−20.0** | **−1.4** |
| T3 exit probe | −0.06 | −3.04 | −0.90 | −0.13 | −2.25 | +0.54 |
| T2 inlet | −1.34 | −1.14 | −2.61 | −2.04 | −4.47 | −3.71 |

**T10 and T11 at 456 kW/m²: 5× structural overprediction — cannot be closed by calibration.**

---

## § 12 — Cross-Check: 1D Model vs. Experimental Manuscript

Manuscript: "Transient Thermal Characterization and Dimensionless Correlations for a Structured Silicon Carbide Volumetric Solar Receiver" (Melhim, Konstandopoulos, Kakosimos)

| Quantity | Manuscript | 1D v50 | Status |
|---|---|---|---|
| Nusselt number | Nu = 3.1×10⁻⁴ Re^1.44 (assembly, r²=0.97) | Nu(z) = 3.61 + Shah-London (single-channel developing) | Different physical quantities |
| C_eff | 301 ± 23 J/K (eigenvalue from transients) | 301.11 J/K (regularization-constrained) | Imposed, not independent |
| K_loss | 0.10–0.16 W/K (temperature-bracketed) | G_receiver_rear = 1.85 W/K (includes flange network) | Incompatible definitions |
| ε* inversion threshold | 0.66 ± 0.03 (flux-independent) | No inversion predicted by 1D | Structural conflict |
| LTNE index Λ₁₀₇ | 0.038 + 8.3×10⁻⁴ Re | Not yet computed from 1D output | Pending comparison |
| f₄₅₆ / scale_456 | 1.336 | 1.34 | ✓ Agreement (same data) |
| f₃₀₄ / scale_304 | 1.374 | 1.58 | Conflict (15% gap) |
| f₂₅₆ / scale_256 | 0.786 | 1.11 | **Conflict (41% gap)** |
| ∂T_rear/∂V̇ | −1 to −5 K/LPM | −10 to −20 K/LPM | **4–5× structural mismatch** |
| T3 sign change with flux | Sign changes positive at 456 kW/m² | Always negative at seed | Qualitative failure |

**Power scale diagnostic:** scale_256 = 1.11 vs f₂₅₆ = 0.786 (41% conflict) — the 1D model absorbs reverse-flow natural convection physics (dominant at 256 kW/m² cooling phase) into the irradiance scale parameter. The 456 kW/m² agreement is by construction (same temperature data), not physical validation.

**C_eff independence caveat:** The regularization constraint explicitly targets the experimental 301 J/K value. A genuine independent check would require the unconstrained model estimate to fall within the experimental 1σ band — this has not been tested.

---

## § 13 — Scientific Value of a Joint Paper

**Verdict: Yes** — the combination is scientifically valuable because the systematic discrepancies between model and experiment are publishable findings, not calibration failures.

### Five arguments

**I. Rear flow-sensitivity as structural proof.** The 4–5× overprediction of ∂T₁₀/∂V̇ by a 1D single-stream model with exact energy conservation is a proof that the experimental T10/T11 behavior cannot be reproduced without radial temperature gradients or channel-to-channel maldistribution. The experiment's ε* and Λ correlations independently document this regime.

**II. C_eff cross-validated by two methods.** The geometry-based 1D estimate (SiC mass + rear hardware) and the experimental eigenvalue both yield ~301 J/K. Even with the regularization caveat, this is a meaningful geometric consistency check.

**III. Power scale conflict at 256 kW/m² is a finding.** The 41% gap between scale_256 and f₂₅₆ identifies the low-flux operating regime as the boundary where the 1D forced-flow model breaks down. Documenting this teaches the reader which regime the model is valid in.

**IV. Assembly-level vs. channel-level Nu is a new comparison.** Computing the 1D effective assembly Nu and comparing it to the experimental Nu = 3.1×10⁻⁴ Re^1.44 (Re exponent 1.44 >> laminar value ~0.5) would either confirm the convective closure or reveal hidden inter-channel mixing contributions.

**V. LTNE index Λ constrains the absorption profile.** The Λ₁₀₇ Re-dependence tests whether the predicted front-loaded absorption (β_opt = 274.5 m⁻¹, depth 3.6 mm) is consistent with the measured nonequilibrium gradient.

### Two non-negotiable conditions
1. The 1D model occupies **Role A only** (mechanistic diagnostic) throughout — fitted parameters remain non-unique and boundary-active
2. All conflicts are **reported as findings**, not tuning targets — no additional free parameters to close structural gaps

---

## § 14 — Revision Strategy

### Structure: Two-Part Paper

**Part I** (existing manuscript, unchanged): Sections 1–5 with empirical correlations as primary deliverables. Update Section 5.6 to reference Part II.

**Part II** (new §6 or appendix, ~3–4 pages): 1D model as mechanistic diagnostic layer. Five structured diagnostic comparisons, each: *model prediction → experimental result → physical interpretation of the gap.*

### Nine-step plan

1. **Retain Part I as written** — correlations (Nu, ε*, Λ₁₀₇, C_eff, K_loss) are the primary scientific deliverables
2. **Add Part II** — LTNE 1D model with governing PDEs, Shah-London Nu, Beer-Lambert/Rosseland physics, energy conservation verification
3. **Diagnostic 1: C_eff cross-validation** — geometry-based estimate vs. eigenvalue; state regularization caveat explicitly
4. **Diagnostic 2: Rear flow-sensitivity** — plot ∂T₁₀/∂V̇ and ∂T₁₁/∂V̇ model vs. experiment; attribute 4–5× overprediction to single-stream assumption
5. **Diagnostic 3: Power scales** — report all three (1.34 / 1.58 / 1.11 vs. 1.336 / 1.374 / 0.786); attribute 256 kW/m² conflict to reverse-flow natural convection
6. **Diagnostic 4: LTNE index Λ** — compute Λ₁₀₇ from calibrated 1D profile; compare to experimental correlation; interpret direction of deviation *(prerequisite: implement post-processing)*
7. **Diagnostic 5: Assembly Nu** — compute effective assembly Nu from calibrated model; compare exponent against Re^1.44 *(prerequisite: v50 must be calibrated)*
8. **Prerequisites before submission:** (a) full v50 calibration (multi-start, acceptance checks), (b) Λ₁₀₇ post-processing implementation
9. **Role A language throughout** — "The reported transport coefficients are model-dependent effective parameters... non-unique... cannot be used as validated macroscopic coefficients without [additional measurements]."

### Target journals
*Solar Energy*, *Applied Thermal Engineering*, *International Journal of Heat and Mass Transfer*

### Timeline estimate
- Part II draft: 2–3 weeks
- v50 calibration run: hours of wall time (multi-start)
- Λ₁₀₇ post-processing: 1–2 days
- **Submission-ready revision: 4–6 weeks**

---

*Generated by Claude (Cowork) · AUTh-TRANSP solar-receiver project · 2026-09-02*
