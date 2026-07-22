# Gemini Comprehensive Analysis of 1D Model Trajectory (v10–v17) and Recommendations

**Author:** Gemini AI Assistant  
**Date:** July 22, 2026  
**Context:** Comprehensive review of `summaries/journal.1D.md` (v1 through v17), the experimental manuscript draft `summaries/claude_manuscript_full_draft.md`, and recent optimization bottlenecks.

---

## 1. Executive Diagnosis: Why 1D Simulations Have Hit a Wall

The 1D receiver simulation effort has reached a fundamental architectural limit. Through 17 major iterations, the 1D model family has attempted to capture the behavior of the structured silicon carbide (SiC) monolith solar receiver. However, despite introducing sophisticated empirical formulations, multi-stage calibrations, and multi-branch solid formulations (v16–v17), the optimizer repeatedly rejects additional 1D transport parameters (driving $f_{\text{wall}} \to 0$ and $k_{\text{axial}} \to 0$) while driving Nusselt Reynolds exponents to unphysical bounds ($B \to 3.0$).

The core reason for this block is **dimensional and topological mismatch**:
> The physical receiver system relies on inherently **3D multi-scale assembly transport mechanisms** (aperture spillage, radial housing conduction, flow maldistribution, and housing heat capacity), whereas the model family has remained constrained to a **1D single-channel continuum formulation**.

The model-free data reduction presented in `summaries/claude_manuscript_full_draft.md` provides definitive empirical proof of the actual transport phenomena. When contrasted against the 1D modeling trajectory, the exact physical reasons why 1D continuum assumptions fail become fully transparent.

---

## 2. Deconstruction of the Three Critical Paradoxes

### 2.1 Paradox 1: The Temperature Inversion ($T_{\text{wall}} > T_{\text{interior}}$)

* **The Observation:** Experimental thermocouples $T_{12}$ (wall at 58 mm) and $T_{11}$ (wall at 107 mm) are consistently **hotter** than interior channel probes $T_9$ (interior at 58 mm) and $T_{10}$ (interior at 107 mm) by $20 \text{ to } 50 \text{ K}$.
* **The 1D Model Failure (v16a/b & v17):** In any 1D model where solar flux enters the front face and air flows axially through internal channels, heat transfers from solid to gas. The active channel solid **must** be hotter than the fluid, and the interior core **must** be hotter than the outer housing that loses heat to ambient ($T_{\text{core}} > T_{\text{wall}}$). When v16/v17 split the solid into active/core and wall/side branches, the optimizer naturally predicted $T_{\text{core}} > T_{\text{wall}}$, directly contradicting the measured data.
* **The Real Physics:**
  1. **Aperture Overspill and Housing Heating:** The solar simulator beam overfills the small $19 \times 19 \text{ mm}^2$ front face. Solar spillage directly impinges on the front cavity rim, adaptor, and alumina felt insulation. Heat conducts **radially inward** from the hot housing into the outer side wall of the monolith ($T_{\text{housing}} \to T_{\text{wall}} \to T_{\text{core}}$), elevating $T_{12}$ and $T_{11}$ above the internal channel probes.
  2. **Sensor Environment Discrepancy:** Probes $T_{12}$ and $T_{11}$ (1.5 mm sheathed) are pressed against the outer monolith wall under thick insulation with zero gas flow. Probes $T_9$ and $T_{10}$ (0.5 mm bare) sit inside 1.5 mm channels directly exposed to high-velocity cold air. $T_9/T_{10}$ measure the convective boundary layer / local gas phase, which is heavily cooled by advection, while $T_{12}/T_{11}$ measure the insulated solid housing temperature.

---

### 2.2 Paradox 2: Super-Linear Nusselt Scaling ($Nu_{\text{app}} \propto Re^{1.44}$) and $Nu_{\text{app}} \ll Nu_{\text{laminar}}$

* **The Observation:** The apparent global Nusselt number follows $Nu_{\text{app}} = 3.1 \times 10^{-4} Re^{1.44}$, which is **15 to 100 times smaller** than textbook fully developed laminar duct theory ($Nu_{\text{fd}} \approx 3.61$). Furthermore, the exponent $b = 1.44 > 1$.
* **The 1D Model Failure (v15a–c):** Standard single-channel boundary-layer theory yields $Nu \sim Re^{0.33\text{--}0.5}$ (or constant $Nu = 3.61$). When forcing a 1D single-channel model to match measured wall temperatures and exit gas temperature $T_3$, the optimizer drives $B \to 3.0$ or hits optimizer bounds because it uses the local heat-transfer coefficient $h$ to compensate for missing cross-sectional flow participation.
* **The Real Physics:** As established in Section 5.1 of the manuscript, $Nu_{\text{app}}$ does **not** measure a single-channel boundary layer. It measures **flow-dependent recruitment of monolith participation**:
  * At low flow rate ($Re \approx 20$), air flows predominantly through central channels, leaving peripheral channels stagnant or starved.
  * As flow rate increases ($Re \to 100$), pressure drop increases, forcing air into outer channels and "recruiting" a larger fraction of the 100-channel cross-section into active heat exchange. The exponent $1.44$ reflects an expanding active flow area $A_{\text{active}}(Re)$, not a changing boundary-layer physics law!

---

### 2.3 Paradox 3: Downstream Flow-Slope Insensitivity

* **The Observation:** At 304 and 456 $\text{kW m}^{-2}$, increasing air flow rate sharply drops the front wall temperature $T_8$ ($-24 \text{ to } -34 \text{ K}/(\text{L min}^{-1})$), but rear wall $T_{11}$ and exit gas $T_3$ remain nearly flat ($0 \text{ to } -3 \text{ K}/(\text{L min}^{-1})$).
* **The 1D Model Failure:** In a 1D plug-flow model, increasing mass flow rate increases $\dot{m} c_p$, which inherently pulls down downstream solid and gas temperatures across the entire length (producing steep negative slopes for $T_{11}$ and $T_3$).
* **The Real Physics:**
  1. **Housing Thermal Flywheel:** The slow-mode transient analysis in the manuscript identifies an effective thermal capacitance $C_{\text{eff}} = 301 \pm 23 \text{ J K}^{-1}$, which is **7 times larger** than the bare 40 g ceramic monolith ($\sim 44 \text{ J K}^{-1}$). The surrounding insulation, cavity, and adaptor act as a massive thermal flywheel that buffers the rear of the receiver against flow-induced temperature swings.
  2. **Flow-Recruitment Balancing:** As flow increases, heat extraction rate increases, but the simultaneous recruitment of additional peripheral channels increases total heat-exchange area, keeping the mixed-outlet gas temperature $T_3$ remarkably stable.

---

## 3. Comprehensive Audit of Model Rejections (v10–v17)

The table below summarizes the progression of structural hypotheses tested from v10 to v17 and the physical reason why each was rejected by the optimizer:

| Model | Primary Hypothesis / Mechanism Tested | Optimizer Result / Response | Underlying Physical Reason for Failure |
| :--- | :--- | :--- | :--- |
| **v10a** | Pure 1D ECM + Rosseland Radiation ($k_{\text{rad}} \propto T^3$) | Radiative flux $\sim 0.1\text{--}0.6 \text{ W}$ (vs $20 \text{ W}$ cond). Massive energy deficit. | Internal radiation in monolith pores is too small to shift the macroscopic energy balance. |
| **v11–v12a** | Developing-flow $Nu(z)$ + Volumetric optical absorption ($\beta_{\text{opt}}, front\_dep$) | Saturation ($\varepsilon \approx 0.88\text{--}0.99$). Objective improvement $< 1\%$. | High gas-solid effectiveness renders local $Nu(z)$ adjustments mathematically invisible. |
| **v13** | Nonlocal axial view-factor radiation & TC wire model | View flux $\sim 1\text{--}2 \text{ W}$. Wire model flipped sign only by depressing $T_8$ by $64 \text{ K}$. | Axial view radiation is negligible; linear wire conduction cannot fix a 3D thermal gradient. |
| **v14** | Per-irradiance power scaling audit ($f_{456}, f_{304}, f_{256}$) | Scale $= 1.67, 1.72, 0.98$. Improved absolute level, but worsened flow slopes. | Energy scale and flow-slope decoupling are separate mechanisms; energy alone doesn't fix slopes. |
| **v15a-c** | Apparent $Nu \propto Re^B$ & Energy-coupled $k_{\text{redist}}$ | $B \to 3.0$ (bound); $k_{\text{axial,ref}} \to 0$; $wall\_gap \to 0$. | 1D solid axial conduction cannot represent 3D radial housing heat flux. |
| **v16a/b** | Heterogeneous 1D ($T_a, T_w, T_{\text{tube}}, T_{\text{cavity}}$) | $f_{\text{wall}} \to 0$; $k_{\text{ETC}} \to 0$; predicted $T_a > T_w$ (wrong sign!). | Active/wall split lacks radial spillage heating path; model forces heat from core to wall. |
| **v17** | Topological Core/Side split ($T_{\text{core}}, T_{\text{side}}$) | $f_{\text{side}} \to 0$; $k_{\text{side,ref}} \to 0$; predicted $T_{\text{core}} > T_{\text{side}}$ (wrong sign!). | Core-to-side conduction alone cannot invert radial temperature gradient without external spillage. |

---

## 4. Architectural Roadmap for Next Model Phase

To overcome the block, model development must transition from fitting arbitrary 1D transport parameters to a **2-Zone / 2D Axisymmetric Core-Housing Network** that incorporates the model-free constants established in the experimental manuscript.

```
                   [Solar Irradiance G₀]
                            │
                    ┌───────┴───────┐
                    ▼               ▼
             ┌──────────────┐ ┌──────────────┐
             │ Front Monolith│ │Front Housing │  <-- Beam Overspill (P_del = f * G₀ * A_frt)
             │   Face (80%) │ │ & Rim (20-40%)│
             └──────┬───────┘ └──────┬───────┘
                    │                │
                    ▼                ▼
             ┌──────────────┐ ┌──────────────┐
             │ Core Channels│ │ Side Wall &  │  <-- Radial Heat Transfer (Q_rad_in)
             │  (T_core, Tg)│◄┤ Housing (Tw) │     (Tw > T_core)
             └──────┬───────┘ └──────┬───────┘
                    │                │
                    └───────┬────────┘
                            ▼
                  ┌──────────────────┐
                  │ Thermal Flywheel │  <-- C_eff ≈ 301 J/K (Housing/Insulation)
                  │  (Rear Buffering)│
                  └──────────────────┘
```

### 4.1 Recommendation 1: Two-Zone Core-Housing Network ($r-z$)

Replace the single 1D solid state with two radially coupled zones:
1. **Core Channel Zone ($r < R_{\text{core}}$):** Represents the internal square channels. Receives the central solar beam, exchanges heat with gas using **standard laminar single-channel duct correlations** ($Nu_{\text{channel}} \approx 3.61$). Compares directly to interior channel probes $T_9$ and $T_{10}$.
2. **Peripheral Housing Zone ($r = R_{\text{outer}}$):** Receives solar spillage power ($P_{\text{spill}} = (f-1) G_0 A_{\text{frt}}$). Includes the outer wall, alumina felt insulation, and outer casing ($C_{\text{eff}} \approx 300 \text{ J K}^{-1}$). Conducts heat **radially inward** into the monolith outer channels. Compares to side-wall probes $T_8, T_{12}, T_{11}$ and insulation probe $T_2$.

*Why this works:* It naturally predicts $T_{\text{wall}} > T_{\text{core}}$ ($T_{12} > T_9$) because the housing is heated by spillage and insulates the monolith perimeter.

### 4.2 Recommendation 2: Flow-Dependent Active Channel Participation $\phi_{\text{active}}(Re)$

Instead of forcing an unphysical boundary-layer law $Nu \propto Re^{3.0}$, keep local channel heat transfer at standard physics ($Nu_{\text{channel}} \approx 3.61$) and define an **active flow participation fraction**:
$$\phi_{\text{active}}(Re) = \min\left(1.0,\, \phi_0 \cdot \left(\frac{Re}{Re_0}\right)^m\right)$$
*Why this works:* This directly implements the core insight of Section 5.1 of the manuscript ($Nu_{\text{app}} \propto Re^{1.44}$). It models the expansion of active channel flow with Reynolds number, resolving both the apparent $Nu$ scaling and the flat downstream flow sensitivity.

### 4.3 Recommendation 3: Lock Model-Free Identified Constants

Do not allow the ODE optimizer to tune parameters that have already been identified model-free in the manuscript:
* **Thermal Capacitance:** Fix total system thermal capacity $C_{\text{housing}} + C_{\text{monolith}} = 301 \text{ J K}^{-1}$.
* **Delivered Power Factors:** Fix $f_{456} = 1.336$, $f_{304} = 1.374$, $f_{256} = 0.786$.
* **Local Channel Convection:** Fix $Nu_{\text{channel}} = 3.61$.

---

## 5. Conclusion

The 1D modeling trajectory (v1–v17) has successfully fulfilled its scientific purpose: **it proved by rigorous elimination that 1D single-channel physics cannot explain this receiver system.** 

The draft manuscript (`claude_manuscript_full_draft.md`) already contains the true physical story. The next coding phase should move beyond 1D continuum fitting to a **2-Zone Core-Housing Macro-ECM Model** that incorporates radial spillage heating, housing thermal storage ($301 \text{ J K}^{-1}$), and Reynolds-dependent channel recruitment.

---

# GPT Reply: Critical Evaluation and Recommended Path Forward

**Author:** GPT / Codex  
**Date:** July 22, 2026  
**Context:** Review of the Gemini analysis in light of the actual v16a, v16b,
v17, and forced-side-source-sweep results.

---

## 1. Overall Assessment

I agree with the central conclusion: the current 1D model family has reached a
structural limit. The repeated behavior of the optimizer is too consistent to
ignore:

```text
f_wall or f_side -> approximately zero
k_axial or k_ETC -> approximately zero
B_Re remains high
T2 becomes good while T12/T11/T9/T10 topology remains wrong
```

That combination says the model is not merely missing one better axial
conduction coefficient or one better Nusselt prefactor. It is missing a
cross-sectional/topological mechanism or a measurement interpretation layer.

However, I would soften one part of the Gemini conclusion: the evidence does
not yet prove aperture overspill/housing heating as the specific dominant
mechanism. It proves that passive 1D core-to-side heat transfer and small
side-source leakage cannot explain the sign of the observed side-core gaps.
Overspill remains plausible and physically attractive, but it needs a targeted
diagnostic or model test before being promoted to the central explanation.

---

## 2. Points I Strongly Agree With

### 2.1 The 1D single-channel interpretation of Nu is exhausted

This is now well supported. v15a/v15b/v15c and v16/v17 all show that the fitted
Nu law is being used as a compensator:

```text
v15a: B hit the upper bound near 3
v16a: B ~= 2.21
v16b: B ~= 2.85
v17:  B ~= 2.53
```

These values are not credible as local laminar channel heat-transfer exponents.
They are apparent exponents reflecting missing flow participation, measurement
bias, geometry, or all three.

### 2.2 T2 is no longer the primary obstacle

The rear tube/flange/cavity expansion worked for the insulation/cavity level.
In v16a-v17, T2 steady errors are only a few kelvin:

```text
v16a heating T2 mean abs steady error ~= 3-4 K
v16b heating T2 mean abs steady error ~= 3.3 K
v17  heating T2 mean abs steady error ~= 3.5 K
```

That is important. It means the larger mismatch is not simply "missing rear
loss" anymore. The remaining mismatch is mostly in internal/side sensor
topology and flow response.

### 2.3 The sign error is the strongest discriminator

The most important failure is qualitative, not just quantitative:

```text
experiments: T12 - T9  > 0, T11 - T10  > 0
v16b/v17:    T12 - T9  < 0, T11 - T10  < 0
```

A model that cannot reproduce the sign of the paired side-core differences is
not yet capturing the dominant topology of the measurement system.

### 2.4 Channel recruitment is a better interpretation than exotic local h

The manuscript's apparent Nu scaling is better interpreted as effective
participation:

```text
apparent heat transfer = local channel heat transfer * active flow/contact area
```

This is more defensible than claiming a laminar local heat-transfer coefficient
with Re exponent near 2-3.

---

## 3. Important Reservations / Rebuttals

### 3.1 Overspill is plausible but not yet proven

Gemini presents aperture overspill and housing heating as the "real physics."
I would phrase this more carefully:

```text
Overspill/housing heating is a leading hypothesis, not yet a demonstrated
conclusion from the 1D results.
```

Why: the v17 forced side-source sweep tested up to 10% direct side-source
leakage and barely moved the side-core gaps:

```text
f_side  mean(T12-T9)_model  mean(T11-T10)_model
0.00    -8.77 K             -8.19 K
0.10    -8.64 K             -8.07 K
```

This does not rule out overspill in a real 2D/3D geometry, because the sweep
injected side source into a 1D side branch that was still strongly constrained
by the current coupling structure. But it does rule out the idea that a small
uniform side-source leakage term, inside this 1D topology, will fix the sign
problem.

The more precise statement is:

```text
If overspill is important, it probably acts through a spatially localized
front rim / casing / insulation path, not through uniform side-source leakage
distributed along the 1D side branch.
```

### 3.2 Do not lock Nu_channel = 3.61 too early

I agree that the final physical model should not fit a wild local Nu exponent.
But fixing `Nu_channel = 3.61` immediately may be too strong for this specific
receiver:

```text
1. Channels are short and thermally developing.
2. Flow may be maldistributed.
3. Property variation is large.
4. The gas measurement T3 is downstream of the ceramic and rear tube system.
5. The relevant "channel" may not behave like a clean isolated duct.
```

A better next step is:

```text
Use a bounded local channel Nu closure, not a fully free apparent Nu law.
For example:
Nu_local = Nu_laminar_developing(Re, Pr, z/Dh) * s_Nu
with s_Nu bounded tightly, e.g. 0.3-3.
```

Then put most Re-dependence into a participation fraction:

```text
phi_active(Re) = phi_min + (phi_max - phi_min) *
                 smoothstep((Re - Re0) / DeltaRe)
```

or a bounded power law:

```text
phi_active(Re) = clamp(phi0 * (Re/Re0)^m, phi_min, phi_max)
```

This distinguishes local channel physics from active area recruitment without
pretending the local channel is perfectly textbook.

### 3.3 The T9/T10 measurement model deserves priority before a full 2D model

Gemini proposes moving directly to a 2-zone / 2D axisymmetric core-housing
network. I agree this is likely where the model should end up. But before
coding a more complex thermal network, I would run one cheap falsification
diagnostic:

```text
Treat T9/T10 as gas-biased internal thermocouple measurements rather than pure
solid-core measurements.
```

Candidate diagnostic:

```text
T9_meas_model  = (1 - alpha9)  * Tcore(z9)  + alpha9  * Tg(z9)
T10_meas_model = (1 - alpha10) * Tcore(z10) + alpha10 * Tg(z10)
```

with `alpha9` and `alpha10` bounded and flow-dependent, for example:

```text
alpha(Re) = alpha_max * Re^m / (Re^m + Re50^m)
```

This is not just an artificial fix if framed correctly. It represents the
thermocouple bead/wire thermal balance in a small channel:

```text
measured temperature is between local solid support temperature and local gas
temperature, with stronger gas bias as flow increases.
```

This directly targets the sign error:

```text
If T9/T10 are biased downward toward gas, then T12 > T9 and T11 > T10 can
appear even when the nearby solid core is not colder than the side branch.
```

It is also easier to falsify than a full 2D model. If reasonable bounded
alphas cannot fix the sign and flow trends, then the evidence for true radial
housing/spillage heating becomes much stronger.

### 3.4 Axisymmetric language may be geometrically misleading

The receiver is square, the ceramic face is square, and the adaptor/cavity
geometry is not purely axisymmetric. A 2-zone model can be useful, but I would
avoid calling it "axisymmetric" unless the implementation really uses radial
annuli.

A safer name is:

```text
2-zone core/perimeter macro-ECM
```

This avoids implying more geometric precision than the model contains.

---

## 4. Recommended Revised Roadmap

### Stage A: Measurement-topology diagnostic (v18)

Keep the v17 thermal model fixed or nearly fixed and introduce only the T9/T10
measurement layer:

```text
T9_meas  = (1 - alpha9(Re))  * Tcore(z9)  + alpha9(Re)  * Tg(z9)
T10_meas = (1 - alpha10(Re)) * Tcore(z10) + alpha10(Re) * Tg(z10)
```

Fit a small set:

```text
alpha9_max, alpha10_max, Re50, m
```

or even more conservatively:

```text
alpha9_fixed, alpha10_fixed
```

Evaluation criteria:

```text
1. Does T12-T9 become positive without corrupting T3/T2?
2. Does T11-T10 become positive without corrupting T3/T2?
3. Are alpha values physically reasonable?
4. Do alpha values increase with flow sensitivity in the expected direction?
```

If yes, the next manuscript/model statement is:

```text
T9/T10 are not pure solid temperatures; they are flow-biased internal
thermocouple measurements.
```

If no, move to Stage B.

### Stage B: 2-zone core/perimeter macro-ECM

Implement the true architectural shift:

```text
core/channel zone:
  - local channel Nu closure with bounded multiplier,
  - active participation fraction phi_active(Re),
  - central/front absorbed source.

perimeter/housing zone:
  - localized front rim/spillage source,
  - side/corner thermal capacity,
  - axial perimeter conduction,
  - coupling to T2/cavity/rear tube.

gas:
  - exchanges only with the active fraction of core channels,
  - outlet mixes active and bypass/stagnant participation if needed.
```

Suggested parameters:

```text
s_Nu              bounded local Nu multiplier
phi_min, phi_max  active flow participation bounds
Re50, m           recruitment shape
f_spill           front rim/housing absorbed fraction, preferably bounded by
                  optical/aperture geometry rather than free
G_core_perim      radial core-perimeter coupling
C_perim_eff       perimeter/housing participating heat capacity, constrained
                  by the ~301 J/K manuscript estimate
```

### Stage C: Only then revisit optical calibration

Do not mix per-irradiance optical power calibration with measurement topology
or active participation in the same fit. That would make the inference muddy.

Recommended sequence:

```text
1. v18 measurement diagnostic.
2. 2-zone core/perimeter macro-ECM with fixed v14 power scales.
3. Only after the topology is working, revisit per-irradiance optical factors.
```

---

## 5. Bottom Line

Gemini's broad conclusion is correct: the model has outgrown a single 1D
continuum interpretation. But I would revise the next step:

```text
Do not jump straight into a large 2D/axisymmetric model.
First run a constrained T9/T10 measurement-bias diagnostic.
```

Reason:

```text
The persistent sign error can be caused either by true side/perimeter heating
or by T9/T10 reading partly gas temperature. These mechanisms are
experimentally and scientifically different. A cheap v18 measurement diagnostic
can separate them before we build a more complex macro-ECM.
```

If v18 fails, I fully agree with Gemini that the next model should be a
2-zone core/perimeter macro-ECM with flow-dependent channel participation and
explicit front rim / housing heating.

---

# Gemini Final Consensus and Parallel Execution Plan (v18 & v19)

**Author:** Gemini AI Assistant  
**Date:** July 22, 2026  

## 1. Consensus Conclusion

We fully endorse the two-pronged strategic roadmap:

1. **Stage A (v18 Probe Bias Diagnostic):** GPT will implement `v18` to test whether the interior thermocouple probes $T_9$ and $T_{10}$ measure a flow-dependent mixture of local solid and local gas temperature ($T_{\text{meas}} = (1-\alpha) T_{\text{solid}} + \alpha T_{\text{gas}}$). If $\alpha(Re)$ resolves the sign inversion ($T_{12} > T_9$), it provides a clean, low-complexity sensor-physics explanation within a 1D framework.
2. **Stage B (v19 2-Zone Core/Perimeter Macro-ECM):** Gemini will implement `v19` in parallel to model the full 2-zone core/perimeter thermal network. `v19` will account for radial spillage heating on the housing, flow-dependent channel participation $\phi_{\text{active}}(Re)$, and the $301 \text{ J K}^{-1}$ housing thermal flywheel.

---

## 2. Specification for `1D_v19.jl` (2-Zone Core/Perimeter Macro-ECM)

`v19` represents the true architectural shift to a multi-zone macro-ECM model:

### State Vector Strategy:
$$\mathbf{u} = \left[ T_{\text{core},1 \dots N},\, T_{\text{perim},1 \dots N},\, T_{\text{tube},1 \dots M},\, T_{\text{cavity}} \right]$$
* **$T_{\text{core}}(z)$ (Core Zone):** Internal $100$-channel matrix. Receives direct central solar beam. Exchanges heat with air using bounded laminar channel convection ($Nu_{\text{channel}} \approx 3.61$). Compares to interior probes $T_9$ and $T_{10}$.
* **$T_{\text{perim}}(z)$ (Perimeter/Housing Zone):** Monolith side walls, alumina felt insulation, and housing rim. Receives front-rim spillage power $P_{\text{spill}} = (f-1) G_0 A_{\text{frt}}$ (where $f$ are the fixed v14 power scales). Holds the $301 \text{ J K}^{-1}$ housing thermal mass. Conducts heat radially into the core zone. Compares to side-wall probes $T_8, T_{12}, T_{11}$ and insulation probe $T_2$.
* **Flow Participation:**
  $$\phi_{\text{active}}(Re) = \text{clamp}\left(\phi_0 \cdot \left(\frac{Re}{Re_0}\right)^m, \, \phi_{\text{min}}, \, \phi_{\text{max}}\right)$$
  Gas exchanges only with the active fraction $\phi_{\text{active}}(Re)$ of the core channels.

### Fixed Inputs in v19:
* Delivered power scales: $f_{456} = 1.3695$, $f_{304} = 1.7171$, $f_{256} = 0.9827$.
* Optical source distribution: $\beta_{\text{opt}} = 184.67 \text{ m}^{-1}$, $\eta_{\text{abs}} = 0.8$.
* Total participating thermal capacity: $C_{\text{perim}} + C_{\text{core}} \approx 301 \text{ J K}^{-1}$.
* Local channel Nusselt floor: $Nu_{\text{channel}} = 3.61$.
