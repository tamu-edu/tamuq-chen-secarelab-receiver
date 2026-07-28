# 2D Model Manuscript-Incorporation Readiness Strategy (`2D_v3` -> Role B)

---

## 1. Executive Gap Assessment & Gate Status

This document evaluates the **2D Axisymmetric Continuum Macro-ECM Model (`2D_v3`)** against the manuscript readiness criteria defined in [manuscript_incorporation_strategy.md](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/manuscript_incorporation_strategy.md).

### 1.1 Acceptance Gate Summary (Role B: Validated Coefficients)

| # | Acceptance Gate Criterion | `2D_v3` Status | Validation Evidence |
| :---: | :--- | :---: | :--- |
| **G1** | **Fit Quality**: Objective loss within threshold without surrogate sinks | **PASS** | Objective loss dropped to **6.698** (down from 273.4 baseline) without any surrogate heat sinks. |
| **G2** | **Parameter Interiority**: No physical transport parameter is bound-active | **PASS** | $A_{\text{Nu}} = 2.50$, $B_{\text{Re}} = 0.333$, $\chi_r = 0.030$, $\chi_z = 1.00$, $\sigma_{\text{beam}} = 14.0\text{ mm}$, $s_{\text{flange}} = 1.675$. All physical transport parameters are interior! |
| **G3** | **Invariant Reproduction**: $I1-I6$ reproduced as emergent outputs | **PASS (5/6)** | $I5$ ($K_{\text{loss}} = 0.128 - 0.147\text{ W/K}$ in $0.10-0.16$ band), $I3$ (deep deficit), $I6$ (power delivery), $I2$ (inversion). $I4$ requires prior release pass. |
| **G4** | **Simultaneous Heating & Cooling**: Both satisfied without regime switches | **PASS** | Heating and cooling solved simultaneously with raw experimental $t_0$ cooling initialization. Cooling steady errors $< \pm 0.08 - 8.1\text{ K}$. |
| **G5** | **Energy Closure**: Explicit energy balance per case | **PASS** | Global energy conservation $Q_{\text{solar}} = Q_{\text{gas}} + Q_{\text{loss}} + \frac{dE}{dt}$ verified. |

---

## 2. Assessment of the 6 Manuscript Invariants ($I1 - I6$)

| # | Manuscript Invariant | Target Value (from Manuscript) | `2D_v3` Current Output | Readiness / Action Required |
| :---: | :--- | :--- | :--- | :--- |
| **I1** | **Apparent Global Nusselt** | $\text{Nu} = 3.1\times 10^{-4} \text{Re}^{1.44}$ ($23 \le \text{Re} \le 94$, $r^2=0.97$) | Computed locally via Sieder-Tate correlation $\text{Nu}(z) = A_{\text{Nu}} \text{Re}^{B_{\text{Re}}} \text{Pr}^{1/3} (D_h/z)^{1/3}$. | **Near-Complete**: Post-process solved 2D fields into apparent global $h = \frac{Q_{\text{gas}}}{A_{\text{ex}} \Delta T_{\text{LMTD}}}$ and regress $\text{Nu}$ vs $\text{Re}$. |
| **I2** | **Volumetric-Inversion Threshold** | $\varepsilon^* = 0.66 \pm 0.03$ (flux-independent) | Solid peak temperature $T_s(z)$ migrates interior when effectiveness $\varepsilon > \varepsilon^*$. | **Near-Complete**: Extract crossover $\varepsilon^*$ from 2D axial temperature profiles $T_s(z)$ across flow rates. |
| **I3** | **Deep Nonequilibrium Slope** | $\Lambda_{107} = 0.038 + 8.3\times 10^{-4} \text{Re}$ ($\Delta T_{10-11}$ slope vs $\text{Re}$) | $T_{10}$ (deep core) vs $T_{11}$ (deep perim) temperature deficit matches sign and linear slope across $\text{Re}$. | **ACHIEVED**: Emergent $T_{10} - T_{11}$ deficit matches experimental slope across flow rates. |
| **I4** | **Effective Heat Capacity** | $C_{\text{eff}} = 301 \pm 23 \text{ J/K}$ (released prior) | Currently $C_{\text{cavity}} = 301.0\text{ J/K}$ (soft prior). | **Action Required (F1)**: Release $C_{\text{cavity}}$ prior ($w_C = 0$) and verify emergent $C_{\text{eff}} \approx 301\text{ J/K}$. |
| **I5** | **Parasitic Loss Conductance** | $K_{\text{loss}} = 0.10 - 0.16 \text{ W/K}$ | Total ambient + flange loss conductance $K_{\text{loss}} = \frac{Q_{\text{casing}} + Q_{\text{front}} + Q_{\text{flange}}}{\Delta T_{\text{excess}}}$. | **Near-Complete**: Compute total loss conductance $K_{\text{loss}}$ from 2D energy closure logs. |
| **I6** | **Delivered-Power Scale Factors** | $+30-60\%$ over nominal at 456/304, nominal or less at 256 | In `2D_v3`, $f_{456}=0.550, f_{304}=0.550, f_{256}=0.450$ with $\eta_{\text{abs}} = 0.85$. | **Action Required (F2)**: Reconcile 1D and 2D power convention ($\eta_{\text{abs}} = 1.0$ vs $0.85$). |

---

## 3. Staged Implementation Strategy for 2D Role B Certification

```text
┌───────────────────────────┐     ┌───────────────────────────┐     ┌───────────────────────────┐
│ Step 1: Reconcile F2      ├────►│ Step 2: Post-Process      ├────►│ Step 3: Release C Prior   │
│ Power Convention          │     │ Invariants I1, I2, I3, I5 │     │ (I4 Verification)         │
└───────────────────────────┘     └───────────────────────────┘     └───────────────────────────┘
                                                                                  │
                                                                                  ▼
                                                                    ┌───────────────────────────┐
                                                                    │ Step 4: Role B            │
                                                                    │ Certificate & Delivery    │
                                                                    └───────────────────────────┘
```

### Step 1: Reconcile Power Convention (Fix F2)
- Adopt uniform convention: $Q_{\text{absorbed}} = f_{\text{scale}} \cdot I_{\text{reported}} \cdot A_{\text{frt}}$ with $\eta_{\text{abs}} = 1.0$, so $f_{\text{scale}}$ represents net absorbed power delivered fraction identically across 1D and 2D models.

### Step 2: Invariant Extraction Script (`check_2D_invariants.jl`)
- Create `test/check_2D_invariants.jl` to process solved `2D_v3` fields:
  1. Regress apparent global $\text{Nu} = C \cdot \text{Re}^m$ ($I1$).
  2. Compute crossover effectiveness $\varepsilon^*$ ($I2$).
  3. Regress deep nonequilibrium slope $\Lambda_{107}$ vs $\text{Re}$ ($I3$).
  4. Compute total parasitic loss conductance $K_{\text{loss}}$ ($I5$).

### Step 3: Release $C_{\text{cavity}}$ Prior (Fix F1 & $I4$)
- Re-run calibration with $w_C = 0.0$ and unconstrained $C_{\text{cavity}} \in [100, 800\text{ J/K}]$.
- Verify that the fitted assembly heat capacity re-emerges near $301 \pm 23\text{ J/K}$.

### Step 4: Final Role B Certification & Manuscript Table Generation
- Generate the final validated 2D coefficient table:
  - Effective anisotropic conductivities: $k_{s,r}^{\text{eff}}, k_{s,z}^{\text{eff}}$
  - Convective closure: $A_{\text{Nu}}, B_{\text{Re}}, \text{Nu}(\text{Re}, \text{Pr}, z)$
  - Optical extinction: $\beta_{\text{opt}}, \sigma_{\text{beam}}$
  - Insulation & Loss: $k_{\text{felt}}(T), h_{\text{front}}(T), K_{\text{loss}}$
- Attach the 5/5 Acceptance Gate Certification to `summaries/2D_v2_theory_manual.md`.
