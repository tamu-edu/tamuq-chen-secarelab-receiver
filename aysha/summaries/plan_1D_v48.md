# Implementation Plan: Entire Converter Model `1D_v48`

Refactor the Entire Converter Model (ECM) into version **`1D_v48`** to resolve all technical critiques from `summaries/critique_1D_v47_2026-09-01.md`, capture the experimental rear-temperature flow-invariance trend via physics-based mechanisms, prove mesh convergence across $N \in \{15, 25, 50\}$ nodes, and perform full Option A calibration and validation.

---

## User Review & Key Design Decisions

> [!NOTE]
> **Key Physical Modeling Decisions for `1D_v48`**:
> 1. **Developing Nusselt Formulation**: Replace the compound double-Graetz $(Gz) \times (Gz)^B$ formula with a standard Shah-London / Hausen laminar developing entry law:
>    $$\text{Nu}(z) = \text{Nu}_\infty + \frac{C_1 \text{Gz}_z}{1 + C_2 \text{Gz}_z^{2/3}}$$
>    with asymptotic theoretical limit $\text{Nu}_\infty = 3.61$ for fully developed square ducts.
> 2. **Rear Axial Conductance Dimensional Fix**: Formulate rear rail inter-node conduction as $\frac{(kA)_{rear, eff}}{\Delta z} (\Delta T)$ in $[\text{W}\cdot\text{m/K}]$ to guarantee dimensional consistency and grid-independence.
> 3. **Decoupled Radiative Mechanisms**:
>    - Solar Band ($0.2\text{--}2.5\ \mu\text{m}$): Optical penetration depth $\beta_{opt}$ with front absorption fraction $f_{front}$.
>    - Thermal IR Band ($2.5\text{--}50\ \mu\text{m}$): Rosseland internal diffusion $\beta_{rad}$ + **Direct Line-of-Sight Cavity Radiation Beaming** ($Q_{rad, LoS}$) from the front hot spot ($1200\text{ K}$) to the rear exit plane ($700\text{ K}$).
> 4. **Downstream Gas Recuperation & Parasitic Drain Reduction**: Reduce static parasitic losses from the rear rail to the cold cavity ($G_{rear-cavity} \to 0$), allowing downstream gas-to-solid reheating ($Q_{gas} < 0$) to buffer and maintain rear solid temperatures at high flow rates.
> 5. **Fixed Water Flange Boundary**: Retain the fixed liquid cooling sink ($T_{water} = 293.15\text{ K}$) without artificial empirical flow-dependent resistance.
> 6. **Grid-Invariant Macroscopic HTC**: Define the bulk macroscopic heat transfer coefficient via LMTD:
>    $$\bar{h}_{macro} = \frac{\dot{Q}_{gas, core}}{P_{exchange} L \cdot \Delta T_{LMTD}}$$
>    and verify its numerical convergence across $N \in \{15, 25, 50\}$ nodes in the pre-calibration smoke test.

---

## Proposed Changes

### Core Physics & Simulation Engine

#### [`1D_v48.jl`](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/1D_v48.jl)
- **Geometry & Material Synchronization**: Consolidate single authoritative geometry parameters and polynomial material properties for SiC and air.
- **Shah-London Nusselt Entry**: Implement single-group entry law with $\text{Nu}_\infty = 3.61$.
- **Decoupled Optical & Thermal Radiation**: Separate $\beta_{opt}$ (solar band) from $\beta_{rad}$ (IR diffusion) and add direct line-of-sight radiative transfer $Q_{rad, LoS} = F_{1 \to N} A_{frt} \sigma_{SB} (T_{core, 1}^4 - T_{core, N}^4)$.
- **Dimensionally Consistent Rear Rail**: Implement $(kA)_{rear, eff} / \Delta z$ conduction.
- **Constrained Parasitic Drains**: Remove unphysical parasitic leakages from the rear rail to the cold cavity.
- **Macroscopic LMTD Extraction**: Implement `macroscopic_htc_v48(result, idx)`.
- **First-Law Energy Audit**: Implement exact `compute_energy_balance_v48` with $< 10^{-13}\text{ W}$ instantaneous residual.

#### [`test/smoke_1D_v48.jl`](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/test/smoke_1D_v48.jl)
- Basic parameter consistency and bounds check.
- **Automated Mesh Sensitivity Suite**: Forward integration of E67 and E72 at $N \in \{15, 25, 50\}$ nodes confirming temperature profile convergence within $< 2\%$ and $\bar{h}_{macro}$ stability within $< 3\%$.
- Exact First-Law energy conservation verification ($|\Delta \dot{E}_{inst}| < 10^{-13}\text{ W}$).
- Cooling case verification for C69, C80 (0 LPM natural cooling), and C81.

#### [`run_1D_v48.jl`](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/run_1D_v48.jl)
- Execution runner with NLopt BOBYQA optimization on the 15 heating runs under Option A.
- Automated generation and export of all CSV data tables and full plot catalogs in `summaries/1D_v48/`.
- Decaying cooling $t_{90}$ calculation fix.

---

### Documentation & Summaries

#### [`summaries/plan_1D_v48.md`](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/plan_1D_v48.md)
- Export of the approved roadmap document in the project repository.

#### [`summaries/theory_1D_v48.md`](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/theory_1D_v48.md)
- Comprehensive mathematical manuscript with LaTeX derivations, tables matching calibrated parameters, mesh invariance proofs, and physical discussion.

#### [`summaries/journal.1D.md`](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/journal.1D.md)
- Document the `1D_v48` formulation, calibration outcomes, mesh convergence metrics, and validation conclusions.

---

## Verification Plan

### Automated Tests
```powershell
# 1. Run Automated Smoke Test (including N=15, 25, 50 mesh convergence and First-Law balance)
& "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. test\smoke_1D_v48.jl

# 2. Run Full Calibration & Post-Processing Pipeline
$env:RECEIVER1D_v48_RUN_CALIBRATION="true"; $env:RECEIVER1D_v48_FIT_SECONDS="600.0"; & "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. run_1D_v48.jl
```

### Verification Criteria
1. Mesh sensitivity test in `test/smoke_1D_v48.jl` passes with $\bar{h}_{macro}$ varying by $< 3\%$ across $N \in \{15, 25, 50\}$.
2. Instantaneous First-Law conservation residual $< 10^{-13}\text{ W}$ across all cases.
3. Flow slope analysis demonstrates significant flattening of rear temperature slopes $(\partial T_{10}/\partial \dot{V}, \partial T_{11}/\partial \dot{V}, \partial T_3/\partial \dot{V})$.
4. All 58 output artifacts and plots generated in `summaries/1D_v48/`.
5. 100% synchronization of formulas, parameters, and tables across code, CSVs, and manuscript.
