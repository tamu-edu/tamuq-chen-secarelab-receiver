# Solar Receiver Validation: Quantitative Metrics Framework

This document defines the quantitative metrics and data extraction strategy used to evaluate simulation runs against experimental data for the solar receiver validation project.

## 1. Primary Physical Metrics

### 1.1 Volumetric Inversion Index ($I_{vol}$)
*   **Definition:** $I_{vol} = T_{9,ss} - T_{8,ss}$
*   **Target:** $I_{vol} > 0$
*   **Physical Meaning:** Quantifies the success of moving the thermal peak from the aperture ($z \approx 0$) deeper into the receiver monolith ($z \approx 50\text{mm}$). A positive value indicates a successful "volumetric effect."
*   **Control Parameters:** `source_y` (Ray focus), `VF_l_f` (Aperture cooling), `k_air_mult` (Inlet quenching).
*   **Data Availability:**
    *   S368–S393: T8 (y=11mm) **not exported** → proxy metric used: $T9(y{=}58) - T10(y{=}107)$.
    *   S396+: T8 export **fixed** (via `data3` → `_T08.csv`). True $I_{vol}$ now computable.

### 1.2 Steady-State Magnitude Error ($\Delta T_{SS}$)
*   **Definition:** $\Delta T_i = T_{i,sim} - T_{i,exp}$ for $i \in \{T_3, T_9, T_2\}$
*   **Success Criteria:** $|\Delta T_i| < 25 	ext{ K}$
*   **Physical Meaning:**
    *   $\Delta T_3$ (Gas Outlet): Overall system energy balance and coupling efficiency.
    *   $\Delta T_9$ (Middle Solid): Peak absorption and solar-to-thermal conversion magnitude.
    *   $\Delta T_2$ (Insulation): Radial heat loss accuracy.
*   **Control Parameters:** `I_f_high/low` (Power scaling), `qlpm_f_all` (Leakage surrogate).

### 1.3 Thermal Lag Indices ($\Gamma$)
*   **Definitions:**
    *   **Solid-core lag:** $\Delta t_{90\%,T9} = t_{90,sim}(T9) - t_{90,exp}(T9)$
    *   **Gas-exit lag:** $\Delta t_{90\%,T3} = t_{90,sim}(T3) - t_{90,exp}(T3)$
*   **Targets:** $\Delta t_{90\%,T9} \approx 0$ and $\Delta t_{90\%,T3} \approx 0$
*   **Physical Meaning:**
    *   $\Delta t_{90\%,T9}$ tracks the solid-core heating timescale.
    *   $\Delta t_{90\%,T3}$ tracks the gas-exit heating timescale and exposes gas-side transport / leakage mismatch more directly.
*   **Control Parameters:** `ins_cp` (Thermal mass), `tc_d` (Contact resistance), plus gas-side surrogates such as `qlpm_f_all` and `k_air_mult` for `T3`.

### 1.4 Longitudinal Gradient Slope ($m_{grad}$)
*   **Definition:** $m_{grad} = (T_{10} - T_9) / \Delta z$
*   **Target:** Match experimental decay slope.
*   **Physical Meaning:** Evaluates the accuracy of internal transport mechanisms.
*   **Control Parameters:** `kz_mult` (Radiative spreading).

### 1.5 Energy Partitioning Ratio ($R_{leak}$)
*   **Definition:** $R_{leak} = (T_3 - T_{amb}) / (T_9 - T_{amb})$
*   **Target:** Derived from experimental averages across flux levels.
*   **Physical Meaning:** Proxy for flow-to-leakage ratio. High values indicate energy is successfully being recovered by the gas stream.
*   **Control Parameters:** `qlpm_f_all` (Leakage factor).

### 1.6 Radial Nonequilibrium Gap Metrics ($\Delta T_{cw}$)
*   **Definitions:**
    *   **Mid-receiver center/wall gap:** $\Delta T_{cw,58} = (T_{12} - T_9)$
    *   **Back-receiver center/wall gap:** $\Delta T_{cw,107} = (T_{11} - T_{10})$
    *   **Gap-growth metric:** $\Delta T_{grow} = \Delta T_{cw,107} - \Delta T_{cw,58}$
*   **Targets:** Match experimental center-vs-wall gaps and, especially, match the increase in that gap from `58 mm` to `107 mm`.
*   **Physical Meaning:**
    *   Quantifies how strongly the receiver remains radially nonequilibrated.
    *   A larger downstream gap indicates that the interior remains thermally distinct from the external wall deeper in the receiver.
    *   This metric is diagnostic of whether the model is over-forcing gas/solid/wall equilibration too early.
*   **Model-side implementation note:** experimental center thermocouples are not assumed to be pure gas measurements. The current script therefore uses the exported **gas-center length profile** sampled at `58 mm` and `107 mm` as the model analog for the center measurement, compared against wall probes `T09` and `T10`.
*   **Reported columns in `analysis_results.csv`:**
    *   `Gap_58_exp_K`, `Gap_107_exp_K`, `Gap_growth_exp_K`
    *   `Gap_58_sim_K`, `Gap_107_sim_K`, `Gap_growth_sim_K`
    *   `dGap_58_K`, `dGap_107_K`, `dGap_growth_K` where `dGap = sim - exp`
*   **Control Parameters:** `k_air_mult` (gas/solid coupling surrogate), receiver transport parameters such as `kz_mult`, and receiver property surrogates such as `k_SiC_scale` / `Cp_SiC_scale`.

---

## 2. Data Extraction Strategy

### 2.1 Simulation Data (`outdata/`)
*   **Format:** Multi-case CSV files (e.g., `S<ID>_T08.csv`).
*   **Logic:**
    1.  Identify case boundaries (E67, E68, E78) using `t=0` markers.
    2.  Extract final value of each case as the Steady-State ($SS$) value.
    3.  Calculate $t_{90\%}$ by finding the first timestamp where $T \ge T_{initial} + 0.9(T_{ss} - T_{initial})$.
    4.  Store explicit lag columns for both `T09` and `T03`, while preserving the legacy `dt90_s` field as the `T09`-based lag alias for backward compatibility.

### 2.2 Experimental Data (`exp/`)
*   **Requirement:** Establish baseline targets for each metric.
*   **Source:** Read corresponding `exp/E<ID>_*.csv` files to calculate the experimental $I_{vol}$, $\Gamma$, and $R_{leak}$ for comparison.

---

## 3. Implementation Workflow

1.  ✅ **Baseline Extraction:** Implemented in `analyze_batch.py` — reads `exp/E##_steady_state.csv` and `exp/E##_transient.csv` for reference values.
2.  ✅ **Simulation Processing:** Implemented in `analyze_batch.py` — parses `outdata/S###_T##.csv` for all discovered runs, parses `S###_solid.csv` / `S###_gas.csv` line exports, computes the legacy metrics plus the new radial nonequilibrium gap metrics per case (E67, E68, E78), stores explicit `dt90_T09_s` and `dt90_T03_s` lag metrics, and outputs `analysis_results.csv`.
3.  ✅ **Sensitivity Mapping:** Implemented in `plot_comparison.py` — generates heatmaps, ΔT bar charts, direct `Δt90%` lag plots for both `T09` and `T03`, transient overlays, radial-gap comparison plots, and parameter-trend plots.
4.  ✅ **Separate scores:** `analyze_batch.py` now reports separate lag-aware scores:
    *   `Score_T09Lag` / `TotalScore_T09Lag` / `Rank_T09Lag`
    *   `Score_T03Lag` / `TotalScore_T03Lag` / `Rank_T03Lag`
    *   The legacy `Score` / `TotalScore` / `Rank` fields are preserved as aliases to the `T09`-lag score path for compatibility with existing scripts.
5.  ✅ **Script Update for S396+:** `analyze_batch.py` updated with corrected CSV mapping (T08=y11mm direct). True I_vol (T9−T8) now computed. `plot_comparison.py` updated with Score-vs-source_y plot.

### Usage
```bash
uv run --with pandas --with numpy analyze_batch.py
uv run --with pandas --with numpy --with matplotlib plot_comparison.py
```
