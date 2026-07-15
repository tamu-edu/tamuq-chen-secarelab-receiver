# COMSOL Model Analysis: Cav_Hex Validation (v7.18)

## 1. Experimental Parameter Cases (`Param_Exp`)
The model correlates three primary flux levels from the experimental campaign. The `I_factor` scales the irradiance (`ray_p`) to match measured outlet temperatures.

| Case ID | Run | Flux $q_{ap}$ [kW/m²] | `I_factor` Configuration |
| :--- | :--- | :--- | :--- |
| **E67-E71** | High | 456 | `I_f_high` (current: 1.15–1.25) |
| **E72-E76** | Med | 304 | `1.0 * (304/456)` |
| **E77-E81** | Low | 256 | `I_f_low` (current: 0.85) |

## 2. Current Calibration Parameters (S396+ Batch)

| Parameter | Symbol | Current Values | Role |
| :--- | :--- | :--- | :--- |
| Source penetration depth | `source_y` | {0.001, 0.002, 0.004} m | Shifts ray focus deeper into monolith |
| Thermal contact thickness | `tc_d` | 0.1 mm | Controls insulation–monolith coupling |
| Air conductivity multiplier | `k_air_mult` | 15.0 | Surface area surrogate for CHT |
| Longitudinal Rosseland factor | `kz_mult` | 20.0 | Long-path radiative spreading |
| Flow scaling factor | `qlpm_f_all` | {0.15, 0.2} | Leakage surrogate |
| High irradiance factor | `I_f_high` | {1.15, 1.25} | Power calibration |
| Low irradiance factor | `I_f_low` | 0.85 | Power calibration |
| Insulation specific heat | `ins_cp` | 2500 J/(kg·K) | Thermal inertia of felt |
| View factor scaling | `VF_l_f` | 4.0 | Cavity reradiation penetration |

## 3. CSV Export Mapping (Corrected from S396)

| Export Node | CSV Filename | COMSOL Probe | y Position | Physical TC |
| :--- | :--- | :--- | :--- | :--- |
| `data2` | `_T02.csv` | cpt10 | y=58mm, z=insulation | T02 (insulation) |
| `data1` | `_T03.csv` | cpt1 | y=137.5mm | T03 (gas exit) |
| `data3` | `_T08.csv` | cpt2 | y=11mm | **T08 (aperture face)** ← NEW |
| `data4` | `_T09.csv` | cpt3 | y=58mm | T09 (monolith middle) |
| `data5` | `_T10.csv` | cpt4 | y=107mm | T10 (monolith back) |

> **Note:** For runs S368–S393, `data4` was exported as `_T08.csv` (y=58mm = Physical T9), causing a label shift. From S396 onward, CSV labels match physical thermocouple positions.

## 4. Key Architectural Interconnections
- `ray_p = P_ap * I_factor / 247646`: Links total power to ray-tracing density.
- `m_tot = rho_air(T_amb) * qlpm`: Converts volumetric flow (L/min) to mass flow (kg/s) based on standard conditions.
- Raytracing study (`std2`) only re-runs when `source_y` changes (outermost loop optimization).
- Flange conduction sink (`hf4`: `-M_k*(T-T_amb)/40mm`) is implemented but currently commented out.

## 5. Summary of Current Status
The model has converged on a narrow parameter window after systematic mapping (S368–S393). The primary remaining challenge is achieving the volumetric inversion ($T9 > T8$), which is being addressed through source-shifting (`source_y > 0`). The corrected T8 export enables direct quantitative tracking of this metric for the first time.
