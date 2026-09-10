# Proposed manuscript edits after the CoolProp property swap

Line numbers are `manuscript_revised_v3_withfigures.md` (the working draft) unless marked *supp*. Values come from the regenerated `outputs/results.json`. Nothing here is applied - these are replacement strings for you to accept or reject.

## Section 3.1, methods - the property source is currently unstated

L76 ends the property sentence at the film-temperature convention. Suggested addition after "...held fixed across every variant in section 5.1.":

> Air properties are evaluated from the Lemmon equation of state and the Lemmon and Jacobsen transport correlations for dry air at 1 atm, as implemented in CoolProp [ref].

That needs one new reference (Bell, Wronski, Quoilin and Lemort, *Industrial & Engineering Chemistry Research* 53 (2014) 2498-2508, doi:10.1021/ie4033999). It would be cited first in section 3.1, which sits ahead of the current [17], so a numbered list ordered by first appearance would have to renumber - your call whether to insert it there or append it as [47].

## Numbers that move

| line | current | replacement |
|---|---|---|
| 191 | `(3.10 \pm 0.12)\times10^{-4}\, Re_{\rm nom}^{1.443}, \qquad r^2 = 0.971` | `(3.14 \pm 0.12)\times10^{-4}\, Re_{\rm nom}^{1.444}, \qquad r^2 = 0.968` |
| 193 | `the exponent is 1.440, 1.473 and 1.522` | `the exponent is 1.440, 1.476 and 1.530` |
| 193 | `pooled regression standard error ±0.069` | `pooled regression standard error ±0.072` |
| 195 | `Absolute values run from 0.028 to 0.213.` | `Absolute values run from 0.028 to 0.216.` |
| 195 | `a factor of 14.5 to 111 below Nu_H2, moving only to 14.0–107 or 17.0–130` | `a factor of 14.3 to 110 below Nu_H2, moving only to 13.8–106 or 16.7–128` |
| 197 | `our 1.443 over Re = 23–94` | `our 1.444 over Re = 23–94` |
| 211 | `The axial Péclet group is at least 1460` | `The axial Péclet group is at least 1481` |
| 211 | `the radiation–conduction number is 1.50 to 5.66` | `the radiation–conduction number is 1.53 to 5.57` |
| 211 | `The Graetz number reaches only 0.704; its reciprocal, x* = 1/Gz_L, runs from 1.42 at highest flow to 5.71 at lowest` | `The Graetz number reaches only 0.715; its reciprocal, x* = 1/Gz_L, runs from 1.40 at highest flow to 5.63 at lowest` |
| 235 | `gives 299 ± 41 J K⁻¹` | `gives 298 ± 41 J K⁻¹` |
| 289 | `(2.74 against 0.99 mbar at 15.3 sL min⁻¹)` | `(2.74 against 1.00 mbar at 15.3 sL min⁻¹)` |
| 305, 319, 369 | `a factor of 28.9` | `a factor of 29.7` |
| 333 | `Re·Pr·L/D_h ≥ 1460` | `Re·Pr·L/D_h ≥ 1481` |
| 333 | `N_rc of 1.50 to 5.66` | `N_rc of 1.53 to 5.57` |
| 363 | `1.470 ± 0.011 with r² = 0.9994` | `1.472 ± 0.011 with r² = 0.9993` |
| 363 | `prefactors 2.55 to 3.20×10⁻⁴` | `prefactors 2.55 to 3.24×10⁻⁴` |
| 363 | `the pooled form gives 1.443 ± 0.069 at r² = 0.971` | `the pooled form gives 1.444 ± 0.072 at r² = 0.968` |
| 363 | `this is 14.5 to 111 below` | `this is 14.3 to 110 below` |
| *supp* 167 | `(2.74 against 0.99 mbar at 15.3 sL min⁻¹)` | `(2.74 against 1.00 mbar at 15.3 sL min⁻¹)` |

## Embedded tables to re-paste from `outputs/`

| manuscript table | source file | what moved |
|---|---|---|
| Table 1, L138-146 | `outputs/table1_envelope.md` | Re, Pr, Gz_L, Nu columns (Pr rises about 2%, Nu about 1.5%); eps, NTU, eta unchanged |
| Table 2, L246 | `outputs/table2_constants.md` | the $Nu_{\rm app}$ prefactor row (3.104 to 3.143 e-4, sd 0.111 to 0.112 e-4) and the exponent row's SE and r-squared; every other row moves in the fourth digit or not at all |
| Table 4 (reference sensitivity), L312 | `outputs/reference_sensitivity.csv` | prefactor and Nu_min/Nu_max columns, up to 2.2% |
| Table S5, *supp* L110-125 | `outputs/tableS5_auxiliary_groups.md` | Pr, Gz_L, x*, Pe, Bi, N_rc columns |
| Tables S3, S4 | `outputs/tableS3_heating_conditionality.md`, `outputs/tableS4_wall_and_falsification.md` | C_eff and K_loss in the fourth digit; the fixed-conductance rms values move by less than 0.03 K |

## Checked and needing no edit

The abstract's transfer-unit exponent (Re_nom^+0.341 +- 0.041), the required -1.00, the 36 to 68 K rms range of the fixed-conductance families, the inversion threshold 0.57 to 0.79, the capacitance span 122 to 281 J K⁻¹, effectiveness 0.573 to 0.781, NTU 0.851 to 1.518, Bi <= 3.4e-5, the 2.31 sensor-selection capacitance ratio, all three crossing fluxes, the 1.4 to 2.8 pressure-drop ratio and the -0.02 against 0.20 mbar low-flow figure are all unchanged to the digits printed.

One judgement call: the abstract says the prefactor spans "a factor of 29", and the span is now 29.7 rather than 28.9. Still 29 if truncated, 30 if rounded - if you want it to read as a round number, "a factor of thirty" is now the honest one.
