# exp_analysis — experimental reduction, figures and tables

Everything the manuscript reports from the experimental campaign is produced by
the scripts in this directory, starting from the raw logger files in `../RAW/`.

**Project rule: every figure and every table ships with the script that
extracted, processed and plotted its data.** A number in the manuscript that no
script here reproduces is a defect, not an oversight. The "Known gaps" section
at the end lists the places where that rule is not yet met.

---

## 1. Quick start

```bash
cd analysis/exp_analysis
python3 exp_analysis.py            # raw -> steady metrics, cooling fits, figs 1-5
python3 dimensionless_analysis.py  # -> dimensionless groups, figs 6-7
python3 eigenvalues_and_power.py   # -> delivered-power closure, eigenvalues
python3 uncertainty_analysis.py    # -> Monte Carlo uncertainties and constants
python3 make_tables.py             # -> manuscript Tables 1 and 2
python3 fig8_model_comparison.py   # -> Figure 7 (needs ../../summaries/1D_v36)
```

Run in that order: each step reads CSVs written by an earlier one. All scripts
resolve paths from their own location (`__file__`), so the tree can be moved or
cloned without editing anything. Outputs are written beside the scripts.

**Environment.** Python 3 with `numpy`, `pandas`, `matplotlib`, `scipy`.
Verified on numpy 2.2, pandas 2.3, matplotlib 3.10, scipy 1.15.

**Runtime.** `exp_analysis.py` and `dimensionless_analysis.py` each read all 22
raw files and take roughly 30–60 s. The others are seconds.

---

## 2. The scripts

### `exp_analysis.py` — raw data to steady-state metrics
The entry point and the only script that touches `../RAW/`.

*Reads* `../RAW/Data_FPT*.csv` (22 files: 15 heating runs E67–E81, 2 unassigned
runs E82–E83, 3 cooling runs C69/C80/C81, plus duplicates of the parent runs).
The column contract follows `import_exp_1D_v2.jl`: time = 1, MFC actual flows =
6:9, T2 = 35, T3 = 36, T8 = 41, T9 = 42, T10 = 43, T11 = 44, T12 = 45. All four
controllers carry air and the flow is their sum.

> **`WTS` is the wall quadrature and is also owned here.** It is computed from
> the probe coordinates z8/z9/z10 and L rather than hardcoded, giving
> 0.2518 / 0.3504 / 0.3978 — the length-average weights for a piecewise-linear
> wall profile through probes at 11/58/107 mm with constant extrapolation to
> the end faces. Every other script imports it. It was previously hardcoded in
> five places as 0.248/0.365/0.387, which are the weights for probes at
> 10/58/110 mm and did not match the stated positions; correcting it shifted
> ε, Nu, ε*, the irradiances and everything downstream by a few tenths of a
> per cent.

*Defines*, and therefore owns for every downstream script: the run dictionaries
`heating` and `cooling`, the frontal area `A_frt`, the length `L`, the sensor
coordinates `z8/z9/z10`, the air property functions `cp_air` and `rho_air`, the
steady-window helper `ss`, the t90 helper, and the delivered-power factors
`F_DEL`.

*Writes* `steady_state_metrics.csv` (steady means of every sensor, the inversion
indicator, interior-to-wall gaps, t90 values, mass flow, gas enthalpy, nominal
and delivered power, efficiencies), `flow_slopes.csv` (per-irradiance
flow-sensitivity regressions), `cooling_time_constants.csv` (per-sensor
exponential fits, full and late window), and `fig1_transients.png`,
`fig2_ss_vs_flow.png`, `fig4_cooling_lin.png`, `fig5_t90.png`.

> **`F_DEL` is the delivered-power basis and is the one number in this directory
> a reader must understand.** It is the group means of the model-free energy
> closure evaluated with the tangent loss conductance from the heating
> eigenvalue (0.119 W/K): `{456: 1.146, 304: 1.343, 256: 0.931}`. The secant
> conductance from cooling (0.097 W/K) gives `{456: 1.052, 304: 1.230,
> 256: 0.845}`, carried as a one-sided systematic band of about −9%. The closure
> is a steady balance and needs the secant conductance at the hot steady state,
> which those two values bracket. No thermal model enters.
>
> Note the two-pass dependency: `F_DEL` lives in `exp_analysis.py`, which runs
> first, but its value comes from `eigenvalues_and_power.py`, which runs third.
> That script prints the factors on each run; keep the constant consistent with
> what it prints. Changing `F_DEL` changes `eta_del` and `PR_del` everywhere and
> requires re-running the whole chain.
>
> `G_DEL` is derived from it: `G_del = f * G0 = {456: 523, 304: 409, 256: 239}`
> kW/m2, the delivered aperture irradiance. **All results are labelled by
> G_del**, and the helper `GLAB(Io)` emits those legend strings; every figure
> and both tables use it. The nominal G0 survives only in the simulator
> characterization, the nominal column of Table 1, the energy closure that
> identifies f, and the statement of the delivered-power discrepancy.

### `dimensionless_analysis.py` — dimensionless groups and correlations
*Imports* `exp_analysis` (which re-runs it, since it has no `__main__` guard).

*Reads* the raw files again through `load`, plus the geometry and properties
from `exp_analysis`.

*Computes* per run: Re, Pr, Graetz number at the outlet, the wall quadrature
`Tw = 0.248·T8 + 0.365·T12 + 0.387·T11`, effectiveness `eps = (T3−Tamb)/(Tw−Tamb)`,
`NTU = −ln(1−eps)`, the apparent heat-transfer coefficient and Nusselt number,
the Hausen reference Nusselt number and their ratio, the Biot number, the
nonequilibrium indices at 58 and 107 mm, the nominal and delivered specific
energy, and the wall inversion indicator.

*Prints* the Nu = a·Re^b fit, the inversion crossings per irradiance group, the
LTNE regressions, and the cooling slow-mode identification.

*Writes* `dimensionless_groups.csv`, `fig6_master_curve.png`,
`fig7_Nu_LTNE_eps.png`.

### `eigenvalues_and_power.py` — slow-mode eigenvalues and the power closure
Added August 2026 to close a reproduction gap: the two files it produces
previously existed with no generating script.

*Reads* the raw files through `exp_analysis.load`, plus
`steady_state_metrics.csv` and `dimensionless_groups.csv`.

*Computes* the model-free delivered-power factor per run,
`f = [Q_gas + K·(Tw − Tamb)] / (G0·A_frt)`, at both ends of the identified
K_loss bracket (0.096 and 0.119 W/K); and the slow-mode eigenvalues, by
regressing `log(T − Tamb)` against time on the late cooling window and
`log(T_ss − T)` against time on the late heating window, averaged over the six
receiver sensors.

*Writes* `delivered_power_check.csv` — **reproduces the committed file exactly**
— and `eigenvalue_verification_recomputed.csv`, deliberately under a separate
name; see "Known gaps".

### `uncertainty_analysis.py` — Monte Carlo uncertainty and identified constants
*Reads* the raw files through `exp_analysis`, plus `eigenvalue_verification.csv`.

*Propagates*, over 4000 realizations: thermocouple bias (σ = 1.1 K, class ±2.2 K)
plus a 0.5 K steady-window drift; per-controller flow bias (σ = 0.025 sL/min on
each of four); property tabulations (±2% on μ and k, ±1% on cp); radiometer
calibration (±3%); the delivered-power factor across its one-sided K bracket
(−9%) with a further ±3% for the closure's flow residual; and each eigenvalue
by its own fit standard error. The whole analysis chain is recomputed per
realization, so nonlinear amplification (notably `dNTU = dε/(1−ε)`) is captured
without linearization.

*Writes* `uncertainty_per_run.csv` (per-run standard deviations) and
`identified_constants_mc.csv` (the Nusselt law, the three inversion thresholds,
C_eff and K_loss from each dataset, and the monolith capacitance, each with a
95% interval). **`identified_constants_mc.csv` is the authoritative source for
manuscript Table 2** — before August 2026 these constants were printed to
stdout only, which is why the table could not be regenerated.

### `make_tables.py` — manuscript Tables 1 and 2
*Reads* `dimensionless_groups.csv`, `uncertainty_per_run.csv`,
`steady_state_metrics.csv`, `eigenvalue_verification.csv`,
`delivered_power_check.csv`, `identified_constants_mc.csv`.

Table 2 takes its values from `identified_constants_mc.csv` and adds an
`independent_recheck` column in which this script re-identifies the same
constant by a different route, so drift between the two appears in the table
rather than hiding. The two agree: Nu prefactor 3.08e-4 against 3.13e-4,
exponent 1.444 against 1.440, ε* identical to three decimals in all three
groups, C_eff 301/288 J/K against 300/287, K_loss 0.097/0.120 against
0.096/0.119.

*Writes* `table1_campaign_envelope.csv`/`.md`,
`table2_identified_constants.csv`/`.md`.

> **Two conventions exist for ε\*** and they differ by 0.01–0.02. The manuscript
> convention (§4.4) regresses the inversion indicator on flow to locate `q*`,
> then evaluates the effectiveness regression there, using all five runs of a
> group. `make_tables.py` reports that value and carries the direct
> ε-interpolation alternative in its note column.

### `fig8_model_comparison.py` — Figure 7, model comparison
The only script that reads outside this directory and `../RAW/`.

*Reads* `../../summaries/1D_v36/steady_results_fitted_energy_accounting_1D_v36.csv`
and `dimensionless_groups.csv`. It runs no simulation: the 1D model output is
an input, produced by the Julia runners in the repository root.

*Applies* the §3.2 reductions identically to model and measurement, and the §4.4
convention to both, so the two are comparable.

*Writes* `fig8_model_comparison.png`, `fig8_model_comparison_data.csv` (every
plotted value) and `fig8_inversion_thresholds.csv` (the per-group thresholds
quoted in §4.9).

---

## 3. Manuscript figures and tables

| Manuscript | File | Script |
|---|---|---|
| Figure 1 | — | apparatus schematic, drawn by hand |
| Figure 2 | `fig2_ss_vs_flow.png` | `exp_analysis.py` |
| Figure 3 | `fig3_inversion_gaps_eta.png` | `dimensionless_analysis.py` |
| Figure 4 | `fig7_Nu_LTNE_eps.png` | `dimensionless_analysis.py` |
| Figure 5 | `fig4_cooling_lin.png` | `exp_analysis.py` |
| Figure 6 | `fig6_master_curve.png` | `dimensionless_analysis.py` |
| Figure 7 | `fig8_model_comparison.png` | `fig8_model_comparison.py` |
| Table 1 | `table1_campaign_envelope.csv`/`.md` | `make_tables.py` |
| Table 2 | `table2_identified_constants.csv`/`.md` | `make_tables.py` |

The file numbering and the manuscript numbering do not match, because
`fig1_transients.png` and `fig5_t90.png` are diagnostics that no longer appear
in the paper. Use this table, not the filenames.

## 4. Derived-data products

| File | Written by | Contents |
|---|---|---|
| `steady_state_metrics.csv` | `exp_analysis.py` | steady means, t90, gaps, efficiencies |
| `flow_slopes.csv` | `exp_analysis.py` | per-group flow-sensitivity regressions |
| `cooling_time_constants.csv` | `exp_analysis.py` | per-sensor τ, full and late window |
| `dimensionless_groups.csv` | `dimensionless_analysis.py` | Re, Gz, ε, NTU, Nu, Λ, I_vol, PR |
| `delivered_power_check.csv` | `eigenvalues_and_power.py` | per-run delivered-power factors |
| `eigenvalue_verification.csv` | *(no generator — see below)* | slow-mode eigenvalues, both phases |
| `uncertainty_per_run.csv` | `uncertainty_analysis.py` | Monte Carlo σ per run |
| `identified_constants_mc.csv` | `uncertainty_analysis.py` | Table 2 constants with intervals |
| `fig8_*.csv` | `fig8_model_comparison.py` | plotted values and thresholds |

## 5. Reproduction test, 31 August 2026

The scripts were copied to an empty directory with only `../RAW/` symlinked, and
run in order. Everything below was regenerated **byte-for-byte identical** to the
committed copies: `steady_state_metrics.csv`, `flow_slopes.csv`,
`cooling_time_constants.csv`, `dimensionless_groups.csv`, `uncertainty_per_run.csv`,
`delivered_power_check.csv`, and figures 1–7 (all file names). The three cooling
eigenvalues and the abscissa of all eighteen rows of `eigenvalue_verification.csv`
also reproduce exactly; its fifteen heating eigenvalues do not (gap 1 below).
Repeat this test after any change to the reduction — the two missing generators
were invisible until it was run, because the files sat in the folder looking
finished.

## 6. Known gaps

1. **Heating eigenvalues: resolved, with a changed method.** The committed
   `eigenvalue_verification.csv` had no generator. Its abscissa convention and
   its three cooling eigenvalues were recovered exactly. Its fifteen heating
   eigenvalues could not be reproduced by the published recipe applied to all
   six sensors, nor by any of 126 combinations of fit window and r² threshold —
   every one returned C_eff between 63 and 129 J/K against the committed 302.
   The cause is sensor selection: the committed values come from the deep and
   outlet probes only. `eigenvalues_and_power.py` now implements that rule
   ({T11, T10, T3}, physically justified in manuscript Section 3.4) and the file
   is regenerated from raw data. The identified constants changed accordingly:
   heating C_eff 304 → 288 J/K, K_loss 0.164 → 0.119 W/K, r² 0.79 → 0.90, and
   the delivered-power factors and every efficiency moved with them. The whole
   chain is now reproducible; there is no remaining orphan data product.

2. **Figure 1** is an apparatus schematic with no generator.
3. **Runs E82 and E83** carry no irradiance level and are excluded from every
   manuscript figure and table. They remain in `steady_state_metrics.csv` and
   `dimensionless_groups.csv` with `Io_kWm2` blank; drop those rows when
   consuming either file.
4. **`exp_analysis.py` has no `__main__` guard**, so importing it re-runs the
   whole reduction. Harmless but slow; add a guard if the reduction grows.
5. **Model-side artefacts** under `../../summaries/1D_v*/` and `2D_v*/` are
   written by the Julia runners in the repository root, which are outside this
   directory. Figure 7 consumes only a saved CSV, so it is reproducible from
   here alone.
