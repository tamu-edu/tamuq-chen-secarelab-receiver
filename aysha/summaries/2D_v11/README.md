# 2D v11 simulation artifacts

Generated on 2026-07-29 with Julia 1.12.6 and the workspace project
environment.

## Reproduce

```powershell
& "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. run_2D_v11.jl
& "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. plot_2D_v11.jl
& "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. test\smoke_2D_v11.jl
& "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. test\mesh_2D_v11.jl
& "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. run_2D_v11_groove_sensitivity.jl
& "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. run_2D_v11_inheritance_sensitivity.jl
& "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. run_2D_v11_betaopt_confirmation.jl
& "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. --threads=6 calibrate_2D_v11_staged.jl
& "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. plot_2D_v11_staged.jl
& "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. --threads=6 test\check_2D_v11_staged.jl
```

## Tables

- `cold_dp1_cases_2D_v11.csv`: selected t0 observations and fitted cold
  predictions.
- `cold_groove_profile_2D_v11.csv`: fixed-`K` groove/resistance
  identifiability profile.
- `model_form_cases_2D_v11.csv`: final sensor, profile, pressure and flow
  metrics for all six variants and 15 experiments.
- `model_form_slopes_2D_v11.csv`: modeled and observed axial-offset slopes
  for each irradiance group.
- `model_form_summary_2D_v11.csv`: aggregate acceptance metrics.
- `axial_profiles_2D_v11.csv`: representative final solid/gas axial fields.
- `ring_profiles_2D_v11.csv`: representative final ring-flow and pressure
  fields.
- `representative_transients_2D_v11.csv`: temporal comparison for E67, E72
  and E77.
- `mesh_cases_2D_v11.csv` and `mesh_comparison_2D_v11.csv`: sensitivity
  versus nominal-mesh confirmation.
- `groove_geometry_fits_2D_v11.csv`,
  `groove_geometry_cases_2D_v11.csv`,
  `groove_geometry_slopes_2D_v11.csv` and
  `groove_geometry_summary_2D_v11.csv`: complete 12/13/14 mm free-diameter
  sensitivity.
- `inheritance_sensitivity_cases_2D_v11.csv` and
  `inheritance_sensitivity_summary_2D_v11.csv`: broad one-at-a-time screen
  of active v9-derived solid/optical parameters and the Graetz scale.
- `betaopt110_confirmation_cases_2D_v11.csv`,
  `betaopt110_confirmation_slopes_2D_v11.csv` and
  `betaopt110_confirmation_summary_2D_v11.csv`: full 15-case confirmation
  of the optical-extinction sensitivity.
- `staged_parameters_2D_v11.csv`, `staged_profile_trace_2D_v11.csv` and
  `staged_power_trace_2D_v11.csv`: selected parameters and all true-model
  fitting evaluations.
- `staged_identifiability_columns_2D_v11.csv` and
  `staged_identifiability_correlation_2D_v11.csv`: local rank diagnostics.
- `staged_sensor_metrics_2D_v11.csv`,
  `staged_case_metrics_2D_v11.csv`,
  `staged_validation_summary_2D_v11.csv` and
  `staged_flow_slopes_2D_v11.csv`: training, held-out heating and cooling
  validation.
- `staged_mesh_cases_2D_v11.csv` and
  `staged_mesh_comparison_2D_v11.csv`: fitted sensitivity/nominal-mesh
  comparison.
- `staged_acceptance_status_2D_v11.txt`: formal coefficient-extraction
  decision.

## Figures

`plots/` contains cold DP1, parity, axial-flow trend, radial-offset,
representative axial, representative temporal and ring-flow figures.

## Decision

The numerical implementation passes conservation, pressure-equality,
asymptotic-Nu and mesh checks. The frozen-v9 parameter point is rejected for
v11 because the literature Graetz variants reverse the measured axial flow
trends and the groove-induced flow maldistribution does not produce the
measured radial solid-temperature offsets. At that stage this was not yet an
independently recalibrated v11-versus-v9 comparison. See `../journal.2D.md` and
`../2D_v11_theory_manual.md`.

The subsequent inheritance screen shows why that qualification matters:
setting optical extinction to 110 rather than the v9-fitted 21.29 1/m
restores all three axial slope signs and reduces axial RMSE to 58.17 K,
although the slopes remain too weak and the radial offsets remain absent.
This supported the controlled v11-specific refit; it did not accept 110 1/m
as a fitted result.

That staged refit is now complete. It restores the three axial flow-slope
signs and reduces all-heating axial RMSE to 44.14 K, but fails held-out
absolute temperatures, transient heating/cooling times, radial ordering and
heating DP1. The fitted vector is rejected for coefficient extraction.
Plots are in `plots/staged/`.
