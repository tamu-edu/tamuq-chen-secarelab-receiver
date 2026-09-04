# Revision 2026-09 — experimental manuscript

`manuscript_revised_v2.md` replaces `../manuscript/manuscript_full_draft.md`.
Fourteen sections consolidated into six; narrative prose throughout; two new
methodological claims in section 5. `manuscript_review.md` is the audit of the
previous draft that motivated the revision.

## Reproducing every number, table and figure

    python receiver_reduction.py --raw ../RAW --out outputs --nmc 4000
    python make_figures.py                      # writes figures/fig2..fig6

Single pass from the raw logger files; no hand-entered constants. Supersedes the
exp_analysis.py / dimensionless_analysis.py / eigenvalues_and_power.py /
uncertainty_analysis.py / make_tables.py chain in ../exp_analysis.

## Corrections relative to that chain

C1  mass flow now uses the controller reference density (1.1996 kg/m3), not
    rho(T_amb); the old form was 0.7-1.9% wrong in mdot, Re, Nu and Q_gas.
C2  gas enthalpy by trapezoidal integration of cp(T) rather than cp(Tmean)*dT.
C3  heating eigenvalues recomputed with an explicit sensor rule plus a window
    sweep. The committed eigenvalue_verification.csv could NOT be reproduced by
    eigenvalues_and_power.py: its docstring documents a three-deep-probe rule
    but the code applied all six sensors, giving C_eff = 125 J/K instead of 287.
    The rule as documented is applied here and gives 286 J/K.
C4  C_eff / K_loss identified from cooling, deep-probe heating, six-probe
    heating and a joint 18-eigenvalue fit; all four reported.
C5  inversion crossings located by local bracketing, with the global regression
    carried as an estimator systematic.
C6  LTNE indices fitted per irradiance group as well as pooled.
C7  T3 systematic-bias band (+-25 K) propagated through every T3-derived result.
C8  reference-probe sensitivity: eps, NTU and Nu recomputed against six choices
    of solid reference temperature.
C9  delivered irradiance reported as an interval bracketed by the radiometric
    flux map and the steady energy closure, with the reconciling T3 bias and the
    Gardon-gauge area-averaging factor computed explicitly.
C10 pressure drop extracted from the raw logs and compared against the laminar
    prediction and the transducer resolution.

## Outputs

outputs/groups.csv                per-run dimensionless groups, 15 runs
outputs/groups_replicates.csv     E82, E83 (archived, excluded from all fits)
outputs/eigenvalues.csv           18 eigenvalues under both heating sensor rules
outputs/reference_sensitivity.csv section 5.1 evidence table
outputs/pressure_drop.csv         section 4.7
outputs/flux_geometry.json        section 4.6 gauge-averaging calculation
outputs/uncertainty.csv           Monte Carlo, 4000 realizations
outputs/results.json              every reported scalar
outputs/table1_envelope.md        Table 1
outputs/table2_constants.md       Table 2
outputs/lit_extract.json          structured extraction of the 23 literature PDFs

## Open items for the authors

Section 4.6 flux-mapping subsection is written from refs [22,23] and should be
extended with the flux-map uncertainty propagated from the instrument chain.
Reference details flagged in the completion note at the end of the manuscript.
