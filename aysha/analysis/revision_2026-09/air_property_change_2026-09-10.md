# Air properties: hand-entered table -> CoolProp

Date 2026-09-10. `receiver_reduction.py`, `outputs/` and `figures/` all regenerated.

## What changed in the code

- The four hand-entered arrays (`_T`, `_CP`, `_MU`, `_K`; 11 points on a 100 K grid, 300-1300 K) are gone. `cp_air`, `mu_air`, `k_air` now interpolate a 1 K grid over 200-1600 K generated from CoolProp `Air` at 1 atm at import (Lemmon et al. 2000 EOS; Lemmon & Jacobsen 2004 transport, the REFPROP formulations). Grid interpolation error is below 1e-6 relative. A query outside 200-1600 K now raises a RuntimeWarning instead of being silently clamped, as the old table's end values were.

- Correction note C11 added to the module docstring. `results.json` carries a new `air_properties` block with the source, CoolProp version (8.0.0), pressure and grid, so an archived run records which property release produced it. CoolProp added to the requirements and software-stack sections of `README_reproduction.md`.

- `fixed_profile_test` gained two deterministic starts per fit (uniform conductance, and the next-coarser fit of the same nested family interpolated onto the finer node set) alongside the 40 random ones, plus a monotonicity check on rms against node count. This was necessary, not cosmetic: with random starts only, the seven-node `shared_h` fit converged to rms 54.2 K under the old properties and 99.5 K under the new ones - same model family, different basin, decided by a 1% shift in the starting scale. Under the old table the warm-started code reproduces every archived rms to five digits, so the fix changes no archived number; it only removes the coin flip.

## Property deltas (CoolProp relative to the old table)

| T [K] | c_p | mu | k |
|---|---|---|---|
| 300 | -0.06% | +0.20% | +0.32% |
| 400 | +0.01% | +0.24% | -1.03% |
| 500 | -0.01% | +0.33% | -1.86% |
| 600 | +0.02% | +0.55% | -1.89% |
| 700 | -0.00% | +1.11% | -1.23% |
| 800 | -0.03% | +1.27% | -0.09% |
| 900 | -0.01% | +1.49% | +0.88% |
| 1000 | -0.00% | +1.83% | +1.46% |
| 1100 | -0.02% | +2.11% | +1.65% |
| 1200 | -0.04% | +2.59% | +1.67% |
| 1300 | -0.06% | +3.06% | +0.47% |

The old table's viscosity ran 1.1 to 3.0% low above 600 K and its conductivity up to 1.9% high between 400 and 700 K; c_p agreed to 0.06%.

## Reported quantities: archived run vs CoolProp run

754 numeric leaves in `results.json`; 495 moved by a non-negligible amount, 64 by more than 1%, 6 by more than 3%, none by more than 6.3% - and that one is a 9 K residual moving 0.6 K. Nothing structural moved.

| quantity | archived | CoolProp | change |
|---|---|---|---|
| NTU exponent (profile-corrected, const_const) | 0.340734 | 0.340812 | +0.02% |
|   its stderr | 0.0404565 | 0.0403707 | -0.21% |
| NTU exponent (uncorrected) | 0.389011 | 0.389104 | +0.02% |
|   its stderr | 0.0452784 | 0.0451763 | -0.23% |
|   r2 | 0.850256 | 0.85089 | +0.07% |
| Nu law exponent | 1.44299 | 1.4445 | +0.10% |
|   stderr | 0.0690127 | 0.0723756 | +4.87% |
| Nu law prefactor | 0.000310126 | 0.000313943 | +1.23% |
|   r2 | 0.971123 | 0.968396 | -0.28% |
| envelope Re lo | 23.3446 | 23.2359 | -0.47% |
| envelope Re hi | 93.9898 | 93.5604 | -0.46% |
| envelope Pr lo | 0.683414 | 0.697954 | +2.13% |
| envelope Pr hi | 0.688945 | 0.699778 | +1.57% |
| envelope Nu lo | 0.027828 | 0.0281645 | +1.21% |
| envelope Nu hi | 0.212616 | 0.216224 | +1.70% |
| envelope NTU lo | 0.851078 | 0.851078 | +0.00% |
| envelope NTU hi | 1.51825 | 1.51825 | +0.00% |
| envelope eps lo | 0.573045 | 0.573045 | +0.00% |
| envelope eps hi | 0.780905 | 0.780905 | +0.00% |
| envelope Bi lo | 4.8079e-06 | 4.80434e-06 | -0.07% |
| envelope Bi hi | 3.38623e-05 | 3.38385e-05 | -0.07% |
| envelope N_rc lo | 1.5034 | 1.53179 | +1.89% |
| envelope N_rc hi | 5.66 | 5.57216 | -1.55% |
| envelope Pe_LD lo | 1459.7 | 1481.21 | +1.47% |
| envelope Pe_LD hi | 5875.71 | 5964.41 | +1.51% |
| envelope Gz_L lo | 0.174987 | 0.177565 | +1.47% |
| envelope Gz_L hi | 0.704371 | 0.715004 | +1.51% |
| envelope eta_nom lo | 0.252204 | 0.252033 | -0.07% |
| envelope eta_nom hi | 1.2236 | 1.22299 | -0.05% |
| Nu/Nu_fd lo | 14.5379 | 14.2953 | -1.67% |
| Nu/Nu_fd hi | 111.075 | 109.748 | -1.19% |
| x* exit lo | 1.41971 | 1.39859 | -1.49% |
| x* exit hi | 5.7147 | 5.63174 | -1.45% |
| C_eff cooling | 298.797 | 298.301 | -0.17% |
| K_loss cooling | 0.0953495 | 0.0951325 | -0.23% |
| C_eff all6 | 121.635 | 121.561 | -0.06% |
| C_eff deep | 280.735 | 280.576 | -0.06% |
| C ratio | 2.30801 | 2.30812 | +0.00% |
| dp pred lo | 0.200674 | 0.201529 | +0.43% |
| dp pred hi | 0.991889 | 0.997066 | +0.52% |
| dp ratio lo | -0.109384 | -0.108919 | +0.42% |
| dp ratio hi | 2.7658 | 2.75144 | -0.52% |

The two headline results are unmoved: the transfer-unit exponent is +0.3408 +- 0.0404 against +0.3407 +- 0.0405 archived (profile-corrected, const_const), +0.3891 +- 0.0452 against +0.3890 +- 0.0453 uncorrected. Effectiveness, NTU and all three inversion crossings are identical to every printed digit, none of them depending on a gas property. The reference-probe prefactor span widens from 28.9 to 29.7 and the fixed-profile rms range stays 36.0 to 67.7 K.

## Manuscript strings that need editing

| value | archived | CoolProp | where (line numbers in manuscript_revised_v3.md) |
|---|---|---|---|
| Pe.L/D lower bound | 1460 | 1481 | 189, 272; supp. S5 table |
| N_rc range | 1.50 to 5.66 | 1.53 to 5.57 | 272 |
| Nu_app/Nu_fd range | 14.5 to 111 | 14.3 to 110 | 173, 302 |
| x* = 1/Gz_L range | 1.42 to 5.71 | 1.40 to 5.63 | 189 |
| Gz_L maximum | 0.704 | 0.715 | 189 |
| Nu_app absolute range | 0.028 to 0.213 | 0.028 to 0.216 | 173 |
| Nu correlation | (3.10 +- 0.12)e-4 Re^1.443, r2 0.971 | (3.14 +- 0.13)e-4 Re^1.445, r2 0.968 | 169 display, 246, 302 |
| per-flux exponents | 1.440, 1.473, 1.522 | 1.440, 1.476, 1.530 | 171 |
| Pr range | 0.683 to 0.689 | 0.698 to 0.700 | Table 1, 138-146; supp. S5 |
| Re range | 23.3 to 94.0 | 23.2 to 93.6 | Table 1; supp. S5 |
| laminar dp prediction at 15.3 sL/min | 0.99 mbar | 1.00 mbar | 235 |
| C_eff cooling | 299 +- 41 J/K | 298 +- 41 J/K | 203 |
| reference-probe prefactor span | 28.9x | 29.7x | section 5.1 |
| Tables 1, 2 and S3-S5 | - | regenerated | embedded copies to be replaced from outputs/*.md |

Effectiveness (0.573 to 0.781), the epsilon* thresholds (0.57 to 0.79), NTU (0.851 to 1.518), Bi (<= 3.4e-5), the C_eff sensor-selection ratio (2.308) and all three inversion crossings need no edit.
