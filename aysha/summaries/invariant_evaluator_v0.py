"""
Stage-0 manuscript-invariant evaluator, prototype (v0).

Applies the manuscript's own reduction (dimensionless_analysis.py) identically to
(a) the measured steady states and (b) a model's predicted steady states, so that
I1 (apparent Nu), I2 (effectiveness / axial inversion) and I3 (deep LTNE) are
compared on a like-for-like basis instead of through temperature RMSE.

Model steady states are reconstructed as  T_model = T_measured + steady_error_K
from summaries/<ver>/staged_sensor_metrics_<ver>.csv, so this runs without
re-executing the solver.

Usage (from the repo root that contains analysis/ and summaries/):
    python summaries/invariant_evaluator_v0.py 2D_v16
"""

import csv
import math
import sys

import numpy as np
from scipy.stats import linregress

VER = sys.argv[1] if len(sys.argv) > 1 else "2D_v16"

MEAS_CSV = "analysis/exp_analysis/delivered_power_check.csv"
MODEL_CSV = f"summaries/{VER}/staged_sensor_metrics_{VER}.csv"

# exact 11/58/107 mm wall quadrature (constant end extrapolation, piecewise-linear wall)
WTS = np.array([0.251825, 0.350365, 0.397810])

TTAB = [300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200]
K_AIR = [0.0263, 0.0338, 0.0407, 0.0469, 0.0524, 0.0573, 0.0620, 0.0667, 0.0715, 0.0763]
CP_AIR = [1005, 1014, 1030, 1051, 1075, 1099, 1121, 1141, 1159, 1175]
k_air = lambda T: np.interp(T, TTAB, K_AIR)
cp_air = lambda T: np.interp(T, TTAB, CP_AIR)

DH, L, N_CH = 1.5e-3, 0.137, 100
P_CH = 4 * DH
SENSORS = ["T2", "T3", "T8", "T9", "T10", "T11", "T12"]


def invariants(T8, T12, T11, T9, T10, T3, Tamb, mdot_gs):
    """Manuscript reduction applied to one steady state."""
    Tw = float(WTS @ np.array([T8, T12, T11]))
    eps = (T3 - Tamb) / (Tw - Tamb)
    Tgm = 0.5 * (Tamb + T3)
    if eps >= 1.0:
        NTU = h = Nu = float("nan")
    else:
        NTU = -math.log(1 - eps)
        mch = mdot_gs * 1e-3 / N_CH
        h = NTU * mch * cp_air(Tgm) / (P_CH * L)
        Nu = h * DH / k_air(0.5 * (Tw + Tgm))
    return dict(
        Tw=Tw, eps=eps, NTU=NTU, Nu=Nu,
        Ivol=T12 - T8,                          # axial inversion metric
        L58=(T12 - T9) / (T12 - Tamb),          # deep LTNE at 58 mm
        L107=(T11 - T10) / (T11 - Tamb),        # deep LTNE at 107 mm
    )


def load():
    meas = {}
    for x in csv.DictReader(open(MEAS_CSV)):
        meas[x["ID"]] = {k: float(x[k]) for k in x if k != "ID"}
    err = {}
    for x in csv.DictReader(open(MODEL_CSV)):
        if not x["phase"].startswith("heating"):
            continue
        err.setdefault(x["simulation_id"], {})[x["sensor"]] = float(x["steady_error_K"])
    return meas, err


def summarize(name, rows):
    Re = np.array([r[0] for r in rows])
    Nu = np.array([r[1]["Nu"] for r in rows])
    ok = np.isfinite(Nu) & (Nu > 0)
    print(f"\n[{name}]")
    if ok.sum() > 2:
        r = linregress(np.log(Re[ok]), np.log(Nu[ok]))
        print(f"  I1  Nu_app   = {math.exp(r.intercept):.3e} Re^{r.slope:.3f}   r2={r.rvalue**2:.3f}  (n={ok.sum()}/{len(Nu)})")
    else:
        print(f"  I1  Nu_app   undefined: eps>=1 in {len(Nu)-ok.sum()}/{len(Nu)} cases")
    eps = np.array([r[1]["eps"] for r in rows])
    Iv = np.array([r[1]["Ivol"] for r in rows])
    print(f"  I2  eps range = {eps.min():.3f} .. {eps.max():.3f}   mid-peak (T12>T8) = {(Iv > 0).sum()}/{len(Iv)}")
    for tag, key in (("L58", "L58"), ("L107", "L107")):
        Lm = np.array([r[1][key] for r in rows])
        rr = linregress(Re, Lm)
        print(f"  I3  {tag:5s}= {rr.intercept:+.4f} {rr.slope:+.2e}*Re   r2={rr.rvalue**2:.3f}   sign>0 = {(Lm > 0).sum()}/{len(Lm)}")


def main():
    meas, err = load()
    hdr = f"{'ID':5}{'Re':>7} | {'eps_m':>6}{'eps_M':>7} | {'Nu_m':>7}{'Nu_M':>8} | {'T12-T8 m':>9}{'M':>8} | {'L107 m':>7}{'M':>8}"
    print(hdr)
    print("-" * len(hdr))
    rm, rM = [], []
    for ID, m in meas.items():
        if ID not in err:
            continue
        M = {s: m[s + "_ss"] + err[ID][s] for s in SENSORS}
        im = invariants(m["T8_ss"], m["T12_ss"], m["T11_ss"], m["T9_ss"], m["T10_ss"],
                        m["T3_ss"], m["Tamb"], m["mdot_gs"])
        iM = invariants(M["T8"], M["T12"], M["T11"], M["T9"], M["T10"],
                        M["T3"], m["Tamb"], m["mdot_gs"])
        rm.append((m["Re"], im))
        rM.append((m["Re"], iM))
        print(f"{ID:5}{m['Re']:7.1f} | {im['eps']:6.3f}{iM['eps']:7.3f} | "
              f"{im['Nu']:7.3f}{iM['Nu']:8.3f} | {im['Ivol']:9.1f}{iM['Ivol']:8.1f} | "
              f"{im['L107']:7.3f}{iM['L107']:8.3f}")
    summarize("MEASURED", rm)
    summarize(f"MODEL {VER}", rM)


if __name__ == "__main__":
    main()
