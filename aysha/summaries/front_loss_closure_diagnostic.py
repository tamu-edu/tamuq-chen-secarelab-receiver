"""
Model-free steady energy-closure diagnostic for the 2D receiver study.

Question: with the delivered power fixed inside a lamp group (which it must be,
physically), what loss law makes the steady energy books close across the flow
sweep?

    Q_abs = Q_gas(measured) + Q_loss(T)
    f     = Q_abs / (Io * A_front)

f MUST be constant within a lamp group. The spread of f across the five flow
settings of a group is therefore a direct, model-free test of the assumed loss
law and of the temperature that drives it.

Result (see summaries/2D_v20_energy_closure_assessment.md):
  * a linear conductance K = 0.096-0.164 W/K (the manuscript slow-mode value,
    and the basis of the power-closure brackets) leaves f varying by 1.25-1.95x
    inside a single group -> rejected;
  * a T^4 loss driven by the FRONT temperature closes f to 0.3-3.4% CV and puts
    every group's f inside its independent closure bracket;
  * the same T^4 loss driven by the rear (T11) temperature is absurd, because
    T11 is nearly clamped across each flow sweep.

Usage:  python summaries/front_loss_closure_diagnostic.py
        (run from the repo root containing analysis/ and summaries/)
"""

import csv

import numpy as np

SIG = 5.67e-8
A_FRONT = 0.019 ** 2                       # 19 x 19 mm monolith face
W = np.array([0.251825, 0.350365, 0.397810])   # exact 11/58/107 mm quadrature
BRACKET = {456: (1.05, 1.34), 304: (1.23, 1.58), 256: (0.84, 1.11)}

TTAB = [300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200]
CP = [1005, 1014, 1030, 1051, 1075, 1099, 1121, 1141, 1159, 1175]
cp_air = lambda T: np.interp(T, TTAB, CP)

SIG_GRID = np.linspace(0.0, 4e-9, 8001)


def load(path="analysis/exp_analysis/delivered_power_check.csv"):
    groups = {}
    for x in csv.DictReader(open(path)):
        Io = float(x["Io_kWm2"])
        Ta = float(x["Tamb"])
        mdot = float(x["mdot_gs"]) * 1e-3
        T8, T12, T11, T3 = (float(x[k + "_ss"]) for k in ("T8", "T12", "T11", "T3"))
        groups.setdefault(Io, []).append(dict(
            ID=x["ID"], Re=float(x["Re"]), Ta=Ta, T8=T8, T12=T12, T11=T11,
            Tw=float(W @ np.array([T8, T12, T11])),
            Qgas=mdot * cp_air(0.5 * (Ta + T3)) * (T3 - Ta),
            Qnom=Io * 1e3 * A_FRONT,
        ))
    return groups


def best_sigma(g, Tdrv):
    """sigma_eff that minimises the coefficient of variation of f in this group."""
    Ta = np.array([c["Ta"] for c in g])
    Q = np.array([c["Qgas"] for c in g])
    Qn = g[0]["Qnom"]
    D = Tdrv ** 4 - Ta ** 4
    cv = [np.std((Q + v * D) / Qn) / np.mean((Q + v * D) / Qn) for v in SIG_GRID]
    i = int(np.argmin(cv))
    return SIG_GRID[i], (Q + SIG_GRID[i] * D) / Qn, cv[i]


def linear_law(groups):
    print("A. Linear loss conductance (the current power-closure convention)\n")
    print(f"{'group':>6} {'K W/K':>7} {'mean f':>7} {'max/min':>8} {'CV %':>6}")
    for Io in sorted(groups, reverse=True):
        g = groups[Io]
        for K in (0.096, 0.164):
            f = np.array([(c["Qgas"] + K * (c["Tw"] - c["Ta"])) / c["Qnom"] for c in g])
            print(f"{Io:6.0f} {K:7.3f} {f.mean():7.2f} {f.max()/f.min():8.2f} "
                  f"{100*f.std()/f.mean():6.1f}")
    print("\n  -> f is not constant within a group. The linear law is rejected.\n")


def driver_scan(groups):
    print("B. Which temperature drives the missing loss?  Qloss = sig*(Tdrv^4 - Ta^4)\n")
    print(f"{'group':>6} {'driver':<14} {'sig_eff':>9} {'eps*A cm2':>10} {'mean f':>7} {'CV %':>6}")
    for Io in sorted(groups, reverse=True):
        g = groups[Io]
        for name, T in (("T8  (11 mm)", np.array([c["T8"] for c in g])),
                        ("T12 (58 mm)", np.array([c["T12"] for c in g])),
                        ("T11 (107 mm)", np.array([c["T11"] for c in g])),
                        ("Tw  (quadr.)", np.array([c["Tw"] for c in g]))):
            s, f, cv = best_sigma(g, T)
            print(f"{Io:6.0f} {name:<14} {s:9.2e} {s/SIG*1e4:10.1f} {f.mean():7.2f} {100*cv:6.1f}")
        print()
    print("  -> the front temperature is the best driver; the rear is nearly clamped.\n")


def front_closure(groups, offsets=(0, 100, 150, 250)):
    print("C. Front-radiation closure, Tfront = T8 + dT\n")
    print(f"{'group':>6} {'dT K':>5} {'sig_eff':>9} {'eps*A cm2':>10} {'mean f':>7} "
          f"{'CV %':>6}  bracket")
    selected = {}
    for Io in sorted(groups, reverse=True):
        g = groups[Io]
        for d in offsets:
            T = np.array([c["T8"] for c in g]) + d
            s, f, cv = best_sigma(g, T)
            lo, hi = BRACKET[Io]
            tag = "INSIDE" if lo <= f.mean() <= hi else ("above" if f.mean() > hi else "below")
            print(f"{Io:6.0f} {d:5d} {s:9.2e} {s/SIG*1e4:10.1f} {f.mean():7.2f} "
                  f"{100*cv:6.1f}  {lo:.2f}-{hi:.2f} {tag}")
            if d == 150:
                selected[Io] = (s, f.mean())
        print()
    return selected


def front_share(groups, selected, offset=150):
    print(f"D. Front-loss share of absorbed power (Tfront = T8 + {offset} K)\n")
    print(f"{'case':<6}{'Re':>7}{'Tfront':>8}{'Qabs W':>9}{'Qfront W':>10}{'% Qabs':>8}")
    for Io in sorted(groups, reverse=True):
        sig, f = selected[Io]
        for c in groups[Io]:
            Tf = c["T8"] + offset
            Qa = f * c["Qnom"]
            Ql = sig * (Tf ** 4 - c["Ta"] ** 4)
            print(f"{c['ID']:<6}{c['Re']:7.1f}{Tf:8.0f}{Qa:9.1f}{Ql:10.1f}{100*Ql/Qa:8.0f}")
    cap = SIG * 1090 ** 4 * A_FRONT
    print(f"\n  model's radiating front area is only {A_FRONT*1e4:.2f} cm2 "
          f"-> blackbody cap at 1090 K = {cap:.1f} W")


if __name__ == "__main__":
    G = load()
    linear_law(G)
    driver_scan(G)
    sel = front_closure(G)
    front_share(G, sel)
