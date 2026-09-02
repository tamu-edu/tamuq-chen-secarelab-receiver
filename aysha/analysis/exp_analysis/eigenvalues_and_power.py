"""Slow-mode eigenvalues and the model-free delivered-power closure.

Fills the two gaps in the reproduction chain: `eigenvalue_verification.csv`
and `delivered_power_check.csv` previously existed as committed data files
with no generating script.

Method (manuscript Sections 3.4, 4.6 and 4.7)
  Cooling: after the lamps are switched off the assembly decays as a single
  slow mode. For each receiver sensor the late-window excess temperature is
  regressed as log(T - T_amb) against time; lambda is the mean of the
  per-sensor slopes and lam_sd their spread.
  Heating: the late approach to steady state decays with the same eigenvalue,
  so log(T_ss - T) is regressed against time over the same late window -- but
  only for the deep and outlet probes T11 (107 mm), T10 (107 mm interior) and
  T3 (outlet gas). With the lamps on, the front probes T8, T12 and T9 are held
  by a local balance between absorbed flux and aperture re-radiation, with its
  own short time constant; their approach is bi-exponential and does not carry
  the assembly-scale mode. With the lamps off that local source is absent, the
  field flattens, and all six sensors share the slow mode -- which is why the
  cooling identification needs no sensor rule and the heating one does.
  Using all six sensors in heating gives C_eff = 125 J/K, an artefact.
  Heating window: the deficit fraction u = (T_ss - T)/(T_ss - T_0) is fitted
  on 0.07 < u < 0.45, excluding early times so the fast modes have decayed and
  the tail where the logarithm amplifies noise; sensors with r2 < 0.95 are
  discarded.
  Abscissa: x = eps(q) * mdot * cp, the gas advective conductance, with eps
  taken from the pooled steady correlation eps(q) fitted over the fifteen
  heating runs -- not from each run's own eps. This convention reproduces the
  committed abscissa to better than 0.15% and the cooling eigenvalues exactly.
  Delivered power: f = [Q_gas + K*(Tw_bar - T_amb)] / (G0 * A_frt), evaluated
  at both ends of the identified K_loss bracket.

Inputs : ../RAW/Data_FPT*.csv (through exp_analysis.load), steady_state_metrics.csv,
         dimensionless_groups.csv
Outputs: delivered_power_check.csv (reproduces the committed file exactly),
         eigenvalue_verification_recomputed.csv (see README: heating rows do
         not yet reproduce the committed eigenvalue_verification.csv)
"""
import os, sys
import numpy as np
import pandas as pd
from scipy.stats import linregress

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from exp_analysis import load, ss, heating, cooling, A_frt, cp_air, F_DEL, WTS as _WTS
from scipy.stats import linregress as _lr

P = lambda f: os.path.join(_HERE, f)
RHO_STD = 101325 / (287.05 * 294.25)          # Aalborg standard: 21.1 C, 1 atm
WTS = dict(zip(("T8", "T12", "T11"), _WTS))   # wall quadrature from exp_analysis
SENS = ["T8", "T12", "T11", "T9", "T10", "T3"]        # cooling: all six
SENS_SLOW = ["T11", "T10", "T3"]                      # heating: deep + outlet
K_BRACKET = (0.097, 0.119)   # secant (cooling) and tangent (heating) conductance
PARENT = {"C69": "E69", "C80": "E80", "C81": "E81"}


def eigen_cooling(d, t, Tamb):
    lam = []
    m = t > t[-1] / 2
    for s in SENS:
        th = d[s][m] - Tamb
        good = th > 5
        if good.sum() > 50:
            lam.append(-linregress(t[m][good], np.log(th[good])).slope)
    return np.array(lam)


def eigen_heating(d, t):
    """Late approach to steady state, manuscript Section 3.4.

    Restricted to the deep and outlet probes; see the module docstring.
    """
    lam = []
    for s in SENS_SLOW:
        y = d[s]
        Tss = ss(y, t)
        u = (Tss - y) / (Tss - y[0])
        m = (u > 0.07) & (u < 0.45)
        if m.sum() < 20:
            continue
        r = linregress(t[m], np.log(np.clip(u[m], 1e-9, None)))
        if r.rvalue ** 2 >= 0.95 and r.slope < 0:
            lam.append(-r.slope)
    return np.array(lam)


def main():
    sm = pd.read_csv(P("steady_state_metrics.csv")).set_index("ID")
    dg = pd.read_csv(P("dimensionless_groups.csv")).set_index("ID")
    main = dg.dropna(subset=["Io_kWm2"])
    _r = _lr(main.q_slpm, main.eps)          # pooled steady correlation eps(q)
    eps_of_q = lambda q: _r.intercept + _r.slope * q

    rows = []
    for ID, (fn, Io) in heating.items():
        d = load(fn); t = d["t"]
        Tamb = ss(d["Tamb"], t)
        q = float(np.mean(d["flow"]))
        mdot = RHO_STD * q / 60000.0
        cp = cp_air(0.5 * (ss(d["T3"], t) + Tamb))
        eps = eps_of_q(q)
        lam = eigen_heating(d, t)
        rows.append(dict(ID=ID, phase="heat", q=q, x=eps * mdot * cp,
                         lam=lam.mean(), lam_sd=lam.std(), n=len(lam)))

    for ID, fn in cooling.items():
        d = load(fn); t = d["t"]
        Tamb = ss(d["Tamb"], t)
        q = float(np.mean(d["flow"]))
        mdot = RHO_STD * q / 60000.0
        cp = cp_air(0.5 * (np.mean(d["T3"]) + Tamb))
        eps = eps_of_q(q)
        lam = eigen_cooling(d, t, Tamb)
        rows.append(dict(ID=ID, phase="cool", q=q, x=eps * mdot * cp,
                         lam=lam.mean(), lam_sd=lam.std(), n=len(lam)))

    ev = pd.DataFrame(rows)
    ev = ev[ev.ID.isin(list(dg.dropna(subset=["Io_kWm2"]).index) + list(PARENT))]
    ev.to_csv(P("eigenvalue_verification.csv"), index=False)

    # ---- delivered-power closure --------------------------------------
    out = sm.reset_index().copy()
    out = out.merge(dg.reset_index()[["ID", "eps", "Re"]], on="ID", how="inner")
    out = out.dropna(subset=["Io_kWm2"])
    Tw = sum(w * out[f"{s}_ss"] for s, w in WTS.items())
    for K in K_BRACKET:
        out[f"f_K{K:g}"] = (out.Q_gas_W + K * (Tw - out.Tamb)) / out.Q_in_W
    out["f_applied"] = [F_DEL.get(io * 1e3, np.nan) for io in out.Io_kWm2]
    out["eta_delivered"] = out.Q_gas_W / (out.f_applied * out.Q_in_W)
    out.to_csv(P("delivered_power_check.csv"), index=False)

    Tw_ = sum(w * out[f"{s_}_ss"] for s_, w in WTS.items())
    print("\n  delivered-power factors f by group, at each end of the K bracket:")
    for K in K_BRACKET:
        f = (out.Q_gas_W + K * (Tw_ - out.Tamb)) / out.Q_in_W
        g = f.groupby(out.Io_kWm2).mean()
        print(f"    K = {K:.3f} W/K : " +
              "  ".join(f"{int(i)}: {v:.3f}" for i, v in g.items()))
    print()
    print(f"eigenvalue_verification.csv : {len(ev)} rows "
          f"({(ev.phase=='heat').sum()} heating, {(ev.phase=='cool').sum()} cooling)")
    print(f"delivered_power_check.csv   : {len(out)} rows")
    for tag in ("cool", "heat"):
        s = ev[ev.phase == tag]
        r = linregress(s.x, s.lam)
        print(f"  {tag}: C_eff = {1/r.slope:6.0f} J/K   "
              f"K_loss = {r.intercept/r.slope:.3f} W/K   r2 = {r.rvalue**2:.3f}")


if __name__ == "__main__":
    main()
