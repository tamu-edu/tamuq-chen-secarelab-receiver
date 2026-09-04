"""Corrected reduction pipeline for the SiC honeycomb volumetric receiver.

Single entry point: raw logger files -> every reported quantity, uncertainty,
figure and table. Supersedes the exp_analysis.py / dimensionless_analysis.py /
eigenvalues_and_power.py / uncertainty_analysis.py / make_tables.py chain,
whose two-pass constant dependency (F_DEL defined in script 1 but produced by
script 3) and unreproducible heating eigenvalues are the defects this replaces.

Corrections relative to the previous chain
------------------------------------------
C1  mdot uses the controller REFERENCE density (21.1 C, 1 atm -> 1.1996 kg/m3),
    not rho(T_amb). The Aalborg units are mass-flow devices reporting standard
    L/min, so the reference density is the correct conversion. Was a 0.7-1.9%
    error in mdot and therefore in Re, Nu, Q_gas and f.
C2  Gas enthalpy by trapezoidal integration of cp(T), matching the stated
    method, rather than cp(T_mean)*dT.
C3  Heating eigenvalues recomputed with an explicit, declared sensor rule and
    a window sweep. Both the six-sensor and the deep-probe subsets are
    reported, because their disagreement (a factor ~2.3 in C_eff) is a result.
C4  C_eff / K_loss identified from a joint fit to all 18 eigenvalues as well as
    the two subsets, with the cooling-only fit reported as the primary value
    and its 1 degree of freedom stated.
C5  Inversion crossings located by interpolation between the bracketing runs,
    with the global linear regression carried as an estimator sensitivity. The
    256 kW/m2 crossing is reported as bounded, not resolved, because the
    receiver was already inverted at every tested flow at that flux.
C6  LTNE indices fitted per irradiance group (common slope, flux-dependent
    offset) as well as pooled.
C7  T3 systematic-bias sensitivity propagated through every T3-derived
    quantity as a declared band, separate from the Monte Carlo interval.
C8  Reference-probe sensitivity: eps and Nu recomputed against the interior
    probes and against single wall probes, to quantify how far an identified
    coefficient moves with the choice of reference sensor.
C9  Delivered aperture irradiance reported as an interval bracketed by the
    radiometric flux map and the steady energy closure, with the T3 bias that
    would reconcile them computed explicitly.
C10 Pressure drop extracted from the raw logs and compared against the laminar
    square-duct prediction and the transducer resolution.

Usage:  python receiver_reduction.py [--raw DIR] [--out DIR]
"""
from __future__ import annotations
import argparse, json, os
import numpy as np
import pandas as pd
from scipy.stats import linregress
from scipy.optimize import brentq, least_squares

# --------------------------------------------------------------------------
# campaign definition
# --------------------------------------------------------------------------
HEATING = {
    "E67": ("Data_FPT0067_231125_161757", 456e3), "E68": ("Data_FPT0068_231126_115725", 456e3),
    "E69": ("Data_FPT0069_231126_140153", 456e3), "E70": ("Data_FPT0070_231127_090339", 456e3),
    "E71": ("Data_FPT0071_231128_102707", 456e3), "E72": ("Data_FPT0072_231129_104140", 304e3),
    "E73": ("Data_FPT0073_231129_132744", 304e3), "E74": ("Data_FPT0074_231130_123228", 304e3),
    "E75": ("Data_FPT0075_231201_162138", 304e3), "E76": ("Data_FPT0076_231203_120521", 304e3),
    "E77": ("Data_FPT0077_231203_161315", 256e3), "E78": ("Data_FPT0078_231204_132252", 256e3),
    "E79": ("Data_FPT0079_231204_172244", 256e3), "E80": ("Data_FPT0080_231205_095122", 256e3),
    "E81": ("Data_FPT0081_231205_135354", 256e3),
}
# replicates: E82 repeats E77, E83 repeats E81. Archived, excluded from all fits.
REPLICATES = {"E82": ("Data_FPT0082_231210_130825", 256e3),
              "E83": ("Data_FPT0083_231211_122053", 256e3)}
COOLING = {"C69": "Data_FPT0069-Cooling_231126_153148",
           "C80": "Data_FPT0080-cooling_231205_112837",
           "C81": "Data_FPT0081-cooling_231205_153409"}
FLUXES = [456, 304, 256]

# --------------------------------------------------------------------------
# geometry, instrument constants
# --------------------------------------------------------------------------
W_CH, T_WEB, N_CH, L_REC = 1.5e-3, 0.4e-3, 100, 0.137
SIDE = 10 * (W_CH + T_WEB)
A_FRT = SIDE ** 2
A_CH = W_CH ** 2
D_H = W_CH
PER = 4 * W_CH
POROSITY = N_CH * A_CH / A_FRT
A_SOLID = A_FRT - N_CH * A_CH
M_MONO = 0.040
K_SIC = 40.0
Z_WALL = {"T8": 0.011, "T12": 0.058, "T11": 0.107}   # side-wall chain
Z_INT = {"T9": 0.058, "T10": 0.107}                  # interior, flow-exposed
RHO_STD = 101325.0 / (287.05 * 294.25)               # C1: Aalborg 21.1 C, 1 atm
DP_FS, DP_ACC = 200.0, 0.001                         # Keller PD33X: +-200 mbar, 0.1% FS
DT3_BAND = 25.0                                      # C7: declared T3 systematic band [K]

# wall quadrature: length average of the piecewise-linear profile through the
# three side-wall probes with constant extrapolation to the end faces
_zw = [Z_WALL[k] for k in ("T8", "T12", "T11")]
_bnd = [0.0] + [0.5 * (_zw[i] + _zw[i + 1]) for i in range(len(_zw) - 1)] + [L_REC]
WTS = {k: (_bnd[i + 1] - _bnd[i]) / L_REC for i, k in enumerate(("T8", "T12", "T11"))}

_T = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300], float)
_CP = np.array([1007, 1014, 1030, 1051, 1075, 1099, 1121, 1141, 1159, 1175, 1189], float)
_MU = np.array([1.85, 2.30, 2.70, 3.06, 3.38, 3.69, 3.98, 4.25, 4.51, 4.75, 4.98]) * 1e-5
_K = np.array([26.3, 33.8, 40.7, 46.9, 52.4, 57.3, 62.0, 66.7, 71.5, 76.3, 82.0]) * 1e-3

cp_air = lambda T: np.interp(T, _T, _CP)
mu_air = lambda T: np.interp(T, _T, _MU)
k_air = lambda T: np.interp(T, _T, _K)


def h_gas(T_lo, T_hi, n=64):
    """C2: specific enthalpy difference by trapezoidal integration of cp(T)."""
    Tg = np.linspace(T_lo, T_hi, n)
    return float(np.trapezoid(cp_air(Tg), Tg))


# --------------------------------------------------------------------------
# raw file handling
# --------------------------------------------------------------------------
COLS = dict(t=1, mfc=(6, 7, 8, 9), dp1=16, dp2=17,
            T2=35, T3=36, T8=41, T9=42, T10=43, T11=44, T12=45, T15=48, T16=49)
RHO_MFC_CASES = (1.0, 0.0)   # fully correlated / independent controllers
# Provenance of each cooling transient: the heating run it decays from. Fixed
# by the experiment log, not inferred by flow proximity -- E81 and E76 differ
# by 0.001 sL/min, so a nearest-flow match is unstable under perturbation.
COOL_PROVENANCE = {"C69": "E69", "C80": "E80", "C81": "E81"}
SENSORS = ["T2", "T3", "T8", "T9", "T10", "T11", "T12"]

# uncertainty model constants (see monte_carlo docstring D1, D2)
# D1b (2026-09-04, revised after review v9): the controller error is TWO-TERM.
# The Aalborg GFC17 specification is +-1.0% of FULL SCALE over 0-100% of range
# with +-0.5% FS repeatability -- an ADDITIVE bound, not a percent of reading.
# The units were calibrated in situ against a bubble flowmeter because the
# factory calibration is gas-specific, which removes the gas-conversion bias
# (a percent-of-reading systematic) but cannot remove the unit's own zero,
# linearity and repeatability floor. Both terms are therefore carried:
#   sigma_unit = sqrt( (MFC_A_FS * MFC_FS)^2 + (MFC_B_REL * reading)^2 )
# with the additive term from the specification and the proportional term from
# the in-house calibration residual. A purely proportional model understates
# the error at low flow and, because the two terms have opposite flow
# dependence, the choice acts on a fitted exponent and not only on a prefactor.
MFC_FS = 10.0           # per-unit full scale [sL/min]; GFC17 max air/N2 range
MFC_A_FS = 0.0025       # additive 1-sigma as a fraction of FS (0.5% FS repeatability, 95% bound)
MFC_B_REL = 0.025       # proportional 1-sigma of reading (bubble-flowmeter calibration residual)
Q_REL_SD = MFC_B_REL    # retained name for the proportional term
TOL_COVERAGE = 2.0      # IEC 60584 class-1 tolerance read as a 95% bound


def mfc_sigma(reading):
    """Per-unit 1-sigma absolute flow error [sL/min] at a given reading."""
    r = np.asarray(reading, float)
    return np.sqrt((MFC_A_FS * MFC_FS) ** 2 + (MFC_B_REL * r) ** 2)


def mfc_rel_perturbation(shares, q_total, z_unit):
    """Relative perturbation of a summed flow from persistent per-unit errors.

    shares: (n, 4) measured share of each controller for each record.
    q_total: (n,) metered total flow [sL/min].
    z_unit: (4,) standard normal draws, persistent across records within one
            Monte Carlo realization (one physical error per controller).
    """
    shares = np.asarray(shares, float)
    q_total = np.asarray(q_total, float)
    reading = shares * q_total[:, None]
    dq = (mfc_sigma(reading) * np.asarray(z_unit, float)[None, :]).sum(axis=1)
    return dq / q_total


def load(raw_dir, fname):
    df = pd.read_csv(os.path.join(raw_dir, fname + ".csv"), header=0,
                     low_memory=False, encoding="latin-1")
    d = df.apply(pd.to_numeric, errors="coerce")
    t = d.iloc[:, COLS["t"]].to_numpy(float)
    keep = ~np.isnan(t)
    d, t = d[keep], t[keep] - t[keep][0]
    K = lambda i: d.iloc[:, i].to_numpy(float) + 273.15
    _mfc = [d.iloc[:, i].to_numpy(float) for i in COLS["mfc"]]
    out = {"t": t, "flow": sum(_mfc), "mfc": _mfc,
           "dp1": d.iloc[:, COLS["dp1"]].to_numpy(float),
           "dp2": d.iloc[:, COLS["dp2"]].to_numpy(float)}
    out.update({s: K(COLS[s]) for s in SENSORS})
    out["Tamb"] = 0.5 * (K(COLS["T15"]) + K(COLS["T16"]))
    return out


def tail_mean(x, t, win=120.0):
    return float(np.nanmean(np.asarray(x)[t >= t[-1] - win]))


# --------------------------------------------------------------------------
# steady-state reduction
# --------------------------------------------------------------------------
def reduce_steady(raw_dir, runs):
    rows = []
    for ID, (fn, Io) in runs.items():
        d = load(raw_dir, fn)
        t = d["t"]
        q = tail_mean(d["flow"], t)
        Tamb = tail_mean(d["Tamb"], t)
        S = {s: tail_mean(d[s], t) for s in SENSORS}
        mdot = RHO_STD * q / 60000.0                                    # C1
        Tw = sum(WTS[k] * S[k] for k in WTS)
        Qgas = mdot * h_gas(Tamb, S["T3"])                              # C2
        # D1: the four measured controller shares for this run. The Monte Carlo
        # draws four PERSISTENT unit errors per realization and combines them
        # with these shares, so one physical unit error propagates to every
        # run (steady and transient) with weights that differ only because the
        # shares do.
        _sh = np.array([tail_mean(m, t) for m in d["mfc"]])
        _f = (_sh / _sh.sum()) if _sh.sum() > 0 else np.full(4, 0.25)
        _rss = float(np.sqrt((_f ** 2).sum()))
        rows.append(dict(ID=ID, Io_kWm2=Io / 1e3, q_slpm=q, Tamb=Tamb, dur_s=t[-1],
                         mfc_rss=_rss,
                         **{f"mfc_f{k+1}": float(_f[k]) for k in range(4)},
                         mdot_gs=mdot * 1e3, Tw_K=Tw, Q_gas_W=Qgas,
                         Q_nom_W=Io * A_FRT, dp1_mbar=tail_mean(d["dp1"], t),
                         dp2_mbar=tail_mean(d["dp2"], t),
                         **{f"{s}_ss": S[s] for s in SENSORS}))
    return pd.DataFrame(rows)


def dimensionless(ss, dT3=0.0):
    """Per-run dimensionless groups. dT3 shifts T3 by a systematic bias (C7)."""
    o = ss.copy()
    T3 = o.T3_ss + dT3
    o["Tg_bar"] = 0.5 * (o.Tamb + T3)
    o["mdot_ch"] = o.mdot_gs * 1e-3 / N_CH
    o["Re"] = o.mdot_ch * D_H / (A_CH * mu_air(o.Tg_bar))
    o["Pr"] = cp_air(o.Tg_bar) * mu_air(o.Tg_bar) / k_air(o.Tg_bar)
    o["Gz_L"] = D_H * o.Re * o.Pr / L_REC
    o["Pe_LD"] = o.Re * o.Pr * L_REC / D_H
    o["eps"] = (T3 - o.Tamb) / (o.Tw_K - o.Tamb)
    o["NTU"] = -np.log(1.0 - o.eps)
    o["h_app"] = o.NTU * o.mdot_ch * cp_air(o.Tg_bar) / (PER * L_REC)
    o["T_film"] = 0.5 * (o.Tw_K + o.Tg_bar)
    o["Nu"] = o.h_app * D_H / k_air(o.T_film)
    o["Bi"] = o.h_app * T_WEB / (2 * K_SIC)
    o["N_rc"] = 4 * 5.670e-8 * o.Tw_K ** 3 * D_H / k_air(o.Tw_K)
    o["Lam58"] = (o.T12_ss - o.T9_ss) / (o.T12_ss - o.Tamb)
    o["Lam107"] = (o.T11_ss - o.T10_ss) / (o.T11_ss - o.Tamb)
    o["I_vol"] = o.T12_ss - o.T8_ss
    o["Q_gas_W"] = o.mdot_gs * 1e-3 * np.array([h_gas(a, b) for a, b in zip(o.Tamb, T3)])
    o["eta_nom"] = o.Q_gas_W / o.Q_nom_W
    return o


# --------------------------------------------------------------------------
# eigenvalue identification
# --------------------------------------------------------------------------
COOL_SENS = ["T8", "T12", "T11", "T9", "T10", "T3"]
DEEP_SENS = ["T11", "T10", "T3"]


def eigen_cooling(d, sensors=COOL_SENS, thresh=5.0):
    t = d["t"]
    Tamb = tail_mean(d["Tamb"], t)
    m0 = t > t[-1] / 2
    lam = []
    for s in sensors:
        th = d[s] - Tamb
        m = m0 & (th > thresh)
        if m.sum() < 20:
            continue
        r = linregress(t[m], np.log(th[m]))
        lam.append(-r.slope)
    return float(np.mean(lam)), float(np.std(lam)), len(lam)


def eigen_heating(d, sensors, u_lo=0.07, u_hi=0.45, r2_min=0.95):
    """log(T_ss - T) regressed against time on the declared deficit window."""
    t = d["t"]
    lam = []
    for s in sensors:
        y = d[s]
        Tss, T0 = tail_mean(y, t), float(np.mean(y[:10]))
        if Tss - T0 < 20:
            continue
        u = (Tss - y) / (Tss - T0)
        m = (u > u_lo) & (u < u_hi) & (Tss - y > 0)
        if m.sum() < 20:
            continue
        r = linregress(t[m], np.log((Tss - y)[m]))
        if r.rvalue ** 2 >= r2_min and r.slope < 0:
            lam.append(-r.slope)
    if not lam:
        return np.nan, np.nan, 0
    return float(np.mean(lam)), float(np.std(lam)), len(lam)


def identify(x, lam):
    """lam = (x + K_loss)/C_eff  with x = eps*mdot*cp -> slope 1/C, intercept K/C."""
    r = linregress(x, lam)
    return 1.0 / r.slope, r.intercept / r.slope, r.rvalue ** 2, r.stderr / r.slope ** 2


# --------------------------------------------------------------------------
# inversion crossing
# --------------------------------------------------------------------------
def crossings(g):
    """Return crossing flow and effectiveness by local bracketing and by a
    global linear regression (C5)."""
    g = g.sort_values("q_slpm")
    q, iv, ep = g.q_slpm.values, g.I_vol.values, g.eps.values
    c = np.polyfit(q, iv, 1)
    q_glob = -c[1] / c[0]
    e_glob = float(np.polyval(np.polyfit(q, ep, 1), q_glob))
    sgn = np.diff(np.sign(iv))
    if not np.any(sgn):
        return dict(q_local=np.nan, eps_local=np.nan, q_global=q_glob,
                    eps_global=e_glob, bracketed=False, q_min=q.min(), q_max=q.max())
    k = int(np.where(sgn)[0][0])
    fr = -iv[k] / (iv[k + 1] - iv[k])
    return dict(q_local=q[k] + fr * (q[k + 1] - q[k]),
                eps_local=ep[k] + fr * (ep[k + 1] - ep[k]),
                q_global=q_glob, eps_global=e_glob,
                bracketed=bool(q.min() <= q_glob <= q.max()),
                q_min=q.min(), q_max=q.max())


# --------------------------------------------------------------------------
# pressure drop
# --------------------------------------------------------------------------
PO_SQ = 56.91   # Darcy friction factor x Re for a square duct


def dp_laminar(mdot_gs, Tg):
    """Monolith-only laminar prediction, all metered flow through N_CH channels."""
    G = (mdot_gs * 1e-3 / N_CH) / A_CH
    rho = 101325.0 / (287.05 * Tg)
    return 0.5 * PO_SQ * mu_air(Tg) * L_REC * (G / rho) / D_H ** 2 / 100.0   # mbar


# --------------------------------------------------------------------------
# delivered power closure and the reconciling T3 bias
# --------------------------------------------------------------------------
def closure(dg, K):
    """f = [Q_gas + K (Tw - Tamb)] / (G0 A_frt), per run."""
    return (dg.Q_gas_W + K * (dg.Tw_K - dg.Tamb)) / dg.Q_nom_W


def reconciling_dT3(dg, K):
    """C9: the uniform T3 bias that would bring the closure to f = 1, i.e. that
    would reconcile the energy closure with the radiometric flux map."""
    out = {}
    for Io, g in dg.groupby("Io_kWm2"):
        need = (g.Q_nom_W - K * (g.Tw_K - g.Tamb) - g.Q_gas_W)     # W of gas enthalpy
        dT = need / (g.mdot_gs * 1e-3 * cp_air(g.Tg_bar))          # K
        out[int(Io)] = (float(dT.mean()), float(dT.min()), float(dT.max()))
    return out


# --------------------------------------------------------------------------
# reference-probe sensitivity  (C8)
# --------------------------------------------------------------------------
def reference_variants(ss):
    """eps / NTU / Nu recomputed against different solid reference temperatures.
    Quantifies how far an 'identified' assembly coefficient moves with the
    choice of reference sensor - the core of the under-instrumentation claim."""
    variants = {
        "wall quadrature (T8,T12,T11)": lambda o: sum(WTS[k] * o[f"{k}_ss"] for k in WTS),
        "interior probes (T9,T10)": lambda o: 0.5 * (o.T9_ss + o.T10_ss),
        "front wall only (T8)": lambda o: o.T8_ss,
        "mid wall only (T12)": lambda o: o.T12_ss,
        "rear wall only (T11)": lambda o: o.T11_ss,
        "wall+interior mean": lambda o: 0.5 * (sum(WTS[k] * o[f"{k}_ss"] for k in WTS)
                                               + 0.5 * (o.T9_ss + o.T10_ss)),
    }
    rows = []
    for name, fn in variants.items():
        o = ss.copy()
        o["Tw_K"] = fn(o)
        dg = dimensionless(o)
        b, a = np.polyfit(np.log(dg.Re), np.log(dg.Nu), 1)
        cr = {int(Io): crossings(g)["eps_local"] for Io, g in dg.groupby("Io_kWm2")}
        rows.append(dict(reference=name, eps_min=dg.eps.min(), eps_max=dg.eps.max(),
                         NTU_min=dg.NTU.min(), NTU_max=dg.NTU.max(),
                         Nu_min=dg.Nu.min(), Nu_max=dg.Nu.max(),
                         Nu_prefactor=np.exp(a), Nu_exponent=b,
                         eps_star_456=cr.get(456), eps_star_304=cr.get(304),
                         eps_star_256=cr.get(256)))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Monte Carlo
# --------------------------------------------------------------------------
def monte_carlo(ss, eig, n=4000, seed=20260902, rho=1.0):
    """Instrument errors: systematic within a run, independent between
    instruments.

    D1 (2026-09-04, revised after review v7): the flow is the SUM OF FOUR
    Aalborg controllers carrying roughly 25/33/28/14% of the total. The 2.5%
    figure is a PER-CONTROLLER accuracy specification and therefore fixes each
    unit's TOTAL variance; it cannot be decomposed into a full-variance common
    component plus a full-variance independent component, which is what the
    previous version did (giving each unit sqrt(2) x 2.5%). Each controller's
    error is instead written as a correlated part of weight sqrt(RHO_MFC) and
    an independent part of weight sqrt(1 - RHO_MFC), so the per-unit total is
    2.5% for any correlation RHO_MFC. Because the manufacturer does not state
    a correlation between units, the pipeline reports the two limiting cases:
    RHO_MFC = 1 (fully correlated, relative error on the summed flow 2.5%,
    purely systematic, displaces a prefactor) and RHO_MFC = 0 (independent,
    relative error 2.5% x sqrt(sum f_i^2) ~ 1.3%, which does act on an
    exponent). The independent part is NOT claimed as run-to-run
    repeatability, which is unspecified for these units; it is the
    between-unit part of the accuracy spec, and it varies between runs only
    because the shares f_i do.

    D2 (2026-09-03/04): the IEC 60584 class-1 tolerance for type N -- the
    greater of 1.5 K and 0.004|t| -- is a maximum permissible error, i.e. a
    bound rather than a standard deviation. We interpret it as a 95% bound and
    divide by TOL_COVERAGE = 2 for the 1-sigma term. The class error is a
    per-sensor CALIBRATION offset and is therefore drawn ONCE PER SENSOR per
    realization and shared across all runs (scaled by each run's temperature-
    dependent tolerance), with an independent 0.5 K per-run drift on top. It
    was previously drawn independently for every run, which is inconsistent
    with calling it a systematic instrument error and mis-states how it
    propagates into a slope.
    """
    rng = np.random.default_rng(seed)
    keep = ["Nu_a", "Nu_b", "Nu_b_grouped", "NTU_corr_b", "eps_star_456",
            "eps_star_304", "eps_star_256", "Lam107_slope", "Lam107_int",
            "C_cool", "K_cool", "C_deep", "K_deep", "C_all", "K_all",
            "C_match", "K_match"]
    acc = {k: [] for k in keep}
    tol = {s: np.maximum(1.5, 0.004 * np.abs(ss[f"{s}_ss"] - 273.15)) / TOL_COVERAGE
           for s in SENSORS}
    for _ in range(n):
        o = ss.copy()
        z_cal_sensor = {}
        for sn in SENSORS:
            z_cal_sensor[sn] = rng.normal(0, 1)           # one calibration draw per sensor
            o[f"{sn}_ss"] = (o[f"{sn}_ss"] + z_cal_sensor[sn] * tol[sn]
                             + rng.normal(0, 0.5, len(o)))
        z_cal_sensor["Tamb"] = rng.normal(0, 1)
        # D1: per-unit accuracy Q_REL_SD split by correlation rho, so the
        # per-unit total variance is Q_REL_SD^2 whatever rho is. Both the
        # common and independent parts are shared with the eigenvalue abscissa
        # below, since they are the same physical controller errors.
        # four PERSISTENT unit errors, drawn once per realization. Each unit's
        # total sd is Q_REL_SD for any rho; rho only sets how much of it the
        # units share. Combined with each record's own measured shares.
        # four PERSISTENT standard-normal unit draws; rho sets how much of the
        # error the units share. The magnitude each draw produces at a given
        # setpoint follows the two-term model mfc_sigma(), so one physical unit
        # error reaches every steady and transient record with the weight its
        # own reading implies.
        z_unit = (np.sqrt(rho) * rng.normal(0, 1)
                  + np.sqrt(1.0 - rho) * rng.normal(0, 1, 4))
        F_ss = ss[[f"mfc_f{k+1}" for k in range(4)]].to_numpy(float)
        o["q_slpm"] = o.q_slpm * (1.0 + mfc_rel_perturbation(F_ss, ss.q_slpm.values, z_unit))
        o["mdot_gs"] = RHO_STD * o.q_slpm / 60000.0 * 1e3
        o["Tw_K"] = sum(WTS[k] * o[f"{k}_ss"] for k in WTS)
        o["Q_gas_W"] = o.mdot_gs * 1e-3 * np.array(
            [h_gas(a2, b2) for a2, b2 in zip(o.Tamb, o.T3_ss)])
        fmu, fk = rng.normal(1, 0.02), rng.normal(1, 0.01)
        dgm = dimensionless(o)
        dgm["Re"] = dgm.Re / fmu
        dgm["Nu"] = dgm.Nu / fk
        try:
            bb, aa = np.polyfit(np.log(dgm.Re), np.log(dgm.Nu), 1)
        except Exception:
            continue
        acc["Nu_a"].append(np.exp(aa)); acc["Nu_b"].append(bb)
        acc["Nu_b_grouped"].append(grouped_powerlaw(dgm, "Re", "Nu")["exponent"])
        try:
            Nc = np.array([_solve_N(r) for _, r in dgm.iterrows()])
            acc["NTU_corr_b"].append(float(np.polyfit(np.log(dgm.Re), np.log(Nc), 1)[0]))
        except Exception:
            pass
        for Io, g in dgm.groupby("Io_kWm2"):
            acc[f"eps_star_{int(Io)}"].append(crossings(g)["eps_local"])
        # common slope with one intercept per irradiance configuration: the
        # pooled single-intercept fit averages the flux-ordered offsets away
        # and inflates the slope standard error ninefold (2026-09-03).
        Dm = pd.get_dummies(dgm.Io_kWm2.astype(int), prefix="I").astype(float)
        Xl = np.column_stack([dgm.Re.values] + [Dm[c].values for c in Dm.columns])
        bl = np.linalg.lstsq(Xl, dgm.Lam107.values, rcond=None)[0]
        acc["Lam107_slope"].append(bl[0]); acc["Lam107_int"].append(float(np.mean(bl[1:])))
        e = eig.copy()
        e["lam"] = e.lam + rng.normal(0, 1, len(e)) * e.lam_sd.fillna(0)
        # D3 (2026-09-03): rebuild the abscissa x = eps*mdot*cp from perturbed
        # flow and temperature realizations, so that its uncertainty is the
        # propagated instrument uncertainty rather than a blanket 1% multiplier.
        # eps comes from the eps(q) fit refitted on THIS realization's steady
        # data, so the correlation between the abscissa and the steady-state
        # perturbations is retained.
        eq = np.polyfit(dgm.q_slpm, dgm.eps, 1)
        # the cooling runs carry the same controller errors; their shares are
        # not logged per cooling run, so the mean steady share factor is used
        # the same four unit errors reach the transient records, weighted by
        # each cooling/heating run's own measured shares
        F_e = e[[f"mfc_f{k+1}" for k in range(4)]].to_numpy(float)
        e["fq"] = 1.0 + mfc_rel_perturbation(F_e, e.q.values, z_unit)
        q_e = e.q.values * e.fq.values
        # D2: the T3 and ambient calibration offsets are the SAME sensors as in
        # the steady frame, so their per-sensor draws z_cal are reused here
        # rather than redrawn, with an independent drift term per row
        dT_e = (z_cal_sensor["Tamb"]
                * np.maximum(1.5, 0.004 * np.abs(e.Tamb_e - 273.15)) / TOL_COVERAGE
                + z_cal_sensor["T3"]
                * np.maximum(1.5, 0.004 * np.abs(e.T3_e - 273.15)) / TOL_COVERAGE
                + rng.normal(0, 0.5, len(e)) + rng.normal(0, 0.5, len(e)))
        e["x"] = (np.polyval(eq, q_e) * (RHO_STD * q_e / 60000.0)
                  * cp_air(0.5 * (e.Tamb_e + e.T3_e) + 0.5 * dT_e))
        for tag, sub in (("cool", e[e.phase == "cool"]),
                         ("deep", e[e.phase == "heat"]), ("all", e)):
            C, K, _, _ = identify(sub.x.values, sub.lam.values)
            acc[f"C_{tag}"].append(C); acc[f"K_{tag}"].append(K)
        # matched-effectiveness abscissa: eps from the heating run this cooling
        # transient actually decays from (COOL_PROVENANCE), c_p from the
        # cooling run's own tail temperatures
        _m = (e.phase.values == "cool")
        ec = e[_m]
        dT_c = np.asarray(dT_e)[_m]
        xm = []
        for (_, cr), dTi in zip(ec.iterrows(), dT_c):
            hh = dgm[dgm.ID == COOL_PROVENANCE[cr.ID]].iloc[0]
            xm.append(float(hh.eps) * (RHO_STD * float(cr.q) * float(cr.fq) / 60000.0)
                      * float(cp_air(0.5 * (cr.Tamb_e + cr.T3_e) + 0.5 * float(dTi))))
        Cm, Km, _, _ = identify(np.array(xm), ec.lam.values)
        acc["C_match"].append(Cm); acc["K_match"].append(Km)
    rows = []
    for k, v in acc.items():
        v = np.array([x for x in v if np.isfinite(x)])
        rows.append(dict(quantity=k, value=float(np.mean(v)), sd=float(np.std(v)),
                         ci_lo=float(np.percentile(v, 2.5)),
                         ci_hi=float(np.percentile(v, 97.5)), n=len(v)))
    return pd.DataFrame(rows)



# --------------------------------------------------------------------------
# grouped regression and profile-corrected transfer units (added 2026-09-03)
# --------------------------------------------------------------------------
def grouped_powerlaw(dg, xcol, ycol, gcol="Io_kWm2"):
    """Common slope in log-log with one intercept per group. The pooled
    single-intercept fit averages group-ordered offsets away and inflates the
    slope standard error ninefold on this dataset."""
    x = np.log(dg[xcol].values); y = np.log(dg[ycol].values)
    Dm = pd.get_dummies(dg[gcol].astype(int), prefix="g").astype(float)
    X = np.column_stack([x] + [Dm[c].values for c in Dm.columns])
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    res = y - X @ beta
    se = float(np.sqrt((res ** 2).sum() / (len(y) - X.shape[1])
                       * np.linalg.inv(X.T @ X)[0, 0]))
    r2 = float(1 - (res ** 2).sum() / ((y - y.mean()) ** 2).sum())
    return dict(exponent=float(beta[0]), stderr=se, r2=r2,
                prefactors={c.replace("g_", ""): float(np.exp(bi))
                            for c, bi in zip(Dm.columns, beta[1:])})


def _Tg_exit(N, T8, T12, T11, Tamb, m=400):
    """Integrate mdot cp dTg/dz = h P (Tw(z) - Tg) with uniform h through the
    measured piecewise-linear wall profile. N = h P L / (mdot cp)."""
    zt = np.linspace(0.0, 1.0, m + 1)
    Tw = np.interp(zt * L_REC, [Z_WALL["T8"], Z_WALL["T12"], Z_WALL["T11"]],
                   [T8, T12, T11])
    Tg, hh = Tamb, 1.0 / m
    for i in range(m):
        k1 = N * (Tw[i] - Tg)
        k2 = N * (Tw[i + 1] - (Tg + hh * k1))
        Tg = Tg + hh * (k1 + k2) / 2.0
    return Tg


def _solve_N(r, dT3=0.0):
    f = lambda N: (_Tg_exit(N, r.T8_ss, r.T12_ss, r.T11_ss, r.Tamb) - (r.T3_ss + dT3))
    return brentq(f, 1e-4, 120.0)


def ntu_profile_corrected(dg, dT3_band=25.0):
    """Transfer-unit count without the isothermal-wall identity.

    NTU_app = -ln(1-eps) equals (1/mdot cp) int h P dz only for an isothermal
    wall. The wall here spans several hundred K and changes its axial shape
    with flow, so instead solve for the N that reproduces the measured outlet
    temperature when integrated through the measured wall profile.
    """
    def solve(frame, dT3=0.0):
        return np.array([_solve_N(r, dT3) for _, r in frame.iterrows()])

    N = solve(dg)
    d = dg.copy(); d["NTU_corr"] = N
    r = linregress(np.log(d.Re), np.log(d.NTU_corr))
    band = sorted(float(linregress(np.log(dg.Re), np.log(solve(dg, sg))).slope)
                  for sg in (-dT3_band, +dT3_band))
    return dict(exponent=float(r.slope), stderr=float(r.stderr),
                r2=float(r.rvalue ** 2),
                grouped=grouped_powerlaw(d, "Re", "NTU_corr"),
                ratio_to_NTU_app=[float((d.NTU_corr / d.NTU).min()),
                                  float((d.NTU_corr / d.NTU).max())],
                exponent_T3_band=band, single_stream_requirement=-1.0,
                NTU_corr=[float(v) for v in N])



# --------------------------------------------------------------------------
# fixed-conductance falsification  (2026-09-04; added to the pipeline in
# response to review v7, which correctly noted it existed only in scratch)
# --------------------------------------------------------------------------
def _wall_z(r, zt, rear="const", front="const"):
    """Wall temperature on the normalized grid zt, with a choice of
    extrapolation over the unmeasured front 11 mm and rear 30 mm."""
    z = zt * L_REC
    zp = [Z_WALL["T8"], Z_WALL["T12"], Z_WALL["T11"]]
    T = np.interp(z, zp, [r.T8_ss, r.T12_ss, r.T11_ss])
    sl_r = (r.T11_ss - r.T12_ss) / (zp[2] - zp[1])
    sl_f = (r.T12_ss - r.T8_ss) / (zp[1] - zp[0])
    if rear == "linear":
        T = np.where(z > zp[2], r.T11_ss + sl_r * (z - zp[2]), T)
    elif rear == "half":
        T = np.where(z > zp[2], r.T11_ss + 0.5 * sl_r * (z - zp[2]), T)
    if front == "linear":
        T = np.where(z < zp[0], r.T8_ss - sl_f * (zp[0] - z), T)
    return T


def _Tg_exit_h(hz_nodes, nodes, r, hscale, rear="const", front="const", m=400):
    """Integrate mdot cp dTg/dz = h(z) P (Tw - Tg) for a piecewise-linear
    h(z) on `nodes` (normalized).  hscale = P L / (mdot_ch cp) for this run,
    so N(z) = h(z) * hscale."""
    zt = np.linspace(0.0, 1.0, m + 1)
    Tw = _wall_z(r, zt, rear, front)
    hz = np.interp(zt, nodes, hz_nodes) * hscale
    Tg, dz = r.Tamb, 1.0 / m
    # The RK2 step is unstable for hz*dz > 2, and the multi-start optimizer
    # does probe log-h values that far. Clip the local transfer-unit density at
    # 1.5/dz, i.e. N ~ 600 per unit zeta: there the gas is already at the wall
    # temperature to machine precision (1 - exp(-600)), so the clip changes no
    # attainable outlet temperature while keeping the integration finite. This
    # removes the optimizer overflow warnings without altering any result.
    hz = np.clip(hz, 0.0, 1.5 / dz)
    for i in range(m):
        k1 = hz[i] * (Tw[i] - Tg)
        k2 = hz[i + 1] * (Tw[i + 1] - (Tg + dz * k1))
        Tg = Tg + dz * (k1 + k2) / 2.0
    return Tg


def fixed_profile_test(dg, n_starts=40, seed=20260904):
    """Can ANY flow-independent axial conductance profile reproduce the fifteen
    measured outlet temperatures?

    Three families are tested, all with non-negative piecewise-linear nodes and
    multi-start least squares in log-h:
      shared_h    - one dimensional h(z) common to all fifteen runs
      per_flux    - one h(z) per irradiance configuration (five nodes), which
                    allows conductance to differ with lamp setting/temperature
      shared_Nu   - one Nu(z) common to all runs, h_i(z) = Nu(z) k_i / D_h,
                    the direct test of the constant-Nusselt model class
    Reported: rms outlet residual, and the residual's regression on ln Re,
    which is the discriminating statistic (a flow-independent profile that
    merely mis-fits gives scattered residuals; one that is the wrong FORM
    gives residuals ordered in flow).
    """
    rng = np.random.default_rng(seed)
    cp = cp_air(dg.Tg_bar.values)
    kk = k_air(dg.Tg_bar.values)
    sc = PER * L_REC / (dg.mdot_ch.values * cp)             # h -> N
    sc_nu = sc * kk / D_H                                    # Nu -> N
    lnRe = np.log(dg.Re.values)

    def _fit(rows_idx, nn, scale):
        nodes = np.linspace(0.0, 1.0, nn)
        sub = dg.iloc[rows_idx]
        # scale-aware starting range: the uniform-conductance solution for these
        # runs sets the order of magnitude. Without this the Nu(z) family, whose
        # values are ~200x smaller than the dimensional h(z) family, is started
        # far from any minimum and the multi-start does not converge.
        c0 = float(np.mean([_solve_N(r) / scale[j]
                            for r, j in zip([sub.iloc[a] for a in range(len(sub))],
                                            rows_idx)]))

        def res(p):
            h = np.exp(p)
            return np.array([_Tg_exit_h(h, nodes, r, scale[j], "const", "const")
                             - r.T3_ss for r, j in zip(
                                 [sub.iloc[a] for a in range(len(sub))], rows_idx)])
        best = None
        for _ in range(n_starts):
            x0 = np.log(c0 * rng.uniform(0.1, 10.0, nn))
            try:
                so = least_squares(res, x0=x0, xtol=1e-13, ftol=1e-13,
                                   max_nfev=3000)
            except Exception:
                continue
            if best is None or so.cost < best.cost:
                best = so
        rr = res(best.x)
        return np.exp(best.x), rr

    out = {"shared_h": {}, "per_flux": {}, "shared_Nu": {}}
    allidx = list(range(len(dg)))
    out["run_ids"] = [str(v) for v in dg.ID.values]
    for nn in (2, 3, 5, 7):
        h, rr = _fit(allidx, nn, sc)
        sl = linregress(lnRe, rr)
        out["shared_h"][str(nn)] = dict(
            nodes=[float(v) for v in h], rms_K=float(np.sqrt((rr ** 2).mean())),
            max_abs_K=float(np.abs(rr).max()), r_lnRe=float(sl.rvalue),
            slope_K_per_lnRe=float(sl.slope),
            residuals_K=[float(v) for v in rr])
    for Io, g in dg.groupby("Io_kWm2"):
        idx = [allidx[i] for i in range(len(dg)) if dg.iloc[i].Io_kWm2 == Io]
        h, rr = _fit(idx, 5, sc)
        sl = linregress(np.log(g.Re.values), rr)
        out["per_flux"][str(int(Io))] = dict(
            nodes=[float(v) for v in h], rms_K=float(np.sqrt((rr ** 2).mean())),
            max_abs_K=float(np.abs(rr).max()), r_lnRe=float(sl.rvalue),
            slope_K_per_lnRe=float(sl.slope), n=int(len(g)),
            run_ids=[str(v) for v in g.ID.values],
            residuals_K=[float(v) for v in rr])
    for nn in (3, 5):
        h, rr = _fit(allidx, nn, sc_nu)
        sl = linregress(lnRe, rr)
        out["shared_Nu"][str(nn)] = dict(
            nodes=[float(v) for v in h], rms_K=float(np.sqrt((rr ** 2).mean())),
            max_abs_K=float(np.abs(rr).max()), r_lnRe=float(sl.rvalue),
            slope_K_per_lnRe=float(sl.slope),
            residuals_K=[float(v) for v in rr])
    return out


def wall_extrapolation_sensitivity(dg):
    """Exponent of the profile-corrected transfer-unit count under each
    admissible reconstruction of the unmeasured front 11 mm and rear 30 mm.
    A full linear rear continuation is inadmissible where it places the exit
    wall below the measured T3, which no single-stream model can produce."""
    zt = np.linspace(0.0, 1.0, 401)
    out = {}
    for rear in ("const", "half", "linear"):
        for front in ("const", "linear"):
            Ns, infeasible = [], []
            for _, r in dg.iterrows():
                Tw = _wall_z(r, zt, rear, front)
                if Tw[-1] <= r.T3_ss:
                    infeasible.append(r.ID); continue

                def f(N, r=r, Tw=Tw):
                    Tg, dz = r.Tamb, 1.0 / (len(zt) - 1)
                    for i in range(len(zt) - 1):
                        k1 = N * (Tw[i] - Tg)
                        k2 = N * (Tw[i + 1] - (Tg + dz * k1))
                        Tg = Tg + dz * (k1 + k2) / 2.0
                    return Tg - r.T3_ss
                Ns.append(brentq(f, 1e-4, 120.0))
            key = f"{rear}_{front}"
            if len(Ns) == len(dg):
                m = dg.iloc[[i for i in range(len(dg))]]
                rg = linregress(np.log(m.Re.values), np.log(np.array(Ns)))
                out[key] = dict(exponent=float(rg.slope), stderr=float(rg.stderr),
                                n=len(Ns), infeasible=[])
            else:
                out[key] = dict(exponent=None, n=len(Ns),
                                infeasible=infeasible,
                                note="exit wall below measured T3 for these runs")
    return out


# --------------------------------------------------------------------------
# tables (generated from the same pass, 2026-09-04: they were previously
# hand-maintained and had drifted from the outputs)
# --------------------------------------------------------------------------
def write_tables(rep, dg, mc, P):
    m = mc.set_index("quantity")
    g = lambda k, c="value": float(m.loc[k, c])
    # units are carried in the header cells rather than on a second header row:
    # a table may have only ONE header line before the delimiter row, and the
    # delimiter must have exactly as many cells as the header (2026-09-04).
    _t1cols = ["Run", "$G_0$ [kW m$^{-2}$]", "$q$ [sL min$^{-1}$]", "$Re_{\\rm nom}$",
               "$Gz_L$", "$\\bar T_w$ [K]", "$T_3$ [K]", "$\\varepsilon$",
               "$NTU_{\\rm app}$", "$N_{\\rm prof}$", "$Nu_{\\rm app}$",
               "$\\Lambda_{58}$", "$\\Lambda_{107}$", "$T_{12}-T_8$ [K]",
               "$\\eta_{\\rm nom}$"]
    T1 = ["| " + " | ".join(_t1cols) + " |",
          "|" + "|".join(["---"] + ["---:"] * (len(_t1cols) - 1)) + "|"]
    npf = rep["ntu_profile"]["NTU_corr"]
    for (_, r), nv in zip(dg.sort_values(["Io_kWm2", "q_slpm"]).iterrows(),
                          [npf[i] for i in dg.sort_values(["Io_kWm2", "q_slpm"]).index]):
        T1.append(f"| {r.ID} | {r.Io_kWm2:.0f} | {r.q_slpm:.2f} | {r.Re:.1f} | {r.Gz_L:.3f} | "
                  f"{r.Tw_K:.0f} | {r.T3_ss:.0f} | {r.eps:.3f} | {r.NTU:.3f} | {nv:.3f} | "
                  f"{r.Nu:.4f} | {r.Lam58:.4f} | {r.Lam107:.4f} | {r.T12_ss - r.T8_ss:+.1f} | {r.eta_nom:.3f} |")
    open(P("table1_envelope.md"), "w").write("\n".join(T1) + "\n")

    nu, pr = rep["nusselt"], rep["ntu_profile"]
    idn = rep["identification"]
    def ci(k, sc=1.0, f="{:.3g}"):
        try:
            lo, hi = float(m.loc[k, "ci_lo"]) * sc, float(m.loc[k, "ci_hi"]) * sc
            return f"[{f.format(lo)}, {f.format(hi)}]"
        except Exception:
            return "\u2014"

    nu, pr = rep["nusselt"], rep["ntu_profile"]
    idn = rep["identification"]
    rows = [
        ("$Nu_{\\rm app}$ prefactor $a$ (pooled)", f"{g('Nu_a')*1e4:.2f}$\\times10^{{-4}}$",
         f"{g('Nu_a','sd')*1e4:.2f}$\\times10^{{-4}}$", ci("Nu_a", 1e4, "{:.2f}"), "$\\times10^{-4}$",
         "15 steady runs"),
        ("$Nu_{\\rm app}$ exponent, pooled", f"{nu['exponent']:.3f}", f"{g('Nu_b','sd'):.3f}",
         ci("Nu_b", 1.0, "{:.3f}"), "\u2013",
         f"instrumental MC; regression SE $\\pm${nu['stderr_exponent']:.3f}, $r^2$={nu['r2']:.3f}"),
        ("$Nu_{\\rm app}$ exponent, grouped (primary)", f"{nu['grouped']['exponent']:.3f}",
         f"{g('Nu_b_grouped','sd'):.3f}", ci("Nu_b_grouped", 1.0, "{:.3f}"), "\u2013",
         f"instrumental MC; regression SE $\\pm${nu['grouped']['stderr']:.3f}, $r^2$={nu['grouped']['r2']:.4f}"),
        ("$N_{\\rm prof}$ exponent (primary)", f"{pr['exponent']:+.3f}", f"{g('NTU_corr_b','sd'):.3f}",
         ci("NTU_corr_b", 1.0, "{:.3f}"), "\u2013",
         f"instrumental MC; regression SE $\\pm${pr['stderr']:.3f}; fixed-$Nu$ requirement "
         f"{rep['ntu_structure']['fixed_Nu_requirement']:+.3f}"),
        ("$NTU_{\\rm app}$ exponent (identity, superseded)", f"{rep['ntu_structure']['exponent']:+.3f}",
         "\u2014", "\u2014", "\u2013", "isothermal-wall identity; retained for comparison only"),
    ]
    for Io in (456, 304, 256):
        rows.append((f"Inversion marker $\\varepsilon^*$, {Io} kW m$^{{-2}}$", f"{g(f'eps_star_{Io}'):.3f}",
                     f"{g(f'eps_star_{Io}','sd'):.3f}", ci(f"eps_star_{Io}", 1.0, "{:.3f}"), "\u2013",
                     "operational marker under the adopted wall convention; see \u00a75.1"))
    rows += [
        ("$\\Lambda_{107}$ slope [$Re^{-1}$]", f"{g('Lam107_slope')*1e4:.2f}$\\times10^{{-4}}$",
         f"{g('Lam107_slope','sd')*1e4:.2f}$\\times10^{{-4}}$", ci("Lam107_slope", 1e4, "{:.2f}"),
         "$\\times10^{-4}$", "common slope, per-flux intercepts"),
        ("$C_{\\rm eff}$, cooling matched-$\\varepsilon$ (primary)", f"{idn['cooling_matched_eps']['C_eff']:.0f}",
         f"{g('C_match','sd'):.0f}", ci("C_match", 1.0, "{:.0f}"), "J K$^{-1}$",
         f"$n$=3, $r^2$={idn['cooling_matched_eps']['r2']:.3f}"),
        ("$C_{\\rm eff}$, cooling pooled-$\\varepsilon$", f"{idn['cooling']['C_eff']:.0f}",
         f"{g('C_cool','sd'):.0f}", ci("C_cool", 1.0, "{:.0f}"), "J K$^{-1}$",
         f"$r^2$={idn['cooling']['r2']:.3f}"),
        ("$C_{\\rm eff}$, joint 18 eigenvalues", f"{idn['joint']['C_eff']:.0f}", f"{g('C_all','sd'):.0f}",
         ci("C_all", 1.0, "{:.0f}"), "J K$^{-1}$", f"$r^2$={idn['joint']['r2']:.3f}"),
        ("$C_{\\rm eff}$, heating deep probes", f"{idn['heating_deep']['C_eff']:.0f}", f"{g('C_deep','sd'):.0f}",
         ci("C_deep", 1.0, "{:.0f}"), "J K$^{-1}$",
         f"{idn['heating_all6']['C_eff']:.0f} J K$^{{-1}}$ if all six probes used"),
        ("$K_{\\rm loss}$, cooling matched-$\\varepsilon$ (primary)", f"{idn['cooling_matched_eps']['K_loss']:.3f}",
         f"{g('K_match','sd'):.3f}", ci("K_match", 1.0, "{:.3f}"), "W K$^{-1}$", "secant conductance"),
        ("$K_{\\rm loss}$, heating deep probes", f"{idn['heating_deep']['K_loss']:.3f}", f"{g('K_deep','sd'):.3f}",
         ci("K_deep", 1.0, "{:.3f}"), "W K$^{-1}$", "tangent conductance"),
        ("Monolith capacitance (measured mass)",
         f"{rep['monolith']['C_monolith_600K']:.1f} \u2013 {rep['monolith']['C_monolith_900K']:.1f}",
         "\u2014", "\u2014", "J K$^{-1}$", "40 g $\\times\\,c_p$(600\u2013900 K)"),
    ]
    T2 = ["| Constant | Value | s.d. | 95% interval | Unit | Notes |",
          "|---|---|---|---|---|---|"]
    T2 += [f"| {r0} | {r1} | {r2} | {r3} | {r4} | {r5} |" for r0, r1, r2, r3, r4, r5 in rows]
    open(P("table2_constants.md"), "w").write("\n".join(T2) + "\n")

# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------
def write_supplementary_tables(rep, dg, P):
    """Markdown tables for supplementary sections S3, S4 and S5.

    Written by the pipeline so that no supplementary table is transcribed by
    hand (2026-09-04, review v9). Values come from `rep`, which is the same
    dict archived as results.json, and from the per-run groups frame.
    """
    hw, ss_ = rep["heating_window_sensitivity"], rep["sensor_selection_swing"]
    T = ["| deficit window $u$ | $C_{\\rm eff}$ [J K$^{-1}$] | "
         "$K_{\\rm loss}$ [W K$^{-1}$] | $r^2$ |", "|---|---:|---:|---:|"]
    T += [f"| ({w['u_lo']:.2f}, {w['u_hi']:.2f}) | {w['C_eff']:.1f} | "
          f"{w['K_loss']:.4f} | {w['r2']:.4f} |" for w in hw]
    T += ["", f"Sensor selection: all six receiver probes give "
              f"{ss_['C_all6']:.1f} J K$^{{-1}}$, the three deep and outlet probes "
              f"{ss_['C_deep']:.1f} J K$^{{-1}}$, a ratio of {ss_['ratio']:.2f}."]
    open(P("tableS3_heating_conditionality.md"), "w").write("\n".join(T) + "\n")

    we = rep["wall_extrapolation"]
    T = ["| Rear / front extrapolation | exponent | s.e. | runs | infeasible |",
         "|---|---:|---:|---:|---|"]
    for k in ("const_const", "const_linear", "half_const", "half_linear",
              "linear_const", "linear_linear"):
        v = we[k]
        e = v.get("exponent")
        lab = k.replace("_", " / ").replace("const", "constant")
        se = v.get("stderr")
        e_txt = "—" if e is None else format(e, "+.4f")
        se_txt = "—" if se is None else format(se, ".4f")
        inf_txt = ", ".join(v.get("infeasible") or []) or "none"
        T.append(f"| {lab} | {e_txt} | {se_txt} | {v['n']} | {inf_txt} |")
    fp = rep["fixed_profile_test"]
    T += ["", "| Family | nodes | RMS residual [K] | max abs [K] | "
              "$r$ vs $\\ln Re$ | slope [K per e-fold] |",
          "|---|---:|---:|---:|---:|---:|"]
    for fam, lab in (("shared_h", "shared dimensional $h(z)$"),
                     ("shared_Nu", "shared $Nu(z)$")):
        for n in sorted(fp[fam], key=int):
            v = fp[fam][n]
            T.append(f"| {lab} | {n} | {v['rms_K']:.1f} | {v['max_abs_K']:.1f} | "
                     f"{v['r_lnRe']:+.3f} | {v['slope_K_per_lnRe']:+.0f} |")
    for grp in sorted(fp["per_flux"], key=int):
        v = fp["per_flux"][grp]
        T.append(f"| per-configuration, {grp} kW m$^{{-2}}$ | 5 | {v['rms_K']:.1f} | "
                 f"{v['max_abs_K']:.1f} | {v['r_lnRe']:+.3f} | "
                 f"{v['slope_K_per_lnRe']:+.0f} |")
    open(P("tableS4_wall_and_falsification.md"), "w").write("\n".join(T) + "\n")

    T = ["| Run | $G_0$ [kW m$^{-2}$] | $q$ [sL min$^{-1}$] | $Re_{\\rm nom}$ | "
         "$Pr$ | $Gz_L$ | $x^*=1/Gz_L$ | $Re\\,Pr\\,L/D_h$ | $Bi$ | $N_{rc}$ |",
         "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for _, r in dg.iterrows():
        T.append(f"| {r.ID} | {r.Io_kWm2:.0f} | {r.q_slpm:.2f} | {r.Re:.1f} | "
                 f"{r.Pr:.4f} | {r.Gz_L:.3f} | {1.0 / r.Gz_L:.2f} | {r.Pe_LD:.0f} | "
                 f"{r.Bi:.2e} | {r.N_rc:.2f} |")
    open(P("tableS5_auxiliary_groups.md"), "w").write("\n".join(T) + "\n")


def main(raw_dir, out_dir, n_mc=4000):
    os.makedirs(out_dir, exist_ok=True)
    P = lambda f: os.path.join(out_dir, f)
    rep = {}

    ss = reduce_steady(raw_dir, HEATING)
    ss_rep = reduce_steady(raw_dir, REPLICATES)
    dg = dimensionless(ss)
    dg_rep = dimensionless(ss_rep)
    dg.to_csv(P("groups.csv"), index=False)
    dg_rep.to_csv(P("groups_replicates.csv"), index=False)

    rep["geometry"] = dict(side_mm=SIDE * 1e3, A_frt_cm2=A_FRT * 1e4, porosity=POROSITY,
                           A_solid_m2=A_SOLID, rho_eff=M_MONO / (A_SOLID * L_REC),
                           wall_weights={k: round(v, 4) for k, v in WTS.items()})

    # --- envelope
    rep["envelope"] = {c: [float(dg[c].min()), float(dg[c].max())] for c in
                       ["Re", "Pr", "Gz_L", "eps", "NTU", "Nu", "Bi", "N_rc",
                        "Lam58", "Lam107", "Pe_LD", "eta_nom"]}

    # --- Nusselt law and the structural NTU result
    b, a = np.polyfit(np.log(dg.Re), np.log(dg.Nu), 1)
    r_nu = linregress(np.log(dg.Re), np.log(dg.Nu))
    b_ntu = linregress(np.log(dg.Re), np.log(dg.NTU))
    # Fully developed laminar square-duct Nusselt numbers, Shah & London
    # (1978), aspect ratio 1: (T) 2.976, (H2) 3.091, (H1) 3.608, fRe = 14.227.
    # (H2) is the primary reference; the previous draft used 3.09 without
    # naming the boundary condition.
    NU_FD = {"T": 2.976, "H2": 3.091, "H1": 3.608}
    Nu_fd = NU_FD["H2"]
    _pf = {}
    for Io, g in dg.groupby("Io_kWm2"):
        rg = linregress(np.log(g.Re), np.log(g.Nu))
        _pf[int(Io)] = dict(exponent=float(rg.slope), prefactor=float(np.exp(rg.intercept)),
                            r2=float(rg.rvalue ** 2))
    rep["nusselt"] = dict(prefactor=float(np.exp(a)), exponent=float(b),
                          r2=float(r_nu.rvalue ** 2), stderr_exponent=float(r_nu.stderr),
                          n=len(dg), fd_reference="H2", Nu_fd=NU_FD,
                          ratio_to_fd=[float(Nu_fd / dg.Nu.max()),
                                       float(Nu_fd / dg.Nu.min())],
                          ratio_to_fd_all={k: [float(v / dg.Nu.max()), float(v / dg.Nu.min())]
                                           for k, v in NU_FD.items()},
                          grouped=grouped_powerlaw(dg, "Re", "Nu"),
                          per_flux=_pf,
                          x_star_exit=[float(1.0 / dg.Gz_L.max()), float(1.0 / dg.Gz_L.min())])
    rep["ntu_profile"] = ntu_profile_corrected(dg)
    # the requirement is -1 only if properties are held fixed; computed from
    # the measured k, c_p and mdot it is -0.9996 for a fixed Nusselt number
    _kk = k_air(dg.Tg_bar.values)
    _req_Nu = float(linregress(np.log(dg.Re), np.log(_kk / (dg.mdot_ch * cp_air(dg.Tg_bar)))).slope)
    _req_h = float(linregress(np.log(dg.Re), -np.log(dg.mdot_ch)).slope)
    rep["ntu_structure"] = dict(exponent=float(b_ntu.slope), r2=float(b_ntu.rvalue ** 2),
                                fixed_Nu_requirement=_req_Nu, fixed_h_requirement=_req_h,
                                stderr=float(b_ntu.stderr),
                                single_stream_requirement=-1.0,
                                gap_in_Re_power=float(b_ntu.slope + 1.0))

    # --- inversion crossings
    rep["crossings"] = {int(Io): {k: (float(v) if isinstance(v, (int, float, np.floating))
                                      else v) for k, v in crossings(g).items()}
                        for Io, g in dg.groupby("Io_kWm2")}

    # --- LTNE
    lt = {}
    for tag, col in (("Lam58", "Lam58"), ("Lam107", "Lam107")):
        pooled = linregress(dg.Re, dg[col])
        per = {}
        for Io, g in dg.groupby("Io_kWm2"):
            r = linregress(g.Re, g[col])
            per[int(Io)] = dict(slope=float(r.slope), intercept=float(r.intercept),
                                r2=float(r.rvalue ** 2),
                                range=[float(g[col].min()), float(g[col].max())])
        lt[tag] = dict(pooled=dict(slope=float(pooled.slope),
                                   intercept=float(pooled.intercept),
                                   r2=float(pooled.rvalue ** 2)), per_flux=per)
    rep["ltne"] = lt

    # --- eigenvalues (C3, C4)
    eps_q = np.polyfit(dg.q_slpm, dg.eps, 1)
    rows = []
    for ID, (fn, Io) in HEATING.items():
        d = load(raw_dir, fn)
        t = d["t"]
        q = tail_mean(d["flow"], t)
        Tamb = tail_mean(d["Tamb"], t)
        mdot = RHO_STD * q / 60000.0
        x = float(np.polyval(eps_q, q)) * mdot * cp_air(0.5 * (Tamb + tail_mean(d["T3"], t)))
        _shh = np.array([tail_mean(m, t) for m in d["mfc"]])
        _fh = (_shh / _shh.sum()) if _shh.sum() > 0 else np.full(4, 0.25)
        for tag, sens in (("heat", DEEP_SENS), ("heat6", COOL_SENS)):
            lam, sd, nn = eigen_heating(d, sens)
            rows.append(dict(ID=ID, phase=tag, q=q, x=x, lam=lam, lam_sd=sd, n=nn,
                             Tamb_e=Tamb, T3_e=tail_mean(d["T3"], t),
                             **{f"mfc_f{k+1}": float(_fh[k]) for k in range(4)}))
    for ID, fn in COOLING.items():
        d = load(raw_dir, fn)
        t = d["t"]
        q = float(np.mean(d["flow"]))
        Tamb = tail_mean(d["Tamb"], t)
        mdot = RHO_STD * q / 60000.0
        x = float(np.polyval(eps_q, q)) * mdot * cp_air(0.5 * (Tamb + tail_mean(d["T3"], t)))
        _shc = np.array([np.mean(m) for m in d["mfc"]])
        _fc = (_shc / _shc.sum()) if _shc.sum() > 0 else np.full(4, 0.25)
        lam, sd, nn = eigen_cooling(d)
        rows.append(dict(ID=ID, phase="cool", q=q, x=x, lam=lam, lam_sd=sd, n=nn,
                         Tamb_e=Tamb, T3_e=tail_mean(d["T3"], t),
                         **{f"mfc_f{k+1}": float(_fc[k]) for k in range(4)}))
    eig = pd.DataFrame(rows)
    eig.to_csv(P("eigenvalues.csv"), index=False)

    ident = {}
    for tag, sub in (("cooling", eig[eig.phase == "cool"]),
                     ("heating_deep", eig[eig.phase == "heat"]),
                     ("heating_all6", eig[eig.phase == "heat6"]),
                     ("joint", eig[eig.phase.isin(["cool", "heat"])])):
        sub = sub.dropna(subset=["lam"])
        C, K, r2, _ = identify(sub.x.values, sub.lam.values)
        ident[tag] = dict(C_eff=float(C), K_loss=float(K), r2=float(r2),
                          n=int(len(sub)), dof=int(len(sub) - 2))
    # Matched-effectiveness abscissa for the cooling runs. The abscissa above
    # uses eps from a pooled eps(q) fit over all fifteen heating runs, which
    # averages over the flux ordering of eps. Each cooling run is instead
    # referred to the heating run it actually decays from, taken from
    # COOL_PROVENANCE rather than by flow proximity (2026-09-04: E81 and E76
    # differ by 0.001 sL/min, so a nearest-flow match is unstable under
    # perturbation and could switch between them inside the Monte Carlo).
    # This is the PRIMARY identification.
    _cool = eig[eig.phase == "cool"].dropna(subset=["lam"]).copy()
    _xm = []
    for _, cr in _cool.iterrows():
        hh = dg[dg.ID == COOL_PROVENANCE[cr.ID]].iloc[0]
        _xm.append(hh.eps * (RHO_STD * cr.q / 60000.0)
                   * cp_air(0.5 * (cr.Tamb_e + cr.T3_e)))
    # archive the matched abscissa per cooling run so figures can plot the
    # PRIMARY identification rather than the pooled-eps one
    eig.loc[_cool.index, "x_matched"] = np.array(_xm)
    eig.to_csv(P("eigenvalues.csv"), index=False)
    Cm, Km, r2m, _ = identify(np.array(_xm), _cool.lam.values)
    ident["cooling_matched_eps"] = dict(C_eff=float(Cm), K_loss=float(Km),
                                       r2=float(r2m), n=int(len(_cool)),
                                       dof=int(len(_cool) - 2))
    rep["identification"] = ident
    rep["sensor_selection_swing"] = dict(
        C_deep=ident["heating_deep"]["C_eff"], C_all6=ident["heating_all6"]["C_eff"],
        ratio=ident["heating_deep"]["C_eff"] / ident["heating_all6"]["C_eff"])

    # heating window sensitivity
    win = []
    for lo, hi in [(0.05, 0.35), (0.07, 0.45), (0.10, 0.50), (0.15, 0.60), (0.20, 0.70)]:
        lam = []
        for ID, (fn, Io) in HEATING.items():
            d = load(raw_dir, fn)
            l_, _, _ = eigen_heating(d, DEEP_SENS, lo, hi)
            lam.append(l_)
        sub = eig[eig.phase == "heat"].copy()
        sub["lam"] = lam
        sub = sub.dropna(subset=["lam"])
        C, K, r2, _ = identify(sub.x.values, sub.lam.values)
        win.append(dict(u_lo=lo, u_hi=hi, C_eff=float(C), K_loss=float(K), r2=float(r2)))
    rep["heating_window_sensitivity"] = win

    rep["monolith"] = dict(C_monolith_600K=M_MONO * 1050.0, C_monolith_900K=M_MONO * 1170.0)

    # --- delivered power (C9)
    # the matched-effectiveness cooling identification is the primary one
    # (better conditioned, r2 0.986 vs 0.964); the closure bracket now spans
    # it and the heating tangent conductance (2026-09-04)
    K_lo, K_hi = ident["cooling_matched_eps"]["K_loss"], ident["heating_deep"]["K_loss"]
    pw = {}
    for Io, g in dg.groupby("Io_kWm2"):
        f_lo, f_hi = closure(g, K_lo).mean(), closure(g, K_hi).mean()
        pw[int(Io)] = dict(f_Klo=float(f_lo), f_Khi=float(f_hi),
                           G_closure_lo=float(Io * f_lo), G_closure_hi=float(Io * f_hi),
                           G_nominal=float(Io),
                           G_interval=[float(min(Io, Io * f_lo, Io * f_hi)),
                                       float(max(Io, Io * f_lo, Io * f_hi))],
                           eta_nom=[float(g.eta_nom.min()), float(g.eta_nom.max())])
    rep["delivered_power"] = dict(K_bracket=[float(K_lo), float(K_hi)], per_flux=pw,
                                  reconciling_dT3=reconciling_dT3(dg, K_hi))

    # --- T3 sensitivity band (C7)
    t3 = {}
    for dT3 in (-DT3_BAND, 0.0, DT3_BAND):
        d2 = dimensionless(ss, dT3=dT3)
        bb, aa = np.polyfit(np.log(d2.Re), np.log(d2.Nu), 1)
        cr = {int(Io): crossings(g)["eps_local"] for Io, g in d2.groupby("Io_kWm2")}
        t3[f"{dT3:+.0f}"] = dict(eps=[float(d2.eps.min()), float(d2.eps.max())],
                                 Nu_prefactor=float(np.exp(aa)), Nu_exponent=float(bb),
                                 NTU_exponent=float(linregress(np.log(d2.Re),
                                                               np.log(d2.NTU)).slope),
                                 eps_star=cr,
                                 eta_nom=[float(d2.eta_nom.min()), float(d2.eta_nom.max())])
        # C7b (2026-09-04, review v9 item 2): the identified constants derive
        # from T3 through the abscissa x = eps*mdot*cp, so the declared band
        # must be carried through the identification as well, not only through
        # the steady-state groups. eps is that of the provenance heating run,
        # recomputed at this offset; c_p uses the cooling run's own tails,
        # shifted by the same offset.
        xm, xd = [], []
        for _, cr_ in eig[eig.phase == "cool"].iterrows():
            hh = d2[d2.ID == COOL_PROVENANCE[cr_.ID]].iloc[0]
            xm.append(float(hh.eps) * (RHO_STD * float(cr_.q) / 60000.0)
                      * float(cp_air(0.5 * (cr_.Tamb_e + cr_.T3_e + dT3))))
        Cm3, Km3, r2m3, _ = identify(np.array(xm),
                                     eig[eig.phase == "cool"].lam.values)
        # the heating-branch abscissa uses the POOLED eps(q) fit, as the stored
        # heating_deep identification does, refitted at this offset so the
        # dT3 = 0 case reproduces the archived point estimate exactly
        eq3 = np.polyfit(d2.q_slpm, d2.eps, 1)
        eh = eig[eig.phase == "heat"].dropna(subset=["lam"])
        for _, hr in eh.iterrows():
            xd.append(float(np.polyval(eq3, hr.q)) * (RHO_STD * float(hr.q) / 60000.0)
                      * float(cp_air(0.5 * (hr.Tamb_e + hr.T3_e + dT3))))
        Cd3, Kd3, r2d3, _ = identify(np.array(xd), eh.lam.values)
        t3[f"{dT3:+.0f}"]["identification"] = dict(
            C_match=float(Cm3), K_match=float(Km3), r2_match=float(r2m3),
            C_deep=float(Cd3), K_deep=float(Kd3), r2_deep=float(r2d3))
    _offsets = list(t3.keys())
    t3["band"] = {}
    for q_ in ("C_match", "K_match", "C_deep", "K_deep"):
        vals = [t3[k]["identification"][q_] for k in _offsets]
        t3["band"][q_] = [float(min(vals)), float(max(vals))]
    rep["T3_sensitivity"] = t3

    # --- reference-probe sensitivity (C8)
    rv = reference_variants(ss)
    rv.to_csv(P("reference_sensitivity.csv"), index=False)
    rep["reference_sensitivity"] = rv.to_dict("records")

    # --- pressure drop (C10)
    dp = dg[["ID", "Io_kWm2", "q_slpm", "mdot_gs", "Tg_bar", "dp1_mbar", "dp2_mbar"]].copy()
    dp["dp_pred_mbar"] = [dp_laminar(m, T) for m, T in zip(dp.mdot_gs, dp.Tg_bar)]
    dp["ratio"] = dp.dp1_mbar / dp.dp_pred_mbar
    dp.to_csv(P("pressure_drop.csv"), index=False)
    rep["pressure_drop"] = dict(resolution_mbar=DP_FS * DP_ACC,
                                pred_range=[float(dp.dp_pred_mbar.min()),
                                            float(dp.dp_pred_mbar.max())],
                                meas_range=[float(dp.dp1_mbar.min()),
                                            float(dp.dp1_mbar.max())],
                                ratio_range=[float(dp.ratio.min()), float(dp.ratio.max())],
                                dp2_range=[float(dp.dp2_mbar.min()), float(dp.dp2_mbar.max())])

    # --- similarity collapse
    # rescaled on the PRIMARY (matched-effectiveness) identification and on
    # each heating run's own measured effectiveness, not the pooled eps(q) fit
    C_ref = ident["cooling_matched_eps"]["C_eff"]
    K_ref = ident["cooling_matched_eps"]["K_loss"]
    half = {"wall": [], "gas": []}
    for ID, (fn, Io) in HEATING.items():
        d = load(raw_dir, fn)
        t = d["t"]
        q = tail_mean(d["flow"], t)
        mdot = RHO_STD * q / 60000.0
        Tamb = tail_mean(d["Tamb"], t)
        eps = float(dg[dg.ID == ID].iloc[0].eps)
        tau = C_ref / (eps * mdot * cp_air(0.5 * (Tamb + tail_mean(d["T3"], t))) + K_ref)
        for key, sig in (("wall", sum(WTS[k] * d[k] for k in WTS)), ("gas", d["T3"])):
            th = (sig - sig[0]) / (tail_mean(sig, t) - sig[0])
            i = int(np.argmax(th >= 0.5))
            half[key].append(t[i] / tau)
    rep["similarity"] = {k: dict(t_half_mean=float(np.mean(v)),
                                 cv_percent=float(100 * np.std(v) / np.mean(v)))
                         for k, v in half.items()}

    # --- Monte Carlo
    # D1: report both controller-correlation limits. The fully correlated case
    # (rho=1) is the headline, being the conservative one for the systematic
    # part; the independent case (rho=0) bounds what reaches an exponent.
    _mcs = {}
    for _rho in RHO_MFC_CASES:
        _m = monte_carlo(ss, eig[eig.phase.isin(["cool", "heat"])], n=n_mc,
                         rho=_rho)
        _m["rho_mfc"] = _rho
        _mcs[_rho] = _m
    mc = _mcs[1.0]
    pd.concat([_mcs[r] for r in RHO_MFC_CASES]).to_csv(P("uncertainty.csv"),
                                                       index=False)
    rep["controller_correlation_bounds"] = {
        f"rho_{r:g}": {q: dict(value=float(v["value"]), sd=float(v["sd"]))
                       for q, v in _mcs[r].set_index("quantity").iterrows()}
        for r in RHO_MFC_CASES}
    rep["fixed_profile_test"] = fixed_profile_test(dg)
    rep["wall_extrapolation"] = wall_extrapolation_sensitivity(dg)
    rep["monte_carlo"] = mc.to_dict("records")
    write_tables(rep, dg, mc, P)
    write_supplementary_tables(rep, dg, P)

    with open(P("results.json"), "w") as fh:
        json.dump(rep, fh, indent=2, default=float)
    return rep, dg, dg_rep, eig, mc, rv, dp


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default="/mnt/d/kkakosim/github/tamuq-chen-secarelab-receiver/"
                                     "aysha/analysis/RAW")
    ap.add_argument("--out", default=".")
    ap.add_argument("--nmc", type=int, default=4000)
    A = ap.parse_args()
    R = main(A.raw, A.out, A.nmc)[0]
    print(json.dumps({k: R[k] for k in
                      ["geometry", "nusselt", "ntu_structure", "identification",
                       "sensor_selection_swing", "heating_window_sensitivity",
                       "crossings", "ltne", "delivered_power", "T3_sensitivity",
                       "pressure_drop", "similarity", "envelope"]},
                     indent=2, default=float))

