"""Figures 2-6 for the SiC volumetric receiver manuscript.

Reads the outputs of receiver_reduction.py plus three trace files it exports,
and writes figures/fig2..fig6. Run receiver_reduction.py first.

Usage:  python make_figures.py [--raw DIR] [--out DIR] [--fig FIGDIR]
"""
from __future__ import annotations
import argparse, json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from scipy.stats import linregress

CS = {456: "#1f4e79", 304: "#c1666b", 256: "#8f9779"}
FL = [456, 304, 256]
CC = {"C69": "#1f4e79", "C80": "#c1666b", "C81": "#8f9779"}
QC = {"C69": 10.5, "C80": 6.6, "C81": 4.5}
REFC = {"front wall only": "#7a3b2e", "mid wall only": "#c1666b",
        "wall quadrature": "#1f4e79", "interior probes": "#8f9779",
        "rear wall only": "#3f5c3a"}
REF_ORDER = ["front wall only", "mid wall only", "wall quadrature",
             "interior probes", "rear wall only"]


def style():
    mpl.rcParams.update({
        "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "figure.dpi": 110, "savefig.dpi": 300, "lines.markeredgewidth": 0,
        "font.family": "sans-serif", "mathtext.default": "regular",
    })


def panel_letter(ax, L):
    ax.text(-0.16, 1.06, L, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="bottom", ha="left")


def export_traces(raw_dir, out_dir):
    """Cooling decays, master curves and per-reference Nu points."""
    import sys
    _here = os.path.dirname(os.path.abspath(__file__))
    if _here not in sys.path:
        sys.path.append(_here)
    import receiver_reduction as rr
    ss = rr.reduce_steady(raw_dir, rr.HEATING)
    dg = rr.dimensionless(ss)
    eps_q = np.polyfit(dg.q_slpm, dg.eps, 1)
    P = lambda f: os.path.join(out_dir, f)

    rows = []
    for ID, fn in rr.COOLING.items():
        d = rr.load(raw_dir, fn); t = d["t"]
        Tamb = rr.tail_mean(d["Tamb"], t); k = max(1, len(t)//400)
        for s in rr.COOL_SENS:
            th = (d[s]-Tamb)/(d[s][0]-Tamb)
            rows.append(pd.DataFrame(dict(ID=ID, sensor=s, t=t[::k], theta=th[::k])))
    pd.concat(rows).to_csv(P("cooling_decays.csv"), index=False)

    C, K = 298.8, 0.0953
    rows = []
    for ID, (fn, Io) in rr.HEATING.items():
        d = rr.load(raw_dir, fn); t = d["t"]
        q = rr.tail_mean(d["flow"], t); Tamb = rr.tail_mean(d["Tamb"], t)
        mdot = rr.RHO_STD*q/60000.0
        Tg = 0.5*(Tamb + rr.tail_mean(d["T3"], t))
        tau = C/(float(np.polyval(eps_q, q))*mdot*rr.cp_air(Tg) + K)
        k = max(1, len(t)//400)
        for key, sig in (("wall", sum(rr.WTS[j]*d[j] for j in rr.WTS)), ("gas", d["T3"])):
            th = (sig-sig[0])/(rr.tail_mean(sig, t)-sig[0])
            rows.append(pd.DataFrame(dict(ID=ID, Io=Io/1e3, signal=key,
                                          tstar=t[::k]/tau, theta=th[::k])))
    pd.concat(rows).to_csv(P("master_curves.csv"), index=False)

    variants = {"wall quadrature": lambda o: sum(rr.WTS[k2]*o[f"{k2}_ss"] for k2 in rr.WTS),
                "interior probes": lambda o: 0.5*(o.T9_ss+o.T10_ss),
                "front wall only": lambda o: o.T8_ss,
                "mid wall only": lambda o: o.T12_ss,
                "rear wall only": lambda o: o.T11_ss}
    out = []
    for name, fn in variants.items():
        o = ss.copy(); o["Tw_K"] = fn(o); d2 = rr.dimensionless(o)
        out.append(pd.DataFrame(dict(reference=name, Re=d2.Re, Nu=d2.Nu,
                                     eps=d2.eps, NTU=d2.NTU, Io=d2.Io_kWm2)))
    pd.concat(out).to_csv(P("reference_points.csv"), index=False)


def main(raw_dir, out_dir, fig_dir):
    style()
    os.makedirs(fig_dir, exist_ok=True)
    P = lambda f: os.path.join(out_dir, f)
    F = lambda f: os.path.join(fig_dir, f)
    if not os.path.exists(P("cooling_decays.csv")):
        export_traces(raw_dir, out_dir)

    R = json.load(open(P("results.json")))
    dg = pd.read_csv(P("groups.csv"))
    eig = pd.read_csv(P("eigenvalues.csv"))
    rv = pd.read_csv(P("reference_sensitivity.csv"))
    cd_ = pd.read_csv(P("cooling_decays.csv"))
    mc_ = pd.read_csv(P("master_curves.csv"))
    rp = pd.read_csv(P("reference_points.csv"))
    sub = lambda f: dg[np.isclose(dg.Io_kWm2, f)].sort_values("q_slpm")
    rr_ = np.linspace(21, 102, 60)

    # ---- Figure 2: steady field vs flow
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.9), sharey=True)
    for ax, f in zip(axs, FL):
        g = sub(f)
        for col, mk in (("T8_ss", "o"), ("T12_ss", "s"), ("T11_ss", "^"), ("T3_ss", "D")):
            ax.plot(g.q_slpm, g[col]-273.15, mk+"-", color=CS[f], ms=3.8, lw=1.1,
                    mfc="white" if col == "T3_ss" else CS[f])
        ax.set_title(f"{f} kW m$^{{-2}}$", loc="left")
        ax.set_xlabel("$q$ [sL min$^{-1}$]")
    axs[0].set_ylabel("Steady temperature [$^\\circ$C]")
    h = [plt.Line2D([], [], color="0.3", marker=m, ls="-", ms=3.8, mfc=c, lw=1.1)
         for m, c in (("o", "0.3"), ("s", "0.3"), ("^", "0.3"), ("D", "white"))]
    axs[0].legend(h, ["wall 11 mm", "wall 58 mm", "wall 107 mm", "gas outlet"],
                  frameon=False, loc="lower left", handlelength=1.5)
    for ax, L in zip(axs, "abc"):
        panel_letter(ax, L)
    fig.tight_layout(); fig.savefig(F("fig2_steady_field.png"), bbox_inches="tight")
    plt.close(fig)

    # ---- Figure 3: assembly-scale limitation
    a, b = R["nusselt"]["prefactor"], R["nusselt"]["exponent"]
    bn = R["ntu_structure"]["exponent"]
    cN = np.polyfit(np.log(dg.Re), np.log(dg.NTU), 1)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.4, 3.1))
    for f in FL:
        g = sub(f); a1.plot(g.Re, g.Nu, "o", color=CS[f], ms=4.5,
                            label=f"{f} kW m$^{{-2}}$")
    # three parallel grouped lines (common slope, per-configuration prefactor)
    # are the primary fit; the pooled single-prefactor line is shown dashed
    GR = R["nusselt"]["grouped"]
    for f in FL:
        g = sub(f); rf = np.linspace(g.Re.min()*0.9, g.Re.max()*1.1, 40)
        a1.plot(rf, GR["prefactors"][str(f)]*rf**GR["exponent"], "-",
                color=CS[f], lw=1.5, zorder=2)
    a1.plot(rr_, a*rr_**b, "--", color="0.35", lw=1.2, zorder=1)
    # Shah & London (1978), square duct: (T) 2.976, (H2) 3.091, (H1) 3.608
    a1.axhspan(2.976, 3.608, color="0.85", zorder=0)
    a1.axhline(3.091, ls="--", color="0.55", lw=1.4)
    a1.text(23, 3.75, "fully developed laminar square duct,\n$Nu_{H2}=3.09$ (band: $Nu_T=2.98$ to $Nu_{H1}=3.61$)",
            color="0.45", va="bottom", fontsize=7)
    a1.text(21.5, 1.55, f"grouped (primary, solid):\n$\\propto Re_{{\\rm nom}}^{{{GR['exponent']:.3f}}}$, "
            f"$r^2$={GR['r2']:.3f}\n"
            f"pooled (dashed):\n${a*1e4:.2f}\\times10^{{-4}}Re_{{\\rm nom}}^{{{b:.2f}}}$, $r^2$=0.971",
            color="0.15", ha="left", va="top", fontsize=7.5)
    a1.set_xscale("log"); a1.set_yscale("log"); a1.set_xlim(20, 110); a1.set_ylim(0.02, 6.5)
    a1.set_xticks([25, 50, 100]); a1.set_xticklabels(["25", "50", "100"])
    a1.xaxis.set_minor_formatter(NullFormatter())
    a1.set_xlabel("Reynolds number, $Re$"); a1.set_ylabel("Apparent Nusselt number, $Nu$")
    a1.set_title("Exchange is limited at the assembly scale", loc="left")
    a1.legend(loc="lower right", frameon=False, handlelength=1.0)
    # profile-corrected transfer units are the primary series (2026-09-03);
    # NTU_app from the isothermal-wall identity is shown open for comparison
    dg["NTU_corr"] = np.array(R["ntu_profile"]["NTU_corr"])
    bnc = R["ntu_profile"]["exponent"]
    anc = np.polyfit(np.log(dg.Re), np.log(dg.NTU_corr), 1)[1]
    for f in FL:
        g = sub(f)
        a2.plot(g.Re, g.NTU, "o", mfc="none", mec=CS[f], ms=4.5, mew=1.0)
        a2.plot(g.Re, g.NTU_corr, "o", color=CS[f], ms=4.5)
    a2.plot(rr_, np.exp(anc)*rr_**bnc, "-", color="0.25", lw=1.6)
    a2.plot(rr_, 1.31*(rr_/72.5)**-1.0, "--", color="0.55", lw=1.4)
    a2.set_xscale("log"); a2.set_yscale("log"); a2.set_xlim(20, 110); a2.set_ylim(0.28, 3.0)
    a2.set_xticks([25, 50, 100]); a2.set_xticklabels(["25", "50", "100"])
    a2.set_yticks([0.5, 1, 2]); a2.set_yticklabels(["0.5", "1", "2"])
    a2.xaxis.set_minor_formatter(NullFormatter()); a2.yaxis.set_minor_formatter(NullFormatter())
    a2.set_xlabel("Reynolds number, $Re$"); a2.set_ylabel("Transfer units, $NTU$")
    a2.set_title("Transfer units rise with flow", loc="left")
    a2.text(102, 2.35, f"measured, $\\propto Re^{{{bnc:+.2f}}}$", color="0.15", ha="right")
    a2.text(21.5, 0.335, "filled: wall-profile integration\nopen: $-\\ln(1-\\varepsilon)$",
            color="0.35", ha="left", va="bottom", fontsize=7)
    a2.text(102, 0.335, "conductance fixed in $z$,\n$\\propto Re^{-1}$", color="0.5", ha="right")
    for ax, L in ((a1, "a"), (a2, "b")):
        panel_letter(ax, L)
    fig.tight_layout(); fig.savefig(F("fig3_assembly_limitation.png"), bbox_inches="tight")
    plt.close(fig)

    # ---- Figure 4: inversion + nonequilibrium
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(7.4, 2.9))
    for f in FL:
        g = sub(f); c = R["crossings"][str(f)]
        a1.plot(g.q_slpm, g.I_vol, "o-", color=CS[f], ms=4, lw=1.1)
        a1.plot([c["q_local"]], [0], "^", color=CS[f], ms=7, zorder=5)
        a2.plot(g.eps, g.I_vol, "o", color=CS[f], ms=4.5, label=f"{f}")
        a2.plot([c["eps_local"]], [0], "^", color=CS[f], ms=7, zorder=5)
        a3.plot(g.Re, g.Lam107, "o-", color=CS[f], ms=4, lw=1.1)
        a3.plot(g.Re, g.Lam58, "s--", color=CS[f], ms=3.4, lw=1.0, mfc="white")
    for ax in (a1, a2):
        ax.axhline(0, color="0.6", lw=0.8, zorder=0)
    a1.set_xlabel("$q$ [sL min$^{-1}$]"); a1.set_ylabel("$T_{12}-T_8$ [K]")
    a1.set_title("Front-to-mid side-wall crossing", loc="left")
    a1.annotate("single negative\npoint only", xy=(5.08, 0), xytext=(6.4, -88), fontsize=7,
                color="0.3", arrowprops=dict(arrowstyle="-", lw=0.7, color="0.45"))
    lo = min(R["crossings"][str(f)]["eps_local"] for f in FL)
    hi = max(R["crossings"][str(f)]["eps_local"] for f in FL)
    a2.axvspan(lo, hi, color="0.85", zorder=0)
    a2.set_xlabel("Effectiveness, $\\varepsilon$")
    a2.set_title("Crossings collapse on $\\varepsilon$", loc="left")
    a2.text(0.5*(lo+hi), 86, "$\\varepsilon^*$", ha="center", color="0.3")
    a2.legend(title="kW m$^{-2}$", loc="lower right", frameon=False,
              handlelength=0.8, fontsize=7, title_fontsize=7)
    a3.set_xlabel("Reynolds number, $Re$")
    a3.set_ylabel("Apparent wall-to-interior deficit, $\\Lambda$")
    a3.set_title("Wall-to-interior deficit grows with $Re$", loc="left")
    a3.text(90, 0.108, "$\\Lambda_{107}$", ha="right", color="0.25")
    a3.text(90, 0.040, "$\\Lambda_{58}$", ha="right", color="0.45")
    for ax, L in ((a1, "a"), (a2, "b"), (a3, "c")):
        panel_letter(ax, L)
    fig.tight_layout(); fig.savefig(F("fig4_inversion_ltne.png"), bbox_inches="tight")
    plt.close(fig)

    # ---- Figure 5: transient identification
    fig, (b1, b2, b3) = plt.subplots(1, 3, figsize=(7.4, 2.9))
    for ID, g in cd_.groupby("ID"):
        for s, gs in g.groupby("sensor"):
            b1.semilogy(gs.t/60, np.clip(gs.theta, 1e-2, None), color=CC[ID], lw=0.7, alpha=0.75)
    b1.set_xlabel("Time [min]"); b1.set_ylabel("Normalised excess, $\\theta/\\theta_0$")
    b1.set_title("One shared slow mode", loc="left"); b1.set_ylim(1e-1, 1.4); b1.set_xlim(0, 105)
    b1.set_yticks([0.1, 0.2, 0.5, 1.0]); b1.set_yticklabels(["0.1", "0.2", "0.5", "1"])
    b1.yaxis.set_minor_formatter(NullFormatter())
    for ID in ("C69", "C80", "C81"):
        b1.plot([], [], color=CC[ID], lw=1.4, label=f"{QC[ID]} sL min$^{{-1}}$")
    b1.legend(loc="lower left", frameon=False, handlelength=1.2, fontsize=7)
    # the cooling series is plotted at the PRIMARY matched-effectiveness
    # abscissa (archived as x_matched); the pooled-eps abscissa is shown open
    for ph, mk, lab, col, key in (
            ("cool", "o", "cooling, matched $\\varepsilon$ (primary)", "#1f4e79",
             "cooling_matched_eps"),
            ("heat", "s", "heating, deep 3", "#c1666b", "heating_deep"),
            ("heat6", "^", "heating, all 6", "#8f9779", "heating_all6")):
        e = eig[eig.phase == ph].dropna(subset=["lam"])
        xcol = e.x_matched if ph == "cool" else e.x
        if ph == "cool":
            b2.plot(e.x, e.lam*1e3, mk, mfc="none", mec=col, ms=4.2, mew=1.0,
                    label="cooling, pooled $\\varepsilon$")
        b2.plot(xcol, e.lam*1e3, mk, color=col, ms=4.2, label=lab)
        xx = np.linspace(0.04, 0.32, 20)
        I = R["identification"][key]
        b2.plot(xx, (xx+I["K_loss"])/I["C_eff"]*1e3, "-", color=col, lw=1.2)
    b2.set_xlabel("$\\varepsilon\\,\\dot m c_p$ [W K$^{-1}$]")
    b2.set_ylabel("$\\lambda$ [$10^{-3}$ s$^{-1}$]")
    b2.set_title("Probe set sets $C_{\\rm eff}$", loc="left")
    b2.legend(loc="upper left", frameon=False, handlelength=1.2, fontsize=7)
    for sig, ls, col in (("wall", "-", "0.35"), ("gas", "--", "#c1666b")):
        for ID, g in mc_[mc_.signal == sig].groupby("ID"):
            b3.plot(g.tstar, g.theta, ls, color=col, lw=0.6, alpha=0.8)
    b3.set_xlim(0, 2.2); b3.set_ylim(0, 1.08)
    b3.set_xlabel("Rescaled time, $t^*$"); b3.set_ylabel("Normalised rise, $\\theta^*$")
    b3.set_title("Transients collapse", loc="left")
    b3.text(1.45, 0.45, "wall", color="0.3")
    b3.text(1.45, 0.26, "gas outlet", color="#c1666b")
    for ax, L in ((b1, "a"), (b2, "b"), (b3, "c")):
        panel_letter(ax, L)
    fig.tight_layout(); fig.savefig(F("fig5_transient_identification.png"), bbox_inches="tight")
    plt.close(fig)

    # ---- Figure 6: under-instrumentation
    fig, (c1, c2) = plt.subplots(1, 2, figsize=(7.4, 3.1))
    for name in REF_ORDER:
        g = rp[rp.reference == name]
        b_, a_ = np.polyfit(np.log(g.Re), np.log(g.Nu), 1)
        c1.plot(g.Re, g.Nu, "o", color=REFC[name], ms=3.2, alpha=0.7)
        c1.plot(rr_, np.exp(a_)*rr_**b_, "-", color=REFC[name], lw=1.5,
                label=f"{name}: $Re^{{{b_:.2f}}}$")
    c1.axhline(3.091, ls="--", color="0.55", lw=1.2)   # Nu_H2, Shah & London
    c1.text(21, 3.4, "laminar duct", color="0.45", va="bottom", fontsize=7)
    c1.set_xscale("log"); c1.set_yscale("log"); c1.set_xlim(20, 112); c1.set_ylim(0.012, 9.0)
    c1.set_xticks([25, 50, 100]); c1.set_xticklabels(["25", "50", "100"])
    c1.xaxis.set_minor_formatter(NullFormatter())
    c1.set_xlabel("Reynolds number, $Re$"); c1.set_ylabel("Apparent Nusselt number, $Nu$")
    c1.set_title("Same data, five reference probes", loc="left")
    c1.legend(loc="lower right", frameon=False, handlelength=1.2, fontsize=6.5)
    keys = ("cooling", "heating_deep", "heating_all6", "joint")
    xlab = ["cooling\n6 probes", "heating\ndeep 3", "heating\nall 6", "joint\n18 eigenvalues"]
    val = [R["identification"][k]["C_eff"] for k in keys]
    colv = ["#1f4e79", "#c1666b", "#8f9779", "0.35"]
    c2.vlines(range(4), 0, val, color=colv, lw=1.2)
    for i, (v, cv) in enumerate(zip(val, colv)):
        c2.plot([i], [v], "o", ms=8, color=cv, zorder=4)
        c2.text(i, v+16, f"{v:.0f}", ha="center", fontsize=7, color="0.2")
    w = R["heating_window_sensitivity"]
    wlo, whi = min(x["C_eff"] for x in w), max(x["C_eff"] for x in w)
    c2.plot([1, 1], [wlo, whi], lw=6, color="#c1666b", alpha=0.25,
            solid_capstyle="butt", zorder=1)
    c2.annotate(f"fit-window\nrange {wlo:.0f}\u2013{whi:.0f}", xy=(1.07, 0.5*(wlo+whi)),
                xytext=(1.38, 150), fontsize=6.5, color="#8a4a4e",
                arrowprops=dict(arrowstyle="-", lw=0.7, color="#c1666b"))
    cm = R["monolith"]
    c2.axhspan(cm["C_monolith_600K"], cm["C_monolith_900K"], color="0.8", zorder=0)
    c2.text(3.45, 58, "bare monolith", ha="right", fontsize=7, color="0.35")
    c2.set_xticks(range(4)); c2.set_xticklabels(xlab, fontsize=7)
    c2.set_xlim(-0.5, 3.5); c2.set_ylim(0, 345)
    c2.set_ylabel("Effective capacitance, $C_{\\rm eff}$ [J K$^{-1}$]")
    c2.set_title("Same runs, four defensible answers", loc="left")
    for ax, L in ((c1, "a"), (c2, "b")):
        panel_letter(ax, L)
    fig.tight_layout(); fig.savefig(F("fig6_under_instrumentation.png"), bbox_inches="tight")
    plt.close(fig)

    print("wrote fig2-fig6 to", fig_dir)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default="../RAW")
    ap.add_argument("--out", default="outputs")
    ap.add_argument("--fig", default="figures")
    A = ap.parse_args()
    main(A.raw, A.out, A.fig)
