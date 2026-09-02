"""
Figure 8 - continuum model comparison (1D_v36) against the measured
effectiveness envelope and the volumetric inversion threshold.

Uses only saved artefacts; runs no simulation.

Inputs
  summaries/1D_v36/steady_results_fitted_energy_accounting_1D_v36.csv
  analysis/exp_analysis/dimensionless_groups.csv   (Tamb per run)

Manuscript definitions (M2):
  Tw_bar = 0.2518 T8 + 0.3504 T12 + 0.3978 T11
  eps    = (T3 - Tamb) / (Tw_bar - Tamb)
  I_vol  = T12 - T8
Model eps uses the identical reduction applied to the model output.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MODEL = os.path.join(ROOT, "summaries", "1D_v36",
                     "steady_results_fitted_energy_accounting_1D_v36.csv")
DIMLESS = os.path.join(HERE, "dimensionless_groups.csv")
OUT = os.path.join(HERE, "fig8_model_comparison.png")

W = dict(T8=0.2518, T12=0.3504, T11=0.3978)   # matches exp_analysis.WTS


def wall_mean(t8, t12, t11):
    return W["T8"] * t8 + W["T12"] * t12 + W["T11"] * t11


def main():
    m = pd.read_csv(MODEL)
    d = pd.read_csv(DIMLESS)[["ID", "Tamb_K", "Re", "eps", "I_vol_wall"]]
    df = m.merge(d, left_on="simulation_id", right_on="ID", how="inner")

    for tag in ("model", "experiment"):
        tw = wall_mean(df[f"T8_{tag}"], df[f"T12_perim_{tag}"],
                       df[f"T11_perim_{tag}"])
        df[f"Tw_{tag}"] = tw
        df[f"eps_{tag}"] = (df[f"T3_{tag}"] - df["Tamb_K"]) / (tw - df["Tamb_K"])
        df[f"NTU_{tag}"] = -np.log(1.0 - df[f"eps_{tag}"])
        df[f"Ivol_{tag}"] = df[f"T12_perim_{tag}"] - df[f"T8_{tag}"]

    df["flux"] = (df["irradiance"] / 1000.0).round().astype(int)
    # label groups by the delivered aperture irradiance G_del (Section 3.5)
    G_DEL = {456: 523, 304: 408, 256: 238}
    df["gdel"] = df["flux"].map(G_DEL)

    # consistency check against the published reduction
    resid = (df["eps_experiment"] - df["eps"]).abs()
    print(f"eps(recomputed) vs eps(published): max |diff| = {resid.max():.4f}")
    print(f"measured eps range   {df.eps_experiment.min():.3f}-{df.eps_experiment.max():.3f}"
          f"   NTU {df.NTU_experiment.min():.2f}-{df.NTU_experiment.max():.2f}")
    print(f"model    eps range   {df.eps_model.min():.3f}-{df.eps_model.max():.3f}"
          f"   NTU {df.NTU_model.min():.2f}-{df.NTU_model.max():.2f}")

    colors = {456: "#B3202E", 304: "#1F6FB4", 256: "#2E8B57"}
    marks = {456: "o", 304: "s", 256: "^"}
    GL = {456: 523, 304: 408, 256: 238}

    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.6))

    # ---- (a) effectiveness parity -------------------------------------
    a = ax[0]
    lo, hi = 0.52, 0.82
    a.plot([lo, hi], [lo, hi], "-", color="0.35", lw=1.0, zorder=1)
    a.fill_between([lo, hi], [lo - 0.03, hi - 0.03], [lo + 0.03, hi + 0.03],
                   color="0.85", alpha=0.55, zorder=0, label="±0.03")
    for f, g in df.groupby("flux"):
        a.scatter(g.eps_experiment, g.eps_model, s=46, marker=marks[f],
                  facecolor=colors[f], edgecolor="k", linewidth=0.5,
                  zorder=3, label=f"{GL[f]} kW m$^{{-2}}$")
    a.set_xlim(lo, hi); a.set_ylim(lo, hi)
    a.set_xlabel(r"measured effectiveness $\varepsilon$")
    a.set_ylabel(r"modelled effectiveness $\varepsilon$")
    a.set_title("(a) effectiveness envelope reproduced", loc="left", fontsize=10)
    a.legend(frameon=False, fontsize=8, loc="upper left")
    a.grid(alpha=0.25, lw=0.5)

    # ---- (b) inversion threshold --------------------------------------
    b = ax[1]
    b.axhline(0.0, color="0.35", lw=1.0, zorder=1)
    b.axvspan(0.63, 0.69, color="0.85", alpha=0.6, zorder=0)
    for f, g in df.groupby("flux"):
        g = g.sort_values("eps_experiment")
        b.plot(g.eps_experiment, g.Ivol_experiment, marks[f] + "-",
               color=colors[f], ms=6, lw=1.2, mec="k", mew=0.4, zorder=3,
               label=f"{GL[f]} measured")
        g = g.sort_values("eps_model")
        b.plot(g.eps_model, g.Ivol_model, marks[f] + "--",
               color=colors[f], ms=6, lw=1.2, mfc="white", mec=colors[f],
               zorder=2, label=f"{GL[f]} model")
    b.set_xlabel(r"effectiveness $\varepsilon$")
    b.set_ylabel(r"inversion indicator $I_{\rm vol}=T_{12}-T_{8}$  [K]")
    b.set_title(r"(b) inversion displaced in $\varepsilon$ and an order of "
                r"magnitude too weak", loc="left", fontsize=10)
    b.legend(frameon=False, fontsize=7, ncol=2, loc="upper left")
    b.grid(alpha=0.25, lw=0.5)
    b.text(0.66, b.get_ylim()[0] * 0.92, r"measured $\varepsilon^*$",
           ha="center", fontsize=8, color="0.3")

    fig.tight_layout()
    fig.savefig(OUT, dpi=220)
    print("wrote", OUT)

    # ---- inversion threshold, manuscript convention (Section 4.4) ------
    # Regress I_vol on flow to locate q*, then evaluate the eps regression
    # at q*. Applied identically to measurement and model so the two are
    # comparable; see make_tables.py, which uses the same reduction.
    def eps_star(q, ivol, eps):
        s_, i_ = np.polyfit(q, ivol, 1)
        qs = -i_ / s_
        se, ie = np.polyfit(q, eps, 1)
        return qs, se * qs + ie

    rows = []
    print("\ninversion threshold (regression convention, Section 4.4)")
    for f, g in df.groupby("flux"):
        q = g.flow_lpm.values
        qe, ee = eps_star(q, g.Ivol_experiment.values, g.eps_experiment.values)
        qm, em = eps_star(q, g.Ivol_model.values, g.eps_model.values)
        rows.append(dict(flux_kW_m2=f, q_star_meas=qe, eps_star_meas=ee,
                         q_star_model=qm, eps_star_model=em,
                         d_eps_star=em - ee,
                         Ivol_max_meas=g.Ivol_experiment.max(),
                         Ivol_max_model=g.Ivol_model.max(),
                         q_range_max=q.max()))
        print(f"  G_del={GL[f]} kW/m2 measured q*={qe:5.2f} sL/min eps*={ee:.3f}"
              f"   model q*={qm:5.2f} eps*={em:.3f}"
              f"   deps*={em - ee:+.3f}"
              f"   max I_vol {g.Ivol_model.max():+6.1f} K model /"
              f" {g.Ivol_experiment.max():+6.1f} K measured")
    cr = pd.DataFrame(rows)
    cr.to_csv(os.path.join(HERE, "fig8_inversion_thresholds.csv"), index=False)
    print("  model q* exceeds the tested flow range in "
          f"{int((cr.q_star_model > cr.q_range_max).sum())} of {len(cr)} groups")

    df.to_csv(os.path.join(HERE, "fig8_model_comparison_data.csv"), index=False)


if __name__ == "__main__":
    main()
