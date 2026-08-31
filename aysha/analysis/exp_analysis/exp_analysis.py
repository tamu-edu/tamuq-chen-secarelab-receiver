"""Experimental data analysis for the SiC honeycomb solar receiver.

Column contract follows import_exp_1D_v2.jl (1-indexed Julia -> 0-indexed here):
time=1, MFC actual flows=6:9, T2=35, T3=36, T8=41, T9=42, T10=43, T11=44,
T12=45, T15=48, T16=49.  All controllers carry air; flow = sum of the four.
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
import os, json

_HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(_HERE, "..", "RAW") + os.sep
OUT = _HERE + os.sep
os.makedirs(OUT, exist_ok=True)

heating = {
 "E67":("Data_FPT0067_231125_161757",456e3),"E68":("Data_FPT0068_231126_115725",456e3),
 "E69":("Data_FPT0069_231126_140153",456e3),"E70":("Data_FPT0070_231127_090339",456e3),
 "E71":("Data_FPT0071_231128_102707",456e3),"E72":("Data_FPT0072_231129_104140",304e3),
 "E73":("Data_FPT0073_231129_132744",304e3),"E74":("Data_FPT0074_231130_123228",304e3),
 "E75":("Data_FPT0075_231201_162138",304e3),"E76":("Data_FPT0076_231203_120521",304e3),
 "E77":("Data_FPT0077_231203_161315",256e3),"E78":("Data_FPT0078_231204_132252",256e3),
 "E79":("Data_FPT0079_231204_172244",256e3),"E80":("Data_FPT0080_231205_095122",256e3),
 "E81":("Data_FPT0081_231205_135354",256e3),
 # not in the Julia import lists -- extra runs, irradiance unknown (flagged)
 "E82":("Data_FPT0082_231210_130825",np.nan),"E83":("Data_FPT0083_231211_122053",np.nan),
}
cooling = {
 "C69":"Data_FPT0069-Cooling_231126_153148",
 "C80":"Data_FPT0080-cooling_231205_112837",
 "C81":"Data_FPT0081-cooling_231205_153409",
}

# geometry (1D_v8b review: w_t=19 mm, 10x10 channels of 1.5 mm, L=137 mm)
A_frt = 0.019**2
# delivered-power calibration factors (ad hoc, per lamp configuration):
# fitted per-irradiance scales from the independent 1D calibration (v14),
# corroborated by the steady energy closure of this analysis.
# Delivered-power factors, model-free steady energy closure (Section 3.5):
#   f = [Q_gas + K_loss*(Tw_bar - Tamb)] / (G0*A_frt),  group means.
# The closure is a steady balance and needs the SECANT loss conductance
# Q_loss/theta at the hot steady state. That quantity is bracketed by the
# cooling eigenvalue (secant conductance at lower temperature, 0.096 W/K) and
# the heating eigenvalue (tangent conductance dQ_loss/dT at the steady state,
# 0.119 W/K). The upper end is applied and the lower carried as a one-sided
# systematic band of about -9%. Values from eigenvalues_and_power.py; that
# script prints them and this constant must be kept consistent with it.
F_DEL = {456e3: 1.147, 304e3: 1.345, 256e3: 0.932}
L = 0.137
z8, z9, z10 = 0.011, 0.058, 0.107

cp_air = lambda T: np.interp(T,[300,400,500,600,700,800,900,1000,1100,1200],
                             [1007,1014,1030,1051,1075,1099,1121,1141,1159,1175])
rho_air = lambda T: 101325.0/(287.05*T)

def load(fname):
    df = pd.read_csv(RAW+fname+".csv", header=0, low_memory=False, encoding="latin-1")
    d = df.apply(pd.to_numeric, errors="coerce")
    t = d.iloc[:,1].to_numpy(float)
    keep = ~np.isnan(t)
    d = d[keep]; t = t[keep]; t = t - t[0]
    g = lambda i: d.iloc[:,i].to_numpy(float) + 273.15
    out = dict(t=t,
        flow=d.iloc[:,6].to_numpy(float)+d.iloc[:,7].to_numpy(float)
            +d.iloc[:,8].to_numpy(float)+d.iloc[:,9].to_numpy(float),
        T2=g(35),T3=g(36),T8=g(41),T9=g(42),T10=g(43),T11=g(44),T12=g(45))
    out["Tamb"] = 0.5*(g(48)+g(49))
    return out

def ss(x, t, win=120.0):
    m = t >= t[-1]-win
    return float(np.mean(x[m]))

def t90(t, x):
    x0, xss = np.mean(x[:10]), ss(x,t)
    if abs(xss-x0) < 5: return np.nan
    thr = x0 + 0.9*(xss-x0)
    idx = np.where(x >= thr)[0] if xss>x0 else np.where(x <= thr)[0]
    return float(t[idx[0]]) if len(idx) else np.nan

rows = []
store = {}
for ID,(fn,Io) in heating.items():
    d = load(fn); store[ID]=d
    t = d["t"]
    q = ss(d["flow"],t); Tamb = ss(d["Tamb"],t)
    S = {k: ss(d[k],t) for k in ["T2","T3","T8","T9","T10","T11","T12"]}
    mdot = rho_air(Tamb)*q/60000.0
    Tm = 0.5*(S["T3"]+Tamb)
    Qgas = mdot*cp_air(Tm)*(S["T3"]-Tamb)
    Qin = Io*A_frt if np.isfinite(Io) else np.nan
    Qdel = Io*F_DEL[Io]*A_frt if np.isfinite(Io) else np.nan
    rows.append(dict(ID=ID, Io_kWm2=Io/1e3 if np.isfinite(Io) else np.nan,
        q_lpm=q, Tamb=Tamb, dur_s=t[-1], **{k+"_ss":v for k,v in S.items()},
        # corrected topology: side-wall axial chain T8(11mm)-T12(58mm)-T11(107mm)
        # T9/T10 are interior, flow-exposed -> biased toward gas temperature
        I_vol_wall=S["T12"]-S["T8"],
        m_grad_wall=(S["T11"]-S["T12"])/(z10-z9),
        gapC58=S["T9"]-S["T12"], gapC107=S["T10"]-S["T11"],
        gap_deepening=(S["T10"]-S["T11"])-(S["T9"]-S["T12"]),
        dT_rad_ins=S["T12"]-S["T2"],
        R_leak=(S["T3"]-Tamb)/(S["T12"]-Tamb),
        t90_T8=t90(t,d["T8"]), t90_T9=t90(t,d["T9"]),
        t90_T10=t90(t,d["T10"]), t90_T3=t90(t,d["T3"]),
        mdot_gs=mdot*1e3, Q_gas_W=Qgas, Q_in_W=Qin, Q_del_W=Qdel,
        eta_gas=Qgas/Qin if np.isfinite(Qin) else np.nan,
        eta_del=Qgas/Qdel if np.isfinite(Qin) else np.nan))
ss_df = pd.DataFrame(rows)
ss_df.to_csv(OUT+"steady_state_metrics.csv", index=False)

# ---- flow-slope regressions per irradiance group -------------------------
slope_rows=[]
for Io,gname in [(456,"456"),(304,"304"),(256,"256")]:
    g = ss_df[np.isclose(ss_df.Io_kWm2,Io)]
    for s in ["T8_ss","T12_ss","T11_ss","T9_ss","T10_ss","T3_ss",
              "I_vol_wall","gapC58","gapC107","eta_gas","eta_del","R_leak"]:
        r = linregress(g.q_lpm, g[s])
        slope_rows.append(dict(Io_kWm2=Io, sensor=s, slope=r.slope,
                               intercept=r.intercept, r2=r.rvalue**2))
pd.DataFrame(slope_rows).to_csv(OUT+"flow_slopes.csv", index=False)

# ---- cooling: exponential time-constant identification -------------------
def expdec(t, dT, tau, Tinf): return Tinf + dT*np.exp(-t/tau)
cool_rows=[]; cool_store={}
for ID,fn in cooling.items():
    d = load(fn); cool_store[ID]=d
    t=d["t"]; Tamb = ss(d["Tamb"],t); q = float(np.mean(d["flow"]))
    for s in ["T8","T12","T11","T9","T10","T3","T2"]:
        y=d[s]
        try:
            p,_ = curve_fit(expdec, t, y, p0=[y[0]-Tamb, 1000.0, Tamb],
                            maxfev=20000)
            # late-time tau (second half) to expose non-single-exponential lag
            m = t> t[-1]/2
            p2,_ = curve_fit(expdec, t[m], y[m], p0=p, maxfev=20000)
            cool_rows.append(dict(ID=ID, q_lpm=q, sensor=s, T0=float(y[0]),
                tau_full_s=p[1], Tinf_full=p[2], tau_late_s=p2[1]))
        except Exception as e:
            cool_rows.append(dict(ID=ID,q_lpm=q,sensor=s,T0=float(y[0]),
                tau_full_s=np.nan,Tinf_full=np.nan,tau_late_s=np.nan))
cool_df = pd.DataFrame(cool_rows)
cool_df.to_csv(OUT+"cooling_time_constants.csv", index=False)

# effective capacitance check from cooling: C_eff ~ tau * (mdot cp eps + K)
# report tau ratio gas vs solid (hysteresis signature)
# ---- figures -------------------------------------------------------------
grp_colors={456:"tab:red",304:"tab:orange",256:"tab:blue"}

# 1. transient overview, one run per irradiance
fig,axs=plt.subplots(1,3,figsize=(15,4),sharey=True)
for ax,ID in zip(axs,["E67","E74","E78"]):
    d=store[ID]
    for s,c in [("T8","tab:red"),("T9","tab:orange"),("T10","tab:green"),
                ("T3","tab:blue"),("T2","tab:gray")]:
        ax.plot(d["t"]/60,d[s]-273.15,c=c,label=s)
    ax.set_title(f"{ID}  (Io={heating[ID][1]/1e3:.0f} kW/m2)")
    ax.set_xlabel("time [min]"); ax.grid(alpha=.3)
axs[0].set_ylabel("T [degC]"); axs[0].legend(fontsize=8)
fig.tight_layout(); fig.savefig(OUT+"fig1_transients.png",dpi=140); plt.close(fig)

# 2. steady state vs flow per irradiance (wall chain + gas out)
fig,axs=plt.subplots(1,4,figsize=(17,4))
for ax,s,ttl in zip(axs,["T8_ss","T12_ss","T11_ss","T3_ss"],
                    ["T8 wall (z=11mm)","T12 wall (z=58mm)",
                     "T11 wall (z=107mm)","T3 gas out"]):
    for Io in [456,304,256]:
        g=ss_df[np.isclose(ss_df.Io_kWm2,Io)].sort_values("q_lpm")
        ax.plot(g.q_lpm,g[s]-273.15,"o-",c=grp_colors[Io],label=f"{Io} kW/m2")
    ax.set_title(ttl); ax.set_xlabel("flow [L/min]"); ax.grid(alpha=.3)
axs[0].set_ylabel("steady T [degC]"); axs[0].legend(fontsize=8)
fig.tight_layout(); fig.savefig(OUT+"fig2_ss_vs_flow.png",dpi=140); plt.close(fig)

# 3. wall volumetric inversion + interior-minus-wall gaps + efficiency
fig,axs=plt.subplots(1,3,figsize=(14,4))
for Io in [456,304,256]:
    g=ss_df[np.isclose(ss_df.Io_kWm2,Io)].sort_values("q_lpm")
    axs[0].plot(g.q_lpm,g.I_vol_wall,"o-",c=grp_colors[Io],label=f"{Io}")
    axs[1].plot(g.q_lpm,g.gapC58,"o-",c=grp_colors[Io])
    axs[1].plot(g.q_lpm,g.gapC107,"s--",c=grp_colors[Io])
    axs[2].plot(g.q_lpm,g.eta_del*100,"o-",c=grp_colors[Io])
    axs[2].plot(g.q_lpm,g.eta_gas*100,"o--",c=grp_colors[Io],mfc="none",alpha=.5)
axs[0].axhline(0,c="k",lw=.8); axs[0].set_title("wall I_vol = T12-T8 [K]")
axs[1].axhline(0,c="k",lw=.8)
axs[1].set_title("interior-wall: o T9-T12 (58mm), s T10-T11 (107mm)")
axs[2].axhline(100,c="k",lw=.8)
axs[2].set_title("gas eff. [%]: solid=delivered, hollow=nominal")
for ax in axs: ax.set_xlabel("flow [L/min]"); ax.grid(alpha=.3)
axs[0].legend(title="kW/m2",fontsize=8)
fig.tight_layout(); fig.savefig(OUT+"fig3_inversion_gaps_eta.png",dpi=140); plt.close(fig)

# 4. cooling decays, log scale vs single-exp fit
fig,axs=plt.subplots(1,3,figsize=(14,4),sharey=True)
for ax,ID in zip(axs,["C69","C80","C81"]):
    d=cool_store[ID]; t=d["t"]; Tamb=ss(d["Tamb"],t)
    for s,c in [("T8","tab:red"),("T9","tab:orange"),("T10","tab:green"),("T3","tab:blue")]:
        theta=(d[s]-Tamb)/(d[s][0]-Tamb)
        ax.semilogy(t/60,np.clip(theta,1e-3,None),c=c,label=s)
    ax.set_title(f"{ID} (flow={np.mean(d['flow']):.1f} lpm)")
    ax.set_xlabel("time [min]"); ax.grid(alpha=.3,which="both")
axs[0].set_ylabel("normalized excess temp"); axs[0].legend(fontsize=8)
fig.tight_layout(); fig.savefig(OUT+"fig4_cooling_lin.png",dpi=140); plt.close(fig)

# 5. t90 lags
fig,ax=plt.subplots(figsize=(7,4))
for Io in [456,304,256]:
    g=ss_df[np.isclose(ss_df.Io_kWm2,Io)].sort_values("q_lpm")
    ax.plot(g.q_lpm,g.t90_T9/60,"o-",c=grp_colors[Io],label=f"{Io} T9")
    ax.plot(g.q_lpm,g.t90_T3/60,"s--",c=grp_colors[Io],label=f"{Io} T3")
ax.set_xlabel("flow [L/min]"); ax.set_ylabel("t90 [min]"); ax.grid(alpha=.3)
ax.legend(fontsize=7,ncol=3); ax.set_title("heating t90: solid core vs gas outlet")
fig.tight_layout(); fig.savefig(OUT+"fig5_t90.png",dpi=140); plt.close(fig)

print(ss_df.round(1).to_string(index=False))
print(); print(pd.DataFrame(slope_rows).round(3).to_string(index=False))
print(); print(cool_df.round(1).to_string(index=False))
