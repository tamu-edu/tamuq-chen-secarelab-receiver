"""Steps 2-4: dimensionless groups, steady correlations, transient collapse.

Corrected topology: wall chain T8(11)-T12(58)-T11(107); T9/T10 interior
flow-exposed. Flow metering: Aalborg GFC17, standard L/min (21.1 C, 1 atm).
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit
import os, sys
_HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,_HERE)
from exp_analysis import load, ss, heating, cooling, A_frt, L, z8, z9, z10, cp_air, F_DEL, G_DEL, GLAB, WTS as _WTS

OUT=_HERE+os.sep
RHO_STD=101325/(287.05*294.25)          # Aalborg std: 21.1 C, 1 atm
w_ch=1.5e-3; Dh=w_ch; n_ch=100; A_ch=w_ch**2; P_ch=4*w_ch
k_SiC=48.8; t_wall=0.4e-3
WTS=np.array(_WTS)   # wall quadrature, computed in exp_analysis from z8/z9/z10

Ttab=[300,400,500,600,700,800,900,1000,1100,1200]
mu_air=lambda T: np.interp(T,Ttab,[1.846e-5,2.301e-5,2.701e-5,3.058e-5,3.388e-5,3.698e-5,3.981e-5,4.244e-5,4.490e-5,4.730e-5])
k_air =lambda T: np.interp(T,Ttab,[0.0263,0.0338,0.0407,0.0469,0.0524,0.0573,0.0620,0.0667,0.0715,0.0763])

rows=[]; store={}
for ID,(fn,Io) in heating.items():
    d=load(fn); store[ID]=d; t=d["t"]
    q=ss(d["flow"],t); Tamb=ss(d["Tamb"],t)
    S={k:ss(d[k],t) for k in ["T2","T3","T8","T9","T10","T11","T12"]}
    mdot=RHO_STD*q/60000; mch=mdot/n_ch
    Tg_m=0.5*(Tamb+S["T3"])                       # mean gas temperature
    Tw=float(WTS@np.array([S["T8"],S["T12"],S["T11"]]))  # energy-wtd wall temp
    Tfilm=0.5*(Tw+Tg_m)
    Re=mch*Dh/(A_ch*mu_air(Tg_m))
    Pr=cp_air(Tg_m)*mu_air(Tg_m)/k_air(Tg_m)
    Gz_L=Dh*Re*Pr/L                                # Graetz at outlet
    eps=(S["T3"]-Tamb)/(Tw-Tamb)                   # effectiveness (R_leak)
    NTU=-np.log(1-eps) if eps<1 else np.nan
    h_g=NTU*mch*cp_air(Tg_m)/(P_ch*L)              # global HTC per channel
    Nu=h_g*Dh/k_air(Tfilm)
    Nu_H=3.66+0.0668*(Dh/L*Re*Pr)/(1+0.04*(Dh/L*Re*Pr)**(2/3))  # Hausen
    Bi=h_g*(t_wall/2)/k_SiC
    # LTNE (interior-wall depression normalized by wall excess)
    ltne58=(S["T12"]-S["T9"])/(S["T12"]-Tamb)
    ltne107=(S["T11"]-S["T10"])/(S["T11"]-Tamb)
    PRn=Io*A_frt/mdot/1e3 if np.isfinite(Io) else np.nan     # nominal kJ/kg
    PRd=PRn*F_DEL[Io] if np.isfinite(Io) else np.nan         # delivered kJ/kg
    rows.append(dict(ID=ID,Io_kWm2=Io/1e3 if np.isfinite(Io) else np.nan,
        q_slpm=q,mdot_gs=mdot*1e3,Tw_K=Tw,T3_K=S["T3"],Tamb_K=Tamb,
        Re=Re,Pr=Pr,Gz_L=Gz_L,eps=eps,NTU=NTU,h_Wm2K=h_g,Nu=Nu,Nu_Hausen=Nu_H,
        Nu_ratio=Nu/Nu_H,Bi=Bi,LTNE58=ltne58,LTNE107=ltne107,
        PR_nom_kJkg=PRn,PR_del_kJkg=PRd,I_vol_wall=S["T12"]-S["T8"]))
dg=pd.DataFrame(rows)
dg.to_csv(OUT+"dimensionless_groups.csv",index=False)

main=dg[np.isfinite(dg.Io_kWm2)]

# ---- Nu = a Re^b fit (Pr ~ const) ---------------------------------------
r=linregress(np.log(main.Re),np.log(main.Nu))
a,b=np.exp(r.intercept),r.slope
print(f"Nu = {a:.4f} Re^{b:.3f}   r2={r.rvalue**2:.3f}  b_stderr={r.stderr:.3f}")
print(f"Nu/Nu_Hausen: mean={main.Nu_ratio.mean():.3f}  range=({main.Nu_ratio.min():.3f},{main.Nu_ratio.max():.3f})")

# ---- inversion criterion: eps at I_vol_wall = 0 -------------------------
print("\ninversion criterion (eps* at wall-inversion crossing):")
for Io in [456,304,256]:
    g=main[np.isclose(main.Io_kWm2,Io)]
    ri=linregress(g.q_slpm,g.I_vol_wall); q0=-ri.intercept/ri.slope
    re_=linregress(g.q_slpm,g.eps); eps0=re_.intercept+re_.slope*q0
    rr=linregress(g.q_slpm,g.Re);  Re0=rr.intercept+rr.slope*q0
    PR0=(Io*1e3)*F_DEL[Io*1e3]*A_frt/(RHO_STD*q0/60000)/1e3  # exact at q0
    print(f"  G_del={G_DEL[Io*1e3]} kW/m2 (nominal {Io}): q0={q0:5.2f} slpm  eps*={eps0:.3f}  Re*={Re0:.1f}  PR_del*={PR0:.0f} kJ/kg")

# ---- LTNE vs Re ---------------------------------------------------------
r107=linregress(main.Re,main.LTNE107); r58=linregress(main.Re,main.LTNE58)
print(f"\nLTNE107 = {r107.intercept:.4f} + {r107.slope:.5f} Re  (r2={r107.rvalue**2:.3f})")
print(f"LTNE58  = {r58.intercept:.4f} + {r58.slope:.5f} Re  (r2={r58.rvalue**2:.3f})")

# ---- transient: C_eff, K_loss with stderr -------------------------------
taus={"C69":None,"C80":None,"C81":None}
cool_late={}
for ID,fn in cooling.items():
    d=load(fn); t=d["t"]; Tamb=ss(d["Tamb"],t)
    q=float(np.mean(d["flow"])); mdotcp=RHO_STD*q/60000*cp_air(0.5*(np.mean(d["T3"])+Tamb))
    lam=[]
    for s in ["T8","T12","T11","T9","T10","T3"]:
        y=d[s]; m=t>t[-1]/2
        th=(y[m]-Tamb)
        good=th>5
        if good.sum()>50:
            rl=linregress(t[m][good],np.log(th[good]))
            lam.append(-rl.slope)
    cool_late[ID]=dict(q=q,mdotcp=mdotcp,lam=np.mean(lam),lam_sd=np.std(lam))
cl=pd.DataFrame(cool_late).T
# The slow mode is lam = (eps*mdot*cp + K)/C, so the abscissa carries eps.
# eps for a cooling run is taken from the pooled steady correlation eps(q),
# the same convention as eigenvalues_and_power.py, so the two scripts return
# the same constants and the master-curve time scale below is the reciprocal
# of the identified eigenvalue rather than an unrelated scaling.
_re=linregress(main.q_slpm,main.eps)
cl["x"]=cl.mdotcp*(_re.intercept+_re.slope*cl.q)
rC=linregress(cl.x,cl.lam)
C_eff=1/rC.slope; K=rC.intercept*C_eff
dC=C_eff*rC.stderr/rC.slope
print(f"\nslow mode: lam = (eps*mdotcp+K)/C ; C_eff={C_eff:.0f}+-{abs(dC):.0f} J/K, K_loss={K:.3f} W/K, r2={rC.rvalue**2:.4f}")
print(cl.round(6).to_string())

# ---- axial apparent diffusivity from cooling tau ordering ---------------
print("\napparent axial diffusivity (pairwise, tau from early exp fits):")
def expdec(t,dT,tau,Ti): return Ti+dT*np.exp(-t/tau)
for ID,fn in cooling.items():
    d=load(fn); t=d["t"]
    tau={}
    for s in ["T8","T12","T11"]:
        p,_=curve_fit(expdec,t,d[s],p0=[d[s][0]-300,600,300],maxfev=20000)
        tau[s]=p[1]
    a1=(z9-z8)**2/(tau["T12"]-tau["T8"])
    a2=(z10-z9)**2/(tau["T11"]-tau["T12"])
    print(f"  {ID}: tau(T8,T12,T11)=({tau['T8']:.0f},{tau['T12']:.0f},{tau['T11']:.0f}) s"
          f"  alpha_app={a1:.2e}, {a2:.2e} m2/s")

# ---- master-curve collapse ----------------------------------------------
fig,axs=plt.subplots(1,2,figsize=(12,4.5))
cols=plt.cm.viridis(np.linspace(0,1,15))
spread={0:[],1:[]}
for i,(ID,(fn,Io)) in enumerate([x for x in heating.items() if np.isfinite(x[1][1])]):
    d=store[ID]; t=d["t"]
    q=ss(d["flow"],t); Tamb=ss(d["Tamb"],t)
    mdotcp=RHO_STD*q/60000*cp_air(0.5*(ss(d["T3"],t)+Tamb))
    eps_run=_re.intercept+_re.slope*q          # pooled steady correlation eps(q)
    tstar=t*(eps_run*mdotcp+K)/C_eff           # = t*lambda, the identified slow mode
    for j,sig in enumerate(["Tw","T3"]):
        y=(WTS[0]*d["T8"]+WTS[1]*d["T12"]+WTS[2]*d["T11"]) if sig=="Tw" else d["T3"]
        y0,yss=np.mean(y[:10]),ss(y,t)
        th=(y-y0)/(yss-y0)
        axs[j].plot(tstar,th,c=cols[i],lw=.8,label=ID if j==0 else None)
        # t* at theta=0.5
        k50=np.argmax(th>=0.5)
        if k50>0: spread[j].append(tstar[k50])
axs[0].set_title("energy-weighted wall temperature")
axs[1].set_title("gas outlet T3")
for ax in axs:
    ax.set_xlabel("t* = t(mcp+K)/C_eff"); ax.set_ylabel("theta")
    ax.set_xlim(0,6); ax.grid(alpha=.3)
axs[0].legend(fontsize=6,ncol=3)
for j,n in [(0,"wall"),(1,"T3")]:
    s=np.array(spread[j]); print(f"master-curve t*(theta=0.5) {n}: mean={s.mean():.2f} CV={s.std()/s.mean()*100:.0f}%")
fig.tight_layout(); fig.savefig(OUT+"fig6_master_curve.png",dpi=140); plt.close(fig)

# ---- Nu-Re figure -------------------------------------------------------
fig,axs=plt.subplots(1,3,figsize=(15,4.2))
gc={456:"tab:red",304:"tab:orange",256:"tab:blue"}
for Io in [456,304,256]:
    g=main[np.isclose(main.Io_kWm2,Io)]
    axs[0].loglog(g.Re,g.Nu,"o",c=gc[Io],label=GLAB(Io))
    axs[1].plot(g.Re,g.LTNE107,"s",c=gc[Io])
    axs[1].plot(g.Re,g.LTNE58,"o",c=gc[Io],mfc="none")
    axs[2].plot(g.q_slpm,g.eps,"o",c=gc[Io])
Rex=np.linspace(main.Re.min(),main.Re.max(),50)
axs[0].loglog(Rex,a*Rex**b,"k-",label=f"fit {a:.2e} Re$^{{{b:.2f}}}$")
axs[0].loglog(Rex,3.66+0.0668*(Dh/L*Rex*0.69)/(1+0.04*(Dh/L*Rex*0.69)**(2/3)),
              "k--",label="Hausen")
axs[0].set_xlabel("Re"); axs[0].set_ylabel("Nu"); axs[0].legend(fontsize=8)
axs[0].set_title("global channel Nusselt number")
axs[1].set_xlabel("Re"); axs[1].set_title("LTNE: o 58mm, s 107mm")
axs[2].set_xlabel("flow [slpm]"); axs[2].set_title("effectiveness eps=(T3-Ta)/(Tw-Ta)")
for ax in axs: ax.grid(alpha=.3,which="both")
fig.tight_layout(); fig.savefig(OUT+"fig7_Nu_LTNE_eps.png",dpi=140); plt.close(fig)

print("\n",dg.round(3).to_string(index=False))
