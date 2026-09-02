"""Monte Carlo uncertainty propagation for the dimensionless receiver analysis.

Sources (all systematic per run, independent between sensors/instruments):
  TC class accuracy: +-2.2 K (type N)  -> normal, sigma = 2.2/2
  MFC: 4 controllers, +-1% FS of 5 slpm each -> sigma = 0.05/2 slpm each
  air properties mu, k: +-2%; cp: +-1%  (correlation/tabulation uncertainty)
  (the radiometric flux map no longer enters any reported quantity: the
   aperture irradiance is determined from the energy closure, Section 3.5)
  steady-state residual drift: sigma = 0.5 K on steady means
N = 4000 draws. Outputs per-run sigma for Re, eps, NTU, Nu, Lambda107, Q_gas,
and Monte Carlo CIs for the Nu-law (a, b), eps*, and C_eff / K_loss.
"""
import numpy as np, pandas as pd, sys, os
from scipy.stats import linregress
_HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,_HERE)
from exp_analysis import load, ss, heating, cooling, A_frt, L, cp_air, F_DEL, WTS as _WTS

rng=np.random.default_rng(1)
N=4000
OUT=_HERE+os.sep
RHO_STD=101325/(287.05*294.25)
w=1.5e-3; Dh=w; A_ch=w*w; P=4*w; n_ch=100
WTS=np.array(_WTS)   # wall quadrature from exp_analysis
Ttab=[300,400,500,600,700,800,900,1000,1100,1200]
MU=[1.846e-5,2.301e-5,2.701e-5,3.058e-5,3.388e-5,3.698e-5,3.981e-5,4.244e-5,4.490e-5,4.730e-5]
KA=[0.0263,0.0338,0.0407,0.0469,0.0524,0.0573,0.0620,0.0667,0.0715,0.0763]

# nominal steady values per run
runs={}
for ID,(fn,Io) in heating.items():
    if not np.isfinite(Io): continue
    d=load(fn); t=d["t"]
    runs[ID]=dict(Io=Io,q=ss(d["flow"],t),Ta=ss(d["Tamb"],t),
        **{k:ss(d[k],t) for k in ["T3","T8","T9","T10","T11","T12"]})

sTC=1.1; sq=0.025; sSS=0.5
per_run={ID:{"Re":[],"eps":[],"NTU":[],"Nu":[],"L107":[],"Qg":[],"eta_d":[]} for ID in runs}
fits={"a":[],"b":[],"eps456":[],"eps304":[],"eps256":[]}

for n in range(N):
    bTC={s:rng.normal(0,sTC)+rng.normal(0,sSS) for s in ["T3","T8","T9","T10","T11","T12","Ta"]}
    fmu=1+rng.normal(0,0.02); fk=1+rng.normal(0,0.02); fcp=1+rng.normal(0,0.01)
    dq=rng.normal(0,sq,4).sum()
    # Delivered-power factor uncertainty: the K_loss bracket dominates. The
    # applied f uses the tangent conductance (0.119 W/K); the secant value
    # from cooling (0.096) lowers f by about 9%. One-sided systematic sampled
    # uniformly across the bracket, plus 3% for the closure's flow residual.
    fdel=(1-0.09*rng.random())*(1+rng.normal(0,0.03))
    X={"Re":[],"Nu":[],"eps":[],"q":[],"Io":[],"Iv":[]}
    for ID,r in runs.items():
        T3=r["T3"]+bTC["T3"]; Ta=r["Ta"]+bTC["Ta"]
        T8=r["T8"]+bTC["T8"]; T12=r["T12"]+bTC["T12"]; T11=r["T11"]+bTC["T11"]
        T9=r["T9"]+bTC["T9"]; T10=r["T10"]+bTC["T10"]
        q=r["q"]+dq
        mdot=RHO_STD*q/60000; mch=mdot/n_ch
        Tg=0.5*(Ta+T3); Tw=float(WTS@[T8,T12,T11]); Tf=0.5*(Tw+Tg)
        mu=np.interp(Tg,Ttab,MU)*fmu; ka=np.interp(Tf,Ttab,KA)*fk
        cp=cp_air(Tg)*fcp
        Re=mch*Dh/(A_ch*mu)
        eps=(T3-Ta)/(Tw-Ta); NTU=-np.log(1-eps)
        Nu=NTU*mch*cp/(P*L)*Dh/ka
        L107=(T11-T10)/(T11-Ta)
        Tgv=np.linspace(Ta,T3,60)
        Qg=mdot*np.trapezoid(cp_air(Tgv)*fcp,Tgv)
        eta_d=Qg/(r["Io"]*F_DEL[r["Io"]]*fdel*A_frt)     # single irradiance basis
        pr=per_run[ID]
        for k,v in [("Re",Re),("eps",eps),("NTU",NTU),("Nu",Nu),("L107",L107),("Qg",Qg),("eta_d",eta_d)]:
            pr[k].append(v)
        X["Re"].append(Re); X["Nu"].append(Nu); X["eps"].append(eps)
        X["q"].append(q); X["Io"].append(r["Io"]); X["Iv"].append(T12-T8)
    # Nu-law fit for this draw
    rr=linregress(np.log(X["Re"]),np.log(X["Nu"]))
    fits["a"].append(np.exp(rr.intercept)); fits["b"].append(rr.slope)
    # eps* per flux group
    for Io,keyn in [(456e3,"eps456"),(304e3,"eps304"),(256e3,"eps256")]:
        m=[i for i,v in enumerate(X["Io"]) if v==Io]
        ri=linregress([X["q"][i] for i in m],[X["Iv"][i] for i in m])
        q0=-ri.intercept/ri.slope
        re_=linregress([X["q"][i] for i in m],[X["eps"][i] for i in m])
        fits[keyn].append(re_.intercept+re_.slope*q0)

tab=[]
for ID,r in per_run.items():
    row={"ID":ID}
    for k in ["Re","eps","NTU","Nu","L107","Qg","eta_d"]:
        a=np.array(r[k]); row[k+"_mean"]=a.mean(); row[k+"_sd"]=a.std()
    tab.append(row)
ur=pd.DataFrame(tab)
ur.to_csv(OUT+"uncertainty_per_run.csv",index=False)
pd.set_option("display.width",250)
print(ur.round(4).to_string(index=False))

print("\nMC parameter CIs (mean +- sd; 2.5/97.5 pct):")
for k in ["a","b","eps456","eps304","eps256"]:
    a=np.array(fits[k])
    print(f"  {k:8s}: {a.mean():.4g} +- {a.std():.2g}  [{np.percentile(a,2.5):.4g},{np.percentile(a,97.5):.4g}]")

# ---- C_eff / K_loss MC: perturb lambda (fit sd), eps model, mdotcp --------
ev=pd.read_csv(OUT+"eigenvalue_verification.csv")
res={}
for tag,sub in [("cool",ev[ev.phase=="cool"]),("heat",ev[ev.phase=="heat"]),("all",ev)]:
    Cs=[];Ks=[]
    for n in range(N):
        lam=rng.normal(sub.lam,sub.lam_sd/np.sqrt(sub.n))
        x=sub.x*(1+rng.normal(0,0.03))*(1+rng.normal(0,0.05))  # mdotcp 3%, eps model 5%
        r=linregress(x,lam)
        if r.slope>0: Cs.append(1/r.slope); Ks.append(r.intercept/r.slope)
    Cs=np.array(Cs);Ks=np.array(Ks)
    res[tag]=(Cs.mean(),Cs.std(),Ks.mean(),Ks.std())
    print(f"C_eff ({tag:4s}): {Cs.mean():.0f} +- {Cs.std():.0f} J/K   K_loss: {Ks.mean():.3f} +- {Ks.std():.3f} W/K")

# monolith capacitance from measured mass
m_meas=0.040
for T,cpS in [(300,680),(600,1050),(900,1170)]:
    print(f"C_monolith(40 g, cp_SiC({T} K)={cpS}) = {m_meas*cpS:.1f} J/K")

# ---- persist the identified constants ------------------------------------
# These are the values quoted in manuscript Table 2. They were previously
# printed only; make_tables.py reads this file so the table cannot drift.
rows=[]
label={"a":"Nu_app prefactor a  (Nu = a Re^b)",
       "b":"Nu_app exponent b",
       "eps456":"eps* at 523 kW/m2",
       "eps304":"eps* at 408 kW/m2",
       "eps256":"eps* at 238 kW/m2"}
for k in ["a","b","eps456","eps304","eps256"]:
    v=np.array(fits[k])
    rows.append(dict(quantity=label[k], value=v.mean(), sd=v.std(),
                     ci_lo=np.percentile(v,2.5), ci_hi=np.percentile(v,97.5),
                     n_draws=len(v), source="uncertainty_analysis.py"))
for tag,name in [("cool","cooling"),("heat","heating"),("all","cooling+heating")]:
    Cm,Cs,Km,Ks=res[tag]
    rows.append(dict(quantity=f"C_eff from {name} [J/K]", value=Cm, sd=Cs,
                     ci_lo=Cm-1.96*Cs, ci_hi=Cm+1.96*Cs, n_draws=N,
                     source="uncertainty_analysis.py"))
    rows.append(dict(quantity=f"K_loss from {name} [W/K]", value=Km, sd=Ks,
                     ci_lo=Km-1.96*Ks, ci_hi=Km+1.96*Ks, n_draws=N,
                     source="uncertainty_analysis.py"))
for T,cpS in [(600,1050),(900,1170)]:
    rows.append(dict(quantity=f"Monolith capacitance at {T} K [J/K]",
                     value=m_meas*cpS, sd=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                     n_draws=np.nan, source="measured mass 40 g x cp_SiC"))
pd.DataFrame(rows).to_csv(OUT+"identified_constants_mc.csv", index=False)
print("\nwrote identified_constants_mc.csv")
