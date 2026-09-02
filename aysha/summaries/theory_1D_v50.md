# 1D_v50: audited repair specification and verification note

**Date:** 2026-09-02  
**Scope:** implementation repair of `1D_v49`; not a new physical closure and not yet a calibrated model release.

## 1. Purpose and status

Version 50 repairs defects that made several advertised v49 features inactive or internally inconsistent. It preserves the v49 plant equations so that the effect of the software repair can be separated from later changes to the receiver physics.

The repaired implementation now:

1. makes parameters 24--28 active in the observation equations;
2. evaluates one unambiguous calibration objective containing the flow-slope penalty;
3. evaluates each heating case once per objective call and reuses those solutions for the signal and slope terms;
4. makes the numerical accessors agree with the optimizer bounds;
5. includes irradiance power scales 5--7 in the full fit;
6. initializes the exit-probe state from the first measured probe temperature instead of changing simulated outputs after integration; and
7. exposes `:plant`, `:observation`, and `:full` calibration stages.

The vector `pnew_v50` is deliberately the saved v49 vector. It is an initialization seed only. Because parameters 24--28 were inactive in the saved v49 run, their inherited numerical values are not v50 estimates.

## 2. Observation equations repaired in v50

### 2.1 Exit thermocouple, T3

The exit thermocouple is a passive, one-way observation state:

$$
C_{TC}\frac{dT_3}{dt}
=h_{TC}(\dot V)A_{TC}(T_{g,out}-T_3)
+w_{rad}\sigma A_{TC}(T_{tube,1}^{4}-T_3^{4})
+G_{stem}(T_{water}-T_3),
$$

where

$$
h_{TC}(\dot V)=h_{TC,15}\left(\frac{\dot V}{15\ \mathrm{L\,min^{-1}}}\right)^{0.60}
$$

for positive forced flow. At zero flow the model uses a small natural-convection coefficient and the tube temperature as the surrounding gas temperature. Parameters 24, 25, and 26 are respectively $h_{TC,15}$, $w_{rad}$, and $G_{stem}$.

This state does not remove energy from the gas or tube balances. That choice is appropriate only if the probe thermal load is negligible relative to the receiver ledger. It must be revisited if a calibrated probe conductance implies a non-negligible heat rate.

For every case, $T_3(t_0)$ is set to the first measured T3 value. The v49 post-solution operation that forced initial cooling outputs to measured values has been removed.

### 2.2 Rear solid observations, T10 and T11

The reported temperatures are now explicit mixtures:

$$
T_{10,obs}=(1-w_{10})T_{core}(z_{10})+w_{10}T_{rear}(z_{10}),
$$

$$
T_{11,obs}=(1-w_{11})T_{perim}(z_{11})+w_{11}T_{water}.
$$

Thus parameters 27 and 28 are active. These equations are empirical observation closures, not derived thermocouple fin models. Estimated weights at their bounds should therefore be treated as evidence that the closure or sensor-location model is inadequate, not as material properties.

## 3. Repaired objective

The v50 heating objective is

$$
\mathcal L = \mathcal L_{signals}+\mathcal L_{slopes}
+\mathcal L_C+\mathcal L_{power}.
$$

`L_signals` retains the v49 normalized transient/endpoint/shape and ordering terms. `L_C` constrains participating heat capacity near 301 J/K. `L_power` penalizes the inherited non-monotonic irradiance power-scale relation.

For each irradiance cluster (256, 304, and 456 kW/m2), an ordinary least-squares endpoint slope against measured flow is computed. The active slope set is T8, T11, T9, T10, T3, and T2. With six signals in three clusters, the normalized penalty is

$$
\mathcal L_{slopes}=\frac{1}{18}\sum_{c=1}^{3}\sum_{s=1}^{6}
\left[\frac{b_{model,c,s}-b_{exp,c,s}}{\sigma_s}\right]^2,
$$

using scales 5 K/LPM for T8, 1 K/LPM for T2, and 3 K/LPM for the other active signals. T12 is reported diagnostically but is not included in this term. This exclusion should be reconsidered only after the role and independence of the two perimeter measurements are established.

There is now exactly one `loss_cases_v50` method. The former duplicate definition and the broken dormant slope path are absent.

## 4. Search domain and fit stages

The active numerical limits for $G_{core-perim}$, flange scale, rear-tube conductance, and suction heat-transfer coefficient now equal their optimizer bounds: 0.1--100 W/(m K), 0.05--20, 0.01--15 W/K, and 0--150 W/(m2 K), respectively.

Fit stages are:

- `:observation`: parameters 24--28;
- `:plant`: parameters 1--2 and 4--23; and
- `:full`: parameters 1--2 and 4--28.

Parameter 3 remains fixed at Nu-infinity = 3.61. Unlike v49, the full and plant stages include power-scale parameters 5--7.

## 5. Verification performed on the inherited seed

The automated v50 smoke suite checks parameter dimensions, bound/accessor agreement at both endpoints, observation-parameter activity, objective composition, forward integration, instantaneous energy conservation, mesh sensitivity, and three cooling cases.

At the inherited seed, the repaired objective components are:

| Component | Value |
|---|---:|
| Signal loss | 0.102844 |
| Flow-slope loss | 8.219941 |
| Capacitance regularization | 0.000000689 |
| Power-scale regularization | 1.690000 |
| Total | 10.012786 |

The nonzero slope term demonstrates that it is active. Its large share of the total also demonstrates why the v49 objective near 1.713 cannot be compared with v50: v49 did not actually optimize this term.

Perturbation tests produce nonzero changes in the intended output for every observation parameter. Representative maximum signal changes for the chosen test perturbations were 84.1, 17.2, 54.3, 44.7, and 188.6 K for parameters 24 through 28, respectively.

The inherited-seed run also gives:

- participating heat capacity: 301.112 J/K;
- maximum instantaneous conservation residual: 5.68e-14 W;
- endpoint flow-slope mean absolute error over all reported signals: 7.357 K/LPM;
- core-gas heat change from 15 to 50 axial nodes: about 0.10%; and
- macroscopic HTC change from 15 to 50 nodes: about 11.2%.

The last value is not strong mesh convergence. A finer-grid study or discretization improvement remains necessary before using the macroscopic HTC as a reported coefficient.

## 6. Performance is not yet acceptable

Running v50 with calibration disabled is a software baseline, not a validation result. The inherited seed produces mean heating RMSE 76.53 K and mean heating endpoint absolute error 83.51 K. T3 dominates the degradation: its mean heating RMSE is 268.16 K and endpoint absolute error is 325.02 K. Mean cooling RMSE is 39.19 K.

The baseline T3 flow slopes are -0.064, -0.903, and -2.246 K/LPM at 256, 304, and 456 kW/m2, versus experimental slopes -3.039, -0.132, and +0.543 K/LPM. Activating a probe model therefore does not, by itself, explain the power-dependent change of sign.

An 80-evaluation observation-only screening run reduced the objective to 6.725 but stopped at the iteration limit. It drove $w_{rad}$ to its lower bound and both T10/T11 weights to their upper bounds. This is evidence of objective responsiveness and observation-model strain; it is not a calibrated parameter set and is not installed as `pnew_v50`.

A separate 160-evaluation full-stage screen reduced the total objective from 10.013 to 2.222 and also stopped at the iteration limit. Its component values were signal 0.189, slope 2.031, capacitance 0.00098, and power 0. The reduction confirms that the full active-index route works, including the irradiance power scales. It is not an acceptable release fit: signal loss increased relative to the inherited seed (0.103 to 0.189), the 456 and 304 kW/m2 power scales reached their lower bounds, and the T11 mixing weight reached its upper bound. The scalar objective is therefore exposing a conflict between absolute/transient fit and flow slopes that the present parameterization cannot yet resolve cleanly.

## 7. Scientific limitations retained from v49

Version 50 does not establish the missing physics behind the systematic flow-rate residuals. In particular, it retains:

- a single prescribed bulk stream with no independently resolved channel-to-channel flow distribution;
- the inherited optical source parameterization, including strong front localization;
- the inherited rear hardware and flange network;
- empirical, one-way sensor closures whose heat removal is absent from the plant ledger;
- a clamped LMTD diagnostic that can conceal an invalid temperature ordering;
- effective geometry and thermal-capacity assumptions whose identifiability is coupled to transport parameters; and
- calibration on the heating data used to form the flow slopes, with cooling as the principal held-out dynamic check.

Accordingly, a successful v50 fit would show that the repaired 1D parameterization can reproduce the data; it would not uniquely identify the cause of the flow bias. The earlier 2D evidence also does not justify assuming a low-flow maldistribution bifurcation or a narrow Gaussian beam without new measurements or a discriminating simulation study.

## 8. Recommended release sequence

1. Run a full v50 calibration from multiple starts and record convergence, parameter bounds, and objective components.
2. Reject fits that merely exchange plant error for observation weights at their bounds.
3. Re-run heating, all three cooling cases, mesh checks, and energy audits with the selected vector.
4. Compare per-power flow-slope residuals, not only aggregate RMSE.
5. Profile or bootstrap the key transport and observation parameters to assess identifiability.
6. If the sign-changing T3 trend and systematic rear residuals persist, test explicit hypotheses one at a time: flow distribution, rear leakage/contact conductance, optical deposition, then coupled probe heat loading.

A full run can be launched in PowerShell with:

```powershell
$env:RECEIVER1D_v50_RUN_CALIBRATION = "true"
$env:RECEIVER1D_v50_FIT_STAGE = "full"
& "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. run_1D_v50.jl
```

Until that workflow converges and the acceptance checks pass, v50 should be described as **repaired and testable**, not **calibrated**, **validated**, or **physically explanatory**.
