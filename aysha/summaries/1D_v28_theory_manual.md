# Theory Manual: 1D_v28 Energy-Accounting Two-Zone Receiver Model

## Manuscript-Level Model Description and Scientific Interpretation

This document describes the mathematical formulation, parameterization, and
interpretation of `1D_v28.jl`, a one-dimensional reduced-order thermal model
for a structured monolithic solar receiver with square channels. The model is
part of a sequence of increasingly constrained macro-scale receiver models
developed to identify and validate effective heat-transfer closures from
transient experimental data.

The specific purpose of v28 is not to introduce another empirical correction
for better fitting. Instead, v28 is a controlled structural test. It removes
the direct distributed rear-core heat sink used in v27, replaces the hard
lamp-on/lamp-off rear-loss switch with a smooth operating-state gate on the
explicit rear tube/flange pathway, and compares the T3 thermocouple to the gas
temperature at 140 mm without a fitted wall-temperature blend. The resulting
fit is intentionally evaluated as a scientific diagnostic, including whether
the model remains identifiable after removing the prior direct rear sink.

The main outcome is that v28 is a useful negative result: the explicit
rear-tube/flange pathway by itself does not reproduce the combined heating and
cooling data with the current two-zone topology. The fit collapses several
thermal-coupling and heat-capacity parameters to lower bounds, indicating that
the removed v27 direct rear sink was standing in for missing rear/adaptor
thermal inventory and loss topology. Therefore, v28 should be interpreted as a
model-form falsification test, not as a preferred predictive baseline.

---

## 1. Physical System and Modeling Objective

The experimental receiver is a silicon carbide monolith with a square frontal
aperture and a 10 x 10 array of square air channels. Concentrated solar flux is
incident on the front face. Air enters the channels, exchanges heat with the
solid matrix, and exits through a downstream tube/adaptor region. The receiver
is embedded in a cavity consisting of alumina insulation and an outer metal
housing, with an additional downstream water-cooled flange environment.

The broader scientific goal of the 1D model family is to obtain effective
macroscopic transport coefficients for:

- convective gas-solid heat exchange inside square monolith channels,
- radiative or volumetric source deposition through the receiver depth,
- conductive and effective axial redistribution in the solid receiver,
- rear/adaptor/flange heat removal and thermal storage.

The v28 model tests whether the downstream cooling and heating behavior can be
explained without a distributed direct receiver-to-flange sink. This is an
important distinction. A distributed sink can produce good fits, but it is hard
to assign scientific credit to a physically specific mechanism. Removing it
forces all downstream heat removal to proceed through explicitly represented
rear tube, cavity, and flange pathways.

---

## 2. State Variables and Domain Decomposition

The receiver is discretized axially into `N` finite volumes of length

```text
dx = L / N
```

where `L = 137 mm` is the monolith length. The model contains four dynamic
thermal domains:

1. Core receiver solid: `T_core,i(t)`, for `i = 1,...,N`.
2. Perimeter/housing receiver branch: `T_perim,i(t)`, for `i = 1,...,N`.
3. Rear alumina tube/adaptor branch: `T_tube,j(t)`, for `j = 1,...,N_rear`.
4. Lumped cavity/housing state: `T_cav(t)`.

The gas phase is not a dynamic state. It is reconstructed quasi-steadily at
each ODE evaluation by marching along the core channels and rear tube using an
exact NTU update. Thus the state vector is

```text
u(t) = [T_core,1 ... T_core,N,
        T_perim,1 ... T_perim,N,
        T_tube,1 ... T_tube,Nrear,
        T_cav]
```

The corresponding ODE system is written in energy form:

```text
C_k(T) dT_k/dt = sum of heat inputs - sum of heat outputs
```

and is integrated using a stiff ODE solver.

---

## 3. Geometric Definitions

The main receiver geometry is inherited from the earlier 1D model definitions:

| Quantity | Symbol | Value |
| :--- | :--- | ---: |
| Receiver frontal width | `w_t` | 19.0 mm |
| Frontal area | `A_frt = w_t^2` | 3.61e-4 m2 |
| Receiver length | `L` | 137.0 mm |
| Number of channels | `n_chnl` | 100 |
| Channel width | `w_chnl` | 1.5 mm |
| Hydraulic diameter | `D_h` | 1.5 mm |
| Total flow area | `A_flow` | `n_chnl w_chnl^2` |
| Gas-solid exchange perimeter | `P_exchange` | `4 w_chnl n_chnl` |

The cavity and downstream geometry in v28 includes:

| Quantity | Symbol | Value |
| :--- | :--- | ---: |
| Cavity outer radius | `R_cav` | 75.0 mm |
| Insulation outer radius | `R_ins` | 57.0 mm |
| Cavity length | `L_cav` | 165.0 mm |
| Adaptor diameter | `D_ad` | 77.6 mm |
| Adaptor length | `L_ad` | 57.0 mm |
| Rear tube length | `L_tube` | 150.0 mm |
| Rear tube gas radius | `r_tube,g` | 6.5 mm |
| Rear tube wall thickness | `t_tube` | 1.5 mm |
| Rear tube cavity-coupled length | `L_tube,cav` | 46.0 mm |

The T3 gas comparison point is fixed at

```text
z_T3 = 140 mm
```

which lies 3 mm downstream of the monolith exit.

---

## 4. Solar Power Accounting

### 4.1 Total Participating Power

For an experimental irradiance level `I`, v28 defines the modeled participating
receiver power as

```text
Q_abs(I) = eta_abs scale(I) I A_frt
```

where

```text
eta_abs = 1
```

and `scale(I)` is a fitted irradiance-level multiplier:

```text
scale(456 kW/m2) = p7
scale(304 kW/m2) = p8
scale(256 kW/m2) = p9
```

These scale factors should not be interpreted as absorptivity alone. They
combine local receiver-average flux, aperture-average flux convention,
calibration uncertainty, and any remaining model-form mismatch.

### 4.2 Core/Perimeter Power Partition

The transverse flux map is represented by a shape-only Gaussian:

```text
F(x,y) = exp[-(x^2 + y^2)/(2 sigma^2)]
```

with `sigma = 30 mm`. Integrating this profile over the circular aperture and
over the receiver square gives

```text
f_receiver = 0.0692
f_spill = 0.9308
```

The fitted spillover-capture fraction `p14` converts part of the off-receiver
Gaussian power into perimeter/housing source:

```text
f_perim = clamp(p14 f_spill, 0, 0.80)
```

Then

```text
Q_core = (1 - f_perim) Q_abs
Q_perim = f_perim Q_abs
```

For the fitted v28 parameter set,

```text
p14 = 0.3518
f_perim = 0.3275
```

### 4.3 Axial Source Deposition

Core power uses the inherited Beer-Lambert-plus-front-deposition weights:

```text
w_core,i proportional to exp(-beta z_left) - exp(-beta z_right)
```

with front deposition fraction `front_dep = 1`. Therefore all core source is
placed in the first axial cell in the current v28 parameterization. The optical
attenuation `beta_opt = 184.7 1/m` is fixed, but it has no practical influence
while `front_dep = 1`.

Perimeter source uses a separate exponential axial distribution:

```text
w_perim,i proportional to exp(-beta_perim z_left) - exp(-beta_perim z_right)
```

where `beta_perim = p15`.

This source structure is one of the most important limitations of v28. Since
the core source remains front-cell dominated, the model cannot independently
test deeper volumetric deposition or radiative redistribution.

---

## 5. Gas-Solid Convective Exchange

The gas is marched through the receiver channels using a quasi-steady plug-flow
NTU update. For cell `i`, the film temperature is

```text
T_film,i = 0.5 (T_core,i + T_g,in,i)
```

Air properties are evaluated at `T_film,i`. The mass flow rate is computed from
standard volumetric flow:

```text
mdot = rho_std Q_lpm / 60000
```

with

```text
rho_std = 101325 / (287.05 * 294.25)
```

The Reynolds number is

```text
Re = mdot D_h / (A_flow mu)
```

and the Prandtl number is

```text
Pr = cp mu / k_f
```

### 5.1 Active Flow Fraction

The model retains the active-flow fraction form:

```text
phi_act = clamp(phi_0 (Re/Re_ref)^m_rec, 0.1, 1.0)
```

In the fitted v28 default,

```text
phi_0 = 1
m_rec = 0
```

so all flow is treated as active in the receiver. This means v28 does not
include bypass or incomplete channel participation.

### 5.2 Nusselt Closure

The local receiver Nusselt number is

```text
Nu_i = max[Nu_floor,
           A_Nu Re_i^B_Re Pr_i^(1/3) (D_h / max(z_i, D_h/2))^(1/3)]
```

with

```text
Nu_floor = 3.61
0 <= B_Re <= 0.6
```

In v28, the upper bound on `A_Nu` was relaxed from 5 to 25 so the optimizer
could test whether the previous upper-bound pressure was artificial. The fitted
result was

```text
A_Nu = 1.915
B_Re = 0.510
```

Thus, after removing the direct rear sink, the fit no longer pushes `A_Nu` to
the old upper bound. This indicates that the prior high `A_Nu` signal was
coupled to the broader energy-path compensation, not only to channel
convection.

### 5.3 Exact NTU Cell Update

The heat-transfer coefficient is

```text
h_i = Nu_i k_f / D_h
```

The cell conductance is

```text
UA_i = h_i P_exchange dx
```

For active gas flow,

```text
epsilon_i = 1 - exp[-UA_i/(mdot_active cp)]
```

and

```text
T_g,out,i = T_g,in,i + epsilon_i (T_core,i - T_g,in,i)
Q_gas,i = mdot_active cp (T_g,out,i - T_g,in,i)
```

The reported mixed receiver gas temperature is

```text
T_g,mix,i = phi_act T_g,active,i + (1 - phi_act) T_in
```

Since `phi_act = 1` in fitted v28, this reduces to the active gas temperature.

---

## 6. Solid Energy Balances

### 6.1 Core Branch

For receiver cell `i`, the core energy balance is

```text
C_core,i dT_core,i/dt =
    Q_core w_core,i
  - Q_gas,i
  + G_cp dx (T_perim,i - T_core,i)
  + axial core conduction
  + receiver-tube coupling at i = N
```

The core heat capacity is

```text
C_core,i = rho_s cp_s(T_core,i) A_solid dx
```

Axial core conduction is represented by

```text
Q_cond,core,i+1/2 =
    k_core_scale k_s,face A_solid (T_core,i - T_core,i+1) / dx
```

where `k_core_scale = p20`.

The fitted v28 value is

```text
k_core_scale = 0
```

so the fitted model turns off effective axial core redistribution.

### 6.2 Perimeter Branch

The perimeter balance is

```text
C_perim,i dT_perim,i/dt =
    Q_perim w_perim,i
  - G_cp dx (T_perim,i - T_core,i)
  - G_cav dx (T_perim,i - T_cav)
  + axial perimeter conduction
  - front losses at i = 1
  + receiver-tube coupling at i = N
```

The fitted perimeter capacity is distributed uniformly:

```text
C_perim,i = C_perim_eff / N
```

and the fitted v28 value is

```text
C_perim_eff = 50 J/K
```

which is exactly the lower bound. This lower-bound result indicates that the
v28 topology is not using the perimeter as a physically meaningful thermal
inventory.

Perimeter axial conduction is

```text
k_perim(T) = k_perim_ref (T/900 K)^3
```

with fitted

```text
k_perim_ref = 0
```

Again, this is a lower-bound collapse rather than a validated material
property.

### 6.3 Core-Perimeter Coupling

The radial/core-perimeter exchange is

```text
Q_cp,i = G_cp dx (T_perim,i - T_core,i)
```

with fitted

```text
G_cp = 0.5 W/m/K
```

which is the lower bound. This is another sign that, after removal of the
distributed rear sink, the optimizer suppresses the two-zone interaction rather
than using it as physical heat redistribution.

---

## 7. Cavity and Front Losses

The perimeter branch couples to a lumped cavity/housing state through an
insulation conductance:

```text
Q_perim-cav,i = G_receiver-cavity dx (T_perim,i - T_cav)
```

The cavity loses heat to ambient by convection and radiation:

```text
Q_cav-amb =
    h_nat A_cav (T_cav - T_amb)
  + epsilon_metal sigma A_cav (T_cav^4 - T_amb^4)
```

The front receiver loss is applied to the first perimeter cell:

```text
Q_front =
    h_front A_frt (T_perim,1 - T_amb)
  + epsilon_front sigma A_frt (T_perim,1^4 - T_amb^4)
```

with fixed values

```text
h_front = 10 W/m2/K
epsilon_front = 0.95
```

---

## 8. Rear Tube, Flange, and Smooth Cooling Gate

### 8.1 Receiver-to-Tube Coupling

At the receiver exit, core and perimeter heat can enter the first rear-tube
solid node:

```text
Q_rec-tube,core =
    G_rec-tube f_core_tube (T_core,N - T_tube,1)

Q_rec-tube,perim =
    G_rec-tube (1 - f_core_tube) (T_perim,N - T_tube,1)
```

with fitted

```text
f_core_tube = 0.9995
```

so almost all receiver-to-tube solid heat is routed through the core exit.

### 8.2 Rear Tube Energy Balance

The rear tube has its own axial finite-volume conduction:

```text
Q_cond,tube,j+1/2 =
    k_tube,face A_tube,face (T_tube,j - T_tube,j+1) / dx_rear
```

It exchanges heat with rear gas:

```text
Q_gas,rear,j =
    mdot cp (T_g,out,j - T_g,in,j)
```

It also exchanges heat with the cavity over the cavity-coupled part of the rear
tube and with the water-flange pathway farther downstream:

```text
Q_tube-cav,j = G_tube-cav,j dx_rear (T_tube,j - T_cav)
Q_tube-flange,j = G_tube-flange,j dx_rear (T_tube,j - T_water)
```

### 8.3 Smooth Cooling Gate

The v27 model used a hard branch: one rear sink value for irradiated cases and
another for zero-irradiance cooling cases. v28 removes that direct sink and
instead lets the explicit flange path strengthen smoothly after lamp shutoff.

The lamp-off gate is

```text
s_off(t,I) =
    [1 / (1 + (I/I_gate)^4)] [1 - exp(-t/tau_cool)]
```

with

```text
I_gate = 1000 W/m2
tau_cool = p19
```

The effective flange multiplier is

```text
flange_scale_eff =
    flange_scale [1 + flange_cool_gain s_off(t,I)]
```

where

```text
flange_scale = p17
flange_cool_gain = p18
```

The fitted values are

```text
flange_scale = 0.1014
flange_cool_gain = 6.319
tau_cool = 116.8 s
```

For heating cases with nonzero irradiance, `s_off` is essentially zero. For
cooling cases, the gain approaches `1 + flange_cool_gain` over the fitted time
scale.

This is more continuous than the v27 branch, but it still represents an
operating-condition surrogate. It should be interpreted as an effective
cooling-side boundary change, not as an independently measured flange
coefficient.

---

## 9. Measurement Mapping

The model compares the following outputs to experiment:

| Experimental channel | v28 model output |
| :--- | :--- |
| T8 | perimeter temperature at 11 mm |
| T12 | perimeter temperature at 58 mm |
| T11 | perimeter temperature at 107 mm |
| T9 | core temperature at 58 mm |
| T10 | core temperature at 107 mm |
| T3 | gas temperature at 140 mm |
| T2 | lumped cavity temperature |

The most important v28 measurement change is:

```text
T3_model = T_gas(140 mm)
```

No wall-temperature fraction, thermocouple lag, or fitted T3 correction is
used. This makes v28 cleaner for scientific interpretation, but it also makes
T3 residuals more diagnostic of missing rear-gas/rear-tube physics.

For cooling cases, the first output row is aligned to the measured initial
temperature:

```text
outputs[1, 1:7] = experiment[1, 1:7]
```

This removes initial-condition mismatch from the cooling objective. It is
appropriate for transient-shape comparison, but it should be disclosed because
the initial cooling residual is zero by construction.

---

## 10. Calibration Objective

The calibration objective combines heating and cooling cases. For each sensor,
the signal loss is

```text
L_signal =
    normalized_mse(model, experiment)
  + w_slope normalized_slope_mse(model, experiment)
  + w_timing [(t90_model - t90_experiment)/duration]^2
```

with

```text
w_slope = 0.25
w_timing = 0.15
```

Heating cases also include an axial ordering penalty based on final-time
temperature offsets:

```text
T12 - T8
T11 - T8
T12 - T9
T11 - T10
```

Cooling cases include a monotonicity/upturn penalty on

```text
T11, T10, T3
```

with a large weight:

```text
w_upturn = 10000
```

The fitted participating heat capacity is weakly regularized toward the
measured assembly value:

```text
C_ref = 301 J/K
w_C = 0.10
```

In v28, however, the fitted participating heat capacity collapses to

```text
C_participating = 122.5 J/K
```

despite this regularization. This is strong evidence of model-form conflict.

---

## 11. Fitted Parameter Set

The full v28 runner completed with:

```text
objective = 117.71188433231664
return_code = MaxTime
```

The fitted parameter set is:

| Parameter | Value | Interpretation |
| :--- | ---: | :--- |
| `A_Nu` | 1.9151 | developing-flow Nu prefactor |
| `B_Re` | 0.5097 | Reynolds exponent |
| `C_Pr` | 0.3333 | fixed Prandtl exponent |
| `phi_0` | 1.0 | fixed full active flow |
| `m_rec` | 0.0 | fixed no recruitment trend |
| `front_dep` | 1.0 | fixed all core source in first cell |
| `scale_456` | 0.7023 | effective 456 kW/m2 power scale |
| `scale_304` | 0.8758 | effective 304 kW/m2 power scale |
| `scale_256` | 0.9516 | effective 256 kW/m2 power scale |
| `G_core_perim` | 0.5 W/m/K | lower bound |
| `C_perim_eff` | 50.0 J/K | lower bound |
| `k_perim_ref` | 0.0 W/m/K | lower bound |
| `beta_opt` | 184.7 1/m | fixed |
| `spill_capture` | 0.3518 | fitted perimeter spill capture |
| `beta_perim` | 6.873 1/m | perimeter source attenuation |
| `f_core_tube` | 0.9995 | near all exit solid heat routed through core |
| `flange_scale` | 0.1014 | near lower bound |
| `flange_cool_gain` | 6.319 | cooling-side flange gain |
| `flange_cool_tau` | 116.8 s | cooling gate time constant |
| `k_core_axial_scale` | 0.0 | lower bound |

Several fitted values at lower bounds should not be interpreted as measured
physical properties. They are signatures of a failed or under-complete model
structure.

---

## 12. Scientific Interpretation of the v28 Result

### 12.1 What v28 Improved Conceptually

v28 removes three major interpretability problems from v27:

1. It eliminates the distributed direct rear-core sink.
2. It replaces the hard lamp-on/lamp-off branch with a continuous gate.
3. It removes the fitted T3 wall-temperature blend and uses pure gas
   temperature.

These choices make the model more transparent and easier to audit.

### 12.2 What v28 Falsified

The full fit shows that the stricter explicit rear-tube/flange pathway does
not, by itself, reproduce the data. The objective is much larger than v27, and
the fitted parameters collapse to bounds:

```text
G_core_perim -> lower bound
C_perim_eff -> lower bound
k_perim_ref -> lower bound
k_core_axial_scale -> lower bound
flange_scale -> near lower bound
```

This pattern means that removing the direct rear sink did not reveal a clean
physical solution. Instead, it revealed that the model lacks another physical
mechanism, most likely rear/adaptor thermal storage and delayed heat exchange.

### 12.3 Implication for Scientific Credit

v28 should be cited as evidence that:

- T3 should be tested as a gas measurement before applying wall or lag
  corrections.
- A hard irradiance-branch rear sink is scientifically undesirable.
- A direct distributed rear sink is too artificial for final coefficient
  identification.
- The current explicit rear tube/flange path is incomplete without additional
  rear/adaptor inventory.

v28 should not be cited as a validated coefficient-identification model.
Specifically, its fitted `A_Nu`, `G_core_perim`, `C_perim_eff`, and
`k_perim_ref` values should not be reported as physical receiver coefficients.

---

## 13. Recommended Next Model Revision

The next defensible revision should not restore the v27 distributed rear sink.
Instead, it should introduce a physically bounded rear/adaptor reservoir:

```text
receiver exit solid -> rear/adaptor reservoir -> explicit rear tube/flange/cavity
```

A minimal v29 state addition would be:

```text
C_rear dT_rear/dt =
    G_core-rear (T_core,N - T_rear)
  + G_perim-rear (T_perim,N - T_rear)
  - G_rear-tube (T_rear - T_tube,1)
  - G_rear-cav (T_rear - T_cav)
```

with `C_rear` bounded by the measured alumina/adaptor mass and heat capacity,
and conductances tied to contact areas or allowed only within physically
plausible ranges.

This revision would preserve the main scientific improvement of v28: heat must
leave through an interpretable topology. It would also restore the delayed rear
thermal inventory that v28 appears to lack.

T3 should remain pure gas for this next revision. A wall blend, thermocouple
lag, or fitted axial T3 position should only be reconsidered after the rear
thermal reservoir is physically represented.

---

## 14. Conclusion

The 1D_v28 model is a manuscript-relevant structural experiment. It improves
interpretability by removing an artificial direct rear-core sink, smoothing the
cooling-side boundary change, widening the Nusselt prefactor bound, and mapping
T3 directly to gas temperature. However, the fitted result is poor and
parameter-bound dominated. The model therefore demonstrates that the v27 direct
rear sink was compensating for missing physical topology, not merely correcting
a numerical inconvenience.

The central scientific conclusion is:

```text
The receiver model needs an explicit rear/adaptor thermal inventory and
physically constrained rear heat pathway. A distributed direct rear sink is not
acceptable as final physics, but removing it without replacing the missing rear
reservoir makes the two-zone 1D model structurally under-complete.
```

