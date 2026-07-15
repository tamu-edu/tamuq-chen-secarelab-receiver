
# A. Theoretical Analysis for Transient Parameter Identification in an End-Heated Ceramic Solar Receiver

## 1. Measurement configuration

Take the irradiated front face as $z=0$ and the rear face as $z=L$, with

$$
L \simeq 137\ \mathrm{mm}.
$$

The available measurements are approximately:

| Position                                  | Measurement interpretation                               |
| ----------------------------------------- | -------------------------------------------------------- |
| $z\simeq10$ mm                            | Ceramic side-wall temperature near the irradiated face   |
| $z\simeq58$ mm, centre                    | Gas-exposed thermocouple or internal ceramic temperature |
| $z\simeq58$ mm, wall                      | Ceramic wall temperature                                 |
| $z\simeq110$ mm, centre                   | Gas-exposed thermocouple or internal ceramic temperature |
| $z\simeq110$ mm, wall                     | Ceramic wall temperature                                 |
| $z\simeq137$ mm                           | Rear-end ceramic temperature                             |
| Near $z\simeq L/2$, 40 mm into insulation | Internal insulation temperature                          |

The experiments contain multiple heating ramps at different:

$$
\dot m,\qquad q''*{\mathrm{sol}}(t),\qquad T*{\mathrm{in}},
$$

but no dedicated cooling transients.

The desired quantities include:

$$
k_s,\quad \rho_sc_{p,s},\quad \alpha_s,
$$

$$
k_{\mathrm{ins}},\quad \rho_{\mathrm{ins}}c_{p,\mathrm{ins}},\quad
\alpha_{\mathrm{ins}},
$$

$$
h_i(\dot m,T),\quad h_{\mathrm{ext}},\quad \epsilon,
$$

$$
\eta_{\mathrm{abs}},\quad R''_{\mathrm{contact}},
$$

and the dynamic response of the temperature sensors.

---

## 2. Recommended hierarchy of models

Three levels should be used.

### Level 1: apparent first-order fits

Each trace is first fitted with one or more exponential ramp responses. This gives apparent gains, delays and time constants without yet assigning them to physical properties.

### Level 2: axial RC-network model

The receiver is divided into axial sections, with ceramic capacitance, axial conduction, gas heat removal and insulation losses. This model is fast enough for global parameter estimation over all experiments.

### Level 3: axisymmetric conjugate model

A two-dimensional $r-z$ ceramic and insulation conduction model, coupled to a flowing-gas energy equation, is used to validate the final parameters and resolve radial gradients.

The 0D model alone cannot explain different response times at 10, 58, 110 and 137 mm. Those delays contain the strongest information about ceramic thermal diffusivity.

---

## 3. Full forward model

### 3.1 Ceramic conduction

In the ceramic,

$$
\rho_s c_{p,s}(T_s)
\frac{\partial T_s}{\partial t}
=

\frac{1}{r}
\frac{\partial}{\partial r}
\left(
r k_s(T_s)\frac{\partial T_s}{\partial r}
\right)
+
\frac{\partial}{\partial z}
\left(
k_s(T_s)\frac{\partial T_s}{\partial z}
\right).
$$

Initially,

$$
T_s(r,z,0)=T_{s,0}(r,z).
$$

The absorbed solar boundary condition at the front face is

$$
-k_s\left.\frac{\partial T_s}{\partial z}\right|_{z=0}
=

\eta_{\mathrm{abs}}q''*{\mathrm{sol}}(t)
-q''*{\mathrm{front,loss}}.
$$

The front loss may include

$$
q''_{\mathrm{front,loss}}
=

h_f(T_s-T_\infty)
+
\epsilon_f\sigma(T_s^4-T_{\mathrm{sur}}^4).
$$

At the rear face,

$$
-k_s\left.\frac{\partial T_s}{\partial z}\right|_{z=L}
=

h_L(T_s-T_\infty)
+
\epsilon_L\sigma(T_s^4-T_{\mathrm{sur}}^4).
$$

### 3.2 Gas model

A practical one-dimensional gas equation is

$$
\rho_f c_{p,f}A_f
\frac{\partial T_g}{\partial t}
+
\dot m c_{p,f}
\frac{\partial T_g}{\partial z}
=

h_iP_i(T_w-T_g),
$$

with

$$
T_g(0,t)=T_{\mathrm{in}}(t).
$$

If gas residence time is much shorter than the ceramic response time, the storage term can be neglected:

$$
\dot m c_{p,f}
\frac{\partial T_g}{\partial z}
=

h_iP_i(T_w-T_g).
$$

The local heat flux from ceramic to gas is

$$
q''_{s\rightarrow g}(z,t)
=

h_i(\dot m,T)
[T_w(z,t)-T_g(z,t)].
$$

### 3.3 Insulation

The insulation should be modeled dynamically because a thermocouple lies 40 mm inside it:

$$
\rho_{\mathrm{ins}}c_{p,\mathrm{ins}}
\frac{\partial T_{\mathrm{ins}}}{\partial t}
=

\frac{1}{r}
\frac{\partial}{\partial r}
\left(
rk_{\mathrm{ins}}
\frac{\partial T_{\mathrm{ins}}}{\partial r}
\right)
+
\frac{\partial}{\partial z}
\left(
k_{\mathrm{ins}}
\frac{\partial T_{\mathrm{ins}}}{\partial z}
\right).
$$

At the ceramic–insulation interface,

$$
-k_s\frac{\partial T_s}{\partial n}
=

\frac{T_s-T_{\mathrm{ins}}}{R''_{\mathrm{contact}}}
-

-k_{\mathrm{ins}}\frac{\partial T_{\mathrm{ins}}}{\partial n}.
$$

At the external insulation surface,

$$
-k_{\mathrm{ins}}\frac{\partial T_{\mathrm{ins}}}{\partial n}
=

h_{\mathrm{ext}}(T_{\mathrm{ins}}-T_\infty)
+
\epsilon_{\mathrm{ext}}\sigma
(T_{\mathrm{ins}}^4-T_{\mathrm{sur}}^4).
$$

The insulation thermocouple is especially valuable because it distinguishes transient storage within the insulation from instantaneous heat loss through a simple thermal resistance.

---

## 4. Reduced axial model for inverse fitting

Divide the receiver into $N$ axial segments. A good initial discretization places nodes near

$$
z=0,\ 10,\ 58,\ 110,\ 137\ \mathrm{mm}.
$$

For ceramic node $i$,

$$
C_{s,i}\frac{dT_{s,i}}{dt}
=

G_{z,i-1}(T_{s,i-1}-T_{s,i})
+
G_{z,i}(T_{s,i+1}-T_{s,i})
-\dot Q_{g,i}
-\dot Q_{\mathrm{loss},i}.
$$

Here,

$$
C_{s,i}
=

\rho_sc_{p,s}A_s\Delta z_i,
$$

and

$$
G_{z,i}
=

\frac{k_sA_s}{\Delta z_{i+1/2}}.
$$

Under quasi-steady gas flow through segment $i$,

$$
T_{g,i+1}
=

T_{s,i}
-

(T_{s,i}-T_{g,i})
\exp\left(
-\frac{h_iP_i\Delta z_i}
{\dot m c_{p,f}}
\right).
$$

Thus,

$$
\dot Q_{g,i}
=

\dot m c_{p,f}(T_{g,i+1}-T_{g,i}).
$$

The side loss is represented either by a conductance,

$$
\dot Q_{\mathrm{loss},i}
=

G_{\mathrm{loss},i}(T_{s,i}-T_\infty),
$$

or by one or more dynamic insulation nodes:

$$
C_{\mathrm{ins},i}\frac{dT_{\mathrm{ins},i}}{dt}
=

G_{s-\mathrm{ins},i}(T_{s,i}-T_{\mathrm{ins},i})
-

G_{\mathrm{out},i}(T_{\mathrm{ins},i}-T_\infty).
$$

The front-node input is

$$
\dot Q_{\mathrm{abs}}(t)
=

\eta_{\mathrm{abs}}
q''*{\mathrm{sol}}(t)A*{\mathrm{irr}}.
$$

This produces a state-space model,

$$
\mathbf C\dot{\mathbf x}
=

\mathbf K(\boldsymbol\theta,\dot m,T)\mathbf x
+
\mathbf B\dot Q_{\mathrm{abs}}(t)
+
\mathbf d,
$$

where $\boldsymbol\theta$ is the parameter vector to be estimated.

---

## 5. Analytical response to a heating ramp

Before fitting the physical model, each temperature trace can be characterized by a modal ramp response.

For a single first-order system,

$$
C\frac{dT}{dt}+K(T-T_0)=\dot Q(t),
$$

define

$$
\tau=\frac CK.
$$

For a linear power ramp

$$
\dot Q(t)=s_Qt,\qquad t\geq0,
$$

the temperature response is

$$
\boxed{
T(t)-T_0
=

\frac{s_Q}{K}
\left[
t-\tau\left(1-e^{-t/\tau}\right)
\right].
}
$$

For a distributed receiver, each sensor generally sees several modes:

$$
\boxed{
T_j(t)-T_{j,0}
=

\sum_{n=1}^{M}
g_{jn}s_Q
\left[
t-\tau_n(1-e^{-t/\tau_n})
\right].
}
$$

Here:

* $\tau_n$ are apparent thermal time constants;
* $g_{jn}$ are sensor-dependent modal gains;
* slower modes normally represent axial conduction and environmental losses;
* faster modes represent local wall, gas or thermocouple response.

If the ramp stops at $t=t_r$, use superposition:

$$
\dot Q(t)
=

s_Qt-s_Q(t-t_r)H(t-t_r).
$$

For each mode,

$$
\Delta T_{j,n}
=

g_{jn}s_Q
\left[
F_n(t)-H(t-t_r)F_n(t-t_r)
\right],
$$

where

$$
F_n(t)=t-\tau_n(1-e^{-t/\tau_n}).
$$

This gives robust initial estimates for the physical inverse model even though cooling curves are unavailable.

---

## 6. What each measurement identifies

### 6.1 Axial ceramic diffusivity

The sequential response at

$$
10,\ 58,\ 110,\ 137\ \mathrm{mm}
$$

is primarily controlled at early time by

$$
\alpha_s=\frac{k_s}{\rho_sc_{p,s}}.
$$

A rough scaling is

$$
t_{\mathrm{propagation}}
\sim
\frac{(\Delta z)^2}{\alpha_s}.
$$

The full inverse model should be used for numerical values, but this scaling provides the initial estimate.

Early-time phase delays and curvature are generally more informative about $\alpha_s$ than absolute temperatures.

### 6.2 Ceramic conductivity versus heat capacity

Transient propagation mainly identifies $\alpha_s$, while absolute temperature gradients identify $k_s$. Then

$$
\rho_sc_{p,s}=\frac{k_s}{\alpha_s}.
$$

Without reliable spatial gradients or an independent value of either $k_s$ or $c_p$, conductivity and volumetric heat capacity can remain correlated.

### 6.3 Internal gas heat-transfer coefficient

At 58 and 110 mm, the difference between wall and gas-exposed measurements is sensitive to

$$
h_i(\dot m,T).
$$

A useful parameterization is

$$
h_i
=

C_h\frac{k_f}{D_h}
Nu_{\mathrm{reference}}(Re,Pr,z/D_h),
$$

where $C_h$ is fitted globally.

Flow-rate variation is essential:

* ceramic properties are common to all runs;
* gas-side heat removal changes strongly with $\dot m$;
* environmental loss changes mainly with temperature, not directly with flow.

If the centre sensors are truly exposed to gas, their response should change significantly with flow rate. If they are embedded in ceramic, they should follow a conduction-dominated response. Both hypotheses can be fitted and compared using residual structure or an information criterion.

### 6.4 Insulation conductivity and diffusivity

The wall temperature near mid-length and the thermocouple 40 mm into the insulation provide:

$$
\alpha_{\mathrm{ins}}
=

\frac{k_{\mathrm{ins}}}
{\rho_{\mathrm{ins}}c_{p,\mathrm{ins}}},
$$

from response delay and attenuation, and information on $k_{\mathrm{ins}}$ from the temperature difference across the insulation.

A contact resistance may be detectable if the wall rises rapidly while the nearby insulation response shows an additional delay:

$$
R''_{\mathrm{contact}}
=

\frac{T_s-T_{\mathrm{ins,surface}}}{q''}.
$$

However, one internal insulation sensor cannot uniquely identify all of

$$
k_{\mathrm{ins}},\quad
\rho_{\mathrm{ins}}c_p,\quad
R''*{\mathrm{contact}},\quad
h*{\mathrm{ext}}
$$

without priors. At least one or two should be supplied from material data or constrained ranges.

### 6.5 Rear-face loss

The thermocouple at approximately 137 mm is particularly sensitive to:

$$
h_LA_L,\qquad
\epsilon_LA_L,
$$

and to axial conduction through the ceramic.

Rear-face loss is most identifiable during the later portion of a long heating ramp, after the thermal front has reached the end.

### 6.6 Absorbed solar fraction

The absorbed fraction enters as

$$
\eta_{\mathrm{abs}}
=

\frac{\dot Q_{\mathrm{absorbed}}}
{q''*{\mathrm{sol}}A*{\mathrm{irr}}}.
$$

It is inferred from the overall energy balance only after storage, gas heat removal and losses are represented correctly:

$$
\dot Q_{\mathrm{abs}}
=

\frac{dU_s}{dt}
+
\frac{dU_{\mathrm{ins}}}{dt}
+
\dot Q_g
+
\dot Q_{\mathrm{loss}}.
$$

If the experiments do not approach a plateau, $\eta_{\mathrm{abs}}$ and thermal capacitance can be strongly correlated.

---

## 7. Sensor dynamics

Each thermocouple should have an observation model:

$$
\tau_{\mathrm{TC},j}
\frac{d\widehat T_j}{dt}
+
\widehat T_j
=

T_{\mathrm{local},j}+b_j.
$$

Here:

* $\widehat T_j$ is the predicted measured temperature;
* $T_{\mathrm{local},j}$ is the model field temperature;
* $\tau_{\mathrm{TC},j}$ is the sensor time constant;
* $b_j$ is a calibration offset.

Ignoring thermocouple lag will generally bias thermal diffusivity downward because sensor delay can be mistaken for conduction delay.

A thermocouple exposed to the gas may not measure the true gas temperature. It may receive radiation from the hot ceramic and conduction through its leads. A more detailed model is

$$
C_{\mathrm{TC}}\frac{dT_{\mathrm{TC}}}{dt}
=

h_{\mathrm{TC}}A_{\mathrm{TC}}(T_g-T_{\mathrm{TC}})
+
\epsilon_{\mathrm{TC}}\sigma A_{\mathrm{TC}}
(T_{\mathrm{sur}}^4-T_{\mathrm{TC}}^4)
+
G_{\mathrm{wire}}(T_{\mathrm{support}}-T_{\mathrm{TC}}).
$$

If these quantities are unknown, the centre-sensor temperature should be treated as an effective gas-related measurement, not automatically as the bulk gas temperature.

---

## 8. Global inverse problem

All experiments should be fitted simultaneously.

Let $e$ denote the experiment and $j$ the sensor. The objective may be written as

$$
J(\boldsymbol\theta)
=

\sum_e\sum_j\sum_k
w_{ejk}
\left[
T^{\mathrm{meas}}_{ej}(t_k)
-

T^{\mathrm{model}}*{ej}(t_k;\boldsymbol\theta)
\right]^2
+
J*{\mathrm{regularization}}.
$$

Parameters common to every experiment should include

$$
k_s(T),\quad
\rho_sc_{p,s}(T),\quad
k_{\mathrm{ins}}(T),\quad
\rho_{\mathrm{ins}}c_{p,\mathrm{ins}},
$$

$$
R''*{\mathrm{contact}},\quad
\epsilon,\quad
h*{\mathrm{ext}},\quad
\eta_{\mathrm{abs}}.
$$

Flow-dependent quantities include

$$
h_i(\dot m,T)
$$

or a common multiplier $C_h$ applied to a selected Nusselt-number model.

Experiment-specific nuisance parameters may include

$$
t_{0,e},\quad
q''*{\mathrm{scale},e},\quad
T*{\infty,e},\quad
T_{\mathrm{in},e},
$$

along with small sensor offsets.

A separate $h_i$ for every experiment should be avoided initially. A shared physical correlation with one or two fitted correction coefficients gives much better identifiability.

---

## 9. Temperature-dependent rather than arbitrary time-dependent properties

Material and transfer properties should first be parameterized as functions of temperature and flow, not as unconstrained functions of time.

For example,

$$
k_s(T)
=

k_{s,0}
\left[
1+a_k(T-T_{\mathrm{ref}})
\right],
$$

$$
c_{p,s}(T)
=

c_{p,0}
\left[
1+a_c(T-T_{\mathrm{ref}})
\right],
$$

and

$$
h_i(\dot m,T)
=

C_h\frac{k_f(T)}{D_h}
Nu[Re(\dot m,T),Pr(T)].
$$

Radiation already introduces a strong temperature dependence:

$$
\dot Q_{\mathrm{rad}}
=

\epsilon\sigma A(T^4-T_{\mathrm{sur}}^4).
$$

An arbitrary fitted $k_s(t)$, $h_i(t)$, or $\eta(t)$ would absorb model errors and would generally not be uniquely identifiable from heating ramps.

Only after a constant or temperature-dependent model has been rejected by systematic residuals should slow temporal variation, such as surface absorptivity degradation, be introduced.

---

## 10. Recommended staged extraction

### Stage A: preprocessing

For every experiment:

1. Synchronize solar-flux, mass-flow and temperature channels.
2. Estimate the true ramp start time.
3. Remove pre-ramp offsets using the initial equilibrium period.
4. Retain the measured solar-flux history rather than replacing it with an ideal ramp.
5. Record $T_{\mathrm{in}}(t)$, $T_\infty(t)$ and $\dot m(t)$.

### Stage B: apparent modal fitting

Fit each trace with one-, two- and three-mode ramp responses. Extract:

$$
\tau_{j,1},\tau_{j,2},\ldots,
$$

the initial delay, and ramp gains.

The increase of apparent time constant with axial position gives an immediate check for conduction-dominated behavior.

### Stage C: ceramic-only early-time fit

Use the early portions of the 10, 58, 110 and 137 mm wall traces to estimate:

$$
\alpha_s,\qquad \tau_{\mathrm{TC},j},
$$

while initially constraining losses and $h_i$ to plausible values.

### Stage D: flow-dependent fit

Fit experiments at different mass flow rates jointly to estimate:

$$
C_h
$$

or the parameters of $h_i(\dot m,T)$.

The common ceramic parameters must remain fixed across flow rates.

### Stage E: insulation and external-loss fit

Use the mid-length ceramic-wall and insulation temperatures to estimate:

$$
\alpha_{\mathrm{ins}},\quad
k_{\mathrm{ins}},\quad
R''*{\mathrm{contact}},\quad
h*{\mathrm{ext}}.
$$

Use prior values for at least one insulation property if only one embedded insulation sensor is available.

### Stage F: high-flux nonlinear fit

Use different heat-flux levels to separate approximately linear losses from radiation:

$$
\dot Q_{\mathrm{loss}}
=

K_{\mathrm{linear}}(T-T_\infty)
+
\epsilon\sigma A(T^4-T_{\mathrm{sur}}^4).
$$

Low-flux tests are more sensitive to linear conductances. High-flux tests provide the information needed for emissivity.

## Stage G: absorbed power and final global fit

Finally fit

$$
\eta_{\mathrm{abs}}
$$

and all remaining parameters globally, followed by uncertainty analysis.

---

## 11. Identifiability with heating ramps only

Heating data are sufficient to estimate many parameters, particularly because several flow rates and flux levels are available. However, the following correlations must be expected.

| Parameter pair                                  | Identifiability issue                                                               |
| ----------------------------------------------- | ----------------------------------------------------------------------------------- |
| $k_s$ and $\rho_sc_p$                           | Early transients identify their ratio $\alpha_s$ more strongly than each separately |
| $\eta_{\mathrm{abs}}$ and $C_s$                 | Both affect heating rate if losses and gradients are not independently constrained  |
| $h_{\mathrm{ext}}$ and $\epsilon$               | Difficult to separate over a narrow temperature range                               |
| $h_i$ and gas-sensor lag                        | Both can delay or attenuate a gas-exposed thermocouple                              |
| $k_{\mathrm{ins}}$ and $R''_{\mathrm{contact}}$ | Both contribute to ceramic-to-insulation temperature difference                     |
| $k_{\mathrm{ins}}$ and $\rho_{\mathrm{ins}}c_p$ | Dynamic data identify diffusivity more strongly than both individual properties     |

Heating-only data are weakest for loss coefficients if none of the ramps is long enough for the temperatures to show substantial curvature toward equilibrium. A perfectly linear temperature rise contains little information about environmental loss.

Cooling curves would improve identifiability, but they are not mathematically required if the heating ramps cover:

* several mass flow rates;
* several heat-flux levels;
* sufficiently long durations;
* a broad temperature range;
* accurately measured input histories.

---

## 12. Minimum recommended fitted parameter set

To avoid over-parameterization, begin with

$$
\boldsymbol\theta_{\min}
=

\left[
\alpha_s,,
k_s,,
\alpha_{\mathrm{ins}},,
k_{\mathrm{ins}},,
C_h,,
G_{\mathrm{ext}},,
\eta_{\mathrm{abs}},,
h_LA_L
\right].
$$

Use known geometry and material density to calculate

$$
\rho_sc_{p,s}=\frac{k_s}{\alpha_s},
\qquad
\rho_{\mathrm{ins}}c_{p,\mathrm{ins}}
=

\frac{k_{\mathrm{ins}}}{\alpha_{\mathrm{ins}}}.
$$

Add only when justified by residuals:

$$
R''*{\mathrm{contact}},\quad
\epsilon,\quad
k_s(T),\quad
c_p(T),\quad
\tau*{\mathrm{TC},j}.
$$

The preferred order is therefore:

$$
\boxed{
\text{modal ramp fits}
\rightarrow
\text{axial RC inverse model}
\rightarrow
\text{axisymmetric validation}.
}
$$

The most robust outputs from the existing measurements should be:

1. ceramic axial thermal diffusivity;
2. effective ceramic conductivity or heat capacity when combined with priors;
3. flow-dependent internal heat-transfer coefficient;
4. insulation diffusivity and effective radial conductance;
5. total environmental loss conductance;
6. effective absorbed solar fraction;
7. sensor-specific lags and possible gas-sensor radiation bias.

# B. 0D single node model
Yes. The appropriate formulation is a **zero-dimensional, single-state model** in which the only dynamic state is the average ceramic temperature,

$$
T_{s,\mathrm{avg}}(t),
$$

and the gas outlet temperature,

$$
T_{g,\mathrm{out}}(t),
$$

is calculated algebraically.

Strictly speaking, it is a **one-node 0D model**. A model with zero thermal nodes would be purely steady-state and could not reproduce the heating ramps.

The spatial thermocouples are used during model development to identify correction coefficients that account for the receiver length and its nonuniform axial temperature profile. Those local temperatures do not need to remain as states in the final control model.

## 1. Definition of the average solid temperature

The physically appropriate average is an energy-weighted temperature:

$$
T_{s,\mathrm{avg}}
=

\frac{
\displaystyle
\int_{V_s}
\rho_s c_{p,s} T_s,dV
}{
\displaystyle
\int_{V_s}
\rho_s c_{p,s},dV
}.
$$

For constant density, heat capacity and cross-sectional area,

$$
T_{s,\mathrm{avg}}
=

\frac{1}{L}
\int_0^L T_s(z),dz.
$$

The local thermocouples can be used to approximate this integral.

If solid-wall temperatures are available at $10$, $58$, $110$, and $137\ \mathrm{mm}$, and the profile is taken as piecewise linear, then approximately

$$
T_{s,\mathrm{avg}}
\approx
0.248,T_{10}
+
0.365,T_{58,w}
+
0.288,T_{110,w}
+
0.099,T_{137,w}.
$$

If the $137\ \mathrm{mm}$ sensor measures gas rather than ceramic temperature, a simpler estimate is

$$
\boxed{
T_{s,\mathrm{avg}}
\approx
0.248,T_{10}
+
0.365,T_{58,w}
+
0.387,T_{110,w}.
}
$$

This assumes that the wall temperature between $110$ and $137\ \mathrm{mm}$ is approximately represented by the $110\ \mathrm{mm}$ wall measurement.

The centre thermocouples should not be included in the solid average unless they are confirmed to be embedded in the ceramic.

## 2. Single-state receiver energy balance

The receiver energy balance is

$$
\boxed{
C_{\mathrm{eff}}
\frac{dT_{s,\mathrm{avg}}}{dt}
=

\dot Q_{\mathrm{abs}}
-

\dot Q_g
-

\dot Q_{\mathrm{loss}}.
}
$$

The effective thermal capacitance is

$$
C_{\mathrm{eff}}
=

\gamma_C \rho_s c_{p,s}V_s,
$$

where

$$
V_s=A_sL.
$$

The coefficient $\gamma_C$ is an active-mass correction factor. Ideally,

$$
\gamma_C=1.
$$

A fitted value below one indicates that only part of the ceramic participates significantly over the duration of the heating ramp.

## 3. Effective absorbed solar power

The absorbed power is written as

$$
\boxed{
\dot Q_{\mathrm{abs}}
=

\eta_{\mathrm{eff}}
q''*{\mathrm{sol}}A*{\mathrm{irr}}.
}
$$

Here $\eta_{\mathrm{eff}}$ includes:

* surface absorptivity;
* aperture losses;
* reflections;
* flux nonuniformity;
* radiation penetration into the receiver;
* any difference between measured incident and absorbed power.

The observed temperature maximum near $58\ \mathrm{mm}$ is not modeled explicitly. Its effect is absorbed into $\eta_{\mathrm{eff}}$ and the gas-coupling correction described below.

## 4. Length-corrected gas outlet relation

For a spatially varying wall temperature, the gas equation is approximately

$$
\dot m c_{p,g}
\frac{dT_g}{dz}
=

U'(z)\left[T_w(z)-T_g(z)\right],
$$

where (U'(z)) is the local wall-to-gas conductance per unit length.

Define

$$
N(z)
=

\int_0^z
\frac{U'(\zeta)}
{\dot m c_{p,g}}
,d\zeta.
$$

The exact outlet solution is

$$
\begin{aligned}
T_{g,\mathrm{out}}
={}&
T_{g,\mathrm{in}}e^{-N(L)}
\\
&+
\int_0^L
e^{-[N(L)-N(z)]}
\frac{U'(z)}
{\dot m c_{p,g}}
T_w(z),dz.
\end{aligned}
$$

This can be written in the usual heat-exchanger form

$$
T_{g,\mathrm{out}}
=

T_{g,\mathrm{in}}
+
\left(1-e^{-N_L}\right)
\left(
T_{w,\mathrm{eff}}-T_{g,\mathrm{in}}
\right),
$$

where

$$
N_L=N(L),
$$

and $T_{w,\mathrm{eff}}$ is a conductance- and flow-weighted wall temperature.

Because the final 0D model does not contain $T_w(z)$, define a length correction coefficient

$$
\boxed{
\chi_L
=

\frac{
T_{w,\mathrm{eff}}-T_{g,\mathrm{in}}
}{
T_{s,\mathrm{avg}}-T_{g,\mathrm{in}}
}.
}
$$

The outlet-temperature model then becomes

$$
\boxed{
T_{g,\mathrm{out}}
=

T_{g,\mathrm{in}}
+
\chi_L
\left(1-e^{-NTU_{\mathrm{eff}}}\right)
\left(
T_{s,\mathrm{avg}}-T_{g,\mathrm{in}}
\right).
}
$$

Define an overall gas effectiveness

$$
\boxed{
\varepsilon_{\mathrm{eff}}
=

\chi_L
\left(1-e^{-NTU_{\mathrm{eff}}}\right).
}
$$

Then

$$
\boxed{
T_{g,\mathrm{out}}
=

T_{g,\mathrm{in}}
+
\varepsilon_{\mathrm{eff}}
\left(
T_{s,\mathrm{avg}}-T_{g,\mathrm{in}}
\right).
}
$$

This is the central algebraic output equation.

The coefficient $\chi_L$ accounts for the axial temperature profile. In particular, because the $58\ \mathrm{mm}$ region is hotter than the average,

$$
\chi_L
$$

may be greater than one. That does not imply an effectiveness greater than the physical wall-to-gas limit; it means the gas-weighted wall temperature is higher than the volume-averaged ceramic temperature.

## 5. Length dependence of the effective $UA$

The effective number of transfer units is

$$
NTU_{\mathrm{eff}}
=

\frac{UA_{\mathrm{eff}}}
{\dot m c_{p,g}}.
$$

The receiver length enters through

$$
UA_{\mathrm{eff}}
=

h_{\mathrm{eff}}P_iL.
$$

A convenient empirical flow-rate relation is

$$
h_{\mathrm{eff}}
=

a_h\dot m^{,n}.
$$

Therefore,

$$
\boxed{
NTU_{\mathrm{eff}}
=

\frac{
a_hP_iL\dot m^{,n}
}{
\dot m c_{p,g}
}
=

\frac{
a_hP_iL
}{
c_{p,g}
}
\dot m^{,n-1}.
}
$$

Thus, the final outlet relation may be written as

$$
\boxed{
T_{g,\mathrm{out}}
=

T_{g,\mathrm{in}}
+
\chi_L
\left[
1-
\exp\left(
-\frac{
a_hP_iL
}{
c_{p,g}
}
\dot m^{,n-1}
\right)
\right]
\left(
T_{s,\mathrm{avg}}-T_{g,\mathrm{in}}
\right).
}
$$

This form explicitly accounts for receiver length.

## 6. Heat transferred to the gas

The gas heat removal is

$$
\dot Q_g
=

\dot m c_{p,g}
\left(
T_{g,\mathrm{out}}-T_{g,\mathrm{in}}
\right).
$$

Using the effective outlet relation,

$$
\boxed{
\dot Q_g
=

K_g
\left(
T_{s,\mathrm{avg}}-T_{g,\mathrm{in}}
\right),
}
$$

where

$$
\boxed{
K_g
=

\dot m c_{p,g}\varepsilon_{\mathrm{eff}}.
}
$$

Therefore,

$$
K_g
=

\dot m c_{p,g}
\chi_L
\left(1-e^{-NTU_{\mathrm{eff}}}\right).
$$

This effective gas conductance includes both the gas heat-transfer coefficient and the axial temperature-shape correction.

## 7. Length-corrected environmental losses

The total environmental loss can be written as

$$
\dot Q_{\mathrm{loss}}
=

\dot Q_{\mathrm{side}}
+
\dot Q_{\mathrm{front}}
+
\dot Q_{\mathrm{rear}}.
$$

For linear losses,

$$
\boxed{
\dot Q_{\mathrm{loss}}
=

K_{\mathrm{loss}}
\left(
T_{s,\mathrm{avg}}-T_\infty
\right).
}
$$

The effective loss conductance may be expressed as

$$
\boxed{
K_{\mathrm{loss}}
=

\psi_{\mathrm{side}}U_{\mathrm{side}}P_oL
+
\psi_FU_FA_F
+
\psi_RU_RA_R.
}
$$

The coefficients

$$
\psi_{\mathrm{side}},\qquad
\psi_F,\qquad
\psi_R
$$

are temperature-shape factors. They convert local surface temperatures into losses based on the average ceramic temperature.

For example,

$$
\psi_{\mathrm{side}}
=

\frac{
\displaystyle
\int_0^L
\left[T_{\mathrm{side}}(z)-T_\infty\right]dz
}{
L\left[T_{s,\mathrm{avg}}-T_\infty\right]
}.
$$

These coefficients can be estimated using the axial wall measurements and the insulation thermocouple.

## 8. Radiation loss

At elevated temperatures, radiation may be retained as

$$
\boxed{
\dot Q_{\mathrm{rad}}
=

\chi_{\mathrm{rad}}
\epsilon_{\mathrm{eff}}
\sigma A_{\mathrm{rad}}
\left(
T_{s,\mathrm{avg}}^4-T_{\mathrm{sur}}^4
\right).
}
$$

The coefficient $\chi_{\mathrm{rad}}$ corrects for the fact that

$$
\overline{T^4}\neq\overline T^{,4}
$$

when the temperature distribution is strongly nonuniform.

The complete nonlinear heat-loss model is then

$$
\begin{aligned}
\dot Q_{\mathrm{loss}}
={}&
K_{\mathrm{lin}}
\left(
T_{s,\mathrm{avg}}-T_\infty
\right)
\\
&+
\chi_{\mathrm{rad}}
\epsilon_{\mathrm{eff}}
\sigma A_{\mathrm{rad}}
\left(
T_{s,\mathrm{avg}}^4-T_{\mathrm{sur}}^4
\right).
\end{aligned}
$$

For a simpler control model, radiation may be linearized around an operating temperature $T_*$:

$$
K_{\mathrm{rad},*}
=

4\chi_{\mathrm{rad}}
\epsilon_{\mathrm{eff}}
\sigma A_{\mathrm{rad}}T_*^3.
$$

The total linearized conductance becomes

$$
K_{\mathrm{loss},*}
=

K_{\mathrm{lin}}+K_{\mathrm{rad},*}.
$$

## 9. Final nonlinear zero-dimensional model

The recommended one-state model is

$$
\boxed{
\begin{aligned}
C_{\mathrm{eff}}
\frac{dT_{s,\mathrm{avg}}}{dt}
={}&
\eta_{\mathrm{eff}}
q''*{\mathrm{sol}}A*{\mathrm{irr}}
\\
&-
\dot m c_{p,g}
\varepsilon_{\mathrm{eff}}
\left(
T_{s,\mathrm{avg}}-T_{g,\mathrm{in}}
\right)
\\
&-
K_{\mathrm{lin}}
\left(
T_{s,\mathrm{avg}}-T_\infty
\right)
\\
&-
\chi_{\mathrm{rad}}
\epsilon_{\mathrm{rad}}
\sigma A_{\mathrm{rad}}
\left(
T_{s,\mathrm{avg}}^4-T_{\mathrm{sur}}^4
\right).
\end{aligned}
}
$$

The gas outlet temperature is

$$
\boxed{
T_{g,\mathrm{out}}
=

T_{g,\mathrm{in}}
+
\varepsilon_{\mathrm{eff}}
\left(
T_{s,\mathrm{avg}}-T_{g,\mathrm{in}}
\right),
}
$$

with

$$
\boxed{
\varepsilon_{\mathrm{eff}}
=

\chi_L
\left[
1-
\exp\left(
-\frac{
a_hP_iL
}{
c_{p,g}
}
\dot m^{,n-1}
\right)
\right].
}
$$

The model therefore has one differential state and two principal outputs:

$$
\boxed{
\mathbf{x}=T_{s,\mathrm{avg}},
}
$$

$$
\boxed{
\mathbf{y}
=

\begin{bmatrix}
T_{s,\mathrm{avg}}\\
T_{g,\mathrm{out}}
\end{bmatrix}.
}
$$

## 10. Linear analytical transient solution

For fixed mass flow, constant inlet and ambient temperatures, and linearized losses, define

$$
K_{\mathrm{tot}}
=

K_g+K_{\mathrm{loss}}.
$$

The governing equation becomes

$$
C_{\mathrm{eff}}
\frac{dT_{s,\mathrm{avg}}}{dt}
+
K_{\mathrm{tot}}T_{s,\mathrm{avg}}
=

\eta_{\mathrm{eff}}
A_{\mathrm{irr}}q''*{\mathrm{sol}}
+
K_gT*{g,\mathrm{in}}
+
K_{\mathrm{loss}}T_\infty.
$$

The time constant is

$$
\boxed{
\tau
=

\frac{C_{\mathrm{eff}}}
{K_g+K_{\mathrm{loss}}}.
}
$$

For a constant solar-flux step, the equilibrium temperature is

$$
\boxed{
T_{s,\mathrm{ss}}
=

\frac{
\eta_{\mathrm{eff}}A_{\mathrm{irr}}q''*{\mathrm{sol}}
+
K_gT*{g,\mathrm{in}}
+
K_{\mathrm{loss}}T_\infty
}{
K_g+K_{\mathrm{loss}}
}.
}
$$

The transient solution is

$$
\boxed{
T_{s,\mathrm{avg}}(t)
=

T_{s,\mathrm{ss}}
+
\left[
T_{s,\mathrm{avg}}(0)-T_{s,\mathrm{ss}}
\right]
e^{-t/\tau}.
}
$$

The outlet temperature is

$$
\boxed{
T_{g,\mathrm{out}}(t)
=

T_{g,\mathrm{in}}
+
\varepsilon_{\mathrm{eff}}
\left[
T_{s,\mathrm{avg}}(t)-T_{g,\mathrm{in}}
\right].
}
$$

## 11. Analytical solution for a solar heating ramp

Suppose

$$
q''_{\mathrm{sol}}(t)
=

q''_0+s_qt,
$$

where $s_q$ is the solar-flux ramp rate.

Define

$$
T_{\mathrm{eq},0}
=

\frac{
\eta_{\mathrm{eff}}A_{\mathrm{irr}}q''*0
+
K_gT*{g,\mathrm{in}}
+
K_{\mathrm{loss}}T_\infty
}{
K_{\mathrm{tot}}
},
$$

and

$$
S_T
=

\frac{
\eta_{\mathrm{eff}}A_{\mathrm{irr}}s_q
}{
K_{\mathrm{tot}}
}.
$$

The analytical response is

$$
\boxed{
\begin{aligned}
T_{s,\mathrm{avg}}(t)
={}&
T_{\mathrm{eq},0}
+
S_T
\left[
t-\tau\left(1-e^{-t/\tau}\right)
\right]
\\
&+
\left[
T_{s,\mathrm{avg}}(0)-T_{\mathrm{eq},0}
\right]
e^{-t/\tau}.
\end{aligned}
}
$$

If the receiver initially starts from the equilibrium corresponding to $q''_0$, then

$$
T_{s,\mathrm{avg}}(0)=T_{\mathrm{eq},0},
$$

and the solution simplifies to

$$
\boxed{
T_{s,\mathrm{avg}}(t)
=

T_{\mathrm{eq},0}
+
S_T
\left[
t-\tau\left(1-e^{-t/\tau}\right)
\right].
}
$$

The corresponding outlet-temperature ramp is

$$
\boxed{
T_{g,\mathrm{out}}(t)
=

T_{g,\mathrm{in}}
+
\varepsilon_{\mathrm{eff}}
\left[
T_{s,\mathrm{avg}}(t)-T_{g,\mathrm{in}}
\right].
}
$$

## 12. How the thermocouples identify the coefficients

The local thermocouples are used offline in four ways.

First, the wall thermocouples reconstruct

$$
T_{s,\mathrm{avg}}(t).
$$

Second, the measured axial profile determines the length-shape coefficient $\chi_L$. Conceptually,

$$
\chi_L
=

\frac{
T_{w,\mathrm{eff}}-T_{g,\mathrm{in}}
}{
T_{s,\mathrm{avg}}-T_{g,\mathrm{in}}
}.
$$

Third, the insulation thermocouple helps determine

$$
K_{\mathrm{loss}}
$$

and its thermal lag.

Fourth, the centre and outlet sensors are used to fit

$$
a_h,\qquad n,\qquad \chi_L,
$$

through the outlet relation.

A practical initial fitted parameter vector is

$$
\boxed{
\boldsymbol{\theta}
=

\begin{bmatrix}
C_{\mathrm{eff}},&
\eta_{\mathrm{eff}},&
K_{\mathrm{lin}},&
\epsilon_{\mathrm{rad}}\chi_{\mathrm{rad}},&
a_h,&
n,&
\chi_L
\end{bmatrix}.
}
$$

For a purely linear first model, use

$$
\boxed{
\boldsymbol{\theta}_{\mathrm{linear}}
=

\begin{bmatrix}
C_{\mathrm{eff}},&
\eta_{\mathrm{eff}},&
K_{\mathrm{loss}},&
a_h,&
n,&
\chi_L
\end{bmatrix}.
}
$$

## 13. Parameter identifiability from heating ramps

The early heating slope primarily identifies

$$
\frac{\eta_{\mathrm{eff}}}{C_{\mathrm{eff}}}.
$$

The curvature of the heating transient identifies

$$
\frac{K_g+K_{\mathrm{loss}}}{C_{\mathrm{eff}}}
=

\frac{1}{\tau}.
$$

Changes with flow rate identify

$$
a_h,\qquad n,\qquad \chi_L.
$$

Changes with heat-flux level help distinguish linear loss from radiation.

The parameters $C_{\mathrm{eff}}$ and $\eta_{\mathrm{eff}}$ may remain correlated unless one of them is constrained using:

* known ceramic mass and heat capacity;
* calorimetric absorbed-power information;
* long ramps showing strong curvature;
* several flux levels.

## 14. Control transfer functions

For fixed flow and linearized heat losses, the transfer function from solar flux to average ceramic temperature is

$$
\boxed{
\frac{
T_{s,\mathrm{avg}}(s)
}{
q''_{\mathrm{sol}}(s)
}
=

\frac{
\eta_{\mathrm{eff}}A_{\mathrm{irr}}
}{
C_{\mathrm{eff}}s+K_{\mathrm{tot}}
}.
}
$$

The transfer function from solar flux to gas outlet temperature is

$$
\boxed{
\frac{
T_{g,\mathrm{out}}(s)
}{
q''_{\mathrm{sol}}(s)
}
=

\frac{
\varepsilon_{\mathrm{eff}}
\eta_{\mathrm{eff}}A_{\mathrm{irr}}
}{
C_{\mathrm{eff}}s+K_{\mathrm{tot}}
}.
}
$$

This is a first-order model with a single dominant pole:

$$
s=-\frac{1}{\tau}.
$$

It is therefore well suited to conventional control, gain scheduling, state estimation, or model-predictive control.

## 15. What this model can and cannot reproduce

The model can predict:

$$
T_{s,\mathrm{avg}}(t),
$$

$$
T_{g,\mathrm{out}}(t),
$$

the dominant heating time constant, absorbed-power response, flow-rate effects, and average thermal losses.

The coefficients

$$
\chi_L,\qquad
\chi_{\mathrm{rad}},\qquad
\eta_{\mathrm{eff}}
$$

encode the effect of the nonuniform receiver length.

However, the model cannot independently predict

$$
T_{10},\qquad
T_{58},\qquad
T_{110},
$$

or the exact hot-spot position. Those observations are used to calibrate the reduced coefficients rather than being retained as control-model outputs.

The 0D model is acceptable when the measured local temperature traces can be represented by approximately the same dominant time constant after scaling. If the front, middle and rear traces exhibit substantially different dominant time constants, the later 1D model should provide a measurable improvement.

The recommended baseline is therefore the one-state model

$$
\boxed{
C_{\mathrm{eff}}\dot T_{s,\mathrm{avg}}
=

\dot Q_{\mathrm{abs}}
-

\dot Q_g
-

\dot Q_{\mathrm{loss}},
}
$$

combined with

$$
\boxed{
T_{g,\mathrm{out}}
=

T_{g,\mathrm{in}}
+
\chi_L
\left(1-e^{-NTU_{\mathrm{eff}}}\right)
\left(
T_{s,\mathrm{avg}}-T_{g,\mathrm{in}}
\right).
}
$$

This preserves a genuinely simple control model while using the thermocouple array to incorporate the otherwise missing axial information.


