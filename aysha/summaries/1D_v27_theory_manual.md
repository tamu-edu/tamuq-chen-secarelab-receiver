# Theory Manual: 1D_v27 Energy-Accounting Two-Zone Receiver Model

## Manuscript-Level Model Description, Parameterization, and Scientific Scope

This document describes the mathematical formulation, fitted parameterization,
and scientific interpretation of `1D_v27.jl`, a one-dimensional reduced-order
thermal model for a structured monolithic solar receiver with square channels.
The model is part of a staged 1D development sequence whose broader purpose is
to obtain and validate effective macroscopic heat-transfer closures for the
receiver, including convective gas-solid exchange, effective solid/perimeter
redistribution, radiative/volumetric source deposition, and rear/adaptor/flange
thermal losses.

Version v27 should be interpreted as an energy-accounting macro model with a
two-zone receiver topology: a core channel matrix and a perimeter/housing
branch. It keeps the corrected incident-power accounting from v25 and the
useful cooling-side corrections introduced in v26, while repairing the most
severe v26 failure mode by separating the rear-core heat-removal strength used
during irradiated heating from the stronger effective heat-removal strength
needed during zero-irradiance cooling. It also restores the T3 comparison point
to 140 mm and introduces a modest fitted gas/tube-wall measurement blend for
T3.

The v27 result is more defensible than v26 as a fitted predictive structure,
but it is not a final coefficient-identification model. Its direct distributed
rear-core sink and fitted T3 wall fraction are effective closures. They should
be credited as diagnostic surrogate terms, not as independently validated
physical coefficients.

---

## 1. Physical System and Modeling Objective

The experimental receiver is a silicon carbide monolithic structure with a
square frontal aperture and a 10 x 10 array of square air channels. Concentrated
solar flux enters through the receiver front face. Air flows axially through the
channels, exchanges heat with the solid matrix, and leaves through a downstream
rear tube/adaptor assembly. The receiver is surrounded by insulation and a metal
cavity/housing that exchanges heat with the laboratory environment.

The 1D model family seeks a reduced-order continuum representation that can be
calibrated against transient thermocouple and gas-temperature measurements.
The scientific target is not merely curve fitting. The intended deliverable is
a set of effective transport closures that can bridge detailed single-channel
physics and an entire-converter macro model.

The v27 model specifically addresses three observations from earlier versions:

1. A single rear heat-removal strength overcooled irradiated heating cases when
   tuned to match cooling.
2. T3 behaved neither like a pure receiver-exit gas temperature nor like a pure
   wall measurement.
3. Some axial redistribution was useful once cooling monotonicity was protected,
   but broad fitted thermal patches could hide missing physics.

---

## 2. Dynamic States and Domain Topology

The receiver is discretized into `N` axial finite volumes of size

```text
dx = L / N
```

where `L = 137 mm`. The downstream rear tube is discretized into `N_rear`
finite volumes of size

```text
dx_rear = L_rear / N_rear
```

The dynamic state vector is

```text
u(t) = [T_core,1 ... T_core,N,
        T_perim,1 ... T_perim,N,
        T_tube,1 ... T_tube,Nrear,
        T_cavity]
```

where:

- `T_core,i` is the temperature of the central channel-matrix/core zone.
- `T_perim,i` is the temperature of the perimeter/housing thermal branch.
- `T_tube,j` is the downstream alumina rear-tube/adaptor wall temperature.
- `T_cavity` is a lumped cavity/housing temperature compared to T2.

The gas phase is reconstructed algebraically at each ODE evaluation by a
quasi-steady plug-flow march. This is appropriate because gas residence time is
much shorter than the solid/cavity transient timescale.

---

## 3. Geometry and Fixed Physical Constants

The inherited receiver geometry is:

| Quantity | Code Symbol | Value |
| :--- | :--- | ---: |
| Receiver frontal width | `w_t` | 19.0 mm |
| Frontal area | `A_frt = w_t^2` | 3.61e-4 m2 |
| Receiver length | `L` | 137.0 mm |
| Number of square channels | `n_chnl` | 100 |
| Channel width | `w_chnl` | 1.5 mm |
| Hydraulic diameter | `Dh` | 1.5 mm |
| Total flow area | `A_flow` | `n_chnl w_chnl^2` |
| Gas-solid exchange perimeter | `P_exchange` | `4 w_chnl n_chnl` |

The v27 cavity and rear domain use:

| Quantity | Code Symbol | Value |
| :--- | :--- | ---: |
| Cavity outer radius | `CAVITY_OUTER_RADIUS_v27` | 75.0 mm |
| Insulation outer radius | `INSULATION_OUTER_RADIUS_v27` | 57.0 mm |
| Cavity length | `CAVITY_LENGTH_v27` | 165.0 mm |
| Adaptor diameter | `ADAPTOR_DIAMETER_v27` | 77.6 mm |
| Adaptor length | `ADAPTOR_LENGTH_v27` | 57.0 mm |
| Rear tube length | `REAR_TUBE_LENGTH_v27` | 150.0 mm |
| Rear tube gas radius | `REAR_TUBE_GAS_RADIUS_v27` | 6.5 mm |
| Rear tube wall thickness | `REAR_TUBE_WALL_THICKNESS_v27` | 1.5 mm |
| Rear-tube cavity-coupled length | `REAR_TUBE_CAVITY_LENGTH_v27` | 46.0 mm |
| T3 comparison position | `T3_SAMPLE_POSITION_v27` | 140.0 mm |

Important fixed thermal constants include:

```text
epsilon_front = 0.95
eta_abs = 1.0
h_front = 10 W/m2/K
Nu_floor = 3.61
T_water_flange = 25 C
T_cooling_room = 17 C
```

The cooling-room temperature is imposed only for cooling simulations, where
both inlet and ambient histories are set to 17 C.

---

## 4. Solar Power Accounting and Source Partition

### 4.1 Total Participating Receiver Power

For irradiance `I`, the modeled participating receiver power is

```text
Q_abs(I) = eta_abs scale(I) I A_frt
```

with

```text
eta_abs = 1.0
```

and three fitted irradiance-level scales:

```text
scale_456 = p7
scale_304 = p8
scale_256 = p9
```

For the fitted v27 parameter set:

```text
scale_456 = 2.0025
scale_304 = 1.9966
scale_256 = 0.9520
```

These scale factors should be interpreted as effective participating flux
multipliers, not as pure absorptivity. They include the relation between
aperture-average reported irradiance and local receiver-average flux, together
with remaining calibration and model-form uncertainty.

### 4.2 Shape-Only Transverse Flux Prior

The transverse flux shape is represented by a Gaussian:

```text
F(x,y) = exp[-(x^2 + y^2)/(2 sigma^2)]
```

with

```text
sigma = 30 mm
```

The model integrates this shape over the circular aperture and the square
receiver face. The resulting fixed fractions are:

```text
flux_receiver_fraction = 0.0692
flux_spillover_fraction = 0.9308
```

The fitted spillover capture factor `p14` determines what fraction of the
off-receiver transverse flux shape participates as perimeter/housing heating:

```text
f_perim = clamp(p14 flux_spillover_fraction, 0, 0.80)
```

For the fitted v27 result:

```text
p14 = 0.5251
f_perim = 0.4888
```

Thus nearly half of the modeled participating power is routed to the perimeter
branch:

```text
Q_core = (1 - f_perim) Q_abs
Q_perim = f_perim Q_abs
```

### 4.3 Axial Source Weights

The core source uses the inherited Beer-Lambert-plus-front-deposition function:

```text
w_core,i proportional to exp(-beta_opt z_left) - exp(-beta_opt z_right)
```

followed by a front-cell deposition fraction:

```text
w_core,1 += front_dep
```

In v27,

```text
front_dep = 1.0
beta_opt = 184.7 1/m
```

so all core-source power is placed in the first axial cell. This is a major
source-distribution limitation.

The perimeter source uses a separate exponential axial attenuation:

```text
w_perim,i proportional to exp(-beta_perim z_left) - exp(-beta_perim z_right)
```

with fitted

```text
beta_perim = 3.629 1/m
```

This means the perimeter source is much more axially distributed than the core
source.

---

## 5. Receiver Gas-Solid Heat Exchange

The receiver gas is modeled as quasi-steady plug flow. For each core cell:

```text
T_film,i = 0.5 (T_core,i + T_g,in,i)
```

Air properties are evaluated at the film temperature. The mass flow rate is
computed from the measured standard volumetric flow:

```text
mdot = rho_std Q_lpm / 60000
rho_std = 101325 / (287.05 * 294.25)
```

The Reynolds and Prandtl numbers are:

```text
Re = mdot Dh / (A_flow mu)
Pr = cp mu / k_f
```

### 5.1 Active Flow Fraction

The general active-flow form is:

```text
phi_act = clamp(phi_0 (Re/Re_ref)^m_rec, 0.1, 1.0)
```

with

```text
Re_ref = 50
```

For fitted v27:

```text
phi_0 = 1.0
m_rec = 0.0
```

so all measured flow is treated as active channel flow. No bypass or partial
participation is active in the fitted v27 baseline.

### 5.2 Nusselt Number Closure

The local receiver Nusselt number is:

```text
Nu_i =
max[Nu_floor,
    A_Nu Re_i^B_Re Pr_i^(1/3) (Dh / max(z_i, Dh/2))^(1/3)]
```

with

```text
Nu_floor = 3.61
0 <= B_Re <= 0.6
```

For fitted v27:

```text
A_Nu = 4.9176
B_Re = 0.5128
```

The old upper bound on `A_Nu` was 5, so v27 places the Nusselt prefactor close
to its upper limit. This should be interpreted cautiously: it indicates strong
effective gas-solid exchange within the v27 energy topology, not a validated
standalone square-channel Nusselt correlation.

The heat-transfer coefficient is

```text
h_i = Nu_i k_f / Dh
```

and the cell NTU/effectiveness update is:

```text
UA_i = h_i P_exchange dx
epsilon_i = 1 - exp[-UA_i/(mdot_active cp)]
T_g,out,i = T_g,in,i + epsilon_i (T_core,i - T_g,in,i)
Q_gas,i = mdot_active cp (T_g,out,i - T_g,in,i)
```

---

## 6. Core and Perimeter Energy Balances

### 6.1 Core Receiver Balance

For cell `i`, the core energy balance is:

```text
C_core,i dT_core,i/dt =
    Q_core w_core,i
  - Q_gas,i
  + G_cp dx (T_perim,i - T_core,i)
  - Q_rear_core_sink,i
  + axial core conduction
  - receiver-to-tube core heat at i = N
```

The core heat capacity is:

```text
C_core,i = rho_s cp_s(T_core,i) A_solid dx
```

The fitted v27 effective axial core conduction term is:

```text
Q_cond,core,i+1/2 =
    k_core_axial_scale k_s,face A_solid (T_core,i - T_core,i+1) / dx
```

with

```text
k_core_axial_scale = 0.0361
```

The value is small but nonzero. Under the cooling monotonicity guard, the fit
uses a limited amount of axial core redistribution.

### 6.2 Perimeter Branch Balance

The perimeter branch obeys:

```text
C_perim,i dT_perim,i/dt =
    Q_perim w_perim,i
  - G_cp dx (T_perim,i - T_core,i)
  - G_cavity dx (T_perim,i - T_cavity)
  - Q_rear_perim_sink,i
  + axial perimeter conduction
  - front loss at i = 1
  - receiver-to-tube perimeter heat at i = N
```

The total perimeter heat capacity is fitted and distributed uniformly:

```text
C_perim,i = C_perim_eff / N
```

For fitted v27:

```text
C_perim_eff = 222.1 J/K
C_core_eff = 72.5 J/K
C_participating = 294.7 J/K
measured reference = 301 J/K
```

The fitted participating thermal inventory is therefore close to the measured
assembly capacitance target.

Perimeter axial conduction uses:

```text
k_perim(T) = k_perim_ref (T/900 K)^3
```

with

```text
k_perim_ref = 7.62 W/m/K
```

This is an effective perimeter/housing axial redistribution term and should not
be reported as a bulk material property without qualification.

### 6.3 Core-Perimeter Coupling

The radial two-zone coupling is:

```text
Q_cp,i = G_core_perim dx (T_perim,i - T_core,i)
```

with fitted:

```text
G_core_perim = 17.86 W/m/K
```

This term allows the perimeter and core branches to exchange heat at each axial
station. It is a macro-scale conductance representing unresolved lateral heat
paths between the central channel matrix and surrounding perimeter/housing
domain.

---

## 7. Direct Rear Sink in v27

The defining v27 repair relative to v26 is the split rear-core sink:

```text
G_rear_core_eff =
    G_rear_core_heat, if irradiance > 0
    G_rear_core_cool, if irradiance = 0
```

with fitted:

```text
G_rear_core_heat = 3.330 W/m/K
G_rear_core_cool = 8.866 W/m/K
```

The rear sink activates only downstream of the T8 station using the axial
weight:

```text
w_rear(z) = 0,                         z <= z_T8
w_rear(z) = [(z - z_T8)/(L - z_T8)]^m, z > z_T8
```

where

```text
m = rear_sink_shape = 0.9998
```

The distributed sink terms are:

```text
Q_rear_core_sink,i =
    G_rear_core_eff dx w_rear(z_i) (T_core,i - T_water)

Q_rear_perim_sink,i =
    G_rear_perim dx w_rear(z_i) (T_perim,i - T_water)
```

For fitted v27:

```text
G_rear_perim = 0
```

so the direct rear sink acts only on the core branch.

This split rear sink is scientifically useful because it confirms the v26
overcooling diagnosis: cooling data require stronger downstream removal than
irradiated heating can tolerate. However, the form is still empirical. The
instantaneous irradiance branch and distributed core-to-water path should not
be interpreted as literal geometry-resolved conduction. It is an effective
closure for unresolved rear/adaptor/flange storage and heat removal.

---

## 8. Explicit Rear Tube and Flange Path

In addition to the distributed rear sink, v27 retains an explicit rear tube
domain. At the receiver exit:

```text
Q_rec-tube,core =
    G_rec-tube f_core_tube (T_core,N - T_tube,1)

Q_rec-tube,perim =
    G_rec-tube (1 - f_core_tube) (T_perim,N - T_tube,1)
```

with fitted:

```text
f_core_tube = 0.9993
```

so nearly all receiver-to-tube solid heat is routed through the core exit.

The rear tube has axial conduction:

```text
Q_cond,tube,j+1/2 =
    k_tube,face A_tube,face (T_tube,j - T_tube,j+1) / dx_rear
```

It exchanges heat with rear gas by an NTU update using standard internal-flow
correlations:

```text
Nu_tube = 3.66,                         Re < 2300
Nu_tube = 0.023 Re^0.8 Pr^0.4,          Re >= 2300
```

It also exchanges heat with the cavity and downstream water-flange path:

```text
Q_tube-cavity,j = G_tube-cavity,j dx_rear (T_tube,j - T_cavity)
Q_tube-flange,j = G_tube-flange,j dx_rear (T_tube,j - T_water)
```

The fitted flange multiplier is:

```text
flange_scale = 0.1014
```

which is near the lower bound. This implies that the explicit flange path is
kept weak while the distributed rear-core sink carries much of the fitted rear
removal. That split is one of the main reasons v27 is not a final physical
topology.

---

## 9. Cavity and Front Losses

The perimeter branch transfers heat to a lumped cavity:

```text
Q_perim-cavity,i =
    G_receiver-cavity dx (T_perim,i - T_cavity)
```

The cavity loses heat to ambient by natural convection and radiation:

```text
Q_cavity-ambient =
    h_nat A_cavity (T_cavity - T_amb)
  + epsilon_metal sigma A_cavity (T_cavity^4 - T_amb^4)
```

The front perimeter cell loses heat to ambient through:

```text
Q_front =
    h_front A_frt (T_perim,1 - T_amb)
  + epsilon_front sigma A_frt (T_perim,1^4 - T_amb^4)
```

with:

```text
h_front = 10 W/m2/K
epsilon_front = 0.95
```

T2 is compared directly to the lumped cavity temperature:

```text
T2_model = T_cavity
```

The fitted v27 model controls T2 well relative to receiver/gas channels,
supporting the use of a single lumped cavity state at this modeling stage.

---

## 10. Experimental Output Mapping

The v27 measurement map is:

| Experimental channel | v27 model output |
| :--- | :--- |
| T8 | perimeter temperature at 11 mm |
| T12 | perimeter temperature at 58 mm |
| T11 | perimeter temperature at 107 mm |
| T9 | core temperature at 58 mm |
| T10 | core temperature at 107 mm |
| T3 | mixed gas/tube-wall output at 140 mm |
| T2 | lumped cavity temperature |

The T3 model is:

```text
T3_model =
    (1 - f_T3_wall) T_gas(140 mm)
  + f_T3_wall T_tube_wall(140 mm)
```

with fitted:

```text
f_T3_wall = 0.1757
```

The fitted value is modest, suggesting that T3 is mostly gas-like in v27, but
with some wall/sensor environmental influence. This is a measurement model, not
a heat-transfer mechanism.

For cooling simulations, the first output row is aligned to the measured
initial values:

```text
outputs[1, 1:7] = experiment[1, 1:7]
```

This is useful for evaluating cooling transient shape, but it should be
reported because the first cooling residual is removed by construction.

---

## 11. Calibration Objective

The calibration objective combines heating and cooling cases. For each sensor,
v27 uses:

```text
L_signal =
    normalized_mse(model, experiment)
  + 0.25 normalized_slope_mse(model, experiment)
  + 0.15 [(t90_model - t90_experiment)/duration]^2
```

Heating cases additionally include an axial-ordering penalty based on final
temperature offsets:

```text
T12 - T8
T11 - T8
T12 - T9
T11 - T10
```

Cooling cases add a large penalty against artificial positive upturns in:

```text
T11, T10, T3
```

with:

```text
w_upturn = 10000
```

The participating heat capacity is regularized toward the measured reference:

```text
C_ref = 301 J/K
w_C = 0.10
```

The fitted value:

```text
C_participating = 294.7 J/K
```

is close to the target.

---

## 12. Fitted v27 Parameter Set

The full v27 runner completed with:

```text
objective = 0.17506355554436148
return_code = MaxTime
```

The fitted parameters are:

| Parameter | Value | Interpretation |
| :--- | ---: | :--- |
| `A_Nu` | 4.9176 | near old upper bound |
| `B_Re` | 0.5128 | fitted Reynolds exponent |
| `C_Pr` | 0.3333 | fixed Prandtl exponent |
| `phi_0` | 1.0 | fixed full active flow |
| `m_rec` | 0.0 | fixed no recruitment trend |
| `front_dep` | 1.0 | fixed all core source in first cell |
| `scale_456` | 2.0025 | effective 456 kW/m2 power scale |
| `scale_304` | 1.9966 | effective 304 kW/m2 power scale |
| `scale_256` | 0.9520 | effective 256 kW/m2 power scale |
| `G_core_perim` | 17.86 W/m/K | fitted two-zone coupling |
| `C_perim_eff` | 222.13 J/K | fitted perimeter capacity |
| `k_perim_ref` | 7.62 W/m/K | effective perimeter axial conductivity |
| `beta_opt` | 184.7 1/m | fixed |
| `spill_capture` | 0.5251 | fitted spillover capture |
| `beta_perim` | 3.629 1/m | fitted perimeter source attenuation |
| `f_core_tube` | 0.9993 | nearly all receiver-tube heat from core |
| `flange_scale` | 0.1014 | near lower bound |
| `G_rear_core_heat` | 3.330 W/m/K | heating distributed rear-core sink |
| `G_rear_perim` | 0.0 W/m/K | perimeter rear sink inactive |
| `rear_sink_shape` | 0.9998 | nearly linear rear-sink growth |
| `k_core_axial_scale` | 0.0361 | small fitted core axial redistribution |
| `G_rear_core_cool` | 8.866 W/m/K | cooling distributed rear-core sink |
| `f_T3_wall` | 0.1757 | T3 wall/environment measurement fraction |

---

## 13. Performance and Diagnostic Interpretation

The full v27 result is a major recovery from the v26 collapse. Representative
mean signed steady errors reported in the model journal are:

```text
Heating:
  T8  =  -49 K
  T12 =  -84 K
  T11 =  -52 K
  T9  =  -79 K
  T10 =  -42 K
  T3  = -117 K
  T2  =   -3 K

Cooling:
  T8  = -11 K
  T12 = -16 K
  T11 = -19 K
  T9  = -19 K
  T10 = -25 K
  T3  = +12 K
  T2  =  -0 K
```

This shows that v27 is much more balanced than v26, especially for cooling and
T2. However, heating remains systematically low for most receiver and gas
outputs, especially T3.

Power-path diagnostics show that the rear-core sink remains a large heating
loss pathway. Representative heating cases include:

```text
E67, 456 kW/m2:
  Q_abs,total = 329.6 W
  Q_receiver_gas = 127 W
  Q_rear_core_sink = 110 W

E80, 256 kW/m2:
  Q_abs,total = 88.0 W
  Q_receiver_gas = 15 W
  Q_rear_core_sink = 41 W
```

Therefore the model should not be described as simply power-limited. A large
part of the fitted power is diverted through the effective rear sink.

---

## 14. Scientific Credit and Limitations

v27 deserves credit for:

- correcting the v26 overcooling failure by separating heating and cooling rear
  removal strengths,
- keeping a conservative finite-volume energy accounting structure,
- preserving a physical rear tube/cavity/flange domain,
- comparing core and perimeter thermocouples with a two-zone macro topology,
- maintaining participating heat capacity near the measured assembly value,
- exposing that T3 is not fully represented by pure receiver-exit gas under the
  current rear-domain topology.

However, v27 has important limitations:

1. The direct distributed rear-core sink is an effective closure, not a literal
   physical conductance.
2. The heating/cooling rear-sink split is a hard irradiance branch.
3. T3 wall mixing is a fitted measurement correction.
4. The core source is still front-cell locked by `front_dep = 1`.
5. `A_Nu` is near its old upper bound and should not be reported as a validated
   Nusselt correlation.
6. The explicit flange path is near its lower bound while the distributed sink
   carries much of the rear removal.

The correct manuscript language is therefore:

```text
v27 is an effective energy-accounting macro model that resolves the v26
rear-overcooling pathology and gives a useful calibrated baseline. Its fitted
rear sink and T3 measurement blend are diagnostic closures, not final physical
coefficients.
```

---

## 15. Recommended Successor Direction

The v27 review and the v28 negative test point to the same next step. The
distributed direct rear sink should not be accepted as final physics, but simply
removing it makes the model structurally under-complete. A better successor
should replace the direct sink with a physically bounded rear/adaptor reservoir:

```text
receiver exit solid -> rear/adaptor reservoir -> explicit rear tube/flange/cavity
```

The rear reservoir should have:

- a heat capacity bounded by measured alumina/adaptor/flange participating mass,
- conductances constrained by contact areas and plausible thermal-contact
  ranges,
- no direct distributed core-to-water sink along the receiver length,
- T3 initially retained as gas temperature until rear topology is represented.

This would preserve the main lesson of v27 while improving scientific
creditability: rear cooling must be stronger than irradiated rear loss, but the
mechanism should be represented by an interpretable rear thermal inventory
rather than a phase-switched distributed sink.

---

## 16. Conclusion

The 1D_v27 model is a meaningful intermediate receiver model. It combines
two-zone receiver energy accounting, fitted irradiance-level participating
power scales, flux-shape-guided perimeter source partitioning, bounded
developing-flow convection, explicit rear tube/cavity dynamics, split rear-core
losses, and a modest T3 gas/wall measurement blend.

Its key scientific contribution is diagnostic: the data require stronger
effective rear removal during cooling than during irradiated heating, and T3
cannot yet be explained by the current gas/rear-tube representation alone.
However, because the rear removal is still implemented as a distributed
phase-dependent sink, v27 should be treated as a calibrated macro-model
baseline rather than a final coefficient-identification model.

