# 2D v11 theory manual

## Purpose and status

`2D_v11.jl` is a model-form experiment for the structured SiC receiver. It
tests whether two physically motivated changes can explain the measured
temperature profiles:

1. a literature square-channel Graetz heat-transfer law; and
2. unequal ring flows caused by the measured rear support groove.

V11 conserves the corrected MFC mass flow. It does not introduce a bypass,
blocked channels, radial gas mixing or a side-weighted irradiance pattern.

The initial E67--E81 matrix rejected the nominal v11 formulation at the
frozen-v9 parameter point. A subsequent v11-specific staged calibration
restored the axial flow-trend signs but failed held-out absolute
temperatures, transient heating/cooling times, radial ordering and heating
DP1. V11 should therefore be read as a documented diagnostic branch, not as
an accepted fitted model.

## 1. Physical representation

The 19 mm by 19 mm, 137 mm long monolith is represented by an
area-equivalent axisymmetric disc:

$$
A_{\mathrm{rec}}=(19\ \mathrm{mm})^2=361\ \mathrm{mm^2},\qquad
R_{\mathrm{eq}}=\sqrt{A_{\mathrm{rec}}/\pi}=10.7196\ \mathrm{mm}.
$$

The receiver contains 100 square channels of width and hydraulic diameter
1.5 mm. Their total flow area is 225 mm2, giving porosity 0.6233. The
remaining area is the SiC solid. Annular numerical cells stand for groups of
parallel channels.

The thermal domains are:

- the porous SiC receiver;
- the surrounding alumina felt;
- the aluminum casing; and
- the rear alumina tube terminating at a water-cooled flange.

The solid model contains axial and radial conduction, effective internal
radiation, volumetric solar deposition, gas/solid exchange, front radiation,
casing convection/radiation, rear-tube gas exchange and flange conduction.
These mechanisms are inherited from v9 and frozen in the v11 no-refit test.

## 2. Flow reported by the MFC

The corrected MFC reading is treated as standard volumetric flow at the
Aalborg reference conditions

$$
p_{\mathrm{std}}=101.4\ \mathrm{kPa},\qquad
T_{\mathrm{std}}=294.25\ \mathrm{K}.
$$

It is converted once to conserved mass flow:

$$
\dot m_{\mathrm{tot}}
=
Q_{\mathrm{std}}
\frac{p_{\mathrm{std}}}
     {R_{\mathrm{air}}T_{\mathrm{std}}}.
$$

The model does not hold actual volumetric flow or velocity fixed. In each
ring and axial cell,

$$
u_i(z,t)=
\frac{\dot m_i(t)}
     {\rho[T_{g,i}(z,t)]A_{\mathrm{flow},i}}.
$$

Thus heating lowers gas density and raises local velocity for the same
conserved mass flow. Viscosity, conductivity and heat capacity are also
evaluated from local temperatures.

## 3. Local channel heat transfer

### 3.1 Nominal v11 law

For each ring and axial cell,

$$
\begin{aligned}
\mathrm{Re}_i
&=
\frac{\dot m_iD_h}
     {A_{\mathrm{flow},i}\mu(T_g)},\\
\mathrm{Pr}_i
&=
\frac{c_p(T_g)\mu(T_g)}{k(T_g)},\\
\mathrm{Gz}_i
&=
\mathrm{Re}_i\mathrm{Pr}_i\frac{D_h}{z_c},\\
\mathrm{Nu}_i
&=
S_h\left[
\mathrm{Nu}_{fd}
+
\frac{0.104\,\mathrm{Gz}_i}
     {1+0.016\,\mathrm{Gz}_i^{0.8}}
\right]
\left(\frac{T_g}{T_w}\right)^n,\\
h_i
&=
\frac{\mathrm{Nu}_ik(T_w)}{D_h}.
\end{aligned}
$$

The cell-center distance `z_c` removes a singular evaluation at the inlet.
The formula approaches the finite fully developed value `Nu_fd` as Graetz
number tends to zero.

The nominal values are:

| Parameter | Value | Status |
| --- | ---: | --- |
| `Nu_fd` | 3.61 | literature alternative |
| temperature exponent `n` | 0 | nominal |
| heat-transfer scale `S_h` | 1 | fixed, not fitted |
| exposed-front coefficient | 0 | v10 mechanism disabled |

The alternatives `Nu_fd=2.98` and `Nu_fd=3.66, n=0.45` are model-form
sensitivities. They are not separately calibrated.

### 3.2 Exact cell gas update

The gas/solid exchange is evaluated through a cell effectiveness:

$$
\epsilon_{ij}
=
1-\exp\left[
-\frac{h_{ij}P_i\Delta z_j}
       {\dot m_i c_{p,ij}}
\right],
$$

$$
T_{g,\mathrm{out}}
=
T_{g,\mathrm{in}}
+
\epsilon_{ij}(T_{w,ij}-T_{g,\mathrm{in}}).
$$

This avoids a numerically unstable explicit gas-temperature step at large
local NTU. The equal and opposite heat rate is applied to the solid energy
equation.

### 3.3 Historical comparison

For diagnosis only, v11 can reproduce the v9 apparent closure:

$$
\mathrm{Nu}
=
A_{\mathrm{Nu}}
\mathrm{Re}^{1.44}
\mathrm{Pr}^{1/3}
\left(\frac{D_h}{z}\right)^{1/3}.
$$

That law is not called a physical single-channel Nusselt correlation. Its
fitted values are much smaller than the literature Graetz values and it
approaches zero downstream.

## 4. Equal-pressure parallel-channel allocation

All channel groups connect the same inlet and outlet pressure levels.
Therefore v11 solves

$$
\Delta p_i(\dot m_i,T_i)
=
\Delta p_{\mathrm{common}},
\qquad
\sum_i\dot m_i=\dot m_{\mathrm{tot}}.
$$

The solver begins with area-proportional flows, marches every ring, updates
flows from their current hydraulic conductances, relaxes the update and
renormalizes the result to exact total mass flow. The nominal relative
pressure-equality tolerance is `1e-5`.

There is no radial mass transfer after the flow split: each channel group
remains axially isolated. Radial interaction occurs only through the solid
thermal model.

## 5. Pressure-drop model and DP1

### 5.1 Channel term

For a square laminar duct, each cell uses

$$
\Delta p_{\mathrm{fric},ij}
=
\chi_h
\frac{56.91}{2}
\frac{\mu_{ij}\Delta z_j u_{ij}}{D_h^2}.
$$

The shared multiplier `chi_h` covers the effective resistance not resolved
by the continuum hydraulic discretization.

### 5.2 Measured rear groove

The rear support leaves approximately 13 mm diameter completely free.
The rest of the monolith is not closed, but its outlet can be obstructed.
For every annulus, v11 calculates the exact fraction `I_i` lying outside
radius 6.5 mm and adds

$$
\Delta p_{\mathrm{groove},i}
=
I_i K_{\mathrm{groove}}
\frac{\rho_{i,\mathrm{out}}u_{i,\mathrm{out}}^2}{2}.
$$

The completely free area is 132.73 mm2, or 36.77% of the receiver face.
The remaining 63.23% is groove-exposed area, not blocked area. Every channel
retains its flow cross-section and heat-transfer perimeter.

### 5.3 DP1 observation

The flush wall tap was inside the water jacket and referenced to atmosphere.
The model observation is

$$
\mathrm{DP1}_{\mathrm{model}}
=
b_{\mathrm{DP1}}
+
\frac{\Delta p_{\mathrm{common}}}{100},
$$

where the division converts pascals to mbar. The sensor offset is allowed
because the cold records show a nonzero zero-flow intercept. No bypass
branch is inferred.

## 6. Solar and solid model assumptions

The v9-fitted optical and solid values are frozen for this model-form test.
Important assumptions are:

- the physical square area is conserved in the equivalent disc;
- the imposed beam is radially Gaussian, not side weighted;
- axial absorption follows finite Beer-Lambert attenuation;
- unabsorbed power can transmit out of the receiver;
- front/rim deposition and external losses remain those of v9;
- gas channels do not exchange mass radially;
- the axisymmetric perimeter temperature is an annular average; and
- side thermocouples are sampled at the perimeter at the verified axial
  positions 5, 58 and 107 mm.

The shallow mounting dips of the side thermocouples are not explicitly
represented. Nor is the four-sided/corner geometry of the square monolith or
local solid contact with the support groove.

## 7. Parameters and what was fitted

### Frozen from v9

Optical delivery and attenuation, solid density and conductivity scales,
radiative extinction, insulation/casing properties, contact resistance,
external losses, rear adaptor/flange parameters and MFC mass-flow scale are
not refitted in v11.

This qualification is essential. The v9-fitted parameters that remain active
in nominal v11 are `k_scale_r`, `k_scale_z`, the three irradiance-group
delivery scales, radiative extinction and optical extinction. They were
estimated together with the v9 apparent heat-transfer closure and may partly
compensate for it. The v9 `A_Nu` and `B_Re` do not control the nominal
Graetz branch, and its heating-only excess pressure coefficient is disabled.
Thus v11 does not inherit every v9 closure error, but its absolute thermal
comparison is still confounded by the shared solid/optical fit.

### Fixed literature tests

`Nu_fd`, the Graetz coefficients and the temperature-correction exponent are
fixed alternatives. `S_h=1`; no heat-transfer scale was fitted because the
profile shapes failed.

### Cold-DP1 parameters

Nine t0 cold cases give:

| Form | `b_DP1` (mbar) | `chi_h` | `K_groove` | RMSE (mbar) |
| --- | ---: | ---: | ---: | ---: |
| no groove | -0.61383 | 1.94402 | 0 | 0.02541 |
| 13 mm groove | -0.54284 | 0.97020 | 184.16 | 0.02452 |

The groove model has worse AICc (`-55.95` versus `-60.11`). Moreover,
profiling `K_groove` shows a continuous tradeoff with `chi_h`. The numerical
value 184.16 therefore identifies one possible decomposition of combined
path resistance; it is not an independently measured groove coefficient.

## 8. Validation results

The full no-refit matrix covers E67--E81.

### 8.1 Axial profiles

The nominal Graetz/groove model gives an axial `T12-T8` RMSE of 176.14 K,
compared with 143.08 K for the historical v9 comparison and the predeclared
v11 target below 33.2 K.

It also predicts the wrong response to flow:

| Irradiance | Measured slope | V11 slope |
| ---: | ---: | ---: |
| 456 kW/m2 | +17.34 K/(L/min) | -14.59 K/(L/min) |
| 304 kW/m2 | +10.48 K/(L/min) | -9.40 K/(L/min) |
| 256 kW/m2 | +6.90 K/(L/min) | -6.29 K/(L/min) |

The 2.98, 3.61 and temperature-corrected 3.66 variants behave almost
identically. This indicates excessive/near-saturated gas-solid coupling at
the literature Nu scale in the present macro model, not a fine choice of
fully developed asymptote.

### 8.2 Radial profiles

The groove creates a mean core-to-edge mass-flux ratio of 2.18 and the ratio
increases with total flow, as expected for a quadratic restriction.
Nevertheless:

| Observable | V11 mean | RMSE |
| --- | ---: | ---: |
| `T12-T9` | -0.25 K | 24.83 K |
| `T11-T10` | -0.70 K | 37.97 K |

The measured offsets are positive and typically tens of kelvin. Strong
hydraulic maldistribution therefore does not survive as a comparable solid
temperature contrast in the current axisymmetric thermal representation.

The complete 15-case branch was also repeated with cold-fitted 12 and 14 mm
free diameters. Their mean core/edge mass-flux ratios are 3.33 and 2.26, yet
the mid-depth modeled radial means remain +0.19 K and -0.25 K. The deep
means remain about -0.63 K and -0.70 K. The missing radial observation is
therefore insensitive to this measured-geometry uncertainty at the inherited
thermal point.

### 8.3 T3 and pressure

Nominal v11 T3 MAE is 99.91 K versus 78.38 K in the historical comparison.
Heating DP1 RMSE is 0.570 mbar with -0.373 mbar bias, versus 0.161 mbar for
the historical v9 branch. The cold-fitted groove does not replace the old
hot-path discrepancy successfully.

## 9. Numerical verification

- v11 smoke tests: 27/27 passed;
- representative nominal-mesh tests: 4/4 passed;
- v9 and v10 regression suites: 34/34 and 21/21 passed;
- exact total-flow normalization is maintained;
- pressure-equality residual is below `8.7e-6` over the full matrix;
- uniform-state energy residual is below `1e-8 W`;
- the Graetz law approaches the requested downstream asymptote; and
- nominal-mesh E67/E72/E77 sensor changes are at most 1.26 K.

The frozen-parameter incompatibility is consequently physical/observational,
not a solver, conservation or reduced-mesh failure. Whether an independently
recalibrated Graetz branch remains incompatible is a separate question.

## 10. Expert interpretation and next model boundary

V11 answers two questions:

1. **Can the measured groove create the proposed flow feedback?** Yes. It
   gives less edge flow at high total flow and weaker maldistribution as
   total flow falls.
2. **Does that feedback explain the measured radial thermocouples in the
   current model?** No. Even a very large flow contrast produces almost no
   modeled side-to-core solid-temperature contrast.

The next bounded test should retain mass conservation, local
temperature-dependent velocity, equal-pressure flow allocation and the
measured groove geometry, but change the observation/solid geometry rather
than the irradiance map. A square-sector or two-solid-population model should
separate internal SiC radial conduction, local holder/groove contact and the
side-thermocouple contact representation.

Before a final transport-form decision, run a v11-specific sensitivity and
identifiability audit, then a constrained refit on the original nine
training heating cases with the six held-out cases untouched. Candidate
active parameters are `S_h`, radial and axial solid conductivity scales,
optical extinction and the three group power scales. They must not all be
released blindly; rank and profile-shape sensitivities should define the
smallest identifiable subset. Cooling remains validation.

The first one-at-a-time inheritance audit has now been completed. Its main
result is that optical extinction is strongly closure-dependent. Changing
the v9-fitted `beta_opt` from 21.29 to 110 1/m reduces the complete
15-case axial RMSE from 176.14 to 58.17 K and changes all three axial
flow-slope signs from negative to positive. The slopes remain much too small
and weakly correlated, so this is not an accepted parameter value. It is
decisive evidence that v11 must be recalibrated before the Graetz form itself
is rejected.

The same audit does not rescue the radial profile. Broad variations of
radial/axial conductivity, optical/radiative extinction, power, or `S_h`
leave the modeled side offsets near zero or negative. At `beta_opt=110`,
mid and deep radial RMSE are 25.61 and 38.19 K. The current rear-groove
redistribution is therefore still insufficient within the axisymmetric
solid/observation representation.

The reporting constraints carried into the staged calibration were:

- do not report v11's literature Nu as an experimentally validated effective
  coefficient;
- do not report `K_groove=184.16` as uniquely identified;
- do not fit `S_h`, optics and radial conductivity together without a new
  rank/identifiability result;
- do not restore the v10 superlinear front coefficient as established
  physics; and
- keep v9 as the better empirical reference while treating its internal Nu
  as an apparent whole-model closure.

## 11. Staged v11-specific calibration

The subsequent calibration used only the nine established heating training
cases. Six heating cases and all three cooling cases remained untouched.

Stage 1 fit `beta_opt`, `S_h` and `k_scale_z` to axial/profile differences
and flow ordering. Stage 2 fit the three group power scales to radially
averaged absolute temperatures. Side/core radial differences were excluded
because the preceding tests identified them as a structural observation
mismatch.

The selected vector is:

```text
beta_opt = 218.8 1/m
S_h = 1.142
k_scale_z = 0.2
scale_456 = 1.25
scale_304 = 1.051748
scale_256 = 0.773518
```

The local six-parameter residual Jacobian is rank 6/6 with condition number
9.37, but `k_scale_z` remains on its lower bound.

This confirms that the frozen-v9 axial failure was not a fair final test.
All three axial flow-slope signs become positive and all-heating axial RMSE
improves to 44.14 K. The slopes remain too weak:

```text
model:    8.02, 4.62, 4.24 K/(L/min)
observed: 17.34, 10.48, 6.90 K/(L/min)
```

The complete model nevertheless fails validation. Held-out heating mean
sensor RMSE is 98.63 K versus 63.92 K for v9; cooling RMSE is 70.99 K versus
34.71 K. Heating and cooling are much too fast, with t90 errors of order
1,200 s. Radial ordering is wrong in all 15 heating cases and held-out
heating DP1 RMSE is 0.366 mbar.

Therefore the order-unity `S_h` cannot be extracted as a validated effective
heat-transfer coefficient. Its numerical identifiability is conditional on
a state/observation model that fails transient, radial and pressure
validation. The next physics must address thermal storage/sensor response
and the square side-wall/support contact before another transport
calibration.
