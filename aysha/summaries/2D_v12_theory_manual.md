# 2D v12 theory manual

## Purpose and status

`2D_v12.jl` tests whether the v11 transient and radial failures arise from
omitted installed-hardware capacitance and incorrect thermocouple observation
physics.  It retains v11 mass-flow conservation, local gas density and
velocity, equal-pressure channel allocation, Graetz heat transfer and the
measured rear groove.

The staged cooling, heating and held-out tests are complete.  V12 improves
cooling and identifies a credible weakly participating adaptor, but fails
held-out heating, axial profiles and every radial-ordering test.  It is
rejected for coefficient extraction and retained as a documented
thermal-inventory model-form result.

## 1. Receiver and enclosure geometry

The receiver is a 19 by 19 mm square, 137 mm long SiC monolith represented by
an area-equivalent axisymmetric domain.  The photographed assembly establishes:

| Quantity | Value |
| --- | ---: |
| Aluminum enclosure outside diameter | 150 mm |
| Aluminum radial thickness | 18 mm |
| Enclosure inner radius | 57 mm |
| Enclosure internal length | 165 mm |
| Rear backplate thickness | 18 mm |
| Alumina adaptor outside diameter | 77.6 mm |
| Adaptor length | 57 mm |
| Adaptor overlap around receiver | 29 mm |
| Adaptor overlap around tube | 28 mm |
| Tube bore diameter | 13 mm |

The felt fills the remaining cavity.  The tube has no direct felt contact.
Its first portion is enclosed by the alumina adaptor and the rest lies outside
the cavity inside the continuously water-cooled flange.

## 2. Hardware capacitances

The solid adaptor volume is

$$
V_a =
\left(\pi R_a^2-A_{\rm rec}\right)L_{\rm overlap}
+
\left(\pi R_a^2-\pi R_{\rm tube}^2\right)L_{\rm tube}.
$$

For $R_a=38.8$ mm, $L_{\rm overlap}=29$ mm and
$L_{\rm tube}=28$ mm,

$$
V_a=255.395\ {\rm mL},\qquad
m_a=\rho_{\rm alumina}V_a=0.99604\ {\rm kg}.
$$

The receiver, felt, adaptor, aluminum sleeve/backplate and rear tube keep
their own heat capacities.  The experimentally observed
$C_{\rm eff}=301\pm23$ J/K is a participating modal capacity, not the sum of
all installed hardware.  Finite contacts determine how much hardware joins
the experimental time scale.

## 3. Assembly states and coupling

In addition to the v11 distributed receiver/felt/casing and rear-tube states,
v12 solves:

$$
C_a(T_a)\frac{dT_a}{dt}
=Q_{\rm rec-a}+Q_{\rm tube-a}-Q_{\rm a-felt},
$$

and

$$
C_{h}\frac{dT_h}{dt}
=Q_{\rm sleeve-h}-Q_{\rm h,loss}.
$$

The historical v11 direct receiver-to-tube adaptor conductance is disabled.
All solid heat transfer between the receiver and tube must pass through the
explicit v12 adaptor/contact topology.

### 3.1 Receiver-adaptor contact

The adaptor loosely contacts the receiver over the final 29 mm:

$$
Q_{\rm rec-a}
=h_{\rm rec-a}A_{\rm overlap}(T_{\rm rec,edge}-T_a).
$$

The coefficient is globally bounded and does not vary by experiment.

### 3.2 Adaptor-tube contact

The first 28 mm of the tube exchanges with the surrounding adaptor:

$$
dQ_{\rm tube-a}
=h_{\rm tube-a}\,2\pi R_{\rm tube,o}\,dz\,(T_{\rm tube}-T_a).
$$

There is no tube-to-felt term.

### 3.3 Adaptor-felt coupling

The felt fills the radial space from adaptor radius $R_a$ to enclosure inner
radius $R_f$.  Its cylindrical conductance is represented as

$$
G_{\rm a-f}
=s_{\rm a-f}\frac{2\pi k_{\rm felt}(T)L_a}
{\ln(R_f/R_a)}.
$$

### 3.4 Rear aluminum housing

The missing 28 mm sleeve extension and 18 mm backplate are one rear housing
state.  It is conductively connected to the existing aluminum sleeve and
loses heat by external natural convection and radiation.  The water-cooled
flange remains the terminal boundary of the rear tube.

## 4. Felt properties

The actual felt grade is unknown.  Density remains fixed at
$140\ {\rm kg\,m^{-3}}$.  The nominal conductivity curve interpolates

| Temperature, degC | 20 | 500 | 800 | 1100 | 1400 | 1600 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| $k$, W/m/K | 0.050 | 0.080 | 0.110 | 0.170 | 0.260 | 0.320 |

and is multiplied by one bounded global scale.  Felt heat capacity has a
separate modest global scale.  These scales represent material-grade and
packing uncertainty and may not vary by experiment.

## 5. Thermocouple observation model

The shallow-dip side probes measure perimeter solid temperature through one
shared response time.  The interior probes T9 and T10 are exposed to channel
flow and do not measure a pure solid state.  Their equilibrium target is

$$
T_{\rm int,target}
=w_{\rm wall}T_{\rm wall}
+(1-w_{\rm wall})T_{\rm gas},
$$

followed by a shared first-order response:

$$
\tau_{\rm int}\frac{dT_{\rm int}}{dt}
=T_{\rm int,target}-T_{\rm int}.
$$

T3 has its own shared outlet-probe response time.  The observation states
carry negligible energy and do not alter the assembly balance.

## 6. Identification sequence

1. Use cooling to constrain felt conductivity/capacity, loose
   receiver-adaptor coupling and sensor response times.
2. Hold those quantities fixed while fitting heating profile and power terms
   on the established nine-case training set.
3. Evaluate the six held-out heating cases and C69/C80/C81 without refitting.
4. Require correct radial ordering, credible transient timing, mass/energy
   closure and compatibility with the measured $301\pm23$ J/K slow mode.

Experiment-specific SiC conductivity or heat capacity is prohibited.

## 7. Initial numerical verification

`test/smoke_2D_v12.jl` passes 31 checks for:

- photographed hardware masses and volumes;
- felt conductivity knots;
- uniform equilibrium;
- adaptor and rear-housing response;
- gas/solid interior-probe leverage;
- conserved standard mass flow; and
- equal-pressure convergence; and
- augmented instantaneous energy closure.

The staged screen-to-nominal mesh test passes its 20 K diagnostic gate; the
largest representative final-sensor difference is 16.04 K at E67.

## 8. Completed staged results and interpretation

Cooling selected global felt scales of 1.20 for conductivity and 1.50 for
heat capacity, loose receiver/adaptor contact
$h_{\rm rec-a}=15$ W/m2/K, and shared response times of 180, 60 and 180 s for
side, interior and outlet probes.  The heat-capacity and response-time values
are bound-active and must be treated as confounded effective quantities.

Heating froze those values and physical SiC properties.  Its training screen
selected

$$
\beta_{\rm opt}=150\ {\rm m^{-1}},\qquad
S_h=1.20,
$$

and power scales 1.25, 1.05 and 0.77 for the three irradiance groups.

Nominal-mesh validation gives:

| Phase | Mean sensor RMSE | Steady MAE | t90 MAE | Axial RMSE |
| --- | ---: | ---: | ---: | ---: |
| heating training | 91.70 K | 89.69 K | 756 s | 53.52 K |
| held-out heating | 90.21 K | 86.15 K | 848 s | 51.80 K |
| cooling | 31.17 K | 30.06 K | 700 s | 9.88 K |

Cooling improves strongly relative to staged v11, supporting missing
participating inventory.  Heating remains worse than the accepted v9
reference, and both radial signs are wrong in all 15 cases.  The calculated
local wall and gas temperatures differ by too little for a thermocouple blend
to repair the radial profiles.

Consequently, the fitted felt scale, optical extinction, Graetz multiplier
and power scales are not validated physical coefficients.  The next model
must separate gas-exposed channel-wall material from the square
side/corner/support-connected population while retaining v12's physical
hardware and finite contacts.
