# 1D Theory Manual, v8 Rear-Domain and Predicted-Cavity Models

This manual documents the v8 family of one-dimensional receiver models:

- `1D_v8.jl`: v8a, rear alumina tube extension with measured T2 as a boundary.
- `1D_v8b.jl`: v8b, rear alumina tube extension with a predicted lumped cavity/T2 state.

Both versions keep the core v5-v7 receiver model, but extend the thermal domain
behind the 137 mm ceramic receiver. The purpose is to test whether the
experimental/numerical mismatch is caused by missing heat storage and heat loss
behind the modeled receiver length.

The target discrepancies are:

- heat transfer deeper in the receiver appearing weaker and less flow-dependent
  than the simple 137 mm model predicts,
- slower deep receiver response,
- difficulty fitting solid and gas steady-state temperatures at the same time,
- overprediction of T2 insulation temperature in higher-dimensional simulations.

v8 does not try to solve all of these by adding arbitrary fitted rear losses.
Instead, it introduces the known rear alumina tube, adaptor, cavity insulation,
metal casing, and water-cooled flange as fixed geometry-driven thermal paths.

## 1. Coordinate System

The axial coordinate follows the previous 1D models:

```text
z = 0        irradiated inlet/front face
z = L        ceramic receiver outlet
z > L        downstream alumina tube / flange region
```

The original ceramic receiver length remains:

```text
L = 137 mm
```

The v8b gas domain extends to:

```text
L_total = 137 mm + 150 mm = 287 mm
```

In v8b, the experimental T3 gas comparison is not taken at the final numerical
outlet. It is sampled from the gas profile at:

```text
z_T3 = 140 mm
```

This places T3 just downstream of the ceramic receiver, consistent with the
user-specified comparison location.

## 2. Baseline Receiver Model Inherited from v7

The ceramic receiver is still represented by axial finite-volume solid
temperatures:

```text
T_s,1 ... T_s,N
```

The gas is not an ODE state. It is quasi-steady and is marched from inlet to
outlet through the current solid temperature field using exact cell NTU updates.

For receiver cell `i`:

```text
UA_i = h_i P_exchange dx
epsilon_i = 1 - exp[-UA_i / (m_dot Cp_f)]
T_g,i+1 = T_g,i + epsilon_i (T_s,i - T_g,i)
Q_gas,i = m_dot Cp_f (T_g,i+1 - T_g,i)
```

The heat-transfer coefficient keeps the reduced v5-v7 form:

```text
Nu = 10^A_Nu Re^B_Re Pr^(10^C_Pr)
h_i = Nu k_f / Dh * s(z_i)
```

with fixed:

```text
B_Re = 1.0
C_Pr = 0.5
```

The axial exchange shape remains:

```text
s(z) = h_floor + (1 - h_floor) exp(-z / L_h)
```

where:

```text
h_floor = downstream residual exchange fraction
L_h     = axial decay length for gas-solid exchange
```

The stagnant/zero-flow behavior follows v7: heat exchange is set to zero, while
the diagnostic gas profile is assigned local solid/tube temperatures instead of
being forced to ambient inlet temperature.

## 3. Solar Source and Front Boundary

The optical deposition model remains the v5-v7 controlled form:

```text
Q_solar = eta_abs f_I Io A_frt
```

where:

```text
eta_abs = 0.80
beta_opt = 50 1/m
front_dep = 0.50
```

The normalized axial weights use Beer-Lambert attenuation plus a fixed front
deposition fraction. The first cell receives the additional front-deposition
source.

Irradiance correction is grouped by nominal irradiance:

```text
f_I = f_I_high    Io >= 3.80e5 W/m2
f_I = f_I_mid     2.80e5 <= Io < 3.80e5 W/m2
f_I = f_I_low     Io < 2.80e5 W/m2
```

The front face loses heat by fixed natural convection and radiation:

```text
Q_front =
    h_front A_frt (T_s,1 - T_amb)
  + eps_front sigma A_frt (T_s,1^4 - T_amb^4)
```

with:

```text
h_front = 10 W/m2/K
eps_front = 0.95
```

## 4. Version v8a: Rear Tube with Measured T2 Boundary

### 4.1 Purpose

v8a was introduced as the first structural test of the rear heat-loss
hypothesis. It keeps T2 as a measured thermal boundary, as in v7, but replaces
the single receiver rear-to-T2 path with an explicit rear alumina tube domain.

v8a asks:

```text
If the receiver is connected to a real downstream tube and cooled flange,
does the added fixed heat path change deep receiver and T3 behavior in the
right direction?
```

### 4.2 Geometry in v8a

The v8a geometry used the first rear-tube estimate:

```text
ceramic receiver length         137 mm
rear alumina tube length        200 mm
tube inside cavity/felt         100 mm
tube inside cooled flange       100 mm
water flange temperature        298.15 K
```

The alumina adaptor was treated as split evenly across receiver and tube:

```text
adaptor length                  57 mm
receiver overlap                28.5 mm
tube overlap                    28.5 mm
```

T2 remains prescribed from experiment:

```text
T2_boundary(t) = measured T2(t)
```

### 4.3 State Vector in v8a

v8a adds dynamic tube wall states:

```text
u = [T_s,1 ... T_s,N, T_tube,1 ... T_tube,M]
```

There is still no dynamic cavity state in v8a.

### 4.4 Receiver Balance in v8a

Each receiver solid cell follows:

```text
C_s,i dT_s,i/dt =
    Q_solar,i
  - Q_gas,i
  - Q_side,T2,i
  + Q_axial,in,i
  - Q_axial,out,i
```

where:

```text
Q_side,T2,i = k_ins_scale G_receiver_to_T2 dx (T_s,i - T2_boundary)
```

The final receiver cell is coupled to the first rear tube cell:

```text
Q_receiver_tube = G_receiver_tube (T_s,N - T_tube,1)
```

The receiver loses this heat and the first tube cell gains it.

### 4.5 Tube Balance in v8a

The rear alumina tube is discretized in axial cells. Each tube cell has:

```text
C_tube,j dT_tube,j/dt =
    Q_axial,in,j
  - Q_axial,out,j
  - Q_gas,tube,j
  - Q_tube,T2,j
  - Q_flange,j
```

The tube gas is marched after the receiver gas:

```text
T_g,receiver_out -> rear tube gas cells -> final outlet
```

Tube-to-gas exchange uses an internal pipe-style heat-transfer coefficient:

```text
Nu = 3.66                              Re < 2300
Nu = 0.023 Re^0.8 Pr^0.4              Re >= 2300
```

The cavity portion of the tube loses radially to measured T2:

```text
Q_tube,T2,j = G_tube_to_T2,j dx_rear (T_tube,j - T2_boundary)
```

The flange portion loses to the water-cooled sink:

```text
Q_flange,j = G_flange,j dx_rear (T_tube,j - T_water)
```

where:

```text
T_water = 298.15 K
```

The flange conductance is based on tube-wall alumina conduction plus an aluminum
radial path. It is geometry-based, not fitted.

### 4.6 v8a Interpretation

v8a is useful as a diagnostic stepping stone, but it cannot predict T2 because
T2 is prescribed. It tests the magnitude of the rear tube/flange heat path while
keeping the cavity temperature externally imposed.

Important v8a diagnostic outputs are:

```text
rear_tube_temperature
receiver_to_tube_heat
tube_to_t2_heat_loss
flange_heat_loss
```

These quantify whether the new downstream domain is extracting energy at a
credible level.

## 5. Version v8b: Predicted Lumped Cavity/T2 State

### 5.1 Purpose

v8b is the main model for testing the T2 hypothesis. It removes measured T2 as
a boundary condition and introduces one predicted cavity thermal state:

```text
T_cavity(t)
```

This state represents the combined thermal mass of:

```text
felt insulation + metal casing + unresolved cavity structure
```

It is not a radial mesh. The model does not solve radial temperature gradients.
Instead, it treats the cavity/felt/metal system as one effective thermal
reservoir whose temperature is compared directly with measured T2.

### 5.2 Corrected Geometry in v8b

The geometry follows the user-supplied schematic and corrections:

```text
ceramic receiver length                 137 mm
rear alumina tube length after receiver 150 mm
total gas domain length                 287 mm
tube length inside cavity                46 mm
tube length inside water-cooled flange  104 mm
T3 comparison position                  140 mm
cavity outer diameter                   150 mm
cavity outer radius                      75 mm
metal thickness                          18 mm
insulation outer radius inside metal     57 mm
alumina adaptor diameter                 77.6 mm
alumina adaptor radius                   38.8 mm
adaptor length                           57 mm
adaptor receiver overlap                 28 mm
adaptor tube overlap                     29 mm
```

The rear alumina tube still uses the inherited gas radius:

```text
r_tube,gas = 6.5 mm
```

The tube wall thickness is currently fixed as:

```text
t_tube = 1.5 mm
```

This should be updated if the physical tube outer diameter is measured.

### 5.3 State Vector in v8b

v8b solves:

```text
u = [
    T_s,1 ... T_s,N,
    T_tube,1 ... T_tube,M,
    T_cavity
]
```

where:

```text
T_s,i       ceramic receiver solid cell temperatures
T_tube,j    alumina rear tube/adaptor wall temperatures
T_cavity    lumped cavity/felt/metal temperature, compared to T2
```

The gas profile is still quasi-steady and is recalculated from the current
solid/tube temperatures during each RHS evaluation.

### 5.4 Receiver-to-Cavity Coupling

In v8b, receiver side heat loss no longer goes to measured T2. It goes to the
predicted cavity state:

```text
Q_side,cavity,i =
    k_ins_scale G_receiver_to_cavity dx (T_s,i - T_cavity)
```

The receiver cell loses this energy:

```text
dE_s,i/dt -= Q_side,cavity,i
```

and the lumped cavity gains it:

```text
dE_cavity/dt += sum_i Q_side,cavity,i
```

The geometry-based conductance per length is:

```text
G_receiver_to_cavity =
    4 pi k_ins / log(r_ins / r_receiver_eq)
```

The factor is intentionally geometry-based and the only scaling comes from the
existing `k_ins_scale` parameter. No new rear/cavity fit parameter is added.

### 5.5 Receiver-to-Tube Coupling

The last receiver cell is connected to the first rear tube cell through the
adaptor/contact path:

```text
Q_receiver_tube =
    G_receiver_tube (T_s,N - T_tube,1)
```

The contact conductance follows the previous 0D/v7 style:

```text
R_contact = 1 / (h_contact A_contact)
G_receiver_tube = 1 / R_contact
```

with:

```text
h_contact = 100 W/m2/K
A_contact = 4 w_t L_overlap
```

For v8b:

```text
L_overlap = 28 mm
```

### 5.6 Tube Axial Conduction

Neighboring rear tube cells exchange heat by axial alumina conduction:

```text
Q_cond,tube,j+1/2 =
    k_face A_face (T_tube,j - T_tube,j+1) / dx_rear
```

with harmonic mean conductivity:

```text
k_face = 2 k_j k_j+1 / (k_j + k_j+1)
```

The alumina conductivity function is inherited from the older 0D material set:

```text
k_alumina(T) = 5.5 + 34.5 exp[-0.0033 (T - 273.15)]
```

The alumina heat capacity is:

```text
Cp_alumina(T) =
    (1.00446 + 1.742e-4 T - 2.796e4 T^-2) * 1000
```

with density:

```text
rho_alumina = 3900 kg/m3
```

### 5.7 Tube-to-Cavity Coupling

Only the first 46 mm of the rear tube is inside the cavity. For tube cells in
that region:

```text
Q_tube,cavity,j =
    G_tube_to_cavity,j dx_rear (T_tube,j - T_cavity)
```

Tube cells outside the cavity have:

```text
G_tube_to_cavity,j = 0
```

The tube loses this heat and the lumped cavity gains it:

```text
dE_tube,j/dt -= Q_tube,cavity,j
dE_cavity/dt += Q_tube,cavity,j
```

### 5.8 Tube-to-Flange Cooling

The remaining 104 mm of tube is inside the water-cooled metal flange. For those
cells:

```text
Q_flange,j =
    G_flange,j dx_rear (T_tube,j - T_water)
```

where:

```text
T_water = 298.15 K
```

The flange conductance combines:

```text
alumina tube-wall resistance
aluminum radial conduction resistance
```

The aluminum conductivity is fixed:

```text
k_aluminum = 205 W/m/K
```

The flange is treated as isothermal because it is water cooled. This is a strong
sink assumption and should be interpreted carefully when comparing gas
temperatures downstream of the receiver.

### 5.9 Cavity Thermal Mass

The cavity state uses a fixed heat capacity:

```text
C_cavity =
    rho_ins Cp_ins V_felt
  + rho_al Cp_al V_metal
```

with:

```text
rho_ins = 140 kg/m3
Cp_ins = 1360 J/kg/K
rho_aluminum = 2700 kg/m3
Cp_aluminum = 900 J/kg/K
```

The v8b quick-run metadata currently reports:

```text
C_cavity = 4025.73 J/K
```

This is the heat capacity of the unresolved cavity/felt/metal lump. It is fixed,
not fitted.

### 5.10 Cavity Energy Balance

The cavity ODE is:

```text
C_cavity dT_cavity/dt =
    sum_i Q_receiver,cavity,i
  + sum_j Q_tube,cavity,j
  - Q_cavity,ambient
```

where:

```text
Q_cavity,ambient =
    h_nat A_outer (T_cavity - T_amb)
  + eps_metal sigma A_outer (T_cavity^4 - T_amb^4)
```

Fixed casing loss values are:

```text
h_nat = 10 W/m2/K
eps_metal = 0.2
```

This means `T_cavity` is a hybrid effective temperature. It is compared to T2,
but it also stands in for the unresolved felt/metal thermal reservoir. This is
the price of avoiding a radial mesh.

## 6. T3 Handling in v8b

Earlier versions compared T3 to the final gas outlet plus an empirical first
order lag:

```text
T3_model = lag(T_g,out, tau_T3)
```

v8b changes this. The user-specified comparison is:

```text
T3_model = T_g(z = 140 mm)
```

Therefore:

- the old `tau_T3` parameter remains in the 10-parameter vector for
  compatibility,
- `tau_T3` is not used in v8b output extraction,
- `tau_T3` is not included in the v8b cooling fit subset.

This makes the T3 comparison a physical location assumption rather than a fitted
sensor-lag correction.

## 7. Experimental Output Channels

v8b compares five measured channels:

```text
T8
T9_pair  = (T9 + T12) / 2
T10_pair = (T10 + T11) / 2
T3       = gas at z = 140 mm
T2       = T_cavity
```

The model output matrix has six columns:

```text
1  T8
2  T9_pair
3  T10_pair
4  T3
5  T2
6  h_effective
```

`h_effective` is the mean receiver gas-solid heat-transfer coefficient over the
receiver cells. It is diagnostic only; it is not compared to experiment.

## 8. Parameter Vector

v8 keeps the 10-parameter vector introduced in v7:

```text
p[1]  gamma_C       receiver solid heat-capacity multiplier
p[2]  k_scale       receiver axial-conductivity multiplier
p[3]  k_ins_scale   receiver/tube-to-cavity conductance multiplier
p[4]  A_Nu          Nusselt multiplier exponent
p[5]  h_floor       downstream exchange floor
p[6]  L_h           axial exchange decay length
p[7]  tau_T3        retained for compatibility; unused in v8b
p[8]  f_I_high      high-irradiance correction factor
p[9]  f_I_mid       middle-irradiance correction factor
p[10] f_I_low       low-irradiance correction factor
```

v8b heating-stage fit:

```text
A_Nu, h_floor, L_h, f_I_high, f_I_mid, f_I_low
```

v8b cooling-stage fit:

```text
gamma_C, k_scale, k_ins_scale
```

No new rear tube, flange, or cavity heat-capacity parameters are fitted.

## 9. Calibration Objective

v8b inherits the v7 signal loss:

```text
loss =
    normalized_mse
  + SLOPE_WEIGHT normalized_slope_mse
  + TIMING_WEIGHT timing_penalty
```

with:

```text
SLOPE_WEIGHT = 0.25
TIMING_WEIGHT = 0.15
```

The case loss averages over five channels:

```text
T8, T9_pair, T10_pair, T3, T2
```

This is important: T2 is not just plotted afterward. It participates in the
objective and can penalize wrong cavity-level predictions.

The loss is still an engineering diagnostic objective, not a strict statistical
likelihood. Its slope and timing terms are intended to discourage fits that
match only final temperatures while missing response shape.

## 10. Heat-Flow Diagnostics

v8b returns several heat-flow histories:

```text
receiver_to_tube_heat
receiver_to_cavity_heat
tube_to_cavity_heat
cavity_ambient_heat_loss
flange_heat_loss
```

These are essential for interpreting whether a calibration is physically
reasonable.

Useful questions:

1. Is the flange heat loss dominating the receiver energy balance?
2. Does the receiver-to-cavity heat flow have the correct sign during heating
   and cooling?
3. Is `T_cavity` being driven mostly by receiver side losses or by the rear tube?
4. Does the cavity cool too quickly because the metal casing loss is too strong?
5. Does the tube gas cool substantially before z = 140 mm?

The v8b model should be judged by these heat-flow histories, not only by RMSE.

## 11. Numerical Solution

The ODE system is solved using:

```text
Rodas5P(autodiff = AutoFiniteDiff())
```

The solver saves only at experimental times:

```text
saveat = time
save_everystep = false
dense = false
```

The simulation also passes:

```text
tstops = time
```

This helps the solver step exactly through piecewise-linear experimental input
history knots.

The gas profile is recomputed during post-processing at every saved time using
the solved receiver/tube/cavity states.

## 12. Quick-Run Results and Interpretation

The v8b quick runner completed and wrote:

```text
summaries/1D_v8b/parameters_1D_v8b.csv
summaries/1D_v8b/analysis_results_1D_v8b.csv
summaries/1D_v8b/plots/
```

The quick-run fitted vector was:

```text
gamma_C   = 1.5875
k_scale   = 1.4229
k_ins     = 0.8711
A_Nu      = 1.8009
h_floor   = 0.2216
L_h       = 0.0444 m
tau_T3    = 20.0     unused in v8b
f_I_high  = 0.9370
f_I_mid   = 1.6000
f_I_low   = 1.1320
```

Quick-run mean metrics were:

```text
heating T8        RMSE 165.1 K    steady -160.9 K
heating T9_pair   RMSE 175.9 K    steady -185.2 K
heating T10_pair  RMSE  99.8 K    steady -101.0 K
heating T3        RMSE  74.0 K    steady  -57.7 K
heating T2        RMSE   6.8 K    steady  -10.4 K

cooling T8        RMSE  13.0 K    steady   11.6 K
cooling T9_pair   RMSE  24.4 K    steady   11.2 K
cooling T10_pair  RMSE  16.5 K    steady   15.7 K
cooling T3        RMSE  26.5 K    steady   20.9 K
cooling T2        RMSE   2.4 K    steady    3.5 K
```

These quick-run numbers should not be treated as final calibration results
because the quick runner uses a tiny optimizer budget. Still, they show useful
structure:

- predicted T2 is in the right range even with a single cavity state,
- cooling behavior is relatively good,
- heating solid temperatures remain underpredicted in the quick fit,
- T3 no longer uses a lag but remains sensitive to the 140 mm comparison
  location and rear-tube gas interaction.

## 13. Scientific Interpretation of v8b

v8b is a targeted compromise:

```text
more physical than measured-T2 boundary,
less flexible than a fitted hidden rear mass,
much cheaper and more identifiable than a radial cavity mesh.
```

If v8b predicts T2 well while still underpredicting receiver heating, the
remaining error is probably not simply "missing cavity heat capacity." It may
instead involve:

- optical/source placement,
- gas bypass or inactive flow fraction,
- excessive modeled gas-solid exchange,
- incorrect T3 measurement location or mixing assumption,
- axial conduction that is too strong for the real porous/channel structure.

If v8b fails badly on T2, then a single lumped cavity is too simple and the next
step should be a small number of axial cavity zones or separate felt/metal
lumps, not a full radial mesh unless absolutely necessary.

## 14. Limitations and Caveats

1. `T_cavity` is not a literal local material temperature.

It represents a combined felt/metal/cavity thermal reservoir. Comparing it to
T2 is useful, but a thermocouple embedded in insulation may not exactly equal
the volume-average cavity lump.

2. The metal casing loss uses `T_cavity` as the outer casing temperature.

This is a simplification. The real metal casing could be cooler or warmer than
the T2 thermocouple location depending on radial gradients.

3. The rear tube wall thickness is assumed.

The current value is:

```text
t_tube = 1.5 mm
```

This affects tube thermal mass, tube-gas exchange, and flange cooling. It should
be replaced by the measured tube OD when available.

4. The water-cooled flange is an isothermal sink.

This likely captures the strong cooling trend, but it may overstate heat
removal if contact to the tube is imperfect or if the tube passes through a
clearance/seal.

5. The adaptor diameter interpretation should remain visible.

v8b uses:

```text
adaptor diameter = 77.6 mm
```

If the intended value was different, the adaptor/cavity conductances and
volumes should be revised.

6. `tau_T3` is inactive in v8b.

It is retained only to preserve the existing vector and runner structure. Full
v8b interpretation should not treat it as a fitted sensor-lag result.

## 15. Recommended Analysis Workflow

1. Run smoke validation:

```powershell
& "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. test\smoke_1D_v8b.jl
```

2. Run quick v8b:

```powershell
$env:RECEIVER1D_V8B_RUNNER_QUICK="true"
& "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. run_1D_v8b.jl
```

3. Inspect T2 first.

If predicted T2 is wrong in level or timing, fix the cavity representation
before judging gas heat transfer.

4. Inspect heat-flow diagnostics.

The key signals are:

```text
receiver_to_cavity_heat
tube_to_cavity_heat
cavity_ambient_heat_loss
flange_heat_loss
receiver_to_tube_heat
```

5. Compare v7, v8a, and v8b with the same parameter vector.

This isolates the structural effect of the rear domain and cavity state from
optimizer compensation.

6. Only then run a full calibration.

Full calibration should report not just objective value but signed residuals,
T2 errors, T3-at-140 mm behavior, and heat-flow magnitudes.

## 16. Open Questions Before Publication-Grade Use

1. What is the measured rear alumina tube outer diameter?
2. Is the flange contact direct, through a seal, or through a clearance?
3. Does T2 physically correspond more closely to felt, metal, or an average
   cavity temperature?
4. Is the T3 thermocouple really best represented by gas at 140 mm, or does it
   have its own thermocouple bead/sleeve thermal inertia?
5. Should the cavity lump be split into two axial zones if T2 level is correct
   but timing remains wrong?

## 17. Summary

v8a and v8b convert the rear of the receiver from a missing boundary condition
into an explicit physical hypothesis.

v8a:

```text
receiver + rear tube + measured T2 boundary + cooled flange
```

v8b:

```text
receiver + rear tube + predicted lumped cavity/T2 + cooled flange
```

The most important v8b advance is that T2 becomes a prediction while the model
still avoids a radial mesh. This keeps the model fast and identifiable while
testing whether the full cavity thermal mass can explain the experimentally
observed slow, flow-weak rear response.
