# 1D Theory Manual, v3-v6

This manual summarizes the theory and modeling changes from `1D_v3.jl` through `1D_v6.jl`. These versions share the same core goal: a fast, calibrated, one-dimensional model of a ceramic solar receiver with axial solid temperatures and quasi-steady gas flow.

## 1. Shared Physical Model

All versions model the receiver along the axial coordinate:

```text
z = 0      irradiated inlet/front
z = L      outlet/rear
```

The solid is discretized into axial finite-volume cells. The gas is quasi-steady and marched cell-by-cell from inlet to outlet.

The base solid energy balance for cell `i` is:

```text
C_i dT_s,i/dt =
    Q_solar,i
  - Q_gas,i
  - Q_side,i
  + Q_cond,in,i
  - Q_cond,out,i
  - boundary losses
```

where:

```text
C_i = rho_s gamma_C Cp_s(T_s,i) A_solid dx
```

The gas update uses an exact NTU relation per cell:

```text
UA_i = h_i P_exchange dx
epsilon_i = 1 - exp[-UA_i / (m_dot Cp_f)]
T_g,i+1 = T_g,i + epsilon_i (T_s,i - T_g,i)
Q_gas,i = m_dot Cp_f (T_g,i+1 - T_g,i)
```

## 2. Shared Geometry and Measurements

The 1D model has no radial degree of freedom. From v4 onward, paired thermocouples at similar axial stations are averaged:

```text
T9_pair  = (T9 + T12) / 2
T10_pair = (T10 + T11) / 2
```

The comparison channels from v4 onward are:

```text
T8, T9_pair, T10_pair, T3
```

This was introduced because comparing a radial 3D/experimental thermocouple layout to a 1D model through only one radial thermocouple can bias calibration.

## 3. Gas Heat-Transfer Coefficient

The later 1D models use the same core A/B/C Nusselt form introduced in the 0D work:

```text
Nu = 10^A Re^B Pr^(10^C)
h = Nu k_f / Dh
```

with:

```text
Re = m_dot Dh / (A_flow mu_f)
Pr = Cp_f mu_f / k_f
```

In v3, `A`, `B`, and `C` were fitted. In v5-v6, `B` and `C` were fixed to reduce identifiability problems, while the scale parameter `A_Nu` remained fitted.

## 4. Optical Deposition

The optical source is the total absorbed power times an axial weight:

```text
Q_solar = eta_abs Io A_frt
Q_solar,i = Q_solar w_i
```

### v3 and v4

The base v3/v4 optical form was Beer-Lambert:

```text
w_i = exp(-beta z_left) - exp(-beta z_right)
```

normalized so:

```text
sum_i w_i = 1
```

### v5 Reset

v5 briefly tested a gamma-like source shape:

```text
w(z) proportional to (z/L)^a exp(-beta z)
```

but this was removed. v5 was reset to the `1D_v1.jl` optical parameterization:

```text
Beer-Lambert cell-integrated weights
+ front-deposition fraction added to the first cell
```

That means:

```text
w_i <- (1 - front_dep) w_i
w_1 <- w_1 + front_dep
```

The optical/input-energy values are currently fixed so incoming energy can be adjusted manually:

```text
eta_abs = 0.80
beta_opt = 50.0 1/m
front_dep = 0.50
```

## 5. Version v3

Files:

```text
1D_v3.jl
run_1D_v3.jl
test/smoke_1D_v3.jl
```

Purpose:

- Build a clean conservative finite-volume model.
- Separate model definitions from the runner.
- Use `Optimization` / `OptimizationNLopt` for calibration.

Core parameters after the A/B/C update:

```text
gamma_C, k_scale, U_side, U_rear,
A_Nu, B_Re, C_Pr, tau_T3,
eta_abs, beta_opt, h_front, eps_front
```

Key lessons:

- Front temperatures could be fit reasonably.
- Deeper temperatures and gas outlet showed stronger simulated flow dependence than experiments.
- Direct fitting of many parameters risked compensation instead of diagnosis.

## 6. Version v4

Files:

```text
1D_v4.jl
run_1D_v4.jl
test/smoke_1D_v4.jl
```

Purpose:

- Test axial weakening of gas-solid exchange.
- Use paired thermocouple averages.

v4 introduced a Graetz-style heat-exchange shape inspired by older `0D_v3.jl` experiments:

```text
Gz = Re Pr Dh / z
shape = (1 - a_Gz Gz^n_Gz) exp(-c_Gz / Gz)
h_i = h_base shape
```

The shape was clamped to avoid negative or runaway exchange.

Outcome:

- The optimizer drove the downstream shape close to the lower clamp.
- This supported the idea that downstream exchange should be weaker.
- But the parameterization was too awkward and encouraged bound-hitting fits.

## 7. Version v5

Files:

```text
1D_v5.jl
run_1D_v5.jl
test/smoke_1D_v5.jl
```

Purpose:

- Reduce the number of fitted parameters.
- Remove optical/input-energy fitting.
- Use a simpler axial exchange shape.

The axial exchange shape became:

```text
s(z) = h_floor + (1 - h_floor) exp(-z / L_h)
h_i = h_base s(z_i)
```

This gives two clear physical controls:

```text
h_floor   residual downstream exchange fraction
L_h       axial length scale for exchange decay
```

v5 uses fixed:

```text
B_Re = 1.0
C_Pr = 0.5
eps_front = 0.95
eta_abs = 0.80
beta_opt = 50.0
front_dep = 0.50
```

The staged calibration evolved to:

```text
1. heating base fit, optical fixed
2. cooling thermophysical fit
3. heating refit, optical fixed
```

Observed issue:

- Cooling profiles improved but still cooled too fast in deeper channels.
- `T10_pair` and `T3` retained timing/shape problems.
- This suggested missing thermal storage near the rear/deep region.

## 8. Version v6

Files:

```text
1D_v6.jl
run_1D_v6.jl
test/smoke_1D_v6.jl
```

Purpose:

- Add a rear/adaptor thermal mass.
- Fit this mass only using cooling data.
- Keep front convection fixed to a literature-style natural-convection value.

The state vector becomes:

```text
u = [T_s,1 ... T_s,N, T_rear]
```

The rear mass equations are:

```text
Q_rear_coupling = K_rear (T_s,N - T_rear)
Q_rear_mass_loss = U_rear_mass (T_rear - T_amb)

C_rear dT_rear/dt = Q_rear_coupling - Q_rear_mass_loss
```

with:

```text
C_rear = C_rear_scale rho_s Cp_s(T_rear) A_solid L
```

The last solid cell loses heat to both the direct rear path and the rear mass:

```text
dE_N/dt -= U_rear (T_s,N - T_amb)
dE_N/dt -= K_rear (T_s,N - T_rear)
```

Fixed in v6:

```text
h_front = 10 W/m2/K
eps_front = 0.95
eta_abs = 0.80
beta_opt = 50.0 1/m
front_dep = 0.50
B_Re = 1.0
C_Pr = 0.5
```

Fitted in heating stages:

```text
A_Nu, h_floor, L_h
```

Fitted in cooling stage:

```text
gamma_C, k_scale, U_side, U_rear,
tau_T3, C_rear_scale, K_rear, U_rear_mass
```

v6 calibration sequence:

```text
1. heating heat-transfer calibration
2. cooling thermophysical/rear-mass calibration
3. heating heat-transfer refit
```

Scientific rationale:

- Heating identifies receiver heat-transfer behavior under irradiation.
- Cooling identifies hard-to-quantify thermal masses, losses, and thermophysical properties.
- Repeating heating after cooling checks whether the heat-transfer parameters remain stable once the missing thermal mass is represented.

## 9. Objective Function Interpretation

The calibration objective is a normalized mean squared error:

```text
NMSE = mean(((model - experiment) / scale)^2)
scale = max(max(experiment) - min(experiment), 20 K)
```

Therefore a value such as:

```text
0.015
```

is not a direct 1.5 percent visual error. The square root is more interpretable:

```text
sqrt(0.015) = 0.122
```

That means the RMS error is about 12 percent of the selected signal scale. A temporal profile can still look visibly wrong if the shape or timing is wrong, even when the normalized squared objective appears numerically small.

## 10. Current Diagnostic Meaning

The current modeling evidence suggests:

1. The downstream flow dependence in experiments is weaker than in simple conjugate heat exchange.
2. Axial exchange weakening is needed, but overly flexible forms create bound-hitting fits.
3. Cooling profiles that decay too fast point to missing rear/deep thermal storage.
4. `T8 < T9_pair` at high irradiance can be reproduced by source/loss placement, but the speed at which this intersection moves with flow remains a key diagnostic.
5. The optical/input-energy terms should remain fixed while testing thermal structure, otherwise they can hide missing physics.

## 11. Recommended Use

Use v3 as the conservative baseline.

Use v4 only as evidence that downstream exchange weakening is needed.

Use v5 as the reduced no-rear-mass model.

Use v6 as the current main structural test:

```powershell
$env:RECEIVER1D_V6_RUNNER_QUICK="true"
julia +1 --project=. run_1D_v6.jl
```

Then, for full calibration:

```powershell
julia +1 --project=. run_1D_v6.jl
```

Inspect cooling profiles first. If v6 improves `T10_pair` and `T3` cooling tails without degrading front sensors, the rear thermal mass hypothesis is supported.
