# 2D v20 theory and identifiability manual

## 1. Purpose and status

V20 is a post-v19 identifiability model. It is not another unconstrained
joint fit and it does not replace the 15-orbit receiver representation.
Its purposes are to:

1. remove the variable-`Cp` gas-energy approximation left in v19;
2. separate the plant fit from the T3 observation;
3. treat the disputed T3 coordinate as discrete model-form uncertainty; and
4. determine whether the existing data identify a receiver integrated
   heat-transfer law after these corrections.

The study objective remains an effective macroscopic receiver coefficient,
not a local single-channel correlation. A v20 parameter is publishable only
if it is interior, transfers to held-out data and mesh, and is not conditional
on an unvalidated T3 operator.

## 2. Architecture retained from v19

V20 retains:

- 15 D4 symmetry orbits representing all 100 physical channels;
- orbit multiplicities, physical channel areas, and total mass flow;
- standard-condition MFC conversion with temperature-dependent local density;
- common-pressure flow allocation and groove obstruction;
- the SiC receiver, full-gap felt, aluminum casing, solid alumina adaptor,
  alumina exit tube, casing/flange contact, and water-cooled flange;
- the centered, exactly normalized source laws;
- the v19 integrated per-channel `Nu/UA` law and normalized axial Graetz
  kernel; and
- the verified side thermocouple coordinates 11, 58, and 107 mm.

No bypass flow, side-weighted irradiance, experiment-specific thermal
properties, independently fitted orbit coefficients, or additional receiver
channels are introduced.

The v19 distributed rear-tube/flange sink is set to zero in the v20
identifiability analysis because its Stage-C optimum was a rejected upper
boundary. The established terminal tube/flange path remains active.

## 3. Exact variable-heat-capacity gas thermodynamics

### 3.1 Air enthalpy

The inherited air heat capacity is

```text
Cp(T) = 1004 [1 + 1.983e-4 T - 4.14e-8 T^2] J/(kg K),
```

with the property temperature clamped below 200 K. V20 integrates this
analytically to a specific enthalpy relative to 200 K:

```text
h(T) = integral(200,T) Cp(theta) dtheta.
```

Below 200 K, the same constant `Cp(200)` continuation used by the property
function is integrated. A safeguarded Newton/bisection inverse recovers
`T(h)` over 100--3500 K.

### 3.2 Receiver-cell exchange

For a cell with constant wall temperature over its local conductance, the gas
outlet is the bounded solution of

```text
UA_cell/mdot = integral(Tin,Tout) Cp(T)/(Twall-T) dT.
```

Because `Cp` is quadratic, the integral has an analytic logarithmic
antiderivative. V20 solves the resulting scalar relation with a safeguarded
Newton iteration. The heat transferred to the gas is then exactly

```text
Qgas = mdot [h(Tout)-h(Tin)].
```

The same update is used in the receiver channels and rear alumina tube.
Heating and cooling remain bounded between the inlet-gas and wall
temperatures.

### 3.3 Orbit mixing

V19 approximated mixing using `Cp(T)T` weights. V20 uses:

```text
hmix = sum_g mdot_g h(Tout,g) / sum_g mdot_g
Tmix = inverse_h(hmix).
```

This closes orbit mixing exactly and makes the receiver plus rear-tube gas
gain telescope to

```text
sum(Qgas) = mdot_total [h(Trear,out)-h(Treceiver,in)].
```

Saved diagnostics separately report receiver, orbit-mixing, rear-tube, and
total-gas enthalpy residuals.

## 4. T3 coordinate and observation branches

The file audit does not support calling the previous global 140 mm coordinate
verified.

Evidence is ranked as follows:

1. **136 mm:** strongest contemporaneous evidence. Early experiment-specific
   scripts explicitly sampled approximately 135.6--136 mm and described it
   as T3 around 136 mm.
2. **137 mm:** the exact geometric receiver exit. This is a clean convention
   for interpreting “at the exit,” but not an independently measured probe
   coordinate.
3. **140 mm:** 3 mm into the rear tube. This is a later inherited model
   convention attributed secondhand to a user specification.

The photograph supports only an outlet-region, off-axis placement and cannot
resolve this 4 mm difference. Therefore v20 implements three discrete
branches:

```text
:receiver_136   global x = 136 mm
:receiver_exit  global x = 137 mm
:rear_003       rear-tube z = 3 mm, global x = 140 mm.
```

At receiver locations, gas is mixed across the 15 orbits by mass and
enthalpy, and the radiative target is the local multiplicity-weighted channel
wall. At 140 mm, the gas and surrounding temperature come from the rear
alumina tube.

The observation alternatives are:

- `P0`: local mixed gas, with no lag or wall/stem exchange;
- `P1`: passive finite-capacity sheath with fixed diameter/emissivity and no
  stem sink; and
- `P2`: the same sheath with a weak lead/stem sink.

The archived evidence does not verify T3 sheath diameter, emissivity, heat
capacity, or lead conductance. The broad capacity/stem grid is therefore an
identifiability profile, not a source of validated probe properties.
Position is never optimized continuously.

## 5. T3-free plant fitting

T3, T9, and T10 are excluded from all plant objective functions. Heating uses
only the side thermocouples T8/T12/T11 and T2. Cooling uses the same four
measurements to constrain felt/assembly behavior.

The corrected cooling initializer is also T3/T9/T10-free:

- T8/T12/T11 initialize the perimeter axial profile;
- the unobserved core is initialized from the perimeter profile;
- T2 initializes the outer assembly; and
- the rear tube is initialized from the back-side T11 temperature.

Measured cooling T3 at `t0` may initialize the one-way probe ODE when T3 is
being scored, but it does not feed back into the plant state.

The integrated law is profiled using the better-conditioned pair
`(Nu50,n)`:

```text
Nu_bar = A_Nu Re^n
A_Nu = Nu50 / 50^n.
```

The nominal-mesh grid is:

```text
n = 1.25, 1.35, 1.45, 1.55, 1.65
Nu50/Nu50_measured = 0.75, 1.35, 1.95.
```

The three heating anchors E72, E69, and E81 span irradiance and flow. The
profile score contains side curves, side final levels, T2 curves, and T2
levels. T3 is reported only as a diagnostic.

## 6. Verification

`test/smoke_2D_v20.jl` checks:

- analytic enthalpy monotonicity and inverse round trips;
- numerical `dh/dT` agreement with inherited `Cp(T)`;
- exact heating and cooling cell-integral closure;
- zero-conductance and zero-flow limits;
- 100 represented channels and exact source normalization;
- a transient solve and all T3-location plumbing;
- exact orbit/rear/total gas enthalpy closure; and
- coupled-flow convergence in the short test.

The smoke suite passes 28/28 tests. The six-case frozen nominal-mesh
comparison gives gas enthalpy residuals below `3.36e-9 W` absolute and
`3.11e-11` relative. The E76 strict flow flag remains false for the same
small capped-iteration pressure residual already documented in v19.

## 7. Frozen v19-to-v20 thermodynamic comparison

The first comparison freezes every v19 plant parameter, including the
rejected distributed tube sink, so only gas thermodynamics changes. Across
E67/E72/E76/E80/C69/C81 on the nominal mesh:

- all-sensor RMS change is 0.030--0.606 K;
- maximum sensor change is 1.45 K;
- side RMS change is 0.013--0.470 K;
- T2 RMS change is below 0.036 K; and
- receiver-exit-gas RMS change is 0.088--1.862 K.

Thus the v19 `Cp*DeltaT` and `Cp*T` treatment was a real conservation
approximation, but it is much too small to explain the 40--180 K model-data
errors or the fitted boundary behavior. Re-fitting unrelated plant
parameters merely because enthalpy accounting changed is not justified.

## 8. Discrete T3 observer result

The nominal-mesh observer profile solves C69/C80/C81 plant trajectories with
the side/T2-only initializer and no T3 in the plant score. It then evaluates
all three locations over:

```text
C'' = 200, 600, 1200, 2400, 3000 J/(m2 K)
G''stem = 0, 20, 60, 120, 200 W/(m2 K).
```

The lowest training score occurs at:

```text
location = 137 mm receiver exit
C'' = 3000 J/(m2 K), upper boundary
G''stem = 0 W/(m2 K), lower boundary.
```

It still gives:

| quantity | result |
|---|---:|
| C69/C80 mean T3 RMSE | 69.93 K |
| worst training-case T3 RMSE | 94.96 K |
| C81 T3 RMSE | 38.28 K |
| C81 final bias | +42.49 K |
| C81 `t90` error | -1800 s |

No coordinate/operator combination passes. The earlier full-sensor
initializer weakly preferred 140 mm, whereas the leakage-free initializer
weakly prefers 137 mm. This branch switching is itself evidence that the
coordinate is not identifiable from these trajectories. It must not be
reported as a T3-location determination.

## 9. T3-free integrated-UA profile

The corrected nominal-mesh profile contains 15 structural cells. Its minimum
is:

```text
n = 1.25, the lower exponent boundary
Nu50/Nu50_measured = 1.35
A_Nu = 8.85657e-4
objective = 4.0651.
```

The corresponding side-history RMSE is 125.96 K and T2 RMSE is 10.81 K.
Four candidates lie within 10% of the minimum and span:

```text
n = 1.25--1.35
Nu50/Nu50_measured = 1.35--1.95.
```

The near-tied candidate at `n=1.25`, `Nu50` ratio 1.95 has a substantially
different final side bias while retaining nearly the same aggregate score.
The surface is therefore a boundary ridge, not an interior coefficient
basin. `A_Nu` and `n` are not identified by side/T2 plant data.

## 10. Bounded nuisance screen

The two near-tied ridge endpoints are carried into a nominal-mesh bounded
screen:

```text
felt k scale  = 0.70, 1.00, 1.30
felt Cp scale = 0.75, 1.00, 1.50
456 power     = 1.05, 1.195, 1.34
304 power     = 1.23, 1.405, 1.58
256 power     = 0.84, 0.975, 1.11.
```

Cooling selects felt from side/T2 only. The three irradiance groups then
select power independently from representative heating cases. This is a
lexicographic bounded screen, not a fully joint covariance calculation.
Every candidate is checked for coupled-flow convergence, and structural,
felt, and power boundary activity is explicit.

The completed numerical result is recorded in
`summaries/2D_v20/plant_nuisance_selected*_2D_v20.csv`.

Both ridge endpoints select the same lower boundaries:

```text
felt k scale = 0.70
felt Cp scale = 0.75
power scales 456/304/256 = 1.05/1.23/0.84.
```

Their plant results are:

| `n`, `Nu50` ratio | heating side | heating T2 | cooling side | cooling T2 |
|---|---:|---:|---:|---:|
| 1.25, 1.35 | 102.47 K | 10.20 K | 66.42 K | 14.55 K |
| 1.25, 1.95 | 79.40 K | 8.47 K | 51.95 K | 14.48 K |

All coupled-flow solves used by the selected rows converge. Nevertheless,
both heating and cooling gates fail, all five nuisance selections are at
their lower bounds, and `n=1.25` is itself at the structural boundary.
`plant_admissible=false` for both endpoints.

The screen therefore stops here. A global joint nuisance fit, T3-conditioned
UA profile, bootstrap coefficient interval, held-out branch selection, and F
mesh confirmation are not authorized because they would quantify a rejected
plant rather than identify a coefficient.

## 11. Interpretation rule

V20 can establish non-identifiability even if it does not find an accepted
model. A coefficient claim requires:

- a T3-free admissible plant;
- an interior connected `(Nu50,n)` basin;
- a T3 location/operator that transfers from C69/C80 to C81;
- no material dependence on initialization convention;
- held-out heating and M/F transfer; and
- conservation and coupled-flow convergence.

If the plant or observer gates fail, T3 cannot rescue the model and no
bootstrap confidence interval for a point coefficient is meaningful. The
defensible result is then a branch-conditional apparent-UA interval or a
formal non-identifiability statement, together with the directly measured
apparent correlation.

V20 reaches that stopping condition. It corrects thermodynamic bookkeeping
and makes the inverse-problem failure cleaner, but it does not produce a
validated `A_Nu`, `n`, felt property, power scale, or T3 probe parameter.

## 12. Files

- model: `2D_v20.jl`
- experiment builder: `run_2D_v20.jl`
- smoke tests: `test/smoke_2D_v20.jl`
- frozen thermodynamic comparison:
  `compare_2D_v19_v20_enthalpy.jl`
- T3 profile: `profile_2D_v20_t3_observer.jl`
- T3-free UA profile: `profile_2D_v20_ua_plant.jl`
- nuisance screen: `screen_2D_v20_plant_nuisance.jl`
- profile aggregation: `aggregate_2D_v20_ua_profile.jl`
- numerical outputs: `summaries/2D_v20/`
