# 2D v17 theory manual

## Purpose

V17 corrects a verified hardware boundary omitted from v12--v16: the rear
face of the aluminum casing/backplate is in contact with the water-cooled
flange.  It tests whether this path explains the systematic T2 overprediction
and the nonphysical v16 felt-conductivity fit.

V17 is a hardware ablation, not a new axial optical or hydraulic model.

## Inherited model

V17 inherits the v16 model with:

- 100 physical channels represented by 15 exact D4 orbits;
- conserved MFC mass flow and temperature-dependent gas properties;
- the exterior SiC skin and centered optical source;
- distributed felt and aluminum casing;
- the solid alumina adaptor and alumina tube;
- terminal alumina-tube cooling by the water flange;
- verified side-TC positions of 11, 58, and 107 mm; and
- the inherited 80% wall / 20% gas T9/T10 observation target.

The receiver, felt, casing, adaptor, tube, and housing masses are unchanged.

## New conservative boundary

The rear aluminum housing state now loses

```text
Q_case_flange =
    G_case_flange (T_housing,rear - T_water,flange)
```

with

```text
T_water,flange = 298.15 K.
```

Its energy balance is

```text
C_housing dT_housing,rear/dt =
    Q_case
  - Q_housing,ambient
  - Q_case_flange.
```

The same `Q_case_flange` is added to the flange heat-loss ledger.  Therefore
the new term is exactly conservative.  It is parallel to, and physically
different from, the receiver--adaptor--tube--flange path.

`G_case_flange` includes real contact area, interface resistance, backplate
spreading, flange spreading, and water-side removal.  It is not the bulk
aluminum conductivity.

## Identification

The v16 fit had three quantities in a largely series-coupled cooling path:
skin/felt contact, felt conductivity, and felt heat capacity.  V17 avoids
repeating this structural confounding by fixing:

```text
skin/felt contact scale = 0.30
felt conductivity scale = 1.00
```

Only two shared quantities were screened on C69/C80:

```text
G_case_flange
felt Cp scale.
```

C81 was held out.  Candidates were reranked on the nominal mesh before
selection.

The nominal selection was:

```text
G_case_flange = 20 W/K
felt Cp scale = 1.00.
```

The high-conductance profile continued to improve only weakly through
500 W/K.  Between 20 and 500 W/K, held-out primary cooling RMSE changed from
22.923 to 22.892 K.  Consequently, cooling temperatures identify a clamped
rear-aluminum regime but not the finite conductance.  The 20 W/K value is a
practical lower-bound/plateau choice and is not a measured contact
coefficient.

## Power treatment

After the hardware correction, one delivered-power factor per lamp group was
re-profiled using the same training/held-out split:

```text
456 group: 1.05
304 group: 1.00
256 group: 0.75.
```

The independent closure brackets are 1.05--1.34, 1.23--1.58, and 0.84--1.11.
Only the 456 factor passes.  The large fall from v16's 1.65/1.80/1.25 proves
that nonphysical felt conductivity and delivered power had formed a
compensating loop.  The two low v17 factors show that the inherited
source/exchange model still has incorrect energy partition.

## Validation results

| Metric | Heating training | Heating held out | Cooling |
|---|---:|---:|---:|
| mean sensor RMSE (K) | 105.20 | 97.24 | 35.73 |
| steady MAE (K) | 100.46 | 90.78 | 21.00 |
| t90 MAE (s) | 1082 | 1196 | 1026 |
| axial RMSE (K) | 50.18 | 44.01 | 5.54 |
| mid-radial RMSE (K) | 27.10 | 27.87 | 4.73 |
| deep-radial RMSE (K) | 44.16 | 44.23 | 9.74 |
| DP1 RMSE (mbar) | 0.502 | 0.337 | 0.031 |

The apparent axial improvement is caused by lower fitted power compressing
the modeled differences.  The model still produces zero 58 mm maxima.

T2 steady behavior improves strongly:

```text
heating mean bias: +35.62 K in v16 -> -3.54 K in v17
heating mean RMSE: 32.46 K in v16 -> 7.98 K in v17.
```

Heating T2 t90 MAE worsens from 653 to 1676 s.  The verified path therefore
explains the steady T2 level but not its dynamics.

## Observation-operator check

The claim that T9/T10 were modeled as pure solid temperatures is incorrect.
V17 retains the v14 calibrated 80% wall / 20% gas target.  Re-observing the
solved fields at wall fractions 0, 0.25, 0.55, 0.75, 0.80, and 1 gives zero
correct radial signs out of 15 at both 58 and 107 mm for every blend.  A
pure-gas observation does not solve the current radial failure.

## Invariants and numerical checks

V17 produces

```text
Nu_app = 6.893e-4 Re^1.4415, r2=0.9949
effectiveness = 0.877--0.937
middle peaks = 0/15
positive Lambda_58 = 0/15
positive Lambda_107 = 0/15.
```

Measured values are

```text
Nu_app = 3.047e-4 Re^1.4449
effectiveness = 0.573--0.781
middle peaks = 10/15
positive Lambda_58 = 15/15
positive Lambda_107 = 15/15.
```

The exponent passes, but exchange magnitude is about 126% high and gas
effectiveness is too high and compressed.

The maximum energy residual is `4.3e-14 W`.  Screen-to-nominal maximum
final-sensor changes are 10.11 K for C69, 5.32 K for C81, and 5.21 K for E72,
passing the pre-declared 20 K transfer threshold.

## Interpretation

The aluminum casing/flange path is real and must remain in subsequent models.
Temperature data do not identify its finite conductance once the rear housing
is effectively clamped.  V17 resolves the T2 steady-bias diagnosis and removes
the v16 felt/power compensation, but it worsens simultaneous heating
prediction and does not solve axial or radial signs.

V17 is rejected for coefficient extraction.  The next screen must combine:

1. current versus reduced global wall/gas exchange magnitude; and
2. single versus conservative two-component axial deposition.

The local exchange exponent remains fixed because the measured exponent
already emerges.  Power factors must be bounded by the independent closure
brackets.  Rear-manifold hydraulics remain deferred until this cheaper
exchange/source screen is complete.
