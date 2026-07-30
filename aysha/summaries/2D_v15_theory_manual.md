# 2D v15 theory manual

## Purpose

V15 tests whether the shallow-dip side thermocouples measure a physical
exterior SiC wall temperature that is different from the orbit-averaged
outer-channel solid temperature used by v14.  It retains the v14 square
channel hydraulics and assembly hardware and adds one exterior perimeter
wall state per axial cell.

This is a structural test.  The new state is formed by partitioning the
existing SiC mass.  No heat capacity is added.

## Inherited model

The following are unchanged from the definitive v14 calculation:

- 100 square channels represented by 15 exact symmetry orbits;
- common-pressure-drop parallel-channel flow;
- local air density, viscosity, heat capacity, and conductivity;
- square-duct friction, rear-groove loss, and Graetz heat transfer;
- the 13 mm completely free rear diameter;
- centered Gaussian/Beer-Lambert heating;
- felt, solid alumina adaptor, alumina exit tube, aluminum enclosure and
  water-cooled flange;
- all interior and outlet thermocouple observation models.

The cold hydraulic parameters are frozen before the v15 thermal fit.

## Exterior-wall inventory

For receiver width `W` and fitted effective skin thickness `delta_skin`,

`A_skin = W^2 - (W - 2 delta_skin)^2`.

For boundary orbit `g`,

`A_skin,g = A_skin Nface,g / sum(Nface)`.

The residual orbit area is

`Acore,g = Ng Amonolith/100 - A_skin,g`.

Thus

`sum(Acore,g) + A_skin = A_SiC`

and the receiver mass is identical to v14.  Temperature-dependent SiC heat
capacity is applied separately at the orbit and skin temperatures while
preserving the same material inventory.

## Energy equations

The residual channel-orbit solid satisfies

`Ccore,g dTg/dt = Qsolar,g + Qax,g + Qlat,g - Qgas,g
                  - Qfront,g - Qg-skin`.

The skin satisfies

`Cskin dTskin/dt = Qsolar,skin + Qax,skin
                   + sum(Qg-skin) - Qskin-felt
                   - Qskin-adaptor - Qfront,skin`.

The v14 exterior-channel/felt and exterior-channel/adaptor terms are removed
before the new paths are applied.  The v14 source, axial-conduction area, and
front radiating area are divided between residual channel material and skin
in proportion to the partitioned SiC area.

Channel-to-skin conduction uses the local SiC conductivity, the physical
receiver perimeter area in the axial cell, a half-pitch spreading distance,
and one global dimensionless spreading multiplier.  The multiplier is the
only unresolved coupling parameter; it is not varied by experiment.

## Thermocouples

The side-wall sensors observe the filtered skin state:

| Sensor | Axial position | V15 state |
|---|---:|---|
| T8 | 5 mm | exterior SiC skin |
| T12 | 58 mm | exterior SiC skin |
| T11 | 107 mm | exterior SiC skin |

T9 and T10 retain the central-channel wall/gas blend at 58 and 107 mm.  T3
remains the rear gas measurement and T2 remains the felt measurement.

## Fitting rule and acceptance

All v14 hydraulic, optical, heat-transfer, loss, and hardware parameters are
initially frozen.  V15 screens only effective skin thickness and the shared
channel-to-skin spreading multiplier.  A subsequent heating refit is allowed
only after a skin configuration demonstrates the correct radial mechanism.

The pre-fit numerical and validation criteria are recorded in
`summaries/journal.2D.md`.  In particular, correct radial signs must occur in
at least 12 of 15 heating cases at both measured depths, held-out radial
RMSE must be below 20 K, and cooling and overall held-out performance must
not materially regress from v14.

## Outcome

V15 was implemented and completed through nominal-mesh validation.  The
numerical model passes mass, equilibrium, hydraulic-invariance, energy, and
mesh-transfer tests.

The structural screen evaluated effective skin thicknesses of 0.05, 0.15,
and 0.30 mm and channel/skin spreading multipliers of 0.01, 0.03, 0.10,
0.30, and 1.00.  All v14 parameters were frozen.  None of the 15
configurations produced the measured radial sign in the three representative
screening experiments.

The selected full-training configuration was:

- effective skin thickness: 0.05 mm, the lower search bound;
- channel/skin conductance multiplier: 1.00, the upper search bound;
- training radial signs: 0/9 at both measured depths.

This combination minimizes the independent effect of the new state and
forces it back toward the v14 outer-channel temperature.

Nominal held-out heating results were:

- mean sensor RMSE: 81.50 K;
- axial RMSE: 47.59 K;
- mid-radial RMSE: 31.00 K;
- deep-radial RMSE: 49.49 K;
- DP1 RMSE: 0.280 mbar.

Cooling mean sensor RMSE was 32.93 K.  Across all 15 heating experiments,
the signs of both `T12 - T9` and `T11 - T10` were wrong.

V15 is therefore **rejected for coefficient extraction**.  Under centered
heating, a separately partitioned exterior wall has no preferential energy
source and retains losses to the felt and adaptor.  Reducing its coupling to
the channel makes it colder; increasing coupling recovers v14.  The new mass
state cannot explain the observed hot side wall.

The recommended next model stage is a cooling-first identification of the
skin/felt contact resistance and the felt conductivity and heat-capacity
scales.  Those values should be fitted globally to C69/C80/C81 and then
frozen before the heating profiles are revisited.  This follows the observed
full but non-firm felt contact and avoids using heating errors to identify
the edge loss path.
