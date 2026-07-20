# 1D_v8b Predicted-Cavity Notes

## Purpose

v8b tests a lumped cavity/T2 prediction without adding a radial mesh. The
receiver and rear tube remain spatially resolved in 1D, while the cavity/felt
and metal casing are represented by one thermal mass:

```text
u = [receiver solid cells, rear alumina tube cells, T_cavity]
```

`T_cavity` is compared directly with measured T2.

## Geometry Used

- Ceramic receiver length: 137 mm.
- Alumina tube length after receiver: 150 mm.
- Total gas domain length: 287 mm.
- Tube inside cavity: 46 mm.
- Tube inside water-cooled flange: 104 mm.
- T3 comparison location: gas temperature at 140 mm.
- Cavity outer diameter: 150 mm.
- Metal thickness: 18 mm.
- Insulation outer radius inside metal: 57 mm.
- Alumina adaptor diameter: 77.6 mm.
- Alumina adaptor length: 57 mm, split as 28 mm receiver overlap and 29 mm tube overlap.

## Remaining Fixed Assumptions

- Rear alumina tube gas radius is inherited as 6.5 mm.
- Rear alumina tube wall thickness is fixed at 1.5 mm until the measured tube OD is confirmed.
- The water-cooled flange is treated as an isothermal 298.15 K sink.
- The old `tau_T3` parameter remains in the vector for compatibility but is not used in v8b;
  T3 is compared directly to gas at 140 mm.

## New Diagnostics

`simulate_v8b` returns:

- `cavity_temperature`
- `receiver_to_cavity_heat`
- `tube_to_cavity_heat`
- `cavity_ambient_heat_loss`
- `flange_heat_loss`
