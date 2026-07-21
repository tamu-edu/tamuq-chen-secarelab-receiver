# 1D_v8a Rear-Tube Extension Notes

## Purpose

v8a tests whether the mismatch between experiments and the 137 mm receiver model
can be explained by a fixed heat-loss/storage path behind the ceramic receiver.
No new rear-component parameters are fitted.

## Geometry Added

- The ceramic receiver remains 137 mm long.
- A 200 mm alumina rear tube is added downstream of the receiver.
- The first 100 mm of the tube is inside the cavity/felt environment.
- The final 100 mm is inside a water-cooled aluminum flange held at 298.15 K.
- The 57 mm alumina adaptor is split evenly:
  - 28.5 mm overlaps the receiver side.
  - 28.5 mm sleeves the rear-tube side.

## Fixed Assumptions

- The existing 6.5 mm adaptor tube radius is used as the rear gas radius.
- Tube wall thickness is fixed at 1.5 mm until the measured OD is available.
- Alumina density, conductivity, and heat capacity follow the older 0D material
  set.
- Aluminum flange conductivity is fixed at 205 W/m/K.
- The flange is treated as an isothermal 25 C sink because it is water cooled.
- T2 remains a measured cavity/felt boundary, not a predicted insulation state.

## New Diagnostics

`simulate_v8` returns:

- `rear_tube_temperature`
- `rear_tube_heat_transfer_coefficient`
- `receiver_to_tube_heat`
- `tube_to_t2_heat_loss`
- `flange_heat_loss`

These are intended to show whether the added rear path has the right magnitude
and timing before adding a predicted insulation/T2 state in a later v8b.
