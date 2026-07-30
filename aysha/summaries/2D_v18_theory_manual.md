# 2D v18 theory and fitting manual

## Purpose and status

V18 is a diagnostic successor to v17.  It asks whether the measured
side-wall, felt, and outlet-air histories can be reconciled by changing:

1. the magnitude, but not the Reynolds-number exponent, of local wall--gas
   exchange;
2. the axial deposition of absorbed lamp power; and
3. modest effective felt conductivity and heat capacity.

It is not permitted to fit the two interior thermocouples T9 and T10.  Their
wall/gas observation operator is not sufficiently known, and v17 showed that
changing that operator cannot repair the radial signs.

## Physical representation

The receiver retains all 100 square channels.  D4 symmetry groups channels
with identical geometry into 15 channel orbits; the orbit multiplicities
recover the full receiver mass, flow, heat-transfer area, and absorbed power.
This is not a 15-channel receiver or a quarter-receiver approximation.

Each orbit contains axial SiC wall states and a quasi-steady gas march.  At
each time step:

- the measured MFC flow is converted to mass flow at the MFC reference
  condition;
- equal pressure drop allocates that mass among the 15 channel orbits,
  including the measured outer-groove obstruction model;
- gas density and thermophysical properties depend on local temperature;
- gas temperature is marched through each channel using local wall--gas
  exchange; and
- the resulting gas heat gain is removed from the solid energy equation.

The surrounding assembly contains the thin side-surface state, filled alumina
silicate felt, aluminum casing, solid alumina adaptor, rear alumina tube,
metal/water-cooled flange, ambient losses, and the verified rear casing-to-
flange contact.  The front side thermocouple is sampled at 11 mm; the other
side thermocouples are at 58 and 107 mm.

## New v18 hypotheses

### Exchange magnitude

V17's apparent Nusselt exponent, 1.442, already agrees with the experimental
reduction, 1.445.  V18 therefore preserves the flow dependence and applies
one multiplier `m_h` to the inherited exchange coefficient:

```text
h_v18(Re,T) = m_h h_v17(Re,T).
```

This is a structural diagnostic.  A value far below unity means that the
v17 local core--gas balance transfers energy much more strongly than the
macroscopic measurements permit; it is not automatically a new validated
Nusselt correlation.

### Axial source

The baseline is the inherited Beer--Lambert deposition.  The alternative is
a conservative near/deep mixture:

```text
w_gj = W_g [(1-f_deep) w_near,j + f_deep w_deep,j(L_deep)].
```

The total for every radial orbit `W_g` is unchanged.  Thus the alternative
cannot create power, alter the centered radial beam, or imitate a power
factor.  Bounds are `0.05 <= f_deep <= 0.80` and
`0.020 <= L_deep <= 0.200 m`.  Because both terms decrease from the front,
this hypothesis can produce a 58 mm temperature maximum only through the
interaction of deposition, axial transport, and end losses.

### Felt

The nominal material properties remain the reference.  After structural
selection, conductivity and heat capacity are screened over
`0.70, 1.00, 1.30` times nominal using C69 and C80 cooling records.  The
felt/receiver contact multiplier remains fixed at 0.30 because contact and
conductivity were structurally confounded in v16.  C81 is untouched
validation.

## Numerical mesh

V18 replaces the older non-nested axial meshes with:

| mesh | felt radial | casing radial | receiver axial | front/rear axial | rear tube |
|---|---:|---:|---:|---:|---:|
| C, screen | 3 | 1 | 24 | 4/20 | 15 |
| M, nominal | 6 | 2 | 48 | 8/40 | 30 |
| F, refined | 12 | 4 | 96 | 16/80 | 60 |

All counts double from C to M to F.  The 14 receiver auxiliary radial entries
remain fixed because physical receiver cross-section is already the 100-
channel/15-orbit topology.

The pre-fit M-to-F test gives:

| case | history RMS (K) | history maximum (K) | final maximum (K) |
|---|---:|---:|---:|
| E72 | 5.51 | 9.36 | 7.24 |
| C69 | 9.53 | 19.93 | 13.15 |
| C81 | 6.04 | 12.16 | 7.53 |

The declared history gates of 10 K RMS and 20 K maximum pass, narrowly for
C69.  Candidate improvements smaller than this numerical scale are not
treated as mechanistic evidence.

## Calibration data and objective

Heating uses interleaved training and validation within every power level:

```text
training: E67/E69/E71, E72/E74/E76, E77/E79/E81
held out: E68/E70, E73/E75, E78/E80
```

Selection observables are T8, T12, T11, T3, and T2.  T9/T10 are diagnostic
only.  Histories use 61 fractional-time samples and exclude the initialized
first point.  Residuals use

```text
rho(x) = 2 [sqrt(1+x^2) - 1]
```

with discrepancy scales 35 K for each side TC, 25 K for T3, and 15 K for T2.
The heating objective is

```text
J = 0.35 J_curve + 0.25 J_final-level
  + 0.20 J_side-shape + 0.20 J_effectiveness.
```

Curve weights are 0.50 side, 0.30 air, and 0.20 felt.  Final-level weights
are 0.45 side, 0.35 air, and 0.20 felt.  Effectiveness is calculated from the
exact 11/58/107 mm wall quadrature and T3.  Cooling property selection uses
0.40 side, 0.30 air, and 0.30 felt curve losses; it does not use a solar
effectiveness term.

One power multiplier for each lamp group is profiled inside every structural
candidate and is hard-bounded by independent power closure:

```text
456 kW/m2: 1.05--1.34
304 kW/m2: 1.23--1.58
256 kW/m2: 0.84--1.11
```

The broad screen runs on C, with all nine training cases used for nested
power profiling.  The four mechanism combinations are transferred to M using
the interleaved E69/E74/E79 representatives, one per lamp group, because
simultaneously solving all nine stiff nominal cases caused severe thread
contention without adding power-group coverage.  All records still enter
final validation.  Only the winner and runner-up are eligible for F
confirmation.

## Interpretation rules

- A power factor at a closure bound is reported as unresolved power
  accounting; it is not silently widened.
- An exchange multiplier at the diagnostic lower bound means the present
  local exchange formulation is rejected, even if aggregate RMSE improves.
- A deeper source is supported only if it improves held-out primary
  observables by more than mesh uncertainty and improves axial peak behavior.
- Felt changes are accepted only if the cooling-training improvement transfers
  to C81 without degrading heating.
- Passing temperature parity alone is insufficient.  V18 must also reproduce
  outlet effectiveness, apparent Nusselt magnitude and exponent, axial
  inversion, mass/power/energy closure, and mesh transfer.

## Final result and decision

The C-mesh marginal screen preferred `m_h=0.05`, the extended lower bound.
The Beer--Lambert source was better than every tested near/deep source.  After
power profiling, all three power factors selected their independent lower
bounds.  The preliminary cooling profile preferred the upper tested felt
conductivity and heat-capacity factors, but its outlet-air loss remained much
larger than its side and felt losses.

Nominal reranking confirms Beer--Lambert/`m_h=0.05` as the scalar-objective
winner (`J=3.087`) and the mild near/deep/`m_h=0.05` candidate as runner-up
(`J=3.135`).  Both use the lower-bound power vector.  The ranking is stable
on F, but the parameter regime is not converged: winner M-to-F history
changes are 12.91--20.12 K RMS and up to 32.86 K maximum.  The declared
10/20 K gates fail.

Untouched M-mesh validation gives:

| phase | primary-five mean RMSE | side RMSE | T3 RMSE | T2 RMSE | axial-inversion RMSE |
|---|---:|---:|---:|---:|---:|
| heating training | 84.66 K | 93.98 K | 132.81 K | 8.56 K | 129.54 K |
| heating held out | 71.50 K | 82.90 K | 101.53 K | 7.26 K | 130.65 K |
| cooling training | 33.25 K | 33.83 K | 50.12 K | 14.64 K | 2.56 K |
| C81 held out | 26.57 K | 29.31 K | 36.21 K | 8.71 K | 6.63 K |

The M-mesh energy residual is below `1.0e-13 W` and source-power error is
zero.  The inherited detailed ledger is intentionally disabled on C/F because
their new axial faces differ from v17; source-power invariants remain exact
on every mesh.  These conservation passes do not compensate for the heating
failures.

The manuscript invariant evaluator gives:

| invariant | measured | v18 selected |
|---|---:|---:|
| apparent-`Nu` prefactor | `3.047e-4` | `6.107e-3` |
| apparent-`Nu` exponent | 1.445 | 0.793 |
| effectiveness range | 0.573--0.781 | 0.712--0.889 |
| effectiveness RMSE | -- | 0.188 |
| middle peaks | 10/15 | 0/15 |

Thus the global multiplier improves one composite objective only by destroying
the previously correct emergent flow exponent.  It is rejected.

The full weak-exchange interaction audit changes the optical interpretation.
The candidate `f_deep=0.60`, `L_deep=50 mm`, `m_h=0.10` creates all three
representative middle peaks with mean inversion error `-2.75 K`, but its
objective is 5.06 versus 3.18 for Beer--Lambert.  A 15-case C-mesh sign audit
finds 7/10 true peaks and five false peaks; heldout axial RMSE is 51.02 K,
side RMSE 124.65 K, and T3 RMSE 130.46 K.  Deeper deposition is therefore
capable of the required shape but is not a transferable complete model.

V18 is rejected for Role B coefficient extraction.  Its fitted exchange,
felt, and power values are diagnostic boundary values, not physical
coefficients.  V19 should retain the 15 D4 orbits but identify normalized
optical shape, integrated wall--gas conductance, and the downstream T3
observation operator as separate invariant-constrained stages.
