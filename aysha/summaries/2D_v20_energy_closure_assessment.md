# v16 -> v20 assessment, and the missing front-aperture loss

Scope: `journal.2D.md` (v17--v20 entries), `2D_v20_theory_manual.md`,
`manuscript_2D_readiness_strategy.md`, `summaries/2D_v20/*`, plus a model-free
steady energy-closure diagnostic (`front_loss_closure_diagnostic.py`).

---

## Part 1 -- Is the model closer to manuscript level?

**No. On every observable that matters it is further away than v16.** What
improved is knowledge, not fit.

| held-out heating | v16 | v17 | v18 | v19 | v20 gate-free |
|---|---:|---:|---:|---:|---:|
| side-wall RMSE (K) | ~62* | ~97 | 82.9 | 94.0 | 94.4 |
| T3 RMSE (K) | -- | -- | 101.5 | 90.0 | 129.8 |
| axial RMSE (K) | 110.5 | 44.0** | 130.7 | 36.4** | -- |
| cooling side RMSE (K) | 33.8 | 35.7 | 26.6 (C81) | 46.0 (C81) | 43.7 (C81) |
| apparent `Nu` prefactor (meas. `3.05e-4`) | `4.91e-4` | `6.89e-4` | `6.11e-3` | -- | -- |
| apparent `Nu` exponent (meas. 1.445) | 1.412 | 1.442 | 0.793 | -- | -- |
| effectiveness range (meas. 0.573--0.781) | 0.732--0.856 | 0.877--0.937 | 0.712--0.889 | -- | -- |
| 58 mm peaks (meas. 10/15) | 0/15 | 0/15 | 0/15 | 7/10 diag. only | 0/15 |

\* v16 mean-sensor 61.16 K includes the easy T2. \*\* v17/v19 axial RMSE fell
only because the fitted power collapsed and flattened the profile; both still
predict 0/15 middle peaks.

The one configuration that came closest to the data is still **v16 with power
1.65/1.80/1.25 and `k_felt` scale 7.2** — the configuration that was rejected
and then explicitly reversed. Part 2 shows why that matters: v16's "nonphysical"
parameters were a crude but *directionally correct* proxy for a real missing
loss, and the corrective actions taken in v17--v20 (enforce the closure
brackets, restore `k_felt` to ~1) both moved away from the physics.

**What v17--v20 did achieve** is a clean, systematic elimination. The following
are now ruled out as the cause of the 40--180 K residual, each with a
boundary-active or null result:

- felt conductivity and heat capacity (v16 bound, v18 upper bound, v20 flat
  plateau at 0.05--0.15 — the felt path is not even rate-limiting);
- casing/backplate-to-flange contact (v17: any value from 20 to 500 W/K is
  equivalent — a clamped plateau);
- global exchange magnitude (v18: lower bound 0.05, and it destroyed the
  emergent `Re^1.44`);
- exchange exponent and integrated `UA` (v19/v20: a boundary ridge, no interior
  basin, and the gate-free search pushed `n` down to 0.25--0.50);
- axial source depth (v18--v20 all select `f_deep -> 0.9`, `L_deep -> 0.12 m`,
  and extending the bounds does not help);
- the T3 probe operator (v19/v20: capacity unbounded to 30000 J/m2/K, stem
  conductance pinned at zero, still 30--95 K RMSE);
- variable-`Cp` enthalpy bookkeeping (v20: worth 0.03--0.6 K);
- mesh, solver tolerance, mass/energy conservation (all verified).

The gate-free stress test is the decisive part: with acceptance gates removed
and only non-negativity retained, the search still leaves 94--116 K heating side
error. **A missing term, not a mis-set bound, is causing this.** Part 2
identifies it.

---

## Part 2 -- What you are missing

### 2.1 The residual has one dominant structure: it is monotonic in flow

From `gatefree_validation_cases_alltrain_power_2D_v20.csv`, final-window side
bias within each lamp group, ordered by decreasing flow:

| 456 kW/m2 | | 304 kW/m2 | | 256 kW/m2 | |
|---|---:|---|---:|---|---:|
| E67 (Re 72) | −37.6 K | E72 (Re 94) | −73.7 K | E77 (Re 79) | −29.1 K |
| E68 (59) | −21.4 | E73 (65) | −89.7 | E78 (56) | −33.9 |
| E69 (50) | +10.5 | E74 (45) | −58.9 | E79 (44) | −27.2 |
| E70 (43) | +11.9 | E75 (35) | −3.6 | E80 (36) | −9.8 |
| E71 (34) | **+97.5** | E76 (23) | **+135.4** | E81 (25) | **+52.6** |

The T3 bias follows the same sweep (E67 −45 K to E71 +328 K; E72 −123 K to E76
+370 K). **Too cold at high flow, too hot at low flow, swinging 130--210 K
inside a single lamp setting.** This is not a spatial-distribution error and not
a level error. It is the global steady energy balance responding wrongly to
flow.

### 2.2 The energy books do not close with the assumed loss law

The delivered power is fixed within a lamp group — the lamps do not know the
flow rate. So

```text
f = [ Q_gas(measured) + Q_loss ] / (Io A_front)
```

must be the *same* for all five cases in a group. It is not:

| loss law | 456: mean f, max/min | 304 | 256 |
|---|---|---|---|
| `K = 0.096 W/K` | 1.05, **1.52x** | 1.22, **1.95x** | 0.84, **1.49x** |
| `K = 0.164 W/K` | 1.34, **1.34x** | 1.58, **1.57x** | 1.11, **1.25x** |
| `sig(Tfront^4 - Ta^4)` | 1.19, **1.09x** | 1.54, **1.13x** | 0.94, **1.07x** |

The manuscript's slow-mode `K_loss = 0.10-0.16 W/K` is a *cooling* quantity,
measured at 300--600 K. **The independent power-closure brackets
(1.05--1.34 / 1.23--1.58 / 0.84--1.11) were built from it, so they inherit a
loss law that the heating data reject.** Those brackets have been enforced as
hard constraints since v18, which is why every power vector since then has sat
on its floor and every heating fit has been too cold at high flow.

### 2.3 The missing loss is at the front, and it is radiative

Fitting `Q_loss = sig_eff (T^4 - Ta^4)` and asking which temperature `T` makes
`f` constant:

| driver | 456 CV(f) | 304 | 256 | implied `eps*A` |
|---|---:|---:|---:|---|
| T8 (11 mm) | 3.4% | 2.0% | 1.0% | 10--31 cm2 |
| T12 (58 mm) | 3.0% | 6.2% | 3.0% | 22--42 cm2 |
| **T11 (107 mm)** | 4.9% | **12.1%** | 5.7% | 167--375 cm2 (absurd) |
| Tw (quadrature) | 2.9% | 5.4% | 2.7% | 30--58 cm2 |
| **front ≈ T8+150 K** | **3.1%** | **1.4%** | **0.6%** | **7--17 cm2** |

The rear is excluded outright, and for a good reason visible in the raw data:
across each group's flow sweep **T11 barely moves (32/70/52 K span) while T8
swings 264/331/194 K.** The rear is nearly clamped by the adaptor/flange; it
cannot be the source of a flow-dependent loss.

The front closure is robust to the assumed front-face offset. For offsets of
0--250 K above T8, `f` closes to 0.3--3.4% CV and lands **inside the independent
closure bracket in 11 of 12 combinations**. The required effective radiating
area is 5.5--31 cm2, i.e. **1.5--8.5x the 19x19 mm monolith face**, equivalent
to a hot front zone of roughly 26--63 mm across.

The implied magnitude:

| case | Re | Tfront | Q_abs | Q_front | share |
|---|---:|---:|---:|---:|---:|
| E67 | 72 | 1090 K | 196 W | 56 W | 28% |
| E71 | 34 | 1354 K | 196 W | 132 W | **67%** |
| E72 | 94 | 887 K | 169 W | 39 W | 23% |
| E76 | 23 | 1218 K | 169 W | 139 W | **82%** |
| E77 | 79 | 715 K | 87 W | 24 W | 28% |
| E81 | 25 | 909 K | 87 W | 64 W | **74%** |

**Your model can radiate at most ~29 W out of the front** (blackbody over
3.61 cm2 at 1090 K). The data require 24--139 W. The model is short by a factor
of 2--5 at high flow and 5--20x at low flow, on a term that scales as `T^4` and
therefore dominates exactly where the residual is largest.

### 2.4 Why this is physically right, and why it was excluded by rule

The illuminated front plane is larger than the 19 mm monolith. The felt face and
frame surrounding the receiver are directly irradiated, run hot, and re-radiate
to the cold surroundings. The receiver's own front face and channel mouths add
to this. The v10 rule set kept "front thermal radiation to the surroundings" but
applied it to the 3.61 cm2 receiver face only, and the readiness strategy's
prohibition on "side-weighted beam / arbitrary side heating" (§8.9) was written
against side heating of the *monolith*. **Neither rule was meant to exclude a
hot front plane wider than the monolith — but together they did, and that plane
is where 25--80% of the energy goes.**

This single term explains every standing symptom:

| symptom | explanation |
|---|---|
| too cold at high flow, too hot at low flow | the missing sink grows as `T^4`, so it bites hardest exactly at low flow |
| power factors run to their bracket floors | with no front sink, the only way to stop low-flow overheating is to remove power — which then starves high flow |
| axial maximum stuck at the entrance, 0/15 | there is no sink at `z=0` to push the maximum downstream |
| fits select `f_deep -> 0.9`, `L_deep -> 0.12 m` | a near-uniform source is the optimizer's substitute for a missing front sink; it moves the peak by starving the front instead of cooling it |
| v16's `k_felt = 7.2` and power 1.65/1.80/1.25 | the least-bad compensation available: buy extra loss on the wrong path. Note 1.65/1.80/1.25 vs the front-closure values 1.19/1.54/0.94 — v16 over-shot power precisely because its extra loss was put in a path with the wrong `T` dependence |
| felt `k` on a flat 0.05--0.15 plateau in v20 | the felt path was never the bottleneck; once `G_case_flange` clamped the casing, felt conductivity stopped mattering at all |
| effectiveness too high, range compressed | with too little front loss, too much of the absorbed power must go into the gas |
| T3 observer capacity runs away to 30000 J/m2/K | an observer lag is being asked to absorb a plant *level* error; it cannot, so it diverges |
| `eps* = 0.66` crossover never emerges | the crossover flow is set by the balance of front re-radiation against entrance convection — remove the first and there is no crossover |

### 2.5 Two supporting sanity checks

- **The mid-peak needs only a few watts.** At E72 the measured 11->58 mm rise
  implies ~6 W conducting forward through the SiC; a 750--890 K front face
  radiating out a 3.6 cm2 aperture supplies about that. You do not need a deep
  optical source to move the peak — you need a front sink.
- **The front face cannot be arbitrarily hot.** SiC axial conductance is
  0.035 W/K, so the `z=0` face can only sit ~100--230 K above T8 before the
  conduction back-flow becomes unphysical. That is why the closure needs an
  effective area of 7--17 cm2, not a 1600 K monolith face: **the extra
  radiating area is the surrounding illuminated front plane, not the monolith.**

---

## Part 3 -- What to do

1. **Confirm the diagnosis in one run (no refit).** For every heating case,
   print the v20 model's `z=0` face temperature, the total front radiative loss,
   and its ratio to absorbed power. If that ratio is below ~15% while the table
   in §2.3 requires 23--82%, the diagnosis is confirmed. This is a reporting
   change, not a model change.

2. **Add exactly one parameter and refit only it plus group power.** An
   effective hot-front radiating area `A_front,rad`, bounded 3.6--40 cm2, at the
   existing emissivity, radiating to ambient from a front-plane temperature tied
   to the receiver front cell. Freeze everything else at v16-era physics
   (Beer--Lambert source, `k_felt` scale ~1, v14 exchange law). Predictions to
   pre-declare:
   - `A_front,rad` lands interior at roughly 7--20 cm2;
   - group power lands inside 1.05--1.34 / 1.23--1.58 / 0.84--1.11;
   - the monotonic flow bias in §2.1 collapses;
   - the axial maximum moves toward 58 mm at high flow only, producing the
     `eps* ≈ 0.66` crossover rather than 0/15 or 15/15;
   - effectiveness range widens toward 0.573--0.781;
   - `f_deep` is no longer driven to its upper bound.
   If `A_front,rad` goes to a boundary or the flow bias survives, the front
   hypothesis is falsified cleanly and cheaply.

3. **Rebuild the power-closure brackets with the radiative loss law.** The
   current brackets are conditional on `K_loss = 0.096-0.164 W/K`, a cooling
   quantity. Recompute them as `f = [Q_gas + K_lin dT + sig(Tfront^4 - Ta^4)]/(Io A_front)`
   and report them as a *band over the plausible `sig`*, not a hard constraint.
   Until then they should not be used as an optimizer bound.

4. **Check the front-plane metrology.** How was `Io` measured — peak or aperture
   average, and over what aperture relative to the 19 mm receiver? How much felt
   face is exposed inside the beam? Those two numbers turn `A_front,rad` from a
   fitted effective area into a geometric one, which is the difference between
   Role A and Role B for this term.

5. **Reinstate v16's power range as admissible.** V16's 1.65/1.80/1.25 was
   rejected for violating brackets that §2.2 shows are built on the wrong loss
   law. It should be re-examined, not treated as a known-bad configuration.

---

## Part 4 -- Manuscript position

The honest current statement is unchanged from v20 and remains defensible:
report the directly measured apparent correlation plus a formal
non-identifiability statement. But the non-identifiability claim should be
qualified — v20 established that the *existing* closure set cannot be
identified, not that the receiver's macroscopic coefficients are
unidentifiable in principle. The front-loss term in §2.3 is a specific,
falsifiable, single-parameter candidate that has never been tested, and it
closes an energy balance the current model cannot. That distinction is worth
one more model version before writing the negative result.

Reproduce §2.2--2.3 with `python summaries/front_loss_closure_diagnostic.py`.
