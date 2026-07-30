# Independent assessment of 2D_v16 and the post-v16 strategy

Scope: `journal.2D.md`, `2D_v16_theory_manual.md`, `manuscript_2D_readiness_strategy.md`,
and a first execution of the Stage-0 invariant evaluator (`invariant_evaluator_v0.py`).

Bottom line: **the v16 rejection is correct, but the stated next-step diagnosis is
only about half right, and the ordering of the staged plan is wrong.** Three
findings below change what v17 should be.

---

## 1. What v16 actually established (agreed)

- The felt-path fit is structurally non-identifiable, not merely unlucky (§3).
- The absolute power level was miscalibrated after the loss path changed; the
  refit was necessary and legitimate (61 K vs 125 K held-out).
- The spatial failure survives both the power refit and the T8 11 mm correction,
  so it is a model-form failure, not a calibration or coordinate artifact.
- Mesh transfer fails (E72: 36.61 K screen-to-nominal).

The v9-v16 audit discipline — pre-declared acceptance, held-out splits, boundary
reporting, explicit rejection — is the strongest part of this work and should be
kept unchanged.

---

## 2. Stage-0 evaluator, run now against v16

I ran the manuscript's own reduction on v16's predicted steady states (model =
measured + `steady_error_K`, so no re-solve was needed). This is the check the
strategy defers to Stage 0; it takes minutes and it reorders everything.

| Invariant | Measured | v16 model | Verdict |
|---|---|---|---|
| I1 apparent `Nu(Re)` | `3.05e-4 Re^1.445`, r2=0.970 | `4.91e-4 Re^1.412`, r2=0.996 | **exponent already passes**; prefactor +61% |
| I2 effectiveness | 0.573 – 0.781 | 0.732 – 0.856 | biased high by 0.07–0.16; **range compressed 2x** |
| I2 axial inversion | mid-peak 10/15, `eps*≈0.66` | 0/15, `T12-T8` = −34 to −194 K | fails, offset ≈ −100 K at matched `eps` |
| I3 `Lambda_58` | `+0.032 + 2.3e-4 Re` | `−0.034 + 8.3e-6 Re` | **sign wrong**, 0/15 |
| I3 `Lambda_107` | `+0.038 + 8.25e-4 Re`, r2=0.90 | `−0.032 − 1.8e-4 Re`, r2=0.47 | **sign and slope wrong**, 0/15 |

Two consequences:

**(a) I1 is not missing physics.** The readiness strategy treats the superlinear
exponent as a signature of flow-dependent solid recruitment that must still
"emerge". It already emerges from the existing v14–v16 network, at
`Re^1.412` against a measured `Re^1.445`, with a *better* correlation than the
data. The recruitment mechanism is present; only its magnitude is 60% too large.
The I1 gate in §8 of the strategy should be restated as a *magnitude* test, and
its current framing overstates how much new physics is required.

**(b) The axial failure decomposes into two independent errors.** Regressing
`T12-T8` on `eps`: the model is offset by roughly −100 K *and* about twice as
sensitive to flow as the data (≈1600 K vs ≈765 K per unit `eps`). Only the offset
is a source-shape problem. The excess slope, the +0.10 `eps` bias, and the +61%
`Nu` prefactor are one error: **the model over-exchanges heat into the gas, so the
gas saturates too early**, which pins the temperature maximum at the entrance in
every case and compresses the flow response. A two-component axial source fixes
the offset and leaves the slope wrong. v17A as specified will therefore improve
axial RMSE and still fail I2.

---

## 3. Why the v16 felt fit ran to its bound (mechanism, not bad luck)

During cooling the receiver is essentially isothermal (`Bi < 1e-5`; 137 mm of SiC
equilibrates fast). A near-isothermal relaxation of a lumped body against ambient
identifies exactly two combinations: a total capacitance and a total series
conductance. V16 fitted three parameters — `S_contact`, `S_k`, `S_C` — of which
`S_contact` and `S_k` appear *in series* in the same `R_series`. That is a
one-dimensional null space by construction, so one of the pair must run to a
boundary and the reported `k_felt` scale of 7.2 carries no information.

This should be stated in the theory manual as a structural result. It also means
the strategy's Stage 3 ("fit global contact/loss quantities to cooling training
cases") will reproduce the same failure unless the felt path is
**re-parameterised to the two identifiable combinations** (`C_total`,
`UA_series`) and material properties are then back-derived with a stated
uncertainty. Fitting three series-coupled scales to a lumped decay is not
recoverable by extending bounds or changing the mesh.

Corollary: v16's cooling axial RMSE of 8.00 K is nearly uninformative. An
isothermal body has no axial signal to get wrong. Cooling data cannot validate
spatial physics and should not be counted as if they do.

---

## 4. The radial signs are almost certainly an observation-model error

This is the finding I would act on first.

The measured "radial" gaps are `T12-T9 ≈ +24 K` at 58 mm and `T11-T10 ≈ +36 K`
(up to +53 K) at 107 mm. Both plans treat these as solid-to-solid radial
temperature differences and propose hydraulic mechanisms to create them.
**A monolithic SiC block cannot sustain them.**

For the 19 mm, 100-channel monolith:

| direction | conducting area | path length | conductance at k=35 W/m/K |
|---|---:|---:|---:|
| axial | 1.360e-4 m2 | 137 mm | 0.0347 W/K |
| radial | 5.480e-4 m2 | 9.5 mm | 2.019 W/K |

The radial conductance is **58x** the axial conductance. So:

- the measured 114 K axial spread costs ≈4 W — entirely affordable out of ~165 W;
- a 50 K radial gap would require **≈101 W of radial conduction, ~61% of total
  input**, flowing *inward* from the felt-cooled skin toward the core.

That is not physically available. The correct statement is that the solid must be
radially uniform to within roughly 5 K, which is exactly what v11 predicted
(0.25–0.70 K) and why v11's rear-groove hydraulics, with a real 2.2–3.3x
core/edge mass-flux contrast, produced no radial solid signal. V16 still runs at
a core/side flow ratio of 2.63–2.69 and still gets 0/15.

Six model versions have now been spent trying to create a solid-phase gradient
that the geometry forbids. The alternative reading is immediate and consistent:

- the manuscript itself labels `Lambda` as **deep LTNE**, i.e. solid-minus-**gas**;
- gas is colder than the wall (`eps < 1`), giving the correct positive sign;
- the deficit grows with flow because the gas saturates less — measured
  `Lambda_107` slope `+8.25e-4 Re`, r2=0.90, is a textbook LTNE signature;
- the model's sign is negative simply because it observes T9/T10 as *interior
  solid*, which sits above the felt-cooled exterior skin.

**Recommended immediate test (post-processing only, no refit, no new physics):**
re-observe T9/T10 in the existing v16 solution as the local gas temperature, or
as a blend `a*Tgas + (1-a)*Tsolid` with a single global `a`, and re-run the
evaluator. If the sign flips and `Lambda_107(Re)` acquires a positive slope, the
entire radial branch of the strategy — including Stage 2 / v17B — is unnecessary.
This costs one afternoon and could retire the largest open failure in the model.

This is the same class of error as the T8 5 mm coordinate (8 versions) and the
pre-v15 side-thermocouple state (6 versions). The recurring root cause in this
project is the observation operator, not the transport physics.

---

## 5. Assessment of the staged plan in `manuscript_2D_readiness_strategy.md`

**Correct and worth keeping:**

- the D4 15-orbit argument (§2) is sound; the 100-channel rewrite is not justified;
- Role A / Role B separation and the nine-item gate (§8);
- nesting power factors as nuisance parameters inside each source candidate;
- refusing to publish a fitted axial deposition law as an optical coefficient;
- the information-assignment table (§7).

**Where I disagree:**

1. **Stage 0 is listed as a prerequisite but was not run before writing the plan.**
   Had it been, §1 would not describe I1 as unachieved and Stage 1 would have been
   specified differently. Run the evaluator against v11, v14, v15 and v16 before
   committing to v17 — the version-over-version trend in `eps` and `Lambda` is the
   cheapest information available and none of it has been extracted.

2. **Mesh convergence is an acceptance item (§8.1) when it must be a
   precondition.** A 36.61 K screen-to-nominal change on E72 is comparable to the
   effects being discriminated. Ranking two axial source laws on a mesh with that
   error is not a valid experiment. Establish a converged baseline (grid
   refinement + Richardson estimate on the axial and radial directions) *before*
   Stage 1, not at Stage 4.

3. **Stage 1 targets the wrong knob alone.** Per §2(b), source-shape redistribution
   cannot fix the `eps` bias or the excess flow sensitivity. Stage 1 should screen a
   2x2: {single Beer–Lambert, two-component source} x {current exchange law,
   exchange law recalibrated to reproduce measured `eps(Re)`}. Rank on the
   invariants, not on RMSE. Without the second axis, a "passing" source law will be
   compensating an over-strong exchange coefficient.

4. **Stage 1 substep 2 (physical axial radiosity) is likely a dead end.** With
   `N_rc = 4 sigma T^3 D_h / k_g` (which reproduces the measured 1.6 at 600 K and
   5.2 at 1000 K to within 4% using `L = D_h = 1.5 mm` — the definition is
   confirmed), the equivalent axial radiative conductance through the open channel
   area at 900 K is ~5.6e-5 W·m/K against ~5.4e-3 W·m/K for axial SiC conduction:
   **about 1%**. In-channel axial radiation is real but two orders of magnitude too
   weak to redistribute the axial profile. `N_rc = O(1-5)` says radiation competes
   with *gas* conduction across the bore, which matters for the wall-to-wall
   coupling and hence the effective local exchange — not for axial transport. This
   argues for folding radiation into the *exchange* law, which is also what §2(b)
   demands, rather than into an axial radiosity network.

5. **Stage 2 (v17B rear manifold) should be deferred, probably cancelled.** v11
   already falsified flow maldistribution as the radial mechanism with a 2.2–3.3x
   contrast, and §4 shows why no hydraulic contrast can work. Do not implement it
   until the T9/T10 observation question is settled.

6. **Sequential freeze-then-advance has now failed three times** (v14 felt
   multiplier frozen into v15; v14 power multipliers frozen into v16; v9 optical fit
   confounding the v11 axial verdict — journal line 1555). The plan repeats the
   pattern. At minimum, carry a bounded loss/contact nuisance parameter through
   Stage 1 rather than freezing v16's rejected felt path into it.

7. **The refitted power factors violate every closure bracket** and the strategy
   does not say so: 1.65 vs 1.05–1.34; 1.80 vs 1.23–1.58; 1.25 vs 0.84–1.11. All
   three are above their upper bounds, two above even the most extreme single-case
   value. This is a closed compensating loop — `k_felt` at 7.2 over-drains the
   assembly during heating, and the power factor is raised to compensate. It is
   direct evidence that §8.3 (parameter interiority) and §8.6 (power accounting)
   fail *together*, from a single cause. A cheap confirmation: rerun the v16 power
   refit with `S_k` forced to a physically plausible value near 1.0. If the group
   factors fall back inside their brackets, the compensation is proven and v17
   inherits a clean baseline.

---

## 6. Recommended revision to the immediate plan

Replace §9 of the readiness strategy with:

1. **Run the evaluator on v11/v14/v15/v16.** (Prototype committed as
   `summaries/invariant_evaluator_v0.py`; extend to `Lambda_58`, `C_eff`, `K_loss`,
   master-curve collapse and energy closure.)
2. **Settle the T9/T10 observation operator** by re-observation of the existing
   v16 fields. Decide gas / blend / solid on evidence, and record it.
3. **Establish mesh convergence** on the current model. Everything downstream is
   conditional on this.
4. **Re-parameterise the felt path** to `(C_total, UA_series)` and re-run the v16
   power refit with a plausible `k_felt`. Confirm or refute the
   `k_felt`/power compensation loop against the closure brackets.
5. **Then v17A**, as a 2x2 screen over {axial source law} x {exchange law},
   ranked on I1 prefactor, `eps(Re)` range, and the `eps*≈0.66` crossover — not on
   RMSE.
6. **Hold v17B** unless step 2 shows the radial gaps really are solid-phase.

Steps 1–4 use no new physics, no new experiments, and no solver changes beyond
re-observation. They are worth doing before any further model version.

---

## 7. Note on the third-party plan

The post-v16 plan from the other model was not present in `analysis/` or
`summaries/` and was not attached, so it is not assessed here. Sections 2–6 apply
to any plan whose next step is an axial source law followed by a rear-manifold
hydraulic network; if that plan differs materially, send it and I will place it
against the same evidence.
