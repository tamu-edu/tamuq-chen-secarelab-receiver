# 2D model code and thermophysical-property audit (v11 -> v20 chain)

Method: read the active inheritance chain (`2D_v11.jl` -> `2D_v12.jl` -> `2D_v14.jl`
-> `2D_v15.jl` -> `2D_v16/17/18/19/20.jl`) and evaluated every property
correlation and loss term analytically against literature and against the
reference tables used by `analysis/exp_analysis/dimensionless_analysis.py`.
Julia is not available in this sandbox, so magnitudes below are computed from
the code's own formulas and default parameters, not from a solver run. The
ledger check in §5 should confirm them.

**Headline: the front-loss diagnosis is confirmed in the code, and three
independent errors all push the same way.** Two parameters that were reported
as pathological fits (`felt k scale = 7.2`, `felt Cp scale = 0.55-0.75`) turn
out to be the optimizer correctly compensating wrong priors.

---

## 1. Front loss — confirmed structurally too small

| # | finding | location |
|---|---|---|
| 1 | Receiver front radiation uses `area = multiplicity * solid_area_per_channel` = **1.36 cm2**, not the 3.61 cm2 front face. The 100 channel mouths (2.25 cm2, **62% of the face**) radiate nothing. | `2D_v14.jl:742-750`, area defined `2D_v14.jl:160-161` |
| 2 | `front_emissivity = 0.85`. At 935 K this gives **5.0 W**. If the mouths radiated at a deep-cavity `eps ~ 1` over 3.61 cm2 it would be 15.5 W. | `2D_v11.jl:120` |
| 3 | `spillage_fraction = 0.10` is **carved out of** `delivered`, not added to it. The rim can therefore never absorb more than the receiver footprint implies. | `2D_v11.jl:132`, `2D_v14.jl:730-734`, deposition `2D_v14.jl:819-829` |
| 4 | `delivered = absorbed_fraction * I(t) * receiver_width^2` with `absorbed_fraction = 1.0`. The illuminated front plane is **pinned at 3.61 cm2**. | `2D_v14.jl:730-732` |

The steady energy balance (`front_loss_closure_diagnostic.py`) requires
**24-139 W** of front-region loss, 23% of absorbed power at high flow rising to
82% at low flow. The code can supply at most ~5-25 W. A 137 mm x 1.5 mm bore
array is close to a blackbody cavity at its mouth; excluding it is not a
conservative approximation, it removes the dominant `T^4` sink.

---

## 2. Felt conductivity — a temperature-SHAPE error, not a magnitude error

`felt_conductivity_v12` (`2D_v12.jl:148-164`) interpolates the RS-3000 sheet:
0.050 / 0.080 / 0.110 / 0.170 / 0.260 / 0.320 W/m/K at 20/500/800/1100/1400/1600 C.

**Ceramic-fibre datasheet conductivity is quoted against MEAN temperature**
(hot face at T, cold face near ambient). The code applies it pointwise at
**local** temperature. That systematically under-predicts hot-face `k`:

| T | v12 prior | v11 law `0.06 + 1.2e-10 T^3` | typical 128-140 kg/m3 blanket |
|---:|---:|---:|---:|
| 400 C | 0.074 | 0.097 | 0.070 |
| 600 C | 0.090 | 0.140 | 0.110 |
| 800 C | 0.110 | **0.208** | **0.200** |
| 1000 C | 0.150 | **0.308** | **0.310** |

The **v11 law was correct** and v12 replaced it with a curve that is right at
~400 C and a factor ~2 low above 800 C. The journal records this as a
"felt-property correction" (`journal.2D.md:1866-1888`); it was a regression.

Consequences, quantified:

- modelled felt lateral conductance (r = 10.7 -> 57 mm, L = 137 mm) is only
  **0.036-0.039 W/K**;
- the steady energy balance requires a total heating-temperature loss
  conductance of **0.21-0.36 W/K**;
- v16's rejected `k_felt scale = 7.2` gives **0.26-0.28 W/K** — *exactly the
  required value*.

**v16's "nonphysical" 7.2 was buying the missing conductance through the only
knob the code exposes.** And a scalar `scale` cannot repair a wrong *shape*: it
inflates the 300 K end (breaking cooling) while still under-supplying 900-1200 K.
That is the textbook cause of the boundary-active pattern — 7.2 at the ceiling
(v16), 1.30 at the ceiling (v18), a flat 0.05-0.15 plateau (v20 stress).

---

## 3. Thermophysical property errors

### 3.1 Air — both errors inflate the apparent Nusselt number

| property | code | error vs the reference table used by the manuscript reduction |
|---|---|---|
| `air_conductivity` | `2.414e-3 * sqrt(T)/(1 + 245.4*10^(-12/T)/T)` | **7.6 - 9.4% LOW** across 300-1200 K |
| `air_heat_capacity` | `1004*(1 + 1.983e-4 T - 4.14e-8 T^2)` | **+5.5 to +6.2% HIGH** at 300-600 K; +0.7% at 1200 K |
| `air_viscosity` | Sutherland | correct, <= 2.2% |

The conductivity coefficient is wrong: the standard correlation is
`k = 2.646e-3 T^1.5/(T + 245.4*10^(-12/T))`. Substituting **2.646e-3 for
2.414e-3 reproduces the reference table to within 1.3% at every temperature.**
`2.414e-3` appears to be a transcription of the water-viscosity constant.

Because apparent `Nu = h Dh/k` and `h ∝ NTU * mdot * cp`, the two errors
compound: **the model's apparent Nu is biased +14 to +17% high at low
temperature from properties alone.** The I1 acceptance gate is "prefactor
within 20%". Most of that budget is consumed before any physics is evaluated.
The 256 group is worst affected — its outlet gas sits at 525-560 K, exactly
where the `cp` error peaks.

### 3.2 SiC

| property | code | vs dense alpha-SiC |
|---|---|---|
| `sic_conductivity` | 116.6 (300 K) ✓ | **+22% at 800 K, +39% at 1000 K, +63% at 1200 K, +84% at 1500 K** |
| `sic_heat_capacity` | 1150 (300 K) | **+67% at 300 K**, +23% at 500 K, +13% at 800 K, +9% above 1000 K |
| `density = 2150 kg/m3` | with nominal 0.4 mm webs | implies **33% porosity** vs dense 3210 |

The `cp` error is largest exactly where the transient is fastest. Receiver
capacity at 300 K is 46 J/K in the model versus ~27.6 J/K correct; at 800 K,
49.2 versus ~43.7 J/K. This is a direct contributor to the persistent
810-893 s `t90` MAE that no version has moved.

The density point matters for conduction, not mass. The model reconciles the
measured 40.06 g correctly, but then applies **full dense-SiC conductivity over
the full nominal web area**. If the monolith really is 33% porous, `k` should
also fall by ~(1-P)^1.5 ~ 0.55; if instead the webs are ~0.26 mm, the
conduction area is 1.5x too large. **Either way the axial conductance is
1.5-1.8x too high.**

### 3.3 In-channel radiative conductivity — beta is ~10x too small

`rad_extinction_coeff = 50.0 1/m` (`2D_v11.jl:95`) feeds
`k_rad = 16 sigma T^3/(3 beta)` (`2D_v11.jl:1035-1036`, `2D_v14.jl:694-701`):

- `k_rad(900 K) = 4.41 W/m/K`, `k_rad(1200 K) = 10.45 W/m/K`;
- `beta = 50 1/m` is a **20 mm radiative mean free path**;
- in a 1.5 mm bore array the photon mean free path is of order the bore, so
  `beta ~ 300-700 1/m` and `k_rad(900 K) ~ 0.55 W/m/K`.

**`k_rad` is therefore ~8x too large.** Combined with §3.2:

| axial conductance at 900 K | value |
|---|---:|
| model `(k_SiC + k_rad,50) * A/L` | **0.0680 W/K** |
| physical `(0.55 k_SiC + k_rad,400) * A/L` | **0.0355 W/K** |

**~1.9x too high.** This smears the axial profile the entire v18-v20 campaign
was trying to fit, keeps the front face artificially cool, and therefore
further suppresses the already-too-small front radiation of §1. It compounds
rather than offsets.

### 3.4 Felt heat capacity — the fitted "lower bound" is the correct value

`felt_heat_capacity = 1360 J/kg/K`, constant (`2D_v11.jl:98`). Alumina-silicate
fibre is 750-850 at 300 K and 1050-1130 at 1000 K: **+20 to +60% high.**

The scales selected repeatedly in v19/v20 (0.75, 0.75, 0.55) give
**748-1020 J/kg/K — the literature value.** The felt-Cp lower-bound result is
not a pathology; it is the optimizer correctly detecting a wrong prior. It
should be read as a property correction, not as a failed identification.

### 3.5 Nusselt asymptote — wrong boundary condition

`fully_developed_nusselt = 3.61` (`2D_v11.jl:114`). 3.61 is the square-duct
**constant-heat-flux (H1)** asymptote. V20's cell integration uses a
**constant wall temperature** (`2D_v20_theory_manual.md §3.2`,
`UA/mdot = int Cp/(Twall - T) dT`), whose square-duct asymptote is
**Nu_T = 2.98**. A **+21% over-statement**, again in the direction of
over-exchange and inflated effectiveness.

### 3.6 Properties that check out

- `alumina_conductivity`: 37.8 / 12.1 / 6.8 W/m/K at 20 / 500 / 1000 C — correct
  for dense 99% alumina. Adaptor mass/volume gives 3900 kg/m3, consistent.
- `alumina_heat_capacity`: 746 (300 K), 1151 (1000 K) — correct.
- Aluminium: `k = 205`, `rho = 2700`, `cp = 900` — all fine.
  `casing_emissivity = 0.20` is defensible for a rough/oxidised sleeve but 2-5x
  high if the sleeve is machined bright; worth one look at the hardware.
- `air_viscosity` (Sutherland) — correct.
- `SIGMA = 5.670374419e-8`, `AIR_SPECIFIC_GAS_CONSTANT = 287.05` — correct.
- `standard_pressure = 101.4 kPa` vs the Aalborg 101.325 kPa — 0.07%, negligible.

---

## 4. Other structural observations

**4.1 The "centered beam" is effectively uniform.** `beam_radius_sigma = 14 mm`
against a receiver half-width of 9.5 mm and a corner radius of 12.1 mm: the
Gaussian only falls to **0.69 at the corner**. Both strategy documents treat
"centered vs side-weighted illumination" as a live hypothesis and forbid
side-weighting, but the implemented radial field is nearly flat and carries
almost no discriminating information. This should be stated explicitly.

**4.2 The baseline source IS strongly front-loaded** —
`front_deposition_fraction = 0.20` plus `beta_opt = 110 1/m` puts roughly
**47% of the power in the first ~3.75 mm** on the C mesh. That is physically
right. The problem is entirely downstream of it: §3.3 conducts it away and §1
fails to radiate it out. The `f_deep -> 0.9`, `L_deep -> 0.12 m` selections in
v18-v20 are the optimizer flattening a correct source to compensate.

**4.3 `radial_conductivity_scale = 0.05`** (`2D_v11.jl:93`). For a square
honeycomb the transverse solid line-fraction is `1 - 10*1.5/19 = 0.211`. 0.05
is ~4x smaller and undocumented. It enters only the half-pitch term of the
receiver/felt series resistance, where it is not rate-limiting, so the impact is
small — but it should be justified or set to the geometric value.

**4.4 `rear_adaptor_conductance = 0.448 W/K`** (`2D_v11.jl:124`, documented as
"100 W/m2/K over channel-wall overlap") is a large fixed rear sink that has
never been identified. Measured T11 is nearly clamped across each flow sweep
(spans of 32 / 70 / 52 K for the 456 / 304 / 256 groups, against T8 spans of
264 / 331 / 194 K), so this parameter is doing substantial unexamined work at
exactly the location where the data are least informative.

---

## 5. Ordered fix list

Corrections 1-4 are unambiguous property/formulation errors and should be made
before any further identification. 5-6 are the structural change.

| # | change | expected effect |
|---|---|---|
| 1 | `air_conductivity`: `2.414e-3 -> 2.646e-3` | removes a 8-9% bias; model and manuscript reductions become comparable |
| 2 | `air_heat_capacity`: refit to the reference table (or use it directly) | removes +5-6% gas-enthalpy bias at 300-600 K |
| 3 | `felt_conductivity`: restore the v11 `0.06 + 1.2e-10 T^3` form, or convert the RS-3000 mean-temperature curve to a local-temperature curve | lateral conductance rises ~2x at the hot face; `k_felt` scale should then sit near 1.0 instead of at a bound |
| 4 | `felt_heat_capacity`: 1360 -> ~1000 J/kg/K (or a T-dependent curve); `sic_heat_capacity`: refit to 690/960/1090/1140/1180; `fully_developed_nusselt`: 3.61 -> 2.98 | felt-Cp and SiC-cp scales become interior; `t90` bias should move |
| 5 | `rad_extinction_coeff`: 50 -> bore-scale (300-700 1/m), and reconcile `density = 2150` with the conduction area/conductivity | axial conductance halves; the axial profile stops being smeared and the front face gets hotter |
| 6 | Front radiation: extend the radiating area to the **full 3.61 cm2 front face** with a deep-cavity effective emissivity for the mouths, **and** add a bounded hot-front-plane area `A_front,rad` (3.6-40 cm2) that is *additional* to `I*A_front`, i.e. decouple `spillage` from `delivered` | supplies the missing 24-139 W `T^4` sink; the monotonic flow bias should collapse |

After 1-5 alone, re-run the v16 configuration. Prediction to pre-declare: the
`k_felt` and `felt_Cp` scales move to the interior near 1.0, and the group
power factors fall from 1.65/1.80/1.25 toward the front-closure values
1.19/1.54/0.94. If they do, the property corrections and the front-loss term
are the same finding seen from two directions, and the v16-v20 boundary-active
history resolves.

## 6. Caveat

The property comparisons are exact (analytic evaluation of the code's own
formulas against literature and against the reference tables in
`dimensionless_analysis.py`). The loss-conductance and axial-conductance
magnitudes are computed from the code's formulas with default parameters and
nominal geometry, not from a solver run. Before acting on §5 item 6, print the
`energy_rate_ledger2D` breakdown (`front_loss`, `casing_loss`, `flange_loss`,
receiver front-face temperature) for each heating case and confirm that
`front_loss / delivered` is below ~15%.
