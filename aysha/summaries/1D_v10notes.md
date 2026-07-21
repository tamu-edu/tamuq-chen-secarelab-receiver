# 1D v10 Notes - Physics Refactor Toward Coefficient Identification

## Executive Position

v10 should be treated as a fundamental model refactor, not as another fitted
variant of v8. The previous versions have been useful because they exposed the
sign pattern of the mismatch:

- the front solid temperatures can be captured more easily than the deeper
  receiver and gas temperatures,
- deeper receiver response is experimentally slower and less flow-dependent
  than the model tends to predict,
- the model has difficulty matching solid and gas steady states at the same
  time,
- T2/cavity behavior improved when the rear domain was represented physically,
  which suggests that part of the previous mismatch was structural rather than
  purely parametric.

The main v10 change should be philosophical: fitted multipliers should no
longer be allowed to stand in for missing physics. Instead, the model should
use fixed geometry and material properties, then identify or validate a small
set of heat-transfer coefficient families:

- gas-solid convective exchange in square monolith channels,
- axial thermal radiation through the structured porous/channel domain,
- solid/conductive heat paths through the monolith, adaptor, insulation, rear
  tube, and cooled flange.

## @contextScopeItemMention - Study Goal

The whole 1D study should be framed as:

```text
Develop a reduced-order, physics-constrained model of a structured monolithic
solar receiver with square channels, using transient thermocouple and gas
temperature data to obtain, test, or validate effective heat-transfer
coefficients for convective gas-solid exchange, axial radiative transfer, and
solid/conductive heat redistribution.
```

This is stronger than saying that the model is being "fit" to the experiments.
The target is not only prediction. The scientific deliverable is a defensible
coefficient closure for this specific structured monolithic receiver.

## What v3-v8 Have Already Taught Us

### v3

v3 established the useful conservative finite-volume structure:

```text
dynamic states: axial solid temperatures
gas: quasi-steady plug-flow march with exact cell NTU updates
solid axial conduction: conservative face fluxes
```

This should remain the core numerical architecture.

### v4

v4 showed that the downstream gas-solid exchange had to weaken, but the chosen
axial exchange shape was too free. The optimizer could suppress exchange until
it effectively became a compensating function rather than a heat-transfer
correlation.

### v5

v5 reduced parameters and introduced the staged logic:

```text
heating -> gas heat transfer and optical/input behavior
cooling -> thermophysical properties and hidden storage/losses
heating -> refit heat transfer
```

The staging was useful while the missing physics was unclear, but it also
encouraged different data subsets to tune different compensating knobs.

### v6

v6 added a rear thermal mass. This correctly targeted the slow deep response,
but the free rear state was too abstract. It could improve transients while
remaining difficult to interpret physically.

### v7

v7 used measured T2 as a boundary and replaced some free rear parameters with a
geometry-based resistance. This was a step toward physical interpretation. The
remaining important result was the signed heating mismatch:

```text
front/mid solid often too cold
gas/T3 can be too hot
```

That sign pattern cannot be fixed by a single energy multiplier. It points to
the coupling among gas-solid exchange, axial redistribution, radiative transfer,
outlet measurement location, and possibly active/bypassed flow.

### v8

v8 added the rear alumina tube, adaptor/flange heat path, and in v8b a predicted
cavity/T2 state. This made T2 a model prediction rather than an imposed thermal
boundary. The v8b quick-run interpretation is especially important:

```text
T2 level can be captured with a simple rear/cavity representation,
cooling is relatively improved,
heating solid temperatures remain underpredicted,
T3 remains sensitive to the gas comparison location and rear-tube interaction.
```

Therefore v10 should not add another empirical rear correction first. The next
missing term to test is axial thermal radiation inside the monolith domain and a
more defensible monolith-channel convective closure.

## 1. Remove k-Scaling and Exchange Factors

### Recommendation

Remove fitted conductivity scaling factors from the v10 parameter vector:

```text
remove: k_scale
remove: k_ins_scale
remove: any f_exchange or similar gas-solid exchange multiplier
```

Use only fixed or temperature-dependent material properties:

```text
solid/substrate conductivity: k_s(T)
alumina adaptor/tube conductivity: k_alumina(T)
insulation conductivity: k_ins(T), or fixed k_ins if no reliable T-dependence
gas conductivity: k_f(T)
gas Cp, viscosity, density: temperature-dependent functions
```

If the current `ks_f(T)` represents dense SiC but the modeled object is a
monolith with voids and square channels, do not multiply dense SiC by a fitted
`k_scale`. Instead define an explicit effective axial conductivity from the
monolith geometry:

```text
k_z,cond,eff(T) = lambda_s k_s(T) + phi k_f(T)
```

where:

```text
phi      = open porosity / channel void fraction
lambda_s = solid substrate volume fraction in the axial direction
```

If washcoat or secondary solids matter:

```text
k_z,cond,eff(T) = lambda_s k_s(T) + xi_w k_w(T) + phi k_f(T)
```

This follows the monolith-modeling logic discussed by Hayes and Cornejo: axial
effective conductivity in a monolith should be directional and structure-based,
not treated as a generic dense-solid multiplier.

### Reason

The fitted `k_scale` has become a mixed surrogate for at least three different
things:

- true solid/substrate conductivity uncertainty,
- monolith anisotropic effective conductivity,
- unresolved axial radiative redistribution.

Those should be separated. v10 should represent conductive and radiative axial
transport as separate terms, otherwise the fitted coefficient cannot be
interpreted.

### Practical Consequence

The axial face flux should become conceptually:

```text
Q_axial_face =
    (k_z,cond,eff(T_face) + k_z,rad,eff(T_face)) A_eff
    (T_left - T_right) / dx
```

not:

```text
Q_axial_face =
    k_scale k_s(T_face) A_solid (T_left - T_right) / dx
```

## 2. Replace the Current Axial Exchange Shape

The current v5-v8 convective form is:

```text
Nu = 10^A_Nu Re^B_Re Pr^(10^C_Pr)
h_i = Nu k_f / Dh * s(z)
s(z) = h_floor + (1 - h_floor) exp(-z / L_h)
```

For v10 this should be retired or kept only as a legacy comparison. The problem
is not the exponential shape itself; the problem is that `h_floor` and `L_h`
are effectively an empirical axial exchange correction. They may fit the data,
but they do not directly validate a square-channel monolith heat-transfer
coefficient.

### Literature-Aligned Alternative

Use a Graetz/entry-length monolith correlation:

```text
Gz = Re Pr Dh / z
Gz_inv = z / (Re Pr Dh)
Nu(z) = Nu_inf [1 + B Gz^n] exp(-C / Gz)
```

For square channels, `Nu_inf` should be tied to the thermal boundary condition:

```text
constant wall temperature, square duct: Nu_inf about 2.98
constant wall heat flux, square duct:  Nu_inf about 3.61
```

The receiver is not exactly either limit, so the two should be treated as
bracketing cases or as two candidate closures. The solar receiver is closer to
a distributed wall/source problem than to an externally imposed constant wall
temperature, but strong axial solid conduction and radiation may move it away
from the simple heat-flux limit.

Cornejo and Hayes emphasize two points that matter directly here:

- variable gas properties can make the entry-region Nu curve non-unique,
- monolith inlet geometry and contraction effects can change the entry region,
  especially for short laboratory monoliths.

Therefore the v10 convective closure should keep temperature-dependent gas
properties and should compute `Re`, `Pr`, and `Gz` locally from the film state.

### Proposed v10 Convective Parameter Set

Start with the least flexible version:

```text
fixed:
  Nu_inf from square-channel boundary-condition choice
  B, C, n from selected literature correlation or initial canonical values

fit/test:
  one bounded Nu correction only if necessary, named C_Nu_model
```

If that fails structurally, move to:

```text
fit:
  B_Gz
  C_Gz
  n_Gz
```

but only with tight bounds and with cross-validation. Do not also fit
`h_floor`, `L_h`, and arbitrary exchange factors in the same version.

## 3. Add Rosseland Axial Radiative Diffusion

### Why This Is the Right v10 Physics Test

Ray tracing of the source showed that incident irradiance does not penetrate
deeply. That is a statement about shortwave solar absorption. It does not rule
out axial thermal radiation after the hot front region re-emits in the infrared
inside the channel/porous structure.

The observed behavior:

```text
front heats rapidly,
internal/deep parts respond more slowly,
deep response is less flow-dependent than predicted,
solid and gas steady states cannot be matched simultaneously
```

is consistent with a missing solid/structure-mediated redistribution term. Axial
thermal radiation is a plausible candidate because it is strongly nonlinear in
temperature and not directly proportional to gas flow.

### Rosseland Form

For an optically thick participating structure, represent axial thermal
radiation as an effective diffusion term:

```text
q_rad,z = - k_rad,R(T) dT_s/dz
```

with:

```text
k_rad,R(T) = 16 sigma n_eff^2 T^3 / (3 beta_R)
```

where:

```text
sigma   = Stefan-Boltzmann constant
n_eff   = effective refractive index, usually near 1 for gas-filled channels
beta_R  = Rosseland or transport extinction coefficient [1/m]
```

If scattering is represented explicitly, use a transport extinction coefficient:

```text
beta_tr = beta_ext (1 - omega g)
```

where:

```text
beta_ext = extinction coefficient
omega    = scattering albedo
g        = asymmetry factor
```

For v10, a practical first implementation is:

```text
k_rad,R(T) = 16 sigma T^3 / (3 beta_tr)
```

and fit or scan `beta_tr`, not a generic radiation multiplier.

### Optical Thickness Check

Rosseland diffusion is not automatically valid. It is mainly appropriate when
the medium is optically thick:

```text
tau_R = beta_tr L_char
```

A working threshold is:

```text
tau_R > 10: Rosseland is defensible
tau_R = 1-10: use as a sensitivity / approximate closure
tau_R < 1: Rosseland is likely not valid; consider two-flux, P1, or view-factor treatment
```

This check should be reported for every fitted `beta_tr`.

### Parameter Bounds

Use the ray-tracing result to define an initial range. If the source penetration
length is `ell_abs`, then:

```text
beta_shortwave about 1 / ell_abs
```

but do not assume this equals the thermal infrared `beta_tr`. The shortwave and
thermal radiation fields see different optical behavior. A sensible v10
optimization parameter is:

```text
ell_rad = 1 / beta_tr
```

because it is easier to interpret as a radiative diffusion length. Bound it by
geometry:

```text
channel hydraulic diameter < ell_rad < receiver length
```

or, more conservatively:

```text
0.5 Dh < ell_rad < L
```

depending on how transparent the channels are expected to be in the IR.

### Finite-Volume Form

At each axial face:

```text
T_face = 0.5 (T_s,i + T_s,i+1)
k_face,total = k_z,cond,eff(T_face) + k_z,rad,eff(T_face)
Q_axial,i+1/2 = k_face,total A_eff (T_s,i - T_s,i+1) / dx
```

where:

```text
k_z,rad,eff(T_face) = chi_rad k_rad,R(T_face)
```

`chi_rad` should not be a free exchange factor in the first v10 attempt. Prefer
deriving it from open area / porosity / line-of-sight geometry. If it must be
fitted, it should be bounded between 0 and the open-area fraction and reported
as a geometry participation factor, not as a heat-transfer coefficient.

### Boundary Treatment

Do not double count front radiation. The existing front-face emission term
already represents radiative loss from the exposed front surface:

```text
Q_front,rad = eps_front sigma A_frt (T_s,1^4 - T_amb^4)
```

The Rosseland term should initially be used only for internal axial
redistribution between receiver cells. Rear radiative coupling to the adaptor or
tube can be added later only if heat-flow diagnostics show that the internal
term is insufficient.

### Literature Alignment

Mendes et al. support using effective conductivity ideas for combined
conduction-radiation in porous structures, but they also warn that Rosseland is
mainly valid for optically thick media and may overpredict effective thermal
conductivity when the optical thickness is marginal. This is exactly why v10
should report `tau_R` and compare:

```text
convective only
convective + Rosseland axial radiation
convective + bounded non-Rosseland sensitivity, if needed
```

Moro Filho and Malalasekera are even more directly relevant to v10 because they
analyze thermal radiation in porous media under local thermal non-equilibrium
(LTNE). Their model uses separate fluid and solid energy equations coupled by a
volumetric/interfacial heat-transfer term:

```text
fluid phase:
  advection = divergence of fluid effective conduction
              + h_v (T_s - T_f)

solid phase:
  divergence of solid effective conduction
  - h_v (T_s - T_f)
  = divergence of radiative heat flux
```

For the Rosseland approximation, they include the radiative term inside the
solid effective conductivity:

```text
K_eff,s =
    (1 - phi) k_s I
  + [16 sigma T_s^3 / (3 K_r)] I
```

where `K_r` is the local Rosseland mean attenuation coefficient. This supports
placing the first v10 Rosseland term in the solid/structure energy balance, not
in the gas enthalpy march.

Their results are also a warning. Rosseland differed substantially from DOM/FVM
radiative-transfer solutions:

```text
channel case: differences up to about +34% and -26%
tube case:    differences up to about 76%
```

The main reason is that the diffusion approximation loses accuracy near
boundaries and in high wall-heat-flux regions, which are exactly relevant for an
irradiated receiver. Therefore v10 should not present Rosseland as the final
radiation truth. It should be the first low-cost radiation closure, reported
with an optical-thickness diagnostic and, if it materially changes the results,
flagged for later comparison against a DOM/P1/two-flux or view-factor style
model.

Moro Filho and Malalasekera also show that the extinction coefficient can affect
the apparent Nusselt number non-monotonically. In their porous tube case, Nu
decreased with increasing extinction coefficient at low beta and increased at
higher beta. That means v10 should not interpret fitted `beta_tr` only as "more
radiation equals more effective heat transfer." The radiation field can
redistribute heat in ways that change both local solid temperatures and the
apparent gas-solid coefficient.

## 4. Refactor the Optimization Stages

The current three-stage sequence was useful diagnostically:

```text
1. heating heat-transfer + irradiance factors
2. cooling thermophysical / conduction scales
3. heating refit
```

For v10, this exact staging is probably no longer scientifically clean because:

- `k_scale` and `k_ins_scale` should be removed,
- axial radiation affects both heating and cooling,
- convection should be represented by a monolith-channel coefficient, not by an
  empirical exchange shape,
- the final fitted vector should be interpretable as coefficient evidence.

### Proposed v10 Workflow

#### Stage 0 - Physics Audit, No Optimization

Fix or document:

```text
geometry:
  receiver length, channel size, hydraulic diameter, porosity, frontal area
  solid area, exchange perimeter, sensor axial/radial positions
  rear tube/adaptor/flange/cavity dimensions

properties:
  k_s(T), Cp_s(T), rho_s
  k_alumina(T), Cp_alumina(T), rho_alumina
  k_ins(T) or fixed k_ins
  gas Cp(T), mu(T), k_f(T), rho(T)

boundary inputs:
  irradiance, flow, Tin, Tamb, T2 if used as observed validation
```

This stage should produce a table of fixed assumptions and not call an
optimizer.

#### Stage 1 - Baseline Forward Runs and Ablation

Run three forward model variants with the same fixed inputs:

```text
A. convection only + fixed conductive properties
B. convection + fixed conductive properties + rear/cavity domain
C. convection + fixed conductive properties + rear/cavity + Rosseland radiation
```

The purpose is to show what radiation changes before any new coefficient is
fitted.

#### Stage 2 - Coefficient Identification

Fit only coefficient-family parameters:

```text
convective:
  C_Nu_model or B_Gz, C_Gz, n_Gz

radiative:
  beta_tr or ell_rad = 1 / beta_tr

optional optical nuisance:
  one global eta_opt or one measured-irradiance correction
```

Avoid fitting separate `f_I_high`, `f_I_mid`, and `f_I_low` unless the
irradiance measurement uncertainty is independently justified. If retained,
they should be reported as input-calibration nuisance parameters, not as
receiver heat-transfer physics.

#### Stage 3 - Validation and Residual Structure

Use holdout experiments or leave-one-case-out validation. Report more than
RMSE:

```text
signed steady error
T8/T9/T10 axial profile error
T3 signed error
T2/cavity error
cooling overshoot
t90 timing error
fitted beta_tr and tau_R
mean h(z), Nu(z), and k_rad(T) ranges
energy partition among gas, side loss, rear/flange, and axial radiation
```

The key test is whether v10 changes the wrong sign pattern:

```text
solid too cold while gas too hot
```

If that sign pattern remains, the problem is probably not axial radiation alone.

## 5. Suggested v10 Parameter Vector

A minimal v10 vector could be:

```text
p[1]  C_Nu_model     bounded correction to selected monolith Nu closure
p[2]  ell_rad        Rosseland radiative mean path, ell_rad = 1 / beta_tr
p[3]  eta_opt        optional global irradiance/input correction
```

If the monolith Nu correlation is fitted rather than corrected:

```text
p[1]  B_Gz
p[2]  C_Gz
p[3]  n_Gz
p[4]  ell_rad
p[5]  eta_opt
```

Optional, only if physically unresolved:

```text
p[6]  h_contact_rear
p[7]  h_flange_contact
```

These should be contact/boundary coefficients, not conductivity multipliers.

Do not include in the main v10 vector:

```text
gamma_C
k_scale
k_ins_scale
h_floor
L_h
tau_T3, unless a physical sensor node is not yet implemented
U_side
U_rear
C_rear_scale
K_rear
U_rear_mass
generic f_exchange factors
```

If `gamma_C` is still needed for numerical sensitivity, it should be renamed
and reported as a thermal-inventory uncertainty, not as a fitted material
property. For the main coefficient-validation study, it is better fixed from
mass, geometry, and Cp.

## 6. T3 Treatment

The old first-order `tau_T3` lag is not ideal for v10. The v8b change to sample
gas at `z = 140 mm` is more physical than fitting a lag, but T3 may still not be
a pure gas mixing-cup temperature.

Recommended v10 treatment:

```text
primary:
  compare T3 to gas temperature at the measured thermocouple location

sensitivity:
  add a small physical sensor/outlet node only if needed
```

A physical node would have:

```text
C_T3 dT_T3/dt =
    h_g,b A_b (T_g - T_T3)
  + h_support A_support (T_support - T_T3)
  + radiation/conduction terms if geometrically justified
```

This is better than a free `tau_T3` because it can be checked against bead size,
support geometry, and local flow.

## 7. Legacy Cleanup List for v10

When code editing starts, remove or isolate:

- `k_scale` and all references to scaled solid axial conductivity.
- `k_ins_scale` and any scaled insulation/adaptor-to-cavity conductance.
- `h_floor` and `L_h` if the monolith/Graetz Nu closure replaces the axial
  exchange shape.
- Unused `tau_T3` compatibility paths, especially if v10 uses a physical T3
  location or sensor node.
- v6 rear-mass parameters and comments:

```text
C_rear_scale
K_rear
U_rear_mass
```

- v5-style empirical side/rear losses:

```text
U_side
U_rear
```

- Any comments saying that the purpose is simply to improve fit quality.
  Replace with coefficient-identification language.
- Version comments that describe obsolete geometry, especially v8a 200 mm rear
  tube assumptions if v8b/v10 uses the corrected 150 mm tube and 46/104 mm
  cavity/flange split.
- Global mutation patterns such as updating `pnew_v*` inside calibration, if a
  clean reproducible v10 runner is being prepared.

Keep:

- conservative finite-volume axial balance,
- exact NTU gas marching,
- temperature-dependent gas properties,
- paired thermocouple comparisons:

```text
T9_pair  = (T9 + T12) / 2
T10_pair = (T10 + T11) / 2
```

- rear tube/flange/cavity heat-flow diagnostics,
- signed residual and timing diagnostics.

## 8. Interpretation Rules for v10 Results

v10 should be judged by whether the fitted coefficients are physically credible,
not only by lower objective value.

### A credible v10 outcome would show:

```text
Nu(z) remains in the expected square-channel monolith range.
beta_tr implies optical thickness consistent with the Rosseland assumption.
k_rad,R(T) is important mainly at high T and does not mimic all conduction.
T2/cavity remains credible without fitted insulation scaling.
solid and gas signed errors improve together.
cooling tails improve without artificial rear heat storage.
```

### A concerning v10 outcome would show:

```text
ell_rad pushed to bounds.
Rosseland tau_R < 1 while the model still relies on diffusion.
Nu correction pushed far above square-channel expectations.
eta_opt pushed above physically credible absorbed-power limits.
solid temperatures improve only by worsening T3.
T2 is recovered only by reintroducing hidden scaling.
```

## 9. Recommended v10 Comparison Matrix

Use a small comparison matrix rather than one all-in fit:

```text
M0: v8b baseline, fixed previous best parameters
M1: fixed properties, no k scales, monolith Nu, no radiation
M2: M1 + Rosseland axial radiation
M3: M2 + fitted beta_tr / ell_rad
M4: M3 + fitted monolith Nu parameters, if necessary
```

For each model report:

```text
heating RMSE and signed steady error by sensor
cooling RMSE and signed steady error by sensor
T3 error at z = 140 mm
T2/cavity error
energy partition
Nu(z) range
k_cond,z range
k_rad,z range
tau_R range
```

This matrix will make it clear whether radiation is doing real explanatory
work or only replacing `k_scale` under a new name.

## 10. Bottom Line

The v10 direction is scientifically justified. The model should now stop
calibrating broad multipliers and move toward coefficient validation:

```text
fixed material/geometry properties
+ square-channel monolith convective Nu closure
+ explicit Rosseland axial radiation term
+ physically represented rear/cavity/flange domain
+ validation by signed transient residuals and heat-flow partition
```

The most important discipline is to avoid adding radiation as another arbitrary
effective conductivity scale. Add it as a separate, interpretable transport
mechanism with its own optical-thickness check.

## 11. v10a Implementation Result

Implemented the first v10 model as a comparison-first ECM baseline:

```text
files:
  1D_v10.jl
  run_1D_v10.jl
  test/smoke_1D_v10.jl

outputs:
  summaries/1D_v10/analysis_results_all_variants_1D_v10.csv
  summaries/1D_v10/plots/*_1D_v10.png
```

The model removes the broad fitted thermal scalars and keeps only coefficient
family parameters:

```text
p[1] C_Nu_model
p[2] ell_rad_m
p[3] eta_opt       fixed to 1.0 for the next stage
p[4] front_dep     fixed to 1.0 for the next stage
p[5] beta_opt
```

The current smoke test passed:

```text
49/49 tests passed
```

The comparison matrix was:

```text
M1_no_rad:        fixed material/geometry properties, monolith Nu, no Rosseland
M2_rad_beta1000:  M1 + Rosseland axial radiation, beta_tr = 1000 1/m
M3_rad_beta300:   Rosseland sensitivity, beta_tr = 300 1/m
M3_rad_beta2700:  Rosseland sensitivity, beta_tr = 2700 1/m
```

Mean heating signed steady errors were approximately:

```text
M1_no_rad:
  T8       -219 K
  T9_pair  -241 K
  T10_pair -156 K
  T3       -118 K
  T2        -11 K

Rosseland variants:
  nearly unchanged over beta_tr = 300-2700 1/m
```

Cooling signed steady errors were much smaller:

```text
T8        +3 K
T9_pair   +2 K
T10_pair  +7 K
T3       +15 K
T2        +3 K
```

### v10a Interpretation

Rosseland axial diffusion does not move the heating residuals enough to explain
the v8b/v9 mismatch. This is a useful negative result: within the tested optical
range, radiation is not acting as the dominant missing transport pathway.

The fixed-property v10a baseline is much too cold during heating, while the
cavity/T2 level remains fairly close. That points away from adding another
hidden cavity or insulation scalar and toward the absorbed-power/source side of
the model:

```text
absorbed input level
shortwave deposition profile
front-boundary loss balance
monolith Nu/source coupling
```

### Recommended v10b

Use the user's fixed optical/source convention first:

```text
eta_opt   = 1.0
front_dep = 1.0
```

Under the current `solar_weights_v5` convention, `front_dep = 1.0` assigns all
shortwave source to the first receiver cell, so `beta_opt` is inactive until
`front_dep` is released again. Per-irradiance absorbed-power calibration should
therefore be deferred to a second stage.

The next v10b fit-ready physics set should be:

```text
Nu A/B/C coefficients
low-dimensional axial Rosseland parameters
```

The next success criterion should be whether Nu and spatial Rosseland can
improve the deeper receiver and T3 residuals without reintroducing `k_scale`,
`k_ins_scale`, `gamma_C`, `h_floor`, `tau_T3`, or `f_exchange`.

After fixing `eta_opt = 1.0` and `front_dep = 1.0`, the regenerated
`M1_no_rad` mean heating signed steady errors are:

```text
T8        -192 K
T9_pair   -248 K
T10_pair  -164 K
T3        -125 K
T2         -11 K
```

Compared with the previous `front_dep = 0.50` matrix, the front thermocouple is
less cold but the deeper receiver and T3 are slightly colder. This confirms that
front-only deposition does not solve the deeper receiver deficit by itself.
