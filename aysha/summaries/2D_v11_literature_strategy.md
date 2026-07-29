# Literature-grounded strategy for 2D_v11

Date: 2026-07-29

## Purpose

This note translates the ceramic volumetric-receiver and monolith literature
in `analysis/literature/` into a falsifiable next step for the 2D receiver
model. It is not a calibration report and it does not accept any new
coefficient.

The primary objectives are to:

1. replace the present apparent entrance heat-transfer law by a physically
   bounded local channel law;
2. use the corrected MFC flow as the conserved total flow;
3. use DP1 to test the resulting temperature-dependent velocity and pressure
   field; and
4. investigate the radial profile without a bypass branch, transverse flow
   through the channel walls, or a side-weighted irradiance pattern.

## Most relevant literature evidence

### Capuano, Fend et al. (2016): ceramic receiver channel model

`analysis/literature/Capuano, Fend, et al (2016) - Numerical models of advanced ceramic absorbers for volumetric solar receivers.pdf`

The reviewed Becker ceramic-receiver model separates the exposed inlet web
from the internal square channel. For developing internal flow it reports the
Kays form

$$
\mathrm{Nu}_{K}(z)
=
3.66+
\frac{0.104\,\mathrm{Gz}(z)}
{1+0.016\,\mathrm{Gz}(z)^{0.8}},
\qquad
\mathrm{Gz}(z)=
\mathrm{Re}(z)\,\mathrm{Pr}(z)\frac{D_h}{z}.
$$

It also reports the variable-temperature correction

$$
\mathrm{Nu}_{K,T}
=
\mathrm{Nu}_{K}
\left(\frac{T_g}{T_w}\right)^{0.45},
\qquad
0.5 < \frac{T_g}{T_w} < 1.5.
$$

The front-web correlation is a Churchill-Bernstein cylinder-crossflow law
based on the web thickness, not on the channel hydraulic diameter. Its stated
range begins at approximately $\mathrm{Re}=100$. In the present receiver,
$b/D_h=0.4/1.5$ and the channel Reynolds range is about 33--135, so
$\mathrm{Re}_b$ is only about 9--36. The reported front correlation is
therefore outside its stated Reynolds range and should not be installed as
the v11 base law.

The paper also separates pressure drop into inlet, outlet and wall-friction
terms and evaluates density at relevant local temperatures.

### Fend (2013): where gas heating occurs

`analysis/literature/Fend (2013) - Numerical investigation offlow and heat transfer in a volumetric solarreceiver.pdf`

The detailed SiC single-channel calculation shows a large wall-to-gas heat
flux at the channel entrance, followed by a rapid decay. The air temperature
rises strongly over approximately the first 20--30 mm; the interphase heat
flow becomes small farther downstream. The continuum model obtains its
volumetric heat-transfer coefficient from this single-channel heat-flow
distribution.

This supports a distributed internal entrance layer. It does not support an
arbitrary instantaneous gas-temperature jump at the exposed front face.

### Cornejo and Hayes (2020): entry geometry and variable properties

`analysis/literature/Cornejo, Hayes, 2020 Nu Number Correlations Monolith Reactors.pdf`

The paper states that a useful monolith correlation must describe the entry
region and approach an analytic fully developed value. Its general form is

$$
\mathrm{Nu}
=
\mathrm{Nu}_{\infty}
\left[
1+B\,\mathrm{Gz}^{n}\exp\left(-C/\mathrm{Gz}\right)
\right].
$$

More importantly for this experiment:

- using inlet rather than local properties can move the apparent Graetz curve
  by 50% or more under strong heating;
- the inlet contraction changes the velocity profile and produces a lower
  entry Nusselt number than idealized straight/reservoir inlets;
- the effect of temperature-dependent viscosity is particularly important;
- at high wall-to-gas temperature ratios, the local Nusselt number can pass
  below the nominal asymptote; and
- a single universal Graetz curve may not represent all heating rates.

Their simulated channel Reynolds range of 50--600 overlaps much of the
present range. Their multi-parameter high-heating correlation should not be
copied directly into v11 because it was derived for a different channel and
thermal boundary condition and the authors themselves recommend parameter
reduction.

### Hayes and Cornejo (2021): local versus average correlations

`analysis/literature/Hayes, 2021 Multi-scale modelling of monolith reactors A 30-year perspective from 1990 to 2020.pdf`

The review distinguishes local Nusselt numbers, which a distributed channel
model needs, from length-averaged Graetz correlations. It warns against
correlations that do not approach a fully developed asymptote. It also notes
that realistic inlet profiles, developing friction and inlet/outlet losses
can matter, particularly for short laboratory monoliths.

The review gives 3.656 and 4.364 as the circular-duct limits for constant wall
temperature and constant wall heat flux and explicitly notes that
non-circular channels have different limits. For an ideal square duct, the
commonly used limits are approximately 2.98 and 3.61 for the corresponding
boundary-condition families. The receiver boundary is neither exactly
uniform-temperature nor uniform-flux, so v11 should test these limits as
fixed model alternatives rather than fit the asymptote.

### Fend, Hoffschmidt et al. (2004): parallel-channel nonuniformity

`analysis/literature/Fend, Hoffschmidt (2004) - Porous materials as open volumetric solar receivers - Experimental determination of thermophysical and heat transfer properties.pdf`

In a cooling comparison, a cordierite foam rapidly compensated a small front
temperature nonuniformity, whereas a parallel-channel catalyst carrier
retained a nonuniform temperature distribution. The authors connect the
difference to architecture and pressure-drop characteristics. They also
report flow instabilities and hot spots as an important receiver limitation.

This supports retaining isolated axial channels and investigating their
parallel hydraulic coupling. It does not justify radial mass exchange through
the walls of a honeycomb.

### Avila-Marin et al. (2019) and Kribus et al. (2014): thermal-hydraulic
feedback

`analysis/literature/Avila-Marin, 2019 review modelling strategies porous structures solar receivers.pdf`

`analysis/literature/Kribus (2014) - The promise and challenge of solar volumetric absorbers.pdf`

The Avila-Marin review recommends LTNE for volumetric receivers and reports
that nonuniform heating can cause hotter regions to have greater flow
resistance, redirecting flow to colder regions and leaving the hot regions
less cooled. It also reports that linear pressure-drop structures are
particularly susceptible to instability. The review says forced convective
loss from the porous inlet is usually small compared with thermal loss.

Kribus likewise uses temperature-dependent density, viscosity, heat capacity
and conductivity and identifies laterally nonuniform velocity and
thermally-induced flow instability as unresolved multidimensional effects.

These observations motivate a conservative equal-pressure parallel-channel
calculation. They do not motivate fitting an unmeasured bypass.

### Evidence used only with caution

The foam and wire-mesh papers are useful for LTNE, internal radiation and
effective-conductivity modelling, but their Darcy-Forchheimer flow and
open-pore transverse mixing should not be transferred directly to a
straight-channel monolith. The 2022 periodic-structure experiment is also a
useful warning that nonuniform flux and high-temperature outlet-flow
measurements complicate interpretation, but its morphology is not the present
square-channel geometry.

## Audit of the current v9/v10 closure

The current internal law is

$$
\mathrm{Nu}_{v10}
=
A_{\mathrm{Nu}}\,
\mathrm{Re}^{1.44}\,
\mathrm{Pr}^{1/3}
\left(\frac{D_h}{z}\right)^{1/3}.
$$

Since

$$
\mathrm{Gz}
=
\mathrm{Re}\,\mathrm{Pr}\frac{D_h}{z},
$$

the same law can be written

$$
\mathrm{Nu}_{v10}
=
A_{\mathrm{Nu}}\,
\mathrm{Re}^{1.44-1/3}\,
\mathrm{Gz}^{1/3}.
$$

It therefore contains an additional $\mathrm{Re}^{1.107}$ dependence beyond
Graetz scaling. With `minimum_nusselt=0`, it also tends toward zero with
depth instead of a square-duct fully developed value.

Using the fitted $A_{\mathrm{Nu}}=0.00123739$, $\mathrm{Pr}=0.70$,
$D_h=1.5$ mm, the verified thermocouple depths 5, 58 and 107 mm, and the
experimental flow range 4.53--18.34 standard L/min gives:

| Flow | Re | Depth | Gz | Current Nu | Kays-form Nu |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 4.53 L/min | 33.3 | 5 mm | 6.98 | 0.114 | 4.335 |
| 4.53 L/min | 33.3 | 58 mm | 0.602 | 0.050 | 3.722 |
| 4.53 L/min | 33.3 | 107 mm | 0.326 | 0.041 | 3.694 |
| 18.34 L/min | 134.6 | 5 mm | 28.27 | 0.856 | 6.047 |
| 18.34 L/min | 134.6 | 58 mm | 2.44 | 0.378 | 3.905 |
| 18.34 L/min | 134.6 | 107 mm | 1.32 | 0.308 | 3.795 |

The comparison is not proof that the Kays expression will fit this receiver.
It does show that the fitted current law is an apparent whole-model
compensator, not an experimentally validated channel Nusselt correlation.
Reintroducing a finite asymptote in v11 is therefore not a reversal of the v8
audit: v11 changes the intended interpretation from an unconstrained
Macro-ECM fit to a channel-grounded LTNE closure whose failure will be
explicitly tested.

The current radial mass split is prescribed by

$$
\psi(r)=1-c_r(r/R)^2,\qquad c_r=0.10.
$$

This cannot respond to the calculated ring temperatures. Moreover, each ring
currently predicts its own pressure drop and the receiver reports an average;
physical parallel channels connected to common upstream and downstream
volumes must have one pressure drop.

## Recommended v11 model sequence

### V11-A: zero-refit local Graetz replacement

Keep the corrected MFC mass flow, geometry, optics, solid properties,
external losses, DP1 zero offset and cold hydraulic scale fixed. Disable the
v10 exposed-front gas jump (`front_coefficient=0`).

Replace the internal apparent law by a local channel model:

$$
\begin{aligned}
\mathrm{Re}_{i,j}
&=
\frac{\dot m_iD_h}{A_{f,i}\mu(T_{g,i,j})},\\
\mathrm{Pr}_{i,j}
&=
\frac{c_p(T_{g,i,j})\mu(T_{g,i,j})}{k_g(T_{g,i,j})},\\
\mathrm{Gz}_{i,j}
&=
\mathrm{Re}_{i,j}\mathrm{Pr}_{i,j}\frac{D_h}{z_j},\\
\mathrm{Nu}_{i,j}
&=
\mathrm{Nu}_{fd}
+
\frac{0.104\,\mathrm{Gz}_{i,j}}
{1+0.016\,\mathrm{Gz}_{i,j}^{0.8}},\\
h_{i,j}
&=
\frac{\mathrm{Nu}_{i,j}k_g}{D_h}.
\end{aligned}
$$

Use cell-center or cell-averaged $z$, not a fitted virtual origin. Test the
following fixed alternatives:

1. `Nu_fd=3.61`, square-duct uniform-flux limit;
2. `Nu_fd=2.98`, square-duct uniform-temperature limit;
3. the Capuano implementation with `Nu_fd=3.66` and its stated temperature
   correction.

No exponent or front coefficient is fitted in this stage. All gas properties
must be local and temperature-dependent, with the exact property-temperature
convention recorded.

### V11-B: thermal-hydraulic parallel-ring coupling

If V11-A improves or materially changes the axial profiles, replace the
prescribed radial flow shape by a conservative equal-pressure solution.
Every ring remains an isolated set of straight channels, so its mass flow is
constant along $z$ and there is no radial mass transfer inside the monolith.

At every thermal evaluation, solve for one common pressure drop and ring
flows:

$$
\Delta p_i(\dot m_i,T_{g,i}(z))
=
\Delta p_{\mathrm{common}},
\qquad
\sum_i\dot m_i=\dot m_{\mathrm{MFC}}.
$$

The pressure model initially retains the t0-calibrated square-channel
resistance scale. The local hot contribution uses the calculated
$\mu(T)$, $\rho(T)$ and

$$
u_i(z)=\frac{\dot m_i}{\rho_i(z)A_{f,i}}.
$$

Thus the total flow remains fixed by the corrected MFC, while hot rings can
receive less flow because their hydraulic resistance is larger. This is a
quantifiable internal redistribution, not a bypass.

The predicted DP1 signal is the common pressure drop plus the independently
established sensor offset/path representation. It must not be a radial
average of mutually different ring pressure drops.

#### Observed rear support-groove restriction

The outer receiver channels are open. The supporting groove at the receiver
outlet leaves an approximately 13 mm diameter region completely free, while
the discharge paths outside that diameter may face partial obstruction.
Therefore:

- do not remove outer channels from the flow geometry;
- do not reduce their channel-wall heat-transfer perimeter;
- do not introduce stagnant/closed-channel regions; and
- add a ring-dependent outlet minor loss only where the groove overlaps the
  outer flow paths.

Represent the groove contribution by

$$
\Delta p_{\mathrm{groove},i}
=
I_{\mathrm{groove},i}\,
K_{\mathrm{groove}}\,
\frac{\rho_{i,\mathrm{out}}u_{\mathrm{groove},i}^{2}}{2},
$$

where $I_{\mathrm{groove},i}$ is the fraction of ring $i$ outside the measured
free radius

$$
r_{\mathrm{free}}=6.5\ \mathrm{mm}.
$$

If the groove leaves a local flow-clearance fraction $\chi_g$,

$$
u_{\mathrm{groove},i}
=
\frac{\dot m_i}
{\rho_{i,\mathrm{out}}\,\chi_g A_{f,i}}.
$$

The outer-ring pressure equation becomes

$$
\Delta p_{\mathrm{common}}
=
\Delta p_{\mathrm{friction},i}
+
\Delta p_{\mathrm{groove},i}
+
\Delta p_{\mathrm{other},i}.
$$

The corrected receiver face area is 361 mm2 and its area-equivalent radius is
10.7196 mm. The measured free area and fractions are therefore

$$
\begin{aligned}
A_{\mathrm{free}}
&=\pi(6.5\ \mathrm{mm})^2=132.73\ \mathrm{mm^2},\\
f_{\mathrm{free}}
&=\frac{A_{\mathrm{free}}}{361\ \mathrm{mm^2}}=0.3677,\\
f_{\mathrm{groove\ exposed}}
&=1-f_{\mathrm{free}}=0.6323.
\end{aligned}
$$

Thus only about 36.8% of the face is guaranteed to have an unobstructed
downstream discharge, while about 63.2% may experience the support-groove
loss. This supersedes the earlier nominal mapping to only the outermost
channel row. It does not imply that 63.2% of the channels are closed; it
defines the region over which a partial quadratic outlet loss may act.

In the axisymmetric grid, calculate $I_{\mathrm{groove},i}$ from exact annular
area overlap with $r>6.5$ mm. Use 12 and 14 mm free diameters as conservative
geometry sensitivities around the gauged 13 mm value unless a more precise
measurement becomes available.

Because the groove term is quadratic in velocity, its relative effect and
the resulting core/edge flow contrast should decrease as total flow
decreases. Temperature-dependent $\mu/\rho$ feedback can oppose or reinforce
that trend, so the net convergence of the radial temperature profiles remains
a prediction rather than an imposed condition.

Use the cold t0 DP1-versus-flow data to identify or bound the single combined
groove strength $K_{\mathrm{groove}}/\chi_g^2$. Do not fit
$K_{\mathrm{groove}}$ and $\chi_g$ independently unless the clearance is
measured. When this physical groove model is enabled, disable the v9 empirical
hot-excess minor-loss coefficient so the same pressure loss is not counted
twice. The cold-identified groove law must then predict hot DP1 using local
outlet density and velocity without a new heating-only hydraulic fit.

### V11-C: inlet and outlet pressure/entry sensitivity

Only after V11-A/B should inlet geometry be added. Use the literature
separation

$$
\Delta p
=
\Delta p_{\mathrm{in}}
+
\Delta p_{\mathrm{friction}}
+
\Delta p_{\mathrm{out}},
$$

and compare:

1. the existing t0-calibrated lumped path;
2. a contraction-aware inlet/outlet model with geometry-fixed terms; and
3. if necessary, a small single-channel CFD surrogate spanning the measured
   Re and wall/gas temperature-ratio ranges.

Do not fit inlet loss, outlet loss, a virtual entry length and a new Nusselt
exponent simultaneously. DP1 at t0 first checks the cold pressure model;
heating DP1 then tests the predicted temperature effect.

### V11-D: at most one heat-transfer scale

If the fixed literature alternatives reproduce the profile shapes but have a
consistent magnitude error, allow one multiplier

$$
\mathrm{Nu}_{\mathrm{eff}}=S_h\,\mathrm{Nu}_{\mathrm{literature}}.
$$

Interpret $S_h$ as an effective-area/model-discrepancy factor, not as a new
fundamental Nusselt exponent. A result far from order unity is a rejection
signal requiring an audit of interfacial area, optical deposition, channel
geometry and sensor representation. It must not be hidden by simultaneous
optical and conductivity refitting.

### V11-E: limited recalibration only after model-form acceptance

Only after the transport model passes the no-refit tests should a small,
predeclared calibration be run. Keep literature exponents and asymptotes
fixed. Separate heating training cases, held-out heating cases and cooling
validation exactly as in v9.

## Required diagnostics

Every run must store:

- local `Re`, `Pr`, `Gz`, `Nu`, `h` and `h*A_v`;
- ring mass flows and their sum;
- groove overlap, groove velocity and groove pressure loss by ring;
- the common receiver pressure drop and DP1 prediction;
- gas temperature at the channel inlet, 5, 20, 30, 58, 107 and 137 mm;
- wall-to-gas heat transferred per axial cell;
- the existing solid sensors and T3;
- `T12-T8`, `T12-T9` and `T11-T10`; and
- energy and mass residuals.

The 5, 58 and 107 mm axial positions are fixed from the verified
thermocouple placement. The side thermocouples remain mapped to the middle of
the side wall at those depths; their shallow mounting dips may affect contact
and response time but do not justify changing the irradiance map.

## Acceptance and rejection gates

V11-A/B advances only if it meets all of the following:

1. exact total mass conservation and equal ring pressure drop;
2. energy-rate closure at the existing numerical tolerance;
3. finite, mesh-convergent entrance heat transfer and
   $\mathrm{Nu}\rightarrow\mathrm{Nu}_{fd}$ downstream;
4. improved axial profile relative to the best v10 sensitivity
   (`T12-T8` RMSE below 33.2 K) without relying on an exposed-front
   exponent;
5. correct signs and improved magnitudes of the three irradiance-group flow
   slopes of `T12-T8`;
6. no degradation of held-out heating DP1 beyond the v9 reference RMSE of
   approximately 0.193 mbar without an explicitly justified pressure-model
   change;
7. improvement in the observed positive radial offsets
   (`T12-T9` and `T11-T10`) when equal-pressure redistribution is enabled;
8. decreasing predicted core/edge flow contrast with decreasing total flow
   when the quadratic groove term dominates, consistent with the observed
   profile convergence;
9. no gross degradation of T3 or the cooling validation; and
10. any fitted $S_h$ remains physically near unity and transfers to held-out
   cases.

Failure of fixed local-Nu alternatives is still useful. It would show that a
channel-scale surrogate, optical/source revision or observation-model audit
is needed. It is not permission to reinstate the fitted
$\mathrm{Re}^{1.44}$ law or the $m_f=1.5$ exposed-front stress test as a
physical coefficient.

## Decision

The next implementation should be `2D_v11` with switchable V11-A and V11-B
physics. The first production run is a no-refit model-form matrix over all 15
heating cases. The recommended nominal branch is:

```text
conserved corrected-MFC total flow
local temperature-dependent gas properties
square-channel Kays/Graetz entry law
Nu_fd = 3.61
no exposed-front gas-temperature jump
no bypass and no radial channel-to-channel mass exchange
equal-pressure parallel-ring flow allocation
quadratic rear support-groove loss over the measured outer belt
existing optical, solid and external-loss parameters frozen
```

This is the shortest literature-supported route to determine whether the
missing flow sensitivity belongs to distributed channel entry heat transfer
and thermal-hydraulic maldistribution rather than to a fitted front-web
coefficient.
