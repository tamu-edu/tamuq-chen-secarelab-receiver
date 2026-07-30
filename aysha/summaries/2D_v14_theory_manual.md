# 2D v14 theory manual

## Purpose

V14 replaces the axisymmetric receiver field and the unsuccessful v13
two-population surrogate with a symmetry-reduced network of the actual
10 by 10 square channels.  Its purpose is to test channel-scale
thermal-hydraulic maldistribution while conserving the corrected MFC total
flow and retaining the validated v12 hardware topology.

## 1. Square-channel orbit reduction

The receiver contains 100 identical 1.5 mm square channels on a 1.9 mm pitch
inside a 19 mm square.  Reflection and diagonal symmetry reduce these to 15
unique channel orbits.  If $n_g$ is an orbit multiplicity,

$$
\sum_{g=1}^{15} n_g=100.
$$

Every channel in an orbit has the same solid temperature, gas profile and
per-channel flow.  The reduction preserves all edge/corner positions and
square-grid neighbor connections.

## 2. Receiver inventory

The SiC area is

$$
A_s=w_{\rm rec}^2-N_{\rm ch}b^2.
$$

Each physical channel cell owns $A_s/100$ of solid.  Orbit $g$ therefore has
axial-cell capacity

$$
C_{g,j}=n_g\rho_{\rm SiC}c_p(T_{g,j})
\frac{A_s}{100}\Delta z_j.
$$

The orbit capacities sum exactly to the measured receiver capacity; no
second mass or artificial heat capacity is added.

## 3. Parallel-channel hydraulics

Each orbit has one isolated gas stream.  At a thermal residual evaluation,
the per-channel mass flows satisfy

$$
\Delta p_g(\dot m_g,T_{g,\rm gas}(z))
=\Delta p_{\rm common},
$$

$$
\sum_g n_g\dot m_g=\dot m_{\rm MFC}.
$$

Density and viscosity are evaluated locally.  The laminar square-channel
friction term is inherited from v11.  Because pressure drop changes with
local temperature, hotter orbits may receive less mass flow.

## 4. Rear support groove

The measured completely free diameter is 13 mm.  V14 calculates the fraction
of each square channel area outside the 6.5 mm radius.  The outlet loss is

$$
\Delta p_{{\rm groove},g}
=I_gK_{\rm groove}\frac{\rho_{g,o}u_{g,o}^2}{2}.
$$

The coefficient remains the cold-DP-calibrated v11/v12 value.  Outer channels
remain open and retain their full heat-transfer perimeter.

## 5. Gas heat transfer

Gas is marched independently through every orbit using the local
temperature-dependent properties and square-channel Graetz correlation:

$$
{\rm Gz}_{g,j}={\rm Re}_{g,j}{\rm Pr}_{g,j}
\frac{D_h}{z_j},
$$

$$
{\rm Nu}_{g,j}=S_h\left[
{\rm Nu}_{fd}
+\frac{0.104{\rm Gz}_{g,j}}
{1+0.016{\rm Gz}_{g,j}^{0.8}}\right].
$$

The heat removed from the orbit solid is exactly the enthalpy gained by its
gas stream.  All streams mix before entering the alumina exit tube.

## 6. Solid conduction

Each orbit conducts axially through its share of SiC area.  Square-grid
neighbors exchange lateral heat through the effective transverse
conductivity inherited from the receiver continuum:

$$
Q_{g-h,j}=N_{gh}k_{\perp,\rm eff}(T)\Delta z_j
(T_{g,j}-T_{h,j}),
$$

where $N_{gh}$ is the exact number of physical neighbor interfaces between
the two symmetry orbits.  Internal exchanges enter with opposite signs and
cancel in the energy ledger.

Exterior channel cells exchange with felt in proportion to the number of
outer square faces they own.  Over the final 29 mm those faces also exchange
with the loose alumina adaptor.

## 7. Optical source and losses

The existing centered Gaussian face field is evaluated at each orbit
centroid and normalized over the square face.  Axial deposition retains the
Beer-Lambert plus front-deposition form.  No side-biased irradiance is
introduced.

Front radiation is evaluated per orbit.  Felt, casing, adaptor, exit tube,
rear housing and water-cooled flange retain the v12 equations and capacities.
There is no tube/felt contact or direct receiver/tube bypass.

## 8. Sensor mapping

The verified axial locations are 5, 58 and 107 mm.

| Sensor | V14 state |
| --- | --- |
| T8, T12, T11 | midpoint exterior-side channel orbit |
| T9, T10 | central channel orbit wall/gas blend |
| T3 | mixed rear gas 3 mm downstream |
| T2 | felt field at 58 mm |

The established shared observation time constants remain fixed for the
initial model-form test.

## 9. Pre-fit gates

V14 must demonstrate:

- 15 orbits with multiplicities summing to 100;
- exact SiC inventory and total-flow closure;
- a common pressure drop within numerical tolerance;
- lower flow through an otherwise identical deliberately hotter orbit;
- correct spatial application of the rear groove;
- uniform equilibrium; and
- instantaneous whole-assembly energy closure.

The initial run freezes v12 hardware, observation, optical, Graetz, power and
hydraulic quantities.  Subsequent training may fit only globally shared,
predeclared quantities.  Held-out heating, per-experiment properties,
side-weighted irradiance and an unmeasured bypass are excluded.

## 10. Verification outcome

The implementation passes 35 smoke checks:

- 15 symmetry orbits and 100 total channels;
- 40 exterior square faces;
- exact receiver mass and flow area;
- uniform equilibrium;
- exact total-flow conservation;
- common-pressure convergence;
- lower flow in a deliberately hotter channel orbit; and
- whole-assembly instantaneous energy closure.

## 11. Cold square-network hydraulics

The annular v11/v12 groove coefficient cannot be transferred directly to
square channel overlaps.  A `t0`-only fit to all 15 heating experiments,
with the established DP1 offset frozen, selected

$$
S_R=0.80,\qquad K_{\rm groove}=335,
$$

with DP1 RMSE 0.02347 mbar.  The best no-groove model gives 0.03359 mbar.
Every thermal calibration and validation result below uses the reidentified
square-network values.

## 12. Staged fit

The topology screen varied transverse conduction, edge/felt contact and the
width of a centered Gaussian beam.  The independent heating stage then
screened optical extinction and the single Graetz multiplier and finally the
three shared power scales.

The definitive selected values are:

| Quantity | Value |
| --- | ---: |
| transverse-conductivity multiplier | 1.00 |
| edge/felt contact multiplier | 3.00 |
| centered beam sigma | 100 mm |
| optical extinction | 100 1/m |
| Graetz multiplier | 1.80 |
| power scales | 1.25, 1.2075, 0.8855 |

Several quantities are bound-active.  A 100 mm sigma is effectively uniform
over the 19 mm square but remains centered, not side weighted.

Weakly connected channel networks can correct the mid radial sign in 7/9
training cases.  Strong isolation corrects four deep signs.  Their heating
and cooling RMSEs rise to approximately 137--142 K and 42--47 K,
respectively, so the independent coefficient refit rejects those
configurations.

## 13. Nominal held-out validation

| Phase | Sensor RMSE | Steady MAE | t90 MAE | Axial RMSE | Mid radial RMSE | Deep radial RMSE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| heating training | 83.10 K | 81.40 K | 630 s | 46.38 K | 28.85 K | 46.56 K |
| held-out heating | 81.57 K | 77.80 K | 717 s | 44.28 K | 29.70 K | 46.73 K |
| cooling | 32.73 K | 36.87 K | 607 s | 12.63 K | 5.47 K | 11.50 K |

Both radial signs remain wrong in all 15 heating cases.  The representative
screen-to-nominal maximum final-sensor change is 6.94 K.

## 14. What the flow hypothesis explains

The common-pressure network predicts a core/side per-channel flow ratio of
about 4.26 at 18.7 L/min and 1.81--1.95 at 4.53 L/min.  This confirms:

1. the rear groove can strongly starve outer channels;
2. temperature-dependent hydraulic feedback is active; and
3. the quadratic flow contrast weakens as total flow falls.

However, the calculated side solid remains 1.8--8.6 K colder than the center,
while the measured side thermocouples are 20--27 K hotter at 58 mm and remain
hotter at 107 mm.  Maldistribution therefore explains the direction of flow
redistribution and profile convergence, but not the measured temperature
ordering.

The interior central-channel gas is generally hotter than the local wall
downstream.  Changing the wall/gas observation fraction cannot create the
measured radial offsets.

## 15. Decision and next boundary

V14 is rejected for coefficient extraction.  The next defensible refinement
is a physical exterior-perimeter wall state representing the shallow-dip side
thermocouple location separately from the orbit-averaged edge channel cell.
It must partition, not add, outer-web mass and must retain the same centered
source, felt/adaptor topology and square-channel flow network.

Exact transverse locations of T9/T10 and an approximate side-dip depth would
substantially reduce ambiguity in that model.
