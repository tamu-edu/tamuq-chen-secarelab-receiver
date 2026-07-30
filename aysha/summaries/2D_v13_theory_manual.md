# 2D v13 theory manual

## Purpose

V13 tests whether the measured side/core temperature separation comes from
two thermally distinct parts of the square monolith: channel walls that
participate strongly in gas heat transfer and side/corner/support-connected
material that participates weakly.  It retains the complete v12 hardware
inventory, loose alumina-adaptor contact, felt model, thermocouple response,
MFC mass flow, local density and velocity, DP1 calculation, equal-pressure
channel allocation and Graetz heat transfer.

V13 does not impose a hotter irradiance pattern at the receiver side and does
not change material properties by experiment.

## 1. Conservative receiver split

At every axial cell, the measured SiC heat capacity is divided into active
and side/support populations:

$$
C_a=(1-f_s)C_{\rm SiC},\qquad
C_s=f_sC_{\rm SiC}.
$$

Thus

$$
C_a+C_s=C_{\rm SiC}.
$$

The new state is one side/support temperature per axial cell.  It adds a
temperature degree of freedom, not mass or heat capacity.

## 2. Solar source

The v12 optical deposition field is unchanged.  A fraction $f_s$ of the
source enters the side population and $1-f_s$ enters the active population.
Both receive exactly the same axial and radial optical weighting:

$$
Q_{{\rm solar},s,j}=f_s Q_{{\rm solar},j}.
$$

This is a uniform inventory split and not side-weighted irradiance.

## 3. Gas heat transfer

V12 calculates the local channel heat transfer and gas-temperature rise using
its density-dependent velocity, equal-pressure channel flow distribution and
Graetz Nusselt model.  V13 preserves the total calculated heat delivered to
the gas but reallocates its solid origin:

$$
Q_{{\rm gas},s,j}=f_{g,s}Q_{{\rm gas},j},
$$

$$
Q_{{\rm gas},a,j}=(1-f_{g,s})Q_{{\rm gas},j},
\qquad 0\le f_{g,s}\le f_s.
$$

When $f_{g,s}<f_s$, the side/support material is less exposed to flowing gas
than its share of receiver mass.  The MFC mass flow itself is unchanged; no
unmeasured bypass branch is introduced.

## 4. Active-to-side exchange

The populations exchange heat through one global conductance per receiver
length:

$$
dQ_{a-s}=G'_{a-s}\,dz\,
(T_{a,\rm mean}-T_s).
$$

The active mean is solid-area weighted.  The exchange enters the two energy
balances with equal and opposite signs.

## 5. Felt and adaptor contacts

The inherited receiver/felt path is reduced in proportion to the active
inventory.  The side population has a separate perimeter contact:

$$
dQ_{s-f}=h_{s-f}(4w_{\rm rec}dz)(T_s-T_f).
$$

Because the photographed alumina adaptor supports the square receiver around
its perimeter, the loose receiver/adaptor contact over the rear 29 mm belongs
to the side/support population:

$$
dQ_{s-ad}=h_{\rm rec-ad}(4w_{\rm rec}dz)
(T_s-T_{\rm ad}).
$$

The v12 adaptor/tube and adaptor/felt paths remain unchanged.  There is no
tube/felt contact.

## 6. Axial transport and front loss

The active population carries $(1-f_s)$ of the inherited axial receiver
transport.  The side population carries the complementary solid-conduction
share.  Front radiation is evaluated separately for the two front
temperatures and areas.  Internal exchange terms cancel exactly in the
whole-assembly energy balance.

## 7. Thermocouple observations

The verified axial positions are 5, 58 and 107 mm from front to back.
T8, T12 and T11 observe the side/support field through the v12 shared side
response time.  T9 and T10 retain the v12 wall/gas observation equation, and
T3 retains the outlet response equation.

## 8. Identification rules

The first v13 test freezes every definitive v12 quantity and fits only:

| Parameter | Meaning |
| --- | --- |
| $f_s$ | side/support share of measured SiC inventory |
| $f_{g,s}$ | absolute gas-heat share assigned to that population |
| $G'_{a-s}$ | active-to-side conductance per length |
| $h_{s-f}$ | side-to-felt contact coefficient |

The nine established heating-training cases determine the profile terms,
with C69/C80/C81 used to penalize loss of the v12 cooling performance.  The
six held-out heating cases are used only after selection.  Acceptance still
requires energy/mass closure, mesh transfer, correct radial ordering,
credible transient and axial behavior, and better held-out performance.

## 9. Numerical verification

`test/smoke_2D_v13.jl` checks:

- exact conservation of receiver heat-capacity inventory;
- uniform thermal equilibrium;
- finite heated response;
- fixed standard mass flow;
- equal-pressure convergence; and
- instantaneous whole-assembly energy closure including the new population,
  felt, front and adaptor exchanges.

Calibration and validation results are appended here only after the staged
screen and nominal-mesh confirmation complete.

## 10. Completed staged fit

The initial bounded population fit selected 7.5% side mass and the maximum
tested active/side conductance.  Expanding the active bounds moved the
selection to:

| Parameter | Selected value |
| --- | ---: |
| $f_s$ | 0.025 |
| $f_{g,s}$ | 0.0125 |
| $G'_{a-s}$ | 100 W/K/m |
| $h_{s-f}$ | 5 W/m2/K |

At 800 K, the full SiC receiver capacity is 49.25 J/K and the selected side
capacity is only 1.23 J/K.  The integrated active/side conductance is
$100(0.137)=13.7$ W/K, corresponding to approximately 0.09 s for the side
capacity.  The fitted second population therefore equilibrates almost
instantaneously and has negligible inventory.  This is parameter collapse
toward v12, not evidence for a distinct side thermal mode.

Moderately separated candidates did correct the $T12-T9$ sign in 7/9
training cases and a strongly separated candidate corrected it in 9/9.
However, neither corrected the deep $T11-T10$ profile; their heating and
cooling errors rose sharply.

## 11. Independent optical and heat-transfer refit

To prevent inherited v12 parameters from determining the verdict, the
population candidates were independently screened over
$\beta_{\rm opt}=75$--400 m$^{-1}$ and $S_h=0.75$--1.80.  The three shared
power-group scales were then refitted on the nine heating-training cases.

The selected point remained the collapsed population and gave

$$
\beta_{\rm opt}=150\ {\rm m^{-1}},\qquad S_h=1.80,
$$

with power scales 1.25, 1.2075 and 0.77.  The Graetz scale and middle power
scale are at the tested upper bounds.  No meaningful population split was
rescued by this refit.

## 12. Nominal-mesh validation

| Phase | Sensor RMSE | Steady MAE | t90 MAE | Axial RMSE | Mid radial RMSE | Deep radial RMSE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| heating training | 88.00 K | 85.74 K | 735 s | 45.35 K | 27.50 K | 43.07 K |
| held-out heating | 82.60 K | 78.15 K | 820 s | 45.59 K | 28.27 K | 43.17 K |
| cooling | 31.04 K | 31.21 K | 688 s | 10.59 K | 4.96 K | 9.77 K |

Both radial signs are wrong in all 15 heating cases.  The largest
representative screen-to-nominal final-sensor change is 7.86 K, so this
failure is not caused by the calibration mesh.

V13 is rejected for coefficient extraction.  Its useful result is narrower:
a spatial distinction between side and interior channel groups can generate
the observed mid-depth direction, but a global co-located side population
cannot reproduce the full axial/radial field.

## 13. Recommended successor

The successor should be an explicit square-sector or channel-group network
with center, edge and corner channels.  Each group should solve flow from the
same pressure drop using local temperature-dependent gas properties and the
measured rear groove obstruction.  Lateral SiC conduction, the v12
felt/adaptor topology, exterior side thermocouple locations and interior
channel-probe locations should remain explicit.

This structure can test whether hotter outer channels acquire higher
hydraulic resistance, receive less flow and become hotter still.  It retains
the measured MFC total flow and does not require an unquantified bypass.
