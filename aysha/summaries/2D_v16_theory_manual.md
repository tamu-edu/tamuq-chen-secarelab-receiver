# 2D v16 theory manual

## Purpose

V16 is the cooling-first identification stage for the explicit exterior SiC
wall introduced in v15.  It tests whether the inherited receiver/felt loss
path, rather than an omitted SiC thermal mass, prevents the model from
reproducing the side-wall temperatures.

## Physical model

The thermal states and conservation equations are the v15 equations.  The
exterior wall is carved from the existing SiC mass and is coupled to:

- the residual exterior-channel solid;
- the alumina-silicate felt; and
- the solid alumina adaptor over the rear overlap.

V16 separates the skin/felt interface multiplier from the historical v14
outer-channel/felt multiplier.  This is required because the measured state
changed in v15, so the earlier fitted multiplier cannot be assumed
transferable.

For axial cell `j`,

`Qskin-felt,j = Scontact Aext,j
                (Tskin,j - Tfelt,j) / Rseries,j`,

where `Scontact` is global and

`Rseries = lpitch/(2 kSiC) + drfelt/(2 kfelt)
           + Rcontact`.

The felt conductivity and heat capacity are

`kfelt(T) = Sk kfelt,reference(T)`

and

`Cfelt = SC rho_felt cp_felt Vfelt`.

`Scontact`, `Sk`, and `SC` are common to every experiment.

## Identification sequence

1. Freeze the v14 hydraulics, centered optics, Graetz law, power scales,
   hardware, and observation filters.
2. Freeze the v15 skin thickness at 0.05 mm and the channel/skin spreading
   multiplier at 1.0.
3. Fit the three global felt-path quantities to C69 and C80 only.
4. Use T8, T12, T11, and T2 as the primary objective.  Use T9, T10, and T3
   only as guard sensors.
5. Evaluate C81 without refitting.
6. Freeze the cooling-derived quantities and rerun all heating cases.

The parameter grid and numerical acceptance criteria were recorded in
`summaries/journal.2D.md` before the calibration was run.

## Important interpretation rule

The observation of full but non-firm felt contact motivates allowing a lower
interface multiplier; it does not predetermine the fit.  If the cooling data
select stronger conductance, that is evidence that the effective macroscopic
path includes more than mechanical contact alone or that another cooling
path remains unresolved.

## Outcome

V16 completed cooling calibration, held-out C81 confirmation, nominal-mesh
re-ranking, frozen-parameter validation, and profile review.

The first screen selected:

- skin/felt contact multiplier: 1.0;
- felt conductivity multiplier: 2.4;
- felt heat-capacity multiplier: 1.5.

It gave primary-sensor RMSE of 23.45 K on C69/C80 and 23.06 K on held-out
C81.  Because conductivity was at the search ceiling and C69 changed by
24.61 K on the nominal mesh, this selection was not accepted.

The conductivity range was extended and the leading candidates were
re-ranked on the nominal mesh.  The nominal selection became:

- skin/felt contact multiplier: 0.30;
- felt conductivity multiplier: 7.20;
- felt heat-capacity multiplier: 1.50.

The reduced interface multiplier agrees qualitatively with non-firm felt
contact.  The conductivity scale does not represent a credible identified
felt property: it remains at the extended upper bound and corresponds to
about 0.36 W/m/K near room temperature.  It is more plausibly compensating
for another cooling path.

The selected model also fails screen-to-nominal transfer, with maximum
final-sensor changes of 28.68 K for C69 and 23.67 K for E72.

Definitive nominal metrics were:

| Metric | Heating held out | Cooling |
|---|---:|---:|
| mean sensor RMSE | 124.91 K | 33.59 K |
| axial RMSE | 83.91 K | 9.25 K |
| mid-radial RMSE | 37.63 K | 5.06 K |
| deep-radial RMSE | 49.87 K | 12.60 K |

The four targeted cooling sensors improve from 26.87 to 23.12 K aggregate
RMSE, but the other cooling states and all heating predictions do not
transfer.  Heating radial signs remain 0/15 at both depths.

V16 is therefore **rejected for coefficient extraction**.  It shows that a
cooling-only felt-path fit can improve selected sensors, but only through a
mesh-sensitive, boundary-seeking effective conductivity that strongly
damages heating.

The next recommended model change is an axial, not radial, optical
refinement.  The measured side temperature peaks at 58 mm in 10 of 15
heating experiments, while v16 peaks at 5 mm in all 15.  A power-conserving
two-component axial source should represent near-entrance wall absorption
plus deeper channel-radiation absorption while retaining the centered radial
beam.  Its shared source split and deep extinction length should be fitted
only to the verified 5/58/107 mm axial profiles before radial validation.

### Power-refit correction

The heating metrics above used the v14 group power multipliers unchanged.
This was a deliberate transfer test, but not a final heating calibration:
the v16 cooling fit changed the assembly loss path, so the absorbed-power
multipliers must subsequently be re-identified.

The pooled parity plot contains cooling overprediction at low temperature
and heating underprediction at high temperature.  It does not contradict
the heating transient plots.  The corrective stage fits one nominal-mesh
power multiplier for each irradiance group using training heating cases,
with all cooling, hydraulic, spatial, and heat-transfer quantities frozen.
The held-out heating cases and all cooling cases are then rerun.  This stage
tests absolute power level only; it cannot validate the axial source law by
itself.

The completed refit selected group multipliers of 1.65, 1.80, and 1.25 for
the 456, 304, and 256 kW/m2 groups, respectively.  Held-out mean sensor RMSE
improved from the frozen-power value of 124.91 K to 61.46 K.  The corrected
parity plot now places heating points on both sides of the 1:1 line.

Spatial errors did not improve:

- held-out axial RMSE: 107.75 K;
- held-out mid-radial RMSE: 43.72 K;
- held-out deep-radial RMSE: 55.70 K;
- correct radial signs: 0/15 at both depths.

For E72, refitted side temperatures are about 688, 648, and 599 K at 5, 58,
and 107 mm, while the measurements are 738, 798, and 734 K.  Refitting total
power corrects much of the temperature level but makes the incorrect
front-loaded axial slope more evident.

Accordingly, the final v16 rejection is not based on unrefitted power.  The
power correction was necessary and successful for absolute level; the
remaining failure is spatial, mesh-sensitive, and consistent with an
inadequate single-component axial absorption law.
