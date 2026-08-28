# Supplementary Note: Model Corroboration

## Scope

The 1D two-zone and 2D axisymmetric LTNE models are used here as Role-A
mechanistic and identifiability studies. They test whether candidate source,
flow, observation, and thermal-inventory closures can reproduce the measured
behavior while conserving energy.

They do not currently provide validated Role-B convective, radiative, or
conductive coefficients. Fitted values are model-dependent effective
parameters and must not be reported as transferable material or channel
properties.

## Findings supported by the model studies

The model history supports the following limited conclusions:

- A uniform single-temperature or single-channel description is inconsistent
  with the measured cross-sectional and gas-solid temperature differences.
- Source placement, peripheral heating, flow distribution, rear hardware, and
  the T3 observation relation materially affect inferred transport parameters.
- Heating and cooling jointly constrain the model more strongly than either
  phase alone. Parameters that improve steady heating often degrade transient
  cooling, exposing compensation among exchange, active flow, and thermal
  inventory.
- Exact energy and enthalpy closure do not by themselves make an inverse result
  identifiable. The v20 2D study conserves energy but retains a boundary ridge
  in `UA(Re)` and an unidentified T3 observer.
- The model-based power factors and the T3-based algebraic closure show a
  similar lamp-group ordering. They use the same campaign and are not
  independent; the factors are provisional configuration-level corrections,
  not validated optical constants.

These results are consistent with assembly-scale nonuniformity and
maldistribution as plausible mechanisms. They do not uniquely determine the
relative contributions of flow bypass, cross-sectional solid-temperature
variation, optical deposition, or downstream probe bias.

## 1D status

v36 is the strongest saved 1D Role-A result. It reproduces the
manuscript-definition effectiveness envelope and approximately recovers the
measured apparent-Nu exponent. It nevertheless misses the measured inversion
threshold and deep flow-slope behavior, and important thermal-capacity
parameters sit on imposed bounds. Its capacitance total is therefore not an
independent confirmation of the experimentally identified `301 +/- 23 J/K`.

v37 is an optical/source-placement sensitivity, v38 is a rejected two-stream
variant within its tested parameterization, and v39 is not citable because its
saved objective and residual artifacts disagree. The current v41 output is
invalid for interpretation because corrected delivered-power factors and
recovered spillover power are applied cumulatively.

The 1D source parameters `scale` and `spill_capture` enter the total absorbed
power mainly through a product. They are not separately identified by the
temperature data. Future work should use one conserved total-power multiplier
and a separate core/perimeter partition rather than interpreting either fitted
source parameter alone.

## 2D status

v20 is the current trustworthy 2D endpoint. It establishes structural
non-identifiability under the available measurements and tested model class:
the T3 observation operator is not identified, the T3-free `UA(Re)` surface is
a boundary ridge, and bounded nuisance parameters do not pass the declared
heating and cooling gates.

The later v22 Phase-4 grid cannot be used as formal falsification evidence. Its
Rosseland coordinate is inactive in the saved sweep, and its fixed macroscopic
settings do not inherit the Phase-3 optimum. The resulting 1500--1800 K errors
primarily diagnose that configuration and wiring, not the limits of every
continuum formulation.

## Manuscript language

The safe manuscript statement is:

> Energy-conserving reduced-order models reproduce selected qualitative
> signatures and expose strong compensation among source placement, flow
> participation, thermal inventory, and outlet-probe physics. Under the
> available observations, the fitted transport coefficients are not uniquely
> identifiable and are retained only as model-dependent effective parameters.

Avoid claims that the models "prove" severe bypass, irreducible
three-dimensionality, or the impossibility of all continuum models. Avoid
claiming that the models independently confirm a capacitance that was imposed
through bounds or priors.

## Measurements needed for Role B

The most valuable additional constraints are:

- an independently shielded or calorimetric bulk outlet-enthalpy measurement;
- an integrated receiver-face flux and spillage map;
- a characterized T3/tube/probe transfer function; and
- pressure-drop or channel-flow-distribution measurements.

Without these data, the experimental apparent correlations remain publishable
assembly-level results, while intrinsic coefficient extraction should remain
an open objective.
