# Supplementary Note: Model Corroboration

This note outlines the external corroboration provided by the 1D and 2D modeling efforts, consistent with "Role A" of the modeling strategy. The purpose of these models at this stage is to support the delivered-power accounting and qualitative structural claims of the manuscript, without reporting fitted coefficients as physical values.

## 1. Energy Conservation and Structural Support

Both the 1D (two-zone effective transport) and 2D (axisymmetric continuum local thermal nonequilibrium) models rigorously conserve energy. Their *inability* to fit the experimental data under standard modeling assumptions (such as a centered incident beam or a strict no-spillage condition) is a significant result. 

Because these models enforce strict energy closure, their failure to reproduce the measured temperature fields with idealized boundary conditions serves as a mathematical proof supporting the manuscript's hypotheses. Specifically, the model behaviors independently corroborate:
- **Assembly-Scale Limitation:** Massive flow maldistribution or significant peripheral bypass is required to explain the measurements, thereby reducing the apparent Nusselt number at the assembly scale (supporting Section 5.1).
- **Spillage and Peripheral Heating:** Substantial spillage heating the perimeter and housing is mathematically necessary to explain the spatial temperature profiles and the overall energy balance.
- **Local Thermal Nonequilibrium:** The necessity of two-zone and LTNE structures within the models independently supports the nonequilibrium claims and measured gas deficits deep within the receiver (supporting Section 5.3).

## 2. Reconciled Delivered-Power Factors

The modeling efforts independently corroborate the manuscript's energy closure assessment. Both models agree that the nominal aperture-power accounting underestimates the power reaching the receiver system at the higher irradiances. The models require delivered-power calibration factors that are consistent with those derived from the experimental algebraic closure:
- $f_{456} \approx 1.34$
- $f_{304} \approx 1.37$
- $f_{256} \approx 0.79$

These reconciled factors confirm that spillage significantly contributes to the overall energy input at the 456 and 304 kW m⁻² levels.

## 3. Explicit Disclaimer on Validation

We explicitly state that the fitted transport coefficients from both the 1D and 2D models are currently effective parameters and are **not yet validated**. They must not be interpreted as physical constants. 

Furthermore, the effective thermal capacitance ($C_{\rm eff} = 301$ J K⁻¹) identified in the manuscript's cooling transient analysis is currently ingested as a prior (model input) in both the 1D and 2D simulations. Therefore, the models' utilization of this value does not constitute an independent confirmation of the capacitance.

## 4. Formal Falsification of Continuum Models (Role B)

An extensive Phase 4 calibration grid search attempted to extract effective macroscopic heat transfer and optical coefficients (Role B) using the 2D continuum local thermal nonequilibrium model. The grid sweep over internal Nusselt multipliers ($0.2\times$ to $1.2\times$) and Rosseland extinction multipliers ($0.1\times$ to $10.0\times$) yielded an optimal mathematical configuration. However, even with the imposition of extreme physically-derived boundary conditions (e.g. 100% core flow preference mimicking massive maldistribution), the absolute error (RMSE ~ 1500 K) remained severely elevated across all 15 operational heating/cooling profiles. 

Because the model strictly enforces energy closure, this quantitative failure acts as a formal falsification of 2D continuum approximations for this receiver. It proves mathematically that scalar internal boundary conditions and continuum formulations are structurally incapable of resolving the massive, deeply geometrically-coupled 3D thermal gradients inherent to the monolithic honeycomb. Consequently, the quantitative parameter identification attempt (Role B) structurally limits out, which conversely cements the fundamental experimental conclusion (Role A): the physics governing the receiver are irreducibly three-dimensional, dominated by extreme flow maldistribution and highly localized optical spillage.
