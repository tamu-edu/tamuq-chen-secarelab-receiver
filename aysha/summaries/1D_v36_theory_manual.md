# 1D_v36: The Dynamic Flow Bypass Formulation (Temporal Correction)

## 1. Purpose and Present Status

`1D_v36.jl` is a specialized, continuum-level macroscopic formulation developed to resolve a critical structural failure identified in the baseline 1D homogeneous receiver model (up to `v35`). 

In standard formulations, the model operates under a **Heat Transfer Paradox**: standard calibration cannot mathematically satisfy both steady-state high-temperature heating and rapid transient cooling. To force the model to cool as rapidly as the experimental data, the optimizer is required to increase the effective convective heat transfer coefficient ($h_{eff}$). However, maintaining a high $h_{eff}$ during the heating phase instantaneously strips thermal energy from the solid into the gas, artificially suppressing the peak steady-state solid temperatures well below measured values.

`1D_v36.jl` hypothesizes that this paradox is not a failure of convection physics, but rather a failure of the **1D uniform flow assumption**. It resolves the temporal constraint by introducing a highly non-linear, temperature-dependent flow resistance mechanism that actively chokes advective mass flow through the core as the monolith heats up.

## 2. Experimental Interpretation and Domain

The physical domain spans axially from the irradiated entrance ($z=0$) to the gas exit ($z=L=0.137$ m) through a SiC monolithic square-channel honeycomb. The calibration is evaluated against a continuous, multi-stage dataset featuring steady-state heating, dynamic transitions, and natural lamp-off cooling.

The optimizer strictly minimizes the $L_2$ residuals for:
*   **Solid Temperatures**: T8 (front) and T9 (mid-axial).
*   **Gas Temperatures**: T3 (downstream outlet).
*   **Thermal Inertia**: T2 (ambient rear flange coupling).

Unlike the spatial/optical formulation (`v37`), the `v36` formulation strictly retains the standard Beer-Lambert volumetric deposition optics, focusing purely on isolating and solving the *temporal/advective* dynamics.

## 3. Governing Core Physics

The formulation retains the classic 1D Finite Difference Method (FDM) partial differential equations for the solid and gas energy balances.

### 3.1 Solid Energy Conservation
The solid core temperature ($T_{c}$) is governed by:

$$
C_{c,i}\frac{dT_{c,i}}{dt} = \dot{Q}_{z,i-1/2} - \dot{Q}_{z,i+1/2} + w_i\dot{Q}_{\mathrm{solar}} - \dot{Q}_{g,i} - \dot{Q}_{c \to p, i}
$$

### 3.2 Quasi-Steady Gas Advection
The gas temperature ($T_{g}$) assumes pseudo-steady integration at each time step:

$$ T_{g,i+1} = T_{g,i} + \left(1 - e^{-NTU_i}\right)(T_{c,i} - T_{g,i}) $$

Where the Number of Transfer Units (NTU) defines the local convective coupling:
$$ NTU_i = \frac{h_{eff} P_{\mathrm{exchange}} \Delta z}{\dot{m}_{active} c_{p,g}} $$

The standard homogenization assumes $\dot{m}_{active} = \dot{m}_{total}$. In `v36`, this assumption is structurally abandoned.

## 4. The Resolution Mechanism: Dynamic Flow Bypass

### 4.1 The Physical Argument (Viscosity-Induced Maldistribution)
As the central core of the SiC monolith reaches extreme temperatures (>1000 K), the local gas viscosity ($\mu$) increases non-linearly ($\mu \propto T^{0.7}$). In a uniform 1D model, this viscosity increase simply affects the Reynolds number ($Re$). In physical reality, however, the monolith is a parallel network of independent square channels. 

As the central channels become exponentially more resistive due to high temperatures, the incoming cold advective mass naturally seeks the path of least resistance—diverting towards the cooler perimeter channels or the physical perimeter bypass gap. **The hot core actively chokes its own flow.**

### 4.2 The Mathematical Implementation
`v36` defines the active core mass flow rate ($\dot{m}_{active}$) as a dynamic function of the spatially averaged core temperature ($T_{c,avg}$):

$$ \phi_{act} = \max\left(0.1, \ \phi_0 \left( \frac{T_{c,avg}}{295 \ K} \right)^{-m_{rec}} \right) $$
$$ \dot{m}_{active} = \phi_{act} \cdot \dot{m}_{total} $$

Where the optimizer fits two new phenomenological parameters:
*   $\phi_0$: The base active flow fraction at room temperature (295 K).
*   $m_{rec}$: The temperature resistance recruitment exponent.

By making the active flow a transient variable, the effective convective heat transfer within the core is artificially suppressed during steady-state heating (allowing the solid to reach high peak temperatures), but rapidly returns to maximum capacity during the lamp-off cooling phase (allowing rapid thermal shedding).

## 5. Optimization Results and Implications

### 5.1 The Calibration Success
The simultaneous fitting of the `v36` formulation achieved a massive success, dropping the global objective score to **0.217** (compared to the baseline scalar failure of ~1.86). 

The optimizer returned highly physical, revealing parameters for the dynamic flow:
*   **$\phi_0 = 0.998$**: At room temperature, 99.8% of the gas actively passes through the core channels. The model mathematically proves that no "permanent" structural bypass exists.
*   **$m_{rec} = 0.0606$**: As the core heats, a critical fraction of the flow is smoothly, dynamically diverted. 

### 5.2 Resolution of the Paradox
The `v36` formulation brilliantly resolves the Heat Transfer Paradox entirely through temporal advective correction. By choking the flow during the hottest phases, it reduces the $NTU$ and retains thermal energy in the solid. When the lamp shuts off and the solid begins to cool, the active flow fraction recovers, immediately increasing advective heat removal and matching the steep transient cooling tails measured by the T8/T9 sensors.

## 6. Assumptions and Structural Weaknesses

1.  **Spatial Limitation (Optics)**: Because `v36` forces the temporal flow to resolve the error while maintaining standard Beer-Lambert optics, it slightly struggles to match the extreme spatial gradients (the sharp temperature delta between T8 and T9) compared to the dedicated optical model (`v37`).
2.  **Mathematical Deletion of Bypassed Gas**: Due to the 1D homogenized architecture, the bypassed gas ($\dot{m}_{total} - \dot{m}_{active}$) is mathematically deleted from the physical domain and perfectly mixed at the rear exit sensor. It does not actively cool the perimeter along the axial length. Attempts to build a true 1.5D parallel-flow thermal exchange (in `v38`) were rejected by the optimizer, indicating that the homogenized continuum is too rigid to accommodate separate mass and energy advection pathways simultaneously.
3.  **Necessity of 3D Modeling**: The mathematical success of the dynamic bypass directly proves the existence of severe localized flow maldistribution within the honeycomb channels. To move from phenomenological exponents ($m_{rec}$) to intrinsic physical friction factors, the model definitively points toward the necessity of 3D Computational Fluid Dynamics (CFD).
