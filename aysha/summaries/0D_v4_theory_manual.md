# Theory Manual: Two-State Lumped Parameter Solar Receiver Model with Hysteresis (v4)

This document provides a comprehensive, technically rigorous description of the mathematical formulation, geometric parameters, boundary conditions, and evaluation metrics of the updated 0D lumped parameter model (`0D_v4.jl`). This version introduces an unresolved axial thermal profile state to correctly reproduce the physical hysteresis observed in cooling experiments.

---

## 1. Physical System Description

The solar receiver setup consists of a high-temperature Silicon Carbide (SiC) honeycomb monolith receiver mounted in a solar simulator cavity. While previous model versions utilized 5 discrete lumped domains, `0D_v4.jl` consolidates the energy balance into an equivalent effective capacitance model ($C_{\mathrm{eff}}$), shifting the analytical focus toward capturing the correct temporal lag between the bulk ceramic temperature and the gas outlet temperature.

Because the fluid residence time in the channels is orders of magnitude shorter than the solid thermal response, the channel fluid is modeled algebraically.

To correctly model the cooling trajectories, the system tracks two ordinary differential equations:
1. $T_{s,\mathrm{avg}}$: The bulk solid average temperature.
2. $a^*$: A hidden axial profile state representing the unresolved thermal gradient that drives gas hysteresis.

---

## 2. Mathematical Formulation (ODE Governing Equations)

The transient thermal response of the system is modeled by the following two-state ODE system.

### 2.1 Node 1: Bulk Solid Average ($T_{s,\mathrm{avg}}$)
The main energy balance on the effective receiver domain is:
$$C_{\mathrm{eff}} \frac{dT_{s,\mathrm{avg}}}{dt} = \eta_{\mathrm{eff}} I_0 A_{\mathrm{irr}} - \dot{Q}_{g} - \dot{Q}_{\mathrm{loss,lin}} - \dot{Q}_{\mathrm{loss,rad}}$$

Where:
*   $C_{\mathrm{eff}} = \gamma_C \cdot \rho_s C_{p,s}(T_s) V_s$ is the tuned effective heat capacity.
*   $\eta_{\mathrm{eff}}$ is the effective absorbed solar fraction.
*   $\dot{Q}_{g} = \dot{m} C_{p,g} (T_{g,\mathrm{out}} - T_{\mathrm{amb}})$ is the enthalpy rise of the gas.
*   $\dot{Q}_{\mathrm{loss,lin}} = K_{\mathrm{lin}} (T_{s,\mathrm{avg}} - T_{\mathrm{amb}})$ is the generalized linear heat loss (conduction/convection).
*   $\dot{Q}_{\mathrm{loss,rad}} = \chi\varepsilon_{\mathrm{rad}} \cdot \sigma A_{\mathrm{rad}} (T_{s,\mathrm{avg}}^4 - T_{\mathrm{amb}}^4)$ is the lumped radiative loss.

### 2.2 Node 2: Unresolved Axial Profile ($a^*$)
To generate physical hysteresis without artificial switches, an internal thermal profile state is tracked:
$$\tau_a \frac{da^*}{dt} + a^* = R_q^* \cdot \eta_{\mathrm{eff}} I_0 A_{\mathrm{irr}} - R_g^* \cdot \dot{Q}_g$$

Where:
*   $\tau_a$ is the characteristic time for axial conduction to redistribute heat.
*   $R_q^*$ governs the profile generation due to incident solar flux.
*   $R_g^*$ governs the profile flattening due to convective heat removal by the gas.

### 2.3 Algebraic NTU Gas Node ($T_{g,\mathrm{out}}$)
The gas outlet temperature is resolved algebraically using the NTU method. However, the gas is exposed to an effective wall temperature that includes the profile state:
$$T_{g,\mathrm{out}} = T_{g,\mathrm{in}} + \varepsilon_{\mathrm{eff}} \left(T_{s,\mathrm{avg}} + a^* - T_{g,\mathrm{in}}\right)$$

Where $\varepsilon_{\mathrm{eff}} = 1 - \exp(-NTU_{\mathrm{eff}})$.

---

## 3. Non-Dimensional Heat Transfer Correlation (Nu-Re-Pr)

The effective gas-side heat transfer coefficient inside the receiver channels is modeled rigorously using a non-dimensional Nusselt number correlation. This allows the model to correctly adapt to massive temperature variations in the fluid properties ($\mu_f, k_f, c_{p,f}$) across the monolith length.

The Nusselt number is parameterized using three fitted exponents ($A, B, C$):
$$
Nu_f(T) = 10^A \cdot Re_f(T)^B \cdot Pr_f(T)^{10^C}
$$

Where the instantaneous fluid Reynolds and Prandtl numbers are dynamically evaluated at the local bulk film temperature:
$$
Re_f(T) = \frac{\dot{m} D_h}{A_{\mathrm{chnl,frt,all}} \mu_f(T)}
$$
$$
Pr_f(T) = \frac{c_{p,f}(T) \mu_f(T)}{k_f(T)}
$$

The effective heat transfer coefficient is then directly mapped to the Number of Transfer Units ($NTU_{\mathrm{eff}}$):
$$
h_{\mathrm{eff}} = \frac{Nu_f(T) \cdot k_f(T)}{D_h}
$$
$$
NTU_{\mathrm{eff}} = \frac{h_{\mathrm{eff}} A_{\mathrm{exchange}}}{\dot{m} c_{p,f}(T)}
$$

This purely physical representation explicitly tracks the geometry ($D_h$, $A_{\mathrm{exchange}}$) and flow scaling without relying on lumped empirical correction factors.

---

## 4. Three-Stage Hysteresis Identification Strategy

Because global optimization algorithms (e.g., `MLSL_LDS`) can struggle to disentangle correlated dynamic parameters if all experiments are lumped together, the 0D model utilizes a strict **3-Stage Sequential Optimization Pipeline** to isolate and identify the hysteresis physics:

### Stage 1: Base Heating Fit
The optimizer strictly uses the **15 heating experiments** to fit the 7 physical base parameters (`γ_C`, `K_lin`, `χε_rad`, `A`, `B`, `C`, `η_eff`). During this stage, the dynamic hysteresis state is explicitly bypassed ($R_q^* = 0$, $R_g^* = 0$). This securely estimates the core properties from a highly robust dataset.

### Stage 2: Cooling Hysteresis Fit
The 7 non-solar parameters from Stage 1 are locked. The optimizer then exclusively uses the **3 cooling experiments** (where $I_0 = 0$) to identify the two hysteresis-specific variables:
* $\tau_a$: The characteristic time for axial profile redistribution.
* $R_g^*$: The gas heat-removal resistance altering the temperature profile.
Because the base properties are locked, the optimizer is strictly constrained to use these parameters solely to explain the dynamic cooling lag.

### Stage 3: Heating Hysteresis Refinement
Using the locked $\tau_a$ and $R_g^*$ from Stage 2, the optimizer is unleashed across the **heating dataset** one last time. It identifies the final solar resistance parameter $R_q^*$ while allowing the 7 base parameters to shift tightly (within ±20% of their Stage 1 optima). This final refinement ensures the core properties adapt to perfectly accommodate the newly active hysteresis physics, providing the best possible fit to the full heating and cooling loops.

---

## 5. Quantitative Thesis Metrics

The model performance is evaluated against transient and steady-state experimental measurements using the following validation metrics:

1.  **Steady-State Magnitude Error ($\Delta T_{ss}$)**
2.  **Thermal Lag ($\Delta t_{90\%}$)**
3.  **Energy Partitioning Ratio ($R_{leak}$):**
    $$R_{leak} = \frac{T_{exit,ss} - T_{amb}}{T_{solid,ss} - T_{amb}}$$
4.  **Gas-to-Solid Gap ($Gap_{ss}$):**
    $$Gap_{ss} = T_{gas,ss} - T_{solid,ss}$$
