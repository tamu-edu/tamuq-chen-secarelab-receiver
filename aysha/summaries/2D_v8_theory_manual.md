# Theoretical Manual & Mathematical Formulation of `2D_v8.jl`: Corrected-Geometry Conservative 2D Axisymmetric Macro-ECM Solar Receiver Model

**Authors**: AGY Advanced Agentic Coding Team  
**Workspace**: `d:\kkakosim\github\tamuq-chen-secarelab-receiver\aysha`  
**Target Codebase**: `2D_v8.jl`  
**Date**: July 29, 2026  

---

## 1. Introduction & Physical System Overview

The `2D_v8.jl` model is an axisymmetric, multi-domain, continuum representation—referred to as the Entire Converter Model (Macro-ECM)—of a high-temperature volumetric solar receiver. The physical core consists of a square-channel monolithic silicon carbide (SiC) honeycomb absorber irradiated by concentrated solar energy entering its front aperture, cooled by a pressurized or atmospheric air stream flowing through $N_{\text{ch}} = 100$ parallel square channels.

```
                  Incident Concentrated Solar Flux I_0(r)
                              ↓↓↓↓↓↓↓↓↓↓↓↓
              ┌──────────────────────────────────────────┐
              │ Front Surface Aperture (z = 0)          │
              ├──────────────────────────────────────────┤  r = R_rec (SiC Disc Boundary)
              │                                          │ ─── Contact Resistance R_contact
              │   Monolithic Porous SiC Receiver Core    │
              │   100 Square Channels (Dh = 1.5 mm)      │  r = R_felt (Alumina Felt Insulation)
              │   Length L_rec = 137 mm                  │
              │                                          │  r = R_case (Aluminum Casing Outer Shell)
              ├──────────────────────────────────────────┤
              │ Rear Exit Face (z = L_rec = 137 mm)      │
              └────────────────────┬─────────────────────┘
                                   │ Mixing Cup Gas Exit
                                   ▼
                      Alumina Rear Exit Tube (L = 150 mm)
                      Thermostatted Flange (298.15 K)
```

### Key Physical & Geometric Invariants in `2D_v8`:
1. **Square-to-Disc Area Equivalence**:
   The physical receiver has a square frontal cross-section of $19.0\text{ mm} \times 19.0\text{ mm}$ ($A_{\text{frt}} = 3.61 \times 10^{-4}\text{ m}^2$). The equivalent 2D disc radius is defined by exact frontal area equality:
   $$R_{\text{rec}} = \sqrt{\frac{W_{\text{rec}}^2}{\pi}} = \sqrt{\frac{(19.0\text{ mm})^2}{\pi}} \approx 10.7196\text{ mm}$$
2. **Exact Porosity & Solid Mass Conservation**:
   For $N_{\text{ch}} = 100$ channels of width $D_h = 1.5\text{ mm}$, the total flow cross-sectional area is $A_{\text{flow}} = 100 \times (1.5\text{ mm})^2 = 2.25 \times 10^{-4}\text{ m}^2$. The receiver porosity is:
   $$\varepsilon = \frac{A_{\text{flow}}}{A_{\text{frt}}} = \frac{2.25 \times 10^{-4}\text{ m}^2}{3.61 \times 10^{-4}\text{ m}^2} = \frac{225}{361} \approx 0.62327 \quad (62.33\%)$$
   The net SiC solid volume is $V_{\text{solid}} = A_{\text{frt}} (1 - \varepsilon) L_{\text{rec}} = 1.863 \times 10^{-5}\text{ m}^3$. With dense SiC solid density $\rho_{\text{SiC}} = 2150\text{ kg/m}^3$, the total active receiver mass is exactly:
   $$m_{\text{rec}} = \rho_{\text{SiC}} V_{\text{solid}} = 2150 \times 1.863 \times 10^{-5} \approx 40.05\text{ g}$$
   restoring the exact physical mass measured experimentally.

---

## 2. Multi-Domain Continuum Governing Equations

The spatial domain is divided into three concentric radial zones across the receiver length ($z \in [0, L_{\text{rec}}]$), plus a 1D rear alumina exit tube domain ($z \in [0, L_{\text{rear}}]$):

1. **SiC Receiver Core Disc**: $r \in [0, R_{\text{rec}}]$ (Porosity $\varepsilon \approx 0.6233$)
2. **Alumina Felt Thermal Insulation**: $r \in [R_{\text{rec}}, R_{\text{felt}}]$ ($R_{\text{felt}} = 57.0\text{ mm}$, Porosity $\varepsilon = 0$)
3. **Aluminum Casing Outer Shell**: $r \in [R_{\text{felt}}, R_{\text{case}}]$ ($R_{\text{case}} = 75.0\text{ mm}$, Porosity $\varepsilon = 0$)
4. **Alumina Rear Exit Tube**: $r \in [r_{i,\text{rear}}, r_{o,\text{rear}}]$ ($r_i = 6.5\text{ mm}, r_o = 8.0\text{ mm}$, Length $L_{\text{rear}} = 150.0\text{ mm}$)

### 2.1 Solid-Phase 2D Heat Conduction Equation
For any cell $(i, j)$ in the solid grid, the transient energy conservation equation is:

$$\rho_i C_{p,i}(T_s) \frac{\partial T_{s,i,j}}{\partial t} = \frac{1}{r} \frac{\partial}{\partial r} \left( r k_{r,i,j}^{\text{eff}} \frac{\partial T_s}{\partial r} \right) + \frac{\partial}{\partial z} \left( k_{z,i,j}^{\text{eff}} \frac{\partial T_s}{\partial z} \right) + q_{\text{solar},i,j}^{\prime\prime\prime} - q_{\text{conv},i,j}^{\prime\prime\prime}$$

where:
- Inside the SiC receiver disc ($i \le N_{r,\text{rec}}$):
  - Solid density $\rho_i = \rho_{\text{SiC}} = 2150\text{ kg/m}^3$.
  - Specific heat capacity $C_{p,i}(T_s)$ is the temperature-dependent SiC heat capacity:
    $$C_{p,\text{SiC}}(T_s) = 1110.0 + 0.15 T_s - \frac{4.2 \times 10^5}{T_s^2} \quad [\text{J/kg}\cdot\text{K}]$$
- Inside Alumina Felt ($N_{r,\text{rec}} < i \le N_{r,\text{rec}} + N_{r,\text{felt}}$):
  - $\rho_{\text{felt}} = 140\text{ kg/m}^3$, $C_{p,\text{felt}} = 1360\text{ J/kg}\cdot\text{K}$.
- Inside Aluminum Casing ($i > N_{r,\text{rec}} + N_{r,\text{felt}}$):
  - $\rho_{\text{case}} = 2700\text{ kg/m}^3$, $C_{p,\text{case}} = 900\text{ J/kg}\cdot\text{K}$, $k_{\text{case}} = 205\text{ W/m}\cdot\text{K}$.

---

### 2.2 Anisotropic Thermal Conductivities & Rosseland Radiative Diffusion

#### Radial Effective Thermal Conductivity $k_{s,r}^{\text{eff}}$
Heat transfer in the radial direction must cross monolith web walls and channel voids. It is modeled using a tortuosity scaling of the temperature-dependent dense SiC thermal conductivity:

$$k_{s,r}^{\text{eff}}(T_s) = \chi_r \cdot k_{\text{SiC}}(T_s)$$

where $\chi_r$ is the fittable radial conductivity multiplier, and $k_{\text{SiC}}(T_s)$ follows the COMSOL alpha-polycrystalline SiC relation:

$$k_{\text{SiC}}(T_s) = 191.9216 - 0.3261784 T_s + 2.739462 \times 10^{-4} T_s^2 - 7.70926 \times 10^{-8} T_s^3 \quad [\text{W/m}\cdot\text{K}]$$

#### Axial Effective Thermal Conductivity $k_{s,z}^{\text{eff}}$ & Rosseland Radiative Diffusion
Along the continuous axial web walls ($z$-direction), heat propagates via solid conduction parallel to the channels, augmented at high temperatures ($T_s > 800\text{ K}$) by thermal infrared radiative diffusion through the open channel voids:

$$k_{s,z}^{\text{eff}}(T_s) = \chi_z \cdot k_{\text{SiC}}(T_s) + k_{\text{rad}}(T_s)$$

The channel radiative conductivity $k_{\text{rad}}(T_s)$ is governed by the Rosseland radiative diffusion approximation:

$$k_{\text{rad}}(T_s) = \frac{16 \sigma T_s^3}{3 \beta_{\text{rad}}}$$

where $\sigma = 5.67037 \times 10^{-8}\text{ W/m}^2\text{K}^4$ is the Stefan-Boltzmann constant, and $\beta_{\text{rad}}$ is the infrared extinction coefficient of the channel array ($\beta_{\text{rad}} \in [20, 1000\text{ m}^{-1}]$).

---

## 3. Local Thermal Non-Equilibrium (LTNE) Fluid-Solid Coupling

Air flows through the $N_{\text{ch}} = 100$ parallel channels under Local Thermal Non-Equilibrium (LTNE) conditions, where the local solid wall temperature $T_{s,i,j}$ and fluid bulk temperature $T_{g,i,j}$ differ.

### 3.1 Gas Mass Flow Rate Allocation & Conservation
The total measured air mass flow rate is:

$$\dot{m}_{\text{total}} = \rho_{\text{air}}(T_{\text{in}}) \cdot \dot{V}_{\text{lpm}} = \frac{1.199 \times \dot{V}_{\text{lpm}}}{60000} \quad [\text{kg/s}]$$

To model non-uniform radial velocity profiles (preferential core flow), the total mass flow is distributed across the concentric SiC radial rings $i \in [1, N_{r,\text{rec}}]$ using a radial weighting function $\psi(r_i)$:

$$\psi(r_i) = 1.0 - c_{\text{radial\_flow}} \left( \frac{r_i}{R_{\text{rec}}} \right)^2$$

The mass flow rate allocated to radial ring $i$ is strictly conservative:

$$\dot{m}_i = \dot{m}_{\text{total}} \cdot \frac{\psi(r_i) A_{\text{flow},i}}{\sum_{k=1}^{N_{r,\text{rec}}} \psi(r_k) A_{\text{flow},k}}$$

ensuring that $\sum_{i=1}^{N_{r,\text{rec}}} \dot{m}_i \equiv \dot{m}_{\text{total}}$ under all operating conditions.

---

### 3.2 Apparent Macro-ECM Convective Heat Transfer Closure
Along each channel ring $i$, the gas temperature $T_{g,i,j}$ marches downstream from the inlet temperature $T_{\text{in}}$ at $z = 0$. In each axial cell $j$ of length $\Delta z_j$:

$$\dot{m}_i C_{p,g}(T_{\text{film}}) \frac{d T_g}{dz} = h_{\text{cell},i,j} a_v A_{\text{frt},i} (T_{s,i,j} - T_{g,i,j})$$

Integrating over cell length $\Delta z_j$ yields the discrete cell effectiveness $\varepsilon_{\text{cell},i,j}$:

$$\varepsilon_{\text{cell},i,j} = 1 - \exp\left( -\frac{h_{\text{cell},i,j} P_{\text{ex},i} \Delta z_j}{\dot{m}_i C_{p,g}(T_{\text{film}})} \right)$$

$$T_{g,i,j+1} = T_{g,i,j} + \varepsilon_{\text{cell},i,j} (T_{s,i,j} - T_{g,i,j})$$

$$Q_{\text{conv},i,j} = \dot{m}_i C_{p,g}(T_{\text{film}}) (T_{g,i,j+1} - T_{g,i,j})$$

where $P_{\text{ex},i} = A_{\text{frt},i} \cdot \left( \frac{4 N_{\text{ch}} D_h}{A_{\text{frt}}} \right)$ is the total channel wetted perimeter in radial ring $i$.

#### Apparent Macro-ECM Nusselt Law
In `2D_v8.jl`, the local heat transfer coefficient $h_{\text{cell},i,j}$ is calculated from an explicit **Apparent Macro-ECM Nusselt Law** that includes thermal entry length enhancement:

$$\text{Nu}_{i,j} = A_{\text{Nu}} \cdot \text{Re}_{i,j}^{B_{\text{Re}}} \cdot \text{Pr}_{i,j}^{C_{\text{Pr}}} \cdot \left( \frac{D_h}{z_j} \right)^{1/3}$$

$$h_{\text{cell},i,j} = \text{Nu}_{i,j} \cdot \frac{k_{\text{air}}(T_{\text{film}})}{D_h}$$

where $\text{Re}_{i,j} = \frac{\dot{m}_i D_h}{A_{\text{flow},i} \mu_{\text{air}}(T_{\text{film}})}$, $\text{Pr}_{i,j} = \frac{C_{p,\text{air}} \mu_{\text{air}}}{k_{\text{air}}}$, $A_{\text{Nu}}$ is the macro-ECM prefactor, $B_{\text{Re}} = 1.440$ is the apparent Reynolds exponent, and $C_{\text{Pr}} = 0.333$. Notice that the arbitrary single-channel laminar floor ($\text{Nu}_{\text{min}} = 3.61$) is eliminated ($\text{Nu}_{\text{min}} = 0$), leaving $A_{\text{Nu}}$ as an unconfounded macroscopic parameter.

---

## 4. Optical Ray Absorption & Non-Renormalized Power Budget

Concentrated solar flux $I_0$ incident on the receiver frontal aperture is distributed radially according to a Gaussian beam profile and absorbed axially following Beer-Lambert optical attenuation.

### 4.1 Gaussian Beam Radial Distribution
The radial beam distribution $I_{\text{rad}}(r_i)$ across the receiver aperture ($r \in [0, R_{\text{rec}}]$) is:

$$I_{\text{rad}}(r_i) = \exp\left( -\frac{r_i^2}{2 \sigma_{\text{beam}}^2} \right)$$

normalized such that $\frac{1}{A_{\text{frt}}} \sum_{i=1}^{N_{r,\text{rec}}} I_{\text{rad}}(r_i) A_{\text{frt},i} = 1.0$.

---

### 4.2 Non-Renormalized Axial Beer-Lambert Attenuation
Unlike previous formulations that forced $100\%$ solar energy absorption within the finite absorber length $L_{\text{rec}}$, `2D_v8.jl` enforces true physical optical transmission:

$$w_z(j) = \exp(-\beta_{\text{opt}} z_{j}) - \exp(-\beta_{\text{opt}} z_{j+1})$$

The total fraction of light absorbed in the depth of the receiver is:

$$\eta_{\text{depth}} = \sum_{j=1}^{N_z} w_z(j) = 1 - \exp(-\beta_{\text{opt}} L_{\text{rec}})$$

The remaining unabsorbed fraction $\exp(-\beta_{\text{opt}} L_{\text{rec}})$ is transmitted through the rear of the channels.

A fraction $f_{\text{front}} = 0.20$ of the optical power is absorbed directly on the front entrance skin ($z = 0$), while the remaining fraction $(1 - f_{\text{front}})$ decays exponentially through the channel depth:

$$w_{z,\text{depth}}(j) = (1 - f_{\text{front}}) w_z(j)$$

$$w_{z,\text{depth}}(1) \leftarrow w_{z,\text{depth}}(1) + f_{\text{front}}$$

---

### 4.3 Delivered Solar Power & Rim Spillage
The total net solar power delivered to the receiver aperture is:

$$P_{\text{delivered}} = \eta_{\text{abs}} \cdot f_{\text{scale}}(I_0) \cdot I_0 \cdot A_{\text{frt}}$$

where $f_{\text{scale}}(I_0)$ is the calibrated flux scaling factor ($f_{456}, f_{304}, f_{256}$).

A fraction $f_{\text{spillage}} = 0.10$ spills outside the SiC receiver disc onto the front surface of the alumina felt insulation ($r \in [R_{\text{rec}}, R_{\text{felt}}]$):

$$P_{\text{spillage}} = f_{\text{spillage}} \cdot P_{\text{delivered}}$$

$$P_{\text{core}} = P_{\text{delivered}} - P_{\text{spillage}}$$

The volumetric heat source $q_{\text{solar},i,j}^{\prime\prime\prime}$ in cell $(i, j)$ is:

$$q_{\text{solar},i,j}^{\prime\prime\prime} = \frac{P_{\text{core}} \cdot I_{\text{rad}}(r_i) \left( \frac{A_{\text{frt},i}}{A_{\text{frt}}} \right) w_{z,\text{depth}}(j)}{A_{\text{solid},i} \Delta z_j}$$

---

## 5. Boundary Conditions, Losses & Exit Tube Dynamics

```
                           Aperture Front Losses (Radiation + Natural Convection)
                                                ↑ ↑ ↑
    -----------------------------------------------------------------------------------------
    SiC Monolith Core (r <= R_rec) | Alumina Felt Insulation (R_rec < r <= R_felt) | Casing (r <= R_case)
    -----------------------------------------------------------------------------------------
                                   | Contact Resistance R_contact                  | Outer Casing Losses
                                                                                   ↓ (Conv + Rad)
```

### 5.1 Front Surface Loss ($z = 0$)
The front surface ($j = 1$) loses heat to ambient air ($T_{\text{amb}}$) via thermal radiation and natural convection:

$$Q_{\text{front},i} = A_{\text{front},i} \cdot \left[ \varepsilon_{\text{front}} \sigma (T_{s,i,1}^4 - T_{\text{amb}}^4) + h_{\text{front}}(T_{s,i,1}, T_{\text{amb}}) (T_{s,i,1} - T_{\text{amb}}) \right]$$

where $h_{\text{front}}$ is calculated from the Churchill-Chu natural convection correlation over a horizontal plate.

---

### 5.2 Outer Casing Surface Loss ($r = R_{\text{case}}$)
The outer aluminum casing shell loses heat radially to the cavity environment:

$$Q_{\text{case},j} = (2 \pi R_{\text{case}} \Delta z_j) \cdot \left[ 10.0 (T_{s,N_r,j} - T_{\text{amb}}) + \varepsilon_{\text{case}} \sigma (T_{s,N_r,j}^4 - T_{\text{amb}}^4) \right]$$

---

### 5.3 Interface Contact Resistance ($r = R_{\text{rec}}$)
Conduction between the SiC receiver disc ($i = N_{r,\text{rec}}$) and the surrounding alumina felt insulation ($i = N_{r,\text{rec}} + 1$) includes a physical areal contact thermal resistance $R_{\text{contact}}$ ($R_{\text{contact}} = 1.0 \times 10^{-3}\text{ m}^2\text{K/W}$):

$$R_{\text{areal},i} = \frac{\Delta r_i}{2 k_i} + \frac{\Delta r_{i+1}}{2 k_{i+1}} + R_{\text{contact}}$$

$$Q_{\text{cond},r,i,j} = (2 \pi r_{\text{face},i+1} \Delta z_j) \cdot \frac{T_{s,i,j} - T_{s,i+1,j}}{R_{\text{areal},i}}$$

---

### 5.4 Rear Exit Tube & Water-Cooled Flange Boundary
Air exiting all channel rings mixes into a single mixing cup temperature $T_{g,\text{rear}}(z=0)$ at the monolith rear exit ($z = L_{\text{rec}}$):

$$T_{g,\text{rear}}(z=0) = \frac{\sum_{i=1}^{N_{r,\text{rec}}} \dot{m}_i C_{p,g}(T_{g,i,N_z+1}) T_{g,i,N_z+1}}{\sum_{i=1}^{N_{r,\text{rec}}} \dot{m}_i C_{p,g}(T_{g,i,N_z+1})}$$

The gas flows through a 1D alumina tube ($r_i = 6.5\text{ mm}, r_o = 8.0\text{ mm}$, length $L_{\text{rear}} = 150.0\text{ mm}$ discretized into $N_{z,\text{rear}} = 30$ cells), exchanging heat with the tube walls via internal forced convection.

The downstream end of the rear tube ($z = L_{\text{rear}}$) is bolted to a water-cooled flange maintained at $T_{\text{flange}} = 298.15\text{ K}$, creating a terminal conductive heat sink:

$$Q_{\text{flange}} = k_{\text{alumina}}(T_{\text{rear},end}) \cdot A_{\text{tube\_wall}} \cdot \frac{T_{\text{rear},end} - T_{\text{flange}}}{0.5 \Delta z_{\text{rear}}}$$

---

## 6. Thermocouple Sensor Mappings & Enforced Initial Conditions

Seven experimental thermocouples record temperatures during testing:

| Sensor Symbol | Physical Location | Radial Position $r$ | Axial Position $z$ | Physical Phase Measured |
| :--- | :--- | :--- | :--- | :--- |
| **$T_8$** | Perimeter skin | $r = R_{\text{rec}} = 10.72\text{ mm}$ | $z = 11.0\text{ mm}$ | SiC Solid |
| **$T_9$** | Core mid-depth | $r = 0.0\text{ mm}$ (Centerline) | $z = 58.0\text{ mm}$ | SiC Solid |
| **$T_{12}$** | Perimeter mid-depth | $r = R_{\text{rec}} = 10.72\text{ mm}$ | $z = 58.0\text{ mm}$ | SiC Solid |
| **$T_{10}$** | Core deep-rear | $r = 0.0\text{ mm}$ (Centerline) | $z = 107.0\text{ mm}$ | SiC Solid |
| **$T_{11}$** | Perimeter deep-rear | $r = R_{\text{rec}} = 10.72\text{ mm}$ | $z = 107.0\text{ mm}$ | SiC Solid |
| **$T_3$** | Exit gas stream | Rear exit tube centerline | $z = 3.0\text{ mm}$ past exit | Rear Gas Stream |
| **$T_2$** | Insulation sleeve | $r = R_{\text{rec}} + 40.0\text{ mm} = 50.72\text{ mm}$ | $z = 58.0\text{ mm}$ | Alumina Felt Insulation |

### Bilinear Spatial Interpolation for Sensor Extraction
Because sensor locations fall between cell centers, sensor predictions in `2D_v8.jl` use exact 2D bilinear spatial interpolation (`_sample_solid` and `_sample_rear_gas`), avoiding nearest-neighbor discretization artifacts.

### Enforced Sample Initial State $u_0$ (`build_initial_state_2D`)
To start transient simulations directly from raw experimental initial readings at $t_0$, `build_initial_state_2D` constructs an initial temperature field that passes **exactly through all measured thermocouple readings** using inverse sample enforcement (`enforce_sample!`).

---

## 7. Optimization Subspace & Parameter Identification

The 12-element parameter vector $\boldsymbol{\theta}$ maps as follows:

$$\boldsymbol{\theta} = \left[ A_{\text{Nu}}, B_{\text{Re}}, \chi_r, \chi_z, \sigma_{\text{beam}}, f_{\text{spillage}}, f_{456}, f_{304}, f_{256}, \beta_{\text{rad}}, R_{\text{contact}}, \beta_{\text{opt}} \right]$$

To prevent ill-conditioned parameter trade-offs, `2D_v8.jl` restricts parameter optimization to the 8-element active subspace `FIT_INDICES_2D = [1, 3, 4, 7, 8, 9, 10, 12]`:

$$\boldsymbol{\theta}_{\text{active}} = \left[ A_{\text{Nu}}, \chi_r, \chi_z, f_{456}, f_{304}, f_{256}, \beta_{\text{rad}}, \beta_{\text{opt}} \right]$$

Structural parameters ($\sigma_{\text{beam}} = 14.0\text{ mm}$, $f_{\text{spillage}} = 0.10$, $R_{\text{contact}} = 1.0 \times 10^{-3}\text{ m}^2\text{K/W}$) are fixed independently, ensuring strict mathematical identifiability.

---

## 8. Summary of Theoretical Breakthroughs in `2D_v8.jl`

```
               2D_v1 - 2D_v7 Legacy                      2D_v8 Refactored Theory
     ┌───────────────────────────────────────┐   ┌───────────────────────────────────────┐
     │ • Receiver Disc Radius = 33.9 mm      │   │ • Receiver Disc Radius = 10.72 mm     │
     │   (Overpredicted Volume & Mass)       │   │   (Exact Area 3.61 cm² & 40.0 g Mass) │
     │ • Re-normalized Optical Absorption    │   │ • Un-renormalized Optical Attenuation │
     │   (100% Absorbed in Finite Mass)      │   │   (True Physical Rear Transmission)   │
     │ • Arbitrary Single-Channel Nu Floor   │   │ • Pure Macro-ECM Apparent Nusselt Law │
     │   (Nu_min = 3.61 Confounded Fit)      │   │   (Nu_min = 0, Clean Unconfounded Fit)│
     └───────────────────────────────────────┘   └───────────────────────────────────────┘
```

1. **Volume & Mass Invariants**: Corrects frontal area to $3.61\text{ cm}^2$ and radius to $10.72\text{ mm}$, guaranteeing exact solid mass ($40.0\text{ g}$) and porosity ($0.6233$).
2. **Mass Flow Conservation**: Strictly conserves total air mass flow across channel rings and into the exit mixing cup.
3. **Physical Optical Transmission**: Retains un-renormalized Beer-Lambert absorption, modeling actual rear optical transmission.
4. **Clean Macro-ECM Closure**: Removes non-identifiable single-channel laminar floors ($\text{Nu}_{\text{min}} = 0$), yielding robust, physically meaningful macroscopic heat transfer parameters.
