# Theoretical Manual & Mathematical Formulation of `2D_v8.jl`: Corrected-Geometry Conservative 2D Axisymmetric Macro-ECM Solar Receiver Model

**Authors**: AGY Advanced Agentic Coding Team  
**Workspace**: `d:\kkakosim\github\tamuq-chen-secarelab-receiver\aysha`  
**Target Codebase**: `2D_v8.jl`  
**Date**: July 29, 2026  

> **2D_v9 implementation addendum:** Section 9 records the implemented
> temperature-dependent velocity and DP1 hydraulic extension that inherits
> the v8 thermal/optical formulation. The cold hydraulic constraint is
> verified; the hot-path pressure closure remains provisional.

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
| **$T_8$** | Perimeter skin | $r = R_{\text{rec}} = 10.72\text{ mm}$ | $z = 5.0\text{ mm}$ | SiC Solid |
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

---

## 9. Implemented `2D_v9` Hydraulic Extension

`2D_v9` inherits the corrected geometry, energy conservation, LTNE coupling,
optics and boundary conditions of `2D_v8`. Its hydraulic objective is to make
the actual local gas velocity explicit and test the prescribed MFC mass flow
against DP1. It does not add an unquantifiable bypass branch.

The equations and numerical values in this section correspond to `2D_v9.jl`
and its DP1 runner. The reference-flow conversion, velocity calculation and
cold hydraulic coefficients are implemented and tested. They are not new
fitted thermal-transfer coefficients.

### 9.1 Verified sensor-coordinate correction

The receiver thermocouple axial planes have been verified as 5, 58 and
107 mm from the illuminated front face. T8 is therefore corrected from the
previously documented 11 mm to 5 mm in the table in Section 6. The other
receiver coordinates remain:

$$
z_{T8}=5\ \mathrm{mm},\qquad
z_{T9}=z_{T12}=58\ \mathrm{mm},\qquad
z_{T10}=z_{T11}=107\ \mathrm{mm}.
$$

The v9 observation map and cooling initialization must use these same verified
coordinates. The change is prescribed experimental geometry, not an
optimization variable.

### 9.2 MFC reference state and conserved mass flow

The Aalborg GFC manual gives the calibration reference state as 14.7 psia
(approximately 101.4 kPa absolute) and 70 degF (21.1 degC or 294.25 K).
The experimental readings were adjusted with the manufacturer's gas factors
and were separately tested with a bubble flow meter. They are treated here as
air-equivalent standard L/min; no second gas-factor correction is applied.

For dry air,

$$
p_{\mathrm{ref}}=101.4\ \mathrm{kPa},\qquad
T_{\mathrm{ref}}=294.25\ \mathrm{K},
$$

$$
\rho_{\mathrm{ref}}
=\frac{p_{\mathrm{ref}}}{R_{\mathrm{air}}T_{\mathrm{ref}}}
\simeq1.200\ \mathrm{kg\,m^{-3}},
$$

and

$$
\dot m_{\mathrm{total}}
=\rho_{\mathrm{ref}}\frac{Q_{\mathrm{std}}}{60000}.
$$

This supersedes the less explicit description in Section 3.1 but is
numerically consistent with its value of 1.199 kg/m3. The corrected MFC flow
remains the primary prescribed flow.

### 9.3 Local density, actual volumetric flow and velocity

The v8 Reynolds-number calculation is already compatible with a prescribed
mass flow:

$$
\mathrm{Re}_{i,j}
=\frac{\rho_{i,j}u_{i,j}D_h}{\mu_{i,j}}
=\frac{\dot m_iD_h}{A_{\mathrm{flow},i}\mu_{i,j}}.
$$

Density cancels from Reynolds number, but it does not cancel from local
velocity or dynamic pressure. V9 therefore evaluates

$$
\rho_{i,j}(t)
=\frac{p_{i,j}(t)}{R_{\mathrm{air}}T_{g,i,j}(t)}
\simeq
\frac{p_{\mathrm{atm}}}{R_{\mathrm{air}}T_{g,i,j}(t)},
$$

$$
\dot V_{i,j,\mathrm{actual}}(t)
=\frac{\dot m_i}{\rho_{i,j}(t)},
$$

$$
u_{i,j}(t)
=\frac{\dot m_i}
{\rho_{i,j}(t)A_{\mathrm{flow},i}}.
$$

The approximation $p_{i,j}\simeq p_{\mathrm{atm}}$ is sufficient initially
because the measured receiver-path pressure difference is only a few
millibar. Thus, standard volumetric flow is fixed by the MFC, mass flow is
conserved, and actual hot-gas volumetric flow and velocity rise as density
falls.

### 9.4 DP1 measurement definition

DP1 is a flush wall static-pressure tap in the gas path just inside the water
jacket. Its reference port is open to atmosphere, as is the receiver front.
It therefore measures static gauge pressure relative to the atmospheric
receiver face and does not contain the directional stagnation-pressure term
that a Pitot-type probe would measure.

### 9.5 Square-channel pressure loss

For fully developed laminar flow in a square channel, the Darcy friction
factor satisfies

$$
f_D\mathrm{Re}=56.91.
$$

The pressure loss through axial cell $j$ of ring $i$ is therefore

$$
\Delta p_{i,j}
=f_{D,i,j}\frac{\Delta z_j}{D_h}
\frac{\rho_{i,j}u_{i,j}^2}{2}
=28.455\frac{\mu_{i,j}u_{i,j}\Delta z_j}{D_h^2}.
$$

Summing the cells gives the ideal receiver-channel contribution:

$$
\Delta p_{\mathrm{square},i}
=\sum_{j=1}^{N_z}\Delta p_{i,j}.
$$

Parallel channel rings experience a common end-to-end pressure difference.
The first v9 implementation retains the conservative v8 ring allocation and
uses this pressure calculation as an aggregate diagnostic. A later
pressure-balanced ring allocation should be considered only after total-flow
validation; one DP1 signal cannot independently identify arbitrary radial
maldistribution.

### 9.6 Cold-$t_0$ DP1 calibration

Near-ambient $t_0$ data from nine heating runs provide the least thermally
confounded normal-configuration pressure-flow set:

```text
E67, E68, E70, E72, E74, E75, E76, E78, E80
```

The calibration uses means from the raw inclusive interval
$0\le t\le10$ s, before the normal 50-sample decimation. Selection requires
both T3 and the mean of T8--T12 to lie within 5 K of ambient.

Their empirical relation is

$$
DP1_{\mathrm{raw}}[\mathrm{mbar}]
=-0.614226
+0.0455545\,Q_{\mathrm{std}}[\mathrm{L\,min^{-1}}],
$$

with

$$
R^2=0.98144,\qquad
\mathrm{RMSE}=0.02768\ \mathrm{mbar}.
$$

The intercept is treated as the cold-dataset DP1 offset unless an independent
zero measurement supersedes it. At the MFC reference temperature of
294.25 K, the ideal fully developed square-channel slope is
0.0233408 mbar/(standard L/min). Relative to this receiver-only pressure loss,
the cold observed slope implies

$$
C_{\mathrm{hyd}}=1.95171\simeq1.95.
$$

This effective multiplier includes the square-channel resistance plus
developing-flow, inlet, outlet and short water-jacket-path effects between the
atmospheric receiver face and the DP1 tap. The provisional observation model,
using the recorded DP1 sign convention, is

$$
DP1_{\mathrm{model}}
=p_{0,\mathrm{DP1}}
+C_{\mathrm{hyd}}\Delta p_{\mathrm{square}},
$$

where

$$
p_{0,\mathrm{DP1}}\simeq-0.614226\ \mathrm{mbar},
\qquad C_{\mathrm{hyd}}\simeq1.95.
$$

No free quadratic minor-loss term is included initially. Such a term may be
added only if full cooling and heating histories show a systematic
flow-dependent residual that the temperature-corrected laminar calculation
cannot explain.

### 9.7 Permitted role of DP1 in inference

DP1 is used in the following hierarchy:

1. validate local density, velocity and temperature-dependent pressure loss;
2. test the corrected MFC standard flow as the prescribed receiver mass flow;
3. constrain, if necessary, one common global flow multiplier;
4. prevent flow-scale error from being absorbed by $A_{\mathrm{Nu}}$ and
   $B_{\mathrm{Re}}$.

DP1 is not used to infer a separate active-flow or bypass fraction for each
experiment. The mass-flow multiplier, hydraulic resistance and DP1 offset
must not all be left unconstrained because they can compensate for one
another, especially in the laminar regime.

### 9.8 Bubble-meter exclusion

The bubble flow-meter tests used to check the manufacturer-factor-adjusted MFC
values were performed **without the receiver**. The bubble meter also imposed
its own pressure restriction. Its measured transfer curve characterizes that
test apparatus, not the installed receiver path, and is excluded from the
receiver-flow/DP1 calibration.

### 9.9 Verification and present acceptance boundary

V9 hydraulic results remain provisional until all of the following are
completed:

1. T8 observation and initialization tests pass at 5 mm.
2. Standard-flow conversion reproduces
   $\rho_{\mathrm{ref}}\simeq1.200\ \mathrm{kg\,m^{-3}}$.
3. Mass-flow conservation remains exact while local actual velocity varies
   with $\rho(T,p)$.
4. Square-channel pressure loss has the correct units, sign, flow scaling and
   temperature response.
5. The cold-$t_0$ DP1 relation is reproduced within its experimental scatter.
6. A common hydraulic closure is validated on full cooling histories and then
   on heating histories.
7. Thermal energy conservation, mesh sensitivity, observation mapping and
   parameter-identifiability tests continue to pass.
8. Only after hydraulic validation are the apparent Nusselt parameters
   recalibrated.

Items 1--5 and the hydraulic part of item 7 were completed on 2026-07-29:

```text
test/smoke_2D_v9.jl          34/34 passed
test/check_2D_v9_physics.jl  46/46 passed
test/check_2D_v9_mesh.jl       2/2 passed
test/smoke_2D_v8.jl          34/34 passed
```

The full 15-heating/3-cooling audit is written to
`summaries/2D_v9/dp1_summary_2D_v9.csv`. With the corrected MFC mass flow fixed
at `mass_flow_scale = 1`, its endpoint statistics are:

| Comparison | Bias, model-data (mbar) | RMSE (mbar) |
| :--- | ---: | ---: |
| 15 heating starts | -0.0159 | 0.0392 |
| 15 heating final points | -0.4507 | 0.5233 |
| 3 cooling initial points | -0.4218 | 0.4811 |
| 3 cooling final points | -0.0336 | 0.0397 |

The mean local channel velocities span 0.3404--1.4701 m/s initially and
0.3395--2.3111 m/s at the final saved states. These are actual local
velocities obtained from conserved standard-flow mass and $\rho(T)$, not
standard volumetric flow divided by channel area.

The cold and cooled endpoints validate the reference-flow conversion and
linear cold resistance. The systematic hot underprediction means item 6 is
not complete. A temperature/velocity-dependent common path loss, such as one
global contraction/expansion coefficient, may be tested using complete
histories after the v9 thermal state is calibrated. It must be calibrated on
a subset and validated on held-out runs. An independent flow or bypass
fraction per experiment remains prohibited, and the global flow scale,
hydraulic multiplier, DP1 offset and minor-loss coefficient must not all be
free simultaneously.

The initial implementation assessment above is superseded by the full
train/validation test in Section 9.10.

### 9.10 Full v9 fitted test and acceptance boundary

#### Calibration/validation separation

Thermal and optical parameters are calibrated on:

```text
E67, E69, E71, E72, E74, E76, E77, E79, E81
```

and evaluated without refitting on:

```text
heating: E68, E70, E73, E75, E78, E80
cooling: C69, C80, C81
```

Every irradiance group supplies three training and two validation heating
cases. The MFC mass-flow multiplier remains fixed at 1.0, and DP1 is excluded
from the thermal objective.

The 120-evaluation derivative-free fit reduced the objective from 14.4361 to
10.3625 but returned `MaxIters`. Its vector is

| Parameter | Fitted value | Status |
| :--- | ---: | :--- |
| $A_{\mathrm{Nu}}$ | 0.00123739 | fitted |
| $B_{\mathrm{Re}}$ | 1.44 | independently fixed |
| $\chi_r$ | 0.180685 | fitted |
| $\chi_z$ | 0.200000 | active lower bound |
| $s_{456}$ | 1.050672 | fitted |
| $s_{304}$ | 1.051748 | fitted |
| $s_{256}$ | 0.621298 | fitted |
| $\beta_{\mathrm{rad}}$ | 102.894 1/m | fitted |
| $\beta_{\mathrm{opt}}$ | 21.2914 1/m | fitted |

#### Incremental hot-path pressure term

The cold empirical multiplier already incorporates the reference-state path
loss. The fitted hot term must therefore vanish at the reference state to
avoid double counting:

$$
DP1_{\mathrm{v9}} =
p_{0,\mathrm{DP1}}+
C_{\mathrm{hyd}}\Delta p_{\mathrm{square}}+
K_{\mathrm{hot}}\left(
\overline{q_{\mathrm{dyn}}}-
\overline{q_{\mathrm{dyn,ref}}}
\right),
$$

$$
q_{\mathrm{dyn}}=\frac{\rho(T_g)u(T_g)^2}{2},
\qquad
q_{\mathrm{dyn,ref}}=
\frac{\rho_{\mathrm{ref}}u_{\mathrm{ref}}^2}{2}.
$$

The ring overbar denotes the same mass-flow-weighted aggregate used for
receiver pressure. Fitting only $K_{\mathrm{hot}}$ to the nine training DP1
histories gives

$$
K_{\mathrm{hot}}=118.479.
$$

This is a lumped coefficient referenced to channel-end dynamic pressure. Its
large magnitude is not the minor-loss coefficient of a known elbow,
contraction or fitting. It includes the unknown tap-to-atmosphere path,
developing flow and thermal-state discrepancy.

| Phase | Base DP1 RMSE (mbar) | With $K_{\mathrm{hot}}$ (mbar) | Bias with $K_{\mathrm{hot}}$ (mbar) |
| :--- | ---: | ---: | ---: |
| Heating training | 0.4398 | 0.1708 | +0.0044 |
| Held-out heating | 0.4106 | 0.1932 | +0.0113 |
| Cooling validation | 0.1648 | 0.0958 | -0.0840 |

The transfer to unseen heating and cooling data supports a common hot-path
correction. It does not identify a bypass or justify changing the MFC mass
flow.

#### Thermal, spatial and numerical results

| Phase | Mean sensor RMSE (K) | Steady MAE (K) | t90 MAE (s) |
| :--- | ---: | ---: | ---: |
| Heating training | 76.21 | 78.72 | 700.65 |
| Held-out heating | 63.92 | 68.46 | 746.43 |
| Cooling validation | 34.71 | 15.75 | 535.71 |

The model reproduces 20/21 flow-slope signs, but it predicts the wrong
core/perimeter sign at both 58 and 107 mm for all 15 heating experiments.
Full-transient sensitivity has rank 7/8 at a relative cutoff of $10^{-3}$,
condition number 21461, and $\chi_z$ is bound-active. Consequently, the
thermal coefficients are not identifiable and validated.

At the fitted point, the E73 maximum sensor change is 5.876 K from coarse to
nominal mesh and 4.909 K from nominal to fine mesh. The nominal-to-fine DP1
change is 0.00943 mbar, and the instantaneous energy-rate residual is
$1.42\times10^{-14}$ W. These results establish numerical convergence and
conservation, not physical coefficient validity.

Final status: **v9 implementation and temperature-dependent velocity PASS;
cold hydraulic constraint PASS; common hot DP1 prediction PROVISIONAL PASS;
thermal coefficient extraction FAIL**.

The completed post-processing set is under `summaries/2D_v9/plots/` and
contains steady-temperature and DP1 parity figures, one four-panel transient
figure for each of the 18 experiments, 15 heating axial profiles, 15 final 2D
temperature fields and the fitted identifiability-correlation heat map.
