# Mathematical and Theoretical Manual: 2D Axisymmetric Continuum Macro-ECM Model (`2D_v2`) for Monolithic Solar Receivers

---

## Executive Summary & Model Purpose

This manual presents the complete mathematical formulation, physical transport phenomena, numerical discretization, boundary conditions, and calibration methodology for the **2D Axisymmetric Continuum Macro-ECM Solar Receiver Model (`2D_v2`)**.

The model represents a structured porous silicon carbide (SiC) monolithic receiver disc ($19 \times 19\text{ mm}^2$ frontal square face, $137\text{ mm}$ length, 100 square channels, open porosity $\varepsilon = 0.65$, porous solid mass $40.0\text{ g}$) enclosed within an alumina felt insulation sleeve ($r \in [33.9, 57.0\text{ mm}]$), an aluminum outer casing ($r \in [57.0, 75.0\text{ mm}]$), a rear metallic adaptor sleeve ($z \in [109, 137\text{ mm}]$), a downstream alumina gas exit tube ($z \in [137, 287\text{ mm}]$), and a water-cooled flange sink ($25.0^\circ\text{C}$).

The primary objective of the 2D continuum model (Entire Converter Model / Macro-ECM) is to extract effective macroscopic transport coefficients—including anisotropic effective thermal conductivities $k_{s,r}^{\text{eff}}$ and $k_{s,z}^{\text{eff}}$, developing-flow Nusselt correlations $\text{Nu}(\text{Re}, \text{Pr}, z)$, Rosseland radiation extinction coefficients $\beta_{\text{opt}}$, front buoyant plume natural convection $h_{\text{front}}(T)$, and temperature-dependent insulation properties $k_{\text{felt}}(T)$—from multi-case transient experimental data.

---

## 1. Domain Geometry & Multi-Domain Radial Discretization

The physical receiver assembly is modeled in cylindrical axisymmetric coordinates $(r, z)$, where $r \in [0, R_{\text{casing}}]$ is the radial coordinate and $z \in [0, L_{\text{rec}} + L_{\text{rear}}]$ is the axial depth coordinate.

```text
  r = 0 mm            r = 33.9 mm              r = 57.0 mm          r = 75.0 mm
  ┌───────────────────────┬────────────────────────┬─────────────────────┐
  │  SiC Monolith Core    │ Alumina Felt           │ Aluminum Outer      │
  │  100 Square Channels  │ Insulation Layer       │ Metallic Casing     │
  │  (ρ = 2150 kg/m³)     │ (ρ = 140 kg/m³)        │ (ρ = 2700 kg/m³)    │
  └───────────────────────┴────────────────────────┴─────────────────────┘
```

### 1.1 Multi-Domain Radial Sub-Regions

1. **SiC Monolithic Channel Core Matrix ($r \in [0, R_{\text{rec}}]$)**:
   - Equivalent disc radius: $R_{\text{rec}} = \sqrt{\frac{A_{\text{frt}}}{\pi}} = \sqrt{\frac{19.0 \times 19.0 \times 10^{-6}}{\pi}} \approx 33.85 \times 10^{-3}\text{ m} \quad (33.9\text{ mm})$.
   - Frontal Area: $A_{\text{frt}} = 3.61 \times 10^{-3}\text{ m}^2$.
   - Channel Geometry: 100 square channels, cell width $w_{\text{ch}} = 1.5\text{ mm}$ (hydraulic diameter $D_h = 1.5\text{ mm}$), web wall thickness $t_w = 0.4\text{ mm}$.
   - Open Flow Porosity: $\varepsilon = \frac{100 \times (1.5\times 10^{-3})^2}{3.61 \times 10^{-3}} \approx 0.623$.
   - Solid Volume Fraction: $1 - \varepsilon \approx 0.377$.
   - Solid Material: Recrystallized SiC, dense solid density $\rho_s = 2150\text{ kg/m}^3$. Measured total porous receiver mass $m_{\text{rec}} = 40.0\text{ g}$.
   - Discretization: $N_{r,\text{rec}} = 10$ radial concentric rings.

2. **Alumina Felt Insulation Ring ($r \in [R_{\text{rec}}, R_{\text{felt}}]$)**:
   - Outer Radius: $R_{\text{felt}} = 57.0\text{ mm}$. Thickness $\Delta r_{\text{felt}} = 23.1\text{ mm}$.
   - Material Properties: Alumina insulation felt, density $\rho_{\text{felt}} = 140\text{ kg/m}^3$, heat capacity $C_{p,\text{felt}} = 1360\text{ J/kg/K}$.
   - Temperature-Dependent Thermal Conductivity: $k_{\text{felt}}(T) = k_0 + \gamma_{\text{rad}} T^3$ (accounting for solid matrix conduction and radiative pore diffusion).
   - Discretization: $N_{r,\text{felt}} = 5$ radial concentric rings.

3. **Aluminum Outer Metallic Casing ($r \in [R_{\text{felt}}, R_{\text{casing}}]$)**:
   - Outer Radius: $R_{\text{casing}} = 75.0\text{ mm}$ ($150\text{ mm}$ outer diameter). Casing wall thickness $\Delta r_{\text{case}} = 18.0\text{ mm}$.
   - Material Properties: Aluminum alloy casing, density $\rho_{\text{Al}} = 2700\text{ kg/m}^3$, conductivity $k_{\text{Al}} = 205\text{ W/m/K}$, heat capacity $C_{p,\text{Al}} = 900\text{ J/kg/K}$.
   - Discretization: $N_{r,\text{case}} = 2$ radial concentric rings.

Total Radial Nodes: $N_r = N_{r,\text{rec}} + N_{r,\text{felt}} + N_{r,\text{case}} = 10 + 5 + 2 = 17$ radial nodes.

---

### 1.2 Axial Sub-Regions & Back-End Assembly

1. **Monolithic Receiver Disc ($z \in [0, L_{\text{rec}}]$)**:
   - Receiver Length: $L_{\text{rec}} = 137\text{ mm}$.
   - Discretization: $N_z = 25$ axial finite-volume cells ($\Delta z = 5.48\text{ mm}$).

2. **Adaptor Sleeve Overlap Region ($z \in [109, 137\text{ mm}]$)**:
   - Adaptor Diameter: $77.6\text{ mm}$ ($r_{\text{adaptor}} = 38.8\text{ mm}$). Overlap length $L_{\text{overlap}} = 28.0\text{ mm}$.
   - Contact Conductance: $G_{\text{adaptor}} = \frac{1}{R_{\text{contact}}}$ couples the rear perimeter nodes to the alumina exit tube wall.

3. **Downstream Alumina Gas Exit Tube ($z \in [L_{\text{rec}}, L_{\text{rec}} + L_{\text{rear}}]$)**:
   - Tube Length: $L_{\text{rear}} = 150.0\text{ mm}$ ($z = 137 \to 287\text{ mm}$).
   - Tube Dimensions: Inner gas radius $r_{\text{inner}} = 6.5\text{ mm}$ ($13.0\text{ mm}$ ID), wall thickness $t_{\text{wall}} = 1.5\text{ mm}$ ($r_{\text{outer}} = 8.0\text{ mm}$).
   - Wall Material Properties: High-purity alumina ($\rho_{\text{alumina}} = 3900\text{ kg/m}^3$, $k_{\text{alumina}}(T) = 5.5 + 34.5 e^{-0.0033(T-273.15)}\text{ W/m/K}$).
   - Discretization: $N_{z,\text{rear}} = 20$ axial finite-volume cells ($\Delta z_{\text{rear}} = 7.50\text{ mm}$).
   - Gas Exit Probe $T_3$: Sampled inside the gas stream at $z = 140\text{ mm}$ ($3\text{ mm}$ behind receiver exit face).

4. **Water-Cooled Flange Thermal Sink**:
   - Temperature: Fixed at $T_{\text{flange}} = 25.0^\circ\text{C} = 298.15\text{ K}$.
   - Conductive Loss: Conduction down the alumina tube wall and aluminum casing into the water-cooled flange.

---

## 2. Governing Continuum Energy Conservation Equations

The model solves Local Thermal Non-Equilibrium (LTNE) energy equations coupling the 2D solid/insulation continuum $T_s(r, z, t)$ to quasi-steady 2D channel gas flow $T_g(r, z, t)$.

```text
  r = 0 mm            r = 33.9 mm              r = 57.0 mm          r = 75.0 mm
  ┌───────────────────────┬────────────────────────┬─────────────────────┐
  │  SiC Monolith Core    │ Alumina Felt           │ Aluminum Outer      │
  │  100 Square Channels  │ Insulation Layer       │ Metallic Casing     │
  │  (ρ = 2150 kg/m³)     │ (ρ = 140 kg/m³)        │ (ρ = 2700 kg/m³)    │
  └───────────────────────┴────────────────────────┴─────────────────────┘
```

### 2.1 Solid Continuum Energy Equation ($T_s(r, z, t)$)

For any cell $(i, j)$ in the 2D domain:

$$\rho_{s,i} \phi_{s,i} C_{p,s}(T_s) \frac{\partial T_s}{\partial t} = \frac{1}{r} \frac{\partial}{\partial r} \left( r k_{s,r}^{\text{eff}}(T_s) \frac{\partial T_s}{\partial r} \right) + \frac{\partial}{\partial z} \left( k_{s,z}^{\text{eff}}(T_s) \frac{\partial T_s}{\partial z} \right) + q_{\text{solar}}^{\prime\prime\prime}(r, z, t) - q_{\text{conv}}^{\prime\prime\prime}(r, z, t)$$

where:
- $\phi_{s,i} = 1 - \varepsilon$ for $i \le N_{r,\text{rec}}$ (porous monolith matrix), and $\phi_{s,i} = 1.0$ for $i > N_{r,\text{rec}}$ (solid insulation and casing).
- $q_{\text{conv}}^{\prime\prime\prime}(r, z, t) = a_v h_c(r, z) (T_s(r, z, t) - T_g(r, z, t))$ is the volumetric gas-solid convective heat exchange.
- Volumetric Interfacial Area Density: $a_v = \frac{4 \cdot N_{\text{ch}} \cdot w_{\text{ch}}}{A_{\text{frt}}} = \frac{4 \times 100 \times 1.5\times 10^{-3}}{3.61\times 10^{-3}} \approx 166.2\text{ m}^2/\text{m}^3$.

---

### 2.2 Anisotropic Effective Thermal Conductivities

1. **Effective Radial Thermal Conductivity ($k_{s,r}^{\text{eff}}$)**:
   In a porous square-channel monolith matrix, radial heat transport must cross web walls and gas void channels, yielding a significantly lower effective thermal conductivity than dense SiC:
   $$k_{s,r}^{\text{eff}}(T_s) = \chi_r \cdot k_{\text{dense, SiC}}(T_s)$$
   where $\chi_r \approx 0.005 - 0.033$ is the calibrated radial scaling factor ($k_{s,r}^{\text{eff}} \approx 0.5 - 3.3\text{ W/m/K}$).

   Dense SiC Thermal Conductivity correlation:
   $$k_{\text{dense, SiC}}(T) = \max\left(5.0, \frac{1}{1.10\times 10^{-4} + 3.42\times 10^{-6} (T - 273.15)}\right) \quad [\text{W/m}\cdot\text{K}]$$

2. **Effective Axial Thermal Conductivity ($k_{s,z}^{\text{eff}}$)**:
   Axial thermal transport combines solid wall conduction and Rosseland solid-phase radiative diffusion:
   $$k_{s,z}^{\text{eff}}(T_s) = \chi_z \cdot k_{\text{dense, SiC}}(T_s) + \frac{16 n^2 \sigma T_s^3}{3 \beta_{\text{opt}}}$$
   where $\sigma = 5.67037\times 10^{-8}\text{ W/m}^2\text{K}^4$, $n \approx 1.0$ is the refractive index, and $\beta_{\text{opt}} \approx 184.67\text{ m}^{-1}$ is the optical extinction coefficient.

3. **Temperature-Dependent Insulation Felt Conductivity ($k_{\text{felt}}(T)$)**:
   For $R_{\text{rec}} < r \le R_{\text{felt}}$, heat transfer through alumina felt insulation combines matrix conduction and high-temperature radiative pore transport:
   $$k_{\text{felt}}(T_s) = \max\left(0.05, k_{\text{felt,0}} + 1.2\times 10^{-10} T_s^3\right) \quad [\text{W/m}\cdot\text{K}]$$
   - At $T = 300\text{ K}$: $k_{\text{felt}} \approx 0.063\text{ W/m/K}$.
   - At $T = 1000\text{ K}$: $k_{\text{felt}} \approx 0.180\text{ W/m/K}$.
   - At $T = 1300\text{ K}$: $k_{\text{felt}} \approx 0.323\text{ W/m/K}$.

---

## 3. Quasi-Steady 2D Channel Gas Flow & Energy Conservation

Because the gas mass residence time in the $137\text{ mm}$ channels ($\tau_{\text{res}} \approx \frac{L}{u_g} \sim 0.01\text{ s}$) is orders of magnitude smaller than the solid thermal response time ($\tau_{\text{solid}} \sim 50 - 200\text{ s}$), gas flow is governed by quasi-steady energy marching.

### 3.1 Ring-Wise Gas Energy Balance

For each radial channel ring $i \in [1, N_{r,\text{rec}}]$:

$$\dot{m}_i C_{p,g}(T_g) \frac{\partial T_g}{\partial z} = h_{c,i}(z) P_{\text{ex},i} (T_s(r_i, z) - T_g(r_i, z))$$

where:
- $P_{\text{ex},i} = A_{\text{frt},i} \cdot a_v$ is the convective heat exchange perimeter for ring $i$.
- $\dot{m}_i$ is the mass flow rate assigned to ring $i$.

---

### 3.2 Exact Cell-by-Cell NTU Effectiveness Integration

Across axial cell $j \in [1, N_z]$ of length $\Delta z$:

$$\text{UA}_{i,j} = h_{c,i,j} \cdot P_{\text{ex},i} \cdot \Delta z$$

$$\text{NTU}_{i,j} = \frac{\text{UA}_{i,j}}{\dot{m}_i C_{p,g}(T_{\text{film}})}$$

$$\varepsilon_{i,j} = 1 - \exp(-\text{NTU}_{i,j})$$

$$T_{g,i,j+1} = T_{g,i,j} + \varepsilon_{i,j} (T_{s,i,j} - T_{g,i,j})$$

$$Q_{\text{gas},i,j} = \dot{m}_i C_{p,g}(T_{\text{film}}) (T_{g,i,j+1} - T_{g,i,j})$$

where film temperature $T_{\text{film}} = \frac{1}{2} (T_{s,i,j} + T_{g,i,j})$.

---

### 3.3 Radial Preferential Flow Distribution & Active Recruitment

1. **Re-Dependent Active Flow Recruitment ($\phi_{\text{act}}(\text{Re})$)**:
   At low flow rates ($\text{Re} \sim 20 - 50$), flow recruitment is partial, scaling active flow advection:
   $$\phi_{\text{act}} = \text{clamp}\left( \phi_0 \left( \frac{\text{Re}}{\text{Re}_0} \right)^m, 0.20, 1.00 \right)$$
   $$\dot{m}_{\text{active}} = \phi_{\text{act}} \cdot \dot{m}_{\text{total}}$$

2. **Preferential Core Channel Flow Distribution**:
   To capture buoyancy and manifold acceleration that direct higher airflow to central channels:
   $$\psi(r_i) = 1.0 - c_{\text{radial\_flow}} \left( \frac{r_i}{R_{\text{rec}}} \right)^2$$
   $$\dot{m}_i = \dot{m}_{\text{active}} \cdot \frac{\psi(r_i) A_{\text{flow},i}}{\sum_{k=1}^{N_{r,\text{rec}}} \psi(r_k) A_{\text{flow},k}}$$

---

### 3.4 Convective Heat Transfer Correlations

Developing laminar flow in square channels is governed by the Sieder-Tate style Nusselt correlation with entry length amplification:

$$\text{Re}_i = \frac{\dot{m}_i D_h}{A_{\text{flow},i} \cdot \mu_g(T_{\text{film}})}$$

$$\text{Pr} = \frac{C_{p,g} \cdot \mu_g}{k_g}$$

$$\text{Nu}_{\text{dev}}(z) = A_{\text{Nu}} \cdot \text{Re}^{B_{\text{Re}}} \cdot \text{Pr}^{1/3} \cdot \left( \frac{D_h}{\max(z, D_h/2)} \right)^{1/3}$$

$$\text{Nu}_{i,j} = \max\left( \text{Nu}_{\text{floor}}, \text{Nu}_{\text{dev}}(z_j) \right) \quad (\text{Nu}_{\text{floor}} = 3.61)$$

$$h_{c,i,j} = \frac{\text{Nu}_{i,j} \cdot k_g(T_{\text{film}})}{D_h}$$

Air thermophysical properties ($T$ in Kelvin):
- $k_g(T) = \frac{2.414\times 10^{-3} \sqrt{T}}{1 + 245.4 \cdot 10^{-12/T} / T} \quad [\text{W/m}\cdot\text{K}]$
- $C_{p,g}(T) = (1 + 1.983\times 10^{-4} T - 4.14\times 10^{-8} T^2) \times 1004.0 \quad [\text{J/kg}\cdot\text{K}]$
- $\mu_g(T) = \frac{1.458\times 10^{-6} T^{1.5}}{T + 110.4} \quad [\text{Pa}\cdot\text{s}]$

---

## 4. Optical Deposition & Gaussian Beam Distribution

Solar energy absorption combines a Gaussian radial beam distribution, a front surface specular absorption fraction, and Beer-Lambert exponential depth attenuation.

```text
  Solar Flux In (Q_solar)               Front Plume Loss (Q_front)
       │                                     ▲
       ▼                                     │
┌───────────────┐     Convection h(Re)  ┌───────────────┐     Housing Loss
│ Front Face    ├──────────────────────►│ Air Flow Stream├───────────────────► Flange Sink
│ SiC Monolith  │     Solid Conduction  │ & Exit Tube T3│     k_felt(T) Radiative
└───────────────┘                       └───────────────┘
```

### 4.1 Delivered Power Accounting

$$Q_{\text{absorbed}} = \eta_{\text{abs}} \cdot (f_{\text{scale}} \cdot I_{\text{reported}}) \cdot A_{\text{frt}}$$

where $\eta_{\text{abs}} = 0.85$, $I_{\text{reported}}$ is the nominal aperture-average irradiance ($456000, 304000, 256000\text{ W/m}^2$), and $f_{\text{scale}} \in \{f_{456}, f_{304}, f_{256}\}$ is the power scaling factor.

Power partition:
- Core Monolith Absorbed Power: $Q_{\text{core}} = (1 - f_{\text{spill}}) Q_{\text{absorbed}}$.
- Rim Spillage Absorbed Power: $Q_{\text{spill}} = f_{\text{spill}} Q_{\text{absorbed}}$ (deposited into insulation ring $r \in [R_{\text{rec}}, R_{\text{felt}}]$).

---

### 4.2 Spatial 2D Deposition Weights

1. **Gaussian Radial Weighting ($w_r(r)$)**:
   $$I_{\text{rad}}(r_i) = \exp\left( -\frac{r_i^2}{2 \sigma_{\text{beam}}^2} \right)$$
   $$w_{r,i} = \frac{I_{\text{rad}}(r_i) \cdot A_{\text{frt},i}}{\sum_{k=1}^{N_{r,\text{rec}}} I_{\text{rad}}(r_k) A_{\text{frt},k}}$$

2. **Axial Depth Weighting ($w_z(z)$)**:
   $$w_{z,\text{depth},j} = \frac{\exp(-\beta_{\text{opt}} z_{j-1}) - \exp(-\beta_{\text{opt}} z_j)}{1 - \exp(-\beta_{\text{opt}} L)}$$
   $$w_{z,j} = (1 - f_{\text{front}}) w_{z,\text{depth},j} + f_{\text{front}} \cdot \delta_{j,1}$$

3. **Combined 2D Volumetric Heat Source ($Q_{\text{solar},i,j}$)**:
   $$Q_{\text{solar},i,j} = Q_{\text{core}} \cdot w_{r,i} \cdot w_{z,j}$$

---

## 5. Boundary Conditions & Heat Loss Models

```text
  Solar Flux In (Q_solar)               Front Plume Loss (Q_front)
       │                                     ▲
       ▼                                     │
┌───────────────┐     Convection h(Re)  ┌───────────────┐     Housing Loss
│ Front Face    ├──────────────────────►│ Air Flow Stream├───────────────────► Flange Sink
│ SiC Monolith  │     Solid Conduction  │ & Exit Tube T3│     k_felt(T) Radiative
└───────────────┘                       └───────────────┘
```

### 5.1 Front Face Boundary ($z = 0$)

At the irradiated front face ($z = 0$), thermal losses include radiative re-emission to ambient and **Churchill-Chu turbulent/laminar buoyant plume natural convection**:

$$Q_{\text{front},i} = A_{\text{solid},i} \left[ \varepsilon_{\text{front}} \sigma (T_s(r_i, 0)^4 - T_{\text{amb}}^4) + h_{\text{front}}(T_s(r_i, 0), T_{\text{amb}}) (T_s(r_i, 0) - T_{\text{amb}}) \right]$$

#### Churchill-Chu Natural Convection Plume Correlation
$$\text{Ra}_L = \frac{g \beta (T_s - T_{\text{amb}}) L_{\text{char}}^3}{\nu \cdot \alpha}$$

$$\text{Nu}_{\text{front}} = \left\{ 0.825 + \frac{0.387 \text{Ra}_L^{1/6}}{\left[ 1 + (0.492 / \text{Pr})^{9/16} \right]^{8/27}} \right\}^2$$

$$h_{\text{front}}(T_s, T_{\text{amb}}) = \max\left( 10.0, \frac{\text{Nu}_{\text{front}} \cdot k_{\text{air}}}{L_{\text{char}}} \right)$$

where $L_{\text{char}} = 2 R_{\text{rec}} = 0.0678\text{ m}$.

---

### 5.2 Outer Radial Housing Boundary ($r = R_{\text{casing}}$)

Heat loss from the outer aluminum casing surface ($r = R_{\text{casing}} = 75.0\text{ mm}$) to ambient room:

$$Q_{\text{casing},j} = (2 \pi R_{\text{casing}} \Delta z) \left[ h_{\text{nat,casing}} (T_s(R_{\text{casing}}, z_j) - T_{\text{amb}}) + \varepsilon_{\text{casing}} \sigma (T_s(R_{\text{casing}}, z_j)^4 - T_{\text{amb}}^4) \right]$$

where $h_{\text{nat,casing}} = 10.0\text{ W/m}^2\text{K}$ and $\varepsilon_{\text{casing}} = 0.20$.

---

### 5.3 Rear Exit Boundary & Flange Sink ($z = L_{\text{rec}} \to L_{\text{rec}} + L_{\text{rear}}$)

1. **Adaptor Sleeve Overlap Conductance**:
   At $z = L_{\text{rec}}$, heat is conducted from receiver rear cells into the alumina exit tube wall:
   $$Q_{\text{rear,exit},i} = G_{\text{adaptor}} \left( \frac{A_{\text{solid},i}}{\sum A_{\text{solid}}} \right) (T_s(r_i, L_{\text{rec}}) - T_{\text{rear},1})$$

2. **Downstream Alumina Exit Tube Wall Conduction**:
   Along the exit tube wall $j \in [1, N_{z,\text{rear}}]$:
   $$Q_{\text{cond,rear},j} = k_{\text{alumina}}(T_{\text{rear},j}) A_{\text{wall,tube}} \frac{T_{\text{rear},j} - T_{\text{rear},j+1}}{\Delta z_{\text{rear}}}$$

3. **Water-Cooled Flange Thermal Sink ($T_{\text{flange}} = 298.15\text{ K}$)**:
   Direct wall conduction to the $25.0^\circ\text{C}$ water-cooled flange assembly:
   $$Q_{\text{flange},j} = s_{\text{flange}} \left( \frac{k_{\text{alumina}} A_{\text{wall,tube}}}{\Delta z_{\text{rear}}} \right) (T_{\text{rear},j} - T_{\text{flange}})$$

4. **Downstream Exit Gas Marching**:
   Gas leaving the monolith channels mixes into the rear tube:
   $$T_{g,\text{rear},1} = \frac{\sum_i \dot{m}_i C_{p,g} T_{g,i,N_z+1}}{\sum_i \dot{m}_i C_{p,g}}$$
   Gas marches through the alumina tube exchanging heat with $T_{\text{rear},j}$ before reaching gas exit thermocouple $T_3$ at $z = 140\text{ mm}$.

---

## 6. Numerical Discretization & SciML DifferentialEquations.jl Implementation

### 6.1 State Vector Architecture

The complete system state vector $u(t) \in \mathbb{R}^{N_{\text{states}}}$ has length:
$$N_{\text{states}} = N_r \cdot N_z + N_{z,\text{rear}} + 1$$

With $N_r = 17$, $N_z = 25$, $N_{z,\text{rear}} = 20$:
$$N_{\text{states}} = 17 \times 25 + 20 + 1 = 446 \text{ dynamic states}$$

State allocation:
1. $u[1 : N_r N_z]$: 2D solid/felt/casing temperatures $T_s(r_i, z_j, t)$.
2. $u[N_r N_z + 1 : N_r N_z + N_{z,\text{rear}}]$: Rear alumina tube wall temperatures $T_{\text{rear}}(z_j, t)$.
3. $u[N_r N_z + N_{z,\text{rear}} + 1]$: Participating housing cavity temperature $T_{\text{cavity}}(t)$.

---

### 6.2 Ordinary Differential Equation (ODE) Problem Formulation

The ODE system is expressed in SciML `OrdinaryDiffEq.jl` standard form:

$$\frac{d u}{d t} = f(u, p, t)$$

Solved using the stiff SDIRK/Rosenbrock ODE solver `Rodas5P(autodiff=AutoFiniteDiff())` with adaptive time stepping (`reltol=1e-6`, `abstol=1e-7`, `dtmax=30.0`).

---

## 7. Thermocouple Probe Spatial Mapping

To evaluate model predictions against experimental thermocouple measurements, simulated 2D temperature fields are mapped to exact physical sensor positions:

| Thermocouple Sensor | Physical Location & Domain | Model Grid Coordinate Mapping |
| :--- | :--- | :--- |
| **$T_8$** | Front Perimeter ($z = 11\text{ mm}, r = R_{\text{rec}}$) | $T_s(r = R_{\text{rec}}, z = 11\text{ mm})$ |
| **$T_{12}$** | Mid-Depth Perimeter ($z = 58\text{ mm}, r = R_{\text{rec}}$) | $T_s(r = R_{\text{rec}}, z = 58\text{ mm})$ |
| **$T_{11}$** | Deep Perimeter ($z = 107\text{ mm}, r = R_{\text{rec}}$) | $T_s(r = R_{\text{rec}}, z = 107\text{ mm})$ |
| **$T_9$** | Mid-Depth Core ($z = 58\text{ mm}, r = 0$) | $T_s(r = 0, z = 58\text{ mm})$ |
| **$T_{10}$** | Deep Core ($z = 107\text{ mm}, r = 0$) | $T_s(r = 0, z = 107\text{ mm})$ |
| **$T_3$** | Exit Gas Stream ($z = 140\text{ mm}$) | $T_{g,\text{rear}}(z = 140\text{ mm})$ |
| **$T_2$** | Housing Insulation ($z = 58\text{ mm}, r = R_{\text{rec}} + 40\text{ mm}$) | $T_s(r = 73.9\text{ mm}, z = 58\text{ mm})$ |

---

## 8. Parameter Calibration Workflow & Validation Metrics

### 8.1 Multi-Objective Calibration Loss Function

The 15-parameter optimization vector $\theta$ is calibrated using `Optimization.jl` and `OptimizationNLopt.jl` (`NLopt.LN_BOBYQA` algorithm) across heating and cooling runs:

$$L(\theta) = \frac{1}{N_{\text{cases}}} \sum_{\text{case}} \sum_{j=1}^7 \left[ \frac{\text{MSE}_j}{\text{scale}_j^2} + w_{\text{shape}} \text{ShapeLoss}_j + w_{\text{t90}} \text{TimeLoss}_j \right] + w_{\text{space}} \text{SpatialOrderPenalty}$$

where:
- $\text{MSE}_j = \frac{1}{N_t} \sum_k (T_{\text{model},j}(t_k) - T_{\text{exp},j}(t_k))^2$.
- $\text{ShapeLoss}_j = \frac{\frac{1}{N_t-1} \sum_k (\Delta T_{\text{model}} - \Delta T_{\text{exp}})^2}{\text{span}_j^2}$.
- $\text{TimeLoss}_j = \left( \frac{t_{90,\text{model}} - t_{90,\text{exp}}}{t_{\text{total}}} \right)^2$.
- $\text{SpatialOrderPenalty} = \frac{1}{100^2} \left[ \max(0, \Delta T_{12-9,\text{model}} - \Delta T_{12-9,\text{exp}})^2 + \max(0, \Delta T_{11-10,\text{model}} - \Delta T_{11-10,\text{exp}})^2 \right]$.

---

### 8.2 Parameter Optimization Bounds ($\text{LB}_{2D} \le \theta \le \text{UB}_{2D}$)

| Parameter Index | Symbol | Description | Lower Bound ($\text{LB}$) | Upper Bound ($\text{UB}$) | Calibrated Value (`2D_v2`) |
| :---: | :--- | :--- | :---: | :---: | :---: |
| $\theta_1$ | $A_{\text{Nu}}$ | Laminar developing Nu prefactor | $0.50$ | $6.00$ | **2.500** |
| $\theta_2$ | $B_{\text{Re}}$ | Reynolds exponent | $0.10$ | $0.60$ | **0.333** |
| $\theta_3$ | $\chi_r$ | Effective radial conductivity scale | $0.001$ | $0.100$ | **0.030** |
| $\theta_4$ | $\chi_z$ | Effective axial conductivity scale | $0.10$ | $3.00$ | **1.000** |
| $\theta_5$ | $\sigma_{\text{beam}}$ | Gaussian beam waist [mm] | $5.0$ | $35.0$ | **14.0 mm** |
| $\theta_6$ | $f_{\text{spill}}$ | Rim spillage power fraction | $0.01$ | $0.40$ | **0.200** |
| $\theta_7$ | $C_{\text{cavity}}$ | Assembly thermal mass [J/K] | $100.0$ | $500.0$ | **301.0 J/K** |
| $\theta_8$ | $f_{456}$ | Power scale factor at 456 kW/m2 | $0.50$ | $2.80$ | **1.336** |
| $\theta_9$ | $f_{304}$ | Power scale factor at 304 kW/m2 | $0.50$ | $2.80$ | **1.374** |
| $\theta_{10}$ | $f_{256}$ | Power scale factor at 256 kW/m2 | $0.40$ | $2.00$ | **0.786** |
| $\theta_{11}$ | $\phi_0$ | Active flow fraction at Re=50 | $0.40$ | $1.00$ | **1.000** |
| $\theta_{12}$ | $m_{\text{rec}}$ | Flow recruitment exponent | $0.00$ | $0.60$ | **0.000** |
| $\theta_{13}$ | $f_{\text{front}}$ | Front surface specular absorption | $0.00$ | $0.50$ | **0.200** |
| $\theta_{14}$ | $s_{\text{flange}}$ | Water flange sink scale factor | $0.10$ | $5.00$ | **1.000** |
| $\theta_{15}$ | $c_{\text{radial\_flow}}$ | Radial flow partition exponent | $0.00$ | $0.50$ | **0.100** |

---

### 8.3 Quantitative Performance Summary Table (`2D_v2`)

```text
Phase/Case        Sensor      2D_v1 Fitted RMSE (K)    2D_v2 Fitted RMSE (K)    2D_v2 Steady Error (K)
Heating (E81)     T8                 171.5                    171.5                    +217.2
Heating (E81)     T12_perim           57.0                     57.0                     +70.8
Heating (E81)     T11_perim           27.5                     25.5                     +25.0
Heating (E81)     T9_core             69.1                     64.6                     +84.1
Heating (E81)     T10_core            46.6                     44.6                     +45.6
Heating (E81)     T3                  30.8                     30.0                     +15.9
Heating (E81)     T2                  53.9                     53.9                     +76.7

Cooling (C81)     T8                 101.9                     98.8                      +1.0
Cooling (C81)     T12_perim           61.3                     58.4                      -1.6
Cooling (C81)     T11_perim           26.7                     25.0                      -5.2
Cooling (C81)     T9_core             60.8                     57.8                      -2.5
Cooling (C81)     T10_core            29.5                     27.8                      -6.2
Cooling (C81)     T3                  21.5                     21.4                      -7.1
Cooling (C81)     T2                  41.3                     42.6                     +10.4
```

---

## Conclusion & Code Structure References

- Model Core: [`2D_v2.jl`](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/2D_v2.jl)
- Calibration & Output Runner: [`run_2D_v2.jl`](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/run_2D_v2.jl)
- Unit Smoke Test Suite: [`test/smoke_2D_v2.jl`](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/test/smoke_2D_v2.jl)
- Model Development Journal: [`summaries/journal.2D.md`](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/journal.2D.md)
- Complete Directory of Calibrated Artifacts: [`summaries/2D_v2/`](file:///d:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/2D_v2/)
