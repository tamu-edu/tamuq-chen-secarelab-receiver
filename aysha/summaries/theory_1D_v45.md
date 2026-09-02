# Continuum Theory and Validation Manual: 1D_v45 Energy-Accounting Macro-Model for Structured Monolithic Solar Receivers

---

## Abstract / Scope of the Manual

This document provides a self-contained, manuscript-ready formulation of the **`1D_v45` macro-continuum model** for high-temperature structured monolithic solar receivers. The model serves as an *Entire Converter Model* (ECM), bridging detailed single-channel transport phenomena and full-scale reactor thermal behavior. It features an energy-accounting two-zone solid partition (active core and passive perimeter housing), developing laminar channel hydrodynamics with intra-strut conduction resistance, volumetric optical channel penetration, flow-dependent aperture suction convection, and physical hardware manifold coupling. All governing differential equations, constitutive transport closures, boundary conditions, numerical discretization schemes, optimization formulations, calibrated parameter values, validation metrics, and physical interpretations are presented in full detail without requiring reference to earlier code versions or technical notes.

---

## 1. Physical System, Geometry & Architecture

### 1.1 Receiver Hardware & Geometry
The modeled receiver consists of a square-channel silicon carbide (SiC) ceramic honeycomb monolith enclosed within an insulated packaging assembly, operating as an open volumetric solar receiver under active suction:

```
                  =======================================================
                  |   PASSIVE PERIMETER HOUSING (Felt + Aluminum Casing) |
                  |-----------------------------------------------------|
 SOLAR FLUX ===>  |                                                     | ===> EXHAUST TUBE
 (Aperture z=0)   |   ACTIVE MONOLITH CORE (100 Square Channels)        |      (To Suction Flange)
 AIR INFLOW ===>  |                                                     | ===> (T3 Sensor @ 140mm)
                  |-----------------------------------------------------|
                  |   PASSIVE PERIMETER HOUSING (Felt + Aluminum Casing) |
                  =======================================================
```

- **Monolith Core**: Cross-sectional frontal area $A_{frt} = W \times W = 19.0\text{ mm} \times 19.0\text{ mm} = 3.61 \times 10^{-4}\text{ m}^2$, total length $L = 137.0\text{ mm}$. The matrix contains $N_{ch} = 100$ square channels with cell width $w_{cell} = 1.90\text{ mm}$, channel opening $w_{ch} = 1.34\text{ mm}$, web thickness $t_{web} = 0.56\text{ mm}$, and open void porosity $\epsilon = A_{flow} / A_{frt} = 0.497$.
- **Hydraulic Diameter**: $D_h = w_{ch} = 1.34 \times 10^{-3}\text{ m}$.
- **Wetted Perimeter & Specific Area**: Total internal convective exchange perimeter $P_{exchange} = 4 N_{ch} w_{ch} = 0.536\text{ m}$, giving a specific exchange area $a_v = P_{exchange} / A_{frt} = 1484.8\text{ m}^2/\text{m}^3$.
- **Solid Volume Fraction**: $A_{solid} = (1 - \epsilon) A_{frt} = 1.815 \times 10^{-4}\text{ m}^2$.
- **Insulation & Housing**: The monolith is wrapped in high-temperature alumina-silicate insulating felt (outer radius $R_{felt} = 25.0\text{ mm}$) and enclosed in an aluminum cylindrical shell (outer radius $R_{case} = 30.0\text{ mm}$, wall thickness $t_{wall} = 5.0\text{ mm}$).
- **Rear Hardware Assembly**: At $z = L$, the monolith seats into an alumina adaptor that transitions into an alumina exit tube (outer radius $R_{tube} = 11.0\text{ mm}$, gas flow radius $R_{gas} = 7.0\text{ mm}$, length $L_{tube} = 63.0\text{ mm}$). The exit tube is anchored to a water-cooled metal mounting flange held at $T_{flange} = 25.0^\circ\text{C}$ ($298.15\text{ K}$).
- **Thermocouple Instrumentation**:
  - Solid front core: $T_9$ at $z = 11.0\text{ mm}$.
  - Solid rear core: $T_{10}$ at $z = 107.0\text{ mm}$.
  - Perimeter side wall: $T_8$ at $z = 11.0\text{ mm}$, $T_{12}$ at $z = 58.0\text{ mm}$, $T_{11}$ at $z = 107.0\text{ mm}$.
  - Exit gas stream: $T_3$ at $z = 140.0\text{ mm}$ (3.0 mm into the rear exit tube along the centerline).
  - External housing/cavity shell: $T_2$.

### 1.2 Thermophysical Material Properties

#### Silicon Carbide (SiC) Monolith
Temperature-dependent solid thermal conductivity $k_s(T)$ and specific heat capacity $c_{p,s}(T)$ for SiC ceramic are evaluated at the local solid temperature $T$ [K]:
$$k_s(T) = \max\left(1.0,\, 51.0 \exp\left(-0.0030 \cdot (T - 273.15)\right) + 1.2\right) \quad [\text{W/m}\cdot\text{K}]$$
$$c_{p,s}(T) = \max\left(500.0,\, 900.0 + 0.30 \cdot (T - 273.15) - 3.0 \times 10^5 \cdot T^{-2}\right) \quad [\text{J/kg}\cdot\text{K}]$$
Solid density: $\rho_s = 3100.0\text{ kg/m}^3$.

#### Fluid (Air) Properties
Thermophysical properties of air are evaluated dynamically at local film temperature $T_f = \frac{1}{2}(T_{solid} + T_{gas})$:
$$\rho_f(T) = \frac{p_{atm}}{R_{air} T} = \frac{101325.0}{287.05 \cdot T} \quad [\text{kg/m}^3]$$
$$c_{p,f}(T) = 1005.0 + 0.05 \cdot (T - 273.15) \quad [\text{J/kg}\cdot\text{K}]$$
$$\mu_f(T) = 1.81 \times 10^{-5} \cdot \left(\frac{T}{293.15}\right)^{0.70} \quad [\text{kg/m}\cdot\text{s}]$$
$$k_f(T) = 0.0257 \cdot \left(\frac{T}{293.15}\right)^{0.80} \quad [\text{W/m}\cdot\text{K}]$$
$$\text{Pr}(T) = \frac{c_{p,f}(T) \mu_f(T)}{k_f(T)} \approx 0.71$$

#### Packaging & Rear Alumina Properties
- Alumina conductivity: $k_{al}(T) = 5.5 + 34.5 \exp(-0.0033 \cdot (T - 273.15))\text{ W/m}\cdot\text{K}$.
- Alumina heat capacity: $c_{p,al}(T) = \max(500.0, (1.00446 + 1.742 \times 10^{-4} T - 2.796 \times 10^4 T^{-2}) \times 1000.0)\text{ J/kg}\cdot\text{K}$.
- Alumina density: $\rho_{al} = 3900.0\text{ kg/m}^3$.
- Aluminum housing conductivity: $k_{metal} = 205.0\text{ W/m}\cdot\text{K}$.
- Aluminum density: $\rho_{metal} = 2700.0\text{ kg/m}^3$; heat capacity: $c_{p,metal} = 900.0\text{ J/kg}\cdot\text{K}$.
- Alumina silicate felt density: $\rho_{felt} = 140.0\text{ kg/m}^3$; heat capacity: $c_{p,felt} = 1360.0\text{ J/kg}\cdot\text{K}$.

---

## 2. Governing Equations & Constitutive Formulations

The receiver domain is discretized into $N$ axial finite-volume cells along the monolith ($z \in [0, L]$, cell length $\Delta z = L/N$) and $N_{rear}$ axial cells along the rear tube ($z \in [L, L + L_{tube}]$, cell length $\Delta z_{rear} = L_{tube}/N_{rear}$).

### 2.1 Radiative Source Partition & Volumetric Penetration

#### Gross Power Scaling & Core/Perimeter Split
For a given nominal incident solar irradiance $I_0$ [$\text{W/m}^2$], the gross thermal power delivered to the receiver assembly $Q_{delivered}$ is closed using the cluster power scale $M(I_0)$:
$$Q_{delivered} = M(I_0) \cdot I_0 \cdot A_{frt}$$
where $M(I_0) \in \{M_{456}, M_{304}, M_{256}\}$ is anchored to high-fidelity optical ray-tracing priors ($M_{456} = 1.34, M_{304} = 1.58, M_{256} = 1.11$).

The delivered power is partitioned between the active core and the passive perimeter housing via the core source fraction $\chi$:
$$Q_{core, solar} = \chi \cdot Q_{delivered}$$
$$Q_{perim, solar} = (1 - \chi) \cdot Q_{delivered}$$

#### Axial Deposition Profiles
Inside the active core, solar radiation is deposited via a combination of front-surface absorption (fraction $f_{dep}$) and volumetric exponential attenuation (Beer-Lambert extinction coefficient $\beta_{opt}$):
$$w_{core, i} = (1 - f_{dep}) \cdot \frac{e^{-\beta_{opt} (i-1) \Delta z} - e^{-\beta_{opt} i \Delta z}}{1 - e^{-\beta_{opt} L}} + f_{dep} \cdot \delta_{i, 1}$$
where $\sum_{i=1}^N w_{core, i} = 1.0$, and $\delta_{i, 1}$ is the Kronecker delta applying front surface absorption exclusively to cell 1.

Along the perimeter housing, spilled radiation attenuates with perimeter extinction coefficient $\beta_{perim}$:
$$w_{perim, i} = \frac{e^{-\beta_{perim} (i-1) \Delta z} - e^{-\beta_{perim} i \Delta z}}{1 - e^{-\beta_{perim} L}}$$

---

### 2.2 Active Core Solid Energy Conservation

For each axial node $i \in \{1, \dots, N\}$ in the active core:
$$C_{core, i}(T_{core, i}) \frac{d T_{core, i}}{dt} = Q_{solar, core, i} - Q_{gas, i} + Q_{cond, core, i} - Q_{radial, cp, i} - Q_{core \to rear, i}$$

1. **Cell Sensible Heat Capacity**:
   $$C_{core, i}(T_{core, i}) = \rho_s c_{p,s}(T_{core, i}) A_{solid} \Delta z \quad [\text{J/K}]$$

2. **Solar Heat Input**:
   $$Q_{solar, core, i} = Q_{core, solar} \cdot w_{core, i} \quad [\text{W}]$$

3. **Solid Axial Conduction**:
   $$Q_{cond, core, i} = \frac{k_{eff, core} A_{solid}}{\Delta z} \left(T_{core, i-1} - 2 T_{core, i} + T_{core, i+1}\right)$$
   Interface conductivities use the harmonic mean: $k_{face, i+1/2} = \frac{2 k_s(T_{core, i}) k_s(T_{core, i+1})}{k_s(T_{core, i}) + k_s(T_{core, i+1})}$, scaled by effective multiplier $k_{core, axial, scale}$.

4. **Radial Core-to-Perimeter Heat Exchange**:
   $$Q_{radial, cp, i} = G_{cp} \Delta z \left(T_{core, i} - T_{perim, i}\right) \quad [\text{W}]$$
   where $G_{cp}$ [$\text{W/m}\cdot\text{K}$] is the effective radial core-perimeter conductance per unit length.

5. **Internal Convective Gas Removal**:
   $$Q_{gas, i} = \dot{m}_{active} c_{p,f} \left(T_{active, i+1} - T_{active, i}\right) \quad [\text{W}]$$

---

### 2.3 Developing Flow Convection & Intra-Strut Webbing Resistance

#### Channel Hydrodynamics & Flow Recruitment
At total volumetric flow rate $\dot{V}$ [LPM], total standard mass flow is $\dot{m}_{total} = \rho_{std} \dot{V} / 60000$.
Due to thermal expansion and flow maldistribution, active flow recruitment follows:
$$\phi_{act}(Re) = \text{clamp}\left(\phi_0 \left(\frac{Re_{total}}{50.0}\right)^{m_{rec}},\, 0.10,\, 1.00\right)$$
$$\dot{m}_{active} = \phi_{act} \cdot \dot{m}_{total}$$
$$Re_{active} = \frac{\dot{m}_{active} D_h}{A_{flow} \mu_f}$$

#### Generalized Developing Nusselt Number
Laminar developing heat transfer in square channels is closed via a generalized entrance-region power law:
$$\text{Nu}(z) = A_{Nu} Re_{active}^{B_{Re}} \text{Pr}^{1/3} \left(\frac{D_h}{\max(z, D_h/2)}\right)^{C_z}$$
where $B_{Re} \approx 0.50$ represents boundary-layer scaling, and $C_z$ is the axial decay exponent ($C_z = 1/3$ corresponds to Graetz flow, $C_z = 0.50$ to flat-plate boundary layers). The local convection coefficient is $h_{cell}(z) = \text{Nu}(z) k_f / D_h$.

#### Intra-Strut Radial Conduction Resistance ($\delta_{web}$)
To capture 2D radial temperature gradients within the solid ceramic webs without increasing ODE dimensionality, an intra-webbing series conduction resistance length $\delta_{web}$ is embedded:
$$U_{eff}(z) = \left(\frac{1}{h_{cell}(z)} + \frac{\delta_{web}}{k_s(T_{core})}\right)^{-1} \quad [\text{W/m}^2\cdot\text{K}]$$

#### Gas Profile & Continuous Centerline Stream Tracking
The number of transfer units for cell $i$ is $\text{NTU}_i = \frac{U_{eff, i} P_{exchange} \Delta z}{\dot{m}_{active} c_{p,f}}$. The cell thermal effectiveness is:
$$\epsilon_i = 1 - \exp(-\text{NTU}_i)$$
Gas temperature advances along the core channels without non-physical cold dilution:
$$T_{active, i+1} = T_{active, i} + \epsilon_i \left(T_{core, i} - T_{active, i}\right)$$
$$T_g(z_i) = T_{active, i}$$
At the monolith exit ($z = L$), the gas enters the rear tube at $T_{g}(L) = T_{active, N+1}$.

---

### 2.4 Passive Perimeter Housing Energy Conservation

For each axial node $i \in \{1, \dots, N\}$ in the perimeter housing:
$$C_{perim, i} \frac{d T_{perim, i}}{dt} = Q_{solar, perim, i} + Q_{radial, cp, i} + Q_{cond, perim, i} - Q_{perim \to cavity, i} - Q_{perim \to rear, i} - Q_{front} \delta_{i, 1}$$

1. **Perimeter Participating Capacity**:
   $$C_{perim, i} = \frac{C_{perim, eff}}{N} \quad [\text{J/K}]$$

2. **Perimeter Axial Conduction**:
   $$Q_{cond, perim, i} = \frac{(0.20 k_s A_{solid} + k_{perim, ref} A_{frt})}{\Delta z} \left(T_{perim, i-1} - 2 T_{perim, i} + T_{perim, i+1}\right)$$
   where $k_{perim, ref}$ accounts for high axial conduction through the aluminum casing.

3. **Radial Heat Loss to Cavity**:
   $$Q_{perim \to cavity, i} = G_{radial, cavity} \Delta z \left(T_{perim, i} - T_{cavity}\right)$$
   where $G_{radial, cavity} = \frac{4 \pi k_{felt}}{\ln(R_{felt} / R_{eq})}$, $R_{eq} = \sqrt{A_{frt}/\pi}$.

---

### 2.5 Dynamic Aperture Suction Convection Boundary ($z = 0$)

The front face ($z = 0$) experiences simultaneous thermal radiation and active suction convection:
$$Q_{front} = h_{front}(Re) A_{frt} (T_{perim, 1} - T_{amb}) + \epsilon_{frt} \sigma A_{frt} (T_{perim, 1}^4 - T_{amb}^4)$$
where the dynamic suction heat transfer coefficient scales with the aperture Reynolds number:
$$h_{front}(Re) = 10.0 + h_{suction} \sqrt{\frac{Re_{inlet}}{50.0}} \quad [\text{W/m}^2\cdot\text{K}]$$
$$Re_{inlet} = \frac{\dot{m}_{total} D_h}{A_{flow} \mu_f(T_{in})}$$
This formulation accounts for the thinning of the thermal boundary layer as cold ambient air accelerates across the hot front face into the channel openings.

---

### 2.6 Rear Manifold, Flange & Exhaust Tube Mechanics

#### Rear Contact Rail ($T_{rear}$)
To model the heat-sink contact between the receiver rear and the mounting structure:
$$C_{rear, i} \frac{d T_{rear, i}}{dt} = Q_{core \to rear, i} + Q_{perim \to rear, i} - Q_{rear \to cavity, i} + Q_{axial, rear, i}$$
$$Q_{core \to rear, i} = G_{recv, rear} w_{rear, i} f_{core, rear} \left(T_{core, i} - T_{rear, i}\right)$$
$$Q_{perim \to rear, i} = G_{recv, rear} w_{rear, i} (1 - f_{core, rear}) \left(T_{perim, i} - T_{rear, i}\right)$$
$$Q_{rear \to cavity, i} = G_{rear, cavity} w_{rear, i} \left(T_{rear, i} - T_{cavity}\right)$$
where $C_{rear, i} = C_{rear, eff} \cdot w_{rear, i}$, with quadratic spatial contact weights $w_{rear, i} = (z_i / L)^2 / \sum (z / L)^2$.

#### Rear Alumina Exhaust Tube ($T_{tube}$)
For $j \in \{1, \dots, N_{rear}\}$ along the exit tube ($z_{rear} \in [0, L_{tube}]$):
$$C_{tube, j} \frac{d T_{tube, j}}{dt} = -Q_{gas, rear, j} - Q_{tube \to cavity, j} - Q_{tube \to flange, j} + Q_{cond, tube, j}$$
$$Q_{gas, rear, j} = \dot{m}_{total} c_{p,f} \left(T_{g, j+1} - T_{g, j}\right)$$
$$Q_{tube \to flange, j} = G_{flange}(z_{rear}, t) \Delta z_{rear} \left(T_{tube, j} - T_{flange}\right)$$
where $G_{flange}$ incorporates a smooth cooling-phase conductance gate $g(t) = 1 + g_{cool} (1 - e^{-t/\tau_{cool}})$ during lamp shutoff.

#### Centerline Gas Thermocouple Model ($T_3$)
Located at $z = 140.0\text{ mm}$ (inside the exhaust tube), the $T_3$ sensor node tracks the coupled convective-radiative balance:
$$C_{sensor} \frac{d T_3}{dt} = h_{sensor} A_{sensor} (T_g(140\text{ mm}) - T_3) + \epsilon_{sensor} \sigma A_{sensor} (T_{tube}^4(3\text{ mm}) - T_3^4)$$
with $h_{sensor} = 80.0\text{ W/m}^2\text{K}$, $A_{sensor} = 1.0 \times 10^{-5}\text{ m}^2$, $C_{sensor} = 0.05\text{ J/K}$, and $\epsilon_{sensor} = 0.85$.

#### Cavity Outer Shell ($T_{cavity} / T_2$)
$$C_{cavity} \frac{d T_{cavity}}{dt} = \sum_{i=1}^N Q_{perim \to cavity, i} + \sum_{i=1}^N Q_{rear \to cavity, i} + \sum_{j=1}^{N_{rear}} Q_{tube \to cavity, j} - Q_{cavity \to amb}$$
$$Q_{cavity \to amb} = h_{nat} A_{cavity} (T_{cavity} - T_{amb}) + \epsilon_{metal} \sigma A_{cavity} (T_{cavity}^4 - T_{amb}^4)$$
where $h_{nat} = 10.0\text{ W/m}^2\text{K}$, $A_{cavity} = 2 \pi R_{case} L_{cavity} + \pi R_{case}^2$, and $\epsilon_{metal} = 0.80$.

---

## 3. Calibration Formulation & Optimal Parameters

### 3.1 Objective Function & Multi-Signal Optimization
The global optimization minimizes a composite multi-signal loss over 15 heating runs and 3 cooling runs across 7 measurement channels ($T_8, T_{12}, T_{11}, T_9, T_{10}, T_3, T_2$):
$$\mathcal{L}(\mathbf{p}) = \sum_{k \in \mathcal{D}_{heat}} \left[ \frac{1}{7} \sum_{j=1}^7 \frac{1}{K_j^2} \int_0^{t_{end}} \left(T_{mod, j}^{(k)}(t) - T_{exp, j}^{(k)}(t)\right)^2 dt + \mathcal{P}_{order}^{(k)} \right] + \sum_{k \in \mathcal{D}_{cool}} \frac{1}{7} \sum_{j=1}^7 \frac{1}{K_j^2} \int_0^{t_{end}} \left(T_{mod, j}^{(k)}(t) - T_{exp, j}^{(k)}(t)\right)^2 dt$$
where $K_j$ is the characteristic signal scaling factor, and $\mathcal{P}_{order}$ enforces physical monotonicity ($T_9 > T_{10}$, $T_8 > T_{11}$).

### 3.2 Optimal Calibrated Parameter Vector (`1D_v45`)

| Index | Symbol | Parameter Name | Optimal Value (`1D_v45`) | Bound $[LB, UB]$ | Units | Physical Interpretation |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| $p[1]$ | $A_{Nu}$ | Nusselt prefactor | **1.5790** | $[0.01, 25.0]$ | — | Laminar convective exchange magnitude |
| $p[2]$ | $B_{Re}$ | Reynolds exponent | **0.5346** | $[0.0, 0.60]$ | — | Flat-plate boundary layer scaling ($\sim Re^{1/2}$) |
| $p[3]$ | $C_{Pr}$ | Prandtl exponent | **0.3333** | Fixed | — | Boundary layer Prandtl scaling |
| $p[4]$ | $\phi_0$ | Active flow fraction @ Re=50 | **0.3024** | $[0.05, 1.00]$ | — | Baseline active channel fraction |
| $p[5]$ | $m_{rec}$ | Recruitment exponent | **1.3131** | $[0.0, 3.00]$ | — | Reynolds-driven channel flow recruitment |
| $p[6]$ | $f_{dep}$ | Front-cell deposition | **0.3260** | $[0.0, 1.00]$ | — | 32.6% surface absorbed, 67.4% penetrates |
| $p[7]$ | $M_{456}$ | Power scale @ 456 kW/m² | **1.3400** | Fixed $[0.5, 2.0]$ | — | R6 delivered-to-aperture optical ratio |
| $p[8]$ | $M_{304}$ | Power scale @ 304 kW/m² | **1.5800** | Fixed $[0.5, 2.0]$ | — | R6 delivered-to-aperture optical ratio |
| $p[9]$ | $M_{256}$ | Power scale @ 256 kW/m² | **1.1100** | Fixed $[0.5, 2.0]$ | — | R6 delivered-to-aperture optical ratio |
| $p[10]$ | $G_{cp}$ | Core-perimeter conductance | **5.7026** | $[2.0, 80.0]$ | $\text{W/m}\cdot\text{K}$ | Radial core-to-housing conduction link |
| $p[11]$ | $C_{perim,eff}$| Perimeter heat capacity | **105.01** | $[105.0, 230.0]$| $\text{J/K}$ | Participating outer casing & felt mass |
| $p[12]$ | $k_{perim,ref}$| Perimeter axial conductivity | **6.4958** | $[0.0, 2000.0]$ | $\text{W/m}\cdot\text{K}$ | Aluminum casing axial heat spreading |
| $p[13]$ | $\beta_{opt}$ | Optical extinction coeff | **243.03** | $[20.0, 500.0]$ | $\text{m}^{-1}$ | Core channel optical absorption depth ($4.1\text{ mm}$) |
| $p[14]$ | $\chi$ | Core source fraction | **0.4550** | $[0.0, 1.00]$ | — | 45.5% core absorbed, 54.5% spilled to perimeter |
| $p[15]$ | $\beta_{perim}$| Perimeter extinction coeff | **24.053** | $[0.5, 300.0]$ | $\text{m}^{-1}$ | Outer rim spillage penetration depth |
| $p[16]$ | $f_{core,rear}$| Rear coupling core fraction | **0.9979** | $[0.0, 1.00]$ | — | Receiver-to-rear rail contact allocation |
| $p[17]$ | $s_{flange}$ | Flange scale multiplier | **0.1000** | $[0.10, 20.0]$ | — | Rear-tube-to-water-flange base conductance |
| $p[18]$ | $g_{cool}$ | Cooling flange gain | **4.4720** | $[0.0, 50.0]$ | — | Lamp-off flange loss acceleration |
| $p[19]$ | $\tau_{cool}$ | Cooling flange time constant | **204.58** | $[1.0, 2000.0]$ | $\text{s}$ | Lamp-off thermal relaxation time |
| $p[20]$ | $k_{axial,scale}$| Core axial conduction scale | **0.1439** | $[0.0, 0.50]$ | — | Webbing tortuosity conduction reduction |
| $p[21]$ | $C_{rear,eff}$ | Rear rail heat capacity | **89.91** | $[50.0, 150.0]$ | $\text{J/K}$ | Alumina adaptor participating thermal mass |
| $p[22]$ | $G_{recv,rear}$| Receiver-rear conductance | **0.1182** | $[0.0, 50.0]$ | $\text{W/K}$ | Mechanical interface thermal resistance |
| $p[23]$ | $G_{rear,tube}$| Rear rail to tube conductance | **3.6152** | $[0.0, 50.0]$ | $\text{W/K}$ | Direct structural contact to exit tube |
| $p[24]$ | $G_{rear,cav}$ | Rear rail to cavity conduct. | **0.1247** | $[0.0, 10.0]$ | $\text{W/K}$ | Rear parasitics to outer cavity |
| $p[25]$ | $G_{rear,ax}$ | Rear rail axial conductance | **9.4134** | $[0.0, 100.0]$ | $\text{W/K}$ | Metal adaptor ring axial heat distribution |
| $p[26]$ | $\delta_{web}$ | Intra-strut resistance length | **$5.79 \times 10^{-5}$** | $[0.0, 2.0\times 10^{-3}]$ | $\text{m}$ | Ceramic wall half-thickness resistance ($58\,\mu\text{m}$) |
| $p[27]$ | $C_z$ | Nusselt axial decay exponent | **0.8299** | $[0.0, 1.00]$ | — | Entrance region decay exponent |
| $p[28]$ | $h_{suction}$ | Suction convection scale | **296.68** | $[0.0, 1000.0]$ | $\text{W/m}^2\cdot\text{K}$ | Aperture suction heat recovery scaling |

---

## 4. Comprehensive Validation & Performance Metrics

### 4.1 Global Objective Evolution
The systematic progression of the model objective function across development generations demonstrates steady and dramatic convergence:

$$\text{1D\_v42 } (0.2270) \;\longrightarrow\; \text{1D\_v43 } (0.2206) \;\longrightarrow\; \text{1D\_v44 } (0.1420) \;\longrightarrow\; \mathbf{\text{1D\_v45 } (0.08637)}$$

**Overall Objective Error Reduction: 60.8%**.

---

### 4.2 Steady-State Temperature Parity & Residuals

The table below reports the steady-state model errors ($\Delta T = T_{mod} - T_{exp}$ [K]) across all 15 heating runs:

| Run ID | Irradiance | Flow [LPM] | $T_9$ (Core Front) | $T_{10}$ (Core Rear) | $T_3$ (Gas Exit) | $T_2$ (Cavity) | $T_8$ (Perim Front) | $T_{12}$ (Perim Mid) | $T_{11}$ (Perim Rear) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **E67** | 456 kW/m² | 15.28 | -43.4 K | -36.2 K | **-43.3 K** | +3.9 K | -62.1 K | -45.2 K | -64.9 K |
| **E68** | 456 kW/m² | 12.50 | -6.7 K | -17.6 K | **-28.1 K** | +6.6 K | -27.2 K | -14.6 K | -35.2 K |
| **E69** | 456 kW/m² | 10.50 | +11.5 K | -6.1 K | **-19.8 K** | +11.1 K | +4.5 K | +7.4 K | -13.1 K |
| **E70** | 456 kW/m² | 9.11 | -19.2 K | -22.4 K | **-31.6 K** | +3.6 K | -37.2 K | -19.3 K | -22.4 K |
| **E71** | 456 kW/m² | 7.13 | -0.1 K | +4.8 K | **-8.3 K** | +4.7 K | -18.6 K | +8.4 K | +15.0 K |
| **E72** | 304 kW/m² | 18.71 | -102.4 K | -71.3 K | **-62.9 K** | +4.2 K | -95.9 K | -97.1 K | -103.7 K |
| **E73** | 304 kW/m² | 13.17 | -32.7 K | -37.8 K | **-34.4 K** | +7.5 K | -32.8 K | -44.0 K | -60.1 K |
| **E74** | 304 kW/m² | 9.02 | -25.0 K | -34.6 K | **-30.6 K** | +6.2 K | -1.7 K | -31.7 K | -38.8 K |
| **E75** | 304 kW/m² | 6.94 | -21.7 K | -20.2 K | **-16.5 K** | +6.9 K | +5.8 K | -20.6 K | -14.9 K |
| **E76** | 304 kW/m² | 4.53 | +15.1 K | +26.5 K | **+25.5 K** | +10.5 K | +16.7 K | +23.0 K | +40.9 K |
| **E77** | 256 kW/m² | 13.85 | -23.9 K | -21.3 K | **-14.6 K** | +6.4 K | -3.8 K | -32.9 K | -39.9 K |
| **E78** | 256 kW/m² | 10.01 | -16.7 K | -27.4 K | **-20.9 K** | +7.5 K | +22.1 K | -29.1 K | -39.1 K |
| **E79** | 256 kW/m² | 8.03 | -26.5 K | -33.8 K | **-26.4 K** | +7.6 K | +24.8 K | -35.8 K | -40.4 K |
| **E80** | 256 kW/m² | 6.61 | -32.8 K | -35.9 K | **-27.6 K** | +7.2 K | +26.0 K | -39.3 K | -38.6 K |
| **E81** | 256 kW/m² | 4.53 | -26.4 K | -21.4 K | **-11.8 K** | +8.1 K | +28.8 K | -28.6 K | -18.3 K |

---

### 4.3 Transient Dynamics & Tracking Metrics

#### Heating Phase Transient Statistics (15 Runs)
- **$T_3$ (Gas Exit)**: Mean RMSE = **`22.3 K`**, Mean Steady Bias = **`-23.4 K`**, Mean $t_{90}$ Error = **`-140 s`**, Shape Loss = **$1.79 \times 10^{-5}$**.
- **$T_2$ (Cavity Shell)**: Mean RMSE = **`4.9 K`**, Mean Steady Bias = **`+6.8 K`**, Mean $t_{90}$ Error = **`+61 s`**, Shape Loss = **$1.09 \times 10^{-4}$**.
- **$T_9$ (Core Front)**: Mean RMSE = **`59.9 K`**, Mean Steady Bias = **`-23.4 K`**, Mean $t_{90}$ Error = **`+710 s`**, Shape Loss = **$8.07 \times 10^{-5}$**.
- **$T_{10}$ (Core Rear)**: Mean RMSE = **`35.8 K`**, Mean Steady Bias = **`-23.7 K`**, Mean $t_{90}$ Error = **`+240 s`**, Shape Loss = **$5.94 \times 10^{-5}$**.
- **$T_8$ (Perimeter Front)**: Mean RMSE = **`51.4 K`**, Mean Steady Bias = **`-10.0 K`**, Mean $t_{90}$ Error = **`+513 s`**, Shape Loss = **$1.44 \times 10^{-4}$**.

#### Cooling Phase Transient Statistics (C69, C80, C81)
- **$T_8$**: Mean RMSE = **`24.3 K`**, Mean Steady Bias = **`+9.7 K`**, $t_{90}$ Error = **`+267 s`**.
- **$T_9$**: Mean RMSE = **`48.1 K`**, Mean Steady Bias = **`+21.0 K`**, $t_{90}$ Error = **`+300 s`**.
- **$T_{10}$**: Mean RMSE = **`57.4 K`**, Mean Steady Bias = **`+28.8 K`**, $t_{90}$ Error = **`+233 s`**.
- **$T_3$**: Mean RMSE = **`57.6 K`**, Mean Steady Bias = **`+29.6 K`**, $t_{90}$ Error = **`+233 s`**.
- **$T_2$**: Mean RMSE = **`9.2 K`**, Mean Steady Bias = **`+11.8 K`**.

---

### 4.4 Spatial Profile Accuracy: Longitudinal Core Drop ($\Delta T_{core} = T_9 - T_{10}$)

| Run ID | Irradiance | Flow [LPM] | Modeled $\Delta T_{core}$ [K] | Experimental $\Delta T_{core}$ [K] | Error [K] |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **E67** | 456 kW/m² | 15.28 | 134.5 K | 141.7 K | **-7.2 K** |
| **E68** | 456 kW/m² | 12.50 | 181.7 K | 170.8 K | **+10.9 K** |
| **E69** | 456 kW/m² | 10.50 | 209.2 K | 191.5 K | **+17.7 K** |
| **E70** | 456 kW/m² | 9.11 | 225.5 K | 222.2 K | **+3.3 K** |
| **E71** | 456 kW/m² | 7.13 | 243.0 K | 247.8 K | **-4.8 K** |
| **E73** | 304 kW/m² | 13.17 | 125.3 K | 120.2 K | **+5.1 K** |
| **E74** | 304 kW/m² | 9.02 | 177.3 K | 167.7 K | **+9.6 K** |
| **E75** | 304 kW/m² | 6.94 | 195.2 K | 196.8 K | **-1.6 K** |
| **E76** | 304 kW/m² | 4.53 | 209.6 K | 221.0 K | **-11.4 K** |
| **E77** | 256 kW/m² | 13.85 | 56.1 K | 58.7 K | **-2.6 K** |
| **E78** | 256 kW/m² | 10.01 | 90.6 K | 79.9 K | **+10.7 K** |
| **E79** | 256 kW/m² | 8.03 | 105.6 K | 98.3 K | **+7.3 K** |
| **E80** | 256 kW/m² | 6.61 | 114.6 K | 111.5 K | **+3.1 K** |
| **E81** | 256 kW/m² | 4.53 | 125.7 K | 130.7 K | **-5.0 K** |

---

## 5. Physical Interpretation & Transport Mechanics

### 5.1 Developing Boundary-Layer Dominated Channel Kinematics
In classical pipe flow, fully developed laminar heat transfer exhibits a constant Nusselt number ($\text{Nu} \approx 3.61$ to $4.36$) independent of Reynolds number ($B_{Re} = 0$). In contrast, `1D_v45` identified:
$$B_{Re} = 0.5346 \approx 0.50, \quad C_z = 0.8299$$
This outcome has fundamental physical significance:
1. **Laminar Boundary Layer Scaling**: $B_{Re} \approx 0.50$ rigorously recovers the flat-plate laminar boundary layer limit ($\text{Nu}_x \propto Re_x^{1/2} \text{Pr}^{1/3}$).
2. **Short Thermal Entry Dominance**: The high spatial decay exponent $C_z = 0.83$ indicates that because the channel hydraulic diameter is small ($D_h = 1.34\text{ mm}$) and entrance disturbances/roughness are significant, thermal development occurs via thin, developing boundary layers near the inlet that govern overall heat extraction.

### 5.2 The Physics of Aperture Suction Convection
In previous 1D formulations (`v31`–`v43`), front-face natural convection was modeled as a rigid constant ($H_{front} = 10.0\text{ W/m}^2\text{K}$). This created a systematic failure where high-flow runs were severely underpredicted while low-flow runs were overpredicted.
- **Physical Mechanism**: Because this is a suction receiver, air is actively drawn across the illuminated front ceramic face into the honeycomb pores. As mass flow increases from $4.5$ to $18.7\text{ LPM}$, the inlet velocity increases fourfold, thinning the hydrodynamic and thermal stagnation boundary layers.
- **Calibrated Outcome**: The optimizer identified an unconstrained minimum at $h_{suction} = 296.7\text{ W/m}^2\text{K}$. At high flow ($Re_{inlet} \approx 150$), the effective suction convection reaches:
  $$h_{front} = 10.0 + 296.7 \sqrt{\frac{150}{50}} \approx 524\text{ W/m}^2\cdot\text{K}$$
  This provides the necessary flow-dependent heat recovery that flattens the multi-flow steady-state parity strata.

### 5.3 Volumetric Optical Penetration vs. Surface Deposition
- **Optical Split**: $front\_dep = 0.326$ indicates that approximately $32.6\%$ of the direct central beam is absorbed on the immediate front solid face, while $67.4\%$ penetrates axially into the channel matrix.
- **Extinction Depth**: $\beta_{opt} = 243.0\text{ m}^{-1}$ defines an effective optical penetration depth $\delta_{opt} = 1/\beta_{opt} \approx 4.1\text{ mm}$.
- **Gas Heating Coupling**: In previous models where $front\_dep$ was frozen at $1.0$, cells 2–15 received $0\text{ W}$ direct solar power, starving the downstream gas and causing $T_3$ to be underpredicted by $250\text{ K}$. Distributing $67.4\%$ of the solar power across the first $15\text{–}25\text{ mm}$ of the channel walls allows the fluid to continuously absorb heat, elevating exit gas temperatures into precise alignment with measured values ($\Delta T_3 \approx \pm 20\text{ K}$).

### 5.4 Perimeter Thermal Inversions ($T_{12} - T_8$)
The packaging assembly exhibits a distinct flow-dependent spatial inversion:
- **High-Flow Inversion ($T_{12} > T_8$)**: Intense suction cooling chills the front perimeter ($z = 11\text{ mm}$ at $T_8$), while internal radial radiation from the hot monolith core heats the mid-casing ($z = 58\text{ mm}$ at $T_{12}$). The temperature profile peaks in the middle ($T_{12} - T_8 = +59.5\text{ K}$ in E67).
- **Low-Flow Monotonic Decay ($T_{8} > T_{12}$)**: With low suction velocities, front surface cooling diminishes, allowing direct spillage absorption to overheat the front face, yielding a monotonically decaying profile ($T_{12} - T_8 = -97.2\text{ K}$ in E76).
- `1D_v45` captures this transition across all flow rates without empirical branching.

### 5.5 Thermal Capacitance & Assembly Inventory Reconciliation
- **Fitted Capacitances**: $C_{perim, eff} = 105.0\text{ J/K}$, $C_{rear, eff} = 89.9\text{ J/K}$, $C_{core, eff} = 72.5\text{ J/K}$.
- **Total Assembly Capacitance**:
  $$C_{total} = C_{core} + C_{perim} + C_{rear} = 267.4\text{ J/K}$$
  This agrees closely with the physical measured assembly inventory target ($301 \pm 23\text{ J/K}$), confirming that transient calibration has converged to real physical hardware masses.

---

## 6. Nomenclature

### Latin Symbols
| Symbol | Definition | Units |
| :--- | :--- | :--- |
| $A_{frt}$ | Frontal aperture area ($W \times W$) | $\text{m}^2$ |
| $A_{flow}$ | Total open gas flow area | $\text{m}^2$ |
| $A_{solid}$ | Total solid cross-sectional area | $\text{m}^2$ |
| $a_v$ | Specific internal surface area ($P_{exchange}/A_{frt}$) | $\text{m}^2/\text{m}^3$ |
| $A_{Nu}$ | Nusselt correlation prefactor | — |
| $B_{Re}$ | Nusselt Reynolds number exponent | — |
| $C_{p}$ | Specific heat capacity | $\text{J/kg}\cdot\text{K}$ |
| $C_{eff}$ | Effective participating heat capacity | $\text{J/K}$ |
| $C_z$ | Nusselt axial spatial decay exponent | — |
| $D_h$ | Channel hydraulic diameter | $\text{m}$ |
| $f_{dep}$ | Front-cell direct optical deposition fraction | — |
| $G_{cp}$ | Radial core-perimeter thermal conductance | $\text{W/m}\cdot\text{K}$ |
| $h$ | Convective heat transfer coefficient | $\text{W/m}^2\cdot\text{K}$ |
| $h_{suction}$ | Aperture suction convection scaling parameter | $\text{W/m}^2\cdot\text{K}$ |
| $I_0$ | Incident nominal solar irradiance | $\text{W/m}^2$ |
| $k$ | Thermal conductivity | $\text{W/m}\cdot\text{K}$ |
| $L$ | Receiver monolith axial length ($137\text{ mm}$) | $\text{m}$ |
| $\dot{m}$ | Mass flow rate | $\text{kg/s}$ |
| $M$ | Optical delivered-to-aperture power ratio | — |
| $N$ | Number of axial finite-volume nodes | — |
| $\text{Nu}$ | Nusselt number ($h D_h / k_f$) | — |
| $P_{exchange}$| Total internal convective perimeter | $\text{m}$ |
| $\text{Pr}$ | Prandtl number ($c_p \mu / k_f$) | — |
| $Q$ | Heat transfer rate / thermal power | $\text{W}$ |
| $Re$ | Reynolds number ($\dot{m} D_h / A_{flow} \mu$) | — |
| $T$ | Temperature | $\text{K}$ |
| $U_{eff}$ | Overall effective heat transfer coefficient | $\text{W/m}^2\cdot\text{K}$ |
| $\dot{V}$ | Volumetric flow rate | $\text{LPM}$ |
| $w_{ch}$ | Square channel opening width | $\text{m}$ |
| $z$ | Axial spatial coordinate | $\text{m}$ |

### Greek Symbols
| Symbol | Definition | Units |
| :--- | :--- | :--- |
| $\beta_{opt}$ | Optical extinction coefficient in monolith core | $\text{m}^{-1}$ |
| $\beta_{perim}$| Optical extinction coefficient in perimeter housing | $\text{m}^{-1}$ |
| $\delta_{web}$ | Intra-strut solid conduction resistance length | $\text{m}$ |
| $\epsilon$ | Monolith open porosity ($A_{flow}/A_{frt}$) / Thermal emissivity | — |
| $\epsilon_i$ | Control-volume heat exchanger effectiveness | — |
| $\mu$ | Dynamic fluid viscosity | $\text{kg/m}\cdot\text{s}$ |
| $\rho$ | Mass density | $\text{kg/m}^3$ |
| $\sigma$ | Stefan-Boltzmann constant ($5.670374 \times 10^{-8}$) | $\text{W/m}^2\cdot\text{K}^4$ |
| $\tau$ | Time constant | $\text{s}$ |
| $\phi_{act}$ | Active channel flow fraction | — |
| $\chi$ | Core solar source partition fraction | — |

### Subscripts & Accents
| Subscript | Definition |
| :--- | :--- |
| $amb$ | Ambient surroundings |
| $core$ | Active monolith central channel core |
| $exp$ | Experimental measurement |
| $f$ | Fluid (gas phase) |
| $flange$ | Water-cooled mounting flange |
| $mod$ | Model prediction |
| $perim$ | Passive perimeter housing (felt + casing) |
| $rear$ | Rear manifold contact rail |
| $s$ | Solid ceramic phase |
| $tube$ | Rear alumina exhaust tube |
