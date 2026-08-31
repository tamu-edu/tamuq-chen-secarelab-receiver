# Continuum Theory Manual: 1D_v42 Receiver Model

## 1. Introduction and Architectural Overview

The `1D_v42` model is a dynamic finite-volume continuum solver built to extract and map effective macroscopic heat transfer coefficients for a monolithically structured square-channel solar receiver. It acts as an "Entire Converter Model" that bridges the high-fidelity geometry of a single channel with the fully-coupled bulk behavior of the reactor hardware. 

The domain consists of a two-zone partition:
1. **The Active Core (`T_core`)**: The central participating solid region representing the illuminated monolith channels carrying the primary working fluid.
2. **The Passive Perimeter (`T_perim`)**: The exterior, un-flowed outer ring of the monolith and its encapsulating housing. 

This core/perimeter arrangement is bounded longitudinally by a **Rear Reservoir (`T_rear`)** (a lumped thermal mass representing the mechanical integration flange) and a **Rear Tube (`T_tube`)** representing the exhaust routing.

## 2. Radiative Source Reparametrization (v42)

The most critical upgrade in the `v42` iteration is the imposition of strict energy conservation on the incoming radiant bounds via a decoupled structural parameterization (`M`, `chi`).

### 2.1 Total Delivered Energy

Rather than fitting a non-identifiable product of unconstrained scales and capture fractions, the total thermal energy deposited into the solid receiver aperture is rigorously anchored to high-fidelity volumetric R6 full-model predictions. 

Given an incident nominal aperture irradiance profile $Q_{aperture}$, the gross delivered power $Q_{delivered}$ is defined by:
```math
Q_{delivered} = M \cdot Q_{aperture}
```
where $M$ is the explicit, irradiance-dependent total delivered-to-aperture ratio. In `1D_v42`, $M$ is pinned firmly to the R6 closures for the respective experimental flux clusters:
- $M_{456} = 1.34$
- $M_{304} = 1.58$
- $M_{256} = 1.11$

### 2.2 Core-Perimeter Partition ($\chi$)

The total absorbed energy is distributed into the active core and the passive housing perimeter based on the scalar partition factor $\chi$ (the `core_source_fraction`).
```math
Q_{core, solar} = \chi \cdot Q_{delivered}
```
```math
Q_{perim, solar} = (1 - \chi) \cdot Q_{delivered}
```

Axially, this energy is deposited throughout the porous network according to standard Beer-Lambert exponential attenuation laws:
```math
q'''_{core}(z) \propto \exp(-\beta_{opt} z)
```
```math
q'''_{perim}(z) \propto \exp(-\beta_{perim} z)
```

## 3. Governing Energy Conservation Equations

The 1D finite-volume discretization splits the domain longitudinally into $N$ control volumes of length $\Delta z$.

### 3.1 Active Core Transport

For a single core control volume $i$, the transient thermal storage is governed by convective heat removal, internal conduction, and external coupling:
```math
C_{core, i} \frac{d T_{core, i}}{dt} = Q_{cell, core} - Q_{gas, i} + Q_{rad, in} + Q_{cond, core, i} - Q_{radial, cp, i}
```
- $C_{core, i} = \rho_s c_{p,s}(T_{core, i}) V_{cell}$: Local solid sensible heat capacity.
- $Q_{cell, core}$: The partitioned volumetric solar source for cell $i$.
- $Q_{gas, i}$: Internal convection to the working fluid, parameterized by a modified Nusselt relation ($A_{Nu}, B_{Re}$).
- $Q_{cond, core, i}$: Axial solid conduction through the structured webbing.
- $Q_{radial, cp, i}$: The lumped inter-zone conductive link moving heat from the core to the perimeter.

### 3.2 Passive Perimeter Transport

The perimeter does not directly interact with the working gas. Its thermal balance is ruled by:
```math
C_{perim, i} \frac{d T_{perim, i}}{dt} = Q_{cell, perim} - Q_{rad, in} - Q_{perim, cavity} + Q_{cond, perim, i} + Q_{radial, cp, i}
```
- $C_{perim, i}$: The localized capacity of the perimeter housing (constrained by an effective global fitting target `C_perim_eff`).
- $Q_{perim, cavity}$: Radiative and natural convective losses to the surrounding domain.

### 3.3 Internal Convection Coupling

Heat extraction to the continuous gas phase occurs via:
```math
Q_{gas, i} = h A_{s} (T_{core, i} - T_{gas, i})
```
The effective heat transfer coefficient $h$ is closed using the fitted laminar flow correlation:
```math
Nu = A_{Nu} Re^{B_{Re}} Pr^{C_{Pr}}
```
Where $A_{Nu}$ and $B_{Re}$ are the primary optimized parameters guiding macroscopic heat extraction, and $Re$ varies locally based on dynamically evaluated fluid properties and fractional active area mapping $\phi_0$.

## 4. Hardware Manifold and Reservoir Coupling

The rear of the continuum domain integrates into a lumped, un-flowed metal structure representing the physical exhaust flange, split into a thermal reservoir and a tube boundary.

```math
C_{rear} \frac{d T_{rear}}{dt} = Q_{core \to rear} + Q_{perim \to rear} - Q_{rear \to cavity} - Q_{rear \to tube}
```
Where the boundary fluxes are mediated by global series conductances ($G_{rear, tube}$, $G_{recv, rear}$) and a geometric split factor ($f_{core\_rear}$) that bridges the final core and perimeter faces directly to the reservoir mass.

By optimizing these macroscopic nodes natively against physical transient operations, the `1D_v42` network ensures the fitted transfer parameters (like $A_{Nu}$) are genuinely isolated from parasitic thermal mass effects or source bookkeeping anomalies.
