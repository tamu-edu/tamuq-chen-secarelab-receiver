# 1D_v1: fast axial solid-gas receiver model

## 1. Purpose and present status

`1D_v1.jl` is a new, conservative one-dimensional model for the 137 mm SiC
honeycomb receiver. It predicts:

- the transient solid temperature along the receiver;
- the quasi-steady gas temperature along every channel and at the outlet;
- the local and length-averaged effective heat-transfer coefficient;
- the instantaneous partition of absorbed power into storage, gas heating, and
  environmental losses.

No legacy source file was modified. The new model deliberately uses a compact
finite-volume formulation rather than the symbolic PDE formulation in
`1D_v1.exp.All.jl`. This makes it fast enough for multi-experiment inverse
fitting and eventual use for long receivers.

The current heat-transfer correlation is a defensible calibration seed, not a
final reported material correlation. The staged calibration code is in
`calibrate_1D_v1.jl`.

## 2. Experimental interpretation

The coordinate is

$$
z=0 \quad \text{at the irradiated air inlet}, \qquad z=L \quad \text{at the rear gas outlet},
$$

with $L=0.137$ m. The manuscript geometry is used by default:

| Quantity | Default |
| --- | ---: |
| Receiver width | 19 mm |
| Channel width | 1.5 mm |
| Wall thickness | 0.4 mm |
| Number of square channels | 100 |
| Channel hydraulic diameter | 1.5 mm |
| Total gas-solid exchange perimeter | 0.600 m |

The solid thermocouples used for primary 1D fitting are TI-8, TI-9, and TI-10
at approximately 11, 58, and 107 mm. TI-3 is the downstream outlet-gas
measurement.

TI-11 and TI-12 are retained as diagnostics. Figure 3 of the manuscript places
TI-11 at 58 mm and TI-12 at 107 mm, whereas both import scripts comment the
opposite mapping. The temperature pairing in the data is more consistent with
TI-12 near TI-9 and TI-11 near TI-10. This conflict must be resolved from the
installation record before TI-11/TI-12 are assigned to physical gas nodes.

## 3. Governing model

### 3.1 Solid finite-volume balance

The monolith is divided into $N$ axial cells of width $\Delta z$. For solid
cell $i$,

$$
C_i\frac{dT_{s,i}}{dt}
=\dot Q_{z,i-1/2}-\dot Q_{z,i+1/2}
+w_i\dot Q_{\mathrm{abs}}
-\dot Q_{g,i}-\dot Q_{\mathrm{side},i}.
$$

The cell heat capacity is

$$
C_i=\rho_s c_{p,s}(T_{s,i})A_s\Delta z,
$$

and the inter-cell conductive heat rate is

$$
\dot Q_{z,i+1/2}
=k_{i+1/2}A_s\frac{T_{s,i}-T_{s,i+1}}{\Delta z},
$$

where $k_{i+1/2}$ is the harmonic mean. Internal conduction is added to one
cell and removed from its neighbor, so it cancels exactly in the receiver-wide
energy balance.

The side loss is currently a calibrated effective insulation conductance,

$$
\dot Q_{\mathrm{side},i}
=U'_{\mathrm{side}}\Delta z(T_{s,i}-T_\infty).
$$

The default $U'_{\mathrm{side}}=0.35$ W m$^{-1}$ K$^{-1}$ is close to the
steady cylindrical resistance estimate for 0.08 W m$^{-1}$ K$^{-1}$ felt. It
does not reproduce insulation storage; a dynamic radial insulation state is a
logical next extension if TI-2 must also be predicted.

### 3.2 Front and rear boundaries

The front loss is

$$
\dot Q_{\mathrm{front}}
=h_fA_{\mathrm{frt}}(T_{s,1}-T_\infty)
+\epsilon_f\sigma A_{\mathrm{frt}}(T_{s,1}^4-T_\infty^4).
$$

The rear/adaptor loss is

$$
\dot Q_{\mathrm{rear}}
=U_{\mathrm{rear}}(T_{s,N}-T_\infty)
+\epsilon_r\sigma A_{\mathrm{frt}}(T_{s,N}^4-T_\infty^4).
$$

The rear term is an effective representation of the adaptor and downstream
hardware. Cooling data should determine it before optical parameters are fit.

### 3.3 Solar deposition

The total absorbed power is

$$
\dot Q_{\mathrm{abs}}(t)
=\eta_{\mathrm{abs}}q''_{\mathrm{solar}}(t)A_{\mathrm{frt}}.
$$

The axial weights use a normalized Beer-Lambert form,

$$
w_i\propto e^{-\beta z_{i-1/2}}-e^{-\beta z_{i+1/2}},
\qquad \sum_iw_i=1.
$$

This permits surface-like or volumetric absorption without changing the total
absorbed power. The optional front-deposition fraction can place a specified
part of the power directly in the first cell. Cooling does not contain
$\eta_{\mathrm{abs}}$ or $\beta$, so these parameters are reserved for the
heating stage.

### 3.4 Quasi-steady gas model

Gas storage is neglected because the channel residence time is much shorter
than the ceramic response. Within cell $i$,

$$
T_{g,i+1}
=T_{g,i}+(1-e^{-NTU_i})(T_{s,i}-T_{g,i}),
$$

where

$$
NTU_i=\frac{h_iP_{\mathrm{exchange}}\Delta z}{\dot m c_{p,g}}.
$$

The solid-to-gas heat rate is evaluated from the same enthalpy rise,

$$
\boxed{
\dot Q_{g,i}
=\dot m c_{p,g}(T_{g,i+1}-T_{g,i}).
}
$$

Therefore gas heat gain and solid heat loss are identical at every axial cell,
not merely after integration.

Standard L/min is converted to mass flow using air density at the measured
inlet temperature. Reynolds number is based on the channel hydraulic diameter:

$$
Re=\frac{\dot mD_h}{A_{\mathrm{flow}}\mu_g}.
$$

### 3.5 Heat-transfer parameterization

The implemented cumulative-average form is

$$
\overline{Nu}(z)
=C_{Nu}Re^mPr^n
\left(\frac{T_{\mathrm{film}}}{T_{\mathrm{ref}}}\right)^r
\left(1+a\frac{D_hRePr}{z}\right)^p.
$$

The local Nusselt number is obtained from

$$
Nu_{\mathrm{local}}(z)
=\frac{d}{dz}\left[z\overline{Nu}(z)\right].
$$

This is important for scale-up: integrating the local coefficient recovers the
length-dependent average correlation. Setting $a=0$ gives a uniform power law.
Setting $C_{Nu}=3.66$, $a=0.095$, and $p=0.45$ reproduces manuscript
correlation 1 in cumulative-average form.

## 4. Provisional effective heat-transfer correlation

A direct steady diagnostic was performed for all 15 heating runs. For each
case, the external wall profile was interpolated through TI-8/TI-9/TI-10 and a
constant $h$ was inverted to reproduce TI-3. The resulting nominal-area
coefficient spans

$$
1.16 \le h_{\mathrm{eff}} \le 8.88\ \mathrm{W\,m^{-2}K^{-1}}.
$$

A joint Reynolds-temperature fit gives

$$
\boxed{
Nu_{\mathrm{eff}}
=6.1\times10^{-4}Re^{1.3835}
\left(\frac{T_{\mathrm{film}}}{600\ \mathrm{K}}\right)^{-0.4883}.
}
$$

At 9 L/min the current default model gives a length-mean value near
$h_{\mathrm{eff}}=4.2$ W m$^{-2}$ K$^{-1}$.

This must be reported as an **effective receiver coefficient**, not yet as the
intrinsic channel film coefficient. A fully developed laminar value $Nu=3.66$
would imply roughly 90-140 W m$^{-2}$ K$^{-1}$ over the relevant air-property
range. The much lower inverse value can include:

- nonuniform flow among the 100 channels;
- only part of the nominal channel perimeter being effectively active;
- a difference between external thermocouple temperature and internal wall
  temperature;
- TI-3 radiation/conduction bias and heat loss between the monolith and TI-3;
- unmodeled adaptor storage.

The discrepancy is a research result to resolve, not a reason to force the
literature coefficient into the data.

## 5. Why cooling must be fitted first

Single-exponential diagnostics already show that the 0D assumption is not
adequate:

| Cooling run | $\tau_{T8}$ (s) | $\tau_{T9}$ (s) | $\tau_{T10}$ (s) | $\tau_{T3}$ (s) |
| --- | ---: | ---: | ---: | ---: |
| C69 | 221 | 351 | 670 | 1058 |
| C80 | 426 | 620 | 1132 | 1616 |
| C81 | 585 | 900 | 1583 | 2099 |

These are apparent, full-record one-mode fits and not final physical time
constants. Their systematic axial ordering is the key result: stored heat is
redistributed along the solid, and the gas outlet responds to the downstream
profile rather than one average solid temperature.

`calibrate_1D_v1.jl` implements this sequence:

1. Cooling fit: conductivity scale, heat-capacity scale, side loss, rear loss,
   effective heat-transfer correlation, and TI-3 response time.
2. Heating fit: absorbed fraction, optical extinction coefficient, front loss,
   and a tightly constrained refinement of the heat-transfer coefficient.
3. Final global fit: should be performed only after the sensor mapping and
   irradiance histories are resolved.

A short workflow-verification fit reduced the normalized cooling objective to
about 0.021, but it pushed the effective heat capacity to roughly 2.6 times the
monolith-only value and the allowed TI-3 response time to its 200 s bound. A
short heating fit also did not converge. Those values are intentionally not
baked into `default_parameters()`: they indicate missing adaptor/downstream
thermal storage or sensor physics rather than a completed identification.

## 6. Usage

```julia
include("1D_v1.jl")
using .Receiver1D

p = default_parameters()
op = OperatingCondition(
    irradiance = 304e3,
    flow_lpm = 9.0,
    inlet_temperature = 295.6,
    ambient_temperature = 295.6,
)

result = simulate(p, op, 0.0:10.0:3600.0;
                  initial_temperature=295.6)

Tgas_out = result.gas_temperature[end, :]
Tsolid_58mm = solid_at(result, 0.058)
h_summary = heat_transfer_summary(result)
```

Measured histories can be supplied using `linear_history`. The calibration
workflow reads the raw CSV files directly and uses the measured sum of the four
flow-controller signals.

Run verification with:

```text
julia --startup-file=no test/test_1D_v1.jl
```

## 7. Problems in the legacy files that are corrected here

The older files remain unchanged. The following issues should be corrected if
they are reused:

1. `1D_v1.exp.All.jl` uses 19 mm (`Lc`) rather than the 1.5 mm channel
   hydraulic diameter in Reynolds number and $h=Nu k/L_c$.
2. Total receiver volumes and total exchange areas appear inside pointwise PDE
   equations. This repeats full-receiver heat exchange at every axial location
   and is dimensionally inconsistent with a distributed equation.
3. The fluid mass flow is based on density evaluated at 1000 K even though the
   controllers report standard L/min.
4. The Stefan-Boltzmann constant is entered as $5.17\times10^{-8}$ rather than
   $5.670374419\times10^{-8}$ W m$^{-2}$ K$^{-4}$.
5. `import_exp2.jl` replaces measured L/min values (about 4.5-18.3) with values
   near 0.6-1.6 while the model still interprets them as L/min.
6. The legacy optimization treats TI-11 and TI-12 as bulk gas temperatures,
   although the manuscript describes them as internal-wall thermocouples and
   their axial mapping is inconsistent.
7. The parameter replacement relies on dictionary value order, which is not a
   safe mapping to symbolic parameter order.
8. `import_exp_0D.jl` applies unexplained irradiance multipliers of 1.15 and
   0.7 to the 304 and 256 kW m$^{-2}$ groups. The new calibration uses the
   manuscript nominal fluxes until measured time histories are available.
9. The experimental manuscript uses 19/1.5/0.4 mm geometry, while `0D_v3.jl`
   uses the COMSOL 20/1.65/0.35 mm geometry. Both may be valid for different
   representations, but they cannot be mixed within one fitted correlation.

## 8. Missing information needed for a publishable coefficient

Before calling the fitted parameter an intrinsic heat-transfer coefficient,
the following are needed:

- resolved TI-11/TI-12 axial locations and whether their junctions contact the
  wall or sample bulk gas;
- the actual irradiance time history for each heating experiment, rather than
  a constant nominal value;
- a TI-3 sensor/lead/downstream heat-loss characterization;
- confirmation of how many channels carry flow and the rear adaptor flow
  distribution;
- one independent prior on SiC effective heat capacity or axial conductivity;
- preferably a no-flow cooling run to separate environmental loss from gas
  heat removal;
- a dynamic insulation/adaptor node if TI-2 and long cooling tails are included
  in the objective.

With those items resolved, the same code can produce a defensible shared
$Nu(Re,Pr,T,L/D_h)$ parameterization and then be tested on receiver lengths not
used for calibration.
