# 2D v19 theory, calibration, and interpretation manual

## Purpose and status

V19 tests a stricter factorization of the inverse problem than v18:

1. normalized side-wall profiles select optical depth shape;
2. measured effectiveness and apparent `Nu(Re)` select an integrated
   wall--gas conductance;
3. cooling records select rear tube/flange topology and a physical T3 probe
   response; and
4. only then are bounded felt properties and lamp-group power factors fitted.

Every stage has an independent transfer gate. A later absolute-temperature
fit cannot compensate for an incorrect earlier signature. T9 and T10 remain
diagnostic and never enter parameter selection.

## Geometry and retained physics

V19 retains the full v18 receiver/assembly:

- all 100 square channels, represented exactly by 15 D4-symmetry orbits;
- common-pressure flow allocation among those orbits;
- the measured outer groove obstruction, with approximately the central
  13 mm diameter unobstructed;
- temperature-dependent air density, viscosity, conductivity, and heat
  capacity;
- SiC receiver, filled alumina-silicate felt, aluminum casing, solid alumina
  adaptor, alumina exit tube, metal/water-cooled flange, and the verified
  rear casing-to-flange contact;
- side thermocouples at 11, 58, and 107 mm; and
- T3 configured 3 mm downstream of the receiver exit in the rear-tube
  coordinate.

The 15 numerical gas histories do not represent 15 physical channels or a
quarter receiver. Orbit multiplicities recover 100 physical channels, their
full flow, heat-transfer area, solid inventory, and absorbed power.

## A. Conservative optical-shape model

The Beer--Lambert baseline is compared with the same centered near/deep basis
used in v18:

```text
w_gj = W_g [(1-f_deep) w_near,j + f_deep w_deep,j(L_deep)].
```

Both axial components are nonnegative and normalized separately. Therefore
every orbit retains the inherited radial power marginal `W_g`, and total
absorbed power is exactly unchanged. This is not a side-weighted irradiance
pattern.

Only the normalized steady side excess profile selects the source:

```text
p = (T8-Tamb, T12-Tamb, T11-Tamb) / sum(Tside-Tamb).
```

The objective also reports the 58-minus-11 mm inversion, the 107-minus-58 mm
rear fall, recovered observed middle peaks, and false peaks. Absolute
temperature, T3, T2, T9, and T10 do not select optical shape.

## B. Integrated wall--gas exchange

V19 removes the v18 local heat-transfer multiplier. It represents the
manuscript's apparent receiver-scale correlation as an equivalent
single-physical-channel conductance:

```text
Nu_bar = A_Nu Re_bar^n (Pr/0.70)^m
UA_ch  = Nu_bar k_film / Dh * P_ch * L.
```

`Re_bar` uses the measured total mass flow divided by 100 channels. Each
physical channel receives the same imposed integrated `UA_ch`; an orbit with
multiplicity `m_g` represents `UA_orbit=m_g UA_ch`, and
`UA_receiver=100 UA_ch`. This deliberately adds no fitted radial `UA`
allocation law and should not be confused with a derived local-channel law.

The axial allocation is

```text
UA_gj = UA_ch I_gj / sum_j(I_gj),
I_gj  = integral_cell phi_Gz(Re_g,Pr_g,z,Tw,Tg) dz.
```

`phi_Gz` is the inherited finite-asymptote square-channel Graetz Nusselt
kernel. Unlike `Re_bar`, its Reynolds number uses the actual per-channel flow
of orbit `g` after common-pressure redistribution. Cell integrals use
three-point Gauss quadrature and are renormalized after every kernel
iteration, so

```text
sum_j UA_gj = UA_ch
```

to machine precision on C, M, and F. Flow maldistribution changes local
Graetz shape and channel NTU, but cannot create or destroy integrated
conductance.

For every cell the gas update is conservative:

```text
epsilon_gj = 1 - exp[-UA_gj/(mdot_g cp_g)]
Tg,out     = Tg,in + epsilon_gj (Tw-Tg,in)
Qgas_gj    = mdot_g cp_g (Tg,out-Tg,in).
```

Exactly the same discretized `Qgas_gj` is removed from the corresponding
solid state, with `cp` frozen at each cell inlet. The gas reference
temperature and common-pressure allocation are iterated together. Thus gas
density and velocity remain temperature-dependent; the input is fixed
standard mass flow from the calibrated MFC convention, not a fixed in-channel
velocity. Rear mixing currently uses a `cp(T)T` weighted approximation rather
than an exact integral-enthalpy inversion; v19 therefore does not claim exact
thermodynamic mixing for temperature-dependent `cp`.

Stage B fits `A_Nu` and `n` only to measured/model effectiveness and apparent
Nusselt number. Its acceptance targets are independent of temperature parity:

```text
training and held-out effectiveness RMSE <= 0.05
training and held-out log(Nu_model/Nu_measured) RMSE <= 0.10
A_Nu within 20% of measured fit
|n - 1.44346| <= 0.10
training maximum effectiveness <= 0.85.
```

There is an important conditional assumption: Stage B compares calculated
gas at the receiver exit/rear-tube inlet with measured T3, while Stage C
argues that T3 is a downstream, cooled, dynamic observation. Consequently
Stage-B `A_Nu,n` can absorb T3 position/probe bias. They remain conditional
apparent parameters unless a constrained B--C iteration or like-for-like
outlet measurement closes this inconsistency.

## C. Receiver outlet versus T3 observation

The `cp`-weighted mixed receiver-outlet approximation and reported T3 are
separate outputs. T3 is modelled as a small sheath/bead at a configured 3 mm
downstream position, driven by:

- cross-flow convection from the calculated local rear-tube gas;
- radiation from the calculated local alumina tube; and
- axial stem conduction to the water-cooled flange.

Per unit exposed sheath area,

```text
C''p dTp/dt =
    hgas (Tgas-Tp)
  + epsilon_p sigma (Ttube^4-Tp^4)
  + G''stem (Twater-Tp).
```

The gas coefficient assumes cylinder cross-flow using local gas properties,
calculated mass flow, and a fixed 1.5 mm diameter. Radiation assumes the
effective unit view used by the equation. `G''stem` is an effective lead/stem
sink to fixed water temperature, not a resolved stem geometry. The probe is
a one-way observation state and does not alter receiver energy. For cooling,
its initial condition is measured T3 at `t0`, not modelled outlet gas.

The 3 mm coordinate is not fully reconciled with the manuscript description
of T3 near 136 mm global axial position. Also, on C the local rear-tube
temperature at 3 mm is represented by the first finite-volume cell. These
coordinate and sampling assumptions must be treated as uncertainty, not as
verified sensor geometry.

The rear alumina tube also has an optional effective distributed flange sink
from 28 mm in the local rear-tube coordinate to its end:

```text
q'_tube-flange =
    h_tube-flange 2 pi r_o (Ttube-Twater).
```

This represents full-circumference tube-to-fixed-water cooling, not an
explicit metal-flange spreading state. When this distributed contact is
active, the standard v19 parameter builder disables the old
terminal-cell-only tube contact to prevent double counting.

Stage C first ablates the distributed contact on C69/C80 using side and felt
cooling histories only; no tube-temperature measurement is fitted. It then
screens global probe parameters against T3 on the discrete grids
`C''p={200,600,1200,2400,3000} J/m2/K` and
`G''stem={0,20,60,120,200} W/m2/K`. C81 is untouched validation. A selected
grid edge is a failed identification.

## D. Felt and power nuisance fit

Stage D is allowed to proceed diagnostically even if A--C fail, but cannot
turn their parameters into validated coefficients.

Felt conductivity and heat capacity are fitted only on cooling C69/C80 over:

```text
k_felt scale  = 0.70, 1.00, 1.30
Cp_felt scale = 0.75, 1.00, 1.50.
```

The best three C-mesh pairs are reranked on M; C81 remains held out. Felt
contact stays fixed at 0.30 and casing/flange conductance at 20 W/K.

Power is then profiled directly on M, independently for the three irradiance
groups:

```text
456 group: 1.05, 1.1225, 1.195, 1.2675, 1.34
304 group: 1.23, 1.3175, 1.405, 1.4925, 1.58
256 group: 0.84, 0.9075, 0.975, 1.0425, 1.11.
```

Only T8/T12/T11, T3, and T2 enter these objectives. A power optimum at a
bracket endpoint is reported as unresolved power accounting and is not
silently extrapolated.

## Numerical and conservation requirements

V19 uses the exact nested v18 C/M/F meshes:

| mesh | felt radial | casing radial | axial | rear tube |
|---|---:|---:|---:|---:|
| C | 3 | 1 | 24 | 15 |
| M | 6 | 2 | 48 | 30 |
| F | 12 | 4 | 96 | 60 |

Smoke tests require exact source power on all meshes, recovery of the
prescribed per-channel and 100-channel integrated `UA`, summed mass-flow
closure, equal-pressure and gas-reference convergence, and exact cooling T3
initialization. The inherited/correction M-mesh ledger is supplemented by an
independent audit that compares the actual capacity-weighted v19 derivative
increment with the separately summed gas-exchange correction and distributed
tube/flange loss. That independent correction audit runs on C, M, and F; the
legacy total-energy ledger remains restricted to M.
Fitted-regime M-to-F primary-history RMS must be below 10 K and maximum below
20 K for all representatives.

## Calibration outcome

Every staged identification fails:

| stage | diagnostic minimum | decisive failure |
|---|---|---|
| A, optical shape | `f_deep=0.90`, `L_deep=120 mm` | 7/10 true peaks and three false peaks |
| B, integrated exchange | `A_Nu=3.80e-4`, `n=1.54` | train/holdout log-Nu RMSE 0.1116/0.1029 and prefactor outside 20% band |
| C, tube/probe | `h=400`, `C''=3000`, `G''stem=0` | all grid edges; T3 RMSE 73.3 K training, 39.3 K C81 |
| D, felt | `k=0.70`, `Cp=0.75` | both lower edges; side/T3 cooling RMSE 49.5/71.0 K |
| D, power | `1.05/1.23/0.84` | all closure floors; heating side/T3 errors remain large |

A four-point constrained Stage-B refinement confirms that its near miss is
not grid resolution. The best admissible point is exactly the upper-right
acceptance corner, `A_Nu=3.693432e-4`, `n=1.54346`, but full-training log-Nu
RMSE worsens to 0.1203. No holdout or M confirmation is authorized.

V19 is therefore rejected regardless of the later parity plots. The full
validation files under `summaries/2D_v19/` quantify errors, numerical
transfer, observation sensitivity, and conservation; they do not convert
these diagnostic boundary values into validated coefficients. The stage gate
status, boundary activity, full sensor errors, mesh transfer, and final model
decision must be read together.

Final M-mesh errors are:

| phase | primary-five mean | side | T3 | T2 | axial inversion |
|---|---:|---:|---:|---:|---:|
| heating training | 101.95 K | 130.87 K | 105.62 K | 11.53 K | 50.03 K |
| heating held out | 76.30 K | 93.96 K | 89.97 K | 9.66 K | 36.35 K |
| cooling training | 44.42 K | 45.56 K | 70.61 K | 14.81 K | 0.42 K |
| C81 held out | 37.26 K | 46.02 K | 39.20 K | 9.02 K | 2.53 K |

M-to-F primary-history RMS is 10.51, 32.17, 26.30, and 16.39 K for heating
E67/E71/E75/E80, so the fitted heating regime fails the 10 K mesh gate.
Cooling C69/C81 gives 9.22/4.68 K and passes. The independent exchange-energy
audit remains below `1.13e-14` relative and the one-case solver-tolerance
change is only 0.024 K RMS; the heating mesh failure is not a time-integrator
artifact.

The complete Stage-C-null branch changes representative heating primary
histories by 14.1--17.9 K RMS and T3 by 31.1--39.5 K RMS. Hence absolute v19
results are non-unique with respect to the rejected outlet-observation
closure. T3 save-point convergence itself is acceptable: final C81 T3 spans
0.027 K and `t90` spans 50 s from 61 to the available 119 record points.

Seventeen of eighteen cases satisfy the internal coupled-flow iteration.
E76 reaches the 24-iteration cap at 112/119 saved points, with maximum
pressure/gas-reference residuals `2.29e-5`/`1.19e-6`. The pressure residual is
small but exceeds the internal target near `1e-5`, so the strict convergence
gate is recorded as failed.

The final output contains 18 transient figures, 15 axial figures, and one
parity figure. Axial plots display the continuous modeled side-skin profile;
measured T8/T12/T11 remain discrete points at 11/58/107 mm.

For all future axial comparison sets, x- and y-axis limits must be fixed
globally and reused in every panel/experiment/version. Per-case autoscaling is
not permitted because it visually changes the apparent magnitude and shape of
profile differences.
