# 1D Receiver Model Theory Manual

This manual documents the theory behind `1D_v3.jl`, a conservative one-dimensional finite-volume model for the solid/gas solar receiver. The model is intended to sit between the very compact 0D lumped receiver models and a full multidimensional heat-transfer simulation.

The main goal is to predict transient axial solid temperatures and outlet gas temperature while keeping the model fast enough for parameter calibration against the experimental runs.

## 1. Physical Picture

The receiver is represented as a straight porous/channelized solid block of length `L`. Solar power enters at the front face, air flows axially through the channels, and heat is exchanged between the hot solid and the gas along the length.

The coordinate system is:

```text
z = 0      irradiated air inlet / front face
z = L      gas outlet / rear face
```

The model resolves only the axial coordinate. Radial and channel-scale details are collapsed into effective geometric and heat-transfer parameters.

The solid is dynamic. The gas is quasi-steady at each time step:

- The solid stores heat and evolves according to an ODE system.
- The gas temperature profile is marched from inlet to outlet using the current solid temperature field.
- The gas has no independent time derivative in this model.

## 2. Main Assumptions

The `1D_v3.jl` model uses the following assumptions:

1. The solid temperature varies only along `z`.
2. Gas flow is one-dimensional plug flow through the channels.
3. Gas thermal capacitance is negligible compared with solid thermal capacitance.
4. Solid axial conduction is conservative between neighboring finite-volume cells.
5. Solar absorption is distributed axially with a Beer-Lambert profile.
6. External side, front, and rear losses are represented by effective fitted coefficients.
7. The outlet thermocouple or downstream gas measurement may lag the true model outlet gas temperature.
8. Temperature-dependent air and solid properties are evaluated from polynomial correlations and clamped to a physically reasonable range.

These assumptions make the model fast and stable while preserving the dominant physics needed for calibration.

## 3. Geometry

The receiver geometry in `1D_v3.jl` is defined by:

| Symbol | Code Name | Meaning |
|---|---|---|
| `L` | `L` | receiver length |
| `w_t` | `w_t` | frontal receiver width |
| `A_frt` | `A_frt` | frontal aperture area |
| `w_chnl` | `w_chnl` | square channel width |
| `n_chnl` | `n_chnl` | number of channels |
| `A_flow` | `A_flow` | total open flow area |
| `A_solid` | `A_solid` | solid cross-sectional area |
| `P_exchange` | `P_exchange` | total gas/solid exchange perimeter per axial length |
| `Dh` | `Dh` | hydraulic diameter |

The current model uses:

```julia
A_frt = w_t^2
A_flow = n_chnl * w_chnl^2
A_solid = A_frt - A_flow
P_exchange = 4.0 * w_chnl * n_chnl
Dh = w_chnl
```

The finite-volume grid divides the length into `N` equal axial cells:

```text
dx = L / N
z_i = (i - 1/2) dx
```

Solid temperatures are stored at cell centers. Gas temperatures are stored at cell faces from inlet face `1` to outlet face `N + 1`.

## 4. State Variables

The dynamic state is the vector of solid temperatures:

```text
T_s(t) = [T_s,1(t), T_s,2(t), ..., T_s,N(t)]
```

Each solid cell represents a finite-volume slab:

```text
V_cell = A_solid dx
```

The thermal capacitance of cell `i` is:

```text
C_i = rho_s gamma_C Cp_s(T_s,i) A_solid dx
```

where:

- `rho_s` is solid density.
- `Cp_s` is temperature-dependent solid heat capacity.
- `gamma_C = p[1]` is a fitted heat-capacity multiplier.

## 5. Parameter Vector

The full calibrated parameter vector has 11 entries:

| Index | Name | Code Meaning |
|---:|---|---|
| `p[1]` | `gamma_C` | effective solid heat-capacity multiplier |
| `p[2]` | `k_scale` | effective axial-conductivity multiplier |
| `p[3]` | `U_side` | side-loss conductance per receiver length |
| `p[4]` | `U_rear` | rear/adaptor conductance |
| `p[5]` | `h_ref` | reference gas/solid heat-transfer coefficient |
| `p[6]` | `n_exp` | flow exponent for heat transfer |
| `p[7]` | `tau_T3` | outlet gas sensor/downstream lag |
| `p[8]` | `eta_abs` | absorbed solar fraction |
| `p[9]` | `beta_opt` | Beer-Lambert extinction coefficient |
| `p[10]` | `h_front` | front convection coefficient |
| `p[11]` | `eps_front` | effective front emissivity |

The model separates calibration into two groups:

```text
p_cool = p[1:7]
p_heat = p[8:11]
```

Cooling data are used first to identify thermal capacity, conduction, gas heat transfer, and losses that are observable without irradiation. Heating data are then used to identify optical and front-loss parameters.

## 6. Material Properties

The model uses temperature-dependent correlations for:

- air density `rho_f_f(T)`
- air viscosity `mu_f_f(T)`
- air thermal conductivity `kf_f(T)`
- air heat capacity `cpf_f(T)`
- solid heat capacity `Cps_f(T)`
- solid thermal conductivity `ks_f(T)`

All property calls pass through:

```julia
property_temperature(T) = clamp(Float64(T), 250.0, 1600.0)
```

The clamp prevents the empirical polynomial correlations from being evaluated far outside their intended range during optimization or transient startup.

## 7. Gas Flow Rate

The gas mass flow rate is computed from the volumetric flow in liters per minute:

```text
m_dot = rho_f(T_in) q_lpm / 60000
```

The factor `60000` converts liters per minute to cubic meters per second.

In code:

```julia
m_dot(flow_lpm, Tin=Tamb) = rho_f_f(Tin) * max(0.0, flow_lpm) / 60000.0
```

The density is evaluated at the inlet temperature. This is a practical approximation that keeps the flow rate tied to the measured flow controller data.

## 8. Quasi-Steady Gas March

At each ODE evaluation, the gas profile is recomputed using the current solid temperatures.

For cell `i`, the gas enters with temperature `T_g,i` and exits with `T_g,i+1`. The local film temperature is:

```text
T_film,i = 0.5 (T_s,i + T_g,i)
```

The effective heat-transfer coefficient is:

```text
h_i = h_ref (q_lpm / q_ref)^n_exp [k_f(T_film,i) / k_f(T_ref)]
```

where:

- `h_ref = p[5]`
- `n_exp = p[6]`
- `q_ref = 10 L/min`
- `T_ref = 600 K`

The cell gas/solid conductance is:

```text
UA_i = h_i P_exchange dx
```

The exact plug-flow NTU update is:

```text
epsilon_i = 1 - exp[-UA_i / (m_dot Cp_f)]
T_g,i+1 = T_g,i + epsilon_i (T_s,i - T_g,i)
```

The heat gained by the gas in the cell is:

```text
Q_gas,i = m_dot Cp_f (T_g,i+1 - T_g,i)
```

The implementation uses:

```julia
effectiveness = -expm1(-UA / (mdot * cp))
```

This is numerically more accurate than `1 - exp(...)` when the NTU is small.

If the flow is zero or nearly zero, the gas profile is set to the inlet temperature and gas heat removal is set to zero.

## 9. Solar Absorption

The total absorbed solar power is:

```text
Q_solar = eta_abs Io A_frt
```

where:

- `eta_abs = p[8]`
- `Io` is the aperture irradiance
- `A_frt` is the receiver frontal area

Axial distribution follows a Beer-Lambert attenuation:

```text
w_i = exp[-beta z_{i-1/2}] - exp[-beta z_{i+1/2}]
```

The weights are normalized so:

```text
sum_i w_i = 1
```

The solar power deposited in cell `i` is:

```text
Q_solar,i = Q_solar w_i
```

If `beta` is essentially zero, the source is distributed uniformly.

Large `beta_opt` concentrates absorption near the front face. Small `beta_opt` distributes absorption deeper into the receiver.

## 10. Solid Axial Conduction

Conduction between neighboring cells is treated as a conservative interface flux.

For the interface between cells `i` and `i + 1`, the face conductivity is the harmonic mean:

```text
k_face = 2 k_i k_{i+1} / (k_i + k_{i+1})
```

where:

```text
k_i = k_scale ks(T_s,i)
```

The conductive heat flow from cell `i` to `i + 1` is:

```text
Q_cond,i->i+1 = k_face A_solid (T_s,i - T_s,i+1) / dx
```

The same flux is subtracted from one cell and added to the neighbor:

```text
dE_i     -= Q_cond
dE_i+1   += Q_cond
```

This guarantees internal conduction conserves energy exactly. Any energy imbalance must come from boundary fluxes, gas exchange, storage, or numerical solve tolerance, not from the discretization itself.

## 11. Heat Losses

The model includes three external loss paths.

### Side Loss

Each cell has an effective side loss to ambient:

```text
Q_side,i = U_side dx (T_s,i - T_amb)
```

where `U_side = p[3]`.

This term absorbs all distributed axial side losses into one fitted conductance per unit length.

### Front Loss

The front face has convection and radiation losses:

```text
Q_front =
    h_front A_frt (T_s,1 - T_amb)
  + eps_front sigma A_frt (T_s,1^4 - T_amb^4)
```

where:

- `h_front = p[10]`
- `eps_front = p[11]`
- `sigma` is the Stefan-Boltzmann constant

The front loss acts only on the first axial cell.

### Rear Loss

The rear loss is represented as:

```text
Q_rear = U_rear (T_s,N - T_amb)
```

where `U_rear = p[4]`.

This term represents heat leakage into rear mounting, adaptors, outlet hardware, and other unmodeled downstream paths.

## 12. Solid Energy Balance

Before division by heat capacity, each cell accumulates a net heat rate:

```text
dE_i/dt =
    Q_solar,i
  - Q_gas,i
  - Q_side,i
  + Q_cond,in,i
  - Q_cond,out,i
  - boundary losses if applicable
```

Then:

```text
dT_s,i/dt = (dE_i/dt) / C_i
```

In vector form, the ODE is:

```text
dT_s/dt = f(T_s, t, p, operating)
```

This is solved with `DifferentialEquations.jl`. The default solver in `1D_v3.jl` is:

```julia
Rodas5P(autodiff=AutoFiniteDiff())
```

This stiff Rosenbrock method is a good default because radiative loss, gas NTU updates, and temperature-dependent properties can create stiffness during fast transients or calibration trials.

## 13. Outlet Gas Sensor Lag

The model outlet gas temperature is the true computed outlet:

```text
T_f,true(t) = T_g,N+1(t)
```

The measured thermocouple response can be delayed by downstream thermal mass or sensor lag. The model applies a first-order discrete lag:

```text
T_f,model[k] =
    T_f,model[k-1]
  + (1 - exp[-dt/tau_T3]) (T_f,true[k] - T_f,model[k-1])
```

where:

```text
tau_T3 = p[7]
```

If `tau_T3` is near zero, no lag is applied.

## 14. Initial Temperature Profile

For experimental cases, the initial axial solid profile is built from the first measured values of thermocouples:

```text
T8  at z = 0.011 m
T9  at z = 0.058 m
T10 at z = 0.107 m
```

The profile is piecewise linear:

- constant at `T8` before the T8 position
- linear from T8 to T9
- linear from T9 to T10
- constant at `T10` after the T10 position

This is more realistic than initializing the whole receiver from one average temperature, especially for cooling runs where axial gradients can be important.

## 15. Experimental Data Interface

The model imports data through:

```julia
include("import_exp_1D_v2.jl")
```

That import provides:

- `measurements`
- `measurements_cooling`
- `simulation_conditions`
- `simulation_conditions_cooling`
- `IDs`
- `IDs_cooling`

For each case, `solve_case_v3(...)` extracts:

- time vector
- measured T8, T9, T10, and gas outlet temperature
- measured flow
- measured or inferred inlet temperature
- measured ambient temperature
- nominal irradiance

These histories become interpolated operating-condition functions:

```julia
operating_condition_v3(
    irradiance = linear_history_v3(time, irradiance),
    flow_lpm = linear_history_v3(time, flow),
    inlet_temperature = linear_history_v3(time, Tin),
    ambient_temperature = linear_history_v3(time, ambient),
)
```

## 16. Model Outputs

The main simulation output from `simulate_v3(...)` is a named tuple containing:

| Field | Meaning |
|---|---|
| `time` | saved time vector |
| `z_solid` | solid cell-center positions |
| `z_gas` | gas face positions |
| `solid_temperature` | `N x Nt` solid temperature matrix |
| `gas_temperature` | `(N + 1) x Nt` gas temperature matrix |
| `heat_transfer_coefficient` | `N x Nt` local gas/solid heat-transfer coefficient |
| `ode_solution` | raw DifferentialEquations solution object |

The experimental comparison output from `extract_outputs_v3(...)` has five columns:

```text
[T8_model, T9_model, T10_model, T3_model, h_mean]
```

The first four columns are compared with the experiment. The fifth column is useful for diagnostics.

## 17. Calibration Objective

For each observed temperature signal, the model uses a normalized mean squared error:

```text
NMSE = mean(((model - experiment) / scale)^2)
```

where:

```text
scale = max(max(experiment) - min(experiment), 20 K)
```

The 20 K lower bound prevents almost-flat signals from becoming excessively influential.

For each case, the loss averages the four measured comparison channels:

```text
loss_case = mean(NMSE_T8, NMSE_T9, NMSE_T10, NMSE_T3)
```

The overall stage loss averages over the selected experimental cases.

## 18. Two-Stage Calibration

Calibration is staged:

### Stage 1: Cooling Calibration

Cooling runs have `Io = 0`, so they primarily identify:

- heat capacity multiplier
- axial conductivity multiplier
- side and rear losses
- gas heat transfer
- gas sensor lag

The fitted cooling vector is:

```text
p_cool = [gamma_C, k_scale, U_side, U_rear, h_ref, n_exp, tau_T3]
```

The current cooling objective also includes a weak regularization:

```text
0.005 [log(gamma_C)^2 + log(k_scale)^2]
```

This discourages non-unique fits from drifting too far from physically natural values.

### Stage 2: Heating Calibration

Heating runs use the calibrated cooling parameters and fit:

```text
p_heat = [eta_abs, beta_opt, h_front, eps_front]
```

These parameters are most visible when solar input is present.

The final parameter vector is:

```julia
pnew = assemble_parameters_v3(p_cool, p_heat)
```

After `calibrate_v3()` runs, the global `pnew` is updated and `calibration_v3` contains the full result.

## 19. Optimization Package Use

`1D_v3.jl` uses SciML optimization tools:

```julia
using Optimization
using OptimizationNLopt
using SciMLBase
```

The calibration helper builds:

```julia
optf = OptimizationFunction(counted_objective, SciMLBase.NoAD())
optprob = Optimization.OptimizationProblem(optf, x0, nothing; lb=lb, ub=ub)
solution = solve(optprob, NLopt.LN_NELDERMEAD(); ...)
```

`SciMLBase.NoAD()` is used because the objective includes ODE solves, data interpolation, and non-smooth guards such as clamps and `Inf` returns. Derivative-free optimization is the safer default.

The default algorithm is:

```julia
NLopt.LN_NELDERMEAD()
```

This is a local derivative-free optimizer. It is slower than gradient methods but tends to be robust for this kind of calibration problem.

## 20. Running the Model Manually

Load the model:

```julia
include("1D_v3.jl")
```

Run a single heating case:

```julia
model, result, experiment = solve_case_v3(pnew, "E74"; nodes=25)
```

Run a cooling case:

```julia
model, result, experiment = solve_case_v3(pnew, "C69"; is_cooling=true, nodes=25)
```

Run a quick calibration:

```julia
calibration_v3 = quick_calibration_v3()
```

Run the full calibration:

```julia
calibration_v3 = calibrate_v3()
```

After calibration:

```julia
pnew == calibration_v3.parameters
```

should be true element-by-element up to normal floating-point identity, because `calibrate_v3()` assigns:

```julia
global pnew = parameters
```

## 21. Running the End-to-End Script

The separate workflow runner is:

```text
run_1D_v3.jl
```

Run the full process:

```powershell
julia +1 --project=. run_1D_v3.jl
```

Run a short sanity process:

```powershell
$env:RECEIVER1D_V3_RUNNER_QUICK="true"
julia +1 --project=. run_1D_v3.jl
```

The runner:

1. includes `1D_v3.jl`
2. loads experimental data through the import script
3. runs calibration
4. saves calibrated parameters
5. computes metrics
6. creates steady-state comparison plots
7. creates transient comparison plots for all heating and cooling cases

Outputs are written under:

```text
summaries/1D_v3/
```

## 22. Plotting and Diagnostics

Transient comparison:

```julia
using StatsPlots
display(plot_case_v3("E74", pnew))
display(plot_case_v3("C69", pnew; is_cooling=true))
```

All heating cases:

```julia
for sm in sim_key_heat
    display(plot_case_v3(sm, pnew))
end
```

All cooling cases:

```julia
for sm in sim_key_cool
    display(plot_case_v3(sm, pnew; is_cooling=true))
end
```

Metrics:

```julia
metrics_v3 = compute_metrics_v3(pnew)
write_metrics_v3("summaries/1D_v3/analysis_results_1D_v3.csv", metrics_v3)
```

Each metric row contains:

- simulation ID
- phase
- sensor
- RMSE
- steady-state error
- t90 timing error

## 23. Energy Balance Check

`energy_rates_v3(...)` recomputes the instantaneous heat-rate terms for a solid profile:

```julia
rates = energy_rates_v3(result.solid_temperature[:, end],
                        result.time[end],
                        pnew,
                        operating)
```

It returns:

| Field | Meaning |
|---|---|
| `absorbed` | absorbed solar power |
| `gas` | heat transferred to gas |
| `side` | distributed side loss |
| `front` | front convection and radiation loss |
| `rear` | rear loss |
| `storage` | solid internal energy rate |
| `residual` | balance residual |

The residual is:

```text
absorbed - gas - side - front - rear - storage
```

For a correctly assembled instantaneous balance, this should be near numerical zero.

## 24. Interpretation of Fitted Parameters

The fitted parameters are effective parameters, not pure material constants.

Examples:

- `gamma_C` can absorb uncertainty in active solid volume, heat capacity, and unresolved thermal mass.
- `k_scale` can absorb axial conduction through the solid, contact paths, and geometry simplification.
- `U_side` and `U_rear` represent complicated environmental and fixture losses.
- `h_ref` is an effective gas/solid exchange coefficient over many channels, not necessarily a direct textbook channel HTC.
- `eta_abs` and `beta_opt` represent optical absorption and axial deposition together.
- `tau_T3` represents sensor response and downstream mixing, not just thermocouple bead dynamics.

This is normal for reduced-order calibrated receiver models. The important requirement is that the fitted values remain physically plausible and predictive across cases.

## 25. Common Failure Modes

### Calibration Gets Stuck at Bounds

If a parameter lands exactly on its lower or upper bound, the model may be missing a loss path, the data may be inconsistent, or the parameter may be weakly identifiable.

Useful checks:

- inspect cooling and heating stage losses separately
- plot cooling cases before heating cases
- reduce the parameter set temporarily
- check whether flow, ambient, and inlet histories are realistic

### Outlet Gas Fits but Solid Sensors Do Not

This can indicate:

- insufficient axial resolution
- incorrect solar deposition depth
- front/rear loss imbalance
- sensor position mismatch
- an overly flexible gas heat-transfer parameter

### Solid Sensors Fit but Gas Does Not

This can indicate:

- `h_ref` or `n_exp` mismatch
- flow-rate interpretation error
- gas sensor lag error
- inlet temperature bias

### Heating Fits but Cooling Does Not

Cooling is the cleaner test for thermal mass, conduction, and non-solar losses. If cooling fails, heating fits may be compensating with optical parameters.

### Cooling Fits but Heating Does Not

This points toward optical deposition, front loss, radiation, or irradiance assumptions.

## 26. Recommended Workflow

The practical workflow is:

1. Run the smoke test to confirm the model loads and one case solves.
2. Run `quick_calibration_v3()` to catch obvious calibration API or data problems.
3. Plot one heating and one cooling case.
4. Run full calibration through `run_1D_v3.jl`.
5. Inspect steady-state comparison first.
6. Inspect transient plots for all cases.
7. Check metrics and energy residuals.
8. Only then adjust parameter bounds or model structure.

This keeps debugging grounded in progressively more expensive checks.

## 27. Relationship to 0D Models

The 0D models compress the receiver to one or two dynamic states. They are useful for fast fitting and qualitative trends.

The 1D model adds:

- axial temperature gradients
- sensor-position-specific predictions
- axial solar deposition
- conservative axial conduction
- local gas heating along the receiver length

The price is more computation and a larger parameter-identifiability problem. This is why the staged cooling/heating calibration is important.

## 28. Relationship to Higher-Fidelity Models

Compared with 2D/3D or CFD-style models, this 1D model omits:

- radial temperature gradients
- individual channel flow variation
- detailed radiation exchange inside pores/channels
- detailed mounting geometry
- gas axial conduction
- gas thermal inertia
- local turbulence and entrance effects

Those effects are folded into effective fitted coefficients. The model is therefore best viewed as a calibrated reduced-order model, not a direct substitute for local high-fidelity heat-transfer analysis.

## 29. File Map

Relevant files:

| File | Role |
|---|---|
| `1D_v3.jl` | model, data interface, calibration functions, plotting helpers |
| `import_exp_1D_v2.jl` | experimental data import contract |
| `run_1D_v3.jl` | end-to-end calibration and plotting workflow |
| `test/smoke_1D_v3.jl` | loads model and solves one representative case |
| `test/test_1D_v3_static.jl` | source-structure regression checks |
| `summaries/1D_v3/` | runner output directory |

## 30. Minimal Example

```julia
include("1D_v3.jl")

# Calibrate.
calibration_v3 = quick_calibration_v3()

# Run a heating case with calibrated parameters.
model, result, experiment = solve_case_v3(pnew, "E74"; nodes=25)

# Plot comparison.
using StatsPlots
display(plot_case_v3("E74", pnew; nodes=25))

# Compute metrics.
metrics_v3 = compute_metrics_v3(pnew; heating_keys=["E74"], cooling_keys=String[])
```

For the full scripted process:

```powershell
julia +1 --project=. run_1D_v3.jl
```
