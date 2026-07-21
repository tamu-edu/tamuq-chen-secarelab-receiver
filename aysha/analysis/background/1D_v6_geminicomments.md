Viewed 1D_v6.jl:1-497

I have thoroughly reviewed the code in `1D_v6.jl` for theoretical consistency, numerical robustness, and proper parameterizations. Overall, the physical basis (energy conservation, convective limits) is very well constructed, but there are a few critical theoretical conflicts and numerical fragility issues that should be addressed for a strict "research-grade" codebase.

Here is the breakdown of the findings:

### 1. Misplaced Loss Regularization (Critical Theoretical Conflict)
There is a major logical mismatch in the objective functions regarding parameters `p[6]` (`h_floor`) and `p[7]` (`L_h`).
- **The Issue:** In `loss_cooling_v6`, you apply a regularization penalty to `p[6]` and `p[7]`. However, during the cooling stage, the optimization subset explicitly excludes them (`fit_cooling_thermal_indices_v6 = [1, 2, 3, 4, 8, 9, 10, 11]`). Therefore, their values never change during this stage, making their regularization a meaningless constant offset. 
- **The Conflict:** Conversely, `p[6]` and `p[7]` are explicitly optimized during the heating stage (`fit_heat_transfer_indices_v6 = [5, 6, 7]`), but `loss_heating_v6` **does not** contain any regularization terms. This means these parameters are entirely unconstrained when they are being fitted, but penalized when they are locked. 
- **Recommendation:** Move the `p[6]` and `p[7]` regularization terms out of `loss_cooling_v6` and into `loss_heating_v6`.

### 2. Fragile Optimization Loop & Error Handling
In `simulate_v6`, the ODE solver's return code is strictly enforced:
```julia
successful_retcode(solution) || error("1D_v6 ODE solve failed with retcode $(solution.retcode)")
```
- **The Issue:** Throwing a hard `error()` completely crashes the Julia script. During Nelder-Mead optimization, the simplex algorithms routinely step into "bad" parameter spaces (e.g., severe stiffness causing `dt <= dtmin`). 
- **The Conflict:** If `simulate_v6` throws an exception, it bypasses your safeguard (`all(isfinite, model) || return Inf`) in `loss_cases_v6` and halts the entire calibration run.
- **Recommendation:** Instead of throwing an error, `simulate_v6` should return `nothing` or an array of `NaN`s upon solver failure. This allows `loss_cases_v6` to catch the failure and gracefully return `Inf`, guiding the optimizer away from unstable parameter regions without crashing.

### 3. Boundary Condition Spatial Discretization (Grid Dependency)
The finite volume discretization defines nodes at cell centers (`z = dx/2, 3dx/2, ...`).
- **The Issue:** The front face convective/radiative losses and the rear direct losses are evaluated using the cell center temperatures (`Ts[1]` and `Ts[end]`), rather than the actual boundary surface temperatures.
```julia
Qfront = H_FRONT_FIXED_V6 * A_frt * (Ts[1] - ambient) + ...
```
- **The Conflict:** Because the true front surface temperature will be notably higher than `Ts[1]` under intense solar flux, using `Ts[1]` mathematically treats the outer `dx/2` of the receiver as having infinite thermal conductivity. This creates an artificial mesh dependency; if you change `nodes`, the distance to `Ts[1]` changes, which inherently alters the model's calibrated parameters.
- **Recommendation:** For a lumped model, this is an acceptable simplification as long as `nodes` remains strictly constant between calibration and validation. For a PDE-strict research model, calculate a boundary surface temperature `T_surf` via a local energy balance (Fourier conduction = external convection + radiation).

### 4. Zero-Flow "Stagnant Gas" Assumptions
In `gas_profile_v6!`, when `flow <= 1e-12`, the internal gas profile is instantly flushed to the inlet temperature (`Tin`):
```julia
fill!(Tg, Tin)
```
- **The Issue:** Physically, if the flow stops (e.g., during a forced cooling test), the trapped stagnant gas would equilibrate to the hot solid temperature `Ts`.
- **The Conflict:** By instantaneously dropping `Tg` to ambient, the predicted gas outlet temperature `Tf_true` instantly flatlines. Consequently, the entire decay of the downstream sensor `T3` is forced onto the numerical lag parameter `p[8]` (`tau_T3`). This means `p[8]` stops representing just the sensor lag, and must artificially absorb the missing thermal inertia of the stagnant hot gas and natural convection physics. 

### 5. Global Variable State Mutation
The function `calibrate_v6()` mutates global variables at its conclusion:
```julia
global pnew_v6 = parameters
```
- **The Issue:** Because the function initializes `p0 = copy(pnew_v6)`, running the calibration script twice in the same REPL session will yield entirely different trajectories (the second run starts from the local minima found by the first). 
- **Recommendation:** While standard for iterative REPL workflows, research-grade code usually isolates state. It's safer to have the function return the optimized parameters without altering global state, enforcing reproducibility for batch processing.

### Positive Theoretical Highlights (What works well)
- **Robust Effectiveness Limit:** Your handling of the `UA / (mdot * cp)` convection via `-expm1(...)` is excellent. It perfectly captures the low-flow asymptotic limit (where exit gas approaches solid temperature) without dividing by zero or risking explicit ODE instabilities.
- **Energy Conservation:** The solid-to-gas energy exchange (`Qgas = mdot * cp * (Tg[i+1] - Tg[i])`) and the subtraction of exactly `Qgas[i]` from the solid states guarantees strict energy conservation. The structural rear thermal mass coupling is also mathematically robust and conservative.
- **Heat Capacity Formulation:** Rewriting $C(T) \frac{dT}{dt} = Q$ into $\frac{dT}{dt} = \frac{Q}{C(T)}$ inside the derivative function is clean and avoids mass-matrix inversion overheads for an explicit representation.

## Codex Review and Revision Disposition

Implemented in both `1D_v6.jl` and `1D_v7.jl`:

1. The heat-exchange-shape regularization is now applied in the heating objective, where `h_floor` and `L_h` are actually fitted. In v6, the constant cooling-stage penalty on `p[6]` and `p[7]` was removed. In v7, the heating objective now regularizes the corresponding `p[5]` and `p[6]`.
2. The calibration loss loops now catch solver failures and return `Inf`, while still rethrowing `InterruptException`. This keeps bad simplex trial points from aborting the full calibration run.
3. The zero-flow gas diagnostic now sets the internal outlet gas profile to the local solid temperatures while keeping `Qgas = 0`. This avoids forcing the outlet gas signal to the inlet temperature during stagnant-flow/cooling cases.

Deferred:

1. The boundary-face temperature correction is physically valid as a mesh-sensitivity concern, but it is not a small patch because it changes how front absorption, conduction, and radiative loss are coupled at the first finite-volume cell. It should be handled as a dedicated mesh/boundary-condition revision.
2. The global mutation of `pnew_v6`/`pnew_v7` is retained for the current notebook/REPL workflow. Batch reproducibility can be improved later by adding a no-mutation calibration wrapper while preserving the interactive behavior.
