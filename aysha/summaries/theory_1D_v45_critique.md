## Overall assessment

`theory_1D_v45.md` is not manuscript-ready in its present form. I would recommend **major revision bordering on rejection as a physical coefficient-identification paper**.

The v45 model is a useful empirical reduced-order temperature model, and it improves several in-sample predictions—especially \(T_3\), \(T_2\), and the two-point \(T_9-T_{10}\) axial drop. However, the current document:

- does not faithfully describe the implementation;
- contains a serious gas-stream mass/enthalpy inconsistency;
- describes calibration data as validation data;
- attributes physical meaning to strongly confounded fitted parameters;
- omits unfavorable invariant and transient results;
- claims optical, heat-transfer, and capacitance validation that the saved evidence contradicts.

The scientifically defensible classification is therefore **Role-A diagnostic/surrogate model**, not a validated Role-B extractor of macroscopic \(Nu\), optical extinction, suction convection, or physical component properties.

No files were modified.

## 1. Critical mismatch between manuscript and implementation

The most immediate problem is reproducibility: the equations and constants in the theory document are not those used to generate the reported results.

| Quantity | Manuscript | Implemented v45 |
|---|---:|---:|
| Channel opening | 1.34 mm | 1.50 mm |
| Web thickness | 0.56 mm | 0.40 mm |
| Porosity | 0.497 | 0.623 |
| \(A_\text{solid}\) | \(1.815\times10^{-4}\,\mathrm{m^2}\) | \(1.36\times10^{-4}\,\mathrm{m^2}\) |
| \(P_\text{exchange}\) | 0.536 m | 0.600 m |
| SiC density | 3100 kg/m³ | 3200 kg/m³ |
| Felt/insulation outer radius | 25 mm | 57 mm |
| Housing outer radius | 30 mm | 75 mm |
| Housing thickness | 5 mm | 18 mm |
| Rear-tube length | 63 mm | 150 mm |
| Rear-tube gas radius | 7 mm | 6.5 mm |
| \(T_9\) position | 11 mm | 58 mm |

The manuscript geometry is stated at [theory_1D_v45.md:28](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/theory_1D_v45.md:28), whereas the inherited geometry actually used by v45 is at [1D_v3.jl:44](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/1D_v3.jl:44). The \(T_9\) location is especially serious: the manuscript calls it a “front core” sensor at 11 mm, but the code and experimental documentation place the interior \(T_9\) probe at 58 mm.

The SiC property equations also differ completely. At 900 K, for example:

- manuscript \(k_s\approx9\ \mathrm{W\,m^{-1}K^{-1}}\);
- implemented polynomial \(k_s\approx64\ \mathrm{W\,m^{-1}K^{-1}}\).

The implemented functions are at [1D_v3.jl:115](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/1D_v3.jl:115), not the formulations printed at [theory_1D_v45.md:43](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/theory_1D_v45.md:43).

Other material discrepancies include:

- manuscript cavity emissivity 0.80 versus code 0.20;
- manuscript \(T_3\) observer \(h=80\ \mathrm{W\,m^{-2}K^{-1}}\), \(A=10^{-5}\,\mathrm{m^2}\), \(\epsilon=0.85\);
- code \(h=150\), \(A=2.83\times10^{-5}\,\mathrm{m^2}\), \(\epsilon=0.80\), at [1D_v45.jl:109](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/1D_v45.jl:109);
- manuscript quadratic rear-contact weights, while the code uses a normalized affine profile \(0.02+0.98s\), at [1D_v45.jl:355](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/1D_v45.jl:355).

Consequently, the printed coefficients and results cannot be reproduced from the printed formulation.

## 2. Gas-flow formulation violates enthalpy conservation

This is the most consequential technical defect.

Inside the monolith, only

\[
\dot m_\text{active}=\phi_\text{act}\dot m_\text{total}
\]

is heated. In the saved cases, \(\phi_\text{act}\) ranges from approximately 0.10 to 0.53. But at the monolith/rear-tube interface, the code passes the active-stream temperature directly into a rear-tube calculation using the **total** mass flow:

- active flow through the receiver: [1D_v45.jl:517](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/1D_v45.jl:517);
- rear tube using total flow at the active-stream temperature: [1D_v45.jl:546](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/1D_v45.jl:546).

If the inactive fraction remains at \(T_\text{in}\), a conservative mixing relation should satisfy approximately

\[
\dot m_\text{total}h(T_\text{mix})
=
\phi\dot m_\text{total}h(T_\text{active})
+
(1-\phi)\dot m_\text{total}h(T_\text{in}).
\]

Instead, the code effectively assigns

\[
\dot m_\text{total}h(T_\text{active})
\]

to the rear gas. The unaccounted enthalpy increase is therefore

\[
(1-\phi)\dot m_\text{total}
\left[h(T_\text{active})-h(T_\text{in})\right].
\]

This is exactly the modification credited with resolving the \(T_3\) deficit. It is not merely “eliminating non-physical dilution,” as claimed at [theory_1D_v45.md:134](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/theory_1D_v45.md:134); it removes a required mass/enthalpy mixing step.

Two physically consistent alternatives exist:

1. Treat \(T_3\) as sampling an unmixed active jet, in which case the rear-tube heat exchange must also use \(\dot m_\text{active}\); or
2. Treat the rear tube as carrying total flow, in which case active and bypass streams must mix conservatively first.

The current model combines the favorable part of both interpretations.

The meaning of \(\phi_\text{act}\) is also internally ambiguous. If it is the fraction of active channels, both active flow area and active perimeter should scale with \(\phi\). The code instead reduces mass flow and Reynolds number but retains the full flow area and full exchange perimeter. It is therefore an empirical heat-exchange modifier, not a rigorous active-channel fraction.

## 3. “Suction heat recovery” is implemented as an unreturned loss

The front-face suction term removes

\[
Q_\text{front}=h_\text{front}A_\text{frt}(T_\text{perim}-T_\text{amb})+\cdots
\]

from the solid at [1D_v45.jl:636](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/1D_v45.jl:636), but this energy is not added to the incoming gas. The gas profile still begins at the measured upstream \(T_\text{in}\).

For an active suction receiver, air heated while crossing the front surface subsequently enters the channels. Calling this process “heat recovery,” as at [theory_1D_v45.md:327](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/theory_1D_v45.md:327), is inconsistent with the equations: numerically it is an external heat sink.

The fitted value \(h_\text{suction}=296.7\ \mathrm{W\,m^{-2}K^{-1}}\) therefore cannot presently be interpreted as a physical convective coefficient. It is a strong flow-dependent correction applied to the front perimeter temperature.

Additional concerns are that:

- internal-channel \(D_h\) and \(A_\text{flow}\) are used for an external aperture boundary layer;
- the term acts over the entire frontal area rather than a defined solid or stagnation area;
- no correlation, characteristic length, or uncertainty is supplied;
- the parameter is strongly confounded with \(\phi_0,m_\text{rec},A_{Nu},B_{Re},C_z,\delta_\text{web}\), and the source distribution.

## 4. The calibration and “validation” descriptions are inaccurate

The manuscript states that the objective simultaneously fits 15 heating and three cooling runs using an integrated squared-error objective with monotonicity penalties at [theory_1D_v45.md:202](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/theory_1D_v45.md:202). The implementation differs materially:

- The runner defaults to **heating-only calibration**, at [run_1D_v45.jl:25](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/run_1D_v45.jl:25).
- The saved summary does not record an override of that dataset.
- The actual loss is a range-normalized MSE plus slope and \(t_{90}\) penalties, not the printed time integral: [1D_v45.jl:1055](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/1D_v45.jl:1055).
- The “ordering” term does not enforce \(T_9>T_{10}\) and \(T_8>T_{11}\). It matches four final temperature offsets to their experimental values.
- The cooling curves are generated after fitting, but the provenance does not demonstrate that they were included in the fitted objective.
- The saved optimization ended with `return_code=MaxTime`, not demonstrated convergence: [optimization_summary_1D_v45.txt:1](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/1D_v45/optimization_summary_1D_v45.txt:1).

All v42–v45 summaries ended at `MaxTime`. The sequence \(0.227\to0.0864\) shows improvement in the best recorded in-sample composite loss, but not “steady and dramatic convergence” to identified minima.

Calling the 60.8% change an “error reduction” is also misleading. It is a reduction in a mixed, dimensionless objective. Because that objective includes MSE, slope, timing, and offset terms, it cannot be converted directly into a percentage reduction in physical temperature error.

There is no:

- held-out irradiance or flow group;
- cross-validation;
- independent experimental validation set;
- multistart convergence study;
- profile likelihood or posterior uncertainty;
- bootstrap interval;
- mesh-convergence study.

The fit is performed on 15 monolith cells, while plots and reported metrics use the inherited 25-cell default. With \(Nu\propto z^{-0.83}\), first-cell exchange is mesh-sensitive, so this difference requires explicit testing.

## 5. Optical interpretation is unsupported

The document calls \(M=1.34/1.58/1.11\) “high-fidelity optical ray-tracing priors.” The project audit explicitly states that these are provisional R6 \(K_\text{heat}\) closures derived from the same experimental campaign and dependent on the unresolved \(T_3\) relation—not independent optical measurements. See [PROJECT_STATUS_2026-08-28.md:140](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/PROJECT_STATUS_2026-08-28.md:140).

Likewise, \(f_\text{dep}\) and \(\beta_\text{opt}\) were fitted solely from the temperature histories. Without an independent axial absorption measurement, they are effective source-shape parameters, not validated optical properties.

The physical description of the fitted source is overstated:

- \(\beta_\text{opt}=243\ \mathrm{m^{-1}}\) gives a 4.1 mm e-folding depth.
- More than 99.7% of its nominal penetrative component is absorbed within 25 mm.
- At the 25-node reporting mesh, approximately 82.2% of total core power is in the first 5.48 mm cell.
- At the 15-node fitting mesh, approximately 92.7% is in the first 9.13 mm cell.

Thus the fit is still strongly front-loaded. It does not heat the gas “continuously along the entire 137 mm” in any substantial sense.

The perimeter exponential is also inconsistent with the documented hardware interpretation that optical spill strikes the front felt/casing rather than penetrating axially through dense packaging, recorded at [journal.1D.md:5337](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/journal.1D.md:5337).

## 6. Saved “energy-accounting” outputs contain stale source semantics

The physical RHS now conservatively partitions

\[
Q_\text{core}+Q_\text{perim}=M I A_\text{frt}
\]

at [1D_v45.jl:609](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/1D_v45.jl:609). However, the output helper still uses the older flux-spillover formula at [1D_v45.jl:317](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/1D_v45.jl:317).

For E67:

- actual RHS delivered input: \(220.6\) W;
- saved `participating_absorbed_W`: \(1570.4\) W;
- saved `core_absorbed_W`: \(220.6\) W, even though the actual core share is only \(\chi Q_\text{delivered}\approx100.4\) W.

This is visible in [steady_results_fitted_energy_accounting_1D_v45.csv:2](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/1D_v45/steady_results_fitted_energy_accounting_1D_v45.csv:2).

Therefore, energy-input columns in the saved artifacts cannot be used for scientific interpretation until regenerated from the same source terms used by the ODE. A true conservation residual is also absent.

## 7. Heat-transfer interpretation is not rigorous

The claim that \(B_{Re}=0.5346\) “rigorously recovers” flat-plate boundary-layer physics is unjustified.

The closure is

\[
Nu=A Re^{B_{Re}}Pr^{1/3}(D_h/z)^{C_z},
\]

with independently fitted \(B_{Re}=0.535\) and \(C_z=0.830\). It is neither a conventional Graetz correlation nor a flat-plate similarity solution. It also has no fully developed asymptote because `NU_FLOOR=0` at [1D_v45.jl:66](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/1D_v45.jl:66). The predicted Nusselt number can consequently decay toward zero downstream, whereas a fully developed square duct has a finite asymptotic value.

The project-status analysis already notes that at the relevant Reynolds numbers the thermal entry length is only a few millimetres compared with a 137 mm receiver; a “long, entrance-dominated channel” interpretation is therefore doubtful. See [PROJECT_STATUS_2026-08-28.md:207](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/PROJECT_STATUS_2026-08-28.md:207).

Most importantly, the saved macroscopic invariants contradict the celebratory interpretation:

- apparent-\(Nu\) exponent: 2.642 versus target 1.44;
- apparent-\(Nu\) prefactor: \(2.08\times10^{-7}\) versus \(3.1\times10^{-4}\);
- \(\Lambda_{107}\) slope: \(4.22\times10^{-4}\) versus \(8.3\times10^{-4}\);
- mean effectiveness range: 0.280–0.677;
- participating capacity: 267.4 J/K, classified as failed.

These are reported in [invariant_summary_1D_v45.csv](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/1D_v45/invariant_summary_1D_v45.csv:1), but omitted from the theory manual.

The local fitted \(B_{Re}\) is therefore not the experimentally observed apparent exponent. Its value is conditional on—and confounded with—the active-flow law, axial decay, series resistance, source depth, and suction term.

## 8. Capacity “confirmation” is false

The manuscript states that

\[
C_\text{total}=267.4\ \mathrm{J/K}
\]

agrees with \(301\pm23\ \mathrm{J/K}\). The stated interval is 278–324 J/K, so 267.4 J/K is outside it. The saved invariant file correctly labels the result `failed`.

Moreover:

- \(C_{\text{perim}}=105.006\ \mathrm{J/K}\) is effectively at its lower bound of 105 J/K;
- the total excludes the modeled cavity capacity, which is approximately 4026 J/K from the implemented geometry;
- it also excludes the distributed rear-tube capacity;
- the printed geometry and properties would not produce the code’s 72.53 J/K core capacity.

The result is an effective subset of participating states, not a reconciliation of the “total assembly inventory.” It should be reported as such.

The fitted \(\delta_\text{web}=58\,\mu\mathrm m\) is similarly mislabelled as the “ceramic wall half-thickness.” Half-thickness is 200 µm using the implemented 0.4 mm web, or 280 µm using the manuscript’s 0.56 mm web. The fitted value can only be described as an equivalent resistance length.

## 9. Actual performance is mixed, not high fidelity

From the saved full-history metrics:

| Signal | Heating mean RMSE | Mean \(t_{90}\) error |
|---|---:|---:|
| \(T_8\) | 51.4 K | +513 s |
| \(T_{12}\) | 69.0 K | +720 s |
| \(T_{11}\) | 52.9 K | +333 s |
| \(T_9\) | 59.9 K | +710 s |
| \(T_{10}\) | 35.8 K | +240 s |
| \(T_3\) | 22.3 K | −140 s |
| \(T_2\) | 4.9 K | +61 s |

The worst heating steady errors exceed 95–103 K for several solid signals. The E67 plot shows the essential transient problem clearly: measured solid temperatures rise rapidly, while the model takes much longer to approach its final level.

The theory’s transient section omits \(T_{12}\) and \(T_{11}\), although both are calibration signals and have among the worst RMSE values. It also omits the cooling \(T_2\) timing result, whose mean is distorted by a worst-case error of 2250 s. This selective reporting is inconsistent with the heading “Comprehensive Validation.”

The claim that \(T_3\) reaches “precise alignment” of approximately ±20 K is also false: steady errors include −43 K for E67 and −63 K for E72.

The axial-drop result \(T_9-T_{10}\) is genuinely encouraging, with errors mostly within about ±18 K. But it is not independent spatial validation:

- both probes were used during calibration;
- there are only two interior axial points;
- the difference cancels their substantial common-mode biases.

Similarly, the model does not capture the \(T_{12}-T_8\) inversion transition for all runs. It gets the sign correct in 12 of 15 cases, failing E74, E79, and E80.

The flow-slope table reveals persistent structural error:

- at 456 kW/m², modeled \(T_{11}\) slope is −10.29 K/LPM versus measured −1.37;
- at 304 kW/m², −11.13 versus −1.93;
- modeled \(T_3\) slopes are −2.97 and −5.41 K/LPM where the measured slopes are approximately +0.54 and −0.13.

Thus deep receiver and outlet flow dependence remains incorrectly structured.

## 10. Observation-model limitations remain unresolved

The 2D journal’s central warning still applies: \(T_3\) is a real measurement, but it is not independent ground truth for bulk outlet enthalpy when the same measurement is also used to determine its downstream observation relation. See [journal.2D.md:3990](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/journal.2D.md:3990).

The archival location evidence favors approximately 136 mm; 137 and 140 mm remain discrete alternatives rather than metrologically established locations, as recorded at [journal.2D.md:4046](D:/kkakosim/github/tamuq-chen-secarelab-receiver/aysha/summaries/journal.2D.md:4046).

The \(T_2\) mapping is also physically coarse. Experimentally, \(T_2\) is a probe embedded roughly 40 mm into the insulation, whereas v45 compares it to a single, perfectly mixed felt-plus-metal cavity state with about 4 kJ/K heat capacity. Its low RMSE is useful, but does not validate that observation mapping.

Finally, the lamp-off flange multiplier is an explicitly phase-dependent correction:

\[
1+g_\text{cool}(1-e^{-t/\tau_\text{cool}}).
\]

It changes the hardware conductance only after irradiation vanishes. Unless supported by an actual temperature-dependent contact mechanism, it is a cooling-specific empirical branch. Under the runner’s heating-only default, \(g_\text{cool}\) and \(\tau_\text{cool}\) are effectively unidentifiable from the fitted objective.

## Defensible conclusions

The following statements are supported:

- The two-zone model can reproduce much of the steady temperature ordering.
- \(T_2\) is predicted accurately in absolute temperature.
- The in-sample \(T_3\) residual is substantially smaller than in v44.
- The \(T_9-T_{10}\) difference is reproduced well.
- The inversion sign is captured for 12 of 15 heating cases.
- The external-loss diagnostic remains near the measured \(K_\text{loss}\) band.
- Reparameterizing the RHS source as conserved \(M,\chi\) was conceptually correct.

The following are not supported:

- validated Nusselt or suction-convection coefficients;
- independently identified optical extinction;
- confirmation of physical wall penetration or active-channel recruitment;
- validated bulk outlet-gas temperature;
- convergence to a unique parameter basin;
- physical assembly-capacitance confirmation;
- comprehensive heating/cooling validation;
- a manuscript-ready, self-contained formulation.

## Recommended research sequence

1. Synchronize the document with the actual geometry, properties, sensor locations, observation constants, and boundary conditions.

2. Repair gas enthalpy conservation. Explicitly choose between active-jet observation and total-flow mixing, and verify the interface enthalpy to numerical tolerance.

3. Couple front suction heat into the inlet gas, or clearly redefine it as an external loss with a physically separate air stream.

4. Remove the stale source-accounting exporter and generate an actual source/storage/gas/loss residual for every run.

5. Record calibration dataset, mesh, code revision, parameter bounds, seed, and return code in every optimization artifact.

6. Replace the heating-only default with a declared training/validation design—ideally hold out an entire irradiance level or selected flow points.

7. Test a physically grounded internal-duct correlation with a finite fully developed \(Nu\) asymptote and perform 15/25/50-node convergence.

8. Treat \(M\), \(f_\text{dep}\), and \(\beta_\text{opt}\) as provisional effective source parameters until independent flux or absorption information exists.

9. Perform multistart/profile-likelihood analysis on the transport block. At minimum, profile \(A_{Nu},B_{Re},\phi_0,m_\text{rec},C_z,\delta_\text{web}\), and \(h_\text{suction}\).

10. Report all seven signal metrics, flow slopes, invariant failures, bound activity, and cooling results—not only the favorable axial-drop and \(T_3\) summaries.

In short: v45 is a promising empirical diagnostic, but its headline improvement is presently inseparable from a non-conservative gas-routing change and several highly flexible compensators. The coefficients should not yet be presented as physical transport properties.