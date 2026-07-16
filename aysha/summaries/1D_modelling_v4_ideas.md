
# First Ideas
The pattern in the two plots indicates a **structural model mismatch**, rather than simply a poor value of $h_{\text{ref}}$ or the flow exponent.

For T8, the model reproduces the temperature range and approximate trend. For T3, the experimental temperature range is strongly compressed: low model temperatures are generally too low and high model temperatures are too high. That is characteristic of a measured temperature that is being buffered by downstream thermal mass, radiation, mixing, or conduction—not merely a delayed gas-temperature signal.

The present model assumes quasi-steady plug flow, a spatially uniform heat-transfer correlation scaled by flow, and only a first-order lag for T3. It also calibrates only T8–T10 and T3 rather than all five solid thermocouples.

## 1. Why the present model necessarily predicts strong flow dependence

Approximately, the gas equation gives

T_g,out ≈ T_s,eff − (T_s,eff − T_in) exp(−NTU)

with

NTU = UA/(ṁ c_p)

and, in the current model,

h ∝ qⁿ  
ṁ ∝ q

Therefore,

NTU ∝ qⁿ⁻¹

For the usual case $n<1$, increasing flow decreases NTU and drives the outlet gas temperature toward the inlet temperature. Thus, the current model is mathematically forced to predict a significant outlet-temperature dependence on flow.

Changing $h_{\text{ref}}$ or $n$ can change the magnitude of this dependence, but it cannot easily produce a nearly flow-independent T3 while simultaneously preserving the correct heat transfer to the solid.

The T3 time constant cannot solve this steady-state problem either. A first-order lag changes the transient response, but at steady state:

T₃,model = T_g,out

independently of the value of τ₃.

## 2. First modification: make the local heat transfer axial-dependent

The fact that T8 is reasonable but downstream solid temperatures are too flow-sensitive strongly suggests that the model applies the influence of flow too uniformly along the receiver.

The current form is essentially

h_i = h_ref (q/q_ref)ⁿ · k_g(T_i)/k_g(T_ref)

with no explicit dependence on axial position.

For developing internal flow, the heat-transfer coefficient is normally largest near the inlet and decreases toward a fully developed value. In a laminar fully developed channel, the Nusselt number—and therefore $h$, apart from property changes—can become almost independent of velocity.

A physically preferable formulation is

h(z,q,T) = k_g(T)/D_h · Nu[Re(q), Pr(T), z/D_h]

A reduced empirical version that would be easier to calibrate is

h(z,q,T) = h_fd(T)  
+ [h_ent(q,T) − h_fd(T)] exp(−z/L_h)

where:

- $h_{\text{ent}}$ represents the flow-dependent entrance-region coefficient;
    
- $h_{\text{fd}}$ represents the weaker-flow-dependent fully developed value;
    
- $L_h$ is an entrance/development length.
    

This would allow:

- strong flow sensitivity near T8;
    
- progressively weaker sensitivity toward T9–T12;
    
- fewer compensating changes in axial conductivity or optical penetration.
    

This is the first change I would make to the receiver model itself.

## 3. Second modification: represent the downstream plenum and T3 measurement

T3 should probably not be modelled as a delayed copy of the channel outlet temperature.

A minimal downstream node could be introduced:

C_p dT_p/dt =  
ṁ c_p (T_g,out − T_p)

- G_sp (T_s,rear − T_p)  
    − G_pa (T_p − T_amb)
    

where $T_p$ is an effective outlet-plenum temperature.

The measured T3 can initially be identified with $T_p$, or a small thermocouple response can be added later:

C_TC dT₃/dt =  
G_g $T_p − T₃$

- G_r (T_rad − T₃)
    
- G_c (T_mount − T₃)
    

This represents three effects that can compress the experimental temperature range:

- convective coupling to the gas;
    
- radiation from the hot rear receiver;
    
- conduction through the thermocouple sheath or mounting.
    

At steady state, the thermocouple behaves approximately as a conductance-weighted temperature:

T₃ ≈  
(G_g T_p + G_r T_rad + G_c T_mount) /  
$G_g + G_r + G_c$

This is preferable to assigning arbitrary statistical weights because the weights arise from identifiable heat-transfer conductances.

Because the cooling data show hysteresis, the plenum/hardware thermal capacitance is particularly relevant. A simple sensor lag alone cannot reproduce hysteresis caused by heat exchange with hot downstream hardware.

## 4. Use T8–T12 to identify where the flow sensitivity becomes incorrect

All five solid thermocouples should be included at their actual axial positions. For each thermocouple, calculate the experimental and modelled steady-state flow sensitivity:

S_i = ∂T_i/∂ln(q)

In practice, estimate it from regressions within groups having similar irradiance.

Plot $S_i$ against axial position. This will distinguish several situations:

- **Experimental sensitivity decreases progressively downstream:** use an axial developing-flow $h(z,q)$.
    
- **Sensitivity changes abruptly near the rear:** add rear hardware/plenum coupling or revise the rear boundary condition.
    
- **All downstream thermocouples move together:** a second solid domain or mounting structure may be required.
    
- **Only T3 has weak sensitivity:** the main problem is the gas-temperature measurement model rather than the receiver core.
    

The axial profiles themselves should also be compared:

T_i − T_in versus z

rather than only separate model-versus-experiment scatterplots.

## 5. Check whether T3 is actually a bulk gas temperature

For every steady experiment, calculate

Q̇_gas,exp = ṁ c_p(T₃ − T_in)

and compare it with the available absorbed solar power.

Warning signs are:

- apparent gas heat removal larger than plausible absorbed power;
    
- large gas heat removal during periods when the solid energy balance cannot support it;
    
- T3 correlating more strongly with T10–T12 than with flow;
    
- different T3 values during heating and cooling at approximately the same irradiance, flow, and solid temperatures.
    

These would indicate that T3 is influenced by radiation, mounting conduction, or downstream stored energy.

Also verify whether the flow measurement is in actual L/min, normal L/min, or standard L/min. The current model computes mass flow using inlet-density times the reported volumetric flow. If the controller reports standard flow, mass flow should instead be based on the specified standard density.

## 6. Recommended calibration sequence

I would avoid fitting the complete model again in one stage.

### Stage A — receiver solid and local heat transfer

Use the cooling experiments and T8–T12, temporarily excluding T3.

Fit only:

- effective thermal capacity;
    
- axial conductivity;
    
- side and rear losses;
    
- the axial $h(z,q)$ parameters.
    

This determines whether the model can reproduce the measured axial-temperature response without using the gas sensor to distort the receiver parameters.

### Stage B — downstream gas/plenum model

Freeze the main solid parameters and fit:

- plenum thermal capacitance;
    
- rear-solid-to-plenum conductance;
    
- plenum-to-ambient conductance;
    
- at most one thermocouple response parameter.
    

Avoid fitting both a highly flexible plenum model and an independent large sensor lag initially, because they will be difficult to identify separately.

### Stage C — heating and optical parameters

Only after cooling and flow dependence are satisfactory, fit:

- absorbed fraction;
    
- optical penetration/deposition;
    
- front convection;
    
- front radiation.
    

Otherwise, $\eta_{\text{abs}}$ and $\beta_{\text{opt}}$ will compensate for errors in convective heat transfer.

### Stage D — cross-validation

Fit using two cooling experiments and predict the third. Then rotate the omitted experiment. The same should be done for the heating cases. This is more informative than obtaining a lower global objective from all runs simultaneously.

## 7. Additional mechanisms to consider only if needed

If the axial $h$ model and plenum node remain insufficient, the next candidates are:

**Flow maldistribution or bypass.** The nominal flow may not all pass through thermally active channels. A two-path model can use a heated branch and a bypass branch, followed by outlet mixing.

**Internal radiative penetration.** The Beer–Lambert source may be too restrictive. A two-component deposition profile could represent front absorption plus deeper channel-wall irradiation.

**Second solid domain.** The receiver core and external shell/mount may have different temperatures and thermal time scales. A shell node coupled to the axial core could explain weak downstream flow sensitivity and cooling hysteresis.

These should not be introduced simultaneously. Each additional mechanism should be justified by a specific residual pattern.

## Recommended next version

The most defensible `1D_v4` would retain the existing axial solid finite-volume model but make three targeted changes:

1. Replace the uniform (h(q)) with a developing-to-fully-developed $h(z,q,T)$.
    
2. Add T11 and T12 explicitly to the objective and initial/profile reconstruction.
    
3. Replace the T3 pure lag with one effective downstream plenum/hardware energy-balance node.
    

That is still a compact reduced-order model, but it directly addresses both observed deficiencies: excessive downstream flow sensitivity and the compressed/hysteretic gas-temperature response.
