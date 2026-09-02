# What this paper adds, against the literature in `analysis/literature`

Plain-language positioning note, 28 August 2026. Written from the abstracts and
introductions of all 21 PDFs in the folder. Raw material for the Introduction
and for the reply to a "what is new here?" referee.

---

## 1. The 21 papers fall into three families

**Family A — "which material or shape is best?" (experiments).**
Fend & Hoffschmidt 2004; Fend, Pitz-Paal et al. 2004; Avila-Marin, De Lara &
Fernandez-Reche 2018; Avila-Marin & Fernandez-Reche 2022; Zaversky et al. 2018;
Wang et al. 2021.

They build absorbers — ceramic foams, wire meshes, screen-printed SiC, diamond
and Voronoi unit cells, single- versus multi-layer, graded porosity — put them
under concentrated light, and report steady efficiency. The output is a
recommendation: use this porosity, this cell density, this layering.

**Family B — "what does the model predict?" (simulation).**
Kribus et al. 2014; Capuano, Fend et al. 2016; Fend et al. 2013; Avila-Marin
et al. 2019 (review of modelling strategies); Avila-Marin, Caliot et al. 2019
(homogeneous equivalent + P1); Smirnova & Fend 2011; Mendes et al. 2014;
Moro Filho & Malalasekera 2020; Salih & Kakosimos 2023.

They treat the absorber as a porous continuum with local thermal
non-equilibrium and a radiation model, and compute where the volumetric effect
should appear and what the best design would be. Local thermal non-equilibrium
here is an *assumption of the model*, not a measured quantity.

**Family C — "what is the right coefficient for a honeycomb channel?"
(catalysis heat transfer).**
Cornejo & Hayes 2020; Hayes et al. 2009; Hayes & Cornejo 2021; Cui & Kær 2018;
Schlereth & Hinrichsen 2014.

They compute channel-level Nusselt numbers, effective conductivities of the
honeycomb solid, and check when a 2D continuum reduction of a monolith
reproduces full 3D CFD. This is the most rigorous heat-transfer work in the
folder, and it is about the channel and the substrate, not about an installed
solar receiver.

## 2. The gap, in one sentence

Every paper in the folder answers either *which absorber to build* or *what a
model of the absorber says*. **None of them measures an installed receiver
while it heats up and cools down, and reads the governing numbers straight off
the measurement.** That is what this paper does.

## 3. What is new, stated simply

**(a) The volumetric effect is shown to be an operating condition, not only a
design property.**
The literature treats the volumetric effect as something you obtain by choosing
the right structure: graded porosity (Avila-Marin 2018), the right cell density
and porosity (Zaversky 2018), multi-layer stacking, specular reflectivity and
simplified geometry (Salih & Kakosimos 2023, this group's own prior work).
We show that on one plain, ungraded SiC honeycomb the hot spot moves from the
front face into the interior once the gas recovers about two thirds of the wall
excess temperature, and that this threshold does not move when the flux changes
by a factor of 1.8. The same receiver is volumetric or not depending on how it
is *run*. That is a design rule anyone can apply without a new material, and it
is a direct experimental complement to the group's 2023 modelling paper, which
sought volumetric behaviour through simplified design rather than complex
porosity grading.

**(b) The limiting resistance is measured, and it is not the channel.**
Solar-side correlations report Nu = a Re^b with exponents around 0.3–0.6, and
the catalysis side (Cornejo & Hayes 2020) gives well-founded channel Nusselt
numbers of order 3–4. Our receiver, reduced at the assembly scale, gives
Nu = 3.1e-4 Re^1.44 — 15 to 100 times *below* fully developed duct theory, with
an exponent above one. This is not a competing channel correlation. It is
evidence that in an installed receiver the bottleneck is how much of the solid
actually participates, and that participation grows with flow. It also explains
why correlations from this class of experiment scatter over two orders of
magnitude: they are assembly-scale quantities being reported as material
properties.

**(c) Local thermal non-equilibrium is measured rather than assumed.**
Families B and C build local thermal non-equilibrium into their equations.
We measure the gas–solid gap deep in the receiver and find it grows linearly
with Reynolds number, and that local thermal equilibrium never holds anywhere
in the tested range. The modelling literature can now be checked against a
number instead of a premise.

**(d) The thermal mass that matters is mostly not the absorber.**
Every paper in the folder models the absorber. From the cooling decay we
identify the effective heat capacity of the installed assembly as 301 J/K,
against 42–47 J/K for the 40 g monolith itself. Six to seven parts in seven of
the participating thermal mass are housing, insulation and mounting. Any
transient or start-up prediction based on absorber properties alone is wrong by
that factor. This is obtainable only from transient data, which is why it is
absent from a steady-state literature.

**(e) All fifteen heating transients collapse onto one curve.**
Under a single time scale built from the identified capacity and loss, the
installed receiver behaves as a one-parameter dynamic system. Nothing in the
folder reports transient dynamics of this class at all.

**(f) The delivered power is reported as a bound, and the reason is stated.**
The measured gas output exceeds the nominal aperture power by 23%. Rather than
absorb this into a calibration factor, we localize it by lamp configuration
with a model-free closure and report every efficiency on a stated basis. Small
apertures under concentrated beams always spill; the literature rarely says
what its efficiency denominators actually contain.

**(g) An identifiability boundary, which no paper in the folder draws.**
Avila-Marin et al. 2019 reviews modelling strategies and Hayes & Cornejo 2021
reviews thirty years of multi-scale monolith modelling; both ask how to model
better. We ask a different question: with wall thermocouples and one outlet gas
probe, what can and cannot be determined at all? Two things cannot. The
receiver conductance UA(Re) is not uniquely identifiable, and the source
magnitude cannot be separated from the captured spillage because only their
product enters the temperature field. We also show which continuum structures
are ruled out by the data. This tells the next experimenter what to instrument,
which is more useful than another fitted coefficient.

## 4. Why our 2D continuum failure is a result, not an embarrassment

Cui & Kær (2018) and Schlereth & Hinrichsen (2014) show that a 2D
pseudo-continuous reduction of a monolith reproduces full 3D CFD very well —
for catalytic reactors with a heated wall and uniform inflow. Our 2D reduction
of the same class of structure fails by 100–200 K while conserving energy to
one part in 10^10. The contrast is the point: what breaks is not the continuum
reduction itself but its use under a concentrated, heavily spilling beam with
an insulated housing and strongly maldistributed flow. Stated this way, the
negative result is a boundary on an established method rather than a failure of
ours.

## 5. What is honestly *not* new

The receiver concept, the material, the idea of the volumetric effect, and
local thermal non-equilibrium as a concept are all long established in this
folder. We add no new absorber and no new model. The contribution is
measurement and identification: a criterion, four constants with intervals, a
similarity collapse, and a statement of what the data cannot determine.

## 6. Suggested one-sentence claim for the Introduction

> Volumetric receivers have been studied largely by building better absorbers
> and by simulating them; this work instead measures an installed SiC honeycomb
> receiver through its transients and extracts, without a model in the loop, the
> operating-point criterion for the volumetric inversion, the assembly-scale
> exchange law, the degree of local thermal non-equilibrium, and the thermal
> mass that actually governs its dynamics.
