# Submission state — read this first

Written 2026-10-11, immediately after submission. It records what was submitted, what is still
unfixed, and which files are current, so the work can be picked up after a long gap.

## 1. Status

The manuscript was submitted on 2026-10-11 as `manuscript_v6.pdf` with `SupplementaryMaterial_v2.pdf`.
`cover_letter.md` is addressed to **Energy Conversion and Management** — confirm against the actual
submission, since the journal was still under discussion when the letter was drafted.

Submitted title:

> What depth-resolved transient measurements on a SiC honeycomb solar receiver imply for
> reduced-order heat transfer closures

Keywords: volumetric solar receiver; silicon carbide honeycomb; local thermal non-equilibrium;
number of transfer units; flow maldistribution; model validation.

Highlights prepared for the submission form (Elsevier caps each at 85 characters):

1. Eighteen transients in a SiC honeycomb receiver, resolved at nine locations
2. Transfer-unit count rises as Re^+0.341; fixed-position closures require Re^-1.00
3. Choice of solid-temperature probe moves the Nusselt prefactor by a factor of 29.7
4. Transient eigenvalue gives 275 +/- 38 J/K, six to seven times the bare monolith
5. Apparent Nusselt number lies one to two orders below the laminar duct value

## 2. Which files are current

| File | Status |
|---|---|
| `manuscript_revised_v3_endnote.md` | **Canonical.** EndNote `{Surname, Year #RecNum}` citations, submitted title and keywords |
| `manuscript_revised_v3_withfigures.md` | Same body text, numeric `[n]` citations and figure embeds, but **one edit behind**: it still carries the old title and the old keyword line |
| `supplementary_material_endnote.md` | Canonical supplementary, EndNote citations |
| `supplementary_material.md` | Numeric-citation counterpart |

A normalized diff of the two manuscript masters shows they differ in exactly two lines, the title
and the keywords; everything else is identical once citation format and figure embeds are set
aside. Syncing those two lines into `manuscript_revised_v3_withfigures.md` is the first thing to do
on returning, and it is left for the author, who applies title edits personally.

Ignore as superseded: `manuscript_revised_v3.md`, all `manuscript_revised_v2*.md`, `README.md` and
`README_reproduction.md` (both describe the retired Python pipeline), `receiver_reduction.py`,
`make_figures.py`. The live pipeline is Julia.

## 3. Outstanding corrections

`submission_check_v6.md` audits the two submitted **PDFs** and lists four must-fix items and about
twenty-five grammar and typo items. Checking each against the markdown masters on 2026-10-11
separates them into two groups, which matters because the fix belongs in a different place for
each.

Already correct in the markdown, so the defect was introduced in the PDF or Word conversion and
needs no edit here:

- **A1, Table 1 decimal points.** `outputs/table_measured_envelope.md` and both masters carry
  `4.53`, `297.8`, `-9.4`; only the PDF lost the separators.
- **A2, the dropped "304" in the S6 closure sentence.** Line 131 of both supplementary files reads
  `1.148–1.323 at 304 and 0.783–0.915 at 256` in full.
- **A3, the undated Fend COMSOL contribution.** Reference [12] in both masters already gives
  2011 and the venue; the gap is in the EndNote library record, not the text.

Still present in the markdown and needing an edit before any resubmission:

- **A4, the duplications.** The S2 "recurring output of" paragraph appears twice in both
  supplementary files, and the sensor-selection result is stated twice in S3.
- **Every item in section B.** All fifteen flagged phrases that were string-searched are still
  present in both manuscript masters.

Separately, both masters still read `<mark>give github</mark>` in the data-availability section
while `submission_check_v6.md` confirms no `<mark>` survived into either PDF, so the markdown and
the submitted file disagree there. The URL it must carry is in section 4 below.

Two further items, neither of them errors, were raised and not resolved:

- The abstract is 357 words. Elsevier journals in this family typically cap abstracts at 250.
- Reviewer suggestions in the cover letter are inherited from the previous submission and match
  that paper's subject (solar thermochemistry, radiative transfer) rather than this one's
  (honeycomb absorber modelling, monolith channel heat transfer).

## 4. Public repository

`https://github.com/auth-transp/transp-receiver`, branch `manuscript`, local clone at
`/mnt/d/kkakosim/github/transp-receiver`. As of 2026-10-11 it is committed and pushed — 53 tracked
files, working tree clean, nothing ahead of `origin/manuscript`. Confirm it is set to public.

Contents, flat at the root: `README.md` (the reproduction README, maintained here as
`README_repository.md`), `receiver_reduction.jl`, `make_figures.jl`, `Project.toml`,
`Manifest.toml`, `.gitignore`, and the directories `RAW/` (20 logger CSVs), `outputs/` (18 archive
files) and `figures/` (8 PNGs). No manuscript or supplementary text is in the repository, by the
author's decision.

`Project.toml` and `Manifest.toml` there are **not** the ones at `aysha/`. They were generated for
the submission by taking the transitive closure of the fourteen packages the two scripts import:
149 packages against the parent project's 456, same pinned versions, referential integrity checked.
Two things remain undone and need a machine with Julia:

- The minimal `Manifest.toml` has no `project_hash`, which cannot be computed without Julia. One
  `Pkg.resolve()` stamps it and confirms the pinned set is self-consistent.
- The README gives no runtime figure for the reduction. The seven minutes quoted in the retired
  `README_reproduction.md` is the Python number and does not transfer.

## 5. If the paper comes back

The journal discussion, if a resubmission is needed: by citation gravity the paper belongs in
**Solar Energy** (7 of 53 references, and the Ávila-Marín, Kribus, Broeske, Becker and Patil work
the argument engages), with **International Journal of Heat and Mass Transfer** (5 references,
including Song 2026 and Kumar 2022) as the equally strong alternative that would judge the
transfer-unit scaling and eigenvalue identification on heat-transfer grounds. **Applied Thermal
Engineering** is the fallback, where the closest experimental analogue (D'Souza 2025) sits. Moving
between Solar Energy and IJHMT needs no restructuring, only a reframing of section 1. Applied
Energy and Energy Conversion and Management carry the highest metrics in the reference set but
reward a demonstrated performance gain or a system-level consequence, neither of which this paper
offers.

The paper's two defensible claims, in the form that survived review of the wording, are worth
restating because every framing decision turned on them. First, the measured transfer-unit count
rises with flow as `Re^+0.341 ± 0.040`, whereas a single-stream formulation whose volumetric
conductance is a fixed function of axial position requires `Re^-1.00`; the disagreement is in sign,
so no recalibration within that functional form removes it. The claim does **not** extend to a
different functional form, and no model was run in this work to test one — the title and abstract
were deliberately worded to avoid implying otherwise. Second, changing only which thermocouple
represents the solid moves the Nusselt prefactor by 29.7× and the identified capacitance from 122
to 281 J K⁻¹, so the conventional validation dataset cannot detect the first failure.

## 6. Open question raised after submission: the radiation–conduction number

Raised 2026-10-13, checked, **not** acted on — the decision is the author's.

`receiver_reduction.jl` line 343 computes `N_rc = 4 σ T̄_w³ D_h / k_air(T̄_w)`, that is with the
**gas** conductivity, 0.045 to 0.070 W m⁻¹ K⁻¹ over the campaign, giving the 1.53 to 5.57 quoted in
§4.2 and S5. The same file uses the solid conductivity correctly everywhere else it appears:
`K_SIC = 40` enters only the web Biot number at line 342, and the monolith capacitance uses the
solid heat capacity.

The difficulty is that both places the manuscript interprets `N_rc` read it as a statement about
the **solid**: §4.2 concludes "the solid redistributes absorbed energy axially by radiation", and
§5.2 calls it "such a lateral conductance" acting as a stabilising mechanism. Recomputing the same
group against `K_SIC = 40` gives 1.7×10⁻³ to 9.8×10⁻³, three orders of magnitude smaller, which
reverses that reading: conduction through the webs would dominate radiation, not the other way
around.

Both readings are defensible as physics; they are answers to different questions. The gas-referenced
group says radiation across the channel exceeds conduction through the gas inside it, which is what
bears on the probe-perturbation argument. The solid-referenced group says axial and lateral
redistribution in the webs is conductive rather than radiative. The §5.2 conclusion — that rising
transfer units occur in the presence of a lateral stabilising conductance — survives either way,
since web conduction is a lateral conductance too; only the attribution to radiation changes.

No edit has been made to the code, the manuscript or the supplementary files.

## 7. Reproducing

From the repository root, with Julia 1.12:

    julia --project=. receiver_reduction.jl --raw RAW --out outputs --nmc 4000
    julia --project=. make_figures.jl       --raw RAW --out outputs --fig figures

In the development tree the scripts live in `revision_2026-09/` with the raw files in
`../RAW/`, which is what their built-in defaults assume, so the flags can be omitted there.
Julia is not installed in the analysis sandbox, so any change to either script has to be run once
on the author's machine before its output can be trusted.
