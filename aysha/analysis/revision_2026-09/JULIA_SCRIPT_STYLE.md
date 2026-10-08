# Julia scientific-script style for this analysis

## Current project-specific conventions

This guide describes the style used for the Julia reduction and plotting
scripts in `analysis/revision_2026-09`. The primary style references are
`0D_v1.jl` and `1D_v1.exp.All.jl`. The goal is a readable, top-to-bottom
scientific calculation in which the inputs, equations, intermediate results,
and exported products remain visible in the source file. Do not use
`1D_v1.jl`, `1D_v46.jl`, or `0D_v4.jl` as the structural model for these
scripts.

The current reference implementations are:

- `receiver_reduction.jl` for data reduction, uncertainty propagation, model
  inversion, and table/data export;
- `make_figures.jl` for manuscript figure generation.

If any older general guidance later in this document conflicts with this
project-specific section, this section takes precedence.

### 1. Prefer a linear scientific script

Organize the file as a sequence of visibly labelled sections. A typical order
is:

1. short file description and example run command;
2. libraries;
3. command-line options and paths;
4. campaign definitions and fixed experimental constants;
5. geometry, instrument uncertainties, and property data;
6. a small set of genuinely reusable numerical or I/O helpers;
7. data loading and preprocessing;
8. calculations in the same order as the physical argument;
9. uncertainty propagation or inversion;
10. table, JSON, or figure export;
11. a guarded executable block that prints progress and output paths.

Use `begin ... end` blocks with descriptive comments to divide those sections
when this improves navigation. Keep the scientific workflow at file scope when
possible. A reader should be able to scroll downward and see how raw data
becomes the reported result without following a chain of orchestration
functions.

### 2. Keep the function budget small

Introduce a function only when at least one of these conditions holds:

- the operation is used more than once;
- a solver, optimizer, interpolation routine, or package API requires a
  callback;
- the operation is a compact scientific equation that benefits from an
  explicit name and can be tested independently;
- file parsing or export would otherwise obscure the calculation;
- the same calculation must be called repeatedly by Monte Carlo sampling.

Functions should be narrow and named after the scientific operation they
perform. Avoid classes, elaborate type hierarchies, manager objects, pipeline
frameworks, and thin wrapper functions that merely call another function.
Avoid putting the whole workflow inside `main()` when doing so hides the
analysis sequence. A few helpers plus a linear top-level calculation is the
preferred balance.

### 3. Separate inclusion from execution

Files may be both run directly and included by another script. Use a simple
include-only flag so including a file exposes its definitions without starting
the full calculation. The reduction script follows this pattern conceptually:

```julia
RECEIVER_REDUCTION_RUN =
    !isdefined(@__MODULE__, :RECEIVER_REDUCTION_INCLUDE_ONLY) ||
    !RECEIVER_REDUCTION_INCLUDE_ONLY

if RECEIVER_REDUCTION_RUN
    # parse options, run the ordered calculation, and write outputs
end
```

Before `make_figures.jl` includes `receiver_reduction.jl`, it temporarily sets
`RECEIVER_REDUCTION_INCLUDE_ONLY = true`, then restores the previous state.
`make_figures.jl` uses the corresponding `MAKE_FIGURES_INCLUDE_ONLY` convention
for its own executable block. Do not silently launch a long reduction merely
because another script needs a shared definition.

When a script is run normally from VS Code or the terminal, it must give
visible progress messages and state where its products were written. A script
that appears to do nothing is not an acceptable user interface for a long
calculation.

### 4. Define campaigns and physical constants once

Keep campaign identifiers such as heat-flux levels and nominal cooling-flow
labels in one visible location and reuse them throughout reduction, tables,
and figures. Do not duplicate those lists in distant functions. Preserve the
distinction between a nominal campaign label and a measured quantity: for
example, read the measured mass-flow-controller signal from the logger rather
than replacing it with the nominal set point.

Write units and physical meaning beside every dimensional constant. Prefer
names that match the manuscript notation or instrument identifier. Make sensor
corrections explicit, for example `T3_offset`, and use the physically neutral
default (normally zero) unless a documented correction is intended. Reject
invalid measurements such as non-positive measured flow instead of inventing
a fallback value.

### 5. Keep equations recognizable

Translate the governing equations directly into Julia expressions, with
variable names close to the notation in the paper. Comments should explain the
physical assumption, unit conversion, or reason for an equation; they should
not narrate elementary syntax. Keep related intermediate quantities visible
instead of compressing a multi-stage derivation into a single opaque
expression.

Use `DataFrame` columns and `NamedTuple` rows for tabular experimental data,
and use nested `Dict` objects only where a hierarchical JSON report is the
natural result. Do not build an object model merely to carry scalar quantities
between consecutive sections.

### 6. Use the project packages for standard numerical work

Use the packages pinned by `Project.toml` and `Manifest.toml` instead of
reimplementing established numerical algorithms:

- `CSV.jl` and `DataFrames.jl` for campaign and result tables;
- `GLM.jl` for regression;
- `Trapz.jl` for numerical integration;
- `Statistics` for means, standard deviations, and quantiles;
- `Roots.jl` for scalar root finding;
- `Optim.jl` and `ForwardDiff.jl` for parameter inversion;
- `Interpolations.jl` for property-table interpolation;
- `JSON3.jl` for the archived machine-readable report;
- `PythonPlot.jl` for the manuscript figures.

The current reduction obtains air properties through the C interface supplied
by `CoolProp_jll` and records the CoolProp version in the output. For speed,
construct a sufficiently resolved property table once and interpolate during
the repeated calculation. Keep the property source, temperature range, and
grid spacing visible and reproducible.

A short custom routine is acceptable when automatic differentiation, exact
reference behaviour, or a small well-defined transformation requires it. It
must not grow into a private replacement for a maintained package.

### 7. Preserve computational parity deliberately

When replacing the Python implementation, match both the physical equations
and the reporting semantics. Check whether a reported value is the nominal
solution, the mean of Monte Carlo samples, a median, or a confidence bound.
These are not interchangeable.

In the present Table 4 export, the reported `Nu_a` value and each reported
`epsilon_star` value follow the Monte Carlo means used by the reference Python
implementation, while the associated uncertainty columns are generated from
the same archived uncertainty results. Generate all rows of a manuscript table
from the computed report; do not manually transcribe selected rows or maintain
a second set of display values.

Use explicit, stable formatting (`@sprintf` where appropriate) so regenerated
Markdown and CSV products produce reviewable diffs.

### 8. Make uncertainty calculations reproducible

Seed stochastic calculations explicitly. Keep the seed, sample count, profile
starts, and relevant analysis options visible as command-line options or named
constants. Store enough nominal and uncertainty information in the JSON output
to reproduce tables and figures without rerunning the expensive analysis.

Use small sample counts only for smoke testing. Production manuscript products
must use the documented production settings; a quick test output must not be
mistaken for an archival result.

### 9. Keep figure construction linear and explicit

In `make_figures.jl`, use one labelled top-level block per manuscript figure.
Shared visual settings and only genuinely repeated plotting operations belong
in small helpers. Each figure block should read like a plotting recipe:
select the data, plot each series, label axes, configure the legend, save, and
print the output path.

Figure labels must identify experimental sensors when that information matters.
For example, the temperature-profile figure labels the thermocouples as well
as their axial positions (`T8`, `T12`, `T11`, and gas thermocouple `T3`). Keep
the legend placement explicit when it is part of the approved layout.

Open or white markers must always have an explicit visible edge. Set the marker
face and edge independently, for example:

```julia
ax.plot(x, y; marker="o", mfc="none", mec=color, mew=1.0)
ax.plot(x, y; marker="s", mfc="white", mec=series_color, mew=1.0)
```

This is required even if the edge colour seems to be inherited, because the
shared plotting style may set the global marker-edge width to zero. The gas
temperature in Figure 3 uses open circles with a solid coloured boundary, and
the white markers in Figure 5c use an explicit series-coloured edge.

### 10. Validate execution and products

Run every Julia command in the repository project environment. This workspace
currently uses Julia 1.12.6:

```powershell
& "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. analysis\revision_2026-09\receiver_reduction.jl
& "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. analysis\revision_2026-09\make_figures.jl
```

Validation should include all of the following, in proportion to the change:

- include-only smoke tests for both scripts;
- a reduced-cost end-to-end reduction in a temporary output directory;
- regeneration of the tables and JSON report;
- execution of the figure script in a temporary figure directory;
- visual inspection of every changed figure, not merely confirmation that a
  PNG exists;
- comparison of key Julia values and exports with the Python reference when
  parity is part of the task;
- a final review of the diff and a whitespace/error check.

Do not overwrite the manuscript's production outputs while performing a smoke
test. Generate test artifacts in a temporary or explicitly separate directory,
then copy only deliberately approved production products.

### 11. Prompt for an AI coding tool

The following wording captures the current style:

> Write this as a linear, literate Julia scientific script modelled on
> `0D_v1.jl` and `1D_v1.exp.All.jl`. Keep campaign inputs, physical constants
> with units, equations, intermediate quantities, and export steps visible in
> top-to-bottom order. Use only a small number of narrow functions for repeated
> calculations, solver callbacks, Monte Carlo kernels, or file I/O; do not
> introduce classes, framework-style pipelines, or a `main()` function that
> hides the workflow. Use the packages in this project's `Project.toml` for
> standard numerical work. Support safe include-only use and normal direct
> execution, print progress and output paths, seed stochastic work, generate
> tables from computed results, preserve the Python implementation's reporting
> semantics, and validate both numerical products and rendered figures. For
> open or white plot markers, set a visible marker edge explicitly.


The canonical style references are `0D_v1.jl` and `1D_v1.exp.All.jl`.
`1D_v1.jl` is not a style reference for this work.  Later architectural
scripts such as `0D_v4.jl` and `1D_v46.jl` are also not the target style.

## Programming structure

Write a single executable scientific script whose order tells the story of the
calculation.  Use this top-to-bottom structure:

1. `begin # libraries`
2. `begin # fixed parameters` or campaign definitions
3. `begin # properties, measurements, and governing equations`
4. `begin # define functions`
5. `begin # reduction or optimization`
6. `begin # plots, tables, and export`
7. a short end-to-end execution block

Keep scientific quantities visible as named Julia variables and keep equations
close to the parameters they use.  Prefer direct formulas, arrays, `Dict`s, and
`DataFrame`s.  A calculation performed only once should normally remain in its
top-to-bottom `begin` block.  Use functions only for operations that repeat,
numerical kernels that need callbacks or recursion, and cohesive algorithms
whose internals would obscure the experimental sequence if written inline.

The intended structure is a **sectioned procedural scientific workflow** or a
**literate/notebook-like Julia script**.  It is not an object-oriented design,
an application framework, or a hierarchy of configuration and model types.

## Naming and presentation

- Use domain names such as `Re`, `Nu`, `NTU`, `Tw_K`, `Q_gas_W`, and
  `C_eff`, including unit suffixes where useful.
- Group constants by physical meaning and add short comments stating units or
  experimental meaning.
- Keep transformations explicit: raw data -> steady values -> dimensionless
  groups -> fitted quantities -> figures and tables.
- Keep the principal workflow at top level; do not hide it in `main`,
  `run_analysis`, `make_figures`, or one function per figure.
- Keep model equations readable in the main script instead of hiding them
  behind generic dispatch layers.
- Use `begin # descriptive section` blocks as visual and executable cells.
- Let the file run from start to finish and produce its scientific artifacts as
  side effects.
- Use the packages already defined by the repository `Project.toml` and
  `Manifest.toml` and run with `julia --project=.`.

## What to avoid for this style

- Do not introduce a module, custom type hierarchy, model class, registry, or
  plugin architecture unless the calculation genuinely requires it.
- Do not add version tags to every function or variable name.
- Do not divide a straightforward scientific calculation into many source
  files merely to enforce software-application layering.
- Do not replace recognizable physical equations with generic configuration
  plumbing.
- Do not use large banner comments and exhaustive architectural narration as
  the main organizing device; use short section labels and local scientific
  comments.

## Copyable instruction for an AI coding tool

> Write this as a top-to-bottom Julia scientific analysis script, following
> `0D_v1.jl` and `1D_v1.exp.All.jl` as the canonical style references.  Organize
> the file with descriptive `begin # ...` sections for libraries, fixed
> parameters, data/properties/equations, helper functions, reduction or
> optimization, and plots/exports.  Keep physical variables and equations
> explicit and near their use.  Prefer direct Julia expressions, arrays,
> dictionaries, and data frames.  Inline every one-use analysis and plotting
> stage; retain functions only for repeated operations or necessary numerical
> kernels.  Do not wrap the workflow in `main`, `run_analysis`, or one function
> per figure.  The result must execute sequentially from start to finish with
> `julia --project=.`.  Do not use `1D_v1.jl`, `0D_v4.jl`, or `1D_v46.jl` as
> style references, and do not redesign the analysis as a module, class/type
> hierarchy, or application framework.
