# Julia scientific-script style for this analysis

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
