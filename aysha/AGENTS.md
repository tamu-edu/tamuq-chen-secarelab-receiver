# Agent Instructions

## Julia Environment Notes

- Always run Julia commands with the project environment:

```powershell
julia --project=. path\to\script.jl
```

- Before choosing a Julia executable, check `Manifest.toml`. This workspace currently has:

```toml
julia_version = "1.12.6"
```

- In this Codex/Windows environment, the plain `julia` launcher may fail with:

```text
The Julia launcher failed to figure out which juliaup channel to use.
```

When that happens, use the explicit juliaup executable matching the manifest, for example:

```powershell
& "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. test\smoke_1D_v7.jl
```

- Do not use Julia 1.11.x with this manifest unless the project is intentionally downgraded. It can produce package precompile errors such as:

```text
UndefVarError: `StaticData` not defined in `Base`
```

- If Julia fails with permission errors under `C:\Users\kkakosim\.julia\compiled\...`, that is usually package-cache access from the sandbox, not a model/test failure. Rerun the same Julia command with elevated permission so Julia can write its normal compiled cache.

- For this receiver project, quick validation commands that have worked are:

```powershell
& "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. test\smoke_1D_v6.jl
& "C:\Users\kkakosim\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe" --project=. test\smoke_1D_v7.jl
```

## Documentation Notes

- Periodically check `summaries/journal.1D.md` when making or evaluating 1D model changes. Update it with meaningful new model versions, calibration outcomes, validation notes, and interpretation changes so the journal remains aligned with the current code and results.
