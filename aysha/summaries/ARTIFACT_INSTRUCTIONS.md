# Instructions: Solar Receiver Analysis Artifacts

**Project:** `tamuq-chen-secarelab-receiver/aysha`
**Artifact URL:** https://claude.ai/code/artifact/2e8296f0-8b10-4dc4-8a61-9ce532664170
**Markdown export pattern:** `summaries/receiver_analysis_YYYY-MM-DD.md`
**HTML source pattern:** `summaries/receiver_analysis_YYYY-MM-DD.html`

---

## What This Artifact Is

A single-page HTML reference document that consolidates the state of the 1D and 2D mathematical modeling work for the SiC volumetric solar receiver into one navigable, theme-aware (light/dark) page. It is organized into numbered sections (§) covering the physical system, model hierarchy, heat transfer physics, model version history, key discoveries, open problems, forward path, and — after the first manuscript comparison — cross-check and revision-strategy sections.

The artifact serves as a **running audit trail** for the study. It is updated every time significant modeling progress is made (new version reaches a milestone, major bug is discovered and corrected, manuscript cross-check is completed) and exported to the repository `summaries/` folder.

---

## Section Structure

| § | Title | When to update |
|---|---|---|
| 1 | Physical System | Only if receiver geometry changes |
| 2 | Modeling Hierarchy | When a new model tier is added |
| 3 | Core Heat Transfer Physics | When governing equations change (new closure, new mechanism) |
| 4 | 1D Model Evolution | After every version milestone |
| 5 | Key Physical Discoveries | When a structural insight is confirmed across versions |
| 6 | 2D Model Evolution | After every 2D version milestone |
| 7 | Flow-Invariance Paradox | When the rear-slope problem is resolved or reframed |
| 8 | Current Calibrated State | After every successful calibration |
| 9 | Open Problems | After each modeling session — add / resolve items |
| 10 | Forward Path | Quarterly review |
| 11 | 1D vN: Repair & Baseline | Add when a new version reaches smoke-suite stage |
| 12 | Cross-Check: Model vs. Manuscript | Add/update when manuscript is revised |
| 13 | Scientific Value Assessment | Add when manuscript combination is being planned |
| 14 | Revision Strategy | Add when manuscript revision is active |

New sections are appended numerically. Never delete a section; mark superseded claims with "**RETRACTED yyyy-mm-dd — reason**" inline.

---

## How to Produce a New Version

### Prerequisites
- Read the relevant version's theory doc (`summaries/theory_1D_vNN.md`) and optimization summary (`summaries/1D_vNN/optimization_summary_1D_vNN.txt`)
- Read the current artifact content: `Artifact(action="read", url="https://claude.ai/code/artifact/2e8296f0-8b10-4dc4-8a61-9ce532664170")`
- The full HTML is saved by the read action to a local file; build all edits from that file

### Editing approach
Use Python to make targeted string replacements on the saved HTML file — never retype the whole document. Specifically:
1. Update header stats (version numbers, RMSE values)
2. Update the TOC `<nav>` with any new sections
3. Update §4 version table with the new version milestone
4. Update §8 (Calibrated State) with the new parameter table and objective breakdown
5. Append new sections (§11+) before the container `</div>` close, using the existing CSS classes (see below)
6. Update the footer version string

### CSS classes available
```html
<!-- Info/warning boxes -->
<div class="warn-box">...</div>     <!-- amber — cautions, not-yet-done -->
<div class="ok-box">...</div>       <!-- green — confirmations, passes -->

<!-- Parameter table rows -->
<div class="param-grid">
  <div class="param-row">
    <span class="param-key">Name</span>
    <span class="param-val">Value</span>
    <span class="param-desc">Description</span>
  </div>
</div>

<!-- Numbered discovery cards -->
<div class="discovery-card">
  <div class="d-num">I</div>
  <div class="d-content"><h4>Title</h4><p>Body</p></div>
</div>

<!-- Numbered step list -->
<div class="step-list">
  <div class="step">
    <div class="step-n">1</div>
    <div class="step-body"><h4>Title</h4><p>Body</p></div>
  </div>
</div>

<!-- Version cards (for §4 model history) -->
<div class="version-grid">
  <div class="version-item">
    <div class="v-tag">vNN</div>
    <div class="v-body"><h4>Title</h4><p>Body</p></div>
  </div>
</div>

<!-- Data table with scrollable container -->
<div style="overflow-x:auto; margin: 24px 0;">
  <table style="width:100%; border-collapse:collapse; font-size:13px; font-family:'IBM Plex Mono',monospace;">
    ...
  </table>
</div>
```

### Status badge colors
```html
<!-- Use these inline spans for status badges in tables -->
<span style="background:var(--ok-bg); color:var(--ok); padding:2px 8px; border-radius:3px; font-size:11px; font-weight:600;">AGREEMENT ✓</span>
<span style="background:var(--warn-bg); color:var(--warn); padding:2px 8px; border-radius:3px; font-size:11px; font-weight:600;">NOT YET COMPARED</span>
<span style="background:var(--accent-bg); color:var(--accent); padding:2px 8px; border-radius:3px; font-size:11px; font-weight:600;">CONFLICT</span>
<span style="background:var(--cool-bg); color:var(--cool); padding:2px 8px; border-radius:3px; font-size:11px; font-weight:600;">DIFFERENT PHYSICS</span>
```

### Publishing
```python
# Publish to same URL (keeps existing link)
Artifact(
    file_path="/tmp/receiver_analysis_vN.html",
    url="https://claude.ai/code/artifact/2e8296f0-8b10-4dc4-8a61-9ce532664170"
)
```

### Exporting to repository
After publishing, save both formats to `summaries/`:
```python
# HTML — commit via device bridge
device_commit_files(files=[{
    "fileUuid": "<uuid from SendUserFile>",
    "devicePath": "<connected_folder>/summaries/receiver_analysis_YYYY-MM-DD.html"
}])

# Markdown — write directly on device
device_bash(f"cp /tmp/receiver_analysis_YYYY-MM-DD.md $HOME/mnt/aysha/summaries/")
```

---

## Markdown Export Guidelines

The Markdown export mirrors the HTML content but uses plain tables and code blocks. It is intended for:
- Reading in any Markdown editor without a browser
- Version-control diff history (HTML diffs are noisy)
- Copy-pasting sections into papers or reports

**Naming convention:** `receiver_analysis_YYYY-MM-DD.md` (date of analysis session, not model version)

**Required sections in every export:**
- All §1–§10 (core)
- Every §11+ that exists in the current artifact
- A footer line: `*Generated by Claude (Cowork) · AUTh-TRANSP solar-receiver project · YYYY-MM-DD*`

---

## What Makes a Good Update Trigger

Update the artifact when **any** of the following occur:

| Event | Which sections to update |
|---|---|
| New 1D version passes smoke suite | §4 (add row), §8 (new params + objective), §9 (resolve or add items), new §11 |
| New 2D version completes a phase | §6 (add row), §9, possibly §7 |
| Bug or retraction found | Inline retraction note in the affected section; §9 add item |
| Manuscript revision | §12 (cross-check), §13, §14 |
| Calibration converges | §8 (replace parameter table + objective), §9 (update open problems) |
| New physical insight confirmed | §5 (add discovery), §7 (if paradox-related) |
| Forward-path priorities change | §10 |

---

## Design Principles

- **Never delete, only amend.** Use inline retraction notes with dates.
- **Status is always honest.** Use the `warn-box` class for anything not yet validated; never imply more certainty than the data supports.
- **Cross-checks report conflicts as findings.** Gaps between model and experiment are scientific results, not errors to hide.
- **Role A language for all model parameters.** Any fitted coefficient must carry the disclaimer that it is non-unique and model-dependent until Role B validation gates are passed.
- **The artifact tracks the study, not the success.** A version that failed is as important to document as one that succeeded.

---

*Instructions maintained as part of the `solar-receiver` project · AUTh-TRANSP · Last updated 2026-09-02*
