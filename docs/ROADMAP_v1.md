# Road to LARMOR 1.0

v0.2.0 is feature-rich. **1.0 is a different claim**: *"you can rely on this for
published work and recommend it to a colleague."* That shifts the priority from
*more features* to *trust, robustness, and sustainability*. This is what 1.0 should
mean and what it takes to get there.

## What "1.0" must guarantee
1. **Trustworthy numbers** — validated against dmfit / literature on real data,
   with honest uncertainties everywhere.
2. **Never loses work, never crashes** — clear errors, recoverable sessions,
   forward-compatible file formats.
3. **Maintainable** — no single 4k-line file; CI keeps it green.
4. **Installable & documented** — a colleague can install it and learn it without
   you.

## Pillars (in priority order)

### A. Trust & validation — *the most important for a scientific tool*
- **Validation report**: reproduce a set of published/known fits (and dmfit
  results) on real datasets, with the agreement tabulated and tolerances stated.
  This is the single most valuable thing for adoption.
- **dmfit round-trip** tested across many real `.fxmla` files (import→refit→export
  agrees).
- **Uncertainties everywhere**: surface the Monte-Carlo option consistently, and
  propagate errors to *every* derived quantity (Cq, P_Q, populations, δ₂).

### B. Robustness & data safety
- **File-format versioning + migration** for the recipe / `.larproj` — an old
  recipe must always open in a newer LARMOR.
- **Project bundles** (`.larproj`) that fully round-trip a session (data refs,
  fits, processing, baselines, batches, figures).
- **Every dialog constructs headlessly** (a smoke test) + **fuzz the loaders**
  (never crash on a malformed file — degrade with a clear message).
- Continue the error-quality work started in v0.2 (`humanize_error`,
  `sanitize_constraints`).

### C. Maintainability
- **Split `app.py`** (4.2k LOC / 216 methods) into cohesive mixins
  (fitting, 2D, session/workspace, reporting, plotting) — no behaviour change, the
  #1 audit recommendation and a prerequisite for sustainable 1.0 upkeep.
- **CI** (GitHub Actions): run the offscreen test suite on every push; track
  coverage; lint.

### D. Packaging & distribution
- **One-click installers** — Windows (PyInstaller spec already exists) and,
  ideally, macOS; handle SmartScreen/signing; an update check on launch.
- A **release pipeline** that builds installers on a version tag.

### E. Documentation & onboarding
- Complete the in-app manuals; add a **getting-started tutorial** with **bundled
  example datasets**.
- A short **README/website** with screenshots and a "LARMOR vs dmfit" page.
- Document the **headless API** (`larmor fit / batchfit / seqfit`, the Python
  entry points).

### F. Feature completeness (close the known gaps)
- Field-use items still open: **#14** labelable peak-pick *assignment table*,
  **#12** TopSpin *drag-to-phase* gesture, **#19** dmfit *amplitude-calibration*.
- Beyond-request roadmap still open: **bootstrapped series errors**, **two-way
  batch table** (row ↔ plot), **figure annotation layer**, **reference/literature
  overlays**, **watch-folder / auto-fit**.
- Quality-of-life: **undo everywhere** (processing & 2D), **keyboard-driven
  fitting**.

## Proposed milestones

| Release | Theme | Headline |
|---|---|---|
| **v0.3** | Trust & safety | validation report · format versioning · project bundles · CI · dialog smoke test |
| **v0.4** | Sustainability & reach | split `app.py` · Windows+macOS installers · docs site + tutorials |
| **v0.5** | Feature completeness | assignment table · drag-to-phase · amp calib · bootstrapped errors · annotation layer |
| **v1.0** | Sign-off | validation signed off · zero known crashes · "recommend to a colleague" |

## If only three things before 1.0
1. **A validation report** (turns "it fits" into "it fits *correctly*").
2. **Project bundles + format versioning** (never lose or strand a fit).
3. **Split `app.py` + CI** (so 1.0 stays maintainable and green).

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
