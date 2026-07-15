# LARMOR

A modern, open successor to [dmfit](https://nmr.cemhti.cnrs-orleans.fr/Dmfit/) for solid-state NMR lineshape fitting and spin-dynamics simulation — built to be usable by students with no NMR background and powerful enough for advanced research work.

Full architecture and workflow design: see the published design document (artifact link shared in the project conversation history). Summary of the plan:

- **Reuse, don't rebuild the physics.** The numerical core comes from [mrsimulator](https://mrsimulator.readthedocs.io), [lmfit](https://lmfit.github.io/lmfit-py/), [mrinversion](https://mrinversion.readthedocs.io), and [csdmpy](https://csdmpy.readthedocs.io) (P. Grandinetti's group, Ohio State). New work here is ingestion, UX, and orchestration.
- **Two-tier UX.** Guided mode for students, Expert mode (with an embedded Python console) for research use — same underlying project file.
- **Bruker/TopSpin import via [nmrglue](https://nmrglue.readthedocs.io)**, strictly read-only against instrument data.
- **Legacy dmfit `.fxmla` import** for migrating existing fits.
- **Advanced simulation as opt-in plug-ins**: SIMPSON (exact density-matrix, e.g. REDOR) and DFT/MD tensor import — not part of the default path.

## Status

**Phase 0 — feasibility: PASSED** (see `notebooks/phase0_feasibility.ipynb`, executed with outputs):

1. `CaAlGlass.fxmla` (dmfit Czjzek fit of a Ca-Al glass, 27Al) reproduced in mrsimulator with
   normalized RMSD 0.027 against the experimental spectrum embedded in the fxmla itself.
   **Key finding:** dmfit's `sCZ_CQ` = 2 × mrsimulator's `CzjzekDistribution` sigma — the
   conversion constant the future `.fxmla` importer must apply to every `CzSimple` line.
2. Bruker EXPNO 1903 (19F Hahn echo) opened with nmrglue end to end; an mtime+size snapshot
   of the whole EXPNO folder proves the read is strictly non-destructive. Lesson for the
   importer: `acqus` MASR (4200 Hz) disagrees with the operator-typed title (35.714 kHz) —
   surface both, trust neither silently.

**Phase 1 — core library: WORKING** (`larmor/` package, `pip install -e .`):

- `larmor.io.fxmla` — full dmfit `.fxmla` parser (both 1D and MQMAS files parse; MQMAS fitting
  itself is Phase 2) + conversion to LARMOR recipes with the σ = sCZ_CQ/2 factor applied.
- `larmor.io.bruker` — read-only TopSpin EXPNO reader with an enforced no-write guarantee
  (mtime+size verification) and metadata-conflict surfacing.
- `larmor.recipe` — the diffable JSON fit format (data referenced by path+SHA-256, never embedded).
- `larmor.engine` — fast fitting: the (Cq, η) Czjzek basis is simulated **once** per
  field/spin-rate/window via mrsimulator, then every fit iteration is a cheap reweighting —
  a full 3-site fit runs in seconds.
- `larmor.fit` — lmfit refinement writing values **and standard errors** back into the recipe.
  On `CaAlGlass.fxmla` the refined fit reaches normalized RMSD 0.0025 (10× better than the
  fixed-parameter replay) — and the uncertainties immediately show that the middle Al site is
  barely determined by the 1D spectrum alone (σ(Cq) = 0.9 ± 6.7 MHz), which is exactly why the
  MQMAS data exists. dmfit never told you that.
- `larmor` CLI — `larmor info <path>` (dmfit file or Bruker EXPNO), `larmor import`, `larmor fit
  recipe.json --plot out.png`. See `examples/`.
- `tests/` — 10 tests pinned against the Phase 0 numbers (`pytest`; `-m "not slow"` to skip the
  full kernel+fit run).

Sites whose center falls outside the fit window are frozen automatically — dmfit's ad-hoc
Gauss/Lor sideband lines (which LARMOR's simulation handles physically) stop wandering into
fake-baseline territory, and the covariance stays well-conditioned.

**Phase 1b — interactive app: first cut** (`larmor app`, then open http://127.0.0.1:8642):

- Load a dmfit `.fxmla` **or** a Bruker EXPNO folder by path; the spectrum plots immediately
  (Plotly, loaded from CDN — needs internet on first page load).
- Live model overlay: edit any site parameter and the simulation redraws in milliseconds
  (the Czjzek kernel is built once per field/spin-rate, then cached).
- One-click **Fit** with per-parameter “± error” shown next to each value, plus the full lmfit
  report (correlations included).
- Add/remove Czjzek and Gauss/Lor sites interactively — works for quadrupolar (27Al Czjzek) and
  spin-1/2 (19F pseudo-Voigt, no kernel needed) alike.
- Saving a recipe **into an instrument data folder is refused** (HTTP 403) — the read-only
  guarantee is enforced server-side, not just promised.

**Phase 1c — constraints (ssNake-inspired, algebraic):**

- Every recipe parameter supports `vary` (fix), `min`/`max` (bounds), and `expr` (links):
  `"expr": "0.29 * s0.amplitude"` locks an amplitude ratio, `"expr": "s0.shift_fwhm_ppm"`
  shares a linewidth. Any lmfit-valid algebra works, with **full error propagation** —
  a linked parameter's stderr is derived, not fitted.
- Bad expressions fail before the fit with a message naming the valid parameters (and come
  back as clean HTTP 422s in the app).
- **At-bounds diagnosis**: parameters that finish a fit pinned at a bound are reported
  (report, app UI ⚠, and a note written into the recipe) — the usual sign that a constraint
  or starting model is fighting the data. Uncertainties are then computed conditional on the
  pinned values instead of silently vanishing.
- In the app: click *constraints ▸* on any site to edit link/min/max per parameter; linked
  parameters grey out, show ⚭, and follow their expression live in the plot.

**Tutorials** (ssNake-style, in [docs/tutorials/](docs/tutorials/)):

1. [Your first fit — ²⁷Al MAS of a glass](docs/tutorials/01-first-fit-27Al-czjzek.md)
2. [Constraining a fit: fix, bound, link](docs/tutorials/02-constraints.md)

Next: Guided-mode layer (plain-language panels, guardrails) on this app, then Phase 2
(MQMAS/2D methods — `CaAlGlassMQ.fxmla` already parses, fitting it comes with the 2D engine).

## Data policy

Raw instrument data (TopSpin EXPNO folders) and legacy `.fxmla` files are **never copied into this repo** and **never written to**. `data/` holds only small reference notes (paths, hashes) — see `data/README.md`.

## Environment

```
conda env create -f environment.yml
conda activate larmor
```
