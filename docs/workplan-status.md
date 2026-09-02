# Workplan status — 2026-09-02 session handoff

Where the 50-item improvement workplan stands, what the last session changed,
and exactly what a continuing session picks up next. The workplan itself
(items, evidence, Sam's keep/defer/drop decisions) lives in the artifact
"LARMOR Workplan"; this file is the ground truth for progress.

## State of the repository

- **Version**: 0.12.1 (`larmor/__init__.py` + `pyproject.toml`, bumped in step).
- **13 local commits ahead of `origin/master`** (v0.10.1 … v0.12.1). **Not
  pushed** — pushing needs Sam's explicit word, after a full-suite green run.
- **Tests**: ~737 collected; the last full run (at v0.11.4 + pool fix) was
  731 passed / 0 failed, and every scoped batch since has been green. The
  terminal now prints a red/green **real-data layer banner** — trust a green
  bar only when the banner says "complete (all 12 datasets present)".
- **Environment**: conda env `larmor` (`C:\Users\samso\xraylarch\envs\larmor`,
  Python 3.11); run tests with that interpreter. Real test data roots at
  `LARMOR_TEST_DATA` (default `C:\Users\samso`); the CaAlGlass fxmla files
  live in `Desktop\larmor_tests\`.
- **Installer**: `dist/LARMOR/` builds and smoke-runs from the pip venv at
  `packaging/.buildenv` (python.org 3.11 — **never** the conda env; the two
  hard-won reasons are in `packaging/README.md`). Still owed: verification on
  a machine with no dev setup, and publishing a release.

## Done — phases 1–4 plus the start of 5

| Phase | Items | Commits |
|---|---|---|
| 1 · Static real | A1 A5 A6 A7 | v0.11.0 |
| 2 · Stalls | B1 B2 B3 B8 C1 C2 C3 | v0.11.1 |
| 3 · Handover | E1 E4 E5 D4 D7 | v0.11.2 |
| 4 · Books + wideline | D1 D2 D3 D5 D6 G3 E7(build) B4 B5 B7 A3 A4 | v0.11.3, v0.11.4, v0.12.0 |
| 5 · (started) | B6, A8 | ec93d39, v0.12.1 |

Highlights a continuing session should know about (details in the commit
bodies, which are the real documentation):

- **Static is a first-class experiment**: `MASR=0` → static (corroborated by
  title/pulse program); seeding probes at the true spin rate; kHz/MHz axis
  display; static physics anchors in `test_physics_validation.py`.
- **Kernels**: axis derived arithmetically (`engine.kernel_axis_ppm`, against
  `Isotope.B0_to_ref_freq`, *not* `larmor_MHz` — 0.08 % apart), float32,
  LRU-bounded in memory, persisted to `%LOCALAPPDATA%\LARMOR\kernels`,
  chunk-built with progress + Stop, pre-warmed in the background on load.
- **dmfit CzSimple interop**: import was 4.1× too tall (hidden by refits);
  one shared constant now serves both directions; round trip is identity and
  both facts are pinned by tests.
- **Wideline acquisition**: `larmor/vocs.py` + Tools ▸ Stitch VOCS spectra;
  `op_wurst_correct` + Process ▸ WURST excitation profile.
- **A8 (last item finished)**: `larmor/staticct.py` reads (C_Q, η, δiso) off
  a static CT pattern from three marker positions — no fit. The features are
  found as peaks of the ideal powder histogram of the exact Freude–Haase
  frequency surface; do **not** "simplify" it back to the classical
  φ = 0°/90° critical-point formulas — they miss the pattern's *dominant*
  horn, which comes from an interior-φ saddle (measured: +604 ppm on an
  ⁸¹Br-like case where the branch formulas put nothing). Powder-averaging the
  surface reproduces `convert.ct_second_order_shift_ppm` exactly; inversion
  round-trips C_Q to 3 %, η to 0.03. UI: Tools ▸ Read static pattern, with
  one-click quad_ct seeding.

## Remaining — Phase 5

In the order I would do them:

1. **F6 — time↔frequency toggle** (M): back to the FID for re-apodization
   without reloading, plus real/imaginary display while phasing. The ops
   (`ift`, `real`, `imag`) already exist in `processing.py`; this is
   workbench state + UI.
2. **F2 — TopSpin drag-to-phase** (M): horizontal drag = p0, vertical = p1
   about the existing pivot line (`plot.py` already has `_pivot` and
   `phase_pivot_frac()`).
3. **F7 — one fit-health panel** (M): `sanity.py`, `identifiability.py`,
   `diagnostics.py` already compute everything; gather into one always-visible
   verdict strip instead of three dialogs.
4. **F8 — command palette** (S/M): fuzzy search over the ~85 menu actions;
   the `_add()` helper in `app.py` is the natural registry to walk.
5. **F3 — sideband auto-detect** (M): find the ±νrot repeat in the data,
   offer the linked manifold or the shifted-copy component in one click
   (build on Decomposition ▸ Add a copy of this spectrum).
6. **F5 — two-way batch table** (M): row ↔ spectrum highlighting in
   `batchfit_dialog.py`.
7. **G4 — three distributions** (L): Gaussian Isotropic Model, correlated
   (δiso, C_Q), two-site exchange. NOTE: every new model must satisfy
   `tests/test_model_tables.py` (four explicit partitions) and the
   `lineshapes.md` backtick gate — the tests name each decision owed.
8. **G1 — complete the project bundle** (L): 2D workspaces, figures, batch
   sessions into `.larproj.json`; bump `PROJECT_BUNDLE_VERSION` and use the
   (new, tested) migration path.
9. **G2 — four tutorials** (M/L): batch fitting, MQMAS, error analysis, and
   a static ⁸¹Br WQCPMG walkthrough (QCPMG dialog → send → Read static
   pattern → fit → infinite-field δiso is the arc; PBi data conventions are
   in Sam's memory notes).
10. **G5 — split app.py** (L, riskiest last): ~4.9k lines. Suggested cuts:
    menus/actions builder, session/workspace persistence, the
    add-line/seeding block, processing glue. Keep `MainWindow` as a facade so
    tests and `_win`-style fixtures keep working; move in small verified
    steps, running `tests/test_desktop.py` between each.

Deferred by Sam (revisit only if asked): A9 C4 C6 E2 E3 F1 F4.
Dropped: A2, C5.

## Conventions that keep this project safe

- Scope test runs to what changed; **full suite + real-data banner before any
  push**, never push without Sam's word.
- Version bumps touch `larmor/__init__.py` **and** `pyproject.toml`.
- Commits are `vX.Y.Z: summary` with a post-mortem body (what was wrong, how
  found, measured effect) — `git log` is the project's bug documentation.
- Adding a lineshape model: the six-tables partition tests + the help gate
  make omissions fail loudly; read `docs/development-notes.md` §6 first.
- `docs/development-notes.md` was updated through v0.11.2; phases 4–5 changes
  are documented in commit bodies but not yet folded into it — worth a pass
  when the workplan closes.
