# LARMOR — full app audit (2026-08-05, updated 2026-08-06)

A pass over every feature area, how they fit together, and where the efficiency
and maintainability risks are. Written after the 2026-08 feature waves; updated
after the 2026-08-06 performance pass (see §6).

## 1 · Scale

| | |
|---|---|
| Python modules | 106 (`larmor/` core + `larmor/desktop/` UI) |
| Lines of code | ~28 k (≈8.6 k core, ≈15 k desktop, rest tests/scripts) |
| Test files | 59 · full suite **green** (433 passed, 16 skipped) |
| Biggest module | `desktop/app.py` — 4.4 k LOC, 227 methods (the main window) |

The core (`larmor/…`) is **Qt-free and testable**; the desktop layer is a thin-ish
GUI over it. This separation is the project's biggest strength: every engine
(fit, batch, seqfit, quantify, diagnostics, methods, chi²) has unit tests that run
headless in seconds.

## 2 · Feature map (what exists, where)

**Ingestion** — `io/bruker.py` (1r/2rr/fid/ser, EXPNO/pdata, pseudo-2D detection),
`io/varian.py` (VnmrJ .fid, validated on real ²⁷Al/²³Na), `io/fxmla.py` (dmfit
round-trip), `io/spectra.py` (CSV/txt), `io/scan.py` (sample auto-identify),
`fourier.py` (States/TPPI/echo-antiecho recombination). One entry: `loader.load_any`.

**Processing** — `processing.py` (EM/GM/SINE/ZF/FT/phase/baseline/LP/…, live &
absolute from an unprocessed base), `baseline.py` (Yon 2020 iterative).

**Models** — `models/` registry (gauss_lor, czjzek, ext_czjzek, quad_ct,
quad_csa, quad_first, csa_mas, amorphous, spectrum, function). One kernel per
(nucleus, field, spin, window), **cached**; analytic models use parameter-level
LRU.

**Fitting** — `fit.py` (lmfit least-squares, constraints/links, bounds, at-bounds
covariance retry, completion threshold, windowed simulation grid,
**live-animation frame_cb**), `multifit.py` (co-fit 1D+2D **genuinely shared**
params via `share=`), `batchfit.py` (one shared-shape model, many 1D spectra
fit **independently** per spectrum via `fit.py` — not a joint optimisation,
since nothing is ever tied across spectra there), `seqfit.py`
(forward-backward warm-start series), `autofit.py` (multi-start + MC errors),
`twod.py` (MQMAS fit/shear/kernels).

**Interpretation & honesty** — `quantify.py` (integral populations ± err),
`sanity.py` (η/width/amplitude/window physical flags), `identifiability.py`
(|r|≥0.95 degenerate pairs), `diagnostics.py` (runs-test + autocorrelation
structure), `chi2map.py` (parameter-pair χ² surface), `recipe_diff.py`
(current vs reference), `convert.py`/`czjzek_dist.py` (Cq↔σ↔P_Q, P(Cq)).

**Experiments** — `satrec.py`/`series.py` (T1/T2, per-site NNLS), `redor.py`,
`qcpmg.py`, `dft.py` (.magres → sites), `vt` (Arrhenius/VFT).

**Output** — `figures.py` (spec-driven 1D/2D/series, journal presets, contour
modes, export png/pdf/svg/tiff), `methods.py` (LaTeX table + methods sentence),
`batch.py` (publication table across fits), `io/export.py` (fit → txt/csv/dmfit).

**UI power tools** — Explorer (procs/fits/pins), Plotting studio (+ built-in file
explorer), batch/sequential dialogs, series-evolution, workspace manager, 2D
phasing, HMQC, constraint library, per-nucleus defaults, timed autosave.

## 3 · How it fits together (data flow)

```
load_any ─▶ (ppm, amp, recipe) ─▶ processing ─▶ workbench
   │                                              │
   ▼                                   place lines (models registry)
 NMRData (2D/fid) ─▶ 2D view / FID dialog          │
                                                   ▼
                                    fit.py ◀── FitWorker (thread, animate, converge)
                                       │            ▲ same engine reused by:
                                       │            ├─ batchfit / seqfit (series)
                                       ▼            ├─ multifit (co-fit)
                        quantify · sanity · identifiability · diagnostics · chi²
                                       │
                                       ▼
                       figures / methods / batch  ─▶ publication bundle
```

Key integration points that make it cohere:
- **One fit engine** (`fit.fit`) is reused by the main window, batch, sequential,
  co-fit, chi²-map and MC — so the completion threshold, constraints and error
  handling behave identically everywhere.
- **One figure engine** (`figures.render` on plain-dict specs) backs the Plotting
  studio, the batch report and the publication bundle; specs are saveable/reloadable.
- **One export dialog** (`export_dialog`) and **one plot right-click menu**
  (`plot_menu`) are attached to every plot (export + send-to-studio); the crash-prone
  native pyqtgraph exporter is disabled globally.
- **QSettings** carries cross-session state: theme, font, recents, pinned folders,
  scroll-nudge, fit threshold, animate flag, per-nucleus site defaults, constraint
  library, batch templates.

## 4 · Efficiency assessment

**Good:**
- **Kernel caching** — the mrsimulator Czjzek basis is built once per
  (nucleus, field, spin, window) and reused; this is why a Czjzek fit runs in
  seconds, not minutes.
- **2D redraw decimation** — contour grids capped at ≤512/axis (`twod_view._decimate`)
  so a 128×2048 map redraws in ~10 ms warm.
- **Fit animation throttle** — the live model is recomputed/redrawn only every ~10
  iterations (`fit.frame_every`), so animation adds ~no cost.
- **Completion threshold** (default 0.1 % ≈ dmfit 1e-3) stops fits when the stdev
  stops moving, cutting needless iterations.
- **Lazy imports** of heavy deps (mrsimulator, matplotlib) keep startup snappy.

**Watch / improve:**
- **`app.py` is a 4.2 k-LOC, 216-method monolith.** It works and is stable, but it
  is the main maintainability risk. Recommendation: extract cohesive mixins
  (fitting, 2D, session/workspace, reporting) — mechanical, low-risk, high payoff.
- **χ² map cost** — `chi2_surface` runs n×n simulations. Fine for gauss_lor and
  cached Czjzek at n≤15; could be slow for a heavy model at n=41. It runs under a
  wait cursor; a coarse-then-refine or a max-time guard would harden it.
- **Undo now stores arrays** for axis-changing ops (calibrate/SR). 1D-only and
  bounded to 60 states, so memory is fine; keep it 1D-scoped.
- ~~`np.trapz` deprecation~~ **done** — both `figures.py` and `measure.py`
  already use `np.trapezoid` when available (`getattr(np, "trapezoid", None)
  or np.trapz`), falling back only on numpy < 2.0. No warning fires on the
  numpy version this app ships against.

## 5 · Recommendations (priority order)

1. **Refactor `app.py`** into mixins (biggest maintainability win; no behaviour change).
2. Harden `chi2_surface` with a max-time / coarse-refine guard.
3. Remaining *ideas-from-field-use* not yet built (all niche/interaction-heavy):
   **#14** labelable peak-pick assignment table, **#12** TopSpin drag-to-phase
   gesture, **#19** dmfit amplitude-calibration UI.
4. Add a smoke test that constructs every dialog headlessly (catches wiring
   regressions cheaply given how much UI now exists).
5. A `pyflakes larmor/` sweep (2026-08-06) turned up only cosmetic debt outside
   today's changes — unused imports in several untouched dialog modules
   (`baseline_dialog.py`, `batch_dialog.py`, `chi2map_dialog.py`,
   `czjzek_dist_dialog.py`, `datasets.py`, `export_dialog.py`, `fid_dialog.py`,
   `figure_dialog.py`, `integrate_dialog.py`, `qcpmg_dialog.py`,
   `twod_dialog.py`, `utilities.py`, `workspaces.py`, `cofit_dialog.py`) and the
   deliberate `__init__.py` re-exports (`larmor/__init__.py`,
   `larmor/models/__init__.py` — real public API, not bugs). Two genuine dead
   locals found and removed in this pass (`io/bruker.py`'s orphaned
   `is2d_expno` lambda, `app.py::edit_mqmas_f1_ref`'s unused `vary` read) —
   both provably inert (an unused local can't change behaviour), verified by
   the existing test suite staying green. The rest is cosmetic and low
   priority — worth a batch cleanup pass sometime, not urgent.

## 6 · 2026-08-06 performance pass

Profiled a real batch fit (7-site 31P, 16 spectra) that took 15 minutes for
amplitude+δiso and grew releasing width too — three root causes found and
fixed, all in the shared `fit.py`/`batchfit.py` engine so every caller
(interactive fit, batch, Monte-Carlo, χ² profile) benefits:

1. **The errorbar-rescue retry was silently doubling cost.** `fit.py::fit` and
   `multifit.py::fit_cofit` reran the *entire* optimization a second time with
   a different algorithm whenever the primary pass couldn't get a valid
   covariance (an ill-conditioned Jacobian — the exact
   `RuntimeWarning: invalid value encountered in sqrt` a user hit). Fine for
   one interactive fit; a needless doubling for batch's initial pass, every
   Monte-Carlo trial, and every χ²-profile scan point — none of which ever
   read the rescued stderr. New `compute_errorbars` flag, default True
   (unchanged everywhere else), set False at those three call sites.
2. **Batch fit was one N-spectrum x M-parameter joint optimisation solving an
   already-fully-decomposable problem.** `batchfit.batch_fit` always called
   `share=()` and only ever set independent per-spectrum bounds — nothing was
   ever tied across spectra, so a 16-spectrum joint problem was mathematically
   identical to 16 independent ones, at real extra cost (every joint-Jacobian
   evaluation re-simulated every OTHER spectrum, and re-did O(spectra x
   parameters) bookkeeping regardless of which one parameter actually
   changed). Now a plain loop over `fit.fit` per spectrum — same answer (all
   pre-existing batchfit tests, which encode the pre-refactor expected values,
   pass unchanged), far less work, and it automatically inherits fixes #1 and
   #3 since they live in `fit.py`.
3. **The model was simulated on the full spectrum, not the fit window.**
   `fit.py::fit` built its simulation context from the whole experimental
   axis, then only used a windowed slice of the result — on every one of
   thousands of Jacobian-probe evaluations. Now restricted to window + a
   margin (6x the widest declared linewidth) for models that tolerate it
   (`engine.grid_restrictable` — a deliberate ALLOWLIST: pointwise closed-form
   lineshapes and the Czjzek kernel family are safe; `quad_ct`/`quad_first`/
   `quad_csa`/`csa_mas`/`csa_czjzek` are excluded because their LRU-cached
   single-site simulation derives its own span directly from the grid and
   convolution-broadens it — restricting there would truncate real signal, not
   just cost). The returned model curve always still spans the full
   experiment regardless — this is an internal optimisation, invisible to
   callers/the UI.

A synthetic benchmark at a comparable scale (9 spectra x 7 sites,
amplitude+width release, full-precision threshold) completed in ~0.6 s / 1634
evaluations post-fix. Real experimental data (noise, genuine parameter
correlation) will behave differently, but the structural redundancy removed —
re-simulating unaffected spectra, re-optimizing unread error bars, simulating
outside the fit window — doesn't depend on the data.

Also written: `docs/GPU_ACCELERATION_PLAN.md` — an opt-in, auto-detected CuPy
path for users with a CUDA GPU, explicitly gated on profiling actual demand
after the fixes above (not built speculatively).

## Verdict

LARMOR is a broad, coherent, well-tested application: a clean Qt-free core, one
fit engine and one figure engine reused everywhere, honest fit diagnostics, and a
publication path from spectrum to LaTeX. The main debt is the size of `app.py`;
the physics/efficiency foundations (kernel caching, decimation, thresholds,
and — as of 2026-08-06 — no redundant work in the fit engine itself) are sound.

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
