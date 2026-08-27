# Architecture

A map of the codebase for anyone who wants to modify it. Line counts are
approximate as of 2026-08 (about 100 Python modules, roughly 30k lines: ~13k
in the core, ~17k in the desktop layer).

## Core/desktop split

Everything under `larmor/` except `larmor/desktop/` is Qt-free. The core does
all the science — reading, processing, simulating, fitting, quantifying,
rendering figures — and every engine has unit tests that run headless in
seconds. The desktop layer (`larmor/desktop/`, PySide6 + pyqtgraph) is a GUI
over that core and holds no physics of its own. If you are adding a
capability, put the logic in the core with tests first and wire the dialog to
it afterwards. The main window, `desktop/app.py`, is the largest module by far
(~4.6k lines); it is stable but is the place to tread most carefully.

## Module map

Ingestion lives in `larmor/io/`: `bruker.py` (1r/2rr/fid/ser, EXPNO/pdata
layout, pseudo-2D detection), `varian.py` (VnmrJ .fid), `fxmla.py` (dmfit
round-trip), `spectra.py` (CSV/txt), `scan.py` (sample auto-identify), and
`export.py` (fit results to txt/csv/dmfit). `loader.load_any` is the single
entry point; callers never dispatch on format themselves. `fourier.py` handles
States/TPPI/echo-antiecho recombination for 2D acquisitions.

Processing is `processing.py` (apodization, zero-fill, FT, phase, baseline,
linear prediction, applied live or absolute from an unprocessed base) plus
`baseline.py` for the iterative baseline corrector.

Lineshape models live in `larmor/models/` behind a registry (`REGISTRY`,
15 models at the time of writing: gauss_lor, voigt, czjzek, ext_czjzek,
quad_ct, quad_csa, quad_first, csa_mas, csa_czjzek, amorphous, sidebands,
jmultiplet, spectrum, function, gl_norm). Models are nucleus-agnostic; spin
and gyromagnetic ratio are resolved from the isotope symbol at simulation
time. `engine.py` builds and caches the simulation kernels (see performance
below) and `recipe.py` defines the `Recipe`/`SiteModel` structures that
everything passes around.

Fitting is one engine plus wrappers. `fit.py` does the lmfit least-squares
work: constraints and links, bounds, a completion threshold, a windowed
simulation grid, and a `frame_cb` hook for the live fit animation.
`batchfit.py` fits one shared-shape model to many 1D spectra, each spectrum
independently through `fit.fit` (nothing is tied across spectra, so a joint
optimization would be pure overhead). `seqfit.py` warm-starts each spectrum in
a series from its fitted neighbour, forward and backward. `multifit.py`
co-fits several datasets with parameters shared through `share=`.
`autofit.py` provides multi-start fitting plus the error machinery
(chi-squared profiles and Monte-Carlo trials). `twod.py` holds the MQMAS
fit, shear, and 2D kernels.

Interpretation sits beside the fitters: `quantify.py` (integral populations
with errors), `sanity.py` (physical-plausibility flags on eta, widths,
amplitudes, window coverage), `identifiability.py` (flags parameter pairs
with |r| >= 0.95), `diagnostics.py` (residual runs-test and autocorrelation),
`chi2map.py` (parameter-pair chi-squared surfaces), and `recipe_diff.py`
(current fit versus a reference).

Output is `figures.py` — `render()` on plain-dict specs for 1D/2D/series
figures with journal presets and png/pdf/svg/tiff export — together with
`methods.py` (LaTeX table and methods sentence) and `batch.py` (publication
table across fits).

## Data flow

```
load_any ──▶ (ppm, amp, recipe) ──▶ processing ──▶ workbench
    │                                                │
    ▼                                     place sites (models registry)
 NMRData (2D/fid) ──▶ 2D view / FID dialog           │
                                                     ▼
                                     fit.py ◀── FitWorker (thread, animate)
                                        │           ▲ same engine reused by
                                        │           batchfit / seqfit / multifit
                                        ▼
                       quantify · sanity · identifiability · diagnostics · chi2
                                        │
                                        ▼
                        figures / methods / batch ──▶ publication bundle
```

A file becomes an axis, an amplitude array, and a `Recipe`; the recipe
accumulates processing steps and fitted sites; the fitted recipe feeds the
diagnostics and the figure/report path. Recipes serialize to `.recipe.json`,
so any fit can be reloaded, diffed, or re-rendered later.

## Integration invariants

Four rules keep the app coherent, and contributions should preserve them.
There is one fit engine: `fit.fit` is called by the main window, batch,
sequential, co-fit, chi-squared map, and Monte-Carlo, so thresholds,
constraints, and error handling behave identically everywhere. There is one
figure engine: `figures.render` on dict specs backs the Plotting studio, the
batch report, and the publication bundle, and specs are saveable. Every plot
gets the same export dialog (`export_dialog`) and right-click menu
(`plot_menu`); the crash-prone native pyqtgraph exporter is disabled
globally. Cross-session UI state (theme, font, recents, pinned folders, fit
threshold, per-nucleus site defaults, constraint library, batch templates)
goes through QSettings, never ad-hoc files.

## Performance design

The expensive object is the Czjzek basis: `engine.build_kernel` simulates the
(Cq, eta) subspectrum grid once with mrsimulator and caches it per process,
keyed by nucleus, Larmor frequency, spinning rate, and spectral window plus
the grid resolution. This is why a Czjzek fit takes seconds rather than
minutes; call `clear_kernel_cache()` if you change kernel settings. Analytic
single-site models use a parameter-level LRU instead.

Heavy imports (mrsimulator, matplotlib) happen lazily inside the functions
that need them, which keeps startup fast — do not hoist them to module level
in the core or the desktop dialogs. During a fit, the model is simulated only
on the fit window plus a margin for models on the `grid_restrictable`
allowlist (an allowlist deliberately, since some cached models derive their
span from the grid), the live-animation redraw is throttled to every
`frame_every` iterations, and error-bar recovery can be skipped via
`compute_errorbars=False` where the stderr is never read (batch initial
passes, Monte-Carlo trials, chi-squared scans). For 2D display,
`twod_view._decimate` caps contour grids at 512 points per axis. Error
analysis (chi-squared profiles, Monte-Carlo) can fan out across a process
pool through `parallel.parallel_map`; it stays sequential by default and for
small workloads, and results are bit-identical either way.
