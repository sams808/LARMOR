# LARMOR — full app audit (2026-08-05, updated 2026-08-06)

A pass over every feature area, how they fit together, and where the efficiency
and maintainability risks are. Written after the 2026-08 feature waves; updated
after the 2026-08-06 performance pass (see §6).

## 1 · Scale

| | |
|---|---|
| Python modules | 106 (`larmor/` core + `larmor/desktop/` UI) |
| Lines of code | ~28 k (≈8.6 k core, ≈15 k desktop, rest tests/scripts) |
| Test files | 60 · full suite **green** (508 passed, 16 skipped, 2026-08-06) |
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

## 7 · 2026-08-06 (evening): publication-plotting & batch-fit generality pass

A code-generality survey (15 registered models, `twod.py`, `figures.py`,
`plotting_studio.py`, real-data test coverage) found the **model registry and
batch-fit engine are already genuinely nucleus-agnostic** — no hardcoded
isotope strings anywhere; every quadrupolar/CSA model resolves gyromagnetic
ratio/spin generically via mrsimulator's `Isotope(symbol=nucleus)`. The real
gap for "any technique" was on the **2D side**, addressed this pass:

- **Batch fit: per-spectrum "Exclude component"** — right-click a cell to
  lock one site's amplitude to exactly zero for that spectrum only
  (`batchfit.is_zeroed_out`); excluded sites are skipped in the CSV export
  and never drawn/legend-listed when rendered, rather than reported as "a
  fitted zero". `shared_table`/`error_table` also now export each site's
  integrated **population %** (`larmor.quantify`, error-method-consistent),
  and a CSV export can auto-save each spectrum's `.recipe.json` alongside it.
- **Plotting studio**: per-component colors + legend visibility, an
  Auto-update toggle (off by default) + Preview button, and — the 2D side —
  a fitted-model contour overlay on a real experimental map plus
  **computed** CS-axis/QIS-axis reference lines (`twod.f1_cs_scale`/
  `qis_slope`, previously only hand-typed).
- **series_grid.recipe_from_csv_rows correctness fix** (found by a
  deliberately mixed-model test — gauss_lor + czjzek + quad_ct in one
  recipe): a site missing a parameter that has a soft render-time default
  (czjzek's `line_fwhm_ppm`) was wrongly rejected as "incomplete" instead of
  filled from the model registry's own default — the exact "worked for one
  dataset, silently wrong for another" class of bug this pass was meant to
  catch. Also fixed: `fit.py::_make_params` crashed on any FIXED parameter
  with `min == max` (lmfit itself rejects degenerate bounds even for
  `vary=False`) — bounds are now cleared for fixed parameters, since they
  have zero effect on optimisation.
- **New real-data regression**: `test_twod.py::test_fit_2d_on_real_mqmas_data_gives_a_sane_fit`
  — every prior `fit_2d` test fit a self-generated synthetic "experiment"
  (a pure round-trip); this fits a real 27Al 3QMAS 2rr and checks for a
  genuinely converged result, closing a real blind spot the survey found.
- **Explicit 2D scope boundary** (unchanged, not a gap introduced here, but
  now documented so it isn't rediscovered): 2D fitting covers 4 models
  (czjzek, ext_czjzek, quad_ct, quad_csa) and MQMAS-family methods only
  (`twod.METHODS`); HMQC/DQ-SQ are visualization-only by design
  (`larmor.correlate` is deliberately unwired, per its own docstring — 2D
  batch fitting, per-dimension 2D lineshape, and the generalized
  multi-experiment correlation engine remain roadmap items, not attempted
  in this pass).
- Also shipped: 4 hidden "aesthetic" themes (Y2K/Dreamcore/Gen X Soft
  Club/Vaporwave — `theme.AESTHETIC_THEMES`, applied live, View ▸ Theme ▸
  More styles…) — a just-for-fun addon, same contrast-floor tests as the
  10 normal presets, kept out of the normal Theme menu.

## 8 · 2026-08-07: exclude-component regression (real-data catch)

Reported against a live 17-spectrum/8-line batch fit (PBi 31P): spectra with
a component excluded ("P0/P1(Bi contact)") showed markedly worse RMSD than
the rest, and the whole batch fit was slower than before the exclude feature
existed — both symptoms traced to one bug in `batchfit.batch_fit()`.

- **Root cause**: the "Exclude component" feature (§7, A1) locks only the
  excluded site's *amplitude* to zero for that spectrum. `batch_fit()`'s
  per-parameter release loop, though, decides `vary` **by parameter name
  across the whole model**, uniformly for every site — it never checked
  whether a given site was excluded before honouring "Release per spectrum"
  for that site's position/width. With δiso and FWHM released (the reported
  configuration), an excluded site's position/width were *still* set
  `vary=True` with relaxed bounds, even though amplitude=0 makes them
  provably inert (zero gradient w.r.t. the residual). Two dead-gradient free
  parameters per excluded spectrum is wasted work at best (slower fits) and,
  in practice, degraded the optimiser's convergence on the real parameters
  too (worse RMSD) — exactly the two symptoms reported.
- **Fix** (`larmor/batchfit.py`, `batch_fit()`): a site whose amplitude is
  `is_zeroed_out` now has *every* one of its other parameters forced
  `vary=False` for that spectrum, unconditionally — before the `released`
  check runs, not after. An excluded site is now fully inert everywhere,
  not just in amplitude.
- **New regression test**
  (`test_batchfit.py::test_excluded_site_params_stay_fixed_even_when_that_param_is_released`):
  reproduces the exact real-world combination (exclusion + a globally
  released position/width), asserts the excluded site's params stay fixed
  while the same params stay released on every other spectrum and on every
  other (non-excluded) site, and asserts the excluded spectrum's RMSD
  exactly matches a control fit where that site never existed in the model
  at all — proving the extra dead parameters were the entire discrepancy.
- Full suite green after the fix (`pytest tests/ -q`).

## 9 · 2026-08-07: CPU-parallel error analysis (Monte-Carlo / χ² profile)

User report: a real 17-spectrum/8-line batch fit's "Compute errors" (χ²
profile, 40 points) had been running 30+ minutes with no end in sight, and
Monte-Carlo (200 trials) was similarly slow. Diagnosis before touching any
code: with "Release per spectrum" on for 2 parameters across 8 sites, that's
24 free parameters/spectrum — a χ² profile is 40 refits PER free parameter,
so one spectrum alone is 960 independent `lmfit` fits, ×17 spectra ≈ 16,000
total; Monte-Carlo is 200 independent refits × 17. Every one of those is a
small, fully independent nonlinear fit (a few thousand points, no shared
state) run **one at a time** in a plain Python loop — the textbook shape for
OS-level process parallelism, and specifically NOT a fit for a GPU (each
unit of work is `lmfit`'s serial per-fit control flow, not a large batched
tensor op the small array sizes here would ever saturate). User chose
"multiprocessing first" when offered the choice explicitly.

- **NEW `larmor/parallel.py`** (Qt-free, tested in isolation first):
  `parallel_map(fn, items, ...)` — sequential by default (`use_processes`
  off, or below `MIN_ITEMS_FOR_PROCESSES=8` items, avoids paying process-pool
  startup cost for nothing); across a `ProcessPoolExecutor` otherwise.
  Preserves original item order in the returned list regardless of
  completion order; `on_result(i, r)` fires per completion for live progress;
  `should_stop()` cancels not-yet-started work (documented as BEST EFFORT —
  `ProcessPoolExecutor` prefetches into its own call queue ahead of a worker
  picking an item up, so exactly how much gets cut off isn't controllable,
  proven by a test using a real per-item delay rather than asserting an
  exact cutoff index); `executor=` lets a caller reuse one already-running
  pool across many calls instead of paying startup cost per call.
- **`larmor/autofit.py`**: `error_profile()` and `monte_carlo_errors()` each
  gained `parallel=False` (default — every existing caller, including every
  test, is unchanged), `max_workers`, `executor`. Both extracted their
  per-item work into a module-level (picklable) worker function
  (`_profile_point_worker`, `_mc_trial_worker`) run via `parallel_map`.
  Monte-Carlo's synthetic noise draws stay a single sequential pass over
  `rng` **before** dispatch, so the trial set — and the whole result, for a
  fixed seed — is identical regardless of how trials get scheduled across
  workers (verified by a dedicated bit-for-bit parallel-vs-sequential test).
- **`larmor/batchfit.py`**: `batch_error_analysis()` gained `parallel=False`,
  `max_workers`. When on, ONE `ProcessPoolExecutor` is created up front and
  reused (via `executor=`) across every spectrum × parameter in the run,
  not a fresh pool per `autofit` call — with a many-spectrum, many-released-
  parameter batch that's potentially hundreds of `error_profile` calls, and
  pool startup has real cost (a fresh interpreter per worker on Windows).
- **GUI opt-in** (library default stays off/unchanged): the batch-fit
  dialog's `_ErrorWorker`, the standalone Monte-Carlo dialog, the standalone
  χ² profile ("Errors Analysis") tool dialog, and the batch-fit-report
  tool's Monte-Carlo call all now pass `parallel=True`.
- **PyInstaller correctness fix, found by reasoning through the packaging
  path, not by reproducing it**: a frozen Windows exe that spawns worker
  processes without `multiprocessing.freeze_support()` at its entry point
  has each worker re-execute the frozen entry point from scratch and launch
  a full second app instance instead of becoming a worker — recursively.
  Added to `packaging/launcher.py`'s `__main__` guard, the only frozen
  entry point per the `.spec` file. A no-op for the normal `larmor desktop`
  /pytest paths (this only matters once LARMOR is rebuilt into an exe).
- **Warning spam fix** (surfaced by the user mid-run: hundreds of lmfit
  `RuntimeWarning: invalid value encountered in sqrt`/`scalar divide` lines).
  Root cause: lmfit computes a covariance-based stderr for every free
  parameter INSIDE `minimize()` itself, unconditionally — with most
  parameters fixed at each χ² scan point (routine there, and common anyway
  for correlated quadrupolar sites), that covariance is often not
  positive-definite, and lmfit `sqrt()`s the diagonal regardless. Harmless
  (this exact stderr is never read in these contexts) but at thousands of
  fits per run, pure noise. Silenced narrowly — scoped inside a
  `warnings.catch_warnings()` block around ONLY the `lmfit.minimize()` call,
  by exact message text, so an unrelated `RuntimeWarning` still surfaces —
  in all three of this codebase's independent `lmfit.minimize()` call sites
  (`fit.py`, `twod.py`, `multifit.py`; found the latter two by re-running
  the full suite after the `fit.py` fix and noticing the warning hadn't
  fully disappeared, rather than assuming one fix covered every call site).
- **New tests**: `tests/test_parallel.py` (8, the module in isolation —
  ordering, exception-as-hole, on_result firing, best-effort should_stop,
  shared-executor reuse, worker-count floor); parallel-vs-sequential
  bit-identical-result tests added to `test_montecarlo.py`, `test_autofit.py`,
  and `test_batchfit.py` (covering both `error_profile` and
  `monte_carlo_errors`, standalone AND through `batch_error_analysis`'s
  shared-pool path). 520 passed, 16 skipped, full suite green.
- **Deliberately not done this pass**: GPU acceleration. Not drafted — the
  user chose "multiprocessing first, GPU only if still needed" when the
  trade-off was raised, and CPU parallelism plausibly closes the entire
  reported gap (a 30+ minute run → minutes) for this workload shape. If a
  genuinely GPU-shaped workload shows up later (e.g. batching the MQMAS 2D
  kernel build's crystallite-orientation sum, or fitting a MUCH larger
  N-spectra batch than seen so far), that would be a separate, from-scratch
  design — a process pool of `lmfit` calls doesn't partially become "a GPU
  version," they're different architectures end to end.

## 10 · 2026-08-07: batch-grid component list undercounted on an excluded panel

Real-data catch (user report + screenshot): the Plotting Studio's "Component
colors / legend" dialog listed 7 of a real 8-site model's components while
the ACTUAL rendered panel correctly drew all 8. Traced to two compounding
issues, both fixed:

- **Root cause, in `series_grid.recipe_from_csv_rows`**: a site excluded for
  one spectrum (batch fit's "Exclude component") was being DROPPED from
  that panel's reconstructed `Recipe.sites` entirely, rather than kept
  present with amplitude locked to zero — unlike a real saved `.recipe.json`
  fit, where an excluded site always stays in the list (that asymmetry was
  introduced two sessions ago while fixing a DIFFERENT bug — "don't treat
  an excluded site as incomplete" — and the chosen fix, silently omitting
  it, created this one). Dropping it doesn't just shrink that one panel's
  count: it shifts every LATER site's **list index** on that panel only,
  while every other (non-excluded) panel keeps the original alignment. Any
  spec keyed by index — `component_colors`, `shade_only`,
  `hide_components`, `legend_hide` — silently means a DIFFERENT physical
  site depending on which panel it's applied to. This is a strictly worse
  failure mode than the visible symptom (an undercounted legend): a color
  picked for "row 5" in the dialog could paint the wrong component once
  applied to any excluded panel. Fixed by reconstructing an excluded site
  as a present-but-zeroed placeholder at its normal sorted position
  (matching `batchfit.is_zeroed_out`'s exact signature), so every panel of
  a batch keeps the identical site list/order regardless of exclusions —
  `render_batch_grid`'s existing `is_zeroed_out` check (already there for
  real saved fits) then hides it for that one panel exactly as before, no
  change needed in the renderer itself.
- **Second layer, in `plotting_studio._grid_detect_sites`**: it only ever
  looked at the FIRST resolvable panel to populate the color/legend dialog,
  on the (now only conditionally true) assumption that any panel's site
  list represents the whole model. Made it union the site list across
  EVERY resolvable panel instead — defense in depth for any other source of
  a genuinely shorter list (e.g. a hand-built/legacy recipe missing a
  TRAILING site), independent of the `series_grid` fix above.
- Both fixed at the correct layer rather than patched at the symptom: the
  studio's own component/legend detection is defense-in-depth, but the
  index-alignment guarantee had to be restored at the reconstruction
  source, since that's the only place that can actually keep every panel's
  site list in sync.
- **New/updated tests**: `test_series_grid.py`'s
  `test_recipe_from_csv_rows_drops_a_site_absent_from_this_scope` renamed
  to `..._zeroes_rather_than_drops_...` and rewritten to assert the site is
  present-and-zeroed, not absent; `test_load_panels_via_csv_with_excluded_
  site_reconstructs_the_rest` extended to assert BOTH panels keep the same
  2-site list and to render the actual batch-grid figure from both;
  `test_plotting_studio.py` gained
  `test_grid_detect_sites_unions_across_panels_not_just_the_first`. 521
  passed, 16 skipped, full suite green.

## Verdict

LARMOR is a broad, coherent, well-tested application: a clean Qt-free core, one
fit engine and one figure engine reused everywhere, honest fit diagnostics, and a
publication path from spectrum to LaTeX. The main debt is the size of `app.py`;
the physics/efficiency foundations (kernel caching, decimation, thresholds,
and — as of 2026-08-06 — no redundant work in the fit engine itself) are sound.

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
