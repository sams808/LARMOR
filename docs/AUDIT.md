# LARMOR — full app audit (2026-08-05)

A pass over every feature area, how they fit together, and where the efficiency
and maintainability risks are. Written after the 2026-08 feature waves.

## 1 · Scale

| | |
|---|---|
| Python modules | ~105 (`larmor/` core + `larmor/desktop/` UI) |
| Lines of code | ~25 k (≈11 k core, ≈14 k desktop) |
| Test files | 56 · full suite **green** |
| Biggest module | `desktop/app.py` — 4.2 k LOC, 216 methods (the main window) |

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
covariance retry, completion threshold, **live-animation frame_cb**),
`multifit.py` (co-fit 1D+2D shared params), `batchfit.py` (one shared model, many
1D), `seqfit.py` (forward-backward warm-start series), `autofit.py` (multi-start +
MC errors), `twod.py` (MQMAS fit/shear/kernels).

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
- **`np.trapz` deprecation** (numpy 2.0 → `np.trapezoid`) in `figures.py` and
  `measure.py` — cosmetic warning, worth a one-line fix.

## 5 · Recommendations (priority order)

1. **Refactor `app.py`** into mixins (biggest maintainability win; no behaviour change).
2. Silence the `np.trapz` deprecation (one line each in `figures.py`, `measure.py`).
3. Harden `chi2_surface` with a max-time / coarse-refine guard.
4. Remaining *ideas-from-field-use* not yet built (all niche/interaction-heavy):
   **#14** labelable peak-pick assignment table, **#12** TopSpin drag-to-phase
   gesture, **#19** dmfit amplitude-calibration UI.
5. Add a smoke test that constructs every dialog headlessly (catches wiring
   regressions cheaply given how much UI now exists).

## Verdict

LARMOR is a broad, coherent, well-tested application: a clean Qt-free core, one
fit engine and one figure engine reused everywhere, honest fit diagnostics, and a
publication path from spectrum to LaTeX. The main debt is the size of `app.py`;
the physics/efficiency foundations (kernel caching, decimation, thresholds) are
sound.

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
