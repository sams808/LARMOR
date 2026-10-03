# Workplan status — 2026-09-23 session handoff

Where the 50-item improvement workplan stands, what the sessions after it
changed, and what a continuing session picks up. The workplan itself (items,
evidence, Sam's keep/defer/drop decisions) lives in the artifact "LARMOR
Workplan"; this file is the ground truth for progress.

## State of the repository

- **Version**: 0.15.0 (`larmor/__init__.py` + `pyproject.toml`, bumped in step;
  committed 2026-09-30 as f627a89, **not pushed** — no push word this session).
- **Pushed** to `origin/master` on 2026-09-23 on Sam's go-ahead ("once done
  make sure the version on github is up to date, and push if needed"):
  v0.13.0 → v0.14.0, the eleven merged branches of the next-ten batch below.
- **Tests**: 1774 collected in 140 files at v0.15.3 (1761 in 139 at v0.15.2, 1728 in 135 at v0.15.0, 1324 in 107 at v0.14.0). Last full run: see the
  line "Full suite" at the end of this section. Trust a green bar only when
  the terminal banner says "real-data layer: complete (all 19 datasets
  present)" — five datasets were added this session (the CaF₂ / NaF magres
  files and their 2026-03 ¹⁹F standards, and the 2026-05 ³¹P EXPNO 3102).
- **Environment**: conda env `larmor` (`C:\Users\samso\xraylarch\envs\larmor`,
  Python 3.11); run tests with that interpreter, `PYTHONIOENCODING=utf-8`
  (several CLI summaries print `→`, which a cp1252 console cannot encode).
  Real test data roots at `LARMOR_TEST_DATA` (default `C:\Users\samso`).
- **Installer**: `dist/LARMOR/` (351 MB) rebuilt at 0.14.0 on 2026-09-23
  from the pip venv at `packaging/.buildenv` (python.org 3.11 — never the
  conda env; reasons in `packaging/README.md`), then wrapped two ways by
  `packaging/build.bat`: **`dist/LARMOR-0.14.0-setup.exe`** (102 MB, Inno
  Setup 6.7.3 from `packaging/larmor.iss`; per-user install without
  administrator rights, Start-menu entry, optional desktop icon, uninstaller
  under Settings ▸ Apps) and `dist/LARMOR-0.14.0-win64.zip` (154 MB, unzip
  and run), both with a plain-language `INSTALL.txt` next to `LARMOR.exe`.
  Verified: the exe smoke-run offscreen 45 s with an empty crash log; the
  archive viewer lists every new module, the 13 manuals and the 7 tutorials;
  the setup installed silently into a scratch folder in 31 s (1952 files),
  the installed app ran 30 s, the uninstall registration read "LARMOR
  0.14.0", and the silent uninstall left no folder and no registry key.
  Still owed: a run on a machine with no development setup (a student's
  laptop), then a GitHub release carrying both files. Rebuilt at **0.14.2**
  on 2026-09-24 (`dist/LARMOR-0.14.2-setup.exe` 105 MB,
  `dist/LARMOR-0.14.2-win64.zip` 151 MB): the frozen exe passes
  `--selftest` (core and desktop fits), the setup installed silently in
  21 s, the installed exe passed the self-test, the Apps entry read
  "LARMOR 0.14.2", and the silent uninstall left no folder and no key. The
  0.14.0 and 0.14.1 files stay in `dist/` for comparison; the student's
  machine should get the 0.14.2 setup and run `LARMOR.exe --selftest` if
  Fit misbehaves again. **0.14.3** rebuilt the same day
  (`dist/LARMOR-0.14.3-setup.exe` 105 MB, `dist/LARMOR-0.14.3-win64.zip`
  151 MB): frozen exe `--selftest` pass, silent install 40 s (1993 files),
  installed exe self-test pass, Apps entry "LARMOR 0.14.3", silent uninstall
  clean. **0.14.4** built after the review fixes
  (`dist/LARMOR-0.14.4-setup.exe` 105 MB, `dist/LARMOR-0.14.4-win64.zip`
  151 MB): frozen exe `--selftest` pass, silent install 35 s (1993 files),
  installed exe self-test pass, Apps entry "LARMOR 0.14.4", silent
  uninstall clean. **0.15.0** built 2026-09-30 after the full suite
  (`dist/LARMOR-0.15.0-setup.exe` 111 MB, `dist/LARMOR-0.15.0-win64.zip`
  159 MB; `packaging/build.bat` run from PowerShell — its version read
  had to be repaired first, see 10eb5d1): frozen exe `--selftest` PASS with
  "solver: ran" on all four fits, also with the completion threshold set to
  the 0.1 % default a fresh install gets (151 / 281 evaluations, RMSD
  0.0524 as at full precision); silent install 38 s (1994 files), Apps
  entry "LARMOR 0.15.0", installed exe self-test PASS, silent uninstall 7 s
  with no folder and no entry left; crash log empty. **0.15.1** built the
  same day after its full suite (`dist/LARMOR-0.15.1-setup.exe` 111 MB,
  `dist/LARMOR-0.15.1-win64.zip` 159 MB): frozen exe `--selftest` PASS
  (core czjzek 557 evaluations with the scaled solver), silent install 39 s
  (1994 files), Apps entry "LARMOR 0.15.1", installed exe self-test PASS,
  silent uninstall clean, crash log empty. `packaging/build.bat` now ends
  with step [6/6], `packaging/desktop_shortcut.ps1`: Sam's standing request
  ("always replace the shortcut to the larmor.exe on the desktop") after a
  shortcut to PyInstaller's intermediate `build\larmor\LARMOR.exe` died with
  "Failed to load Python DLL …\build\larmor\_internal\python311.dll" — the
  step rewrites `Desktop\LARMOR.lnk` to `dist\LARMOR\LARMOR.exe`. **0.15.2**
  built 2026-10-02 after its full suite (`dist/LARMOR-0.15.2-setup.exe`
  111 MB, `dist/LARMOR-0.15.2-win64.zip` 159 MB; build.bat's step [6/6]
  refreshed the desktop shortcut): frozen exe `--selftest` PASS, silent
  install 33 s (1994 files), Apps entry "LARMOR 0.15.2", installed exe
  self-test PASS, silent uninstall 6 s clean, the desktop shortcut restored
  afterwards (the setup's own desktop icon shares its name), crash log
  empty. **Hand out the 0.15.2 setup**; every earlier file in `dist/` is
  superseded.
- **Full suite**: at v0.14.0 (0de34c7), **1324 passed / 0 failed** in 20 min 39 s,
  real-data layer complete (all 19 datasets present). v0.14.1 (faab8f9) adds
  `larmor/selftest.py`, the quick overlays and their tests (scoped runs
  green; the full suite was not re-run for it). **v0.14.2 (cc8b1ca): 1366
  passed / 0 failed** in 20 min 21 s, real-data layer complete (19 datasets).
  **v0.14.3: 1561 passed / 0 failed** in 34 min 12 s (1561 tests in 116
  files), real-data layer complete. **v0.14.4: 1570 passed / 0 failed**,
  real-data layer complete — the run reports 8 h 53 min of WALL time
  because the machine slept from 23:34 to 08:00 (Kernel-Power 42/107);
  the compute was ~27 min, and a before/after benchmark of the plot paths
  the review touched (set_experiment, set_overlays, setXRange, zoom_full,
  set_model, Y-mode round trip) shows no regression: every one is within
  a millisecond of 0.14.3 except zoom_full, 0.97 → 2.15 ms.
  **v0.15.0 (f627a89, 2026-09-30): 1728 passed / 0 failed** in 23 min 31 s,
  run alone after the four worktree merges, real-data layer complete (all 19
  datasets present). **v0.15.1 (a6d2251): 1740 passed / 1 failed** in
  41 min 03 s, real-data layer complete — the one failure was the fixture
  precondition of `test_error_analysis_covariance_refits_when_the_fast_fit_had_none`
  ("all six stderr None"), which the Jacobian-scaled solver no longer
  satisfies for every spectrum of the degenerate pair (2 of 6 None); the
  assertion now requires at least one spectrum without error bars, and the
  test passes. **v0.15.2 (8bbc39d, 2026-10-02): 1761 passed / 0 failed** in
  52 min 08 s, run alone, real-data layer complete (all 19 datasets
  present). **v0.15.3 (2026-10-02): 1774 passed / 0 failed** in 26 min 56 s,
  run alone, real-data layer complete (all 19 datasets present).
- **0.14.1 (2026-09-23, same day)**: a student's frozen 0.14.0 opened spectra
  but "fitting did not work". The fit path was verified in the exe's own
  package set — console-less Python of `packaging/.buildenv`, then the frozen
  exe itself — on synthetic and on a real 32k-point ²⁷Al spectrum (Gauss/
  Lorentz 1.7 s, Czjzek 18 s with the kernel build), so the cause is on that
  machine. **Resolved in 0.15.0 — it was not that machine**: the completion
  threshold's callback abort (b70e9a7, in every 0.14.x) ended every
  Fit-button fit after ~5 evaluations at the default 0.1 %. The verification
  and the self-test below judged a fit by "a parameter changed by more than
  1e-9", which the amplitude pre-scale and the first Jacobian probe satisfy
  even when the solver is aborted at once — so they could not see it; since
  0.15.0 the self-test also requires that the lmfit result was not aborted
  and ran past its first Jacobian sweep (`selftest.solver_verdict`). See
  the 2026-09-30 section. Shipped
  then: **`LARMOR.exe --selftest [spectrum]`** fits
  through the core and through the window's own Fit path and writes
  `%USERPROFILE%\LARMOR_selftest.log` (every step, timing, traceback, dialogs
  the user may not have seen); the spec now collects every larmor / lmfit /
  scipy.optimize / scipy.stats submodule and the dist-info metadata. Lesson
  recorded in `packaging/README.md`: a build cannot replace `dist\LARMOR`
  while a copy runs — PyInstaller ends with the OLD exe in place; check the
  timestamp, or build with `--distpath dist_test` and wrap it with
  `/DSourceDir=..\dist_test\LARMOR`. Also in 0.14.1, on Sam's request: quick
  overlays in the main window (Shift + drop, File ▸ Overlay a spectrum…
  Ctrl+Shift+A, Explorer right-click, View ▸ Overlays Ctrl+Shift+V / Clear
  overlays, "match height" in the Datasets dock).
- **0.14.2 (2026-09-24)**: Sam's fourteen field-use items after the first
  student session, built on three worktree branches and merged (WB 4afa9a2,
  WC 446fcad, WA 2c92c46; release commit cc8b1ca). Fixed: the Baseline
  iterative dialog crashed on open (`PlotItem.getPlotItem`) and its Apply
  crashed after closing (`dlg.Accepted` on the instance); the spectrum view
  "drifted" to thousands of ppm because pyqtgraph's auto-range unioned the
  model drawn on the Czjzek kernel axis (≥ 150 kHz) *and* the "fitting —
  iter n" label re-pinned to the view corner every frame (a runaway loop);
  every model-side item is now `ignoreBounds` and masked to the data range,
  and a fit never moves the view. Speed: a czjzek_corr fit with five linked
  sidebands went from a 40-minute stall (1 evaluation) to 5.7 s (84
  evaluations, error bars) — scipy's first trust-region step scaled by the
  amplitude (~3·10⁶) sent dCS to 3·10⁶ ppm and a 2.4-million-point Gaussian
  filter, and the error-bar rescue dithered at lmfit's default `epsfcn`;
  per-evaluation renders shared between linked copies (czjzek 37.6 → 2.9 ms);
  Stop returns within a second even inside a kernel build (cancel checked
  every ≤ 96 tensors and between site renders). Adding a line takes 15–50 ms;
  the "seconds" were the cold kernel build in the SimWorker, now announced
  with a busy cursor and built once for both threads. Linked sidebands carry
  `SiteModel.sideband` (parent, order) and their exprs are recomputed when
  the spin rate changes. Fit health: red only for impossible values (sanity
  rules incl. C_Q > 120 MHz), amber "check" for residual structure, bounds,
  degenerate pairs, no chip at all for an unknown T1 or tip angle; lb of the
  Czjzek family held by default (`ParamDef.default_fixed`, dmfit's greyed
  Lb). UI: the "＋ Add line" split button with grouped models
  (`desktop/addline_toolbar.py`), a line right-click menu (sidebands,
  duplicate, fix/free, remove the selection) with multi-select Remove, the
  menu bar reorganised into File / Edit / Process / Fit / Series / Tools /
  View / Plotting / Help (214 golden rows, every doc path rewritten, a Menu
  map in getting-started.md), the Processing panel in five titled groups
  that fit 1080 px without scrolling, Explorer ▸ Rename… (alias kept by
  LARMOR in `%LOCALAPPDATA%/LARMOR/aliases.json`, or the folder renamed on
  disk after confirmation, logged to `rename_log.jsonl`). Follow-ups:
  nmrglue `read_jcamp` loops forever on a truncated JCAMP array (a corrupt
  acqus would hang the app — guard it); the `LinesTable` embedded in
  `panels.SiteCard` / the sequential-fit dialog is not wired for multi-select
  Remove.
- **0.14.3 (2026-09-24, later the same day)**: Sam's second field-use batch,
  four branches (WD 490f31b…38bbf93, WE 0a13a3a…59d1d0b, WF 150e103…41ac5f5,
  WG on master 76d1c6b/e93b364; release commit 37ece2d). **Help and
  viewers**: manuals and tutorials open as non-modal tool windows
  (`desktop/windowtray.show_tool_window`) with an A− / A+ / reset bar
  (Ctrl + wheel, size in QSettings) — `QTextDocument.toHtml()` had been
  pinning the manuals at the 9-pt app font all along; thirteen read-only
  viewers and calculators (NMR table, conversions, correlations, Czjzek
  distribution, compare fits, χ² map, integrals, dataset info, About,
  comparability, series plot…) open the same way; value-returning dialogs
  and the analysis tools that push results from a snapshot stay modal, on
  purpose. **Open-windows bar** at the bottom of the left sidebar: every
  tool window listed oldest-first with a count, minimised ones hidden from
  the desktop and dimmed in the bar (no floating title bar), right-click
  Restore / Minimise / Close / Close all. **Plot**: the Full view is
  computed from the data (x extent, y with room for the residual strip),
  a range requested entirely outside the data or wider than ~3 spans snaps
  back, pan limits ±1 span; **View ▸ Y axis**: raw / normalised to maximum /
  to area / to the area of a region, ONE display transform in SpectrumView
  (`larmor/display.py` Qt-free) that scales experiment, model, components,
  residual, frames, paddles (drags map back to raw), anchors, overlays (own
  max / area) and the y label, while recipe, fit and exports stay raw.
  **Datasets dock**: per-overlay ×scale, shift (ppm) and ↑offset rows,
  "match height" writes the factor, right-click Reset / Make active /
  Remove; the three keys travel with the session and the project bundle.
  **QCPMG**: Save as dataset writes `<base>_sumecho.csv` and
  `<base>_spikelets.csv` (own axes, unit maximum, `intensity_scale`,
  `twin_file`, `spectrum_kind` in the headers; `qcpmg.dataset_pair_paths`,
  `write_dataset_pair`), a "normalise to max" checkbox, phasing collapsed
  under "Phasing (optional)", magnitude the documented route (manual +
  tutorial 7). **Themes**: canonical Solarized and Nord palettes, distinct
  identities for Sepia / Slate / Ocean / Y2K (now light silver + lime) /
  Dreamcore / Gen X Soft Club, a `scheme` on every Theme, CIE Lab
  discernibility test over all 91 pairs (window ΔE ≥ 6 or accent ΔE ≥ 25),
  swatch icons and scheme tooltips in the menus, the plot follows the theme.
  **Axis prefix**: `desktop/axes.plain_units()` on every physical axis of
  every dialog plot (23 creation sites) — and the actual trigger of the
  "6 … −6 ppm" report: pyqtgraph recomputes the SI prefix inside
  `enableAutoSIPrefix(False)` and never again once off, so re-applying it
  from `apply_theme()` while a ±6000 ppm spectrum was shown froze a "k"
  prefix on the main plot; `plain_units` also resets the frozen scale.
  Follow-ups: `measure.centre_of_mass` raises on an empty region mask (an
  Integrals region dragged to zero width); no per-theme site palette for
  Sepia / Ocean / Slate; the magnitude checkbox is not ticked by default
  (a pinned test expects absorption on load) — a two-line change if wanted.
- **0.14.4 (2026-09-24, the review of 0.14.3)**: because four agents wrote
  0.14.3 in parallel and only scoped tests had seen it, the whole diff went
  through a six-lens adversarial review (Qt lifetime, the display transform,
  the view range, QCPMG data integrity, persistence, merge damage and tests
  that cannot fail), each finding verified by a skeptic told to refute it:
  **29 raised, 15 confirmed, 14 refuted**, deduplicating to seven defects —
  four lenses had independently found the first one. All fixed, each with a
  test checked to fail against the unfixed code:
  **(A)** "match height" stored a ratio of RAW maxima while
  `display.overlay_display` had already normalised each overlay by its own
  trace, so under *Normalise to maximum* an overlay was drawn
  max(active)/max(overlay) times too tall (20× measured); the stored factor
  is now display-space (`display.match_scale_display`), so ×1 always means
  "as tall as the active spectrum". **(G)** the same factor went stale after
  *Make active*, a new spectrum or a mode change; a ticked box now re-derives
  in `_refresh_overlays` (the one path they all funnel through) and an untick,
  a typed value or Reset stops it. **(E, high)** `set_experiment` assigned
  `_y_scale` directly, so the baseline anchors and paddles already on the plot
  kept the old display units: an anchor clicked at raw 2.0 read back 3.96
  after loading a spectrum twice as tall, and the manual baseline subtracted
  that. `_set_y_scale()` is now the only assignment and rescales both.
  **(D)** *View all* framed only the active trace, dropping the compared
  spectra that the same batch had made scalable and shiftable (before the
  batch it auto-ranged and included them); `full_extents(..., extra=)` unions
  the visible overlays as drawn, still excluding the model and far paddles.
  **(C)** the Datasets ×-box clamped to [0.01, 1000] while the dict and the
  plot used the real factor, so it showed a different number and one arrow
  click destroyed it; the range is now 0 … 1e9 with adaptive decimals and
  `_set_quiet` widens the box for anything outside it. **(F)** a Y-mode change
  with the FID on screen lost the spectrum's saved zoom (the stored range was
  in the old display units). **(B)** the intensity axis was the one physical
  axis left auto-SI-prefixed: under *Normalise to area* a 0 … 6e-8 axis read
  0 … 60 with an "n" prefix, and under *maximum* a 0 … 1 axis read micro —
  `apply_theme` passed `axes=("bottom",)` where the default is both.
  Deliberately NOT acted on: the fourteen refuted claims, several of them
  plausible (that the QCPMG twin files lose raw units; that the Plotting
  studio is handed normalised intensities; that the pan-limit envelope hides
  the residual strip). One consequence recorded in the manual: a project saved
  under a normalisation mode stores the DISPLAYED match factor, so re-ticking
  the box after reopening in another mode restores the right one.

## Done — the workplan is closed

All 41 kept items are implemented. Phase 5, finished on 2026-09-22:

| Item | What | Commit |
|---|---|---|
| B6 | kernel pre-build on load | ec93d39 |
| A8 | read (C_Q, η, δiso) off a static pattern, three markers | c44475f |
| E6 | pyflakes held at zero by `tests/test_pyflakes_gate.py` | 12eb6a6 |
| F5 | two-way batch table (rows ↔ spectrum cells) | edc2775 |
| F8 | command palette, Ctrl+Shift+P | e682a91 |
| G4 | `czjzek_d`, `czjzek_corr`, `exchange2` + the Czjzek σ read-out fix | 1e9f8b9, 4ca5d2e |
| G2 | tutorials 04–07 + Help ▸ Tutorials | d1a52fd |
| F7 | fit-health strip under the workbench | 2d258d3 |
| G1 | project bundle v2: 2D maps by reference, figures, batch sessions | c2ff861 |
| F3 | sideband auto-detect with a one-click linked manifold | 60f8e49 |
| F6 | FID ⇄ spectrum toggle, real / imaginary / magnitude channels | a938b48 |
| F2 | TopSpin drag-to-phase (Ctrl+P) | 0acf234 |
| G5 | `app.py` split: 387-line facade over eleven mixins + `workers.py` | b2de803 … 130b5e1 |

Beyond the workplan, from Sam's requests and from what the work uncovered
(2026-09-22):

- Component names on the plot (View ▸ Component labels, or hover), Delete /
  Ctrl+D / Ctrl+H on the parameter table, sidebands linked to their parent
  line, the ¹⁹F crystalline fluoride ladder in the literature overlay
  (e5602c5); arrow-key nudging of parameter cells (282cc2a).
- Herzfeld–Berger sideband analysis (Tools), CSA tensor conventions in the
  Conversion tools, `csa_mas` sideband count sized from the static span
  (a969585).
- File ▸ Watch the source file: auto-reload on change, fit kept (5178d7b).
- Decomposition ▸ Label lines from literature ranges (adbb2bf).
- Referencing audit (Tools + `larmor srcheck`): every SR of a session against
  its ¹H adamantane reference by indirect Ξ referencing, TopSpin-ready `sr`
  list, append-only log of old and new values (e193587).
- Defects found and fixed: the χ² profile's 1σ/2σ levels were in raw
  intensity units and collapsed every interval to zero width (bf9939f); a
  positive MASR under a static pulse program loaded as MAS with no warning
  (2fe4800); the Experiment dialog moved the axis the wrong way on an SR
  change (e193587); Czjzek read-outs (P(C_Q) dialog, √⟨P_Q²⟩, Methods text)
  described a distribution half as wide as the one fitted, and the kernel
  clipped 21 % of the mass beyond 5σ (1e9f8b9); worker threads outlived closed
  windows and could abort the interpreter; the frozen exe crashed at start
  (`faulthandler.enable()` with no stderr, 099c09f/510f342).

## The next-ten batch (2026-09-23)

After the workplan closed, Sam asked for a judged shortlist of ten further
items useful to the group, plus a thorough re-check of the QCPMG multi-field
extrapolation. Each item was designed by two independent drafts and a judge
(the designs live in the session scratchpad, `designs/Rank_*.md` and
`QF_qcpmg_multifield_fixes.md`), then implemented on its own branch in a git
worktree and merged into master with a scoped regression. Merge commits, in
order:

| Item | What | Merge |
|---|---|---|
| N1 | Publication bundle for a series: `io/bundle.py` (manifest, README, curves CSV, software versions), `--curves` on the CLI, a "Publication bundle…" button in Batch fit and Sequential fit | f6a5dfd |
| N10 | DFT tensor import from `.magres` (CASTEP / QE-GIPAW): spin-aware seeding, equivalent atoms merged per site, `shiftcal.py` multi-reference σ → δ calibration line with covariance, `larmor shiftcal`, Help ▸ DFT tensors | bc4f4c9 |
| N9 | MAS rate from three sources (acqus MASR, title, NMRFAM booking sidecar `experiment_addenda.xml`) by majority rule; the Experiment dialog lists them with Use / Measure; a per-session confirmation store (`masrate.py`) | 6719054 |
| N2 | Quantitativity chips in the fit-health strip: recycle delay against the sibling T₁ (TopSpin `ct1t2.txt`), flip angle against the (I+½)·θ limit, the tail of every line outside the window; the 90° pulse remembered per nucleus / probe / power (`quantitativity.py`) | aa07f54 |
| QF | QCPMG multi-field, 19 verified fixes: C_Q bounds and an honest σ policy, lever-arm and ν₀ > 0 checks, nucleus check per spectrum, spikelet-comb seeding, MAS-aware CG windows with a convergence σ, one CG estimator, processing mode carried into the report, width split over every field, χ² for ≥ 3 fields, duplicate names fitted separately, per-field provenance, P_Q with its η-free error, a whole-chain test on simulated ³⁵Cl patterns, tutorial 7 regenerated | fedd782 |
| N5 | Parameter status markers † fixed · ‡ at a bound · § linked, derived from the recipe (`paramstatus.py`), in every CSV / LaTeX / Markdown table, the Methods sentence and the batch report; the at-bound note no longer accumulates across refits | e78a009 |
| N8 | Comparability check before a batch or sequential fit: acqus / procs / auditp of every member against the series majority (`comparability.py`), an amber verdict line, a Details table, "Reprocess all from fid" with one common chain, `larmor compare` | 6244e8f |
| N7 | Session inventory (Tools, `larmor inventory`): every EXPNO of a month folder as a sample × nucleus grid with the production pick per block; sample identity from the folder name (`scan.sample_name`), the title kept in the provenance | e3c400d |
| N6 | Series table (`series_table.py`, one dialog for Batch fit and Sequential fit): names, replicate groups and order of the series, composition columns joined from a CSV, the Series plot on a numeric x axis with replicate averaging and an OLS line, Species bar, DUST export, `_series.csv` sidecar | de906b4 |
| N4 | Provenance block: the acquisition record and the SHA-256 of the data file in every recipe, a software stamp on every fit, the full Experimental paragraph behind Copy methods, Table S1 across a series (Tools ▸ Experimental section…, `larmor acqtable`), reload checks of hash and SR | 5644efb |
| N3 | Site families: `SiteModel.family`, summed populations and named ratios (N₄, ⟨CN⟩ Al, ⟨n⟩) with errors on three stated bases — amplitude covariance, per-trial Monte-Carlo integrals, or a flagged independent fallback (`families.py`); family column and menu in the table, rows in the Report and the batch tables | 4c4e242 |

Measured on the group's data while implementing (each number is pinned in a
test or a tutorial block):

- The Tutorial-4 ¹¹B series was processed with LB = 0 / 100 / 0 / 100 / 100 Hz
  and acquired with D1 = 12.5 … 36 s; six of the ten 2026-01 glasses sit at
  2.5–3.3 T₁ for their slowest boron site (N2, N8).
- The 2026-05 ³¹P series is cut at TDeff 384 of TD 9590 with ABSG 0 or 5, and
  the three P5-Bi8-12 repeats were acquired at D1/NS = 300 s/512, 300/140 and
  90 s/64 (N7, N8).
- The acqus MASR of almost every NMRFAM EXPNO is a stale 4200 Hz; title and
  booking agree on 405 EXPNOs and disagree three ways on the 2026-05 series
  (title 20 kHz, booking 22 kHz) (N9).
- Derived parameter status reproduces the fit's own at-bound notes on all 32
  Final2 recipes (246 of 504 free parameters at a bound, 190 fixed) (N5).
- The CaF₂ / NaF GIPAW calibration line is δ = −0.687 σ + 51–53 ppm; the
  2026-03 fluoride standards carry a +1.67 ppm session offset (N10).
- Base0Ca ¹¹B, a three-line demonstration model: N₄ = 0.5347 ± 0.0010
  (covariance) ± 0.0010 (Monte-Carlo, 200 trials) ± 0.0020 (independent) —
  the two overlapping BO₃ lines are anticorrelated, so the family sum is known
  ~3.7× better than quadrature says and the ratio ~2× worse (N3).
- Whole-chain QCPMG test: static ³⁵Cl C_Q 3.0 / 4.5 / 5.5 MHz recovered to
  0.1 ppm / 0.02 MHz; under MAS the first-minima window drifts 21 ppm at
  78 MHz for C_Q 5.5 and the convergence flag fires; the whole-manifold route
  with all sidebands recovers within 1.5 ppm / 0.1 MHz (QF).

Open observations from the batch (not fixed; for Sam):

- LARMOR's kernel-based `czjzek` model gives a whole-manifold centre of
  gravity biased −21 ppm at 78 MHz for σ = 1.0 MHz (−4 ppm at 108 MHz, under
  0.1 ppm at σ = 0.6) against the first-moment theorem, with zero intensity
  at the axis edges — consistent with the 150 kHz kernel window aliasing the
  large-C_Q rows. Direct mrsimulator sites are exact; the QCPMG tests use
  those. Worth a separate look at `engine.build_kernel`'s window.
- `redor.analyze` uses the spin-½ second moment for ¹⁹F{²⁷Al}; whether the
  S(S+1) factor should enter is Sam's call.
- A site frozen outside the fit window derives ‡ at min 0 while the fit's own
  at-bound list omits it (development-notes §8); one line in `fit()` setting
  the frozen site's `Param.vary = False` would align them.
- The shipped `examples/pCABS2-4_11B.recipe.json` fit (2 Amorphous + 1
  Gauss/Lorentz, 19 free parameters) converges without a covariance matrix, so
  its family block reports the independent basis; Monte-Carlo gives N₄ =
  0.5704 ± 0.0003.

## Knowledge that must not be re-learned

- **Czjzek σ convention.** mrsimulator's `sigma` is HALF of Czjzek's σ_Cz
  (`sigma_ = 2*sigma` in its density): dmfit sCZ_CQ = σ_Cz = 2σ. Fits and
  recipes were always right; only derived read-outs were wrong. Kernels are
  requested to 10σ (`CZJZEK_KERNEL_HEADROOM`), σ max is 40 MHz.
- **The kernel axis is not ppm-by-γB₀**: `engine.kernel_axis_ppm` against
  `Isotope.B0_to_ref_freq` (0.08 % off γB₀ for ²⁷Al).
- **Static detection**: `MASR=0` means static, corroborated by title or
  pulse program; a POSITIVE MASR under a wcpmg/wurst pulse program is flagged
  (the MagLab ⁸¹Br set records 14 kHz on a static probe).
- **MAS rate by majority**: acqus MASR, the title's `MASR … kHz` and the
  booking sidecar are three candidates; two agreeing settle it, a three-way
  disagreement takes the highest and flags; a confirmation is remembered per
  session folder / rotor / nucleus. The acqus value alone is a leftover on the
  NMRFAM spectrometer.
- **Sample identity** comes from the sample folder (`scan.sample_name`:
  date / rotor / operator tokens stripped), never from the title's pulse
  note; the title's first line and the raw folder live in the provenance.
- **Herzfeld–Berger intensities**: average over the tensor azimuth γ and β;
  the rotation about the rotor axis is a pure time shift. The second moment
  Σ I_N (Nν_r)² = (ζν₀)²(1+η²/3)/5 is the check.
- **The χ² profile levels** are χ²_min(1 + 1.00/dof) and χ²_min(1 + 3.84/dof):
  the fit is unweighted, so Δχ² must be scaled by the residual variance.
- **Indirect referencing**: SF_X = SF_¹H · Ξ_X/Ξ_¹H with mrsimulator's ratios
  reproduces TopSpin `xiref` to 0.14 Hz; adamantane ¹H is at 1.82 ppm in
  this group; a session is one month folder; the audit log is append-only.
- **Two-field extrapolation**: `qcpmg_fields.infinite_field_diso` carries the
  σ policy (equal weights when a σ is missing, the 0.1 ppm floor, ± scaled by
  √χ²_red never down) inline; `weighted_line` (the DFT calibration's line)
  shares its algebra with the pairwise denominator Σ w_i w_j (x_i − x_j)²;
  a test pins the two agree to rounding.
- **Parameter status is derived, never stored**: fixed = `vary False`, linked
  = `expr`, at bound = the fit's own 1e-3 relative rule against the effective
  bound (recipe min/max, else the registry's). No recipe field, so it is
  retroactive on every saved recipe.
- **Family errors need integrals, not amplitudes**: amplitude is peak height
  for every model but `gl_norm`, so the Monte-Carlo basis re-integrates each
  site per trial; the covariance basis uses k_i = integral/amplitude at the
  best fit over lmfit's `uvars`. The basis is always printed.
- **Threads**: every QThread the window owns is joined in `closeEvent`;
  `tests/conftest.py` sets `LARMOR_NO_KERNEL_WARM=1` so fixture windows never
  start a kernel pre-build (a QThread destroyed while running aborts Python).
- **Shortcuts**: the palette owns Ctrl+Shift+P; Save project is Ctrl+Alt+S;
  Ctrl+P drag-to-phase; Ctrl+T FID toggle; Ctrl+I channel; Ctrl+Shift+D
  sideband detect. `tests/test_app_split.py` and the palette test assert
  uniqueness; the menu tree golden there is regenerated whenever a Tools row
  is added (N7 and N4 each added one).
- **Adding a lineshape model** fails loudly if under-declared:
  `tests/test_model_tables.py` (five partitions, the fifth is the DFT seed
  keys) plus the `lineshapes.md` backtick gate; parameter names must be unique
  across models (the `function` model owns `d`, hence `czjzek_d`).
- **Writing files with backslashes from this tooling**: never through a
  shell heredoc (it collapsed doubled backslashes into control bytes once, in
  a display equation, and a heredoc with LaTeX-heavy Python failed to parse
  once more this session); write a Python script file or use an editor tool.
- **Worktree testing**: the conda env's editable install resolves `larmor` to
  the main checkout; a worktree is tested with `PYTHONPATH=<worktree>`.
- **A test must never restyle the QApplication.** `QApplication.setStyleSheet` re-polishes every live widget (0.05 ms
  each); the suite never destroys its windows, so in its last fifth one
  such call takes minutes. Six parametrised theme cases blew the
  10-minute `faulthandler_timeout` three times and killed a run at 83 %.
  A test that needs only the active theme calls `theme.apply(None, name)`.
- **A green pytest bar is not enough**: read the real-data banner.

## 2026-09-30 session — Sam's third field-use batch → v0.15.0

Eleven items from field use, all delivered; the two "broken" ones were real
defects with the same signature Sam described.

**Root causes found (each with a regression test):**

- *"The fit button is broken, it's visibly not fitting, where the autofit
  button is working."* — the completion threshold (b70e9a7, shipped in
  0.14.0, hence also the student's report): `_emit_progress` compared the
  residual of successive `iter_cb` calls, but lmfit calls `iter_cb` on every
  residual EVALUATION, and the Jacobian probes differ by 1e-8…1e-4, so at the
  default 0.1 % every Fit-button fit aborted after 5 evaluations with the
  positions untouched. Auto fit passes no threshold, hence worked. Fixed:
  the callback aborts only on Stop / Cancel; the threshold is the solver's
  `ftol` alone (a matched `xtol` stopped a synthetic line 1.1 ppm short once
  the 1e6-scale amplitude settled). Probe on the shipped 27Al example at
  0.1 %: nfev 5 / RMSD 0.082 before, nfev 106 / RMSD 0.05363 after (0.05362
  at full precision). `tests/test_fit_threshold.py`.
- *"autophasing isn't working"* — two halves. (1) A Bruker 1r was loaded
  real-only, so a phase rotation was `y·cos φ`: Autophase changed nothing
  and a typed p0 only scaled the line. `bruker.read_imag` now reads
  TopSpin's 1i (a 90° p0 turns the real channel into exactly −1i on 3616);
  CSV / dmfit sources get `hilbert` inserted before their first phase step
  (`_complex_for_phase`). (2) An opaque `autophase` op was invisible to the
  processing panel, whose next live tick re-emitted the widgets' phase on
  the raw spectrum — the line snapped back on the first nudge. Now
  `resolve_autophase` writes the angles into the Phase controls as a plain
  phase step (like `apk` filling PHC0 / PHC1), `apply_processing` syncs the
  panel to the recorded chain after every apply, and the Process menu
  entries append to the chain instead of replacing it
  (`append_processing_step`). `tests/test_autophase_channel.py`,
  `tests/test_autophase_resolution.py`.

**UI items:** checked push buttons take the accent fill (`QPushButton:checked`
— Pick anchors / Pick 2 points / Drag to phase now look on); the
spinning-sideband offer on load is OFF by default under a NEW settings key
(`SSB_OFFER_KEY`; the old key had been persisted True on every install);
the NMR table fills every cell with Youngman 2018's glass-NMR feasibility
(five categories read cell by cell from Fig. 1, `nuclei.FEASIBILITY`, legend +
tooltips); File ▸ "Auto-reload when the file changes" says what the watch is
for (process in TopSpin, fit here); "Autophase (ACME)" is "Autophase" (the
op has always run the p0-sweep search; ACME is the recorded op's `method`).

**Features (four worktrees, merged one at a time with scoped regressions):**

- **WS** — Workspaces dock: Ctrl/Shift multi-select and a right-click menu
  (Switch / Open, Rename, Save, Close, Close others, Send to Plotting studio
  with the data + fit + components of every selected document, Fit parameter
  table, Overlay on the active spectrum). `larmor/fittable.py` (PARAM_COLUMNS
  moved to the core; FitEntry from a recipe or a saved fit file; wide / long
  tables; CSV / TSV) + `desktop/fittable_dialog.py`; reached from the
  Workspaces menu, the Explorer's fit rows (multi-select) and Fit ▸ Fit
  parameter table.
- **INV** — Session inventory remedies: right-click a row / Fix… → apply the
  audited SR when the EXPNO opens in LARMOR (the per-EXPNO override store
  `%LOCALAPPDATA%/LARMOR/sr_overrides.json`, applied by `loader` while the
  file's SR still matches, logged in the referencing log, undoable), all
  flagged picks at once, copy the TopSpin `sr` command, open the ¹H
  reference, rename (alias / on disk), process an unprocessed fid in
  LARMOR, pick anyway, dataset info. `inventory.fixes_for` decides; the
  dialog renders.
- **PYB** — Process ▸ Baseline ▸ pybaselines: 13 methods (arPLS, asLS,
  airPLS, iarPLS, drPLS, asPLS; ModPoly, IModPoly, penalized poly; SNIP,
  rolling ball, mor, mormol) behind a live-preview dialog, recorded as the
  replayable `pybaseline` op, λ scaled with the spectrum length; dependency
  `pybaselines>=1.1` in pyproject / environment.yml / the PyInstaller spec /
  `packaging/.buildenv`. Also corrected a manual that called `op_baseline`
  arPLS.
- **PP** — Tools ▸ Pulse program (and the Explorer's EXPNO menu): the
  `pulseprogram` file parsed (`larmor/pulseprog.py`: definitions, `(center …)`
  blocks, loops, mc macros, phase programs, legend), values resolved from
  acqus with a safe expression evaluator, a TopSpin-like not-to-scale timing
  diagram (channels with nuclei, pulses / shaped pulses / delays / decoupling
  bars / acquisition glyph / loop brackets, PNG + SVG export), the text and a
  parameter table; help page "Reading a pulse program" in ? ▸ User manuals.

**Docs:** `development-notes.md` §1 module map and §5 rows for the fit
threshold, the 1i channel and the settings key, §8 items 18–22; the manuals'
menu map and the Autophase wording; `validation.md`.

**Full suite:** see the State-of-the-repository line for v0.15.0.

**Follow-ups noted, not done:** the batch-fit dialog's baseline choices have
no pybaselines entry; asPLS is slow and poor on 64k-point spectra; the SR
override applies only through `loader._load_any` (not Open FID / 2D); the
panel re-emits a chain recorded as `[baseline, phase]` as `[phase, baseline]`;
`Baseline_Corrector.m` (the Yon 2020 MATLAB source) sits untracked in the
repo root and is not part of the package.

## 2026-09-30, later — "check all the error calculation when there is only one line" → v0.15.1

Sam could not calculate errors on a single-line Czjzek fit of ²³Na. Reproduced
on a real ²³Na MAS spectrum (a CEMHTI NaAlSiO glass, zg, 132 MHz) with one
Czjzek site seeded the way the app seeds it: every tool *returned*, and four
things were wrong.

- **The χ² profile dialog looked dead.** It used the process pool
  unconditionally; the pool's start-up (a spawned interpreter per worker
  importing mrsimulator, ~15 s here, longer in the frozen exe) blocked the
  modal window with a fixed "scanning…" label, while the whole scan was 2 s
  of work. Now `parallel="auto"` (the first point is timed in-process; the
  pool only when the rest would outlast `autofit.POOL_BREAKEVEN_S`), a
  `heartbeat` that pumps the event loop while a pool starts, a progress bar,
  Stop, and span / points controls in the dialog; Monte-Carlo the same. A
  pool that breaks falls back to sequential work with a warning
  (`parallel.parallel_map`).
- **The profile disagreed with the table in silence.** It refits everything
  else at each scan point and found a χ² 35–59 % below the fit's own: the
  fit on screen had stopped on a slope. The profile now reports
  `fit_at_minimum` and says "the fit had not converged — Fit again, then
  rescan"; Monte-Carlo names the parameters its refit moved by more than σ.
- **The fit really had stopped on a slope** — the root cause. scipy's
  trust region and `xtol` step test work in raw parameter units, where the
  amplitude (1e5–1e7 counts) dominates: a step moving a position by a few
  hundredths of a ppm already "converged". `fit.X_SCALE = "jac"` scales every
  parameter by its Jacobian column. Measured: the ²³Na line goes from RMSD
  0.116 at a nonsensical +52.7 ppm to 0.073 at −9.8 ppm (σ(C_Q) 0.48 MHz,
  dCS 38 ppm); the shipped ²⁷Al example from its 0.0468 to 0.0458 — the
  minimum a 12-start Auto fit finds too (4 of 5 restarts), and whose
  positions 62.8 / 28.6 / −0.2 ppm agree with the 2D MQMAS fit (62.2 / 29.7 /
  −1.1) better than the old 63.8 / 30.9 / −2.0; the ¹¹B example is unchanged
  (0.00372); the 127 fit-related tests pass either way. The shipped ²⁷Al
  recipe was re-saved at the deeper minimum (`larmor fit … --window 150
  -80`), the README figure re-rendered, the two tests that pin its numbers
  and Tutorials 1, 5 and 6 re-measured.
- **The Report gave one line a population of 100 ± 18 %.** The fraction
  error was `100·σ_i/T`, the integral's error dressed up as the error of a
  share that cannot vary; `quantify.fraction_err_pct` now differentiates
  `I_i/ΣI` through the normalisation (a single line: exactly 0).

Tests: `tests/test_error_tools_single_line.py` (the single-line population,
the propagated fraction error against a numerical derivative, the
not-at-minimum notes of both tools, the interval-narrower-than-step note,
the "auto" decision, the pool fallback, the heartbeat through a real pool,
the shipped ²⁷Al deeper minimum, the dialog with progress and Stop).

## 2026-10-02 — three items and two frozen windows → v0.15.2

- **"Why I cannot open more than 4-5 spectra … they are not staying open in
  the workspace"** — `_register_ws` reused the active 1D workspace whenever
  it carried no fit, so browsing a sample folder replaced the previous
  spectrum each time and only fitted ones accumulated. Only a workspace that
  holds no spectrum is reused now (Reset to original / Make active still
  replace in place). Six double-clicks are six workspaces
  (`tests/test_workspace_policy.py`).
- **"need in the File Menu a open recipe or open fit to use on currently
  open data"** — File ▸ Open a fit on this spectrum… (the browse half of
  Edit ▸ Apply recipe), and the Explorer's fit rows offer "Apply to the open
  spectrum".
- **"Rethink the sequential fit … really rethink and rework it completely"**
  (with two screenshots of the old dialog "(Not Responding)") — replaced by
  a **series mode of the main workbench** (worktree wip/SEQ, merged as
  b87cd4f): Series ▸ Sequential fit… (Explorer selection, inventory picks,
  or the open spectra) opens every member as its own workspace and shows a
  series bar above the spectrum — member buttons with status dots (grey
  unfitted, green fitted, amber edited since fit, red failed), ◀ ▶, Carry ▾
  (seed the next spectrum when I move; which parameter groups; copy or
  replace the model everywhere), Fit → next (the window's own Fit, in its
  thread, then move and seed), Auto sweep ▾ (passes / start / smooth; the
  `seqfit.run_sequential` sweep in a `SeqWorker` thread with Stop), Table…,
  Plot…, Acquisition…, Save ▾ (all fits, publication bundle), the
  comparability chip, End series. Because a member is an ordinary workspace,
  every tool applies per spectrum: add / remove lines, the full lines table,
  right-click menus, processing, paddles, constraints, Auto fit, the
  fit-health strip, undo. Carry rule: an empty member gets a copy of the
  model with amplitudes scaled by the intensity ratio; a member with its own
  lines keeps them and receives only the carried values (matching index and
  name, clipped to its bounds). The series survives a project save / open.
  Found on the way: the old dialog's "Auto ⇄ forward–backward fit" button
  had been broken since its worker stored the start choice as `self.start`,
  shadowing `QThread.start()`.
- The Session inventory also showed "(Not Responding)": its scan and the
  referencing audit ticked progress per sample folder; they tick per EXPNO
  now.

Tests: `tests/test_seriesmode.py` (10), `tests/test_series_mode_ui.py`
(14, through the FitWorker and the SeqWorker), `tests/test_workspace_policy.py`
(4); `seqfit_dialog.py` and its tests removed; golden menu 234 rows.

## 2026-10-02 (later) — the series mode's carry rule, label pairing, the bar's black box → v0.15.3

- **"by using the sequential fit, it is never well fitting the previous one
  … keep the fit of the sample 1 good when changing to sample 2, so when I'm
  going back to sample 1 it is still good"** — the first series mode
  re-seeded the spectrum you landed on from the one you left, fitted or not,
  so walking back ruined finished fits. Now **moving never changes a spectrum
  that has lines**: only an empty member takes a copy of the model you leave;
  the value carry is explicit (Carry ▾ *Seed this spectrum from the previous /
  next one*, or the member's right-click menu) and lands only on lines that
  share a **label**; **Keep this fit (🔒)** protects a member from every carry
  and from the sweep, which still starts its neighbours from it
  (`run_sequential(fixed=…)`). The member menu also renames and removes a
  spectrum from the series.
- **"multiple 19F spectra … evolution of the components along a series with
  different components"** — lines are paired across spectra by label
  (`larmor/components.py`: `pair_sites`, `component_map`), by index only
  where a label cannot identify its line: the carry, the sweep's warm start
  and smoothing, and the Series plot (components read `(n/N)` in its list;
  a gap in the shape parameters and 0 % population where a spectrum's model
  has no such line; the Species bar follows).
- **The screenshot's black box over the plot and the nameless current
  member** — one cause: the strip body's bare `background: transparent`
  style sheet cascaded onto every member button (the checked one lost its
  accent background: white name on white) and onto their tool tips (a
  transparent top-level window, painted black; only the colour-emoji ⚠
  survived). Every rule in the bar now carries a selector; the strip's
  scrollbar is hidden (wheel / ◀ ▶ scroll it) so it never overlays the
  buttons.
- Found on the way: `SeqWorker` turned an all-unticked Carry list (`()`)
  into "carry everything".

Tests: `tests/test_components.py` (3), `tests/test_seqfit.py` (+3),
`tests/test_seriesmode.py` (+2), `tests/test_series_mode_ui.py` (the carry
tests rewritten, +4: a fitted member untouched by moves, the explicit seed by
label, keep-this-fit through seeding / copy / sweep / project, member-menu
removal, the bar's rendering), `tests/test_series_export.py` (+1 series whose
members differ).

## Remaining

1. **E7** — `dist/LARMOR-0.14.0-setup.exe` and the zip are built and tested
   on the build machine; run the setup on a machine with no development
   setup (a student's laptop: install, open a Bruker folder, fit, save a
   recipe, uninstall), then publish a GitHub release carrying both
   (`gh release create v0.14.0 dist/LARMOR-0.14.0-setup.exe
   dist/LARMOR-0.14.0-win64.zip`).
2. The open observations of the next-ten batch above (kernel-window CG bias
   of the Czjzek model, REDOR S(S+1), the frozen-site marker).
3. Follow-ups noted earlier, none blocking:
   - `_update_sn` re-runs sideband detection and rebuilds the literature
     overlay on every live processing apply (about 15 ms at 64k points);
     trim it if drag-to-phase feels sluggish on wideline data.
   - `czjzek_corr` and `exchange2` export to dmfit as the commented
     Gauss/Lorentz envelope; no 2D (MQMAS) path for the three new models.
   - Manual baseline is not undoable; `save_recipe` shows its status twice.
   - Deferred by their designs: CLI `--series` for the series table, a
     QThread port of the synchronous reprocess loop, the Rician correction of
     magnitude-mode QCPMG centres of gravity, batch grids marking spectra
     whose MAS rate is unconfirmed.
4. Deferred by Sam (revisit only if asked): A9 C4 C6 E2 E3 F1 F4.
   Dropped: A2, C5.

## Conventions that keep this project safe

- Scope test runs to what changed; **full suite + real-data banner before any
  push**, never push without Sam's word.
- Version bumps touch `larmor/__init__.py` **and** `pyproject.toml`.
- Commits are `vX.Y.Z: summary` or `<item>: summary` with a post-mortem body
  (what was wrong, how found, measured effect) — `git log` is the project's
  bug documentation.
- Public docs: neutral, professional voice, no first person.
- Parallel items go on `wip/<ID>` branches in `../LARMOR_wt_<ID>` worktrees
  and are merged one at a time with a scoped regression; shared files are
  edited additively so the merges stay textual.
- `docs/development-notes.md` is current through this session (§1 module map,
  §3 conventions table and §8 defects index carry the batch's additions).
