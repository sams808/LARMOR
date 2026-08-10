# Road to LARMOR 1.0

v0.7.1 is feature-rich — arguably feature-*complete* for its headline use case
(Czjzek/quadrupolar fitting, MQMAS, batch workflows, publication figures). **1.0
is a different claim**: *"you can rely on this for published work and recommend
it to a colleague."* This file is a living status check against that bar, not a
feature backlog — re-audited 2026-08-07 against the actual repo state (the
previous version of this file, from 2026-08-05, predated ~15 phases of work and
had drifted from reality in both directions: some things it asked for now exist,
some things it assumed were fine have not been touched since).

**Status key:** ✅ done · 🟡 partial/exists but thin · ❌ not started

## What "1.0" must guarantee

1. **Trustworthy numbers** — validated against dmfit / literature on real data,
   with honest uncertainties everywhere. → 🟡 (see Pillar A)
2. **Never loses work, never crashes** — clear errors, recoverable sessions,
   forward-compatible file formats. → 🟡 (see Pillar B)
3. **Maintainable** — no single 4k-line file; CI keeps it green. → ❌ (see Pillar C)
4. **Installable & documented** — a colleague can install it and learn it
   without you. → 🟡, and the weakest link right now (see Pillar E)

## Pillars — current status

### A. Trust & validation

- ✅ **Validation report, now with a real multi-sample check** (2026-08-09):
  `docs/LARMOR_VALIDATION_REPORT.md` §5.7 imports **20 real dmfit fits** from
  an actual peer-reviewed paper (Soudani *et al.* 2024, $^{11}$B pressure
  series) — zero warnings, exact parameter agreement against the paper's own
  published Table 2 — closing the "one real-data comparison isn't many
  published fits" gap flagged below. Found and fixed a real crash
  (`_lorentz_convolve` kernel-length bug) in the process. Still predates
  batch fit/sequential fit/exclusion/population%/CPU-parallel error analysis
  validation-wise (those remain un-cross-checked against any real published
  dataset).
- ✅ **dmfit round-trip against many real files**: as of the §5.7 pass, 20
  real `.fxml` files (not just `CaAlGlass.fxmla`) parse and reproduce their
  source paper's values — the specific gap this line used to describe.
- ✅ **Uncertainties**: covariance / Monte-Carlo / χ² profile all present,
  now CPU-parallelized (phase 2aj), propagated into `quantify()` populations
  and batch CSV exports (phase 2ah's `population_pct` rows). This pillar's
  actual mechanics are in good shape.
- ❌ **`humanize_error`**: named in the original roadmap as "continue the
  work started in v0.2" — it does not exist anywhere in the codebase today
  (`sanitize_constraints` does, in `constraints_util.py`). Either it was
  renamed away at some point or the claim was already stale when written;
  either way, raw exception text is what surfaces in several dialogs today
  (spot-check: `try/except Exception as exc: self.status.setText(f"...{exc}")`
  is a common pattern across `desktop/*.py` — functional, not colleague-friendly).

### B. Robustness & data safety

- 🟡 **Project bundles (`.larproj`) — correction, 2026-08-09: this item was
  WRONG in the previous version of this doc.** A real, working, versioned
  bundle format already existed (`app.py::save_project`/`open_project`,
  `"larmor_project_version": 1`, `tests/test_project.py`) — every open 1D
  workspace's recipe + processed spectrum round-trips through one file. The
  earlier audit pass concluded "doesn't exist" from seeing `.larproj.json`
  only as a `QFileDialog` filter string, without reading the functions
  behind it — a real methodology gap in that pass, corrected here.
  **Genuine remaining gaps**, found and partly fixed this session:
  - ✅ **Fixed**: overlays (comparison spectra) were tracked in the live
    snapshot but never written into the saved file — silently dropped on
    every reopen. Also fixed, in the same code path: `add_overlay_dialog`
    unpacked `load_any()`'s return tuple in the wrong order, so adding a
    `.recipe.json`/`.fxmla`/`.csv` file as an overlay always failed (only a
    raw Bruker path worked, "by accident", via the fallback branch).
  - ✅ **Fixed**: `larmor.batch.load_entries` (the batch-report tool) claimed
    in its own docstring/help text to accept `.larproj` files but had no
    code path for their actual multi-workspace shape — every project file
    silently failed to load with "could not load" and was skipped. Now
    expands one entry per workspace.
  - ❌ **Still not done**: 2D workspaces are explicitly excluded from
    `save_project` (`Data2D`/view state isn't JSON-serializable as-is —
    real, non-trivial engineering, not attempted this pass). Figures and
    live batch-fit-dialog sessions are still not part of a project bundle.
  - **Known, accepted design tradeoff** (not a bug): a saved project embeds
    the processed `exp_ppm`/`exp_amp` arrays directly rather than only a
    `source_path` reference + replayed processing — unlike `Recipe`'s own
    "never embed" convention. Keeps a project openable even if the source
    file moves, at the cost of file size; not changed this pass.
- 🟡 **File-format versioning**: `Recipe.from_dict` is forward-compatible in
  practice (new optional fields default sanely, old files still load — this
  was exercised repeatedly this session, e.g. `recipe_from_csv_rows`'s
  registry-default fallback) but there's no explicit schema version field or
  migration path if a FUTURE breaking change is ever needed.
- 🟡 **Dialog smoke tests**: `test_desktop.py` is large (29 top-level test
  functions, many multi-assertion) and a `win` fixture builds the full
  `MainWindow` headlessly for most of them — real, substantial coverage —
  but there's no single dedicated "construct every dialog class with
  minimal/empty input and assert no crash" sweep; coverage is a byproduct of
  feature tests, not a systematic guarantee.
- 🟡 **Fuzzing malformed files**: `test_processing_v2.py` has some malformed-
  input cases; the universal reader (`larmor/io/bruker.py` etc.) has organic
  robustness from real-data testing, but there's no dedicated "throw garbage
  bytes at every loader, assert a clean error not a crash" suite.
- ✅ **Undo**: real coverage — `_capture_state`/`_restore_state`,
  `snapshot(with_axis=True)`, wired into calibrate/SR/constraint edits per
  phase 2ad. Not audited exhaustively for every processing op or the 2D
  path, but substantially more than "nothing."
- ❌ **Error-quality work**: see `humanize_error` above.

### C. Maintainability

- ❌ **Split `app.py`**: **not done, and it has grown, not shrunk** —
  4,426 lines now (vs. "4.2k" when this was first flagged as the #1
  recommendation). Every phase since has added more methods to the same
  file. This is the single largest deferred item in the whole roadmap.
- ❌ **CI**: no `.github/workflows/` directory exists at all. Every "full
  suite green" claim in this project's history has been a manually-run
  local `pytest` — real, but never enforced automatically on push, and
  never run on a clean machine/environment other than this dev conda env.

### D. Packaging & distribution

- 🟡 **Windows installer**: `packaging/larmor.spec` + `launcher.py` exist
  and are actively maintained (this session added `multiprocessing.
  freeze_support()` to it) — but there's no evidence in this repo of a
  recently *built and run* `.exe` being verified end-to-end; the spec being
  correct and an actual working installer being confirmed are different
  claims.
- ❌ **macOS installer**: not attempted — no spec, no notes, and this
  project has only ever been developed/tested on Windows.
- ❌ **Release pipeline**: no automated "tag → build installer" anything;
  releases so far are `git tag` + a manually-run local build when needed,
  and even `gh release create` is blocked in this environment (noted in
  memory) — GitHub Release *pages* for v0.3.0/v0.3.1 were never actually
  created, only the tags.

### E. Documentation & onboarding — **the weakest pillar right now**

- ❌ **README.md is badly stale** and actively misrepresents the app's
  current state:
  - Describes the project as being at "Phase 0/1/1b/1c/1d/1e" — i.e., an
    early proof of concept — with **zero mention** of the native desktop app
    being primary (has been since Phase 1f, ~40 phases ago), 2D MQMAS
    interactive fitting, batch/sequential fitting, Monte Carlo/χ² error
    analysis, satrec/T1/REDOR/QCPMG, DFT import, SIMPSON, the Plotting
    Studio, or hidden themes — i.e. most of what the app actually does today.
  - Contains a **broken file reference**: `See ROADMAP.md for what's next`
    — that file does not exist (only `docs/ROADMAP_v1.md`, this one, does).
  - Presents the FastAPI web app (`larmor/app.py`) as "a browser variant...
    secondary" with no caveat — that file **hasn't been touched since
    2026-07-15** (Phase 1e) and has **zero test coverage**; it is almost
    certainly broken against the current recipe/model schema and should
    either be fixed, tested, or explicitly marked unmaintained, not
    presented as a working alternative.
  - Says a standalone installer "is on the roadmap (Phase 5)" — the
    PyInstaller spec has existed since Phase 1k; this line is simply wrong.
- 🟡 **In-app manuals**: `larmor/help/*.md` is substantial and has been kept
  reasonably current (multi-dataset.md especially, updated again this
  session) — genuinely one of the better-maintained parts of the docs.
- 🟡 **Tutorials**: 3 exist (`docs/tutorials/01-03`), all from Phase 1d/1e —
  none cover anything from the desktop-native era onward (no batch fit
  tutorial, no MQMAS tutorial, no error-analysis tutorial).
- ❌ **Bundled example datasets**: `examples/` has exactly **one** dataset
  (`CaAlGlass.recipe.json` + a png) — no 2D/MQMAS example, no batch-series
  example, no relaxation/QCPMG example, despite the app supporting all of
  these as headline features.
- 🟡 **Headless API docs**: CLI subcommands exist and are reasonably
  numerous (`info/import/fit/app/desktop/satrec/redor/magres/multifit/
  seqfit/batchfit`) but aren't documented anywhere outside `--help` text and
  scattered README mentions.
- ⏳ **Two other planning docs exist but were never reconciled with this
  one**: `docs/GPU_ACCELERATION_PLAN.md` (written 2026-08-06, thorough,
  explicitly gated on "profile again after CPU fixes land, decide if still
  worth it") and `docs/ideas_beyond_request.md`/`docs/ideas-from-field-use.md`
  (mostly implemented per the memory log, but the field-use one has **no
  status markers at all** — every item still reads as open even though most
  are done). This doc, that plan, and both ideas lists should really be one
  coherent picture; right now a reader has to already know the project's
  history to know which parts of which document are current.

### F. Feature completeness

Everything the *previous* version of this roadmap listed as "still open" —
**confirmed still open, zero code found for any of them**:
- ❌ **#14** labelable peak-pick assignment table
- ❌ **#12** TopSpin drag-to-phase gesture
- ❌ **#19** dmfit amplitude-calibration UI
- ❌ **Bootstrapped series errors**
- ❌ **Two-way batch table** (row ↔ plot selection)
- ❌ **Reference/literature overlays** — when built, seed the data from
  Edén 2023's compiled tables (δiso ranges per coordination for ²⁷Al/¹¹B/
  ¹⁷O/²³Na/²⁵Mg and P_Q ranges, his Tables 3–5 — PDF at `Desktop\bib\`),
  drawn as labeled shaded δiso spans for the current nucleus
- ❌ **Watch-folder / auto-fit**
- 🟡 **Figure annotation layer**: a spec-level primitive exists
  (`{"x","y","text"}` → `ax.text(...)`, `render_1d` only) but there's no
  Studio UI to place/drag/edit one — usable only by hand-editing a JSON spec.
- 🟡 **Undo everywhere**: see Pillar B — real but not exhaustively audited.
- ❌ **Keyboard-driven fitting**: no evidence found.

None of these are small effort individually, but none of them are blocking
either — they're genuine "nice to have," not "can't trust the numbers."

## If only three things before 1.0 (revised 2026-08-09)

1. ✅ **A real validation report** — substantially strengthened 2026-08-09
   (20-fit real published-paper cross-check, §5.7). Remaining gap: the
   features added since 2026-08-01 (batch fit, exclusion, population%,
   sequential fit) still have no real-dataset cross-check of their own.
2. 🟡 **Project bundles** — turned out to already exist (`save_project`/
   `open_project`, corrected above); its two concrete bugs (overlays
   dropped, `add_overlay_dialog` mis-reading its loader's return tuple)
   are fixed 2026-08-09. Still missing 2D workspaces/figures/batch state —
   real remaining work, just narrower than previously described.
3. ❌ **Split `app.py` + CI** — neither started; `app.py` has grown since
   this was first flagged, not shrunk, and there is no automated test gate
   on push at all. **This is now the single largest genuinely-untouched
   item on this whole roadmap.**

Already done as of this pass: **`README.md` rewritten** (2026-08-09) — was
describing a different, much earlier program, linked to a file that didn't
exist, and recommended a dormant, untested component as if current.

**Process note for future roadmap passes**: the project-bundle correction
above happened because a *test file* (`tests/test_project.py`) was noticed
in passing and read fully, not because the earlier grep-based audit found
it — that audit searched for the string `"larproj"` and stopped once it
found file-dialog filter strings, without checking whether the FUNCTIONS
behind those dialogs (`save_project`/`open_project`) actually did anything.
Grepping for a feature's *name* finds where it's mentioned; it does not
confirm whether it's implemented — checking the actual function bodies
(or the test suite, which will not pass for a feature that doesn't exist)
is the only reliable way to answer "does X exist," and is worth doing
before writing "not started" into a doc like this one again.

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
*Last audited: 2026-08-07 (v0.7.1), against the actual repo state, not the
previous version of this document.*
