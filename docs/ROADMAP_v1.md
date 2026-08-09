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

- ❌ **Project bundles (`.larproj`)**: this is the single most overstated
  item in the previous version of this doc. `.larproj.json` exists **only
  as a file-extension string** in a few `QFileDialog` filters
  (`app.py`, `batch_dialog.py`) — opening one just loads it as an ordinary
  single-fit recipe. There is no bundle format, no multi-fit/session
  container, no round-trip of "data refs + fits + processing + baselines +
  batches + figures" as one file. A session today is reconstructed from
  scattered `.recipe.json` files + QSettings-remembered paths, not a single
  portable file.
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
- ❌ **Reference/literature overlays**
- ❌ **Watch-folder / auto-fit**
- 🟡 **Figure annotation layer**: a spec-level primitive exists
  (`{"x","y","text"}` → `ax.text(...)`, `render_1d` only) but there's no
  Studio UI to place/drag/edit one — usable only by hand-editing a JSON spec.
- 🟡 **Undo everywhere**: see Pillar B — real but not exhaustively audited.
- ❌ **Keyboard-driven fitting**: no evidence found.

None of these are small effort individually, but none of them are blocking
either — they're genuine "nice to have," not "can't trust the numbers."

## If only three things before 1.0 (revised)

The original three (validation report, project bundles, split app.py + CI)
are **still** the right three — re-confirmed, not superseded, by this pass:

1. **A real validation report** — ✅ substantially strengthened 2026-08-09
   (20-fit real published-paper cross-check, §5.7). Remaining gap: the
   features added since 2026-08-01 (batch fit, exclusion, population%,
   sequential fit) still have no real-dataset cross-check of their own.
2. **Project bundles + format versioning** — ❌ genuinely not started; the
   `.larproj` extension in file dialogs is a false signal that this exists.
3. **Split `app.py` + CI** — ❌ neither started; `app.py` has grown since
   this was first flagged, not shrunk, and there is no automated test gate
   on push at all.

**New fourth item, found by this pass, arguably higher priority than any
single item above because it costs almost nothing to fix and actively
misleads anyone who looks**: **fix `README.md`.** It describes a different,
much earlier program than the one in this repository, links to a file that
doesn't exist, and recommends a component (`larmor/app.py`, the web app)
that's been dormant and untested for three weeks. A colleague opening this
repo today would be actively misinformed about what LARMOR is before they
even install it — this is the fastest, cheapest fix on this entire list
relative to its impact on Pillar E ("a colleague can install and learn it
without you").

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
*Last audited: 2026-08-07 (v0.7.1), against the actual repo state, not the
previous version of this document.*
