# Development notes

What a contributor needs to know that reading the code does not tell you: the
conventions the numbers depend on, the bugs that have already been fixed and
must not come back, and the places where a plausible-looking change breaks
something silently.

Every claim here was checked against the source at the version in
`larmor/__init__.py`. File references are `path:line` at the time of writing;
treat the line numbers as hints and the file names as reliable.

---

## 1 · Orientation

**What it is.** A desktop application for fitting solid-state NMR lineshapes,
an open successor to dmfit. The physics comes from mrsimulator, the
optimisation from lmfit; LARMOR adds ingestion, the interactive UI, batch and
series workflows, uncertainties, and reproducible figures.

**Size.** 112 Python modules, ~33k lines under `larmor/`: 58 modules / ~13.8k
lines of Qt-free core, 41 modules / ~18.5k lines of desktop, plus
`larmor/xfact/` (13 modules, an easter egg). Tests: 65 files, 678 collected.

**The split.** Everything outside `larmor/desktop/` and `larmor/xfact/` is
Qt-free — verified by importing all 71 core modules and finding no PySide6 in
`sys.modules`. Nothing enforces this; it is convention, and there is no test
asserting it. The dependency direction is clean: no core module imports from
`larmor.desktop`.

**The biggest files**, where most trouble lives: `desktop/app.py` (4,677
lines — by far the largest and the main window's monolith),
`desktop/plotting_studio.py` (1,432), `desktop/qcpmg_dialog.py` (1,408),
`desktop/batchfit_dialog.py` (1,323), `figures.py` (1,134), `qcpmg.py` (992).

### Module map

| Group | Modules |
|---|---|
| Ingestion | `io/bruker.py` (1r/2rr/fid/ser, EXPNO or pdata, self-identifies 1D/2D and raw/processed), `io/varian.py`, `io/fxmla.py` (dmfit), `io/spectra.py` (CSV with a metadata header), `io/scan.py`, `io/export.py`; `loader.py` is the single entry point (`load_any`, `apply_processing`), `fourier.py` handles States/TPPI/echo-antiecho |
| Model | `recipe.py` — `Param` / `SiteModel` / `Recipe`, the diffable JSON format; data referenced by path + SHA-256, never inlined |
| Processing | `processing.py` (the replayable op pipeline), `baseline.py`, `qcpmg.py`, `qcpmg_fields.py` |
| Simulation | `models/` (the registry), `engine.py` (Czjzek kernel + `simulate`), `twod.py` (MQMAS), `estimate.py` (starting values measured from data) |
| Fitting | `fit.py`, `batchfit.py`, `seqfit.py`, `multifit.py`, `autofit.py`, `parallel.py` |
| Interpretation | `quantify.py`, `sanity.py`, `identifiability.py`, `diagnostics.py`, `chi2map.py`, `czjzek_dist.py`, `convert.py`, `nuclei.py`, `refranges.py` |
| Output | `figures.py` (spec-driven renderers), `methods.py` (auto-written Methods text), `series_grid.py` |

**Data flow.** `loader.load_any(path)` → `(ppm, amp, recipe, meta, warnings)`
→ `engine.make_context(recipe, exp_ppm)` builds a `SimContext` →
`engine.simulate_site` dispatches through the registry → `fit.fit()` returns a
`FitResult` carrying the updated `Recipe` → `Recipe.save()` → `figures.render(spec)`.

---

## 2 · Running it

```
larmor desktop                     # the app
pytest -q                          # the suite: ~663 passed, 15 skipped, ~2.5 min
pytest tests/test_qcpmg.py -q      # one file
python -m pyflakes larmor tests    # the de-facto linter
```

There is **no linter config** (no ruff/flake8/black/pre-commit) and **no CI**.
pyflakes currently reports ~55 cosmetic findings (mostly unused imports); some
are deliberate re-exports in `larmor/__init__.py`, so "fix all pyflakes" would
break the public surface.

**Environment**: Python 3.11.15, mrsimulator 1.0.0, lmfit 1.3.4, numpy 2.4.6,
scipy 1.17.1, PySide6 6.11.1, pyqtgraph 0.14.0. `pyproject.toml` pins
`>=3.10,<3.13`; mrsimulator has no wheels above that.

### Testing traps

- **`QT_QPA_PLATFORM=offscreen` and `LARMOR_NO_SESSION=1`** are set by each Qt
  test module with `os.environ.setdefault`, so they do not need exporting —
  but `setdefault` means an **inherited value wins**. If either is already
  exported in your shell, tests will use it and may open real windows or read
  your real session.
- **`LARMOR_NO_SESSION` guards five write paths**: session restore
  (`app.py:4536`), `paths.remembered_dir`/`remember_dir`,
  `_remember_site_defaults`, and the QCPMG dialog geometry. It does **not**
  protect anything a test writes to `QSettings` directly.
- **QSettings hygiene is the recurring test bug.** Copy the pattern in
  `test_desktop.py:847`, `test_ui_extras.py:766` or `test_qcpmg_dialog.py:288`:
  read the old value, write, restore in a `finally`, and `remove()` the key if
  it was absent. A test that writes `QSettings` without restoring destroys real
  user state — `test_fit_diagnostics.py::test_per_nucleus_seed` did exactly
  that to the learned `siteDefaults` until it was wrapped.

### Test data lives outside the repo

`tests/conftest.py` hard-codes eleven absolute paths under `C:\Users\samso\`,
and `tests/test_qcpmg.py:478` adds the MagLab 35Cl set. `require()` **skips**
rather than fails when they are missing, so:

- on the development machine: 663 passed, 15 skipped
- on a machine with none of that data: 639 passed, 39 skipped, **still green**

That means the entire real-data integration layer — including the
published-value acceptance tests — vanishes silently elsewhere. Any CI added
later validates far less than a local run. When judging whether a change is
safe, check that the real-data tests actually ran.

---

## 3 · Conventions the numbers depend on

These are decided once and relied on everywhere. Changing one without changing
every consumer produces plausible, wrong numbers.

| Quantity | Convention | Where |
|---|---|---|
| ppm ↔ Hz | ν[Hz] = δ[ppm] · SFO[MHz] | `convert.ppm_to_Hz` |
| Quadrupolar product | P_Q = C_Q·√(1+η²/3) | `convert.pq_from_cq_eta` |
| CT second-order shift | δ₂ = −(3/40)·[I(I+1)−¾]/[I²(2I−1)²]·(P_Q/ν₀)²·10⁶, **always negative** | `convert.ct_second_order_shift_ppm` |
| Czjzek width | LARMOR stores **σ**; dmfit's `sCZ_CQ` = **2σ**; dmfit's displayed C_Q = **4σ**; √⟨P_Q²⟩ = **√5·σ** | `desktop/table.py:74` (`CZJZEK_DISPLAYS`) |
| EFG → C_Q | C_Q[MHz] = 234.9647·Q[barn]·V_zz[a.u.] | `convert.cq_from_efg` |
| Axis | IUPAC δ, increasing to the left | `figures.py`, plot widgets |

The Czjzek factor of two is the single most dangerous number in the project:
a value copied from dmfit and stored without dividing by two makes every
lineshape twice as broad, and the fit will happily absorb it elsewhere.

**`fwhm_hz` has a units contract**: it returns `ppm_span * sfo_MHz`. The QCPMG
dialog passes the real Larmor frequency and gets Hz; both infinite-field
dialogs deliberately pass `1.0` to get a width in **ppm**, which is what
Sandland Eq. 2 needs. Changing that contract silently corrupts the width split.

---

## 4 · Validation anchors

`docs/validation.md` is the evidence document; only three of its checks are
automated, in `tests/test_physics_validation.py`:

1. **CT centroid vs analytic δ₂** — a simulated `quad_ct` ²⁷Al site's
   intensity-weighted centroid must match `ct_second_order_shift_ppm` to
   **0.5 ppm** for (C_Q, η) = (4.0, 0.3), (2.5, 0.0), (3.0, 0.8). (The document
   claims < 0.03 ppm from a finer sweep; the test is deliberately looser.)
2. **Czjzek convention relations** — mode of the C_Q marginal within
   1.7σ–2.1σ, and √⟨P_Q²⟩ = √5·σ exactly.
3. **Physical constants** — ¹H–¹H dipolar at 1.5 Å ≈ 35.6 kHz; the EFG
   constant 234.9647.

**Documented but not automated** — reproduced by hand with
`docs/figures/make_report_figures.py`, and easy to invalidate without noticing:
CT centroid vs δ₂ to < 0.03 ppm; Czjzek kernel vs a direct ensemble ≤ 0.64 %
RMSD; MQMAS δiso recovery 59.8 vs 60.0 ppm; two-field extrapolation recovering
δiso 57.8 ± 0.5 (true 58.0) and C_Q 4.16 (true 4.2); pCABS ²⁷Al 3QMAS
62.2 / 29.7 / −1.1 vs dmfit 62.7 / 30 / −0.35.

**The strongest regression net is real data.** `tests/test_qcpmg.py` reproduces
published ssNake T₂ values for **12 ³⁵Cl samples within 7 %**, with period 293
and echo top 147 on every one. `examples/pCABS2-4/` ships a real ²⁷Al and ¹¹B
dataset whose fits must stay at RMSD 0.04676 and 0.00372. Run these before
believing any change to the engine, the kernel or the fit is harmless.

---

## 5 · Do not reintroduce

Each of these was a real bug that produced confident wrong numbers. The
comments in the source explain them; this is the index.

| Trap | Where | What happens |
|---|---|---|
| **Whole-echo transform** — zero-fill must be **mid-array**, split at `m − top`, and apodization must be **circularly symmetric** about the top | `qcpmg.py:871` ("the single place the convention lives") | Appending zeros made the "absorption" spectrum 60 % dispersive; a one-sided window discarded 44 % of the echo |
| **`coadd_spectrum` must delegate to `whole_echo_ft`** | `qcpmg.py:398` | Its own copy disagreed by +53 % on a known width |
| **Two time axes** — ssNake's Split copies the acquisition SW onto the echo dimension: `T2_physical = T2_ssNake × points_per_echo` | `qcpmg.py` module docstring | A T₂ copied between programs is wrong by the echo length; never write a bare `T2` |
| **`fit_t2` needs `t_s` explicitly** when points are excluded | `qcpmg.py:685`, `qcpmg_dialog.py:794` | Survivors get re-timed onto 0, τ, 2τ…; excluding one echo shifted T₂ by −36 % |
| **`detect_period`** must blank the autocorrelation lobe by its **own measured width** | `qcpmg.py:39` | A true period of 293 came back as 8 |
| **`noise_floor` from the two edges**, never the global median | `qcpmg.py:164` | Once a wide pattern fills the window the median ran 6.2× the true floor and moved δ_CG by 44 ppm |
| **Kernel window is a floor, Cq ceiling is a ladder** | `engine.py:113–130` | A fixed 150 kHz window is only 696 ppm at 216 MHz, so a wide site renders as **all zeros**; a fixed 25 MHz grid makes Czjzek saturate |
| **`SIMULATED_DIFF_STEP = 1e-3`** for every model not in `_ANALYTIC_MODELS` | `fit.py:52–73` | At scipy's ~1.5e-8 step the Jacobian of a grid model is quantisation noise: the optimiser "converges" without moving a shape parameter |
| **`_lorentz_convolve` kernel clamp** | `quadrupolar.py:62` | `np.convolve(mode="same")` returns `max(len(a), len(v))`; an unbounded FWHM grew the array and crashed downstream |
| **`fit_2d` must refuse an all-zero model** | `twod.py:669` | Zero gradient → instant "perfect" fit of nothing with β pinned |
| **Constraint expressions must be remapped on site delete** | `constraints_util.py`, `app.py:2951` | `E = D + 5.3` silently becomes `E = E + 5.3` and recurses forever |
| **Excluded batch sites are zeroed placeholders, not omitted** | `series_grid.py:116` | Omitting shifts every later site index; index-based figure specs then colour the wrong component |
| **`setMenuEnabled(False)` destroys the pyqtgraph menu** | `plot.py:269`, `batchfit_dialog.py:559` | Custom context-menu items vanish after the first toggle; re-attach every time |
| **`keyboardTracking(False)` on every spinbox** in the QCPMG dialog | `qcpmg_dialog.py:161` | Typing "293" acts on "29" and clamps the echo top; one point of top moved T₂ by up to 7400 % |
| **`PARAM_COLUMNS` must keep its automatic fallback column** | `table.py:102` | Without it a model's parameters become fitted-but-invisible, as Amorphous ΔC_Q was |
| **`load_any` returns `(ppm, amp, recipe, meta, warnings)`** | `app.py:2331` | It was once unpacked as `(recipe, ppm, amp, …)`; every overlay format except raw Bruker silently failed |

---

## 6 · Adding a lineshape model

The registry makes this look like a one-file change. It is not — there are
**six hand-maintained tables** outside the registry, none of which cross-check
each other, and omission from each degrades behaviour silently.

1. Write `def _render(v: dict, ctx: SimContext) -> np.ndarray` and
   `register(Model(name=…, params=(ParamDef(…),…), render=_render))` in one of
   `models/{analytic,quadrupolar,csa,external}.py` — or a new module added to
   the import list in `models/__init__.py:5`, or it never registers.
2. `engine._GRID_RESTRICTABLE` (`engine.py:182`) — an **allowlist**; omission
   means full-grid simulation on every Jacobian probe (correct, just slow).
3. `fit._ANALYTIC_MODELS` (`fit.py:52`) — omission means the coarse
   `diff_step` is applied. Adding a *grid-based* model here is the dangerous
   direction: it returns the optimiser to measuring quantisation noise.
4. `estimate._WIDTH_KEY` (`estimate.py:23`) — omission means no data-driven
   starting values.
5. `desktop/table.py:PARAM_COLUMNS` — for column order (the fallback covers
   omission, but the ordering will be arbitrary).
6. `constraints_util._PEAK_FWHM_MODELS` (`:144`) — only for peak-FWHM-aware
   constraints.

**There is a hard test gate**:
`tests/test_lineshapes_help.py::test_lineshapes_manual_covers_every_model`
fails unless the model's name appears in backticks in
`larmor/help/lineshapes.md`.

---

## 7 · Limits enforced in code

Worth knowing before concluding "the model cannot fit this".

- **Kernel**: `KERNEL_MIN_SW_HZ = 150000` (a floor), `KERNEL_SPAN_MARGIN = 1.25`,
  `CQ_MAX_LADDER = (25, 50, 100, 200, 400)` MHz, npts hard-capped at 16384 in
  four places. 2D MQMAS is separate: `twod.MQMAS_SETTINGS` cq_max 16 MHz, 40×6.
- **Parameter bounds**: `CQ_MAX_MHZ = 120` (quad_ct, quad_first, quad_csa,
  ext_czjzek); `AMORPH_CQ_MAX = 6.0`; CSA ζ ±1000 ppm; ≤32 sidebands per side;
  J-multiplicity ≤12. `czjzek.sigma_Cq_MHz` has **no upper bound** and can run
  away to absurd values, degenerating into a plain Gaussian.
- **Interop**: `io/fxmla.py` converts only three dmfit line models (CzSimple,
  Gaus/Lor, Amorphous) and skips the rest with a warning.
  `refranges.py` covers exactly 8 nuclei and gives a status hint, never a guess,
  for anything else.
- **Save protection**: `desktop/paths.py:20` blocks only *replacing* an
  existing acquired file by exact name. New files anywhere, including inside
  EXPNO folders, are allowed by design.

---

## 8 · Known defects and loose ends

Ordered by how likely they are to mislead someone.

1. **Two Computing-parameters controls do nothing.** "Cq max (MHz)"
   (`desktop/dialogs.py:227`) writes `KERNEL_SETTINGS["cq_max_MHz"]`, whose
   only reader is `make_context` — which then *discards* that kernel and keeps
   just its axis. The real render paths compute their own ceiling via
   `kernel_cq_max`. "eta steps" (`:229`) is never passed to the Czjzek render
   at all, so `build_kernel`'s default `n_eta = 11` always wins. Both look
   like working controls.
2. **`make_context` builds a full kernel just to get an axis**, so a
   cold-cache Czjzek fit pays for at least two kernel builds.
3. **QCPMG provenance is dropped on "Send to fit".** The dialog emits 21
   `qcpmg_*` keys; the only receiver (`app._fid_to_workbench`) reads five of
   them and `Recipe` has no field to hold the rest. The processing record
   survives only through "Copy CSV" — which itself omits `p2_deg`, the split
   offset and the realign flag.
4. **Two dmfit amplitude constants point in opposite directions** and were
   calibrated on different single fits: import uses 3.55
   (`io/fxmla.py:284`), export uses 3.92 (`io/export.py:131`). No test
   exercises an import → export round trip, and the file the export constant
   was calibrated on is not in the repo and no longer on the machine.
5. **Version fields are write-only.** `RECIPE_VERSION` is written and then
   discarded on load; the project bundle's version is a bare literal never
   read. There is no migration hook for a future schema change.
6. **`AMORPH_CQ_MAX` is duplicated as a literal** in `io/fxmla.py:385`;
   raising the model's bound would leave imported dmfit Amorphous lines
   truncated at 6 MHz.
7. **η = 0.7 is written three times** with no shared constant
   (`qcpmg_fields.py:39`, both field dialogs). The dialogs pass their spinbox
   explicitly, so behaviour follows the UI — but the three can desynchronise.
8. **`figures.py` does not close figure handles**, triggering matplotlib's
   ">20 figures" warning in the studio tests; a memory-growth risk in long
   sessions.
9. **`docs/validation.md:452`** still advises raising `cq_max` for large σ in
   1D — stale since the ladder replaced the fixed ceiling.
10. **`README.md:111`** says "about 600 tests"; the collected count is 678.

Genuinely open work is in `docs/roadmap.md`. The largest structural items:
`app.py` is a 4.7k-line monolith, and there is no CI — so every "suite green"
claim is one machine, one environment, with the real-data tests present.

---

## 9 · The QCPMG / wide-line subsystem

The most actively developed area, and where most of the recent bug fixes came
from. Three core modules — `qcpmg.py` (992 lines), `qcpmg_fields.py` (253) —
and three dialogs: `qcpmg_dialog.py` (six stages, 1,408 lines),
`qcpmg_fields_dialog.py`, `qcpmg_batch_dialog.py`.

**The chain**: `echo_period_from_meta` → (`find_period_by_correlation` if the
period was guessed) → `centre_offset` → `split_echoes` → `echo_top_point` →
`echo_decay` → `fit_t2` → `sum_echoes` → `whole_echo_ft` → `autophase_best` →
`phase_spectrum` → `cg_window` → `centre_of_gravity` / `fwhm_hz`.

**What the dialog decides for you on load**, and how to override each:

| Decision | Rule | Override |
|---|---|---|
| Period | `CNST7 → CNST8 → CNST11 → CNST15 → CNST14 → MASR`, each range-checked to 8…n/2. If the source was MASR or none, the period is **measured** from echo-repeat correlation and replaces the guess when the data disagrees | type it, or **Find period** |
| Split offset | `centre_offset` returns the points to skip so each block holds one *centred* echo; 0 when the top is already within 25 % of centre | the offset field, or **Centre echo** |
| Echo top | `argmax\|mean(echoes)\|` — the **coherent** average, not the mean of magnitudes | drag the marker, or **Auto** |
| Phase | `autophase_best` fits p0/p1, then p0/p1/p2, keeping the quadratic only if it cuts the negative area by >25 % | type p0/p1/p2, or tick **magnitude (mc)** |
| δ_CG window | re-seeded from `cg_window` on every recompute **until you drag it**, after which the latch holds for that dataset | drag the region; **Auto window** re-seeds without clearing the latch |

**A trap in that table**: changing the period resets the split offset to 0,
clears exclusions and resets the echo counts. Because **Find period** goes
through `period.setValue()`, pressing it also wipes the offset — press
**Centre echo** afterwards.

**Measured vs fitted.** Measured (no model): period, echo top, split offset,
δ_CG and its σ (spread over nine jittered window-edge combinations), FWHM,
noise floor. Fitted: T₂ (`C + B·exp(−t/T₂)`), the phase, and the
infinite-field line. `fit_t2.ok` means **usable**, not "curve_fit returned" —
it is False for too few points, no dynamic range, a pinned bound, r² ≤ 0.5, or
err ≥ 0.5·T₂.

**What is defensible to publish**, per the help file: for a broad,
distribution-dominated pattern at a *single field*, a lineshape fit cannot
separate δiso from the second-order quadrupolar shift, so quote **δ_CG and the
central-band width** as primary and treat a fit as supporting information. For
the multi-field route, C_Q depends on the **assumed η** — quote P_Q when η is
unknown.

**Two cross-checks that make results trustworthy**: the sum-echo envelope must
trace the spikelet tops (stage 5 overlays both), and the magnitude and
p2-phased δ_CG must agree (−310 vs −314 ppm on the real ⁸¹Br WCPMG dataset).

**Documentation gaps** in `larmor/help/qcpmg.md` (accurate on all substantive
numerics, but incomplete): it does not mention the stage-2 "realign echo tops"
checkbox, stage 6's "→ infinite-field δiso" button (the primary route into the
extrapolation), the "Export figure package" button, the automatic autophase on
load, or that changing the period resets the offset. Its "all echoes overlaid"
plots actually draw the first 60.

---

## 10 · Working conventions

- **Version lives in two places** — `larmor/__init__.py` and `pyproject.toml`.
  Bump both.
- **Commits** are `vX.Y.Z: <summary>` with a body that reads as a post-mortem:
  what was wrong, how it was found, what the measured effect was. That history
  is the project's real bug documentation — `git log` is worth reading before
  touching an area.
- **Push** to `origin master` (github.com/sams808/LARMOR).
- **Before believing a change is safe**: run the full suite, and separately
  confirm the real-data tests ran rather than skipped. For anything touching
  the engine, kernel, fit or QCPMG paths, re-check the 12-sample ³⁵Cl
  acceptance test and the two shipped example fits.
