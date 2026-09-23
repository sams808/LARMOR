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

**Size.** 140 Python modules, ~43k lines under `larmor/`: 67 modules / ~17.9k
lines of Qt-free core, 60 modules / ~23.9k lines of desktop, plus
`larmor/xfact/` (13 modules, an easter egg). Tests: 89 files, ~1000 collected
(counts as of 0.13.0).

**The split.** Everything outside `larmor/desktop/` and `larmor/xfact/` is
Qt-free — `tests/test_core_qt_free.py` imports every core module and fails if
PySide6 or pyqtgraph appears in `sys.modules`. The dependency direction is
clean: no core module imports from `larmor.desktop`.

**The biggest files**, where most trouble lives: `desktop/batchfit_dialog.py`
(1,765), `desktop/plotting_studio.py` (1,449), `desktop/qcpmg_dialog.py`
(1,411), `figures.py` (1,134), `qcpmg.py` (992). The main window is no
longer one file: `desktop/app.py` is a 387-line facade and its behaviour
lives in eleven `desktop/mw_*.py` mixin modules of 207–920 lines each plus
`desktop/workers.py` (see the Desktop shell row below and §11).

### Module map

| Group | Modules |
|---|---|
| Ingestion | `io/bruker.py` (1r/2rr/fid/ser, EXPNO or pdata, self-identifies 1D/2D and raw/processed), `io/varian.py`, `io/fxmla.py` (dmfit), `io/spectra.py` (CSV with a metadata header), `io/scan.py`, `io/export.py`, `io/bundle.py` (series publication bundle); `loader.py` is the single entry point (`load_any`, `apply_processing`), `fourier.py` handles States/TPPI/echo-antiecho |
| Model | `recipe.py` — `Param` / `SiteModel` / `Recipe`, the diffable JSON format; data referenced by path + SHA-256, never inlined. `project.py` — the `.larproj.json` bundle: schema version, the v1 → v2 migration (through `recipe.run_migrations`), the per-kind entry builders (1D embedded; 2D maps and batch spectra by reference), path relocation |
| Processing | `processing.py` (the replayable op pipeline), `baseline.py`, `qcpmg.py`, `qcpmg_fields.py`, `sidebands.py` (autocorrelation νrot / sideband-manifold detector; mirrors qcpmg's period finder), `phasedrag.py` (drag-to-phase gesture arithmetic + pivot re-expression, Qt-free) |
| Simulation | `models/` (the registry), `engine.py` (Czjzek kernel + `simulate`), `twod.py` (MQMAS), `estimate.py` (starting values measured from data) |
| Fitting | `fit.py`, `batchfit.py`, `seqfit.py`, `multifit.py`, `autofit.py`, `parallel.py` |
| Interpretation | `quantify.py`, `sanity.py`, `identifiability.py`, `diagnostics.py`, `fithealth.py` (one verdict from the previous four, rendered by the desktop strip), `chi2map.py`, `czjzek_dist.py`, `convert.py`, `nuclei.py`, `refranges.py`, `dft.py` (magres → sites, the `[calculation]` header, equivalent-atom grouping, the `_SEED_KEYS` partition), `shiftcal.py` (σ → δ calibration line with covariance) |
| Output | `figures.py` (spec-driven renderers), `methods.py` (auto-written Methods text), `series_grid.py` |
| Desktop shell | `desktop/app.py` is the `MainWindow` facade (construction, the two Qt event overrides, `main()`); every other method is defined on one mixin in `desktop/mw_*.py` — `mw_menus`, `mw_chrome`, `mw_files`, `mw_session`, `mw_overlays`, `mw_editing`, `mw_sidebands`, `mw_fitting`, `mw_processing`, `mw_cofit`, `mw_tools` — and the QThreads are in `desktop/workers.py`. The layout and its rules are in §11 |

**Data flow.** `loader.load_any(path)` → `(ppm, amp, recipe, meta, warnings)`
→ `engine.make_context(recipe, exp_ppm)` builds a `SimContext` →
`engine.simulate_site` dispatches through the registry → `fit.fit()` returns a
`FitResult` carrying the updated `Recipe` → `Recipe.save()` → `figures.render(spec)`.

---

## 2 · Running it

```
larmor desktop                     # the app
pytest -q                          # the suite: ~1000 tests, ~15 min with the real data
pytest tests/test_qcpmg.py -q      # one file
python -m pyflakes larmor tests    # the linter; must print nothing
```

There is **no linter config** (no ruff/flake8/black/pre-commit) and **no CI**.
pyflakes must report nothing: `tests/test_pyflakes_gate.py` runs it over both
trees and fails on the first finding (the sixty cosmetic findings were curated
once in 0.12.x, and two of them were real). Deliberate re-exports are declared
through `__all__`, side-effect imports go through `importlib.import_module`.
A test that runs longer than 10 minutes dumps every thread's stack
(`faulthandler_timeout` in `pyproject.toml`).

**Environment**: Python 3.11.15, mrsimulator 1.0.0, lmfit 1.3.4, numpy 2.4.6,
scipy 1.17.1, PySide6 6.11.1, pyqtgraph 0.14.0. `pyproject.toml` pins
`>=3.10,<3.13`; mrsimulator has no wheels above that.

### Testing traps

- **`QT_QPA_PLATFORM=offscreen` and `LARMOR_NO_SESSION=1`** are set by each Qt
  test module with `os.environ.setdefault`, so they do not need exporting —
  but `setdefault` means an **inherited value wins**. If either is already
  exported in your shell, tests will use it and may open real windows or read
  your real session.
- **`LARMOR_NO_SESSION` guards the session read AND write paths**
  (`_restore_session`, `_persist_session`/`_flush_session` — the write
  side was unguarded until 0.11.1, so every test that called `snapshot()`
  polluted the developer's real session), plus
  `paths.remembered_dir`/`remember_dir`, `_remember_site_defaults`, and
  the QCPMG dialog geometry. It does **not** protect anything a test
  writes to `QSettings` directly. The session itself now lives in
  `%LOCALAPPDATA%/LARMOR/session.json` (debounced 5 s), not the registry.
- **QSettings hygiene is the recurring test bug.** Copy the pattern in
  `test_desktop.py:847`, `test_ui_extras.py:766` or `test_qcpmg_dialog.py:288`:
  read the old value, write, restore in a `finally`, and `remove()` the key if
  it was absent. A test that writes `QSettings` without restoring destroys real
  user state — `test_fit_diagnostics.py::test_per_nucleus_seed` did exactly
  that to the learned `siteDefaults` until it was wrapped.
- **Patch a module global on the module whose code reads it.** After the
  G5 split a `MainWindow` method reads its globals from the `mw_*.py`
  module that defines it, not from `larmor.desktop.app`:
  `test_kernel_warm_worker_spin_gate_and_dedup` patches `KernelWarmWorker`
  on `larmor.desktop.mw_files` (the module of `_warm_kernel`). A patch on
  `larmor.desktop.app` is silently ignored and the real worker thread runs.

### Test data lives outside the repo

`tests/conftest.py` roots every real dataset at **`LARMOR_TEST_DATA`**
(default: the original dev machine's home) and keeps the full manifest in
`ALL_DATASETS`; the MagLab 35Cl set rides the same root. `require()` **skips**
rather than fails when data is missing, so:

- on the development machine: 663 passed, 15 skipped
- on a machine with none of that data: 639 passed, 39 skipped, **still green**

The real-data integration layer — including the published-value
acceptance tests — vanishes elsewhere, but no longer silently: a
`pytest_terminal_summary` hook prints a red **"LARMOR real-data layer:
INCOMPLETE"** banner naming each missing dataset (and a green "complete" line
when all are present). The banner earned its keep the day it landed: both
CaAlGlass fxmla acceptance files had been moved to `Desktop/larmor_tests/`
and their tests had been skipping silently on the dev machine itself. When
judging whether a change is safe, check the banner.

---

## 3 · Conventions the numbers depend on

These are decided once and relied on everywhere. Changing one without changing
every consumer produces plausible, wrong numbers.

| Quantity | Convention | Where |
|---|---|---|
| ppm ↔ Hz | ν[Hz] = δ[ppm] · SFO[MHz] | `convert.ppm_to_Hz` |
| Quadrupolar product | P_Q = C_Q·√(1+η²/3) | `convert.pq_from_cq_eta` |
| CT second-order shift | δ₂ = −(3/40)·[I(I+1)−¾]/[I²(2I−1)²]·(P_Q/ν₀)²·10⁶, **always negative** | `convert.ct_second_order_shift_ppm` |
| Czjzek width | LARMOR stores **σ** (mrsimulator's: the std of each EFG component); the width in Czjzek's formula is **σ_Cz = dmfit `sCZ_CQ` = 2σ**; dmfit's displayed C_Q = **4σ** ≈ the mode of \|C_Q\| (3.73σ, exactly the mode of P_Q); √⟨P_Q²⟩ = **2√5·σ** = √5·σ_Cz (general d: √d·σ_Cz) | `czjzek_dist.py`, `desktop/table.py` (`CZJZEK_DISPLAYS`); pinned to mrsimulator's weights in `tests/test_physics_validation.py` |
| EFG → C_Q | C_Q[MHz] = 234.9647·Q[barn]·V_zz[a.u.] | `convert.cq_from_efg` |
| Axis | IUPAC δ, increasing to the left | `figures.py`, plot widgets |
| Processing replay source | raw fid only when a time-domain op precedes the first `ift` (`processing.chain_start_domain`, `loader.apply_processing`); the Processing panel emits the ABSOLUTE chain from its widgets, so it is synced from `recipe["processing"]` (`ProcessingPanel.sync_from_ops`, unrepresentable frequency-domain steps carried) before a forced re-apply | `processing.py`, `loader.py`, `desktop/panels.py` |

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
   3.5σ–4.0σ (≈ 3.73σ), √⟨P_Q²⟩ = 2√5·σ exactly, and
   `czjzek_dist.czjzek_weights(σ, d = 5)` ≡ `CzjzekKernel.weights(σ)` on the
   kernel grid to 1e-6 (the pin that ties every read-out to the distribution
   the fit actually uses).
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
| **Constraint expressions must be remapped on site delete** | `constraints_util.py`, `desktop/mw_editing.py` `_remap_exprs_after_delete` | `E = D + 5.3` silently becomes `E = E + 5.3` and recurses forever |
| **Excluded batch sites are zeroed placeholders, not omitted** | `series_grid.py:116` | Omitting shifts every later site index; index-based figure specs then colour the wrong component |
| **`setMenuEnabled(False)` destroys the pyqtgraph menu** | `plot.py:269`, `batchfit_dialog.py:559` | Custom context-menu items vanish after the first toggle; re-attach every time |
| **`keyboardTracking(False)` on every spinbox** in the QCPMG dialog | `qcpmg_dialog.py:161` | Typing "293" acts on "29" and clamps the echo top; one point of top moved T₂ by up to 7400 % |
| **`PARAM_COLUMNS` must keep its automatic fallback column** | `table.py:102` | Without it a model's parameters become fitted-but-invisible, as Amorphous ΔC_Q was |
| **`load_any` returns `(ppm, amp, recipe, meta, warnings)`** | `desktop/mw_files.py` `_load_source_body`, `mw_overlays.py` `_read_overlay_source`, `mw_editing.py` `add_background_spectrum`, `mw_cofit.py` `_cofit_load` | It was once unpacked as `(recipe, ppm, amp, …)`; every overlay format except raw Bruker silently failed. The same slip survived in `add_background_spectrum` until it was found in use: the recipe dict landed in `amp`, `np.asarray(..., float)` raised a `TypeError` out of the Qt slot, and no background/reference spectrum could ever be added. Grep every `_load_any(` call site when this shape changes |
| **Never embed 2D arrays in a project bundle** — a 2D entry is a source path (+ project-relative path) plus the view's op log | `project.py` (`entry_2d`, `view2d_persisted`), `twod_view.py` (`_ops`) | A 2rr is 10⁵–10⁶ points, four quadrants when hypercomplex; the bundle would dwarf the data it references, and the recorded ops replay exactly (`twod.replay_ops` dispatches to the same functions the view calls live) |
| **Every `Contour2DView` method that replaces `_orig` / `_committed` must append to `_ops`**; the phase-Reset rule lives in one place (`_ops_after_reset`: drop phases after the last rebasing op, keep shifts) | `twod_view.py` | `replay_ops` raises on an unknown op, but an operation that never logs itself reopens as a map that silently differs from what was saved. `_shift_axes` relabels `_orig` and `_committed` without rebasing, so a phase before a calibrate is still undone by Reset while the shift survives |
| **v2 2D / figure / batch bundle entries carry no `exp_ppm`** | `batch.py` (`load_entries`), `project.py` | The batch report's bundle reader skips them through its existing `not ppm` guard; giving them arrays would turn them into phantom fits. `load_bundle` refuses a JSON object without a `workspaces` list so a recipe file never migrates into an empty project |
| **Phase pivot fraction is computed on the frequency axis** (`SpectrumView._freq_x`), never on the displayed trace | `plot.py` `phase_pivot_frac` / `show_phase_pivot` / `_snap_peak` | With the FID or the imaginary channel shown, `_exp.xData` is milliseconds or the wrong channel: every live tick would phase about a garbage pivot |
| **`op_ift` derives `sw_Hz` and parks the axis in `x_ppm_hold`; `op_ft` anchors on `hold[n//2]`** (the fftshift zero bin), not the mid-span | `processing.py` `op_ift` / `op_ft` | A `from_processed` spectrum has `sw_Hz = 0` (every window then divides by zero, silently — numpy gives inf/NaN) and the ift→ft round trip landed at 0 ppm; the mid-span centre is half a bin off for even n |
| **Capturing a pipeline stage needs `dataclasses.replace(s, y=s.y.copy())`** | `desktop/mw_processing.py` `apply_processing` | Ops mutate the `Spectrum1D` in place and `op_ft` reassigns `y` / `x_ppm` / `domain` on the same object — a plain reference to the "FID" becomes the spectrum |
| **`hilbert` must precede `ift` before any window** on a real-only spectrum | `desktop/panels.py` `_emit` (forced and locked in re-apodize mode) | The IFT of a real spectrum is two-sided (hermitian); a one-sided EM window damps the mirrored half, loses ~half the signal and distorts the line — the whole-echo trap in a new guise |
| **`FidDialog` puts `ft` BEFORE the `phase` step** | `fid_dialog.py` `_chain` | `fourier.ft1d` appends `ft` at the END when none is given; with the phase op ahead of it every non-zero p0 / p1 made the preview fail on time-domain data |
| **A reopened `.json` recipe seeds `_proc_base` from `load_any(path, replay=False)`** | `desktop/mw_files.py` `_load_source_body` | The exp arrays arrive ALREADY replayed; seeding the live pipeline's base from them compounds the recorded chain on the first panel touch (p0 40 became 80) |
| **The magres import routes tensor components by ROLE through `dft._SEED_KEYS`**, never by the shared parameter name `eta` | `dft.py` `to_site_dict` | The quadrupolar η landed in `csa_mas`'s shielding `eta` and a spin-1/2 nucleus got C_Q = 3 MHz; guarded by `tests/test_dft_simpson.py::test_to_site_dict_routes_shielding_and_quadrupolar_eta_by_model` and the fifth partition test |

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

**There are hard test gates**:
`tests/test_lineshapes_help.py::test_lineshapes_manual_covers_every_model`
fails unless the model's name appears in backticks in
`larmor/help/lineshapes.md` — and since 0.11.2, `tests/test_model_tables.py`
fails unless the model is explicitly placed in BOTH sides of four of the
tables above (`fit._ANALYTIC_MODELS`/`_SIMULATED_MODELS`,
`engine._GRID_RESTRICTABLE`/`_GRID_FULL_REQUIRED`,
`estimate._WIDTH_KEY`/`_NO_WIDTH_SEED`,
`constraints_util._PEAK_FWHM_MODELS`/`_NOT_PEAK_FWHM_MODELS`), so the silent
omissions in items 2–4 and 6 of the list above are no longer possible. A
fifth partition, `dft._SEED_KEYS` / `dft._NOT_SEEDABLE`, decides what the
DFT tensor import seeds — a new model must declare which of its parameters
take the shielding ζ / η_CS and C_Q / η_Q (or that a computed tensor cannot
seed it) before Tools ▸ Import DFT tensors offers it.

7. **Model-name tuples NOT covered by the partition tests** — each is a
   hand-maintained allowlist that degrades silently when a new model is
   left out: `desktop/mw_editing.py` `_seed_nucleus_defaults` (per-nucleus σ/dCS seeds),
   `_sim_busy_on` (kernel-build wait message), `show_czjzek_dist` and the
   `czjzek_dist_dialog.py` site filter (P(C_Q) dialog), `_MODELS_2D` (2D
   fit refusal list); `engine._KERNEL_AXIS_MODELS` (`needs_kernel`);
   `methods.py` `_MODEL_PHRASE` / `_COLS` / the Czjzek-sentence gate;
   `batch._site_columns`; `multifit.DEFAULT_SHARE`; `io/fxmla.py` +
   `io/export.py` (dmfit mapping); `panels.PARAM_LABELS`; cofit `_SHORT` and
   `app._COFIT_LABEL`; `help/spectra-1d.md` "Which lineshape?" table;
   `docs/validation.md` §6. Grep for an existing model name (`"ext_czjzek"`)
   before declaring a new model finished.

**Parameter names must be unique across models.** `table.PARAM_COLUMNS`,
`panels.PARAM_LABELS`, `multifit.DEFAULT_SHARE` and the cofit tie bar are
keyed by parameter NAME across all models, so two models using the same name
for different quantities share a column, a label and a tie. The `function`
model owns `a`, `b`, `c`, `d` — which is why the general-d Czjzek parameter is
`czjzek_d` (lmfit key `d` is fine: keys are prefixed per site).

---

## 7 · Limits enforced in code

Worth knowing before concluding "the model cannot fit this".

- **Kernel**: `KERNEL_MIN_SW_HZ = 150000` (a floor), `KERNEL_SPAN_MARGIN = 1.25`,
  `CQ_MAX_LADDER = (25, 50, 100, 200, 400)` MHz, npts hard-capped at 16384 in
  four places. The Czjzek family requests the ladder step covering
  `CZJZEK_KERNEL_HEADROOM = 10` × σ (`models/quadrupolar.py`; 21 % of the
  d = 5 mass lies beyond 5σ, 4.5e-5 beyond 10σ) — the app prewarm uses the
  same rule. 2D MQMAS is separate: `twod.MQMAS_SETTINGS` cq_max 16 MHz, 40×6.
- **Parameter bounds**: `CQ_MAX_MHZ = 120` (quad_ct, quad_first, quad_csa,
  ext_czjzek); `AMORPH_CQ_MAX = 6.0`; CSA ζ ±1000 ppm; ≤32 sidebands per side;
  J-multiplicity ≤12. `czjzek.sigma_Cq_MHz` ≤ **40 MHz** (10σ against the
  400 MHz ladder top; `tests/test_models_v2.py` pins the product for every
  model carrying a `sigma_Cq_MHz`) — at the bound the at-bounds diagnosis
  fires instead of the lineshape saturating into a plain Gaussian.
  `czjzek_d`/`czjzek_corr`: σ ≤ 40 MHz likewise; `czjzek_d` ∈ [1, 5]
  (pinned by default); `shift_slope_ppm_per_MHz` ±50 ppm/MHz; `exchange2`
  `k_ex_hz` ∈ [0, 1e8] s⁻¹, `pop_a` ∈ [0.01, 0.99], `split_ppm` ≥ 0.
  `czjzek_d`, `czjzek_corr` and `exchange2` have no 2D (MQMAS)
  implementation — `app._MODELS_2D` refuses them with the supported list.
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
2. ~~`make_context` builds a full kernel just to get an axis~~ **fixed in
   0.11.1**: `engine.kernel_axis_ppm` derives the axis arithmetically
   (against `Isotope.B0_to_ref_freq`, NOT `larmor_MHz` — 0.08 % apart),
   and a test asserts `make_context` never calls `build_kernel`.
3. **QCPMG provenance is dropped on "Send to fit".** The dialog emits 21
   `qcpmg_*` keys; the only receiver (`mw_files._fid_to_workbench`) reads five of
   them and `Recipe` has no field to hold the rest. The processing record
   survives only through "Copy CSV" — which itself omits `p2_deg`, the split
   offset and the realign flag.
4. **Two dmfit amplitude constants point in opposite directions** and were
   calibrated on different single fits: import uses 3.55
   (`io/fxmla.py:284`), export uses 3.92 (`io/export.py:131`). No test
   exercises an import → export round trip, and the file the export constant
   was calibrated on is not in the repo and no longer on the machine.
5. ~~Version fields are write-only~~ **fixed in 0.11.2**:
   `recipe._MIGRATIONS` runs on load with a note per hop. The migrations
   dict is no longer empty: the project bundle's v1 → v2 hop
   (`project._MIGRATIONS`, `kind` defaults on every entry) is the first
   real migration, and both tables run through `recipe.run_migrations`;
   the recipe schema itself is still v1. A bundle newer than
   `project.PROJECT_BUNDLE_VERSION` opens with a note in the "Project
   opened" box rather than a refusal.
6. **`AMORPH_CQ_MAX` is duplicated as a literal** in `io/fxmla.py:385`;
   raising the model's bound would leave imported dmfit Amorphous lines
   truncated at 6 MHz.
7. **η = 0.7 is written three times** with no shared constant
   (`qcpmg_fields.py:39`, both field dialogs). The dialogs pass their spinbox
   explicitly, so behaviour follows the UI — but the three can desynchronise.
8. **`figures.py` does not close figure handles**, triggering matplotlib's
   ">20 figures" warning in the studio tests; a memory-growth risk in long
   sessions.
9. ~~`docs/validation.md:452` still advises raising `cq_max` for large σ in
   1D~~ **fixed with G4**: the row now states the 1D rule (ladder at 10σ,
   σ ≤ 40 MHz).
10. ~~README test count stale~~ kept current (README says ~700; 712
    collected at 0.11.2).
11. ~~**Czjzek read-outs were a factor 2 low.**~~ **fixed with G4** (first
    commit). The kernel weights are mrsimulator's, whose formula carries
    `sigma_ = 2*sigma`: the stored σ is HALF the width in Czjzek's formula
    (σ_Cz = dmfit sCZ_CQ = 2σ). `czjzek_dist.py`, the P(C_Q) dialog
    ("mode 2σ"), `batch.py`'s √⟨P_Q²⟩ column, `methods.py`, the table's
    P_Q display (√5·σ) and `docs/validation.md` all described a distribution
    half as wide as the one the fit used. Measured against mrsimulator's own
    density: marginal mode 3.73σ (not 2σ), √⟨P_Q²⟩ = 4.472σ = 2√5·σ (not
    √5·σ), 21.3 % of the mass beyond 5σ, 4.5e-5 beyond 10σ; for d = 2..5 the
    rms equals √d·σ_Cz exactly (chi law). Fits, stored σ, recipes and the
    dmfit relations (sCZ_CQ = 2σ, CQ = 4σ) were right throughout; only the
    interpretation changed. `czjzek_weights(σ, 5)` is now pinned to
    `CzjzekKernel.weights(σ)` to 1e-6.
12. ~~**Czjzek kernel headroom truncated the distribution.**~~ **fixed with
    G4** (first commit). `_render_czjzek` requested a (C_Q, η) grid to 5σ,
    beyond which 21.3 % of the d = 5 mass lies; whenever 5σ sat just under a
    ladder step (σ ≈ 4–5, 8–10, 16–20 MHz — the heavy-halide range) a fifth
    of the distribution was dropped and renormalised away, biasing σ upward
    there. Requests are now 10σ (`CZJZEK_KERNEL_HEADROOM`), σ max 80 → 40,
    and the app prewarm follows the same rule. The pCABS anchors
    (²⁷Al σ ≈ 1.5, ¹¹B σ ≈ 1: 10σ ≤ 25) keep their kernel and their RMSDs.
13. ~~`_last_lmfit` survived load / workspace switch / 2D~~ **fixed with the
    fit-health strip (F7)**: the covariance of the previous fit was only ever
    written by `_fit_done`, so Decomposition ▸ Parameter correlations could
    show another spectrum's matrix. `_health_reset()` now drops verdict and
    covariance together in `_update_sn` (the active-1D-document funnel),
    `_show_2d` and `_apply_doc`'s 2D branch; workspace snapshots carry them
    in memory so a switch back restores the right ones.
14. ~~Auto Fit bypassed the diagnostics~~ **fixed with the fit-health strip**:
    `run_auto_fit` wrote only an RMSD; the winning `FitResult`
    (`AutoFitResult.result`) now goes through `_health_from_result` like a
    plain fit, so it gets the same verdict, chi text and correlations.
15. ~~Residual diagnostics never ran for kernel-model fits~~ **fixed with the
    fit-health strip**: `_fit_done` subtracted `y_fit` (on the model axis,
    16k points for a Czjzek recipe) from `y_exp` (the data axis) and swallowed
    the shape error, so the residual/noise ratio and the runs test were
    silently absent from every Czjzek / Amorphous fit's Report header.
    `fithealth.assess` interpolates the model onto the data axis first
    (`x_fit` + `ppm`). `MainWindow._residual_noise_ratio` is now a two-line
    delegate to `fithealth.residual_noise_ratio`; it stayed on
    `desktop/mw_fitting.py` through the G5 split because tests call it on the
    class, and removing it is a separate small item.
16. **`MainWindow._data2d` (what `run_fit_2d` fits) is the as-loaded map**
    and does not follow the view's phase / shear / calibrate; the 2D page's
    footer Compute is a no-op (`_simulate_now` returns early for the 2D
    view) and only `_fit2d_done` draws a 2D model. This is also why a
    reopened project restores a 2D map's fit parameters but not its model
    overlay — the MQMAS kernel build takes seconds and is memory-cached
    only, so the notes box says to run Fit.

Genuinely open work is in `docs/roadmap.md`. The largest structural item left
is that there is no CI — so every "suite green" claim is one machine, one
environment, with the real-data tests present. (The 6.7k-line `app.py`
monolith was split in G5; see §11.)

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

---

## 11 · Main-window layout (after G5)

`MainWindow` is assembled from mixins:
`class MainWindow(_MenusMixin, _ChromeMixin, _FilesMixin, _SessionMixin,
_OverlaysMixin, _EditingMixin, _SidebandsMixin, _FittingMixin,
_ProcessingMixin, _CofitMixin, _ToolsMixin, QMainWindow)`. Every method was
moved verbatim from the pre-split `app.py`; nothing user-visible changed, and
`win.<method>` / `MainWindow.<staticmethod>` / `from larmor.desktop.app import
…` all resolve as before.

| Module (lines) | Owns |
|---|---|
| `desktop/app.py` (387) | the facade: `MainWindow.__init__` (window state and build order), `keyPressEvent`, `closeEvent`, `ClickableLabel`, `asset_path`, `main()`; re-exports through `__all__` |
| `desktop/workers.py` (226) | `FitWorker`, `Fit2DWorker`, `SimWorker`, `KernelWarmWorker`, `_emit_progress`, `humanize_error`, `_fit_tol` |
| `desktop/mw_menus.py` (920) | menu bar, toolbar, sidebar, `_add()`, recent / apply-recipe rebuilders, theme / text-size / axis-unit / Czjzek-display submenus and setters, `_update_enabled`, Help slots, `TUTORIALS` |
| `desktop/mw_chrome.py` (392) | dock builders, `_fit_to_screen`, progress bar with Stop / Cancel, experiment and MAS status labels, zoom / autoscale / toggles, copy and save plot |
| `desktop/mw_files.py` (772) | open dialogs, `load_source`, explorer and 2D navigation, File > Watch, `_warm_kernel`, save recipe / fit / spectrum / figure; module functions `_load_any`, `_nmrdata_to_data2d` |
| `desktop/mw_session.py` (597) | workspaces, `.larproj.json` save / open, figure and batch session rows, the crash-recovery session file |
| `desktop/mw_overlays.py` (207) | Datasets-dock overlays, 1D-on-2D projections, `PROJ_COLOR` |
| `desktop/mw_editing.py` (704) | undo / redo, add-line mode, `_NUCLEUS_START` seeding, every `add_*`, site structure edits, paddles, static-CT and Herzfeld-Berger seeds |
| `desktop/mw_sidebands.py` (423) | manual `add_sidebands` and the F3 auto-detect offer (`_ssb_*`, `detect_sidebands`, the banner) |
| `desktop/mw_fitting.py` (552) | simulation debounce, 1D / 2D fit runs and completions, `_health_*`, quantify, CSV / LaTeX / Methods, publication bundle, auto fit |
| `desktop/mw_processing.py` (737) | `apply_processing`, FID / channel display refresh, drag-to-phase, calibrate / measure / reset, baseline tools, zones, experiment parameters, WURST, subtract |
| `desktop/mw_cofit.py` (471) | the co-fit page and its state machine |
| `desktop/mw_tools.py` (624) | dialog launchers of the Tools / Analysis menus and the small recipe edits they hand back (`_magres_add_sites` applies a DFT import: sites, note, `provenance['dft_import']`); the two-tab magres dialog itself is `desktop/magres_dialog.py` |

**Where a new `MainWindow` method goes.** In the mixin that owns the state it
reads or writes (the module docstrings list the owned attributes); a slot that
only opens a dialog goes in `mw_tools.py`; a new QThread goes in
`workers.py`. `app.py` receives nothing but construction — a "quick fix"
landed there fails `test_facade_holds_only_construction_and_qt_overrides`.

**Mixin rules** (PySide6 multiple inheritance), enforced by
`tests/test_app_split.py`: a mixin defines no `__init__`, no `Signal` and no
Qt event override (`keyPressEvent` / `closeEvent` stay in `app.py`); mixin
method names are pairwise disjoint (a duplicate would be shadowed silently by
MRO order); a mixin never imports `larmor.desktop.app` (app imports the mixins
at module level); module-level names live where they are used and `app.py`
re-exports them through `__all__` (`_load_any` in `mw_files.py`, `TUTORIALS`
in `mw_menus.py`, the workers). The same test file freezes the member
surface, checks that every `self.X(` call resolves, and pins the menu tree,
the shortcut set and the toolbar action count of an offscreen window. When a
label or shortcut changes on purpose, paste the fresh `GOLDEN_MENU` /
`SHORTCUTS` that the failure message prints.

**Noticed during the split, deliberately not changed** (every G5 commit is a
provable move): `apply_manual_baseline` takes no `snapshot()` while
`apply_twopoint_bg` does, so a manual baseline is not undoable; `save_recipe`
emits two status-bar messages back to back; the new-site `params` literal
(`{p.name: {"value": p.default, "stderr": None, "vary": p.vary, …} for p in
model_registry.get(name).params}`) is written out six times in
`mw_editing.py` although `mw_sidebands.py` already has it as
`_fresh_params(model_name)`.
