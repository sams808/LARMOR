# Workplan status — 2026-09-23 session handoff

Where the 50-item improvement workplan stands, what the sessions after it
changed, and what a continuing session picks up. The workplan itself (items,
evidence, Sam's keep/defer/drop decisions) lives in the artifact "LARMOR
Workplan"; this file is the ground truth for progress.

## State of the repository

- **Version**: 0.14.0 (`larmor/__init__.py` + `pyproject.toml`, bumped in step).
- **Pushed** to `origin/master` on 2026-09-23 on Sam's go-ahead ("once done
  make sure the version on github is up to date, and push if needed"):
  v0.13.0 → v0.14.0, the eleven merged branches of the next-ten batch below.
- **Tests**: 1324 collected in 107 files. Last full run, at v0.14.0: see the
  line "Full suite" at the end of this section. Trust a green bar only when
  the terminal banner says "real-data layer: complete (all 19 datasets
  present)" — five datasets were added this session (the CaF₂ / NaF magres
  files and their 2026-03 ¹⁹F standards, and the 2026-05 ³¹P EXPNO 3102).
- **Environment**: conda env `larmor` (`C:\Users\samso\xraylarch\envs\larmor`,
  Python 3.11); run tests with that interpreter, `PYTHONIOENCODING=utf-8`
  (several CLI summaries print `→`, which a cp1252 console cannot encode).
  Real test data roots at `LARMOR_TEST_DATA` (default `C:\Users\samso`).
- **Installer**: `dist/LARMOR/` (350 MB) was rebuilt at 0.13.0 from the pip
  venv at `packaging/.buildenv` (python.org 3.11 — never the conda env;
  reasons in `packaging/README.md`); not rebuilt at 0.14.0. Still owed:
  verification on a machine with no development setup, and publishing a
  release (zip `dist/LARMOR/` or wrap it with Inno Setup).
- **Full suite**: FULL_SUITE_PLACEHOLDER

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
- **A green pytest bar is not enough**: read the real-data banner.

## Remaining

1. **E7** — `dist/LARMOR/` is at 0.13.0; rebuild it at 0.14.0
   (`packaging/README.md`), verify it on a machine with no development setup
   (a student's laptop: unzip, double-click `LARMOR.exe`, open a Bruker
   folder, fit, save a recipe), then publish the release.
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
