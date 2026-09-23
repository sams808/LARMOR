# Workplan status — 2026-09-22 session handoff

Where the 50-item improvement workplan stands, what this session changed,
and what a continuing session picks up. The workplan itself (items,
evidence, Sam's keep/defer/drop decisions) lives in the artifact "LARMOR
Workplan"; this file is the ground truth for progress.

## State of the repository

- **Version**: 0.13.0 (`larmor/__init__.py` + `pyproject.toml`, bumped in step).
- **Pushed** to `origin/master` on 2026-09-22 on Sam's go-ahead (62 commits,
  v0.10.1 … v0.13.0). The three commits after the full-suite run touch only
  the frozen-exe start-up guard (`_install_faulthandler`) and this file.
- **Tests**: 1004 collected in 89 files. Last full run, at v0.13.0
  (7ec83e7): **1004 passed / 0 failed** in 9 min 26 s with the real-data
  layer complete. Trust a green bar only when the terminal banner says
  "real-data layer: complete (all 14 datasets present)"
  — two datasets were added this session (the MagLab ⁸¹Br WCPMG set and the
  LAW ¹¹B series used by the tutorials).
- **Environment**: conda env `larmor` (`C:\Users\samso\xraylarch\envs\larmor`,
  Python 3.11); run tests with that interpreter. Real test data roots at
  `LARMOR_TEST_DATA` (default `C:\Users\samso`); the CaAlGlass fxmla files
  live in `Desktop\larmor_tests\`; the ⁸¹Br static set is
  `Desktop\WSU_work\NMR\MagLab\DATA\81Br_2026-08\{30..34}`.
- **Installer**: `dist/LARMOR/` (350 MB) rebuilt at 0.13.0 from the pip venv
  at `packaging/.buildenv` (python.org 3.11 — never the conda env; reasons in
  `packaging/README.md`) and smoke-run offscreen for 40 s with empty crash
  logs; the archive viewer confirms every `desktop/mw_*.py`, `workers.py` and
  the new core modules are in the bundle. The first rebuild caught a
  start-up crash present in every frozen build since 62bbd6b
  (`faulthandler.enable()` with no stderr; fixed in 099c09f/510f342). Still
  owed: verification on a machine with no development setup, and publishing
  a release (zip `dist/LARMOR/` or wrap it with Inno Setup).

## Done — the workplan is closed

All 41 kept items are implemented. Phase 5, finished this session:

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

Beyond the workplan, from Sam's requests and from what the work uncovered:

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
  windows and could abort the interpreter (this session's last commit).

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
- **Herzfeld–Berger intensities**: average over the tensor azimuth γ and β;
  the rotation about the rotor axis is a pure time shift. The second moment
  Σ I_N (Nν_r)² = (ζν₀)²(1+η²/3)/5 is the check.
- **The χ² profile levels** are χ²_min(1 + 1.00/dof) and χ²_min(1 + 3.84/dof):
  the fit is unweighted, so Δχ² must be scaled by the residual variance.
- **Indirect referencing**: SF_X = SF_¹H · Ξ_X/Ξ_¹H with mrsimulator's ratios
  reproduces TopSpin `xiref` to 0.14 Hz; adamantane ¹H is at 1.82 ppm in
  this group; a session is one month folder; the audit log is append-only.
- **Threads**: every QThread the window owns is joined in `closeEvent`;
  `tests/conftest.py` sets `LARMOR_NO_KERNEL_WARM=1` so fixture windows never
  start a kernel pre-build (a QThread destroyed while running aborts Python).
- **Shortcuts**: the palette owns Ctrl+Shift+P; Save project is Ctrl+Alt+S;
  Ctrl+P drag-to-phase; Ctrl+T FID toggle; Ctrl+I channel; Ctrl+Shift+D
  sideband detect. `tests/test_app_split.py` and the palette test assert
  uniqueness.
- **Adding a lineshape model** fails loudly if under-declared:
  `tests/test_model_tables.py` (four partitions) plus the `lineshapes.md`
  backtick gate; parameter names must be unique across models (the `function`
  model owns `d`, hence `czjzek_d`).
- **Writing files with backslashes from this tooling**: never through a
  shell heredoc (it collapsed doubled backslashes into control bytes once, in
  a display equation); write a Python script or use an editor tool.
- **A green pytest bar is not enough**: read the real-data banner.

## Remaining

1. **E7** — `dist/LARMOR/` is rebuilt and smoke-run at 0.13.0; verify it on a
   machine with no development setup (a student's laptop: unzip, double-click
   `LARMOR.exe`, open a Bruker folder, fit, save a recipe), then publish the
   release.
2. Follow-ups noted during the session, none blocking:
   - `_update_sn` re-runs sideband detection and rebuilds the literature
     overlay on every live processing apply (about 15 ms at 64k points);
     trim it if drag-to-phase feels sluggish on wideline data.
   - `czjzek_corr` and `exchange2` export to dmfit as the commented
     Gauss/Lorentz envelope; no 2D (MQMAS) path for the three new models.
   - Manual baseline is not undoable; `save_recipe` shows its status twice.
   - Two test runs stalled at `tests/test_fit_health_ui.py::test_fit_done_feeds_strip…`
     under heavy parallel load before the thread-join fix; if it recurs the
     new `faulthandler_timeout` prints the stack after 10 minutes.
3. Deferred by Sam (revisit only if asked): A9 C4 C6 E2 E3 F1 F4.
   Dropped: A2, C5.

## Conventions that keep this project safe

- Scope test runs to what changed; **full suite + real-data banner before any
  push**, never push without Sam's word.
- Version bumps touch `larmor/__init__.py` **and** `pyproject.toml`.
- Commits are `vX.Y.Z: summary` or `<item>: summary` with a post-mortem body
  (what was wrong, how found, measured effect) — `git log` is the project's
  bug documentation.
- Public docs: neutral, professional voice, no first person.
- `docs/development-notes.md` is current through this session (§1 module map
  and §11 describe the split; §8 lists the closed defects).
