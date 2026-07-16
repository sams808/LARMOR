# LARMOR — Roadmap to a reference solid-state NMR fitting software

Goal: the tool a solid-state NMR lab opens by default — dmfit's fluency, ssNake's
processing depth, the Grandinetti stack's physics, plus what none of them have:
uncertainties everywhere, reproducible recipes, and multi-dataset fits.

Feature sources, verified against the real things:
- **dmfit 20230120** (menus captured live: Decomposition → Zones / Auto Fit /
  Errors Analysis / Report / MAS spinning; paddles; Fit-Parameters table) —
  Massiot et al., *Magn. Reson. Chem.* 40, 70 (2002).
- **ssNake source** (local copy; ops inventoried from `spectrum.py`) —
  van Meerten et al., *J. Magn. Reson.* 301, 56 (2019).
- **TopSpin** processing parameters (WDW/LB/GB/SSB, TDeff, SI, SR, PHC0/1, FCOR).
- **mrsimulator / mrinversion / csdmpy** — Srivastava et al., *J. Chem. Phys.*
  161, 212501 (2024); Grandinetti PANACEA school notes 2024.
- **SIMPSON** — Bak, Rasmussen, Nielsen, *J. Magn. Reson.* 147, 296 (2000).

Conventions: every feature ships with tests; instrument data stays read-only;
every result carries its uncertainty; recipes/specs are diffable JSON.

---

## Done (v0.2, 2026-07-15)

- Model registry: gauss_lor, czjzek (σ = dmfit sCZ_CQ/2, validated on
  CaAlGlass.fxmla, RMSD 0.0025), **ext_czjzek**, quad_ct, **quad_first**
  (satellites + 64 sidebands), **quad_csa**, csa_mas. Shared LRU single-site
  simulator.
- Native desktop app (PySide6 + pyqtgraph), dmfit-faithful: paddles
  (pos+amp square handle, width side handles), Fit-Parameters spreadsheet
  (pin = fixed), menus, zones (union-of-regions fit), interactive PCHIP
  anchor baseline, experiment dialog (νrot / Larmor / nucleus), Fusion light
  palette (Windows dark-mode proof).
- Constraints: fix / bounds / algebraic links with error propagation;
  dependent positions via dialogs in **ppm or Hz**; at-bounds diagnosis with
  pinned-covariance recovery.
- Processing: EM/GM/SINE/QSINE/TRAF, TDeff, ZF (factor or SI), FCOR,
  shift_fid, FT, phase, robust autophase (scan+NM; ACME option), **Hilbert**
  (rephase pdata), SR, magnitude, auto+manual baseline.
- **satrec module**: ser+vdlist → per-slice processing → integrals → auto T1
  (stretched option); validated r > 0.98 against TopSpin t1ints.txt on real
  19F data; desktop dialog + CLI.
- **multifit engine**: simultaneous multi-dataset fits with cross-dataset
  parameter sharing (multi-field 1D lifts the Cq/δiso degeneracy — proven in
  tests at 9.4 T + 17.6 T); CLI `larmor multifit`.
- Quantification table (% ± err), figure studio (1D/2D/series, styles,
  png+svg+pdf), fxmla import, read-only Bruker I/O, recipes, tutorials 1-3,
  52 tests.

---

## v0.3 — Multi-dataset UI & fit power tools

*The engine exists; give it the cockpit.*

- [ ] **Multi-dataset workspace in the desktop app**: open several datasets as
  tabs/overlays; a shared-parameters panel (checkbox per parameter name, per
  site) driving `larmor.multifit`; per-dataset RMSD; combined report.
  *Accept: refit the two-field synthetic case entirely from the GUI.*
- [ ] **MQMAS ∥ 1D parallel fit**: fit an MQMAS isotropic projection together
  with the 1D MAS spectrum (shared δiso/Cq per site) as a first
  multi-experiment case (dmfit's workflow, automated).
- [ ] **Auto Fit** (dmfit Decomposition > Auto Fit): multi-start heuristic —
  N random restarts within bounds, keep the best χ²; progress in the UI.
- [ ] **Errors Analysis** (dmfit): 1D χ² profile per parameter around the
  minimum (fix-and-refit scan), plotted, with the 1σ/2σ crossings marked —
  honest errors when the covariance is unreliable.
- [ ] **Undo for processing** (op-pipeline history with step removal, ssNake
  style) instead of reset-to-original only.
- [ ] Per-line notes + site colors editable; report includes the full recipe.

## v0.4 — True 2D: MQMAS fitting end-to-end

*CaAlGlassMQ.fxmla finally refits natively.*

- [ ] **2D data model**: CSDM-backed 2D spectra in the app (contour view =
  figure-studio renderer made interactive: zoom, levels slider, projections).
- [ ] **Shearing / referencing** (ssNake `shear`): 3Q→isotropic shear with
  the standard ratios per spin; F1 referencing conventions (Amoureux).
- [ ] **2D MQMAS simulation kernel**: mrsimulator `ThreeQ_VAS` per (Cq, η)
  grid point → 2D kernel; czjzek/ext_czjzek weights reweight it exactly like
  1D (same architecture, one more dimension).
- [ ] **2D residual fit** with the same registry/constraints/quantify stack;
  fxmla MQMAS import already parses — wire `FitMode=MQMAS` to this engine.
  *Accept: CaAlGlassMQ.fxmla parameters reproduced within uncertainties.*
- [ ] 2D figure export via the existing figure studio (overlay fitted
  contours on data).

## v0.5 — Processing parity (ssNake + TopSpin complete)

- [ ] Remaining apodizations: Hamming, Kaiser, shifted Gaussian with echo
  center (ssNake `wholeEcho`), JMOD.
- [ ] **Linear prediction** (LPSVD forward/backward — ssNake `lpsvd`):
  truncated-echo reconstruction, first-points repair.
- [ ] **Whole-echo processing** (swapEcho + wholeEcho apodization + magnitude)
  for Hahn-echo data like EXPNO 1903.
- [ ] 2D-ready processing: States / States-TPPI / echo-antiecho recombination
  (ssNake), F1 FT, t1-noise ridge suppression (median).
- [ ] Region tools: extract/truncate region (ssNake `extract`), align spectra
  (`align`), regrid, spectra algebra (add/subtract/scale — dmfit Dual,
  background subtraction workflow from the NMRVEW notebooks).
- [ ] Peak picking (threshold + parabolic interpolation) feeding "add line at
  every peak".
- [ ] Processing presets saved in the recipe (`processing:` block) and
  replayed on load — the fit and its processing become one reproducible unit.
- [ ] Baseline: polynomial degree UI for auto mode; sine/cosine (Bernstein)
  bases; per-zone baseline.

## v0.6 — Relaxation & series suite

- [ ] Generalize satrec into a **series engine**: inversion recovery
  (vdlist, I(t)=I0(1-2f·exp(-t/T1))), T2 CPMG (vclist), T1ρ, DOSY-style decay;
  auto-detection from the pulse program (satrec/invrec/cpmg patterns).
- [ ] **Per-site series fitting**: instead of integrating a window, fit the
  1D model to every slice with shared lineshapes and per-slice amplitudes
  (the dmfit "spectra series" mode) → T1 per SITE, not per window.
  *Accept: two overlapping sites with different T1 separated on synthetic data.*
- [ ] Arrhenius/VFT helper for variable-temperature series (plots + fits).
- [ ] mrinversion integration UI: T1/T2 **distribution** inversion
  (regularized ILT) as the model-free alternative, with the L-curve/CV shown.
- [ ] REDOR module: S0/S pairs → ΔS/S0 with the universal short-time curve
  and dipolar second moment (Bertmer & Eckert), errors included.

## v0.7 — Advanced simulation & import

- [ ] **SIMPSON bridge**: export the current spin system to SIMPSON input,
  run locally if installed, overlay/fit the returned curve (REDOR, RFDR...).
- [ ] **DFT/MD tensor import**: CASTEP `.magres`, Gaussian/ORCA EFG+CSA →
  predicted sites in the table (the design-doc feature, mrsimulator
  cartesian_tensor does the math).
- [ ] Second-field prediction: one click "what would this look like at X T?"
  from the current model (teaching + planning tool).
- [ ] Isotope database panel (γ, Q, abundance, reference compounds).

## v0.8 — Packaging & distribution

- [ ] **PyInstaller one-folder + installer** (Dataapp pattern): no conda, no
  Python — download, run. Kernel caches in %LOCALAPPDATA%.
- [ ] Auto-update check against GitHub releases.
- [ ] Tutorials 4-8 (multi-field, MQMAS, satrec, processing, SIMPSON) with
  bundled example data; help menu opens them.
- [ ] Crash reporter (local log + "copy diagnostics").

## v1.0 — Reference release

- [ ] **Validation corpus**: ≥10 published fits (dmfit examples, literature
  glasses/crystals) reproduced within stated uncertainties; results table in
  the docs.
- [ ] Performance: 1M-point spectra fluid; kernel build < 5 s warm path;
  full-suite CI on GitHub Actions (fast tests, offscreen Qt).
- [ ] Paper draft (software announcement, JMR/SSNMR style) with the
  validation corpus.
- [ ] Public issue templates, versioned file-format spec (recipe/figure).

---

### Standing engineering rules

1. Physics from mrsimulator/lmfit — never hand-rolled lineshapes.
2. A feature without a test does not merge; UI features get offscreen tests
   plus a real-display screenshot check (Windows dark mode!).
3. Instrument folders are read-only; writes into them are refused in code.
4. Every fitted number is printed with its uncertainty or with the reason
   there is none.
