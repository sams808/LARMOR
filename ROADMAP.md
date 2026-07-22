# LARMOR — Roadmap & status

Goal: the tool a solid-state NMR lab opens by default — dmfit's fluency,
ssNake's processing depth, the mrsimulator physics stack, plus what none of them
combine: uncertainties everywhere, reproducible recipes, multi-dataset fits, and
a modern desktop UX.

Feature sources (verified against the real programs, screenshots captured live):
- **dmfit** — Massiot et al., *Magn. Reson. Chem.* 40, 70 (2002). Model list,
  Decomposition menu, Computing-parameters, 2D operations, copy-to-clipboard.
- **ssNake** — van Meerten et al., *J. Magn. Reson.* 301, 56 (2019). Tools /
  Fitting / Utilities menus, QCPMG tutorial, NMR table.
- **TopSpin** processing (WDW/LB/GB/SSB, TDeff, SI, SR, PHC0/1, FCOR, t1t2).
- **mrsimulator / mrinversion / csdmpy** — Srivastava et al. 2024.
- **Literature**: Czjzek (Czjzek 1981), Extended Czjzek (Le Caër 2010),
  QCPMG (Larsen/Jakobsen/Nielsen 1997/98), Herzfeld–Berger sidebands (1980).

Conventions: every feature ships with tests; instrument folders stay read-only;
every fitted number carries its uncertainty; recipes/spectra are diffable.

---

## Verified working — audit 2026-07-21 (188 tests green)

**Import / data.** Universal read-only Bruker reader (1r/2rr/fid/ser, pdata,
EXPNO; 1D/2D, time/freq; pseudo-2D detection; hypercomplex quadrants). dmfit
`.fxmla` import (1D + MQMAS). LARMOR recipe `.json`. Plain **CSV/TXT** spectra
with metadata header. "Open any type" — nothing is rejected.

**Models (registry).** gauss_lor (pseudo-Voigt), **voigt (true)**, czjzek,
ext_czjzek, quad_ct, quad_first (satellites+sidebands), quad_csa, csa_mas (with
physical sidebands), spectrum (external background component). All carry
uncertainties; constraints = fix / bounds / algebraic links (ppm or Hz) with
error propagation.

**Processing.** EM/GM/SINE/QSINE/TRAF, TDeff, ZF, FCOR, shift_fid, FT, phase,
autophase (ACME), Hilbert, SR, magnitude, swap_echo, echo_apodize, LP, extract,
scale/offset/normalize, auto + manual (PCHIP) baseline. **Live** processing
panel with ±90/180° phase steps.

**1D fitting.** dmfit-style paddles, Fit-Parameters spreadsheet (letters, links,
bounds, pin), zones, quantification (% ± err), Auto Fit (multi-start), Errors
Analysis (χ² profile), S/N, calibrate-by-click, two-cursor measure.

**2D.** Interactive contour view (decimated, ~10 ms warm), positive/negative/both
contours, **exact hypercomplex phasing** (pick 1–2 peaks → row/column traces),
manual shear, projections, trace-to-workbench. **MQMAS fit in-app** (kernel +
czjzek/ext_czjzek/quad, fitted contours overlaid). **HMQC**: overlay 1D on
projections (picked from the Explorer, colour-coded) + un-correlated-feature
difference.

**Relaxation / QCPMG.** Guided T1/T2 (satrec/invrec/cpmg/t1rho, zones, outlier
exclusion, slice-range selection, TopSpin-parity T1). **QCPMG**: full ssNake
sum-echo workflow (split → T2 → matched-filter weighting → whole-echo →
absorption to fit) or spikelets; Help manual.

**Co-fit.** Mixed 1D + 2D (MQMAS) simultaneous fit with shared parameters
(multi-field, or MQMAS ∥ 1D), from a GUI dialog.

**Workspaces.** TopSpin-style: open a 2D / extract a trace → a new switchable
workspace (list dock, close, save); resource-light snapshots.

**Utilities.** Background subtraction (+ save reusable CSV), figure studio,
REDOR, `.magres` DFT tensor import, experiment params incl. SR.

---

## Documentation — extensive user manuals (parallel workstream ⭐)

Every feature gets a written, in-app manual (Markdown opened by a **Help**
button, like `larmor/help/qcpmg.md` already does). Manuals are **organized by
data type**, not by feature, so a user opening the manual for what's in front of
them gets everything relevant — opening, processing, fitting, and the tools for
that experiment — in one place. Shared processing content (apodization, phasing,
baseline, referencing) is **deliberately duplicated** into each data-type manual
so nobody has to hunt across documents.

Planned manuals (`larmor/help/*.md`), each with worked examples on real data:
- [ ] `getting-started.md` — the workspace model, the Explorer, opening any type.
- [ ] `spectra-1d.md` — 1D processing (WDW/LB/GB, phasing, baseline, SR/calibrate)
  and fitting (models, paddles, zones, constraints, quantify, S/N).
- [x] `qcpmg.md` — wide-line QCPMG / echo trains (done; expand as features land).
- [ ] `relaxation.md` — guided T1/T2, slice selection, T2-weighting, buildups.
- [ ] `mqmas.md` — 2D MQMAS: contours, shear/referencing, in-app fit, projections.
- [ ] `correlation-hmqc.md` — HMQC/DQ-SQ, projection overlays, un-correlated
  extraction, and the multi-experiment correlation decomposition (P0.5).
- [ ] `2d-processing.md` — 2D phasing (hypercomplex), contours, measure/calibrate.
- [ ] `multi-dataset.md` — overlays/compare, co-fit (multi-field, 1D + MQMAS).
- [ ] `models-reference.md` — every lineshape, its parameters, and when to use it.
- [ ] `processing-reference.md` — the full op list (shared source for the
  per-data-type sections above).

Standing rule addition: **a feature merges with its manual section.** A short
docs index lists them; the Help button in each tool opens the matching manual.

---

## Priority 0 — Lineshape rework  ⭐ (the user's headline)

The model catalogue works but is not yet at dmfit breadth, and the width /
amplitude conventions deserve a cleanup. This is the next big block.

**Architecture.**
- [x] **Width entry in Hz *or* ppm per line**: the cell accepts `300Hz` /
  `1.5kHz` / `2ppm` (and on MHz params, `3000kHz`). *(done)*
- [x] **Amplitude = area** via the **GL-Norm** model (amplitude is the integral).
  *(done)* — [ ] a per-line peak/area toggle on any model is still open.
- [ ] **Per-dimension lineshape for 2D** (dmfit "GL/F1"): independent F1/F2
  widths and shapes on a fitted MQMAS site.
- [ ] A **line base** so a site declares `lineshape ∈ {gauss, lorentz, voigt,
  pseudo-voigt}` instead of separate models where it's just a shape choice.

**New analytic / physical models (dmfit parity, ranked by usefulness).**
- [x] **Voigt (true)** — Gaussian ⊗ Lorentzian, independent widths. *(done)*
- [x] **GL Norm** — area-normalized Gauss/Lorentz for clean quantification.
  *(done)*
- [x] **J-multiplet** — binomial (n+1)-line pattern split by J (Hz). *(done)*
- [x] **Spinning sidebands (generic "ss band")** — centre + sidebands at ±k·νrot
  with a geometric intensity ratio (`sidebands` model). *(done)*
- [ ] **J-multiplet** (dmfit "Jmultiplet") — n equivalent couplings → binomial
  multiplet with a J (Hz) and a per-component lineshape; and the **J-dispersion**
  variants (residual dipolar / distribution of J).
- [x] **Function fit** (ssNake) — a user-typed `y(x; a,b,c,d)` expression, safely
  evaluated (Decomposition ▸ Add function line). *(done)*
- [ ] **Shape-from-file / user shape** — generalize the current `spectrum`
  component into a first-class shifted+scaled+broadened basis line (dmfit "W
  Shape file" / "user shape").

**Disorder / distribution physics (beyond dmfit).**
- [x] **CSA Czjzek** (disordered shielding): CSA pattern averaged over a Gaussian
  ζ distribution (`csa_czjzek`). *(done)*
- [ ] **Correlated δiso–Cq distribution** for glasses (a 2D Gaussian on the
  (δiso, Cq) plane) — closer to real amorphous lineshapes.
- [ ] **Extended Czjzek in 1D** already exists; expose ρ/ε clearly and validate
  against Le Caër's figures.

*Accept for P0: reproduce a published glass fit (e.g. CaAlGlass variants) with
Voigt + GL-Norm + sidebands, and a J-multiplet crystal, within uncertainties.*

---

## Priority 0.5 — HMQC & multi-experiment correlation  ⭐⭐ (a very big deal for us)

Heteronuclear-correlation processing is central to our work; the current overlay
+ un-correlated-difference is only the start.

**Improve HMQC processing (near-term).**
- [x] **Robust projection ↔ 1D scaling**: least-squares "fit scale" over the
  visible range, not only peak-match. *(done)* — [ ] still: per-peak / region UI.
- [ ] **Projection choice**: skyline / sum / integral, and an **external
  projection** (use a separately-acquired 1D as the projection reference) on
  either F1 or F2.
- [ ] **Indirect-dimension processing**: F1 apodization / phasing / zero-fill and
  **F1 referencing** (the indirect nucleus' SR), plus **t1-noise ridge
  suppression** (median) — HMQC F1 is usually the ugly axis.
- [ ] **Un-correlated extraction polish**: keep/normalize options, error bars on
  the difference, and send *both* the correlated and un-correlated parts to
  separate workspaces.
- [ ] Manual: a dedicated **HMQC / correlation** guide (see Documentation).

The **generalized multi-experiment correlation decomposition** that grows out of
this is deliberately scheduled **last** (see "Priority 9", below) — the engine
architecture may be built ahead of time, but it is not wired into the app yet.

---

## Priority 1 — Utilities & tools (ssNake Utilities menu)

Small, self-contained, high daily value.
- [x] **NMR table** — interactive periodic table of Larmor frequencies at a
  settable B0 (or set the ¹H frequency of *your* magnet). Double-click an element
  → isotopes with spin, abundance, γ, Q, receptivity. *(done)*
- [x] **Chemical-shift / Cq / dipolar-distance** conversion tools *(done)*;
  [ ] **MQMAS-parameter extraction** (δ1/δ2 → δiso, PQ) still pending
  (convention-sensitive; ships when validated against a literature example).
- [ ] **Temperature-calibration** helper (Pb(NO3)2 / MeOH / etc.).
- [ ] **Reference manager** (ssNake): named references (Set / Save / Load /
  Apply) reused across spectra — SR presets per nucleus.

## Priority 2 — Processing depth (ssNake Tools parity)

- [ ] **Reference deconvolution** (ssNake) — divide out a reference lineshape to
  remove field inhomogeneity.
- [ ] **LPSVD** linear prediction (forward/backward) — first-point repair,
  truncated-FID extension (the MQMAS F1 truncation the user hit). (An
  autocorrelation LP op already exists; LPSVD is the upgrade.)
- [x] **Subtract averages / scale SW / scale car-ref (SR) / inverse FT /
  real·imag·conj** (ssNake Tools) as pipeline ops. *(done)*
- [x] **Processing history with per-step undo** (Process ▸ Processing steps):
  list the applied ops, remove any, re-apply the rest. *(done)*
- [ ] Apodizations: Hamming, Kaiser, shifted Gaussian (whole-echo), JMOD.
- [x] **Peak picking** (threshold + parabolic) → "add a line at every peak"
  (Decomposition ▸ Add a line at every peak). *(done)*

## Priority 3 — 2D depth (dmfit 2D menu)

- [x] **2D operations**: transpose, reverse F1/F2, extract diagonal, symmetrize,
  save F1/F2 projection to CSV, send projections/rows to fit (contour "2D ops ▾").
  *(done)* — [ ] still: save individual row/col, "make all 1Ds", diff-by-row.
- [ ] **DQ/SQ and Make-MQ**: build the double-quantum axis; sum/projection
  combinations (dmfit `Sum F2/F1`, `Proj`), 2D↔3D handling.
- [ ] **Referencing/shear conventions** per method (3Q/5Q ratios, Amoureux F1
  referencing) exposed and validated.
- [ ] Full-size interactive MQMAS fit (drag sites on the map; per-site F1/F2
  widths) — ties into the P0 per-dimension lineshape work.

## Priority 4 — Series, simulation, distribution

- [x] **Per-site relaxation** (Tools ▸ Per-site relaxation): decompose each slice
  on the current fit's lineshapes → T1/T2 per site. *(done)*
- [x] **Arrhenius/VFT** for VT series (Tools ▸ Variable temperature). *(done)* —
  [ ] mrinversion ILT (T1/T2 distributions, L-curve) still open.
- [ ] **SIMPSON bridge** (export spin system, run, overlay/fit) for REDOR/RFDR
  and sequences mrsimulator doesn't cover.
- [x] **Second-field prediction** ("what at X T?") from the current model
  (Decomposition ▸ Predict at another field). *(done)*

## Priority 5 — Packaging, validation, release

- [ ] PyInstaller one-folder + installer; auto-update check; crash reporter.
- [ ] Tutorials with bundled example data; Help opens them.
- [ ] **Validation corpus**: ≥10 published fits reproduced within uncertainties.
- [ ] CI (offscreen Qt, fast tests) on GitHub Actions; software paper draft.

---

## dmfit / ssNake / TopSpin parity — remaining items, tiered by likely use

Ranked by the probability they're needed for *this* work (solid-state NMR of
glasses/materials: quadrupolar nuclei, Czjzek, MQMAS, QCPMG, relaxation, figures
for papers/talks). The "unlikely for this scope" tier (DOSY, 3D data, SpinEv /
Shape Hole / Rector, Scroll F1 / Shuffle Imag Power, Print, File Manager) is left
off intentionally; revisit only if the work scope shifts.

### P6 — very likely, do first

- [x] **Integration tool + integral table** (Tools ▸ Integrals & measurements):
  drag regions → integral, %, centre, FWHM; Copy CSV. *(done)*
- [x] **Copy plot / export image** "with all lines" to clipboard and png/svg.
  *(done)*
- [x] **Export residual (Diff)** + per-component columns (export_text now writes
  ppm, experiment, model, residual, per-line). *(done)*
- [x] **Recent files** (File ▸ Open recent, last 12). *(done)*
- [ ] **Dual / compare display**: aligned side-by-side or overlaid two (or N)
  datasets for composition series (e.g. LAW3Cl0→4Ca, Na series), beyond the
  current overlay cockpit.

### P7 — likely, regular use

- [x] **FWHM & Centre-of-Mass readout**: in the Integrals & measurements table.
  *(done)*
- [x] **Computing-parameters dialog** (dmfit): Czjzek/MQMAS kernel resolution —
  computed points, Cq max, (Cq, η) step counts — via Decomposition ▸ Computing
  parameters. *(done)*
- [x] **MAS-spinning / sideband control**: the `sidebands` model (rate from the
  experiment params). *(done)*
- [x] **dmfit interop**: already via `.fxmla` read/write (the dmfit file format);
  a separate XML round-trip isn't needed.
- [ ] **Toggle Time ↔ Frequency** (inverse FT back to the FID): re-apodize /
  reprocess without reloading.
- [ ] **Real / Imag / Abs component views**: inspect the imaginary channel while
  phasing 2D and echoes.

### P8 — occasionally, when the case arises

- [x] **2D Symmetrize** (diagonal). *(done)* — [ ] diagonal-slope helper open.
- [ ] **Multi-curve / per-site relaxation UI**: the engine
  (`series.analyze_per_site`) exists — decompose each slice on fixed lineshapes
  for a per-site T1/T2; just needs the dialog.
- [ ] **Edge processing**: Check Eta (validity flag), Open Imaginary, Complex
  Conjugate (spectral reversal), nBandesMax (max-lines guard).
- [ ] **Quasar** model (a Czjzek-distribution variant) — low incremental value
  since `czjzek`/`ext_czjzek` already cover the amorphous quadrupolar case.

---

## Priority 9 (last, before "unlikely") — Multi-experiment correlation

Our own idea, generalizing the HMQC (1D − projection) difference. **Scheduled
last on purpose; the engine architecture may be built ahead but must NOT appear
in the app until everything above is done.**

- [ ] Given N datasets that **share a nucleus/dimension** — any mix of 1D, MQMAS,
  HMQC, REDOR — align their axes/projections and decompose features into
  **correlated vs un-correlated across arbitrary combinations**:
  - *1D + HMQC*: what correlates heteronuclearly; the rest is un-correlated.
  - *1D + MQMAS*: sites resolved by the isotropic axis vs lumped in the 1D.
  - *1D + HMQC + REDOR*: species that correlate **and** dipolar-dephase → assign,
    subtract to isolate the rest.
  - *1D + MQMAS + HMQC*: cross-check assignments across three experiments.
- [ ] Engine: an "experiment-set" model (each dataset declares its shared axis +
  a projection/observable), a scaling/alignment step, and set algebra
  (intersection = correlated, difference = specific-to-one) → each result a
  fittable workspace.
  *Accept: reproduce the HMQC 1D−projection result as the two-dataset special
  case, then a 1D+HMQC+REDOR three-way isolation on real data.*

---

## To think about (not scheduled)

- Multiplet/J-coupling in the indirect dimension; homonuclear recoupling shapes.
- A scripting/macro layer (ssNake macros) for batch processing.
- Cloud/collaborative recipes; a public recipe library for common materials.

### Standing engineering rules

1. Physics from mrsimulator/lmfit where possible — hand-rolled shapes only when
   they are genuinely analytic (Gauss/Lorentz/Voigt) and tested against theory.
2. A feature without a test does not merge; UI features get an offscreen test
   plus a real-display screenshot check (Windows dark mode).
3. Instrument folders are read-only; writes into them are refused in code.
4. Every fitted number is printed with its uncertainty, or with the reason there
   is none.
5. Manage resources: one reused widget set, lightweight snapshots, shared
   caches; free on close.
6. A feature merges with its **manual section** (Markdown in `larmor/help/`,
   organized by data type) and a Help button that opens it.
