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

## Priority 0 — Lineshape rework  ⭐ (the user's headline)

The model catalogue works but is not yet at dmfit breadth, and the width /
amplitude conventions deserve a cleanup. This is the next big block.

**Architecture.**
- [ ] **Width entry in Hz *or* ppm per line** (dmfit toggles this): store the
  canonical value, let the cell accept `300Hz` / `2ppm` and display either.
- [ ] **Amplitude = peak or area** toggle (dmfit "GL Norm" is the
  area-normalized twin): quantification should be able to report either
  consistently; expose an `amp_mode` on each line.
- [ ] **Per-dimension lineshape for 2D** (dmfit "GL/F1"): independent F1/F2
  widths and shapes on a fitted MQMAS site.
- [ ] A **line base** so a site declares `lineshape ∈ {gauss, lorentz, voigt,
  pseudo-voigt}` instead of separate models where it's just a shape choice.

**New analytic / physical models (dmfit parity, ranked by usefulness).**
- [x] **Voigt (true)** — Gaussian ⊗ Lorentzian, independent widths. *(done
  2026-07-21)*
- [ ] **GL Norm** — area-normalized Gauss/Lorentz for clean quantification.
- [ ] **Spinning sidebands (generic "ss band")** — a centre + Herzfeld–Berger
  sideband manifold at arbitrary intensity ratios, independent of CSA (dmfit's
  manual sideband lines), for cases csa_mas doesn't cover.
- [ ] **J-multiplet** (dmfit "Jmultiplet") — n equivalent couplings → binomial
  multiplet with a J (Hz) and a per-component lineshape; and the **J-dispersion**
  variants (residual dipolar / distribution of J).
- [ ] **Function fit** (ssNake) — a user-typed `y(x; a, b, …)` expression with
  free parameters, safely evaluated, for the odd empirical case.
- [ ] **Shape-from-file / user shape** — generalize the current `spectrum`
  component into a first-class shifted+scaled+broadened basis line (dmfit "W
  Shape file" / "user shape").

**Disorder / distribution physics (beyond dmfit).**
- [ ] **CSA Czjzek / Gaussian-isotropic** distribution (disordered shielding),
  the CSA analogue of the quad Czjzek.
- [ ] **Correlated δiso–Cq distribution** for glasses (a 2D Gaussian on the
  (δiso, Cq) plane) — closer to real amorphous lineshapes.
- [ ] **Extended Czjzek in 1D** already exists; expose ρ/ε clearly and validate
  against Le Caër's figures.

*Accept for P0: reproduce a published glass fit (e.g. CaAlGlass variants) with
Voigt + GL-Norm + sidebands, and a J-multiplet crystal, within uncertainties.*

---

## Priority 1 — Utilities & tools (ssNake Utilities menu)

Small, self-contained, high daily value.
- [ ] **NMR table** — interactive periodic table of Larmor frequencies at a
  settable B0 (or set the ¹H frequency of *your* magnet and read every nucleus).
  Double-click an element → isotopes with spin, natural abundance, γ, Q,
  receptivity, reference compound. (ssNake `nmrTable`.)
- [ ] **Chemical-shift / Cq / dipolar-distance / MQMAS-parameter** conversion
  tools (ssNake): the standard algebra (δiso↔ν, Cq↔νQ↔PQ, r↔D, MQMAS δ1/δ2 → δiso,
  PQ) as dialogs, with copy-to-clipboard.
- [ ] **Temperature-calibration** helper (Pb(NO3)2 / MeOH / etc.).
- [ ] **Reference manager** (ssNake): named references (Set / Save / Load /
  Apply) reused across spectra — SR presets per nucleus.

## Priority 2 — Processing depth (ssNake Tools parity)

- [ ] **Reference deconvolution** (ssNake) — divide out a reference lineshape to
  remove field inhomogeneity.
- [ ] **LPSVD** linear prediction (forward/backward) — first-point repair,
  truncated-FID extension (the MQMAS F1 truncation the user hit).
- [ ] **Subtract averages / offset correction / scale SW / scale car-ref**
  (ssNake Tools) as one-click ops.
- [ ] **Processing history with per-step undo** (ssNake pipeline) instead of
  reset-to-original only; the pipeline already lives in the recipe.
- [ ] Apodizations: Hamming, Kaiser, shifted Gaussian (whole-echo), JMOD.
- [ ] Peak picking (threshold + parabolic) → "add a line at every peak".

## Priority 3 — 2D depth (dmfit 2D menu)

- [ ] **2D operations**: transpose, reverse F1/F2, save row/col/projection,
  extract diagonal, "make all 1Ds", diff-by-row, add-sidebands (dmfit 2D
  right-click).
- [ ] **DQ/SQ and Make-MQ**: build the double-quantum axis; sum/projection
  combinations (dmfit `Sum F2/F1`, `Proj`), 2D↔3D handling.
- [ ] **Referencing/shear conventions** per method (3Q/5Q ratios, Amoureux F1
  referencing) exposed and validated.
- [ ] Full-size interactive MQMAS fit (drag sites on the map; per-site F1/F2
  widths) — ties into the P0 per-dimension lineshape work.

## Priority 4 — Series, simulation, distribution

- [ ] Per-site relaxation (decompose each slice on fixed lineshapes → T1/T2 per
  site, not per window).
- [ ] Arrhenius/VFT for VT series; mrinversion ILT (T1/T2 distributions,
  L-curve).
- [ ] **SIMPSON bridge** (export spin system, run, overlay/fit) for REDOR/RFDR
  and sequences mrsimulator doesn't cover.
- [ ] Second-field prediction ("what at X T?") from the current model.

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

- [ ] **Integration tool + integral table** (ssNake Integrals / TopSpin
  integration): drag regions → integrals with errors, exportable — quantify site
  fractions without a full fit.
- [ ] **Copy plot / export image** (dmfit copy-to-clipboard): as-presented,
  spec-only, and **"with all lines"** (each component) to clipboard and to
  file (png/svg) — figures for papers and slides.
- [ ] **Export residual (Diff)** and **per-component export** ("with all lines"):
  data + total + each line as columns / traces for publication figures.
- [ ] **Recent files / "Open Last fit"**: reopen the last datasets/recipes in a
  click.
- [ ] **Dual / compare display**: aligned side-by-side or overlaid two (or N)
  datasets for composition series (e.g. LAW3Cl0→4Ca, Na series), beyond the
  current overlay cockpit.

### P7 — likely, regular use

- [ ] **FWHM & Centre-of-Mass readout** (ssNake): measure linewidth / centroid
  over a region without fitting.
- [ ] **Computing-parameters dialog** (dmfit): expose the Czjzek/MQMAS kernel
  resolution — computed size (2ⁿ), (Cq, η) step counts, sweep/Gauss multipliers,
  ssb max, distribution threshold — accuracy vs speed.
- [ ] **MAS-spinning / explicit sideband control**: set the rate and sideband
  order for CSA/quad sideband manifolds (slow-MAS cases).
- [ ] **dmfit-XML export/import**: round-trip fits with dmfit for cross-checking
  against the reference program.
- [ ] **Toggle Time ↔ Frequency** (inverse FT back to the FID): re-apodize /
  reprocess without reloading.
- [ ] **Real / Imag / Abs component views**: inspect the imaginary channel while
  phasing 2D and echoes.

### P8 — occasionally, when the case arises

- [ ] **MQMAS post-processing**: Diagonal-slope / shear helpers and 2D
  Symmetrize.
- [ ] **Multi-curve / per-site relaxation**: fit several sites' decays together
  (beyond the guided single-site T1/T2).
- [ ] **Edge processing**: Check Eta (validity flag), Open Imaginary, Complex
  Conjugate (spectral reversal), nBandesMax (max-lines guard).
- [ ] **Quasar** model (a Czjzek-distribution variant) — low incremental value
  since `czjzek`/`ext_czjzek` already cover the amorphous quadrupolar case.

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
