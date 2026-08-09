# LARMOR

A modern, open successor to [dmfit](https://nmr.cemhti.cnrs-orleans.fr/Dmfit/) for solid-state
NMR lineshape fitting and analysis — a native desktop application (PySide6 + pyqtgraph) built on
the [mrsimulator](https://mrsimulator.readthedocs.io) / [lmfit](https://lmfit.github.io/lmfit-py/)
/ [csdmpy](https://csdmpy.readthedocs.io) / [nmrglue](https://nmrglue.readthedocs.io) stack, with
one thing dmfit doesn't give you: **an uncertainty on every fitted number**, and fully
reproducible fits.

## Install (quick start)

**Full step-by-step instructions — Windows, macOS, Linux — are in
[INSTALL.md](INSTALL.md)**, including troubleshooting for the usual snags.

- **Windows, one-click:** double-click **`install.bat`** in this folder. It sets
  everything up and puts a **`LARMOR.bat` on your Desktop** to launch the app.
  Update later by double-clicking **`update.bat`**.

- **Manual (any platform),** from inside the repository folder:
  ```
  # Recommended (Conda — handles the compiled packages for you):
  conda env create -f environment.yml
  conda activate larmor

  # …or with pip, into a Python 3.11 virtual environment:
  pip install -r requirements.txt
  ```
  **Launch:** `larmor desktop`.

> Use **Python 3.11** (3.10–3.12 fine, **not 3.13** yet). If `larmor desktop`
> says a package is missing or `mrsimulator` won't install, see
> [INSTALL.md → Troubleshooting](INSTALL.md#troubleshooting).

## Capabilities

**Fitting** — dmfit-style paddles (drag position + amplitude, side handles for width), a
spreadsheet parameter table with pin-to-fix, live re-simulation, and:
- **15 lineshape models** in a self-describing registry (`larmor/models/`) — Gauss/Lorentz,
  Voigt, J-multiplet, sidebands, Czjzek, extended Czjzek, dmfit-style "Amorphous"
  (Gaussian Cq/η disorder), discrete 2nd-order quadrupolar CT, 1st-order quadrupolar
  (satellites), quad+CSA, CSA powder, external spectrum/function backgrounds — every
  quadrupolar/CSA model resolves gyromagnetic ratio and spin generically, so a new nucleus
  needs no code changes. Fast cached (Cq, η) kernel where applicable — a full multi-site fit
  runs in seconds, not minutes.
- Constraints: fix, bounds, algebraic links; **dependent positions in ppm _or_ Hz** and
  amplitude/width ratios via dialogs (no expression writing); full error propagation.
- Fit **zones** (dmfit-style union of regions), editable νrot / Larmor / nucleus, quantification
  table (% ± error, `larmor.quantify`), CSV export.
- **Auto Fit** (multi-start, escapes local minima).

**Uncertainty, always** — this is LARMOR's actual headline difference from dmfit:
- **Covariance**, **Monte-Carlo** (parametric bootstrap), and **χ² profile** (real 1σ/2σ
  confidence intervals, not just the covariance) error methods everywhere a fit happens —
  single spectrum, batch, or sequential. Monte-Carlo and χ² profile run **across all your CPU
  cores** (one left free for the app), not one refit at a time.
- **At-bounds diagnosis**: a parameter that finishes pinned at a bound is flagged, with
  uncertainties computed conditional on that pin instead of silently vanishing.
- **Identifiability**: unidentifiable parameter pairs (from the fit covariance) are flagged with
  a warning, not left to be discovered later.

**Batch & multi-spectrum workflows**
- **Batch fit**: one shared model fit to many spectra at once, amplitudes free per spectrum,
  optional per-parameter "release" (let δiso/width drift a little, per spectrum); per-cell
  **exclude a component** for samples where a line doesn't apply; exports a long-format CSV
  (value + error + integrated population %) with each spectrum's saved fit alongside it.
- **Sequential fit**: forward/backward warm-started fitting across a series (e.g. a pressure or
  composition series) — each spectrum starts from its converged neighbour, 1/2/4/8/16 passes.
- **Multi-field / multi-dataset** simultaneous fits (`larmor multifit`) — lifts the Cq/δiso
  degeneracy a single field can't resolve.

**2D / MQMAS** (`larmor.twod`) — interactive fitting (click to place a site, fitted model
overlaid as a dashed contour), exact hypercomplex 2D phasing, manual shear, contour figures with
noise-measured levels and projections.

**Processing** (TopSpin/ssNake parity) — EM/GM/SINE/QSINE/TRAF windows, TDeff, ZF (factor or SI),
FCOR, FT, manual/ACME phase, SR, magnitude, Hilbert reconstruction, linear prediction
(forward/backward), whole-echo, polynomial/iterative/interactive-anchor baselines, region
extract, spectra algebra, align, peak picking, 2-point background subtraction. Pipelines are
stored in the recipe and **replayed on load**, including inside every batch/publication figure.

**Relaxation & recoupling** — automatic **T1/T2** from arrayed EXPNOs (satrec/invrec/CPMG/T1ρ,
window- or **per-site** via NNLS decomposition), **REDOR** dipolar couplings and distances
(model-free short-time or full pair curve), **QCPMG** (echo-train sum, spikelet or sum-echo
absorption spectra), two-field infinite-field extrapolation.

**Import** (ssNake-style, universal) — point at almost anything: a legacy dmfit `.fxmla`/`.fxml`,
a LARMOR recipe, or **any Bruker path** — a processed `1r`/`2rr` file, a raw `fid`/`ser`, a
`pdata/N` folder, or an EXPNO folder — plus Varian/Agilent `.fid`. The reader figures out 1D vs
2D, raw vs processed, and a real spectroscopic 2D vs a pseudo-2D arrayed experiment. **Open
FID…** processes the raw fid/ser before the Fourier transform (windowing, zero-fill, phase, and
for 2D the indirect quadrature mode). Everything is strictly read-only against instrument
folders.

**Publication figures** — the Plotting Studio: 1D overlay/stack, 2D contour (with a fitted-model
overlay and computed CS/QIS reference lines), deconvolution/composition grids from a whole batch
fit, species-distribution bars, and relaxation/REDOR series, all from a declarative JSON spec
(savable, re-renderable identically later). Per-component colors and legend control, journal
style presets (`article`, `article-wide`, `presentation`, `thesis`), png + svg + pdf export.

**Advanced** — **DFT tensor import** (CASTEP/QE `.magres` → fittable sites, EFG→Cq validated
against first principles); **SIMPSON** bridge for exact density-matrix recoupling simulations.

Reuse-first design: the physics comes from mrsimulator + lmfit; LARMOR adds ingestion, the
dmfit-faithful UX, orchestration, uncertainties, and reproducibility. Instrument folders are
always read-only.

## LARMOR vs dmfit

dmfit is the tool this project is modeled on and still cites throughout (docs, Methods text,
`.fxml`/`.fxmla` import). What's different:

| | dmfit | LARMOR |
|---|---|---|
| Uncertainty on a fitted value | Optional, separate Monte-Carlo tool | Covariance/MC/χ² profile everywhere, parallelized across CPU cores |
| Fit file format | Binary/XML, dmfit-only | JSON recipe, diffable, forward-compatible, data referenced by path+hash never embedded |
| Batch / series fitting | Manual, one spectrum at a time | Built-in shared-model batch fit + sequential warm-start fit |
| Scripting / automation | None | Python library + `larmor` CLI (`fit`/`batchfit`/`seqfit`/`multifit`/…) |
| Publication figures | External plotting | Built-in Plotting Studio, spec-driven, reproducible |
| Platform | Windows-native (Delphi) | Cross-platform Python/Qt |

## Trust & validation

`docs/LARMOR_VALIDATION_REPORT.md` cross-checks LARMOR's physics against analytic theory,
direct ensemble simulation, and **real published data** — including reproducing 20 fits from a
peer-reviewed paper (Soudani *et al.*, *J. Non-Cryst. Solids* **638**, 123085 (2024)) with exact
parameter agreement. `docs/ROADMAP_v1.md` tracks, honestly, what's still needed before a 1.0
sign-off (splitting the 4000+-line main window module and setting up CI are the biggest open
items — neither affects the correctness of a fit today).

## Status

LARMOR's native desktop app is the primary, actively-developed interface and has been for most
of this project's history — this is not an early prototype. The `larmor` Python package
(`pip install -e .`) is a fully working library on its own (see `larmor.engine`/`larmor.fit` for
the fitting core, `larmor.io` for readers, `larmor.figures` for the publication-plotting engine),
and every capability above is backed by an automated test (500+ tests, `pytest tests/`).

A browser-based variant (`larmor app`, FastAPI + Plotly) also exists in the codebase but is
**not actively maintained** — the desktop app absorbed all feature development early on. Treat
it as experimental; if you need a browser-based/shared-lab-server interface, open an issue
rather than assuming it's current.

**Tutorials** (in [docs/tutorials/](docs/tutorials/)) currently cover the fitting basics; a
refresh covering batch fitting, MQMAS, and error analysis is planned (`docs/ROADMAP_v1.md`
Pillar E).

## Data policy

Raw instrument data (TopSpin EXPNO folders) and legacy `.fxmla` files are **never copied into this repo** and **never written to**. `data/` holds only small reference notes (paths, hashes) — see `data/README.md`.

## Installation & launching

**See [INSTALL.md](INSTALL.md) for the complete, cross-platform guide with
troubleshooting.** In brief, from inside this folder:

```
conda env create -f environment.yml     # recommended (or: pip install -r requirements.txt)
conda activate larmor
larmor desktop                           # …or double-click LARMOR.bat on Windows
```

- The **native desktop app** (PySide6 + pyqtgraph — instant zoom/pan/drag, no
  browser) is the primary interface.
- CLI without any GUI: `larmor info <path>`, `larmor import <fxmla>`, `larmor fit <recipe>`,
  plus `satrec`/`redor`/`magres`/`multifit`/`batchfit`/`seqfit` for their respective workflows.
- A packaged Windows build (`packaging/larmor.spec`, PyInstaller) exists for a Conda-free
  install, but is not yet part of a signed, tested release pipeline — the Conda/pip route above
  is the reliable path today.
- If an existing env predates a feature, refresh it with
  `conda env update -f environment.yml`.
