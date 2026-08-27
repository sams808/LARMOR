# Roadmap

Where LARMOR stands and where it is going, checked against the code. There are
no dates on any of this; items move up when they block real work.

## Where it stands

The fitting engine is validated on real published work: 20 dmfit fits from a peer-reviewed ¹¹B pressure series (Soudani et al. 2024) import with zero warnings and match the paper's published Table 2 parameters to its own rounding. Every fitted number carries an uncertainty from one of three estimators (least-squares covariance, Monte Carlo, or a χ² profile); the heavier two run across a process pool, and the errors propagate into site populations, the batch CSV columns, and the error bars on series plots. The newer batch and sequential-fit paths have not yet had a cross-check of their own against a published dataset.

File formats are versioned. A recipe writes a `larmor_recipe_version` field, and when a file written by a newer version is opened, unrecognized fields are dropped with a note rather than a refusal to load. A project bundle (`.larproj.json`) saves every open 1D workspace, its processed spectrum, and its overlays in one file. 2D workspaces, figures, and live batch-fit sessions are not part of a bundle yet.

The weak points are engineering and onboarding. The main window lives in a single module of about 4,600 lines. There is no continuous integration; the test suite is substantial but runs only when started by hand, in one development environment. Releases so far are git tags: no installer executable has been built, verified end to end, and published, and nothing is signed. Development and testing happen on Windows only. The three tutorials predate the desktop app and cover none of the batch, MQMAS, or error-analysis workflows, and the repository bundles a single example (a fitted recipe and figure for a ²⁷Al glass — no spectrum data ships with it).

## Next

- Split the main-window module and put continuous integration in place, so the test suite runs on every push in a clean environment.
- Build the Windows installer from the existing PyInstaller spec, verify it end to end on a machine without a development setup, and publish it as a release.
- Extend project bundles to cover 2D workspaces, figures, and batch-fit sessions.
- Bring the documentation up to the app: tutorials for batch fitting, MQMAS, and error analysis, and example datasets beyond the single glass — a 2D set, a composition series, a relaxation set.
- A peak-pick assignment table: pickable, labelable peaks (AlIV / AlV / AlVI and the like) that export and double as fit starting points.

## Further out

Lineshapes and physics. Per-dimension lineshapes for 2D fits, so an MQMAS site can carry independent F1 and F2 widths and shapes. A correlated δiso–Cq distribution (a 2D Gaussian on the (δiso, Cq) plane) for lineshapes closer to real amorphous materials.

Processing depth. LPSVD linear prediction — an autoregressive LP op already exists; the SVD variant is the upgrade for truncated FIDs. Reference deconvolution to divide out field inhomogeneity. A time↔frequency toggle that returns to the FID for re-apodization without reloading, and real/imaginary component views for inspecting the imaginary channel while phasing.

2D. DQ/SQ handling: building the double-quantum axis and the associated sum and projection combinations.

Series and batch. Inverse Laplace transforms of relaxation series (T1/T2 distributions via mrinversion). A dual/compare display that aligns two or more datasets for composition series, beyond the current overlays. A two-way batch table where selecting a row highlights the corresponding spectrum. An explicit monotonic-series penalty for sequential fits; the current between-pass trajectory smoothing regularizes toward smooth series but does not enforce monotonicity. A fit-history timeline, and a watch folder that queues new spectra for automatic fitting.

Interaction. TopSpin's drag-to-phase gesture (horizontal drag for zero order, vertical for first order, about the existing pivot line). Keyboard-driven fitting: shortcuts to nudge the selected parameter, cycle sites, and start a fit. Automatic spinning-sideband detection that offers to add a linked sideband manifold in one click. An interactive dmfit amplitude calibration — today the import and export amplitude factors are fixed constants calibrated on one known fit, and pasting a single dmfit amplitude for a known line would pin the conversion exactly.

Utilities and platform. A temperature-calibration helper (Pb(NO₃)₂, methanol, and similar standards). A macOS build. Last of all, a multi-experiment correlation engine that aligns any set of experiments sharing a nucleus (1D, MQMAS, HMQC, REDOR) and decomposes features into correlated and un-correlated parts; the core of it exists in `larmor/correlate.py` with tests, but wiring it into the app comes after everything above.

## Left manual on purpose

A few things that look automatable stay with the scientist, because they are physical judgments. Soft priors or penalty terms from a known crystal structure depend on how far that structure can be trusted for the sample at hand. Whether an extra site is justified depends on the system and on the question being asked, so no AIC/BIC score or F-test decides it. What counts as structured residual likewise depends on the model: the app prints a note when the residual exceeds the noise level and flags parameter pairs with |r| ≥ 0.95 as unidentifiable, but it does not grade the fit.

## Not planned

GPU acceleration was evaluated in detail and CPU parallelism was chosen instead; the expensive error estimators already run across a process pool.

The browser variant (`larmor/app.py`) is unmaintained. It has no test coverage and has not tracked the current recipe and model schema; all development goes into the desktop app.
