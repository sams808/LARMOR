# LARMOR

A modern, open successor to [dmfit](https://nmr.cemhti.cnrs-orleans.fr/Dmfit/) for solid-state NMR lineshape fitting and spin-dynamics simulation — built to be usable by students with no NMR background and powerful enough for advanced research work.

Full architecture and workflow design: see the published design document (artifact link shared in the project conversation history). Summary of the plan:

- **Reuse, don't rebuild the physics.** The numerical core comes from [mrsimulator](https://mrsimulator.readthedocs.io), [lmfit](https://lmfit.github.io/lmfit-py/), [mrinversion](https://mrinversion.readthedocs.io), and [csdmpy](https://csdmpy.readthedocs.io) (P. Grandinetti's group, Ohio State). New work here is ingestion, UX, and orchestration.
- **Two-tier UX.** Guided mode for students, Expert mode (with an embedded Python console) for research use — same underlying project file.
- **Bruker/TopSpin import via [nmrglue](https://nmrglue.readthedocs.io)**, strictly read-only against instrument data.
- **Legacy dmfit `.fxmla` import** for migrating existing fits.
- **Advanced simulation as opt-in plug-ins**: SIMPSON (exact density-matrix, e.g. REDOR) and DFT/MD tensor import — not part of the default path.

## Status

**Phase 0 — feasibility: PASSED** (see `notebooks/phase0_feasibility.ipynb`, executed with outputs):

1. `CaAlGlass.fxmla` (dmfit Czjzek fit of a Ca-Al glass, 27Al) reproduced in mrsimulator with
   normalized RMSD 0.027 against the experimental spectrum embedded in the fxmla itself.
   **Key finding:** dmfit's `sCZ_CQ` = 2 × mrsimulator's `CzjzekDistribution` sigma — the
   conversion constant the future `.fxmla` importer must apply to every `CzSimple` line.
2. Bruker EXPNO 1903 (19F Hahn echo) opened with nmrglue end to end; an mtime+size snapshot
   of the whole EXPNO folder proves the read is strictly non-destructive. Lesson for the
   importer: `acqus` MASR (4200 Hz) disagrees with the operator-typed title (35.714 kHz) —
   surface both, trust neither silently.

Next: **Phase 1 — core workbench** (Guided/Expert-mode app over 1D MAS/static spectra, Bruker +
dmfit import, lmfit uncertainty).

The MQMAS twin (`CaAlGlassMQ.fxmla`) is deferred to Phase 2 with the other 2D methods.

## Data policy

Raw instrument data (TopSpin EXPNO folders) and legacy `.fxmla` files are **never copied into this repo** and **never written to**. `data/` holds only small reference notes (paths, hashes) — see `data/README.md`.

## Environment

```
conda env create -f environment.yml
conda activate larmor
```
