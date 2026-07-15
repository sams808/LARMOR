# LARMOR

A modern, open successor to [dmfit](https://nmr.cemhti.cnrs-orleans.fr/Dmfit/) for solid-state NMR lineshape fitting and spin-dynamics simulation — built to be usable by students with no NMR background and powerful enough for advanced research work.

Full architecture and workflow design: see the published design document (artifact link shared in the project conversation history). Summary of the plan:

- **Reuse, don't rebuild the physics.** The numerical core comes from [mrsimulator](https://mrsimulator.readthedocs.io), [lmfit](https://lmfit.github.io/lmfit-py/), [mrinversion](https://mrinversion.readthedocs.io), and [csdmpy](https://csdmpy.readthedocs.io) (P. Grandinetti's group, Ohio State). New work here is ingestion, UX, and orchestration.
- **Two-tier UX.** Guided mode for students, Expert mode (with an embedded Python console) for research use — same underlying project file.
- **Bruker/TopSpin import via [nmrglue](https://nmrglue.readthedocs.io)**, strictly read-only against instrument data.
- **Legacy dmfit `.fxmla` import** for migrating existing fits.
- **Advanced simulation as opt-in plug-ins**: SIMPSON (exact density-matrix, e.g. REDOR) and DFT/MD tensor import — not part of the default path.

## Status

**Phase 0 — feasibility.** Confirming the reuse thesis before building any UI:

1. Reproduce `CaAlGlass.fxmla` and `CaAlGlassMQ.fxmla` (existing dmfit Czjzek fits of a Ca-Al glass) in mrsimulator + lmfit.
2. Open one real Bruker TopSpin EXPNO from the `NMRFAM/DATA/2026-07` dataset with nmrglue, end to end, read-only.

See `notebooks/phase0_feasibility.ipynb`.

## Data policy

Raw instrument data (TopSpin EXPNO folders) and legacy `.fxmla` files are **never copied into this repo** and **never written to**. `data/` holds only small reference notes (paths, hashes) — see `data/README.md`.

## Environment

```
conda env create -f environment.yml
conda activate larmor
```
