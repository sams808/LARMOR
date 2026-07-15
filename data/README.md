# data/

This folder holds **references to** raw data, never copies of it.

Raw Bruker TopSpin experiments live under `Desktop\WSU_work\NMR\NMRFAM\DATA\...` and legacy dmfit fits live wherever they were saved (e.g. `Desktop\CaAlGlass.fxmla`). LARMOR notebooks and code read those paths directly and read-only — nothing under `data/` in this repo should ever be the instrument data itself, only small text notes (source paths, content hashes, what an EXPNO contains) that are safe to commit.
