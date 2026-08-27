# examples/

A small, real dataset that every tutorial and CLI example can run against:
`pCABS2-4/` holds three Bruker acquisitions of one Ca-aluminoborosilicate
glass (the pCABS2-4 sample, pure Ca, quenched at 2 GPa), trimmed to the files
LARMOR needs.

- `pCABS2-4/3616` — ²⁷Al MAS, single pulse, 26 kHz. The dmfit fit made for
  it (`pdata/1/1r.fxml`) is included, so the dmfit import path can be tried on
  a real file.
- `pCABS2-4/1118` — ¹¹B MAS, single pulse, 20 kHz, with its dmfit fit.
- `pCABS2-4/3620` — ²⁷Al 3QMAS (mp3qdfsz), the 2D dataset used by the MQMAS
  workflow.

`pCABS2-4_27Al.recipe.json` and `pCABS2-4_11B.recipe.json` are LARMOR fits of
the two 1D spectra, started from the dmfit models. `pCABS2-4_fits.figure.json`
is the figure spec that renders both deconvolutions side by side
(`pCABS2-4_fits.png`, the figure shown in the README). Recipes reference their
source data by a path relative to the repository root and a SHA-256 hash, so
run commands from the repository root:

```
larmor info examples/pCABS2-4/3616
larmor fit examples/pCABS2-4_27Al.recipe.json --window 150 -80 --plot fit.png
```

LARMOR never modifies the acquired files themselves; fits saved next to
them are ordinary new files.
