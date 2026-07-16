# Tutorial 3 — The figure studio: publication figures from any experiment

*Time: ~20 minutes. Works with everything LARMOR can read: 1D spectra (zg,
echo), fits, 2D maps (MQMAS, HMQC, SQ-DQ), and series (saturation recovery,
REDOR).*

A LARMOR figure is a small JSON document — a **spec** — that describes what to
draw. Like the fit recipe, it references data by path and can be saved and
re-rendered identically later: when a reviewer asks for the same figure with
different axis limits eight months after submission, you edit two numbers and
re-export. No mouse archaeology.

## 1. Fastest path: templates in the app

```
larmor app
```

Load your data, then in the **Figure studio** panel click *Templates for
loaded data*. LARMOR inspects the source and offers what applies: a plain 1D
for a spectrum, `2d` if the EXPNO has an indirect dimension, `satrec`/`redor`
if TopSpin analysis files (`t1ints.txt`, `redor.txt`) are present. Click one,
edit the JSON, **Preview**, then **Export** (png + svg + pdf at 600 dpi, ready
for the journal). Exports into instrument data folders are refused — figures
go in your own folders.

## 2. Anatomy of a 1D spec

```json
{
  "kind": "1d",
  "style": "article-wide",
  "xlim": [150, -80],
  "traces": [
    {"path": "C:/data/CaAlGlass.fxmla", "label": "experiment", "color": "black"},
    {"recipe": "CaAlGlass.recipe.json", "part": "total", "label": "fit", "color": "crimson"},
    {"recipe": "CaAlGlass.recipe.json", "part": "site", "site": 0, "label": "AlO$_4$", "linestyle": "--"},
    {"recipe": "CaAlGlass.recipe.json", "part": "residual", "label": "residual", "offset": -260}
  ]
}
```

- **style** — `article` (single column), `article-wide`, `presentation`,
  `thesis`. Each is a complete font/linewidth/size bundle; switching styles
  regenerates the same figure for a different medium.
- **traces** pull from a spectrum file (`path`), a saved fit (`recipe` +
  `part`: `total`, `site`, `residual`), or inline arrays (`data`).
- Every trace takes `scale`, `offset` (for stacked comparisons and offset
  residuals), and `normalize` (peak-normalize inside a ppm window — for
  comparing spectra acquired with different gains, exactly like NMRVEW's
  `norm_0_to_1`).
- Labels accept LaTeX: `"AlO$_4$"`, `"$^{27}$Al"`.

## 3. 2D maps (MQMAS, HMQC, SQ-DQ, ...)

```json
{
  "kind": "2d",
  "style": "thesis",
  "path": "C:/data/my2Dexpno",
  "xlabel": "$^{27}$Al NMR shift (ppm)",
  "ylabel": "Isotropic dimension (ppm)",
  "levels": {"mode": "log", "n": 12},
  "proj_top": true, "proj_right": true,
  "overlay_top": [
    {"path": "C:/data/my1Dexpno", "label": "zg", "normalize": true}
  ],
  "subproj": [
    {"f1": [55, 90], "label": "AlO$_4$", "scale": 1.0},
    {"f1": [20, 55], "label": "AlO$_5$", "scale": 0.5}
  ],
  "slopes": [{"slope": 1.0, "intercept": 0.0}],
  "annotation": "my glass, 3QMAS"
}
```

The moves that make a 2D figure publishable, all declarative:

- **Contour floor from the noise.** If you don't set `levels.min_frac`, LARMOR
  measures the noise in the matrix edges and puts the lowest contour at ~8σ —
  no more contouring the noise floor. `"mode": "log"` spaces the levels
  logarithmically (both weak and strong features visible).
- **`overlay_top`** draws any external 1D (e.g. the quantitative zg) over the
  2D's top projection — the classic "is my MQMAS projection representative?"
  panel.
- **`subproj`** integrates an F1 band and plots it as an extra top trace —
  per-site projections.
- **`slopes`** draws diagonal/CS/QIS reference lines. `"negative": true` adds
  dashed red negative contours (phase-sensitive experiments, HMQC artifacts).

## 4. Series: saturation recovery and REDOR

```json
{"kind": "series", "mode": "satrec", "path": "C:/data/expno1901", "stretched": true}
```

Reads TopSpin's own `t1ints.txt`, plots the integrals on a log time axis, fits
a (stretched) exponential recovery, and annotates T₁ ± error (and β). If the
error is huge, believe it — it usually means the longest delay didn't reach
the plateau.

```json
{"kind": "series", "mode": "redor", "path": "C:/data/expno1901"}
```

Reads `redor.txt` (S₀/S pairs), converts to ΔS/S₀ against recoupling time from
the spinning speed in the file. Inline `data` with `x`/`y`/`yerr` works too,
for curves computed elsewhere (e.g. a SIMPSON simulation to overlay — that
hook is deliberate).

## 5. From Python, for batch work

```python
from larmor import figures
spec = {...}                                  # same JSON as the app
figures.export(spec, "fig/CaAlGlass_fit")     # writes .png, .svg, .pdf
```

Because specs are plain dicts, a composition series is a for-loop that edits
one path per iteration. Save each spec next to its figure — that pair is the
reproducibility unit.
