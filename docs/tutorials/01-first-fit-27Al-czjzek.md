# Tutorial 1 — Your first fit: ²⁷Al MAS spectrum of a glass

*Time: ~15 minutes. You need: the `larmor` conda environment, and a dmfit
`.fxmla` file with an embedded spectrum (this tutorial uses `CaAlGlass.fxmla`,
a ²⁷Al MAS spectrum of a Ca-aluminosilicate glass fitted with Czjzek
distributions).*

In this tutorial you will load a spectrum, inspect the model, run a fit, and —
the part dmfit never gave you — read the uncertainties.

## 1. Look before you fit

Open a terminal, activate the environment, and ask LARMOR what's in the file:

```
conda activate larmor
larmor info C:\path\to\CaAlGlass.fxmla
```

You should see something like:

```
dmfit fit file (version 20110208, mode 'Fit 1D')
dimension F2: 27Al at 195.483 MHz, 9 lines
  [0] CzSimple   pos=66.18 ppm  sCZ_CQ=4549 kHz
  [1] CzSimple   pos=37.20 ppm  sCZ_CQ=4407 kHz
  [2] Gaus/Lor   pos=236.27 ppm  wid=20.25 ppm
  ...
embedded spectrum: 1D, NP=8192
```

Things to notice:

- Three `CzSimple` lines — dmfit's simple Czjzek distribution, the standard
  model for a distribution of quadrupolar couplings in a disordered material.
  For ²⁷Al in a glass, the three positions (~66, ~37, ~14 ppm) are the classic
  AlO₄ / AlO₅ / AlO₆ coordination assignment.
- Two `Gaus/Lor` lines at 236 and 208 ppm. These are **not** extra aluminum
  sites — at 33.3 kHz MAS and 195.5 MHz, one rotor frequency is 170 ppm, so
  236 ≈ 66 + 170 and 208 ≈ 37 + 170: they are the +1 spinning sidebands, which
  dmfit models as separate ad-hoc lines. LARMOR simulates sidebands physically,
  so these lines get frozen automatically during the fit.

## 2. Convert to a LARMOR recipe

```
larmor import C:\path\to\CaAlGlass.fxmla -o CaAlGlass.recipe.json
```

The recipe is a small, readable JSON file: the model *only*, with the data
referenced by path and SHA-256 hash instead of copied. Open it in any text
editor. Every parameter looks like:

```json
"sigma_Cq_MHz": { "value": 2.274, "stderr": null, "vary": true, "min": 0.05, "max": null, "expr": null }
```

The Czjzek width was converted with `sigma = sCZ_CQ / 2` — dmfit and
mrsimulator use conventions that differ by exactly a factor of two (LARMOR
established this against this very file).

## 3. Fit

```
larmor fit CaAlGlass.recipe.json --window 150 -80 --plot fit.png
```

The window `150 -80` (high ppm, low ppm) covers the central-transition region
and excludes the sidebands. The first run builds the Czjzek simulation kernel
(~15 s); after that, iterations are milliseconds.

You should see, at the end of the report:

```
normalized RMSD: 0.0025
```

and `fit.png` shows the experiment (black), total fit (red), and each site's
contribution (dashed).

## 4. Read the uncertainties — the important part

Open the updated `CaAlGlass.recipe.json` and look at each site's `stderr`
fields, or read them from the fit report. On this dataset:

- The AlO₄ site is well determined: δiso = 65.1 ± 0.4 ppm, σ(Cq) = 1.84 ± 0.13 MHz.
- The middle (AlO₅) site is **not**: its σ(Cq) error is several times its
  value. The 1D lineshape simply does not contain enough information to pin
  three overlapping Czjzek sites independently.

That second bullet is not a failure of the fit — it is the fit telling you a
scientific truth: to quantify the AlO₅ site you need more data (an MQMAS
spectrum, a second field, or a constraint from chemistry). dmfit reports the
same parameter values with no error bars, so this degeneracy stays invisible.

Tutorial 2 shows how to add exactly those constraints.

## 5. Same thing, interactively

```
larmor app
```

Open http://127.0.0.1:8642, paste the path to the `.fxmla` file, press Load.
The spectrum appears with the model overlaid; every parameter edit redraws the
model immediately. `Fit` runs the same engine and shows `± error` next to each
parameter. (You can also load a Bruker EXPNO folder the same way — LARMOR
never writes into it.)
