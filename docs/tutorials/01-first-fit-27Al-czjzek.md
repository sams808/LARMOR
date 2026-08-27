# Tutorial 1 — A first fit: ²⁷Al MAS spectrum of a glass

*Time: ~15 minutes. Everything used here ships with the repository: a ²⁷Al MAS
spectrum of a Ca-aluminoborosilicate glass (`examples/pCABS2-4/3616`, a Bruker
EXPNO) together with the dmfit fit that was originally made for it. Run the
commands from the repository root.*

This tutorial loads a spectrum, inspects an existing dmfit model, runs a fit,
and reads the uncertainties, which dmfit does not report.

## 1. Look before you fit

Open a terminal and activate the environment:

```
conda activate larmor
larmor info examples/pCABS2-4/3616
```

```
EXPNO: examples\pCABS2-4\3616
nucleus: 27Al   SFO1: 130.3230 MHz
pulse program: zg   TD: 2048   SW: 100000 Hz
MASR (acqus): 26000.0 Hz
```

A single-pulse ²⁷Al acquisition at 26 kHz MAS. The dmfit fit made for this
spectrum sits next to the processed data; LARMOR reads it directly:

```
larmor info examples/pCABS2-4/3616/pdata/1/1r.fxml
```

```
dmfit fit file (version 20200306, mode 'Fit 1D')
dimension F2: 27Al at 130.318 MHz, 3 lines
  [0] CzSimple   pos=64.10 ppm  sCZ_CQ=3508 kHz
  [1] CzSimple   pos=33.80 ppm  sCZ_CQ=3965 kHz
  [2] CzSimple   pos=3.88 ppm  sCZ_CQ=3020 kHz
```

Three `CzSimple` lines — dmfit's simple Czjzek distribution, the standard
model for a distribution of quadrupolar couplings in a disordered material.
The three positions (~64, ~34, ~4 ppm) are the classic AlO₄ / AlO₅ / AlO₆
coordination assignment.

## 2. The recipe

A LARMOR fit lives in a small, readable JSON "recipe": the model only, with
the data referenced by path and SHA-256 hash instead of copied. The shipped
`examples/pCABS2-4_27Al.recipe.json` was produced by importing the dmfit file
above (`larmor import <file.fxml> -o out.recipe.json` does the conversion).
Open it in any text editor; every parameter looks like:

```json
"sigma_Cq_MHz": { "value": 1.754, "stderr": null, "vary": true, "min": 0.05, "max": null, "expr": null }
```

The Czjzek width was converted with `sigma = sCZ_CQ / 2`: dmfit and
mrsimulator use conventions that differ by exactly a factor of two (see the
validation report for the numerical check).

## 3. Fit

```
larmor fit examples/pCABS2-4_27Al.recipe.json --window 150 -80 --plot fit.png
```

The window `150 -80` (high ppm, low ppm) covers the central-transition region.
The first run builds the Czjzek simulation kernel (~15 s); after that,
iterations are milliseconds. At the end of the report:

```
normalized RMSD: 0.0462
```

and `fit.png` shows the experiment (black), total fit (red), and each site's
contribution (dashed).

## 4. Read the uncertainties

Look at each site's `stderr` fields in the updated recipe, or read them from
the fit report. On this dataset:

- The AlO₄ site is well determined: δiso = 63.8 ± 0.6 ppm,
  σ(Cq) = 1.74 ± 0.14 MHz.
- The AlO₅ site is marginal: δiso = 31.3 ± 6.0 ppm.
- The AlO₆ site is not determined at all: δiso = −2 ± 39 ppm, and its σ(Cq)
  error is ten times its value. The 1D lineshape simply does not contain
  enough information to pin three overlapping Czjzek sites independently.

That last bullet isn't a failure of the fit. To quantify the minor sites you
need more data — an MQMAS spectrum (one ships in `examples/pCABS2-4/3620`), a
second field, or a constraint from chemistry. dmfit reports the same parameter
values with no error bars, so this degeneracy stays invisible there.

Tutorial 2 shows how to add exactly those constraints.

## 5. Same thing, interactively

```
larmor desktop
```

In the desktop application, use **File > Open EXPNO / folder…** and pick
`examples/pCABS2-4/3616` — the explorer also lists the `.fxml` fit under the
proc, and double-clicking it loads the model over the spectrum. Every
parameter edit redraws the model immediately; **Fit** (F5) runs the same
engine and shows `± error` next to each parameter. Fits can be saved next
to the data — LARMOR never modifies the acquired files themselves.
