# Tutorial 4 — Batch fitting a composition series

*Time: ~25 minutes. Data: the five ¹¹B MAS spectra of a calcium composition
series — EXPNO 24 of the five 2026-01 sample folders under
`<NMR>/NMRFAM/DATA`, where `<NMR>` stands for the instrument-data root of the
development machine (`…/Desktop/WSU_work/NMR`). These are the developer's
instrument files: described here, not part of the repository. §6 is
self-contained (a synthetic three-spectrum series) and runs anywhere. Run the
commands from the repository root.*

A composition series is measured the same way for every sample, and the
question is how the site populations move with composition. LARMOR's batch
fit applies **one model** to every spectrum and lets only the amplitudes (and,
optionally, a few named parameters) change from spectrum to spectrum. This
tutorial builds the model on the first glass, batch-fits the five, releases
one parameter, computes per-spectrum errors, and does the same from the
command line.

## 1. Batch versus sequential

Two tools in LARMOR fit a series; they answer different questions.

- **Batch fit** (**Tools > Batch fit spectra…**, or the Explorer's
  **Batch fit selected…**): one shared model, every lineshape parameter held
  at its recipe value, amplitudes free per spectrum. Parameters that are
  known to drift can be *released* within a ±fraction. Use it when the sites
  are the same in every sample and only their populations change.
- **Sequential fit** (**Tools > Sequential fit…**, `larmor seqfit`):
  independent fits, each warm-started from its fitted neighbour, in forward
  and backward passes. Use it when positions or widths march along the
  series and a shared value would be wrong for every member.

The *Multi-dataset & co-fitting* manual, §6 and §8, describes both in
detail.

## 2. Look before you fit

```
larmor info <NMR>/NMRFAM/DATA/2026-01/01192026_SR31649_Base0Ca_SS_ALP/24
```

```
EXPNO: <NMR>\NMRFAM\DATA\2026-01\01192026_SR31649_Base0Ca_SS_ALP\24
nucleus: 11B   SFO1: 192.4307 MHz
pulse program: zg   TD: 7988   SW: 100000 Hz
MASR (acqus): 4200.0 Hz
title: 11B with short tip angle
CONFLICT: MAS rate: acqus says 4200 Hz but the title says 35714 Hz -- confirm before fitting
```
<!-- measured v0.12.1, 2026-09-22: `larmor info` on the EXPNO (the first line echoes the path as given) -->

The acqus file records a 4.2 kHz spinning rate while the third line of the
operator's title says `MASR 35.714 kHz`. The two disagree by far more than
2 %, so LARMOR prints the `CONFLICT` line and, when the spectrum is opened in
the app, resolves the rate to the *higher* of the two sources and flags it:
the indicator at the bottom right reads `⚠ MAS 35714 Hz — check!` until the
rate is confirmed in **Process > Experiment parameters…**. The MAS rate
enters every quadrupolar lineshape, so confirm it before fitting — here the
title is right, the probe was spinning at 35.7 kHz. All five EXPNOs of the
series carry the same acquisition and the same conflict.

## 3. Build the shared model on the first glass

The shipped ¹¹B recipe, `examples/pCABS2-4_11B.recipe.json`, has the model
this series needs: two *Amorphous (Gaussian Cq/eta dist.)* BO₃ sites and one
*Gauss/Lorentz* BO₄ line. It was fitted at 160 MHz and 20 kHz; the series is
at 192 MHz and 35.7 kHz, and the recipe's positions and quadrupolar
parameters carry over unchanged because they are physical quantities.

```
larmor desktop
```

1. **File > Open…** (Ctrl+O) the shipped recipe. It loads with its own
   spectrum (`examples/pCABS2-4/1118`).
2. **File > Open EXPNO / folder…** (Ctrl+Shift+O) and pick
   `…/01192026_SR31649_Base0Ca_SS_ALP/24`. Because a fit is open and the new
   spectrum brings none, LARMOR asks *Keep fit parameters?* — "A fit is
   already open. Keep the current lines and fit them against the new
   spectrum? Yes = keep the lines · No = start empty". Answer **Yes**: the
   three lines stay, and the nucleus, Larmor frequency and spin rate are
   taken from the new EXPNO. (The alternative route is **Decomposition >
   Apply recipe > Browse for recipe…** after opening the EXPNO first.)
3. Confirm the 35 714 Hz rate in **Process > Experiment parameters…**; the
   red indicator clears.
4. **Fit** (F5), then read the results strip (`RMSD … · χ²ᵣ …`, with the
   flags described in Tutorial 6, §2) and the `± error` column.
5. **File > Save recipe** (Ctrl+S) as `base0Ca_11B.recipe.json` next to the
   data. This file is the shared model for the batch.

No fitted values are quoted for this step: the fit was not re-run for this
tutorial. The result to expect is the two BO₃ lines near 16 and 11 ppm and
the BO₄ line near 0 ppm, with the BO₄ fraction read from **Decomposition >
Report (quantify)** (F6).

## 4. Batch fit in the app

1. In the *Explorer* dock press **Browse…** and choose the `2026-01` folder;
   expand the five sample folders and Ctrl-click their EXPNO `24` rows
   (Shift-click selects a range). Selecting fewer than two and pressing the
   button leaves the status bar saying
   `Ctrl/Shift-select at least two spectra in the Explorer first`.
2. Press **Batch fit selected…** (**Tools > Batch fit spectra…** is the same
   action). The dialog *Batch fit — one shared model, amplitudes per
   spectrum* opens with one panel per spectrum. The current fit is
   pre-loaded as the model; **Model from recipe…** loads
   `base0Ca_11B.recipe.json` instead.
3. Options worth knowing before the first run: `components` overlays each
   site's curve on every panel; `shared scale` puts all panels on one axis
   scale; **Fit baseline…** estimates and subtracts a baseline from every
   spectrum (**Reset** restores the raw data); **Save setup…** /
   **Load setup…** keep the release set and baseline choice for the next
   series.
4. **Fit**. Every lineshape parameter is held at the model value and each
   spectrum's amplitudes are fitted; the status line then reads
   `batch fit: 5 spectra, N shared parameters · mean RMSD …` and each panel
   shows its own RMSD (a panel whose RMSD stands out is flagged with ⚠ and a
   tooltip saying why). **Cancel** stops and discards the run; **Stop**
   stops and keeps the latest values.

The five-spectrum run was not re-run for this tutorial, so no RMSD or
population numbers are quoted; §6 shows the same outputs on a series that
runs anywhere.

## 5. Release, exclude, and compute errors

- **Release per spectrum (else held fixed):** the row of checkboxes names
  the parameters the model holds — `δiso (ppm)`, `Cq (MHz)`, `η`,
  `FWHM (ppm)` and so on. Tick `δiso (ppm)` and set the `±` box to 5 % to let
  each spectrum's positions move within ±5 % of the model value, then **Fit**
  again. The status line gains `· released isotropic_chemical_shift_ppm
  (±5%)`. A released parameter that ends at the edge of its ±range is
  telling you the model value does not fit that member of the series (§6
  shows exactly this).
- Right-click a panel and use **Exclude component** for a species that is
  absent in one glass: its amplitude is locked at zero in that spectrum and
  it disappears from the tables.
- **Error calculation:** choose `Covariance (from the fit)`, `Monte-Carlo
  (synthetic-noise refits)` (with `trials`) or `χ² profile (error analysis)`
  (with `points`) and press **Compute errors**; the three estimators are
  those of Tutorial 6, applied to every spectrum. **Export CSV…** writes one
  row per parameter with the columns `scope, site, label, param, value,
  stderr, sigma_pct, ci68_lo, ci68_hi, error_method, model, source_path`,
  plus a `population_pct` row per site and spectrum.
- **Save table…** writes `batch_table.csv` (shared values once, then the
  per-spectrum amplitudes, released parameters and populations); the
  checkbox `also save individual fits (.recipe.json) next to the CSV` is on
  by default. **Save individual fits…** writes one recipe per spectrum.
- **Publication bundle…** writes the whole batch to a folder you choose: the
  table, one recipe and one `_curves.csv` per spectrum (experiment exactly
  as fitted, model, residual, components), `manifest.csv` (source, SHA-256,
  EXPNO, NS, D1, SF/SR, window, RMSD, versions) and `README.txt` with the
  Methods paragraph — see the manual, §6.
- **Series plot…** opens *Series evolution*: pick one or more lines and a
  parameter (BO₄ population against sample, for instance), choose the error
  bars from the computed estimators, and **Export figure…** or
  **Export parameters (CSV)…**.

Two more tools read what the batch wrote. **Plotting > Plotting studio…**
has templates for a series — `Stacked series`, `Deconvolution grid`,
`Composition series (shaded component)`, `Composition trend` and
`Species distribution` — that take the saved table and the recipes next to
it. **Tools > Batch fit report…** (*Batch fit report — publication table,
plots & report*) takes any set of saved fits (**Add fits…**), re-fits them
for fresh errors and writes `table.csv`, `table.tex`, `report.md` and one
overlay figure per fit to a folder.

## 6. Same thing from the command line

This section runs anywhere: it writes three noiseless spectra of one
Gaussian line marching from 14 to 16 ppm, plus a deliberately wrong starting
model at 12 ppm. Save the following as `make_series.py` in an empty folder
and run it there:

```python
import numpy as np
from larmor.recipe import Recipe, SiteModel, Param
from larmor import engine

x = np.linspace(-20, 60, 400)
for k, pos in enumerate([14.0, 15.0, 16.0]):
    truth = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0,
                   sites=[SiteModel(model="gauss_lor", label="A", params={
                       "isotropic_chemical_shift_ppm": Param(pos),
                       "shift_fwhm_ppm": Param(6.0), "amplitude": Param(100.0),
                       "gl": Param(1.0, vary=False)})])
    _, y, _ = engine.simulate(truth, exp_ppm=x)
    with open(f"s{k}.csv", "w") as f:
        f.write("# nucleus = 11B\n# larmor_MHz = 160\n")
        f.writelines(f"{xi:.4f} {yi:.4f}\n" for xi, yi in zip(x, y))
model = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0,
               sites=[SiteModel(model="gauss_lor", label="A", params={
                   "isotropic_chemical_shift_ppm": Param(12.0, min=0, max=30),
                   "shift_fwhm_ppm": Param(5.0, min=0.1),
                   "amplitude": Param(80.0, min=0),
                   "gl": Param(1.0, vary=False)})])
model.save("model.recipe.json")
```

The `# key = value` header lines are how a plain CSV carries its nucleus and
Larmor frequency into LARMOR. Now the batch fit:

```
larmor batchfit s0.csv s1.csv s2.csv --model model.recipe.json -o out
```

```
batch fit: 3 spectra, 3 shared parameters · mean RMSD 0.1735
wrote 3 recipe(s) + batch_table.csv to out
```
<!-- measured v0.12.1, 2026-09-22: the three commands of this section, run in a scratch folder -->

Adding `--curves` writes the publication bundle next to the recipes and the
table:

```
larmor batchfit s0.csv s1.csv s2.csv --model model.recipe.json -o out --curves
```

and prints one more line after the two above:

```
wrote 3 curve file(s) + manifest.csv + README.txt to out
```

`out/s0_batch_curves.csv` holds `ppm`, `experiment` (the fitted array, written
exactly), `model`, `residual` and `s0_A` on the experimental axis;
`out/manifest.csv` records per spectrum the source file and its SHA-256, the
fit window, the RMSD and the software versions, and `out/README.txt` the
Methods paragraph. The table itself now also carries `model` and
`source_path` columns after the six above.

The mean RMSD is poor (0.17) because the shared model insists on 12 ppm for
lines that sit at 14–16 ppm; only the amplitudes could move, and they
absorbed the misfit by shrinking (`out/batch_table.csv` lists 90.6, 72.2 and
52.5 for a true amplitude of 100). Releasing the position within ±5 %:

```
larmor batchfit s0.csv s1.csv s2.csv --model model.recipe.json --release isotropic_chemical_shift_ppm --release-frac 0.05 -o out_rel
```

```
batch fit: 3 spectra, 2 shared parameters · released isotropic_chemical_shift_ppm (±5%) · mean RMSD 0.1485
wrote 3 recipe(s) + batch_table.csv to out_rel
```

All three positions end at 12.6 ppm — exactly 12 × 1.05, the edge of the
released range — and the RMSD barely improves. That is the at-the-edge
signal of §5: the model is wrong for the whole series, not just drifting.
The sequential fit, which lets every spectrum find its own position, is the
right tool here:

```
larmor seqfit s0.csv s1.csv s2.csv --model model.recipe.json --passes 2 -v -o sout
```

```
  pass 1 · spectrum 1: RMSD 0.0001
  pass 1 · spectrum 2: RMSD 0.0001
  pass 1 · spectrum 3: RMSD 0.0002
  pass 2 · spectrum 3: RMSD 0.0002
  pass 2 · spectrum 2: RMSD 0.0001
  pass 2 · spectrum 1: RMSD 0.0001
sequential fit: 3 spectra, 2 passes · mean RMSD 0.0001044 → 0.0001044
wrote 3 recipe(s) + seq_table.csv to sout
```

`sout/seq_table.csv` recovers 14.000, 15.000 and 16.000 ppm with widths of
6.000 ppm and amplitudes of 100.0 — the truth, to the fourth digit, since the
spectra are noiseless. The second pass changes nothing here; on real data it
is the backward sweep that catches a member the forward pass fitted badly.

Both commands accept EXPNO folders as well as CSVs, so the five-glass series
of §4 is fitted headlessly with

```
larmor batchfit <the five EXPNO 24 paths> --model base0Ca_11B.recipe.json --window 30 -30 -o out
larmor seqfit <the five EXPNO 24 paths> --model base0Ca_11B.recipe.json --window 30 -30 --passes 2 -v -o sout --curves
```

(`--window` is high ppm then low ppm; the model recipe's own window is used
when it is omitted). The outputs are the same `batch_table.csv` /
`seq_table.csv` and one recipe per spectrum, ready for **Tools > Batch fit
report…** and the Plotting studio.

## 7. What a shared model can and cannot claim

A batch fit with everything but the amplitudes held is a **hypothesis test**:
"the same sites, in the same environments, with changing populations". When
it fits every member of the series to within the noise, the populations are
comparable across the series precisely because nothing else moved. When one
member's RMSD stands out, or a released parameter ends at the edge of its
range, the hypothesis failed for that member — and the honest answers are a
sequential fit (positions and widths free, populations no longer strictly
comparable) or a co-fit that ties only the parameters physics says are
shared (**Decomposition > Advanced > Co-fit datasets…**, Tutorial 5, §10).

## Where to read more

**Help > User manuals**: *Multi-dataset & co-fitting*, §5 (batch report),
§6 (batch fit), §7 (series figures) and §8 (sequential fit).
