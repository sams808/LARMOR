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
MAS (booking sidecar): 35714 Hz
title: 11B with short tip angle
CONFLICT: MAS rate: acqus says 4200 Hz but the title says 35714 Hz and the booking sidecar agrees -- acqus outvoted, using 35714 Hz
```
<!-- measured v0.13.0 (N9 worktree), 2026-09-23: `larmor info` on the EXPNO (the first line echoes the path as given) -->

The acqus file records a 4.2 kHz spinning rate — a leftover of that
spectrometer's rotor controller, the same 4200 Hz in almost every NMRFAM
EXPNO — while the third line of the operator's title says `MASR 35.714 kHz`
and the NMRFAM booking sidecar (`experiment_addenda.xml`, written when the
instrument was reserved) says 35.714 kHz too. Two independent sources agree,
so LARMOR settles the rate at 35 714 Hz, prints the `CONFLICT` line with
acqus reported as *outvoted*, and shows no red indicator when the spectrum
is opened in the app. All five EXPNOs of the series carry the same
acquisition and the same outcome. The MAS rate enters every quadrupolar
lineshape, so it is worth the glance: the case where the indicator *does*
appear is a three-way disagreement, such as the 2026-05 ³¹P series where
the title says 20 kHz and the booking 22 kHz — there the highest is taken,
`⚠ MAS 22000 Hz — check!` appears, and **Process > Experiment parameters…**
lists the sources side by side; confirming once with *Remember for this
session* covers every EXPNO of that session, rotor and nucleus with the same
three source values.

The title also says `short tip angle` without a number. EXPNO 23 of every
folder is the ¹¹B saturation-recovery measurement, and its TopSpin result
(`pdata/1/ct1t2.txt`) gives T₁ = 4.64 s for the BO₃ region and 3.90 s for
the BO₄ region of Base1Ca against a recycle delay D1 = 14 s. After the fit
the fit-health strip therefore shows `D1 = 3.0–3.6 T1 → 95–97 % (90°
assumed)`; typing the ¹¹B 90° pulse in **Process > Experiment
parameters…** turns it into the steady-state value for the short pulse
actually used. Six of the ten glasses of this session (Base1Ca, Base2Ca and
four of the fluorinated ones) were acquired with their slowest boron site
at 2.5–3.3 T₁ (D1 = 9 … 18 s), the two sites relaxing differently in every
glass; the other four sit at 4.7–7.3 T₁. A caveat to carry into the Methods
text, not something a fit can repair.

The five EXPNOs share the pulse program and NS = 256, but not their
processing or their recycle delay:

```
larmor compare <NMR>/NMRFAM/DATA/2026-01/01192026_SR31649_Base0Ca_SS_ALP/24 <NMR>/NMRFAM/DATA/2026-01/01202026_SR31649_Base1Ca_SS_ALP/24 <NMR>/NMRFAM/DATA/2026-01/01202026_SR31649_Base2Ca_SS_ALP/24 <NMR>/NMRFAM/DATA/2026-01/01202026_SR31648_Base3Ca_SS_ALP/24 <NMR>/NMRFAM/DATA/2026-01/01202026_SR31649_Base4Ca_SS_ALP/24
```

```
⚠ processed differently: LB 0 / 100 Hz (EM) · acquired differently: D1 varies 12.5–36 s — populations from one shared model are not strictly comparable
group        key  label    level  majority  01192026_SR31649_Base0Ca_SS_ALP  01202026_SR31649_Base1Ca_SS_ALP  01202026_SR31649_Base2Ca_SS_ALP  01202026_SR31648_Base3Ca_SS_ALP  01202026_SR31649_Base4Ca_SS_ALP
processing   LB   LB (Hz)  check  100       0 *                              0 *                              100                              100                              100
acquisition  D1   D1 (s)   check  12.5      12.5                             14                               18                               29                               36
* differs from the series majority
```
<!-- measured v0.13.0, 2026-09-23: `larmor compare` on the five EXPNO 24 folders -->

Base0Ca and Base1Ca were processed with LB = 0 Hz and the other three with a
100 Hz exponential window — the majority, but the worse choice for glasses
(*Glass fitting*, §2) — and the recycle delay runs from 12.5 to 36 s. The
Batch fit dialog shows the same in an amber line under the panels;
**Reprocess all from fid…** with the window set to none (or GM, the glass
rule) removes the LB difference for every member, while the D1 difference
cannot be reprocessed away and is recorded as a note in each saved fit.

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
3. Open **Process > Experiment parameters…** to see the three sources side
   by side (acqus 4 200 Hz outvoted by the title and the booking sidecar at
   35 714 Hz); nothing to confirm here, the rate is already settled.
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
   series. Under the panels, the amber comparability line says how the five
   differ (LB and D1, as §2 measured); **Details…** opens the parameter table
   and **Reprocess all from fid…** rebuilds all five from their fids with one
   window and TDeff before fitting.
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
  stderr, sigma_pct, ci68_lo, ci68_hi, error_method, model, source_path,
  vary, min, max, expr, at_bound`, plus a `population_pct` row per site and
  spectrum. The last five are the row's status: a shared row reads
  `vary = False` because the batch holds it, a released row carries its
  ± window in `min` / `max`, and `at_bound` reads `min` or `max` when the
  value stopped at that edge — the same cell shows ‡ in the results table,
  with the remedy in its tooltip.
- **Save table…** writes `batch_table.csv` (shared values once, then the
  per-spectrum amplitudes, released parameters and populations, with the
  same `vary, min, max, expr, at_bound` columns); the
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
<!-- measured v0.13.0 (N5 worktree), 2026-09-23: the four commands of this section, run in a scratch folder; transcripts unchanged from v0.12.1, at_bound column read from out_rel/batch_table.csv -->

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
released range — and the RMSD barely improves. The table says so in words:
`out_rel/batch_table.csv` reads `max` in its `at_bound` column for all
three position rows, with the released window `11.4` / `12.6` in `min` /
`max` (the shared width and `gl` rows read `vary = False`, held by the
batch). That is the at-the-edge
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
