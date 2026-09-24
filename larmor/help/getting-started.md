# Getting started with LARMOR

> **LARMOR** — *Lineshape Analysis & Refinement for Magnetic-resonance Of solids*.
> An open desktop successor to dmfit for solid-state NMR: dmfit's fitting
> fluency, ssNake's processing depth, and the **mrsimulator** physics stack, with
> uncertainties on every number and reproducible, diffable recipes.

This guide is the front door. It explains the ideas the whole program is built
on — the **workspace model**, the **Explorer**, and the *open-anything* import —
and then points you to the manual for whatever experiment is in front of you.

---

## The three ideas

| Idea | What it means for you |
|---|---|
| **Open anything** | 1D, 2D, a raw FID or `ser`, a dmfit `.fxmla`, a plain CSV — LARMOR reads it, shows it, and *then* offers the right processing. Nothing is rejected at the door. |
| **Workspaces** | Every dataset you open lives in its own switchable workspace (like TopSpin windows). Extract a trace from a 2D and it becomes a new workspace; the map stays where it was. |
| **Reproducible recipes** | Processing and fitting are stored as a **recipe** replayed from the raw data. A saved recipe *is* the analysis — it reopens exactly, and it diffs cleanly in version control. |

Every fitted number is reported **with its uncertainty**, or with the reason
there is none. Instrument files are opened **read-only** — LARMOR never modifies
into your Bruker data.

---

## 1 · Opening data

**File ▸ Open**, the **Explorer** dock on the left, or **drag-and-drop** onto the
plot. LARMOR's universal reader accepts:

| Source | Notes |
|---|---|
| Bruker `1r` / `2rr` | processed 1D / 2D spectrum (or point at the `pdata` folder, or the EXPNO) |
| Bruker `fid` / `ser` | raw time-domain — opens with a processing preview so you apodize/phase, then transform |
| dmfit `.fxmla` | a dmfit model + data (1D and MQMAS) imported directly |
| LARMOR `.json` | a saved recipe (data reference + full processing + model) |
| `.csv` / `.txt` / `.dat` | a two-column *(shift, intensity)* spectrum, with an optional `# key=value` metadata header |

The reader detects **1D vs 2D**, **time vs frequency** domain, real vs
**hypercomplex** (it loads the `2rr`/`2ri`/`2ir`/`2ii` quadrants when present),
and it distinguishes a genuine spectroscopic 2D from a **pseudo-2D** relaxation
series (a `vdlist`/`vclist` with a placeholder F1 axis). What you opened decides
where it lands:

- a **1D frequency** spectrum → the workbench (process & fit);
- a **2D** dataset → the contour view;
- a raw **FID** → a magnitude-FT preview (EM 100 Hz) whose pipeline is recorded —
  *FID ⇄ spectrum* (Ctrl+T) shows the FID and re-apodizes live; untick *magnitude*
  and phase in the Processing panel, or use *File ▸ Open FID* for the full dialog;
- a raw **`ser`** → a 2D preview (or the guided relaxation tool if it's a series).

**Several processed datasets.** If you open a bare EXPNO that has more than one
`pdata/N` (several TopSpin processings), LARMOR asks **which one** (labelled
1D/2D + title). Re-open the sample to work with a different processing — it opens
in its own workspace, so you can keep two processings of one experiment side by
side.

**MAS-rate check.** The MAS rate is read from three sources: `acqus` MASR, the
operator's title line (`MASR 20 kHz`) and, on NMRFAM data, the booking sidecar
`experiment_addenda.xml` written when the instrument was reserved. Two sources
agreeing within 2 % settle the rate and the third is reported as outvoted in
the load message (on that spectrometer `acqus` MASR is a constant leftover, so
it is usually the one outvoted). When the sources all disagree the highest is
taken and a **red "⚠ MAS … — check!"** badge appears bottom-right; a title rate
outside 1–150 kHz (`35741 kHz`, `35.714 Hz`) is re-read in the other unit and
flagged as a probable typo; nothing found at all falls back to 35 714 Hz,
flagged. **Process ▸ Experiment parameters…** (or double-click the strip or the
badge) lists the sources side by side with a **Use** button each, can
**Measure** νrot from the spinning sidebands of the spectrum on screen, and
offers **Remember for this session**: the confirmed rate then applies to every
later spectrum of the same session folder, rotor and nucleus that shows the
same three source values, without the badge (**Forget** undoes it).

## 2 · The Explorer

The Explorer browses a sample folder and **auto-identifies** each experiment —
nucleus, 1D/2D, and kind (single pulse, MQMAS, satrec, QCPMG…) — by reading
`acqus` and the pulse program, so you can see at a glance what an EXPNO holds
before opening it. **File ▸ Open sample** points it at a folder.

**Session inventory (Series ▸ Session inventory…).** One month folder read into
a sample × nucleus grid. Every EXPNO gets a role — production, candidate,
short, setup, failed, arrayed, 2D, reference, unprocessed — and the production
spectrum of each block is pre-picked by the operator's rule: the highest EXPNO
that has a `pdata/1/1r`, demoted when its NS is a small fraction of the block's
maximum (the fraction is adjustable, 0.25 by default) or when its title says
power check / optimisation / test / failed. The sample name is the normalised
folder — date, rotor ID and the `_SS_ALP` operator suffix removed, so
`01192026_SR31649_Base0Ca_SS_ALP` is `Base0Ca` — or, for sets that keep one
EXPNO per sample, the title's `Sample …` line; never a pulse note such as
"11B with short tip angle". The same rule names every Bruker spectrum opened in
the workbench (plot title, workspace row, default save names, batch scopes);
the title's first line is kept in the recipe's provenance. Flags mark a title
whose rotor ID or sample name contradicts the folder, a leading nucleus that is
not NUC1, `zg` in the title of a non-zg pulse program, two folders sharing one
sample name, and an unprocessed EXPNO with more scans than the pick; flags
never move a pick. A tick in the detail table overrides a pick, and **Batch fit
picks…** (or **Sequential fit picks…**) hands the picks of the chosen nucleus
over in sample order. The Explorer's right-click on a month or sample folder
opens the same window; `larmor inventory <month>` prints the grid and the
demoted or flagged rows, `--picks --nucleus 31P` the picks' paths only.

## 3 · Workspaces

The **Workspaces** dock lists everything you have open. Each entry carries an
icon for its kind: bare spectrum (∿), spectrum-with-fit (⤳), 2D map (▦), and
the session rows a project keeps alongside them — a figure kept from the
Plotting studio (◫) and a batch-fit session (☷). Click to switch to a
document, **Close** to free it, **Save** to write its recipe; a figure or batch
row opens in its dialog on **Enter**, double-click or **Save**, and is updated
when that dialog closes. Opening a 2D, or pulling a trace off a map, spawns a
**new** workspace so your existing fit is never disturbed; **Back to 2D map**
(Ctrl + 2) returns to the parent map. Snapshots are lightweight (arrays + a
recipe reference), and closing a workspace frees it. **File ▸ Save project…**
(Ctrl + Shift + P) keeps the whole dock as one `.larproj.json` file; a project
opens from **File ▸ Open project…** or plain **File ▸ Open…**.

## 4 · Processing, fitting, tools

Once your data is open, the workflow is the same shape everywhere — process,
then fit, then measure/export — but the details depend on the experiment. The
**Process** panel is live (edits re-apply as you type); the **Fit-Parameters**
table at the bottom is a dmfit-style spreadsheet with paddles on the plot; the
**Edit** menu builds the model, **Fit** runs and checks it, **Series** and
**Tools** hold the experiment-specific machinery (see the menu map below).
Under the spectrum, a **fit-health strip** states whether the last fit's
numbers can be read as they stand — unphysical values, degenerate pairs,
structured residual, parameters at bounds, missing error bars, a line cut by
the integration window, an under-relaxed recycle delay, a too-long pulse — and
each chip opens the matching detail (**F7** lists them all).

**Appearance.** **View ▸ Theme** offers ten colour presets (five light, five dark
— Light, Sepia, Solarized, Nord, Slate, Ocean, High-Contrast, …). The choice is
remembered between sessions and re-colours everything — widgets, plot, curves,
site markers — while keeping text, buttons and markers readable (every preset is
contrast-checked). Pick whichever is easiest on your eyes.

**Finding a command.** **? ▸ Command palette…** (Ctrl + Shift + P) lists every
menu and toolbar entry — menu path, label and shortcut — and filters as you type.
Fragments and ASCII spellings work (`czj`, `theme dark`, `F5`, `chi2` for χ²); a
check mark shows the current state of toggles; greyed entries are not available
yet (usually: open a spectrum first). Enter runs the highlighted entry; with
nothing typed the list is a complete map of the menus.

Each experiment has its own manual with worked steps and the science behind it:

| If you have… | Read |
|---|---|
| a normal 1D spectrum | **1D spectra — processing & fitting** |
| a raw FID / echo you must transform | **1D spectra** + **Processing reference** |
| a QCPMG echo train (broad lines) | **QCPMG** |
| a T₁/T₂ relaxation series | **Relaxation (T1/T2)** |
| a 2D MQMAS map | **MQMAS (2D)** |
| an HMQC / DQ-SQ correlation | **HMQC & correlation** |
| any 2D to phase / measure | **2D processing** |
| several datasets to compare or co-fit | **Multi-dataset & co-fitting** |
| a question about a lineshape model | **Lineshapes — models & physics** |
| a question about a processing step | **Processing reference** |
| a step-by-step walkthrough on real data | **Help ▸ Tutorials** (seven tutorials) |

All are under **? ▸ User manuals**; the reference documents are direct **?**-menu
items. Every tool with a **Help** button opens the matching section.

## 5 · Menu map

Nine menus, each split into groups; every row shows its explanation as a
tooltip, and **? ▸ Command palette…** (Ctrl + Shift + P) searches all of them.

| Menu | Groups |
|---|---|
| **File** | Open (spectrum, sample, EXPNO, FID, Varian, recent) · Overlay a spectrum, Watch the source file · Open / Save project · Save fit, Save fit as, **Export ▸** (spectrum as CSV, figure, plot image, copy plot, report table as CSV, LaTeX table, publication bundle) · Quit |
| **Edit** | Undo, Redo · New fit · **Add line ▸** (Simple · Quadrupolar · Disordered · CSA · Other), Add a line at every peak, Add function line, Add background spectrum, **Apply recipe ▸** · **Spinning sidebands ▸** (detect, add manifold, add a shifted copy, offer on load) · **Constraints ▸** (label from literature ranges, restrict around current values, save / apply a constraint set) · Add fit zone, Clear zones |
| **Process** | Experiment parameters, Processing steps, Show processing panel · **Phase ▸** (autophase, drag to phase) · **Baseline ▸** (polynomial, iterative, 2-point, subtract averages) · **Reference ▸** (calibrate axis, measure Δ, referencing audit) · **Region / algebra ▸** (integrals & measurements, subtract a spectrum, WURST profile, stitch VOCS) · FID ⇄ spectrum, **Display channel ▸** · Reset to original |
| **Fit** | Simulate (F9), Fit (F5), Auto fit · Report (F6), Fit health details (F7), **Errors ▸** (χ² profile, Monte-Carlo, parameter correlations, χ² map), Compare with a saved fit · Co-fit datasets, Predict at another field, **MQMAS ▸** (2D viewer / fit, F1 reference), **Fit settings ▸** (computing parameters, completion threshold, animate fits) |
| **Series** | Batch fit spectra, Sequential fit, Batch fit report · Session inventory, Experimental section, Compare acquisition parameters |
| **Tools** | NMR table, Conversion tools · Herzfeld–Berger, Read static pattern, Czjzek distribution · Import DFT tensors, **Relaxation ▸** (T1/T2 series, per-site, variable temperature), **QCPMG ▸** (echo train → spectrum, infinite-field δiso, batch), REDOR |
| **View** | Residual, Components, Component labels, Paddles, Overlays, Clear overlays, Literature shift ranges, Scroll nudges fit values · **Zoom ▸** (full, sites, back to 2D map) · **Axis unit ▸**, **Czjzek width display ▸** · **Panels ▸**, **Theme ▸**, **Text size ▸** |
| **Plotting** | Plotting studio, Plot current spectrum, New 2D contour plot |
| **Help (?)** | Command palette · **User manuals ▸**, **Tutorials ▸** · About LARMOR, More… |

---

## Conventions used everywhere

- **Chemical shift** increases to the **left** (decreasing frequency, IUPAC δ
  scale). In 2D, F2 (direct) high-shift is left, F1 (indirect) high-shift is top.
- **Quadrupolar coupling** is reported as $C_Q = e^2qQ/h$ (MHz) and asymmetry
  $\eta_Q\in[0,1]$; for disordered sites the **quadrupolar product**
  $P_Q = C_Q\sqrt{1+\eta_Q^2/3}$ is the invariant that keeps its meaning.
- **Shielding anisotropy** uses the **Haeberlen** convention (δ_iso, ζ, η_CS).
- **Widths** may be entered in **ppm or Hz** (`300Hz`, `1.5kHz`, `2ppm`);
  LARMOR converts with the Larmor frequency.

---

## References

- D. Massiot, F. Fayon, M. Capron, I. King, S. Le Calvé, B. Alonso, J.-O. Durand,
  B. Bujoli, Z. Gan, G. Hoatson, "Modelling one- and two-dimensional solid-state
  NMR spectra", *Magn. Reson. Chem.* **40**, 70 (2002). *(dmfit)*
- S. G. J. van Meerten, W. M. J. Franssen, A. P. M. Kentgens, "ssNake: A
  cross-platform open-source NMR data processing and fitting application",
  *J. Magn. Reson.* **301**, 56 (2019). *(ssNake)*
- D. J. Srivastava, P. J. Grandinetti *et al.*, **mrsimulator**,
  github.com/deepanshs/mrsimulator. *(the physics engine)*
- M. H. Levitt, *Spin Dynamics: Basics of Nuclear Magnetic Resonance*, 2nd ed.,
  Wiley (2008). *(general reference)*

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
