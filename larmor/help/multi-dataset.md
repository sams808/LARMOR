# Multi-dataset — compare & co-fit

> Real assignments rarely come from one spectrum. LARMOR lets you **overlay**
> datasets for comparison and **co-fit** several measurements of the same sample
> at once — a 1D MAS spectrum with its MQMAS map, or the same nucleus at two
> magnetic fields — sharing the parameters that *must* agree while decorrelating
> those that need not.

---

## 1 · Overlay & compare (the cockpit)

The **Datasets** dock overlays spectra behind the active one for visual
comparison — a composition series (e.g. LAW3Cl0→4Ca), before/after processing,
or a reference. **＋ Add spectrum to compare…** draws each overlay in its own
colour with per-overlay **visible / colour / remove** controls and a global
**stack offset**. **Compare acquisition…** (beside the add button, live once
an overlay with a Bruker source exists) opens the comparability table of §6
over the active spectrum and every overlay — `acqus` and `procs` side by
side — and an overlay acquired or processed unlike the others carries a ⚠
marker with the differing parameter under its name. Promote an overlay to
**active** (it becomes the fit target)
and the previous active demotes back to an overlay. The active spectrum is always
the single object the 1D fitter works on, so overlays never disturb a fit.

**Quick overlays, to check something against the spectrum on screen.** Three
ways add a compared spectrum without leaving the main window, and none of
them touches the active spectrum or its fit: **Shift + drop** one or several
files onto the plot (a plain drop still *opens* the file), **File ▸ Overlay a
spectrum…** (Ctrl + Shift + A), or right-click an EXPNO in the Explorer and
choose **Overlay on the current spectrum**. **View ▸ Overlays** (Ctrl + Shift + V)
hides and shows every compared spectrum at once without removing it — adding
a new one turns the display back on — and **View ▸ Clear overlays** drops them
all. In the Datasets dock, **match height** scales each overlay so its maximum
equals the active spectrum's, so shapes compare when intensities do not
(display only: the stored data, the fit and every export are unchanged);
**stack offset** spreads them out instead. The status bar names what was
overlaid; a file that cannot be read as a 1D spectrum is refused in the
status bar, never with a dialog.

**File ▸ Save project…** (Ctrl + Shift + P) captures the whole session as one
reopenable `.larproj.json` file — every row of the **Workspaces** dock, in
order: 1D spectra with their processing, fit *and overlays* (overlays by
reference); 2D maps **by reference** — the source path (absolute and
project-relative) plus the recorded phase / shear / transpose / reverse /
symmetrize / calibrate operations, the contour settings and the projection
overlays — replayed on the freshly loaded data when the project reopens;
figures kept from the Plotting / Figure studio; and batch-fit sessions with
their spectra, model, release / baseline / error settings and fitted results.
**File ▸ Open project…** (or **File ▸ Open…** on a `.larproj.json`) restores it
all; opening into a session that already holds a fit offers **Replace / Add /
Cancel**. Anything relocated or missing — a moved 2D source, an overlay or a
batch spectrum that cannot be found — is listed once in a "Project opened"
box rather than silently dropped; the fitted 2D model overlay is redrawn by
**Fit** rather than stored. The co-fit page (§3) is not part of a project.

## 2 · Background subtraction (a related but different tool)

**Process ▸ Subtract a spectrum (background)** removes a *measured* background
(an empty rotor, a probe/impurity signal) from the sample:

$$S_\text{clean}(\nu) = S_\text{sample}(\nu) - k\,S_\text{bg}(\nu),$$

with $k$ from a least-squares fit over the view (or set by hand) and an optional
shift. The result **replaces** the workbench spectrum and can be saved as a
reopenable CSV. *(Distinct from the `spectrum` fit component, which instead fits
a background's amplitude as one term of a model — see the Lineshapes reference.)*

---

## 3 · Co-fitting

**Decomposition ▸ Co-fit** puts the workspace into a **split view**: the 1D panel
and the 2D panel side by side, each with **its own parameter table**. You add the
second dataset with **＋ Add / replace dataset** (a 1D or a 2D/MQMAS file).

### Why decorrelate?

Two experiments on the *same* sample are still *different measurements*. A 1D MAS
spectrum and an MQMAS map of the same aluminoborosilicate glass do not, in
practice, share identical fitted parameters — resolution, second-order effects
and referencing differ — so forcing every parameter to be common biases both
fits. LARMOR therefore keeps a **separate recipe per dataset** and lets you
choose, parameter by parameter, what is **tied**:

- The **tie bar** lists only the parameters that actually influence the selected
  lineshape (a Czjzek co-fit shows δ_iso / σ(Cq) / dCS / line — not every possible
  parameter). Tick one to tie it across both datasets; untick to let each fit it
  independently.
- A **tied** parameter is optimised as one shared value (lmfit `expr` linkage);
  an **untied** one is free in each dataset.
- Sites keep the **same colour** across the 1D and 2D panels, so Al⁽⁴⁾ is the same
  colour everywhere.

### Two common cases

- **1D + MQMAS of one sample.** Tie $\delta_\text{iso}$ and the Czjzek widths you
  trust from the high-resolution MQMAS, but let the 1D keep its own line
  broadening. This anchors the isotropic shifts while respecting that the 1D MAS
  envelope is broader.
- **Multi-field (same nucleus, two B₀).** The **quadrupolar** parameters ($C_Q$,
  $\eta$, $P_Q$) are field-independent and *should* be tied; the second-order
  quadrupolar shift scales as $1/\nu_0^2$, so tying $C_Q$ across fields is a
  powerful constraint that separates it from the chemical shift.

### Running it

- **Preview** simulates both panels at the **current** values (no optimisation),
  auto-scaling each overlay and auto-aligning the 2D F1 reference — so you can
  hand-tune δ_iso / σ / dCS / F1-ref and *watch* before committing.
- **Run co-fit** minimises the joint residual across datasets. Amplitudes are
  **pre-scaled** per dataset (a 1D in raw counts and a normalised 2D map start on
  very different scales), and a **progress bar** shows iteration and RMS. Fitted
  values, with errors, are written back into each panel's table.
- **Close** adopts the 1D recipe as the main-window model.

Show/hide, bounds, `Fix`, links and paddles all work in the co-fit tables exactly
as in the normal Fit-Parameters spreadsheet (see **1D spectra**).

---

## 4 · Predicting another field

**Decomposition ▸ Predict at another field** re-simulates the current model at a
target ¹H frequency into a new workspace — useful to plan an experiment, or to
sanity-check that a fit's $C_Q$/$\delta_\text{iso}$ split behaves correctly when
you change $B_0$ (the quadrupolar shift moving as $1/\nu_0^2$).

---

## 5 · Batch fit report (publication table + plots)

Once you have a **set of finished fits** — a whole glass series fitted the same
way — **Tools ▸ Batch fit report** turns them into a paper-ready package in one
pass, so you never hand-copy numbers.

1. **Add fits.** Point it at the saved fits: LARMOR `.recipe.json`, dmfit
   `.fxmla`, or `.larproj`. Each carries its own data, so nothing else is
   needed. Ideally they share a **nucleus and acquisition** — the tool flags
   mixed nuclei or fields, because a table is only comparable within them.
2. **Choose the errors.** *Covariance* (fast — the lmfit standard errors) or
   *Monte-Carlo* (the parametric bootstrap of the Errors tools, slower but
   honest for correlated/glassy fits). LARMOR **re-fits every dataset** so the
   errors are fresh, not whatever was saved.
3. **Pick a folder and Generate.** It writes:
   - **`table.csv`** — every site's δ_iso, C_Q (or σ, with the derived
     `C_Q = 2σ` and `√⟨P_Q²⟩`), η, width and **population % — each with its
     error** — one row per site, machine-readable. A `… flag` column beside
     each value reads `fixed`, `linked`, `at_min` or `at_max` (blank for a
     free value), so a filter finds every held or pinned number.
   - **`table.tex`** — the same table as a LaTeX `tabular` for direct inclusion;
     cells are marked † (held fixed), ‡ (finished at a bound — the value
     shown IS the bound) and § (linked), with one legend row.
   - **`report.md`** — a Markdown report: the table (same markers), the
     nucleus/field summary, the error method, a **per-fit overlay**
     (experiment + model + components + residual), and a *Constraints*
     section listing, per sample, the fixed parameters, the at-bound ones
     with their bound and the linked ones with their expression — the
     paragraph a referee asks for.
   - **`figures/*.png`** — the individual overlays.
   - **`acquisition.csv`** / **`acquisition.tex`** (the *acquisition table
     (Table S1)* checkbox, on by default) — the acquisition parameters of every
     fit's spectrum from its recipe's acquisition block, one row per EXPNO, with
     a `# varies` row / bold cells for every column that differs across the
     series; `report.md` gains an **`## Acquisition`** section with the
     Experimental paragraph (varying parameters as ranges) and the table.
     Nothing is written for CSV / dmfit sources, which carry no block.

The columns are model-aware (Czjzek sites report σ and the field-independent
`√⟨P_Q²⟩`; discrete/Amorphous sites report C_Q and η) and populations come from
the same integral-over-the-window quantification as **Report** (F6). It is the
fastest route from "a folder of fits" to a table you can paste into a manuscript.

## 6 · Batch fit — one shared model, many spectra (1D)

For a **series measured the same way** (a composition series, a time course),
you often want *one* set of lineshape/position parameters describing every
spectrum, with only the **amplitudes** free — the sites are the same, their
populations change. **Ctrl/Shift-click** the spectra in the Explorer and press
**Batch fit selected…** (or *Tools ▸ Batch fit spectra*). Tutorial 4
(Help ▸ Tutorials) runs this on a five-glass composition series and on a
synthetic series from the command line.

**Comparable before shared.** A shared model assumes the series was
*measured and processed alike*, and TopSpin never says otherwise. When the
dialog opens it reads every Bruker member's `acqus`, `pdata/N/procs` and
TopSpin's own command trail (`auditp.txt`) and compares the series against
its **majority** — the most frequent value; an exact tie takes the first
spectrum's value and says so; a parameter where every member differs reports
once as *varies min–max* without singling anyone out. Compared and flagged:
the window and its LB / GB / SSB, TDeff, SI, FCOR, the phase mode and the
`abs` baseline order (ABSG, BC_mod), the referencing (SF, in ppm against the
Referencing audit's tolerance) and, on the acquisition side, the pulse
program, SW, TD, D1 (5 %), P1, PLW1 and the probe. NS, RG, the date, the
PHC0 / PHC1 values and the audit command line are reported but never flagged —
they are legitimately per spectrum (NS and RG are still named in the verdict
when they differ). Mixed nuclei or fields are *bad* (red): no shared model can
describe them. The verdict is one line under the panels — dim when the series
is alike, amber when it is not (`⚠ processed differently: LB 0 / 100 Hz (EM) ·
acquired differently: D1 varies 12.5–36 s`); each deviating panel's name ends
with an amber ⚠ whose tooltip says what differs from the majority, the results
table gains a sortable **comparability** column (✓, the difference, or
*reprocessed*) and the spotlight status line repeats the reason. **Details…**
opens the full parameter × spectrum table (majority bold, deviants amber,
every processing cell tooltipped with that spectrum's TopSpin command line,
*only differences* on by default) with **Copy as text** and **Save CSV…**; the
same table comes from `larmor compare <spectra…> [--all] [--csv out.csv]`, and
`larmor batchfit` / `larmor seqfit` print one `warning:` line for a differing
series. Apodisation is not a detail for a shared-width model: an EM of 100 Hz
adds 0.5 ppm of Lorentzian width at 192 MHz, so a member apodised differently
is fitted with the wrong width by construction. **Reprocess all from fid…**
removes that: every member with a raw fid is rebuilt with one common chain —
TDeff (in TopSpin's real points), the window (none, EM, GM or QSINE with LB /
GB), zero-fill to SI, an FT referenced by each spectrum's *own* SF, then each
spectrum's own TopSpin phase (PHC0 / PHC1) or **Autophase** — pre-filled from
the series majority or from one chosen spectrum, and editable. The chain is
recorded on every saved fit as its `processing` with `processing_from_raw` and
the EXPNO as source, so **Save individual fits…**, the batch table's automatic
recipes, the Plotting studio and a reopened project replay it from the
instrument fid (the instrument folders are never written to). Members without
a fid keep their TopSpin spectrum and stay flagged; a previous fit is cleared
(**Fit** again) and baselines reset — **Fit baseline…** still applies
afterwards, and it is the tool for what TopSpin's `abs` did, which is reported
but not replayed. **Use TopSpin processing** reverts to the 1r files.
Acquisition differences (D1, NS, pulse) cannot be reprocessed away; a member
that still differs from the majority carries the fact as a note in its recipe.
The glass protocol (*Glass fitting*, §2) asks for a Gaussian window (GM) rather
than EM — the reprocess form is the place to apply it to the whole series at
once.

**Which EXPNO is the spectrum?** For a month of instrument folders,
**Tools ▸ Session inventory…** reads every EXPNO into a sample × nucleus grid,
pre-picks the production spectrum of each block (the highest EXPNO with a
`pdata/1/1r`, demoted for a short NS or a setup / failed title) and
**Batch fit picks…** hands the picks over in sample order, named by sample
folder (`Base0Ca` … `Base4Ca`), so the panels, the table scopes and the
auto-named recipes carry the sample and not a shared title. Two spectra of one
glass never share a scope: the folder's date token (or the EXPNO) is appended
— `P5-Bi8-12 (04272026)` and `P5-Bi8-12 (05082026)`. The *Getting started*
manual, §2, describes the rule and the flags.

1. **One model, applied.** The batch uses a single model for all spectra — your
   current fit, or a recipe you load in the dialog. The recipe is treated as the
   **answer for lineshape**: every parameter is **held fixed at its recipe value
   except the amplitude**, which is always free per spectrum (and may fall to
   **zero** where a line is absent). This holds regardless of the recipe's own
   pin/vary flags. To let a parameter adapt across the series, tick it under
   **Release per spectrum** (see step 5) — nothing else moves.
2. **See them all — as spectra.** The spectra show in a **3×3 grid** with tabs
   (page through 10–15 at a time). Each cell is a real NMR plot: **sample name**
   top-left, ppm running **high→low**, and you can **drag to zoom** (right-click ▸
   *View All* to reset) exactly like the main window. Toggles above the buttons:
   **components** overlays each site's curve on every fit; **shared scale** puts
   all the plots on one common x/y range for honest comparison (off = each
   auto-scales). Each cell reports its RMSD, updating live after the fit.
   **…and as a table.** A results table sits under the grid, one row per
   spectrum: sample, S/N, RMSD and — after a fit — each site's amplitude, any
   released parameter and its integrated population % (the same numbers
   **Save table…** / **Export CSV…** write; an excluded component is blank).
   **Click a row** to spotlight that spectrum in the grid: its frame takes the
   accent colour, the other cells dim, and its page comes forward. **Click a
   spectrum** to find its row; the arrow keys then step along the series with
   the spotlight following. **Sort any column** (click its header) — RMSD
   descending puts the worst fit on top, so an outlier in a series is one
   click from the data behind it; a flagged spectrum (RMSD outlier or low
   S/N) shows its RMSD in red with the reason as a tooltip. **Esc** clears
   the spotlight; untick **table** to hide it, or drag the divider to trade
   grid height for table height.
3. **Name and order the series.** Each panel is titled by its **sample** as
   the loader derives it (*Getting started*, §2): the sample folder without
   its date, rotor and operator tokens (`01192026_SR31649_Base0Ca_SS_ALP` →
   `Base0Ca`), the title's `Sample …` line for an EXPNO-per-sample layout,
   and the folder's date token when one glass was measured twice
   (`P5-Bi8-12 (04272026)`). The TopSpin title line is not the name — the
   five 2026-01 glasses all read `11B with short tip angle` — it stays in the
   panel tooltip and in the `title` column. **Series table…** opens the one
   editable table behind the series: `#` | `name` | `group` | `folder` |
   `title` | one column per numeric metadata column. The Explorer hands
   spectra over in tree order, which is alphabetical by folder and therefore
   by rotor (`SR31648` sorts before `SR31649`: 0Ca, 3Ca, 1Ca, 2Ca, 4Ca); the
   ▲ ▼ arrows and **Sort by…** (natural name order — `0Ca, 1Ca, …, 10Ca` — or
   any numeric column) set the series order, and OK applies it everywhere at
   once: the grid, the results table, the CSV rows, the saved-recipe names
   and the Series plot. A fitted result follows its spectra (paired by
   source path), so reordering never asks for a refit. Names are kept
   **unique** — a typed collision gets ` (2)` — because a name is a CSV scope
   and a recipe file stem; replicates are expressed through the **group**
   column instead (two folders of one glass share the base name as their
   group, which is what *average replicates* in the Series plot collapses).
   **Join CSV…** brings composition columns in: a wide table with one row
   per sample (an EPMA export — delimiter and BOM are detected, the sample
   column is found by name), or another batch's long `batch_table*.csv`,
   recognised by its `scope, site, param, value` header and pivoted to one
   `<site> <parameter>` column per scope (a 31P batch joined to a 27Al
   series, for instance). The proposed row mapping is shown first — exact
   name or group, then a normalised key (`P5-Bi8-12` = `P5Bi8-12`), an amber
   `(none)` where nothing matched, a blank choice with the candidates listed
   where several rows match loosely — with a checklist of the numeric
   columns (a `_sd` / `_err` / `_u` / ` ±` partner is paired as the column's
   ± automatically; oxide columns are ticked by default) and a tag,
   *analysed* / *nominal* / none, that prefixes the labels. Nothing is joined
   on Cancel. **Add column…** types a column by hand; **Save table…** /
   **Load table…** keep the whole table as `<name>.series.json` (rows pair by
   source path, else by name, and the saved order is followed). **Save
   table…** and **Export CSV…** (step 7) also write `<stem>_series.csv`
   beside the long CSV — position, name, group, folder, title, source path,
   every column and the RMSD — the wide table a notebook otherwise rebuilds
   by hand. The batch session and the project bundle carry the table; the
   Sequential fit (§8) shares the same dialog.
4. **Baseline, per spectrum.** **Fit baseline…** estimates and subtracts a
   baseline from every spectrum *independently* before fitting — **Polynomial**
   (robust asymmetric, choose the order), **Iterative** (Yon 2020), or a flat
   edge-median level. **Reset** restores the raw spectra. For a spectrum that
   needs its own manual correction, **right-click its cell ▸ Add 2-point linear
   baseline** — click two points (one each side of the peaks); each is
   **draggable** with a live preview of the line, so a bad click is fixed by
   dragging rather than starting over. Nothing is applied until you
   **right-click again** to **Apply this baseline** (or **Cancel**) — placing
   the second point never silently commits. The baseline menu is **always
   there on right-click**, including right after applying one, so you can add
   another correction on top (they compose) or clear it. **Clear this
   spectrum's baseline** (same menu) restores just that one spectrum to raw.
   This correction is **recorded on the spectrum's own fit**
   (its `.recipe.json` carries the two points and the source file), so
   **Save individual fits…** exports it faithfully and reopening that fit later
   reproduces the exact corrected spectrum — re-running the global **Fit
   baseline…** afterward recomputes from raw and clears any per-spectrum
   corrections layered on top.
   The same right-click menu has **Exclude component ▸**, for a line that
   only belongs in *some* of the spectra (e.g. a Bi-contact line only real
   for Bi-loaded glasses): pick it for the spectra where it doesn't apply,
   and its amplitude is **locked to exactly zero for that spectrum only**,
   instead of fit. An excluded component never draws, never gets a legend entry, and is
   **left out of the exported table/CSV and any plot built from it**
   entirely — not reported as "a fitted zero". The panel title grows an
   "(excluded: …)" note so it's never mistaken for a fit that simply found
   nothing there. An excluded component's position/width/shape are also held
   fixed for that spectrum, even if "Release per spectrum" is ticked for
   them elsewhere — a line that isn't there has nothing to release, and
   letting it drift would only hand the fit useless free parameters.
5. **One Fit button; choose what may move.** **Fit** refines **only the
   amplitudes** per spectrum — everything else is held at the recipe. Whichever
   parameters you tick under **Release per spectrum** are additionally fit,
   **independently per spectrum**, allowed to drift by **±X %** around their
   recipe value (a relaxation — e.g. let δ_iso wander ±5 % across the series while
   widths stay pinned). You choose **parameter by parameter**; anything unticked
   does not move. A **completion threshold** (Δσ %) stops the fit once the
   residual stdev stops improving. Interrupt any time with **Cancel** (discard,
   revert) or **Stop** (keep the latest iteration) — the same two modes as the
   main fitter.
6. **Error calculation.** After the fit, choose how the per-spectrum errors are
   estimated from the **Error calculation** menu, then **Compute errors**:
   * **Covariance** — the least-squares covariance stderr. The batch fit's own
     pass skips the (potentially costly) errorbar-rescue step for speed, so
     **Compute errors** here does one quick confirming re-fit per spectrum to
     get real numbers, rather than reporting nothing when that first pass
     couldn't get a clean covariance (common for several overlapping,
     correlated released parameters) — usually still fast, since each
     spectrum starts from its already-converged values.
   * **Monte-Carlo** — refit *N* synthetic noisy copies of each spectrum and take
     the spread; captures correlations and non-linearity the covariance misses.
   * **χ² profile (error analysis)** — scan each fitted parameter, refit the rest,
     and read a real 1σ interval off the χ² curve.
   These are the same estimators as the single-fit **Errors** tools, run for
   every spectrum. **Export CSV…** writes one row per fitted parameter — value,
   error, %-error, and (for the χ² profile) the 1σ interval — tagged with the
   selected method (it computes that method first if you have not yet). Switching
   the menu never loses a method you already computed. Both this CSV and
   **Save table…** end with `vary`, `min`, `max`, `expr`, `at_bound`: a shared
   row reads `vary = False` because the batch holds it, and a released value
   that stopped at the edge of its ± range reads `min` / `max` there — the
   same cell shows ‡ in the table under the grid, with the remedy in its
   tooltip (widen the release %, or fit the series sequentially).
   Monte-Carlo and χ² profile are each hundreds to thousands of independent
   refits (every trial, or every scan point of every released parameter of
   every spectrum), so both run across all of your CPU cores (one left free
   for the app itself) instead of one refit at a time — the same numbers,
   much less waiting. **Stop** still works: it finishes whatever's already in
   flight on each core rather than cutting off instantly, the same "keeps
   the current work" behaviour as everywhere else in LARMOR that can be
   interrupted.
7. **Save.** **Save individual fits…** writes one LARMOR `.recipe.json` per
   spectrum, named **automatically**
   (`sample_nucleus_recipe_batch_YYYYMMDD_HHMM`) or with a **name you type for
   each**
   (it prompts in turn, showing the sample and proc number) — each carries the
   errors from the last error-calculation you ran. **Save table…**
   writes a `batch_table.csv` of the shared and per-spectrum values — the
   same numbers as the on-screen table, in long form and with the shared
   parameters, each row with its `vary` / `min` / `max` / `expr` /
   `at_bound` status — **plus each site's integrated population %** (a
   `population_pct` row per site,
   same integral-over-the-window quantification as Report/§5) — the exact
   column the Plotting studio's species-distribution chart wants, without a
   separate export step. An excluded component (above) is left out entirely,
   not reported at 0%. Lines tagged with a **family** in the model (the
   family column of the Fit-parameters table) add a `family_pct` row per
   family (site ids `f0`, `f1`, …; label = the family) and a `ratio` row per
   named ratio of the nucleus (`r0`, …; label `N4`, `⟨CN⟩ Al`, …) in the same
   long format — the headers do not change. Their error is the exact
   per-spectrum value of the last error calculation when one was run
   (Monte-Carlo: the spread of the per-trial re-integrated sums;
   covariance: propagated over the amplitude covariance), else the flagged
   *independent* estimate; the results table under the grid shows them as
   trailing `Σ BO4 / family %` and `N4 / ratio` columns in a neutral colour
   (blank where every line of the family is excluded), and the tooltip
   names the basis. The "also save individual fits… next to the CSV"
   checkbox (on by default) auto-writes each spectrum's `.recipe.json`
   alongside the table too, so the Plotting studio's batch-grid finds the
   real saved fits automatically (bounds, `vary`, baseline processing
   included) instead of only having the CSV's bare values to work from.
   Every saved individual fit carries the source kind, the SHA-256 of its
   data file and the acquisition block of its spectrum. **Acquisition
   table…** needs no fit: it opens the same window as **Tools ▸ Experimental
   section…** over the loaded spectra — Table S1 with every parameter that
   varies across the series highlighted, the Experimental paragraph with
   those parameters as ranges, **Copy paragraph** / **Copy LaTeX table** /
   **Save CSV + LaTeX…**.
   **Publication bundle…** writes everything about the batch to one folder
   of your choice: `batch_table.csv` (plus the error table for the selected
   error-calculation method when it has been computed), one `.recipe.json`
   and one `_curves.csv` per spectrum — named `01_<sample>`, `02_<sample>`…
   in series order — `manifest.csv` and `README.txt`. A `_curves.csv` carries
   the experiment **exactly as fitted** (after the recipe's processing steps
   and any per-spectrum baseline; when a baseline was applied the pre-baseline
   trace is included as `experiment_raw`), the total model, the residual and
   one `s<i>_<label>` column per component, all on the experimental ppm axis,
   and reopens in LARMOR as a spectrum through **File ▸ Open**. An excluded
   component stays as a zero column, named in the file's header and in
   `manifest.csv`, so every file of a series has the same layout.
   `manifest.csv` records per spectrum the source file and its SHA-256, EXPNO
   and procno, NS, D1, BF1/SF/SR, the processing steps, the fit window, the
   RMSD (normalised by the window maximum, as everywhere in LARMOR), the
   excluded sites, the error method and the software versions; `README.txt`
   carries the Methods paragraph, the versions, the conventions (RMSD,
   populations, number precision) and a file glossary. Read the CSVs with
   `pandas.read_csv(path, comment="#")`. If the chosen folder already holds a
   `manifest.csv`, the tool asks before replacing it. The same bundle comes
   from the command line with `larmor batchfit … --curves`.
   **Series plot…** charts how any parameter (δ_iso, width, C_Q, η, or population %)
   evolves along the series. With tagged lines the list also offers **Σ BO4**-style
   family items (drawn on the population subplot) and the named ratios (their own
   *named ratio* subplot), each with the error bars of the chosen method and
   exported as `family:BO4` / `ratio:N4` columns. Its **Error bars** menu chooses which computed error
   to draw and export — *covariance*, *Monte-Carlo*, or *χ² profile* (whichever
   you ran in step 6), or *none*. The **integrated population %** carries an
   error too — first-order from the amplitude's error under the chosen method
   (the other sites' amplitude errors, which also shift the total, are
   neglected — the same approximation the Report table uses). Export the
   numbers (the ± column is labelled with the chosen method) or the figure, and
   **Send to Plotting studio** carries the points *and their error bars* into
   the studio, where the axes, limits,
   ticks, legend and fonts are fully customisable. The **x axis** menu
   replaces the series order by any numeric column of the Series table
   (step 3): points sit at the column's value, sorted by it, its ± becomes
   x error bars and the column label the axis title — population against
   analysed P2O5, or a 27Al shift against the 31P bonded-P fraction joined
   from another batch. **Average replicates** (shown when a group has more
   than one member, on by default then) collapses each group to its mean;
   the bar is the sample standard deviation (ddof = 1) — a singleton keeps
   its own fit error — and the exported CSV carries `n`, `x` and `x ±`.
   **Fit line** adds an ordinary-least-squares line through three or more
   points, labelled with the slope ± its error, Pearson r and n; it travels
   with the figure export and into the studio. **Species bar…** opens the
   Plotting studio on a 100 %-stacked bar of every site's integral
   population, one bar per name (or per group), in the current x order.
   **Export DUST CSV…** writes a composition file DUST imports directly:
   `Sample`, the Series table's oxide columns under DUST's canonical names
   (`P2O5_mol`, `Bi2O3 (mol%)` → `P2O5`, `Bi2O3`; non-oxide columns are
   skipped and named), then `N4_measured` and `N4_measured_err` — the
   selected sites' summed integral population as a fraction, with its
   error. DUST's import dialog maps the two N4 columns to *Ignore*; join
   them against DUST's Results CSV by `Sample` to compare measured and
   predicted N4. (For a fuller publication table
   across independent fits, see the **Batch fit report** tool.)

It builds on the same co-fit engine (§3), so the shared parameters carry full
uncertainties. The completion threshold is global (set it under **Decomposition ▸
Advanced ▸ Fit completion threshold**) and honoured by every fit in LARMOR.

## 7 · Plotting studio — publication figures from a batch

The **Plotting studio** (any *Send to Plotting studio*, or *Tools ▸ Plotting
studio*) builds a figure as a plain, reloadable spec — style, labels, ticks,
legend and size are shared across every plot kind. Two kinds are purpose-built
for a **whole series at once**, so a batch fit (§6) becomes a submission-ready
figure without hand-assembling panels:

1. **Templates.** The **Template** picker at the top is a set of named,
   nucleus-generic starting points copied from common published NMR figure
   styles — *Stacked series*, *Deconvolution grid*, *Composition series
   (shaded component)*, *Composition trend*, *Species distribution*, *2D
   correlation*. Picking one sets the plot kind and sensible layout/style
   defaults; every field underneath is still yours to change, and it combines
   freely with any journal **Style** (Nature, ACS, RSC, …).
2. **Batch grid** — one panel per spectrum, small-multiples: experiment + total
   fit + components, laid out in a grid. **Load CSV…** takes a
   `batch_table*.csv` (Batch fit's **Save table…** / **Export CSV…**, §6) and
   auto-matches its rows to the saved `.recipe.json` fits — by sample name next
   to the CSV, and by each row's own **`source_path`** column (every CSV
   exported from Batch fit now carries it, so the studio finds the right
   spectrum even without **Save individual fits…** too). **Load folder…**
   takes a folder of saved fits directly. Check/reorder/remove panels in the
   list; **Components** chooses fill / dashed-outline / total-only, **Shade
   only** highlights one component per panel (the "composition series"
   style), and **Peak labels** adds position, letter, or position + integrated
   population %.
   - **Older CSVs and missing files.** A CSV written before this feature (no
     `source_path` column), or any sample neither method above can place, is
     never just dropped — it's flagged **⚠ locate data…** in the panel list,
     and the studio asks directly (one file dialog per sample, at load time)
     for that spectrum's data — pick its dmfit fit, or its EXPNO/`pdata`
     folder (a Bruker `1r`). Double-click the row any time afterward to try
     again if you cancelled, and double-click a resolved row to rename its
     panel title.
   - **CSV-only reconstruction.** When no saved `.recipe.json` matches, the
     studio rebuilds a full fit straight from the CSV's own rows (needs the
     `model` column, on every export since §6) — a site the CSV never gave
     every parameter for (an unreleased/held one, or one a model defaults
     when omitted) fills in from that model's own registry default, the same
     value a freshly-added site of it would start from. Works across mixed
     models in one recipe (e.g. a Gauss/Lorentz line next to a Czjzek site)
     — each site keeps its own model and parameter set. An excluded
     component (§6) is simply absent from the reconstructed recipe, exactly
     as if it were never part of that spectrum's model.
   - **Component colors / legend…** — a color swatch and an "in legend"
     checkbox per detected component (from the first resolved panel's fit).
     **Hide** (next to **Shade only**) drops a component entirely — no line,
     fill, or legend entry, in every panel — while an unchecked legend box
     keeps the line but drops just its label (for a component that's obvious
     from position/color and doesn't need one competing for space).
3. **Species distribution** — a 100%-stacked bar of species/oxygen population
   vs. composition. Type the category × species table directly, or **Load
   from batch CSV…** to pivot one parameter (e.g. `amplitude`) out of a
   `batch_table*.csv` automatically, one row per sample — each bar normalizes
   to 100% on its own, so raw amplitudes work without pre-converting to %.
   From the Series plot, **Species bar…** builds the same chart directly,
   its categories following the Series table's names and order.
   From a batch whose model carried family tags (§6), pick `family_pct` for
   a bar stacked by family (Σ BO3 / Σ BO4 per glass) instead of by line;
   `ratio` also appears in that picker but is not a stackable quantity —
   N4 per glass belongs in a *Composition trend* or the Series plot.
4. **Auto update / Preview.** Auto update is **off by default** — a batch
   grid with many panels (each a full reconstruction + population-%
   integral) can be slow to redo on every tweak. **Preview** renders on
   demand regardless of the toggle; turning Auto update on immediately
   re-renders once so it never shows a stale preview.
5. **2D publication figures** (the "2D contour" kind). **Nucleus** and
   **Larmor (MHz)** drive axis labels *and* the computed reference lines
   below. **Fit overlay** takes a saved 2D fit (`Decomposition ▸ Fit` on an
   MQMAS map, then save the recipe) and draws it as a dashed contour over the
   experimental one, at the **MQMAS method** you fit it with (not stored on
   the recipe itself — pick the one you used). **Add iso/quad line…** now
   offers **Compute** for the two lines an MQMAS figure actually needs: the
   **CS axis** (the diagonal a pure-chemical-shift site would sit on) and the
   **QIS axis** (the direction a site moves as C_Q grows, drawn from its own
   δ_iso) — both use `larmor.twod`'s own physics for the chosen
   nucleus/method rather than a hand-typed slope, and stay fully editable
   afterward.

## 8 · Sequential fit — forward / backward series sweep (1D)

The batch tool (§6) assumes one *shared* model. Some series don't work that way:
the lineshape **evolves smoothly** from one end-member to the other (a
composition or temperature series), and each spectrum deserves its own fit — just
one that starts from where its neighbour ended. **Tools ▸ Sequential fit** does
exactly that. Ctrl/Shift-select the series in the Explorer in any order —
**Series table…** sets the sweep order and the names (§6, step 3) — open it,
and you get a **one-spectrum-at-a-time** workbench:

1. **Precise, per-spectrum control.** The current spectrum shows with its model
   and components, and its **full fit-parameters table** — set values, bounds,
   fixes and links exactly as in the main window. **Fit current** fits just this
   one.
2. **Carry it forward.** **◀ Prev / Next ▶** move along the series; when you move,
   the spectrum you land on is **seeded from the one you left** (tick which
   parameters carry — positions/widths/quadrupolar by default, amplitudes always
   re-fit fresh). **Fit → seed next ▶** fits the current spectrum and steps on.
   This is the manual forward (or backward) chain.
3. **Automate it.** **Auto ⇄ forward–backward fit** runs the whole sweep itself:
   choose the number of **passes** (1, 2, 4, 8, 16 — each pass sweeps one
   direction, so 2 = forward then back), which end to **start** from, and an
   optional **smoothing** window that gently smooths each parameter's trajectory
   *between* passes so the series doesn't jitter. A live plot shows the **RMSD of
   every spectrum** updating and the **mean RMSD falling** pass over pass, plus a
   **trajectory plot** of any chosen parameter across the series. **Cancel**
   reverts; **Stop** keeps what's done.
4. **Save.** **Save individual fits…** (auto `sample_nucleus_seq_YYYYMMDD_HHMM` or
   a name per fit) and **Series plot…** (parameter/population evolution, with
   export) — as in the batch tool; every saved recipe carries its source
   path, the SHA-256 of its data file and its acquisition block, and
   **Acquisition table…** opens the Experimental-section window over the
   series (no fit needed). **Publication bundle…** is the batch tool's
   bundle (§6, step 7) for the series — `seq_table.csv` instead of
   `batch_table.csv`, otherwise the same files — and works after manual **Fit
   current** steps as well as after an auto sweep (a member never fitted gets a
   manifest row marked *not fitted*); the saved recipes carry their source
   path. `larmor seqfit … --curves` writes the same from the command line.
   The comparability line of §6 (with **Details…**, without the reprocess —
   every member gets its own model here) sits under **◀ Prev / Next ▶**, and
   the current spectrum's title carries the amber ⚠ when it was acquired or
   processed unlike the rest of the series.

Use §6 when the sites are genuinely the *same* everywhere and only populations
change; use §8 when the sites themselves **evolve** along the series.

## References

- D. Massiot *et al.*, "Modelling one- and two-dimensional solid-state NMR
  spectra", *Magn. Reson. Chem.* **40**, 70 (2002). *(joint 1D/2D fitting; the
  quadrupolar product as the field-independent invariant)*
- M. Newville *et al.*, **lmfit** (constrained/linked least squares with
  uncertainties), doi:10.5281/zenodo.11813 (2014).
- G. Czjzek *et al.*, *Phys. Rev. B* **23**, 2513 (1981); J.-B. d'Espinose de
  Lacaillerie, C. Fretigny, D. Massiot, *J. Magn. Reson.* **192**, 244 (2008).
  *(the Czjzek widths shared in a glass co-fit)*

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
