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
**stack offset**. Promote an overlay to **active** (it becomes the fit target)
and the previous active demotes back to an overlay. The active spectrum is always
the single object the 1D fitter works on, so overlays never disturb a fit.

**File ▸ Save project…** captures every open 1D workspace — spectrum,
processing, fit, *and its overlays* — as one reopenable `.larproj.json`
file; **File ▸ Open project…** restores it all, relocating each overlay
from its saved source path (an overlay whose file has moved is reported,
not silently dropped). 2D workspaces are not yet included in a project —
reopen those from their source.

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
     error** — one row per site, machine-readable.
   - **`table.tex`** — the same table as a LaTeX `tabular` for direct inclusion.
   - **`report.md`** — a Markdown report: the table, the nucleus/field summary,
     the error method, and a **per-fit overlay** (experiment + model +
     components + residual).
   - **`figures/*.png`** — the individual overlays.

The columns are model-aware (Czjzek sites report σ and the field-independent
`√⟨P_Q²⟩`; discrete/Amorphous sites report C_Q and η) and populations come from
the same integral-over-the-window quantification as **Report** (F6). It is the
fastest route from "a folder of fits" to a table you can paste into a manuscript.

## 6 · Batch fit — one shared model, many spectra (1D)

For a **series measured the same way** (a composition series, a time course),
you often want *one* set of lineshape/position parameters describing every
spectrum, with only the **amplitudes** free — the sites are the same, their
populations change. **Ctrl/Shift-click** the spectra in the Explorer and press
**Batch fit selected…** (or *Tools ▸ Batch fit spectra*).

1. **One model, applied.** The batch uses a single model for all spectra — your
   current fit, or a recipe you load in the dialog. The recipe is treated as the
   **answer for lineshape**: every parameter is **held fixed at its recipe value
   except the amplitude**, which is always free per spectrum (and may fall to
   **zero** where a line is absent). This holds regardless of the recipe's own
   pin/vary flags. To let a parameter adapt across the series, tick it under
   **Release per spectrum** (see step 4) — nothing else moves.
2. **See them all — as spectra.** The spectra show in a **3×3 grid** with tabs
   (page through 10–15 at a time). Each cell is a real NMR plot: **sample name**
   top-left, ppm running **high→low**, and you can **drag to zoom** (right-click ▸
   *View All* to reset) exactly like the main window. Toggles above the buttons:
   **components** overlays each site's curve on every fit; **shared scale** puts
   all the plots on one common x/y range for honest comparison (off = each
   auto-scales). Each cell reports its RMSD, updating live after the fit.
3. **Baseline, per spectrum.** **Fit baseline…** estimates and subtracts a
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
4. **One Fit button; choose what may move.** **Fit** refines **only the
   amplitudes** per spectrum — everything else is held at the recipe. Whichever
   parameters you tick under **Release per spectrum** are additionally fit,
   **independently per spectrum**, allowed to drift by **±X %** around their
   recipe value (a relaxation — e.g. let δ_iso wander ±5 % across the series while
   widths stay pinned). You choose **parameter by parameter**; anything unticked
   does not move. A **completion threshold** (Δσ %) stops the fit once the
   residual stdev stops improving. Interrupt any time with **Cancel** (discard,
   revert) or **Stop** (keep the latest iteration) — the same two modes as the
   main fitter.
5. **Error calculation.** After the fit, choose how the per-spectrum errors are
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
   the menu never loses a method you already computed.
   Monte-Carlo and χ² profile are each hundreds to thousands of independent
   refits (every trial, or every scan point of every released parameter of
   every spectrum), so both run across all of your CPU cores (one left free
   for the app itself) instead of one refit at a time — the same numbers,
   much less waiting. **Stop** still works: it finishes whatever's already in
   flight on each core rather than cutting off instantly, the same "keeps
   the current work" behaviour as everywhere else in LARMOR that can be
   interrupted.
6. **Save.** **Save individual fits…** writes one LARMOR `.recipe.json` per
   spectrum, named **automatically**
   (`sample_nucleus_recipe_batch_YYYYMMDD_HHMM`) or with a **name you type for
   each**
   (it prompts in turn, showing the sample and proc number) — each carries the
   errors from the last error-calculation you ran. **Save table…**
   writes a `batch_table.csv` of the shared and per-spectrum values, **plus
   each site's integrated population %** (a `population_pct` row per site,
   same integral-over-the-window quantification as Report/§5) — the exact
   column the Plotting studio's species-distribution chart wants, without a
   separate export step. An excluded component (above) is left out entirely,
   not reported at 0%. The "also save individual fits… next to the CSV"
   checkbox (on by default) auto-writes each spectrum's `.recipe.json`
   alongside the table too, so the Plotting studio's batch-grid finds the
   real saved fits automatically (bounds, `vary`, baseline processing
   included) instead of only having the CSV's bare values to work from.
   **Series plot…** charts how any parameter (δ_iso, width, C_Q, η, or population %)
   evolves along the series. Its **Error bars** menu chooses which computed error
   to draw and export — *covariance*, *Monte-Carlo*, or *χ² profile* (whichever
   you ran in step 5), or *none*. The **integrated population %** carries an
   error too — first-order from the amplitude's error under the chosen method
   (the other sites' amplitude errors, which also shift the total, are
   neglected — the same approximation the Report table uses). Export the
   numbers (the ± column is labelled with the chosen method) or the figure, and
   **Send to Plotting studio** carries the points *and their error bars* into
   the studio, where the axes, limits,
   ticks, legend and fonts are fully customisable. (For a fuller publication table
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
exactly that. Ctrl/Shift-select the series in the Explorer (in order), open it,
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
   export) — as in the batch tool.

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
