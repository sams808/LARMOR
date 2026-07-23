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
