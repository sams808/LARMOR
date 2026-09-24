# Tutorial 5 — MQMAS: opening, checking and fitting a 2D map

*Time: ~25 minutes. Everything used here ships with the repository:
`examples/pCABS2-4/3620` is the ²⁷Al 3QMAS map of the glass fitted in
Tutorial 1 (TopSpin-processed `2rr` with its `2ri`/`2ir`/`2ii` quadrants and
the raw `ser`), and `examples/pCABS2-4/3616` is the quantitative single-pulse
spectrum of the same sample. Run the commands from the repository root. The
first 2D fit builds the MQMAS kernel, which takes noticeably longer than a
1D kernel.*

Tutorial 1 ended on a degenerate 1D fit: three overlapping Czjzek sites whose
minor components the 1D lineshape could not pin down. An MQMAS spectrum
separates the sites along a second, isotropic dimension. This tutorial opens
the map, checks that it is sheared and phased, places the three sites on their
ridges, fits them in 2D, and finally co-fits the 2D map with the quantitative
1D spectrum.

## 1. Look before you fit

```
larmor info examples/pCABS2-4/3620
```

```
EXPNO: examples\pCABS2-4\3620
nucleus: 27Al   SFO1: 130.3230 MHz
pulse program: mp3qdfsz   TD: 1024   SW: 119048 Hz
MASR (acqus): 26000.0 Hz
title: 02/12/21 : conditions OK
```
<!-- measured v0.12.1, 2026-09-22: `larmor info examples/pCABS2-4/3620` -->

`mp3qdfsz` is Bruker's split-t₁ 3QMAS sequence with a z-filter; TD and SW are
those of the direct dimension. The title continues with `2D MAS=26 kHz`, which
agrees with the acqus rate, so no `CONFLICT` line is printed (compare
Tutorial 4, where the two sources disagree). The command line describes a 2D
EXPNO but fits only 1D data and dmfit files; the 2D work below happens in the
desktop application.

## 2. Open the map

```
larmor desktop
```

**File > Open EXPNO / folder…** (Ctrl+Shift+O) and pick
`examples/pCABS2-4/3620` — or **File > Open…** (Ctrl+O) on its
`pdata/1/2rr`. A processed 2D is not fittable as a 1D spectrum, so LARMOR
routes it to the contour view, titled `2D spectrum — 3620`, with the status
bar reading:

```
MQMAS 2D displayed — Fit ▸ MQMAS to fit · Tools ▸ Relaxation for a series · or send a 1D trace to fitting below
```

Conventions: F2 (horizontal) is the MAS dimension in ppm, high ppm to the
left as usual; F1 (vertical) is the isotropic dimension. The control bar
above the map holds, from left to right: the `contours` sign (positive,
negative, both), the `display` mode (contour, density, filled,
contour+values), the colour map, the number of `levels`, the contour
`floor ×σ`, the `shear` factor with **Apply shear**, **Phase 2D**,
**Calibrate**, **Measure**, **CS/QIS axes**, **Overlay 1D ▾**, **HMQC** and
**2D ops ▾**.

Press **CS/QIS axes** first: it draws the chemical-shift (CS) and
quadrupole-induced-shift (QIS) reference lines for ²⁷Al at this field and
method. Sites of one coordination lie along a QIS line; their spread along
the CS line is the distribution of isotropic shifts.

## 3. Is the map sheared?

The shipped `2rr` was sheared in TopSpin, so the three ridges run horizontal
(constant F1). The `shear` spin box and **Apply shear** exist for data
processed without the shear step: a positive or negative factor tilts the map
until the ridges are horizontal (the *MQMAS (2D)* manual, §1, gives the
²⁷Al 3Q value). Leave it at 0 here — applying a shear to an already sheared
map tilts the ridges the other way.

## 4. Is the map phased?

TopSpin's `2rr` is normally phased. When it is not, **Phase 2D** opens the
phasing bar: choose `rows (F2)` or `columns (F1)`, click one or two reference
peaks on the map, press **Phase ▶**, adjust `p0`/`p1` and **Apply**. Because
the four quadrants ship, the correction is exact, not an approximation from
the real part alone (the *2D processing* manual, §3).

## 5. Is the projection representative?

MQMAS intensities are weighted by the triple-quantum excitation and
conversion efficiency, which depends on each site's quadrupolar coupling: the
F2 projection of a 3QMAS map is not quantitative. Check how far it deviates
before reading any population off the 2D:

1. **Overlay 1D ▾** → `From file → F2 (direct)…` → `examples/pCABS2-4/3616`.
   The single-pulse spectrum is drawn over the F2 projection.
2. The **HMQC** toggle shows a second bar with a `scale` box and **fit scale**
   (a least-squares scale of the projection to the 1D) when the two need to
   be brought to a common height.

Whatever the AlO₅ and AlO₆ shoulders look like in the projection compared
with the 1D, quantify from `3616`, not from the map (§10).

## 6. Place and fit the three sites

1. **Edit > Add line > Czjzek (quad. distribution)**. The action
   stays checked: every click on the map now adds a site.
2. Click the map once on each ridge: AlO₄ (F1 near the 1D position of the
   main site), AlO₅, then AlO₆. Each click seeds the site's isotropic shift
   from the clicked F1 position and reports:

   ```
   added czjzek at F1≈63.8 ppm — click to add more, or click czjzek again (Esc) to stop; Fit ▸ Fit to fit the 2D
   ```

   (the number is the F1 value under the cursor). Press Esc, or click the
   model entry again, to leave add mode. The new sites appear in the *Fit
   parameters* dock as `Czjzek-0`, `Czjzek-1`, `Czjzek-2`; rename them to
   AlO₄ / AlO₅ / AlO₆ there.
3. **Fit** (F5). The status bar reads
   `fitting the 2D (building the MQMAS kernel, first time is slow)…`; the
   kernel is a grid of 2D lineshapes over C_Q and η that is built once per
   nucleus, field and window and then cached. Only the models with a 2D
   implementation can take part: Czjzek, ext. Czjzek, Quad CT (2nd order) and
   Quad CT + CSA — a fit containing any other line is refused with a message
   listing them.
4. When the fit finishes the results strip reads
   `2D MQMAS fit · RMSD … · F1 ref +… ppm` and the *Report* dock starts with
   `MQMAS F1 isotropic-axis reference offset: … ppm (auto-fitted)`, followed
   by the usual lmfit report with a standard error per parameter. The fitted
   contours are drawn over the data in the colours of the table; toggle
   `contours` between positive and both to see the residual sign.

No fitted values are quoted here: the 2D fit was not re-run for this
tutorial. Two references exist for comparison. Tutorial 1's 1D fit gave
δiso = 63.8 / 30.9 / −2.0 ppm for AlO₄ / AlO₅ / AlO₆ (with a ±24 ppm error on
the last). The validation report (`docs/validation.md`) records that a fit of
this map reproduced the dmfit MQMAS analysis to about 0.5 ppm:
62.2 / 29.7 / −1.1 ppm against dmfit's 62.7 / 30 / −0.35 ppm. A fit started
from the three ridge clicks should land within a few ppm of those values for
the two strong sites; the AlO₆ position is set far better by the 2D than by
the 1D because the ridge is resolved in F1.

The kernel resolution is under **Fit > Fit settings > Computing
parameters…**: the *MQMAS kernel* block defaults to 192 F2 points, 96 F1
points, 40 C_Q steps, 6 η steps and a C_Q maximum of 16 MHz. Larger values
cost time on every first fit; smaller ones make the ridges blocky.

## 7. The F1 reference

The kernel's F1 axis follows the isotropic-shift convention (chemical shift on
the diagonal); an experimental F1 axis differs by a referencing offset that
depends on how the map was processed. The fit determines this offset as one
extra parameter, which is why the strip reports `F1 ref`. To hold it at a
known value instead (dmfit style), use **Fit > MQMAS > MQMAS F1
reference…**: entering a value and pressing OK fixes it — the status bar reads
`MQMAS F1 reference fixed at … ppm — Fit to apply` — while Cancel returns it to
auto-fitting. A fixed reference is reported as `(held fixed)` in the *Report*
dock.

## 8. Constrain what the 2D cannot resolve

Three Czjzek sites on a 2D map are far better determined than on a 1D
lineshape, but not perfectly: the AlO₅ ridge overlaps both neighbours and its
σ(C_Q) trades against its amplitude. The constraint tools of Tutorial 2 apply
unchanged — open a site's **⚙** card to bound its δiso to the ridge, or fix
σ(C_Q) when the tail is under-determined. The *MQMAS (2D)* manual's
*Constraining overlapping sites* section discusses which parameters a 2D map
determines and which it does not.

## 9. From 2D to 1D and back

Any 1D trace can be pulled out of the map into a new fitting workspace:
**2D ops ▾** → `F2 skyline → fit` or `F2 sum → fit`, or the buttons under
*Send a 1D trace to fitting* (**F2 skyline →**, **F2 sum →**,
**row at cursor →** after dragging the dashed line onto a row). The trace
opens as a new workspace — status
`F2 skyline — new workspace; add lines and Fit` — and the map stays open in
its own. **View > Zoom > Back to 2D map** (Ctrl+2) returns to it; the *Workspaces*
dock lists both. A row through one ridge is the cleanest way to read that
site's C_Q distribution in isolation.

## 10. Co-fit with the quantitative 1D

The 2D map fixes positions and quadrupolar parameters; the single-pulse
spectrum carries the populations. **Fit > Co-fit
datasets…** fits both at once with one physical model:

1. The dialog *Co-fit datasets (shared model)* starts from the current fit's
   sites. Press **＋ Add dataset…** and pick `examples/pCABS2-4/3616`.
2. Under *Tie across datasets*, tick the parameters that must be identical in
   both experiments — the isotropic shift, σ(C_Q) and the shift-distribution
   width are ticked by default; amplitudes are always independent, because
   the 3QMAS efficiency scales them differently.
3. **Run co-fit**. The report lists `shared: …`, one RMSD per dataset (the 2D
   one with its `F1 ref`), and every site's parameters with their standard
   errors. **Preview (no fit)** draws the current model over both datasets
   without fitting.
4. **Apply shared params to current fit** writes the tied values back into
   the open workspace.

> 3QMAS volumes are not populations. Quantify (**Report (quantify)**, F6) on
> the `3616` workspace, where the amplitudes are those of a single-pulse
> experiment; the co-fit has already constrained its positions and widths
> with the information from the 2D.

## 11. Save and draw

**File > Save fit** (Ctrl+S) writes the 2D fit — including the F1
reference and the kernel method — as an ordinary recipe next to the data.
Tutorial 3, §3, shows the 2D figure spec for this dataset (path
`examples/pCABS2-4/3620`, the `3616` overlay on the top projection, per-site
sub-projections); in **Plotting > Plotting studio…** the same spec is built
interactively, and a saved 2D fit recipe can be given as the fit overlay.

## 12. Optional: the raw route

The `ser` file ships too. **File > Open FID…** (Ctrl+F), **Open fid / ser…**,
pick `examples/pCABS2-4/3620/ser`; the F1 mode for this acquisition is
`States-TPPI` (acqu2s FnMODE 4). **Transform preview** shows the processed
map; **Use this spectrum →** opens it in the *2D MQMAS* viewer (levels,
`floor ×σ`, shear and **Apply shear**), which displays but does not fit — the
same viewer as **Fit > MQMAS > 2D MQMAS viewer / fit…**. Fitting is done on a
processed `2rr` in the main window, as above.

## Where to read more

**Help > User manuals**: *MQMAS (2D)* (physics, shear factor, F1 referencing,
constraining sites), *2D processing* (phasing, 2D operations) and
*Multi-dataset & co-fitting*, §3.
