# MQMAS (2D) — display, process, fit

## Open
Open the processed `2rr` (or its pdata folder). It lands on the **contour map**
in its own workspace. F2 is the MAS (direct) dimension, F1 the isotropic
(indirect) one.

## Display
- **contours**: positive / negative / both; **levels** and **floor ×σ**.
- The **projections** (top = F2, left = F1) track the map. **Measure** and
  **Calibrate** work on the 2D as on a 1D.
- Axes follow the standard NMR/dmfit convention: F2 high-ppm **left**, F1
  high-ppm **top** (diagonal runs top-left → bottom-right). **2D ops ▾ ▸ Flip
  F2/F1 axis** flips either direction (view only, data untouched) to match
  another program or a figure.
- **2D ops ▾**: Transpose F1↔F2, Reverse F1/F2, Diagonal → fit, projections →
  fit.

## Superpose a 1D spectrum
**Overlay 1D ▾** puts a 1D spectrum on the F2 (direct) or F1 (indirect)
projection for direct comparison: **Current 1D** uses the spectrum already open,
**From file…** loads any 1D. Tune **scale** / **fit scale** in the overlay bar;
**Clear 1D overlays** removes them.

## Phase (if needed)
**Phase 2D** → click 1 or 2 reference peaks → the selected rows/columns show as
full-width traces with the pivot marked → p0/p1 + ±90/180° → Apply. With the
hypercomplex quadrants present the correction is exact.

## Shear / reference
**Apply shear** for data processed without xfshear. (mrsimulator MQMAS is already
sheared.) Use the F1 axis referencing conventions for your spin/method.

## Isotropic-axis (F1) referencing
mrsimulator's kernel puts a pure-CS site on the diagonal (F1 = δiso), but an
experimental Bruker/dmfit F1 axis is referenced by its own convention — often
tens of ppm away. The fit therefore carries a single **F1 reference offset** that
slides the whole model along F1 to line up with your data (δiso still sets the
diagonal position, so the two don't fight). It is **auto-fitted** and reported
after every 2D fit. To pin it (dmfit-style manual referencing) use
**Decomposition ▸ MQMAS F1 reference…**; re-run with Cancel to return to auto.
If a Czjzek 2D fit looks well-shaped but sits off in F1, this is the knob.

## Fit
1. Pick **Czjzek** (glasses) or **ext. Czjzek / Quad CT** (crystalline) from
   Models, then **click the contour** to place a site (its δiso starts at the
   clicked F1).
2. **Fit** — the MQMAS kernel builds once, then the fitted contours overlay in
   orange with the RMSD and per-site δiso/Cq/η/σ.
3. Or pull the **isotropic projection** to the workbench and fit it as a 1D, or
   **co-fit** it with a 1D MAS spectrum (Decomposition ▸ Co-fit datasets). The
   co-fit window now plots every dataset with its shared-model fit (1D as an
   overlay, 2D as experiment/model contours) beside the numeric report.

### Glass broadening (dmfit CzSimple parity)
The **Czjzek** model separates the two glass broadenings that dmfit calls
`dCS` and `wid`:
- **dCS** (`isotropic-shift distribution FWHM`) — the spread of isotropic
  chemical shifts. In 2D it elongates the peak **along the diagonal** (δiso moves
  F2 and F1 together), the hallmark of a disordered site.
- **line** (`round point/line broadening`) — an isotropic (round) linewidth.
  Leave it at 0 to reproduce older fits; set ~2–5 ppm for the residual width.
- **σ (Cq)** is the Czjzek quadrupolar width (half of dmfit's `sCZ_CQ`) and sets
  the **off-diagonal** quadrupolar tail.

Parameter extraction (δiso, PQ from the two centres of gravity) is a Utilities
tool on the roadmap.
