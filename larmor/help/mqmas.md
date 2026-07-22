# MQMAS (2D) — display, process, fit

## Open
Open the processed `2rr` (or its pdata folder). It lands on the **contour map**
in its own workspace. F2 is the MAS (direct) dimension, F1 the isotropic
(indirect) one.

## Display
- **contours**: positive / negative / both; **levels** and **floor ×σ**.
- The **projections** (top = F2, left = F1) track the map. **Measure** and
  **Calibrate** work on the 2D as on a 1D.
- **2D ops ▾**: Transpose F1↔F2, Reverse F1/F2, Diagonal → fit, projections →
  fit.

## Phase (if needed)
**Phase 2D** → click 1 or 2 reference peaks → the selected rows/columns show as
full-width traces with the pivot marked → p0/p1 + ±90/180° → Apply. With the
hypercomplex quadrants present the correction is exact.

## Shear / reference
**Apply shear** for data processed without xfshear. (mrsimulator MQMAS is already
sheared.) Use the F1 axis referencing conventions for your spin/method.

## Fit
1. Pick **Czjzek** (glasses) or **ext. Czjzek / Quad CT** (crystalline) from
   Models, then **click the contour** to place a site (its δiso starts at the
   clicked F1).
2. **Fit** — the MQMAS kernel builds once, then the fitted contours overlay in
   orange with the RMSD and per-site δiso/Cq/η/σ.
3. Or pull the **isotropic projection** to the workbench and fit it as a 1D, or
   **co-fit** it with a 1D MAS spectrum (Decomposition ▸ Co-fit datasets).

Parameter extraction (δiso, PQ from the two centres of gravity) is a Utilities
tool on the roadmap.
