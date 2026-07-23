# 2D processing — transform, phase, contour, measure

> The processing common to **every** 2D dataset — MQMAS, HMQC, DQ-SQ, HETCOR — is
> here: the second Fourier transform and its **quadrature** modes, **hypercomplex
> phasing** done exactly, the contour display, and measuring/calibrating on the
> map. Experiment-specific analysis lives in the **MQMAS** and **HMQC &
> correlation** manuals.

---

## 1 · From `ser` to a spectrum

A 2D is a stack of FIDs acquired at incremented $t_1$. Making it a spectrum means
two transforms — one in $t_2$ (direct), one in $t_1$ (indirect) — plus a rule for
recovering **sign** in the indirect dimension, since a single amplitude-modulated
$t_1$ series cannot tell $+\nu_1$ from $-\nu_1$. LARMOR reads the acquisition mode
and recombines accordingly:

| $t_1$ quadrature mode | Recombination |
|---|---|
| **States** (Haberkorn–Ruben) | cosine + i·sine → phase-sensitive, pure absorption |
| **TPPI** | time-proportional phase incrementation → real FT, folded reference |
| **States-TPPI** | States with alternating sign to move axial peaks to the edge |
| **Echo–Antiecho** | gradient/phase-encoded P- and N-type coherence combined |
| **QF** | magnitude (no sign) — a last resort |

Opening a raw `ser` gives a 2D preview so you can process it; the reader loads the
**hypercomplex quadrants** (`2rr`, `2ri`, `2ir`, `2ii`) when they exist, which is
what makes exact phasing possible (§3). Apodization, zero-fill and FT are applied
per dimension (see the **Processing reference**); the indirect dimension is
usually short and benefits most from a window and zero-fill.

---

## 2 · Contour display

- **contours**: positive / negative / **both** (negatives red); set the number of
  **levels** and the **floor ×σ** above the noise.
- **display mode**: `contour`, `density` (viridis heat-map), `filled`, or
  `contour + values` (levels labelled with %); pick a **colormap**.
- **projections**: F2 (top) and F1 (left) skyline/sum projections track the map;
  pull either to the workbench to treat it as a 1D.
- **cursor readout**: live F2 / F1 / z under the mouse.
- performance: the grid is decimated for drawing, so even a 128 × 2048 map pans
  and zooms smoothly.

---

## 3 · Phasing — exact, hypercomplex

**Phase 2D** → click **1 or 2** reference peaks → the rows through them (F2) and
columns (F1) are shown **full-width**, stacked, with a red pivot line → adjust
**p0 / p1** with the sliders and **−90 / +90 / 180°** buttons → **Apply to all
rows/cols**; **Re-pick** or **Reset** to start over.

$$S_\text{corr}(\nu) = S(\nu)\,e^{i(\phi_0+\phi_1(\nu-\nu_\text{pivot})/\text{SW})}$$

is applied independently in each dimension. **With the four quadrants present the
correction is exact** — LARMOR rotates `2rr`/`2ri`/`2ir`/`2ii` together, so a
mis-set 50° recovers to ~10⁻⁹. When only a real `2rr` is available it
reconstructs the dispersive part by a **Hilbert transform** (Kramers–Kronig),
which is exact for *re*-phasing already-absorptive data and good enough
interactively.

---

## 4 · Measure & calibrate

- **Measure** — drop two markers → Δ in **ppm and Hz** on each axis, shown in the
  cursor label.
- **Calibrate** — click a peak → enter its known F2 and F1 ppm → the axes shift
  rigidly to reference the map. As always, the raw data is not modified.

---

## 5 · 2D operations (`2D ops ▾`)

| Op | Does |
|---|---|
| **Transpose F1↔F2** | swap the two axes |
| **Reverse F1 / F2** | flip an axis (spectral reversal) |
| **Flip F2 / F1 axis** | change display direction only (data untouched) |
| **Symmetrize** | average the map with its transpose (square grids) |
| **Diagonal → fit** | extract the diagonal trace to a workbench |
| **Projection → fit / CSV** | send or save an F1/F2 projection |
| **Apply shear** | shear an unsheared MQMAS dataset (see the MQMAS manual) |

---

## References

- D. J. States, R. A. Haberkorn, D. J. Ruben, "A two-dimensional nuclear Overhauser
  experiment with pure absorption phase in four quadrants", *J. Magn. Reson.*
  **48**, 286 (1982). *(States quadrature)*
- D. Marion, K. Wüthrich, "Application of phase-sensitive 2D NMR by time-
  proportional phase incrementation (TPPI)", *Biochem. Biophys. Res. Commun.*
  **113**, 967 (1983).
- A. G. Palmer III, J. Cavanagh, P. E. Wright, M. Rance, "Sensitivity improvement
  in proton-detected two-dimensional heteronuclear correlation NMR" (echo–
  antiecho), *J. Magn. Reson.* **93**, 151 (1991).
- R. R. Ernst, G. Bodenhausen, A. Wokaun, *Principles of NMR in One and Two
  Dimensions*, Oxford (1987). *(2D transforms, quadrature, phasing)*
- J. Keeler, *Understanding NMR Spectroscopy*, 2nd ed., Wiley (2010). *(a very
  readable account of 2D quadrature detection)*

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
