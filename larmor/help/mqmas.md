# MQMAS (2D) — display, process & fit

> **Multiple-Quantum Magic-Angle Spinning** (Frydman & Harwood 1995) is *the*
> high-resolution experiment for half-integer quadrupolar nuclei (²⁷Al, ¹¹B,
> ²³Na, ¹⁷O…). It refocuses the anisotropic second-order quadrupolar broadening
> into a narrow **isotropic** dimension, resolving crystallographically distinct
> sites that overlap hopelessly in a 1D MAS spectrum. This manual covers opening,
> displaying, phasing, **referencing** (the part everyone gets wrong), and
> **fitting** an MQMAS map in LARMOR.

---

## 1 · The physics in one page

For a half-integer spin $I$, the central transition is unshifted to first order
but carries a **second-order quadrupolar** broadening that MAS cannot remove.
MQMAS correlates a symmetric multiple-quantum coherence (order $p$, e.g. the
triple-quantum $\pm3Q$) evolving in $t_1$ with the single-quantum CT in $t_2$.
Because the second-order broadening of the $p$Q and $1Q$ coherences are
**proportional**, a suitable combination refocuses it — leaving, along F1, only
the **isotropic** shift.

The centre of gravity of each site in the two dimensions is set by the isotropic
chemical shift $\delta_\text{iso}$ and the second-order **quadrupolar isotropic
shift**

$$\delta_\text{QIS} = -\frac{3}{40}\left[\frac{I(I+1)-\frac{3}{4}}{I^2(2I-1)^2}\right]\left(\frac{P_Q}{\nu_0}\right)^2\times 10^6\ \text{(ppm)},\qquad P_Q = C_Q\sqrt{1+\eta_Q^2/3}$$

with Larmor frequency $\nu_0$ (Samoson, Kundla & Lippmaa 1982). The direct
dimension F2 shows the full MAS powder lineshape centred at
$\delta_\text{iso}+\delta_\text{QIS}$; the indirect dimension F1, after the
**shearing** transformation, shows the isotropic position — so a site's F1–F2
offset is a *direct read-out of $P_Q$*, and its position along the isotropic axis
gives $\delta_\text{iso}$.

### The shear, and why F1 needs a convention

Raw MQMAS data lie on a sheared grid: the ridge from an anisotropic site runs at
a spin-dependent slope, not along the diagonal. A **shearing transformation**
(Massiot *et al.* 1996) rotates the isotropic axis to be horizontal/pure. The
scaling of a pure-shift site in the unsheared representation is
$\delta_\text{F1} = c\,\delta_\text{iso}$ with

$$c = -\frac{p-R}{1+R},\qquad R = \left|\frac{C_4(I,\,p/2)}{C_4(I,\,1/2)}\right|,\qquad C_4(I,m) = m\left[18\,I(I+1) - 34m^2 - 5\right].$$

For $I=\frac{5}{2}$ triple-quantum, $R = \frac{19}{12}$ and $c = -\frac{17}{31}$.
**mrsimulator's `ThreeQ_VAS` returns F1 in exactly this sheared convention.**
LARMOR rescales the kernel's F1 axis by $1/c$ at build time so that its internal
convention is the **δ₁-isotropic** one (a pure-CS site sits on the diagonal),
which is what a Bruker/dmfit `2rr` axis already approximates. Only a small
residual **referencing offset** then remains — see §4.

---

## 2 · Open & display

Open the processed `2rr` (or its `pdata` folder). It lands on the **contour map**
in its own workspace. **F2** is the MAS (direct) dimension, **F1** the isotropic
(indirect) one.

- **Contours** — positive / negative / **both** (negatives drawn red); set the
  number of **levels** and the **floor ×σ**.
- **Display mode** — `contour`, `density` (a viridis heat-map), `filled`, or
  `contour + values` (each level labelled with its %); pick a **colormap**
  (viridis / magma / inferno / cividis / gray). Inspired by common
  matplotlib `contourf` figure styles.
- **Projections** — the F2 (top) and F1 (left) skyline/sum projections track the
  map. **Measure** (two markers → Δ in ppm & Hz on each axis) and **Calibrate**
  (click → set ppm) work as on a 1D.
- **Axes** follow the standard NMR/dmfit convention: F2 high-ppm **left**, F1
  high-ppm **top**, so the diagonal runs top-left → bottom-right. **2D ops ▾ ▸
  Flip F2/F1 axis** flips either direction *for display only* (data untouched) to
  match another program or a published figure.
- **2D ops ▾** — Transpose F1↔F2, Reverse F1/F2, extract the **Diagonal → fit**,
  send **projections → fit**, **symmetrize**, save a projection to CSV.

### Overlay a 1D spectrum

**Overlay 1D ▾** puts a 1D on the F2 (direct) or F1 (indirect) projection for
direct comparison: **Current 1D** uses the open spectrum, **From file…** loads
any 1D. Tune **scale** / **fit scale** in the overlay bar; **Clear 1D overlays**
removes them. (For heteronuclear correlations see the **HMQC & correlation**
manual.)

---

## 3 · Phase (if the data isn't already `2rr`-clean)

**Phase 2D** → click 1 or 2 reference peaks → the selected rows (F2) and columns
(F1) appear as full-width traces with the pivot marked → adjust **p0 / p1** with
the ±90/180° steps → **Apply to all rows/cols**. When the hypercomplex quadrants
(`2ri`/`2ir`/`2ii`) are present LARMOR rotates all four exactly, so the
correction is **exact**; with only a real `2rr` it reconstructs the dispersive
part by Hilbert transform (good for re-phasing absorptive data). Full details in
the **2D processing** manual.

**Apply shear** is available for data processed *without* `xfshear` on the
spectrometer. mrsimulator-simulated MQMAS is already in the sheared convention.

---

## 4 · Isotropic-axis (F1) referencing — the one knob to know

Even after the kernel is rescaled to δ₁-isotropic (§1), the experimental F1 axis
is referenced by its own convention (3Q/5Q ratio, transmitter offset, the choice
of δ₁ vs the Amoureux/Massiot scales) and typically sits a few — sometimes tens —
of ppm away from the model. LARMOR therefore carries a **single global F1
reference offset** that slides the whole model along F1 to line up with the data.
It is **independent** of $\delta_\text{iso}$ (which also moves F2, so F2 anchors
it) and it is **auto-fitted** and reported after every 2D fit.

> **If a Czjzek 2D fit is well-shaped but sits off in F1, this is the knob.**
> To pin it (dmfit-style manual referencing): **Decomposition ▸ MQMAS F1
> reference…**; re-run with Cancel to return to auto-fit. The offset is held when
> you pin it and reported in the fit summary.

Extracting δ_iso and $P_Q$ from the two centres of gravity (F1, F2) is a
Utilities tool on the roadmap; it is convention-sensitive and ships validated
against a literature example.

---

## 5 · Fit

1. Pick **Czjzek** (glasses) or **ext. Czjzek / Quad CT** (crystalline) from
   Models, then **click the contour** to place a site (its $\delta_\text{iso}$
   starts at the clicked F1). Placing is sticky — drop as many sites as needed.
2. **Fit.** The MQMAS kernel builds once (cached), then the fitted model overlays
   as per-site contours **in the same colour the site has in the 1D table**, with
   the RMSD and each site's $\delta_\text{iso}$ / $C_Q$ / $\eta$ / σ ± error.
3. Or pull the **isotropic projection** to the workbench and fit it as a 1D, or
   **co-fit** the map with a 1D MAS spectrum of the same sample (see the
   **Multi-dataset & co-fitting** manual).

### Glass broadening (dmfit *CzSimple* parity)

A disordered (glassy) quadrupolar site has **two** independent broadenings, which
LARMOR keeps separate exactly as dmfit's `dCS` and `wid`:

| Parameter | Meaning | Effect on the 2D map |
|---|---|---|
| **σ (Cq)** | Czjzek quadrupolar width; half of dmfit's `sCZ_CQ` | sets the **off-diagonal** quadrupolar tail |
| **dCS** (`shift_fwhm_ppm`) | spread of **isotropic chemical shifts** | elongates the peak **along the diagonal** (δ_iso moves F2 and F1 together) — the hallmark of disorder |
| **line** (`line_fwhm_ppm`) | residual **round** (isotropic) linewidth | a symmetric broadening; leave 0 to reproduce older fits, ~2–5 ppm otherwise |

The Czjzek distribution itself (the $d=5$ EFG statistics behind σ) is documented
in the **Lineshapes** reference.

### Constraining overlapping sites

If a weak site drifts under a strong one (e.g. Al⁽⁵⁾ under Al⁽⁴⁾), constrain it
in the parameter table: cap its δ_iso with a bound (`[.. 40]`) or **pin** it, and
fix σ if the 2D tail is under-determined. Co-fitting with a clean 1D of the same
sample also stabilises the shared shifts.

---

## References

- L. Frydman, J. S. Harwood, "Isotropic spectra of half-integer quadrupolar spins
  from bidimensional MAS NMR", *J. Am. Chem. Soc.* **117**, 5367 (1995).
- A. Medek, J. S. Harwood, L. Frydman, "Multiple-quantum MAS NMR: a new method for
  the study of quadrupolar nuclei in solids", *J. Am. Chem. Soc.* **117**, 12779
  (1995).
- D. Massiot, B. Touzo, D. Trumeau, J. P. Coutures, J. Virlet, P. Florian,
  P. J. Grandinetti, "Two-dimensional magic-angle spinning isotropic
  reconstruction sequences for quadrupolar nuclei", *Solid State Nucl. Magn.
  Reson.* **6**, 73 (1996). *(shearing & F1 referencing)*
- J.-P. Amoureux, C. Fernandez, S. Steuernagel, "Z-filtering in MQMAS NMR",
  *J. Magn. Reson. A* **123**, 116 (1996).
- A. Samoson, E. Kundla, E. Lippmaa, "High resolution MAS-NMR of quadrupolar
  nuclei in powders", *J. Magn. Reson.* **49**, 350 (1982). *(2nd-order shift)*
- Y. Millot, P. P. Man, "Procedures for labeling the high-resolution axis of
  MQMAS NMR spectra", *Solid State Nucl. Magn. Reson.* **21**, 21 (2002).
- S. E. Ashbrook, S. Wimperis, "High-resolution NMR of quadrupolar nuclei in
  solids: the satellite-transition MAS and MQMAS experiments", *Prog. NMR
  Spectrosc.* **45**, 53 (2004). *(review)*
- D. Massiot *et al.*, dmfit, *Magn. Reson. Chem.* **40**, 70 (2002).

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
