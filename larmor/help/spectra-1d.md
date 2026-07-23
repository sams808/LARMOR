# 1D spectra — processing & fitting

> Everything you do to a one-dimensional spectrum: **open** it, **process** it
> from the FID if need be, and **fit** it to a physical model with uncertainties.
> Shared processing steps (apodization, phasing, baseline, referencing) are
> documented here in full and repeated in the other data-type manuals so you
> never have to hunt across documents.

---

## At a glance

```
FID ──▶ apodize ──▶ zero-fill ──▶ Fourier transform ──▶ phase ──▶ baseline ──▶ reference
                                                                                   │
                                                                                   ▼
                              add lineshapes ──▶ fit ──▶ quantify (% ± err) ──▶ export
```

Every step is stored in the **recipe** and replayed from the raw data, so a fit
is reproducible and diffable.

---

## 1 · Open

**File ▸ Open**, the **Explorer**, or drag-and-drop. LARMOR reads a Bruker `1r`,
a raw `fid`, a dmfit `.fxmla`, a LARMOR `.json` recipe, or a two-column
`.csv`/`.txt`. A raw `fid` opens with the processing preview so you can apodize
and phase **before** the transform. (See **Getting started** for the full import
matrix.)

---

## 2 · Process (the live right-hand panel)

The panel is **live**: every control re-applies the whole pipeline on a short
debounce, rebuilding from an unprocessed baseline each time so steps never
compound. Below, $t$ is acquisition time, $T$ the acquisition length, and
$\nu$ frequency.

### Apodization (window functions — WDW/LB/GB/SSB)

Multiplying the FID by a window $w(t)$ trades resolution against sensitivity
(Lindon & Ferrige 1980). LARMOR offers the TopSpin/ssNake set:

- **Exponential (EM)** — $w(t) = e^{-\pi\,\text{LB}\,t}$. A line-broadening of
  `LB` Hz; the classic S/N-boosting **matched filter** when `LB` equals the
  natural linewidth.
- **Gaussian (GM)** — $w(t) = \exp\!\big[\pi\,\text{LB}\,t - (\pi\,\text{LB}\,t)^2/(2\cdot\text{GB}\,T)\big]$,
  resolution enhancement (negative `LB`, `GB` sets the Gaussian maximum).
- **Sine / Sine² (SINE/QSINE)** — $w(t) = \sin\!\big(\pi(1-\phi)\,t/T + \phi\big)$
  with a shift $\phi$ (`SSB`); $\phi = \pi/2$ is a pure cosine.
- **TRAF** — the Traficante window (Traficante & Nemeth 1987), an optimal
  S/N-vs-resolution compromise.

### Zero-filling & the transform

**ZF** pads the FID with zeros before the FFT: one zero-fill (doubling the
points) is *information-preserving* — it interpolates the real part from the
imaginary and is always worth doing (Bartholdi & Ernst 1973). LARMOR then applies
**FCOR** (first-point scaling: the $t=0$ point is halved to remove a DC
baseline offset) and the FFT. **TDeff** truncates the FID to a chosen number of
points first (drop a noisy tail).

### Phasing (p0 / p1)

A complex spectrum is rotated to pure absorption by
$S_\text{corr}(\nu) = S(\nu)\,e^{i(\phi_0 + \phi_1\,\nu/\text{SW})}$: a
frequency-independent **p0** and a linear **p1** (pivoted). The panel gives
sliders, an exact type-in, and **−90 / +90 / 180°** quick steps (two +90 equals
one 180). **Autophase** minimises the spectral entropy of the real part (the
**ACME** criterion of Chen *et al.* 2002), which finds p0 **and** p1 robustly
even on crowded spectra.

### Baseline

- **Automatic** — an asymmetrically reweighted penalized least-squares baseline
  (**arPLS**, Baek *et al.* 2015): it iteratively fits a smooth curve
  ($\lambda \approx 10^7$) that follows the baseline but not the peaks. Robust
  for rolling baselines under broad lines.
- **Manual (PCHIP)** — drop anchor points and LARMOR interpolates a
  shape-preserving monotone cubic through them; drag anchors live.

### Referencing (SR / Calibrate)

Type a spectral-reference **SR** (Hz), or **Process ▸ Calibrate**, click a peak
(it snaps to the local maximum), and set its known ppm — LARMOR reports the
resulting SR. Double-click the **experiment strip** to edit nucleus / field /
νrot / SR, or copy the SR from another spectrum. All referencing is a rigid ppm
shift of the axis; the raw data is untouched.

> **Processing history.** *Process ▸ Processing steps* lists every applied op;
> remove any one and LARMOR re-applies the reduced pipeline. The full op list is
> in the **Processing reference** manual.

---

## 3 · Fit

1. **Pick a model** from the **Models** menu (or toolbar) and **click the
   spectrum** to drop lines. The mode is *sticky* — place as many as you like;
   click the model again or press **Esc** to stop.
2. **Adjust.** Drag the **paddles** on the plot (square = position + amplitude,
   round side-handles = width) or edit the **Fit-Parameters** spreadsheet at the
   bottom. Each site is a lettered row (`A`, `B`, …); a cell accepts
   - a **value** (in ppm, or `300Hz` / `1.5kHz` for a width),
   - a **bound** `[0..100]` (bounded cells get a teal border),
   - a **link** to another line by its letter — `A`, `A+20`, `A+20kHz`, `0.5B` —
     with error propagation, or
   - a **pin** ☑ to fix the parameter.
3. **Fit** (F5). It minimises $\chi^2=\sum_i w_i\,[y_i - f(x_i)]^2$ by
   Levenberg–Marquardt (via **lmfit**), and every fitted value comes back with a
   **standard error** from the covariance matrix. Read **RMSD** and **S/N** next
   to the buttons.
4. **Report** (F6) gives the **quantification**: each site's integrated area as a
   population **% ± error**. Use the area-normalised `gl_norm` model (or a Czjzek
   site) so an amplitude *is* a population.

**Zones** restrict the fit to chosen spectral regions (union of intervals) — fit
only where the model is valid and let peaks outside float frozen. **Auto Fit**
does a multi-start search to escape local minima. **Errors Analysis** profiles
$\chi^2$ around a chosen parameter to show its true (possibly asymmetric)
confidence interval, beyond the linear covariance estimate.

### Which lineshape?

| Model | Use it for |
|---|---|
| **Gauss/Lorentz** (pseudo-Voigt) | symmetric lines, quick fits |
| **Gauss/Lorentz (area)** | as above, amplitude = integral (quantification) |
| **Voigt (true)** | separable Gaussian ⊗ Lorentzian broadening |
| **J-multiplet** | scalar J splitting to *n* equivalent spins |
| **Czjzek / ext. Czjzek** | amorphous quadrupolar sites (glasses) |
| **Quad CT / 1st / +CSA** | crystalline quadrupolar sites |
| **CSA powder** | spin-½ shielding anisotropy (+ sidebands) |
| **Spectrum (background)** | fit a measured impurity/phase's amplitude & shift |

The full physics, equations and literature for each are in the **Lineshapes —
models & physics** reference (**? ▸ Lineshapes**).

---

## 4 · Measure & export

- **Tools ▸ Integrals & measurements** — drag regions → integral, %, centre of
  mass, FWHM; **Copy CSV**.
- **File ▸ Copy plot** (with all fitted lines) / **Save plot image** (PNG/SVG)
  for slides and papers.
- **File ▸ Save fit as** — `txt` / `csv` / `json` / **dmfit `.fxmla`**; writes the
  data, model, residual and every component (the dmfit export round-trips the
  Czjzek σ ↔ `sCZ_CQ = 2σ` relation).
- **File ▸ Save spectrum as** — a reopenable CSV with a metadata header.

---

## References

- J. C. Lindon, A. G. Ferrige, "Digitisation and data processing in Fourier
  transform NMR", *Prog. NMR Spectrosc.* **14**, 27 (1980). *(apodization)*
- D. D. Traficante, G. A. Nemeth, "The TRAF window function", *J. Magn. Reson.*
  **71**, 237 (1987).
- E. Bartholdi, R. R. Ernst, "Fourier spectroscopy and the causality principle",
  *J. Magn. Reson.* **11**, 9 (1973). *(zero-filling)*
- L. Chen, Z. Weng, L. Goh, M. Garland, "An efficient algorithm for automatic
  phase correction of NMR spectra based on entropy minimization" (ACME),
  *J. Magn. Reson.* **158**, 164 (2002).
- S.-J. Baek, A. Park, Y.-J. Ahn, J. Choo, "Baseline correction using
  asymmetrically reweighted penalized least squares smoothing" (arPLS),
  *Analyst* **140**, 250 (2015).
- M. Newville *et al.*, **lmfit**: non-linear least-squares minimization for
  Python, doi:10.5281/zenodo.11813 (2014). *(the optimiser + uncertainties)*
- R. R. Ernst, G. Bodenhausen, A. Wokaun, *Principles of NMR in One and Two
  Dimensions*, Oxford (1987). *(general reference)*

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
