# Processing reference — the full operation list

> The complete catalogue of LARMOR's processing operations, the shared source for
> the per-experiment manuals. Each op is a step in the **recipe**: applied live,
> listed in *Process ▸ Processing steps*, removable individually, and replayed
> from the raw data so the whole analysis is reproducible. Notation: $s(t)$ is the
> complex time-domain signal, $S(\nu)$ its spectrum, $T$ the acquisition time.

---

## 1 · Time-domain (before the transform)

### Apodization / window functions

A window $w(t)$ trades resolution for sensitivity (Lindon & Ferrige 1980):

| Op | $w(t)$ | Use |
|---|---|---|
| **EM** (exponential) | $e^{-\pi\,\text{LB}\,t}$ | S/N; matched filter at `LB` = linewidth |
| **GM** (Gaussian) | $\exp[\pi\,\text{LB}\,t-(\pi\,\text{LB}\,t)^2/(2\,\text{GB}\,T)]$ | resolution enhancement |
| **SINE / QSINE** | $\sin/\sin^2\!\left[\pi(1{-}\phi)t/T+\phi\right]$ | sine bells, shift `SSB` = $\phi$ |
| **TRAF** | Traficante window | optimal S/N-vs-resolution (Traficante 1987) |
| **Shifted Gaussian** | Gaussian centred on the echo top | whole-echo data |

### Other time-domain ops

- **TDeff** — truncate the FID to a chosen number of points (drop a noisy tail).
- **shift_fid** — circular/linear shift; **swap echo** — move a whole echo's top
  to $t=0$ so its transform is pure absorption.
- **Linear prediction (LP)** — extend a truncated FID or repair corrupted first
  points by an autoregressive model fitted to the signal (Barkhuijsen *et al.*
  1985). Especially useful for the short, truncated indirect dimension of a 2D.
  *(An LPSVD upgrade is on the roadmap.)*
- **Hilbert** — reconstruct the imaginary channel from a real-only signal
  (Kramers–Kronig), so a real spectrum can be re-phased.

---

## 2 · The transform

- **ZF (zero-fill)** — pad with zeros before the FFT; one zero-fill is
  information-preserving and always worth doing (Bartholdi & Ernst 1973).
- **FCOR** — halve the first FID point ($s(0)\to s(0)/2$) to suppress the DC
  baseline offset it would otherwise create.
- **FT / IFT** — forward and inverse Fourier transform (`ift` round-trips `ft` to
  ~$10^{-9}$), letting you return to the FID, re-apodize, and transform again.
  On a processed spectrum (TopSpin 1r, CSV) `ift` derives the spectral width from
  the ppm axis (`SW = Δδ · n · SFO1`, so the Larmor frequency must be set) and
  parks the axis so the next `ft` restores it exactly; zero-filling in between
  keeps the zero-frequency bin, so peaks stay at their ppm. A uniform ppm grid is
  assumed. *Process ▸ FID ⇄ spectrum* (**Ctrl+T**) shows the FID the transform
  sees and re-applies the window live; for a spectrum that did not come from a
  raw fid the panel's **re-apodize this spectrum** checkbox builds
  `hilbert → ift → window → ft` (Hilbert first is mandatory: the inverse
  transform of a real-only spectrum is two-sided, and a one-sided window on it
  loses half the signal).
- **2D quadrature recombination** — States / TPPI / States-TPPI / Echo–Antiecho /
  QF for the indirect dimension (see **2D processing**).

---

## 3 · Phase

$$S_\text{corr}(\nu) = S(\nu)\,e^{i(\phi_0+\phi_1(\nu-\nu_\text{pivot})/\text{SW})}$$

- **p0 / p1** — zero- and first-order phase; sliders, exact entry, **±90 / 180°**
  quick steps.
- **Pivot** ($\nu_\text{pivot}$) — a draggable line (default: the tallest peak,
  shown while the Processing panel is open) about which p1 rotates, TopSpin-style,
  so the peak under it stays in phase.
- **Drag to phase** — *Drag to phase* in the panel or *Process ▸ Phase ▸ Drag to phase*
  (**Ctrl+P**): a horizontal drag on the spectrum changes p0 (0.25° per pixel), a
  vertical drag changes p1 (1° per pixel, up = positive) about the pivot; each
  drag locks to the direction it starts in; **Shift** = fine (×0.1), **Ctrl+drag**
  pans, **Esc** stops. One Undo reverts a whole drag. In pdata mode *Hilbert
  first* is switched on automatically, since the displayed spectrum carries no
  imaginary part. While the FID is displayed (Ctrl+T) the gesture is suspended
  and a drag pans.
- Moving the pivot while phasing re-expresses the zero order,
  $\phi_0' = \phi_0 + \phi_1\,\Delta\nu_\text{pivot}/\text{SW}$, so the spectrum
  on screen does not change; the recipe records the new pivot.
- **Autophase** — finds p0 and p1 automatically: a fine p0 sweep maximising the
  positive real signal with a penalty on negative excursions, then a
  Nelder–Mead refinement of (p0, p1) — the criterion that holds up on wide
  all-positive powder patterns (QCPMG included), where the entropy criterion
  (ACME, Chen *et al.* 2002; available as the recorded op's `method: acme`)
  can prefer a dispersive-looking answer. The angles it finds are **written
  into the Phase controls** as an ordinary p0 / p1 step (like TopSpin's `apk`
  filling PHC0 / PHC1), so a nudge of the sliders continues from the
  autophased state. Phasing needs the imaginary channel: TopSpin's `1i` is
  read next to the `1r`; a source without one (CSV, dmfit, magnitude data)
  gets it reconstructed by a Hilbert transform first — *Signal ▸ Hilbert
  first* is ticked for you.

---

## 4 · Baseline & referencing

- **Automatic baseline (polynomial)** — *Process ▸ Baseline ▸ Polynomial*: a
  polynomial (order 3 from the menu; the panel sets the order) fitted to every
  point, the points sticking up above it clipped, refitted until stable — the
  standard automatic baseline of most NMR software (`baseline`). For a smooth
  baseline that is not a polynomial, or one under broad peaks, use the
  **pybaselines** methods below (arPLS and its relatives).
- **Manual baseline (PCHIP, dmfit-style anchors)** — click **Pick anchors**, drop
  as many points as needed (drag to shape); a shape-preserving monotone cubic
  previews live and is **subtracted automatically when you exit anchor mode**.
- **Iterative baseline (dead-time; Yon et al. 2020)** — *Process ▸ Baseline ▸
  Iterative*. For the **rolling baseline caused by receiver dead time** in
  pulse-acquire (single-pulse) MAS, where a polynomial or arPLS fails. It builds
  the baseline **iteratively**: at each pass a **histogram filter** picks the
  noise/baseline points automatically (the noise band around the histogram mode —
  positive peaks *and* negative spikes excluded), a smoothing spline is fit to
  them, and the result is subtracted; repeat until the correction falls below the
  noise. Optionally the estimate is restricted to **broad** components by keeping
  only the first time-domain points (the **dead-time window**, ≈ 2·DE/DW) — this
  is what separates a genuine dead-time roll from real (narrow) peaks. It opens an
  **interactive dialog** with a live preview (spectrum + estimated baseline, and
  the corrected result below) and three controls — dead-time points (≈ 2·DE/DW; 0 =
  off), spline **smoothness**, and the histogram **threshold** — so you can tune and
  see the effect before **Apply**. It leaves an already-flat spectrum essentially
  untouched, so it is safe to try.
- **SR / calibrate** — reference the axis: type an SR (Hz), or click a peak and set
  its ppm. A rigid ppm shift; the raw data is untouched.
- **scale SW / car-ref** — stretch the ppm axis about its centre (correct a
  spectral-width/referencing mismatch between datasets).

### pybaselines — *Process ▸ Baseline ▸ pybaselines*

Thirteen baseline algorithms of the **pybaselines** library (Erb 2022), each
recorded as one `pybaseline` step (`method` plus its parameters) so a saved fit
replays it exactly. The submenu opens the dialog on one method (*arPLS…*,
*asLS…*, *airPLS…*, *SNIP…*, *ModPoly…*, *Rolling ball…*, *Morphological…*);
*All methods…* lists every one, grouped by family. The dialog shows the
spectrum with the estimated baseline (top) and the corrected result (bottom)
and updates them as you move a control; the step is appended to the recorded
processing chain, so an earlier phase or window is kept. A spectrum longer
than 8192 points is previewed on every n-th point with λ and the windows
rescaled accordingly (the status line says so) and **Apply** runs on all
points. The baseline is subtracted from the real part; the imaginary part is
kept, so phasing afterwards still works.

**Which method.**

| Family | Methods | Suits | What to set |
|---|---|---|---|
| Whittaker smoothers | arPLS, asLS, airPLS, IarPLS, drPLS, asPLS | a smooth, slowly curving baseline under peaks of any width — the usual case in solid-state NMR (probe background, a residual dead-time roll, broad humps) | λ, the smoothness |
| Iterative polynomials | ModPoly, IModPoly, penalized poly | a gently curved baseline that a low-order polynomial describes; IModPoly for a noisy spectrum | the polynomial order (2–4) |
| Peak clipping and morphology | SNIP, rolling ball, Mor, MorMol | narrow lines on a wide rolling baseline, where a smoother would start to follow the peaks | the half-window in points, wider than the widest peak's half-width |

Start with **arPLS**: its weights come from the noise, so λ is the only choice.
asLS adds a fixed asymmetry *p* (0.001–0.05); airPLS discounts the peak regions
harder at each pass (dense, overlapping peaks); IarPLS is gentler on weak,
broad features; drPLS follows a strong curvature; asPLS lets λ vary along the
axis (the slowest — prefer it on spectra up to ~16k points). SNIP and the
morphological methods need the half-window above the widest peak's half-width
or the baseline climbs into the peak; MorMol removes the flat steps of a plain
opening. A baseline that dips under a broad, weak feature is the sign that the
smoothing is too weak (λ too low, or for asLS *p* too high, or a window too
narrow).

**λ rule of thumb.** λ is written in points, so it depends on how densely the
spectrum is sampled: the same curve sampled *f* times more finely has every
second difference *f*² smaller and costs *f*⁴ more penalty. LARMOR's default
is 10⁵ (arPLS) at 2048 points and scales as (*n* / 2048)⁴, rounded to a decade
— 10⁶ at 4096, 10⁹ at 16k, 10¹¹ at 64k points (`pybaseline.suggest_lam`). From
there move by decades: one down when the baseline misses a real curvature, one
up when it starts to follow the peaks. The window defaults are 1/32 of the
spectrum (`pybaseline.suggest_half_window`).

**Citing.** D. Erb, "pybaselines: A Python library of algorithms for the
baseline correction of experimental data", Zenodo, doi:10.5281/zenodo.5608581
— together with the method's own paper (arPLS: Baek *et al.* 2015; asLS: Eilers
& Boelens 2005; airPLS: Zhang *et al.* 2010; ModPoly: Lieber & Mahadevan-Jansen
2003; IModPoly: Zhao *et al.* 2007; penalized poly: Mazet *et al.* 2005; SNIP:
Ryan *et al.* 1988; rolling ball: Kneen & Annegarn 1996; MorMol: Koch *et al.*
2017; the IarPLS, drPLS and asPLS papers are listed in the library's
documentation). Every algorithm's parameters and reference:
https://pybaselines.readthedocs.io.

---

## 5 · Amplitude, components & algebra

- **scale / offset / normalize** — multiply, add a constant, or normalise
  (max or area) — e.g. to put datasets on a common scale before overlay.
- **magnitude** — $|S| = \sqrt{\text{Re}^2+\text{Im}^2}$ (phase-insensitive; a
  destructive step — the **|S|** display channel shows the same without altering
  the data).
- **real / imag / conj** — take a single channel or complex-conjugate (spectral
  reversal). These are destructive pipeline steps; to merely *look at* the
  imaginary channel while phasing use the non-destructive **Display channel**
  (Real / Imaginary / |S|, **Ctrl+I**, the radios in the Processing panel), which
  changes nothing in the pipeline.
- **extract** — keep a spectral region.
- **combine / align / subtract averages** — algebra between spectra: add, align by
  cross-correlation, or subtract an average reference (ssNake *Subtract
  Averages*). Full background subtraction with least-squares scaling is
  **Process ▸ Region / algebra ▸ Subtract a spectrum** (see **Multi-dataset**).
- `wurst_correct` — divide out a WURST/chirp sweep's excitation profile
  (centre, sweep width, order N, floor); **Process ▸ Region / algebra ▸ WURST excitation
  profile…**. The amplitude half of swept-pulse physics — autophase's p2
  handles the phase half.
- Patterns wider than one acquisition window: **Process ▸ Region / algebra ▸ Stitch
  frequency-stepped (VOCS) spectra…** combines sub-spectra acquired at
  stepped transmitter offsets (skyline or coverage-weighted average, with a
  per-edge trim for the excitation roll-off) into one spectrum, and records
  the sources in the recipe's provenance.
- **pick peaks** — threshold + parabolic interpolation → peak list, and
  *Edit ▸ Add a line at every peak*.

---

## 6 · Reproducibility

Every op above is recorded in the recipe with its parameters. *Process ▸
Processing steps* shows the applied sequence; remove any step and LARMOR
re-applies the reduced pipeline from the raw data (steps never silently
compound). A saved `.json` recipe reopens to the identical processed state, and a
dmfit `.fxmla` or CSV export carries the result out.

Where a chain replays from is decided by its first domain-restricted step
(`processing.chain_start_domain`): a chain beginning with a time-domain op (`em`,
`zf`, `ft`, …) needs the raw fid; a chain beginning with a frequency-domain op and
then `ift` (a re-apodized 1r or CSV: `hilbert, ift, em, ft`) replays from the
processed data. *File ▸ Open FID* records the chain behind the spectrum it hands
over (window, zero-fill, `ft`, then the typed phase or `autophase`), and a raw
`fid` opened directly records its preview chain. The FID / spectrum display and
the display channel are view state, not recipe content. A batch reprocessed
from its fids (*Multi-dataset*, §6) records the common chain — `tdeff` (TopSpin's
TDeff halved: it counts real points, the fid is complex), the window, `zf` to
SI, `ft` with the spectrum's own referencing offset, then the spectrum's own
`phase` (p0 = −PHC0, p1 = PHC1, pivot at the high-frequency end) or
`autophase` — on every recipe with `processing_from_raw` and the EXPNO as
source, replayed from the instrument fid on reopen.

A saved recipe also records where its numbers came from: `source_sha256`, the
SHA-256 of the exact data file the fit saw (`pdata/<procno>/1r`, or the `fid`
behind a *File ▸ Open FID* fit); `acquisition`, the flat block read from acqus,
`pdata/<procno>/procs`, the title, `uxnmr.info` and the booking sidecar (with
the procno the fit was made on, so a pdata/2 fit reopens pdata/2); and
`software`, the software that produced the fit (LARMOR version and git commit
when run from a checkout, mrsimulator, lmfit, numpy, scipy, nmrglue, Python,
the fitting time). Reopening a recipe re-hashes the data file and re-reads the
SR, and reports in the status bar — never a dialog — when the data file
differs from the fitted one (*source data changed*: TopSpin rewrote the 1r),
when the axis was re-referenced since the fit (*SR a → b Hz*), and when a
recorded read-out change (e.g. the 0.13.0 Czjzek read-outs, for recipes with a
Czjzek-family site) lies between the fitting version and the current one. A
recipe without a hash (written before this record existed) stays silent.

---

## References

- J. C. Lindon, A. G. Ferrige, *Prog. NMR Spectrosc.* **14**, 27 (1980).
  *(digitisation & apodization)*
- D. D. Traficante, G. A. Nemeth, *J. Magn. Reson.* **71**, 237 (1987). *(TRAF)*
- E. Bartholdi, R. R. Ernst, *J. Magn. Reson.* **11**, 9 (1973). *(zero-filling)*
- H. Barkhuijsen, R. de Beer, W. M. M. J. Bovée, D. van Ormondt, "Retrieval of
  frequencies, amplitudes, damping factors and phases from time-domain signals"
  (LPSVD), *J. Magn. Reson.* **61**, 465 (1985).
- L. Chen, Z. Weng, L. Goh, M. Garland, *J. Magn. Reson.* **158**, 164 (2002).
  *(ACME autophase)*
- M. Yon, F. Fayon, D. Massiot, V. Sarou-Kanian, "Iterative baseline correction
  algorithm for dead time truncated one-dimensional solid-state MAS NMR spectra",
  *Solid State Nucl. Magn. Reson.* **110**, 101699 (2020),
  doi:10.1016/j.ssnmr.2020.101699; github.com/maximeYon/Baseline_Corrector.
  *(iterative dead-time baseline)*
- S.-J. Baek, A. Park, Y.-J. Ahn, J. Choo, *Analyst* **140**, 250 (2015). *(arPLS
  baseline)*
- D. Erb, "pybaselines: A Python library of algorithms for the baseline
  correction of experimental data", Zenodo, doi:10.5281/zenodo.5608581;
  https://pybaselines.readthedocs.io. *(pybaselines)*
- P. H. C. Eilers, H. F. M. Boelens, "Baseline correction with asymmetric least
  squares smoothing", Leiden University Medical Centre report (2005). *(asLS)*
- Z.-M. Zhang, S. Chen, Y.-Z. Liang, *Analyst* **135**, 1138 (2010). *(airPLS)*
- C. A. Lieber, A. Mahadevan-Jansen, *Appl. Spectrosc.* **57**, 1363 (2003).
  *(ModPoly)*
- J. Zhao, H. Lui, D. I. McLean, H. Zeng, *Appl. Spectrosc.* **61**, 1225 (2007).
  *(IModPoly)*
- V. Mazet, C. Carteret, D. Brie, J. Idier, B. Humbert, *Chemom. Intell. Lab.
  Syst.* **76**, 121 (2005). *(penalized polynomial)*
- C. G. Ryan, E. Clayton, W. L. Griffin, S. H. Sie, D. R. Cousens, *Nucl.
  Instrum. Methods B* **34**, 396 (1988). *(SNIP)*
- M. Kneen, H. Annegarn, *Nucl. Instrum. Methods B* **109–110**, 209 (1996).
  *(rolling ball)*
- M. Koch, C. Suhr, B. Roth, M. Meinhardt-Wollweber, *J. Raman Spectrosc.* **48**,
  336 (2017). *(MorMol)*
- R. R. Ernst, G. Bodenhausen, A. Wokaun, *Principles of NMR in One and Two
  Dimensions*, Oxford (1987).

*LARMOR — Sam Soudani, NOME Laboratory (McCloy group), Washington State University.*
