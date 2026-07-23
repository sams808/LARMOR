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
- **2D quadrature recombination** — States / TPPI / States-TPPI / Echo–Antiecho /
  QF for the indirect dimension (see **2D processing**).

---

## 3 · Phase

$$S_\text{corr}(\nu) = S(\nu)\,e^{i(\phi_0+\phi_1(\nu-\nu_\text{pivot})/\text{SW})}$$

- **p0 / p1** — zero- and first-order phase; sliders, exact entry, **±90 / 180°**
  quick steps, adjustable step size.
- **Autophase (ACME)** — minimises the entropy of the real spectrum to find p0 and
  p1 automatically (Chen *et al.* 2002); a negative-area criterion is used for
  all-positive powder patterns (QCPMG).

---

## 4 · Baseline & referencing

- **Automatic baseline (arPLS)** — asymmetrically reweighted penalized
  least-squares (Baek *et al.* 2015, $\lambda\approx10^7$): a smooth curve that
  follows the baseline but not the peaks.
- **Manual baseline (PCHIP)** — shape-preserving monotone cubic through
  drag-placed anchors.
- **SR / calibrate** — reference the axis: type an SR (Hz), or click a peak and set
  its ppm. A rigid ppm shift; the raw data is untouched.
- **scale SW / car-ref** — stretch the ppm axis about its centre (correct a
  spectral-width/referencing mismatch between datasets).

---

## 5 · Amplitude, components & algebra

- **scale / offset / normalize** — multiply, add a constant, or normalise
  (max or area) — e.g. to put datasets on a common scale before overlay.
- **magnitude** — $|S| = \sqrt{\text{Re}^2+\text{Im}^2}$ (phase-insensitive
  display).
- **real / imag / conj** — take a single channel or complex-conjugate (spectral
  reversal); inspect the imaginary channel while phasing.
- **extract** — keep a spectral region.
- **combine / align / subtract averages** — algebra between spectra: add, align by
  cross-correlation, or subtract an average reference (ssNake *Subtract
  Averages*). Full background subtraction with least-squares scaling is
  **Process ▸ Subtract a spectrum** (see **Multi-dataset**).
- **pick peaks** — threshold + parabolic interpolation → peak list, and
  *Decomposition ▸ Add a line at every peak*.

---

## 6 · Reproducibility

Every op above is recorded in the recipe with its parameters. *Process ▸
Processing steps* shows the applied sequence; remove any step and LARMOR
re-applies the reduced pipeline from the raw data (steps never silently
compound). A saved `.json` recipe reopens to the identical processed state, and a
dmfit `.fxmla` or CSV export carries the result out.

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
- S.-J. Baek, A. Park, Y.-J. Ahn, J. Choo, *Analyst* **140**, 250 (2015). *(arPLS
  baseline)*
- R. R. Ernst, G. Bodenhausen, A. Wokaun, *Principles of NMR in One and Two
  Dimensions*, Oxford (1987).

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
