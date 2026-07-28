# QCPMG — wide-line echo trains

> **QCPMG** (Quadrupolar Carr–Purcell–Meiboom–Gill) records a *train of echoes*
> after a single excitation. It is the method of choice for **broad quadrupolar
> lineshapes** (³⁵Cl, ⁸¹Br, low-field ²⁷Al…) where the signal dephases far faster
> than it relaxes — capturing many echoes in one shot multiplies the
> signal-to-noise. There are two ways to turn the train into a spectrum, and
> **which one you fit matters.**

---

## At a glance

| Representation | How | S/N | Fit it? |
|---|---|---|---|
| **Sum echo** (absorption) | split the train, coadd the echoes, process one echo | high | ✅ **yes** — a continuous powder lineshape |
| **Spikelets** | Fourier-transform the whole train | spectacular | ❌ no — a comb with no lineshape *between* spikes |

The spikelet **maxima** trace the powder pattern, but a smooth model cannot fit
the comb — fit the **sum echo** and use spikelets only to inspect the manifold.

---

## 1 · The two spectra

**Sum echo (absorption) — fit this.** Split the train into individual echoes, add
them together, and process the resulting single echo. You get a normal,
continuous powder lineshape at high S/N — send it to the workbench and model it
with `ext Czjzek`, `Quad CT`, etc.

**Spikelets.** Fourier-transform the whole train without splitting. You get a
manifold of sharp spikelets spaced by $1/\tau_\text{echo}$, whose maxima trace the
powder pattern. All the signal is concentrated in the spikes, so S/N is
spectacular — but **the area between spikelets carries no lineshape
information**, and fitting a smooth model onto the comb is meaningless.

---

## 2 · The controls, step by step (mirrors the ssNake tutorial)

**period (points per echo)** — the number of data points between two consecutive
echo tops. Auto-detected from the autocorrelation of the magnitude FID; it equals
$\text{SW}/(\text{spikelet spacing})$ (e.g. $200\,\text{kHz}/2000\,\text{Hz} =
100$ points). It is the split length; if the echoes don't align, nudge it.

**echo top** — the index of the echo maximum *within one period*. Used for
**whole-echo processing**: the summed echo is circularly shifted so its top sits
at $t=0$, which makes the transform come out in **pure absorption** (ssNake's
*swap echo* step). Auto-detected from a clean mid-train echo.

**T₂ (echo-top decay → the evolution time)** — the intensity at the echo top,
echo by echo, decays as

$$I_k = I_0\,e^{-t/T_2},\qquad t = k\,\tau_\text{echo},\qquad \tau_\text{echo} = \frac{\text{period}}{\text{SW}}.$$

LARMOR fits this and reports **T₂ in ms** — your transverse *evolution time*.
Beside it is the **matched apodization**

$$\text{LB} = \frac{1}{\pi\,T_2}\ \text{(Hz)},$$

the Lorentzian that weights the echoes by how much signal they still carry.

**T₂ weighting** — when ticked, echo $k$ is scaled by
$e^{-k\,\tau_\text{echo}/T_2}$ before summing: the **matched filter**. Echoes that
have mostly decayed (noise-dominated) count less, giving the best S/N — exactly
ssNake's "apply a Lorentzian $\text{LB}=1/(\pi T_2)$ along the echo dimension".

**GB (Hz)** — Gaussian broadening of the summed echo before the transform;
smooths the lineshape.

**p0 / p1 and step** — zero- and first-order phase. The **step** sets how much the
arrows / wheel move the phase per click — lower it for fine control. **Autophase**
optimises p0 **and** p1 by minimising the negative area of the real spectrum (a
robust criterion for an all-positive powder pattern). Fine-tune p1 by hand if the
wings aren't flat.

---

## 3 · Recommended workflow

1. **Open the raw fid** of the QCPMG experiment (the echo train).
2. Check **period** and **echo top** on the echo-train plot (red lines mark the
   period). The auto values are usually right.
3. Read **T₂** off the decay plot; optionally tick **T₂ weighting**.
4. Choose **sum echo (absorption)** mode.
5. **Autophase**, then nudge **p1** (small step) until the lineshape is clean
   absorption.
6. Adjust **GB** to taste.
7. **Send to fit →**, add a quadrupolar model (`ext Czjzek`, `Quad CT`, …) and fit
   as usual (see the **Lineshapes** reference).

Use **spikelets** mode only to inspect the manifold or confirm the spikelet
spacing; do not fit it directly.

---

## 4 · Infinite-field δiso from two (or more) fields

For a half-integer quadrupolar nucleus the **central-transition centre of
gravity** carries a second-order quadrupolar shift that scales as $1/\nu_0^2$.
Measuring δcg at several fields and extrapolating to $1/\nu_0^2 \to 0$ removes it,
giving the true isotropic chemical shift **δiso** and the quadrupolar coupling
$C_Q$ (Sandland *et al.* 2004, Eq. 1; Baasner *et al.* 2014, Fig. 6):

$$\delta_\text{cg} = \delta_\text{iso} - \frac{10^6}{40}\,\frac{C_Q^2(3+\eta^2)}{\nu_0^2\,I^2(2I-1)^2}\left(I(I+1)-\frac{3}{4}\right)$$

so a plot of δcg (ppm) vs $1/\nu_0^2$ is a straight line: the **intercept is
δiso**, and the **slope gives $C_Q$** (with an assumed η, conventionally 0.7 —
two centres of gravity cannot determine η).

**Tools ▸ QCPMG: infinite-field δiso** opens the extrapolation. Enter each
field's Larmor frequency and its δcg (type it, or **grab it from the open
spectrum's visible range** — zoom to the CT band first), tick whether that field
was CT-selective, set η, and **Compute**. It plots δcg vs $1/\nu_0^2$ with the
fit line and reports δiso, $C_Q$, and $P_Q$ with propagated uncertainties.

> **Selective vs non-selective pulses.** Equation (1) is the shift of the
> *central-transition* centre of gravity. In the **large-$C_Q$ limit**
> ($C_Q \gtrsim 1.5$ MHz) only the ½ ↔ −½ transition is excited even by a
> non-selective (hard) pulse — the satellites are too broad — so the measured
> centroid is the CT centroid at **both** fields regardless of pulse
> selectivity, and a non-selective field can be combined with a CT-selective one
> (Baasner *et al.* 2014). Take the centroid over the **CT band only** at each
> field. The CT-selective flag is recorded for provenance; it does not change
> Eq. (1) in this limit.

A companion **two-field width separation** (Sandland Eq. 2) splits the linewidth
into a quadrupolar part $W_q \propto 1/B_0^2$ and a chemical-shift-distribution
part $W_\text{csd} \propto B_0$; ask if you want that wired in too.

---

## 5 · Background

The echo train trades acquisition time for S/N; you recover the **true**
lineshape by coadding the echoes (sum echo), while the spikelet view is a
convenient, high-S/N but lineshape-sparse alternative.

## References

- F. H. Larsen, H. J. Jakobsen, P. D. Ellis, N. C. Nielsen,
  "Sensitivity-enhanced quadrupolar-echo NMR of half-integer quadrupolar nuclei",
  *J. Phys. Chem. A* **101**, 8597 (1997). *(the original QCPMG)*
- F. H. Larsen, H. J. Jakobsen, P. D. Ellis, N. C. Nielsen, "QCPMG-MAS NMR of
  half-integer quadrupolar nuclei", *J. Magn. Reson.* **131**, 144 (1998).
- H. Y. Carr, E. M. Purcell, *Phys. Rev.* **94**, 630 (1954); S. Meiboom, D. Gill,
  *Rev. Sci. Instrum.* **29**, 688 (1958). *(the CPMG echo train)*
- S. G. J. van Meerten *et al.*, ssNake, *J. Magn. Reson.* **301**, 56 (2019), and
  its QCPMG tutorial — on which this workflow is modelled.
- T. O. Sandland, L.-S. Du, J. F. Stebbins, J. D. Webster, "Structure of Cl-
  containing silicate and aluminosilicate glasses: A ³⁵Cl MAS-NMR study",
  *Geochim. Cosmochim. Acta* **68**, 5059 (2004). *(infinite-field δiso, Eq. 1–2)*
- J. F. Stebbins, L.-S. Du, "Chloride ion sites in silicate and aluminosilicate
  glasses: A preliminary study by ³⁵Cl solid-state NMR", *Am. Mineral.* **87**,
  359 (2002).
- A. Baasner *et al.*, "The behavior of chlorine in aluminosilicate glasses",
  *Geochim. Cosmochim. Acta* (2014); δcg vs $1/\nu_0^2$ extrapolation, Fig. 6.
- H. Schmidt *et al.* (2000) and D. Freude, J. Haase, "Quadrupole effects in
  solid-state NMR", *NMR Basic Principles and Progress* **29** (1993). *(the CT
  second-order shift and width formulas)*

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
