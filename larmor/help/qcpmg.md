# QCPMG processing in LARMOR

**QCPMG** (Quadrupolar Carr–Purcell–Meiboom–Gill) records a *train of echoes*
after a single excitation. It is the method of choice for **broad quadrupolar
lineshapes** (e.g. ³⁵Cl, ⁸¹Br, ²⁷Al at low field), where the signal dephases far
faster than it relaxes (T₂). Because many echoes are captured in one shot, the
signal-to-noise ratio is much higher than a single echo.

There are two ways to turn the echo train into a spectrum, and LARMOR gives you
both. **Which one you fit matters** — see below.

---

## 1. The two spectra

### Sum echo (absorption) — *fit this one*
Split the train into individual echoes, add them together, and process the
resulting single echo. You get a **normal, continuous powder lineshape** at high
S/N. This is what you send to the fit workbench and model with `ext Czjzek`,
`Quad CT`, etc.

### Spikelets
Fourier-transform the **whole train** without splitting. You get a manifold of
sharp **spikelets** spaced by `1/τ_echo`, whose *maxima* trace the powder
pattern. All the signal is concentrated in the spikes, so S/N is spectacular —
but **the area between spikelets carries no lineshape information**. A smooth
model cannot fit the comb: if you try to fit `ext Czjzek` directly onto the
spikelet spectrum the result is meaningless. Use spikelets for display / peak
picking, and fit the **sum echo** instead.

---

## 2. The controls, step by step (this mirrors the ssNake tutorial)

**period (points per echo).**
The number of data points between two consecutive echo tops. LARMOR auto-detects
it from the autocorrelation of the magnitude FID; it equals `SW / spikelet
spacing` (e.g. 200 kHz / 2000 Hz = 100 points). It is the number the train is
split by. If the split looks wrong (echoes not aligned), nudge it.

**echo top.**
The point index of the echo maximum *within one period*. It is used for
**whole-echo processing**: the summed echo is circularly shifted so its top sits
at t = 0, which makes the Fourier transform come out in pure absorption (this is
ssNake's *Swap echo* step). Auto-detected from a clean mid-train echo.

**T₂ (echo-top decay → the evolution time).**
The intensity at the echo top, echo by echo, decays as `exp(-t/T₂)` with
`t = k·τ_echo` (`τ_echo = period / SW`). LARMOR fits this and reports **T₂ in ms**
— your transverse *evolution time*. Right next to it is the **matched
apodization** `LB = 1/(π·T₂)` in Hz: apply this Lorentzian to weight the echoes
by how much signal they still carry (see *T₂ weighting*).

**T₂ weighting.**
When ticked, echo *k* is scaled by `exp(-k·τ_echo/T₂)` before summing. This is
the *matched filter*: echoes that have mostly decayed (noise-dominated) count
less, giving the best S/N. It is exactly ssNake's "apply a Lorentzian
`LB = 1/(πT₂)` along the echo dimension".

**GB (Hz).**
Gaussian line broadening applied to the summed echo before the transform —
smooths the lineshape.

**p0 / p1 and step.**
Zero- and first-order phase. The **step** box sets how much the arrows / mouse
wheel move the phase per click — *lower it for fine control* (the sliders felt
too fast; now you type an exact value or step by e.g. 0.1°). **Autophase**
optimizes p0 **and** p1 by minimising the negative area of the real spectrum (a
robust criterion for an all-positive powder pattern). Fine-tune p1 by hand
afterwards if the wings aren't flat.

---

## 3. Recommended workflow

1. **Open the raw fid** of the QCPMG experiment (the echo train).
2. Check the **period** and **echo top** on the echo-train plot (red lines mark
   the period). The auto values are usually right.
3. Read **T₂** off the decay plot. Optionally tick **T₂ weighting**.
4. Choose **sum echo (absorption)** mode.
5. **Autophase**, then nudge **p1** (with a small step) until the lineshape is
   clean absorption.
6. Adjust **GB** to taste.
7. **Send to fit →**, then add a quadrupolar model (`ext Czjzek`, `Quad CT`, …)
   and fit as usual.

Use **spikelets** mode only to inspect the manifold or to confirm the spikelet
spacing; do not fit it directly.

---

## 4. Background

- E. Larsen, H. J. Jakobsen, P. D. Ellis, N. C. Nielsen,
  *Sensitivity-Enhanced Quadrupolar-Echo NMR of Half-Integer Quadrupolar
  Nuclei*, **J. Phys. Chem. A** 101 (1997) 8597 — the original QCPMG.
- F. H. Larsen, H. J. Jakobsen, P. D. Ellis, N. C. Nielsen,
  *QCPMG-MAS NMR of half-integer quadrupolar nuclei*, **J. Magn. Reson.** 131
  (1998) 144.
- ssNake reference manual and the QCPMG tutorial (Radboud University), on which
  this workflow is modelled.

The key idea in all of them: the echo train lets you trade acquisition time for
S/N, and you recover the true lineshape by **coadding the echoes** (sum echo),
while the **spikelet** view is a convenient, high-S/N but lineshape-sparse
alternative.
