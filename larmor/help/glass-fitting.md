# Fitting glasses for publication — a working protocol

This page turns the fitting-methodology guidance of **M. Edén,
"Probing oxide-based glass structures by solid-state NMR: Opportunities and
limitations", *J. Magn. Reson. Open* 16–17, 100112 (2023)** — especially its
§8, *NMR spectra deconvolution* — into a concrete LARMOR workflow. Where the
review names a requirement, the LARMOR tool that implements it is given in
**bold**.

Edén's opening warning deserves quoting, because it describes the failure
mode this whole page exists to prevent: a *"worrying recent trend … 
conclusions are drawn from a multi-component deconvolution of an essentially
featureless (near-Gaussian) NMR peak shape, without justifying the physical
and statistical relevance of the best-fit chemical shifts and peak widths."*

## 1 · Choose the components from chemistry, not from the residual

Pick the *number* of peaks and their approximate shifts from what the
composition allows (which Qⁿ species can exist, which coordinations are
plausible — his §5 compiles the shift systematics per nucleus), **before**
looking at fit quality. Adding a component because the residual improves is
exactly the trend the review criticises: with enough Gaussians, any smooth
peak fits perfectly and means nothing.

**View ▸ Literature shift ranges** draws this guidance on the spectrum:
shaded, labeled δiso spans of the typical species for the current nucleus
(²⁷Al⁽⁴/⁵/⁶⁾, ¹¹B BO₃/BO₄, ²⁹Si Qⁿ/Si⁽⁵/⁶⁾, ³¹P Qⁿ, ¹⁷O BO/NBO, ²³Na,
²⁵Mg, ¹⁹F species), with each label/tooltip carrying the corresponding
**P_Q or C_Q range** — a width is not a shift-axis quantity, but it is what
separates e.g. BO₃ (C_Q ≈ 2.4–2.8 MHz) from BO₄ (≈ 0.2–0.8 MHz) where their
shifts approach. The spans are typical multi-component-oxide ranges, an
assignment *guide*, not bounds — e.g. ²⁷Al shifts run significantly lower
in aluminophosphates (the review's own caveat). Every number's source (and
the reference-compound traps for ²³Na/¹⁹F) is documented in
**Help ▸ Literature shift ranges — data & sources**.

## 2 · Peak shape: start free, then justify Gaussian

For spin-½ nuclides (²⁹Si, ³¹P, ¹¹B(IV), ¹³C) in diamagnetic glasses at
B₀ ≥ 5 T, the physical Lorentzian (T₂) contribution is ≲10 % (≲0.5 ppm) of
each peak's FWHM — the rest is a Gaussian shift distribution (§8.3). So:

1. start with a **free Gaussian:Lorentzian ratio** — the `gl` parameter of
   **Gauss/Lorentz**, or the true **Voigt** model;
2. confirm the fitted Lorentzian fraction is small;
3. only then continue with (near-)purely Gaussian shapes.

Related processing rule from the same section: apodize glass spectra with a
**Gaussian window (Process ▸ GM)**, not exponential/EM — Lorentzian
broadening measurably degrades deconvolution accuracy.

## 3 · Restricted ranges — never fixed values

The review's central practical recommendation (§8.3): after one
unconstrained benchmark fit, all further fitting should use **restricted
parameter ranges**, and it stresses that restricted-range fitting *"is to
our knowledge not supported by most currently available public
deconvolution software."* It is first-class in LARMOR — **type `[lo..hi]`
into any fit-table cell**, or use
**Decomposition ▸ Advanced ▸ Restrict around current values (glass
protocol)…** to apply the review's default ranges to every site in one
step:

- **FWHM** — lower bound ≈ **4 ppm** for any amorphous-phase component
  (below that a "peak" in a glass is not physically credible), upper bound
  from the nucleus/Qⁿ systematics;
- **δiso** — confined to about **±3 ppm** around its expected/current
  value;
- **populations (amplitudes)** — always **free**.

Two cautions from the same passage: a constrained fit will normally sit
slightly *above* the unconstrained rms minimum — that is the point, not a
failure; and restricted is not fixed — fixing shifts/widths outright is
called out as something *"never to be utilized"*.

## 4 · Verify starting-point independence

A best fit must not depend on where the optimiser started (§8.3). That is
**Decomposition ▸ Auto Fit** — multi-start refits from randomized starting
points inside the bounds; if restarts land on different answers with
comparable rms, the decomposition is under-determined and the honest next
step is fewer components or better constraints, not a nicer-looking single
run.

## 5 · Statistical relevance, before publishing

The review demands justification of the *statistical* relevance of a
deconvolution. LARMOR's tools map one-to-one:

- **Errors Analysis (χ² profile)** — real 1σ/2σ intervals per parameter,
  not covariance guesses;
- **Monte-Carlo errors** — the parametric bootstrap, robust for correlated
  parameters;
- **identifiability warnings** — parameter pairs the data cannot separate
  are flagged after every fit: the fit-health strip under the spectrum turns
  red with a `degenerate` chip, and a click on it opens the correlation
  matrix;
- **residual diagnostics** — the runs-test flags structured residuals even
  when the rms looks small (a systematically wrong model, e.g. too few or
  wrongly-placed components); the strip shows an amber `structured residual`
  or `residual N× noise` chip, and a click turns the residual trace on.

If a component's population has a 100 % relative uncertainty, the data do
not support that component. Report that, or remove it. The strip shows an
amber `population ±≥100 %` chip for exactly that case, naming the component.

The integral is taken over the window in view, so a broad tail reaching
beyond it is cut off and that population is biased low. The Report table's
*outside window (%)* column gives the share of each line's simulated area
the window leaves out; above 2 % the strip shows an amber `tail outside
window` chip naming the line, and a click widens the window to contain
99.5 % of every line and re-integrates (F5 then refits over the wider
window).

## 6 · Quantitative intensities from quadrupolar nuclei

Peak areas of half-integer quadrupolar nuclei (²⁷Al, ¹¹B, ²³Na, ¹⁷O …) are
only proportional to populations under specific excitation conditions
(§3.6). Either:

- **CT-selective excitation**, ν₁ ≪ νQ for *every* site in the sample
  (typically ν₁ < 10 kHz — beware off-resonance effects); or
- the **short-pulse linear regime** (his Eq. 38):

  τ_pulse ≤ 1 / [12 (I + ½) ν₁]   ⇔   flip angle ≤ 30°/(I + ½)

  i.e. ≤ 15° for I = 3/2 and ≤ 10° for I = 5/2, at high ν₁ (~100 kHz).

If your `zg` used a longer pulse, sites with different C_Q were excited
with different efficiency and the fitted populations are biased — no fit
can repair that afterwards. Also verify **full relaxation** (recycle delay
vs the *longest* T₁ — measure it with **Tools ▸ Relaxation**, don't assume
it) before treating any integral as quantitative.

**What LARMOR checks for you.** After a fit of a Bruker dataset the
fit-health strip reads D1, AQ, P1/PLW1 and the title of the EXPNO
(read-only), finds the saturation-recovery EXPNO of the same nucleus in the
sample folder and takes the T₁ per integral region from TopSpin's
`ct1t2.txt`, matched to each site by its δiso (the nearest region when the
δiso falls outside every integral, marked as such). It reports the
steady-state recovery $(1-E)/(1-E\cos\theta)$ with $E = e^{-(D_1+AQ)/T_1}$,
where θ is the flip angle — 90° when none is known — and $(I+1/2)\,\theta$
replaces θ for a half-integer quadrupolar nucleus (the central-transition
nutation of a short pulse); the chip turns amber below 99 % (≈ 5 T₁ at
90°). The flip angle itself comes from the title (`P1(90)=3.750; 30 deg
tip`, `11 degree tip`) or from the 90° pulse typed once in **Process ▸
Experiment parameters…** (remembered per nucleus, probe and power), and is
compared with 30°/(I + ½) for single-pulse spectra; echo and CPMG
excitation is not judged. Grey chips name what could not be checked (no T₁
measurement in the folder, a TopSpin fit that did not converge, an unknown
flip angle); the **F7** menu's *Measure T1 per site* runs LARMOR's own
per-site decomposition of the relaxation series on the fitted lineshapes
instead. Nothing is corrected — the chips are the caveats to carry into
the Methods text.

## 7 · MQMAS: resolution yes, populations no

3QMAS excitation and 3Q→1Q conversion efficiencies depend on C_Q (§6.2.2,
§6.3.2), so:

- **3QMAS peak volumes are not populations** — quantify from the 1D MAS
  fit, use the 2D for resolving/positioning sites;
- 3QMAS-derived δiso and P_Q values run systematically **lower** than
  1D-fit values — don't mix the two sources in one table without saying so;
- a quick sanity check the review recommends: compare the 3QMAS 1Q ("MAS")
  projection against the quantitative single-pulse spectrum — significant
  disagreement means the 2D under-represents high-C_Q sites.

## 8 · What to report

- B₀ (or Larmor frequency), MAS rate, pulse length/flip angle and recycle
  delay (the quantitativity conditions of §6 above);
- per site: δiso, FWHM, population with **uncertainties** (and how they
  were obtained: covariance / MC / χ² profile);
- for Czjzek sites: **σ and/or √⟨P_Q²⟩ = 2√5·σ, with the convention stated**
  (σ is the width LARMOR stores and fits, half of Czjzek's own σ_Cz; dmfit's
  sCZ_CQ = 2σ and its displayed CQ = 4σ) — see *Lineshapes ▸ Czjzek width
  conventions*;
- the parameter ranges used in the restricted fit (one sentence);
- the **Report dock's "Copy methods"** sentence includes the essentials
  automatically.
