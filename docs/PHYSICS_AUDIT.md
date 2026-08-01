# LARMOR — Physics audit (viability & trustworthiness for publication)

*Deep review of the physics behind LARMOR, 2026-08. Every formula was checked
against first principles / the primary literature, and the simulation path was
validated numerically against mrsimulator ground truth and known analytic
relations. This document states what is trustworthy, what is an approximation,
and what to check per use-case before publishing numbers from LARMOR.*

---

## Verdict

**LARMOR's physics is sound and its fitted parameters are trustworthy for
publication**, provided the documented caveats (integration window, kernel
resolution/`cq_max`, discrete-2D grid, CSA sign convention) are respected. The
reason is architectural: **the hard spin physics is not re-implemented — it is
delegated to `mrsimulator`**, a peer-reviewed, community-validated spin-dynamics
engine (Srivastava & Grandinetti). LARMOR adds ingestion, the fitting
orchestration, unit conventions, distributions, and processing — and each of
those layers was verified here against theory and numerical ground truth.

Two independent cross-validations anchor this verdict:

| Check | Result |
|---|---|
| Simulated `quad_ct` CT centroid vs the analytic 2nd-order shift formula (`convert.ct_second_order_shift_ppm`) | **agree to 0.01 ppm** |
| LARMOR Czjzek kernel-reweighting vs a *direct* mrsimulator Czjzek ensemble simulation | **0.5 % peak-normalised RMSD** |

and, on real data, LARMOR's fitted MQMAS δiso reproduced dmfit to **~0.5 ppm**
(pCABS2-4 ²⁷Al: 62.2/29.7/−1.1 vs dmfit 62.7/30/−0.35).

---

## 1 · Trust architecture

- **Delegated to mrsimulator (trusted):** all quadrupolar CT/satellite powder
  patterns (`BlochDecayCTSpectrum`/`BlochDecaySpectrum`), the Czjzek and extended
  Czjzek distributions (`CzjzekDistribution`, `ExtCzjzekDistribution`), CSA powder
  patterns and spinning sidebands, and the MQMAS 2D methods (`ThreeQ_VAS`,
  `FiveQ_VAS`, `ST1_VAS`). B₀ is set from ν₀/γ; Cq is passed in Hz; the
  central-transition-only vs full-manifold choice is correct per model.
- **Hand-rolled but verified (this audit):** the unit conversions
  (`convert.py`), the analytic lineshapes (`models/analytic.py`), the
  quantification integral (`quantify.py`), the Fourier transform & phasing
  (`processing.py`/`fourier.py`), the relaxation model forms (`series.py`), the
  QCPMG workflow (`qcpmg.py`), and the two-field δiso/width extrapolation
  (`qcpmg_fields.py`).

The **only** places a physics error could originate are the hand-rolled layer
and the conventions/rescalings LARMOR applies on top of mrsimulator. Those are
exactly what the rest of this document verifies.

---

## 2 · Conventions & unit relations (all verified)

| Quantity | LARMOR formula | Status |
|---|---|---|
| ppm ↔ Hz | Hz = ppm · SFO(MHz) | ✓ |
| Quadrupolar product | P_Q = C_Q·√(1+η²/3) (SOQE) | ✓ standard |
| First-order quad freq. | ν_Q = 3C_Q/[2I(2I−1)] | ✓ standard |
| CT 2nd-order iso shift | δ₂ = −(3/40)·[I(I+1)−¾]/[I²(2I−1)²]·(P_Q/ν₀)²·10⁶ | ✓ = Samoson 1982 / Sandland Eq. 1 |
| EFG → C_Q | C_Q[MHz] = 234.9647·Q[barn]·V_zz[a.u.] | ✓ constant exact |
| Dipolar coupling | d = (μ₀/4π)γ₁γ₂ℏ/r³ | ✓ (¹H–¹H @1.5 Å = 35.6 kHz) |
| Czjzek width ↔ dmfit | dmfit `sCZ_CQ` = 2σ; mode of \|C_Q\| = 2σ | ✓ exact |
| Czjzek invariant | √⟨P_Q²⟩ = √5·σ (Eq. 7) | ✓ exact |
| MQMAS F1 shear factor | c = −(p−R)/(1+R), R=\|C₄(I,p/2)/C₄(I,½)\|; −17/31 for ²⁷Al 3Q | ✓ **measured**, not hardcoded |

**The δ₂ formula is used in two publication-critical places** — MQMAS δiso/PQ
extraction and the QCPMG two-field extrapolation — and it agrees with a full
mrsimulator CT simulation to 0.01 ppm (§1). This is the single most important
formula to get right, and it is right.

---

## 3 · Lineshape models

- **Czjzek / ext-Czjzek** (`models/quadrupolar.py`, `engine.py`): a
  precomputed (C_Q, η) basis of single-site CT subspectra is reweighted by the
  mrsimulator Czjzek PDF — an *exact discretisation* of the Czjzek integral, not
  an approximation of the physics. Validated to 0.5 % RMSD against a direct
  simulation. The isotropic-shift distribution **dCS** is a Gaussian convolution
  (correct), and dCS + round `line` add **in quadrature** in 1D (correct for two
  independent Gaussian broadenings).
- **quad_ct / quad_first / quad_csa** (`models/_singlesite.py`): exact per-site
  mrsimulator simulations (cached by rounded parameters); CT-only vs full
  manifold selected correctly; combined Cq+CSA handled by mrsimulator.
- **CSA (`csa_mas`, `csa_czjzek`)**: mrsimulator shielding (Haeberlen ζ, η).
- **Analytic** (`models/analytic.py`): pseudo-Voigt (correct FWHM→σ and →HWHM
  mappings, peak-normalised), area-normalised `gl_norm` (correct unit-area G and
  L), J-multiplet (binomial weights `comb(n,i)/2ⁿ`, split by J/ν₀), true Voigt
  (`scipy.special.voigt_profile` with correct σ and γ). All verified.
- **`sidebands`**: honestly labelled **empirical** (geometric intensity ratio,
  not Herzfeld–Berger). Use `csa_mas` for physical CSA sidebands; do **not**
  read CSA parameters from the `sidebands` model.

---

## 4 · MQMAS 2D (the subtle part — sound)

- The isotropic (F1) convention is handled by **measuring** the shear/scaling
  factor `c` (and the QIS-axis slope) from reference mrsimulator simulations
  rather than hardcoding — so it is automatically correct for any nucleus and
  coherence order. This is the safest possible design and avoids the classic
  MQMAS-referencing sign errors.
- The kernel F1 axis is rescaled to the **δ1-isotropic** convention, so a
  pure-CS site sits on the diagonal (F1 = δiso). A quadrupolar site's intensity
  then moves along the **QIS axis** (verified: ΔF1 = −0.58·ΔF2 for ²⁷Al 3Q,
  self-consistent with the measured slope). **The fitted δiso is therefore the
  true isotropic chemical shift** — the publication quantity — while C_Q/P_Q come
  from the quadrupolar kernel. dCS smears **along the diagonal** (correct
  shift-disorder physics); the round `line` is isotropic.
- A single fitted **F1 reference offset β** (`mqmas_f1_ref_ppm`) absorbs the
  residual experiment-vs-model F1 referencing; it should refine to a *small*
  value — if it is large, check the F1 referencing of the data.

---

## 5 · Quantification (populations) — correct, with one caveat

`quantify.py` integrates the **actual simulated lineshape** over the fit window
(`np.trapezoid`), *not* the amplitude parameter — so the peak-normalisation of
the models does not bias populations (a broad quadrupolar site correctly counts
more area per unit peak). The reported error is a **stated first-order
approximation** (amplitude covariance only; shape-parameter covariance
neglected), which is printed with the table.

> **Caveat for publication:** the integral is taken over the fit window. A very
> broad quadrupolar tail extending *outside* the window is truncated, biasing
> that site's population low. Set the window (or the full spectrum) wide enough
> to contain the tails before quoting populations. (dmfit has the same
> requirement.)

---

## 6 · Processing

- **Fourier transform** (`processing.op_ft`): `fftshift(fft)` paired with
  `fftshift(fftfreq)` — the handedness is **empirically validated** against real
  Bruker `1r` data (a raw-FID FT peaks at the same ppm as TopSpin's own processed
  spectrum, up to the SR offset). Other vendors' quadrature conventions may need
  a Reverse/Flip (available in the UI).
- **Phasing**: `S·exp(i(φ₀+φ₁(ν−ν_pivot)/SW))` with a user-settable pivot
  (TopSpin-style) — standard and correct. ACME autophase (Chen 2002).
- **Apodization / zero-fill / FCOR**: standard window functions and
  information-preserving zero-fill (see the Processing manual).
- **Referencing (SR/calibrate)**: a rigid ppm shift; re-referencing now also
  moves the fitted sites so the model tracks the peaks.

---

## 7 · Relaxation & QCPMG

- **Relaxation** (`series.py`): saturation recovery `I₀(1−e^(−(t/T₁)^β))`,
  inversion recovery `I₀(1−2f·e^(−(t/T₁)^β))`, CPMG/T₁ρ exponential decay — the
  standard TopSpin `t1/t2` forms; stretched β (KWW) optional. Validated against
  TopSpin (T₁ within ~10 %, r>0.98 on real ²⁷Al satrec).
- **QCPMG** (`qcpmg.py`): the ssNake sum-echo workflow (split → echo-top T₂'
  decay `a·e^(−t/T₂)` → matched filter LB=1/πT₂ → whole-echo → absorption). Fit
  the coadded sum echo, **not** the spikelet comb (documented).
- **Two-field δiso / width** (`qcpmg_fields.py`): Sandland Eq. 1 (δcg vs 1/ν₀²
  → δiso, C_Q) and Eq. 2 (W_q ∝ 1/ν₀² vs W_csd, field-independent). Forward/inverse
  roundtrip exact; propagated errors reproduce Baasner 2014 (±12 ppm δiso, ±0.3
  MHz C_Q). Selective/non-selective fields combine in the large-C_Q limit.

---

## 8 · Approximations & their impact (know these before publishing)

| Approximation | Where | Impact / mitigation |
|---|---|---|
| (C_Q, η) kernel grid (80×11 1D; 40×6 2D) | Czjzek fits | 0.5 % lineshape RMSD; raise in *Computing parameters* for demanding cases |
| `cq_max` (25 MHz 1D / 16 MHz 2D) | Czjzek | truncates the distribution tail for **large σ** — raise `cq_max` if σ is big (mode ≈ 2σ) |
| Discrete quad in 2D snaps to the grid | `quad_ct`/`quad_csa` MQMAS | crystalline 2D C_Q/η is grid-limited (~0.4 MHz); use 1D `quad_ct` (exact) for precise crystalline C_Q, or a finer grid |
| Parameter rounding for the sim cache | all mrsimulator paths | C_Q to 1 kHz, η to 0.001 — negligible vs experimental error |
| 5-point Gaussian for dCS / σζ | Czjzek 2D, csa_czjzek | coarse but adequate for a smooth Gaussian |
| Quantification error = amplitude error only | `quantify.py` | first-order, **stated in output**; use Errors Analysis for the full χ² profile |
| `sidebands` intensities are empirical | `sidebands` model | not for CSA extraction — use `csa_mas` |
| fxmla Czjzek amp ×3.92 to dmfit | export only | calibrated to one ²⁷Al@195 MHz glass; verify for other nuclei/fields (does not affect LARMOR's own fit) |
| CSA ζ uses mrsimulator's *shielding* sign | `csa_mas` | confirm sign vs your convention if you quote ζ |

None of these affect the **core** trustworthiness of a Czjzek or MQMAS glass fit
(the main use case), which is validated end-to-end.

---

## 9 · Recommendations before quoting numbers

1. **Populations:** integrate over a window wide enough to contain every tail;
   report the stated first-order errors, and cross-check important ones with
   *Errors Analysis* (χ² profile) rather than the covariance alone.
2. **Czjzek glasses:** report **σ / √⟨P_Q²⟩ / dCS**, not a fictitious single C_Q;
   check `cq_max` covers the distribution if σ is large.
3. **MQMAS:** confirm the fitted **β (F1 reference)** is small; the fitted δiso is
   the CS shift and C_Q/P_Q the quadrupolar product.
4. **Two-field QCPMG δiso:** state the assumed η (0.7) and that both fields are in
   the large-C_Q (CT-only) limit; the ±errors already reflect the two-field,
   unknown-η uncertainty.
5. **dmfit interop:** LARMOR's own fits and reported numbers are independent of
   the fxmla export factor; only verify the ×3.92 Czjzek amp if you rely on the
   exported file's absolute amplitudes in dmfit.

*Validated by: numerical comparison to mrsimulator ground truth, analytic
first-principles checks, and real-data cross-validation against dmfit and
TopSpin. — LARMOR, Sam Soudani, McCloy group, Washington State University.*
