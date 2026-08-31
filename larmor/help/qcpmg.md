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
| **Sum echo, magnitude (mc)** | as above, then \|spectrum\| | high | ❌ no — but the only option when the pattern cannot be phased |
| **Both, overlaid** | Stage 5, tick both | — | the validation: the envelope must trace the spikelet tops |

The spikelet **maxima** trace the powder pattern, but a smooth model cannot fit
the comb — fit the **sum echo** and use spikelets only to inspect the manifold.

---

## 1 · The two spectra

**Sum echo.** Split the train into its echoes, add them (optionally weighted by
the T₂ decay), and process the single resulting echo with *whole-echo*
processing: the echo top is moved to t = 0 so the transform is **pure
absorption**. This is the spectrum you fit with the usual quadrupolar models.

**Spikelets.** Fourier-transform the whole train untouched. The manifold of
sharp lines, spaced by 1/τ_echo, traces the powder pattern with superb S/N —
but there is no information *between* the spikes, so no smooth lineshape model
can fit it. Use it to check that the sum-echo envelope really does follow the
intensity distribution: **Stage 5 draws both on one axes** for exactly that
comparison.

---

## 2 · The six stages

Each stage shows one plot and gives you one number. The headline readout and
**Send to fit →** stay visible from every stage.

> **Open the raw `fid`**, not the EXPNO folder — an EXPNO resolves to the
> *processed* data, which is no longer an echo train. LARMOR picks the `fid`
> for you if you hand it a folder.

### 1 · Train & split

The echo period is **read from the pulse program** — not guessed. Sequences
disagree on where they put it, so LARMOR checks the lot: Bruker's `CNST7`
(spikelet spacing, Hz) and `CNST8` (points), then the NMRFAM/Perras
`qcpmg.av4.nmrfam` trio `CNST11` (Hz), `CNST15` (echo period, µs) and
`CNST14` (points). The readout names the one it used. Failing all of those it
falls back to a rotor-synchronisation guess or the autocorrelation, and says
so — **Find period** then measures it from the data itself (echo-repeat
correlation), which is the reliable answer when nothing was recorded.

You can type the period as **points or Hz** — they stay in sync. Two health
scores sit in the readout: **echo-repeat** (does the train actually repeat at
this period?) and **alignment**. Both collapse when the period is wrong.

#### Split offset

Acquisition does not always begin half an echo before the first top. The
NMRFAM sequence starts recording **at** a top, so the natural blocks each
hold the *right half of one echo and the left half of the next* — two
different echoes, of different amplitude, glued into a fake one. LARMOR
detects this on load and skips the right number of points so the echo lands
in the middle of its block; **Centre echo** recomputes it, and the field is
editable.

The cost of getting this wrong is not cosmetic. On a real ³⁵Cl train the
straddled split gave **T₂ = 4.0 ms instead of 10.3 ms**, demanded p1 = 493°
and p2 = 331° to phase (against 0° and 0° once centred), and inflated the
FWHM by 24 %.

### 2 · Echo & top

All echoes overlaid, plus **first vs last** — the ssNake validation that the
split is right: their features (and the flat tail) must line up.

**Drag the vertical line onto the echo top** (this is ssNake's "Pos *N*"); the
numeric field follows, and vice versa. A dotted line marks the block centre as
a **tolerance marker**: for whole-echo processing the top should be *within a
point or two* of it.

⚠️ **Put the marker on the measured echo maximum, not on the centre line.** For
an odd echo length the true top is typically `centre + 1`, and this is the most
expensive mistake available in this dialog — moving the top by a single point
changed the fitted T₂ by 4 % to 7400 % (median 39 %) across a real 12-sample
set. The **Auto** button uses the coherent average of all echoes, which found
the published top on all twelve; trust it unless the marker visibly misses.

That one point drives three things — where the decay is sampled, where the
whole-echo swap happens, and therefore the phase of the final spectrum.

### 3 · Decay & T₂

The echo-top intensity versus echo number, fitted with `C + B·exp(-t/T₂)`.
**Click any point to exclude it**; it refits immediately.

- **Signed real** (default) is what ssNake samples and what reproduces
  published values. *Magnitude* has a rectified noise floor that biases the
  tail upward.
- **Fit a constant offset** (default on): the decay sits on a noise floor.
  Across a real 12-sample set `B/C` ranged 10–77 (median 18), and dropping
  `C` moved T₂ by 6–67 % (median 26 %).
- You get **T₂ ± uncertainty and R²** — so a decay that is not actually
  measurable says so instead of returning a confident number. A sample whose
  signal dies within one or two echoes will show a large error and a warning;
  that is a real limit of the experiment, not a fit you should tune.

⚠️ **Two time axes.** ssNake's `Split` copies the *acquisition* sweep width
onto the echo dimension, so a T₂ read there is in "one dwell per echo" units,
**not seconds**:

> `T₂(physical) = T₂(ssNake D1) × points-per-echo`

For a 293-point echo, ssNake's `1.446e-05 s` **is** `4.24 ms`, and its matched
`22013 Hz` is physically `75 Hz`. Both are self-consistent; only the physical
pair transfers to another program. LARMOR shows **both**, always labelled, and
the CSV export names them `T2_physical_s` and `T2_ssnake_D1_s` — never a bare
`T2`.

### 4 · Apodization

**Matched filter**: one click applies `LB = 1/(π·T₂)`, with the value on the
button. Weighting each echo by `exp(-t/T₂)` before summing is the matched
filter that maximises S/N — echoes that are mostly noise contribute
proportionately less.

You see the **apodized echoes** and the **weighting curve along the echo
dimension** (ssNake's "Apodised echoes" / "Apodised D1"), plus how many
echoes effectively survive. Toggling the weighting no longer rescales the
spectrum, so before/after are directly comparable.

### 5 · Spectrum

Sum echo and spikelets, **independently toggleable** — tick both for the
overlay. A correct whole-echo transform is already near-pure absorption, so
**p0 alone is normally enough**.

A large p1 is a *diagnostic*, not a nuisance. An echo top that falls between
two samples needs a first-order phase to compensate it exactly, and up to
±180° is legitimate (a −0.76-dwell top on a real ⁸¹Br train needed +148°).
Beyond about ±200°, the period or the top is genuinely wrong: go back to
stage 1 (**Find period**) and stage 2 (**Auto**).

#### When p0/p1 is not enough: the second-order phase

A **frequency-swept refocusing pulse** — WURST or chirp, as in WURST-CPMG
(`WCPMG`) — imprints a *quadratic* phase across the swept band. That is
invisible to p0 and p1, so a spectrum that looks hopeless under ordinary
phasing is usually not: it needs a **p2** term.

On a real ⁸¹Br WCPMG dataset (2 MHz sweep, pattern ~400 kHz wide) the best
p0/p1 phasing still left −47 % negative dips and a meaningless δ_CG of
−94 ppm. Adding p2 takes the dips to **−3.8 %** and gives δ_CG = −314 ppm.
**Autophase fits p2 automatically** and keeps it only when it genuinely
helps, so an ordinary echo still reports p2 = 0.

#### Magnitude (mc)

Ticking **magnitude (mc)** plots |spectrum| — TopSpin's `mc` — which is
phase-independent by construction. Use it when a phase error cannot be
written as a polynomial at all, or as a cross-check: on the ⁸¹Br sample above
magnitude gives δ_CG = −310 ppm against the p2-phased −314 ppm, which is the
agreement that makes both trustworthy.

Two things to know, and one common myth to drop:

- **Magnitude does *not* cost you ×√3 in width here.** That factor applies to
  a *causal* (one-sided) FID, where absorption and dispersion are Hilbert
  partners. A whole echo is symmetric about t = 0, so its ideal spectrum is
  **real** and |spectrum| simply recovers the absorption lineshape — measured
  widening 1.000 through LARMOR's own transform, against 1.70 for the same
  linewidth from a one-sided FID.
- **The rectified noise floor is subtracted**, estimated from the two *edges*
  of the spectrum (never the median of the whole trace, which stops being a
  floor once a wide pattern fills the window — that mistake moved a real
  δ_CG by 44 ppm).
- **The sign is what you lose.** |spectrum| ≥ 0 before the floor subtraction,
  so a genuinely negative feature folds upward and cannot be recognised.

The ppm axis uses the **processing reference** (`SF` from `procs`), i.e. the
same zero TopSpin puts on the axis. If `procs` is missing, the readout warns
that the carrier fell back to `O1/BF1` — on a referenced dataset the two can
differ by tens of ppm, which would shift every shift you report.

### 6 · Measure — δ_CG and FWHM

For a broad, distribution-dominated pattern measured at a **single field**, a
lineshape fit cannot separate δ_iso from the second-order quadrupolar shift
(both are distributed and correlated). The defensible numbers are the
**centre of gravity** of the central band and its **width**.

Drag the window onto the first intensity minima either side of the central
peak. LARMOR reports **δ_CG ± σ**, where σ is the spread obtained by jittering
each window edge — a deterministic replacement for integrating three times by
hand. Read σ as a **quality flag**: a few ppm means the window is well
defined; tens of ppm means the edges are running down a tail and should be
placed by hand (the dialog says so).

---

## 3 · Recommended workflow

1. **Open the raw `fid`** of the QCPMG EXPNO.
2. **Stage 1** — confirm the period was read from the pulse program and the
   markers sit on the echoes.
3. **Stage 2** — check first-vs-last overlap; put the top marker on the echo
   maximum (it is usually already there).
4. **Stage 3** — read T₂ ± error and R². Exclude obvious outliers by clicking.
5. **Stage 4** — click **Use matched LB**.
6. **Stage 5** — **Autophase**, then tick *spikelets* to confirm the envelope
   traces the spikelet tops.
7. **Stage 6** — place the window, read δ_CG ± σ and FWHM.
8. **Copy CSV** for the lab book, then **Send to fit →** to model the
   lineshape (see the *Lineshapes* manual — for a glass, `czjzek` or
   `ext_czjzek`).

> Fitting a single-field QCPMG spectrum gives correlated δ_iso/C_Q. Quote
> δ_CG and the central-band width as the primary numbers, and treat the fit as
> supporting information — or measure at a second field (§4).

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

**Two-field width split (Sandland Eq. 2).** Fill the **FWHM (ppm)** column at
both fields and press **Split W_q / W_csd**: it separates the CT linewidth into a
**quadrupolar** part $W_q \propto 1/\nu_0^2$ (broader at low field) and a
**chemical-shift-distribution** part $W_\text{csd}$ (field-independent in ppm),

$$\text{FWHM}_1^2 = W_q^2 + W_\text{csd}^2, \qquad \text{FWHM}_2^2 = W_q^2\left(\frac{\nu_1}{\nu_2}\right)^4 + W_\text{csd}^2$$

so $W_\text{csd}$ reports the intrinsic shift disorder of the site independent of
the quadrupolar broadening.

---


### Many samples at once

**Tools ▸ QCPMG: batch infinite-field δiso…** opens a grid: set how many
samples and how many fields, then **drop the processed spectra straight onto
the cells** (the `.csv` files stage 5 writes with *Save as dataset…*), or
double-click a cell to browse. Dropping several files at once fills a row
from that column onwards.

Each cell is measured exactly as the single-sample dialog does it — δcg over
an automatic window — and the column header learns its field from the files
themselves, warning if the frequencies in one column disagree. **Select any
cell** to see its spectrum with a draggable band and supervise that one
measurement; the fit is invalidated whenever you move a window, so a stale
result can never be exported.

**Compute all** extrapolates every sample. Then:

- **Export report…** writes a plain-text record: every input point, every
  fitted δiso, C_Q and P_Q with uncertainties, the W_q/W_csd split where two
  fields allow it, and the assumptions (η, spin) spelled out. A sample that
  could not be fitted is listed as such rather than silently dropped.
- **Export figures…** writes the **merged** figure — every sample on one
  δcg vs 1/ν₀² axes, each with its extrapolation and a starred intercept —
  **and one figure per sample**, all as `.png` + `.svg` + `.pdf` at 600 dpi.

The single-field-pair dialog has the same two export buttons for one sample.
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
