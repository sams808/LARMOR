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

> ⚠️ **Changing the period resets the split offset to 0**, clears the echo
> exclusions and resets the echo counts — including when **Find period** sets
> it. Press **Centre echo** afterwards to re-derive the offset.

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

All echoes overlaid (the first **60**, to keep the plot readable — a long
train's remaining echoes only re-trace the same shape), plus **first vs
last** — the ssNake validation that the split is right: their features (and
the flat tail) must line up.

**realign echo tops before summing** (checkbox): each echo is shifted so its
own maximum lands on the shared top before the coherent sum. Use it when a
drifting field or timing jitter visibly staggers the tops in this overlay;
leave it off otherwise — realigning also lines up noise, which slightly
biases the summed echo.

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

**Phasing runs once automatically on load**: on an unphased spectrum the
tallest feature is a noise sliver, and stage 6 would happily report a δ_CG
from it. The automatic pass is exactly the **Autophase** button (p0/p1, then
p0/p1/p2, keeping the quadratic only if it cuts the negative area by more
than 25 %); type into the phase fields or press Autophase again to override
it at any time.

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

**Under MAS the first-minima window is safe only when the pattern is
narrower than the rotor rate.** A QCPMG train acquired while spinning keeps
its spinning sidebands in the sum-echo spectrum (at ±ν_r/ν₀ — 204 ppm at
16 kHz and 78 MHz). For a pattern narrower than ν_r the first minima bracket
a genuine centreband and the window is fine. For a distribution-broadened
pattern wider than ν_r the "first minimum" is a valley *between* sidebands:
the window catches the centreband plus one sideband on one side only, with
opposite sign at the two fields, and δ_CG moves by tens of ppm (a simulated
35Cl Czjzek glass: −129.7 ppm with the first-minima window against −145.6
for the whole manifold at 78 MHz). The multi-field tools (§4) know the rotor
rate of a saved dataset, tick the sideband positions on the supervision plot
and flag *window catches one sideband only*; use the **whole manifold**
window there.

**The window must reach the noise on both sides.** The jitter σ is a local
sensitivity: it does not see a tail that the window cuts. The tools also
re-measure δ_CG with the window at 1.5×, 2× and 3× its width (each edge
alone as well) and report the drift |CG(2w) − CG(w)|; the quoted ± is the
larger of σ and that drift, and *CG not converged — window cuts the pattern*
is flagged when the drift exceeds max(2σ, 5 ppm). The out-of-window median
is subtracted first, so a raw magnitude pedestal cannot fake convergence.

---

**→ infinite-field δiso…** (stage 6) sends this dataset's (field, δ_CG,
FWHM) straight into the multi-field extrapolation dialog — the primary route
into §4 below: process each field's dataset, press this button on each, and
the extrapolation dialog collects them.

**Export figure package…** writes the assembled six-stage figure as
`.png`/`.svg`/`.pdf` next to a base name you choose — the same composite the
right-click *Export figure / Send to studio* menu offers per panel.

**Reading C_Q and η without a fit**: for a static CT pattern with visible
structure, **Tools ▸ Read static pattern (C_Q, η)…** puts three draggable
markers on the spectrum — the two divergent horns and the step-like outer
limit — and reports (C_Q, η, P_Q, δiso) from their positions alone, with the
feature mismatch as its honesty figure. That reading can seed a `quad_ct`
line in one click. It is the classic singularity measurement: defensible at
a single field where a full lineshape fit is not, and a good cross-check on
δ_CG.

**Intensity across the sweep** (WCPMG): the swept pulse also imprints an
*amplitude* weighting — near-flat mid-band, rolling off at the sweep edges.
After sending the spectrum to the workbench, **Process ▸ WURST excitation
profile…** divides the computed WURST-N weighting out (enter the sweep width,
e.g. 2000 kHz; the correction is clamped at a floor so edge noise is not
amplified). Do this before quantifying anything across a wide pattern.

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
   `ext_czjzek`). Tutorial 7 (Help ▸ Tutorials) walks this end to end on a
   static ⁸¹Br WURST-CPMG dataset, through *Read static pattern* and the
   infinite-field δiso reasoning.

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
two centres of gravity cannot determine η). The slope must be **negative**:
the second-order shift lowers δcg more at the lower field. A slope that is
positive, or negative but within 2σ of zero, cannot come from Eq. (1) — the
two δcg are then not the same observable (different processing mode, window
or referencing) — and the tool reports only a 2σ **upper bound** on $C_Q$
instead of a value, together with the slope and its σ.

**Tools ▸ QCPMG: infinite-field δiso** opens the extrapolation. Enter each
field's Larmor frequency and its δcg (type it, or **grab it from the open
spectrum's visible range** — zoom to the CT band first), set η, and
**Compute**. The *CT-selective (declared)* box starts **unknown** and is not
read from the data: selectivity is the operator's judgement (ν_rf against
ν_Q), recorded for provenance only — it does not enter the fit. It plots δcg vs $1/\nu_0^2$ with the
fit line and reports δiso, $C_Q$, and $P_Q$ with propagated uncertainties.

**Add from datasets…** expects the **sum-echo dataset** written by *Save as
dataset…* (its header carries the Larmor frequency, the nucleus, the
processing mode and the rotor rate under `qcpmg_rotor_Hz` — the recipe's
`spin_rate_Hz` stays 0 because a sum-echo spectrum is not re-modelled with
sidebands). A TopSpin `1r` of a QCPMG EXPNO can be opened too, but it
is a **spikelet comb**: the first-minima window would stop at the first
spikelet gap (measured: δcg −105.1 ppm against −112.8 ppm from the sum echo
of the same EXPNO, 8.6 ppm on δiso through the low-field lever). The comb is
detected from the data, the window is seeded from the envelope over one
spikelet period, and the row is flagged — process the `fid` in *Tools ▸
QCPMG* and use *Save as dataset…* / *→ infinite-field δiso…* instead. A
`1r` processed with `mc` (procs `PH_mod = 2`) is recorded as magnitude.

**Window modes (MAS data).** Each dataset row or batch cell has a window
mode: *first minima / dragged band* (the default — supervise it), *whole
manifold* and *centreband*. When the rotor rate is known the supervision
plot ticks the sideband positions of the tallest peak and flags a window
that catches **one sideband only**. The **whole manifold** integrates the
full axis with the edge-noise floor subtracted and repeats the centre of
gravity with the trace cut at 5 / 2 / 1 / 0 % of the peak (their spread is
its σ): by the first-moment theorem it is exact for any distribution of
sites *provided the whole sideband manifold is in the spectrum*. Rectified
noise pulls a full-axis centroid towards the axis centre, so use it on the
absorption sum echo; the report flags *magnitude + whole manifold*. The
**centreband** window (peak ± ν_r/2) is accepted only when it holds ≥ 80 %
of the intensity and its FWHM is below ν_r/2 — a pattern narrower than
ν_r — and is refused with the measured fraction otherwise, because for a
distribution it weights every site by its centreband share, which falls
with P_Q (simulated glass: −89.1 ppm against −70 true).

**Convergence.** Every measurement re-integrates the window at 1.5×, 2× and
3× its width and quotes ± = max(jitter σ, |CG(2w) − CG(w)|); the report
prints the sequence CG(w, 1.5w, 2w, 3w) per field and flags *CG not
converged — window cuts the pattern* (drift beyond max(2σ, 5 ppm) or 10 %
of the separation between the two fields' δcg). On the tutorial's LAW4Ca
dataset at 78 MHz the sequence runs −73.4, −68.1, −57.3, −40.1 ppm while
the jitter σ was 14.5 ppm.

**What the slope measures for a distribution.** For a Czjzek (or any)
distribution of sites the slope gives $\sqrt{\langle P_Q^2\rangle}$ — in
LARMOR's σ convention $2\sqrt{5}\,\sigma$ (= $\sqrt{5}\,\sigma_\text{Cz}$ in
the d'Espinose convention), `larmor.czjzek_dist.rms_pq` — not a single
$C_Q$. Because $C_Q^2(3+\eta^2) = 3P_Q^2$ the slope needs no η at all; only
the $C_Q$ line does. The window must reach the noise on **both** sides: a
wide-window $P_Q$ is a lower bound when the tail runs into the noise (a
Czjzek lineshape fit recovers the tail weight).

> **Selective vs non-selective pulses.** Equation (1) is the shift of the
> *central-transition* centre of gravity. In the **large-$C_Q$ limit**
> ($C_Q \gtrsim 1.5$ MHz) only the ½ ↔ −½ transition is excited even by a
> non-selective (hard) pulse — the satellites are too broad — so the measured
> centroid is the CT centroid at **both** fields regardless of pulse
> selectivity, and a non-selective field can be combined with a CT-selective one
> (Baasner *et al.* 2014). Take the centroid over the **CT band only** at each
> field. The CT-selective flag is recorded for provenance; it does not change
> Eq. (1) in this limit.

**Width split (Sandland Eq. 2).** Fill the **FWHM (ppm)** column at two or
more fields and press **Split W_q / W_csd**: it separates the CT linewidth
into a **quadrupolar** part $W_q \propto 1/\nu_0^2$ (broader at low field) and
a **field-independent** part $W_\text{csd}$ (constant in ppm),

$$\text{FWHM}_i^2 = W_q^2\left(\frac{\nu_\text{ref}}{\nu_i}\right)^4 + W_\text{csd}^2$$

fitted over every field entered (the two-field case is the closed form).
$W_\text{csd}$ collects **everything constant in ppm** — the distribution of
isotropic shifts *and* the chemical-shift anisotropy — so it is an upper
bound on shift disorder unless the CSA is known to be small; δcg itself is
CSA-invariant, so δiso, $C_Q$ and $P_Q$ are unaffected. Eq. 2 assumes
near-Gaussian, disorder-broadened lines: on a pure second-order CT pattern
it returns a spurious $W_\text{csd}$ of 10–17 ppm (static) and a CSA span
comparable to $W_q$ inflates both widths, so the tool flags a split whose
FWHM ratio follows $(\nu_\text{lo}/\nu_\text{hi})^2$ within 15 % or whose
$W_\text{csd}$ is below $0.3\,W_q$. A processing line broadening is removed
in quadrature when it is known. The same separation on the **second
moment** (`qcpmg.second_moment_ppm`, `qcpmg_fields.second_moment_split`)
needs no Gaussian assumption — variances add exactly under convolution —
and is reported as a Gaussian-equivalent FWHM; it is window-sensitive
(the window must hold the whole band and exclude spinning sidebands, which
are fixed in Hz).

---


### Many samples at once

**Tools ▸ QCPMG: batch infinite-field δiso…** opens a grid: set how many
samples and how many fields, then **drop the processed spectra straight onto
the cells** (the `.csv` files stage 5 writes with *Save as dataset…*), or
double-click a cell to browse. Dropping several files at once fills a row
from that column onwards.

Each cell is measured with the same functions and the same automatic window
the single-sample dialog proposes for *Add from datasets…* (the stage-6 send
instead uses the band placed there, so supervise each batch cell), and the
column header learns its field from the files themselves, warning if the
frequencies in one column disagree. A cell whose jitter σ exceeds 8 ppm, whose
centre of gravity drifts as the window widens, or whose window catches one
sideband only shows a **!** flag in the cell and in the report. **Select any
cell** to see its spectrum with a draggable band and supervise that one
measurement; the fit is invalidated whenever you move a window, so a stale
result can never be exported.

**Compute all** extrapolates every sample. Then:

- **Export report…** writes a plain-text record: every input point, every
  fitted δiso, C_Q and P_Q with uncertainties, the W_q/W_csd split over all
  the fields that carry a FWHM, and the assumptions (η, spin) spelled out. A sample that
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
