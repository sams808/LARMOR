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
- **Gaussian (GM)** — $w(t) = \exp\!\left[\pi\,\text{LB}\,t - (\pi\,\text{LB}\,t)^2/(2\cdot\text{GB}\,T)\right]$,
  resolution enhancement (negative `LB`, `GB` sets the Gaussian maximum).
- **Sine / Sine² (SINE/QSINE)** — $w(t) = \sin\!\left(\pi(1-\phi)\,t/T + \phi\right)$
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
$S_\text{corr}(\nu) = S(\nu)\,e^{i(\phi_0 + \phi_1(\nu-\nu_\text{pivot})/\text{SW})}$:
a frequency-independent **p0** and a linear **p1**. The panel gives sliders, an
exact type-in, and **−90 / +90 / 180°** quick steps (two +90 equals one 180).
**TopSpin-style pivot:** while the Processing panel is open a **draggable purple
pivot line** appears (default: the tallest peak); p1 rotates *about* the pivot,
so the peak under it stays in phase and only the wings twist — drag the pivot
onto whichever peak you want to hold. **Autophase** minimises the spectral
entropy of the real part (the **ACME** criterion of Chen *et al.* 2002), finding
p0 **and** p1 robustly even on crowded spectra.

**Drag to phase (TopSpin gesture):** press *Drag to phase* in the panel (or
**Ctrl+P**, *Process ▸ Drag to phase*) and drag on the spectrum — left/right
turns p0 (0.25°/px), up/down turns p1 (1°/px, up = positive) about the pivot; a
drag locks to whichever direction it starts in, **Shift** makes it fine (×0.1),
**Ctrl+drag** (or the middle button) still pans, and **Esc** leaves the mode.
Start the drag on empty canvas: the pivot, paddles, rulers, anchors and zones
keep their own drags. Moving the pivot while phasing does not change the
spectrum (p0 is re-expressed), and one Undo takes back a whole drag. *Hilbert
first* is switched on automatically in pdata mode: without the reconstructed
imaginary part a phase change would only scale the spectrum. While the FID is
displayed (Ctrl+T) the gesture is suspended and a drag pans.

### Back to the FID (time ↔ frequency) and the imaginary channel

**FID ⇄ spectrum** (*Process ▸ FID ⇄ spectrum*, **Ctrl+T**, the sidebar **FID**
button, or the button in the panel's **Display** row) flips the canvas to the
**windowed, zero-filled FID the transform sees** — the state of the pipeline just
before its last `ft` — on a non-inverted time axis in ms. Model, components,
residual, paddles, zones, overlays, pivot and literature ranges are hidden until
you return; clicks on the canvas place nothing, the cursor reads ms, and **Fit** is
refused until the spectrum is back (its window is read from the frequency axis).
Every **WDW / LB / GB / SSB / TDeff / ZF** change redraws the FID live, so an echo
train or a truncated signal is re-apodized while watching the decay; **Ctrl+T**
again returns to the spectrum with the new lineshape and the previous zoom.
Nothing is reloaded: the chain is re-applied from the unprocessed data.

**Display channel.** The **real / imag / |S|** radios in the same row (*Process ▸
Display channel*, **Ctrl+I** cycles) draw one channel of the complex result on the
same ppm axis: inspect the **imaginary channel while phasing** — under the pivot
the dispersion signal should vanish when p0 / p1 are right — or the magnitude.
The channel is display only: the fit, S/N, save and overlays always use the real
frequency-domain spectrum, and **|S|** is not the destructive `magnitude`
pipeline step (that is the checkbox next to SR). A spectrum that came in
real-only (a TopSpin 1r, a CSV) has an identically zero imaginary channel, so
selecting imag or |S| ticks **Hilbert first** automatically (the real part is
unchanged by it).

**Re-apodizing a spectrum that did not come from a raw fid.** The raw-fid mode
re-applies the window to the instrument file and is exact. For a TopSpin 1r, a
CSV, or a spectrum sent from the FID / QCPMG / VOCS dialogs, the toggle ticks
**re-apodize this spectrum** in the window block, which builds
`Hilbert → IFT → window → FT` ahead of the phase steps. Hilbert first is mandatory
and stays locked: the inverse transform of a real-only spectrum is two-sided
(hermitian), and a one-sided window on it loses half the signal and distorts the
line — the whole-echo trap in a new guise. The result is a **reconstruction**: it
assumes a well-phased, flat-baseline spectrum, and the new window compounds with
the one already applied in TopSpin, so for a genuine re-processing use the
raw-fid mode or *File ▸ Open FID*. Open FID records its own chain in the recipe,
so after *Use this spectrum* the first **Ctrl+T** shows the TRUE windowed FID from
the instrument file; a raw `fid` opened directly does the same with its preview
chain (`fcor`, EM 100 Hz, `ft`, `magnitude` — untick *magnitude* and phase to turn
the preview into a working spectrum).

### Baseline

- **Automatic** — an asymmetrically reweighted penalized least-squares baseline
  (**arPLS**, Baek *et al.* 2015): it iteratively fits a smooth curve
  ($\lambda \approx 10^7$) that follows the baseline but not the peaks. Robust
  for rolling baselines under broad lines.
- **Manual anchors (dmfit-style)** — click **Pick anchors**, then click the
  spectrum to drop as many baseline points as you like (drag to shape); a live
  PCHIP curve previews the baseline, and it is **subtracted automatically when
  you turn Pick anchors back off** — no separate Subtract click.
- **Iterative (dead-time; Yon et al. 2020)** — for the rolling baseline from
  receiver **dead time** in pulse-acquire MAS, where the above fail. Iterative,
  histogram-thresholded smoothing baseline, optionally restricted to broad
  (dead-time) components. See the **Processing reference** manual.

### Referencing (SR / Calibrate)

Type a spectral-reference **SR** (Hz), or **Process ▸ Calibrate**, click a peak
(it snaps to the local maximum), and set its known ppm — LARMOR reports the
resulting SR. Double-click the **experiment strip** to edit nucleus / field /
νrot / SR, or copy the SR from another spectrum. For a spectrum read from an
EXPNO the same dialog shows where νrot came from — `acqus` MASR, the title, the
NMRFAM booking sidecar `experiment_addenda.xml`, or a rate measured from the
spinning sidebands — with a **Use** button per source, and **Remember for this
session** applies the confirmed rate to the session's other spectra with the
same source values (see *Getting started*, MAS-rate check). All referencing is
a rigid ppm shift of the axis; the raw data is untouched.

**Referencing audit (Tools ▸ Referencing audit…).** Checks a whole session at
once. Point it at one month folder of the data tree (`…/DATA/2026-05`): it lists
the ¹H spectra that carry a reference (adamantane at 1.82 ppm by default;
the tallest line of each is checked against that value), then compares every
other acquisition's stored SR with the value indirect referencing gives,
SF_X = SF_¹H · Ξ_X / Ξ_¹H (IUPAC unified scale, the same rule as TopSpin's
`xiref`). Verdicts: **ok**, **off** (a stale SR carried over from another day),
**unreferenced** (SR = 0, `xiref` was never applied), **no reference** (no
referenced ¹H on that magnet in the session). References are never borrowed
from another month. The window exports a CSV plus a TopSpin-ready list of `sr`
values per EXPNO, appends every old and new value to a permanent log
(`%LOCALAPPDATA%\LARMOR\referencing_log.jsonl`, never truncated) so a
correction typed at the spectrometer can be reversed later, and can
re-reference the spectrum open in the workbench (the recipe keeps the old SR
in its provenance). It is meant to be run once on a new dataset; a repeat
run says when the session was audited before. The same check runs from the
command line: `larmor srcheck <session folder> --csv audit.csv`.

> **Processing history.** *Process ▸ Processing steps* lists every applied op;
> remove any one and LARMOR re-applies the reduced pipeline. Removing `ift` alone
> breaks a re-apodize chain (`hilbert, ift, em, ft`) — remove the window step
> instead. The full op list is in the **Processing reference** manual.

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
   - a **pin** ☑ to fix the parameter (a pinned value reads dimmed, like
     dmfit's greyed parameters; the Czjzek family's `lb` starts pinned at its
     model default and is freed by unticking it).
   **Right-click a row** (its letter, model or family cell, or any value) for
   the line actions: **Add spinning sidebands…** (the Decomposition dialog
   with this line preselected — linked copies at ±k·νrot), **Duplicate**,
   **Rename…**, **Hide / Show on plot**, **Fix all / Free all parameters of
   the line**, **Move up / down**, **Remove**, and the **Family ▸** presets.
   Select several rows (Ctrl- or Shift-click their letter cells) and
   **Remove N selected lines** — or the **Delete** key — takes them all in
   one undo step; links between the remaining lines are renumbered and links
   to a removed line dropped, never left dangling. The value cells keep
   their own menus (Constrain min / max…, links, Unlink).
3. **Fit** (F5). It minimises $\chi^2=\sum_i w_i\,[y_i - f(x_i)]^2$ by
   Levenberg–Marquardt (via **lmfit**), and every fitted value comes back with a
   **standard error** from the covariance matrix. Read **RMSD** and **S/N** next
   to the buttons.
4. **Report** (F6) gives the **quantification**: each site's integrated area as a
   population **% ± error**. Use the area-normalised `gl_norm` model (or a Czjzek
   site) so an amplitude *is* a population.
   Lines that belong to one structural species can be tagged in the **family**
   column of the Fit-Parameters table (the last column; right-click a row ▸
   **Family** for the presets of the nucleus — BO3/BO4, Al(IV)/Al(V)/Al(VI),
   Q0–Q4 — or type any name; **Decomposition ▸ Label lines from literature
   ranges** fills it for ²⁷Al and ¹¹B). The Report then adds a bold **Σ** row
   per family and the named ratios of the nucleus (N4 = BO4/(BO3+BO4), ⟨CN⟩ Al,
   ⟨n⟩); their uncertainty propagates the covariance between the line
   amplitudes after a fit, or the spread of the Monte-Carlo trials after
   **Use as fit errors**, and the row tooltip states which basis was used —
   *covariance*, *Monte-Carlo* or, when neither is available, *independent*
   (see **Fitting glasses for publication** §5).
   ⚠ For **quadrupolar** nuclei, integrated areas only reflect populations if
   the *acquisition* was quantitative — CT-selective excitation, or the
   short-pulse regime (flip angle ≤ 30°/(I+½), i.e. ≤ 10° for I = 5/2), plus
   full relaxation. No fit can repair a non-quantitative acquisition; see
   **Fitting glasses for publication** §6 for the exact conditions.

The Report dock's **Copy methods** gives the full Experimental paragraph —
spectrometer and field, probe, MAS rate, pulse and flip angle, recycle delay,
scans, referencing, TopSpin and LARMOR processing, the fit and the software
versions — built from acqus / procs / title (a CSV or dmfit source gets the fit
sentence alone); *Process ▸ Experiment parameters…* shows, read-only, what was
read from those files.

### Reading the fit-health strip

Under the spectrum, a one-line strip states whether the last fit's numbers can
be read as they stand. It is filled after every **Fit** and **Auto Fit** and
re-evaluated as values change: a coloured pill, the RMSD and $\chi^2_r$ of the
fit, then one chip per flag. A click on a chip opens the matching detail;
**F7** (Decomposition ▸ Fit health details…) or a click on the pill lists every
flag together with the analysis tools. The colours classify the *kind* of
flag, not the quality of the fit — LARMOR does not grade a fit, and a fit
that is merely not accounting for sidebands or a phasing error is never
called unphysical.

| Pill | Meaning |
|---|---|
| `no fit yet — F5 fits the lines` | the model has not been fitted; only the values are checked while lines are placed |
| `✓ Fit OK` (green) | every check passed — the pill's tooltip lists them, and names any acquisition fact that could not be judged |
| `⚠ Fit: check` (amber) | something to look at: a structured residual, a parameter at a bound, a degenerate pair, a population the data do not support, a line outside the window, a short recycle delay or a long pulse. The values are readable as they stand |
| `✗ Fit: not physical` (red) | a value that cannot be physical — a negative amplitude, η outside [0, 1], a width of zero, a C_Q above 120 MHz. Only these turn the pill red |

| Chip | What it measures | A click opens |
|---|---|---|
| `residual N.N× noise[, structured]` / `residual: structure left` | residual RMS in the signal region ÷ the noise of the spectrum edges, at or above 1.5×, and/or a Wald–Wolfowitz runs test with $\vert z\vert > 3$ or a lag-1 autocorrelation above 0.4. The tooltip names the usual causes — spinning sidebands not modelled, a phasing or baseline error, one component too few | the residual trace (View ▸ Residual) |
| `unphysical ×n` (red) | η or the Gauss/Lorentz mix outside [0, 1], a width ≤ 0, a negative amplitude or C_Q, a C_Q above 120 MHz | the offending cell of the Fit-Parameters table |
| `δiso outside the fit window: …` | a line whose centre the fit cannot see — widen the window or move the line | its δiso cell |
| `degenerate pair: a↔b (r)` / `degenerate ×n: …` | a parameter pair with $\vert r\vert \geq 0.95$ in the covariance: the values are fine to read, their error bars are not | Parameter correlations |
| `at a bound: A lb at its lower bound 0` | a free parameter that finished at its min/max bound, named as the table names it (letter + column). A parameter pinned since the fit, or held at its model default (the Czjzek lb), is never listed | its cell; the value is marked ‡ in every export |
| `no error bars (no covariance)` | the fit returned no covariance matrix | Errors Analysis (χ² profile) |
| `population ±≥100 %: …` | a population whose relative error reaches 100 % (see *Fitting glasses for publication* §5) | the Report |
| `tail outside window: …` | a line whose simulated area lies more than 2 % outside the integration window — its population is biased low (*Fitting glasses for publication* §5); the Report table's *outside window (%)* column carries the number | widens the window to 99.5 % of every line and re-integrates (F5 then refits over it) |
| `D1 = n T1 → m %` | the recycle delay D1 + AQ against a KNOWN T1 of each site — TopSpin's `ct1t2.txt` of the same-nucleus saturation-recovery EXPNO in the sample folder, or a T1 typed in Experiment parameters (steady-state recovery; 90° assumed until a flip angle is known) — below 99 % | Tools ▸ Relaxation on that EXPNO |
| `flip n° > m° limit (I = …)` | a KNOWN flip angle (from `P1(90)=` or `n deg tip` in the title, or the 90° pulse typed in Experiment parameters) above 30°/(I + ½) for single-pulse spectra of half-integer quadrupolar nuclei | Process ▸ Experiment parameters… |
| `frozen: …` | a site the fit held because its centre lies outside the window | its δiso cell |

An **unknown** acquisition fact — no T1 measurement in the folder, a TopSpin
fit that did not converge, no 90° pulse anywhere — produces **no chip**: it is
not a problem with the fit. The pill's tooltip and the **F7** list still name
it (`recycle 14 s — T1 unknown — not judged`, `flip angle unknown (I = 3/2) —
not judged`), and **Process ▸ Experiment parameters…** is where the fact can
be supplied; once it is, the check runs and a chip appears only if it fails.

After an edit the pill reads **Model** instead of **Fit**, with `edited since
fit`: the residual and value chips follow the live model, while the
covariance-based chips (degenerate, at a bound, no error bars, population, tail
outside window) are dimmed — they describe the last fit until the next **F5**.
The recycle-delay and flip-angle chips describe the acquisition rather than the
fit, so they stay live and appear as soon as lines are placed on a Bruker
dataset. Re-processing the spectrum drops the verdict altogether, since the
covariance no longer describes the data. **View ▸ Panels ▸ Fit health strip**
hides the strip.

**Zones** restrict the fit to chosen spectral regions (union of intervals) — fit
only where the model is valid and let peaks outside float frozen. **Auto Fit**
does a multi-start search to escape local minima.

**Three ways to get errors** (Decomposition ▸ Analyze), in increasing rigour and
cost — the strip's `no error bars` chip opens the second one, the rescue when
the covariance failed:

- **Covariance** — the standard error printed next to every fitted value after a
  **Fit**. Instant, but assumes a locally quadratic, well-conditioned $\chi^2$ —
  it under-reports for strongly correlated or non-linear parameters.
- **Errors Analysis (χ² profile)** — fixes one parameter, re-fits everything else
  at each scan point, and reads the true (possibly asymmetric) confidence interval
  off the real $\chi^2$ curve. One parameter at a time.
- **Monte-Carlo** — a parametric bootstrap (dmfit ▸ Errors ▸ Monte Carlo): adds
  synthetic Gaussian noise at the residual level to the best-fit model and re-fits
  $N$ times (default 200); every free parameter's error is the spread of its
  fitted values, reported as **mean ± σ** with a percentage and a distribution
  **histogram**. It captures correlations *and* non-linearity for all parameters
  at once. Press **Use as fit errors** to write the σ's into the table. It is the
  most honest error bar for a glass/overlapping fit — and the slowest (N refits).

### Which lineshape?

| Model | Use it for |
|---|---|
| **Gauss/Lorentz** (pseudo-Voigt) | symmetric lines, quick fits |
| **Gauss/Lorentz (area)** | as above, amplitude = integral (quantification) |
| **Voigt (true)** | separable Gaussian ⊗ Lorentzian broadening |
| **J-multiplet** | scalar J splitting to *n* equivalent spins |
| **Czjzek / ext. Czjzek** | amorphous quadrupolar sites (glasses) |
| **Czjzek, general d** | dmfit-style d fits; testing whether d = 5 is justified (the GIM ≡ Czjzek at d = 5) |
| **Czjzek + δiso–C_Q correlation** | glasses whose isotropic shift tracks C_Q (tilted MQMAS ridges) |
| **Two-site exchange** | chemical exchange between two isotropic sites (VT series → Arrhenius) |
| **Amorphous** | BO₃ in ¹¹B; a well-defined C_Q with modest Gaussian disorder |
| **Quad CT / 1st / +CSA** | crystalline quadrupolar sites |
| **CSA powder** | spin-½ shielding anisotropy (+ sidebands) |
| **Spectrum (background)** | fit a measured impurity/phase's amplitude & shift |

The full physics, equations and literature for each are in the **Lineshapes —
models & physics** reference (**? ▸ Lineshapes**).

**Spectrum components have two entry points**, both on the **Decomposition**
menu:

- **Add background spectrum…** — a spectrum from disk (an empty rotor, an
  impurity, a separately measured phase).
- **Add a copy of this spectrum…** — the trace **currently on screen**,
  processed exactly as displayed, offered shifted by ±νrot. That is how a
  satellite-transition or spinning-sideband manifold is taken out: the manifold
  repeats the whole pattern at the MAS rate, so a shifted copy of the measured
  spectrum models it with no lineshape to parameterise. The shift is held fixed
  by default (a copy of the same data with a free amplitude *and* a free shift
  is degenerate with the rest of the model), and the copy is a snapshot —
  re-processing the workbench afterwards does not update it.

### Spinning sidebands — detected for you

LARMOR looks for the **repeat of the spectrum at ±ν_rot** whenever a 1D
spectrum lands on the workbench (open, workspace switch, re-processing, a
processed FID, background subtraction, and after the Experiment or Calibrate
dialogs). The spectrum is correlated with a shifted copy of itself; the lag at
which it repeats is the sideband spacing $\Delta\delta = \nu_\mathrm{rot} / \nu_0$
(in ppm), searched within ±2 % of the recorded rate — a rate that is slightly
off is re-measured from the data. The spacing is then checked against the
peaks: the tallest band is the centreband and the orders +1 **and** −1 must
both sit on a resolved peak, the one cheap test that tells a manifold from two
unrelated lines or a two-horn static pattern. When it holds, a small **banner**
appears over the plot ("Spinning sidebands: ν_rot 20 000 Hz (124.6 ppm at
160.46 MHz) · 2 orders each side · confidence 0.90") together with dotted
**guide lines** labelled *centre*, *+1*, *−1*, … on the very peaks it claims,
so the claim can be checked before anything is clicked. Two one-click actions,
each a single undo step:

- **Add linked manifold** (or **Return**) — one centreband line (an existing
  line sitting on it is reused) plus one line per detected order whose position
  is the constraint `A+124.6` / `A-124.6`, whose width and shape are tied to the
  centreband (`A`), and whose **amplitude stays free**: sideband intensities
  follow the shielding or quadrupolar tensor (Herzfeld & Berger), not a fixed
  ratio, so they are seeded from the measured teeth and left to the fit. The
  paddles of linked lines are not draggable — drag the centreband and the
  manifold follows. The button's drop-down offers the same manifold **as one
  `sidebands` model line** instead (geometric ratio *r* and *n_ssb*, both
  seeded from the teeth).
- **Add shifted copy** — one *spectrum* component per detected side, shifted by
  ±Δδ with the shift held and scaled to the first sideband: the copy carries the
  whole pattern, higher orders included (the "Add a copy of this spectrum…"
  route, without the dialog).

LARMOR stays silent on a confirmed-static dataset (ν_rot = 0), when nothing
repeats within ±2 % of a certain recorded rate (so a J-multiplet or two sites
split by about ν_rot are never mistaken for sidebands in the automatic offer),
when the model already carries a manifold (a `sidebands`, CSA-MAS, Czjzek-CSA
or first-order quadrupolar site, linked positions, or a shifted copy), and for
a trace dismissed with **✕** / **Esc** until the data or the rate changes. A
rate recorded as *uncertain* (the red MAS pill) is a question the data can
answer: the search then also scans 1–80 kHz, and accepting an offer **writes
the measured ν_rot** into the recipe and clears the warning. A certain rate is
never written — when the measured and recorded rates disagree by more than
0.5 % the banner shows **⚠ acquisition says …** and the experiment strip is
the place to fix it — except for the `sidebands` model line, which has no
spacing of its own and needs the rate to sit on the teeth.
**Decomposition ▸ Detect spinning sidebands** (Ctrl+Shift+D) runs the same
search on demand, scanning when the recorded rate yields nothing, and reports
the reason in the status bar when it finds no repeat. **View ▸ Offer
spinning-sideband detection on load** turns the automatic offer off.

---

## 4 · Measure & export

- **Tools ▸ Integrals & measurements** — drag regions → integral, %, centre of
  mass, FWHM; **Copy CSV**.
- **File ▸ Copy plot** (with all fitted lines) / **Save plot image** (PNG/SVG)
  for slides and papers.
- **File ▸ Save fit as** — `txt` / `csv` / `json` / **dmfit `.fxmla`**; writes the
  data, model, residual and every component (the dmfit export round-trips the
  Czjzek σ ↔ `sCZ_CQ = 2σ` relation). The `csv` parameter table lists, per
  parameter, `min`, `max`, `link`, `vary` and `at_bound` (`min` / `max` when a
  free value sits at its bound, blank otherwise), so a held or pinned value
  is never mistaken for a fitted one.
- The Report dock's **Copy LaTeX** table and the **Publication bundle…**'s
  `table.tex` mark every held value † (printed without its `± 0.00`), every
  value that finished at a bound ‡ and every linked value §, with one
  footnote line naming each bound and each expression. The status is read
  from the recipe itself (`vary`, `expr`, value against bound), so it is
  right for a saved fit reopened years later and clears itself on the next
  edit — nothing to switch on. The status bar reports the count
  (`marked: 1 fixed · 2 at a bound`).
- **File ▸ Save spectrum as** — a reopenable CSV with a metadata header.
- **Provenance.** A saved recipe carries `acquisition` (the block read from
  acqus / procs / title, with the procno), `software` (the versions that
  produced the fit, stamped at every Fit) and `source_sha256` (the hash of the
  exact data file). On reopen the status bar says *source data changed* when
  the file no longer hashes the same (TopSpin reprocessed it — refit before
  publishing) and *re-referenced since the fit* when the SR moved (the ppm
  axis is no longer the fitted one); a recipe without a hash stays silent.

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
- J. Herzfeld, A. E. Berger, "Sideband intensities in NMR spectra of samples
  spinning at the magic angle", *J. Chem. Phys.* **73**, 6021 (1980).
  *(sideband intensities)*
- M. M. Maricq, J. S. Waugh, "NMR in rotating solids", *J. Chem. Phys.* **70**,
  3300 (1979). *(spinning sidebands)*
- R. R. Ernst, G. Bodenhausen, A. Wokaun, *Principles of NMR in One and Two
  Dimensions*, Oxford (1987). *(general reference)*

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
