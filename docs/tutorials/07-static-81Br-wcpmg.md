# Tutorial 7 — A static wideline: ⁸¹Br WURST-CPMG from echo train to δiso

*Time: ~40 minutes. Data: one static ⁸¹Br WURST-CPMG acquisition,
`<NMR>/MagLab/DATA/81Br_2026-08/30` (the 0 Ca member of a calcium bromide
glass series; ⁸¹Br at 216.007 MHz, pulse program `WCPMG.jk`, 2 MHz sweep,
TD 58000, SW 2.5 MHz, 512 000 scans; EXPNOs 31–34 are the 1–4 Ca glasses),
and, for §12, a ³⁵Cl sample measured at two fields. `<NMR>` stands for the
instrument-data root of the development machine (`…/Desktop/WSU_work/NMR`).
These are the developer's instrument files: described here,
not part of the repository; any WCPMG or QCPMG `fid` can stand in for
§2–§10. Run the commands from the repository root.*

A static half-integer quadrupolar pattern several hundred kHz wide is the
opposite of the narrow MAS lines of Tutorials 1–6: it is acquired as an echo
train, its shape is a second-order powder pattern smeared by disorder, and a
single-field fit cannot separate the isotropic shift from the quadrupolar
shift. This tutorial processes the train, reads the pattern, fits it, and
ends with what one field can and cannot say about δiso — and how a second
field settles it.

## 1. Why not the processed spectrum

```
larmor info <NMR>/MagLab/DATA/81Br_2026-08/30
```

```
EXPNO: <NMR>\MagLab\DATA\81Br_2026-08\30
nucleus: 81Br   SFO1: 216.0067 MHz
pulse program: WCPMG.jk   TD: 58000   SW: 2500000 Hz
MASR (acqus): 14000.0 Hz
title: 08/10/2026
```
<!-- measured v0.12.1, 2026-09-22: `larmor info` on the EXPNO (the first line echoes the path as given) -->

Two things are wrong with taking this EXPNO's `pdata/1/1r` at face value.
First, TopSpin processed it as a magnitude spikelet spectrum (procs
`PH_mod = 2`, `WDW = 0`): a comb of spikelets 24 ppm apart, not a lineshape
to fit. Second, `MASR = 14000` is stale — the probe was static, the rotor
controller of a static probe is not connected, and the parameter is whatever
the previous MAS session left behind. No `CONFLICT` line is printed because
the title states no rate, but a positive MASR under a WURST/CPMG pulse
program is flagged as a disagreement: opening the EXPNO in the app shows the
red indicator `⚠ MAS 14000 Hz — check!` until the rate is confirmed or set
to 0 in **Process > Experiment parameters…**. The route below sidesteps
both problems — the QCPMG tool works from the raw `fid` and forwards a spin
rate of 0 Hz.

## 2. Stage 1 — train and split

```
larmor desktop
```

**Tools > QCPMG > QCPMG (echo train → spectrum)…** opens *QCPMG processing — guided
sum-echo workflow*, a non-modal window that stays open beside the main one.
**Open echo-train FID…** and pick `30/fid` (with the EXPNO already loaded,
its `fid` is chosen automatically). The six tabs are the six stages.

Tab `1 · Train & split`. The echo period should come from the pulse
program's constants; `WCPMG.jk` records none of the usual ones (CNST7, 8, 11,
14, 15 are all 1 or 0), so the first guess falls back to the stale MASR
(2.5 MHz / 14 kHz = 179 points) — and the data contradict it, so LARMOR
measures the period from the echo-repeat correlation instead. The readout:

```
61 echoes × 475 pts · τecho = 190.0 µs · spikelet spacing 5,263.2 Hz (24.37 ppm) · period MEASURED from the data (echo-repeat correlation) · echo-repeat 0.98 · alignment 0.60 · signal above 3σ through echo 60 (all are kept — the fitted constant absorbs the noise floor)
```
<!-- measured v0.12.1, 2026-09-22: the dialog's stage-1 computation reproduced with larmor.qcpmg (echo_period_from_meta, find_period_by_correlation, split_echoes, split_alignment, period_correlation, n_usable_echoes) on 30/fid -->

**Find period** repeats that measurement; the dotted markers on the train
must sit between echoes. Changing the period resets the split offset, so do
any period work first.

## 3. Stage 2 — echo and top

Tab `2 · Echo & top` overlays the echoes and marks the top and the block
centre. With the period found above the readout is

```
echo top = 283 · block centre = 237 · offset +46 pt — the top is far from the centre; the transform will not be pure absorption
```

The acquisition did not start half an echo before the first top, so each
475-point block holds the tail of one echo and the head of the next.
**Centre echo** (stage 1) corrects this automatically only when the top is
more than a quarter of a period off centre; 46 points is inside that
tolerance, so set the split `offset` box in stage 1 to 46 by hand. The train
re-splits into 60 whole echoes and stage 2 reads

```
echo top = 237 · block centre = 237 · offset +0 pt — whole-echo condition met
```
<!-- measured v0.12.1, 2026-09-22: split_echoes(fid, 475, first=46), echo_top_point, echo_centre -->

**Auto (mean of all echoes)** re-derives the top from the coherent mean of
the echoes; leave `realign echo tops before summing` off unless a fractional
period makes the tops drift across a long train. The first-vs-last overlay
below shows whether the echo shape survives to the end of the train — here
it does not, which stage 3 quantifies.

## 4. Stage 3 — decay and T₂

Tab `3 · Decay & T₂` plots the echo-top intensity against time and fits
`C + B·exp(−t/T₂)` (ssNake's model; `signed real (ssNake)` and `fit a
constant offset` are the defaults):

```
T₂ = 0.282 ± 0.018 ms (R² = 0.9539, 60 of 60 echoes, decays in 1.5 echoes) · matched apodization LB = 1,127.0 Hz = 1/πT₂
```
<!-- measured v0.12.1, 2026-09-22: qcpmg.fit_t2 on the signed-real echo-top decay of the 60 centred echoes -->

The signal is gone after two echoes: this WURST-CPMG train carries almost
all of its information in the first echo, and the other 58 add noise. Click
a point to exclude an outlier (it turns into a cross and the fit is redone
with the remaining points at their true times).

## 5. Stage 4 — apodization

Tab `4 · Apodization`. `T₂-weight the echo sum (matched filter)` is on by
default; press **Use matched LB = 1,127.0 Hz** (the button carries the value
once T₂ is fitted) to set the Lorentzian broadening to 1/(πT₂). The readout
`matched filter ON · LB 1,127.0 Hz · GB 0.0 Hz · effective echoes ≈ 2.0 of
60` says what the weighting does: the sum is dominated by the first two
echoes.

## 6. Stage 5 — the spectrum

Tab `5 · Spectrum` shows the whole-echo transform of the weighted echo sum,
zero-filled ×16. **Autophase** fits p0/p1 first and tries p0/p1/p2 second,
keeping the quadratic term only when it cuts the residual negative area by
more than a quarter — a frequency-swept refocusing pulse can imprint such a
phase. On this train the linear phase suffices:

```
sum echo · scale unit max · zero-fill ×16 · p0 -101.6° p1 141.5° · carrier -0.04 ppm
```
<!-- measured v0.12.1, 2026-09-22: qcpmg.autophase_best on the summed echo; the residual negative area is 0.45 % of full scale with p0/p1, so p2 stays 0 -->

Tick `spikelets` to overlay the spikelet spectrum (the comb TopSpin
produced): its envelope must trace the sum-echo lineshape. Tick
`magnitude (mc)` for the phase-independent cross-check of §7, then untick it
— the absorption spectrum is the one to fit. **Save as dataset…** writes the
processed spectrum as a LARMOR `.csv` with the Larmor frequency in its
header; this is the file the batch route of §12 reads.

## 7. Stage 6 — measure

Tab `6 · Measure`. **Auto window (first minima)** places the integration
window at the first intensity minima either side of the tallest point; the
edges are then jittered by 10 % of the width to give the sensitivity σ:

```
δCG = -238.7 ± 10.7 ppm · FWHM = 95,825 Hz (443.6 ppm) · window 610.3 … -1189.6 ppm
the window is sensitive — drag its edges onto the first minima either side of the central band
```
<!-- measured v0.12.1, 2026-09-22: qcpmg.cg_window, centre_of_gravity (jitter 10 %) and fwhm_hz on the phased sum echo -->

The warning fires because σ exceeds 8 ppm: the pattern has a long tail to
low ppm and no clean minimum on that side, so the automatic edge runs down
the tail. Drag the window's edges by hand and watch δCG; the spread you see
is the honest uncertainty of the centre of gravity. On the magnitude
spectrum the same window logic gives δCG = −237.3 ± 21.0 ppm — the same
number within its (wider) sensitivity, which confirms the phasing. The
headline at the bottom of the window summarises the session
(`T₂ 0.28 ms · LB 1,127 Hz · δCG -238.7 ± 10.7 ppm · FWHM 95,825 Hz`);
**Copy CSV** copies every processing value and the decay points for the lab
book, and **Export figure package…** writes the six-panel figure.

## 8. Send to the workbench, correct the sweep

**Send to fit →** loads the absorption spectrum into the main window (a
warning asks for confirmation if `magnitude (mc)` was left on). The
workspace is titled `08/10/2026 (QCPMG)`, the status bar reads `spectrum
from processed FID loaded — add lines and fit`, and the recipe now carries a
spin rate of 0 Hz — no stale MASR — plus every `qcpmg_*` processing value in
its provenance, so a saved fit records how the spectrum was made.

**View > Axis unit > kHz** shows the axis as an offset in kHz, the natural
unit for a pattern this wide; every internal value stays in ppm.

**Process > Region / algebra > WURST excitation profile…** divides out the amplitude envelope of
the WURST sweep. Set `sweep centre (carrier)` to the carrier, `sweep width`
to 2000 kHz (the pulse program's 2 MHz sweep), `WURST order N` 80 and
`correction floor` 10 %; the status bar confirms `WURST profile divided out
(2000 kHz sweep, N=80) — File ▸ Export ▸ Save spectrum as CSV… to keep it`. On this
dataset the correction changes nothing visible: a WURST-80 profile is flat
over about 90 % of its sweep, and the 2 MHz sweep spans ±4 600 ppm while
the pattern occupies about 1 800 ppm around the carrier. It matters for a
pattern that reaches the sweep edges, or a narrower sweep. The step is
undoable with Ctrl+Z, which restores the spectrum as well as the recipe.

Finally zoom to the pattern (drag on the plot): the fit of §10 uses the
visible range as its window.

## 9. Read the static pattern

**Tools > Read static pattern (C_Q, η)…** opens *Read static pattern —
81Br* with three draggable markers over the spectrum: two solid `horn`
lines and a dashed `edge` line. Put the solid lines on the pattern's two
divergent horns and the dashed line on its step-like outer limit — the side
away from the taller horn; the other limit coincides with a horn and adds
nothing, and placing the edge there is refused with the message
`that edge sits on a horn and adds nothing — mark the opposite (step) limit
of the pattern`. The readout updates live:

```
C_Q = … MHz   η = …   P_Q = … MHz   δiso = … ppm   (feature mismatch … ppm)
```

The inversion runs on position ratios, so it is independent of the axis
reference and scale; C_Q follows from the horn separation and δiso from the
absolute positions; the mismatch is how far the three markers sit from the
positions a single crystalline site with that C_Q and η would give. On a
simulated single site the round trip recovers C_Q to about 3 % and η to
about 0.03. **Add as quad_ct line** seeds a *Quad CT (2nd order)* site
named `read-0` from the reading — status `added quad_ct from the reading:
C_Q … MHz, η …, δiso … ppm — refine with Fit if the lineshape supports it` —
and the step is undoable with Ctrl+Z.

No reading is quoted for this glass, and its pattern shows why: the
WURST-corrected spectrum has a single broad maximum near −184 ppm, falls to
half height at about +45 and −403 ppm and to 20 % at +304 and −821 ppm, with
no resolved horns. A distribution of C_Q and η in the glass has smeared the
crystalline features the reader expects. Placing the horn markers on the
shoulders and the edge on the low-ppm step gives an *effective* single-site
reading, useful as a starting point for a fit and nothing more; the feature
mismatch tells you how poorly a single site describes the shape.
<!-- measured v0.12.1, 2026-09-22: local maxima and outer half-height / 20 % crossings of the smoothed WURST-corrected spectrum -->

## 10. Fit

**Fit** (F5) refines the seeded *Quad CT (2nd order)* line in the visible
window: δiso, C_Q, η, a Gaussian broadening and the amplitude. A static
pattern half a megahertz wide is a heavy simulation, so the first fit builds
its kernel before iterating; **Stop** keeps the latest values. Read the
results strip (`RMSD … · χ²ᵣ …` with the flags of Tutorial 6, §2), the
`± error` column, **Fit > Errors > Parameter correlations…** for the
δiso ↔ C_Q trade-off that a single field never fully breaks, and
**Fit > Report** (F6).

No fitted values are quoted: the fit was not run for this tutorial. On a
glass like this one expect a structured residual under the single-site
model; the physically right description is **Edit > Add line >
Czjzek (quad. distribution)** — the QCPMG manual recommends it for glasses —
and the reading of §9 gives the scale of C_Q to start it from.

## 11. δiso at one field — three estimates, one honest answer

Three numbers now claim to be a shift:

1. **δCG** from stage 6 (−238.7 ± 10.7 ppm here). It is a model-free
   observable, but it contains the second-order quadrupolar shift of the
   centre of gravity, which is negative and scales as C_Q²/ν₀².
2. **δiso from the reading** (§9) and 3. **δiso from the fit** (§10). Both
   remove the quadrupolar shift *through a model* — a single crystalline site
   — so they are only as good as that model is for a glass.

The difference δCG − δiso should equal the second-order shift of the
centre of gravity. **Tools > Conversion tools…** computes it: in the
*Quadrupole (spin I)* group set `spin I` to 1.5, `Larmor ν₀ (MHz)` to
216.007 and type the fitted `Cq (MHz)` and `η`; the line
`PQ = … MHz · νQ = … kHz · CT 2nd-order shift = … ppm` gives the number to
compare with. This is arithmetic on the fitted values, not an independent
measurement. **Fit > Predict at another field…**
(*Target ¹H frequency*: the dialog proposes the present field, 800 MHz;
enter 1100) opens a new workspace with the same model simulated at the
higher field: the pattern narrows and its centre of gravity moves towards
δiso, which is what the next section exploits with real data.

> At one field, quote δCG and the central-band width as the primary results
> and treat any single-field δiso — read or fitted — as supporting
> information conditional on the lineshape model. This is the QCPMG
> manual's recommendation, and for a distribution-broadened glass pattern it
> is the only defensible one.

## 12. With a second field — the model-free route

Measuring δCG at two fields removes the quadrupolar shift without any
lineshape model: δCG = δiso + slope/ν₀², so a straight line through the two
points, plotted against 1/ν₀², has intercept δiso and a slope that gives
C_Q for an assumed η (0.7 by convention; two centres of gravity cannot
determine η). The ⁸¹Br series exists at one field only, so the
demonstration uses a ³⁵Cl chloride glass measured at 78.354 MHz
(`<NMR>/MagLab/DATA/35Cl_2025-12/1`) and at 107.811 MHz
(`<NMR>/NMRFAM/DATA/2026-06_35Cl/06102026_RS40175_LAW3Cl0Ca_SS_ALP/3`); both
are MAS-QCPMG trains (16 kHz and 20 kHz — the NMRFAM `acqus` still says
`MASR = 5000`, a stale controller value the title and `CNST31` correct, and
LARMOR flags the conflict). MAS is fine for the extrapolation only when the
window is a genuine centreband narrower than ν_r, as it is for this 0 Ca
glass (171 and 135 ppm windows against sideband spacings of 204 and 186
ppm), or when the window is the **whole sideband manifold**. For the
distribution-broadened 2–4 Ca glasses the pattern is wider than ν_r, the
first-minima window is one-sided (it catches the centreband plus one
sideband, with opposite sign at the two fields) and δiso shifts by tens of
ppm; the batch grid flags those cells and offers the whole-manifold window.

**From the processing dialog.** Run each EXPNO's `fid` through §2–§7 and
press **→ infinite-field δiso…** in stage 6. The first press opens *QCPMG —
infinite-field δiso (2+ fields)* and adds that field's row (field, δcg ± σ,
FWHM) below two empty rows kept for typed entries; the second adds the other
row; the window stays open and accumulates. On the two EXPNOs the automatic
processing (matched LB 75 Hz on the 78 MHz train; on the higher-field train
the T₂ fit does not converge — `T₂ fit did not converge — no matched filter
offered` — so its sum is unweighted) and the automatic first-minima windows
give stage-6 readings of δCG = −113.1 ± 1.5 ppm and −92.1 ± 1.1 ppm. The
rows the dialog fills read δCG = −113.08 ± 2.7 ppm at 78.354 MHz and
δCG = −92.07 ± 1.1 ppm at 107.811 MHz: the ± of a dataset row is the larger
of the jitter σ and the drift of δCG when the window is doubled (2.7 ppm at
78 MHz — the row's tooltip lists CG(w, 1.5w, 2w, 3w) = −113.0, −115.4,
−115.7, −120.8), so it is never smaller than what the window placement can
do. **Compute δiso** then reads

```
δiso = -68.55 +- 3.82 ppm  (intercept, 1/ν₀²→0)   ·   P_Q = 3.307 +- 0.229 MHz (η-independent)   ·   C_Q = 3.066 +- 0.213 MHz   (η = 0.7 assumed; 2.864-3.307 over η 0-1)
```
<!-- measured 2026-09-23 (wip/QF): both EXPNOs through the stage 1…6 functions as above, sent with the stage-6 button, the dialog's cells as quoted; the headline is fmt_result_lines on those cells -->

together with a warning that the 78 MHz centre of gravity drifts 2.7 ppm
when its window is doubled — 13 % of the 21 ppm separation between the two
fields, i.e. the window uncertainty is a sizeable part of the lever the
slope is built from. The whole-manifold window is no remedy here: the
107.8 MHz sweep integrated in full is dominated by noise (σ 27 ppm), and
for this narrow pattern the first-minima centreband window is the right
one. P_Q carries the η-free uncertainty; C_Q needs the assumed η and its
range over η = 0–1 is a systematic the ± does not contain.

Rows can also be typed in (**＋ Add field**), taken from the open workspace
(**δcg from open spectrum (visible range)**) or read from saved datasets
(**Add from datasets…**). **Split W_q / W_csd** separates the FWHM values of
all the fields entered into a quadrupolar width and a field-independent
width (shift distribution *and* CSA — an upper bound on shift disorder);
**Export report…** and **Export figure…** write the record (the figure with
a `.json` sidecar listing every point's window, mode and source).

**A whole series at once.** **Tools > QCPMG > Batch infinite-field δiso…**
opens *QCPMG — batch infinite-field δiso*: set the grid to five samples and
two fields and drop the ten `.csv` files written earlier with **Save as
dataset…** onto the cells (`<NMR>/MagLab/DATA/LAW{0-4}Ca-3Cl_850_MHz.csv` —
the `850` in these file names is a misnomer, the header's `larmor_MHz =
78.354` is authoritative — and
`<NMR>/NMRFAM/DATA/2026-06_35Cl/LAW{0-4}Ca-3Cl_1p1GHz.csv`). The row names
are prefilled from the file stems (`LAW0Ca-3Cl` …); each cell is measured
with the same functions and the automatic first-minima window (every
dataset here is a magnitude spectrum, which the cells show as *(mc)* and
the report states). **Compute all** reports `5 of 5 samples extrapolated`
and **Export report…** writes, for the first sample:

```
--- LAW0Ca-3Cl --------------------------------------------------
    nu0 (MHz)      dcg (ppm)   +- err    resid   mode        CT-selective (declared)
       78.3541      -112.76     2.25     0.00  magnitude   ?
                 window -206.8 ... -35.4 ppm (minima) · magnitude (mc) · dataset LAW0Ca-3Cl_850_MHz.csv
                 CG(w, 1.5w, 2w, 3w) = -112.7, -114.9, -115.0, -119.8
      107.8113       -95.78     0.51     0.00  magnitude   ?
                 window -168.7 ... -33.7 ppm (minima) · magnitude (mc) · dataset LAW0Ca-3Cl_1p1GHz.csv
                 CG(w, 1.5w, 2w, 3w) = -95.7, -96.0, -95.8, -97.9
    all points measured on magnitude (mc) spectra -- dcg is the |spectrum| centroid, not the absorption one;
    rectified noise pulls it toward the window centre, growing with window width and 1/(S/N); compare with the
    phased spectrum on the same window before quoting P_Q to better than a few %
    delta_iso      = -76.77 +- 2.74 ppm
    P_Q            = 2.973 +- 0.202 MHz (eta-independent)
    C_Q            = 2.756 +- 0.187 MHz   (eta = 0.7 assumed; 2.575-2.973 over eta 0-1)
    slope          = -220953 +- 30021 ppm.MHz^2
    chi2/dof       = exact (2 points, no redundancy)
    ! CG drift 2.2 ppm at 78.354 MHz is 13 % of the 17.0 ppm separation between the fields -- the window cuts the pattern
    width split over 2 fields (78.4, 107.8 MHz; FWHM, Sandland Eq. 2)
    W_q            =     48.2 ppm (low field) / 25.5 ppm (high field)
    W_csd          =     34.7 ppm (field-independent: shift distribution + CSA)
```
<!-- measured 2026-09-23 (wip/QF): the ten CSVs through the batch grid (read_field_spectrum / measure_cg / fit_samples / report_text), spin 1.5, η 0.7; the two dcg rows are the grid's cells, which the grid fits as printed -->

followed by the other four samples and a summary table; δiso runs from
-76.8 ± 2.7 ppm (0 Ca) to +63.7 ± 25.9 ppm (4 Ca) across the series, with
P_Q from 2.97 to 5.95 MHz and C_Q from 2.76 to 5.51 MHz at η = 0.7. The
1–4 Ca cells carry flags the 0 Ca cells do not: the 1 Ca window at 78 MHz
cuts the pattern (δCG drifts 8 ppm when it is doubled), the 2–4 Ca windows
are wider than the 204 ppm sideband spacing and *sensitive* (σ 9–16 ppm),
and the 4 Ca window at 78 MHz *catches one sideband only* — the
whole-manifold window mode is made for those cells, on the absorption sum
echo rather than these magnitude datasets.

Compare the two routes on the same 0 Ca sample: −68.6 ± 3.8 ppm from the
freshly processed EXPNOs against -76.8 ± 2.7 ppm from the saved magnitude
datasets, whose higher-field δcg is −95.8 instead of −92.1 ppm. That
3.7 ppm difference comes from a different EXPNO (2 against 3), line
broadening and window, not from the processing mode — the rectified-noise
mechanism of a magnitude spectrum predicts +0.2 ppm at this S/N. The lever
arm from 1/107.8² to zero is short, so a 4 ppm difference in one centre of
gravity moves the intercept by about 8 ppm. The fit's ± is only the
propagation of the two δCG uncertainties (window jitter and convergence
drift); window placement, processing mode / EXPNO / LB and the assumed η
add systematic errors it does not contain — the −68.6 against −76.8 gap is
about 2σ, amplified 2.1× by the short lever arm — and a two-point line has
no degrees of freedom to reveal them. Process both fields the same way
(both phased, or both magnitude — the dialog and the report say *NOT
COMPARABLE* when they are mixed), place the windows deliberately, and read
the flags and the CG(w, 1.5w, 2w, 3w) sequence in the report before quoting
a number.

## 13. Save

**File > Save fit** (Ctrl+S) stores the fit with its QCPMG provenance
next to the data; **File > Save project…** keeps every open workspace
(the QCPMG spectrum, the WURST-corrected copy, the predicted-field
simulation) in one file. The other four bromide glasses (EXPNOs 31–34)
repeat §2–§10, or their saved `.csv` datasets go through the batch fit of
Tutorial 4 with the model built here.

## Where to read more

**Help > User manuals**: *QCPMG* (the six stages, the infinite-field
extrapolation, Sandland's equations), *Lineshapes — models & physics* (the
`quad_ct` and Czjzek models) and *Processing reference*, §5 (the WURST
correction).
