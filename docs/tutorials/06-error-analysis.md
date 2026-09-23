# Tutorial 6 — Error analysis: covariance, Monte-Carlo and the χ² profile

*Time: ~20 minutes plus the Monte-Carlo run. Data: the two shipped recipes,
`examples/pCABS2-4_27Al.recipe.json` (the degenerate three-site ²⁷Al fit of
Tutorial 1) and `examples/pCABS2-4_11B.recipe.json` (a well-converged ¹¹B
fit). Work on copies placed in the repository root, so the shipped files stay
pristine and the recipes' relative data paths resolve. Run the commands from
the repository root.*

Every fit in LARMOR comes with an uncertainty per parameter, which dmfit does
not report. This tutorial shows the three ways LARMOR estimates that
uncertainty, what each assumes, and — on a fit that is genuinely degenerate —
how differently they can answer.

## 1. Three estimators

| Estimator | Where | Cost | Assumes | Fails when |
|---|---|---|---|---|
| Covariance | printed after every **Fit** (F5), `± error` column | none | χ² is locally quadratic and well conditioned | parameters are strongly correlated or pinned at a bound; sometimes unavailable |
| χ² profile | **Decomposition > Errors Analysis (χ² profile)…** | 15 refits per parameter | nothing about the shape; one parameter at a time | the scan range is wrong for the parameter |
| Monte-Carlo | **Decomposition > Monte-Carlo errors…** | 200 refits | the model is right and the residual is white noise at its measured level | the residual is structured, or a bound truncates the trials |

dmfit's *Errors Analysis* is the χ² profile; its *Monte Carlo* is the third
row. Two further views support them: **Decomposition > Parameter
correlations…** shows which pairs the data cannot separate, and
**Decomposition > χ² map (parameter pair)…** draws χ² over a plane of two
parameters.

## 2. Covariance — what every fit prints

Copy the ²⁷Al recipe to `copy_27Al.recipe.json` in the repository root and
fit it:

```
larmor fit copy_27Al.recipe.json --window 150 -80
```

The report ends with the fitted parameters, their standard errors and the
strongest correlations (abridged):

```
[[Variables]]
    s0_pos:    63.8426905 +/- 0.58285702 (0.91%) (init = 63.80425)
    s0_sigma:  1.74518403 +/- 0.15053457 (8.63%) (init = 1.744092)
    s0_dCS:    12.9583804 +/- 0.82633519 (6.38%) (init = 13.03105)
    s0_amp:    3898915.49 +/- 50182.4137 (1.29%) (init = 3911614)
    s1_pos:    30.8997341 +/- 7.45155915 (24.12%) (init = 31.34656)
    s1_sigma:  1.31839293 +/- 0.73063415 (55.42%) (init = 1.359326)
    s1_dCS:    23.0806229 +/- 8.80905266 (38.17%) (init = 21.45238)
    s1_amp:    736477.871 +/- 257925.217 (35.02%) (init = 736477.9)
    s2_pos:   -1.99292032 +/- 23.8299558 (1195.73%) (init = -2.185331)
    s2_sigma:  0.58338151 +/- 3.24852434 (556.84%) (init = 0.514625)
    s2_dCS:    17.6512577 +/- 14.1136381 (79.96%) (init = 18.30487)
    s2_amp:    357377.292 +/- 163660.390 (45.79%) (init = 357377.3)
[[Correlations]] (unreported correlations are < 0.100)
    C(s2_pos, s2_sigma)   = +0.9968
    C(s0_sigma, s1_amp)   = -0.9886
    C(s1_pos, s1_sigma)   = +0.9708
    C(s0_pos, s0_sigma)   = +0.9577
    …

normalized RMSD: 0.0468
```
<!-- measured v0.12.1, 2026-09-22: `larmor fit copy_27Al.recipe.json --window 150 -80` on a fresh copy of the shipped recipe -->

These are the covariance errors: the inverse of the curvature of χ² at the
minimum, propagated through the correlations. They are honest about the AlO₆
site (δiso = −2 ± 24 ppm) but they say nothing about *why* — for that, §3.
Note also that the report is that of the first refit of an already converged
recipe: the values differ from the shipped ones in the third digit (63.80 →
63.84 ppm for the AlO₄ position), an amount well inside the ±0.58 ppm error.

In the app, load the same copy (**File > Open…**, Ctrl+O) and press **Fit**
(F5). Each parameter shows its `± error` in the *Fit parameters* dock, and the
results strip above the report reads, for this fit:

```
RMSD 0.0468 · χ²ᵣ 33549623772.66   ·   ⚠ residual 8.7× noise (structure left)   ·   ⚠ structured residual   ·   ⚠ 4 unidentifiable pairs (see Correlations)
```
<!-- measured v0.12.1, 2026-09-22: the same fit result passed through the functions the strip calls (MainWindow._residual_noise_ratio, diagnostics.residual_structure, sanity.check_recipe, identifiability.unidentifiable_pairs) -->

Read it as a checklist. χ²ᵣ is in the spectrum's intensity² units (residuals
are not divided by a noise estimate), so only its change between fits of the
same spectrum means anything. `residual 8.7× noise` compares the residual in
the signal region with the noise at the edges of the window: the model leaves
structure behind. `structured residual` is a runs test on the residual's
sign. `4 unidentifiable pairs` counts parameter pairs with |r| ≥ 0.95. Two
more flags can appear here: `⚠ at bounds: …` when a parameter finished at
its min/max (Tutorial 2, §4 — the errors are then conditional on the pinned
value, and the recipe carries a note saying so), and `⚠ N physical warnings`
for unphysical values such as η outside [0, 1].

## 3. Correlations — where the covariance comes from

**Decomposition > Parameter correlations…** opens the correlation matrix of
the last fit (red +1, blue −1). Below it:

```
Strongest: s2.pos↔s2.sigma (+1.00)  ·  s0.sigma↔s1.amp (-0.99)  ·  s1.pos↔s1.sigma (+0.97)  ·  s0.pos↔s0.sigma (+0.96)  ·  s1.sigma↔s2.amp (-0.94)  ·  s0.pos↔s1.amp (-0.93)
```

and, because four pairs exceed |r| = 0.95, the red banner:

```
⚠ Unidentifiable (|r| ≥ 0.95): s2.pos↔s2.sigma (+1.00)  ·  s0.sigma↔s1.amp (-0.99)  ·  s1.pos↔s1.sigma (+0.97)  ·  s0.pos↔s0.sigma (+0.96). These pairs are degenerate — fix or link one, or add a constraint; don't trust their separate values.
```
<!-- measured v0.12.1, 2026-09-22: correlation matrix of the fit above, formatted as correlation_dialog.py does -->

The AlO₆ position and its Czjzek width move together with r = +0.997: a
larger σ(C_Q) shifts the simulated centre of gravity to lower ppm, which a
higher δiso undoes. The ±24 ppm error on δiso is the length of that valley,
not a measurement of the position. This is the case Tutorial 2 resolves with
a constraint from independent information.

## 4. The χ² profile

**Decomposition > Errors Analysis (χ² profile)…** opens *Errors Analysis —
χ² profile*. Choose site `s2 — AlO$_6$` and parameter
`isotropic_chemical_shift_ppm`, then **Scan**. The parameter is fixed at 15
values spanning ±3 standard errors around the fitted value (here
−73 … +70 ppm) and *every other free parameter is refitted at each point*,
so the profile absorbs the correlations of §3 instead of assuming them away.
The curve is χ² against the scanned value, with dashed lines at χ²min + 1.00
(1σ) and χ²min + 3.84 (2σ). On this parameter the scan gives, relative to its
minimum:

```
scanned δiso (ppm):  -73.5  -63.3  -53.1  -42.8  -32.6  -22.4  -12.2   -2.0    8.2   18.4   28.6   38.9   49.1   59.3   69.5
χ² / χ²min:          1.119  1.119  1.124  1.123  1.126  1.120  1.081  1.000  1.078  1.124  1.077  1.094  1.129  1.143  1.128
```
<!-- measured v0.12.1, 2026-09-22: autofit.error_profile(recipe, ppm, amp, site=2, param="isotropic_chemical_shift_ppm", window_ppm=(150, -80)) on the fit of §2 -->

χ² rises by only 8–14 % over a 140 ppm excursion of the AlO₆ position, and
not monotonically — the other eleven parameters keep finding almost equally
good fits wherever this one is put. That flat, ragged profile is the
signature of an undetermined parameter; a determined one gives a parabola
(§7).

The readout under the plot for this scan is

```
s2.isotropic_chemical_shift_ppm = -1.993 [1σ: -1.993 … -1.993]
```

a zero-width interval, which is not a measurement either. The +1.00 and
+3.84 levels are in the same unweighted units as χ² itself (here χ²min is
about 1.9 × 10¹³), so on real-intensity data they sit indistinguishably
close to the minimum and the interpolated crossings collapse onto the best
value. With residuals of unknown scale the standard remedy is to estimate
the noise variance from the fit, χ²min/ν with ν = 614 − 12 = 602 degrees of
freedom, which puts the 1σ level 0.17 % above χ²min. Read the curve with
that in mind rather than the printed interval: the nearest scanned points
(±10 ppm) are already 8 % above the minimum, so the scan range — set from the
±24 ppm covariance error — is far too wide to resolve the interval it was
meant to bracket, and the honest statement remains the one of §3.

## 5. Monte-Carlo errors

**Decomposition > Monte-Carlo errors…** opens *Monte-Carlo errors —
synthetic-noise refits*: `trials` 200 and `seed` 0 by default. **Run** first
refits the recipe to fix the best fit, measures the noise as the standard
deviation of the residual over the window, then refits 200 synthetic spectra
(best-fit model + Gaussian noise at that level), each starting from the best
fit. The table lists every free parameter with its trial mean and standard
deviation; the `histogram` selector draws the distribution of any one of
them. **Use as fit errors** writes the standard deviations into the table's
`± error` column; **Copy report** copies the text below, which the run on the
²⁷Al fit produced:

```
Monte-Carlo errors from 200/200 synthetic refits · noise σ = 1.255e+05

s0.isotropic_chemical_shift_ppm       63.3163 ± 0.3812   (0.60%)
s0.sigma_Cq_MHz                       1.60846 ± 0.08943   (5.56%)
s0.shift_fwhm_ppm                     13.5566 ± 0.5245   (3.87%)
s0.amplitude                      3.84248e+06 ± 9.069e+04   (2.36%)
s1.isotropic_chemical_shift_ppm       29.4921 ± 3.944   (13.37%)
s1.sigma_Cq_MHz                      0.907177 ± 0.3922   (43.23%)
s1.shift_fwhm_ppm                     27.8836 ± 4.47   (16.03%)
s1.amplitude                           996768 ± 1.575e+05   (15.80%)
s2.isotropic_chemical_shift_ppm       -1.5664 ± 2.341   (149.44%)
s2.sigma_Cq_MHz                      0.541492 ± 0.3624   (66.93%)
s2.shift_fwhm_ppm                     17.1736 ± 3.031   (17.65%)
s2.amplitude                           509763 ± 5.875e+04   (11.53%)
```
<!-- measured v0.12.1, 2026-09-22: autofit.monte_carlo_errors(recipe, ppm, amp, window_ppm=(150, -80), n_trials=200, seed=0) on the fit of §2; the run took several minutes on the development machine -->

Compare with §2 site by site:

- AlO₄: δiso 63.3 ± 0.4 ppm (Monte-Carlo) against 63.8 ± 0.6 ppm
  (covariance). Same scale; the two estimators agree on a well-determined
  site. The mean itself sits 0.5 ppm below the fitted value because the
  initial refit inside the Monte-Carlo run settled at a slightly different
  point of the shallow minimum — again within the quoted error.
- AlO₅: ±3.9 ppm against ±7.5 ppm.
- AlO₆: ±2.3 ppm against ±24 ppm — ten times smaller. The trials for this
  site ranged from −5.6 to +4.0 ppm, and its σ(C_Q) trials ran down to the
  0.05 MHz lower bound in some refits: the bound truncates the distribution,
  and every trial starts from the best fit and is fitted with the same
  tolerance, so the optimizer rarely walks far along the flat valley of §4.
  The Monte-Carlo σ of a bounded, degenerate parameter is not a confidence
  interval; look at its histogram, which is skewed and far from Gaussian
  here, before quoting it.

The noise level also deserves a look: σ = 1.26 × 10⁵ is the residual over the
window, which the strip in §2 flagged as 8.7 times the edge noise. The
Monte-Carlo run treats that misfit as if it were random noise — it answers
"how would the parameters scatter if the model were right", not "is the model
right".

## 6. The χ² map

**Decomposition > χ² map (parameter pair)…** draws χ² over a 15 × 15 grid
(`grid:` box) spanning ±25 % of each of two parameters' fitted values; pick
them as `X:` and `Y:` and press **Compute**, **Export figure…** saves the
image. Two maps of the fit above:

- `s2.isotropic_chemical_shift_ppm` against `s2.sigma_Cq_MHz`: χ² varies by
  less than 1 % over the whole grid. Because δiso is −2 ppm, ±25 % is only
  ±0.5 ppm — a tiny patch at the bottom of a valley that is tens of ppm long,
  so the map is flat and shows nothing beyond that flatness.
- `s0.isotropic_chemical_shift_ppm` against `s0.sigma_Cq_MHz`: a closed
  basin. Along δiso (±16 ppm) χ² climbs to 30 times its minimum, along
  σ(C_Q) (±0.44 MHz) to 6 times, and the ellipse is visibly tilted along the
  +0.96 correlation of §3.
<!-- measured v0.12.1, 2026-09-22: chi2map.chi2_surface on the fit of §2 for both pairs, n=15, span 0.25 -->

The map is most useful for a pair the correlation table flags: a tilted,
closed ellipse means the pair is correlated but determined; an open valley
running off the grid means it is not.

## 7. A control case: the ¹¹B fit

Copy `examples/pCABS2-4_11B.recipe.json` to `copy_11B.recipe.json` and run

```
larmor fit copy_11B.recipe.json --window 30 -30
```

```
[[Fit Statistics]]
    # fitting method   = least_squares
    # function evals   = 40
    # data points      = 663
    # variables        = 19
    …
##  Warning: uncertainties could not be estimated:
    s0_pos:   at initial value
    s0_deta:  at initial value
    …
normalized RMSD: 0.0037
```
<!-- measured v0.12.1, 2026-09-22: `larmor fit copy_11B.recipe.json --window 30 -30` on a fresh copy of the shipped recipe -->

The fit is excellent (RMSD 0.0037, twelve times lower than the ²⁷Al one) and
the recipe was already converged, so the optimizer barely moves — and lmfit
cannot form a covariance from a Jacobian at the starting point. The
`± error` column stays empty, and **Decomposition > Parameter correlations…**
answers `No covariance available — run a fit first (a fit pinned at bounds or
with fixed parameters may not report one).` LARMOR retries the covariance
with a second algorithm before giving up; when that fails too, the two other
estimators still work.

The χ² profile does: **Decomposition > Errors Analysis (χ² profile)…**, site
`s2 — BO$_4$`, parameter `isotropic_chemical_shift_ppm`, **Scan**. With no
standard error to set the range, the scan spans ±25 % of the value
(−0.54 … −0.08 ppm), and the profile is a clean parabola that rises to 40
times its minimum at the edges:

```
scanned δiso (ppm):  -0.537  -0.504  -0.471  -0.438  -0.405  -0.373  -0.340  -0.307  -0.274  -0.241  -0.208  -0.175  -0.142  -0.110  -0.077
χ² / χ²min:          40.16   29.64   21.37   13.92    8.13    4.15    1.80    1.00    1.77    4.06    7.86   13.18   20.03   28.34   38.12
```
<!-- measured v0.12.1, 2026-09-22: autofit.error_profile(recipe, ppm, amp, site=2, param="isotropic_chemical_shift_ppm", window_ppm=(30, -30)) on the fit above -->

Using the χ²min/ν rule of §4 (ν = 663 − 19 = 644), the 1σ level lies 0.16 %
above the minimum: from the parabola's curvature the BO₄ position is
determined to a few thousandths of a ppm — a well-conditioned parameter,
even though the covariance could not say so. The Monte-Carlo estimate for this
fit was not run for this tutorial; on a well-determined, unbounded parameter
it agrees with the profile.

## 8. Constraints, bounds and populations

- A linked parameter (Tutorial 2) is not fitted; its error is propagated
  from the parameters it depends on, in all three estimators.
- A parameter that finishes at a bound is reported with the flag of §2 and a
  note in the recipe; every error on that fit is conditional on the pinned
  value. Loosen the bound or fix the parameter deliberately, then refit.
- **Decomposition > Report (quantify)** (F6) integrates each site over the
  window and reports populations `± error`; those errors are first-order
  propagation of the amplitude covariance and inherit its limitations.
- The batch-fit dialog of Tutorial 4 has an *Error calculation:* selector
  with the same three estimators — `Covariance (from the fit)`,
  `Monte-Carlo (synthetic-noise refits)`, `χ² profile (error analysis)` —
  applied to every spectrum of a series at once (**Compute errors**, then
  **Export CSV…**).

## 9. From Python

No command-line subcommand exists for the profile, the Monte-Carlo run or the
map; the functions the dialogs call are public:

```python
from larmor import autofit, chi2map
from larmor.loader import load_any
from larmor.recipe import Recipe

ppm, amp, rec, meta, warnings = load_any("copy_27Al.recipe.json")
r = Recipe.from_dict(rec)
prof = autofit.error_profile(r, ppm, amp, site=2,
                             param="isotropic_chemical_shift_ppm",
                             window_ppm=(150, -80), parallel=True)
print(prof.summary)                       # the readout of §4
print(prof.values, prof.chi2 / prof.chi2_min)
mc = autofit.monte_carlo_errors(r, ppm, amp, window_ppm=(150, -80),
                                n_trials=200, seed=0, parallel=True)
print(mc.report())                        # the text of §5
A, B, Z, (a0, b0) = chi2map.chi2_surface(rec, ppm, amp, (150, -80),
                                         (2, "isotropic_chemical_shift_ppm"),
                                         (2, "sigma_Cq_MHz"), n=15)
```

`parallel=True` spreads the refits over a process pool, as the dialogs do.
The same `seed` reproduces the same Monte-Carlo trial set regardless of how
the work is scheduled.

## Where to read more

**Help > User manuals**: *1D spectra — processing & fitting*, §3 (the three
estimators and the results strip) and *Fitting glasses for publication
(Edén 2023)* for what to report.
