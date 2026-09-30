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
| χ² profile | **Fit > Errors > χ² profile…** | 15 refits per parameter | nothing about the shape; one parameter at a time | the scan range is wrong for the parameter |
| Monte-Carlo | **Fit > Errors > Monte-Carlo errors…** | 200 refits | the model is right and the residual is white noise at its measured level | the residual is structured, or a bound truncates the trials |

dmfit's *Errors Analysis* is the χ² profile; its *Monte Carlo* is the third
row. Two further views support them: **Fit > Errors > Parameter
correlations…** shows which pairs the data cannot separate, and
**Fit > Errors > χ² map (parameter pair)…** draws χ² over a plane of two
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
    s0_pos:    62.8210630 +/- 0.99261661 (1.58%) (init = 62.81532)
    s0_sigma:  1.47181970 +/- 0.27242975 (18.51%) (init = 1.469871)
    s0_dCS:    14.0125910 +/- 0.86658410 (6.18%) (init = 14.01425)
    s0_amp:    3754949.40 +/- 392220.860 (10.45%) (init = 3751916)
    s1_pos:    28.1114260 +/- 112.626550 (400.64%) (init = 28.59157)
    s1_sigma:  0.58777647 +/- 14.1532210 (2407.93%) (init = 0.6415538)
    s1_dCS:    33.2876870 +/- 29.3828970 (88.27%) (init = 33.25384)
    s1_amp:    1261200.60 +/- 510545.430 (40.48%) (init = 1264975)
    s2_pos:   -0.18019479 +/- 4.33778690 (2407.28%) (init = -0.1825506)
    s2_sigma:  0.80709882 +/- 0.61714862 (76.47%) (init = 0.8088475)
    s2_dCS:    17.7327010 +/- 8.10558130 (45.71%) (init = 17.68899)
    s2_amp:    575527.580 +/- 260626.200 (45.28%) (init = 573659.8)
[[Correlations]] (unreported correlations are < 0.100)
    C(s1_pos, s1_sigma)   = +0.9993
    C(s0_sigma, s1_amp)   = -0.9944
    C(s0_pos, s0_sigma)   = +0.9704
    C(s0_pos, s1_amp)     = -0.9469
    …

normalized RMSD: 0.0458
```
<!-- measured v0.15.1, 2026-09-30: `larmor fit copy_27Al.recipe.json --window 150 -80` on a fresh copy of the shipped recipe, which since 0.15.1 is saved at the deeper minimum the Jacobian-scaled solver reaches (0.13.0 quoted RMSD 0.0468, a shallower minimum) -->

These are the covariance errors: the inverse of the curvature of χ² at the
minimum, propagated through the correlations. They are honest about the AlO₅
site (δiso = 28 ± 113 ppm, σ(C_Q) = 0.6 ± 14 MHz) but they say nothing about
*why* — for that, §3. Note also that the report is that of the first refit of
an already converged recipe: the values differ from the shipped ones in the
third digit (62.815 → 62.821 ppm for the AlO₄ position), an amount well
inside the ±1.0 ppm error — while the degenerate AlO₅ errors move a lot
between two such refits (±213 ppm stored in the recipe, ±113 here): the
curvature of a flat valley is not a stable number.

In the app, load the same copy (**File > Open…**, Ctrl+O) and press **Fit**
(F5). Each parameter shows its `± error` in the *Fit parameters* dock, and the
fit-health strip under the spectrum reads, for this fit:

```
⚠ Fit: check    RMSD 0.0458 · χ²ᵣ 32177719412.09    degenerate ×3: s1.pos↔s1.sigma (+1.00)    residual: structure left
```
<!-- measured v0.15.1, 2026-09-30: fithealth.assess on the fit result above (chi_text, pair count); the Report dock header reads "RMSD 0.0458 · χ²ᵣ 32177719412.09 · residual within noise · ⚠ structured residual · ⚠ 3 unidentifiable pairs (see Correlations)" -->

Read it left to right. The pill is the verdict: `✓ Fit OK`, `⚠ Fit: check`
(something to look at — the values are readable as they stand) or `✗ Fit: not
physical` (an impossible value, and only that). χ²ᵣ is in the spectrum's
intensity² units (residuals are not divided by a noise estimate), so only its
change between fits of the same spectrum means anything. `degenerate ×4`
counts parameter pairs with |r| ≥ 0.95 and names the worst one; clicking the
chip opens the correlation matrix of §3. `residual: structure left` is a runs
test on the residual's sign: the misfit is not random noise — its tooltip
names the usual causes (sidebands not modelled, phasing, baseline). The same
chip reads `residual N× noise` when the residual in the signal region exceeds
the noise at the edges of the window. `degenerate ×3` here: it read ×4 on the
shallower minimum earlier versions fitted into. Other chips appear when they apply:
`at a bound: A lb at its lower bound 0` when a free parameter finished at its
min/max (Tutorial 2, §4 — the errors are then conditional on the pinned
value, and the recipe carries a note saying so), `frozen` for sites the fit
held outside the window, and the red `unphysical` chip for values such as η
outside [0, 1]. Editing any parameter after the fit dims the strip and
re-labels it "Model … edited since fit" until the next Fit.

## 3. Correlations — where the covariance comes from

**Fit > Errors > Parameter correlations…** opens the correlation matrix of
the last fit (red +1, blue −1). Below it:

```
Strongest: s1.pos↔s1.sigma (+1.00)  ·  s0.sigma↔s1.amp (-0.99)  ·  s0.pos↔s0.sigma (+0.97)  ·  s0.pos↔s1.amp (-0.95)  ·  s1.pos↔s2.amp (-0.93)  ·  s0.amp↔s1.amp (-0.93)
```

and, because three pairs exceed |r| = 0.95, the red banner:

```
⚠ Unidentifiable (|r| ≥ 0.95): s1.pos↔s1.sigma (+1.00)  ·  s0.sigma↔s1.amp (-0.99)  ·  s0.pos↔s0.sigma (+0.97). These pairs are degenerate — fix or link one, or add a constraint; don't trust their separate values.
```
<!-- measured v0.15.1, 2026-09-30: correlation matrix of the fit above, formatted as correlation_dialog.py does -->

The AlO₅ position and its Czjzek width move together with r = +0.999: a
larger σ(C_Q) shifts the simulated centre of gravity to lower ppm, which a
higher δiso undoes. The ±113 ppm error on δiso is the length of that valley,
not a measurement of the position. This is the case Tutorial 2 resolves with
a constraint from independent information.

## 4. The χ² profile

**Fit > Errors > χ² profile…** opens *Errors Analysis —
χ² profile*. Choose site `s1 — AlO$_5$` and parameter
`isotropic_chemical_shift_ppm`, then **Scan**. The parameter is fixed at
`points` values (15) spanning `span` standard errors (3) each way around the
fitted value — here ±338 ppm, from the ±113 ppm error of §2 — and *every
other free parameter is refitted at each point*, so the profile absorbs the
correlations of §3 instead of assuming them away. The progress bar counts the
refits and **Stop** keeps what has been scanned. The curve is χ² against the
scanned value, with dashed lines at the 1σ and 2σ levels and a dotted
vertical at the fitted value. The fit minimises an unweighted sum of squares,
so those levels are χ²min·(1 + 1.00/ν) and χ²min·(1 + 3.84/ν) with
ν = 614 − 12 = 602 degrees of freedom — the residual variance χ²min/ν stands
in for the unknown noise, exactly as it does in the covariance errors of §2;
here the 1σ level sits 0.17 % above the minimum. On this parameter the first
scan gives, relative to its minimum:

```
scanned δiso (ppm):  -310   -262   -213   -165   -117    -68    -20     28     76    125    173    221    269    318    366
χ² / χ²min:          1.119  1.119  1.119  1.119  1.119  1.122  1.119  1.000  2.138  1.119  1.119  1.119  1.119  1.119  1.119
```
<!-- measured v0.15.1, 2026-09-30: autofit.error_profile(recipe, ppm, amp, site=1, param="isotropic_chemical_shift_ppm", window_ppm=(150, -80), n_points=15, span=3.0) on the fit of §2 -->

Every point but the fitted one sits on a plateau 12 % above the minimum:
pushed anywhere from −310 to +366 ppm, the AlO₅ line simply leaves the
central band and the other two sites absorb what they can (at +76 ppm the
refit even lost its way, 2.1× the minimum). The readout under the plot says
what is wrong with this scan:

```
s1.isotropic_chemical_shift_ppm = 28.11 [1σ: 27.44 … 28.18]   ·   1σ interval narrower than the scan step — rerun with a smaller `span` to resolve it
```

The ±0.4 ppm interval is not a measurement: the crossings are interpolated
inside the first 48 ppm step. The scan range — set from the ±113 ppm
covariance error — is far too wide to resolve the interval it was meant to
bracket. Set `span` to 0.1 (±11 ppm) and **Scan** again:

```
scanned δiso (ppm):   16.8   18.5   20.1   21.7   23.3   24.9   26.5   28.1   29.7   31.3   32.9   34.5   36.2   37.8   39.4
χ² / χ²min:          1.080  1.056  1.035  1.017  1.005  1.001  1.000  1.000  1.000  1.000  1.001  1.003  1.006  1.010  1.016

s1.isotropic_chemical_shift_ppm = 26.5 [1σ: 24.53 … 33.53]
```
<!-- measured v0.15.1, 2026-09-30: the same call with span=0.1; span=0.3 (±34 ppm) gives a ragged profile whose left end, at −5.7 ppm, is within 0.4 % of the minimum: the site has swapped roles with AlO₆ -->

Now the bottom of the valley is visible: flat to within 0.1 % from 26.5 to
31.3 ppm, so the honest interval for the AlO₅ position is about 25–34 ppm
(the minimum is reported at 26.5 rather than 28.1 because the bottom is
flat), nine ppm wide rather than the ±113 ppm of §2 or the ±0.4 ppm of the
first scan. At `span` 0.3 (±34 ppm) the profile turns ragged and its left end,
at −5.7 ppm, comes back within 0.4 % of the minimum: there the AlO₅ line has
taken the place of AlO₆. A flat or ragged profile is the signature of an
undetermined parameter; a determined one gives a parabola (§7) — the AlO₄
position, scanned the same way, gives `62.82 [1σ: 62.39 … 63.19]` at the
default span, in line with its ±1.0 ppm covariance error.

If a scan ever reports *"the scan reached a χ² N % below the fit's own … the
fit had not converged"*, the fit on screen had stopped short of the minimum
the refits find: press **Fit** (or **Auto fit**) again before quoting any of
its numbers.

## 5. Monte-Carlo errors

**Fit > Errors > Monte-Carlo errors…** opens *Monte-Carlo errors —
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
Monte-Carlo errors from 200/200 synthetic refits · noise σ = 1.247e+05

s0.isotropic_chemical_shift_ppm       62.7104 ± 0.5333   (0.85%)
s0.sigma_Cq_MHz                       1.44914 ± 0.1206   (8.33%)
s0.shift_fwhm_ppm                     13.9607 ± 0.5675   (4.06%)
s0.amplitude                      3.67701e+06 ± 1.811e+05   (4.92%)
s1.isotropic_chemical_shift_ppm       29.4495 ± 3.231   (10.97%)
s1.sigma_Cq_MHz                      0.616652 ± 0.3146   (51.02%)
s1.shift_fwhm_ppm                     33.6551 ± 4.987   (14.82%)
s1.amplitude                      1.30532e+06 ± 2.165e+05   (16.59%)
s2.isotropic_chemical_shift_ppm      -1.20396 ± 2.311   (191.99%)
s2.sigma_Cq_MHz                      0.640149 ± 0.3332   (52.06%)
s2.shift_fwhm_ppm                     17.4384 ± 2.938   (16.85%)
s2.amplitude                           556053 ± 5.904e+04   (10.62%)
```
<!-- measured v0.15.1, 2026-09-30: autofit.monte_carlo_errors(recipe, ppm, amp, window_ppm=(150, -80), n_trials=200, seed=0) on the fit of §2; 92 s across the process pool (parallel="auto" chose it: the first trial timed the rest at well over the pool's start-up) -->

Compare with §2 site by site:

- AlO₄: δiso 62.7 ± 0.5 ppm (Monte-Carlo) against 62.8 ± 1.0 ppm
  (covariance). Same scale; the two estimators agree on a well-determined
  site, the covariance figure a little wider because it carries the +0.97
  correlation with σ(C_Q) (§3).
- AlO₅: ±3.2 ppm against ±113 ppm — thirty times smaller. The trials for this
  site ranged from 23.6 to 42.8 ppm: every trial starts from the best fit
  and is fitted with the same tolerance, so it settles in the flat bottom of
  the valley §4 mapped (25–34 ppm) and rarely walks along it. The
  Monte-Carlo σ of a degenerate parameter is not a confidence interval;
  look at its histogram, which is skewed and far from Gaussian here, before
  quoting it.
- AlO₆: ±2.3 ppm against ±4.3 ppm; the trials ranged from −6.1 to +3.3 ppm.

Had the fit on screen not been at the minimum, the run's status line would
have said so — *"the refit moved … by more than σ: the fit on screen had not
converged"* — because the trials scatter around the run's own refit, not
around the numbers in the table.

The noise level also deserves a look: σ = 1.25 × 10⁵ is the residual over the
window, which the strip in §2 flagged as structured (a runs test on its
sign). The Monte-Carlo run treats that misfit as if it were random noise —
it answers "how would the parameters scatter if the model were right", not
"is the model right".

## 6. The χ² map

**Fit > Errors > χ² map (parameter pair)…** draws χ² over a 15 × 15 grid
(`grid:` box) spanning ±25 % of each of two parameters' fitted values; pick
them as `X:` and `Y:` and press **Compute**, **Export figure…** saves the
image. Two maps of the fit above:

- `s1.isotropic_chemical_shift_ppm` against `s1.sigma_Cq_MHz` (±7 ppm,
  ±0.15 MHz): along σ(C_Q) at the fitted δiso, χ² changes by 2 % over the
  whole range; along δiso at the fitted σ(C_Q) it climbs to 1.65× — yet the
  profile of §4 was flat over the same ±7 ppm. The two agree: the map holds
  σ(C_Q) fixed while the profile lets it follow, and the r = +0.999 valley
  runs diagonally across the map, off the grid at both ends.
- `s0.isotropic_chemical_shift_ppm` against `s0.sigma_Cq_MHz`: a closed
  basin. Along δiso (±16 ppm) χ² climbs to 31 times its minimum, along
  σ(C_Q) (±0.37 MHz) to 4.5 times, and the ellipse is visibly tilted along
  the +0.97 correlation of §3.
<!-- measured v0.15.1, 2026-09-30: chi2map.chi2_surface on the fit of §2 for both pairs, n=15, span 0.25 -->

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
<!-- measured v0.13.0, 2026-09-22: `larmor fit copy_11B.recipe.json --window 30 -30` on a fresh copy of the shipped recipe -->

The fit is excellent (RMSD 0.0037, twelve times lower than the ²⁷Al one) and
the recipe was already converged, so the optimizer barely moves — and lmfit
cannot form a covariance from a Jacobian at the starting point. The
`± error` column stays empty, and **Fit > Errors > Parameter correlations…**
answers `No covariance available — run a fit first (a fit pinned at bounds or
with fixed parameters may not report one).` LARMOR retries the covariance
with a second algorithm before giving up; when that fails too, the two other
estimators still work.

The χ² profile does: **Fit > Errors > χ² profile…**, site
`s2 — BO$_4$`, parameter `isotropic_chemical_shift_ppm`, **Scan**. With no
standard error to set the range, the scan spans ±25 % of the value
(−0.54 … −0.08 ppm), and the profile is a clean parabola that rises to 40
times its minimum at the edges:

```
scanned δiso (ppm):  -0.537  -0.504  -0.471  -0.438  -0.405  -0.373  -0.340  -0.307  -0.274  -0.241  -0.208  -0.175  -0.142  -0.110  -0.077
χ² / χ²min:          40.16   29.64   21.37   13.92    8.13    4.15    1.80    1.00    1.77    4.06    7.86   13.18   20.03   28.34   38.12
```
<!-- measured v0.13.0, 2026-09-22: autofit.error_profile(recipe, ppm, amp, site=2, param="isotropic_chemical_shift_ppm", window_ppm=(30, -30)) on the fit above -->

With ν = 663 − 19 = 644 the 1σ level lies 0.16 % above the minimum, and the
readout is `-0.3068 [1σ: -0.3069 … -0.3068]` with the same "narrower than the
scan step" note as in §4 — this time for the opposite reason: the parabola is
so steep that the interval is a few thousandths of a ppm wide. Rerun with
`span` 0.02 to trace it. The BO₄ position is a well-conditioned parameter,
even though the covariance could not say so. The Monte-Carlo estimate for this
fit was not run for this tutorial; on a well-determined, unbounded parameter
it agrees with the profile.

## 8. Constraints, bounds and populations

- A linked parameter (Tutorial 2) is not fitted; its error is propagated
  from the parameters it depends on, in all three estimators.
- A parameter that finishes at a bound is reported with the flag of §2 and a
  note in the recipe; every error on that fit is conditional on the pinned
  value, and every exported table marks the value with ‡ (fixed values †,
  linked values §) with a footnote naming the bound. Loosen the bound or fix
  the parameter deliberately, then refit.
- **Fit > Report** (F6) integrates each site over the
  window and reports populations `± error`; those errors are first-order
  propagation of the amplitude covariance and inherit its limitations.
- The Report's family sums and named ratios (the Σ rows and N₄ of lines
  tagged in the `family` column) take their error from the Monte-Carlo
  per-trial integrals once **Use as fit errors** is pressed (the dialog
  lists `s0.population_pct`, `family.BO4`, `ratio.N4` with histograms), from
  the amplitude covariance (lmfit's correlated values) right after a plain
  fit, and are flagged *independent* otherwise — the row tooltip, Copy CSV
  and the Methods sentence state which basis was used.
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
