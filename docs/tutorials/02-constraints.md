# Tutorial 2 — Constraining a fit: fix, bound, link

*Time: ~15 minutes. Follows on from Tutorial 1 (you have
`CaAlGlass.recipe.json` and know the three-site ²⁷Al model is degenerate).*

LARMOR supports the three kinds of constraint you know from ssNake, in a more
general algebraic form:

| Constraint | ssNake | LARMOR recipe |
|---|---|---|
| Fix a value | tick the checkbox | `"vary": false` |
| Bound a range | — | `"min": 0.05, "max": 12` |
| Link parameters | `('Amplitude', 0, 0.5, 0)` tuple | `"expr": "0.5 * s0.amplitude"` |

An `expr` is any algebraic expression. Site parameters are addressed as
`s<index>.<parameter>` — site indices are shown in the app (`s0:`, `s1:` …)
and correspond to the order in the recipe file. A linked parameter is no
longer varied independently: its value *and its standard error* are derived
from the parameters it depends on (full error propagation).

## 1. Why constrain?

In Tutorial 1 the fit told you the AlO₅ site was undetermined — the optimizer
can trade its width, position, and amplitude against the neighboring sites
almost freely. Suppose you know from an MQMAS spectrum (or a composition
argument) that the AlO₅/AlO₄ population ratio is ~0.29 and that both sites,
being in the same glass, should share the same chemical-shift-distribution
width. Those two facts are exactly two constraints.

## 2. Add constraints in the recipe file

Open `CaAlGlass.recipe.json` and edit site 1 (the AlO₅ site):

```json
"amplitude":      { "value": 465, "vary": true, "min": 0.0, "expr": "0.29 * s0.amplitude" },
"shift_fwhm_ppm": { "value": 27,  "vary": true, "min": 0.1, "expr": "s0.shift_fwhm_ppm" }
```

(Only the `expr` field matters; `value` becomes the starting point and is then
derived.) Re-run:

```
larmor fit CaAlGlass.recipe.json --window 150 -80 --plot fit_constrained.png
```

You should see:

- `s1_amp` reported as `== '0.29 * s0_amp'` in the fit report, with a stderr
  that is exactly 0.29 × the stderr of `s0_amp` — that's error propagation,
  not a coincidence.
- The RMSD rises slightly (constraints remove freedom; ~0.004 vs 0.0025 free).
  A small rise is the price of a physically meaningful model. A large rise
  means your constraint is wrong.

## 3. Or add them in the app

```
larmor app
```

Load the file, then click **constraints ▸** on any site. Each parameter gains
three fields:

- **link** — type the expression, e.g. `0.29 * s0.amplitude`. The value box
  greys out and shows a ⚭ symbol; the parameter now follows its expression
  live in the plot.
- **min / max** — box bounds for the fit.
- The plain checkbox next to each value still fixes it outright.

## 4. When a constraint fights the data

Try linking the amplitude with a deliberately wrong ratio, e.g.
`0.5 * s0.amplitude`, and fit again. The report now warns:

```
⚠ at bounds (constraint or start value fighting the data?): s4.amplitude
```

What happened: forced to hold a 2:1 ratio the data doesn't support, the
optimizer pushed other parameters to the edges of their allowed ranges — the
third site's amplitude collapsed to zero. LARMOR detects parameters that
finish at a bound, warns you, and reports the remaining uncertainties
*conditional* on those pinned values (the note is also written into the
recipe, so the caveat travels with the result).

The rule of thumb: **parameters at bounds after a constrained fit mean the
constraint and the data disagree** — revisit one of them.

## 5. Constraint cookbook

| Goal | expr on which parameter | expression |
|---|---|---|
| Population ratio from MQMAS/chemistry | site j `amplitude` | `0.29 * s0.amplitude` |
| Shared Gaussian width across sites of one phase | site j `shift_fwhm_ppm` | `s0.shift_fwhm_ppm` |
| Fixed shift difference (e.g. crystallographic pair) | site j `isotropic_chemical_shift_ppm` | `s0.isotropic_chemical_shift_ppm - 12.5` |
| Two sites, equal populations | site j `amplitude` | `s0.amplitude` |
| Keep Cq width physical | site j `sigma_Cq_MHz` | set `min`/`max` instead of expr |

Anything lmfit accepts is valid — `sin`, `exp`, ratios of other parameters —
but if you find yourself writing something elaborate, consider whether the
model itself should change instead.
