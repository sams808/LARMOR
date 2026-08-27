# Tutorial 2 — Constraining a fit: fix, bound, link

*Time: ~15 minutes. Follows on from Tutorial 1 (you have fitted
`examples/pCABS2-4_27Al.recipe.json` and know the three-site ²⁷Al model is
degenerate). Work on a copy of the recipe so the shipped example stays
pristine.*

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

In Tutorial 1 the fit told you the AlO₅ and AlO₆ sites were poorly
determined — the optimizer can trade their widths, positions, and amplitudes
against the neighboring sites almost freely. Suppose independent information
— an MQMAS spectrum, or a composition argument — gives you an AlO₅/AlO₄
population ratio of about 0.19, and says the two sites, being in the same
glass, share the same chemical-shift-distribution width. Those two facts are
exactly two constraints.

## 2. Add constraints in the recipe file

Copy the recipe, then edit site 1 (the AlO₅ site) in the copy:

```json
"amplitude":      { "value": 740000, "vary": true, "min": 0.0, "expr": "0.19 * s0.amplitude" },
"shift_fwhm_ppm": { "value": 21,     "vary": true, "min": 0.0, "expr": "s0.shift_fwhm_ppm" }
```

(Only the `expr` field matters; `value` becomes the starting point and is then
derived.) Re-run:

```
larmor fit my_constrained.recipe.json --window 150 -80 --plot fit_constrained.png
```

You should see:

- `s1_amp` reported as linked in the fit report, with a stderr that is
  exactly 0.19 × the stderr of `s0_amp` — that is error propagation through
  the link, not a coincidence of rounding.
- The AlO₅ position error drops from ±6.0 to ±2.0 ppm: the constraints
  removed exactly the freedom that made the site undetermined.
- The RMSD is essentially unchanged (≈ 0.047), consistent with the
  constraint. If it had risen sharply, the constraint would be fighting the
  data.
- A warning: `s2.sigma_Cq_MHz` finished at a bound. Constraining two sites
  changed what the optimizer could do with the third, and the weakest site's
  Czjzek width collapsed to its lower bound — which brings us to §4.

## 3. Or add them in the app

```
larmor desktop
```

Open the file, then click the **⚙** button on a site's card (its tooltip
reads "constraints: link expression / min / max"). Each parameter gains a
constraint row with three fields:

- **link** — type the expression, e.g. `0.29 * s0.amplitude`. The value box
  greys out and the parameter label gains a ⚭ mark; the parameter now follows
  its expression live in the plot.
- **min / max** — box bounds for the fit.
- The plain checkbox next to each value still fixes it outright (checked =
  fitted, unchecked = fixed).

## 4. When a constraint fights the data — and when it hides

When a parameter finishes at a bound, the fit flags it and a note is written
into the recipe so the caveat travels with the result:

```
parameters finished at a bound (check constraints/starting model; uncertainties are conditional on them): s2.sigma_Cq_MHz
```

That is what happened in §2: constraining sites 0 and 1 changed what the
optimizer could do with site 2, and its Czjzek width ran to the edge of its
allowed range. LARMOR reports the remaining uncertainties *conditional* on
the pinned value. Parameters at bounds after a constrained fit mean the
constraints, the starting model, and the data disagree somewhere — revisit
one of them.

The opposite failure is quieter. Refit with a deliberately wrong ratio,
`0.5 * s0.amplitude`: on this dataset the RMSD barely moves and nothing
finishes at a bound — the undetermined AlO₆ site simply absorbs the error by
shifting its own position and width. A degenerate fit can hide a wrong
constraint completely, which is why the ratio should come from independent
information (MQMAS, composition), not from trying values until the RMSD looks
good.

## 5. Constraint cookbook

| Goal | expr on which parameter | expression |
|---|---|---|
| Population ratio from MQMAS/chemistry | site j `amplitude` | `0.19 * s0.amplitude` |
| Shared Gaussian width across sites of one phase | site j `shift_fwhm_ppm` | `s0.shift_fwhm_ppm` |
| Fixed shift difference (e.g. crystallographic pair) | site j `isotropic_chemical_shift_ppm` | `s0.isotropic_chemical_shift_ppm - 12.5` |
| Two sites, equal populations | site j `amplitude` | `s0.amplitude` |
| Keep Cq width physical | site j `sigma_Cq_MHz` | set `min`/`max` instead of expr |

Anything lmfit accepts is valid — `sin`, `exp`, ratios of other parameters —
but if you find yourself writing something elaborate, consider whether the
model itself should change instead.
