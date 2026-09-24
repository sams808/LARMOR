# DFT tensors — import & shift calibration

> **Tools ▸ Import DFT tensors (.magres)** turns the shielding and EFG tensors
> of a GIPAW calculation (Quantum ESPRESSO, CASTEP) into starting sites for a
> fit: one site per crystallographic position, the model chosen for the
> nucleus's spin, and — the part that decides whether the result means
> anything — the computed shieldings converted to chemical shifts through a
> calibration line **fitted over reference compounds**, with its uncertainty
> carried into every predicted shift and into the Methods text.

---

## 1 · What a .magres file carries

A `.magres` file (the magres 1.0 text format written by CASTEP and by QE-GIPAW)
holds, per atom, the **magnetic shielding tensor** `ms` (3 × 3, ppm) and, when
requested, the **electric field gradient** `efg` (3 × 3, atomic units). LARMOR
reads both from the `[magres]` block and derives, in the Haeberlen convention,
the isotropic shielding $\sigma_{\rm iso}$, the anisotropy $\zeta_\sigma$ and
the asymmetry $\eta_{\rm CS}$; from the EFG it derives $C_Q$ and $\eta_Q$ with
the isotope's quadrupole moment (mrsimulator's table) — but only for a nucleus
with $I > \frac{1}{2}$. A spin-½ nucleus (¹⁹F, ³¹P, ²⁹Si, ¹³C, ¹H) has no
quadrupole moment: an `efg` record for it describes the crystal, not the
spectrum, and is ignored (the table says so).

The `[calculation]` block records how the numbers were obtained — code and
version, exchange–correlation functional, wave-function and density cutoffs,
one pseudopotential per element, the k-point grid. LARMOR shows a one-line
summary ("QE-GIPAW 7.5 · PBE · 60/480 Ry") and uses the block to decide
whether two files sit on the **same shielding scale** (§2).

**Symmetry-equivalent atoms.** A unit cell lists every atom, so CaF₂ has eight
fluorine records that are the same crystallographic site. Their tensors agree
to numerical noise (their *orientations* differ, their invariants do not).
The import collapses them into one site per position, keeps the first atom's
label (F5), records the members and the **multiplicity** (×8), and averages
the principal values — never the raw 3 × 3 components, which would shrink the
anisotropy of atoms related by a mirror or a rotation. The tolerance
(default 0.05 ppm on $\sigma_{\rm iso}$ and $\zeta$, 0.02 on $\eta$, 0.01 MHz
on $C_Q$) is far above the spread in a well-converged calculation (CaF₂'s Ca
atoms differ by 0.009 ppm, NaF's F by 0.001 ppm); a low-symmetry cell whose
"equivalent" atoms genuinely differ can lower it or switch grouping off.

---

## 2 · Shielding is not shift

A calculated shielding $\sigma$ is an absolute number on a scale fixed by the
calculation. A measured shift $\delta$ is relative to a reference compound.
The two are related linearly,

$$\delta = a\,\sigma + b,$$

but **neither $a$ nor $b$ is universal**. The textbook relation
$\delta = \sigma_{\rm ref} - \sigma$ assumes $a = -1$; in practice the slope
fitted over a set of reference compounds is $-0.7$ to $-1.1$ for ¹⁹F, ²⁹Si and
³¹P at the PBE level (for the CaF₂/NaF pair computed with QE-GIPAW 7.5, PBE,
60/480 Ry, $a = -0.687$), because the functional compresses the shielding
range and the pseudopotential and cutoffs shift its origin. Using a single
$\sigma_{\rm ref}$ therefore predicts the *spread* of the shifts wrongly by
that factor — the twelve kogarkoite lines land at the wrong separations even
when one of them is placed right.

**Why two references are not enough.** A straight line has two parameters. Two
points fix it exactly: the residuals are zero *by construction*, and the fit
has **zero residual degrees of freedom**. The ± on $a$ and $b$ that LARMOR still
reports for a two-point line is propagated from the errors of the two measured
shifts only; it says nothing about whether a third compound would fall on the
line. LARMOR flags this in amber wherever the line is used and writes it into
the Methods sentence. Three or more references give residuals with meaning and
a root-mean-square deviation (RMSE) that is the honest uncertainty of a
predicted shift; the (a, b) covariance is then inflated by the reduced χ² when
the residuals exceed the reference errors.

**Which references may share one line.** Only files computed with the same
code and version, functional, both cutoffs and the same pseudopotential **for
the observed element**. The k-point grid and the pseudopotentials of the other
elements are cell-specific and legitimately differ (CaF₂ carries a Ca
pseudopotential, NaF a Na one); they are recorded, not compared. LARMOR
refuses to fit a line over mismatched references and refuses to apply a line
to a file it does not match, naming the differing keys; an explicit
*Apply anyway* is offered and **recorded in the provenance**, so a researcher
who mixed cutoffs after a convergence test has a path, and the paper knows.
Whether the geometry was relaxed or taken as published is not in the header:
it is a methodological choice to state in the paper, not something LARMOR can
check.

---

## 3 · Workflow — building the calibration line

1. **Measure the references** under the *same referencing* as the sample
   (the same ¹H adamantane reference at 1.82 ppm through the indirect Ξ, the
   same SR; see **Process ▸ Reference ▸ Referencing audit**). Literature shifts are a
   fallback, not a substitute: the 2026-03 CaF₂ and NaF standards both read
   +1.7 ppm relative to the literature values, an offset that the measured
   line absorbs into $b$ and a literature line does not.
2. **Fit each reference spectrum** in LARMOR with one `gl_norm` line (the
   Gaussian fraction free). Its position is $\delta_{\rm exp}$, its standard
   error is the weight of that reference. Save the fit.
3. Open **Tools ▸ Import DFT tensors ▸ Calibration**. *Add reference
   (.magres)…* adds one row per equivalent group of the chosen isotope (the
   name comes from `calc_prefix`); *δ from a saved fit…* fills the selected
   row from the recipe (the largest-amplitude site is proposed, a combo picks
   another), *δ from the current fit* takes the open fit, or type $\delta$ and
   ± into the cells. A row without ± makes the whole fit unweighted, and the
   result says so.
4. **Fit line.** Read $a \pm \sigma_a$, $b \pm \sigma_b$, cov(a, b), n, the
   degrees of freedom, the RMSE and the per-point residuals over the plot of
   $\delta_{\rm exp}$ against $\sigma_{\rm iso}$ (the band is the ±1σ
   prediction when dof ≥ 1). A reference that falls off the line by more
   than its error is a wrong assignment, a wrong referencing or a different
   structure — not a reason to drop it silently.
5. **Save calibration…** writes `<isotope>.shiftcal.json` (the line, its
   covariance, every reference with its source and the `[calculation]`
   header). **Use this line** switches to the Sites tab with the line
   active; *Load…* on the Sites tab reuses a saved one.

The command line does the same: `larmor shiftcal --ref NaF NaF.nmr.magres
-225.0 0.2 --ref CaF2 CaF2.nmr.magres -108.82 0.2 --ref cryolite
cryolite.magres -189.85 --isotope 19F -o f19.shiftcal.json`.

---

## 4 · Import — seeding the sites

**Isotope and spin.** The isotope combo lists what the file contains (the
element's default NMR isotope) with its spin, pre-selected to the open
spectrum's nucleus; when the file lacks that nucleus the dialog says so and
lets the choice stand.

**Model.** Only models the tensor can seed for that spin are offered, the
recommended one first. For $I = \frac{1}{2}$: `gl_norm` (recommended under
MAS — its amplitude is an **area**, so multiplicity locks are exact
populations), `csa_mas` (recommended for a static spectrum; $\zeta$ and
$\eta_{\rm CS}$ seeded from the tensor), `csa_czjzek`, then `gauss_lor` and
`voigt` (position only). For $I > \frac{1}{2}$: `quad_ct` first, then
`quad_csa` (both tensors), `quad_first`, `ext_czjzek`, `amorphous`, the CSA
models, and the Czjzek family (position only — a distribution is not one
tensor). A quadrupolar model is never offered for a spin-½ nucleus; the
command line exits with a one-line reason.

**Sign convention.** LARMOR's `csa_mas` stores the **shielding** anisotropy
$\zeta$ exactly as mrsimulator takes it, with $\delta_{\rm aniso} = -\zeta$
(the same sign **Tools ▸ Herzfeld–Berger** reports). With $\sigma_{\rm ref}$
($a = -1$) the seeded $\zeta$ equals the computed $\zeta_\sigma$; with a
calibration slope $a$ the shift tensor is $a\,\sigma + b$, its anisotropy
$a\,\zeta_\sigma$, and the seeded value is $-a\,\zeta_\sigma$. The dialog
prints this line under the table so the convention is never implicit.

**Multiplicity locks.** *Lock relative populations to multiplicity* links the
amplitude of every site to the first one as $m_k / m_0$ ("0.5A" in the lines
table), and *One shared linewidth* links the width parameter. For `gl_norm`
the lock is on by default and exact; for height-amplitude models (`csa_mas`,
`gauss_lor`, `voigt`) a height ratio equals a population ratio only for
equal-shape lines, so the lock defaults off and its label says *heights —
approximate* — read populations from the Report. Quadrupolar models share no
width: their breadth parameter is $C_Q$.

**Add.** The button is disabled, with the reason shown, until a conversion is
set, a spectrum is open and the line's `[calculation]` fingerprint matches the
file. Adding is one undo step.

---

## 5 · What is recorded

- **A recipe note**: file, header summary, model, the conversion used and one
  entry per site ("s3 F5: ×8 (F5, F6, …), δ_pred −108.8 ± 0.4 ppm, amplitude
  locked to s3").
- **`provenance['dft_import']`**: file path and SHA-256, isotope, model, the
  `[calculation]` record, every grouped site with members, $\sigma_{\rm iso}$,
  predicted shift and its uncertainty, the full calibration record (or
  $\sigma_{\rm ref}$), and which fingerprint keys were overridden. The per-site
  uncertainty lives here, never in a parameter's stderr — stderr belongs to the
  fit.
- **The Methods sentence** (Report ▸ Copy methods, the publication bundle)
  gains: "Starting values came from computed ¹⁹F shielding tensors (QE-GIPAW
  7.5 · PBE · 60/480 Ry); shieldings were converted to chemical shifts through
  δ = a·σ + b fitted over 4 reference compounds (NaF, CaF₂, cryolite,
  sulphohalite; a = −0.698 ± 0.012, b = 51.2 ± 3.9 ppm, RMSE 0.9 ppm)." — or,
  for two references, "the line is exact by construction, no residual degrees
  of freedom".

**Command line.** `larmor magres kog.nmr.magres --isotope 19F --calibration
f19.shiftcal.json --lock-multiplicity --larmor 564.27 -o kog.recipe.json`
writes the grouped, locked recipe with the same provenance; `--reference 560`
keeps the single-σ_ref form, `--no-group` writes one site per atom,
`--calc-override` applies a mismatched line and records it.

---

## References

- U. Haeberlen, *High Resolution NMR in Solids: Selective Averaging*, Academic
  Press (1976). *(the ζ, η convention)*
- J. Herzfeld, A. E. Berger, "Sideband intensities in NMR spectra of samples
  spinning at the magic angle", *J. Chem. Phys.* **73**, 6021 (1980).
- C. J. Pickard, F. Mauri, "All-electron magnetic response with
  pseudopotentials: NMR chemical shifts", *Phys. Rev. B* **63**, 245101 (2001).
  *(GIPAW)*
- C. Bonhomme, C. Gervais, F. Babonneau, C. Coelho, F. Pourpoint, T. Azaïs,
  S. E. Ashbrook, J. M. Griffin, J. R. Yates, F. Mauri, C. J. Pickard,
  "First-principles calculation of NMR parameters using the gauge including
  projector augmented wave method: a chemist's point of view", *Chem. Rev.*
  **112**, 5733 (2012). *(why δ = a·σ + b is fitted, and over what)*
- A. Sadoc, M. Body, C. Legein, M. Biswal, F. Fayon, X. Rocquefelte,
  F. Boucher, "NMR parameters in alkali, alkaline earth and rare earth
  fluorides from first principle calculations", *Phys. Chem. Chem. Phys.*
  **13**, 18539 (2011). *(¹⁹F shielding–shift lines over fluoride standards)*
- J. S. McCloy *et al.*, "Fluorine speciation in alkali/alkaline-earth
  aluminosilicate glasses", *Inorg. Chem.* **63**, 4669 (2024). *(the fluoride
  standards and their measured shifts)*
- S. Sturniolo, T. F. G. Green, R. M. Hanson, M. Zilka, K. Refson, P. Hodgkinson,
  S. P. Brown, J. R. Yates, "Visualization and processing of computed
  solid-state NMR parameters: MagresView and MagresPython", *Solid State Nucl.
  Magn. Reson.* **78**, 64 (2016). *(the magres file format)*
- D. J. Srivastava *et al.*, mrsimulator, *J. Chem. Phys.* **160**, 234110
  (2024). *(the EFG → C_Q conversion and the CSA lineshape)*

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
