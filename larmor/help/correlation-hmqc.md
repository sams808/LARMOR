# HMQC & correlation (2D)

> Heteronuclear- and homonuclear-**correlation** experiments answer *which
> nucleus is near which*. An HMQC cross-peak says two spins are coupled — through
> a bond (J) or through space (dipolar) — and a DQ-SQ map reports spatial
> proximity within one nuclear species. LARMOR treats a correlation map as a
> **2D whose projections you compare against separately-measured 1D spectra**, and
> it can subtract the correlated part to isolate what does *not* correlate.

---

## 1 · What the experiment tells you

| Experiment | Coupling used | Cross-peak means |
|---|---|---|
| **J-HMQC** | scalar $J$ (through-bond) | the two nuclei are **bonded** |
| **D-HMQC** | dipolar (through-space) | the two nuclei are **spatially close** (≲ a few Å) |
| **DQ-SQ** (e.g. ¹¹B, ¹H) | homonuclear dipolar | the two like-spins are **proximate**; the F1 (DQ) coordinate is the **sum** of the two SQ shifts |

In an HMQC, the **indirect (F1) dimension carries the indirectly-detected
nucleus** (often ¹⁴N, ²⁷Al, ³⁵Cl) reconstructed through the sensitive spin, which
is why F1 is usually the ugly, low-resolution axis (Gan 2000; Cavadini 2010). In a
DQ-SQ map the diagonal is where a spin correlates with *itself*; off-diagonal
pairs are distinct-site proximities, and an on-diagonal DQ peak reports a like
pair (Brown 2007).

---

## 2 · Open & overlay a 1D on a projection

Open the `2rr` on the **contour map**. Press **HMQC**, then **F2 1D…** /
**F1 1D…**: click the spectrum in the left **Explorer** (or Browse to one
elsewhere). It is overlaid on that projection and highlighted in the axis colour
(**F2 = orange, F1 = purple**). This lets you check, directly, whether the peaks
in a single-pulse 1D of that nucleus all appear in the correlation — i.e. whether
*everything* is correlated, or something is missing.

## 3 · Scale the overlay

The overlay is peak-matched automatically. Press **fit scale** for a robust
least-squares scaling of the projection onto the 1D over the visible range, or
type a manual multiplier. Correct scaling matters for the next step, because the
subtraction is only meaningful once the correlated intensity is matched.

## 4 · Extract the un-correlated features

**uncorrelated F2 →** / **uncorrelated F1 →** sends

$$r(\nu) = S_\text{1D}(\nu) - k\,P_\text{proj}(\nu)$$

(the measured 1D minus the scaled projection $P_\text{proj}$) to a **new
workspace**. The projection is the *correlated* signal, so the residual $r$ is
exactly what does **not** produce a cross-peak — an unbonded / isolated species,
an impurity, a site with no near neighbour of the partner nucleus. Fit $r$ like
any 1D. This is the practical core of correlation-editing: use the correlation to
label the connected signal, and read the difference to isolate the rest.

## 5 · Navigate

Extracting a trace or a difference opens a **new workspace**; the map stays in
its own. Switch via the **Workspaces** dock or **Back to 2D map** (Ctrl + 2).
General 2D display, phasing and measurement are in the **2D processing** manual.

---

## 6 · Where this is going (multi-experiment correlation)

The (1D − projection) difference above is the two-dataset special case of a more
general idea we care about a lot: given **N** datasets that share a nucleus or a
dimension — any mix of 1D, MQMAS, HMQC, REDOR — align their axes and decompose
features into **correlated vs un-correlated across arbitrary combinations**:

- *1D + HMQC* — what correlates heteronuclearly; the rest is isolated.
- *1D + MQMAS* — sites resolved by the isotropic axis vs lumped in the 1D.
- *1D + HMQC + REDOR* — species that correlate **and** dipolar-dephase → assign,
  then subtract to isolate the remainder.

The engine for this is scheduled deliberately **last** on the roadmap; the
architecture may be built ahead of time but it is not wired into the app yet.

---

## References

- A. Lesage, D. Sakellariou, S. Steuernagel, L. Emsley, "Carbon–proton chemical
  shift correlation in solid-state NMR by through-bond multiple-quantum
  spectroscopy", *J. Am. Chem. Soc.* **120**, 13194 (1998).
- Z. Gan, "Isotropic NMR spectra of half-integer quadrupolar nuclei using inverse
  high-resolution detection", *J. Am. Chem. Soc.* **122**, 3242 (2000). *(¹⁴N HMQC)*
- N. Cavadini, "Indirect detection of nitrogen-14 in solid-state NMR
  spectroscopy", *Prog. NMR Spectrosc.* **56**, 46 (2010). *(D-HMQC review)*
- J. Trébosc, B. Hu, J. P. Amoureux, Z. Gan, "Through-space ¹H–X heteronuclear
  correlation with the D-HMQC sequence", *J. Magn. Reson.* **186**, 220 (2007).
- S. P. Brown, "Applications of high-resolution ¹H solid-state NMR", *Solid State
  Nucl. Magn. Reson.* **41**, 1 (2012). *(DQ-SQ)*
- T. Gullion, J. Schaefer, "Rotational-echo double-resonance NMR" (REDOR),
  *J. Magn. Reson.* **81**, 196 (1989).

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
