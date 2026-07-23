# Relaxation series (T₁ / T₂)

> A relaxation experiment is a **pseudo-2D**: the same spectrum, re-acquired at a
> list of delays (`vdlist`/`vclist`), so an integral builds up or decays with the
> delay. **Tools ▸ Relaxation** opens a guided, TopSpin-style workflow — process
> the slices, drag integration zones, and fit the build-up to a relaxation model
> with an uncertainty on every τ.

---

## 1 · What is being measured

| Kind | Delay is… | Model | Reads out |
|---|---|---|---|
| **Saturation recovery** | recovery time after saturation | $I(\tau)=I_0\big(1-e^{-\tau/T_1}\big)$ | longitudinal $T_1$ |
| **Inversion recovery** | recovery after a 180° | $I(\tau)=I_0\big(1-A\,e^{-\tau/T_1}\big)$ | longitudinal $T_1$ |
| **CPMG / echo** | echo evolution time | $I(t)=I_0\,e^{-t/T_2}$ | transverse $T_2$ |
| **T₁ρ** | spin-lock time | $I(t)=I_0\,e^{-t/T_{1\rho}}$ | rotating-frame $T_{1\rho}$ |

The kind is **auto-detected from the pulse program**. $T_1$ (Vold *et al.* 1968)
governs how long you must wait between scans; $T_2$ (from the CPMG train of Carr
& Purcell 1954 / Meiboom & Gill 1958) governs linewidth and the usable echo
length; $T_{1\rho}$ probes mid-kHz motion.

### The fitting form

LARMOR uses the robust **3-parameter** build-up that TopSpin's `t1/t2` module
uses,

$$I(\tau) = I_0\Big(1 - \text{sat}\cdot e^{-(\tau/T_1)^{\beta}}\Big),$$

where **sat** floats to absorb imperfect saturation/inversion and a baseline
offset (rather than forcing the curve through zero), and the **stretch**
$\beta\le 1$ is fixed at 1 (mono-exponential) by default. Tick **stretched β**
for a **distribution** of relaxation times — the Kohlrausch–Williams–Watts form
(Williams & Watts 1970), appropriate for glasses and other disordered solids
where a single $T_1$ is a fiction. The fitted formula is shown above the curve.

---

## 2 · The guided workflow

1. **Process slices** — every slice is EM/ZF/FT'd and phased on the
   **maximum-signal** slice (not the last, which in a `ser` is often a trailing
   dummy row); the relaxed spectrum is shown. LARMOR removes the group delay
   (GRPDLY), applies FCOR, and uses an ascending frequency axis so the build-up
   is clean and monotonic.
2. **Slices [first]…[last]** — restrict the range to trim early-stopped or
   already-relaxed acquisitions (useful when fewer real points were acquired than
   the `TD1` claims).
3. **Integration zones** — drag one or more regions over your peak(s); each
   becomes an independent build-up curve, so overlapping sites can be separated
   spatially before fitting.
4. **Fit T₁/T₂** — mono-exponential by default; tick **stretched β** for a
   distribution. Click a point to **exclude an outlier** and it refits live.
   Toggle **log delay** and use **Fit view** to rescale the axis to the points.
   Each zone reports **τ ± error**.
5. **Copy CSV** — the delays, integrals and fitted τ for your notes or a figure.

---

## 3 · Beyond a single series

- **Per-site relaxation** (Tools ▸ Per-site relaxation) — instead of integrating a
  region, decompose *every slice* on the **current fit's lineshapes** (NNLS on the
  fixed model) and follow each **site's** amplitude → a $T_1$/$T_2$ per resolved
  site, even where the sites overlap.
- **Variable temperature** (Tools ▸ Variable temperature) — feed a τ(T) series and
  fit an **Arrhenius** $\tau = \tau_0\,e^{E_a/RT}$ or **Vogel–Fulcher–Tammann**
  $\tau = \tau_0\,e^{B/(T-T_0)}$ law to extract an activation energy $E_a$ (VFT for
  the super-Arrhenius slowing near a glass transition).

---

## References

- H. Y. Carr, E. M. Purcell, "Effects of diffusion on free precession in NMR
  experiments", *Phys. Rev.* **94**, 630 (1954). *(CPMG)*
- S. Meiboom, D. Gill, "Modified spin-echo method for measuring nuclear
  relaxation times", *Rev. Sci. Instrum.* **29**, 688 (1958).
- R. L. Vold, J. S. Waugh, M. P. Klein, D. E. Phelps, "Measurement of spin
  relaxation in complex systems", *J. Chem. Phys.* **48**, 3831 (1968).
  *(inversion recovery)*
- J. L. Markley, W. J. Horsley, M. P. Klein, "Spin-lattice relaxation
  measurements in slowly relaxing complex spectra", *J. Chem. Phys.* **55**, 3604
  (1971).
- G. Williams, D. C. Watts, "Non-symmetrical dielectric relaxation behaviour…",
  *Trans. Faraday Soc.* **66**, 80 (1970). *(stretched-exponential / KWW)*
- H. Vogel, *Phys. Z.* **22**, 645 (1921); G. S. Fulcher, *J. Am. Ceram. Soc.*
  **8**, 339 (1925). *(VFT law)*

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
