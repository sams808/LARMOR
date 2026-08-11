# Literature shift ranges — data & sources

This page documents every number behind **View ▸ Literature shift ranges**:
the shaded assignment-guide spans drawn on the spectrum for the current
nucleus, and where each one comes from. The overlay is a *guide* built from
typical values in multi-component oxide glasses — never bounds, and never a
substitute for the primary literature of your specific system.

**Reference conventions matter.** Every range below assumes the standard
reference of its nucleus. Two classic traps, called out by the sources
themselves: ²³Na values quoted against NaCl(s) sit ≈ +7 ppm from
NaCl(aq)-referenced ones (Edén 2023), and ¹⁹F work often uses 1 M NaF(aq)
as a secondary reference at −120 ppm vs the primary CFCl₃ scale
(Stebbins & Zheng 2000) — always check which zero a paper uses.

## ²⁷Al (aluminate/aluminosilicate glasses)

| Species | δiso (ppm) | P_Q (MHz) |
|---|---|---|
| Al[4] | 60 – 75 | 5.9 – 10.9 |
| Al[5] | 35 – 40 (up to 33–44 across systems) | 5.3 – 10.3 |
| Al[6] | 0 – 12 | 4.0 – 8.4 |

δiso: Edén 2023 §5.2; P_Q aggregated across the glass systems of his
Table 4. **Caveat from the same source**: ²⁷Al shifts run *significantly
lower* in aluminophosphate glasses.

## ¹¹B (borate/borosilicate glasses)

| Species | δiso (ppm) | C_Q (MHz) |
|---|---|---|
| B[3] (BO₃) | ≈ 10 – 20 | 2.4 – 2.8 (η ≈ 0 for 0/3 NBO; 0.4–0.8 for 1–2 NBO) |
| B[4] (BO₄) | ≈ −3 – 3 | 0.2 – 0.8 (typically 0.3–0.5) |

Edén 2023 §5.2 (BO₄ "around 0 ppm", BO₃ 15–20 ppm above) and §6.3.1
(C_Q/η systematics; BO₃ shows a second-order lineshape, BO₄ near-Gaussian).

## ²⁹Si / ³¹P (Qⁿ speciation)

| Species | δiso (ppm) |
|---|---|
| Si[4] Qⁿ window | −120 … −60 (Q⁴ ≈ −110 in v-SiO₂; +7–12 ppm per Qⁿ→Qⁿ⁻¹) |
| Si[5] | −150 ± 10 (rare) |
| Si[6] | −200 ± 20 (rare) |
| P Qⁿ window | −50 … +5 (Q³ ≈ −47 in v-P₂O₅; 15–20 ppm per Q³→Q²→Q¹, 7–10 for Q¹→Q⁰) |

Edén 2023 §5.2–5.3 (his Fig. 12 compiles the same spans graphically).

## ¹⁷O (bridging / non-bridging oxygen)

| Species | δiso (ppm) | C_Q or P_Q (MHz) | Source |
|---|---|---|---|
| Si–O–Si (BO) | 40 – 60 (glass: 49 ± 1) | C_Q 4.4 – 5.5; P_Q 4.9 – 5.8 | Dirken 1997; Edén Table 5 |
| Si–O–Al (BO) | 27 – 40 (glass: 33 ± 1) | C_Q ≈ 3.2 – 4.2; P_Q 3.3 – 4.3 | Dirken 1997; Edén Table 5 |
| NBO (Si–O–M) | 30 – 75 — **strongly modifier-dependent**: ≈38 (Na), ≈42 (Li), ≈71 (K) | P_Q 1.9 – 2.3 (alkali borosilicates) | Du & Stebbins 2003 |

C_Q, not δiso, is the workhorse discriminator between BO types
(Dirken 1997). Species with well-established **P_Q but no drawable δiso
span here** (from Edén 2023 Table 5, so the overlay does not draw them):
Al–O–Al 2.0–2.8 · B–O–B 4.8–5.6 · B–O–Si 5.3–5.6 · P–O–P 7.0–8.0 ·
Al–NBO 1.5–2.3 · P–NBO 3.8–5.5 MHz.

## ²³Na (network modifier)

One broad distribution, not discrete sites:

| | δiso (ppm) | C_Q (MHz) |
|---|---|---|
| Na in silicate/aluminosilicate glasses | −20 … +10 | ≈ 1.8 – 1.9 (± 1) |

Anchors: soda-lime silicate +4.7 ± 2, Ca-Al-Na silicate −6.9 ± 2
(Gambuzzi 2014, MQMAS + GIPAW). δiso rises with NBO content and falls
with mean Na–O distance (Stebbins 1998) — the *position within* the band
is structural information, not noise.

## ²⁵Mg (network modifier)

| Species | δCS (ppm) | P_Q (MHz) |
|---|---|---|
| Mg[6] (incl. strongly distorted) | 0 – 25 | 3.3 – 7.3 (MgSiO₃ glass sites) |
| Mg[4]/Mg[5] | 30 – 55 | ≈ 2.8 – 3 |

Shimoda 2007 (3QMAS at 21.8 T), Table 1: MgSiO₃ glass sites at 1–24 ppm
(concluded to be *distorted MgO₆*, not low coordination); K₂MgSi₂O₆ glass
37 ppm; crystalline anchors diopside (ⱽᴵMg) 6 ppm, åkermanite (ᴵⱽMg)
49 ppm, grandidierite (ⱽMg) 55 ppm. Edén 2023 §7.4 adds that recent glass
fits consistently give C_Qη ≈ 6.5–8.5 MHz (earlier 3–5 MHz estimates were
biased low) — expect large widths.

## ¹⁹F (vs CFCl₃; spin-½, no quadrupolar note)

The sources report **peak positions**, not ranges — the overlay draws each
as a ± 10 ppm guide band:

| Species | position (ppm) | Anchor crystals |
|---|---|---|
| F–Ca(n) | ≈ −113 | CaF₂ (F–Ca(4)) −108 … −112 |
| Si–F–Na(n) / Al–F–Ca(n) | ≈ −146 (overlapping) | |
| Al–F–Al | ≈ −168 | |
| Al–F–Na(n) | ≈ −188 | cryolite-type Na–Al–F phases nearby |
| F–Na(n) | ≈ −225 | NaF −225 |

Baasner 2014 (peralkaline Na/Ca aluminosilicate glasses); crystal anchors
from Kiczenski & Stebbins 2002. In F-rich multicomponent glasses expect
crystalline NaF/CaF₂/cryolite contributions on top of the glassy species —
sharp lines at the anchor positions are a crystallization warning, not a
new glass species (see e.g. McCloy et al. 2024 on F-loaded
aluminoborosilicates, where CaF₂/NaF/Na₃AlF₆ appear from ≈2 mol% F).

## Sources

- M. Edén, *J. Magn. Reson. Open* **16–17**, 100112 (2023) — §5.2–5.3,
  Tables 4–5, §6.3.1, §7.4.
- P. J. Dirken, S. C. Kohn, M. E. Smith, E. R. H. van Eck, *Chem. Phys.
  Lett.* **266**, 568 (1997) — ¹⁷O MQMAS, Si–O–Si vs Si–O–Al.
- L.-S. Du, J. F. Stebbins, *J. Phys. Chem. B* **107**, 10063 (2003) —
  ¹⁷O NBO δiso/P_Q per alkali in borosilicates.
- E. Gambuzzi *et al.*, ²³Na MQMAS + GIPAW of (Ca,Na)/(Ca,Al,Na) silicate
  glasses (2014).
- J. F. Stebbins, *Solid State Ionics* **112**, 137 (1998) — ²³Na shift ↔
  Na–O distance / NBO correlations.
- K. Shimoda *et al.*, *Am. Mineral.* **92**, 695 (2007) — ²⁵Mg 3QMAS at
  21.8 T.
- J. Baasner *et al.*, *J. Non-Cryst. Solids* (2014) — ¹⁹F speciation in
  Na/Ca aluminosilicate glasses.
- T. J. Kiczenski, J. F. Stebbins, *J. Non-Cryst. Solids* **306**, 160
  (2002) — ¹⁹F in crystalline fluorides/oxyfluorides.
- J. F. Stebbins, Q. Zheng (2000) — ¹⁹F in silicate glasses; the
  NaF(aq)/CFCl₃ referencing note.
- J. S. McCloy *et al.*, *Inorg. Chem.* (2024) — F speciation and
  crystallization in complex aluminoborosilicate glasses.

*To extend*: add an entry to `larmor/refranges.py` (`REF_RANGES` + a
citation in `REFS`) — one dict per species, and keep the source honest.
