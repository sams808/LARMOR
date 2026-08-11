"""Literature isotropic-shift ranges (with quadrupolar-product notes) for
oxide-glass species — the data behind View ▸ "Literature shift ranges".

Every number below carries its source (per-entry ``ref`` key, resolved in
``REFS``), so the overlay is a citable assignment guide rather than
folklore. The full tables, conventions and caveats are in the manual:
Help ▸ User manuals ▸ "Literature shift ranges — data & sources".

Conventions the numbers assume (state them when comparing!):
- 27Al/11B/29Si/31P/17O/23Na/25Mg: standard aqueous references
  (Al(NO3)3, BF3·Et2O, TMS, 85% H3PO4, H2O, NaCl(aq), Mg(2+)(aq));
  NB 23Na values quoted against NaCl(s) differ by ≈ +7 ppm (Edén 2023);
- 19F: CFCl3 (values vs 1 M NaF(aq) differ by the −120 ppm secondary
  reference — Stebbins & Zheng 2000 describe exactly this trap);
- a "P_Q"/"C_Q" note rides along as label/tooltip text — a quadrupolar
  width is not a shift-axis quantity, so it cannot be drawn as a span,
  but it is what disambiguates species with overlapping shifts;
- 19F species bands are REPORTED PEAK POSITIONS ± 10 ppm (a stated guide
  width), because the primary sources quote positions, not ranges.

Qt-free and testable; the desktop overlay is only a consumer.
"""
from __future__ import annotations

#: citation keys — every REF_RANGES entry's ``ref`` resolves here
REFS: dict[str, str] = {
    "eden2023": ("Edén 2023, J. Magn. Reson. Open 16–17, 100112 "
                 "(§5.2–5.3, Tables 4–5, §6.3.1, §7.4)"),
    "dirken1997": ("Dirken et al. 1997, Chem. Phys. Lett. 266, 568 "
                   "(17O 3QMAS, Si–O–Si / Si–O–Al)"),
    "du2003": ("Du & Stebbins 2003, J. Phys. Chem. B 107, 10063 "
               "(17O NBO δiso/P_Q per alkali)"),
    "gambuzzi2014": ("Gambuzzi et al. 2014, (Ca,Na)/(Ca,Al,Na) silicate "
                     "glasses, 23Na MQMAS + GIPAW"),
    "stebbins1998": ("Stebbins 1998, Solid State Ionics 112, 137 "
                     "(23Na shift vs Na–O distance / NBO correlation)"),
    "shimoda2007": ("Shimoda et al. 2007, Am. Mineral. 92, 695 "
                    "(25Mg 3QMAS at 21.8 T, glasses + crystals)"),
    "baasner2014": ("Baasner et al. 2014, J. Non-Cryst. Solids "
                    "(19F speciation in Na/Ca aluminosilicate glasses)"),
    "kiczenski2002": ("Kiczenski & Stebbins 2002, J. Non-Cryst. Solids 306, "
                      "160 (19F in crystalline fluorides/oxyfluorides)"),
    "stebbins2000f": ("Stebbins & Zheng 2000, 19F in silicate glasses"),
}

CITATION = REFS["eden2023"]          # module default (27Al/11B/29Si/31P data)

#: nucleus -> list of {label, lo_ppm, hi_ppm, quad, note, ref}
REF_RANGES: dict[str, list[dict]] = {
    "27Al": [
        # δiso: Edén §5.2 "typically observed ... in aluminosilicate
        # glasses"; P_Q: min–max across the glass systems of his Table 4
        {"label": "Al[4]", "lo_ppm": 60.0, "hi_ppm": 75.0,
         "quad": "P_Q 5.9–10.9 MHz (system-dependent, Table 4)",
         "note": "AlO4; δ[4] depends mainly on n_Al/n_Si", "ref": "eden2023"},
        {"label": "Al[5]", "lo_ppm": 35.0, "hi_ppm": 40.0,
         "quad": "P_Q 5.3–10.3 MHz (Table 4)",
         "note": "AlO5 (Table 4 spans reach 33–44 ppm across systems)",
         "ref": "eden2023"},
        {"label": "Al[6]", "lo_ppm": 0.0, "hi_ppm": 12.0,
         "quad": "P_Q 4.0–8.4 MHz (Table 4)",
         "note": "AlO6", "ref": "eden2023"},
    ],
    "11B": [
        # §5.2: BO4 "resonates around 0 ppm", BO3 separated by 15–20 ppm;
        # C_Q ranges from §6.3.1 (crystalline + (dis)ordered borate phases)
        {"label": "B[3] (BO3)", "lo_ppm": 10.0, "hi_ppm": 20.0,
         "quad": "C_Q 2.4–2.8 MHz; η ≈ 0 (0 or 3 NBO) / 0.4–0.8 (1–2 NBO)",
         "note": "trigonal boron; 2nd-order quadrupolar lineshape",
         "ref": "eden2023"},
        {"label": "B[4] (BO4)", "lo_ppm": -3.0, "hi_ppm": 3.0,
         "quad": "C_Q 0.2–0.8 MHz (typ. 0.3–0.5) — near-Gaussian peak",
         "note": "tetrahedral boron, ≈0 ppm", "ref": "eden2023"},
    ],
    "29Si": [
        # §5.2/§5.3: Q4(SiO2) ≈ −110 ppm, +7–12 ppm per Qn→Qn−1;
        # Si[5] −150 ± 10; Si[6] −200 ± 20
        {"label": "Si[4] Qn", "lo_ppm": -120.0, "hi_ppm": -60.0,
         "quad": "",
         "note": "Q4 ≈ −110 (v-SiO2); each Qn→Qn−1 deshields by ≈7–12 ppm "
                 "(Q4→Q3 10–15, Q1→Q0 5–8)", "ref": "eden2023"},
        {"label": "Si[5]", "lo_ppm": -160.0, "hi_ppm": -140.0,
         "quad": "", "note": "SiO5, −150 ± 10 ppm (rare)", "ref": "eden2023"},
        {"label": "Si[6]", "lo_ppm": -220.0, "hi_ppm": -180.0,
         "quad": "", "note": "SiO6, −200 ± 20 ppm (rare)", "ref": "eden2023"},
    ],
    "31P": [
        # §5.3: Q3 (v-P2O5) ≈ −47 ppm; Q3→Q2→Q1 steps ≈15–20 ppm;
        # Q1→Q0 ≈7–10 ppm; Q0/Q1 with Na/Ca reach ≈0–3 ppm
        {"label": "P Qn", "lo_ppm": -50.0, "hi_ppm": 5.0,
         "quad": "",
         "note": "Q3 ≈ −47 (v-P2O5); Q3→Q2→Q1 steps ≈15–20 ppm, "
                 "Q1→Q0 ≈7–10 ppm; Q0/Q1 (Na/Ca) up to ≈0–3 ppm",
         "ref": "eden2023"},
    ],
    "17O": [
        # bridging Si–O–Si: Dirken 1997 (glass δiso 49±1, C_Q 5.1±0.3;
        # lit. C_Q 4.4–5.5); P_Q 4.9–5.8 aggregated in Edén Table 5
        {"label": "Si–O–Si (BO)", "lo_ppm": 40.0, "hi_ppm": 60.0,
         "quad": "C_Q 4.4–5.5 MHz (P_Q 4.9–5.8, Edén Table 5)",
         "note": "bridging oxygen between two SiO4", "ref": "dirken1997"},
        {"label": "Si–O–Al (BO)", "lo_ppm": 27.0, "hi_ppm": 40.0,
         "quad": "C_Q ≈ 3.2–4.2 MHz (P_Q 3.3–4.3, Edén Table 5)",
         "note": "glass δiso 33 ± 1 (Dirken); lower C_Q than Si–O–Si "
                 "is the discriminator", "ref": "dirken1997"},
        {"label": "NBO (Si–O–M)", "lo_ppm": 30.0, "hi_ppm": 75.0,
         "quad": "P_Q ≈ 1.9–2.3 MHz (borosilicates); Edén Table 5: "
                 "Si–NBO 1.7–2.8",
         "note": "STRONGLY modifier-dependent: δiso(NBO) ≈ 38 (Na), "
                 "42 (Li), 71 (K) in alkali borosilicates", "ref": "du2003"},
    ],
    "23Na": [
        # Gambuzzi 2014: soda-lime silicate δiso +4.7±2 / Ca-Al-Na −6.9±2,
        # C_Q 1.8–1.9 MHz; Stebbins 1998: δiso rises with NBO content and
        # falls with mean Na–O distance — one broad band, not per-species
        {"label": "Na (modifier)", "lo_ppm": -20.0, "hi_ppm": 10.0,
         "quad": "C_Q ≈ 1.8–1.9 MHz (±1) in (Ca,Al)Na silicate glasses",
         "note": "one distribution, not discrete sites; δiso rises with "
                 "NBO content / shorter Na–O. Vs NaCl(aq); NaCl(s)-"
                 "referenced values differ by ≈ +7 ppm",
         "ref": "gambuzzi2014"},
    ],
    "25Mg": [
        # Shimoda 2007 Table 1 (3QMAS, 21.8 T): MgSiO3 glass sites δCS
        # 1–24 ppm P_Q 3.3–7.3 (distorted MgO6); K2MgSi2O6 δCS 37 P_Q 2.8;
        # diopside (VIMg) δCS 6; åkermanite (IVMg) δCS 49; grandidierite
        # (VMg) δCS 55. Edén §7.4: recent glass fits give C_Qη 6.5–8.5 MHz.
        {"label": "Mg[6]", "lo_ppm": 0.0, "hi_ppm": 25.0,
         "quad": "P_Q 3.3–7.3 MHz (MgSiO3 glass); recent glass fits "
                 "C_Qη 6.5–8.5 MHz (Edén §7.4)",
         "note": "(distorted) octahedral Mg; diopside VIMg at 6 ppm",
         "ref": "shimoda2007"},
        {"label": "Mg[4]/Mg[5]", "lo_ppm": 30.0, "hi_ppm": 55.0,
         "quad": "P_Q ≈ 2.8–3 MHz (K2MgSi2O6 glass, åkermanite)",
         "note": "åkermanite IVMg 49 ppm; grandidierite VMg 55 ppm; "
                 "K2MgSi2O6 glass 37 ppm", "ref": "shimoda2007"},
    ],
    "19F": [
        # Baasner 2014 (peralkaline Na/Ca aluminosilicate glasses):
        # reported POSITIONS −113 / −146 / −168 / −188 / −225 ppm vs CFCl3
        # — drawn here as ±10 ppm guide bands (see module docstring).
        # Crystals for anchoring: CaF2 −108/−112, NaF −225 (Kiczenski 2002).
        {"label": "F–Ca(n)", "lo_ppm": -123.0, "hi_ppm": -103.0,
         "quad": "",
         "note": "≈ −113 ppm in glass; CaF2 (F–Ca(4)) at −108…−112",
         "ref": "baasner2014"},
        {"label": "Si–F–Na(n) / Al–F–Ca(n)", "lo_ppm": -156.0,
         "hi_ppm": -136.0, "quad": "",
         "note": "≈ −146 ppm (the two species overlap here)",
         "ref": "baasner2014"},
        {"label": "Al–F–Al", "lo_ppm": -178.0, "hi_ppm": -158.0,
         "quad": "", "note": "≈ −168 ppm", "ref": "baasner2014"},
        {"label": "Al–F–Na(n)", "lo_ppm": -198.0, "hi_ppm": -178.0,
         "quad": "", "note": "≈ −188 ppm", "ref": "baasner2014"},
        {"label": "F–Na(n)", "lo_ppm": -235.0, "hi_ppm": -215.0,
         "quad": "",
         "note": "≈ −225 ppm in glass; crystalline NaF at −225",
         "ref": "baasner2014"},
    ],
}


def ranges_for(nucleus: str | None) -> list[dict]:
    """The literature ranges for a nucleus string as LARMOR stores it
    ('27Al', '11B', …); [] when none are compiled — callers show nothing
    rather than guessing."""
    return list(REF_RANGES.get((nucleus or "").strip(), ()))


def citation_for(entry: dict) -> str:
    """The full citation string behind one REF_RANGES entry."""
    return REFS.get(entry.get("ref", ""), CITATION)
