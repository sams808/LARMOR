"""Literature isotropic-shift ranges (with quadrupolar-product notes) for
oxide-glass species — the data behind View ▸ "Literature shift ranges".

Every number below is sourced from ONE review, so the overlay is a citable
assignment guide rather than folklore:

    M. Edén, "Probing oxide-based glass structures by solid-state NMR:
    Opportunities and limitations", J. Magn. Reson. Open 16–17, 100112
    (2023) — §5.2/§5.3 (shift systematics, Fig. 12), Table 4 (27Al P_Q per
    glass system), §6.3.1 (11B C_Q ranges).

Scope and honesty notes:
- ranges are TYPICAL spans for multi-component oxide glasses (O the sole
  anion); outliers exist — e.g. 27Al[p] shifts run significantly lower in
  aluminophosphate glasses (Edén's own caveat to his Fig. 12);
- 17O and 25Mg are omitted for now: the review gives their quadrupolar
  products (his Table 5; §7.4) but no drawable δiso spans in the text;
- a "P_Q"/"C_Q" note rides along as label/tooltip text — a quadrupolar
  width is not a shift-axis quantity, so it cannot be drawn as a span, but
  it is exactly what disambiguates e.g. BO3 (C_Q ≈ 2.4–2.8 MHz) from BO4
  (C_Q ≈ 0.2–0.8 MHz) once a peak sits in an overlapping shift region.

Qt-free and testable; the desktop overlay is only a consumer.
"""
from __future__ import annotations

CITATION = ("Edén 2023, J. Magn. Reson. Open 16–17, 100112 "
            "(§5.2–5.3, Table 4, §6.3.1)")

#: nucleus -> list of {label, lo_ppm, hi_ppm, quad, note}
REF_RANGES: dict[str, list[dict]] = {
    "27Al": [
        # δiso: §5.2 "typically observed ... in aluminosilicate glasses";
        # P_Q: min–max across the glass systems compiled in Table 4
        {"label": "Al[4]", "lo_ppm": 60.0, "hi_ppm": 75.0,
         "quad": "P_Q 5.9–10.9 MHz (system-dependent, Table 4)",
         "note": "AlO4; δ[4] depends mainly on n_Al/n_Si"},
        {"label": "Al[5]", "lo_ppm": 35.0, "hi_ppm": 40.0,
         "quad": "P_Q 5.3–10.3 MHz (Table 4)",
         "note": "AlO5 (Table 4 spans reach 33–44 ppm across systems)"},
        {"label": "Al[6]", "lo_ppm": 0.0, "hi_ppm": 12.0,
         "quad": "P_Q 4.0–8.4 MHz (Table 4)",
         "note": "AlO6"},
    ],
    "11B": [
        # §5.2: BO4 "resonates around 0 ppm", BO3 separated by 15–20 ppm;
        # C_Q ranges from §6.3.1 (crystalline + (dis)ordered borate phases)
        {"label": "B[3] (BO3)", "lo_ppm": 10.0, "hi_ppm": 20.0,
         "quad": "C_Q 2.4–2.8 MHz; η ≈ 0 (0 or 3 NBO) / 0.4–0.8 (1–2 NBO)",
         "note": "trigonal boron; 2nd-order quadrupolar lineshape"},
        {"label": "B[4] (BO4)", "lo_ppm": -3.0, "hi_ppm": 3.0,
         "quad": "C_Q 0.2–0.8 MHz (typ. 0.3–0.5) — near-Gaussian peak",
         "note": "tetrahedral boron, ≈0 ppm"},
    ],
    "29Si": [
        # §5.2/§5.3: Q4(SiO2) ≈ −110 ppm, +7–12 ppm per Qn→Qn−1;
        # Si[5] −150 ± 10; Si[6] −200 ± 20
        {"label": "Si[4] Qn", "lo_ppm": -120.0, "hi_ppm": -60.0,
         "quad": "",
         "note": "Q4 ≈ −110 (v-SiO2); each Qn→Qn−1 deshields by ≈7–12 ppm "
                 "(Q4→Q3 10–15, Q1→Q0 5–8)"},
        {"label": "Si[5]", "lo_ppm": -160.0, "hi_ppm": -140.0,
         "quad": "", "note": "SiO5, −150 ± 10 ppm (rare)"},
        {"label": "Si[6]", "lo_ppm": -220.0, "hi_ppm": -180.0,
         "quad": "", "note": "SiO6, −200 ± 20 ppm (rare)"},
    ],
    "31P": [
        # §5.3: Q3 (v-P2O5) ≈ −47 ppm; Q3→Q2→Q1 steps ≈15–20 ppm;
        # Q1→Q0 ≈7–10 ppm; Q0/Q1 with Na/Ca reach ≈0–3 ppm
        {"label": "P Qn", "lo_ppm": -50.0, "hi_ppm": 5.0,
         "quad": "",
         "note": "Q3 ≈ −47 (v-P2O5); Q3→Q2→Q1 steps ≈15–20 ppm, "
                 "Q1→Q0 ≈7–10 ppm; Q0/Q1 (Na/Ca) up to ≈0–3 ppm"},
    ],
}


def ranges_for(nucleus: str | None) -> list[dict]:
    """The literature ranges for a nucleus string as LARMOR stores it
    ('27Al', '11B', …); [] when none are compiled — callers show nothing
    rather than guessing."""
    return list(REF_RANGES.get((nucleus or "").strip(), ()))
