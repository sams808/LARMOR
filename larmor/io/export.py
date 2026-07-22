"""Export a fit in several formats: text, CSV, JSON recipe, and dmfit .fxmla.

- text  : columns ppm, experiment, model, one column per component — for
          plotting elsewhere (Origin, matplotlib, gnuplot).
- csv   : the parameter table (line, model, parameter, value, ± error,
          min, max, link) — for a paper's supporting information.
- json  : the LARMOR recipe (via Recipe.save) — the reproducible unit.
- fxmla : a dmfit-compatible file, so a fit made in LARMOR opens in dmfit
          (the σ = sCZ_CQ/2 convention is inverted on the way out).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from larmor.recipe import Recipe

SCZ_FROM_SIGMA = 2.0            # dmfit sCZ_CQ = 2 × mrsimulator σ (Phase 0)


# --------------------------------------------------------------------------
def export_text(recipe: Recipe, exp_ppm: np.ndarray, exp_amp: np.ndarray,
                path: str | Path) -> str:
    """Whitespace columns: ppm, experiment, model, per-component curves."""
    from larmor import engine

    x, total, per_site = engine.simulate(recipe, exp_ppm=exp_ppm)
    exp_on_x = np.interp(x, exp_ppm, exp_amp)
    residual = exp_on_x - total
    cols = [x, exp_on_x, total, residual, *per_site]
    labels = ["ppm", "experiment", "model", "residual"] + \
        [s.label or f"site{i}" for i, s in enumerate(recipe.sites)]
    lines = ["# " + "\t".join(labels)]
    for row in zip(*cols):
        lines.append("\t".join(f"{v:.6g}" for v in row))
    text = "\n".join(lines) + "\n"
    Path(path).write_text(text, encoding="utf-8")
    return text


def export_csv_params(recipe: Recipe, path: str | Path) -> str:
    """The parameter table as CSV (one row per parameter)."""
    rows = ["line,model,parameter,value,stderr,min,max,link"]
    for i, site in enumerate(recipe.sites):
        letter = _letter(i)
        for pname, p in site.params.items():
            rows.append(",".join([
                letter, site.model, pname,
                _num(p.value), _num(p.stderr),
                _num(p.min), _num(p.max),
                (p.expr or "").replace(",", ";"),
            ]))
    if recipe.fit_rmsd is not None:
        rows.append(f"# RMSD,{recipe.fit_rmsd:.6g}")
    text = "\n".join(rows) + "\n"
    Path(path).write_text(text, encoding="utf-8")
    return text


def export_json(recipe: Recipe, path: str | Path) -> None:
    recipe.save(path)


def _letter(i: int) -> str:
    from larmor.cellparse import index_to_letter

    return index_to_letter(i)


def _num(v) -> str:
    return "" if v is None else f"{v:.8g}"


# --------------------------------------------------------------------------
def export_fxmla(recipe: Recipe, exp_ppm: np.ndarray, exp_amp: np.ndarray,
                 path: str | Path) -> str:
    """Write a dmfit-compatible .fxmla.

    Czjzek → CzSimple (sCZ_CQ = 2σ, CQ = 2·sCZ_CQ), Gauss/Lorentz → Gaus/Lor.
    Other models are written as Gauss/Lorentz envelopes with a note so the
    file still opens. The experimental spectrum is embedded as a dmfit SIMP
    block on an ascending-Hz grid.
    """
    freq = recipe.larmor_frequency_MHz or 1.0
    lines_xml = []
    for i, site in enumerate(recipe.sites):
        lines_xml.append(_line_xml(site, i, freq))
    n_lines = len(lines_xml)

    header = (
        '<?xml version="1.0" encoding="utf-8" ?>\n'
        "<NMRFit>\n"
        "<FitParameters>\n"
        "\t<DMFitVersion>LARMOR-export</DMFitVersion>\n"
        "\t<FitMode>1</FitMode>\n"
        "\t<FitModeAsc>Fit 1D</FitModeAsc>\n"
        "\t<ComputeParam>\n"
        "\t\t<ShapeSize>32768</ShapeSize>\n"
        "\t\t<npasab>64</npasab>\n"
        "\t</ComputeParam>\n"
        f"\t<Dimension>F2\n"
        f"\t\t<nucleus>{recipe.nucleus}</nucleus>\n"
        f"\t\t<frequency>{freq:.6f}</frequency>\n"
        f"\t\t<nb_lines>{n_lines}</nb_lines>\n"
        f"\t\t<spinrate>{(recipe.spin_rate_Hz or 0.0) / 1000.0:.4f}</spinrate>\n"
    )
    body = "".join(lines_xml)
    footer = "\t</Dimension>\n</FitParameters>\n"
    expdata = _simp_block(recipe, exp_ppm, exp_amp, freq)
    text = header + body + footer + expdata + "</NMRFit>\n"
    Path(path).write_text(text, encoding="utf-8")
    return text


def _line_xml(site, i: int, freq_MHz: float) -> str:
    p = {k: v.value for k, v in site.params.items()}
    pos = p.get("isotropic_chemical_shift_ppm", 0.0)
    amp = p.get("amplitude", 1.0)
    name = site.label or f"line{i}"
    if site.model == "czjzek":
        sigma = p.get("sigma_Cq_MHz", 1.0)
        scz_khz = sigma * SCZ_FROM_SIGMA * 1000.0
        cq_khz = 2.0 * scz_khz
        dcs = p.get("shift_fwhm_ppm", 10.0)
        return (
            "\t\t<line>\n"
            "\t\t\t<ModelName>CzSimple</ModelName>\n"
            "\t\t\t<ModelNb>46</ModelNb>\n"
            f"\t\t\t<Name>{name}</Name>\n"
            "\t\t\t<GaussLor>\n"
            f"\t\t\t\t<amp>{amp:.6f}</amp>\n"
            f'\t\t\t\t<pos Unit="ppm">{pos:.6f}</pos>\n'
            f'\t\t\t\t<wid Unit="ppm">{dcs:.6f}</wid>\n'
            "\t\t\t\t<gl>0.5</gl>\n"
            f'\t\t\t\t<dCS Unit="ppm">{dcs:.6f}</dCS>\n'
            "\t\t\t</GaussLor>\n"
            "\t\t\t<QUAD>\n"
            f'\t\t\t\t<CQ Unit="KHz">{cq_khz:.6f}</CQ>\n'
            "\t\t\t\t<d>5</d>\n"
            f'\t\t\t\t<sCZ_CQ Unit="KHz">{scz_khz:.6f}</sCZ_CQ>\n'
            f'\t\t\t\t<CQ_max Unit="KHz">{cq_khz:.6f}</CQ_max>\n'
            "\t\t\t</QUAD>\n"
            "\t\t</line>\n"
        )
    # everything else -> a Gauss/Lorentz envelope (opens in dmfit)
    wid = p.get("shift_fwhm_ppm", 10.0)
    gl = p.get("gl", 0.5)
    note = "" if site.model == "gauss_lor" else \
        f"  <!-- LARMOR {site.model} exported as Gaus/Lor -->"
    return (
        "\t\t<line>" + note + "\n"
        "\t\t\t<ModelName>Gaus/Lor</ModelName>\n"
        "\t\t\t<ModelNb>1</ModelNb>\n"
        f"\t\t\t<Name>{name}</Name>\n"
        "\t\t\t<GaussLor>\n"
        f"\t\t\t\t<amp>{amp:.6f}</amp>\n"
        f'\t\t\t\t<pos Unit="ppm">{pos:.6f}</pos>\n'
        f'\t\t\t\t<wid Unit="ppm">{wid:.6f}</wid>\n'
        f"\t\t\t\t<gl>{gl:.4f}</gl>\n"
        "\t\t\t</GaussLor>\n"
        "\t\t</line>\n"
    )


def _simp_block(recipe: Recipe, exp_ppm: np.ndarray, exp_amp: np.ndarray,
                freq_MHz: float) -> str:
    ppm = np.asarray(exp_ppm, float)
    amp = np.asarray(exp_amp, float)
    order = np.argsort(ppm)[::-1]         # dmfit SIMP: descending ppm
    ppm, amp = ppm[order], amp[order]
    n = ppm.size
    hz = ppm * freq_MHz                   # ppm -> Hz on the transmitter ref
    x0 = float(hz[0])
    dx = float((hz[-1] - hz[0]) / (n - 1)) if n > 1 else 1.0
    sw = abs(dx) * n
    rows = "\n".join(f"{a:.6g}\t0" for a in amp)
    sample = recipe.sample or "exported from LARMOR"
    return (
        "<ExpData>\n<Data>\nSIMP\n"
        f"Comment={sample}\n"
        f"NP={n}\nX0={x0:.6g}\ndX={dx:.6g}\nSW={sw:.6g}\n"
        f"Sf={freq_MHz:.6f}\nSr=0\nTYPE=SPE\nDATA\n"
        f"{rows}\n</Data>\n</ExpData>\n"
    )


# --------------------------------------------------------------------------
FORMATS = {
    "text (.txt)": ("txt", export_text),
    "parameters (.csv)": ("csv", export_csv_params),
    "LARMOR recipe (.json)": ("json", export_json),
    "dmfit (.fxmla)": ("fxmla", export_fxmla),
}
