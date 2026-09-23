"""Export a fit in several formats: text, CSV, JSON recipe, and dmfit .fxmla.

- text  : columns ppm, experiment, model, one column per component — for
          plotting elsewhere (Origin, matplotlib, gnuplot).
- csv   : the parameter table (line, model, parameter, value, ± error,
          min, max, link, vary, at_bound) — for a paper's supporting
          information. ``vary`` is True/False as stored; ``at_bound`` reads
          ``min`` / ``max`` when a free value sits at its effective bound
          (larmor.paramstatus, the fit's own rule) and is blank otherwise.
- json  : the LARMOR recipe (via Recipe.save) — the reproducible unit.
- fxmla : a dmfit-compatible file, so a fit made in LARMOR opens in dmfit
          (the σ = sCZ_CQ/2 convention is inverted on the way out).
- curves: the publication-bundle CSV (``export_curves_csv``, not in FORMATS):
          the experiment exactly as fitted on its own ppm axis, model,
          residual and one column per component, with a ``# key=value``
          header that reopens in LARMOR as a spectrum.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from larmor import paramstatus
from larmor.recipe import Recipe

SCZ_FROM_SIGMA = 2.0            # dmfit sCZ_CQ = 2 × mrsimulator σ (Phase 0)


# --------------------------------------------------------------------------
def curve_table(recipe: Recipe, exp_ppm: np.ndarray, exp_amp: np.ndarray, *,
                on_experiment_axis: bool = False, site_prefix: bool = False,
                ) -> tuple[list[str], list[np.ndarray]]:
    """The curves of a fit as ``(labels, columns)``: ppm, experiment, model,
    residual (= experiment - model), then one column per site.

    ``on_experiment_axis=False`` (export_text's layout): the MODEL axis
    (``engine.simulate``'s, a Czjzek kernel axis covering the data) with the
    experiment interpolated onto it. ``True``: the EXPERIMENTAL axis sorted
    ascending (``argsort(kind="stable")``), the experiment never interpolated
    or otherwise altered, and the model/components brought onto it with the
    same ``np.interp`` call ``fit.fit`` uses for its RMSD -- so a file of
    these columns reproduces ``Recipe.fit_rmsd``. ``site_prefix`` names the
    components ``s<i>_<label>`` (pandas-friendly, joins onto the batch table's
    ``site`` column) instead of the bare label."""
    from larmor import engine

    exp_ppm = np.asarray(exp_ppm, float)
    exp_amp = np.asarray(exp_amp, float)
    x_model, total, per_site = engine.simulate(recipe, exp_ppm=exp_ppm)
    names = [s.label or f"site{i}" for i, s in enumerate(recipe.sites)]
    if site_prefix:
        names = [f"s{i}_{n}" for i, n in enumerate(names)]
    if on_experiment_axis:
        order = np.argsort(exp_ppm, kind="stable")
        x = exp_ppm[order]
        experiment = exp_amp[order]
        model = np.interp(x, x_model, total)
        comps = [np.interp(x, x_model, np.asarray(y, float)) for y in per_site]
    else:
        x = x_model
        experiment = np.interp(x, exp_ppm, exp_amp)
        model = total
        comps = list(per_site)
    residual = experiment - model
    return (["ppm", "experiment", "model", "residual", *names],
            [x, experiment, model, residual, *comps])


def export_text(recipe: Recipe, exp_ppm: np.ndarray, exp_amp: np.ndarray,
                path: str | Path) -> str:
    """Whitespace columns: ppm, experiment, model, per-component curves."""
    labels, cols = curve_table(recipe, exp_ppm, exp_amp)
    lines = ["# " + "\t".join(labels)]
    for row in zip(*cols):
        lines.append("\t".join(f"{v:.6g}" for v in row))
    text = "\n".join(lines) + "\n"
    Path(path).write_text(text, encoding="utf-8")
    return text


def export_curves_csv(recipe: Recipe, exp_ppm: np.ndarray, exp_amp: np.ndarray,
                      path: str | Path, *, raw: np.ndarray | None = None,
                      header: dict | None = None) -> dict:
    """The publication-bundle curves file: ``ppm, experiment,
    [experiment_raw,] model, residual, s0_<label>, s1_<label> …`` on the
    experimental axis (ascending). ``ppm``, ``experiment`` and
    ``experiment_raw`` are written with ``repr`` so they round-trip
    bit-exactly (the experiment is the array the fit saw, never re-derived);
    model, residual and components with ``.9g``. ``raw`` (the pre-baseline
    trace, same axis) adds the ``experiment_raw`` column. Header lines follow
    ``spectra.read_csv``'s ``# key=value`` convention -- sample, nucleus,
    larmor_MHz, spin_rate_Hz from the recipe plus every entry of ``header``
    (source_path, fit_window_ppm, processing, excluded…) -- so File ▸ Open
    reopens the file as the fitted spectrum. Returns the arrays written
    (``ppm``, ``experiment``, ``experiment_raw``, ``model``, ``residual``,
    ``components``, ``columns``) so a caller can compute its RMSD from
    exactly what is on disk."""
    import csv

    labels, cols = curve_table(recipe, exp_ppm, exp_amp,
                               on_experiment_axis=True, site_prefix=True)
    x, experiment, model, residual = cols[:4]
    comps = cols[4:]
    raw_col = None
    if raw is not None:
        order = np.argsort(np.asarray(exp_ppm, float), kind="stable")
        raw_col = np.asarray(raw, float)[order]
    meta = {"sample": recipe.sample or "", "nucleus": recipe.nucleus or "",
            "larmor_MHz": recipe.larmor_frequency_MHz,
            "spin_rate_Hz": recipe.spin_rate_Hz}
    meta.update(header or {})
    exact = [x, experiment] + ([raw_col] if raw_col is not None else [])
    approx = [model, residual, *comps]
    columns = (labels[:2] + (["experiment_raw"] if raw_col is not None else [])
               + labels[2:])
    with open(path, "w", newline="", encoding="utf-8") as f:
        f.write("# LARMOR curves\n")
        for k, v in meta.items():
            f.write(f"# {k}={v}\n")
        w = csv.writer(f, lineterminator="\n")
        w.writerow(columns)
        for i in range(x.size):
            w.writerow([repr(float(c[i])) for c in exact]
                       + [format(float(c[i]), ".9g") for c in approx])
    return {"columns": columns, "ppm": x, "experiment": experiment,
            "experiment_raw": raw_col, "model": model, "residual": residual,
            "components": comps}


def export_csv_params(recipe: Recipe, path: str | Path) -> str:
    """The parameter table as CSV (one row per parameter). The trailing
    ``vary`` / ``at_bound`` columns say which values were held and which
    finished pinned at a bound (``min`` / ``max``); stderr is written as
    stored (a held value's 0.0 stays 0.0 -- ``vary`` explains it)."""
    rows = ["line,model,parameter,value,stderr,min,max,link,vary,at_bound"]
    for i, site in enumerate(recipe.sites):
        letter = _letter(i)
        for pname, p in site.params.items():
            rows.append(",".join([
                letter, site.model, pname,
                _num(p.value), _num(p.stderr),
                _num(p.min), _num(p.max),
                (p.expr or "").replace(",", ";"),
                str(bool(p.vary)),
                paramstatus.param_status(site.model, pname, p).side,
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
    # Amplitude convention. dmfit's CzSimple <amp> is the *seed* (pre-Czjzek-
    # broadening) Gauss/Lor amplitude, whereas LARMOR peak-normalises the fully
    # broadened component (rendered peak = amplitude). Writing LARMOR's peak
    # amplitude straight out therefore makes dmfit render Czjzek lines far too
    # short, by a width-dependent factor. Convert peak -> seed amplitude via the
    # area-preserving Czjzek relation: seed_amp = area / (Gauss-Lor area/peak at
    # dCS). Gauss/Lorentz lines need no conversion (no quadrupolar broadening).
    amp_scale = _czjzek_amp_scales(recipe, exp_ppm)
    lines_xml = []
    for i, site in enumerate(recipe.sites):
        lines_xml.append(_line_xml(site, i, freq, amp_scale.get(i, 1.0)))
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


#: LARMOR peak-amp → dmfit CzSimple <amp>: the SAME constant the importer
#: divides by (larmor.io.fxmla.DMFIT_CZSIMPLE_AMP_RATIO), so the round trip
#: is the identity — the two sides used to carry independently-calibrated
#: numbers (3.92 out, nothing in). Gauss/Lorentz lines need no conversion
#: (ratio ≈ 1). NOTE: calibrated on one ²⁷Al glass at 195 MHz — verify if a
#: different system exports mis-scaled.
from larmor.io.fxmla import DMFIT_CZSIMPLE_AMP_RATIO as _DMFIT_CZJZEK_AMP_FACTOR  # noqa: E402,E501


def _czjzek_amp_scales(recipe, exp_ppm) -> dict:
    """Per-site Czjzek amplitude multiplier for dmfit export (see the factor
    above). exp_ppm is unused now (kept for signature stability)."""
    return {i: _DMFIT_CZJZEK_AMP_FACTOR
            for i, s in enumerate(recipe.sites)
            if s.model in ("czjzek", "czjzek_d", "ext_czjzek")}


def _line_xml(site, i: int, freq_MHz: float, amp_scale: float = 1.0) -> str:
    p = {k: v.value for k, v in site.params.items()}
    pos = p.get("isotropic_chemical_shift_ppm", 0.0)
    amp = p.get("amplitude", 1.0) * amp_scale
    name = site.label or f"line{i}"
    if site.model in ("czjzek", "czjzek_d"):
        sigma = p.get("sigma_Cq_MHz", 1.0)
        scz_khz = sigma * SCZ_FROM_SIGMA * 1000.0
        cq_khz = 2.0 * scz_khz
        dcs = p.get("shift_fwhm_ppm", 10.0)
        d = float(p.get("czjzek_d", 5.0))        # dmfit's <d>; 5 for czjzek
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
            f"\t\t\t\t<d>{d:g}</d>\n"
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
