"""Copy-ready outputs: a LaTeX results table and a short methods sentence from a
finished fit — the last mile from "a fit on screen" to "text in a manuscript".

Qt-free and testable; the desktop layer just puts the strings on the clipboard.
"""
from __future__ import annotations

_MODEL_PHRASE = {
    "gauss_lor": "Gauss/Lorentz lines",
    "czjzek": "a Czjzek distribution of quadrupolar parameters",
    "ext_czjzek": "an extended (Gaussian-isotropic) Czjzek distribution",
    "quad_ct": "second-order quadrupolar central-transition lineshapes",
    "quad_csa": "combined quadrupolar + CSA lineshapes",
    "quad_first": "first-order quadrupolar lineshapes with spinning sidebands",
    "csa_mas": "CSA (Herzfeld–Berger) lineshapes",
    "amorphous": "dmfit 'Amorphous' distributions",
    "spectrum": "an experimental background component",
}

#: the site parameters that go into the table, with a display header + format
_COLS = [
    ("isotropic_chemical_shift_ppm", "δiso (ppm)", "{:.2f}"),
    ("Cq_MHz", "C_Q (MHz)", "{:.2f}"),
    ("sigma_Cq_MHz", "σ(C_Q) (MHz)", "{:.2f}"),
    ("eta", "η", "{:.2f}"),
    ("shift_fwhm_ppm", "FWHM (ppm)", "{:.1f}"),
    ("line_fwhm_ppm", "FWHM (ppm)", "{:.1f}"),
]


def _fmt(v, err, fmt):
    s = fmt.format(v)
    if err is not None:
        s += " ± " + fmt.format(err)
    return s


def latex_table(recipe: dict, quant: dict | None = None,
                caption: str = "", label: str = "tab:fit") -> str:
    """A LaTeX ``tabular`` (booktabs) of the fitted sites: the model-relevant
    parameters plus the integral population %. ``quant`` is a Report result
    (``run_quantify``) for the populations."""
    sites = recipe.get("sites", [])
    # which parameter columns actually appear in this model
    present = []
    for key, head, fmt in _COLS:
        if any(key in s.get("params", {}) for s in sites):
            present.append((key, head, fmt))
    pops = {}
    if quant:
        for r in quant.get("rows", []):
            pops[r.get("label")] = (r.get("fraction_pct"), r.get("fraction_err_pct"))

    ncol = 1 + len(present) + 1
    lines = [r"\begin{table}[h]", r"\centering",
             r"\begin{tabular}{l" + "c" * (ncol - 1) + "}", r"\toprule"]
    header = ["site"] + [h for _, h, _ in present] + ["pop. (\\%)"]
    lines.append(" & ".join(header) + r" \\")
    lines.append(r"\midrule")
    for s in sites:
        params = s.get("params", {})
        label_s = s.get("label") or s.get("model", "")
        row = [label_s]
        for key, _h, fmt in present:
            p = params.get(key)
            if p is None:
                row.append("--")
            else:
                v = p.get("value") if isinstance(p, dict) else p
                e = p.get("stderr") if isinstance(p, dict) else None
                row.append(_fmt(float(v), e, fmt))
        frac, ferr = pops.get(label_s, (None, None))
        row.append("--" if frac is None else _fmt(frac, ferr, "{:.1f}"))
        lines.append(" & ".join(str(c) for c in row) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    if caption:
        lines.append(r"\caption{" + caption + "}")
    lines.append(r"\label{" + label + "}")
    lines.append(r"\end{table}")
    return "\n".join(lines)


def methods_sentence(recipe: dict, error_method: str = "covariance") -> str:
    """A short, paper-ready methods sentence describing the fit."""
    sites = recipe.get("sites", [])
    nucleus = recipe.get("nucleus", "") or "the"
    field = recipe.get("larmor_frequency_MHz", 0.0) or 0.0
    models = []
    for s in sites:
        phrase = _MODEL_PHRASE.get(s.get("model", ""), s.get("model", ""))
        if phrase and phrase not in models:
            models.append(phrase)
    model_txt = models[0] if len(models) == 1 else \
        (" and ".join([", ".join(models[:-1]), models[-1]]) if len(models) > 1
         else "the fitted lineshapes")
    err_txt = ("the least-squares covariance" if error_method == "covariance"
               else "a Monte-Carlo (parametric bootstrap) analysis")
    field_txt = f" (Larmor frequency {field:.1f} MHz)" if field else ""
    # a Czjzek width is quoted in four incompatible conventions across the
    # literature (σ / 2σ / dmfit's displayed CQ = 4σ / P_Q = √5σ) — a paper
    # that names its convention costs one sentence and saves every reader
    # a factor-of-4 ambiguity, so the generated Methods text always does
    czjzek_txt = ""
    if any(s.get("model") in ("czjzek", "ext_czjzek", "csa_czjzek")
           for s in sites):
        czjzek_txt = (
            " Czjzek widths are reported as the distribution parameter σ "
            "(mrsimulator convention) together with the rms quadrupolar "
            "product P_Q = √5·σ; for comparison, dmfit's displayed CQ for "
            "the same fit corresponds to 4σ (2 × sCZ_CQ)."
        )
    return (
        f"The {nucleus} MAS NMR spectra{field_txt} were deconvoluted into "
        f"{len(sites)} site{'s' if len(sites) != 1 else ''} using {model_txt} in "
        f"LARMOR (an open dmfit-successor built on mrsimulator and lmfit). "
        f"Isotropic chemical shifts, quadrupolar parameters and relative "
        f"populations (integrated over the fit window) are reported with "
        f"uncertainties from {err_txt}." + czjzek_txt
    )
