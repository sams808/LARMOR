"""Copy-ready outputs: a LaTeX results table, a short methods sentence and the
full Experimental paragraph from a finished fit — the last mile from "a fit
on screen" to "text in a manuscript".

Also supplies the software block (``software_versions``) that exports record
next to their numbers -- the series publication bundle's README and, through
``larmor.provenance.software_stamp``, every fitted recipe.

Qt-free and testable; the desktop layer just puts the strings on the clipboard.
"""
from __future__ import annotations

from larmor import paramstatus

_MODEL_PHRASE = {
    "gauss_lor": "Gauss/Lorentz lines",
    "gl_norm": "area-normalised Gauss/Lorentz lines",
    "czjzek": "a Czjzek distribution of quadrupolar parameters",
    "czjzek_d": "a generalised (d-parameter) Czjzek distribution",
    "czjzek_corr": "a Czjzek distribution with a correlated (δiso, C_Q) "
                   "shift dependence",
    "ext_czjzek": "an extended (Gaussian-isotropic) Czjzek distribution",
    "exchange2": "two-site Bloch–McConnell exchange lineshapes",
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
    ("czjzek_d", "d", "{:.2f}"),
    ("shift_slope_ppm_per_MHz", "dδ/dC_Q (ppm/MHz)", "{:.2f}"),
    ("eta", "η", "{:.2f}"),
    ("shift_fwhm_ppm", "FWHM (ppm)", "{:.1f}"),
    ("line_fwhm_ppm", "FWHM (ppm)", "{:.1f}"),
    ("split_ppm", "Δδ (ppm)", "{:.2f}"),
    ("pop_a", "p_A", "{:.2f}"),
    ("k_ex_hz", "k_ex (s⁻¹)", "{:.3g}"),
    ("lorentz_fwhm_ppm", "Lorentz FWHM (ppm)", "{:.1f}"),
]


def _fmt(v, err, fmt, status=None):
    """``value ± err`` with the status marker appended. A HELD value prints
    without its error (saved recipes carry stderr 0.0 on fixed parameters,
    which read as a fitted ``1.00 ± 0.00``); a value that finished at a bound
    has no covariance error at all after the fit's retry, so ‡ is its
    explanation."""
    s = fmt.format(v)
    fixed = status is not None and status.kind == "fixed"
    if err is not None and not fixed:
        s += " ± " + fmt.format(err)
    if status is not None:
        s += status.latex
    return s


def _tex_name(s: str) -> str:
    """'N4' -> 'N$_{4}$', 'BO4/(BO3+BO4)' -> 'BO$_{4}$/(BO$_{3}$+BO$_{4}$)',
    '⟨CN⟩ Al' -> '$\\langle$CN$\\rangle$ Al': the family / ratio names of
    larmor.families set in LaTeX, digits glued to letters as subscripts."""
    import re

    s = str(s).replace("⟨", r"$\langle$").replace("⟩", r"$\rangle$")
    return re.sub(r"(?<=[A-Za-z])(\d+)(?!\w)", r"$_{\1}$", s)


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
    # † fixed · ‡ at a bound · § linked, derived per cell; the footnote lists
    # only the kinds that occurred, so an all-free recipe prints exactly the
    # plain table
    entries: list[tuple[str, paramstatus.ParamStatus]] = []
    for s in sites:
        params = s.get("params", {})
        label_s = s.get("label") or s.get("model", "")
        model_s = s.get("model")
        row = [label_s]
        for key, head, fmt in present:
            p = params.get(key)
            if p is None:
                row.append("--")
            else:
                v = p.get("value") if isinstance(p, dict) else p
                e = p.get("stderr") if isinstance(p, dict) else None
                st = paramstatus.param_status(model_s, key, p)
                if st.kind != "free":
                    entries.append((f"{label_s} {head}", st))
                row.append(_fmt(float(v), e, fmt, st))
        frac, ferr = pops.get(label_s, (None, None))
        # the population inherits the amplitude's status: a held (or zeroed)
        # amplitude is not a fitted population
        amp_st = paramstatus.param_status(model_s, "amplitude", params.get("amplitude"))
        if frac is not None and amp_st.kind != "free":
            entries.append((f"{label_s} pop. (\\%)", amp_st))
        row.append("--" if frac is None else
                   _fmt(frac, ferr, "{:.1f}", amp_st if amp_st.kind != "free" else None))
        lines.append(" & ".join(str(c) for c in row) + r" \\")
    # Σ family rows and named ratios (N3) -- only when a line is tagged, so
    # an untagged table is byte-identical to before
    fams = (quant or {}).get("families") or []
    if fams:
        lines.append(r"\midrule")
        for f in fams:
            row = [f"Σ {f['family']}"] + ["--"] * len(present)
            row.append(_fmt(float(f["fraction_pct"]), f.get("fraction_err_pct"),
                            "{:.1f}"))
            lines.append(" & ".join(row) + r" \\")
        for r in quant.get("ratios") or []:
            if not r.get("defined") or r.get("value") is None:
                continue
            txt = _tex_name(r["name"])
            if r.get("description"):
                txt += " = " + _tex_name(r["description"])
            txt += " = " + _fmt(float(r["value"]), r.get("err"), r.get("fmt", "{:.3f}"))
            lines.append(r"\multicolumn{" + str(ncol) + r"}{l}{" + txt + r"} \\")
        lines.append(r"\multicolumn{" + str(ncol) + r"}{l}{\footnotesize "
                     + "family/ratio errors: " + str(quant.get("family_basis", ""))
                     + r"} \\")
    lines.append(r"\bottomrule")
    if entries:
        lines.append(r"\multicolumn{" + str(ncol) + r"}{l}{\footnotesize "
                     + paramstatus.footnote(entries, "latex") + r"} \\")
    lines.append(r"\end{tabular}")
    if caption:
        lines.append(r"\caption{" + caption + "}")
    lines.append(r"\label{" + label + "}")
    lines.append(r"\end{table}")
    return "\n".join(lines)


def _dft_clause(recipe: dict) -> str:
    """One sentence on where the starting values came from when the sites
    were seeded from computed (GIPAW) shielding tensors -- and, above all,
    HOW the shieldings became shifts: the fitted line with its reference
    compounds and uncertainties, or the single sigma_ref. Empty when the
    recipe carries no ``provenance['dft_import']``."""
    d = (recipe.get("provenance") or {}).get("dft_import")
    if not d:
        return ""
    from larmor import dft
    from larmor.shiftcal import ShiftCalibration

    iso = d.get("isotope") or recipe.get("nucleus") or ""
    calc = dft.MagresCalc.from_dict(d.get("calculation") or {})
    calc_txt = calc.describe()
    head = f" Starting values came from computed {iso} shielding tensors"
    if calc_txt != "no calculation header":
        head += f" ({calc_txt})"
    cal_d = d.get("calibration")
    if not cal_d:
        return (head + "; the isotropic shieldings were used unconverted as "
                "starting positions.")
    cal = ShiftCalibration.from_dict(cal_d)
    return head + "; " + cal.methods_clause()


def _family_clause(quant: dict | None) -> str:
    """One sentence on the summed families and named ratios (N3) with the
    basis their uncertainty was propagated on; '' for an untagged fit so the
    untagged Methods text is byte-identical."""
    fams = (quant or {}).get("families") or []
    if not fams:
        return ""
    names = ", ".join(f["family"] for f in fams)
    ratios = [r for r in (quant.get("ratios") or [])
              if r.get("defined") and r.get("value") is not None]
    txt = f" Site populations were summed into structural families ({names})"
    if ratios:
        txt += ("; the ratio" + ("s " if len(ratios) > 1 else " ")
                + ", ".join(r["name"] + (f" ({r['description']})" if r.get("description")
                                          else "") for r in ratios)
                + (" are" if len(ratios) > 1 else " is") + " reported")
    basis = quant.get("family_basis", "")
    if basis == "covariance":
        txt += (" with uncertainties that propagate the covariance between "
                "line amplitudes.")
    elif basis == "montecarlo":
        txt += (" with uncertainties from the spread of per-trial sums over "
                "the Monte-Carlo refits.")
    else:
        txt += (" with uncertainties propagated to first order with the lines "
                "treated as independent.")
    return txt


def methods_sentence(recipe: dict, error_method: str = "covariance",
                     software: dict | None = None, quant: dict | None = None) -> str:
    """A short, paper-ready methods sentence describing the fit. With
    ``software`` (a recipe's ``software`` block or ``provenance.
    software_stamp()``) the LARMOR clause names the versions; the default
    keeps the version-free wording. ``quant`` (a quantify() result) adds the
    family/ratio sentence when lines are tagged; without it, or untagged, the
    text is unchanged."""
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
    err_txt = {"covariance": "the least-squares covariance",
               "profile": "one-parameter χ² profiles (1σ intervals)",
               }.get(error_method, "a Monte-Carlo (parametric bootstrap) analysis")
    field_txt = f" (Larmor frequency {field:.1f} MHz)" if field else ""
    # a Czjzek width is quoted in four incompatible conventions across the
    # literature (σ / σ_Cz = sCZ_CQ = 2σ / dmfit's displayed CQ = 4σ /
    # P_Q = 2√5·σ) — a paper that names its convention costs one sentence and
    # saves every reader a factor-of-4 ambiguity, so the generated Methods
    # text always does
    czjzek_txt = ""
    if any(s.get("model") in ("czjzek", "czjzek_d", "czjzek_corr",
                              "ext_czjzek", "csa_czjzek")
           for s in sites):
        czjzek_txt = (
            " Czjzek widths are reported as the distribution parameter σ "
            "(mrsimulator convention; the Czjzek-paper width is σ_Cz = 2σ) "
            "together with the rms quadrupolar product P_Q = 2√5·σ "
            "(= √5·σ_Cz); for comparison, dmfit's displayed CQ for the same "
            "fit corresponds to 4σ (2 × sCZ_CQ)."
        )
    larmor_txt = "LARMOR (an open dmfit-successor built on mrsimulator and lmfit)"
    if software:
        v = software.get("larmor") or ""
        commit = str(software.get("git_commit") or "")[:7]
        mrs = software.get("mrsimulator") or ""
        lm = software.get("lmfit") or ""
        larmor_txt = (f"LARMOR{' ' + v if v else ''}{', commit ' + commit if commit else ''}"
                      f" (an open dmfit-successor built on mrsimulator{' ' + mrs if mrs else ''}"
                      f" and lmfit{' ' + lm if lm else ''})")
    return (
        f"The {nucleus} MAS NMR spectra{field_txt} were deconvoluted into "
        f"{len(sites)} site{'s' if len(sites) != 1 else ''} using {model_txt} in "
        f"{larmor_txt}. "
        f"Isotropic chemical shifts, quadrupolar parameters and relative "
        f"populations (integrated over the fit window) are reported with "
        f"uncertainties from {err_txt}." + czjzek_txt + _dft_clause(recipe)
        + _family_clause(quant)
    )


def methods_paragraph(recipe: dict, error_method: str = "covariance", *,
                      block: dict | None = None,
                      referencing_text: str | None = None,
                      quant: dict | None = None) -> str:
    """The full Experimental paragraph of one fit: the acquisition and
    processing sentences from the recipe's acquisition block (``larmor.
    acquisition.sentences`` -- spectrometer and field, probe, the confirmed
    MAS rate or a bracketed candidate list, pulse program, P1 and the flip
    angle when the title supports one, PLW1, D1, NS, SW/TD, referencing with
    evidence or a bracket, TopSpin processing then the recipe's own steps)
    followed by :func:`methods_sentence` naming the software versions. For a
    CSV / dmfit / Varian source (no block) the paragraph is the fit sentence
    alone. ``block`` overrides the stored / re-read block;
    ``referencing_text`` overrides the evidence search."""
    from larmor import acquisition, provenance

    recipe = recipe or {}
    if block is None:
        block = acquisition.block_for_recipe(recipe)
    sents: list[str] = []
    if block:
        if referencing_text is None:
            referencing_text = acquisition.referencing_evidence([block], recipe=recipe)
        sents = acquisition.sentences(
            block, spin_rate_Hz=recipe.get("spin_rate_Hz"),
            mas_uncertain=recipe.get("mas_uncertain"),
            mas_rate=(recipe.get("provenance") or {}).get("mas_rate"),
            sr_hz=recipe.get("sr_hz"), referencing_text=referencing_text,
            larmor_ops=recipe.get("processing"))
    software = recipe.get("software") or provenance.software_stamp()
    return " ".join(sents + [methods_sentence(recipe, error_method,
                                              software=software, quant=quant)])


# ---------------------------------------------------------------- software
def software_versions() -> dict:
    """The software block an export records next to its numbers: LARMOR,
    mrsimulator, lmfit and numpy versions (importlib.metadata, so nothing
    heavy is imported for a string), the Python version and, when LARMOR
    runs from a git checkout, the short commit. Never raises; a version that
    cannot be determined is ``""``."""
    import platform

    import larmor

    return {"larmor": larmor.__version__,
            "mrsimulator": _dist_version("mrsimulator"),
            "lmfit": _dist_version("lmfit"),
            "numpy": _dist_version("numpy"),
            "python": platform.python_version(),
            "git_commit": _git_commit()}


def _dist_version(name: str) -> str:
    from importlib.metadata import PackageNotFoundError, version

    try:
        return version(name)
    except PackageNotFoundError:
        return ""


def _git_commit() -> str:
    """``git rev-parse --short=12 HEAD`` of the checkout larmor/ lives in --
    only when a ``.git`` sits next to the package (a worktree's ``.git`` is a
    file, hence ``exists`` rather than ``is_dir``); ``""`` for an installed
    copy, a missing git, a timeout or any other failure."""
    import subprocess
    from pathlib import Path

    import larmor

    root = Path(larmor.__file__).resolve().parents[1]
    if not (root / ".git").exists():
        return ""
    try:
        kw = {}
        if hasattr(subprocess, "CREATE_NO_WINDOW"):      # no console flash on Windows
            kw["creationflags"] = subprocess.CREATE_NO_WINDOW
        out = subprocess.run(["git", "-C", str(root), "rev-parse", "--short=12", "HEAD"],
                             capture_output=True, text=True, timeout=2, **kw)
        return out.stdout.strip() if out.returncode == 0 else ""
    except Exception:
        return ""
