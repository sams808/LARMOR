"""Copy-ready outputs: LaTeX results table + methods sentence."""
from larmor.recipe import Recipe, SiteModel, Param
from larmor import methods, quantify


def _czjzek_recipe():
    return Recipe(nucleus="27Al", larmor_frequency_MHz=130.3, sites=[
        SiteModel(model="czjzek", label="Al(IV)", params={
            "isotropic_chemical_shift_ppm": Param(60.0),
            "sigma_Cq_MHz": Param(3.0), "line_fwhm_ppm": Param(8.0),
            "amplitude": Param(100.0)}),
        SiteModel(model="czjzek", label="Al(VI)", params={
            "isotropic_chemical_shift_ppm": Param(5.0),
            "sigma_Cq_MHz": Param(2.0), "line_fwhm_ppm": Param(6.0),
            "amplitude": Param(30.0)})])


def test_latex_table_has_sites_and_populations():
    rec = _czjzek_recipe()
    q = quantify.quantify(rec, (120.0, -40.0))
    tex = methods.latex_table(rec.to_dict(), q)
    assert r"\begin{tabular}" in tex and r"\toprule" in tex
    assert "Al(IV)" in tex and "Al(VI)" in tex
    assert "δiso (ppm)" in tex
    # a population column with numbers
    assert "pop." in tex


def test_methods_sentence_names_nucleus_model_and_errors():
    s = methods.methods_sentence(_czjzek_recipe().to_dict())
    assert "27Al" in s and "130.3 MHz" in s
    assert "Czjzek" in s and "2 sites" in s
    assert "covariance" in s
    mc = methods.methods_sentence(_czjzek_recipe().to_dict(), error_method="montecarlo")
    assert "Monte-Carlo" in mc


def test_latex_table_and_methods_sentence_carry_families():
    """N3: tagged sites add Σ family rows after a \\midrule, one
    \\multicolumn line per defined ratio and a 'family/ratio errors: <basis>'
    line; the Methods sentence names the families and the basis. Untagged
    output is byte-identical to before."""
    rec = _czjzek_recipe()
    rec.sites[0].family = "Al(IV)"
    rec.sites[1].family = "Al(VI)"
    q = quantify.quantify(rec, (120.0, -40.0))
    tex = methods.latex_table(rec.to_dict(), q)
    lines = tex.splitlines()
    i4 = next(i for i, ln in enumerate(lines) if ln.startswith("Σ Al(IV)"))
    i6 = next(i for i, ln in enumerate(lines) if ln.startswith("Σ Al(VI)"))
    assert lines[i4 - 1] == r"\midrule" and i6 == i4 + 1
    assert lines[i4].count("&") == lines[1 + lines.index(r"\midrule")].count("&")
    assert " -- & -- & -- & " in lines[i4]                  # no parameter cells
    multi = [ln for ln in lines if ln.startswith(r"\multicolumn")]
    assert any(r"$\langle$CN$\rangle$ Al" in ln and "mean Al coordination" in ln
               for ln in multi)
    assert any("family/ratio errors: independent" in ln for ln in multi)
    assert lines.index(r"\bottomrule") > lines.index(multi[-1])
    # N4 is typeset with subscripts
    assert methods._tex_name("N4") == "N$_{4}$"
    assert methods._tex_name("BO4/(BO3+BO4)") == "BO$_{4}$/(BO$_{3}$+BO$_{4}$)"
    assert methods._tex_name("mean Al coordination number") == "mean Al coordination number"

    s = methods.methods_sentence(rec.to_dict(), quant=q)
    assert "structural families (Al(IV), Al(VI))" in s
    assert "⟨CN⟩ Al (mean Al coordination number) is reported" in s
    assert "treated as independent" in s
    assert s.startswith(methods.methods_sentence(rec.to_dict()))     # appended only
    q["family_basis"] = "covariance"
    assert "covariance between line amplitudes" in methods.methods_sentence(
        rec.to_dict(), quant=q)
    q["family_basis"] = "montecarlo"
    assert "per-trial sums" in methods.methods_sentence(rec.to_dict(), quant=q)

    plain = _czjzek_recipe()
    qp = quantify.quantify(plain, (120.0, -40.0))
    assert "Σ" not in methods.latex_table(plain.to_dict(), qp)
    assert "family/ratio" not in methods.latex_table(plain.to_dict(), qp)
    base = methods.methods_sentence(plain.to_dict())
    assert methods.methods_sentence(plain.to_dict(), quant=qp) == base
    assert "structural families" not in base


def test_latex_handles_gauss_lor_without_quad_columns():
    rec = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
        SiteModel(model="gauss_lor", label="B(3)", params={
            "isotropic_chemical_shift_ppm": Param(15.0),
            "shift_fwhm_ppm": Param(6.0), "amplitude": Param(100.0),
            "gl": Param(1.0, vary=False)})])
    tex = methods.latex_table(rec.to_dict())
    assert "C_Q" not in tex          # no quadrupolar columns for a gl-only model
    assert "FWHM (ppm)" in tex


def test_latex_marks_fixed_linked_and_at_bound_with_a_footnote():
    """† held fixed (value WITHOUT '± 0.00'), ‡ finished at a bound, § linked,
    and one \\multicolumn footnote row after \\bottomrule naming the bound
    and the expression -- so a table never presents a held value as fitted."""
    rec = _czjzek_recipe()
    al4, al6 = rec.sites
    al6.params["isotropic_chemical_shift_ppm"] = Param(5.0, stderr=0.0, vary=False)
    al4.params["sigma_Cq_MHz"] = Param(1.0, min=1.0)              # at its floor
    al6.params["line_fwhm_ppm"] = Param(6.0, stderr=0.3, expr="s0.line_fwhm_ppm")
    tex = methods.latex_table(rec.to_dict())
    lines = tex.splitlines()
    row4 = next(ln for ln in lines if ln.startswith("Al(IV)"))
    row6 = next(ln for ln in lines if ln.startswith("Al(VI)"))
    assert r"5.00$^{\dagger}$" in row6 and "5.00 ± " not in row6
    assert r"1.00$^{\ddagger}$" in row4 and r"\ddagger" not in row6
    assert r"6.0 ± 0.3$^{\S}$" in row6 and r"\S" not in row4
    foot = [ln for ln in lines if ln.startswith(r"\multicolumn")]
    assert len(foot) == 1
    i_bottom, i_foot, i_end = (lines.index(r"\bottomrule"), lines.index(foot[0]),
                               lines.index(r"\end{tabular}"))
    assert i_bottom < i_foot < i_end
    assert foot[0].startswith(r"\multicolumn{5}{l}{\footnotesize ")   # site + 3 cols + pop
    assert "held fixed" in foot[0] and "lower bound 1" in foot[0]
    assert "s0.line_fwhm_ppm" in foot[0]
    assert "Al(IV) σ(C_Q) (MHz) (lower bound 1)" in foot[0]
    assert "Al(VI) FWHM (ppm) = s0.line_fwhm_ppm" in foot[0]
    assert foot[0].index(r"$\dagger$") < foot[0].index(r"$\ddagger$") < foot[0].index(r"$\S$")
    assert foot[0].endswith(r"} \\")


def test_latex_population_inherits_the_amplitude_status():
    """A fully held site (the Final2 'excluded' signature: amplitude value =
    min = max = 0, vary False) prints its population with † instead of
    looking fitted; a free amplitude leaves the population unmarked."""
    rec = _czjzek_recipe()
    rec.sites[1].params["amplitude"] = Param(0.0, stderr=0.0, vary=False, min=0.0, max=0.0)
    q = quantify.quantify(rec, (120.0, -40.0))
    tex = methods.latex_table(rec.to_dict(), q)
    row6 = next(ln for ln in tex.splitlines() if ln.startswith("Al(VI)"))
    row4 = next(ln for ln in tex.splitlines() if ln.startswith("Al(IV)"))
    assert row6.rstrip(r" \\").endswith(r"0.0$^{\dagger}$")
    assert r"\dagger" not in row4
    # fixed entries are not enumerated in the footnote (the cell glyph says
    # which); only the legend line appears, once
    assert tex.count(r"\multicolumn") == 1 and r"$\dagger$ held fixed." in tex


def test_latex_is_unchanged_for_an_all_free_recipe():
    """No marker and no footnote row unless a status occurred -- byte-identical
    to the plain table for an unconstrained fit."""
    rec = _czjzek_recipe()
    q = quantify.quantify(rec, (120.0, -40.0))
    tex = methods.latex_table(rec.to_dict(), q)
    for bad in ("dagger", "ddagger", r"\S}", r"\multicolumn"):
        assert bad not in tex
    assert tex.splitlines()[-4:-3] == [r"\bottomrule"] or r"\bottomrule" in tex
    assert tex.index(r"\bottomrule") < tex.index(r"\end{tabular}")
    # a fixed gl (not a tabulated column) does not trigger a footnote either
    gl_only = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, sites=[
        SiteModel(model="gauss_lor", label="B(3)", params={
            "isotropic_chemical_shift_ppm": Param(15.0),
            "shift_fwhm_ppm": Param(6.0), "amplitude": Param(100.0),
            "gl": Param(1.0, vary=False)})])
    tex2 = methods.latex_table(gl_only.to_dict())
    for bad in ("dagger", "ddagger", r"\S}", r"\multicolumn"):
        assert bad not in tex2


def test_methods_sentence_states_the_czjzek_width_convention():
    """A Czjzek fit's Methods text must name the width convention -- the
    same fitted width is quoted as sigma / 2 sigma (dmfit sCZ_CQ) /
    4 sigma (dmfit's CQ box) / sqrt5 sigma (P_Q) across the literature,
    so an unlabeled number is unusable by the next reader."""
    s = methods.methods_sentence(_czjzek_recipe().to_dict())
    assert "P_Q" in s and "4σ" in s and "sCZ_CQ" in s

    # a fit with no Czjzek-family site gets no convention sentence
    plain = {"nucleus": "11B", "larmor_frequency_MHz": 160.0, "sites": [
        {"model": "gauss_lor", "label": "A", "params": {}}]}
    assert "P_Q" not in methods.methods_sentence(plain)


def test_methods_paragraph_includes_acquisition_and_software():
    """The Experimental paragraph: acquisition sentences from the block (the
    recipe's confirmed MAS rate, the TopSpin and LARMOR processing), then
    the fit sentence naming the versions; methods_sentence without the
    software kw is byte-identical to before."""
    from test_acquisition import _block_3102_like

    rec = _czjzek_recipe().to_dict()
    rec.update({"nucleus": "31P", "acquisition": _block_3102_like(), "spin_rate_Hz": 20000.0,
                "mas_uncertain": False, "sr_hz": -210.86,
                "processing": [{"op": "twopoint_bg"}],
                "software": {"larmor": "0.13.0", "git_commit": "2fe4800f" + "0" * 32,
                             "mrsimulator": "1.0.0", "lmfit": "1.3.4"}})
    para = methods.methods_paragraph(rec, referencing_text=None)
    for must in ("242.79 MHz", "20.0 kHz", "300 s", "14 transients", "two-point linear baseline",
                 "LARMOR 0.13.0, commit 2fe4800", "mrsimulator 1.0.0", "lmfit 1.3.4",
                 "[state the reference standard]"):
        assert must in para, must
    assert para.rstrip().endswith("covariance.") or "covariance" in para.split(". ")[-2]
    assert para.index("14 transients") < para.index("deconvoluted")
    plain = methods.methods_sentence(rec)
    assert "LARMOR (an open dmfit-successor built on mrsimulator and lmfit)" in plain
    assert "0.13.0" not in plain
    # the recipe's own mas_rate block (confirmed) removes the bracket even
    # when the stored flag says uncertain
    rec2 = dict(rec, mas_uncertain=True,
                provenance={"mas_rate": {"acqus_Hz": 4200.0, "title_Hz": 20000.0,
                                         "booking_Hz": 22000.0, "confirmed": "2026-09-23"}})
    assert "[confirm" not in methods.methods_paragraph(rec2, referencing_text=None)
    rec3 = dict(rec, mas_uncertain=True, provenance={})
    assert "[confirm MAS rate: acqus 4.2 kHz, title 20 kHz, booking 22 kHz]" in \
        methods.methods_paragraph(rec3, referencing_text=None)


def test_methods_paragraph_without_acquisition_is_the_sentence_plus_software(monkeypatch):
    from larmor import provenance

    stamp = {"larmor": "0.13.0", "git_commit": "", "mrsimulator": "1.0.0", "lmfit": "1.3.4",
             "numpy": "2.4.6", "fitted": "2026-09-23T10:00:00"}
    monkeypatch.setattr(provenance, "software_stamp", lambda: stamp)
    rec = _czjzek_recipe().to_dict()
    rec["source_kind"] = "csv"
    para = methods.methods_paragraph(rec)
    assert para == methods.methods_sentence(rec, software=stamp)
    assert "LARMOR 0.13.0 (an open dmfit-successor built on mrsimulator 1.0.0 and lmfit 1.3.4)" in para
    for never in ("probe", "recycle", "MHz probe", "commit"):
        assert never not in para


def test_methods_sentence_names_the_profile_estimator_and_software_versions(monkeypatch):
    """'profile' used to fall into the Monte-Carlo wording; the software
    block never raises and names every field an export README records."""
    import re
    import subprocess
    import time

    import larmor

    rec = _czjzek_recipe().to_dict()
    prof = methods.methods_sentence(rec, "profile")
    assert "χ² profile" in prof and "Monte-Carlo" not in prof
    assert "covariance" in methods.methods_sentence(rec, "covariance")
    assert "Monte-Carlo" in methods.methods_sentence(rec, "montecarlo")

    t0 = time.perf_counter()
    v = methods.software_versions()
    assert time.perf_counter() - t0 < 3.0
    assert set(v) == {"larmor", "mrsimulator", "lmfit", "numpy", "python", "git_commit"}
    assert v["larmor"] == larmor.__version__
    assert v["mrsimulator"] and v["lmfit"] and v["numpy"] and v["python"]
    assert v["git_commit"] == "" or re.fullmatch(r"[0-9a-f]{12}", v["git_commit"])

    def boom(*a, **k):
        raise OSError("no git here")

    monkeypatch.setattr(subprocess, "run", boom)
    assert methods.software_versions()["git_commit"] == ""
    assert methods._dist_version("no-such-distribution-xyz") == ""


def _f19_recipe():
    return {"nucleus": "19F", "larmor_frequency_MHz": 564.27, "sites": [
        {"model": "gl_norm", "label": "F5", "params": {}},
        {"model": "gl_norm", "label": "F7", "params": {}}]}


def test_methods_sentence_states_the_dft_shift_conversion():
    """Sites seeded from GIPAW tensors: the sentence must say how the
    shieldings became shifts -- the fitted line, its references and its
    uncertainties -- or that a single sigma_ref was used."""
    base = _f19_recipe()
    before = methods.methods_sentence(base)
    rec = dict(base)
    rec["provenance"] = {"dft_import": {
        "isotope": "19F",
        "calculation": {"code": "QE-GIPAW", "version": "7.5", "xc": "PBE",
                        "cutoff_wfc_Ry": 60.0, "cutoff_rho_Ry": 480.0},
        "calibration": {"kind": "fit", "slope": -0.698, "intercept": 56.5,
                        "cov": [[1.44e-4, -0.047], [-0.047, 15.2]],
                        "n": 4, "dof": 2, "rmse_ppm": 0.9,
                        "points": ["NaF", "CaF2", "cryolite", "sulphohalite"]}}}
    s = methods.methods_sentence(rec)
    assert s.startswith(before)                    # appended, nothing altered
    assert "δ = a·σ + b" in s and "-0.698 ± 0.012" in s
    assert "4 reference compounds" in s and "QE-GIPAW 7.5" in s
    assert "NaF" in s and "sulphohalite" in s and "RMSE 0.9 ppm" in s
    # two references: no RMSE claim, the honesty clause instead
    rec["provenance"]["dft_import"]["calibration"].update(
        {"n": 2, "dof": 0, "points": ["NaF", "CaF2"]})
    s2 = methods.methods_sentence(rec)
    assert "exact by construction" in s2 and "RMSE" not in s2
    # sigma_ref
    rec["provenance"]["dft_import"]["calibration"] = {
        "kind": "sigma_ref", "slope": -1.0, "intercept": 560.0}
    s3 = methods.methods_sentence(rec)
    assert "σ_ref = 560.0 ppm" in s3 and "δ = σ_ref − σ" in s3
    # no conversion recorded at all
    rec["provenance"]["dft_import"]["calibration"] = None
    assert "unconverted" in methods.methods_sentence(rec)
    # a recipe without dft_import is byte-identical to before
    assert methods.methods_sentence(_f19_recipe()) == before
    assert methods.methods_sentence(dict(base, provenance={"referencing": {}})) == before


def test_methods_phrase_for_gl_norm_is_not_the_raw_name():
    s = methods.methods_sentence(_f19_recipe())
    assert "area-normalised" in s and "gl_norm" not in s


def test_latex_does_not_footnote_an_lb_held_at_its_model_default():
    """A fresh Czjzek site holds lb at the registry default (dmfit's greyed
    Lb): the table prints it as a plain held value -- no '± 0.00', no dagger,
    no footnote row -- because a model default is not a user constraint. A
    dagger appears only once the user pins lb at a value of their own."""
    from larmor import models

    rec = _czjzek_recipe()
    fresh = models.get("czjzek").defaults()["line_fwhm_ppm"]
    fresh.stderr = 0.0                              # as a saved recipe stores it
    rec.sites[0].params["line_fwhm_ppm"] = fresh
    tex = methods.latex_table(rec.to_dict())
    row4 = next(ln for ln in tex.splitlines() if ln.startswith("Al(IV)"))
    assert "0.0 &" in row4 and "0.0 ±" not in row4 and r"\dagger" not in row4
    assert r"\multicolumn" not in tex
    rec.sites[0].params["line_fwhm_ppm"] = Param(2.0, stderr=0.0, vary=False)
    tex2 = methods.latex_table(rec.to_dict())
    row4 = next(ln for ln in tex2.splitlines() if ln.startswith("Al(IV)"))
    assert r"2.0$^{\dagger}$" in row4 and r"$\dagger$ held fixed." in tex2
