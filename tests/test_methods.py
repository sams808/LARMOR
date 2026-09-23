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
