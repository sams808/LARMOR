"""larmor.fittable: one parameter table from several fits (Qt-free).

Pinned against the shipped example fits (two LARMOR recipes, two dmfit
files), a synthetic mixed-model set for the column union, the paramstatus
glyphs and the CSV / TSV outputs."""
import csv
from pathlib import Path

import pytest

from larmor import fittable as F
from larmor import paramstatus

ROOT = Path(__file__).resolve().parents[1]
REC_AL = ROOT / "examples" / "pCABS2-4_27Al.recipe.json"
REC_B = ROOT / "examples" / "pCABS2-4_11B.recipe.json"
FXML_B = ROOT / "examples" / "pCABS2-4" / "1118" / "pdata" / "1" / "25_fit2021.fxml"
FXML_AL = ROOT / "examples" / "pCABS2-4" / "3616" / "pdata" / "1" / "1r.fxml"


def _param(value, **kw):
    p = {"value": value, "stderr": None, "vary": True, "min": None, "max": None,
         "expr": None}
    p.update(kw)
    return p


def _mixed_entries():
    """A gauss_lor fit and a czjzek fit: their parameter sets overlap only
    in position / width / amplitude, so the wide table's column union and
    its empty cells are both exercised."""
    gl = {"nucleus": "27Al", "larmor_frequency_MHz": 130.3, "spin_rate_Hz": 20000.0,
          "fit_rmsd": 0.02, "source_path": "src_gl",
          "sites": [{"model": "gauss_lor", "label": "A", "params": {
              "isotropic_chemical_shift_ppm": _param(60.0, stderr=0.3),
              "shift_fwhm_ppm": _param(10.0, stderr=0.4),
              "gl": _param(0.5, vary=False, stderr=0.0),       # fixed: †, no error bar
              "amplitude": _param(1.0, stderr=0.01)}}]}
    cz = {"nucleus": "27Al", "larmor_frequency_MHz": 130.3, "spin_rate_Hz": 20000.0,
          "fit_rmsd": 0.03, "source_path": "src_cz",
          "sites": [{"model": "czjzek", "label": "", "params": {
              "isotropic_chemical_shift_ppm": _param(64.0, stderr=0.5),
              "sigma_Cq_MHz": _param(0.05, stderr=None, min=0.05),   # at its bound: ‡
              "shift_fwhm_ppm": _param(12.0, stderr=0.9),
              "amplitude": _param(2.0, stderr=0.05, expr="2 * s0.amplitude")}},
                    {"model": "czjzek", "label": "B", "params": {
              "isotropic_chemical_shift_ppm": _param(30.0, stderr=0.7),
              "sigma_Cq_MHz": _param(1.2, stderr=0.1, min=0.05),
              "shift_fwhm_ppm": _param(15.0, stderr=1.1),
              "amplitude": _param(0.5, stderr=0.02)}}]}
    return [F.entry_from_recipe("glass A", gl), F.entry_from_recipe("glass B", cz)]


# ------------------------------------------------------------------ loading
def test_load_the_shipped_fits():
    al = F.load_fit_file(REC_AL)
    assert (al.name, al.nucleus) == ("pCABS2-4_27Al", "27Al")
    assert al.larmor_MHz == pytest.approx(130.3175616)
    assert al.spin_rate_Hz == 26000.0
    assert al.rmsd == pytest.approx(0.04677, abs=1e-4)      # the README anchor
    assert al.source == str(REC_AL)
    assert [s["model"] for s in al.sites] == ["czjzek"] * 3

    b = F.load_fit_file(REC_B)
    assert b.rmsd == pytest.approx(0.00372, abs=1e-4)
    assert [s["model"] for s in b.sites] == ["amorphous", "amorphous", "gauss_lor"]

    # dmfit files carry no RMSD; the nucleus / field / rate come from the file
    fx = F.load_fit_file(FXML_B)
    assert fx.name == "25_fit2021" and fx.rmsd is None
    assert (fx.nucleus, fx.spin_rate_Hz) == ("11B", 20000.0)
    assert len(fx.sites) == 3
    fx2 = F.load_fit_file(FXML_AL)
    assert fx2.name == "1r" and fx2.nucleus == "27Al" and len(fx2.sites) == 3


def test_load_rejects_what_neither_reader_accepts(tmp_path):
    bad = tmp_path / "x.recipe.json"
    bad.write_text("not json", encoding="utf-8")
    with pytest.raises(Exception):
        F.load_fit_file(bad)


def test_fit_name_strips_only_a_fit_suffix():
    assert F.fit_name("a/b/pCABS2-4_27Al.recipe.json") == "pCABS2-4_27Al"
    assert F.fit_name("1r.fxml") == "1r"
    assert F.fit_name("x.fxmla") == "x"
    assert F.fit_name("plain.json") == "plain"
    assert F.fit_name("notes.txt") == "notes.txt"


def test_entry_from_recipe_reads_the_facts_and_defaults_the_source():
    e = F.entry_from_recipe("ws", {"nucleus": "11B", "larmor_frequency_MHz": 160.0,
                                   "spin_rate_Hz": 0.0, "fit_rmsd": "nan",
                                   "source_path": "d/1r", "sites": []})
    assert e.rmsd is None and e.source == "d/1r" and e.nucleus == "11B"
    assert F.entry_from_recipe("ws", {}, source="given").source == "given"
    assert F.entry_from_recipe("", None).name == "fit"


# --------------------------------------------------------------------- wide
def test_wide_columns_are_the_union_in_curated_order():
    headers, rows = F.build_wide(_mixed_entries())
    assert headers[:3] == ["fit", "site", "model"]
    assert headers[-5:] == ["RMSD", "nucleus", "ν0 (MHz)", "νrot (Hz)", "source"]
    params = headers[3:-5]
    # PARAM_COLUMNS order: amplitude, position, width, gl, σ(Cq); no η, no Cq
    assert params == ["Amplitude", "Position (ppm)", "Width (ppm)", "xG/(1-x)L",
                      "σ(Cq) (MHz)"]
    assert len(rows) == 3                       # one row per (fit, site)
    assert [r[0] for r in rows] == ["glass A", "glass B", "glass B"]
    assert [r[1] for r in rows] == ["A", "A", "B"]   # an unlabelled site: its letter
    assert [r[2] for r in rows] == ["gauss_lor", "czjzek", "czjzek"]
    gl_col, sig_col = headers.index("xG/(1-x)L"), headers.index("σ(Cq) (MHz)")
    assert rows[0][sig_col] == ""               # gauss_lor has no σ(Cq)
    assert rows[1][gl_col] == ""                # czjzek has no gl
    rmsd, nuc, nu0, rot, src = rows[0][-5:]
    assert (rmsd, nuc, nu0, rot, src) == ("0.02", "27Al", "130.300", "20000", "src_gl")


def test_wide_cells_carry_value_error_and_the_status_glyphs():
    headers, rows = F.build_wide(_mixed_entries())
    pos = headers.index("Position (ppm)")
    gl_col = headers.index("xG/(1-x)L")
    sig = headers.index("σ(Cq) (MHz)")
    amp = headers.index("Amplitude")
    assert rows[0][pos] == "60 ± 0.3"
    assert rows[0][gl_col] == "0.5 " + paramstatus.MARK["fixed"]      # held: no ± 0
    assert rows[1][sig] == "0.05 " + paramstatus.MARK["at_bound"]     # at min 0.05
    assert rows[1][amp] == "2 ± 0.05 " + paramstatus.MARK["linked"]
    assert rows[2][sig] == "1.2 ± 0.1"                                # free, unmarked


def test_wide_on_the_shipped_fits_and_unknown_keys_get_a_header():
    entries = [F.load_fit_file(REC_AL), F.load_fit_file(REC_B), F.load_fit_file(FXML_AL)]
    headers, rows = F.build_wide(entries)
    assert len(rows) == 9
    assert "Cq (MHz)" in headers and "σ(Cq) (MHz)" in headers
    assert headers.index("σ(Cq) (MHz)") < headers.index("Cq (MHz)")   # curated order
    assert rows[8][0] == "1r" and rows[8][headers.index("RMSD")] == ""
    # a key the curated list does not know is appended, headed from the model
    e = F.entry_from_recipe("odd", {"sites": [{"model": "gauss_lor", "params": {
        "amplitude": _param(1.0), "mystery_key": _param(3.0)}}]})
    h2, r2 = F.build_wide([e])
    assert h2[3:-5] == ["Amplitude", "mystery_key"]
    assert r2[0][4] == "3"
    assert F.auto_header("czjzek", "sigma_Cq_MHz") == "sigma\n(MHz)"
    assert F.auto_header("no_such_model", "k") == "k"
    assert F.header_text("Position\n(ppm)") == "Position (ppm)"


def test_format_value_rules():
    st_free = paramstatus.ParamStatus("free")
    assert F.format_value(_param(1.23456789, stderr=0.012), st_free) == "1.2346 ± 0.012"
    assert F.format_value(_param(1.0, stderr=0.0), st_free) == "1"     # no zero error bar
    assert F.format_value(_param(1.0, stderr=0.2), paramstatus.ParamStatus("fixed")) == "1 †"
    assert F.format_value(_param(1.0, stderr=0.2), paramstatus.ParamStatus("default")) == "1"
    assert F.format_value(_param(None), st_free) == ""
    assert F.format_value("junk", st_free) == ""


# --------------------------------------------------------------------- long
def test_long_rows_one_per_parameter_with_the_csv_vocabulary():
    headers, rows = F.build_long(_mixed_entries())
    assert headers == F.LONG_HEADERS
    assert len(rows) == 4 + 4 + 4
    by = {(r[0], r[1], r[3]): r for r in rows}
    fixed = by[("glass A", "A", "gl")]
    assert fixed[4] == 0.5 and fixed[5] is None and fixed[6] == "fixed"
    bound = by[("glass B", "A", "sigma_Cq_MHz")]
    assert bound[6] == "at_min" and bound[7] == 0.05 and bound[8] is None
    linked = by[("glass B", "A", "amplitude")]
    assert linked[6] == "linked" and linked[9] == "2 * s0.amplitude"
    free = by[("glass B", "B", "sigma_Cq_MHz")]
    assert free[4] == 1.2 and free[5] == 0.1 and free[6] == ""


# ------------------------------------------------------------------ outputs
def test_csv_round_trip_wide_and_long(tmp_path):
    entries = _mixed_entries() + [F.load_fit_file(REC_B)]
    for build in (F.build_wide, F.build_long):
        headers, rows = build(entries)
        out = tmp_path / f"{build.__name__}.csv"
        F.to_csv(headers, rows, out)
        with open(out, encoding="utf-8", newline="") as f:
            back = list(csv.reader(f))
        assert back[0] == list(headers)
        assert len(back) == len(rows) + 1
        for got, want in zip(back[1:], rows):
            assert got == [F.cell_text(c) for c in want]
        # the glyphs and Greek headers survive the UTF-8 file
    text = (tmp_path / "build_wide.csv").read_text(encoding="utf-8")
    assert "σ(Cq) (MHz)" in text and paramstatus.MARK["fixed"] in text


def test_tsv_has_headers_and_one_line_per_row():
    headers, rows = F.build_wide(_mixed_entries())
    tsv = F.to_tsv(headers, rows)
    lines = tsv.split("\n")
    assert lines[0].split("\t") == headers
    assert len(lines) == len(rows) + 1
    assert lines[1].split("\t")[0] == "glass A"
    h2, r2 = F.build_long(_mixed_entries())
    l2 = F.to_tsv(h2, r2).split("\n")
    assert l2[0].split("\t") == list(h2)
    assert l2[1].split("\t")[4] == "60"          # a float cell as %.8g, None as ''
    assert l2[1].split("\t")[7] == ""


def test_param_columns_is_shared_with_the_desktop_table():
    """The desktop lines table imports the curated list from here (no second
    copy to drift); the automatic-fallback contract is documented on it."""
    pytest.importorskip("PySide6")
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from larmor.desktop import table
    assert table.PARAM_COLUMNS is F.PARAM_COLUMNS
    assert table._auto_header is F.auto_header
