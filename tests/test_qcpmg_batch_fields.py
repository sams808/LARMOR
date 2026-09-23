"""Batch infinite-field: many samples x several fields in one grid."""
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

pytest.importorskip("PySide6")
from PySide6.QtWidgets import QApplication  # noqa: E402

from larmor.io import spectra  # noqa: E402
from larmor.qcpmg_fields import (FieldPoint, InfiniteFieldResult,  # noqa: E402
                                 fit_samples, report_text)


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def _write(tmp, name, nu, centre, width=40.0):
    x = np.linspace(centre - 400, centre + 400, 3000)
    y = np.exp(-((x - centre) / width) ** 2)
    p = tmp / name
    spectra.write_csv(p, x, y, {"nucleus": "35Cl", "larmor_MHz": nu})
    return str(p)


def test_fit_samples_groups_and_survives_a_bad_sample():
    rows = [("A", FieldPoint(78.35, -123.0, 0.5)),
            ("A", FieldPoint(107.8, -99.2, 1.0)),
            ("B", FieldPoint(78.35, -140.0, 0.6)),
            ("B", FieldPoint(107.8, -112.0, 0.9)),
            ("only-one", FieldPoint(78.35, -100.0, 1.0))]
    res = fit_samples(rows, spin=1.5, eta=0.7)
    assert list(res) == ["A", "B", "only-one"]          # order preserved
    assert isinstance(res["A"], InfiniteFieldResult)
    assert res["A"].delta_iso_ppm == pytest.approx(-72.5, abs=0.5)
    assert isinstance(res["only-one"], str)             # reported, not raised
    txt = report_text(res, 1.5, 0.7, "35Cl")
    assert "delta_iso" in txt and "NOT FITTED" in txt
    assert "eta (ASSUMED)" in txt                       # the caveat is stated


def test_batch_dialog_loads_a_grid_and_extrapolates(qapp, tmp_path):
    from larmor.desktop.qcpmg_batch_dialog import QcpmgBatchFieldsDialog

    d = QcpmgBatchFieldsDialog(None, "35Cl")
    d.nSamples.setValue(2)
    d.nFields.setValue(2)
    # two samples, each at two fields, with a real 1/nu0^2 dependence
    for si, (iso, cq) in enumerate(((-70.0, 3.2), (-80.0, 3.6))):
        for fi, nu in enumerate((78.3541, 107.811)):
            centre = iso - 1e6 * cq ** 2 * 0.02 / nu ** 2
            d._drop_files(si, 1 + fi,
                          [_write(tmp_path, f"s{si}f{fi}.csv", nu, centre)])
    assert len(d.cells) == 4
    d.table.item(0, 0).setText("LAW0Ca")
    d.table.item(1, 0).setText("LAW4Ca")
    d._compute()
    assert "2 of 2" in d.msg.text()
    assert set(d._results) == {"LAW0Ca", "LAW4Ca"}
    assert d._results["LAW0Ca"].delta_iso_ppm == pytest.approx(-70.0, abs=3.0)
    assert d.btnReport.isEnabled() and d.btnFig.isEnabled()
    # the column header learns the field it holds
    assert "78.354" in d.table.horizontalHeaderItem(1).text()
    d.close()


def test_batch_dialog_supervision_and_dirty_state(qapp, tmp_path):
    from larmor.desktop.qcpmg_batch_dialog import QcpmgBatchFieldsDialog

    d = QcpmgBatchFieldsDialog(None, "35Cl")
    d.nSamples.setValue(1)
    d._drop_files(0, 1, [_write(tmp_path, "a.csv", 78.3541, -120.0)])
    d._drop_files(0, 2, [_write(tmp_path, "b.csv", 107.811, -100.0)])
    d._compute()
    assert d.btnFig.isEnabled()
    before = d.cells[(0, 1)]["cg"]

    d.table.setCurrentCell(0, 1)                  # select -> supervision view
    assert d._region is not None
    d._region.setRegion((-160.0, -80.0))
    d._region_moved()
    assert d.cells[(0, 1)]["window"] == (-160.0, -80.0)
    assert d.cells[(0, 1)]["cg"] == pytest.approx(before, abs=2.0)
    # a changed window invalidates the previous fit rather than leaving it
    assert not d.btnFig.isEnabled()
    d.close()


def test_batch_figure_spec_has_one_entry_per_sample(qapp, tmp_path):
    import matplotlib
    matplotlib.use("Agg")
    from larmor import figures
    from larmor.desktop.qcpmg_batch_dialog import QcpmgBatchFieldsDialog

    d = QcpmgBatchFieldsDialog(None, "35Cl")
    d.nSamples.setValue(2)
    for si in range(2):
        for fi, nu in enumerate((78.3541, 107.811)):
            d._drop_files(si, 1 + fi,
                          [_write(tmp_path, f"x{si}{fi}.csv", nu,
                                  -120.0 + 10 * si + 20 * fi)])
    d._compute()
    merged = d._figure_spec()
    assert len(merged["samples"]) == 2
    one = d._figure_spec(only=list(d._results)[0])
    assert len(one["samples"]) == 1
    fig = figures.render(merged)          # and it actually draws
    assert len(fig.axes) == 1
    d.close()


# ---------------------------------------------------------------- fix 5
def test_batch_grid_checks_the_files_nucleus(qapp, tmp_path, monkeypatch):
    """A 35Cl CSV dropped on a 27Al grid used to be fitted with I = 5/2 (C_Q
    2x too large, report titled 27Al). The file's nucleus is now read from
    load_any's recipe and a mismatch is refused into the error list."""
    from PySide6.QtWidgets import QMessageBox

    from larmor.desktop.qcpmg_batch_dialog import (QcpmgBatchFieldsDialog,
                                                   _spin_of)

    assert _spin_of("27Al") == 2.5 and _spin_of("35Cl") == 1.5
    assert _spin_of("93Nb") == 4.5 and _spin_of("bogus") == 1.5

    d = QcpmgBatchFieldsDialog(None, "27Al")
    assert d.spin.value() == 2.5
    p = _write(tmp_path, "cl.csv", 78.3541, -120.0)          # nucleus 35Cl
    with pytest.raises(ValueError) as ei:
        d._load_cell(0, 1, p)
    assert "35Cl" in str(ei.value) and "27Al" in str(ei.value)
    shown = []
    monkeypatch.setattr(QMessageBox, "warning",
                        staticmethod(lambda *a, **k: shown.append(a[2])))
    d._drop_files(0, 1, [p])
    assert shown and "35Cl" in shown[0] and "27Al" in shown[0]
    assert not d.cells and d.spin.value() == 2.5                 # unchanged
    d.close()

    # an anonymous grid adopts the first file's nucleus, then holds it
    d = QcpmgBatchFieldsDialog(None, "")
    d._load_cell(0, 1, p)
    assert d._nucleus == "35Cl" and d.spin.value() == 1.5
    assert "35Cl" in d.lblNuc.text()
    x = np.linspace(-400, 400, 3000)
    spectra.write_csv(tmp_path / "al.csv", x, np.exp(-(x / 40.0) ** 2),
                      {"nucleus": "27Al", "larmor_MHz": 130.3})
    with pytest.raises(ValueError, match="27Al"):
        d._load_cell(0, 2, str(tmp_path / "al.csv"))
    d._load_cell(0, 2, _write(tmp_path, "cl2.csv", 107.811, -100.0))
    d._compute()
    res = d._results[list(d._results)[0]]
    from larmor.qcpmg_fields import cq_from_slope
    assert res.cq_MHz == pytest.approx(cq_from_slope(res.slope, 1.5, 0.7))
    d.close()


# ---------------------------------------------------------------- fix 6
def test_batch_cell_seeds_a_comb_from_its_envelope(qapp, tmp_path):
    from larmor import qcpmg
    from larmor.desktop.qcpmg_batch_dialog import QcpmgBatchFieldsDialog

    x = np.linspace(-300.0, 100.0, 4001)
    env = np.where(x < -105.0, np.exp(-(((x + 105.0) / 80.0) ** 2)),
                   np.exp(-(((x + 105.0) / 30.0) ** 2)))
    comb = np.zeros_like(x)
    for c in np.arange(-300.0, 100.0, 4.0):
        comb += np.exp(-((x - c) / 0.3) ** 2)
    p = tmp_path / "comb.csv"
    spectra.write_csv(p, x, env * comb, {"nucleus": "35Cl", "larmor_MHz": 78.354})
    d = QcpmgBatchFieldsDialog(None, "35Cl")
    d._load_cell(0, 1, str(p))
    cell = d.cells[(0, 1)]
    assert cell["comb"] and "spikelet" in cell["seed_note"]
    assert "spikelet" in d.table.item(0, 1).text()
    cg_env = qcpmg.centre_of_gravity(x, env, cell["window"])[0]
    assert cell["cg"] == pytest.approx(cg_env, abs=1.0)
    assert abs(qcpmg.centre_of_gravity(x, env * comb)[0] - cg_env) > 8.0
    d.close()


# ---------------------------------------------------------------- fix 7 / 8
def _write_manifold(tmp, name, nu, rot, centre, side=0.4, tail=False, **hdr):
    """A centreband with +-1 sidebands at +-nu_r/nu0, optionally with a long
    low-frequency tail; the header carries the rotor rate like a dataset
    written by Save as dataset..."""
    x = np.linspace(centre - 700.0, centre + 500.0, 6001)
    step = rot / nu
    y = np.exp(-(((x - centre) / 15.0) ** 2))
    for k in (-1, 1):
        y += side * np.exp(-(((x - centre - k * step) / 15.0) ** 2))
    if tail:
        y += 0.5 * np.where(x < centre, np.exp(-(((x - centre) / 150.0) ** 2)), 0.0)
    p = tmp / name
    meta = {"nucleus": "35Cl", "larmor_MHz": nu, "qcpmg_rotor_Hz": rot,
            "spectrum_mode": "absorption"}
    meta.update(hdr)
    spectra.write_csv(p, x, y, meta)
    return str(p)


def test_batch_cell_carries_the_rotor_rate_and_flags_the_window(qapp, tmp_path):
    """The saved dataset's qcpmg_rotor_Hz reaches the cell; the supervision
    view ticks the sidebands; a window that catches one sideband, a sloppy
    sigma (> 8 ppm) and a cut tail are flagged in the cell text and the
    report; the cell's sigma is max(jitter, drift)."""
    from larmor import qcpmg
    from larmor.desktop.qcpmg_batch_dialog import QcpmgBatchFieldsDialog

    d = QcpmgBatchFieldsDialog(None, "35Cl")
    d.nSamples.setValue(1)
    lo = _write_manifold(tmp_path, "lo.csv", 78.354, 16000.0, -100.0, tail=True)
    hi = _write_manifold(tmp_path, "hi.csv", 107.811, 20000.0, -80.0)
    d._drop_files(0, 1, [lo, hi])
    c = d.cells[(0, 1)]
    assert c["rotor_Hz"] == 16000.0 and c["magnitude"] is False
    m = c["meas"]
    assert m.ticks_ppm                                     # sidebands located
    assert c["sigma"] == pytest.approx(max(m.jitter_ppm, m.drift_ppm, 0.1))
    # drag the window over the centreband and the +1 sideband only
    d.table.setCurrentCell(0, 1)
    d._region.setRegion((-160.0, 160.0))
    d._region_moved()
    c = d.cells[(0, 1)]
    assert c["mode"] == "manual"
    assert "! window catches one sideband only" in c["meas"].flags
    assert "! window catches one sideband only" in d.table.item(0, 1).text()
    assert "!" in d.plot.plotItem.titleLabel.text
    d._compute()
    assert "! window catches one sideband only" in d.report.toPlainText()
    assert "CG(w, 1.5w, 2w, 3w) =" in d.report.toPlainText()
    # a window that cuts the long tail: the drift flag, and sigma = drift
    d._region.setRegion((-250.0, -20.0))
    d._region_moved()
    c = d.cells[(0, 1)]
    assert "! CG not converged -- window cuts the pattern" in c["meas"].flags
    assert "! CG not converged" in d.table.item(0, 1).text()
    assert c["sigma"] == pytest.approx(c["meas"].drift_ppm)
    d._compute()
    assert "! CG not converged" in d.report.toPlainText()
    # whole-manifold mode: the window becomes the full axis
    d.table.setCurrentCell(0, 2)
    d.winMode.setCurrentIndex(1)
    c2 = d.cells[(0, 2)]
    assert c2["mode"] == "manifold"
    assert c2["window"][0] == pytest.approx(c2["ppm"].min())
    assert c2["window"][1] == pytest.approx(c2["ppm"].max())
    assert c2["cg"] == pytest.approx(qcpmg.manifold_cg(c2["ppm"], c2["amp"])[0])
    d.close()


def test_batch_sigma_above_8ppm_is_flagged_in_cell_and_report(qapp, tmp_path):
    from larmor.desktop.qcpmg_batch_dialog import QcpmgBatchFieldsDialog

    d = QcpmgBatchFieldsDialog(None, "35Cl")
    d.nSamples.setValue(1)
    # a broad Gaussian (sigma 200 ppm) with the window cutting it at
    # half height on both sides: jitter 11.5 ppm, drift 0
    x = np.linspace(-900.0, 500.0, 7001)
    y = np.exp(-(((x + 150.0) / 200.0) ** 2))
    p = tmp_path / "broad.csv"
    spectra.write_csv(p, x, y, {"nucleus": "35Cl", "larmor_MHz": 78.354})
    d._drop_files(0, 1, [str(p)])
    d._drop_files(0, 2, [_write(tmp_path, "b2.csv", 107.811, -100.0)])
    d.table.setCurrentCell(0, 1)
    d._region.setRegion((-350.0, 50.0))
    d._region_moved()
    c = d.cells[(0, 1)]
    assert c["meas"].jitter_ppm > 8.0
    assert "! window sensitive" in d.table.item(0, 1).text()
    d._compute()
    assert "! window sensitive" in d.report.toPlainText()
    d.close()


# ---------------------------------------------------------------- fix 10
def test_batch_grid_knows_the_processing_mode(qapp, tmp_path):
    """All ten real datasets are magnitude spectra but the batch grid had no
    idea: the mode is read from the header (spectrum_mode, or the legacy
    ', magnitude)' sample suffix), shown as '(mc)' in the cell and printed
    in the report; a magnitude/absorption mix is NOT COMPARABLE."""
    from larmor.desktop.qcpmg_batch_dialog import QcpmgBatchFieldsDialog

    def write(name, nu, centre, **hdr):
        x = np.linspace(centre - 400, centre + 400, 3000)
        y = np.exp(-((x - centre) / 40.0) ** 2)
        meta = {"nucleus": "35Cl", "larmor_MHz": nu}
        meta.update(hdr)
        p = tmp_path / name
        spectra.write_csv(p, x, y, meta)
        return str(p)

    d = QcpmgBatchFieldsDialog(None, "35Cl")
    d.nSamples.setValue(2)
    d._drop_files(0, 1, [write("m1.csv", 78.354, -120.0,
                               spectrum_mode="magnitude(mc)"),
                         write("m2.csv", 107.811, -100.0,
                               sample="x · QCPMG sum echo (LB 51 Hz, magnitude)")])
    assert d.cells[(0, 1)]["magnitude"] is True
    assert d.cells[(0, 2)]["magnitude"] is True            # legacy header
    assert "(mc)" in d.table.item(0, 1).text()
    d._drop_files(1, 1, [write("a1.csv", 78.354, -120.0,
                               spectrum_mode="magnitude(mc)"),
                         write("a2.csv", 107.811, -100.0,
                               spectrum_mode="absorption")])
    assert d.cells[(1, 2)]["magnitude"] is False
    assert "(mc)" not in d.table.item(1, 2).text()
    d.table.item(0, 0).setText("both-mc")
    d.table.item(1, 0).setText("mixed")
    d._compute()
    txt = d.report.toPlainText()
    assert "magnitude" in txt and "all points measured on magnitude" in txt
    assert "NOT COMPARABLE" in txt
    assert "NOT COMPARABLE" in d.msg.text() and "mixed" in d.msg.text()
    assert "both-mc" not in d.msg.text().split("NOT COMPARABLE")[1]
    d.close()


# ---------------------------------------------------------------- fix 14
def test_identically_named_rows_are_fitted_separately(qapp, tmp_path):
    """Two rows both typed 'LAW0Ca' became one pooled 4-point fit (and the
    widths loop pooled their FWHMs). They are two fits, labelled apart, and
    the status line says so; the widths keys match the result keys."""
    from larmor.desktop.qcpmg_batch_dialog import QcpmgBatchFieldsDialog

    d = QcpmgBatchFieldsDialog(None, "35Cl")
    d.nSamples.setValue(3)
    for si, (iso, cq) in enumerate(((-70.0, 3.2), (-80.0, 3.6), (-60.0, 2.8))):
        for fi, nu in enumerate((78.3541, 107.811)):
            centre = iso - 1e6 * cq ** 2 * 0.02 / nu ** 2
            d._drop_files(si, 1 + fi, [_write(tmp_path, f"d{si}{fi}.csv", nu, centre)])
    d.table.item(0, 0).setText("LAW0Ca")
    d.table.item(1, 0).setText("LAW4Ca")
    d.table.item(2, 0).setText("LAW0Ca")                  # the duplicate
    d._compute()
    assert len(d._results) == 3                          # not 2 (pooled)
    assert set(d._results) == {"LAW0Ca", "LAW4Ca", "LAW0Ca [row 3]"}
    assert d.report.toPlainText().count("--- ") == 3
    assert all(len(r.points) == 2 for r in d._results.values())
    assert d._results["LAW0Ca"].delta_iso_ppm == pytest.approx(-70.0, abs=3.0)
    assert d._results["LAW0Ca [row 3]"].delta_iso_ppm == pytest.approx(-60.0, abs=3.0)
    assert set(d._widths) <= set(d._results) and len(d._widths) == 3
    assert "rows 1 and 3 share the name LAW0Ca" in d.msg.text()
    assert "fitted separately" in d.msg.text()
    # the figure spec follows the same labels
    assert [s_["label"] for s_ in d._figure_spec()["samples"]] == list(d._results)
    d.close()


# ---------------------------------------------------------------- fix 15
def test_batch_report_records_each_cells_path_window_and_referencing(qapp, tmp_path,
                                                                     monkeypatch):
    """The exported record must let a reader reproduce and audit every
    field: window, mode, LB, SR and the file it came from; an unreferenced
    axis (SR 0) is flagged while a referenced one is not."""
    from PySide6.QtWidgets import QFileDialog

    from larmor.desktop.qcpmg_batch_dialog import QcpmgBatchFieldsDialog

    def write(name, nu, centre, sr):
        x = np.linspace(centre - 400, centre + 400, 3000)
        y = np.exp(-((x - centre) / 40.0) ** 2)
        p = tmp_path / name
        spectra.write_csv(p, x, y, {"nucleus": "35Cl", "larmor_MHz": nu,
                                    "spectrum_mode": "absorption", "lb_Hz": 75.0,
                                    "sr_hz": sr, "referenced": abs(sr) > 0.5,
                                    "qcpmg_rotor_Hz": 16000.0 if nu < 100 else 20000.0})
        return str(p)

    d = QcpmgBatchFieldsDialog(None, "35Cl")
    d.nSamples.setValue(1)
    lo = write("LAW0Ca-3Cl_850_MHz.csv", 78.3541, -112.8, 3982.9)
    hi = write("LAW0Ca-3Cl_1p1GHz.csv", 107.811, -95.8, 0.0)
    d._drop_files(0, 1, [lo, hi])
    d.table.item(0, 0).setText("LAW0Ca")
    d._compute()
    out = tmp_path / "report.txt"
    monkeypatch.setattr(QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(out), "")))
    d._export_report()
    txt = out.read_text(encoding="utf-8")
    for name in ("LAW0Ca-3Cl_850_MHz.csv", "LAW0Ca-3Cl_1p1GHz.csv"):
        assert name in txt
    for (r, c), cell in d.cells.items():
        lo_w, hi_w = cell["window"]
        assert f"window {lo_w:.1f} ... {hi_w:.1f} ppm" in txt
    assert "LB 75 Hz" in txt and "SR +3982.9 Hz" in txt and "MAS 16000 Hz" in txt
    assert "! 107.8110 MHz: SR = +0.0 Hz -- unreferenced" in txt
    assert "! 78.3541 MHz" not in txt
    # the figure sidecar carries the windows too
    base = tmp_path / "fig"
    monkeypatch.setattr(QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(base), "")))
    d._export_figures()
    import json
    rec = json.loads((tmp_path / "fig.json").read_text(encoding="utf-8"))
    assert rec["samples"][0]["provenance"][0]["window"] == pytest.approx(
        list(d.cells[(0, 1)]["window"]))
    assert "LAW0Ca-3Cl_850_MHz.csv" in rec["samples"][0]["provenance"][0]["source"]
    d.close()
