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
