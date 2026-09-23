"""Offscreen test of the Experimental-section window (AcquisitionTableDialog)
over a synthetic five-EXPNO 11B session: rows, the varying-column marking,
the paragraph's ranges, the clipboard copies, the CSV + TeX export, the
duplicate rule, the folder scan and the confirmed-rate override."""
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

pytest.importorskip("PySide6")
from PySide6.QtWidgets import QApplication  # noqa: E402

from test_acquisition import make_expno, series_11b  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def _col(dlg, header: str) -> int:
    for j in range(dlg.table.columnCount()):
        if dlg.table.horizontalHeaderItem(j).text().startswith(header):
            return j
    raise AssertionError(header)


def test_dialog_lists_series_highlights_varying_and_exports(qapp, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QFileDialog
    from larmor.desktop.acquisition_dialog import AcquisitionTableDialog

    root = tmp_path / "2026-01"
    expnos = series_11b(root)
    dlg = AcquisitionTableDialog(None, [str(p) for p in expnos[:4]] + [str(expnos[4] / "pdata" / "1" / "1r")])
    assert dlg.table.rowCount() == 5 and len(dlg.blocks) == 5
    d1 = dlg.table.horizontalHeaderItem(_col(dlg, "D1 (s)"))
    assert d1.text().endswith("▲") and "varies" in d1.toolTip() and "12.5 – 36" in d1.toolTip()
    ns = dlg.table.horizontalHeaderItem(_col(dlg, "NS"))
    assert not ns.text().endswith("▲") and ns.toolTip() == ""
    cell = dlg.table.item(0, _col(dlg, "D1 (s)"))
    assert cell.background().color().name().lower() != "#000000" and "varies" in cell.toolTip()
    assert dlg.table.item(0, _col(dlg, "NS")).toolTip() == ""
    assert "varying: Rotor, D1 (s), LB (Hz)" in dlg.status.text()
    assert "5 spectra · 1 nucleus" in dlg.status.text() and "MAS unconfirmed for 0" in dlg.status.text()
    para = dlg.para.toPlainText()
    assert "recycle delays of 12.5–36 s (Table S1)" in para and "LB = 0–100 Hz" in para
    # the rate is confirmed (title + booking agree): no warning glyph
    rot = dlg.table.item(0, _col(dlg, "νrot"))
    assert rot.text() == "35.7" and "⚠" not in rot.text()

    dlg.btnCopyPara.click()
    assert QApplication.clipboard().text() == para
    dlg.btnCopyTex.click()
    tex = QApplication.clipboard().text()
    assert r"\begin{tabular}" in tex and max(ord(c) for c in tex) < 0x80
    assert "Sample & EXPNO & Nucleus & Rotor & D1 (s) & LB (Hz)" in tex
    dlg.chkCompact.setChecked(False)
    dlg.btnCopyTex.click()
    assert "NS &" in QApplication.clipboard().text()

    target = tmp_path / "out" / "acq.csv"
    target.parent.mkdir()
    monkeypatch.setattr(QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(target), "")))
    dlg.btnSave.click()
    assert target.exists() and target.with_suffix(".tex").exists()
    assert target.read_text(encoding="utf-8").count("\n") == 7
    assert "wrote acq.csv and acq.tex" in dlg.status.text()
    # never over an instrument file
    monkeypatch.setattr(QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(expnos[0] / "acqus"), "")))
    before = (expnos[0] / "acqus").read_bytes()
    dlg.btnSave.click()
    assert (expnos[0] / "acqus").read_bytes() == before and "refused" in dlg.status.text()

    # a duplicate (same EXPNO through its 1r) adds no row; a CSV is skipped
    dlg.add_paths([str(expnos[0] / "pdata" / "1" / "1r"), str(tmp_path / "nothing.csv")])
    assert dlg.table.rowCount() == 5 and "skipped" in dlg.status.text()
    # remove one, then the folder scan adds it back (and nothing else twice)
    dlg.list.setCurrentRow(4)
    dlg.list.item(4).setSelected(True)
    dlg._remove()
    assert dlg.table.rowCount() == 4
    dlg.scan_folder(root)
    assert dlg.table.rowCount() == 5
    # a sample folder adds only its own EXPNOs
    other = make_expno(root, "01222026_SR31649_Base5Ca_SS_ALP", 24, nuc="11B",
                       bf1=192.430693, sf=192.431528, ns=256, d1=40.0, p1=0.425,
                       plw1=100.0, lb=100, tdeff=1024, si=32768, masr=4200,
                       title="11B with short tip angle\nRotor SR31648\nMASR 35.714 kHz")
    dlg.scan_folder(other.parent)
    assert dlg.table.rowCount() == 6 and "12.5–40 s" in dlg.para.toPlainText()
    dlg.close()


def test_dialog_spin_rates_override_and_uncertain_badge(qapp, tmp_path):
    from larmor.desktop.acquisition_dialog import AcquisitionTableDialog

    # acqus 4200 vs title 20 kHz vs booking 22 kHz: unconfirmed -> the badge
    e = make_expno(tmp_path / "2026-05", "04272026_P5-Bi8-12_SS_ALP", 3102, nuc="31P",
                   bf1=242.792156, sf=242.79194514, masr=4200, booking_khz=22,
                   title="P1(90)=3.750; 30 deg tip (pi/6) = 1.25\nRotor RS2427418\nMASR 20 kHz")
    dlg = AcquisitionTableDialog(None, [str(e)])
    rot = dlg.table.item(0, _col(dlg, "νrot"))
    assert rot.text().startswith("⚠ 22.0") and "acqus 4.2" in rot.text() and "booking 22" in rot.text()
    assert "MAS unconfirmed for 1" in dlg.status.text()
    assert "[confirm MAS rate: acqus 4.2 kHz, title 20 kHz, booking 22 kHz]" in dlg.para.toPlainText()
    dlg.close()
    # the fitted recipe confirmed 20 kHz (keyed by the path given): no badge
    dlg2 = AcquisitionTableDialog(None, [str(e / "pdata" / "1" / "1r")],
                                  spin_rates={str(e / "pdata" / "1" / "1r"): (20000.0, False)})
    rot = dlg2.table.item(0, _col(dlg2, "νrot"))
    assert rot.text() == "20.0" and "MAS unconfirmed for 0" in dlg2.status.text()
    assert "spinning rate was 20.0 kHz" in dlg2.para.toPlainText()
    assert "[confirm" not in dlg2.para.toPlainText()
    # an empty dialog is inert
    dlg3 = AcquisitionTableDialog(None, [])
    assert dlg3.table.rowCount() == 0 and not dlg3.btnSave.isEnabled()
    assert dlg3.para.toPlainText() == ""
    dlg2.close(); dlg3.close()
