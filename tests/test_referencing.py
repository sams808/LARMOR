"""Referencing audit: indirect (Xi) referencing arithmetic, the session
scan, verdicts, the TopSpin list, the permanent log, the window, the CLI and
the app's re-referencing -- plus one real NMRFAM session as an anchor."""
import json
import os
from pathlib import Path

import numpy as np
import pytest

from conftest import DATA_ROOT, require
from larmor import referencing as R

ADA = 1.82


def _jcamp(path: Path, **params):
    lines = ["##TITLE= Parameter file, TopSpin 4.1", "##JCAMPDX= 5.0",
             "##DATATYPE= Parameter Values", "##ORIGIN= Bruker BioSpin GmbH",
             "##OWNER= nmr"]
    for k, v in params.items():
        lines.append(f"##${k}= {v}")
    lines.append("##END=")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _expno(root, sample, expno, nuc, bf1, sf, date, title="", probe="Z1 3.2mm"):
    e = root / sample / str(expno)
    _jcamp(e / "acqus", NUC1=f"<{nuc}>", BF1=bf1, SFO1=bf1, O1=0, DATE=date,
           PROBHD=f"<{probe}>", PULPROG="<zg>")
    if sf is not None:
        _jcamp(e / "pdata" / "1" / "procs", SF=sf, SI=1024, SR=(sf - bf1) * 1e6)
        (e / "pdata" / "1" / "title").write_text(title, encoding="utf-8")
    return e


SF_H = 599.771479328723          # the 2026-05 session's referenced 1H SF
BF1_H = 599.772
BF1_AL, BF1_P = 156.281744, 242.792156


def _session(tmp_path):
    root = tmp_path / "2026-05"
    t0 = 1_777_400_000
    _expno(root, "04272026_RS40339_P5", 1, "1H", BF1_H, SF_H, t0,
           "1H spectrum for later referencing")
    al_ok = R.expected_sf_MHz(SF_H, "27Al")
    _expno(root, "04272026_P1-Bi1-12", 2701, "27Al", BF1_AL, al_ok, t0 + 3600, "27Al ok")
    _expno(root, "04272026_P1-Bi1-12", 2702, "27Al", BF1_AL, BF1_AL, t0 + 7200,
           "27Al forgot xiref")
    _expno(root, "04272026_P1-Bi1-12", 3101, "31P", BF1_P,
           R.expected_sf_MHz(SF_H, "31P") + 2.0e-4, t0 + 9000, "31P stale SR")
    _expno(root, "04272026_P1-Bi1-12", 3102, "31P", BF1_P, None, t0 + 9500, "unprocessed")
    # a second magnet (1.1 GHz): its 35Cl has no 1H reference in this session
    _expno(root, "06052026_LAW", 5, "35Cl", 107.85, 107.8503, t0 + 20000, "35Cl")
    return root


def test_xi_arithmetic_reproduces_topspins_xiref():
    """The 04272026 session: xiref from the 1H gave SF(27Al) 156.281609 and
    SF(31P) 242.791945; the Xi recomputation lands on the same values."""
    assert R.xi_ratio("27Al") == pytest.approx(0.26056859, abs=2e-8)
    assert R.xi_ratio("31P") == pytest.approx(0.40480742, abs=2e-8)
    assert R.expected_sf_MHz(SF_H, "27Al") == pytest.approx(156.281608690994, abs=1e-6)
    assert R.expected_sf_MHz(SF_H, "31P") == pytest.approx(242.791945, abs=1e-6)
    assert R.expected_sr_hz(SF_H, "27Al", BF1_AL) == pytest.approx(-135.31, abs=0.02)
    assert R.expected_sr_hz(SF_H, "31P", BF1_P) == pytest.approx(-210.86, abs=0.2)


def test_axis_shift_follows_topspin_direction():
    """SR from 0 to -135.31 Hz at 156.28 MHz moves the axis by +0.866 ppm (a
    forgotten xiref leaves peaks too LOW in ppm by that amount)."""
    assert R.axis_shift_ppm(0.0, -135.31, 156.281744) == pytest.approx(+0.8658, abs=1e-3)
    assert R.axis_shift_ppm(-135.31, 0.0, 156.281744) == pytest.approx(-0.8658, abs=1e-3)
    assert R.sf_from_peak(SF_H, 1.85, 1.82) == pytest.approx(SF_H * (1 + 1.85e-6) / (1 + 1.82e-6))


def test_scan_audit_list_and_log(tmp_path, monkeypatch):
    root = _session(tmp_path)
    monkeypatch.setenv("LARMOR_REF_LOG", str(tmp_path / "log.jsonl"))
    acqs = R.scan_session(root)
    assert [a.label for a in acqs] == ["04272026_RS40339_P5/1", "04272026_P1-Bi1-12/2701",
                                       "04272026_P1-Bi1-12/2702", "04272026_P1-Bi1-12/3101",
                                       "04272026_P1-Bi1-12/3102", "06052026_LAW/5"]
    refs = R.find_references(acqs)
    assert [r.label for r in refs] == ["04272026_RS40339_P5/1"]
    assert refs[0].sr_hz == pytest.approx((SF_H - BF1_H) * 1e6, abs=0.01)

    rows = {r.acq.label: r for r in R.audit(acqs, tol_ppm=0.05, ada_ppm=ADA)}
    assert rows["04272026_RS40339_P5/1"].verdict == "1H reference"
    ok = rows["04272026_P1-Bi1-12/2701"]
    assert ok.verdict == "ok" and abs(ok.delta_ppm) < 1e-3 and ok.ref.expno == 1
    forgot = rows["04272026_P1-Bi1-12/2702"]
    assert forgot.verdict == "unreferenced" and forgot.flagged
    assert forgot.expected_sr_hz == pytest.approx(-135.31, abs=0.02)
    assert forgot.delta_ppm == pytest.approx(+0.866, abs=2e-3)
    stale = rows["04272026_P1-Bi1-12/3101"]
    assert stale.verdict == "off" and stale.delta_ppm == pytest.approx(0.824, abs=0.01)
    assert rows["04272026_P1-Bi1-12/3102"].verdict == "unprocessed"
    assert rows["06052026_LAW/5"].verdict == "no reference"    # other magnet
    # the partly referenced 27Al pair is noted independently of the 1H
    assert any("27Al" in n for n in ok.notes)
    counts = R.summary(list(rows.values()))
    assert counts == {"1H reference": 1, "ok": 1, "unreferenced": 1, "off": 1,
                      "unprocessed": 1, "no reference": 1}

    txt = R.topspin_list(list(rows.values()), session=str(root))
    assert "04272026_P1-Bi1-12/2702" in txt and "sr    -135.31" in txt
    assert "was 0.00" in txt and "2701" not in txt             # ok rows not listed
    assert "3101" in txt and "sr    -210.86" in txt
    full = R.topspin_list(list(rows.values()), only_flagged=False)
    assert "2701" in full

    csv_path = R.to_csv(list(rows.values()), tmp_path / "audit.csv")
    body = csv_path.read_text(encoding="utf-8")
    assert "unreferenced" in body and "expected_sr_hz" in body

    assert R.previous_audits(root) == []
    log = R.append_log(list(rows.values()), root, ada_ppm=ADA, tol_ppm=0.05, action="export")
    assert log == tmp_path / "log.jsonl"
    rec = json.loads(log.read_text(encoding="utf-8").splitlines()[-1])
    assert rec["action"] == "export" and rec["ada_ppm"] == ADA
    labels = {(e["sample"], e["expno"]) for e in rec["entries"]}
    assert labels == {("04272026_P1-Bi1-12", 2702), ("04272026_P1-Bi1-12", 3101)}
    e2702 = next(e for e in rec["entries"] if e["expno"] == 2702)
    assert e2702["old_sr_hz"] == pytest.approx(0.0, abs=0.01)
    assert e2702["new_sr_hz"] == pytest.approx(-135.31, abs=0.02)
    assert e2702["reference"].endswith("1")
    # append-only: a second run adds a record, the first survives
    R.append_log(list(rows.values()), root, ada_ppm=ADA, tol_ppm=0.05, action="apply")
    assert len(R.previous_audits(root)) == 2
    assert len(log.read_text(encoding="utf-8").splitlines()) == 2


def test_session_root_and_reference_check(tmp_path):
    root = _session(tmp_path)
    e = root / "04272026_P1-Bi1-12" / "2701"
    assert R.session_root(e) == root
    assert R.session_root(e / "pdata" / "1") == root
    assert R.session_root(root / "04272026_P1-Bi1-12") == root
    assert R.session_root(root) == root
    ref = R.find_references(R.scan_session(root))[0]
    x = np.linspace(-10.0, 10.0, 2001)

    def load_ok(_p):
        return x, np.exp(-((x - 1.83) / 0.05) ** 2)

    def load_tms(_p):
        return x, np.exp(-((x - 0.01) / 0.05) ** 2)

    c = R.check_reference(ref, ADA, load=load_ok)
    assert c.ok and c.peak_ppm == pytest.approx(1.83, abs=0.02)
    c2 = R.check_reference(ref, ADA, load=load_tms)
    assert not c2.ok and "0 ppm" in c2.status


def test_real_session_2026_05_is_consistently_referenced():
    """The user's 04272026 session: every 27Al/31P EXPNO of P1-Bi1-12 was
    indirectly referenced from RS40339/1 -- the audit must say so."""
    root = require(DATA_ROOT / "Desktop/WSU_work/NMR/NMRFAM/DATA/2026-05")
    acqs = R.scan_session(root)
    refs = R.find_references(acqs)
    assert any(r.label == "04272026_RS40339_P5-Bi1-12/1" for r in refs)
    rows = R.audit(acqs, tol_ppm=0.05, ada_ppm=ADA)
    sample = [r for r in rows if r.acq.sample == "04272026_P1-Bi1-12_SS_ALP"]
    assert sample
    for r in sample:
        assert r.verdict == "ok", (r.acq.label, r.verdict, r.delta_ppm)
        assert abs(r.delta_ppm) < 0.002
        assert r.ref.label == "04272026_RS40339_P5-Bi1-12/1"


def test_cli_srcheck(tmp_path, monkeypatch, capsys):
    from larmor import cli

    root = _session(tmp_path)
    monkeypatch.setenv("LARMOR_REF_LOG", str(tmp_path / "log.jsonl"))
    out_csv = tmp_path / "out.csv"
    rc = cli.main(["srcheck", str(root / "04272026_P1-Bi1-12" / "2701"),
                   "--csv", str(out_csv)])
    assert rc == 0
    text = capsys.readouterr().out
    assert "2702" in text and "unreferenced" in text and "summary:" in text
    assert out_csv.exists() and (tmp_path / "out_topspin.txt").exists()
    assert (tmp_path / "log.jsonl").exists()


# ---------------------------------------------------------------- Qt
def test_dialog_scans_and_exports(tmp_path, monkeypatch):
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication, QFileDialog
    QApplication.instance() or QApplication([])
    from larmor.desktop.referencing_dialog import ReferencingAuditDialog

    root = _session(tmp_path)
    monkeypatch.setenv("LARMOR_REF_LOG", str(tmp_path / "log.jsonl"))
    x = np.linspace(-10.0, 10.0, 2001)
    monkeypatch.setattr(R, "check_reference",
                        lambda acq, ada, load=None: R.ReferenceCheck(
                            acq, 1.82, "adamantane line at 1.82 ppm — reference confirmed",
                            True, x, np.exp(-((x - 1.82) / 0.05) ** 2)))
    d = ReferencingAuditDialog(None, str(root), str(root / "04272026_P1-Bi1-12" / "2702"))
    try:
        d.scan()
        assert d.refs.rowCount() == 1
        assert d.table.rowCount() == 4          # 2 flagged + no reference + unprocessed
        d.chkAll.setChecked(True); d._fill_rows()
        assert d.table.rowCount() == 6
        assert "unreferenced" in d.windowTitle()
        assert d.btnApply.isEnabled()           # the open spectrum is a flagged row
        got = []
        d.apply_requested.connect(lambda *a: got.append(a))
        d._apply()
        assert got and got[0][1] == pytest.approx(-135.31, abs=0.02) and got[0][2] == 0.0
        out = tmp_path / "audit.csv"
        monkeypatch.setattr(QFileDialog, "getSaveFileName",
                            staticmethod(lambda *a, **k: (str(out), "CSV (*.csv)")))
        d.export()
        assert out.exists() and (tmp_path / "audit_topspin.txt").exists()
        assert len(R.previous_audits(root)) == 2          # apply + export
        d.scan()
        assert "audited before" in d.status.text()
    finally:
        d.close()


def test_app_re_references_the_open_spectrum():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    os.environ["LARMOR_NO_SESSION"] = "1"
    from PySide6.QtWidgets import QApplication, QMenu
    QApplication.instance() or QApplication([])
    from larmor.desktop.app import MainWindow

    win = MainWindow()
    try:
        labels = [a.text() for m in win.menuBar().findChildren(QMenu) for a in m.actions()]
        assert any("Referencing" in t for t in labels)
        x = np.linspace(-50.0, 150.0, 801)
        win._display_1d(x, np.exp(-((x - 60.0) / 5.0) ** 2), "27Al", 156.281744,
                        20000.0, "t", "x")
        win.source_path = r"C:\data\2026-05\S\2701"
        before = win.exp_ppm.copy()
        win._apply_sr_correction(r"C:\data\2026-05\S\2701", -135.31, 0.0, "test")
        assert np.allclose(win.exp_ppm, before + 0.8658, atol=1e-3)   # TopSpin direction
        assert win.recipe["sr_hz"] == pytest.approx(-135.31)
        assert win.recipe["provenance"]["referencing"]["old_sr_hz"] == 0.0
        win.undo()
        assert np.allclose(win.exp_ppm, before)
        # a different EXPNO than the open one is refused
        win._apply_sr_correction(r"C:\data\2026-05\S\2702", -135.31, 0.0, "")
        assert np.allclose(win.exp_ppm, before)
        # the Experiment dialog's SR edit moves the axis the same way
        from PySide6.QtWidgets import QDialog
        orig = QDialog.exec

        def fake_exec(dlg_self):            # the user typed SR = -135.31 and accepted
            win.recipe["sr_hz"] = -135.31
            return True

        QDialog.exec = fake_exec
        try:
            win.recipe["sr_hz"] = 0.0
            win.edit_experiment()
            assert np.allclose(win.exp_ppm, before + 0.8658, atol=1e-3)
        finally:
            QDialog.exec = orig
    finally:
        win.close()
