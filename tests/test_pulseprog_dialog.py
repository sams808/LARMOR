"""Offscreen tests of the pulse program viewer (larmor.desktop.pulseprog_dialog)
on the two shipped EXPNOs, plus the Explorer hook and the Tools slot: the
dialog builds, the diagram paints into a QImage, Export writes PNG and SVG,
the Parameters tab lists p1 with its acqus value, the "?" button opens the
manual, an EXPNO without a pulse program opens calmly."""
import os
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")

pytest.importorskip("PySide6")
from PySide6.QtWidgets import QApplication  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ZG = ROOT / "examples" / "pCABS2-4" / "1118"
MQ = ROOT / "examples" / "pCABS2-4" / "3620"

pytestmark = pytest.mark.skipif(not (ZG / "pulseprogram").exists(),
                                reason="shipped example data missing")


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def _rows(table):
    out = {}
    for r in range(table.rowCount()):
        out[table.item(r, 0).text()] = tuple(table.item(r, c).text() for c in range(1, 4))
    return out


def _distinct_colours(img, step=7):
    seen = set()
    for y in range(0, img.height(), step):
        for x in range(0, img.width(), step):
            seen.add(img.pixel(x, y))
    return seen


def test_dialog_builds_for_zg_with_three_tabs_text_and_parameters(qapp):
    from larmor.desktop.pulseprog_dialog import PulseProgramDialog

    dlg = PulseProgramDialog(None, str(ZG))
    assert "zg" in dlg.windowTitle() and "1118" in dlg.windowTitle()
    assert [dlg.tabs.tabText(i) for i in range(dlg.tabs.count())] == ["Diagram", "Text", "Parameters"]
    assert "zg" in dlg.header.text() and "f1 11B" in dlg.header.text()
    assert dlg.summary.text() == "zg: d1 – p1 (ph1) – acquire go=2 (ph31); phase cycle of 8 steps"
    text = dlg.text.toPlainText()
    assert "  go=2 ph31" in text and ";p1 : f1 channel -  high power pulse" in text
    assert '# 1 "' not in text                      # cpp markers removed, legend kept
    rows = _rows(dlg.table)
    assert rows["p1"] == ("2.1 µs", "µs", "f1 channel -  high power pulse")
    assert rows["d1"] == ("5 s", "s", "relaxation delay; 1-5 * T1")
    assert rows["ns"][0] == "512" and rows["ns"][1] == ""
    assert rows["MCWRK"][0] == "10 ms" and rows["MCWRK"][2].startswith("= 0.333333*30m")
    assert rows["ph31"] == ("0 2 2 0 1 3 3 1", "× 90°", "8 steps")
    assert list(rows)[0] == "p1"                    # pulses first, then delays
    assert dlg.table.item(0, 0).toolTip() == "from acqus"
    tl = dlg.diagram.timeline()
    assert tl is not None and tl.channels == ["f1", "acq"]
    assert dlg.notes.text() == ""
    assert dlg.btnExport.isEnabled()
    dlg.close()


def test_diagram_paints_into_an_image_and_exports_png_and_svg(qapp, tmp_path):
    from PySide6.QtGui import QImage

    from larmor.desktop.pulseprog_dialog import PulseProgramDialog

    dlg = PulseProgramDialog(None, str(ZG))
    w, h = dlg.diagram.sizeHint().width(), dlg.diagram.sizeHint().height()
    assert w > 500 and h > 150
    img = dlg.diagram.image()
    assert not img.isNull() and img.width() == 2 * w and img.height() == 2 * h
    assert len(_distinct_colours(img)) > 3           # ink, accent and background at least
    # painting the live widget (paintEvent) must not raise either
    dlg.diagram.resize(w, h)
    dlg.diagram.grab()
    png = tmp_path / "zg_diagram.png"
    assert dlg.export_to(str(png)) == str(png)
    assert png.stat().st_size > 2000
    back = QImage(str(png))
    assert back.width() == 2 * w and back.height() == 2 * h
    svg = tmp_path / "zg_diagram.svg"
    out = dlg.export_to(str(svg))
    assert Path(out).exists()
    if out.endswith(".svg"):                          # QtSvg present
        head = svg.read_text(encoding="utf-8", errors="replace")[:600]
        assert "<svg" in head
        assert "not to scale" in svg.read_text(encoding="utf-8", errors="replace")
    else:                                             # fallback when QtSvg is missing
        assert out.endswith(".png")
    dlg.close()


def test_mq_dialog_shows_shaped_pulse_loops_and_definition_values(qapp):
    from larmor.desktop.pulseprog_dialog import PulseProgramDialog

    dlg = PulseProgramDialog(None, str(MQ))
    assert "mp3qdfsz" in dlg.windowTitle()
    assert "f1 27Al" in dlg.header.text() and "4 pulse z-filtered 3QMAS with DFS" in dlg.header.text()
    tl = dlg.diagram.timeline()
    pulses = [e for e in tl.events if e.kind in ("pulse", "shaped_pulse")]
    assert len(pulses) == 4 and sum(e.kind == "shaped_pulse" for e in pulses) == 1
    assert len(tl.loops) == 4 and {lp.text for lp in tl.loops} >= {"× ns = 3600", "× 2", "× ST1CNT = 11"}
    rows = _rows(dlg.table)
    assert rows["p2"][0] == "9.615 µs" and rows["p2"][2].startswith("= 1s/(cnst31*cnst0)")
    assert "duration of sweep" in rows["p2"][2]
    assert rows["sp1"] == ("220 W", "W", "-23.4 dB · shape dfs — power level for frequency sweep")
    assert rows["pl21"][0] == "7.5 W" and "-8.8 dB" in rows["pl21"][2]
    assert rows["td1"] == ("22", "", "number of experiments")
    assert rows["ST1CNT"][0] == "11"
    assert rows["ph1"][1] == "× 30°"
    img = dlg.diagram.image(1.0)
    assert len(_distinct_colours(img)) > 3
    dlg.close()


def test_help_button_and_context_menu_open_the_manual(qapp, monkeypatch):
    from larmor.desktop import help_dialog
    from larmor.desktop.pulseprog_dialog import HELP_PAGE, PulseProgramDialog

    calls = []
    monkeypatch.setattr(help_dialog, "show_help",
                        lambda parent, name, title="Help", kind="manual": calls.append((name, title)))
    dlg = PulseProgramDialog(None, str(ZG))
    dlg.btnHelp.click()
    assert calls == [HELP_PAGE] == [("pulse-programs", "Reading a pulse program")]
    assert help_dialog.help_path("pulse-programs") is not None
    dlg._copy_image()
    assert not QApplication.clipboard().image().isNull()
    assert "copied" in dlg.notes.text()
    dlg.close()


def test_expno_without_pulse_program_opens_calmly(qapp, tmp_path):
    from larmor.desktop.pulseprog_dialog import PulseProgramDialog

    e = tmp_path / "12"
    e.mkdir()
    (e / "acqus").write_text("##TITLE= Parameter file\n##$PULPROG= <zg>\n##$NUC1= <7Li>\n##END=\n")
    dlg = PulseProgramDialog(None, str(e))
    assert "zg" in dlg.windowTitle()
    assert "no pulse program file" in dlg.summary.text()
    assert dlg.diagram.timeline() is None
    assert not dlg.btnExport.isEnabled()
    assert dlg.table.rowCount() == 0
    dlg.diagram.grab()                                # the placeholder paints
    dlg.close()
    empty = PulseProgramDialog(None, "")
    assert "no EXPNO" in empty.summary.text()
    empty.close()


def test_explorer_signal_and_dataset_info_line(qapp):
    from PySide6.QtCore import SignalInstance

    from larmor.desktop.explorer import ExplorerPanel

    panel = ExplorerPanel()
    assert isinstance(panel.pulseprog_requested, SignalInstance)
    got = []
    panel.pulseprog_requested.connect(got.append)
    panel.pulseprog_requested.emit(str(ZG))
    assert got == [str(ZG)]
    text = ExplorerPanel.dataset_info_text(str(ZG))
    assert "pulse program: pulseprogram (" in text and "lines)" in text
    assert "Pulse program…" in text
    panel.close()


def test_tools_slot_uses_the_open_expno_or_says_what_to_do(qapp):
    from PySide6.QtWidgets import QStatusBar, QWidget

    from larmor.desktop.mw_tools import _ToolsMixin
    from larmor.desktop.pulseprog_dialog import PulseProgramDialog

    class Host(QWidget):
        def __init__(self):
            super().__init__()
            self.source_path = None
            self._sb = QStatusBar(self)

        def statusBar(self):
            return self._sb

    host = Host()
    _ToolsMixin.open_pulse_program(host, False)              # a QAction's bool
    assert "open a Bruker EXPNO first" in host._sb.currentMessage()
    assert not getattr(host, "_tool_windows", [])
    host.source_path = str(ZG / "pdata" / "1" / "1r")
    _ToolsMixin.open_pulse_program(host, False)
    wins = [w for w in host._tool_windows if isinstance(w, PulseProgramDialog)]
    assert len(wins) == 1 and "zg" in wins[0].windowTitle() and not wins[0].isModal()
    _ToolsMixin.open_pulse_program(host, str(MQ))              # an Explorer path wins
    titles = [w.windowTitle() for w in host._tool_windows if isinstance(w, PulseProgramDialog)]
    assert any("mp3qdfsz" in t for t in titles)
    for w in list(host._tool_windows):
        w.close()
    host.close()
