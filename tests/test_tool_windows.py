"""Read-only viewers and calculators open as NON-modal tool windows.

Item 1 of the 0.14.2 field-use batch (wip/WD): dialogs that only display
information or compute something for the user -- NMR table, Conversion tools,
Variable temperature, Parameter correlations, Czjzek distribution, Compare
fits, chi-square map, Integrals & measurements, the Explorer's dataset info,
About -- used to run through ``exec()`` and froze the workbench behind them.
They now go through ``windowtray.show_tool_window``: a real window with
minimise / maximise / close buttons, kept alive on the main window, freed on
close, listed in the open-windows bar. Dialogs whose caller reads a result
back (``if dlg.exec(): ...``) and the analysis tools that push results into
the workbench are deliberately left modal and are not covered here.
"""
import os
import types

import numpy as np
import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QEvent, Qt  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def win(qapp, monkeypatch):
    monkeypatch.setenv("LARMOR_NO_SESSION", "1")
    from larmor.desktop.app import MainWindow

    w = MainWindow()
    w.show()
    qapp.processEvents()
    yield w
    w._window_bar.close_all()
    w.close()


def _recipe():
    return {"nucleus": "27Al", "larmor_frequency_MHz": 130.3, "spin_rate_Hz": 12000.0,
            "sites": [
                {"model": "gauss_lor", "label": "a", "params": {
                    "isotropic_chemical_shift_ppm": {"value": 10.0},
                    "shift_fwhm_ppm": {"value": 5.0}, "amplitude": {"value": 1.0},
                    "gl": {"value": 0.5}}},
                {"model": "czjzek", "label": "b", "params": {
                    "isotropic_chemical_shift_ppm": {"value": 60.0},
                    "sigma_Cq_MHz": {"value": 1.2}, "amplitude": {"value": 0.7}}}]}


def _spectrum(win):
    # a dense grid: the Integrals dialog seeds its first region from the
    # plot's view range, which offscreen (never painted) is still (0, 1) ppm
    x = np.linspace(200.0, -100.0, 6001)
    y = np.exp(-((x - 10.0) ** 2) / 50.0) + 0.7 * np.exp(-((x - 60.0) ** 2) / 400.0)
    win.exp_ppm, win.exp_amp = x, y
    win.recipe = _recipe()


def _opened(win, qapp, slot, *args, **kw):
    """Run a launcher slot; return the tool window it showed."""
    n = len(getattr(win, "_tool_windows", []))
    slot(*args, **kw)
    qapp.processEvents()
    keep = win._tool_windows
    assert len(keep) == n + 1, "the slot did not open a kept-alive tool window"
    d = keep[-1]
    assert not d.isModal() and d.isVisible()
    assert d.windowFlags() & Qt.Window
    assert d.windowFlags() & Qt.WindowMinMaxButtonsHint
    assert d.testAttribute(Qt.WA_DeleteOnClose)
    assert d in win._window_bar.entries()
    # the main window keeps taking actions while the viewer is open
    before = win.actResid.isChecked()
    win.actResid.trigger()
    assert win.actResid.isChecked() != before
    return d


def _flush(qapp):
    QApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    qapp.processEvents()


def test_nmr_table_and_its_isotope_details(win, qapp):
    from larmor.desktop.utilities import IsotopeDetailsDialog, NmrTableDialog

    d = _opened(win, qapp, win.open_nmr_table)
    assert isinstance(d, NmrTableDialog)
    d._details("Al")
    qapp.processEvents()
    det = d._tool_windows[-1]                                 # owned by the table
    assert isinstance(det, IsotopeDetailsDialog) and not det.isModal()
    assert det.isVisible() and "27Al" in det.windowTitle() or "Al" in det.windowTitle()
    d.close()
    _flush(qapp)
    assert d not in win._tool_windows


def test_conversion_tools_and_variable_temperature(win, qapp):
    from larmor.desktop.utilities import ConvertDialog
    from larmor.desktop.vt_dialog import VtDialog

    assert isinstance(_opened(win, qapp, win.open_convert), ConvertDialog)
    assert isinstance(_opened(win, qapp, win.open_vt), VtDialog)
    assert len(win._tool_windows) == 2 and win._window_bar.count() == 2


def test_parameter_correlations(win, qapp):
    from larmor.desktop.correlation_dialog import CorrelationDialog

    win._last_lmfit = types.SimpleNamespace(
        var_names=["a", "b"], covar=np.array([[1.0, 0.9], [0.9, 1.0]]))
    assert isinstance(_opened(win, qapp, win.show_correlations), CorrelationDialog)


def test_czjzek_distribution(win, qapp):
    from larmor.desktop.czjzek_dist_dialog import CzjzekDistDialog

    _spectrum(win)
    assert isinstance(_opened(win, qapp, win.show_czjzek_dist), CzjzekDistDialog)


def test_compare_fits_diff_table(win, qapp, monkeypatch):
    from larmor.desktop import mw_tools
    from larmor.desktop.diff_dialog import RecipeDiffDialog

    _spectrum(win)
    monkeypatch.setattr(mw_tools.QFileDialog, "getOpenFileName",
                        staticmethod(lambda *a, **k: ("C:/x/ref.json", "")))
    ref = _recipe()
    ref["sites"][0]["params"]["isotropic_chemical_shift_ppm"]["value"] = 12.0
    monkeypatch.setattr(win, "_recipe_model_from", lambda path: ref)
    d = _opened(win, qapp, win.compare_with_saved_fit)
    assert isinstance(d, RecipeDiffDialog)


def test_chi2_map_and_integrals(win, qapp):
    from larmor.desktop.chi2map_dialog import Chi2MapDialog
    from larmor.desktop.integrate_dialog import IntegralsDialog

    _spectrum(win)
    assert isinstance(_opened(win, qapp, win.show_chi2_map), Chi2MapDialog)
    assert isinstance(_opened(win, qapp, win.open_integrals), IntegralsDialog)


def test_explorer_dataset_info_and_about(win, qapp, tmp_path):
    d = _opened(win, qapp, win.explorer._dataset_info, str(tmp_path))
    assert "dataset info" in d.windowTitle()
    about = _opened(win, qapp, win._about)
    assert about.windowTitle() == "About LARMOR"
    # Close through the dialog's own button plumbing frees the window
    about.accept()
    _flush(qapp)
    assert about not in win._tool_windows
    assert win._window_bar.count() == 1


def test_value_returning_dialogs_stay_modal():
    """The guard for the audit: the launchers that read a result back keep
    ``exec()`` -- none of them was converted by mistake."""
    import pathlib

    src = (pathlib.Path(__file__).resolve().parents[1] / "larmor" / "desktop"
           / "mw_tools.py").read_text(encoding="utf-8")
    for kept in ("ComputingParamsDialog(self).exec()", "dlg.exec() and dlg.result",
                 "SatrecDialog(self, expno).exec()", "RedorDialog(self, expno).exec()",
                 "TwoDDialog(self, expno).exec()", "BatchReportDialog(self, start).exec()",
                 "SeqFitDialog(self, paths, model).exec()"):
        assert kept in src, kept
    for gone in ("NmrTableDialog(self, h1).exec()", "ConvertDialog(self, sfo).exec()",
                 "VtDialog(self).exec()", "CorrelationDialog(self, lm).exec()",
                 "CzjzekDistDialog(self, self.recipe).exec()"):
        assert gone not in src, gone
