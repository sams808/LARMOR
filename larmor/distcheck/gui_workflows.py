"""Desktop sweep, part 3: the dialog TOOLS of the desktop -- each one built by
its own ``MainWindow`` launcher, driven through its widgets, its computation
run (small spans, few points), its exports written to the scratch folder,
then closed. Where ``gui_menus`` proves every menu action can be triggered
and ``gui_tools`` the workbench workflows, these stages prove that the
error estimators, the batch tools, the line / literature helpers, the
relaxation and wideline tools and the utility windows actually COMPUTE and
EXPORT in this installation -- the frozen build's failures hide in exactly
those lazy imports (matplotlib backends, scipy sub-packages, pybaselines,
the help renderer).

Each stage is ``(name, fn(say, ctx))`` under the ``gui_common`` plumbing:
one offscreen QApplication, the user's settings put back on exit
(``SettingsGuard``), exceptions raised inside Qt slots trapped and turned
into a stage failure (``ExceptionTrap``), every modal box a workflow
provokes rejected within ~100 ms (``ModalDismisser``). A dialog whose
launcher ``exec()``s it is NOT exec'd here (the dismisser would close it):
its ``exec`` is replaced, for the duration of the launcher call only, by a
function that shows the dialog, drives it the way a user would, and returns
the code the launcher expects (``_captured``). So the launcher's own
construction -- the arguments it passes -- is what is tested, and the code
after the modal loop (``_remember_batch``, ``_magres_add_sites``, the
``exec() and dlg.result`` branches) runs on a dialog that really did its
work. Non-modal tool windows (``windowtray.show_tool_window``) are taken
from the window's ``_tool_windows`` list. File dialogs, QInputDialog and
the figure-export options dialog are replaced by answers for the one call
that needs them and restored at once (``_static``).

Example data: the bundled ``examples/pCABS2-4`` (3616 = 1D 27Al MAS, 1118 =
1D 11B MAS, 3620 = 27Al 3QMAS 2D) and the two shipped recipes; without them
a stage falls back to synthetic spectra and says so. Quick mode keeps the
real data but holds the shipped recipe's shape parameters (a free 12-
parameter Czjzek fit takes ~45 s on a laptop; δiso + amplitude take 0.1 s),
scans 5 profile points, runs the minimum 10 Monte-Carlo trials, 2 auto-fit
restarts and a 7×7 χ² grid.
"""
from __future__ import annotations

import contextlib
import json
import os
import time
from pathlib import Path

import numpy as np

from larmor.distcheck.gui_common import (ExceptionTrap, ModalDismisser, SettingsGuard,
                                         close_tool_windows, ensure_app, pump)

__all__ = ["stages"]


def stages() -> list:
    return [("dialog-error-tools", error_tools),
            ("dialog-batch", batch),
            ("dialog-lines-and-literature", lines_and_literature),
            ("dialog-relaxation-wideline", relaxation_wideline),
            ("dialog-tools-misc", tools_misc)]


# ====================================================================== plumbing
_MISSING = object()


@contextlib.contextmanager
def _patched(target, name, value, *, static=False):
    """``setattr`` for the duration, restoring the raw class-dict entry (the
    descriptor itself, so a PySide static or instance method comes back
    exactly). ``static=True`` wraps a function as a staticmethod -- for
    ``QFileDialog.getOpenFileName``, ``QInputDialog.getDouble`` and the
    like, which the application calls on the class."""
    d = getattr(target, "__dict__", {})
    old = d[name] if name in d else _MISSING
    if static and callable(value) and not isinstance(value, staticmethod):
        value = staticmethod(value)
    setattr(target, name, value)
    try:
        yield
    finally:
        if old is _MISSING:
            try:
                delattr(target, name)
            except AttributeError:
                pass
        else:
            setattr(target, name, old)


def _static(cls, name, fn):
    return _patched(cls, name, fn, static=True)


@contextlib.contextmanager
def _captured(cls, drive=None, code: int = 0):
    """Replace ``cls.exec`` while a launcher runs: the dialog is shown
    non-modal, ``drive(dialog)`` works it, and ``exec`` returns ``code``
    (or what ``drive`` returned when not None). Yields the list of dialogs
    the launcher built, so the stage can assert on them afterwards."""
    got: list = []

    def fake_exec(self):
        from PySide6.QtCore import Qt
        got.append(self)
        try:
            self.setModal(False)
            self.setWindowModality(Qt.NonModal)
        except Exception:                                 # noqa: BLE001
            pass
        self.show()
        pump(60)
        rc = drive(self) if drive is not None else None
        pump(30)
        return code if rc is None else rc

    with _patched(cls, "exec", fake_exec):
        yield got


@contextlib.contextmanager
def _env(**kv):
    old = {k: os.environ.get(k) for k in kv}
    os.environ.update({k: str(v) for k, v in kv.items()})
    try:
        yield
    finally:
        for k, v in old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def _wait(cond, timeout: float = 90.0, what: str = "condition") -> None:
    t0 = time.perf_counter()
    while not cond():
        if time.perf_counter() - t0 > timeout:
            raise AssertionError(f"timed out after {timeout:.0f} s waiting for {what}")
        pump(20)
    pump(20)


def _running(worker) -> bool:
    try:
        return worker is not None and hasattr(worker, "isRunning") and worker.isRunning()
    except RuntimeError:                                  # deleted on the C++ side
        return False


def _fit_idle(win) -> bool:
    return (not _running(getattr(win, "_fit_worker", None))
            and not _running(getattr(win, "_fit2d_worker", None))
            and getattr(win, "_active_fit_worker", None) is None)


def _settle(win, timeout: float = 30.0) -> None:
    """Let the debounced simulation land and its thread finish."""
    t0 = time.perf_counter()
    while time.perf_counter() - t0 < timeout:
        pump(10)
        sw = win._sim_worker
        if not (win._sim_timer.isActive() or win._sim_pending or _running(sw)):
            break


def _status(win) -> str:
    return win.statusBar().currentMessage()


@contextlib.contextmanager
def _desktop(say, ctx):
    """One stage's desktop: the main window with the stubs the UI tests use
    (no recent-files writes, no project-notes box, open mode "replace"),
    shown offscreen, under the three guards. On exit every tool window and
    kept-alive dialog is closed, the workers joined, the window closed, and
    the trapped slot exceptions / Qt criticals fail the stage."""
    ensure_app()
    with SettingsGuard(), ExceptionTrap() as trap, ModalDismisser() as dis:
        from larmor.desktop.app import MainWindow

        win = MainWindow()
        win._add_recent = lambda p: None
        win._report_project_notes = lambda notes: None
        win._confirm_open_mode = lambda: "replace"
        win.show()
        pump(100)
        try:
            yield win, trap, dis
        finally:
            try:
                _teardown(win)
            finally:
                if dis.seen:
                    say("auto-dismissed modal(s):", "; ".join(dis.seen))
                if trap.warnings:
                    say(f"{len(trap.warnings)} Qt warning(s), first: {trap.warnings[0][:160]!r}")
                for e in trap.errors:
                    say("TRAPPED:", e.strip())
        if trap.errors:
            raise AssertionError(
                f"{len(trap.errors)} exception(s) escaped into Qt slots / Qt criticals: "
                + " || ".join(e.strip().splitlines()[-1] for e in trap.errors))


def _teardown(win) -> None:
    from PySide6.QtCore import QEvent
    from PySide6.QtWidgets import QApplication

    try:
        bar = getattr(win, "_window_bar", None)
        if bar is not None and hasattr(bar, "close_all"):
            bar.close_all()
    except Exception:                                     # noqa: BLE001
        pass
    for name in ("_qcpmg_dlg", "_qcpmg_batch_dlg", "_ref_audit", "_inventory_dlg",
                 "_acq_dlg"):
        d = getattr(win, name, None)
        if d is not None:
            try:
                d.close()
            except RuntimeError:
                pass
    try:
        from larmor.desktop import qcpmg_fields_dialog as qfd
        if qfd._shared is not None:
            try:
                qfd._shared.close()
            except RuntimeError:
                pass
            qfd._shared = None
    except Exception:                                     # noqa: BLE001
        pass
    close_tool_windows(keep=[win])
    _settle(win)
    _wait(lambda: _fit_idle(win), 90, "the fit worker to finish")
    win.close()                       # joins the workers, shuts the shared pool
    pump(100)
    left = [n for n in ("_warm_worker", "_sim_worker", "_fit_worker", "_fit2d_worker",
                        "_seq_worker") if _running(getattr(win, n, None))]
    QApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    pump(50)
    if left:
        raise AssertionError(f"QThread(s) still running after the window closed: {left}")


def _tool_window(win, cls):
    """The newest kept-alive tool window of the main window, type-checked."""
    keep = getattr(win, "_tool_windows", None) or []
    if not keep:
        raise AssertionError(f"no tool window opened (expected {cls.__name__})")
    d = keep[-1]
    if not isinstance(d, cls):
        raise AssertionError(f"last tool window is {type(d).__name__}, expected {cls.__name__}")
    if not d.isVisible():
        raise AssertionError(f"{cls.__name__} opened but is not visible")
    return d


def _load(win, path) -> None:
    win.load_source(str(path), keep_fit=False)
    pump(100)
    _settle(win)
    if not (win.exp_ppm.size or win._data2d is not None):
        raise AssertionError(f"nothing loaded from {path}: {_status(win)}")


def _pin_range(win, hi: float, lo: float) -> None:
    """A never-shown offscreen view keeps a (0, 1) range: pin the window the
    fit / error tools read from the view box, as a user zoomed to the data."""
    vb = win.view.getPlotItem().getViewBox()
    vb.setXRange(min(hi, lo), max(hi, lo), padding=0)
    pump(30)


def _examples(ctx) -> dict | None:
    ex = ctx.get("examples")
    if not ex:
        return None
    ex = Path(ex)
    d = {"root": ex, "al": ex / "pCABS2-4" / "3616", "b": ex / "pCABS2-4" / "1118",
         "twod": ex / "pCABS2-4" / "3620", "rec_al": ex / "pCABS2-4_27Al.recipe.json",
         "rec_b": ex / "pCABS2-4_11B.recipe.json"}
    return d if all(p.exists() for p in d.values()) else None


def _scratch(ctx, name: str) -> Path:
    out = Path(ctx["out"]) / name
    out.mkdir(parents=True, exist_ok=True)
    return out


# ================================================================ synthetic data
def _gauss(x, centre, fwhm, height):
    return height * np.exp(-4.0 * np.log(2.0) * ((x - centre) / fwhm) ** 2)


def _write_spectrum(path, x, y, nucleus, larmor_MHz, spin_rate_Hz=0.0, sample="") -> str:
    from larmor.io import spectra

    spectra.write_csv(path, np.asarray(x, float), np.asarray(y, float),
                      {"nucleus": nucleus, "larmor_MHz": larmor_MHz,
                       "spin_rate_Hz": spin_rate_Hz, "sample": sample or Path(path).stem})
    return str(path)


def _fresh_site(model: str, label: str, **values) -> dict:
    """A site dict with the registry defaults, ``values`` overriding."""
    from larmor import models as model_registry

    m = model_registry.get(model)
    params = {p.name: {"value": p.default, "stderr": None, "vary": p.vary,
                       "min": p.min, "max": p.max, "expr": None} for p in m.params}
    for k, v in values.items():
        params[k]["value"] = float(v)
    return {"model": model, "label": label, "params": params}


def _synthetic_27al(win):
    """Two Gauss/Lorentz lines with 1 % noise on the workbench, with the
    two-line model that fits them -- the stand-in when the bundled
    example data is absent."""
    x = np.linspace(-40.0, 110.0, 1200)
    y = _gauss(x, 60.0, 14.0, 100.0) + _gauss(x, 8.0, 12.0, 45.0)
    y = y + np.random.default_rng(1).normal(0.0, 1.0, x.size)
    win._display_1d(x, y, "27Al", 130.3, 12000.0, "synthetic 27Al", "")
    win.recipe["sites"] = [
        _fresh_site("gauss_lor", "AlO4", isotropic_chemical_shift_ppm=58.0,
                    shift_fwhm_ppm=12.0, amplitude=90.0),
        _fresh_site("gauss_lor", "AlO6", isotropic_chemical_shift_ppm=10.0,
                    shift_fwhm_ppm=10.0, amplitude=40.0)]
    for s in win.recipe["sites"]:
        s["params"]["gl"]["vary"] = False
    win.on_structure_changed()
    return (100.0, -30.0)


def _series_csvs(folder: Path):
    """Three 11B gauss_lor spectra at 13 / 15 / 17 ppm with different heights
    (the series-mode UI tests' helper) and the one-line model to batch-fit
    them with."""
    from larmor import engine
    from larmor.recipe import Param, Recipe, SiteModel

    x = np.linspace(-20.0, 60.0, 500)
    paths = []
    for k, (pos, amp) in enumerate(zip((13.0, 15.0, 17.0), (100.0, 50.0, 80.0))):
        tr = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0,
                    sites=[SiteModel(model="gauss_lor", label="A", params={
                        "isotropic_chemical_shift_ppm": Param(pos),
                        "shift_fwhm_ppm": Param(6.0), "amplitude": Param(amp),
                        "gl": Param(1.0, vary=False)})])
        _, m, _ = engine.simulate(tr, exp_ppm=x)
        m = m + np.random.default_rng(k).normal(0.0, 0.4, x.size)
        paths.append(_write_spectrum(folder / f"s{k}.csv", x, m, "11B", 160.0, sample=f"s{k}"))
    model = {"nucleus": "11B", "larmor_frequency_MHz": 160.0, "spin_rate_Hz": 0.0,
             "sites": [{"model": "gauss_lor", "label": "A", "params": {
                 "isotropic_chemical_shift_ppm": {"value": 14.0, "min": 0, "max": 30},
                 "shift_fwhm_ppm": {"value": 5.0, "min": 0.1},
                 "amplitude": {"value": 80.0, "min": 0},
                 "gl": {"value": 1.0, "vary": False}}}]}
    return paths, model


def _manifold_csv(path: Path, rate_hz: float = 20000.0, ratio: float = 0.3, n: int = 2) -> str:
    """An 11B MAS manifold at 160.46 MHz: a centreband at 0 ppm and ±k·νrot
    Gaussians of height ratio**k (the desktop tests' sideband fixture)."""
    lar = 160.46
    x = np.linspace(-400.0, 400.0, 4096)
    y = np.zeros_like(x)
    spacing = rate_hz / lar
    for k in range(-n, n + 1):
        y += ratio ** abs(k) * _gauss(x, k * spacing, 4.0, 1.0)
    return _write_spectrum(path, x, y, "11B", lar, rate_hz, "manifold")


def _magres_text() -> str:
    """A magres 1.0 file with a QE-GIPAW header: six PAIRS of 19F at distinct
    shieldings (pairwise equal to 1e-4 ppm) plus two Na -- the DFT tests'
    multiplicity fixture, inlined."""
    recs = [("Na", 1, (520.0, 521.0, 522.0)), ("Na", 2, (520.0, 521.0, 522.0))]
    idx = 3
    for s in (250.0, 262.0, 275.0, 291.0, 304.0, 318.0):
        for rep in range(2):
            eps = 1e-4 * rep
            recs.append(("F", idx, (s - 12.0 + eps, s + eps, s + 12.0 + eps)))
            idx += 1
    head = ["#$magres-abinitio-v1.0", "[calculation]", "calc_code QE-GIPAW",
            "calc_code_version 7.5", "calc_prefix kog", "calc_xcfunctional PBE",
            "calc_cutoffenergy 60.00 Ry", "calc_cutoffenergy_rho 480.00 Ry",
            "calc_pspot Na.pbe-spn-kjpaw_psl.1.0.0.UPF",
            "calc_pspot F.pbe-n-kjpaw_psl.1.0.0.UPF", "calc_kpoint_mp_grid 4 4 4",
            "[/calculation]", "[atoms]", "units lattice Angstrom", "units atom Angstrom"]
    head += [f"atom {el} {el} {i} 0.0 0.0 0.0" for el, i, _ in recs]
    head += ["[/atoms]", "[magres]", "units ms ppm"]
    body = []
    for el, i, (a, b, c) in recs:
        vals = [a, 0.0, 0.0, 0.0, b, 0.0, 0.0, 0.0, c]
        body.append(f"ms {el} {i} " + " ".join(f"{v:.6f}" for v in vals))
    return "\n".join(head + body + ["[/magres]", ""])


def _redor_expno(folder: Path, d_hz: float = 350.0, masr: float = 10000.0) -> Path:
    """A TopSpin-style ``pdata/1/redor.txt`` of an isolated 13C–15N pair with
    dipolar coupling ``d_hz`` (larmor.redor's own forward curve)."""
    from larmor import redor

    n = np.arange(2, 30, 2)
    ds = redor.redor_pair_curve(d_hz, n / masr)
    s0 = -1.0e10
    lines = ["Dataset :", str(folder), "", f"Spinning speed : {masr:.6f}", "",
             "Peak 1 (xy -10.000 ppm)", "",
             "    number integral(S0) integral(S*) intensity(S0) intensity(S*)", ""]
    for k, v in zip(n, ds):
        lines.append(f"{int(k):10d} {s0:14.4e} {s0 * (1.0 - v):14.4e} 0 0")
    p = folder / "pdata" / "1"
    p.mkdir(parents=True, exist_ok=True)
    (p / "redor.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return folder


def _hold_shapes(recipe: dict, keep=("isotropic_chemical_shift_ppm", "amplitude")) -> int:
    """Hold every parameter but ``keep`` (vary=False); returns how many."""
    n = 0
    for s in recipe.get("sites", []):
        for name, p in s.get("params", {}).items():
            if name not in keep and p.get("vary", True) and not p.get("expr"):
                p["vary"] = False
                n += 1
    return n


def _rows(table) -> int:
    return int(table.rowCount())


def _assert_file(p: Path, what: str, min_bytes: int = 20) -> None:
    if not Path(p).is_file():
        raise AssertionError(f"{what}: {p} was not written")
    size = Path(p).stat().st_size
    if size < min_bytes:
        raise AssertionError(f"{what}: {p} is {size} bytes")


def _export_opts(parent, formats, **kw):
    """The answer to the shared figure-export options dialog."""
    fmt = "PNG" if "PNG" in formats else formats[0]
    return {"format": fmt, "dpi": 80, "width_cm": 10.0, "height_cm": 8.0}


# ============================================================ 1. error tools
def error_tools(say, ctx):
    """The example 27Al 1r with the shipped recipe fitted, then every
    estimator and viewer of the Fit menu: χ² profile (Errors analysis),
    Monte-Carlo, χ² map with its figure export, Parameter correlations,
    Czjzek distribution, Compare with a saved fit, Integrals, Auto fit."""
    from PySide6.QtWidgets import QFileDialog, QInputDialog

    quick = bool(ctx.get("quick"))
    out = _scratch(ctx, "error_tools")
    ex = _examples(ctx)
    with _desktop(say, ctx) as (win, trap, dis):
        if ex:
            _load(win, ex["al"])
            win.apply_recipe(str(ex["rec_al"]))
            _settle(win)
            n_sites = len(win.recipe["sites"])
            if n_sites < 3:
                raise AssertionError(f"apply_recipe left {n_sites} site(s): {_status(win)}")
            say(f"example 27Al 1r: {win.exp_ppm.size} points; shipped recipe applied: "
                f"{n_sites} Czjzek sites")
            window = tuple(json.loads(ex["rec_al"].read_text(encoding="utf-8"))
                           .get("fit_window_ppm") or (150.0, -80.0))
            if quick:
                held = _hold_shapes(win.recipe)
                win.lines_table.rebuild(win.recipe, win.hidden)
                say(f"quick: {held} shape parameter(s) held (σ_Cq, widths); δiso and "
                    "amplitude of each site free")
        else:
            window = _synthetic_27al(win)
            _settle(win)
            say("no bundled example data: synthetic two-line 27Al spectrum")
        _pin_range(win, *window)

        # -- the fit through the window's FitWorker
        t0 = time.perf_counter()
        win.run_fit()
        if win._fit_worker is None:
            raise AssertionError(f"run_fit started no worker: {_status(win)}")
        _wait(lambda: _fit_idle(win), 180, "the 1D fit")
        _settle(win)
        if win._last_lmfit is None:
            raise AssertionError(f"the fit produced no lmfit result: {_status(win)}")
        say(f"fit done in {time.perf_counter() - t0:.1f} s · {_status(win)}")

        # -- Errors analysis: the χ² profile of one parameter
        from larmor.desktop.tool_dialogs import ErrorsDialog

        def drive_errors(d):
            d.points.setValue(5)
            d.span.setValue(2.0)
            if d.param.count() == 0:
                raise AssertionError("the Errors dialog lists no free parameter")
            t1 = time.perf_counter()
            d._run()
            if d.prof is None:
                raise AssertionError(f"χ² profile failed: {d.res.text()}")
            if d.prog.value() != 5:
                raise AssertionError(f"progress bar ended at {d.prog.value()} of 5")
            say(f"χ² profile of {d.param.currentText()} ({d.points.value()} points, "
                f"{time.perf_counter() - t1:.1f} s): {d.res.text()}")
            if d.note.text():
                say("  note:", d.note.text())

        with _captured(ErrorsDialog, drive_errors) as got:
            win.run_errors_analysis()
        if not got:
            raise AssertionError("run_errors_analysis built no ErrorsDialog")

        # -- Monte-Carlo (the spinbox floor is 10 trials)
        from larmor.desktop.montecarlo_dialog import MonteCarloDialog

        def drive_mc(d):
            d.n.setValue(10)
            d.seed.setValue(1)
            t1 = time.perf_counter()
            d._run()
            if d._result is None:
                raise AssertionError(f"Monte-Carlo failed: {d.status.text()}")
            if _rows(d.table) < 1:
                raise AssertionError("Monte-Carlo table is empty")
            d._draw_hist()
            if len(d.plot.plotItem.items) < 2:
                raise AssertionError("Monte-Carlo histogram drew nothing")
            d._copy()
            say(f"Monte-Carlo {d.n.value()} trials ({time.perf_counter() - t1:.1f} s): "
                f"{_rows(d.table)} table rows · {d._result.summary}")

        with _captured(MonteCarloDialog, drive_mc) as got:
            win.run_monte_carlo()
        if not got:
            raise AssertionError("run_monte_carlo built no MonteCarloDialog")

        # -- χ² map (tool window) on a 7×7 grid, exported as PNG
        from larmor.desktop import export_dialog
        from larmor.desktop.chi2map_dialog import Chi2MapDialog

        win.show_chi2_map()
        pump(50)
        d = _tool_window(win, Chi2MapDialog)
        d.n.setValue(7)
        t1 = time.perf_counter()
        d._compute()
        if not d.btnExp.isEnabled() or d.img.pixmap() is None or d.img.pixmap().isNull():
            raise AssertionError(f"χ² map computed no image: {d.img.text()}")
        target = out / "chi2_map.png"
        with _patched(export_dialog, "choose", _export_opts), \
                _static(QFileDialog, "getSaveFileName", lambda *a, **k: (str(target), "")):
            d._export()
        _assert_file(target, "χ² map export", 1000)
        say(f"χ² map {d.cbA.currentText()} × {d.cbB.currentText()} 7×7 "
            f"({time.perf_counter() - t1:.1f} s) exported: {target.name}")
        d.close()

        # -- Parameter correlations, Czjzek distribution
        from PySide6.QtWidgets import QTableWidget

        from larmor.desktop.correlation_dialog import CorrelationDialog
        from larmor.desktop.czjzek_dist_dialog import CzjzekDistDialog

        win.show_correlations()
        pump(50)
        d = _tool_window(win, CorrelationDialog)
        tables = d.findChildren(QTableWidget)
        if not tables or _rows(tables[0]) < 1:
            raise AssertionError("the correlation matrix is empty")
        say(f"parameter correlations: {_rows(tables[0])}×{tables[0].columnCount()} matrix")
        d.close()
        if any(s.get("model", "").startswith(("czjzek", "ext_czjzek"))
               for s in win.recipe["sites"]):
            import pyqtgraph as pg

            win.show_czjzek_dist()
            pump(50)
            d = _tool_window(win, CzjzekDistDialog)
            plots = d.findChildren(pg.PlotWidget)
            drawn = sum(len(p.plotItem.items) for p in plots)
            if not plots or drawn == 0:
                raise AssertionError("the Czjzek distribution dialog drew nothing")
            say(f"Czjzek distribution: {len(plots)} plot(s), {drawn} item(s)")
            d.close()
        else:
            win.show_czjzek_dist()
            say("no Czjzek site in the synthetic model:", _status(win))

        # -- Compare with a saved fit (the shipped recipe, or the model itself)
        from larmor.desktop.diff_dialog import RecipeDiffDialog

        if ex:
            ref_path = str(ex["rec_al"])
        else:
            ref_path = str(out / "reference.recipe.json")
            Path(ref_path).write_text(json.dumps(win.recipe), encoding="utf-8")
        with _static(QFileDialog, "getOpenFileName", lambda *a, **k: (ref_path, "")):
            win.compare_with_saved_fit()
        pump(50)
        d = _tool_window(win, RecipeDiffDialog)
        tables = d.findChildren(QTableWidget)
        if not tables or _rows(tables[0]) < 1:
            raise AssertionError("the fit comparison table is empty")
        say(f"compare with {Path(ref_path).name}: {_rows(tables[0])} parameter rows")
        d.close()

        # -- Integrals & measurements
        from larmor.desktop.integrate_dialog import IntegralsDialog

        win.open_integrals()
        pump(50)
        d = _tool_window(win, IntegralsDialog)
        d._add_region()
        pump(30)
        if _rows(d.table) < 1:
            raise AssertionError("the Integrals dialog has no region row")
        d._csv()
        say(f"integrals: {_rows(d.table)} region(s), CSV copied")
        d.close()

        # -- Auto fit (synchronous in the slot; 2 restarts)
        t1 = time.perf_counter()
        with _static(QInputDialog, "getInt", lambda *a, **k: (2, True)):
            win.run_auto_fit()
        _settle(win)
        msg = _status(win)
        if "auto" not in msg.lower() and "rmsd" not in msg.lower():
            raise AssertionError(f"auto fit left no verdict: {msg!r}")
        say(f"auto fit, 2 restarts ({time.perf_counter() - t1:.1f} s): {msg}")
        close_tool_windows(keep=[win])


# ================================================================== 2. batch
def batch(say, ctx):
    """Series ▸ Batch fit on three synthetic 11B spectra through the window's
    launcher: the dialog filled, the fit run in its worker thread, then its
    outputs -- error CSV, Save table, individual fits, the Series plot with
    its CSV and figure exports and the Species bar (Plotting studio), the
    publication bundle -- and Series ▸ Batch fit report on the saved fits."""
    from PySide6.QtWidgets import QFileDialog

    out = _scratch(ctx, "batch")
    paths, model = _series_csvs(out)
    say("three synthetic 11B gauss_lor spectra (13 / 15 / 17 ppm); the two example "
        "1D spectra are different nuclei and cannot share one model")
    with _desktop(say, ctx) as (win, trap, dis):
        _load(win, paths[0])
        win.recipe["sites"] = json.loads(json.dumps(model["sites"]))
        win.on_structure_changed()
        _settle(win)

        from larmor.desktop import export_dialog
        from larmor.desktop.batchfit_dialog import BatchFitDialog
        from larmor.desktop.plotting_studio import PlottingStudio
        from larmor.desktop.series_plot import SeriesPlotDialog

        def drive_batch(d):
            if len(d._data) != 3:
                raise AssertionError(f"the batch dialog loaded {len(d._data)} of 3 spectra")
            if not d.btnFit.isEnabled():
                raise AssertionError(f"Fit is disabled: {d.status.text()}")
            chk = d._rel_checks.get("isotropic_chemical_shift_ppm")
            if chk is not None:
                chk.setChecked(True)                   # the shift moves along the series
            t0 = time.perf_counter()
            d._run()
            if d._worker is None:
                raise AssertionError(f"the batch Fit started no worker: {d.status.text()}")
            _wait(lambda: not _running(d._worker), 90, "the batch worker")
            pump(100)
            if d._result is None:
                raise AssertionError(f"batch fit gave no result: {d.status.text()}")
            say(f"batch fit ({time.perf_counter() - t0:.1f} s): {d.status.text()}")
            for name in ("btnSave", "btnTable", "btnBundle", "btnSeries", "btnErr", "btnErrCsv"):
                if not getattr(d, name).isEnabled():
                    raise AssertionError(f"{name} stayed disabled after the fit")
            # the covariance error CSV
            d.errCombo.setCurrentIndex(0)
            errs = out / "batch_errors.csv"
            d._write_err_csv(str(errs))
            _assert_file(errs, "batch error CSV")
            # Save table…
            table = out / "batch_table.csv"
            with _static(QFileDialog, "getSaveFileName", lambda *a, **k: (str(table), "CSV (*.csv)")):
                d._save_table()
            _assert_file(table, "batch table")
            # the individual fits (what the Batch fit report reads)
            fits = out / "fits"
            fits.mkdir(exist_ok=True)
            n = d._save_all_recipes_to(fits)
            if n != 3 or len(list(fits.glob("*.recipe.json"))) != 3:
                raise AssertionError(f"{n} individual fit(s) saved, 3 expected")
            # Series plot… (a tool window owned by the batch dialog)
            d._series_plot()
            pump(100)
            sp = (getattr(d, "_tool_windows", None) or [None])[-1]
            if not isinstance(sp, SeriesPlotDialog):
                raise AssertionError("Series plot… opened no SeriesPlotDialog")
            sp.list.selectAll()
            pump(50)
            if not sp._selected():
                raise AssertionError("the Series plot lists no component")
            series_csv = out / "series.csv"
            with _static(QFileDialog, "getSaveFileName", lambda *a, **k: (str(series_csv), "")):
                sp._export_csv()
            _assert_file(series_csv, "series CSV")
            series_png = out / "series.png"
            with _patched(export_dialog, "choose", _export_opts), \
                    _static(QFileDialog, "getSaveFileName", lambda *a, **k: (str(series_png), "")):
                sp._export_fig()
            _assert_file(series_png, "series figure", 1000)
            with _captured(PlottingStudio) as studios:
                sp._species_bar()
            if not studios:
                raise AssertionError("Species bar… opened no Plotting studio")
            spec = studios[0]._spec()
            if spec.get("kind") != "species_bar":
                raise AssertionError(f"the Species bar spec is {spec.get('kind')!r}")
            say(f"series plot: {len(sp._selected())} component(s) · {series_csv.name}, "
                f"{series_png.name} · species bar with {len(spec.get('categories', []))} categories")
            sp.close()
            # Publication bundle…
            bundle = out / "bundle"
            bundle.mkdir(exist_ok=True)
            with _static(QFileDialog, "getExistingDirectory", lambda *a, **k: str(bundle)):
                d._export_bundle()
            names = {p.name for p in bundle.iterdir()}
            need = {"manifest.csv", "README.txt", "batch_table.csv"}
            if not need <= names:
                raise AssertionError(f"bundle lacks {need - names}: {d.status.text()}")
            if len(list(bundle.glob("*.recipe.json"))) != 3 or len(list(bundle.glob("*_curves.csv"))) != 3:
                raise AssertionError(f"bundle holds {sorted(names)}")
            say(f"publication bundle: {len(names)} files · {d.status.text()}")

        with _captured(BatchFitDialog, drive_batch) as got:
            win.run_batch_fit(paths)
        if not got:
            raise AssertionError("run_batch_fit built no BatchFitDialog")
        if not any(ws["kind"] == "batch" for ws in win.workspaces):
            raise AssertionError("the batch session was not kept in the Workspaces dock")
        say("batch session kept in the Workspaces dock:", _status(win))

        # -- Series ▸ Batch fit report… on the saved individual fits
        from larmor.desktop.batch_dialog import BatchReportDialog

        recs = sorted(str(p) for p in (out / "fits").glob("*.recipe.json"))
        report = out / "report"

        def drive_report(d):
            with _static(QFileDialog, "getOpenFileNames", lambda *a, **k: (recs, "")):
                d._add()
            if d.list.count() != 3:
                raise AssertionError(f"the report lists {d.list.count()} fits")
            d.out.setText(str(report))
            d.method.setCurrentIndex(0)                       # covariance
            t0 = time.perf_counter()
            d._generate()
            if d.status.text().startswith("failed") or not d.btnOpen.isEnabled():
                raise AssertionError(f"batch report: {d.status.text()}")
            files = sorted(p.relative_to(report).as_posix() for p in report.rglob("*") if p.is_file())
            if not files:
                raise AssertionError("the batch report wrote nothing")
            say(f"batch report ({time.perf_counter() - t0:.1f} s): {d.status.text()} · "
                f"{len(files)} files: {', '.join(files[:8])}{' …' if len(files) > 8 else ''}")

        with _captured(BatchReportDialog, drive_report) as got:
            win.run_batch_report()
        if not got:
            raise AssertionError("run_batch_report built no BatchReportDialog")
        close_tool_windows(keep=[win])


# ==================================================== 3. lines and literature
def lines_and_literature(say, ctx):
    """The Edit-menu helpers that add or constrain lines through small
    dialogs: peak picking, literature labels, the field prediction, the
    function line, a measured background spectrum as a component, the
    glass protocol, saved constraint sets, the fit threshold, the MQMAS F1
    reference, the computing parameters; then the sideband detector's
    offer and the Add spinning sidebands dialog on a synthetic manifold."""
    from larmor.desktop.prefs import settings as prefs_settings
    from PySide6.QtWidgets import QDialog, QFileDialog, QInputDialog

    out = _scratch(ctx, "lines")
    ex = _examples(ctx)
    with _desktop(say, ctx) as (win, trap, dis):
        if ex:
            _load(win, ex["al"])
            win.new_fit()
            _settle(win)
            say(f"example 27Al 1r: {win.exp_ppm.size} points, empty model")
            bg_path = str(ex["b"] / "pdata" / "1" / "1r")
        else:
            _synthetic_27al(win)
            win.new_fit()
            _settle(win)
            x = np.linspace(-40.0, 110.0, 800)
            bg_path = _write_spectrum(out / "background.csv", x, _gauss(x, -10.0, 30.0, 1.0),
                                      "27Al", 130.3)
            say("no bundled example data: synthetic 27Al spectrum")
        _pin_range(win, float(win.exp_ppm.max()), float(win.exp_ppm.min()))

        # -- Add a line at every peak
        with _static(QInputDialog, "getDouble", lambda *a, **k: (5.0, True)):
            win.autopick_lines()
        _settle(win)
        n_pk = len(win.recipe["sites"])
        if n_pk < 1:
            raise AssertionError(f"autopick added no line: {_status(win)}")
        say(f"autopick: {_status(win)}")

        # -- Label lines from literature ranges (27Al: AlO4 / AlO5 / AlO6 bands)
        win.label_from_literature()
        msg = _status(win)
        if "no literature ranges" in msg:
            raise AssertionError(msg)
        say(f"literature labels: {msg}")

        # -- Predict at another field: a NEW workspace with the simulated model
        home, n_ws = win.active_ws, len(win.workspaces)
        with _static(QInputDialog, "getDouble", lambda *a, **k: (600.0, True)):
            win.predict_at_field()
        _settle(win)
        if len(win.workspaces) != n_ws + 1 or "predicted" not in _status(win):
            raise AssertionError(f"predict_at_field: {_status(win)} ({len(win.workspaces)} workspaces)")
        say(f"predict at field: {_status(win)}")
        win.switch_workspace(home)
        _settle(win)
        if len(win.recipe["sites"]) != n_pk:
            raise AssertionError("switching back lost the picked lines")

        # -- Add function line…
        expr = "a * exp(-((x - b) / c)**2) + d"
        with _static(QInputDialog, "getText", lambda *a, **k: (expr, True)):
            win.add_function_line()
        _settle(win)
        if win.recipe["sites"][-1].get("model") != "function":
            raise AssertionError(f"no function line added: {_status(win)}")
        say(f"function line: {_status(win)}")

        # -- Add background spectrum… (a measured trace as a fit component)
        with _static(QFileDialog, "getOpenFileName", lambda *a, **k: (bg_path, "")):
            win.add_background_spectrum()
        _settle(win)
        last = win.recipe["sites"][-1]
        if last.get("model") != "spectrum" or len(last.get("ref", {}).get("ppm", [])) < 100:
            raise AssertionError(f"no background spectrum component: {_status(win)}")
        say(f"background spectrum: {_status(win)} ({len(last['ref']['ppm'])} points)")

        # -- Constraints: restrict around current values, save / apply a set
        with _static(QInputDialog, "getDouble", lambda *a, **k: (3.0, True)):
            win.restrict_glass_protocol()
        if "restricted" not in _status(win):
            raise AssertionError(f"restrict_glass_protocol: {_status(win)}")
        say(f"glass protocol: {_status(win)}")
        name = "distcheck set"
        with _static(QInputDialog, "getText", lambda *a, **k: (name, True)):
            win.save_constraint_set()
        if "saved constraint set" not in _status(win):
            raise AssertionError(f"save_constraint_set: {_status(win)}")
        lib = json.loads(prefs_settings().value("constraintLibrary", "{}") or "{}")
        if name not in lib:
            raise AssertionError("the constraint set did not reach the library")
        (out / "constraint_set.json").write_text(json.dumps(lib[name], indent=1), encoding="utf-8")

        def pick_item(parent, title, label, items, *a, **k):
            return next(i for i in items if i.startswith(name)), True

        with _static(QInputDialog, "getItem", pick_item):
            win.apply_constraint_set()
        _settle(win)
        if "applied" not in _status(win):
            raise AssertionError(f"apply_constraint_set: {_status(win)}")
        say(f"constraint set (library in QSettings, copy in {out.name}/constraint_set.json): "
            f"{_status(win)}")

        # -- Fit settings: completion threshold, MQMAS F1 reference, computing parameters
        with _static(QInputDialog, "getDouble", lambda *a, **k: (0.1, True)):
            win.edit_fit_tol()
        tol = prefs_settings().value("fitStdevPct")
        if tol is None or abs(float(tol) - 0.1) > 1e-9:
            raise AssertionError(f"fit threshold not stored: {tol!r} · {_status(win)}")
        say(f"fit threshold: {_status(win)}")
        with _static(QInputDialog, "getDouble", lambda *a, **k: (0.5, True)):
            win.edit_mqmas_f1_ref()
        if win.recipe.get("mqmas_f1_ref_ppm") != 0.5 or win.recipe.get("mqmas_f1_ref_vary"):
            raise AssertionError(f"MQMAS F1 reference not set: {_status(win)}")
        say(f"MQMAS F1 reference: {_status(win)}")

        from larmor import engine, twod
        from larmor.desktop.dialogs import ComputingParamsDialog

        kernel_before = dict(engine.KERNEL_SETTINGS)
        mqmas_before = dict(twod.MQMAS_SETTINGS)
        try:
            with _captured(ComputingParamsDialog, lambda d: (d._accept(), 1)[1]) as got:
                win.edit_computing_params()
            if not got or "computing parameters" not in _status(win):
                raise AssertionError(f"edit_computing_params: {_status(win)}")
            say(f"computing parameters: {_status(win)}")
        finally:
            engine.KERNEL_SETTINGS.clear()
            engine.KERNEL_SETTINGS.update(kernel_before)
            twod.MQMAS_SETTINGS.clear()
            twod.MQMAS_SETTINGS.update(mqmas_before)
        _settle(win)

        # -- Spinning sidebands on a synthetic MAS manifold: detect, offer, add
        _load(win, _manifold_csv(out / "manifold.csv"))
        win.detect_sidebands()
        pump(50)
        banner = getattr(win, "ssb_banner", None)
        if banner is None or banner.isHidden():
            raise AssertionError(f"detect_sidebands offered nothing: {_status(win)}")
        say(f"sideband detection: {_status(win)}")
        win._apply_sideband_choice("linked")
        _settle(win)
        n_lines = len(win.recipe["sites"])
        if n_lines < 3 or not any("sideband" in s for s in win.recipe["sites"]):
            raise AssertionError(f"the linked manifold was not added: {_status(win)}")
        say(f"linked manifold: {_status(win)}")
        with _patched(QDialog, "exec", lambda self: QDialog.Accepted):
            win.add_sidebands_for_line(0)
        _settle(win)
        if len(win.recipe["sites"]) <= n_lines:
            raise AssertionError(f"Add spinning sidebands… added nothing: {_status(win)}")
        say(f"add sidebands dialog: {_status(win)}")
        win.undo()
        _settle(win)
        close_tool_windows(keep=[win])


# ============================================ 4. relaxation and wideline tools
def relaxation_wideline(say, ctx):
    """Tools ▸ Relaxation and the wideline tools: guided T1 (synthetic
    build-up through the dialog's slice API), REDOR (a synthetic
    ``redor.txt``), Variable temperature, VOCS stitching, per-site
    relaxation (launcher), QCPMG processing on a synthetic echo train with
    its hand-offs to the workbench and to the infinite-field dialog, the
    batch infinite-field grid, Read static pattern, Herzfeld–Berger, WURST
    correction and Subtract a spectrum -- each on small synthetic data,
    exports to the scratch folder."""
    from PySide6.QtWidgets import QDialog, QFileDialog

    quick = bool(ctx.get("quick"))
    out = _scratch(ctx, "relaxation")
    ex = _examples(ctx)
    with _desktop(say, ctx) as (win, trap, dis):
        if ex:
            _load(win, ex["al"])
        else:
            _synthetic_27al(win)
            _settle(win)
            say("no bundled example data: synthetic 27Al spectrum")

        # -- Guided relaxation (T1): no ser in the example -> synthetic build-up
        from larmor.desktop.satrec_dialog import SatrecDialog

        def drive_satrec(d):
            x = np.linspace(-60.0, 60.0, 600)
            t = np.array([0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0])
            T1 = 1.5
            build = 1.0 - np.exp(-t / T1)
            rng = np.random.default_rng(3)
            slices = (build[:, None] * _gauss(x, 10.0, 6.0, 100.0)[None, :]
                      + rng.normal(0.0, 0.3, (t.size, x.size)))
            d.x_ppm, d._all_delays, d._all_slices = x, t, slices
            d._src = "synthetic vdlist"
            for sb in (d.firstSlice, d.lastSlice):
                sb.blockSignals(True)
                sb.setRange(1, t.size)
            d.firstSlice.setValue(1)
            d.lastSlice.setValue(t.size)
            for sb in (d.firstSlice, d.lastSlice):
                sb.blockSignals(False)
            d._apply_slice_range()
            if not d.zones:
                raise AssertionError("no integration zone was placed on the relaxed slice")
            d.kind.setCurrentText("satrec")
            d._fit()
            if not d.results:
                raise AssertionError(f"T1 fit failed: {d.res.text()}")
            tau = float(d.results[0]["tau"])
            if abs(tau - T1) / T1 > 0.25:
                raise AssertionError(f"T1 {tau:.3g} s off the synthetic {T1} s")
            d._csv()
            say(f"guided T1 on a synthetic build-up: T1 = {tau:.3g} s (truth {T1}), CSV copied")

        with _captured(SatrecDialog, drive_satrec) as got:
            win.open_satrec()
        if not got:
            raise AssertionError("open_satrec built no SatrecDialog")

        # -- REDOR on a synthetic TopSpin redor.txt
        from larmor.desktop.tool_dialogs import RedorDialog

        redor_dir = _redor_expno(out / "redor_expno")

        def drive_redor(d):
            d.expno = str(redor_dir)
            d.lbl.setText(str(redor_dir))
            d.iso1.setText("13C")
            d.iso2.setText("15N")
            d.regime.setCurrentText("pair")
            d._run()
            if d.res.text().startswith("failed") or not d.res.text():
                raise AssertionError(f"REDOR: {d.res.text()!r}")
            if len(d.plot.plotItem.items) < 2:
                raise AssertionError("REDOR drew no points / curve")
            say(f"REDOR: {d.res.text()[:160]}")

        with _captured(RedorDialog, drive_redor) as got:
            win.open_redor()
        if not got:
            raise AssertionError("open_redor built no RedorDialog")

        # -- Variable temperature (Arrhenius then VFT)
        import re

        from larmor.desktop.vt_dialog import VtDialog

        win.open_vt()
        pump(50)
        d = _tool_window(win, VtDialog)
        d.table.setRowCount(0)
        for T in (300.0, 320.0, 340.0, 360.0, 380.0):
            d._add_row(T, f"{1e8 * np.exp(-50e3 / (8.314 * T)):.6g}")
        d.mode.setCurrentText("Arrhenius")
        d._fit()
        m = re.search(r"Ea = ([\d.]+)", d.res.text())
        if not m or abs(float(m.group(1)) - 50.0) > 1.0:
            raise AssertionError(f"Arrhenius fit: {d.res.text()}")
        d.mode.setCurrentText("VFT")
        d._fit()
        if "T0" not in d.res.text():
            raise AssertionError(f"VFT fit: {d.res.text()}")
        say(f"variable temperature: Ea = {m.group(1)} kJ/mol (truth 50) · VFT: {d.res.text()}")
        d.close()

        # -- VOCS: three offset sub-spectra stitched onto the workbench
        from larmor.desktop.vocs_dialog import VocsDialog

        pieces = []
        for i, centre in enumerate((-400.0, 0.0, 400.0)):
            x = np.linspace(centre - 320.0, centre + 320.0, 641)
            y = _gauss(x, -300.0, 150.0, 1.0) + _gauss(x, 250.0, 120.0, 0.7)
            pieces.append(_write_spectrum(out / f"vocs{i}.csv", x, y, "81Br", 216.0))

        def drive_vocs(d):
            with _static(QFileDialog, "getOpenFileNames", lambda *a, **k: (pieces, "")):
                d._add()
            if not d.btnApply.isEnabled():
                raise AssertionError("VOCS: → Workbench stayed disabled after adding 3 pieces")
            d._apply()

        with _captured(VocsDialog, drive_vocs) as got:
            win.open_vocs()
        _settle(win)
        prov = (win.recipe or {}).get("provenance") or {}
        if not got or len(prov.get("vocs_sources", [])) != 3 or win.exp_ppm.size < 100:
            raise AssertionError(f"VOCS stitch did not reach the workbench: {_status(win)}")
        say(f"VOCS: {win.exp_ppm.size} points from 3 pieces · {_status(win)}")

        # -- Per-site relaxation: the launcher up to its EXPNO prompt (needs a ser)
        win.recipe["sites"] = [_fresh_site("gauss_lor", "A", isotropic_chemical_shift_ppm=-300.0,
                                           shift_fwhm_ppm=150.0, amplitude=1.0)]
        win.on_structure_changed()
        _settle(win)
        with _static(QFileDialog, "getExistingDirectory", lambda *a, **k: ""):
            win.open_per_site_relaxation()
        say("per-site relaxation: launcher reached its EXPNO prompt (a relaxation ser is "
            "needed to run the decomposition; none is bundled)")

        # -- QCPMG: the kept-alive dialog, a synthetic echo train, both hand-offs
        from larmor.desktop import qcpmg_fields_dialog as qfd
        from larmor.desktop.qcpmg_dialog import QcpmgDialog

        if ex:
            _load(win, ex["al"])                 # its fid is the dialog's source
        qfd._shared = None
        win.open_qcpmg()
        pump(100)
        d = getattr(win, "_qcpmg_dlg", None)
        if not isinstance(d, QcpmgDialog) or not d.isVisible():
            raise AssertionError("open_qcpmg kept no visible QcpmgDialog")
        say(f"QCPMG opened on {Path(d.source).name if d.source else 'no source'}: "
            f"{d.res.text()[:120]}")
        period, top, n_ech, sw = 128, 64, 24, 50000.0
        tt = np.arange(period) - top
        echo = np.exp(-np.abs(tt) / 9.0) * np.exp(2j * np.pi * 0.06 * tt)
        fid = np.concatenate([echo * np.exp(-k / 6.0) for k in range(n_ech)]).astype(complex)
        d._loading = True
        d.fid = fid
        d.meta = {"sw_Hz": sw, "larmor_MHz": 78.0, "nucleus": "35Cl",
                  "title": "synthetic", "expno": "1"}
        d._carrier, d._referenced = 0.0, True
        d._period_src = "manual"
        d.period.setValue(period)
        d.periodHz.setValue(sw / period)
        d.nEch.setMaximum(fid.size // period)
        d.nEch.setValue(fid.size // period)
        d.top.setMaximum(period - 1)
        d.top.setValue(top)
        d._loading = False
        d._recompute()
        pump(50)
        for w in (d.btnCsv, d.btnSend):
            w.setEnabled(True)
        if d._spec is None or d._ppm is None:
            raise AssertionError(f"QCPMG processed no spectrum from the synthetic train: {d.res.text()}")
        d._autophase()
        pump(30)
        d._send_infinite()
        if "sent to infinite-field" not in d.res.text():
            raise AssertionError(f"QCPMG → infinite-field: {d.res.text()}")
        d._copy_csv()
        if quick:
            say("quick: the QCPMG figure package (15 files at 600 dpi) is left to the full run")
        else:
            base = out / "qcpmg_fig"
            with _static(QFileDialog, "getSaveFileName", lambda *a, **k: (str(base), "")):
                d._export_package()
            _assert_file(base.with_suffix(".png"), "QCPMG figure package", 1000)
            say(f"QCPMG figure package: {d.res.text()}")
        n_pts = int(d._ppm.size)
        d._send()
        _settle(win)
        if win.exp_ppm.size != n_pts or "qcpmg_period_pts" not in json.dumps(win.recipe):
            raise AssertionError(f"Send to fit did not reach the workbench: {_status(win)}")
        say(f"QCPMG: {n_ech} synthetic echoes → {n_pts}-point sum-echo spectrum on the "
            f"workbench with its qcpmg_* provenance; T2 {d.t2.T2_s:.3g} s"
            if d.t2 is not None and d.t2.ok else
            f"QCPMG: {n_ech} synthetic echoes → {n_pts}-point spectrum on the workbench")

        # -- Infinite-field δiso (the shared dialog that received the row)
        win.open_qcpmg_fields()
        pump(50)
        f = qfd._shared
        if f is None or not f.isVisible() or _rows(f.table) < 1:
            raise AssertionError("open_qcpmg_fields showed no shared dialog with the sent row")
        x = np.linspace(-300.0, 100.0, 2001)
        f.add_dataset_spectrum(160.0, x, np.exp(-(((x + 101.0) / 15.0) ** 2)))
        f._compute()
        if "iso" not in f.result.text():
            raise AssertionError(f"infinite-field compute: {f.result.text()[:200]}")
        rep = out / "infinite_field_report.txt"
        with _static(QFileDialog, "getSaveFileName", lambda *a, **k: (str(rep), "")):
            f._export_report()
        _assert_file(rep, "infinite-field report")
        figb = out / "infinite_field"
        with _static(QFileDialog, "getSaveFileName", lambda *a, **k: (str(figb), "")):
            f._export_figure()
        _assert_file(figb.with_suffix(".png"), "infinite-field figure", 1000)
        n_ds = sum(1 for r in range(_rows(f.table)) if f._row_ds_id(r) is not None)
        if n_ds != 2:
            raise AssertionError(f"the infinite-field dialog holds {n_ds} dataset rows, 2 expected")
        say(f"infinite-field δiso: {n_ds} dataset rows (QCPMG hand-off + added spectrum) · "
            f"{rep.name}, {figb.name}.png/.svg/.pdf")
        f.close()
        qfd._shared = None

        # -- Batch infinite-field: 2 samples × 2 fields
        from larmor.desktop.qcpmg_batch_dialog import QcpmgBatchFieldsDialog

        win.open_qcpmg_batch_fields()
        pump(50)
        b = getattr(win, "_qcpmg_batch_dlg", None)
        if not isinstance(b, QcpmgBatchFieldsDialog) or not b.isVisible():
            raise AssertionError("open_qcpmg_batch_fields kept no visible dialog")
        b.nSamples.setValue(2)
        b.nFields.setValue(2)
        for si, (iso, cq) in enumerate(((-70.0, 3.2), (-80.0, 3.6))):
            for fi, nu in enumerate((78.3541, 107.811)):
                centre = iso - 1e6 * cq ** 2 * 0.02 / nu ** 2
                xx = np.linspace(centre - 400.0, centre + 400.0, 3000)
                p = _write_spectrum(out / f"batch_s{si}f{fi}.csv", xx,
                                    np.exp(-((xx - centre) / 40.0) ** 2), "35Cl", nu)
                b._drop_files(si, 1 + fi, [p])
        if len(b.cells) != 4:
            raise AssertionError(f"the batch grid holds {len(b.cells)} of 4 cells")
        b.table.item(0, 0).setText("sampleA")
        b.table.item(1, 0).setText("sampleB")
        b._compute()
        verdict = b.msg.text()
        if "2 of 2" not in verdict:
            raise AssertionError(f"batch infinite-field: {verdict}")
        brep = out / "batch_infinite_field.txt"
        with _static(QFileDialog, "getSaveFileName", lambda *a, **k: (str(brep), "")):
            b._export_report()
        _assert_file(brep, "batch infinite-field report")
        bfig = out / "batch_infinite_field"
        with _static(QFileDialog, "getSaveFileName", lambda *a, **k: (str(bfig), "")):
            b._export_figures()
        _assert_file(bfig.with_suffix(".png"), "batch infinite-field figure", 1000)
        say(f"batch infinite-field: {verdict} · {b.msg.text()}")
        b.close()

        # -- Read static pattern (C_Q, η) on a simulated static 81Br CT pattern
        from larmor import engine
        from larmor.desktop.staticct_dialog import StaticCtDialog
        from larmor.recipe import Param, Recipe, SiteModel
        from larmor.staticct import ct_static_features

        cq, eta, diso, lar = 30.0, 0.5, -300.0, 216.0
        r = Recipe(nucleus="81Br", larmor_frequency_MHz=lar, spin_rate_Hz=0.0, sites=[
            SiteModel(model="quad_ct", label="s", params={
                "isotropic_chemical_shift_ppm": Param(diso), "Cq_MHz": Param(cq),
                "eta": Param(eta), "shift_fwhm_ppm": Param(10.0), "amplitude": Param(1.0)})])
        gx, gy, _ = engine.simulate(r, exp_ppm=np.linspace(-5200.0, 3000.0, 4001))
        win._display_1d(gx, np.clip(gy, 0, None), "81Br", lar, 0.0, "static 81Br", "")
        _settle(win)

        def drive_static(d):
            feat = ct_static_features(cq, eta, 1.5, lar)
            dists = [min(abs(e - h) for h in feat.horns) for e in feat.edges]
            d.horn_a.setValue(feat.horns[0] + diso)
            d.horn_b.setValue(feat.horns[1] + diso)
            d.edge.setValue(feat.edges[int(np.argmax(dists))] + diso)
            pump(30)
            if d.reading is None or not d.reading.ok:
                raise AssertionError("the static reading is not ok")
            if abs(d.reading.Cq_MHz - cq) / cq > 0.1:
                raise AssertionError(f"static C_Q {d.reading.Cq_MHz:.3g} vs {cq}")
            d._seed()
            say(f"static pattern: C_Q {d.reading.Cq_MHz:.3g} MHz, η {d.reading.eta:.2f} "
                f"(truth {cq}, {eta}) → quad_ct line seeded")

        with _captured(StaticCtDialog, drive_static) as got:
            win.open_staticct()
        _settle(win)
        if not got or not win.recipe["sites"] or win.recipe["sites"][-1]["model"] != "quad_ct":
            raise AssertionError(f"no quad_ct line seeded: {_status(win)}")

        # -- Herzfeld–Berger on a synthetic 31P sideband manifold
        from larmor import herzfeld_berger as hb
        from larmor.desktop.herzfeld_berger_dialog import HerzfeldBergerDialog

        nu0, nur, zeta, eta_cs, centre = 162.0, 5000.0, 60.0, 0.4, -12.0
        xh = np.linspace(-200.0, 200.0, 6001)
        inten = hb.sideband_intensities(zeta, eta_cs, nu0, nur, n_max=6)
        pos = hb.sideband_positions_ppm(centre, nu0, nur, inten)
        sig = 1.5 / (2 * np.sqrt(2 * np.log(2)))
        yh = np.full_like(xh, 0.02)
        for n, v in inten.items():
            yh += v / (sig * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((xh - pos[n]) / sig) ** 2)
        win._display_1d(xh, yh, "31P", nu0, nur, "31P manifold", "")
        _settle(win)

        def drive_hb(d):
            d.centre.setValue(centre)
            d.n_side.setValue(4)
            d._fit()
            if d.fit is None or not d.btnSeed.isEnabled():
                raise AssertionError(f"Herzfeld–Berger fit failed: {d.res.text()}")
            if abs(d.fit.zeta_ppm - zeta) > 3.0:
                raise AssertionError(f"ζ {d.fit.zeta_ppm:.3g} vs {zeta}")
            d._seed()
            say(f"Herzfeld–Berger: ζ {d.fit.zeta_ppm:.3g} ppm, η {d.fit.eta:.2f} "
                f"(truth {zeta}, {eta_cs}) → csa_mas line seeded")

        with _captured(HerzfeldBergerDialog, drive_hb) as got:
            win.open_herzfeld_berger()
        _settle(win)
        if not got or win.recipe["sites"][-1]["model"] != "csa_mas":
            raise AssertionError(f"no csa_mas line seeded: {_status(win)}")

        # -- WURST excitation profile (its own form dialog, accepted) and undo
        before = win.exp_amp.copy()
        with _patched(QDialog, "exec", lambda self: QDialog.Accepted):
            win.open_wurst_correct()
        prov = (win.recipe.get("provenance") or {}).get("wurst_correct")
        if not prov or np.allclose(win.exp_amp, before):
            raise AssertionError(f"WURST correction changed nothing: {_status(win)}")
        win.undo()
        _settle(win)
        if not np.allclose(win.exp_amp, before):
            raise AssertionError("undo did not restore the spectrum after the WURST division")
        say(f"WURST: sweep {prov['sweep_kHz']:g} kHz, N {prov['n']:g}, floor {prov['floor']:g} "
            "divided out and undone")

        # -- Subtract a spectrum (background): load, auto-scale, apply
        from larmor.desktop.subtract_dialog import SubtractDialog

        bg = _write_spectrum(out / "subtract_bg.csv", xh, 0.02 + _gauss(xh, 50.0, 120.0, 0.3),
                             "31P", nu0)

        def drive_subtract(d):
            with _static(QFileDialog, "getOpenFileName", lambda *a, **k: (bg, "")):
                d._load()
            if not d.btnApply.isEnabled():
                raise AssertionError("Subtract: Apply stayed disabled after loading the background")
            d._auto()
            d._apply()

        with _captured(SubtractDialog, drive_subtract) as got:
            win.open_subtract()
        _settle(win)
        if not got or "background" not in (win.recipe.get("sample") or ""):
            raise AssertionError(f"the subtraction did not reach the workbench: {_status(win)}")
        say(f"subtract: {_status(win)}")
        close_tool_windows(keep=[win])


# ============================================================ 5. misc tools
def tools_misc(say, ctx):
    """The utility windows: Conversion tools, the NMR table (every element's
    feasibility fill and tooltip, one isotope detail), the DFT tensor import
    on a synthetic .magres through ``_magres_add_sites``, the 2D MQMAS
    viewer, Co-fit (the example 1D + 2D) and its close, the Figure studio
    with preview and export, the Experimental section (Table S1) with its
    CSV + LaTeX export, the pulse program viewer on both shipped
    sequences, the Session inventory (scan, Fix menu, export) and the
    Referencing audit on the examples folder, the Fit parameter table from
    the open workspaces and from the shipped recipe files."""
    from PySide6.QtWidgets import QFileDialog

    out = _scratch(ctx, "tools")
    ex = _examples(ctx)
    store = out / "store"
    store.mkdir(exist_ok=True)
    with _env(LARMOR_SR_OVERRIDES=store / "sr_overrides.json",
              LARMOR_REF_LOG=store / "referencing_log.jsonl",
              LARMOR_ALIASES=store / "aliases.json",
              LARMOR_RENAME_LOG=store / "rename_log.jsonl"), \
            _desktop(say, ctx) as (win, trap, dis):
        if ex:
            _load(win, ex["al"])
            win.apply_recipe(str(ex["rec_al"]))
            _settle(win)
            say(f"example 27Al 1r with the shipped recipe ({len(win.recipe['sites'])} sites) "
                "as the fitted workspace")
        else:
            _synthetic_27al(win)
            _settle(win)
            say("no bundled example data: synthetic 27Al spectrum")
        fitted_ws = win.active_ws

        # -- Conversion tools: the CSA conventions group
        from larmor.desktop.utilities import ConvertDialog, IsotopeDetailsDialog, NmrTableDialog

        win.open_convert()
        pump(50)
        d = _tool_window(win, ConvertDialog)
        d.cs_iso.setValue(10.0)
        d.cs_zeta.setValue(50.0)
        d.cs_eta.setValue(0.3)
        pump(20)
        if d.cs_d11.value() == 0.0 or d.cs_span.value() <= 0.0:
            raise AssertionError("the CSA conventions did not propagate")
        d.cq.setValue(3.0)
        d.eta.setValue(0.5)
        pump(20)
        say(f"conversion tools: δ11 {d.cs_d11.value():.2f}, span {d.cs_span.value():.2f} ppm; "
            f"quadrupolar: {d.qout.text()[:80]!r}")
        d.close()

        # -- NMR table: every element's fill + tooltip, one detail window
        win.open_nmr_table()
        pump(50)
        d = _tool_window(win, NmrTableDialog)
        fills: dict = {}
        bare = []
        for el, btn in d._buttons.items():
            if not btn.toolTip():
                bare.append(el)
            fills[d.cell_fill(el)] = fills.get(d.cell_fill(el), 0) + 1
        if bare:
            raise AssertionError(f"{len(bare)} element(s) without a tooltip: {bare[:10]}")
        d._details("Al")
        pump(50)
        det = (getattr(d, "_tool_windows", None) or [None])[-1]
        if not isinstance(det, IsotopeDetailsDialog):
            raise AssertionError("the Al cell opened no isotope details window")
        say(f"NMR table: {len(d._buttons)} elements, fills {fills} · details: {det.windowTitle()}")
        det.close()
        d.close()

        # -- Import DFT tensors (.magres) through the launcher's exec-and-apply branch
        from larmor.desktop.magres_dialog import MagresDialog

        magres = out / "kog.nmr.magres"
        magres.write_text(_magres_text(), encoding="utf-8")
        xf = np.linspace(-300.0, -50.0, 2001)
        yf = sum(_gauss(xf, c, 2.0, 1.0) for c in (-250.0, -238.0, -225.0, -209.0, -196.0, -182.0))
        win._display_1d(xf, yf, "19F", 564.27, 35714.0, "kog", "")
        _settle(win)

        def drive_magres(d):
            d._load(str(magres))
            d.ref.setValue(560.0)                 # σ_ref conversion
            pump(30)
            if not d.btnAdd.isEnabled():
                raise AssertionError(f"magres Add stayed disabled: {d.warn.text()}")
            d._accept()
            if d.result is None:
                raise AssertionError("the magres dialog accepted without a result")
            return 1

        with _captured(MagresDialog, drive_magres, code=1) as got:
            win.open_magres()
        _settle(win)
        if not got or len(win.recipe["sites"]) != 6 \
                or "dft_import" not in (win.recipe.get("provenance") or {}):
            raise AssertionError(f"the DFT import did not land: {_status(win)}")
        say(f"magres import: {_status(win)}")

        # -- 2D MQMAS viewer on the example 2D
        from larmor.desktop.twod_dialog import TwoDDialog

        def drive_twod(d):
            if d.data is None and ex:
                d._load(str(ex["twod"]))
            if ex and d.data is None:
                raise AssertionError(f"the 2D viewer read nothing: {d.res.text()}")
            if d.data is not None:
                d.nlevels.setValue(8)
                pump(30)
                say(f"2D viewer: {d.res.text()} · {len(d.p_main.items)} contour items")
            else:
                say("2D viewer built without data (no example 2D)")

        if ex:
            _load(win, ex["twod"])
            if win._data2d is None:
                raise AssertionError(f"the example 2D did not open as a map: {_status(win)}")
        with _captured(TwoDDialog, drive_twod) as got:
            win.open_twod()
        if not got:
            raise AssertionError("open_twod built no TwoDDialog")

        # -- Co-fit: the fitted 1D plus the example 2D, simulated, closed
        win.switch_workspace(fitted_ws)
        _settle(win)
        if ex:
            t0 = time.perf_counter()
            with _static(QFileDialog, "getOpenFileName", lambda *a, **k: ("", "")), \
                    _static(QFileDialog, "getExistingDirectory", lambda *a, **k: str(ex["twod"])):
                win.open_cofit()
            st = win._cofit
            if not st or not st.get("d1") or not st.get("d2"):
                raise AssertionError(f"co-fit did not take both datasets: {_status(win)}")
            win._cofit_simulate_now()
            pump(50)
            if win.cofit_view1d._model is None or win.cofit_view2d._model_sites is None:
                raise AssertionError("the co-fit panels show no model overlay")
            say(f"co-fit: 1D {st['d1'][2]} + 2D {st['d2'][1]} simulated "
                f"({time.perf_counter() - t0:.1f} s, MQMAS kernel) · {_status(win)}")
            win.close_cofit()
            if win.central_stack.currentWidget() is win.cofit_page:
                raise AssertionError("close_cofit left the co-fit page showing")
            _settle(win)
        else:
            say("co-fit needs the example 2D: skipped")

        # -- Figure studio: template, preview, export, kept in the dock
        from larmor.desktop.figure_dialog import FigureDialog

        fig_base = out / "figure"

        def drive_figure(d):
            d.load_templates()
            pump(30)
            buttons = [d.tpl_box.itemAt(i).widget() for i in range(d.tpl_box.count())]
            buttons = [b for b in buttons if b is not None]
            if buttons:
                buttons[0].click()
            else:
                src = str(ex["al"]) if ex else ""
                traces = ([{"path": src, "label": "experiment"}] if src else
                          [{"data": {"x": win.exp_ppm.tolist(), "y": win.exp_amp.tolist()},
                            "label": "experiment"}])
                d.spec.setPlainText(json.dumps({"kind": "1d", "style": "article-wide",
                                                "traces": traces}))
            d.preview()
            pix = d.preview_label.pixmap()
            if pix is None or pix.isNull():
                raise AssertionError("the Figure studio preview is empty")
            d.fmt_pdf.setChecked(False)
            with _static(QFileDialog, "getSaveFileName", lambda *a, **k: (str(fig_base), "")):
                d.export()
            _assert_file(fig_base.with_suffix(".png"), "figure export", 1000)
            _assert_file(fig_base.with_suffix(".svg"), "figure export (svg)", 1000)
            say(f"figure studio: {'template' if buttons else 'inline spec'}, preview "
                f"{pix.width()}×{pix.height()}, exported {fig_base.name}.png/.svg")

        n_ws = len(win.workspaces)
        with _captured(FigureDialog, drive_figure) as got:
            win.open_figure_dialog()
        if not got or not any(ws["kind"] == "figure" for ws in win.workspaces[n_ws:]):
            raise AssertionError(f"the figure was not kept in the Workspaces dock: {_status(win)}")

        # -- Experimental section (Table S1) on the example EXPNOs
        from larmor.desktop.acquisition_dialog import AcquisitionTableDialog

        if ex:
            win.open_acquisition_table([str(ex["al"]), str(ex["b"]), str(ex["twod"])])
            pump(100)
            d = getattr(win, "_acq_dlg", None)
            if not isinstance(d, AcquisitionTableDialog) or _rows(d.table) != 3:
                raise AssertionError("the acquisition table does not list the 3 EXPNOs")
            d.btnCopyPara.click()
            para = d.para.toPlainText()
            if len(para) < 40:
                raise AssertionError(f"the Experimental paragraph is empty: {para!r}")
            acq = out / "acq.csv"
            with _static(QFileDialog, "getSaveFileName", lambda *a, **k: (str(acq), "")):
                d.btnSave.click()
            _assert_file(acq, "acquisition CSV")
            _assert_file(acq.with_suffix(".tex"), "acquisition LaTeX")
            say(f"experimental section: 3 rows · {d.status.text()} · paragraph {len(para)} chars")
            d.close()
        else:
            win.open_acquisition_table([])
            say("acquisition table: no EXPNO bundled, opened empty")

        # -- Pulse program viewer on both shipped sequences
        from larmor.desktop.pulseprog_dialog import PulseProgramDialog

        if ex:
            for key in ("b", "twod"):
                win.open_pulse_program(str(ex[key]))
                pump(80)
                d = _tool_window(win, PulseProgramDialog)
                tl = d.diagram.timeline()
                if tl is None or not tl.channels:
                    raise AssertionError(f"no timeline drawn for {ex[key].name}: {d.notes.text()}")
                png = out / f"pulseprog_{ex[key].name}.png"
                written = d.export_to(str(png))
                _assert_file(Path(written), "pulse program diagram", 1000)
                say(f"pulse program {ex[key].name}: channels {tl.channels} · {d.summary.text()[:100]}")
                d.close()
        else:
            win.open_pulse_program()
            say("pulse program:", _status(win))

        # -- Session inventory on the examples folder: scan, Fix menu, export
        from PySide6.QtWidgets import QMenu

        from larmor.desktop.inventory_dialog import SessionInventoryDialog

        if ex:
            win.open_session_inventory(str(ex["root"]))
            pump(50)
            d = getattr(win, "_inventory_dlg", None)
            if not isinstance(d, SessionInventoryDialog):
                raise AssertionError("open_session_inventory kept no dialog")
            d.scan()
            pump(50)
            if _rows(d.grid) < 1:
                raise AssertionError(f"the inventory scan found no sample: {d.status.text()}")
            d.grid.setCurrentCell(0, 0)
            d._grid_selected()
            pump(30)
            if _rows(d.detail) < 1:
                raise AssertionError("the inventory detail table is empty")
            d.detail.setCurrentCell(0, 1)
            row = d.current_row()
            menu = d.fix_menu(row)
            if not isinstance(menu, QMenu) or not menu.actions():
                raise AssertionError("the Fix menu is empty")
            inv = out / "inventory.csv"
            with _static(QFileDialog, "getSaveFileName", lambda *a, **k: (str(inv), "")):
                d.export()
            _assert_file(inv, "inventory CSV")
            say(f"session inventory: {d.status.text()[:140]} · Fix menu {len(menu.actions())} "
                f"actions · {inv.name}")
            d.close()

            # -- Referencing audit on the same folder
            from larmor.desktop.referencing_dialog import ReferencingAuditDialog

            win.open_referencing_audit()
            pump(50)
            d = getattr(win, "_ref_audit", None)
            if not isinstance(d, ReferencingAuditDialog):
                raise AssertionError("open_referencing_audit kept no dialog")
            d.folder.setText(str(ex["root"]))
            d.scan()
            pump(50)
            d.chkAll.setChecked(True)
            d._fill_rows()
            if _rows(d.table) < 1:
                raise AssertionError(f"the referencing audit lists no EXPNO: {d.status.text()}")
            audit = out / "audit.csv"
            with _static(QFileDialog, "getSaveFileName", lambda *a, **k: (str(audit), "CSV (*.csv)")):
                d.export()
            _assert_file(audit, "referencing audit CSV")
            say(f"referencing audit: {_rows(d.table)} rows · {d.status.text()[:120]} · {audit.name}")
            d.close()
        else:
            say("session inventory / referencing audit need the examples folder: skipped")

        # -- Fit parameter table: the open fitted workspaces, then the recipe files
        from larmor.desktop.fittable_dialog import FitTableDialog

        win.switch_workspace(fitted_ws)
        _settle(win)
        win.open_fit_table()
        pump(50)
        d = _tool_window(win, FitTableDialog)
        if _rows(d.table) < 1:
            raise AssertionError(f"the fit table lists no fitted workspace: {_status(win)}")
        ft = out / "fit_table.csv"
        if not d.export_csv(ft):
            raise AssertionError("fit table CSV export returned False")
        _assert_file(ft, "fit table CSV")
        tsv = d.copy_tsv()
        say(f"fit table (workspaces): {_rows(d.table)} rows · {ft.name} · "
            f"{len(tsv.splitlines())} TSV lines")
        d.close()
        if ex:
            win.open_fit_table_files([str(ex["rec_al"]), str(ex["rec_b"])])
            pump(50)
            d = _tool_window(win, FitTableDialog)
            if _rows(d.table) != 6:
                raise AssertionError(f"the two shipped recipes give {_rows(d.table)} rows, 6 expected")
            d.chkLong.setChecked(True)
            pump(20)
            ftl = out / "fit_table_long.csv"
            if not d.export_csv(ftl):
                raise AssertionError("long fit table CSV export returned False")
            _assert_file(ftl, "long fit table CSV")
            say(f"fit table (files): 6 wide rows, {_rows(d.table)} long rows · {ftl.name}")
            d.close()
        close_tool_windows(keep=[win])
