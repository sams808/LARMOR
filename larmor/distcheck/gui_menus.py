"""Desktop sweep, part 1: the main window built offscreen, every menu
action triggered with its dialogs auto-dismissed -- in three window states --
every theme and text size applied, every manual and tutorial rendered, every
View toggle flipped. Each stage is ``(name, fn(say, ctx))``; plumbing in
``gui_common``.

What a stage proves is that the DESKTOP RUNS HERE: no module missing from the
bundle behind a menu entry, no exception in any slot (PySide6 only prints
those; ``ExceptionTrap`` collects them), every page and theme rendering. The
spectra the sweeps act on are copies of the bundled examples inside
``ctx["out"]`` (a synthetic CSV when the examples do not ship), so an action
that writes next to its data writes into the scratch folder; the
installation's own example files, the user's home folder and the current
folder are checked untouched afterwards, and the user's saved settings come
back exactly as they were (``SettingsGuard``). Nothing reaches the network:
``webbrowser.open`` / ``os.startfile`` are recorders for the duration and the
? ▸ More card runs on its offline pack.

The menu is walked through ``larmor.desktop.palette_dialog.collect_commands``
-- the same leaf list the command palette offers, in menu order -- so the
sweep covers what the palette covers, and never through ``QAction.menu()``
(which re-parents the window's own menu wrappers in PySide6).
"""
from __future__ import annotations

import collections
import contextlib
import gc
import os
import re
import shutil
import time
import traceback
from pathlib import Path

import numpy as np

from larmor.distcheck.gui_common import (ExceptionTrap, ModalDismisser, SettingsGuard,
                                         close_tool_windows, ensure_app, pump)

__all__ = ["stages", "EXPECTED_LEAF_ACTIONS", "SKIP"]

#: the leaf rows (no submenu) of tests/test_app_split.GOLDEN_MENU when this
#: was written: a sweep must find at least this many actions outside the two
#: menus made of the user's own files (Open recent, Apply recipe), or the menu
#: bar was not really walked. tests/test_distcheck_gui_menus.py checks the
#: number against the golden tree, so it cannot go stale silently.
EXPECTED_LEAF_ACTIONS = 198

#: actions the sweep never triggers, keyed by (menu path, label) the way the
#: command palette spells them (mnemonics stripped), with the reason
SKIP = {
    (("File",), "Quit"): "closes the window the sweep is driving",
}

#: submenus whose entries are the user's own files, never triggered (the
#: file dialog of Apply recipe ▸ Browse for recipe… is -- and dismissed)
_USER_FILE_MENUS = {
    ("File", "Open recent"): "a recent file of this user",
    ("Edit", "Apply recipe"): "a recent recipe of this user",
}

#: how long a fit started by a swept action may run before the sweep asks it
#: to stop (keeping the values so far) -- the point is the slot, not the fit
_FIT_STOP_S = {"quick": 4.0, "full": 30.0}
#: the hard bound on waiting for the window to go idle after one action
_IDLE_TIMEOUT_S = 60.0


# ------------------------------------------------------------------ data
def _write_synthetic(folder: Path) -> tuple[Path, Path]:
    """A two-line 11B gauss_lor spectrum as a LARMOR CSV plus its recipe --
    the stand-in when the bundled examples are not in this installation."""
    from larmor import engine
    from larmor.recipe import Param, Recipe, SiteModel

    folder.mkdir(parents=True, exist_ok=True)
    x = np.linspace(-20, 60, 600)
    rec = Recipe(nucleus="11B", larmor_frequency_MHz=160.0, spin_rate_Hz=0.0, sites=[
        SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(15.0), "shift_fwhm_ppm": Param(6.0),
            "amplitude": Param(100.0), "gl": Param(1.0, vary=False)}),
        SiteModel(model="gauss_lor", label="B", params={
            "isotropic_chemical_shift_ppm": Param(2.0), "shift_fwhm_ppm": Param(4.0),
            "amplitude": Param(40.0), "gl": Param(1.0, vary=False)})])
    _, model, _ = engine.simulate(rec, exp_ppm=x)
    y = model + np.random.default_rng(0).normal(0.0, 0.5, model.size)
    csv = folder / "synthetic_11B.csv"
    with open(csv, "w", encoding="utf-8") as f:
        f.write("# nucleus = 11B\n# larmor_MHz = 160\n")
        for xi, yi in zip(x, y):
            f.write(f"{xi:.4f} {yi:.4f}\n")
    recipe = folder / "synthetic_11B.recipe.json"
    rec.save(recipe)
    return csv, recipe


class _Data:
    """The spectra and the recipe the sweeps act on, copied under
    ``ctx["out"]``: the bundled 27Al 1r (EXPNO 3616) with its recipe and the
    MQMAS 2rr (EXPNO 3620) when the examples ship, the synthetic CSV
    otherwise. ``written()`` lists what an action left next to them."""

    def __init__(self, ctx):
        self.root = Path(ctx["out"]) / "menus" / "data"
        self.root.mkdir(parents=True, exist_ok=True)
        examples = ctx.get("examples")
        self.synthetic = examples is None
        self.spectrum_2d = None
        if examples is not None:
            tree = self.root / "pCABS2-4"
            if not tree.exists():
                shutil.copytree(Path(examples) / "pCABS2-4", tree)
            for name in ("pCABS2-4_27Al.recipe.json", "pCABS2-4_11B.recipe.json"):
                src = Path(examples) / name
                if src.is_file() and not (self.root / name).exists():
                    shutil.copy(src, self.root / name)
            self.spectrum = tree / "3616" / "pdata" / "1" / "1r"
            self.recipe = self.root / "pCABS2-4_27Al.recipe.json"
            two_d = tree / "3620" / "pdata" / "1" / "2rr"
            if two_d.is_file():
                self.spectrum_2d = two_d
            if not (self.spectrum.is_file() and self.recipe.is_file()):
                raise AssertionError(
                    f"the bundled examples are incomplete: {self.spectrum} / {self.recipe}")
        else:
            self.spectrum, self.recipe = _write_synthetic(self.root)
        self._initial = _tree_signature(self.root)

    def describe(self) -> str:
        if self.synthetic:
            return f"synthetic 11B CSV ({self.spectrum.name}) + {self.recipe.name}"
        return (f"examples copy: {self.spectrum.relative_to(self.root)} + {self.recipe.name}"
                + (f" + 2D {self.spectrum_2d.relative_to(self.root)}" if self.spectrum_2d else ""))

    def written(self) -> list[str]:
        now = _tree_signature(self.root)
        return sorted(p for p in now if p not in self._initial or now[p] != self._initial[p])


def _tree_signature(root) -> dict:
    out = {}
    root = Path(root)
    for p in root.rglob("*"):
        if p.is_file():
            try:
                st = p.stat()
            except OSError:
                continue
            out[str(p.relative_to(root))] = (st.st_size, int(st.st_mtime))
    return out


class _Outside:
    """What a swept action must never do without a dialog: write outside the
    scratch folder. New entries in the user's home folder and in the current
    folder, and any change under the installation's own examples, between
    construction and ``check()``."""

    def __init__(self, ctx):
        from larmor.distcheck import log_path

        lp = log_path()
        self._ignore = {lp.name, lp.with_suffix(".json").name, "larmor_crash.log"}
        self._roots = []
        for r in (Path.home(), Path.cwd()):
            if r not in self._roots:
                self._roots.append(r)
        self._before = {r: self._listing(r) for r in self._roots}
        self._examples = ctx.get("examples")
        self._ex_before = (_tree_signature(self._examples)
                           if self._examples is not None else None)

    def _listing(self, root: Path) -> set:
        try:
            return {n for n in os.listdir(root)
                    if n not in self._ignore and not n.startswith(".")
                    and n != "__pycache__"}
        except OSError:
            return set()

    def check(self) -> list[str]:
        out = []
        for r in self._roots:
            new = self._listing(r) - self._before[r]
            if new:
                out.append(f"new entries in {r} during the sweep: {sorted(new)}")
        if self._ex_before is not None and _tree_signature(self._examples) != self._ex_before:
            out.append(f"the installation's examples changed under {self._examples}")
        return out


# --------------------------------------------------------------- threads
def _thread_ptr(t) -> int:
    import shiboken6
    return shiboken6.getCppPointer(t)[0]


def _running_threads() -> list:
    """Every Python-owned QThread still running (FitWorker, SimWorker,
    SeqWorker, the More card's FactWorker...), the main thread excluded."""
    import shiboken6
    from PySide6.QtCore import QThread

    main = _thread_ptr(QThread.currentThread())
    out = []
    for obj in gc.get_objects():
        if isinstance(obj, QThread):
            try:
                if shiboken6.isValid(obj) and _thread_ptr(obj) != main and obj.isRunning():
                    out.append(obj)
            except RuntimeError:            # deleted on the C++ side meanwhile
                pass
    return out


def _join_threads(timeout: float = 15.0) -> list[str]:
    """Wait (pumping) for every running QThread; the names of those that did
    not finish in ``timeout`` s each."""
    left = []
    for t in _running_threads():
        end = time.perf_counter() + timeout
        try:
            while t.isRunning() and time.perf_counter() < end:
                pump(20)
            if t.isRunning():
                left.append(type(t).__name__)
        except RuntimeError:
            pass
    return left


# --------------------------------------------------------------- rendering
def _paints(widget) -> tuple[bool, str]:
    """Grab ``widget`` to an image (offscreen, never shown) and say whether
    something was drawn: more than two distinct colours over a sampling grid."""
    img = widget.grab().toImage()
    if img.isNull() or img.width() < 8 or img.height() < 8:
        return False, f"empty image ({img.width()}×{img.height()})"
    colours = set()
    for y in range(0, img.height(), max(1, img.height() // 24)):
        for x in range(0, img.width(), max(1, img.width() // 32)):
            colours.add(img.pixel(x, y))
        if len(colours) > 12:
            break
    n = f"{len(colours)}{'+' if len(colours) > 12 else ''}"
    return len(colours) > 2, f"{img.width()}×{img.height()} px, {n} colours sampled"


# --------------------------------------------------------------- watchdog
class _Watchdog:
    """A 1 s QTimer alongside the sweep: closes the popup menus an action
    opens (``ModalDismisser`` sees only modal widgets, a ``QMenu.exec()``
    would block forever), and names in the log the action a nested event
    loop has been sitting in for too long, with the windows on screen -- so a
    hang in the field is attributable from the log alone."""

    def __init__(self, say, warn_after: float = 20.0):
        from PySide6.QtCore import QTimer

        self.say = say
        self.warn_after = warn_after
        self.current = ""
        self.since = 0.0
        self.popups: list = []
        self._warned = False
        self._timer = QTimer()
        self._timer.setInterval(1000)
        self._timer.timeout.connect(self._tick)

    def __enter__(self):
        self._timer.start()
        return self

    def __exit__(self, *exc):
        self._timer.stop()
        return False

    def begin(self, name: str):
        self.current, self.since, self._warned = name, time.perf_counter(), False

    def end(self):
        self.current = ""

    def _tick(self):
        from PySide6.QtWidgets import QApplication

        p = QApplication.activePopupWidget()
        if p is not None:
            self.popups.append(f"{type(p).__name__} ({self.current or 'between actions'})")
            try:
                p.close()
            except Exception:                           # noqa: BLE001
                pass
        if self.current and not self._warned and time.perf_counter() - self.since > self.warn_after:
            self._warned = True
            modal = QApplication.activeModalWidget()
            shown = [f"{type(w).__name__} {w.windowTitle()!r}"
                     for w in QApplication.topLevelWidgets() if w.isVisible()]
            self.say(f"    still inside {self.current} after {self.warn_after:.0f} s; "
                     f"active modal: {type(modal).__name__ if modal is not None else None}; "
                     f"windows: {shown}")


# --------------------------------------------------------------- harness
class _Harness:
    """Everything a desktop stage needs around its window: the settings
    guard, the exception trap, the modal dismisser, the stubs that keep the
    sweep offline, the data copy, the write watch, and a teardown that closes
    every window and joins every thread. Use as a context manager; ``build()``
    makes the window, ``load()`` puts the spectrum and the recipe on it."""

    def __init__(self, say, ctx):
        self.say = say
        self.ctx = ctx
        self.quick = bool(ctx.get("quick"))
        self.win = None
        self.urls: list = []
        self.opened: list = []
        self.leftover_threads: list = []

    # ------------------------------------------------------------ enter/exit
    def __enter__(self):
        ensure_app()
        self._stack = contextlib.ExitStack()
        self._stack.__enter__()
        self._stack.enter_context(SettingsGuard())
        self.trap = self._stack.enter_context(ExceptionTrap())
        self.modals = self._stack.enter_context(ModalDismisser())
        self.watchdog = self._stack.enter_context(_Watchdog(self.say))
        self.trace = bool(os.environ.get("LARMOR_DISTCHECK_TRACE"))
        self._install_stubs()
        self._stack.callback(self._restore_stubs)
        self.data = _Data(self.ctx)
        self.outside = _Outside(self.ctx)
        return self

    def __exit__(self, exc_type, exc, tb):
        try:
            if self.win is not None:
                self.leftover_threads = self.close_window()
        finally:
            return self._stack.__exit__(exc_type, exc, tb)

    def _install_stubs(self):
        import webbrowser
        self._wb_open = webbrowser.open
        webbrowser.open = lambda url, *a, **k: (self.urls.append(str(url)), True)[1]
        self._startfile = getattr(os, "startfile", None)
        if self._startfile is not None:
            os.startfile = lambda p, *a, **k: self.opened.append(str(p))
        self._http = None
        try:
            from larmor.xfact.core import http
            self._http = (http.get_json, http.get_bytes)
            http.get_json = lambda *a, **k: None       # the card falls back to
            http.get_bytes = lambda *a, **k: None      # its offline pack
        except Exception:                               # noqa: BLE001
            pass

    def _restore_stubs(self):
        import webbrowser
        webbrowser.open = self._wb_open
        if self._startfile is not None:
            os.startfile = self._startfile
        if self._http is not None:
            from larmor.xfact.core import http
            http.get_json, http.get_bytes = self._http

    # ------------------------------------------------------------ the window
    def build(self):
        from larmor.desktop.app import MainWindow

        win = MainWindow()
        win._add_recent = lambda p: None               # the user's recent list stays
        win._report_project_notes = lambda notes: None  # never a modal box here
        win._confirm_open_mode = lambda: "replace"
        self.win = win
        pump(50)
        return win

    def load(self):
        """The 1D spectrum on the workbench with the recipe's lines on it,
        simulated."""
        self.win.load_source(str(self.data.spectrum), keep_fit=False)
        pump(50)
        self.win.apply_recipe(str(self.data.recipe))
        self.settle()

    def busy(self) -> bool:
        win = self.win
        try:
            sw = win._sim_worker
            if win._sim_timer.isActive() or win._sim_pending or (sw is not None and sw.isRunning()):
                return True
            if win._active_fit_worker is not None:
                return True
            for name in ("_fit_worker", "_fit2d_worker", "_seq_worker"):
                w = getattr(win, name, None)
                if w is not None and w.isRunning():
                    return True
            t = getattr(win, "_cofit_timer", None)
            if t is not None and t.isActive():
                return True
        except RuntimeError:                            # a worker deleted meanwhile
            return False
        return False

    def settle(self, timeout: float = _IDLE_TIMEOUT_S, stop_after: float | None = None) -> bool:
        """Pump until the window is idle (no simulation, fit or sweep in
        flight). A fit still running after ``stop_after`` s is asked to stop
        -- it keeps its values, the slot is what the sweep tests. Returns
        False when the window is still busy at ``timeout``."""
        t0 = time.perf_counter()
        stopped = False
        while self.busy():
            dt = time.perf_counter() - t0
            if dt > timeout:
                return False
            if (stop_after is not None and dt > stop_after and not stopped
                    and self.win._active_fit_worker is not None):
                self.win._interrupt_fit("stop")
                stopped = True
            pump(20)
        pump(20)
        return True

    def cancel_workers(self):
        try:
            self.win._interrupt_fit("cancel")
        except Exception:                               # noqa: BLE001
            pass
        pump(200)

    def wait_more_card(self, timeout: float = 10.0):
        """The ? ▸ More card fetches in a QThread parented to the dialog:
        join it before anything closes the dialog."""
        dlg = getattr(self.win, "_more_dlg", None)
        w = getattr(dlg, "_worker", None) if dlg is not None else None
        if w is None:
            return
        end = time.perf_counter() + timeout
        try:
            while w.isRunning() and time.perf_counter() < end:
                pump(20)
        except RuntimeError:
            pass

    def close_window(self) -> list[str]:
        """Close every tool window, join every thread, close and delete the
        main window; the names of the threads still running afterwards."""
        win = self.win
        if win is None:
            return []
        try:
            self.settle(timeout=_IDLE_TIMEOUT_S, stop_after=0.0)
        except RuntimeError:
            pass
        self.wait_more_card()
        close_tool_windows(keep=[win])
        left = _join_threads()
        self.win = None
        try:
            win.close()
            pump(100)
            win.deleteLater()
            pump(100)
        except RuntimeError:
            pass
        gc.collect()
        return left

    # ------------------------------------------------------------ the states
    def _leave_modes(self):
        """Out of every transient mode an action may have left on: the series
        bar, the FID view, the add-line mode, the co-fit page."""
        win = self.win
        if getattr(win, "_series", None) is not None:
            win.end_series()
        if getattr(win.view, "domain", "freq") == "time":
            win._toggle_time_domain(False)
        win._set_add_mode(None)
        if win.central_stack.currentWidget() is win.cofit_page:
            win.close_cofit()

    def restore_loaded(self):
        """The 1D workbench showing the spectrum with the recipe's lines."""
        win, d = self.win, self.data
        self.settle(stop_after=0.0)
        self._leave_modes()
        if (win.central_stack.currentWidget() is not win.view or not len(win.exp_ppm)
                or str(win.source_path) != str(d.spectrum)):
            win.load_source(str(d.spectrum), keep_fit=False)
            pump(30)
        if not win.recipe or not win.recipe.get("sites"):
            win.apply_recipe(str(d.recipe))
        self.settle()

    def restore_empty(self):
        """The 1D workbench showing the spectrum with NO lines."""
        win, d = self.win, self.data
        self.settle(stop_after=0.0)
        self._leave_modes()
        if (win.central_stack.currentWidget() is not win.view or not len(win.exp_ppm)
                or str(win.source_path) != str(d.spectrum)):
            win.load_source(str(d.spectrum), keep_fit=False)
            pump(30)
        if win.recipe and win.recipe.get("sites"):
            win.new_fit()
        self.settle()

    def restore_2d(self):
        """The 2D contour map of the example MQMAS on screen."""
        win, d = self.win, self.data
        self.settle(stop_after=0.0)
        self._leave_modes()
        if (win.central_stack.currentWidget() is not win.view2d
                or str(win.source_path) != str(d.spectrum_2d)):
            win.load_source(str(d.spectrum_2d), keep_fit=False)
            pump(30)
        self.settle()

    # ------------------------------------------------------------ the walk
    def commands(self) -> list:
        """The leaf actions of the menu bar, in menu order (the toolbar's
        entries left out), as the command palette lists them."""
        from larmor.desktop.palette_dialog import collect_commands

        return [c for c in collect_commands(self.win) if c.path and c.path[0] != "Toolbar"]

    def sweep(self, label: str, restore) -> list[dict]:
        """Trigger every menu action once in the current window state, pump
        and wait for its workers, dismiss its dialogs, close its tool windows,
        then put the state back with ``restore``. One record per entry."""
        import shiboken6

        win = self.win
        cmds = self.commands()
        rows = [(c.path, c.label) for c in cmds]
        by_row = {(c.path, c.label): c for c in cmds}
        stop_after = _FIT_STOP_S["quick" if self.quick else "full"]
        records = []
        for path, lab in rows:
            rec = {"pass": label, "path": path, "label": lab, "status": "", "reason": "",
                   "modals": [], "errors": [], "seconds": 0.0,
                   "golden": path[:2] not in _USER_FILE_MENUS}
            records.append(rec)
            cmd = by_row.get((path, lab))
            reason = _skip_reason(cmd)
            if reason:
                rec["status"], rec["reason"] = "skipped", reason
                continue
            action = cmd.action
            if not shiboken6.isValid(action):          # a rebuilt dynamic menu
                fresh = {(c.path, c.label): c for c in self.commands()}
                cmd = fresh.get((path, lab))
                if cmd is None:
                    rec["status"], rec["reason"] = "gone", "the entry disappeared (menu rebuilt)"
                    continue
                action = cmd.action
            if not action.isEnabled():
                rec["status"] = "disabled"
                continue
            n_modal = len(self.modals.seen)
            t0 = time.perf_counter()
            if self.trace:
                self.say(f"    → {_name(rec)}")
            self.watchdog.begin(_name(rec))
            try:
                action.trigger()
            except Exception:                           # noqa: BLE001
                self.trap.errors.append(traceback.format_exc())
            pump(40)
            if not self.settle(stop_after=stop_after):
                self.trap.errors.append(
                    f"{_name(rec)}: the window was still busy {_IDLE_TIMEOUT_S:.0f} s after "
                    "the action (a worker did not finish)")
                self.cancel_workers()
            self.wait_more_card()
            close_tool_windows(keep=[win])
            try:
                restore()
            except Exception:                           # noqa: BLE001
                self.trap.errors.append(
                    f"{_name(rec)}: putting the window state back failed:\n"
                    + traceback.format_exc())
            self.watchdog.end()
            rec["seconds"] = time.perf_counter() - t0
            rec["modals"] = list(self.modals.seen[n_modal:])
            rec["errors"] = self.trap.take()
            rec["status"] = "error" if rec["errors"] else "ok"
            if self.trace:
                self.say(f"      {rec['status']} in {rec['seconds']:.1f} s"
                         + (f", dismissed {rec['modals']}" if rec["modals"] else ""))
        return records


def _skip_reason(cmd) -> str:
    if cmd is None:
        return ""
    r = SKIP.get((tuple(cmd.path), cmd.label))
    if r:
        return r
    menu = tuple(cmd.path[:2])
    if menu in _USER_FILE_MENUS and not (menu == ("Edit", "Apply recipe")
                                         and cmd.label.startswith("Browse")):
        return _USER_FILE_MENUS[menu]
    return ""


def _name(rec: dict) -> str:
    return " ▸ ".join(tuple(rec["path"]) + (rec["label"],))


def _summarize(say, records: list[dict]) -> list[str]:
    """Log one pass of the sweep (counts, per menu, dialogs, the slowest
    entries, the failures in full) and return its problems."""
    label = records[0]["pass"] if records else "?"
    ok = [r for r in records if r["status"] == "ok"]
    err = [r for r in records if r["status"] == "error"]
    skipped = [r for r in records if r["status"] == "skipped"]
    disabled = [r for r in records if r["status"] == "disabled"]
    gone = [r for r in records if r["status"] == "gone"]
    golden = [r for r in records if r["golden"]]
    total_s = sum(r["seconds"] for r in records)
    say(f"[{label}] {len(records)} menu entries: {len(ok) + len(err)} triggered "
        f"({len(err)} failing), {len(skipped)} skipped, {len(disabled)} disabled, "
        f"{len(gone)} gone, {len(records) - len(golden)} in the user-file menus; "
        f"{total_s:.0f} s")
    for menu, group in collections.OrderedDict(
            (r["path"][0], None) for r in records).items():
        rs = [r for r in records if r["path"][0] == menu]
        n_modals = sum(len(r["modals"]) for r in rs)
        say(f"    {menu}: {sum(r['status'] == 'ok' for r in rs)} ok, "
            f"{sum(r['status'] == 'error' for r in rs)} failing, "
            f"{sum(r['status'] == 'skipped' for r in rs)} skipped, "
            f"{sum(r['status'] == 'disabled' for r in rs)} disabled, "
            f"{n_modals} dialogs dismissed, {sum(r['seconds'] for r in rs):.1f} s")
    for r in skipped:
        say(f"    skipped {_name(r)}: {r['reason']}")
    if disabled:
        say("    disabled (Qt ignores trigger()): " + ", ".join(_name(r) for r in disabled))
    if gone:
        say("    gone: " + ", ".join(_name(r) for r in gone))
    titles = collections.Counter(t for r in records for t in r["modals"])
    if titles:
        say(f"    {sum(titles.values())} dialogs dismissed: "
            + ", ".join(f"{t} ×{n}" if n > 1 else t for t, n in titles.most_common()))
    slow = sorted(records, key=lambda r: -r["seconds"])[:5]
    say("    slowest: " + ", ".join(f"{_name(r)} {r['seconds']:.1f} s" for r in slow))
    problems = []
    if len(golden) < EXPECTED_LEAF_ACTIONS:
        problems.append(
            f"[{label}] only {len(golden)} leaf actions found outside the user-file menus; "
            f"the pinned menu tree has {EXPECTED_LEAF_ACTIONS}: the menu bar was not fully walked")
    expected_triggered = len(golden) - sum(1 for r in golden if r["status"] in
                                           ("skipped", "disabled", "gone"))
    if sum(1 for r in golden if r["status"] in ("ok", "error")) < expected_triggered:
        problems.append(f"[{label}] fewer actions triggered than found")
    if label.startswith("loaded") and len(disabled) > 4:
        problems.append(f"[{label}] {len(disabled)} entries disabled with a spectrum and a "
                        f"fit on screen: {', '.join(_name(r) for r in disabled)}")
    for r in err:
        head = r["errors"][0].strip().splitlines()[-1] if r["errors"] else "?"
        problems.append(f"[{label}] {_name(r)}: {head}")
        say(f"    FAIL {_name(r)}:")
        for e in r["errors"]:
            say(e.rstrip())
    return problems


def _finish(h: _Harness, problems: list[str]):
    """Common tail of every stage: close the window, join the threads, check
    the write watch, raise once with every problem named."""
    problems = list(problems)
    left = h.close_window()
    if left:
        problems.append(f"QThread(s) still running after the window closed: {left}")
    problems += [e.strip().splitlines()[-1] for e in h.trap.take()]
    problems += h.outside.check()
    if h.urls or h.opened:
        h.say(f"external openers intercepted: urls {h.urls}, files {h.opened}")
    if h.watchdog.popups:
        h.say(f"popup menus closed by the watchdog: {h.watchdog.popups}")
    written = h.data.written()
    if written:
        h.say(f"written next to the data (inside the scratch folder): {written}")
    if h.trap.warnings:
        uniq = collections.Counter(w.split("\n")[0][:120] for w in h.trap.warnings)
        h.say(f"{len(h.trap.warnings)} Qt warnings (not failures), most frequent: "
              + "; ".join(f"{w} ×{n}" for w, n in uniq.most_common(3)))
    if problems:
        raise AssertionError(f"{len(problems)} problem(s):\n  " + "\n  ".join(problems))


# ================================================================ stages
def menu_window(say, ctx):
    """The main window builds offscreen with its central views, docks, lines
    table, Workspaces panel and status bar; the example spectrum loads, the
    shipped recipe applies, the model simulates and the table fills."""
    from PySide6.QtWidgets import QDockWidget, QToolBar

    with _Harness(say, ctx) as h:
        say("data:", h.data.describe())
        t0 = time.perf_counter()
        win = h.build()
        say(f"MainWindow built in {time.perf_counter() - t0:.1f} s")
        problems = []
        menus = [a.text().replace("&", "") for a in win.menuBar().actions()]
        say("menus:", ", ".join(menus))
        if menus != ["File", "Edit", "Process", "Fit", "Series", "Tools", "View",
                     "Plotting", "Help"]:
            problems.append(f"unexpected menu bar: {menus}")
        for name in ("view", "view2d", "central_stack", "cofit_page", "health_strip",
                     "lines_table", "ws_panel", "explorer", "proc_panel", "exp_label",
                     "pos_label"):
            if getattr(win, name, None) is None:
                problems.append(f"window member missing: {name}")
        docks = []
        for name in ("explorer_dock", "datasets_dock", "ws_dock", "lines_dock",
                     "results_dock", "proc_dock"):
            d = getattr(win, name, None)
            if not isinstance(d, QDockWidget):
                problems.append(f"dock missing: {name}")
            else:
                docks.append(d.windowTitle())
        say("docks:", ", ".join(docks))
        toolbars = [tb.windowTitle() for tb in win.findChildren(QToolBar) if tb.parent() is win]
        say("toolbars:", ", ".join(toolbars))
        if len(toolbars) < 2:
            problems.append(f"expected the main toolbar and the sidebar, found {toolbars}")
        if win.statusBar() is None or not win.statusBar().currentMessage():
            problems.append("no status bar message after construction")
        if win.central_stack.currentWidget() is not win.view:
            problems.append("the 1D view is not the initial central page")

        t0 = time.perf_counter()
        win.load_source(str(h.data.spectrum), keep_fit=False)
        pump(50)
        n_pts = len(win.exp_ppm)
        say(f"loaded {h.data.spectrum.name}: {n_pts} points, nucleus "
            f"{(win.recipe or {}).get('nucleus')}, {(win.recipe or {}).get('larmor_frequency_MHz')} "
            f"MHz in {time.perf_counter() - t0:.1f} s; experiment strip: {win.exp_label.text()!r}")
        if n_pts < 100 or not win.recipe:
            problems.append(f"the spectrum did not load ({n_pts} points, recipe {bool(win.recipe)})")
        if win.central_stack.currentWidget() is not win.view:
            problems.append("a 1D spectrum did not land on the 1D workbench")
        if not win.exp_label.text():
            problems.append("the experiment strip in the status bar is empty")

        t0 = time.perf_counter()
        win.apply_recipe(str(h.data.recipe))
        idle = h.settle()
        sites = (win.recipe or {}).get("sites") or []
        rows = win.lines_table.table.rowCount()
        model = win._last_model
        say(f"recipe applied: {len(sites)} line(s), {rows} table row(s), simulated in "
            f"{time.perf_counter() - t0:.1f} s" + ("" if idle else " (still busy)"))
        if not idle:
            problems.append("the simulation did not finish within the bound")
        if not sites or rows != len(sites):
            problems.append(f"lines table shows {rows} rows for {len(sites)} lines")
        if model is None:
            problems.append("no model curve after the simulation")
        else:
            x, total = model
            mx, my = win.view._model.getData()
            if len(x) < 10 or not np.all(np.isfinite(total)) or float(np.max(total)) <= 0:
                problems.append("the simulated model is empty or not finite")
            # the view draws the model inside the experiment's range (plus a
            # margin), so it holds a subset of the simulation grid
            if mx is None or len(mx) < 10 or not np.all(np.isfinite(np.asarray(my, float))):
                problems.append("the model curve on the plot does not carry the simulation")
            else:
                say(f"model curve: {len(mx)} of {len(x)} simulated points drawn, "
                    f"max {float(np.max(total)):.3g}")
        if win._health is None:
            problems.append("no fit-health verdict after the simulation")
        else:
            hv = win._health
            say(f"fit health: fitted={getattr(hv, 'fitted', '?')}, "
                f"{len(getattr(hv, 'flags', []))} flag(s), "
                f"{len(getattr(hv, 'unchecked', []))} fact(s) not judged")
        ok, detail = _paints(win.centralWidget())
        say("central area renders:", detail)
        if not ok:
            problems.append(f"the central area did not render ({detail})")
        if h.modals.seen:
            problems.append(f"unexpected dialogs while loading: {h.modals.seen}")
        _finish(h, problems)


def menu_sweep(say, ctx):
    """Every leaf action of the menu bar triggered -- with the spectrum and
    its fit on screen, then with no lines, then on the 2D map -- its dialogs
    dismissed, its workers waited for, its exceptions trapped."""
    with _Harness(say, ctx) as h:
        say("data:", h.data.describe())
        win = h.build()
        h.load()
        n_cmds = len(h.commands())
        say(f"{n_cmds} leaf entries under the menu bar (the pinned tree has "
            f"{EXPECTED_LEAF_ACTIONS} outside Open recent / Apply recipe)")
        passes = [("loaded: the spectrum with its fit on screen", h.restore_loaded),
                  ("empty: the spectrum with no lines", h.restore_empty)]
        if h.data.spectrum_2d is not None:
            passes.append(("2D: the MQMAS map on screen", h.restore_2d))
        else:
            say("no 2D example in this installation: the 2D pass is not run")
        problems = []
        for label, restore in passes:
            restore()
            t0 = time.perf_counter()
            records = h.sweep(label, restore)
            problems += _summarize(say, records)
            say(f"[{label}] done in {time.perf_counter() - t0:.0f} s")
        if win.actWatch.isChecked():
            win.actWatch.setChecked(False)
        _finish(h, problems)


def menu_themes(say, ctx):
    """Every theme of View ▸ Theme, every hidden style of More styles…, every
    text size, applied through the window's own setters with a spectrum and
    a fit on screen; the central area must render under each."""
    from PySide6.QtWidgets import QApplication

    from larmor.desktop import theme

    with _Harness(say, ctx) as h:
        win = h.build()
        h.load()
        app = QApplication.instance()
        font0 = app.font()
        theme0 = theme.active().name
        problems = []

        def check(kind, name):
            pump(50)
            if not h.settle():
                problems.append(f"{kind} {name}: the re-simulation did not finish")
            ok, detail = _paints(win.centralWidget())
            say(f"{kind} {name}: {detail}")
            if not ok:
                problems.append(f"{kind} {name}: the central area did not render ({detail})")
            errs = h.trap.take()
            for e in errs:
                problems.append(f"{kind} {name}: {e.strip().splitlines()[-1]}")
                say(e.rstrip())

        for name in theme.names():
            t = theme.get(name)
            win._set_theme(name)
            if theme.active().name != name:
                problems.append(f"theme {name}: active theme is {theme.active().name}")
            if theme.swatch_icon(t).isNull():
                problems.append(f"theme {name}: empty swatch icon")
            check("theme", name)
        for name in list(theme.aesthetic_names()) + [""]:
            win._set_aesthetic_override(name)
            check("style", name or "Normal")
        for pt in (8, 9, 11, 13):
            win._set_text_size(pt)
            if app.font().pointSize() != pt:
                problems.append(f"text size {pt} pt: application font is {app.font().pointSize()} pt")
            check("text size", f"{pt} pt")
        # back to where the sweep started (the settings guard restores the keys)
        app.setFont(font0)
        win._apply_theme_live(theme0)
        h.settle()
        say(f"{len(theme.names())} themes, {len(theme.aesthetic_names())} hidden styles "
            f"+ Normal, 4 text sizes")
        _finish(h, problems)


def _count_math(md: str) -> int:
    """How many ``$...$`` / ``$$...$$`` snippets mdrender will try to typeset."""
    from larmor.desktop import mdrender
    try:
        md, _blocks = mdrender._protect_code(md)
    except AttributeError:
        pass
    n_display = len(re.findall(r"\$\$(.+?)\$\$", md, flags=re.DOTALL))
    md = re.sub(r"\$\$(.+?)\$\$", " ", md, flags=re.DOTALL)
    n_inline = len(re.findall(r"(?<!\$)\$(?!\$)(.+?)(?<!\$)\$(?!\$)", md, flags=re.DOTALL))
    return n_display + n_inline


def _check_page(hw, kind: str, name: str) -> tuple[list[str], str]:
    """One rendered help window: non-trivial text, no 'not found', every
    equation typeset to an image, every linked page shipping."""
    from larmor.desktop.help_dialog import help_path, tutorial_path

    probs = []
    md, html = hw.markdown, hw.html
    text = hw.tb.toPlainText()
    if len(html) < 1000 or len(text.strip()) < 200:
        probs.append(f"{kind} {name}: only {len(text)} characters rendered")
    if "not found" in text[:200].lower():
        probs.append(f"{kind} {name}: the page was not found")
    heading = re.search(r"^#\s+(.+)$", md, flags=re.M)
    if heading:
        first = re.sub(r"[`*_$]", "", heading.group(1)).split()
        if first and first[0] not in text:
            probs.append(f"{kind} {name}: the title {heading.group(1)!r} is not in the rendered text")
    n_math = _count_math(md)
    n_img = html.count('<img src="data:image/png;base64,')
    if n_img != n_math:
        probs.append(f"{kind} {name}: {n_math} equations in the Markdown, {n_img} typeset images")
    links = 0
    for target in re.findall(r"\]\(([^)\s]+)\)", md):
        if target.startswith(("http://", "https://", "mailto:")):
            continue
        links += 1
        page = target.split("#")[0]
        if page.endswith(".md"):
            stem = Path(page).stem
            if help_path(stem) is None and tutorial_path(stem) is None:
                probs.append(f"{kind} {name}: links to {target}, which does not ship")
        elif page:
            probs.append(f"{kind} {name}: link target {target!r} is neither a page nor a URL")
    return probs, f"{len(text)} chars, {n_math} equations, {links} page links"


def menu_help(say, ctx):
    """Every manual and tutorial opened from the Help menu and rendered (text,
    typeset equations, links), plus the welcome hint, About, the command
    palette and the ? ▸ More card (offline)."""
    from PySide6.QtCore import QEvent

    from larmor.desktop.prefs import settings as prefs_settings
    from PySide6.QtWidgets import QApplication, QTextBrowser

    import larmor
    from larmor.desktop.help_dialog import help_path, tutorial_path
    from larmor.desktop.mw_menus import TUTORIALS

    with _Harness(say, ctx) as h:
        win = h.build()
        problems = []

        # the first-run hint on the empty canvas (a user with no recent files)
        prefs_settings().setValue("recent", [])
        win._maybe_show_welcome()
        pump(30)
        hint = getattr(win.view, "_placeholder", None)
        if hint is None or hint.isHidden() or "Open a spectrum" not in hint.text():
            problems.append("the first-run hint did not appear on the empty canvas")
        else:
            say("first-run hint shown:", hint.text().splitlines()[0])
        h.load()
        if hint is not None and not hint.isHidden():
            problems.append("the first-run hint stayed on after a spectrum loaded")

        def flush():
            pump(30)
            QApplication.sendPostedEvents(None, QEvent.DeferredDelete)
            pump(10)

        def take_page():
            reg = getattr(win, "_help_windows", {}) or {}
            if len(reg) != 1:
                return None, None
            (kind, name), hw = next(iter(reg.items()))
            return (kind, name), hw

        opened = {"manual": [], "tutorial": []}
        manuals = [c for c in h.commands() if tuple(c.path) == ("Help", "User manuals")]
        for cmd in manuals:
            cmd.action.trigger()
            pump(30)
            key, hw = take_page()
            if hw is None:
                problems.append(f"manual entry {cmd.label!r}: no help window registered")
                close_tool_windows(keep=[win])
                continue
            probs, detail = _check_page(hw, *key)
            problems += probs
            opened[key[0]].append(key[1])
            say(f"manual {key[1]} ({cmd.label}): {detail}")
            hw.close()
            flush()
        for stem, title in TUTORIALS:
            win._open_tutorial(stem, title)
            pump(30)
            key, hw = take_page()
            if hw is None:
                problems.append(f"tutorial {stem}: no help window registered")
                close_tool_windows(keep=[win])
                continue
            probs, detail = _check_page(hw, *key)
            problems += probs
            opened[key[0]].append(key[1])
            say(f"tutorial {stem}: {detail}")
            hw.close()
            flush()
        # every page that ships has a menu entry, and the other way round
        help_dir = help_path("getting-started")
        tut_dir = tutorial_path(TUTORIALS[0][0])
        if help_dir is None or tut_dir is None:
            problems.append("the help / tutorial folders are not where the viewer looks")
        else:
            shipped = {p.stem for p in help_dir.parent.glob("*.md")}
            if shipped != set(opened["manual"]):
                problems.append(f"manual pages vs menu entries differ: "
                                f"{sorted(shipped ^ set(opened['manual']))}")
            shipped_t = {p.stem for p in tut_dir.parent.glob("*.md")}
            if shipped_t != set(opened["tutorial"]):
                problems.append(f"tutorial pages vs TUTORIALS differ: "
                                f"{sorted(shipped_t ^ set(opened['tutorial']))}")
        say(f"{len(opened['manual'])} manuals and {len(opened['tutorial'])} tutorials rendered")
        problems += [e.strip().splitlines()[-1] for e in h.trap.take()]

        # About
        win._about()
        pump(50)
        about = [w for w in QApplication.topLevelWidgets()
                 if w.windowTitle() == "About LARMOR" and w.isVisible()]
        if not about:
            problems.append("About LARMOR did not open")
        else:
            tb = about[0].findChild(QTextBrowser)
            text = tb.toPlainText() if tb is not None else ""
            if "LARMOR" not in text or larmor.__version__ not in text:
                problems.append("the About text lacks the name or the version")
            say(f"About LARMOR: v{larmor.__version__}, {len(text)} chars")
            about[0].close()
            flush()

        # the command palette: the live list, a query, a miss; then the real
        # modal through the menu slot (dismissed like any dialog)
        from larmor.desktop.palette_dialog import PLACEHOLDER, CommandPalette, collect_commands
        cmds = collect_commands(win)
        if len(cmds) < EXPECTED_LEAF_ACTIONS:
            problems.append(f"the command palette lists {len(cmds)} commands, fewer than "
                            f"the {EXPECTED_LEAF_ACTIONS} menu leaves")
        dlg = CommandPalette(win, cmds)
        dlg.edit.setText("czjzek")
        pump(20)
        n_hits = dlg.tree.topLevelItemCount()
        top = dlg.tree.topLevelItem(0).text(0) if n_hits else ""
        if n_hits < 1 or "zjzek" not in top:
            problems.append(f"palette query 'czjzek': {n_hits} rows, first {top!r}")
        dlg.edit.setText("qqqzzzxxx")
        pump(20)
        miss = dlg.tree.topLevelItem(0).text(0) if dlg.tree.topLevelItemCount() else ""
        if miss != PLACEHOLDER:
            problems.append(f"palette miss shows {miss!r}")
        dlg.close()
        dlg.deleteLater()
        flush()
        n_modal = len(h.modals.seen)
        chosen = win.open_command_palette()
        pump(50)
        if chosen is not None or "Command palette" not in h.modals.seen[n_modal:]:
            problems.append(f"the palette modal was not dismissed ({h.modals.seen[n_modal:]})")
        say(f"command palette: {len(cmds)} commands, 'czjzek' -> {n_hits} rows "
            f"(first: {top}), miss -> placeholder, modal dismissed")

        # ? ▸ More: the card on its offline pack (http stubbed to None)
        win._show_more()
        card = getattr(win, "_more_dlg", None)
        if card is None:
            problems.append("? ▸ More opened no card")
        else:
            h.wait_more_card()
            pump(50)
            name = card.name_lab.text() if hasattr(card, "name_lab") else ""
            if not name:
                problems.append("the More card showed no fact (offline pack)")
            say(f"More card (offline): {name!r}")
            card.close()
            flush()
        _finish(h, problems)


def menu_panels(say, ctx):
    """Every checkable View entry flipped twice (the plain toggles, the
    docks, the fit-health strip), every axis unit, every Y-axis mode and
    every Czjzek width convention switched through -- the window must still
    render and each state must take."""
    from larmor.desktop import table as _table

    with _Harness(say, ctx) as h:
        win = h.build()
        h.load()
        problems = []

        def rendered(what):
            pump(30)
            h.settle()
            ok, detail = _paints(win.centralWidget())
            if not ok:
                problems.append(f"{what}: the central area did not render ({detail})")
            for e in h.trap.take():
                problems.append(f"{what}: {e.strip().splitlines()[-1]}")
                say(e.rstrip())

        toggles = [("Residual", win.actResid), ("Components", win.actComp),
                   ("Component labels", win.actLabels), ("Paddles", win.actPaddles),
                   ("Overlays", win.actOverlaysVisible),
                   ("Literature shift ranges", win.actRefRanges),
                   ("Scroll nudges fit values", win.actScrollNudge),
                   ("Fit health strip", win.actHealthStrip)]
        for label, act in toggles:
            before = act.isChecked()
            act.trigger()
            rendered(f"View ▸ {label} {'off' if before else 'on'}")
            if act.isChecked() == before:
                problems.append(f"View ▸ {label}: trigger() did not flip the check")
            act.trigger()
            rendered(f"View ▸ {label} back")
            if act.isChecked() != before:
                problems.append(f"View ▸ {label}: not back to its initial state")
        say(f"{len(toggles)} view toggles flipped twice")

        docks = [win.explorer_dock, win.datasets_dock, win.ws_dock, win.lines_dock,
                 win.results_dock, win.proc_dock]
        for d in docks:
            act = d.toggleViewAction()
            hidden0 = d.isHidden()
            act.trigger()
            rendered(f"View ▸ Panels ▸ {d.windowTitle()}")
            if d.isHidden() == hidden0:
                problems.append(f"Panels ▸ {d.windowTitle()}: the dock did not toggle")
            act.trigger()
            rendered(f"View ▸ Panels ▸ {d.windowTitle()} back")
            if d.isHidden() != hidden0:
                problems.append(f"Panels ▸ {d.windowTitle()}: not back to its initial state")
        say(f"{len(docks)} docks toggled twice: "
            + ", ".join(f"{d.windowTitle()} ({'hidden' if d.isHidden() else 'shown'})"
                        for d in docks))

        unit0 = win._axis_unit
        for unit, act in win._axis_unit_actions.items():
            act.trigger()
            rendered(f"Axis unit {unit}")
            if win._axis_unit != unit or not act.isChecked():
                problems.append(f"Axis unit {unit}: the window reports {win._axis_unit}")
        win._axis_unit_actions[unit0].trigger()
        rendered("Axis unit back")
        say(f"axis units: {', '.join(win._axis_unit_actions)}; readout "
            f"{win._format_x(float(win.exp_ppm[len(win.exp_ppm) // 2]))!r}")

        mode0 = win.view.y_mode()[0]
        n_modal = len(h.modals.seen)
        for key, act in win._y_axis_actions.items():
            act.trigger()
            rendered(f"Y axis {key}")
            mode = win.view.y_mode()[0]
            if key == "region":
                # its dialog is dismissed: the mode stays, the radio re-syncs
                if mode == "region" or act.isChecked():
                    problems.append("Y axis region: a cancelled dialog changed the mode")
            elif mode != key or not act.isChecked():
                problems.append(f"Y axis {key}: the view reports {mode}")
        win._y_axis_actions[mode0].trigger()
        rendered("Y axis back")
        say(f"Y-axis modes: {', '.join(win._y_axis_actions)}; dialogs "
            f"{h.modals.seen[n_modal:]}")

        cz0 = _table.czjzek_display_mode()
        keys = list(_table.CZJZEK_DISPLAYS)
        acts = win._czjzek_display_group.actions()
        if len(acts) != len(keys):
            problems.append(f"Czjzek width display: {len(acts)} entries for {len(keys)} conventions")
        for key, act in zip(keys, acts):
            act.trigger()
            rendered(f"Czjzek width display {key}")
            if _table.czjzek_display_mode() != key:
                problems.append(f"Czjzek width display {key}: table module reports "
                                f"{_table.czjzek_display_mode()}")
            t = win.lines_table.table
            headers = [t.horizontalHeaderItem(i).text() for i in range(t.columnCount())
                       if t.horizontalHeaderItem(i) is not None]
            label = _table.CZJZEK_DISPLAYS[key][0]
            if (win.recipe and any(s.get("model", "").startswith("czjzek")
                                   for s in win.recipe.get("sites", []))
                    and not any(label in hd for hd in headers)):
                problems.append(f"Czjzek width display {key}: no table header reads {label!r}")
        acts[keys.index(cz0)].trigger()
        rendered("Czjzek width display back")
        say(f"Czjzek width conventions: {', '.join(keys)}")
        _finish(h, problems)


def stages() -> list:
    return [("menu-window", menu_window), ("menu-sweep", menu_sweep),
            ("menu-themes", menu_themes), ("menu-help", menu_help),
            ("menu-panels", menu_panels)]
