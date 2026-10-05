"""Desktop sweep, part 2: the main workflows of the desktop run end to end,
offscreen, on the bundled example data (synthetic CSV spectra when the
examples are not shipped) -- load, apply a recipe, fit, every export
format, the Plotting studio and the live plot's exporters, overlays, the
series mode, the 2D view, the processing panel, project save / open, the
Explorer and Workspaces docks. Each stage is ``(name, fn(say, ctx))``;
plumbing in ``gui_common``.

Every stage builds its own ``MainWindow`` (never shown), drives the same
methods the menus and buttons call, and tears the window down with its
worker threads joined. Where a workflow asks for a file path or a Yes, the
static dialog method is replaced for the duration of that call
(``_patched``) and the path points into the scratch folder; a dialog whose
result the workflow reads back (``exec``) is replaced by a stand-in that
records the dialog and answers at once; any OTHER modal the sweep provokes
is rejected by ``ModalDismisser`` and listed in the log. Exceptions raised
inside Qt slots are caught by ``ExceptionTrap`` and fail the stage; a
worker thread still running after the window closed fails it too. The
user's saved settings are restored by ``SettingsGuard``.

The bug this sweep exists for: the frozen 0.15.3 exe raised "No module
named 'matplotlib.backends.backend_svg'" on the Plotting studio's SVG
export -- nothing had ever exported a figure from the frozen build.
"""
from __future__ import annotations

import contextlib
import csv
import json
import os
import time
from pathlib import Path

import numpy as np

from larmor.distcheck.gui_common import (ExceptionTrap, ModalDismisser, SettingsGuard,
                                          close_tool_windows, ensure_app, pump)

__all__ = ["stages"]

#: the bundled example data under examples/pCABS2-4: the 1D 27Al MAS, the
#: 1D 11B MAS and the 27Al 3QMAS 2D, plus the two shipped recipes and the
#: two-panel figure spec next to them
_EXPNO_27AL, _EXPNO_11B, _EXPNO_2D = "3616", "1118", "3620"

#: the window's worker threads (MainWindow.closeEvent joins them)
_WORKERS = ("_warm_worker", "_sim_worker", "_fit_worker", "_fit2d_worker", "_seq_worker")

_FIT_TIMEOUT = 90.0          # a fit / sweep / kernel build, bounded

_MAGIC = {"png": (b"\x89PNG",), "pdf": (b"%PDF",), "eps": (b"%!PS",),
          "tiff": (b"II*\x00", b"MM\x00*"), "jpg": (b"\xff\xd8",), "svg": (b"<svg",)}


# ====================================================================== plumbing
def _examples(ctx) -> dict | None:
    """The bundled example files, or None when any of them is missing."""
    ex = ctx.get("examples")
    if not ex:
        return None
    ex = Path(ex)
    s = ex / "pCABS2-4"
    d = {"root": ex, "sample": s,
         "al_expno": s / _EXPNO_27AL, "al_1r": s / _EXPNO_27AL / "pdata" / "1" / "1r",
         "al_fid": s / _EXPNO_27AL / "fid",
         "b_expno": s / _EXPNO_11B, "b_1r": s / _EXPNO_11B / "pdata" / "1" / "1r",
         "twod_expno": s / _EXPNO_2D, "twod_2rr": s / _EXPNO_2D / "pdata" / "1" / "2rr",
         "recipe_al": ex / "pCABS2-4_27Al.recipe.json",
         "recipe_b": ex / "pCABS2-4_11B.recipe.json",
         "figure": ex / "pCABS2-4_fits.figure.json"}
    return d if all(p.exists() for p in d.values()) else None


def _write_synthetic(out: Path, name: str, pos: float, amp: float, nucleus: str = "11B",
                     larmor_MHz: float = 160.0) -> str:
    """One gauss_lor line (6 ppm wide, a little noise) at ``pos`` on a
    -20..60 ppm grid, written as a LARMOR CSV spectrum. Returns the path."""
    from larmor import engine
    from larmor.io import spectra
    from larmor.recipe import Param, Recipe, SiteModel

    x = np.linspace(-20.0, 60.0, 500)
    rec = Recipe(nucleus=nucleus, larmor_frequency_MHz=larmor_MHz, spin_rate_Hz=0.0,
                 sites=[SiteModel(model="gauss_lor", label="A", params={
                     "isotropic_chemical_shift_ppm": Param(pos),
                     "shift_fwhm_ppm": Param(6.0), "amplitude": Param(amp),
                     "gl": Param(1.0, vary=False)})])
    _, model, _ = engine.simulate(rec, exp_ppm=x)
    y = model + 0.002 * amp * np.random.default_rng(7).standard_normal(x.size)
    out.mkdir(parents=True, exist_ok=True)
    return str(spectra.write_csv(out / name, x, y, {"nucleus": nucleus, "larmor_MHz": larmor_MHz,
                                                    "sample": Path(name).stem}))


def _synthetic_series(out: Path) -> tuple[list, dict]:
    """Three 11B spectra at 13 / 15 / 17 ppm with different heights, and the
    one-line model a user would have on screen (the series tests' fixture)."""
    paths = [_write_synthetic(out, f"s{k}.csv", pos, amp)
             for k, (pos, amp) in enumerate(zip((13.0, 15.0, 17.0), (100.0, 50.0, 80.0)))]
    model = {"nucleus": "11B", "larmor_frequency_MHz": 160.0, "spin_rate_Hz": 0.0,
             "sites": [{"model": "gauss_lor", "label": "A", "params": {
                 "isotropic_chemical_shift_ppm": {"value": 12.0, "min": 0, "max": 30},
                 "shift_fwhm_ppm": {"value": 5.0, "min": 0.1},
                 "amplitude": {"value": 80.0, "min": 0},
                 "gl": {"value": 1.0, "vary": False}}}]}
    return paths, model


def _new_window(wins: list):
    """A MainWindow the way the UI tests build one: never shown, its three
    user-facing prompts stubbed (Open recent stays the user's, project notes
    go nowhere, Open project replaces)."""
    from larmor.desktop.app import MainWindow

    win = MainWindow()
    win._add_recent = lambda p: None
    win._report_project_notes = lambda notes: None
    win._confirm_open_mode = lambda: "replace"
    wins.append(win)
    pump(50)
    return win


def _running(w) -> bool:
    try:
        return w is not None and hasattr(w, "isRunning") and w.isRunning()
    except RuntimeError:                       # deleted on the C++ side
        return False


def _wait(cond, timeout: float, what: str) -> None:
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance()
    t0 = time.perf_counter()
    while not cond():
        if time.perf_counter() - t0 > timeout:
            raise AssertionError(f"timed out after {timeout:.0f} s waiting for {what}")
        app.processEvents()
        time.sleep(0.01)
    app.processEvents()


def _settle(win, timeout: float = 60.0) -> None:
    """Let the debounced simulation land and its thread finish (a first
    Czjzek / amorphous simulation builds its kernel: seconds)."""
    def idle():
        sw = win._sim_worker
        return not (win._sim_timer.isActive() or win._sim_pending
                    or (sw is not None and _running(sw)))
    _wait(idle, timeout, "the simulation to settle")


def _fit_idle(win) -> bool:
    return not _running(win._fit_worker) and win._active_fit_worker is None


def _wait_fit(win, what: str = "the fit", timeout: float = _FIT_TIMEOUT) -> None:
    _wait(lambda: _fit_idle(win), timeout, what)
    _settle(win)


def _status(win) -> str:
    return win.statusBar().currentMessage()


def _teardown(win) -> list:
    """Close the window the way the user would, every tool window first;
    returns the names of the worker threads still running afterwards."""
    close_tool_windows(keep=[win])
    bar = getattr(win, "_window_bar", None)
    if bar is not None:
        try:
            bar.close_all()
        except Exception:                      # noqa: BLE001
            pass
    try:
        _settle(win, 30.0)
    except AssertionError:
        pass
    win.close()                                # closeEvent joins the workers
    pump(100)
    return [n for n in _WORKERS if _running(getattr(win, n, None))]


def _run(body, say, ctx) -> None:
    """One stage: the app, the settings guard, the exception trap and the
    modal dismisser around ``body(say, ctx, wins, dismiss)``; every window
    the body registered in ``wins`` is torn down afterwards, whatever
    happened, and the trap's catch fails the stage."""
    ensure_app()
    wins: list = []
    with SettingsGuard(), ExceptionTrap() as trap, ModalDismisser() as dismiss:
        err = None
        left: list = []
        try:
            body(say, ctx, wins, dismiss)
        except BaseException as exc:           # noqa: BLE001 -- reported below
            err = exc
        finally:
            for w in wins:
                try:
                    left += _teardown(w)
                except Exception as exc:       # noqa: BLE001
                    left.append(f"teardown failed: {exc!r}")
            close_tool_windows()
            for w in wins:
                try:
                    w.deleteLater()
                except RuntimeError:
                    pass
            wins.clear()
            pump(150)
        if dismiss.seen:
            say("modal dialogs auto-dismissed:", ", ".join(dismiss.seen))
        notes = []
        if trap.errors:
            notes.append(f"{len(trap.errors)} exception(s) escaped Qt slots:\n"
                         + "\n".join(trap.errors))
        if left:
            notes.append(f"worker thread(s) still running after close: {left}")
        # a style sheet Qt refuses is a widget that silently lost its look
        # (the fit table's held-value italics, found this way)
        bad_css = [w for w in trap.warnings if "Could not parse stylesheet" in w]
        if bad_css:
            notes.append(f"{len(bad_css)} style sheet(s) Qt could not parse: {bad_css[0][:120]!r}")
        if err is not None:
            if notes:
                raise AssertionError(f"{err!r}; also: " + "; ".join(notes)) from err
            raise err
        if notes:
            raise AssertionError("; ".join(notes))
        if trap.warnings:
            say(f"{len(trap.warnings)} Qt warning(s), first: {trap.warnings[0][:160]!r}")


@contextlib.contextmanager
def _patched(target, name: str, value):
    """``target.name = value`` for the block, the original put back after
    (an attribute the class only inherited is removed again)."""
    own = name in vars(target)
    old = getattr(target, name)
    setattr(target, name, value)
    try:
        yield
    finally:
        if own:
            setattr(target, name, old)
        else:
            delattr(target, name)


def _save_dialog(path, chosen_filter: str = ""):
    from PySide6.QtWidgets import QFileDialog
    return _patched(QFileDialog, "getSaveFileName",
                    staticmethod(lambda *a, **k: (str(path), chosen_filter)))


def _open_dialog(path):
    from PySide6.QtWidgets import QFileDialog
    return _patched(QFileDialog, "getOpenFileName",
                    staticmethod(lambda *a, **k: (str(path), "")))


def _dir_dialog(path):
    from PySide6.QtWidgets import QFileDialog
    Path(path).mkdir(parents=True, exist_ok=True)
    return _patched(QFileDialog, "getExistingDirectory",
                    staticmethod(lambda *a, **k: str(path)))


def _answer_yes():
    from PySide6.QtWidgets import QMessageBox
    return _patched(QMessageBox, "question",
                    staticmethod(lambda *a, **k: QMessageBox.Yes))


def _export_choice(fmt: str, dpi: int = 150):
    """The shared Export figure dialog answered with ``fmt`` at 12 x 9 cm."""
    from larmor.desktop import export_dialog
    return _patched(export_dialog, "choose",
                    lambda *a, **k: {"format": fmt, "dpi": dpi,
                                     "width_cm": 12.0, "height_cm": 9.0})


def _export_path(path):
    """The exporters' Save figure file dialog answered with ``path``."""
    from larmor.desktop import export_dialog
    return _patched(export_dialog, "_ask_path", lambda *a, **k: str(path))


def _exec_captured(cls, created: list, result: int = 0, before=None):
    """``cls.exec`` replaced: the dialog is recorded in ``created``, ``before``
    (if any) acts on it, and ``result`` comes back at once -- the way to
    drive a modal dialog's widgets from here, or to accept it with its
    defaults, without a nested event loop."""
    def fake_exec(self, *a, **k):
        created.append(self)
        if before is not None:
            before(self)
        return result
    return _patched(cls, "exec", fake_exec)


def _check_file(path, kind: str | None = None, min_bytes: int = 200) -> int:
    """The file exists, is not suspiciously small and starts the way its
    format must (``_MAGIC``); returns its size."""
    p = Path(path)
    if not p.is_file():
        raise AssertionError(f"{p.name} was not written")
    size = p.stat().st_size
    if size < min_bytes:
        raise AssertionError(f"{p.name} is suspiciously small ({size} bytes)")
    if kind:
        head = p.read_bytes()[:4096]
        magic = _MAGIC[kind]
        ok = (b"<svg" in head) if kind == "svg" else any(head.startswith(m) for m in magic)
        if not ok:
            raise AssertionError(f"{p.name} does not look like {kind.upper()}: {head[:24]!r}")
    return size


def _clip_text() -> str:
    from PySide6.QtWidgets import QApplication
    return QApplication.clipboard().text()


def _dump(recipe) -> str:
    return json.dumps(recipe, sort_keys=True, default=str)


def _absolute_figure_spec(ex: dict, out: Path) -> Path:
    """The shipped two-panel figure spec, its recipes copied into ``out``
    with ABSOLUTE source paths: the originals name their data relative to
    the repository root, which is not where a user's shell stands."""
    spec = json.loads(ex["figure"].read_text(encoding="utf-8"))
    out.mkdir(parents=True, exist_ok=True)
    for panel in spec.get("panels", []):
        name = Path(panel["recipe"]).name
        rec = json.loads((ex["root"] / name).read_text(encoding="utf-8"))
        sp = Path(str(rec.get("source_path") or ""))
        if not sp.is_absolute() and sp.parts and sp.parts[0] == "examples":
            sp = ex["root"].joinpath(*sp.parts[1:])
        rec["source_path"] = str(sp.resolve())
        copy = out / name
        copy.write_text(json.dumps(rec, indent=1), encoding="utf-8")
        panel["recipe"] = str(copy)
    target = out / ex["figure"].name
    target.write_text(json.dumps(spec, indent=1), encoding="utf-8")
    return target


def _load_1d(win, path, what: str, nucleus: str | None = None) -> None:
    """``File ▸ Open`` on ``path``; it must land on the 1D workbench."""
    win.load_source(str(path), keep_fit=False)
    pump(50)
    _settle(win)
    if win.central_stack.currentWidget() is not win.view or not win.exp_ppm.size:
        raise AssertionError(f"{what} did not open on the 1D workbench: {_status(win)!r}")
    if nucleus and (win.recipe or {}).get("nucleus") != nucleus:
        raise AssertionError(f"{what} opened as {(win.recipe or {}).get('nucleus')!r}, "
                             f"expected {nucleus}")


# ====================================================================== stages
# ------------------------------------------------------------ fit workflow
def _fit_workflow(say, ctx, wins, dismiss):
    """Open the example 27Al, apply the shipped recipe, Fit (the FitWorker
    thread), the verdict and the quantify table, the Report dock's three
    clipboard buttons, Save recipe, Save fit as… in every format, the
    publication bundle, Save spectrum, undo / redo / New fit."""
    from PySide6.QtWidgets import QFileDialog

    from larmor.io import export as io_export
    from larmor.io import spectra
    from larmor.recipe import Recipe
    from larmor.seriesmode import rmsd_of

    out = Path(ctx["out"]) / "fit"
    out.mkdir(parents=True, exist_ok=True)
    ex = _examples(ctx)
    win = _new_window(wins)
    if ex:
        src, recipe_path, nucleus = ex["al_1r"], ex["recipe_al"], "27Al"
        what = f"{_EXPNO_27AL}/pdata/1/1r"
    else:
        say("no bundled examples: a synthetic 11B spectrum and a one-line model stand in")
        paths, model = _synthetic_series(out / "synthetic")
        src, nucleus, what = paths[0], "11B", "s0.csv"
        recipe_path = out / "synthetic" / "model.recipe.json"
        Recipe.from_dict(model).save(recipe_path)
    t0 = time.perf_counter()
    _load_1d(win, src, what, nucleus)
    say(f"opened {what}: {win.exp_ppm.size} points, {nucleus}, {time.perf_counter() - t0:.1f} s")

    # Fit ▸ Apply recipe: the shipped lines onto the data; the simulation
    # that follows builds the lineshape kernel (once per process)
    t0 = time.perf_counter()
    win.apply_recipe(str(recipe_path))
    _settle(win, _FIT_TIMEOUT)
    n_sites = len(win.recipe["sites"])
    if not n_sites:
        raise AssertionError(f"apply_recipe put no lines on the spectrum: {_status(win)!r}")
    if win._last_model is None:
        raise AssertionError("the model was not simulated after apply_recipe")
    say(f"applied {Path(recipe_path).name}: {n_sites} line(s) "
        f"({', '.join(s['model'] for s in win.recipe['sites'])}) simulated in "
        f"{time.perf_counter() - t0:.1f} s")

    # Fit (F5): the worker thread, then the verdict
    before = _dump(win.recipe)
    t0 = time.perf_counter()
    win.run_fit()
    if win._fit_worker is None:
        raise AssertionError(f"run_fit started no worker: {_status(win)!r}")
    _wait_fit(win)
    msg = _status(win)
    if "fit done" not in msg:
        raise AssertionError(f"the fit did not finish normally: {msg!r}"
                             + (f" (dismissed: {dismiss.seen})" if dismiss.seen else ""))
    rmsd = rmsd_of(win.recipe)
    if rmsd is None or not np.isfinite(rmsd):
        raise AssertionError("the fitted recipe carries no fit_rmsd")
    h = win._health
    if h is None or not h.fitted:
        raise AssertionError("no fit verdict after the fit")
    if not win.health_strip.isVisibleTo(win):
        raise AssertionError("the fit-health strip stayed hidden after the fit")
    rows = win.qtable.rowCount()
    if rows < n_sites:
        raise AssertionError(f"the quantify table has {rows} rows for {n_sites} lines")
    if not win.lines_table.chi2.text().strip():
        raise AssertionError("the fit table's chi2 label is empty after the fit")
    if not win.report.toPlainText().strip():
        raise AssertionError("the Report is empty after the fit")
    after = _dump(win.recipe)
    say(f"fit done in {time.perf_counter() - t0:.1f} s: RMSD {rmsd:.4g}, verdict "
        f"{win.health_strip.pill.text()!r}, {rows} quantify rows, chi2 "
        f"{win.lines_table.chi2.text()!r}")

    # the Report dock's clipboard buttons
    win.copy_csv()
    text = _clip_text()
    if not text.startswith("line,model,position_ppm") or len(text.splitlines()) < n_sites + 1:
        raise AssertionError(f"Copy CSV put an unexpected text on the clipboard: {text[:120]!r}")
    win.copy_latex()
    tex = _clip_text()
    if "tabular" not in tex:
        raise AssertionError(f"Copy LaTeX put no table on the clipboard: {tex[:120]!r}")
    win.copy_methods()
    meth = _clip_text()
    if len(meth) < 40:
        raise AssertionError(f"Copy methods put nothing useful on the clipboard: {meth!r}")
    say(f"clipboard: CSV {len(text)} chars, LaTeX {len(tex)} chars, methods {len(meth)} chars")

    # File ▸ Save recipe…, then Save fit as… in every format the dialog offers
    rec_path = out / "fit.recipe.json"
    with _save_dialog(rec_path):
        win.save_recipe()
    _check_file(rec_path, min_bytes=100)
    saved = Recipe.load(rec_path)
    if len(saved.sites) != n_sites or rmsd_of(saved) is None:
        raise AssertionError("the saved recipe lost its lines or its RMSD")
    sizes = {}
    for name, (ext, _fn) in io_export.FORMATS.items():
        target = out / f"fit.{ext}"
        with _patched(QFileDialog, "getSaveFileName",
                      staticmethod(lambda *a, _t=target, _n=name, **k: (str(_t), _n))):
            win.save_fit_as()
        if "fit exported" not in _status(win):
            raise AssertionError(f"Save fit as {name}: {_status(win)!r}")
        sizes[ext] = _check_file(target, min_bytes=40)
    with open(out / "fit.json", encoding="utf-8") as f:
        if len(json.load(f).get("sites", [])) != n_sites:
            raise AssertionError("fit.json does not carry the lines")
    if not (out / "fit.csv").read_text(encoding="utf-8").splitlines()[0]:
        raise AssertionError("fit.csv has no header line")
    if not (out / "fit.fxmla").read_text(encoding="utf-8", errors="replace").lstrip().startswith("<"):
        raise AssertionError("fit.fxmla is not XML")
    say("Save recipe + Save fit as: " + ", ".join(f"{e} {s} B" for e, s in sizes.items()))

    # the Report dock's Publication bundle… (figure in three formats, table,
    # methods, report.md) into a folder
    bundle = out / "bundle"
    with _dir_dialog(bundle):
        win.export_publication_bundle()
    if "publication bundle written" not in _status(win):
        raise AssertionError(f"publication bundle: {_status(win)!r}")
    for fname, kind in (("figure.png", "png"), ("figure.pdf", "pdf"), ("figure.svg", "svg"),
                        ("methods.txt", None), ("table.tex", None), ("report.md", None)):
        _check_file(bundle / fname, kind, min_bytes=40)
    say("publication bundle: " + ", ".join(sorted(p.name for p in bundle.iterdir())))

    # File ▸ Export ▸ Save spectrum as CSV…, and it reads back
    spec_csv = out / "spectrum.csv"
    with _save_dialog(spec_csv):
        win.save_spectrum()
    _check_file(spec_csv)
    ppm, amp, meta = spectra.read_csv(spec_csv)
    if ppm.size != win.exp_ppm.size or meta.get("nucleus") != nucleus:
        raise AssertionError(f"the saved spectrum reads back as {ppm.size} points, {meta}")
    say(f"spectrum saved and read back: {ppm.size} points, {meta.get('nucleus')}")

    # Edit ▸ Undo / Redo around the fit, New fit, and undo again
    if not win.undo_stack:
        raise AssertionError("the fit left nothing on the undo stack")
    win.undo()
    _settle(win)
    if _dump(win.recipe) != before:
        raise AssertionError("undo did not restore the pre-fit recipe")
    win.redo()
    _settle(win)
    if _dump(win.recipe) != after:
        raise AssertionError("redo did not bring the fit back")
    win.new_fit()
    _settle(win)
    if win.recipe["sites"]:
        raise AssertionError("New fit left lines behind")
    win.undo()
    _settle(win)
    if len(win.recipe["sites"]) != n_sites:
        raise AssertionError("undo after New fit did not bring the lines back")
    say("undo / redo / New fit / undo: the lines go and come back exactly")


# ------------------------------------------------------------ plotting studio
def _plotting_studio(say, ctx, wins, dismiss):
    """Plotting ▸ Plot current spectrum seeds the studio; export in every
    format the Export dialog offers (the frozen-build bug), every style
    preset, every plot kind and template at least previews, Save spec… /
    Load spec… round-trip, the shipped two-panel figure spec loads."""
    from larmor import figures
    from larmor.desktop import export_dialog
    from larmor.desktop.plotting_studio import PlottingStudio

    out = Path(ctx["out"]) / "studio"
    out.mkdir(parents=True, exist_ok=True)
    ex = _examples(ctx)
    win = _new_window(wins)
    if ex:
        _load_1d(win, ex["al_1r"], f"{_EXPNO_27AL}/pdata/1/1r", "27Al")
    else:
        say("no bundled examples: a synthetic 11B spectrum stands in")
        _load_1d(win, _write_synthetic(out, "s0.csv", 13.0, 100.0), "s0.csv")
    created: list = []
    with _exec_captured(PlottingStudio, created):
        win.plot_current_spectrum()
    if not created:
        raise AssertionError("Plot current spectrum opened no Plotting studio")
    st = created[-1]
    if len(st._traces) != 1 or len(st._traces[0]["data"]["x"]) != win.exp_ppm.size:
        raise AssertionError(f"the studio was seeded with {len(st._traces)} trace(s)")
    st._refresh()
    if st.msg.text() or not st._canvas.figure.axes:
        raise AssertionError(f"the preview did not render: {st.msg.text()!r}")
    say(f"studio seeded with the spectrum ({win.exp_ppm.size} points); preview rendered")

    # Export… in every format the shared dialog offers
    ext_of = {**export_dialog._RASTER, **export_dialog._VECTOR}
    sizes = {}
    for fmt, ext in ext_of.items():
        target = out / f"studio.{ext}"
        with _export_choice(fmt), _export_path(target):
            st._export()
        if not st.msg.text().startswith("exported"):
            raise AssertionError(f"export {fmt}: {st.msg.text()!r}")
        sizes[ext] = _check_file(target, kind=ext)
    say("exported: " + ", ".join(f"{e} {s} B" for e, s in sizes.items()))

    # every style preset renders the same figure
    for name in figures.STYLES:
        st.style.setCurrentText(name)
        st._refresh()
        if st.msg.text():
            raise AssertionError(f"style {name!r}: {st.msg.text()!r}")
    st.style.setCurrentText("article")
    say(f"{len(figures.STYLES)} style presets render: {', '.join(figures.STYLES)}")

    # every plot kind: 1D (above), 2D contour on the example map, the batch
    # grid from the shipped figure spec (Load spec…), a species bar from a
    # small table; the Series kind needs a relaxation series nothing ships
    kinds = ["1d"]
    if ex:
        st.kind.setCurrentIndex(1)
        st.path2d.setText(str(ex["twod_expno"]))
        st.nuc2d.setText("27Al")
        st._refresh()
        if st.msg.text():
            raise AssertionError(f"2D contour of {_EXPNO_2D}: {st.msg.text()!r}")
        kinds.append("2d")
        spec_path = _absolute_figure_spec(ex, out / "shipped")
        t0 = time.perf_counter()
        with _open_dialog(spec_path):
            st._load_spec()
        if st.kind.currentIndex() != 3 or len(st._panels) != 2:
            raise AssertionError(f"{ex['figure'].name} did not load as a two-panel batch grid")
        if st.msg.text():
            raise AssertionError(f"batch grid from {ex['figure'].name}: {st.msg.text()!r}")
        kinds.append("batch_grid")
        say(f"{ex['figure'].name} loaded through Load spec… and rendered "
            f"({time.perf_counter() - t0:.1f} s, both recipes simulated on their data)")
    else:
        for k in (1, 3):
            st.kind.setCurrentIndex(k)
            if "cannot render" in st.msg.text():
                raise AssertionError(f"kind {st.kind.currentText()!r}: {st.msg.text()!r}")
        say("no bundled examples: the 2D contour and batch grid kinds show their empty-spec hints")
    st._apply_spec({"kind": "species_bar", "categories": ["a", "b"],
                    "series": [{"label": "s0", "values": [60.0, 40.0]},
                               {"label": "s1", "values": [40.0, 60.0]}]})
    st._kind_changed()
    if st.msg.text():
        raise AssertionError(f"species bar: {st.msg.text()!r}")
    kinds.append("species_bar")
    st.kind.setCurrentIndex(2)
    if "cannot render" in st.msg.text():
        raise AssertionError(f"Series kind: {st.msg.text()!r}")
    say(f"plot kinds rendered: {', '.join(kinds)}; Series shows its hint "
        f"({st.msg.text()!r}) -- no relaxation series is bundled")

    # every template applies (kind + defaults) and previews on what is loaded
    names = []
    for i in range(1, st.template.count()):
        name = st.template.itemData(i)
        st.template.setCurrentIndex(i)
        if "cannot render" in st.msg.text():
            raise AssertionError(f"template {name!r}: {st.msg.text()!r}")
        names.append(name)
    st.template.setCurrentIndex(0)
    say(f"{len(names)} templates applied: {', '.join(names)}")

    # Save spec… / Load spec… round trip into a fresh studio
    st.kind.setCurrentIndex(0)
    st.title.setText("round trip")
    spec1 = st._spec()
    spec_file = out / "studio_spec.json"
    with _save_dialog(spec_file):
        st._save_spec()
    _check_file(spec_file, min_bytes=100)
    st2 = PlottingStudio(win)
    with _open_dialog(spec_file):
        st2._load_spec()
    spec2 = st2._spec()
    if (spec2.get("kind") != spec1.get("kind") or spec2.get("title") != "round trip"
            or len(spec2.get("traces", [])) != len(spec1.get("traces", []))
            or spec2.get("style") != spec1.get("style")):
        raise AssertionError("the spec did not survive Save spec… / Load spec…")
    if st2.msg.text():
        raise AssertionError(f"the reloaded spec did not render: {st2.msg.text()!r}")
    say(f"Save spec… / Load spec… round trip: {spec_file.stat().st_size} B, "
        f"{len(spec2['traces'])} trace(s), style {spec2.get('style')!r}")
    st2.close()
    st.close()


# ------------------------------------------------------------ pyqtgraph export
def _pyqtgraph_export(say, ctx, wins, dismiss):
    """The live plot's exporters: View ▸ Save plot image… in PNG / TIFF /
    JPG / SVG (pyqtgraph's ImageExporter and SVGExporter, lazy imports a
    frozen build can miss), the plot's own right-click Export figure… and
    Send to Plotting studio, Edit ▸ Copy plot to the clipboard."""
    from PySide6.QtWidgets import QApplication

    from larmor.desktop.plotting_studio import PlottingStudio

    out = Path(ctx["out"]) / "plot"
    out.mkdir(parents=True, exist_ok=True)
    ex = _examples(ctx)
    win = _new_window(wins)
    if ex:
        _load_1d(win, ex["al_1r"], f"{_EXPNO_27AL}/pdata/1/1r", "27Al")
    else:
        say("no bundled examples: a synthetic 11B spectrum stands in")
        _load_1d(win, _write_synthetic(out, "s0.csv", 13.0, 100.0), "s0.csv")
    sizes = {}
    for fmt, ext in (("PNG", "png"), ("SVG", "svg"), ("TIFF", "tiff"), ("JPG", "jpg")):
        target = out / f"plot.{ext}"
        with _export_choice(fmt), _export_path(target):
            win.save_plot_image()
        if "plot saved" not in _status(win):
            raise AssertionError(f"Save plot image {fmt}: {_status(win)!r}")
        sizes[ext] = _check_file(target, kind=ext)
    say("Save plot image: " + ", ".join(f"{e} {s} B" for e, s in sizes.items()))

    # the plot's right-click menu (plot_menu.attach_plot_menu)
    menu = win.view.getPlotItem().getViewBox().menu
    acts = {a.text(): a for a in menu.actions() if a.text()}
    export_act = next((a for t, a in acts.items() if t.startswith("Export figure")), None)
    studio_act = next((a for t, a in acts.items() if t.startswith("Send to Plotting studio")), None)
    if export_act is None or studio_act is None:
        raise AssertionError(f"the plot's right-click menu lacks its entries: {list(acts)}")
    target = out / "menu.svg"
    with _export_choice("SVG"), _export_path(target):
        export_act.trigger()
    _check_file(target, kind="svg")
    created: list = []
    with _exec_captured(PlottingStudio, created):
        studio_act.trigger()
    if not created or not created[-1]._traces:
        raise AssertionError("Send to Plotting studio opened no studio with the plotted curve")
    created[-1].close()
    say(f"right-click ▸ Export figure… wrote {target.name} ({target.stat().st_size} B); "
        f"Send to Plotting studio carried {len(created[-1]._traces)} trace(s)")

    # Edit ▸ Copy plot: an image on the clipboard
    clip = QApplication.clipboard()
    clip.clear()
    win.copy_plot()
    md = clip.mimeData()
    if md is None or not md.hasImage():
        raise AssertionError(f"Copy plot left no image on the clipboard: {_status(win)!r}")
    pix = clip.pixmap()
    say(f"Copy plot: a {pix.width()}x{pix.height()} image on the clipboard")


# ------------------------------------------------------------ overlays
def _overlays(say, ctx, wins, dismiss):
    """The Datasets dock: Add… (the file dialog), colour, match height, the
    per-row scale / shift / offset and Reset, visibility, Shift-drop, the
    acquisition comparison, Remove, Make active, Remove all."""
    from larmor import display

    out = Path(ctx["out"]) / "overlays"
    out.mkdir(parents=True, exist_ok=True)
    ex = _examples(ctx)
    win = _new_window(wins)
    if ex:
        active, other, nuc_other = ex["al_1r"], ex["b_1r"], "11B"
        _load_1d(win, active, f"{_EXPNO_27AL}/pdata/1/1r", "27Al")
    else:
        say("no bundled examples: two synthetic 11B spectra stand in")
        paths, _model = _synthetic_series(out)
        active, other, nuc_other = Path(paths[0]), Path(paths[1]), "11B"
        _load_1d(win, active, "s0.csv")
    with _open_dialog(other):
        win.add_overlay_dialog()
    if len(win._overlays) != 1:
        raise AssertionError(f"Add… did not add the overlay: {_status(win)!r}")
    ov = win._overlays[0]
    if len(win.datasets_panel._row_widgets) != 1:
        raise AssertionError("the Datasets dock shows no row for the overlay")
    say(f"overlay added through the file dialog: {ov['label']!r}, {ov.get('npts')} points, "
        f"{ov.get('nucleus')!r}")

    win.overlay_set_color(0, "#ff0000")
    if win._overlays[0]["color"] != "#ff0000":
        raise AssertionError("the overlay colour did not change")
    win.datasets_panel.match.setChecked(True)          # a standing request
    pump(20)
    scale = win._overlays[0]["scale"]
    if not (np.isfinite(scale) and scale > 0):
        raise AssertionError(f"match height produced the scale {scale!r}")
    win.datasets_panel.match.setChecked(False)
    pump(20)
    if win._overlays[0]["scale"] != 1.0:
        raise AssertionError("unticking match height did not put the scale back to x1")
    win.overlay_set_scale(0, 2.0)
    win.overlay_set_shift(0, 1.5)
    win.overlay_set_yoff(0, 0.1)
    o = win._overlays[0]
    if (o["scale"], o["shift"], o["yoff"]) != (2.0, 1.5, 0.1):
        raise AssertionError("the row's scale / shift / offset did not land")
    win.overlay_reset(0)
    if any(win._overlays[0][k] != v for k, v in display.OVERLAY_DEFAULTS.items()):
        raise AssertionError("Reset did not restore the display defaults")
    win.overlay_visibility(0, False)
    pump(20)
    if len(win.view._overlay_items) != 0:
        raise AssertionError("a hidden overlay is still drawn")
    win.overlay_visibility(0, True)
    pump(20)
    if len(win.view._overlay_items) != 1:
        raise AssertionError("a re-shown overlay is not drawn")
    win._toggle_overlays_visible(False)
    win._toggle_overlays_visible(True)
    say(f"colour, match height (x{scale:.3g}), scale / shift / offset, Reset, visibility ok")

    # Shift + drop of a file onto the plot: the view's signal
    win.view.files_dropped_overlay.emit([str(active)])
    pump(50)
    if len(win._overlays) != 2:
        raise AssertionError(f"Shift-drop did not overlay the file: {_status(win)!r}")
    n_tools = len(getattr(win, "_tool_windows", []))
    win.compare_overlays()                             # Datasets ▸ Compare acquisition…
    pump(100)
    if ex and len(getattr(win, "_tool_windows", [])) <= n_tools:
        raise AssertionError(f"Compare acquisition… opened no table: {_status(win)!r}")
    say("Shift-drop overlaid the active spectrum's own file; "
        + ("the acquisition comparison table opened" if ex else
           "no Bruker acquisition to compare (CSV spectra)"))
    win.overlay_remove(1)
    if len(win._overlays) != 1:
        raise AssertionError("Remove left the overlay in place")

    # Make active: the overlay becomes the spectrum on the workbench (in
    # place, no new workspace) and the former active becomes an overlay
    prev = win.source_path
    win.overlay_make_active(0)
    pump(50)
    _settle(win)
    if win.source_path != str(other) or (win.recipe or {}).get("nucleus") != nuc_other:
        raise AssertionError(f"Make active did not load the overlay's file: {_status(win)!r}")
    if len(win._overlays) != 1 or win._overlays[0].get("source") != prev:
        raise AssertionError("the former active spectrum did not become an overlay")
    if len(win.workspaces) != 1:
        raise AssertionError("Make active opened a second workspace instead of swapping in place")
    win.clear_overlays()
    if win._overlays:
        raise AssertionError("Remove all left overlays behind")
    say(f"Make active swapped the two spectra in place (active now {nuc_other}); Remove all ok")


# ------------------------------------------------------------ series mode
def _series_mode(say, ctx, wins, dismiss):
    """Three synthetic spectra as a series: Copy model, Fit → next (the
    FitWorker chain), Seed, Keep this fit, the auto sweep (SeqWorker), the
    walk, Plot… with its CSV and figure exports, Table…, Acquisition…,
    Save all fits…, the publication bundle, End series."""
    from PySide6.QtWidgets import QDialog

    from larmor.desktop.acquisition_dialog import AcquisitionTableDialog
    from larmor.desktop.series_plot import SeriesPlotDialog
    from larmor.desktop.series_table_dialog import SeriesTableDialog
    from larmor.seriesmode import rmsd_of

    out = Path(ctx["out"]) / "series"
    out.mkdir(parents=True, exist_ok=True)
    win = _new_window(wins)
    paths, model = _synthetic_series(out / "spectra")
    say("three synthetic 11B spectra, one gauss_lor line at 13 / 15 / 17 ppm")
    win.start_series(paths, model)
    _settle(win)
    spec = win._series
    if spec is None or spec.n != 3 or len(win.workspaces) != 3 or win.active_ws != 0:
        raise AssertionError(f"start_series did not tag three workspaces: {_status(win)!r}")
    if not win.series_bar.isVisibleTo(win) or len(win.series_bar.member_buttons()) != 3:
        raise AssertionError("the series bar did not show its three members")
    if len(win.recipe["sites"]) != 1:
        raise AssertionError("the model on screen did not land on the first member")
    say(f"series of {spec.names()}: bar up, model on member 1")

    win.series_copy_model()                            # Carry ▾ Copy this model to every spectrum
    if "copied into 2 spectra" not in _status(win):
        raise AssertionError(f"Copy model: {_status(win)!r}")
    t0 = time.perf_counter()
    win.series_fit_then_next()                         # Fit → next
    if win._fit_worker is None:
        raise AssertionError(f"Fit → next started no fit: {_status(win)!r}")
    _wait(lambda: _fit_idle(win) and win._series_chain is None and win.active_ws == 1,
          _FIT_TIMEOUT, "Fit → next to land on member 2")
    _settle(win)
    r0 = rmsd_of(win.workspaces[0]["snap"]["recipe"])
    if r0 is None or win.series_bar.statuses()[0] != "fitted":
        raise AssertionError(f"member 1 is not fitted after Fit → next: {_status(win)!r}")
    say(f"Fit → next: member 1 fitted (RMSD {r0:.4g}) in {time.perf_counter() - t0:.1f} s, "
        "now on member 2")

    win.series_seed(2, 0)                              # Carry ▾ Seed … from spectrum 1
    msg = _status(win)
    if "seeded" not in msg and "copied" not in msg:
        raise AssertionError(f"Seed: {msg!r}")
    win.series_set_locked(0, True)                     # Keep this fit
    if not win._series.members[0].locked or win.series_bar.locked() != [True, False, False]:
        raise AssertionError("Keep this fit did not lock member 1")
    win.series_set_locked(0, False)
    if win._series.members[0].locked:
        raise AssertionError("the lock did not release")

    t0 = time.perf_counter()
    win.series_auto_sweep(2, "first", 0)               # Auto sweep: 2 passes
    if win._seq_worker is None:
        raise AssertionError(f"the auto sweep did not start: {_status(win)!r}")
    _wait(lambda: not _running(win._seq_worker) and win._active_fit_worker is None,
          _FIT_TIMEOUT, "the auto sweep")
    _settle(win)
    statuses = win.series_bar.statuses()
    if statuses != ["fitted"] * 3:
        raise AssertionError(f"after the sweep the members read {statuses}: {_status(win)!r}")
    pos = [win._series_recipe_at(i)["sites"][0]["params"]["isotropic_chemical_shift_ppm"]["value"]
           for i in range(3)]
    if not np.allclose(pos, [13.0, 15.0, 17.0], atol=0.5):
        raise AssertionError(f"the sweep landed the lines at {pos}, expected 13 / 15 / 17 ppm")
    say(f"auto sweep (2 passes) in {time.perf_counter() - t0:.1f} s: positions "
        f"{', '.join(f'{p:.2f}' for p in pos)} ppm")
    win.series_go(0)
    _settle(win)
    win.series_next()
    _settle(win)
    if win.active_ws != 1 or win.series_bar.current() != 1:
        raise AssertionError("the walk (click, ▶) did not move to member 2")

    # Plot…: a tool window; the line selected, its CSV and figure exported
    win.series_plot()
    pump(100)
    dlg = next((w for w in getattr(win, "_tool_windows", []) if isinstance(w, SeriesPlotDialog)),
               None)
    if dlg is None:
        raise AssertionError(f"Plot… opened no Series plot: {_status(win)!r}")
    if dlg.list.count() == 0:
        raise AssertionError("the Series plot lists no line to pick")
    dlg.list.item(0).setSelected(True)
    dlg._draw()
    pump(50)
    if not dlg._selected():
        raise AssertionError("the Series plot's line could not be selected")
    csv_path = out / "series.csv"
    with _save_dialog(csv_path):
        dlg._export_csv()
    _check_file(csv_path, min_bytes=50)
    with open(csv_path, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    if len(rows) != 4 or rows[0][0] != "spectrum":
        raise AssertionError(f"series.csv has {len(rows)} rows, header {rows[:1]}")
    fig_path = out / "series.png"
    with _export_choice("PNG"), _export_path(fig_path):
        dlg._export_fig()
    _check_file(fig_path, kind="png")
    say(f"Series plot: {len(rows) - 1} rows x {len(rows[0])} columns to {csv_path.name}, "
        f"figure {fig_path.name} ({fig_path.stat().st_size} B)")
    dlg.close()

    with _patched(SeriesTableDialog, "exec", lambda self, *a, **k: QDialog.Accepted):
        win.series_table()                             # Table… OK, nothing changed
    if "series table applied" not in _status(win):
        raise AssertionError(f"Table…: {_status(win)!r}")
    win.series_acquisition_table()                     # Acquisition…
    pump(100)
    acq = next((w for w in getattr(win, "_tool_windows", [])
                if isinstance(w, AcquisitionTableDialog)), None)
    if acq is None:
        raise AssertionError(f"Acquisition… opened no table: {_status(win)!r}")
    say(f"Table… applied; Acquisition… built with {len(acq.blocks)} Bruker block(s) "
        "(CSV spectra carry no acqus)")
    acq.close()

    fits = out / "fits"
    with _dir_dialog(fits), _answer_yes():
        win.series_save_all()                          # Save ▾ Save all fits…
    files = sorted(fits.glob("*.recipe.json"))
    if len(files) != 3:
        raise AssertionError(f"Save all fits wrote {len(files)} recipes: {_status(win)!r}")
    bundle = out / "bundle"
    with _dir_dialog(bundle):
        win.series_bundle()                            # Save ▾ Publication bundle…
    names = {p.name for p in bundle.iterdir()}
    if not {"manifest.csv", "seq_table.csv", "README.txt"} <= names:
        raise AssertionError(f"the series bundle lacks its tables: {sorted(names)}")
    if len(list(bundle.glob("*.recipe.json"))) != 3 or len(list(bundle.glob("*_curves.csv"))) != 3:
        raise AssertionError(f"the series bundle lacks recipes / curves: {sorted(names)}")
    with open(bundle / "manifest.csv", newline="", encoding="utf-8") as f:
        manifest = list(csv.DictReader(f))
    if len(manifest) != 3 or any(r.get("note") for r in manifest):
        raise AssertionError(f"manifest.csv: {manifest}")
    say(f"Save all fits: {len(files)} recipes; bundle: {len(names)} files, manifest of "
        f"{len(manifest)} fitted rows")

    win.end_series()
    if win._series is not None or win.series_bar.isVisibleTo(win) or len(win.workspaces) != 3:
        raise AssertionError(f"End series: {_status(win)!r}")
    say("End series: the bar is gone, the three spectra stay open")


# ------------------------------------------------------------ 2D
def _twod(say, ctx, wins, dismiss):
    """The example 3QMAS map in the 2D view: F2 skyline and a row to the
    workbench and back, a 2D site, 1D overlays on both projections (the
    current 1D, a file, an Explorer pick), the 2D MQMAS dialog, the contour
    export; the 2D fit outside quick mode."""
    from larmor.desktop.twod_dialog import TwoDDialog

    out = Path(ctx["out"]) / "twod"
    out.mkdir(parents=True, exist_ok=True)
    ex = _examples(ctx)
    win = _new_window(wins)
    if ex:
        win.load_source(str(ex["twod_expno"]), keep_fit=False)      # Open EXPNO… on 3620
        pump(100)
        name = _EXPNO_2D
        src1d = ex["al_1r"]
    else:
        from larmor import twod
        say("no bundled examples: a synthetic 2D map and a synthetic 1D stand in")
        f2 = np.linspace(-50.0, 50.0, 120)
        f1 = np.linspace(-30.0, 30.0, 60)
        z = (np.exp(-((f2[None, :] - 10) / 4) ** 2) * np.exp(-((f1[:, None] - 5) / 4) ** 2))
        win._show_2d(twod.Data2D(f2_ppm=f2, f1_ppm=f1, z=z, nucleus="27Al",
                                 larmor_MHz=130.3), "synthetic", "2D")
        pump(100)
        name = "synthetic map"
        src1d = Path(_write_synthetic(out, "s0.csv", 13.0, 100.0, nucleus="27Al",
                                      larmor_MHz=130.3))
    if win.central_stack.currentWidget() is not win.view2d:
        raise AssertionError(f"{name} did not open in the 2D view: {_status(win)!r}")
    d = win.view2d.data
    if d is None or d.z.ndim != 2 or not win.workspaces or win.workspaces[0]["kind"] != "2d":
        raise AssertionError("the 2D view holds no map")
    say(f"{name} opened: {d.z.shape[0]} x {d.z.shape[1]} points, {d.nucleus or '?'}, "
        f"fittable {bool(getattr(win, '_data2d_fittable', False))}")

    # traces to the workbench (each a new workspace), and Back to 2D
    win.view2d._emit_projection("skyline")             # F2 skyline → fit
    pump(50)
    _settle(win)
    if win.central_stack.currentWidget() is not win.view or len(win.workspaces) != 2:
        raise AssertionError(f"the F2 skyline did not become a 1D workspace: {_status(win)!r}")
    if win.exp_ppm.size != d.z.shape[1]:
        raise AssertionError(f"the skyline has {win.exp_ppm.size} points for {d.z.shape[1]} columns")
    win.back_to_2d()
    pump(50)
    if win.central_stack.currentWidget() is not win.view2d or win.active_ws != 0:
        raise AssertionError("Back to 2D did not return to the map")
    if win.view2d._slice_line is None:
        raise AssertionError("the 2D view has no row cursor after drawing the map")
    win.view2d._emit_row()                             # row at cursor → fit
    pump(50)
    _settle(win)
    if win.central_stack.currentWidget() is not win.view or len(win.workspaces) != 3:
        raise AssertionError(f"the row did not become a 1D workspace: {_status(win)!r}")
    win.back_to_2d()
    pump(50)
    say("F2 skyline and the row at the cursor each became a 1D workspace; Back to 2D restores the map")

    # a 2D site placed by a click, with the czjzek model armed
    if getattr(win, "_data2d_fittable", False):
        win._model_actions["czjzek"].setChecked(True)
        win.add_site_2d(float(np.mean(d.f2_ppm)), float(np.mean(d.f1_ppm)))
        pump(20)
        win._set_add_mode(None)
        if len(win.recipe["sites"]) != 1 or win.recipe["sites"][0]["model"] != "czjzek":
            raise AssertionError(f"add_site_2d placed no czjzek site: {_status(win)!r}")
        say("a czjzek site placed on the map by a click")

    # 1D overlays on the projections: the current 1D, a file, an Explorer pick
    win.overlay_1d_on_2d("f2", "current")
    pump(20)
    if win.view2d._hmqc["f2"] is None:
        raise AssertionError(f"Overlay 1D ▸ current on F2: {_status(win)!r}")
    with _open_dialog(src1d):
        win.overlay_1d_on_2d("f1", "file")
    pump(20)
    if win.view2d._hmqc["f1"] is None:
        raise AssertionError(f"Overlay 1D ▸ from file on F1: {_status(win)!r}")
    win.view2d.clear_projection_1d("f1")
    win.load_projection_1d("f1")                       # HMQC: arm the Explorer pick
    if win._proj_pick_axis != "f1":
        raise AssertionError("load_projection_1d did not arm the pick")
    win._explorer_open(str(src1d))                     # the Explorer click
    pump(20)
    if win._proj_pick_axis is not None or win.view2d._hmqc["f1"] is None:
        raise AssertionError(f"the Explorer pick did not land on F1: {_status(win)!r}")
    say("1D overlays on F2 (current) and F1 (file, then an Explorer pick)")

    # Fit ▸ 2D MQMAS…: the dialog on the open EXPNO
    created: list = []
    with _exec_captured(TwoDDialog, created):
        win.open_twod()
    if not created:
        raise AssertionError("open_twod built no dialog")
    dlg = created[-1]
    if ex and dlg.data is None:
        raise AssertionError(f"the 2D MQMAS dialog did not load {name}: {dlg.res.text()!r}")
    say(f"2D MQMAS dialog built: {dlg.res.text()!r}")
    dlg.close()

    # the contour map through the shared exporter (View ▸ Save plot image…)
    target = out / "map.png"
    with _export_choice("PNG"), _export_path(target):
        win.save_plot_image()
    _check_file(target, kind="png")
    say(f"the 2D map exported to {target.name} ({target.stat().st_size} B)")

    # the 2D fit: a kernel build plus the fit, tens of seconds -- full mode only
    if ctx.get("quick"):
        say("quick mode: the 2D MQMAS fit (kernel build + fit) is skipped")
    elif getattr(win, "_data2d_fittable", False) and win.recipe.get("sites"):
        t0 = time.perf_counter()
        win.run_fit_2d()
        if win._fit2d_worker is None:
            raise AssertionError(f"run_fit_2d started no worker: {_status(win)!r}")
        _wait(lambda: win._active_fit_worker is None and not _running(win._fit2d_worker),
              240.0, "the 2D fit")
        msg = _status(win)
        if "2D fit done" not in msg:
            raise AssertionError(f"the 2D fit did not finish normally: {msg!r}")
        if win.view2d._model is None:
            raise AssertionError("the fitted 2D model overlay was not drawn")
        say(f"2D fit done in {time.perf_counter() - t0:.1f} s: {msg}")
    else:
        say("no fittable map with a site: the 2D fit is not exercised")


# ------------------------------------------------------------ processing
def _processing(say, ctx, wins, dismiss):
    """The Processing panel on the example 27Al: phase, Autophase, the
    polynomial / subtract-averages / iterative / pybaselines baselines, the
    steps editor, Reset to original, the raw-fid window functions and
    zero-fill, drag-to-phase, manual anchors, the two-point background,
    FID ⇄ spectrum and the display channels, Open FID…, Experiment
    parameters…, Calibrate, Measure."""
    from PySide6.QtWidgets import QDialog, QInputDialog

    from larmor.desktop.baseline_dialog import BaselineDialog
    from larmor.desktop.dialogs import ExperimentDialog, ProcessingStepsDialog
    from larmor.desktop.fid_dialog import FidDialog
    from larmor.desktop.pybaseline_dialog import PybaselineDialog
    from larmor.phasedrag import wrap_p0

    out = Path(ctx["out"]) / "processing"
    out.mkdir(parents=True, exist_ok=True)
    ex = _examples(ctx)
    win = _new_window(wins)
    if ex:
        _load_1d(win, ex["al_1r"], f"{_EXPNO_27AL}/pdata/1/1r", "27Al")
    else:
        say("no bundled examples: a synthetic 11B spectrum stands in (no raw fid)")
        _load_1d(win, _write_synthetic(out, "s0.csv", 13.0, 100.0), "s0.csv")
    pp = win.proc_panel
    original = win.exp_amp.copy()

    def ops() -> list:
        return [o["op"] for o in (win.recipe.get("processing") or [])]

    def applied(fn, what: str):
        n0 = win._proc_apply_count
        fn()
        pump(50)
        if win._proc_apply_count == n0:
            raise AssertionError(f"{what}: the pipeline did not run ({_status(win)!r})")
        if "failed" in _status(win):
            raise AssertionError(f"{what}: {_status(win)!r}")

    applied(lambda: pp.set_phase(15.0, 0.0, apply=True), "phase p0 15°")
    if "phase" not in ops() or np.allclose(win.exp_amp, original):
        raise AssertionError("a 15° phase left the spectrum or the chain unchanged")
    applied(lambda: win.append_processing_step({"op": "autophase"}), "Autophase")
    if "autophase" not in _status(win):
        raise AssertionError(f"Autophase: {_status(win)!r}")
    applied(lambda: win.append_processing_step({"op": "baseline", "order": 3}),
            "polynomial baseline")
    applied(lambda: win.append_processing_step({"op": "subtract_avg"}), "subtract averages")
    with _patched(BaselineDialog, "exec", lambda self, *a, **k: QDialog.Accepted):
        applied(win.apply_iterbaseline, "iterative baseline")
    with _patched(PybaselineDialog, "exec", lambda self, *a, **k: QDialog.Accepted):
        applied(lambda: win.apply_pybaseline("arpls"), "pybaselines arpls")
    pyb = [o for o in win.recipe["processing"] if o["op"] == "pybaseline"]
    if not {"baseline", "subtract_avg", "iterbaseline"} <= set(ops()) or not pyb \
            or pyb[-1].get("method") != "arpls":
        raise AssertionError(f"the recorded chain lacks its baseline steps: {ops()}")
    say(f"recorded chain after the Process menu: {ops()}")
    with _patched(ProcessingStepsDialog, "exec", lambda self, *a, **k: QDialog.Accepted):
        applied(win.edit_processing_steps, "Processing steps… OK")
    n_ws = len(win.workspaces)
    win.reset_processing()                             # Reset to original
    pump(50)
    _settle(win)
    if len(win.workspaces) != n_ws:
        raise AssertionError("Reset to original opened a new workspace")
    if not np.allclose(win.exp_amp, original):
        raise AssertionError("Reset to original did not bring the loaded spectrum back")
    if ops():
        raise AssertionError(f"Reset to original kept the processing record {ops()} while "
                             "the data went back to the unprocessed source")
    if pp.phase_values() != (0.0, 0.0):
        raise AssertionError(f"Reset to original left the panel at p0/p1 {pp.phase_values()}")
    say("Processing steps… re-applied; Reset to original reloaded in place, the record and "
        "the panel cleared")

    # Source ▸ raw fid: the window functions and the zero-fill, applied live
    if ex:
        pp.rb_raw.setChecked(True)
        pp._live_timer.stop()
        for wdw in ("EM", "GM", "QSINE"):
            pp.wdw.setCurrentText(wdw)
            pp._live_timer.stop()
            applied(pp.btnApply.click, f"raw fid · {wdw}")
            if not win.recipe.get("processing_from_raw") or win._proc_fid is None:
                raise AssertionError(f"raw fid · {wdw}: no fid captured")
        pp.zf.setValue(4)
        pp._live_timer.stop()
        applied(pp.btnApply.click, "raw fid · ZF 4")
        n_fid = win._proc_fid.y.size
        say(f"raw-fid chain {ops()}: FID of {n_fid} points, spectrum of {win.exp_ppm.size}")
        pp.rb_pdata.setChecked(True)
        pp._live_timer.stop()
        applied(pp.btnApply.click, "back to pdata")
    else:
        say("CSV source: the raw-fid window functions are not exercised")

    # Process ▸ Phase ▸ Drag to phase: the gesture's slots
    win.start_phase_drag()
    pump(20)
    if not pp.btnDrag.isChecked() or not win.view.phase_drag_active():
        raise AssertionError("Drag to phase did not arm")
    p0_before = pp.phase_values()[0]
    win.on_phase_dragged(12.0, 0.0)
    pump(20)
    win.on_phase_drag_released()
    pump(20)
    if abs(pp.phase_values()[0] - wrap_p0(p0_before + 12.0)) > 1e-6:
        raise AssertionError(f"the drag did not move p0 by 12°: {pp.phase_values()}")
    win.start_phase_drag()
    pump(20)
    if pp.btnDrag.isChecked():
        raise AssertionError("Drag to phase did not disarm")
    say(f"drag to phase: p0 {p0_before:+.1f}° → {pp.phase_values()[0]:+.1f}°, one gesture")

    # manual baseline anchors (Pick anchors; releasing the button applies)
    x0, x1 = float(win.exp_ppm.min()), float(win.exp_ppm.max())
    pp.btnBlPick.setChecked(True)
    pump(10)
    for frac in (0.1, 0.5, 0.9):
        win.view._add_baseline_anchor(x0 + frac * (x1 - x0), 0.0)
    pp.btnBlPick.setChecked(False)
    pump(20)
    if "manual baseline subtracted" not in _status(win):
        raise AssertionError(f"manual baseline: {_status(win)!r}")
    win.start_twopoint_bg()                            # Process ▸ Baseline ▸ 2-point background…
    pump(10)
    if not pp.btnTpPick.isChecked():
        raise AssertionError("2-point background did not arm the picker")
    win.view._add_baseline_anchor(x0 + 0.05 * (x1 - x0), 0.0)
    win.view._add_baseline_anchor(x1 - 0.05 * (x1 - x0), 0.0)
    pp.btnTpPick.setChecked(False)
    pump(20)
    if "2-point" not in _status(win):
        raise AssertionError(f"2-point background: {_status(win)!r}")
    say("manual anchors and the 2-point background subtracted")

    # FID ⇄ spectrum (Ctrl+T) and the display channels
    win._toggle_time_domain(True)
    pump(50)
    if win.view.domain != "time":
        raise AssertionError(f"FID ⇄ spectrum did not show the FID: {_status(win)!r}")
    win._toggle_time_domain(False)
    pump(50)
    if win.view.domain != "freq":
        raise AssertionError("FID ⇄ spectrum did not return to the spectrum")
    for channel in ("imag", "magnitude", "real"):
        win._set_channel(channel)
        pump(50)
        if pp.view_state() != ("freq", channel):
            raise AssertionError(f"channel {channel}: the panel reads {pp.view_state()} "
                                 f"({_status(win)!r})")
    say("FID ⇄ spectrum and the real / imag / |S| channels")

    # File ▸ Open FID… on the example fid, accepted with the dialog's defaults
    if ex:
        created: list = []
        with _exec_captured(FidDialog, created, result=QDialog.Accepted,
                            before=lambda dlg: dlg._accept()):
            win.open_fid_path(str(ex["al_fid"]))
        pump(50)
        _settle(win)
        if not created or created[-1].data is None:
            raise AssertionError("the Open FID dialog did not load the fid")
        if "processed FID loaded" not in _status(win) or not win.recipe.get("processing_from_raw"):
            raise AssertionError(f"Open FID: {_status(win)!r}")
        say(f"Open FID: {_EXPNO_27AL}/fid processed with the defaults → {win.exp_ppm.size} points, "
            f"chain {ops()}")

    # Process ▸ Experiment parameters… (OK with the values as they are; the
    # per-session MAS confirmation store is pointed into the scratch folder)
    created = []
    env_before = os.environ.get("LARMOR_MAS_LOG")
    os.environ["LARMOR_MAS_LOG"] = str(out / "mas_confirmations.jsonl")
    try:
        with _exec_captured(ExperimentDialog, created, result=1, before=lambda dlg: dlg._accept()):
            win.edit_experiment()
    finally:
        if env_before is None:
            os.environ.pop("LARMOR_MAS_LOG", None)
        else:
            os.environ["LARMOR_MAS_LOG"] = env_before
    if not created or "experiment updated" not in _status(win):
        raise AssertionError(f"Experiment parameters…: {_status(win)!r}")
    _settle(win)
    say(f"Experiment parameters… accepted: {_status(win)[:90]!r}")

    # Process ▸ Reference ▸ Calibrate axis…: the picked peak is set 1 ppm higher
    peak = float(win.exp_ppm[int(np.argmax(win.exp_amp))])
    with _patched(QInputDialog, "getDouble", staticmethod(lambda *a, **k: (peak + 1.0, True))):
        win.start_calibrate()
        win.on_calibrate_picked(peak)
    pump(20)
    moved = float(win.exp_ppm[int(np.argmax(win.exp_amp))])
    if abs(moved - (peak + 1.0)) > 1e-6 or "axis shifted by +1.000 ppm" not in _status(win):
        raise AssertionError(f"calibrate moved the peak from {peak:.3f} to {moved:.3f}: "
                             f"{_status(win)!r}")
    win.undo()                                         # calibration is undoable
    _settle(win)
    win.toggle_measure(True)
    win.on_measure_changed(10.0, 20.0)
    if "Δ = 10.000 ppm" not in _status(win):
        raise AssertionError(f"Measure: {_status(win)!r}")
    win.toggle_measure(False)
    say(f"calibrate ({peak:.2f} → {moved:.2f} ppm, undone) and Measure")


# ------------------------------------------------------------ project round trip
def _project_roundtrip(say, ctx, wins, dismiss):
    """A fitted spectrum with an overlay, a three-member series and the 2D
    map: File ▸ Save project…, then a FRESH window opens it (replace mode)
    and everything is back -- workspaces, titles, the fit, the overlay by
    reference, the series tags and bar, the 2D by reference."""
    from larmor import project
    from larmor.recipe import Recipe
    from larmor.seriesmode import rmsd_of

    out = Path(ctx["out"]) / "project"
    out.mkdir(parents=True, exist_ok=True)
    ex = _examples(ctx)
    win = _new_window(wins)
    if ex:
        _load_1d(win, ex["al_1r"], f"{_EXPNO_27AL}/pdata/1/1r", "27Al")
        win.apply_recipe(str(ex["recipe_al"]))        # the shipped, already fitted, recipe
        _settle(win, _FIT_TIMEOUT)
        win.add_overlay_path(str(ex["b_1r"]))
        pump(20)
    else:
        say("no bundled examples: synthetic spectra stand in (no 2D map)")
        paths, model = _synthetic_series(out / "single")
        _load_1d(win, paths[0], "s0.csv")
        rec = out / "single" / "model.recipe.json"
        Recipe.from_dict(model).save(rec)
        win.apply_recipe(str(rec))
        _settle(win)
        win.add_overlay_path(paths[1])
        pump(20)
    n_sites = len(win.recipe["sites"])
    if not n_sites or len(win._overlays) != 1:
        raise AssertionError(f"the fitted spectrum with an overlay was not set up: {_status(win)!r}")
    paths, model = _synthetic_series(out / "spectra")
    win.start_series(paths, model)
    _settle(win)
    win.series_copy_model()
    win.series_fit_then_next()                         # member 1 fitted: an RMSD to find again
    if win._fit_worker is None:
        raise AssertionError(f"Fit → next started no fit: {_status(win)!r}")
    _wait(lambda: _fit_idle(win) and win._series_chain is None
          and win._series_current_index() == 1, _FIT_TIMEOUT, "Fit → next to land on member 2")
    _settle(win)
    sid = win._series.id
    member0 = win._series_rows()[0][0]
    r_member = rmsd_of(win.workspaces[member0]["snap"]["recipe"])
    if r_member is None:
        raise AssertionError(f"the first series member carries no RMSD after its fit: "
                             f"{_status(win)!r}")
    if ex:
        win.load_source(str(ex["twod_expno"]), keep_fit=False)
        pump(100)
        if win.central_stack.currentWidget() is not win.view2d:
            raise AssertionError(f"{_EXPNO_2D} did not open in the 2D view: {_status(win)!r}")
    expected = 1 + 3 + (1 if ex else 0)
    if len(win.workspaces) != expected:
        raise AssertionError(f"{len(win.workspaces)} workspaces open, expected {expected}")
    titles = [ws["title"] for ws in win.workspaces]
    kinds = [ws["kind"] for ws in win.workspaces]

    proj = out / "session.larproj.json"
    with _save_dialog(proj):
        win.save_project()
    if "project saved" not in _status(win):
        raise AssertionError(f"Save project…: {_status(win)!r}")
    _check_file(proj, min_bytes=1000)
    bundle = json.loads(proj.read_text(encoding="utf-8"))
    if [w["kind"] for w in bundle["workspaces"]] != kinds:
        raise AssertionError(f"the bundle holds {[w['kind'] for w in bundle['workspaces']]}, "
                             f"the session {kinds}")
    say(f"saved {proj.name} ({proj.stat().st_size} B): {project.summary(bundle)}")

    win2 = _new_window(wins)
    t0 = time.perf_counter()
    win2._open_project_path(str(proj))
    pump(100)
    _settle(win2, _FIT_TIMEOUT)
    if len(win2.workspaces) != expected:
        raise AssertionError(f"{len(win2.workspaces)} workspaces came back, expected {expected}: "
                             f"{_status(win2)!r}")
    if [ws["title"] for ws in win2.workspaces] != titles:
        raise AssertionError(f"titles came back as {[ws['title'] for ws in win2.workspaces]}, "
                             f"saved {titles}")
    if [ws["kind"] for ws in win2.workspaces] != kinds:
        raise AssertionError("the workspace kinds did not come back")
    rec0 = win2.workspaces[0]["snap"]["recipe"]
    if len(rec0.get("sites", [])) != n_sites or not win2.workspaces[0]["has_fit"]:
        raise AssertionError(f"the {n_sites} lines of workspace 1 did not come back")
    rows2 = win2._series_rows() if win2._series is not None else []
    r_back = rmsd_of(win2.workspaces[rows2[0][0]]["snap"]["recipe"]) if rows2 else None
    if r_back is None or abs(r_back - r_member) > 1e-12:
        raise AssertionError(f"the fitted series member came back with RMSD {r_back!r}, "
                             f"saved {r_member!r}")
    ovs = win2.workspaces[0]["snap"].get("overlays") or []
    if len(ovs) != 1:
        raise AssertionError(f"the overlay did not come back by reference ({len(ovs)})")
    spec2 = win2._series
    if spec2 is None or spec2.id != sid or spec2.n != 3 or spec2.names() != ["s0", "s1", "s2"]:
        raise AssertionError("the series did not come back from its tags")
    if not all(len(ws["snap"]["recipe"]["sites"]) == 1 for ws in win2.workspaces[1:4]):
        raise AssertionError("the series members lost their lines")
    if len(win2.series_bar.member_buttons()) != 3:
        raise AssertionError("the series bar did not come back with its three members")
    if ex:
        w2d = win2.workspaces[4]
        if w2d["kind"] != "2d" or w2d["snap"].get("data2d") is None:
            raise AssertionError("the 2D map did not come back by reference")
        if win2.central_stack.currentWidget() is not win2.view2d:
            raise AssertionError("the active (2D) workspace is not the one on screen")
        # the bar belongs to the 1D page: over the 2D map it stays hidden
        if win2.series_bar.isVisibleTo(win2):
            raise AssertionError("the series bar shows over the 2D map")
    # walking to a member brings the 1D page, and the bar, back on it
    win2.series_go(0, seed=False)
    pump(50)
    _settle(win2, _FIT_TIMEOUT)
    if not win2.series_bar.isVisibleTo(win2) or win2.series_bar.current() != 0:
        raise AssertionError("switching to the first series member did not bring the bar back "
                             f"on it: {_status(win2)!r}")
    say(f"reopened in a fresh window in {time.perf_counter() - t0:.1f} s: {len(win2.workspaces)} "
        f"workspaces {titles}, the {n_sites} lines of the first, {len(ovs)} overlay, the series "
        f"{spec2.names()} with its bar and its fitted member (RMSD {r_back:.4g})"
        + (", the 2D map by reference" if ex else ""))


# ------------------------------------------------------------ explorer + workspaces
def _explorer_and_workspaces(say, ctx, wins, dismiss):
    """The Explorer dock scanned on the example sample folder (experiments,
    their saved fits, selection, the saved-fit menu and its Fit parameter
    table, Dataset info…, Rename… cancelled), opening from it, and the
    Workspaces dock's menu with its slots (studio, fit table, overlay,
    rename dismissed, close)."""
    from PySide6.QtWidgets import QPlainTextEdit

    from larmor.desktop import explorer as expl
    from larmor.desktop.fittable_dialog import FitTableDialog
    from larmor.desktop.plotting_studio import PlottingStudio

    out = Path(ctx["out"]) / "explorer"
    out.mkdir(parents=True, exist_ok=True)
    ex = _examples(ctx)
    win = _new_window(wins)
    if ex:
        folder = ex["sample"]
        win.explorer.load_sample(str(folder))         # File ▸ Open sample folder…
        pump(50)
        items = list(win.explorer._iter_items())
        exps = [it for it in items if it.data(0, expl._ROLE_KIND) == "exp"]
        if len(exps) != 3:
            raise AssertionError(f"the Explorer lists {len(exps)} experiments in {folder.name}, "
                                 "expected 3")
        for it in exps:
            it.setExpanded(True)                      # the single-proc fits appear
            pump(10)
        fits = [it for it in win.explorer._iter_items() if it.data(0, expl._ROLE_KIND) == "fit"]
        if len(fits) < 2:
            raise AssertionError(f"the Explorer lists {len(fits)} saved fits, expected the two "
                                 "dmfit fits")
        for it in exps:
            it.setSelected(True)
        spectra = win.explorer.selected_spectra()
        if len(spectra) != 3:
            raise AssertionError(f"selected_spectra returned {len(spectra)} paths for 3 rows")
        say(f"Explorer: {len(exps)} experiments, {len(fits)} saved fits "
            f"({', '.join(Path(it.data(0, expl._ROLE_OPEN)).name for it in fits)}); "
            f"{len(spectra)} spectra selected")

        # the saved-fit menu, built, and its Fit parameter table (a tool window)
        fit_paths = [it.data(0, expl._ROLE_OPEN) for it in fits]
        menu = win.explorer.fit_menu(fit_paths[0], fit_paths)
        acts = {a.text(): a for a in menu.actions() if a.text()}
        table_act = next((a for t, a in acts.items() if t.startswith("Fit parameter table")), None)
        if table_act is None or "Open" not in acts or "Apply to the open spectrum" not in acts:
            raise AssertionError(f"the saved-fit menu reads {list(acts)}")
        table_act.trigger()
        pump(100)
        tbl = next((w for w in getattr(win, "_tool_windows", []) if isinstance(w, FitTableDialog)),
                   None)
        if tbl is None or tbl.table.rowCount() < 1:
            raise AssertionError(f"the saved fits gave no parameter table: {_status(win)!r}")
        csv_path = out / "explorer_fits.csv"
        with _save_dialog(csv_path):
            tbl._export_csv()
        _check_file(csv_path, min_bytes=20)
        say(f"saved-fit menu ▸ Fit parameter table: {tbl.table.rowCount()} rows, "
            f"exported {csv_path.name}")
        tbl.close()

        # Dataset info… (a tool window) and Rename… (its dialog dismissed)
        info = expl.show_dataset_info(win, str(ex["al_expno"]))
        pump(50)
        text = info.findChild(QPlainTextEdit).toPlainText()
        if "27Al" not in text:
            raise AssertionError(f"dataset info of {_EXPNO_27AL} does not name 27Al: {text[:120]!r}")
        info.close()
        seen0 = len(dismiss.seen)
        got = expl.rename_flow(win, str(ex["al_expno"]), open_paths=win._open_source_paths())
        if got is not None:
            raise AssertionError("rename_flow renamed although its dialog was dismissed")
        if len(dismiss.seen) == seen0:
            raise AssertionError("the Rename dialog never showed")
        say("Dataset info… names the acquisition; Rename… cancelled, nothing renamed")

        # open from the Explorer (a double-click), apply the fit saved next to it
        first = sorted(spectra)[0]                   # 1118: the 11B 1D
        win._explorer_open(first)
        pump(50)
        _settle(win)
        if win.central_stack.currentWidget() is not win.view:
            raise AssertionError(f"{Path(first).name} did not open from the Explorer: {_status(win)!r}")
        expno = Path(win.source_path).parents[2] if Path(win.source_path).name == "1r" \
            else Path(win.source_path)
        own = [p for p in fit_paths if Path(p).is_relative_to(expno)]
        if own:
            win.explorer.fit_menu(own[0], own).actions()
            win.apply_recipe(own[0])                  # "Apply to the open spectrum"
            _settle(win, _FIT_TIMEOUT)
            if not win.recipe["sites"]:
                raise AssertionError(f"the saved dmfit fit put no lines on its spectrum: "
                                     f"{_status(win)!r}")
            say(f"{Path(first).parent.parent.parent.name} opened from the Explorer; its dmfit fit "
                f"{Path(own[0]).name} applied: {len(win.recipe['sites'])} line(s)")
        second = sorted(spectra)[1]                   # 3616: the 27Al 1D
        win.load_source(second, keep_fit=False)
        pump(50)
        _settle(win)
    else:
        say("no bundled examples: synthetic spectra stand in; the Explorer scan is not exercised")
        paths, _model = _synthetic_series(out / "spectra")
        _load_1d(win, paths[0], "s0.csv")
        win._model_actions["gauss_lor"].setChecked(True)
        win.add_site_at(13.0, 100.0)
        win._set_add_mode(None)
        _settle(win)
        win.load_source(paths[1], keep_fit=False)
        pump(50)
        _settle(win)
    if len(win.workspaces) != 2:
        raise AssertionError(f"{len(win.workspaces)} workspaces open, expected 2")

    # the Workspaces dock: the menu for one row and for a multi-selection
    panel = win.ws_panel
    one = [a.text() for a in panel.build_menu(1, [1]).actions() if a.text()]
    both = [a.text() for a in panel.build_menu(1, [0, 1]).actions() if a.text()]
    if "Switch to" not in one or "Rename…" not in one or "Close" not in one:
        raise AssertionError(f"the single-row menu reads {one}")
    if "Send to Plotting studio  (2 spectra)" not in both or "Close 2 selected" not in both:
        raise AssertionError(f"the multi-row menu reads {both}")
    say(f"Workspaces menu: {len(one)} entries for one row, {len(both)} for two")

    created: list = []
    with _exec_captured(PlottingStudio, created):
        win.send_workspaces_to_studio([0, 1])         # Send to Plotting studio
    if not created or len(created[-1]._traces) < 2:
        raise AssertionError(f"Send to Plotting studio: {_status(win)!r}")
    say(f"Send to Plotting studio: {len(created[-1]._traces)} trace(s) from 2 workspaces")
    created[-1].close()

    fitted = [i for i, ws in enumerate(win.workspaces) if ws["has_fit"]]
    win.open_fit_table(fitted or None)                # Fit parameter table…
    pump(100)
    tbl2 = next((w for w in reversed(getattr(win, "_tool_windows", []))
                 if isinstance(w, FitTableDialog)), None)
    if tbl2 is None:
        raise AssertionError(f"Fit parameter table… opened no table: {_status(win)!r}")
    if fitted and tbl2.table.rowCount() < 1:
        raise AssertionError("the fit table of the fitted workspace is empty")
    if fitted and not tbl2.export_csv(out / "workspace_fits.csv"):
        raise AssertionError(f"the fit table did not export: {tbl2.status.text()!r}")
    say(f"Fit parameter table: {tbl2.table.rowCount()} rows from {len(fitted)} fitted workspace(s)"
        + (", exported" if fitted else ""))
    tbl2.close()

    win.overlay_workspaces([0])                       # Overlay on the active spectrum
    pump(20)
    if len(win._overlays) != 1:
        raise AssertionError(f"Overlay on the active spectrum: {_status(win)!r}")
    title = win.workspaces[0]["title"]
    seen0 = len(dismiss.seen)
    win.rename_workspace(0)                           # Rename… (dismissed)
    pump(20)
    if win.workspaces[0]["title"] != title or len(dismiss.seen) == seen0:
        raise AssertionError("Rename… did not show its dialog, or renamed although dismissed")
    win.close_workspaces([0])                         # Close
    pump(20)
    if len(win.workspaces) != 1:
        raise AssertionError("Close did not close the workspace")
    say("overlay from the dock, Rename… dismissed, Close: ok")


# ====================================================================== registry
def _stage(name: str, body):
    def run(say, ctx):
        _run(body, say, ctx)
    run.__name__ = name.replace("-", "_")
    run.__doc__ = body.__doc__
    return name, run


def stages() -> list:
    return [_stage("tool-fit-workflow", _fit_workflow),
            _stage("tool-plotting-studio", _plotting_studio),
            _stage("tool-pyqtgraph-export", _pyqtgraph_export),
            _stage("tool-overlays", _overlays),
            _stage("tool-series-mode", _series_mode),
            _stage("tool-2d", _twod),
            _stage("tool-processing", _processing),
            _stage("tool-project-roundtrip", _project_roundtrip),
            _stage("tool-explorer-and-workspaces", _explorer_and_workspaces)]
