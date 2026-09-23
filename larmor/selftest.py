"""Fit-engine self-test: does this installation fit at all?

Runs the two fits a first-time user tries -- a Gauss/Lorentz line and a
Czjzek line on a synthetic 27Al spectrum -- first through the Qt-free core
(``larmor.fit.fit``) and then through the desktop window exactly as a click
on Fit does (add a line, run the FitWorker thread, wait for its result), and
writes every step, timing and traceback to ``~/LARMOR_selftest.log``. Meant
for a machine where the frozen app opens but fitting "does nothing": the
log says which stage fails and why.

    LARMOR.exe --selftest            (the frozen app)
    python -m larmor.selftest        (a development install)

Exit code 0 when every stage passed, 1 otherwise. The log is appended to,
never overwritten, so several attempts stay visible.
"""
from __future__ import annotations

import datetime
import os
import platform
import sys
import time
import traceback
from pathlib import Path

__all__ = ["run", "log_path", "core_fits", "desktop_fits"]

NUCLEUS, LARMOR_MHZ, MAS_HZ = "27Al", 130.3, 20000.0


def log_path() -> Path:
    return Path(os.environ.get("LARMOR_SELFTEST_LOG") or (Path.home() / "LARMOR_selftest.log"))


class _Log:
    def __init__(self, path: Path):
        self.path = path
        self.f = open(path, "a", encoding="utf-8")
        self.failures = 0

    def __call__(self, *parts):
        line = " ".join(str(p) for p in parts)
        self.f.write(line + "\n")
        self.f.flush()
        if sys.stdout is not None:
            try:
                print(line, flush=True)
            except Exception:                            # noqa: BLE001
                pass

    def fail(self, stage: str, exc: BaseException | None = None):
        self.failures += 1
        self("FAIL", stage)
        if exc is not None:
            self(traceback.format_exc())

    def close(self):
        self.f.close()


def _environment(say):
    import numpy

    say(f"=== LARMOR self-test {datetime.datetime.now().isoformat(timespec='seconds')}")
    say("python", sys.version.split()[0], "| frozen:", bool(getattr(sys, "frozen", False)),
        "| executable:", sys.executable)
    say("platform:", platform.platform(), "| cpu count:", os.cpu_count(),
        "| stdout/stderr present:", sys.stdout is not None, sys.stderr is not None)
    import larmor
    say("larmor", larmor.__version__, "from", Path(larmor.__file__).parent)
    for name in ("numpy", "scipy", "lmfit", "mrsimulator", "uncertainties", "asteval", "csdmpy", "nmrglue"):
        try:
            mod = __import__(name)
            say(f"  {name} {getattr(mod, '__version__', '?')}")
        except Exception as exc:                          # noqa: BLE001
            say(f"  {name}: IMPORT FAILED: {exc!r}")
    try:
        blas = numpy.__config__.CONFIG.get("Build Dependencies", {}).get("blas", {})
        say("  numpy BLAS:", blas.get("name", "?"), blas.get("version", ""))
    except Exception:                                    # noqa: BLE001
        pass
    for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "QT_QPA_PLATFORM"):
        if os.environ.get(var):
            say(f"  env {var}={os.environ[var]}")


def _synthetic():
    """A two-line 27Al spectrum: a broad Czjzek site and a narrow Gaussian."""
    import numpy as np

    from larmor import engine
    from larmor.recipe import Param, Recipe, SiteModel

    x = np.linspace(150.0, -100.0, 2048)
    truth = Recipe(nucleus=NUCLEUS, larmor_frequency_MHz=LARMOR_MHZ, spin_rate_Hz=MAS_HZ, sites=[
        SiteModel(model="czjzek", label="Al4", params={
            "isotropic_chemical_shift_ppm": Param(62.0), "sigma_Cq_MHz": Param(1.6),
            "shift_fwhm_ppm": Param(6.0), "line_fwhm_ppm": Param(0.5), "amplitude": Param(100.0)}),
        SiteModel(model="gauss_lor", label="Al6", params={
            "isotropic_chemical_shift_ppm": Param(5.0), "shift_fwhm_ppm": Param(12.0),
            "amplitude": Param(40.0), "gl": Param(0.3, vary=False)}),
    ])
    _, y, _ = engine.simulate(truth, exp_ppm=x)
    y = np.asarray(y, float) + np.random.default_rng(1).normal(0.0, 0.4, x.size)
    return x, y


def _start_recipe(model: str):
    from larmor.recipe import Param, Recipe, SiteModel

    if model == "gauss_lor":
        sites = [SiteModel(model="gauss_lor", label="A", params={
            "isotropic_chemical_shift_ppm": Param(58.0), "shift_fwhm_ppm": Param(20.0),
            "amplitude": Param(60.0), "gl": Param(0.0, vary=False)})]
    else:
        sites = [SiteModel(model="czjzek", label="A", params={
            "isotropic_chemical_shift_ppm": Param(58.0), "sigma_Cq_MHz": Param(1.0),
            "shift_fwhm_ppm": Param(8.0), "line_fwhm_ppm": Param(0.5), "amplitude": Param(60.0)})]
    return Recipe(nucleus=NUCLEUS, larmor_frequency_MHz=LARMOR_MHZ, spin_rate_Hz=MAS_HZ,
                  sites=sites, sample="selftest")


def core_fits(say) -> bool:
    """larmor.fit.fit on the synthetic spectrum, both models. True when both
    converge with moved parameters and a finite RMSD."""
    ok = True
    try:
        t0 = time.time()
        x, y = _synthetic()
        say(f"synthetic spectrum: {x.size} points, simulated in {time.time() - t0:.1f} s")
    except Exception as exc:                              # noqa: BLE001
        say.fail("simulate the synthetic spectrum", exc)
        return False
    from larmor import fit as fitmod
    for model in ("gauss_lor", "czjzek"):
        try:
            rec = _start_recipe(model)
            before = {k: p.value for k, p in rec.sites[0].params.items()}
            calls = []
            t0 = time.time()
            res = fitmod.fit(rec, x, y, window_ppm=(150.0, -100.0),
                             iter_cb=lambda *a, **k: calls.append(1))
            after = {k: p.value for k, p in res.recipe.sites[0].params.items()}
            moved = any(abs(after[k] - before[k]) > 1e-9 for k in before)
            stderr = res.recipe.sites[0].params["amplitude"].stderr
            say(f"core {model}: {len(calls)} evaluations, {time.time() - t0:.1f} s, "
                f"rmsd {res.rmsd:.4g}, amplitude {after['amplitude']:.3g} ± {stderr}, "
                f"shift {after['isotropic_chemical_shift_ppm']:.2f} ppm, moved: {moved}")
            if not moved or not (res.rmsd == res.rmsd):
                say.fail(f"core {model}: the fit did not change the parameters")
                ok = False
        except Exception as exc:                          # noqa: BLE001
            say.fail(f"core {model}", exc)
            ok = False
    return ok


def desktop_fits(say, timeout_s: float = 300.0, spectrum: str | None = None) -> bool:
    """The window's own path: display the spectrum, add a line at the tallest
    point, run the FitWorker thread, wait, read back what the window shows.
    ``spectrum`` (a 1r / EXPNO / recipe path) replaces the synthetic data --
    `LARMOR.exe --selftest <path>` runs the test on the user's own file."""
    ok = True
    try:
        os.environ.setdefault("LARMOR_NO_SESSION", "1")
        os.environ.setdefault("LARMOR_NO_KERNEL_WARM", "1")
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        import numpy as np
        from PySide6.QtWidgets import QApplication, QMessageBox

        app = QApplication.instance() or QApplication([])
        dialogs: list = []
        for kind in ("warning", "critical", "information"):
            setattr(QMessageBox, kind, staticmethod(
                lambda *a, _k=kind, **kw: dialogs.append(
                    f"{_k}: {' | '.join(str(v) for v in a[1:3])}") or QMessageBox.Ok))
        from larmor.desktop.app import MainWindow

        t0 = time.time()
        win = MainWindow()
        win.show()
        app.processEvents()
        say(f"desktop: window built in {time.time() - t0:.1f} s")
        if spectrum:
            t0 = time.time()
            win.load_source(spectrum)
            for _ in range(50):
                app.processEvents()
                time.sleep(0.02)
            x, y = np.asarray(win.exp_ppm), np.asarray(win.exp_amp)
            say(f"desktop: loaded {spectrum} in {time.time() - t0:.1f} s: "
                f"{win.recipe.get('nucleus') if win.recipe else '?'} at "
                f"{(win.recipe or {}).get('larmor_frequency_MHz')} MHz, {x.size} points, "
                f"view domain {win.view.domain!r}, status {win.statusBar().currentMessage()!r}")
            if x.size == 0 or win.view.domain != "freq":
                say.fail("desktop: the file did not open as a 1D spectrum (Fit needs the "
                         "frequency-domain view; a raw fid must be processed first)")
                return False
        else:
            x, y = _synthetic()
            win._display_1d(x, y, NUCLEUS, LARMOR_MHZ, MAS_HZ, "selftest", "selftest")
        app.processEvents()
        i = int(np.argmax(y))
        for model in ("gauss_lor", "czjzek"):
            try:
                win.recipe["sites"] = []
                win.on_structure_changed()
                for name, act in win._model_actions.items():
                    act.setChecked(name == model)
                win.add_site_at(float(x[i]), float(y[i]))
                app.processEvents()
                before = [dict((k, v["value"]) for k, v in s["params"].items())
                          for s in win.recipe["sites"]]
                t0 = time.time()
                win.run_fit()
                running = bool(win._fit_worker and win._fit_worker.isRunning())
                say(f"desktop {model}: Fit pressed, worker running: {running}, "
                    f"status: {win.statusBar().currentMessage()!r}")
                while (win._fit_worker is not None and win._fit_worker.isRunning()
                       and time.time() - t0 < timeout_s):
                    app.processEvents()
                    time.sleep(0.05)
                for _ in range(20):
                    app.processEvents()
                    time.sleep(0.02)
                if win._fit_worker is not None and win._fit_worker.isRunning():
                    say.fail(f"desktop {model}: the fit thread is still running after {timeout_s:.0f} s")
                    ok = False
                    continue
                after = [dict((k, v["value"]) for k, v in s["params"].items())
                         for s in win.recipe["sites"]]
                moved = any(abs(a[k] - b[k]) > 1e-9 for a, b in zip(after, before) for k in a)
                health = getattr(win, "_health", None)
                say(f"desktop {model}: finished in {time.time() - t0:.1f} s, status "
                    f"{win.statusBar().currentMessage()!r}, parameters moved: {moved}, "
                    f"model curve: {getattr(win, '_last_model', None) is not None}, "
                    f"health: {health.summary() if health is not None else None}")
                if dialogs:
                    say("desktop dialogs:", *dialogs)
                if dialogs or not moved:
                    say.fail(f"desktop {model}: a dialog appeared or nothing changed")
                    ok = False
                dialogs.clear()
            except Exception as exc:                      # noqa: BLE001
                say.fail(f"desktop {model}", exc)
                ok = False
        win.close()
        app.processEvents()
    except Exception as exc:                              # noqa: BLE001
        say.fail("desktop window", exc)
        ok = False
    return ok


def run(gui: bool = True, spectrum: str | None = None) -> int:
    path = log_path()
    say = _Log(path)
    try:
        _environment(say)
        core_ok = core_fits(say)
        gui_ok = desktop_fits(say, spectrum=spectrum) if gui else True
        verdict = "PASS" if (core_ok and gui_ok) else "FAIL"
        say(f"=== {verdict}: core fits {'ok' if core_ok else 'FAILED'}, "
            f"desktop fits {'ok' if gui_ok else 'FAILED'}  (log: {path})")
        return 0 if verdict == "PASS" else 1
    except Exception as exc:                              # noqa: BLE001
        say.fail("self-test", exc)
        return 1
    finally:
        say.close()


def spectrum_arg(argv) -> str | None:
    """The first argument that is not a flag: the user's own spectrum."""
    return next((a for a in argv if not a.startswith("--")), None)


if __name__ == "__main__":                                 # python -m larmor.selftest [path]
    raise SystemExit(run(gui="--no-gui" not in sys.argv,
                         spectrum=spectrum_arg(sys.argv[1:])))
