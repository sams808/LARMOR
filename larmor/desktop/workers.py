"""Background threads of the main window, and the lmfit callback glue.

``FitWorker`` / ``Fit2DWorker`` run a 1D / 2D fit off the GUI thread and
report progress through Qt signals; ``SimWorker`` runs one simulation;
``KernelWarmWorker`` pre-builds the Czjzek kernel for a freshly loaded
dataset. ``_emit_progress`` is the lmfit ``iter_cb`` that feeds the progress
bar and honours Stop / Cancel, ``humanize_error`` turns an lmfit / numpy
message into one a user can act on, and ``_fit_tol`` reads the saved
tolerance. ``MainWindow`` (``larmor.desktop.app``) re-exports every public
name here, so ``from larmor.desktop.app import FitWorker`` keeps working.
"""
from __future__ import annotations

import numpy as np
from PySide6.QtCore import QThread, Signal

from larmor.recipe import Recipe


def _emit_progress(sig, should_stop=None, converge_frac=None):
    """An lmfit iter_cb that reports (iteration, residual stdev) to a Qt signal.

    Returning True aborts the minimisation. It halts if ``should_stop()`` is
    truthy, or — dmfit-style — once the residual stdev changes by less than
    ``converge_frac`` (a fraction, e.g. 1e-3 = 0.1 %) between iterations, keeping
    the latest parameters."""
    state = {"n": 0, "prev": None}

    def cb(params, it, resid, *args, **kws):
        state["n"] += 1
        try:
            rms = float(np.sqrt(np.mean(np.asarray(resid, float) ** 2)))
        except Exception:
            rms = float("nan")
        sig.emit(state["n"], rms)
        if should_stop is not None and should_stop():
            return True
        prev = state["prev"]
        if (converge_frac and prev is not None and prev > 0 and rms == rms
                and abs(rms - prev) / prev < converge_frac):
            return True                      # sdev not changing more than threshold
        state["prev"] = rms
        return None
    return cb


def humanize_error(msg) -> str:
    """Turn a cryptic engine/Qt exception into a plain-language explanation the
    user can act on (keeps the original text as a tail for reference)."""
    m = str(msg)
    low = m.lower()
    if "recursion" in low:
        return ("A line's parameter is linked to itself (or two lines link to "
                "each other in a loop), so the fit can't evaluate it. Check your "
                "position/width links — the “+/−” entries must reference ANOTHER "
                "line, not their own line. LARMOR removes such links before "
                "fitting; if this persists, re-open the recipe and re-check the "
                "links.")
    if "singular" in low or ("nan" in low and "residual" in low):
        return ("The fit became numerically unstable (a singular/NaN residual). "
                "Usually a parameter hit a bound or a line collapsed to zero "
                "amplitude — widen the bounds or remove a redundant line.\n\n"
                f"(details: {m})")
    if "at least two" in low:
        return m
    if "no experimental data" in low:
        return ("This dmfit file holds a model but no spectrum — open it onto a "
                "spectrum (File ▸ Open a spectrum first, then Decomposition ▸ "
                "Apply recipe) or pick data when prompted.")
    return m


def _fit_tol():
    """The user's global completion threshold (% change in the residual stdev at
    which a fit is considered done). Default 0.1 % ≈ dmfit's 1.0e-3; set to 0 for
    the full-precision solver default. Honoured by every fit button in the app."""
    from PySide6.QtCore import QSettings
    try:
        return float(QSettings("LARMOR", "app").value("fitStdevPct", 0.1) or 0.0)
    except (TypeError, ValueError):
        return 0.1


class _StoppableFit:
    """Mixin: a stop flag + mode ('stop' keep / 'cancel' revert) for a fit thread."""

    def _init_stop(self):
        self._stop = False
        self._stop_mode = ""

    def request_stop(self, mode: str):
        self._stop_mode = mode          # "stop" (keep) | "cancel" (revert)
        self._stop = True


class FitWorker(QThread, _StoppableFit):
    done = Signal(object, str)          # (result, stop_mode)
    failed = Signal(str)
    progress = Signal(int, float)                # (iteration, residual rms)
    frame = Signal(object, object, int)          # (x_ppm, y_model, iteration)
    kernel_progress = Signal(int, int)           # (chunk done, chunks total)

    def __init__(self, recipe_dict, ppm, amp, window, animate=False):
        super().__init__()
        self.recipe_dict, self.ppm, self.amp, self.window = \
            recipe_dict, ppm, amp, window
        self.animate = animate
        self._init_stop()

    def run(self):
        try:
            from larmor import engine
            from larmor import fit as fitmod

            recipe = Recipe.from_dict(self.recipe_dict)
            frame_cb = ((lambda x, y, it: self.frame.emit(x, y, it))
                        if self.animate else None)
            # a cold wideline kernel build (~23 s) now reports its chunks
            # and honours the same Stop button as the fit itself
            with engine.kernel_build_feedback(
                    cb=lambda i, n: self.kernel_progress.emit(i, n),
                    cancel=lambda: self._stop):
                result = fitmod.fit(recipe, self.ppm, self.amp,
                                    window_ppm=self.window, tol=_fit_tol(),
                                    iter_cb=_emit_progress(
                                        self.progress, lambda: self._stop,
                                        converge_frac=(_fit_tol() / 100.0)
                                        or None),
                                    frame_cb=frame_cb)
            self.done.emit(result, self._stop_mode)
        except engine.KernelBuildCancelled:
            self.failed.emit("cancelled while building the lineshape kernel")
        except Exception as exc:
            self.failed.emit(str(exc))


class Fit2DWorker(QThread, _StoppableFit):
    done = Signal(object, str)
    failed = Signal(str)
    progress = Signal(int, float)

    def __init__(self, recipe_dict, data2d, method):
        super().__init__()
        self.recipe_dict, self.data2d, self.method = \
            recipe_dict, data2d, method
        self._init_stop()

    def run(self):
        try:
            from larmor import twod

            recipe = Recipe.from_dict(self.recipe_dict)
            result = twod.fit_2d(recipe, self.data2d, method=self.method,
                                 tol=_fit_tol(),
                                 iter_cb=_emit_progress(
                                     self.progress, lambda: self._stop,
                                     converge_frac=(_fit_tol() / 100.0) or None))
            self.done.emit(result, self._stop_mode)
        except Exception as exc:
            self.failed.emit(str(exc))


class KernelWarmWorker(QThread):
    """Pre-build the default Czjzek kernel in the background while the user
    is still placing lines. The nucleus, field, spin rate and window are all
    known the moment a spectrum opens; waiting until the first Simulate
    meant eating the whole build (up to ~23 s wideline, cold) right at the
    click. Failures are silent -- this is an optimisation, never a feature.
    Worst case is a race with an explicit early Simulate building the same
    kernel twice; both land in the cache, last write wins, nothing wrong."""

    def __init__(self, nucleus, larmor_MHz, spin_rate_Hz, exp_ppm):
        super().__init__()
        self._args = (nucleus, float(larmor_MHz), float(spin_rate_Hz),
                      np.asarray(exp_ppm, float))

    def run(self):
        try:
            from larmor import engine
            from larmor.nuclei import all_isotopes

            nucleus, lar, spin_hz, exp_ppm = self._args
            iso = next((i for i in all_isotopes() if i.symbol == nucleus),
                       None)
            if iso is None or iso.spin < 1.0 or lar <= 0:
                return                     # spin-1/2 (or unknown): no kernel
            sw, ref = engine.kernel_window(exp_ppm, lar)
            npts = min(int(engine.KERNEL_SETTINGS["npts"]
                           * max(1.0, sw / engine.KERNEL_MIN_SW_HZ)), 16384)
            # exactly the kernel a fresh czjzek site (default sigma 2 MHz,
            # 10 sigma headroom -> ladder step 25) asks for on its first render
            from larmor.models.quadrupolar import CZJZEK_KERNEL_HEADROOM
            engine.build_kernel(
                nucleus, lar, spin_hz, sw_Hz=sw, npts=npts,
                ref_offset_ppm=ref,
                cq_max_MHz=engine.kernel_cq_max(CZJZEK_KERNEL_HEADROOM * 2.0),
                n_cq=engine.KERNEL_SETTINGS["n_cq"],
                n_eta=int(engine.KERNEL_SETTINGS["n_eta"]))
        except Exception:                                  # noqa: BLE001
            pass


class SimWorker(QThread):
    done = Signal(object, object, object)
    failed = Signal(str)
    kernel_progress = Signal(int, int)           # (chunk done, chunks total)

    def __init__(self, recipe_dict, exp_ppm):
        super().__init__()
        self.recipe_dict, self.exp_ppm = recipe_dict, exp_ppm

    def run(self):
        try:
            from larmor import engine
            from larmor import fit as fitmod

            recipe = Recipe.from_dict(self.recipe_dict)
            params = fitmod._make_params(recipe)
            fitmod._apply_params(recipe, params)
            with engine.kernel_build_feedback(
                    cb=lambda i, n: self.kernel_progress.emit(i, n)):
                x, total, per_site = engine.simulate(recipe,
                                                     exp_ppm=self.exp_ppm)
            self.done.emit(x, total, per_site)
        except Exception as exc:
            self.failed.emit(str(exc))
