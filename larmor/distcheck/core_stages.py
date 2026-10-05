"""The core stages of the distribution check (Qt-free): environment,
imports from the manifest, bundled resources, figure export formats, the
process pool, the command line."""
from __future__ import annotations

import contextlib
import importlib
import io
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np

from larmor.distcheck.manifest import load_manifest, resource_root

__all__ = ["STAGES"]


# ------------------------------------------------------------ environment
def environment(say, ctx):
    import larmor
    say("python", sys.version.split()[0], "| executable", sys.executable)
    say("platform", platform.platform(), "| cpu count", os.cpu_count())
    say("larmor", larmor.__version__, "from", Path(larmor.__file__).parent)
    say("resources from", resource_root())
    if ctx["frozen"]:
        say("_MEIPASS", getattr(sys, "_MEIPASS", ""))
    for name in ("numpy", "scipy", "lmfit", "mrsimulator", "csdmpy", "nmrglue",
                 "matplotlib", "pyqtgraph", "PySide6", "PIL", "pybaselines",
                 "uncertainties", "asteval", "requests"):
        try:
            mod = importlib.import_module(name)
            say(f"  {name} {getattr(mod, '__version__', '?')}")
        except Exception as exc:                          # noqa: BLE001
            say(f"  {name}: IMPORT FAILED: {exc!r}")
    say("stdout present", sys.stdout is not None, "| stderr present", sys.stderr is not None)


# ------------------------------------------------------------ imports
def imports(say, ctx):
    """Import every module the manifest names; a missing one is the bug the
    frozen build hides until a user reaches that feature."""
    m = load_manifest()
    optional = m.get("optional", {})
    groups = [("larmor", m["imports"]["larmor"]), ("external", m["imports"]["external"]),
              ("stdlib", m["imports"]["stdlib"])]
    for feature, names in m.get("dynamic", {}).items():
        groups.append((f"dynamic · {feature}", names))
    missing: list = []
    broken: list = []
    skipped: list = []
    n_ok = 0
    for group, names in groups:
        bad_here = []
        for name in names:
            try:
                importlib.import_module(name)
                n_ok += 1
            except ImportError as exc:
                # an accepted absence: the optional module itself, or a module
                # whose import fails because an optional one is missing
                # (larmor.app raises its own ImportError naming the extra)
                culprit = getattr(exc, "name", None)
                if name in optional or culprit in optional:
                    skipped.append(f"{name} ({optional.get(name) or optional.get(culprit)})")
                    continue
                if (sys.platform != "win32" and name in
                        ("multiprocessing.popen_spawn_win32", "msvcrt", "winreg")):
                    skipped.append(f"{name} (Windows only)")
                    continue
                bad_here.append(f"{name}: {exc}")
                (missing if isinstance(exc, ModuleNotFoundError) else broken).append(name)
            except Exception as exc:                      # noqa: BLE001
                bad_here.append(f"{name}: {type(exc).__name__}: {exc}")
                broken.append(name)
        say(f"{group}: {len(names)} named, {len(bad_here)} failing")
        for line in bad_here:
            say("   ", line)
    if skipped:
        say("accepted absences:", "; ".join(skipped))
    say(f"imported {n_ok} modules")
    if missing or broken:
        raise AssertionError(
            f"{len(missing)} module(s) missing from this installation: {missing}"
            + (f"; {len(broken)} failing to import: {broken}" if broken else ""))


# ------------------------------------------------------------ resources
def resources(say, ctx):
    """Every file the app opens at run time is where the frozen build
    expects it (help pages, tutorials, assets, static files, examples)."""
    m = load_manifest()
    root = resource_root()
    import larmor
    pkg = Path(larmor.__file__).resolve().parent
    bases = {"help": pkg / "help", "tutorials": root / "docs" / "tutorials",
             "assets": root / "assets", "static": pkg / "static",
             "examples": root / "examples", "xfact": pkg / "xfact" / "packs" / "birds" / "assets",
             "docs": root}
    missing = []
    for kind, files in m["resources"].items():
        base = bases[kind]
        lost = [f for f in files if not (base / f).is_file()]
        say(f"{kind}: {len(files) - len(lost)} of {len(files)} present under {base}")
        missing += [f"{kind}/{f}" for f in lost]
    # the help viewer's own locators must agree
    try:
        from larmor.desktop.help_dialog import help_path, tutorial_path
    except Exception as exc:                              # noqa: BLE001
        say("help locators not importable:", repr(exc))
    else:
        for f in m["resources"]["help"]:
            if help_path(Path(f).stem) is None:
                missing.append(f"help locator: {f}")
        for f in m["resources"]["tutorials"]:
            if tutorial_path(Path(f).stem) is None:
                missing.append(f"tutorial locator: {f}")
    if missing:
        raise AssertionError(f"{len(missing)} resource(s) missing: {missing}")


# ------------------------------------------------------------ exports
def _spec_1d() -> dict:
    x = np.linspace(-20, 60, 400)
    y = 100 * np.exp(-0.5 * ((x - 15) / 4) ** 2) + 30 * np.exp(-0.5 * ((x - 2) / 3) ** 2)
    return {"kind": "1d", "x_is_ppm": True, "xlabel": "shift (ppm)", "ylabel": "a.u.",
            "traces": [{"data": {"x": x.tolist(), "y": y.tolist()}, "label": "synthetic",
                        "color": "#222222"}]}


def exports(say, ctx):
    """Every figure format the export dialogs offer, through the same calls:
    PNG / SVG / PDF / EPS via larmor.figures.export, TIFF / JPG through
    Agg + Pillow, PNG bytes for the help renderer."""
    from larmor import figures
    out = Path(ctx["out"]) / "exports"
    out.mkdir(parents=True, exist_ok=True)
    saved = figures.export(_spec_1d(), out / "fig", formats=("png", "svg", "pdf", "eps", "json"),
                           dpi=120)
    for p in saved:
        size = Path(p).stat().st_size
        say(f"{Path(p).suffix[1:]}: {size} bytes")
        if size < 200:
            raise AssertionError(f"{p} is suspiciously small ({size} bytes)")
    fig = figures.render(_spec_1d())
    for fmt in ("tiff", "jpg"):
        target = out / f"fig.{fmt}"
        fig.savefig(target, dpi=100, bbox_inches="tight")
        say(f"{fmt}: {target.stat().st_size} bytes")
    import matplotlib.pyplot as plt
    plt.close(fig)
    png = figures.render_png_bytes({"kind": "species_bar", "categories": ["a", "b"],
                                    "series": [{"label": "s0", "values": [60.0, 40.0]},
                                               {"label": "s1", "values": [40.0, 60.0]}]})
    say(f"species bar PNG bytes: {len(png)}")
    if not png.startswith(b"\x89PNG"):
        raise AssertionError("render_png_bytes did not return a PNG")


# ------------------------------------------------------------ pool
def _square(x):
    return int(x) * int(x), os.getpid()


def pool(say, ctx):
    """The process pool of the error tools runs here, in real worker
    processes: a frozen exe must re-enter its launcher as a worker (not as
    a second LARMOR), and ``parallel_map``'s sequential fallback for a
    broken pool must not be what passed this stage -- the worker pids are
    checked against this process's."""
    import warnings
    from larmor import parallel
    # above parallel.MIN_ITEMS_FOR_PROCESSES (8): fewer items run in-process
    # by design, which is not what this stage is about
    items = list(range(12))
    t0 = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        got = parallel.parallel_map(_square, items, max_workers=2)
    dt = time.perf_counter() - t0
    try:
        parallel.shutdown_shared_pool()
    except Exception:                                     # noqa: BLE001
        pass
    notes = [str(w.message) for w in caught if issubclass(w.category, RuntimeWarning)]
    for n in notes:
        say("runtime warning:", n)
    if any(r is None for r in got):
        raise AssertionError(f"pool left holes in the results: {got}")
    values = [r[0] for r in got]
    pids = sorted({r[1] for r in got})
    say(f"{len(items)} items through a 2-worker pool in {dt:.1f} s: {values}; worker pids "
        f"{pids}, this process {os.getpid()}")
    if values != [i * i for i in items]:
        raise AssertionError(f"pool results wrong: {values}")
    if os.getpid() in pids:
        raise AssertionError("the pool ran in THIS process: the worker processes did not "
                             "start (sequential fallback)" + (f": {notes}" if notes else ""))


# ------------------------------------------------------------ cli
def cli(say, ctx):
    """The command line, in-process: --help of every subcommand, then info /
    compare / acqtable / inventory on the bundled example data."""
    from larmor import cli as climod
    subs = ["info", "srcheck", "compare", "inventory", "acqtable", "import", "fit",
            "satrec", "redor", "magres", "shiftcal", "multifit", "batchfit", "seqfit"]
    for sub in subs:
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
                climod.main([sub, "--help"])
        except SystemExit as exc:
            if exc.code not in (0, None):
                raise AssertionError(f"larmor {sub} --help exited {exc.code}: {buf.getvalue()[-300:]}")
        if "usage" not in buf.getvalue().lower():
            raise AssertionError(f"larmor {sub} --help printed no usage: {buf.getvalue()[:200]!r}")
    say(f"--help of {len(subs)} subcommands ok")
    ex = ctx.get("examples")
    if not ex:
        say("no bundled example data: the data commands are not exercised here")
        return
    expnos = sorted(p for p in (ex / "pCABS2-4").iterdir() if p.is_dir() and p.name.isdigit())
    one = str(expnos[0])
    csv_out = str(Path(ctx["out"]) / "acquisition.csv")   # acqtable's default is the cwd
    for argv in (["info", one], ["acqtable", *map(str, expnos[:2]), "-o", csv_out],
                 ["compare", *map(str, expnos[:2])], ["inventory", str(ex / "pCABS2-4")]):
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
                rc = climod.main(argv)
        except SystemExit as exc:
            rc = exc.code
        text = buf.getvalue()
        say(f"larmor {' '.join(Path(a).name if os.sep in a else a for a in argv)}: "
            f"exit {rc}, {len(text.splitlines())} lines")
        # `larmor compare` exits 2 when the spectra were NOT acquired alike
        # (a lint-style verdict, not an error) -- the two example EXPNOs differ
        allowed = (0, None, 2) if argv[0] == "compare" else (0, None)
        if rc not in allowed:
            raise AssertionError(f"larmor {argv[0]} exited {rc}: {text[-400:]}")


STAGES = [("environment", environment), ("imports", imports), ("resources", resources),
          ("exports", exports), ("pool", pool), ("cli", cli)]
