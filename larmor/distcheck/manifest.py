"""The import manifest: every module LARMOR reaches and every file it opens,
written from the source tree and checked inside the frozen app.

A windowed PyInstaller build fails in one way the test suite never sees: a
module imported lazily -- inside a function, through importlib or a
registry string, by matplotlib's ``savefig(format=...)`` -- that the
dependency scan did not follow, surfacing on the user's machine as
"No module named ..." (0.15.3: ``matplotlib.backends.backend_svg`` on the
first SVG export from the Plotting studio). This module turns that class
of bug into a check that runs in the frozen app itself:

* :func:`scan_source` walks ``larmor/**/*.py`` with ``ast`` and records
  every module of the package and every import statement anywhere in it
  (module level or not, resolved to real modules), plus :data:`DYNAMIC`
  (modules reached only through strings, listed by the feature that needs
  them) and the resource files the app opens at run time (help pages,
  tutorials, assets, static files, the example data).
* :func:`write_manifest` stores that as ``larmor/static/import_manifest.json``,
  which the spec bundles; ``tests/test_distcheck.py`` fails when the source
  moved on without it (the message says how to regenerate).
* :func:`load_manifest` reads it back -- in the frozen app too -- for the
  imports and resources stages of :mod:`larmor.distcheck`.

Qt-free; ``python -m larmor.distcheck.manifest --write`` regenerates.
"""
from __future__ import annotations

import ast
import importlib
import importlib.util
import json
import sys
from pathlib import Path

__all__ = ["MANIFEST_NAME", "DYNAMIC", "OPTIONAL", "repo_root", "resource_root",
           "scan_source", "write_manifest", "load_manifest", "manifest_path"]

MANIFEST_NAME = "import_manifest.json"

#: modules reached only through strings (importlib, matplotlib's backend
#: registry, pyqtgraph's exporter registry, multiprocessing's spawn), by the
#: feature that needs them -- the dependency scan cannot see these
DYNAMIC = {
    "figure export: PNG / TIFF / JPG through Agg + Pillow, SVG, PDF, EPS": [
        "matplotlib.backends.backend_agg", "matplotlib.backends.backend_svg",
        "matplotlib.backends.backend_pdf", "matplotlib.backends.backend_ps",
        "PIL.Image", "PIL.PngImagePlugin", "PIL.TiffImagePlugin", "PIL.JpegImagePlugin"],
    "the Plotting studio's live canvas": ["matplotlib.backends.backend_qtagg"],
    "plot right-click > Save image / SVG (pyqtgraph exporters)": [
        "pyqtgraph.exporters.ImageExporter", "pyqtgraph.exporters.SVGExporter",
        "pyqtgraph.exporters.CSVExporter"],
    "process pool of the error tools and batch fits (Windows spawn)": [
        "multiprocessing.spawn", "multiprocessing.popen_spawn_win32",
        "multiprocessing.reduction", "concurrent.futures.process"],
    "lmfit / uncertainties internals reached through their registries": [
        "lmfit.models", "lmfit.printfuncs", "lmfit.confidence", "asteval"],
    "scipy paths the fit and the error tools reach lazily": [
        "scipy.optimize._lsq.trf", "scipy.optimize._lsq.least_squares",
        "scipy.stats._continuous_distns", "scipy.special._ufuncs",
        "scipy.interpolate._cubic", "scipy.signal._savitzky_golay",
        "scipy.ndimage._filters", "scipy.linalg._decomp_svd"],
    "mrsimulator model registry": [
        "mrsimulator.models.czjzek", "mrsimulator.spin_system.isotope"],
}

#: modules whose absence from the frozen DESKTOP build is accepted, with the
#: reason the distribution check prints instead of a failure
OPTIONAL = {
    "larmor.app": "the FastAPI / Plotly web mode (`larmor app`) is a server feature, "
                  "not part of the desktop build",
    "fastapi": "web mode only", "uvicorn": "web mode only", "plotly": "web mode only",
    "starlette": "web mode only", "pydantic": "web mode only",
    "multiprocessing.popen_spawn_win32": "Windows only",
    "msvcrt": "Windows only", "winreg": "Windows only",
}

_SKIP_DIRS = {"__pycache__"}


def repo_root() -> Path:
    """The source checkout (``larmor/..``) -- meaningful outside a frozen app."""
    import larmor
    return Path(larmor.__file__).resolve().parents[1]


def resource_root() -> Path:
    """Where the bundled resources live: PyInstaller's ``_MEIPASS`` in the
    frozen app (``larmor/help``, ``docs/tutorials``, ``assets``, ``examples``
    keep their source-tree layout under it), the checkout otherwise."""
    frozen = getattr(sys, "_MEIPASS", "")
    return Path(frozen) if frozen else repo_root()


def manifest_path() -> Path:
    import larmor
    return Path(larmor.__file__).resolve().parent / "static" / MANIFEST_NAME


# ------------------------------------------------------------------ scanning
def _module_name(pkg_root: Path, path: Path) -> str:
    rel = path.relative_to(pkg_root.parent).with_suffix("")
    parts = list(rel.parts)
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join(parts)


def _import_candidates(path: Path) -> set:
    """Every module name an import statement in ``path`` could mean, at any
    depth of the file: ``import a.b`` -> ``a.b``; ``from a import b`` ->
    ``a`` and the candidate ``a.b`` (a submodule or an attribute -- the
    scan resolves which)."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    out: set = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                out.add(a.name)
        elif isinstance(node, ast.ImportFrom):
            if node.level or not node.module:
                continue                          # relative: inside the package
            out.add(node.module)
            for a in node.names:
                if a.name != "*":
                    out.add(f"{node.module}.{a.name}")
    return out


def _is_module(name: str) -> bool | None:
    """True for an importable module, False for an attribute of its parent,
    None when the top-level package is not installed here."""
    top = name.split(".")[0]
    try:
        if importlib.util.find_spec(top) is None:
            return None
    except (ImportError, ValueError):
        return None
    try:
        spec = importlib.util.find_spec(name)
    except (ImportError, ValueError, AttributeError):
        spec = None
    if spec is not None:
        return True
    parent, _, leaf = name.rpartition(".")
    if not parent:
        return False
    try:
        mod = importlib.import_module(parent)
    except Exception:                                 # noqa: BLE001
        return None
    return False if hasattr(mod, leaf) else None


def _category(name: str) -> str:
    top = name.split(".")[0]
    if top == "larmor":
        return "larmor"
    if top in getattr(sys, "stdlib_module_names", ()):
        return "stdlib"
    return "external"


def _files_under(base: Path, exclude_dirs=()) -> list:
    if not base.is_dir():
        return []
    out = []
    for p in sorted(base.rglob("*")):
        if p.is_file() and not any(part in exclude_dirs or part in _SKIP_DIRS
                                   for part in p.relative_to(base).parts):
            out.append(p.relative_to(base).as_posix())
    return out


def scan_source(root: Path | None = None) -> dict:
    """Build the manifest from the checkout at ``root`` (default: this one)."""
    root = Path(root) if root else repo_root()
    pkg = root / "larmor"
    files = sorted(p for p in pkg.rglob("*.py")
                   if not any(part in _SKIP_DIRS for part in p.parts))
    modules = sorted(_module_name(pkg, p) for p in files)
    imported_by: dict = {}
    for p in files:
        mod = _module_name(pkg, p)
        for cand in _import_candidates(p):
            imported_by.setdefault(cand, set()).add(mod)
    imports: dict = {"larmor": [], "stdlib": [], "external": []}
    unavailable: dict = {}
    where: dict = {}
    for name in sorted(imported_by):
        if name.startswith("larmor"):
            kind = True if name in modules else _is_module(name)
        else:
            kind = _is_module(name)
        if kind is False:
            continue                                  # an attribute, not a module
        if kind is None:
            unavailable[name] = sorted(imported_by[name])
            continue
        imports[_category(name)].append(name)
        where[name] = sorted(imported_by[name])
    for name in modules:
        if name not in imports["larmor"]:
            imports["larmor"].append(name)
    imports["larmor"] = sorted(set(imports["larmor"]))
    # a DYNAMIC name that does not exist in this environment is a typo (or a
    # library that renamed a private module): say so here, not in the frozen app
    bad = [n for names in DYNAMIC.values() for n in names
           if n not in OPTIONAL and _is_module(n) is not True]
    if bad:
        raise ValueError(f"DYNAMIC names that are not modules here: {bad}")
    dynamic = {k: sorted(v) for k, v in DYNAMIC.items()}
    resources = {
        "help": _files_under(pkg / "help"),
        "tutorials": _files_under(root / "docs" / "tutorials"),
        "assets": _files_under(root / "assets"),
        "static": [f for f in _files_under(pkg / "static") if f != MANIFEST_NAME],
        "examples": _files_under(root / "examples", exclude_dirs=("bib",)),
        "xfact": _files_under(pkg / "xfact" / "packs" / "birds" / "assets"),
        "docs": [d for d in ("README.md",) if (root / d).is_file()],
    }
    return {"modules": modules, "imports": imports, "imported_by": where,
            "unavailable_in_dev": unavailable, "dynamic": dynamic,
            "optional": dict(OPTIONAL), "resources": resources}


def write_manifest(root: Path | None = None, path: Path | None = None) -> Path:
    data = scan_source(root)
    path = Path(path) if path else manifest_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    return path


def load_manifest(path: Path | None = None) -> dict:
    path = Path(path) if path else manifest_path()
    return json.loads(path.read_text(encoding="utf-8"))


def _main(argv=None) -> int:
    import argparse

    ap = argparse.ArgumentParser(prog="python -m larmor.distcheck.manifest",
                                 description="write or show the import manifest")
    ap.add_argument("--write", action="store_true", help="regenerate larmor/static/"
                    + MANIFEST_NAME + " from this checkout")
    args = ap.parse_args(argv)
    if args.write:
        p = write_manifest()
        m = load_manifest(p)
        print(f"wrote {p}: {len(m['modules'])} larmor modules, "
              f"{len(m['imports']['external'])} external and "
              f"{len(m['imports']['stdlib'])} stdlib imports, "
              f"{sum(len(v) for v in m['dynamic'].values())} dynamic, "
              f"{sum(len(v) for v in m['resources'].values())} resource files"
              + (f"; not installed here: {sorted(m['unavailable_in_dev'])}"
                 if m["unavailable_in_dev"] else ""))
        return 0
    m = scan_source()
    print(json.dumps({k: (len(v) if isinstance(v, (list, dict)) else v)
                      for k, v in m.items()}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
