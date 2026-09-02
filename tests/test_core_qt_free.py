"""The Qt-free core is a load-bearing convention -- now enforced.

Everything outside larmor/desktop/ and larmor/xfact/ must import without
pulling PySide6 (or any Qt binding): batch fitting, the CLI, and worker
processes all rely on it, and one careless `from larmor.desktop...` import in
a core module would break headless use with no test noticing. Verified the
way the development notes describe: import every core module in one clean
subprocess and look at sys.modules.
"""
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

#: modules excused from the sweep: the unmaintained browser variant needs
#: fastapi (an optional extra), nothing else may fail to import at all
_OPTIONAL = {"larmor.app"}

_PROBE = r"""
import importlib, json, sys

mods, failed = json.loads(sys.argv[1]), {}
for m in mods:
    try:
        importlib.import_module(m)
    except ImportError as exc:
        failed[m] = str(exc)
qt = sorted(m for m in sys.modules
            if m.split(".")[0] in ("PySide6", "PyQt5", "PyQt6", "shiboken6"))
print(json.dumps({"qt": qt, "failed": failed}))
"""


def _core_modules() -> list[str]:
    mods = []
    for py in sorted((REPO / "larmor").rglob("*.py")):
        rel = py.relative_to(REPO)
        parts = rel.with_suffix("").parts
        if "desktop" in parts or "xfact" in parts or "__pycache__" in parts:
            continue
        name = ".".join(parts[:-1] if parts[-1] == "__init__" else parts)
        mods.append(name)
    return mods


def test_core_imports_are_qt_free():
    mods = _core_modules()
    assert len(mods) > 50, f"module sweep looks wrong: {len(mods)} found"
    out = subprocess.run(
        [sys.executable, "-c", _PROBE, json.dumps(mods)],
        capture_output=True, text=True, cwd=str(REPO), timeout=300)
    assert out.returncode == 0, out.stderr[-2000:]
    res = json.loads(out.stdout.strip().splitlines()[-1])
    hard_failures = {m: e for m, e in res["failed"].items()
                     if m not in _OPTIONAL}
    assert not hard_failures, f"core modules failed to import: {hard_failures}"
    assert not res["qt"], (
        "importing the Qt-free core loaded Qt modules -- some core module "
        f"now imports desktop code: {res['qt']}")
