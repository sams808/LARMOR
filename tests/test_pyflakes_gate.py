"""pyflakes must report nothing on the package and on the suite (E6).

The forty-odd cosmetic findings (unused imports, f-strings without
placeholders, results assigned and dropped) were curated once so that a REAL
finding -- an undefined name, a shadowed import, a fit result silently thrown
away -- is visible the moment it appears instead of buried in noise. Deliberate
re-exports are declared through ``__all__``; imports made only for their side
effects go through ``importlib.import_module``.
"""
import io
from pathlib import Path

import pytest

pytest.importorskip(
    "pyflakes", reason="pyflakes is a dev dependency: pip install larmor[dev]")

ROOT = Path(__file__).resolve().parent.parent
TREES = ("larmor", "tests")


def _sources():
    for tree in TREES:
        for p in sorted((ROOT / tree).rglob("*.py")):
            if "__pycache__" not in p.parts:
                yield p


def test_pyflakes_is_clean():
    from pyflakes.api import checkPath
    from pyflakes.reporter import Reporter

    out, err = io.StringIO(), io.StringIO()
    reporter = Reporter(out, err)
    count = sum(checkPath(str(p), reporter) for p in _sources())
    assert count == 0, (
        f"pyflakes reported {count} finding(s); fix them or declare the "
        f"re-export in __all__:\n{out.getvalue()}{err.getvalue()}")


def test_gate_actually_scans_the_trees():
    """A gate that scans nothing passes for the wrong reason."""
    files = list(_sources())
    assert len(files) > 100, len(files)
    assert any(p.name == "app.py" for p in files)
    assert any(p.name == "conftest.py" for p in files)
