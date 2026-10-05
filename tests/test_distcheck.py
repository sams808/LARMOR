"""The distribution check (larmor.distcheck): the import manifest is current,
the help pages link to pages that exist, and the core stages pass in a
development install (the frozen app runs the same stages through
LARMOR.exe --distcheck)."""
import json
import os
import re
from pathlib import Path

import pytest

from larmor.distcheck import manifest as mf

ROOT = Path(__file__).resolve().parents[1]


def test_manifest_is_current():
    """larmor/static/import_manifest.json must describe THIS source tree --
    a new module, a new import or a new help page without a regenerated
    manifest would leave the frozen check blind to it."""
    fresh = mf.scan_source(ROOT)
    stored = mf.load_manifest()
    for key in ("modules", "imports", "dynamic", "resources", "optional"):
        assert stored.get(key) == fresh[key], (
            f"the import manifest is stale ({key} differs) -- regenerate it with:\n"
            "    python -m larmor.distcheck.manifest --write")
    assert len(fresh["modules"]) >= 170
    assert "larmor.distcheck.core_stages" in fresh["modules"]
    # the lazy imports the 0.15.3 exe missed are named, by feature
    assert "matplotlib.backends.backend_svg" in sum(fresh["dynamic"].values(), [])
    assert fresh["resources"]["help"] and fresh["resources"]["tutorials"]
    assert "pCABS2-4/3616/pdata/1/1r" in fresh["resources"]["examples"]


def test_help_pages_link_to_pages_that_ship():
    pages = {p.stem for p in (ROOT / "larmor" / "help").glob("*.md")}
    tutorials = {p.stem for p in (ROOT / "docs" / "tutorials").glob("*.md")}
    bad = []
    for p in sorted((ROOT / "larmor" / "help").glob("*.md")) + \
            sorted((ROOT / "docs" / "tutorials").glob("*.md")):
        for target in re.findall(r"\]\(([^)#\s]+\.md)", p.read_text(encoding="utf-8")):
            stem = Path(target).stem
            if stem not in pages and stem not in tutorials:
                bad.append(f"{p.name} -> {target}")
    assert not bad, bad


def test_core_stages_pass_in_a_development_install(tmp_path, monkeypatch):
    from larmor import distcheck

    log = tmp_path / "distcheck.log"
    monkeypatch.setenv("LARMOR_DISTCHECK_LOG", str(log))
    rc = distcheck.run(["--no-gui", "--quick", "--out", str(tmp_path / "out"),
                        "--only", "environment,imports,resources,exports,pool,cli"])
    text = log.read_text(encoding="utf-8")
    assert rc == 0, text[-3000:]
    report = json.loads(log.with_suffix(".json").read_text(encoding="utf-8"))
    assert [s["name"] for s in report["stages"]] == [
        "environment", "imports", "resources", "exports", "pool", "cli"]
    assert all(s["ok"] for s in report["stages"])
    assert "=== PASS" in text
    exported = {p.suffix for p in (tmp_path / "out" / "exports").iterdir()}
    assert {".png", ".svg", ".pdf", ".eps", ".tiff", ".jpg", ".json"} <= exported


def test_list_and_unknown_stage_names():
    from larmor import distcheck

    names = [n for n, _ in distcheck.all_stages(gui=False)]
    assert names[:6] == ["environment", "imports", "resources", "exports", "pool", "cli"]
    assert len(names) == len(set(names))
    log = os.environ.get("LARMOR_DISTCHECK_LOG")
    assert log is None or Path(log).parent.exists()
