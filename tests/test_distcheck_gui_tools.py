"""The desktop workflow stages of the distribution check
(larmor.distcheck.gui_tools) pass in a development install: load, apply a
recipe, fit, every export format, the Plotting studio and the live plot's
exporters, overlays, the series mode, the 2D view, the processing panel,
project save / open, the Explorer and Workspaces docks -- offscreen, in
quick mode, through the same ``larmor.distcheck.run`` the frozen exe calls
(``LARMOR.exe --distcheck --gui --quick``). The environment is set before
Qt loads so no window shows and nothing of the user's session is read."""
import json
import os

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")
os.environ.setdefault("LARMOR_NO_KERNEL_WARM", "1")

from larmor.distcheck import gui_tools  # noqa: E402

NAMES = [name for name, _fn in gui_tools.stages()]


def test_stage_names_are_tool_prefixed_kebab_case_and_unique():
    assert len(NAMES) == 9 and len(NAMES) == len(set(NAMES))
    assert all(n.startswith("tool-") and n == n.lower() and " " not in n and "_" not in n
               for n in NAMES)
    assert "tool-fit-workflow" in NAMES and "tool-plotting-studio" in NAMES


def test_tool_stages_pass_in_a_development_install(tmp_path, monkeypatch):
    from larmor import distcheck

    log = tmp_path / "distcheck.log"
    monkeypatch.setenv("LARMOR_DISTCHECK_LOG", str(log))
    rc = distcheck.run(["--gui", "--quick", "--out", str(tmp_path / "out"),
                        "--only", ",".join(NAMES)])
    text = log.read_text(encoding="utf-8")
    if rc != 0:
        print(text[-12000:])
    assert rc == 0, text[-6000:]
    report = json.loads(log.with_suffix(".json").read_text(encoding="utf-8"))
    assert [s["name"] for s in report["stages"]] == NAMES
    assert all(s["ok"] for s in report["stages"]), text[-6000:]
    assert "=== PASS" in text
    # the exports the frozen 0.15.3 build could not make are on disk (the
    # GUI groups run in a child process of their own, under <out>/tool_stages)
    out = tmp_path / "out"
    if (out / "tool_stages").is_dir():
        out = out / "tool_stages"
    studio = out / "studio"
    assert {p.suffix for p in studio.glob("studio.*")} >= {".png", ".pdf", ".svg", ".eps",
                                                           ".tiff", ".jpg"}
    assert (out / "plot" / "plot.svg").stat().st_size > 200
