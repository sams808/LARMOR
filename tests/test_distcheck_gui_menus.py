"""The menu sweeps of the distribution check (larmor.distcheck.gui_menus)
pass in a development install, offscreen: the window builds and loads the
example, every menu action runs in three window states, every theme, help
page and View toggle renders. The frozen app runs the same stages through
``LARMOR.exe --distcheck --gui``.

The whole set is one run (the stages share nothing but the scratch folder);
the log tail is the failure message. The pinned leaf count the sweep checks
itself against is held to the golden menu tree of tests/test_app_split.py.
"""
import json
import os
from pathlib import Path

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")
os.environ.setdefault("LARMOR_NO_KERNEL_WARM", "1")

ROOT = Path(__file__).resolve().parents[1]

STAGES = ["menu-window", "menu-sweep", "menu-themes", "menu-help", "menu-panels"]


def test_stage_names_and_the_pinned_leaf_count():
    from larmor.distcheck import all_stages, gui_menus
    from test_app_split import GOLDEN_MENU

    names = [n for n, _ in gui_menus.stages()]
    assert names == STAGES
    assert all(n.startswith("menu-") for n in names)
    every = [n for n, _ in all_stages(gui=True)]
    assert len(every) == len(set(every)), "stage names must be unique"
    assert set(names) <= set(every)
    # the sweep's own sanity floor is the golden tree's leaf count: a menu
    # change that regenerates GOLDEN_MENU must update the constant too
    leaves = sum(1 for row in GOLDEN_MENU if not row[4])
    assert gui_menus.EXPECTED_LEAF_ACTIONS == leaves, (
        f"gui_menus.EXPECTED_LEAF_ACTIONS is {gui_menus.EXPECTED_LEAF_ACTIONS}, the golden "
        f"menu tree has {leaves} leaf actions -- update the constant")


def test_menu_stages_pass_in_a_development_install(tmp_path, monkeypatch):
    from larmor import distcheck

    log = tmp_path / "distcheck.log"
    monkeypatch.setenv("LARMOR_DISTCHECK_LOG", str(log))
    rc = distcheck.run(["--gui", "--quick", "--out", str(tmp_path / "out"),
                        "--only", ",".join(STAGES)])
    text = log.read_text(encoding="utf-8")
    if rc != 0:
        print(text[-12000:])
    assert rc == 0, text[-6000:]
    report = json.loads(log.with_suffix(".json").read_text(encoding="utf-8"))
    assert [s["name"] for s in report["stages"]] == STAGES
    assert all(s["ok"] for s in report["stages"]), [s["name"] for s in report["stages"]
                                                    if not s["ok"]]
    assert "=== PASS" in text
    sweep = next(s for s in report["stages"] if s["name"] == "menu-sweep")
    notes = "\n".join(sweep["notes"])
    # the menu bar was really walked, in more than one window state
    assert "[loaded: " in notes and "[empty: " in notes
    assert "triggered (0 failing)" in notes
    assert "skipped File ▸ Quit" in notes
    print("menu stages: " + ", ".join(f"{s['name']} {s['seconds']:.0f} s"
                                      for s in report["stages"]))
