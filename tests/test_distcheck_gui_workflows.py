"""The dialog-tool stages of the distribution check (larmor.distcheck.
gui_workflows) pass in a development install, offscreen, in quick mode.

The orchestrator does not list these stages yet, so this reproduces its loop
(``larmor.distcheck.run``): a scratch folder, a ``say`` that collects the
notes, the same ``ctx`` (the bundled examples when present), every stage
run in turn, the failures collected and the whole log printed when one
fails. The environment is set BEFORE PySide6 is imported, as ``run`` does:
no window ever shows, native dialogs become widgets the sweep can dismiss,
nothing of the user's session is read, no kernel pre-build thread starts.
"""
import os
import time
import traceback
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("LARMOR_NO_SESSION", "1")
os.environ.setdefault("LARMOR_NO_KERNEL_WARM", "1")

pytest.importorskip("PySide6")

ROOT = Path(__file__).resolve().parents[1]


def test_stage_names_are_dialog_prefixed_and_unique():
    from larmor.distcheck import gui_workflows

    names = [n for n, _fn in gui_workflows.stages()]
    assert names == ["dialog-error-tools", "dialog-batch", "dialog-lines-and-literature",
                     "dialog-relaxation-wideline", "dialog-tools-misc"]
    assert len(names) == len(set(names))
    assert all(n.startswith("dialog-") and n == n.lower() and " " not in n for n in names)
    assert all(callable(fn) for _n, fn in gui_workflows.stages())


def test_dialog_stages_pass_in_a_development_install(tmp_path):
    from larmor.distcheck import gui_workflows

    examples = ROOT / "examples"
    ctx = {"out": tmp_path / "out", "examples": examples if (examples / "pCABS2-4").is_dir() else None,
           "quick": True, "gui": True, "frozen": False}
    ctx["out"].mkdir()
    log: list = []
    failed: list = []
    t_all = time.perf_counter()
    for name, fn in gui_workflows.stages():
        log.append(f"--- {name}")
        t0 = time.perf_counter()
        try:
            fn(lambda *parts: log.append("    " + " ".join(str(p) for p in parts)), ctx)
        except Exception:                                  # noqa: BLE001
            log.append(traceback.format_exc())
            failed.append(name)
            log.append(f"FAIL {name} ({time.perf_counter() - t0:.1f} s)")
        else:
            log.append(f"ok   {name} ({time.perf_counter() - t0:.1f} s)")
    log.append(f"=== {len(failed)} of {len(gui_workflows.stages())} stages failed in "
               f"{time.perf_counter() - t_all:.0f} s")
    text = "\n".join(log)
    print(text)
    assert not failed, text
    # the exports the stages promise are where the ctx said
    out = ctx["out"]
    assert (out / "error_tools" / "chi2_map.png").stat().st_size > 1000
    assert (out / "batch" / "bundle" / "manifest.csv").is_file()
    assert (out / "batch" / "series.csv").is_file()
    assert (out / "relaxation" / "infinite_field_report.txt").is_file()
    assert (out / "tools" / "fit_table.csv").is_file()
