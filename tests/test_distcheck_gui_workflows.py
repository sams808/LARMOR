"""The dialog-tool stages of the distribution check (larmor.distcheck.
gui_workflows) pass in a development install, offscreen, in quick mode,
through ``larmor.distcheck.run`` -- the path ``LARMOR.exe --distcheck --gui``
takes, with its sandbox for the user's stores. The environment is set
BEFORE PySide6 is imported, as ``run`` does: no window ever shows, native
dialogs become widgets the sweep can dismiss, nothing of the user's session
is read, no kernel pre-build thread starts.
"""
import os
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


def test_dialog_stages_pass_in_a_development_install(tmp_path, monkeypatch):
    import json

    from larmor import distcheck
    from larmor.distcheck import gui_workflows

    names = [n for n, _fn in gui_workflows.stages()]
    log = tmp_path / "distcheck.log"
    monkeypatch.setenv("LARMOR_DISTCHECK_LOG", str(log))
    out = tmp_path / "out"
    rc = distcheck.run(["--gui", "--quick", "--out", str(out), "--only", ",".join(names)])
    text = log.read_text(encoding="utf-8")
    print(text[-8000:])
    assert rc == 0, text[-8000:]
    report = json.loads(log.with_suffix(".json").read_text(encoding="utf-8"))
    assert [s["name"] for s in report["stages"]] == names
    assert all(s["ok"] for s in report["stages"]), text[-8000:]
    # the exports the stages promise are where the ctx said (the GUI groups
    # run in a child process of their own, under <out>/dialog_stages)
    if (out / "dialog_stages").is_dir():
        out = out / "dialog_stages"
    assert (out / "error_tools" / "chi2_map.png").stat().st_size > 1000
    assert (out / "batch" / "bundle" / "manifest.csv").is_file()
    assert (out / "batch" / "series.csv").is_file()
    assert (out / "relaxation" / "infinite_field_report.txt").is_file()
    assert (out / "tools" / "fit_table.csv").is_file()
