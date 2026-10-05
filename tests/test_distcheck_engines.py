"""The engine stages of the distribution check (larmor.distcheck.engines):
every Qt-free capability on synthetic / bundled data, run the way
``LARMOR.exe --distcheck --quick`` runs them, in a development install."""
import json
import re
from pathlib import Path


from larmor.distcheck import engines


def _names() -> list:
    return [n for n, _ in engines.stages()]


def test_stage_names_are_engine_prefixed_and_unique():
    names = _names()
    assert len(names) >= 6, names
    assert len(names) == len(set(names)), names
    for n in names:
        assert n.startswith("engine-"), n
        assert re.fullmatch(r"engine(-[a-z0-9]+)+", n), n
    # every stage is a plain callable with the (say, ctx) contract
    for _n, fn in engines.stages():
        assert callable(fn)


def test_engine_stages_pass_in_a_development_install(tmp_path, monkeypatch):
    from larmor import distcheck

    names = _names()
    log = tmp_path / "distcheck.log"
    monkeypatch.setenv("LARMOR_DISTCHECK_LOG", str(log))
    rc = distcheck.run(["--no-gui", "--quick", "--out", str(tmp_path / "out"),
                        "--only", ",".join(names)])
    text = log.read_text(encoding="utf-8")
    assert rc == 0, text[-6000:]
    report = json.loads(log.with_suffix(".json").read_text(encoding="utf-8"))
    assert [s["name"] for s in report["stages"]] == names
    failed = [s["name"] for s in report["stages"] if not s["ok"]]
    assert not failed, (failed, text[-6000:])
    assert "=== PASS" in text
    # nothing was written outside the scratch folder the check was given
    assert (tmp_path / "out").is_dir()
    assert not list(Path.home().glob("larmor_distcheck_*"))
