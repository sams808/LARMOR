"""The engine stages of the distribution check (larmor.distcheck.engines):
every Qt-free capability on synthetic / bundled data, run the way
``LARMOR.exe --distcheck --quick`` runs them, in a development install."""
import json
import re
from pathlib import Path


from larmor.distcheck import engines


#: the capability groups the engine set must keep covering
GROUPS = ("engine-readers", "engine-processing", "engine-models", "engine-physics",
          "engine-fitting", "engine-error-tools", "engine-constraints", "engine-multi",
          "engine-records", "engine-specialty", "engine-figures")

#: the groups that read the bundled examples when present and must fall back
#: to synthetic data, saying so, when they are not
FALLBACK = ("engine-records", "engine-specialty", "engine-figures")


def _names() -> list:
    return [n for n, _ in engines.stages()]


def test_stage_names_are_engine_prefixed_and_unique():
    names = _names()
    assert set(GROUPS) <= set(names), sorted(set(GROUPS) - set(names))
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
    # the figure stage exported every renderer in the three formats the
    # export dialogs offer, as complete files
    from larmor import figures
    figs = tmp_path / "out" / "figures"
    for kind in figures.RENDERERS:
        for fmt, head in (("png", b"\x89PNG"), ("svg", b"<?xml"), ("pdf", b"%PDF")):
            data = (figs / f"fig_{kind}.{fmt}").read_bytes()
            assert data.startswith(head) and len(data) > 1000, (kind, fmt, data[:16])


def test_engine_stages_fall_back_to_synthetic_data_without_the_examples(tmp_path, monkeypatch):
    """An installation without examples/ (a stripped or development copy):
    the stages that read the bundled data run on synthetic data instead,
    say so in the log, and still pass."""
    from larmor import distcheck

    monkeypatch.setattr(distcheck, "_examples_dir", lambda: None)
    log = tmp_path / "distcheck.log"
    monkeypatch.setenv("LARMOR_DISTCHECK_LOG", str(log))
    rc = distcheck.run(["--no-gui", "--quick", "--out", str(tmp_path / "out"),
                        "--only", ",".join(FALLBACK)])
    text = log.read_text(encoding="utf-8")
    assert rc == 0, text[-6000:]
    report = json.loads(log.with_suffix(".json").read_text(encoding="utf-8"))
    assert [s["name"] for s in report["stages"]] == list(FALLBACK)
    for s in report["stages"]:
        assert s["ok"], (s["name"], s["error"][-3000:])
        # "no bundled example data: ..." / "(no bundled examples: ...)"
        assert any("no bundled example" in n for n in s["notes"]), (s["name"], s["notes"])
