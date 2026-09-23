"""larmor.selftest: the fit-engine self-test behind `LARMOR.exe --selftest`.

The core stage must pass in the development environment (it is what a
student's machine is compared against), write an appendable log where
LARMOR_SELFTEST_LOG points, and report a failing stage without raising.
"""
import os

from larmor import selftest


class _Say:
    def __init__(self):
        self.lines, self.fails = [], []

    def __call__(self, *parts):
        self.lines.append(" ".join(str(p) for p in parts))

    def fail(self, stage, exc=None):
        self.fails.append(stage)


def test_core_fits_pass_here_and_report_both_models():
    say = _Say()
    assert selftest.core_fits(say) is True
    assert say.fails == []
    joined = "\n".join(say.lines)
    assert "core gauss_lor:" in joined and "core czjzek:" in joined
    assert "moved: True" in joined and "evaluations" in joined


def test_run_writes_an_appendable_log_and_exit_code(tmp_path, monkeypatch):
    log = tmp_path / "selftest.log"
    monkeypatch.setenv("LARMOR_SELFTEST_LOG", str(log))
    assert selftest.log_path() == log
    # a broken stage is reported, never raised
    monkeypatch.setattr(selftest, "core_fits", lambda say: (_ for _ in ()).throw(RuntimeError("boom")))
    assert selftest.run(gui=False) == 1
    text = log.read_text(encoding="utf-8")
    assert "=== LARMOR self-test" in text and "FAIL self-test" in text and "boom" in text
    assert "larmor 0.14" in text or "larmor " in text
    monkeypatch.setattr(selftest, "core_fits", lambda say: True)
    assert selftest.run(gui=False) == 0
    text2 = log.read_text(encoding="utf-8")
    assert text2.startswith(text) and text2.count("=== LARMOR self-test") == 2   # appended
    assert "=== PASS" in text2
    assert os.environ.get("LARMOR_SELFTEST_LOG") == str(log)
