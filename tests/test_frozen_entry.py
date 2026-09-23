"""The frozen (PyInstaller, windowed) entry point has no stderr. Everything
main() does before the window appears must survive ``sys.stderr is None``:
the 0.13.0 installer smoke test found faulthandler.enable() raising there,
which killed every frozen build at start-up."""
import faulthandler
import sys

import pytest

pytest.importorskip("PySide6")


def test_faulthandler_install_survives_a_windowed_exe(monkeypatch, tmp_path):
    from larmor.desktop import app as appmod

    was_enabled = faulthandler.is_enabled()
    try:
        # console run: stderr present -> dumps there
        assert appmod._install_faulthandler() == "stderr"
        assert faulthandler.is_enabled()

        # frozen windowed exe: no stderr -> dumps to the crash log, no exception
        faulthandler.disable()
        monkeypatch.setattr(sys, "stderr", None)
        log = tmp_path / "larmor_crash.log"
        monkeypatch.setattr(appmod, "_crash_log_path", lambda: str(log))
        assert appmod._install_faulthandler() == str(log)
        assert faulthandler.is_enabled()
        assert log.exists()
    finally:
        faulthandler.disable()
        if was_enabled:
            faulthandler.enable()
