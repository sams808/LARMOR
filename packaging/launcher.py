"""Frozen-app entry point.

Kept tiny and import-light so PyInstaller's dependency graph starts clean.
Installs a last-resort crash handler that writes a diagnostics file next to
the user, because a windowed frozen app otherwise dies silently.
"""
import sys
import traceback
from pathlib import Path


def _crash_log(exc_type, exc, tb):
    import datetime

    try:
        target = Path.home() / "LARMOR_crash.log"
        with open(target, "a", encoding="utf-8") as f:
            f.write(f"\n=== {datetime.datetime.now().isoformat()} ===\n")
            traceback.print_exception(exc_type, exc, tb, file=f)
        # also surface a dialog if Qt is already up
        try:
            from PySide6.QtWidgets import QApplication, QMessageBox

            if QApplication.instance():
                QMessageBox.critical(
                    None, "LARMOR error",
                    f"{exc}\n\nDetails written to {target}")
        except Exception:
            pass
    except Exception:
        traceback.print_exception(exc_type, exc, tb)


def main() -> int:
    sys.excepthook = _crash_log
    from larmor.desktop.app import main as desktop_main

    return desktop_main()


if __name__ == "__main__":
    raise SystemExit(main())
