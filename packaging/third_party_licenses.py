"""Regenerate packaging/THIRD_PARTY_LICENSES.txt for the build that is about
to be made: the hand-written header (scope, the Qt / PySide6 LGPL terms, the
bundled data files with their own terms) with the version filled in, then
one block per installed package of the build environment from pip-licenses.

    packaging\\.buildenv\\Scripts\\python.exe packaging\\third_party_licenses.py

Run by packaging/build.bat after larmor is reinstalled into the build
environment, so the list names the versions that actually ship. Needs
``pip-licenses`` in the build environment (installed on first use).
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "THIRD_PARTY_LICENSES.txt"
SEPARATOR = "=" * 68


def header(version: str, python: str) -> str:
    title = f"THIRD-PARTY LICENCES -- LARMOR {version}, Windows build"
    return f"""{title}
{"=" * len(title)}

Generated from the frozen-build environment (packaging/.buildenv, Python
{python}) by packaging/third_party_licenses.py, with

    pip-licenses --format=plain-vertical --with-license-file --no-license-path
                 --with-urls --order=name

One block per package: name, version, licence, URL, then the licence text that
the installed wheel carries. Blocks are in alphabetical order after this header.

Scope
-----
The list covers every package installed in the build environment, which is a
superset of what the distributed LARMOR folder contains. Build-time only, NOT
present in the distributed folder: pyinstaller, pyinstaller-hooks-contrib,
altgraph, pefile, pywin32-ctypes, pip-licenses, prettytable, wcwidth, pytest,
pluggy, iniconfig. "larmor" itself is MIT (see LICENSE next to this file).

PyInstaller is GPL-2.0-or-later WITH Bootloader-exception: the bootloader
linked into LARMOR.exe may be distributed under LARMOR's own licence (MIT).

Qt, PySide6 and shiboken6 (The Qt Company)
------------------------------------------
PySide6, PySide6_Essentials, PySide6_Addons and shiboken6 are used under the
GNU Lesser General Public License version 3 (LGPL-3.0-only), one of the
licence options their metadata offers. The only licence file their wheels carry
is the Qt commercial-licence notice reproduced in their blocks below; that
notice does not apply to this build. The Qt libraries ship unmodified, as
separate DLLs under _internal\\PySide6, and can be replaced by the user with
another build of the same Qt major version, as the LGPL requires. The Qt
modules LARMOR does not use -- among them Qt Virtual Keyboard, which is
GPL-3.0-only -- are left out of the build. Sources and terms:

    Qt source           https://code.qt.io
    PySide6 source      https://code.qt.io/cgit/pyside/pyside-setup.git
    PySide6 licensing   https://doc.qt.io/qtforpython-6/licenses.html
    LGPL-3.0 text       LICENSES\\LGPL-3.0.txt next to this file
    GPL-3.0 text        LICENSES\\GPL-3.0.txt next to this file
                        (the LGPL-3.0 is a set of permissions on top of it)

Other bundled data files with their own terms
---------------------------------------------
matplotlib fonts: DejaVu (Bitstream Vera licence) and STIX (SIL OFL), licence
files shipped under _internal\\matplotlib\\mpl-data\\fonts\\ttf.
pyqtgraph CET colour maps: CC BY 4.0, licence file shipped under
_internal\\pyqtgraph\\colors\\maps.
Bird photographs of Help > More... (offline pack): photographer and Creative
Commons licence of each in _internal\\larmor\\xfact\\packs\\birds\\assets\\offline\\CREDITS.txt.

{SEPARATOR}

"""


def main() -> int:
    py = sys.executable
    try:
        import piplicenses  # noqa: F401
    except ImportError:
        subprocess.run([py, "-m", "pip", "install", "--quiet", "pip-licenses"], check=True)
    import larmor
    blocks = subprocess.run(
        [py, "-m", "piplicenses", "--format=plain-vertical", "--with-license-file",
         "--no-license-path", "--with-urls", "--order=name"],
        check=True, capture_output=True, text=True, encoding="utf-8").stdout
    if f"\nlarmor\n{larmor.__version__}\n" not in "\n" + blocks:
        raise SystemExit(f"pip-licenses did not list larmor {larmor.__version__}: "
                         "reinstall it into the build environment first")
    OUT.write_text(header(larmor.__version__, sys.version.split()[0]) + blocks.lstrip(),
                   encoding="utf-8")
    print(f"wrote {OUT.name}: LARMOR {larmor.__version__}, {blocks.count(chr(10))} lines of blocks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
