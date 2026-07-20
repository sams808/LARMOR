# PyInstaller spec for LARMOR -- build with:
#     conda activate larmor
#     pip install pyinstaller
#     pyinstaller packaging/larmor.spec --noconfirm
#
# Produces dist/LARMOR/LARMOR.exe: a standalone desktop app with no conda and
# no Python required on the target machine.
#
# Notes learned the hard way with this stack:
#  * mrsimulator ships compiled extensions AND data files; collect both.
#  * lmfit pulls asteval/uncertainties dynamically -> hidden imports.
#  * matplotlib/PySide6 bring huge optional backends we do not use -> exclude
#    them to keep the folder near ~500 MB instead of well over 1 GB.
import sys
from pathlib import Path

from PyInstaller.utils.hooks import (
    collect_data_files, collect_dynamic_libs, collect_submodules,
)

block_cipher = None
ROOT = Path(SPECPATH).parent

datas = []
binaries = []
hiddenimports = []

# --- scientific stack that PyInstaller cannot see through -------------------
for pkg in ("mrsimulator", "csdmpy", "nmrglue"):
    datas += collect_data_files(pkg)
    binaries += collect_dynamic_libs(pkg)
    hiddenimports += collect_submodules(pkg)

hiddenimports += [
    "lmfit", "asteval", "uncertainties",
    "scipy.special._cdflib", "scipy._lib.messagestream",
    "pkg_resources.py2_warn",
]

# --- LARMOR's own resources -------------------------------------------------
datas += [(str(ROOT / "larmor" / "static"), "larmor/static")]
if (ROOT / "assets").exists():
    datas += [(str(ROOT / "assets"), "assets")]
for doc in ("README.md", "ROADMAP.md"):
    if (ROOT / doc).exists():
        datas.append((str(ROOT / doc), "."))
if (ROOT / "docs" / "tutorials").exists():
    datas.append((str(ROOT / "docs" / "tutorials"), "docs/tutorials"))

a = Analysis(
    [str(ROOT / "packaging" / "launcher.py")],
    pathex=[str(ROOT)],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    runtime_hooks=[],
    excludes=[
        # never used by the desktop app; each is tens of MB
        "tkinter", "PyQt5", "PyQt6", "PySide2",
        "matplotlib.backends._backend_tk", "matplotlib.backends.backend_qt5agg",
        "IPython", "jupyter", "notebook", "jupyterlab", "nbconvert",
        "pytest", "sphinx", "setuptools._distutils",
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="LARMOR",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,          # windowed app; errors go to the crash log
    disable_windowed_traceback=False,
    icon=str(ROOT / "assets" / "larmor.ico")
        if (ROOT / "assets" / "larmor.ico").exists() else None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name="LARMOR",
)
