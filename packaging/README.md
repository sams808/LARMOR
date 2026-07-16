# Packaging LARMOR as a standalone Windows app

Produces `dist/LARMOR/LARMOR.exe` — no conda, no Python needed on the target
machine. Students double-click the exe; that is the whole install.

## Build

```
conda activate larmor
pip install pyinstaller
pyinstaller packaging/larmor.spec --noconfirm
```

The result is `dist/LARMOR/` (~500 MB, dominated by scipy/mrsimulator/PySide6).
Zip that folder to distribute, or wrap it with an installer (Inno Setup /
NSIS) for a Start-menu entry.

## First run

- Kernel caches (the one-time Czjzek simulation per field/spin rate) are held
  in memory per session; nothing is written outside the user's own folders.
- If the app fails to start, a `LARMOR_crash.log` appears in the user's home
  directory with the full traceback (the launcher installs this handler).

## Notes

- `larmor.spec` collects mrsimulator/csdmpy/nmrglue data files and compiled
  libraries, and lists lmfit's dynamic dependencies (asteval, uncertainties)
  as hidden imports — PyInstaller cannot infer these.
- Jupyter, tkinter, the other Qt bindings and test tooling are excluded to
  keep the folder from ballooning past 1 GB.
- To add an icon, drop `packaging/larmor.ico`; the spec picks it up
  automatically.

## CI (future, roadmap v1.0)

A GitHub Actions job on `windows-latest` can run this spec and attach
`LARMOR.zip` to each tagged release, so `larmor.bat` / conda become optional
and the download is the reference distribution.
