@echo off
rem Build the standalone Windows distribution of LARMOR in one go:
rem   1. reinstall the package into the pip build venv (never the conda env)
rem   2. PyInstaller  -> dist\LARMOR\LARMOR.exe (+ _internal\)
rem   3. INSTALL.txt next to the exe
rem   4. Inno Setup   -> dist\LARMOR-<version>-setup.exe
rem   5. zip          -> dist\LARMOR-<version>-win64.zip
rem Run from anywhere; it works in the repository root. See README.md here.
setlocal
cd /d "%~dp0.."
set "PY=packaging\.buildenv\Scripts\python.exe"
if not exist "%PY%" (
    echo Build venv missing. Create it first:
    echo     py -3.11 -m venv packaging\.buildenv
    echo     packaging\.buildenv\Scripts\python -m pip install ".[desktop]" pyinstaller
    exit /b 1
)

echo [1/5] reinstalling larmor into the build venv
"%PY%" -m pip install --no-deps --quiet ".[desktop]" || exit /b 1
rem The version through a file, not a for /f: cmd strips the outer quotes of
rem a for /f command that starts with a quoted path, so '"%PY%" -c "import
rem larmor; ..."' became "import' is not recognized" when this script ran
rem under PowerShell (cmd /c), VER stayed empty and the build produced
rem LARMOR--setup.exe with no version in the installer metadata (0.15.0).
set "VERFILE=%TEMP%\larmor_build_version.txt"
"%PY%" -c "import larmor, sys; sys.stdout.write(larmor.__version__)" > "%VERFILE%" || exit /b 1
set /p VER=<"%VERFILE%"
del "%VERFILE%" >nul 2>&1
if "%VER%"=="" (
    echo could not read larmor.__version__ from the build venv
    exit /b 1
)
echo       version %VER%
rem the third-party list names the versions that actually ship
"%PY%" packaging\third_party_licenses.py || exit /b 1

echo [2/5] PyInstaller
"%PY%" -m PyInstaller packaging\larmor.spec --noconfirm --clean --log-level WARN || exit /b 1

echo [3/5] INSTALL.txt, LICENSE, THIRD_PARTY_LICENSES.txt, LICENSES\
copy /y packaging\INSTALL.txt dist\LARMOR\INSTALL.txt >nul || exit /b 1
rem a public binary must carry its own licence and the LGPL texts of Qt /
rem PySide6 (LGPL-3.0 section 4): the wizard only DISPLAYS LICENSE, the zip
rem had nothing
copy /y LICENSE dist\LARMOR\LICENSE >nul || exit /b 1
copy /y packaging\THIRD_PARTY_LICENSES.txt dist\LARMOR\THIRD_PARTY_LICENSES.txt >nul || exit /b 1
xcopy /y /i /q packaging\LICENSES dist\LARMOR\LICENSES\ >nul || exit /b 1

echo [4/5] Inno Setup
set "ISCC=%LOCALAPPDATA%\Programs\Inno Setup 6\ISCC.exe"
if not exist "%ISCC%" set "ISCC=%ProgramFiles(x86)%\Inno Setup 6\ISCC.exe"
if not exist "%ISCC%" (
    echo Inno Setup not found; install it with
    echo     winget install --id JRSoftware.InnoSetup -e --scope user
    exit /b 1
)
"%ISCC%" /Q /DMyAppVersion=%VER% packaging\larmor.iss || exit /b 1

echo [5/6] zip
powershell -NoProfile -Command "Compress-Archive -Path 'dist\LARMOR' -DestinationPath 'dist\LARMOR-%VER%-win64.zip' -CompressionLevel Optimal -Force" || exit /b 1

echo [6/6] desktop shortcut
rem Sam's standing request (2026-09-30): every build replaces the desktop
rem shortcut so it points at dist\LARMOR\LARMOR.exe -- the runnable copy.
rem A shortcut to build\larmor\LARMOR.exe (PyInstaller's intermediate exe,
rem no _internal beside it) dies with "Failed to load Python DLL
rem ...\build\larmor\_internal\python311.dll".
powershell -NoProfile -ExecutionPolicy Bypass -File packaging\desktop_shortcut.ps1 -Version %VER% || echo       (shortcut not refreshed)

echo.
echo done:
dir /b dist\LARMOR-%VER%-setup.exe dist\LARMOR-%VER%-win64.zip
echo smoke test: set QT_QPA_PLATFORM=offscreen ^&^& dist\LARMOR\LARMOR.exe  (crash log: %%USERPROFILE%%\LARMOR_crash.log must stay empty)
endlocal
