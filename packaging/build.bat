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
for /f %%v in ('"%PY%" -c "import larmor; print(larmor.__version__)"') do set "VER=%%v"
echo       version %VER%

echo [2/5] PyInstaller
"%PY%" -m PyInstaller packaging\larmor.spec --noconfirm --clean --log-level WARN || exit /b 1

echo [3/5] INSTALL.txt
copy /y packaging\INSTALL.txt dist\LARMOR\INSTALL.txt >nul || exit /b 1

echo [4/5] Inno Setup
set "ISCC=%LOCALAPPDATA%\Programs\Inno Setup 6\ISCC.exe"
if not exist "%ISCC%" set "ISCC=%ProgramFiles(x86)%\Inno Setup 6\ISCC.exe"
if not exist "%ISCC%" (
    echo Inno Setup not found; install it with
    echo     winget install --id JRSoftware.InnoSetup -e --scope user
    exit /b 1
)
"%ISCC%" /Q /DMyAppVersion=%VER% packaging\larmor.iss || exit /b 1

echo [5/5] zip
powershell -NoProfile -Command "Compress-Archive -Path 'dist\LARMOR' -DestinationPath 'dist\LARMOR-%VER%-win64.zip' -CompressionLevel Optimal -Force" || exit /b 1

echo.
echo done:
dir /b dist\LARMOR-%VER%-setup.exe dist\LARMOR-%VER%-win64.zip
echo smoke test: set QT_QPA_PLATFORM=offscreen ^&^& dist\LARMOR\LARMOR.exe  (crash log: %%USERPROFILE%%\LARMOR_crash.log must stay empty)
endlocal
