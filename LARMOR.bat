@echo off
REM ============================================================
REM  LARMOR launcher -- double-click to start the app.
REM  One-time setup (in an Anaconda/Miniconda prompt, from this
REM  folder):   conda env create -f environment.yml
REM ============================================================
setlocal
set "PORT=8642"
cd /d "%~dp0"

REM 1) explicit override: set LARMOR_PYTHON to a python.exe if you know it
if defined LARMOR_PYTHON if exist "%LARMOR_PYTHON%" goto :run

REM 2) look in the usual conda env locations
for %%P in (
  "%USERPROFILE%\miniconda3\envs\larmor\python.exe"
  "%USERPROFILE%\anaconda3\envs\larmor\python.exe"
  "%USERPROFILE%\mambaforge\envs\larmor\python.exe"
  "%USERPROFILE%\miniforge3\envs\larmor\python.exe"
  "%USERPROFILE%\xraylarch\envs\larmor\python.exe"
  "%LOCALAPPDATA%\miniconda3\envs\larmor\python.exe"
  "%LOCALAPPDATA%\anaconda3\envs\larmor\python.exe"
  "C:\ProgramData\miniconda3\envs\larmor\python.exe"
  "C:\ProgramData\anaconda3\envs\larmor\python.exe"
) do if exist "%%~P" (
  set "LARMOR_PYTHON=%%~P"
  goto :run
)

REM 3) fall back to conda on PATH
where conda >nul 2>nul
if %errorlevel%==0 (
  echo Starting LARMOR via conda...
  conda run -n larmor python -m larmor.cli app --port %PORT% --open
  goto :end
)

echo.
echo Could not find the 'larmor' conda environment.
echo.
echo One-time setup: open an Anaconda/Miniconda prompt in this folder and run
echo     conda env create -f environment.yml
echo then double-click this file again.
echo.
pause
exit /b 1

:run
echo Starting LARMOR (%LARMOR_PYTHON%)...
echo The app opens in your browser; keep this window open while working.
"%LARMOR_PYTHON%" -m larmor.cli app --port %PORT% --open

:end
echo.
echo LARMOR stopped.
pause
