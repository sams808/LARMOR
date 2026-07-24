# =====================================================================
#  LARMOR one-click installer (Windows).
#  Double-click install.bat, or right-click this file -> Run with PowerShell.
#
#  It creates the 'larmor' Conda environment (or a pip virtual env if you
#  don't have Conda), installs LARMOR, and puts a LARMOR shortcut on your
#  Desktop that launches the app. Safe to re-run: it updates in place.
# =====================================================================
$ErrorActionPreference = "Stop"
$repo = $PSScriptRoot
Write-Host ""
Write-Host "==== LARMOR installer ====" -ForegroundColor Cyan
Write-Host "Repository: $repo"
Write-Host ""

function Find-Conda {
    $cmd = Get-Command conda -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }
    $roots = @(
        "$env:USERPROFILE\miniconda3", "$env:USERPROFILE\anaconda3",
        "$env:USERPROFILE\miniforge3", "$env:USERPROFILE\mambaforge",
        "$env:USERPROFILE\xraylarch",
        "$env:LOCALAPPDATA\miniconda3", "$env:LOCALAPPDATA\anaconda3",
        "C:\ProgramData\miniconda3", "C:\ProgramData\anaconda3"
    )
    foreach ($r in $roots) {
        foreach ($sub in @("Scripts\conda.exe", "condabin\conda.bat")) {
            $p = Join-Path $r $sub
            if (Test-Path $p) { return $p }
        }
    }
    return $null
}

$python = $null
$conda = Find-Conda

if ($conda) {
    Write-Host "Found Conda: $conda" -ForegroundColor Green
    $envList = & $conda env list 2>$null
    $hasLarmor = $false
    foreach ($line in $envList) {
        if ($line -match '^\s*larmor\s' -or $line -match '[\\/]envs[\\/]larmor\s*$') {
            $hasLarmor = $true
        }
    }
    if ($hasLarmor) {
        Write-Host "Updating the existing 'larmor' environment..." -ForegroundColor Cyan
        & $conda env update -n larmor -f "$repo\environment.yml"
    } else {
        Write-Host "Creating the 'larmor' environment (a few minutes)..." -ForegroundColor Cyan
        & $conda env create -f "$repo\environment.yml"
    }
    # locate the environment's python.exe. Parsing the prefix out of
    # `conda env list` is reliable across conda layouts (Scripts / condabin /
    # Library\bin); fall back to `conda run` if the listing looks unusual.
    $python = $null
    foreach ($line in (& $conda env list 2>$null)) {
        if ($line -match '^\s*larmor\s+\*?\s*([A-Za-z]:[\\/].+?)\s*$') {
            $python = Join-Path $matches[1] "python.exe"
        } elseif ($line -match '([A-Za-z]:[\\/].*[\\/]envs[\\/]larmor)\s*$') {
            $python = Join-Path $matches[1] "python.exe"
        }
    }
    if (-not $python -or -not (Test-Path $python)) {
        $out = & $conda run -n larmor python -c "import sys;print(sys.executable)" 2>$null
        $cand = ($out | Where-Object { $_ -match 'python\.exe\s*$' } | Select-Object -Last 1)
        if ($cand) { $python = $cand.Trim() }
    }
} else {
    Write-Host "Conda not found - falling back to Python + pip." -ForegroundColor Yellow
    $py = Get-Command python -ErrorAction SilentlyContinue
    if (-not $py) { $py = Get-Command python3 -ErrorAction SilentlyContinue }
    if (-not $py) {
        Write-Host "No Python found either." -ForegroundColor Red
        Write-Host "Please install Miniconda (recommended) or Python 3.11, then re-run."
        Write-Host "  https://docs.conda.io/en/latest/miniconda.html"
        Write-Host "See INSTALL.md for details."
        Read-Host "Press Enter to close"; exit 1
    }
    $ver = (& $py.Source -c "import sys;print('.'.join(map(str, sys.version_info[:2])))").Trim()
    Write-Host "Found Python $ver at $($py.Source)"
    if ($ver -eq "3.13" -or $ver -eq "3.14" -or $ver -eq "3.15") {
        Write-Host ""
        Write-Host "Python $ver is too new for one of LARMOR's dependencies (mrsimulator)." -ForegroundColor Red
        Write-Host "Install Miniconda (recommended - handles this for you), or Python 3.11,"
        Write-Host "then re-run this installer. See INSTALL.md."
        Read-Host "Press Enter to close"; exit 1
    }
    $venv = Join-Path $repo ".venv"
    if (-not (Test-Path (Join-Path $venv "Scripts\python.exe"))) {
        Write-Host "Creating a virtual environment in .venv ..." -ForegroundColor Cyan
        & $py.Source -m venv $venv
    }
    $python = Join-Path $venv "Scripts\python.exe"
    Write-Host "Installing LARMOR and its dependencies (a few minutes)..." -ForegroundColor Cyan
    & $python -m pip install --upgrade pip
    & $python -m pip install -r "$repo\requirements.txt"
}

if (-not $python -or -not (Test-Path $python)) {
    Write-Host "Installation did not produce a working Python. See INSTALL.md." -ForegroundColor Red
    Read-Host "Press Enter to close"; exit 1
}
Write-Host ""
Write-Host "LARMOR Python: $python" -ForegroundColor Green

# sanity check: LARMOR + the GUI toolkit import
& $python -c "import larmor, PySide6, pyqtgraph" 2>$null
if ($LASTEXITCODE -ne 0) {
    Write-Host "Warning: LARMOR or the GUI toolkit did not import cleanly." -ForegroundColor Yellow
    Write-Host "The Desktop shortcut was still created; see INSTALL.md if it fails to launch."
} else {
    Write-Host "Import check passed." -ForegroundColor Green
}

# create the Desktop launcher (paths baked in, so no activation needed)
$desktop = [Environment]::GetFolderPath("Desktop")
$bat = Join-Path $desktop "LARMOR.bat"
$batLines = @(
    '@echo off',
    'title LARMOR',
    "cd /d `"$repo`"",
    "`"$python`" -m larmor.cli desktop",
    'if errorlevel 1 pause'
)
Set-Content -Path $bat -Value $batLines -Encoding ASCII
Write-Host ""
Write-Host "==== Done ====" -ForegroundColor Green
Write-Host "A 'LARMOR.bat' is on your Desktop - double-click it to launch the app."
Write-Host "(Or run 'larmor desktop' after activating the environment.)"
Read-Host "Press Enter to close"
