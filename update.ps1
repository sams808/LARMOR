# =====================================================================
#  LARMOR updater (Windows).
#  Double-click update.bat, or right-click this file -> Run with PowerShell.
#
#  Pulls the latest code (if this is a git checkout), then re-runs the
#  installer, which updates the environment in place and refreshes the
#  Desktop shortcut.
# =====================================================================
$ErrorActionPreference = "Stop"
$repo = $PSScriptRoot
Write-Host ""
Write-Host "==== LARMOR updater ====" -ForegroundColor Cyan

if (Test-Path (Join-Path $repo ".git")) {
    $git = Get-Command git -ErrorAction SilentlyContinue
    if ($git) {
        Write-Host "Pulling the latest code..." -ForegroundColor Cyan
        & git -C $repo pull
    } else {
        Write-Host "git is not installed - skipping the code update." -ForegroundColor Yellow
        Write-Host "(Download the latest ZIP from GitHub if you want the newest code.)"
    }
} else {
    Write-Host "Not a git checkout - updating the environment only." -ForegroundColor Yellow
    Write-Host "(Re-download the ZIP from GitHub for newer code, then run this again.)"
}

Write-Host "Updating the environment..." -ForegroundColor Cyan
& (Join-Path $repo "install.ps1")
