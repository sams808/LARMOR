@echo off
REM Double-click this to install LARMOR. It runs install.ps1 with the
REM PowerShell execution policy bypassed (so no policy tweaking is needed).
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0install.ps1"
