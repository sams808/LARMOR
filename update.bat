@echo off
REM Double-click this to update LARMOR. It runs update.ps1 with the
REM PowerShell execution policy bypassed.
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0update.ps1"
