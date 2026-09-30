# Replace the desktop shortcut "LARMOR.lnk" so it launches the runnable
# build, dist\LARMOR\LARMOR.exe. Called by packaging\build.bat as its last
# step; safe to run by hand from anywhere:
#     powershell -NoProfile -ExecutionPolicy Bypass -File packaging\desktop_shortcut.ps1
# Only dist\LARMOR\LARMOR.exe is runnable: PyInstaller also leaves an
# intermediate build\larmor\LARMOR.exe with no _internal folder beside it,
# and a shortcut to that one fails with "Failed to load Python DLL".
param([string]$Version = "")

$repo = Split-Path -Parent $PSScriptRoot
$exe = Join-Path $repo "dist\LARMOR\LARMOR.exe"
if (-not (Test-Path $exe)) {
    Write-Host "      no dist\LARMOR\LARMOR.exe to point at"
    exit 1
}
$desk = [Environment]::GetFolderPath("Desktop")
$lnk = Join-Path $desk "LARMOR.lnk"
$sh = New-Object -ComObject WScript.Shell
$s = $sh.CreateShortcut($lnk)
$s.TargetPath = $exe
$s.WorkingDirectory = Split-Path $exe
$s.IconLocation = "$exe,0"
$s.Description = if ($Version) { "LARMOR $Version (dist build)" } else { "LARMOR (dist build)" }
$s.Save()
Write-Host "      $lnk -> $exe"
exit 0
