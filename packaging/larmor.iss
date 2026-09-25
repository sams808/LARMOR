; Inno Setup script for LARMOR -- wraps the PyInstaller folder dist\LARMOR
; into dist\LARMOR-<version>-setup.exe (Start-menu entry, optional desktop
; icon, uninstaller). Compile after PyInstaller, from the repository root:
;
;     "%LOCALAPPDATA%\Programs\Inno Setup 6\ISCC.exe" /DMyAppVersion=0.14.0 packaging\larmor.iss
;
; or run packaging\build.bat, which does the whole chain. The version is
; passed on the command line so this file never carries a stale number.
;
; Per-user by default (no administrator rights needed on a student's laptop:
; the app lands in %LOCALAPPDATA%\Programs\LARMOR); the wizard offers the
; all-users install when the account can elevate.

#ifndef MyAppVersion
  #define MyAppVersion "0.0.0"
#endif
#define MyAppName "LARMOR"
#define MyAppPublisher "McCloy group, Washington State University"
#define MyAppURL "https://github.com/sams808/LARMOR"
#define MyAppExeName "LARMOR.exe"
#ifndef SourceDir
  ; the PyInstaller folder to wrap; /DSourceDir=..\dist_test\LARMOR builds
  ; from a side folder when dist\LARMOR is locked by a running instance
  #define SourceDir "..\dist\LARMOR"
#endif

[Setup]
AppId={{A5B3F6C2-7D41-4E0B-9F3A-2C6E8D1B4F70}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppVerName={#MyAppName} {#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppPublisherURL={#MyAppURL}
AppSupportURL={#MyAppURL}
AppUpdatesURL={#MyAppURL}
DefaultDirName={autopf}\{#MyAppName}
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
OutputDir=..\dist
OutputBaseFilename={#MyAppName}-{#MyAppVersion}-setup
UninstallDisplayIcon={app}\{#MyAppExeName}
UninstallDisplayName={#MyAppName} {#MyAppVersion}
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
CloseApplications=yes
#ifexist "..\assets\larmor.ico"
SetupIconFile=..\assets\larmor.ico
#endif
#ifexist "..\LICENSE"
LicenseFile=..\LICENSE
#endif
#ifexist "INSTALL.txt"
InfoBeforeFile=INSTALL.txt
#endif

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"

[InstallDelete]
; wipe the previous version's _internal before laying down the new one:
; an overwrite install leaves files the new version no longer ships, and a
; stale compiled module there is exactly the kind of fault that is
; undebuggable from another machine. Everything under _internal belongs to
; LARMOR, and the user's recipes, projects, settings and logs live
; elsewhere, so nothing of theirs is at risk.
Type: filesandordirs; Name: "{app}\_internal"

[Files]
; the whole PyInstaller folder: LARMOR.exe, _internal\, INSTALL.txt
Source: "{#SourceDir}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"
Name: "{group}\{#MyAppName} install notes"; Filename: "{app}\INSTALL.txt"
Name: "{group}\{cm:UninstallProgram,{#MyAppName}}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "{cm:LaunchProgram,{#StringChange(MyAppName, '&', '&&')}}"; Flags: nowait postinstall skipifsilent

[UninstallDelete]
; PyInstaller leaves nothing else behind; recipes, projects and figures live
; where the user saved them and are never touched by the uninstaller
Type: filesandordirs; Name: "{app}\_internal"
