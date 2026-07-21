; Inno Setup script for XFit. Install Inno Setup, then compile this file.
#define MyAppName "XFit WAXD Batch Analyzer"
#define MyAppVersion "2.0.0"
#define MyAppPublisher "XFit Team"
#define MyAppExeName "XFit.exe"

[Setup]
AppId={{B7E82310-CE12-4B95-AE8D-XFITWAXD2000}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
DefaultDirName={autopf}\XFit
DefaultGroupName=XFit
OutputDir=dist\installer
OutputBaseFilename=XFitSetup-2.0.0
Compression=lzma
SolidCompression=yes
WizardStyle=modern
ArchitecturesAllowed=x64
ArchitecturesInstallIn64BitMode=x64

[Files]
Source: "..\dist\XFit\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\XFit"; Filename: "{app}\{#MyAppExeName}"
Name: "{autodesktop}\XFit"; Filename: "{app}\{#MyAppExeName}"; Tasks: desktopicon

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; GroupDescription: "Additional icons:"

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "Launch XFit"; Flags: nowait postinstall skipifsilent
