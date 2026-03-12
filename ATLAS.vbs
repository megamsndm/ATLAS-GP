Set WshShell = CreateObject("WScript.Shell")
Set fso = CreateObject("Scripting.FileSystemObject")

' Project root = folder containing this .vbs
Dim root
root = fso.GetParentFolderName(WScript.ScriptFullName)

Dim bin
bin = root & "\build\atlas.exe"

If Not fso.FileExists(bin) Then
    MsgBox "atlas.exe not found." & Chr(13) & Chr(13) & _
           "Please run ATLAS.bat first to build the project." & Chr(13) & _
           "(" & bin & ")", vbExclamation, "ATLAS-GP"
    WScript.Quit 1
End If

Dim gui
gui = root & "\gui\atlas_gui.py"

If Not fso.FileExists(gui) Then
    MsgBox "atlas_gui.py not found." & Chr(13) & _
           "(" & gui & ")", vbExclamation, "ATLAS-GP"
    WScript.Quit 1
End If

' Find Python
Dim python
python = ""
Dim candidates(3)
candidates(0) = "python"
candidates(1) = "python3"
candidates(2) = "py"
candidates(3) = "python.exe"

Dim i
For i = 0 To 3
    On Error Resume Next
    Dim ret
    ret = WshShell.Run("cmd /c " & candidates(i) & " --version", 0, True)
    If Err.Number = 0 And ret = 0 Then
        python = candidates(i)
        Exit For
    End If
    On Error GoTo 0
Next

If python = "" Then
    MsgBox "Python not found." & Chr(13) & Chr(13) & _
           "Install Python 3.10+ from https://python.org", vbExclamation, "ATLAS-GP"
    WScript.Quit 1
End If

' Launch GUI silently
Dim cmd
cmd = python & " """ & gui & """ --binary """ & bin & """"
WshShell.Run cmd, 0, False
WScript.Quit 0
