@echo off
setlocal EnableDelayedExpansion
title ATLAS-GP  ^|  MACE-MP-0 Installer

echo.
echo =====================================================================
echo   ATLAS-GP  ^|  MACE-MP-0 Neural Network Potential Installer
echo =====================================================================
echo.
echo This will install:
echo   - PyTorch (CPU version, ~800MB)
echo   - mace-torch (~50MB)
echo   - ase (Atomic Simulation Environment, ~10MB)
echo   - numpy, scipy (if not present)
echo.
echo Total download: ~900MB
echo Time: 5-15 minutes depending on internet speed
echo.

:: ── Find Python ──────────────────────────────────────────────────────────────
set "PYTHON="

:: Check python launcher first (most reliable on Windows)
py --version >nul 2>&1
if not errorlevel 1 (
    set "PYTHON=py"
    goto :found_python
)

:: Check python3
python3 --version >nul 2>&1
if not errorlevel 1 (
    set "PYTHON=python3"
    goto :found_python
)

:: Check python
python --version >nul 2>&1
if not errorlevel 1 (
    set "PYTHON=python"
    goto :found_python
)

:: Check common install locations
for %%p in (
    "%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python311\python.exe"
    "%LOCALAPPDATA%\Programs\Python\Python310\python.exe"
    "C:\Python312\python.exe"
    "C:\Python311\python.exe"
    "C:\Python310\python.exe"
    "%USERPROFILE%\AppData\Local\Programs\Python\Python312\python.exe"
    "%USERPROFILE%\AppData\Local\Programs\Python\Python311\python.exe"
) do (
    if exist %%p (
        set "PYTHON=%%~p"
        goto :found_python
    )
)

echo [ERROR] Python not found.
echo.
echo Please install Python 3.10 or newer from:
echo   https://www.python.org/downloads/
echo.
echo IMPORTANT: During install, check "Add Python to PATH"
echo.
pause
exit /b 1

:found_python
echo [OK]   Python found: %PYTHON%

:: ── Check Python version ─────────────────────────────────────────────────────
for /f "tokens=2" %%v in ('%PYTHON% --version 2^>^&1') do set PY_VER=%%v
echo [OK]   Python version: %PY_VER%

:: Extract major.minor
for /f "tokens=1,2 delims=." %%a in ("%PY_VER%") do (
    set PY_MAJOR=%%a
    set PY_MINOR=%%b
)

if %PY_MAJOR% LSS 3 (
    echo [ERROR] Python 3.10+ required. Found Python %PY_VER%
    pause & exit /b 1
)
if %PY_MAJOR% EQU 3 if %PY_MINOR% LSS 10 (
    echo [ERROR] Python 3.10+ required. Found Python %PY_VER%
    echo         Please install Python 3.10, 3.11, or 3.12
    pause & exit /b 1
)

:: ── Check if already installed ───────────────────────────────────────────────
echo.
echo [..] Checking existing installation...
%PYTHON% -c "import mace; print('MACE version:', mace.__version__)" >nul 2>&1
if not errorlevel 1 (
    echo [OK]   MACE is already installed!
    %PYTHON% -c "import mace; print('       Version:', mace.__version__)"
    echo.
    echo Nothing to do. MACE is ready to use.
    echo You can now select MACE in the ATLAS-GP GUI.
    echo.
    pause
    exit /b 0
)

:: ── Upgrade pip first ────────────────────────────────────────────────────────
echo [..] Upgrading pip...
%PYTHON% -m pip install --upgrade pip --quiet
if errorlevel 1 (
    echo [WARN] pip upgrade failed, continuing anyway...
)
echo [OK]   pip ready.

:: ── Install PyTorch CPU ──────────────────────────────────────────────────────
echo.
echo [1/3] Installing PyTorch (CPU version)...
echo       This is the largest download (~800MB), please wait...
echo.
%PYTHON% -m pip install torch --index-url https://download.pytorch.org/whl/cpu --quiet
if errorlevel 1 (
    echo [WARN] PyTorch CPU install failed, trying standard install...
    %PYTHON% -m pip install torch --quiet
    if errorlevel 1 (
        echo [ERROR] PyTorch installation failed.
        echo         Check your internet connection and try again.
        pause & exit /b 1
    )
)
echo [OK]   PyTorch installed.

:: ── Install ASE ──────────────────────────────────────────────────────────────
echo.
echo [2/3] Installing ASE (Atomic Simulation Environment)...
%PYTHON% -m pip install ase --quiet
if errorlevel 1 (
    echo [ERROR] ASE installation failed.
    pause & exit /b 1
)
echo [OK]   ASE installed.

:: ── Install mace-torch ───────────────────────────────────────────────────────
echo.
echo [3/3] Installing mace-torch...
%PYTHON% -m pip install mace-torch --quiet
if errorlevel 1 (
    echo [ERROR] mace-torch installation failed.
    echo.
    echo Possible causes:
    echo   - No internet connection
    echo   - pip version too old (upgrade with: python -m pip install --upgrade pip)
    echo   - Python version incompatible (need 3.10-3.12)
    echo.
    pause & exit /b 1
)
echo [OK]   mace-torch installed.

:: ── Verify installation ──────────────────────────────────────────────────────
echo.
echo [..] Verifying installation...
%PYTHON% -c "import mace; print('[OK]   MACE version:', mace.__version__)"
if errorlevel 1 (
    echo [ERROR] MACE installed but cannot be imported.
    echo         Try restarting this installer.
    pause & exit /b 1
)

%PYTHON% -c "import torch; print('[OK]   PyTorch version:', torch.__version__)"
%PYTHON% -c "import ase; print('[OK]   ASE version:', ase.__version__)"

:: ── Write python path for ATLAS ──────────────────────────────────────────────
set "ROOT=%~dp0"
set "ROOT=%ROOT:~0,-1%"
echo %PYTHON% > "%ROOT%\ml\python_path.txt"
echo [OK]   Python path saved for ATLAS.

echo.
echo =====================================================================
echo   MACE-MP-0 installation complete!
echo =====================================================================
echo.
echo You can now:
echo   1. Run ATLAS.bat to launch the GUI
echo   2. Select "MACE" from the potential dropdown
echo   3. Or use --pot MACE from the command line
echo.
echo Note: First run will download the MACE-MP-0 model weights (~50MB)
echo       This happens automatically when MACE is first used.
echo.
pause
exit /b 0
