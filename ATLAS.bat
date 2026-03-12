@echo off
setlocal EnableDelayedExpansion

title ATLAS-GP  ^|  Adaptive Topology-driven Landscape-Aware Search

echo.
echo =====================================================================
echo   ATLAS-GP  ^|  Adaptive Topology-driven Landscape-Aware Search
echo   Universal cluster optimizer: all 103 elements, MACE-MP-0 ready
echo =====================================================================
echo.

:: ── Project root = folder containing this .bat ───────────────────────────────
set "ROOT=%~dp0"
set "ROOT=%ROOT:~0,-1%"
set "BUILD=%ROOT%\build"
set "BIN=%BUILD%\atlas.exe"
set "MACE=%ROOT%\ml\mace_server.py"
set "GUI=%ROOT%\gui\atlas_gui.py"

if not exist "%BUILD%" mkdir "%BUILD%"
set "LOG=%BUILD%\setup_log.txt"
echo ATLAS-GP setup log - %DATE% %TIME% > "%LOG%"

echo [INFO] Project : %ROOT%

:: ── Parse flags ───────────────────────────────────────────────────────────────
set REBUILD=0
set INSTALL_MACE=0
for %%a in (%*) do (
    if /I "%%a"=="--rebuild"      set REBUILD=1
    if /I "%%a"=="--install-mace" set INSTALL_MACE=1
)

:: ─────────────────────────────────────────────────────────────────────────────
:: FAST PATH  —  binary already exists, jump straight to Python + launch
:: ─────────────────────────────────────────────────────────────────────────────
if exist "%BIN%" (
    if "%REBUILD%"=="0" (
        echo [OK]   Binary already built — skipping compilation.
        goto :phase_python
    )
    echo [INFO] --rebuild requested. Recompiling.
)

:: =============================================================================
:: PHASE 1  —  MSYS2  (auto-download + install if missing)
:: =============================================================================
echo.
echo [1/4] C++ build toolchain ^(MSYS2 + GCC + CMake^)
echo ─────────────────────────────────────────────────

set "MSYS2="
for %%p in ("C:\msys64" "C:\msys2" "D:\msys64" "%LOCALAPPDATA%\msys64" "%ProgramFiles%\msys64") do (
    if exist "%%~p\usr\bin\bash.exe" (
        set "MSYS2=%%~p"
        goto :msys2_found
    )
)

:: ── MSYS2 not found — download and silent-install ───────────────────────────
echo [..] MSYS2 not found. Downloading installer (~90 MB^)...
echo [..] This is a one-time operation.

set "MSYS2_EXE=%TEMP%\msys2-setup.exe"
set "MSYS2_URL=https://github.com/msys2/msys2-installer/releases/download/2024-01-13/msys2-x86_64-20240113.exe"

powershell -NoProfile -Command ^
  "[Net.ServicePointManager]::SecurityProtocol='Tls12';" ^
  "(New-Object Net.WebClient).DownloadFile('%MSYS2_URL%','%MSYS2_EXE%')" ^
  2>>"%LOG%"

if not exist "%MSYS2_EXE%" (
    echo [ERROR] Download failed. Check internet connection.
    echo         Install MSYS2 manually from https://www.msys2.org then re-run.
    pause & exit /b 1
)

echo [OK]   Downloaded. Installing to C:\msys64 ^(~1 min^)...
"%MSYS2_EXE%" install --root C:\msys64 --confirm-command --accept-messages 2>>"%LOG%"

if not exist "C:\msys64\usr\bin\bash.exe" (
    echo [ERROR] MSYS2 installation failed. See %LOG%
    echo         Install manually from https://www.msys2.org then re-run.
    pause & exit /b 1
)
set "MSYS2=C:\msys64"
echo [OK]   MSYS2 installed to C:\msys64

:: First-run package database update (required on fresh MSYS2)
echo [..] Updating package database...
"%MSYS2%\usr\bin\bash.exe" -lc "pacman -Sy --noconfirm" >>"%LOG%" 2>&1

:msys2_found
echo [OK]   MSYS2 : %MSYS2%
echo [OK] MSYS2 : %MSYS2% >> "%LOG%"

set "MINGW=%MSYS2%\mingw64\bin"
set "PATH=%MINGW%;%MSYS2%\usr\bin;%PATH%"

:: ── Install compiler packages (pacman skips already-installed ones) ──────────
echo [..] Installing / verifying compiler packages...
"%MSYS2%\usr\bin\bash.exe" -lc ^
  "pacman -S --noconfirm --needed mingw-w64-x86_64-gcc mingw-w64-x86_64-cmake mingw-w64-x86_64-make mingw-w64-x86_64-openmp git" ^
  >>"%LOG%" 2>&1

"%MINGW%\gcc.exe" --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] GCC not found after install. See %LOG%
    pause & exit /b 1
)
echo [OK]   Toolchain ready.

:: =============================================================================
:: PHASE 2  —  LBFGSpp  (clone if missing)
:: =============================================================================
echo.
echo [2/4] LBFGSpp dependency
echo ─────────────────────────────────────────────────

if exist "%ROOT%\external\LBFGSpp\include\LBFGSB.h" (
    echo [OK]   LBFGSpp already present.
    goto :phase_build
)

echo [..] Cloning LBFGSpp from GitHub...
"%MSYS2%\usr\bin\bash.exe" -lc ^
  "cd '%ROOT:\=/%' && git clone --depth=1 https://github.com/yixuan/LBFGSpp external/LBFGSpp" ^
  >>"%LOG%" 2>&1

if not exist "%ROOT%\external\LBFGSpp\include\LBFGSB.h" (
    echo [ERROR] LBFGSpp clone failed. Check internet connection.
    echo         Log: %LOG%
    pause & exit /b 1
)
echo [OK]   LBFGSpp cloned.

:: =============================================================================
:: PHASE 3  —  CMake configure + compile
:: =============================================================================
:phase_build
echo.
echo [3/4] Compiling ATLAS-GP  ^(first time: 2-4 min^)
echo ─────────────────────────────────────────────────

cd /d "%BUILD%"

:: Clear stale CMake cache when rebuilding
if "%REBUILD%"=="1" (
    if exist "%BUILD%\CMakeCache.txt" del /Q "%BUILD%\CMakeCache.txt"
    echo [INFO] CMake cache cleared.
)

cmake "%ROOT%" -G "MinGW Makefiles" ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DCMAKE_C_COMPILER="%MINGW%\gcc.exe" ^
    -DCMAKE_CXX_COMPILER="%MINGW%\g++.exe" ^
    -DCMAKE_MAKE_PROGRAM="%MINGW%\mingw32-make.exe" ^
    -DCMAKE_VERBOSE_MAKEFILE=OFF ^
    >>"%LOG%" 2>&1

if errorlevel 1 (
    echo [ERROR] CMake configuration failed.
    echo.
    echo Last 20 lines of build log:
    echo ─────────────────────────────────────────────────
    powershell -NoProfile -Command "Get-Content '%LOG%' -Tail 20"
    echo ─────────────────────────────────────────────────
    echo Full log: %LOG%
    pause & exit /b 1
)

echo [..] Compiling with %NUMBER_OF_PROCESSORS% threads...
"%MINGW%\mingw32-make.exe" -j%NUMBER_OF_PROCESSORS% >>"%LOG%" 2>&1

if not exist "%BIN%" (
    echo [ERROR] Build failed — atlas.exe was not produced.
    echo.
    echo Last 30 lines of build log:
    echo ─────────────────────────────────────────────────
    powershell -NoProfile -Command "Get-Content '%LOG%' -Tail 30"
    echo ─────────────────────────────────────────────────
    echo Full log: %LOG%
    pause & exit /b 1
)
echo [OK]   Build successful: %BIN%

:: ── Copy runtime DLLs ────────────────────────────────────────────────────────
echo [..] Copying runtime DLLs...
for %%d in (libgcc_s_seh-1.dll libstdc++-6.dll libwinpthread-1.dll libgomp-1.dll) do (
    if not exist "%BUILD%\%%d" (
        if exist "%MINGW%\%%d" (
            copy /Y "%MINGW%\%%d" "%BUILD%\" >nul
            echo [OK]   Copied %%d
        )
    )
)

:: =============================================================================
:: PHASE 4  —  Python  (auto-download + install if missing)
:: =============================================================================
:phase_python
echo.
echo [4/4] Python environment
echo ─────────────────────────────────────────────────

:: ── Find existing Python (must have tkinter) ────────────────────────────────
set "PYTHON="
for %%p in (python py python3) do (
    if not defined PYTHON (
        %%p --version >nul 2>&1 && (
            %%p -c "import tkinter" >nul 2>&1 && set "PYTHON=%%p"
        )
    )
)
:: Fallback: accept any Python even without tkinter (warn later)
if not defined PYTHON (
    for %%p in (python py python3) do (
        if not defined PYTHON (
            %%p --version >nul 2>&1 && set "PYTHON=%%p"
        )
    )
)

:: ── Auto-install Python 3.12 if not found ────────────────────────────────────
if not defined PYTHON (
    echo [..] Python not found. Downloading Python 3.12 installer...

    set "PY_EXE=%TEMP%\python-3.12-setup.exe"
    set "PY_URL=https://www.python.org/ftp/python/3.12.3/python-3.12.3-amd64.exe"

    powershell -NoProfile -Command ^
      "[Net.ServicePointManager]::SecurityProtocol='Tls12';" ^
      "(New-Object Net.WebClient).DownloadFile('%PY_URL%','%PY_EXE%')" ^
      2>>"%LOG%"

    if not exist "!PY_EXE!" (
        echo [ERROR] Could not download Python.
        echo         Install Python 3.10+ manually from https://www.python.org
        echo         Tick "Add Python to PATH" during install, then re-run.
        pause & exit /b 1
    )

    echo [OK]   Downloaded. Installing Python silently ^(adds to PATH^)...
    "!PY_EXE!" /quiet InstallAllUsers=0 PrependPath=1 Include_test=0 >>"%LOG%" 2>&1

    :: Reload user PATH so Python is visible in this session
    for /f "tokens=*" %%i in ('powershell -NoProfile -Command ^
      "[Environment]::GetEnvironmentVariable(\"PATH\",\"User\")"') do (
        set "PATH=%%i;%PATH%"
    )

    :: Try again after install
    for %%p in (python py python3) do (
        if not defined PYTHON (
            %%p --version >nul 2>&1 && set "PYTHON=%%p"
        )
    )

    if not defined PYTHON (
        echo [OK]   Python installed. Please close this window and re-run ATLAS.bat.
        echo        ^(Windows needs one restart to see the new PATH.^)
        pause & exit /b 0
    )
    echo [OK]   Python installed successfully.
)

for /f "tokens=*" %%v in ('%PYTHON% --version 2^>^&1') do echo [OK]   %PYTHON% : %%v

:: ── Verify tkinter ───────────────────────────────────────────────────────────
%PYTHON% -c "import tkinter" >nul 2>&1
if errorlevel 1 (
    echo [WARN] tkinter not available in this Python build.
    echo        GUI may not start. If using a custom Python, reinstall from python.org.
)

:: ── Upgrade pip silently ─────────────────────────────────────────────────────
%PYTHON% -m pip install --upgrade pip --quiet >>"%LOG%" 2>&1

:: ── Optional: install MACE-MP-0 ML potential (~2 GB) ─────────────────────────
if "%INSTALL_MACE%"=="0" goto :launch

echo.
echo [MACE] Installing MACE-MP-0 ML potential ^(~2 GB download, several minutes^)...
%PYTHON% -c "import mace" >nul 2>&1
if errorlevel 1 (
    echo [..] Installing PyTorch ^(CPU build^)...
    %PYTHON% -m pip install torch --index-url https://download.pytorch.org/whl/cpu --quiet >>"%LOG%" 2>&1
    echo [..] Installing ASE...
    %PYTHON% -m pip install ase --quiet >>"%LOG%" 2>&1
    echo [..] Installing mace-torch...
    %PYTHON% -m pip install mace-torch --quiet >>"%LOG%" 2>&1
    %PYTHON% -c "import mace" >nul 2>&1
    if errorlevel 1 (
        echo [WARN] MACE install failed. ATLAS will use built-in potentials ^(LJ/GUPTA/UFF^).
        echo        See %LOG% for details.
    ) else (
        for /f "tokens=*" %%v in ('%PYTHON% -c "import mace;print(mace.__version__)"') do (
            echo [OK]   MACE-MP-0 %%v installed.
        )
    )
) else (
    for /f "tokens=*" %%v in ('%PYTHON% -c "import mace;print(mace.__version__)"') do (
        echo [OK]   MACE-MP-0 %%v already installed.
    )
)

:: =============================================================================
:: LAUNCH  —  Start the GUI
:: =============================================================================
:launch
echo.
echo ─────────────────────────────────────────────────
echo  All checks passed. Launching ATLAS-GP GUI...
echo ─────────────────────────────────────────────────
echo  Log saved to: %LOG%
echo.

if not exist "%BIN%" (
    echo [ERROR] atlas.exe not found: %BIN%
    echo         See %LOG% for build errors.
    pause & exit /b 1
)
if not exist "%GUI%" (
    echo [ERROR] atlas_gui.py not found: %GUI%
    pause & exit /b 1
)
if not defined PYTHON (
    echo [ERROR] Python not found. Install from https://www.python.org
    pause & exit /b 1
)

start "" %PYTHON% "%GUI%" --binary "%BIN%" --mace "%MACE%"
exit /b 0
