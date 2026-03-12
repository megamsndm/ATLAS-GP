#!/usr/bin/env python3
"""
install_mace.py
Called by the ATLAS GUI "Install MACE" button.
Reports progress via stdout JSON lines.

Strategy:
  1. If mace is already importable in the current Python → done.
  2. If current Python >= 3.13, look for py -3.11 or py -3.12 which have
     working matscipy wheels, and install + use that instead.
  3. Otherwise install into current Python normally.

The mace_server.py is always invoked via the discovered compatible Python,
so the GUI passes the right interpreter path back via the "python" key in
the success JSON message.
"""
import sys
import subprocess
import json
import shutil

def report(status, message, **extra):
    d = {"status": status, "message": message}
    d.update(extra)
    print(json.dumps(d), flush=True)

def run(cmd, description):
    """Run a command, stream nothing, return (ok, stderr_tail)."""
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        tail = (result.stdout + result.stderr)[-600:].strip()
        report("error", f"{description} failed:\n{tail}")
        return False
    report("ok", f"{description} done")
    return True

def can_import(python_exe, package):
    try:
        r = subprocess.run(
            [python_exe, "-c",
             f"import {package}; print({package}.__version__)"],
            capture_output=True, text=True, timeout=15
        )
        return r.returncode == 0, r.stdout.strip()
    except Exception:
        return False, ""

def find_compatible_python():
    """
    Return the path to a Python 3.10-3.12 interpreter that can install
    mace-torch (matscipy has no wheels for 3.13+).
    Preference order: py -3.11, py -3.12, py -3.10
    """
    py_launcher = shutil.which("py")
    if py_launcher:
        for ver in ("3.11", "3.12", "3.10"):
            try:
                r = subprocess.run(
                    ["py", f"-{ver}", "--version"],
                    capture_output=True, text=True, timeout=5
                )
                if r.returncode == 0:
                    # Get the actual executable path
                    r2 = subprocess.run(
                        ["py", f"-{ver}", "-c",
                         "import sys; print(sys.executable)"],
                        capture_output=True, text=True, timeout=5
                    )
                    if r2.returncode == 0:
                        exe = r2.stdout.strip()
                        report("info",
                               f"Found compatible Python {ver}: {exe}")
                        return exe
            except Exception:
                continue
    return None

def install_mace(python_exe):
    """Install mace-torch into the given Python interpreter."""
    pip = [python_exe, "-m", "pip"]

    # Upgrade pip
    report("running", "Upgrading pip")
    subprocess.run(pip + ["install", "--upgrade", "pip", "--quiet"],
                   capture_output=True)
    report("ok", "Upgrading pip done")

    # PyTorch CPU
    report("running", "Installing PyTorch CPU (~800MB download, please wait)...")
    ok = run(pip + ["install", "torch",
                    "--index-url", "https://download.pytorch.org/whl/cpu",
                    "--quiet"],
             "Installing PyTorch (CPU)")
    if not ok:
        report("warn", "Trying standard PyTorch...")
        if not run(pip + ["install", "torch", "--quiet"],
                   "Installing PyTorch (standard)"):
            report("error", "PyTorch installation failed.")
            sys.exit(1)

    # ASE
    if not run(pip + ["install", "ase", "--quiet"], "Installing ASE"):
        report("error", "ASE installation failed.")
        sys.exit(1)

    # mace-torch
    if not run(pip + ["install", "mace-torch", "--quiet"],
               "Installing mace-torch"):
        report("error", "mace-torch installation failed.")
        sys.exit(1)

def main():
    report("start",
           f"Using Python: {sys.executable} "
           f"{sys.version.split()[0]} "
           f"({sys.version})")

    if sys.version_info < (3, 10):
        report("error",
               f"Python 3.10+ required. You have {sys.version}.\n"
               "Install Python 3.11 from python.org and re-try.")
        sys.exit(1)

    # ── Case 1: mace already installed in current Python ─────────────────────
    ok, ver = can_import(sys.executable, "mace")
    if ok:
        report("already_installed",
               f"MACE {ver} is already installed in {sys.executable}")
        sys.exit(0)

    # ── Case 2: current Python is 3.13+ — matscipy has no wheels yet ─────────
    target_exe = sys.executable   # default: install into current Python

    if sys.version_info >= (3, 13):
        report("warn",
               f"Python {sys.version_info.major}.{sys.version_info.minor} "
               "detected. matscipy (required by mace-torch) does not yet "
               "have wheels for this version. Looking for Python 3.11/3.12...")

        compat = find_compatible_python()
        if compat is None:
            report("error",
                   "Python 3.11 or 3.12 not found.\n\n"
                   "Please install Python 3.11 from:\n"
                   "  https://www.python.org/ftp/python/3.11.9/python-3.11.9-amd64.exe\n"
                   "During install: UNCHECK 'Add to PATH'.\n"
                   "Then click Install MACE again.")
            sys.exit(1)

        # Check if mace is already installed there
        ok, ver = can_import(compat, "mace")
        if ok:
            report("already_installed",
                   f"MACE {ver} already installed in {compat}",
                   python=compat)
            sys.exit(0)

        target_exe = compat
        report("info",
               f"Installing mace-torch into {target_exe} ...")

    # ── Install ───────────────────────────────────────────────────────────────
    install_mace(target_exe)

    # ── Verify ────────────────────────────────────────────────────────────────
    ok, ver = can_import(target_exe, "mace")
    if not ok:
        report("error",
               "Installation completed but mace cannot be imported.\n"
               "Try restarting ATLAS.")
        sys.exit(1)

    _, torch_ver = can_import(target_exe, "torch")
    _, ase_ver   = can_import(target_exe, "ase")

    report("success",
           f"MACE {ver} installed successfully!\n"
           f"Python:  {target_exe}\n"
           f"PyTorch: {torch_ver}\n"
           f"ASE:     {ase_ver}\n\n"
           f"Select 'MACE' from the potential dropdown to use it.\n"
           f"First run downloads model weights (~50 MB automatically).",
           python=target_exe)
    sys.exit(0)

if __name__ == "__main__":
    main()
