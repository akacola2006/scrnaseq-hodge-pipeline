"""
scRNAseq Hodge Pipeline — Windows Launcher
============================================
Double-click this (or the .exe) to start the Streamlit app.
Automatically finds system Python, checks dependencies, and opens the browser.
"""
import os
import shutil
import subprocess
import sys
import time
import webbrowser
from pathlib import Path


# When running from PyInstaller exe, __file__ points to a temp folder.
# Use the exe's actual location instead.
if getattr(sys, "frozen", False):
    PIPELINE_ROOT = Path(sys.executable).resolve().parent
else:
    PIPELINE_ROOT = Path(__file__).resolve().parent
APP_PY = PIPELINE_ROOT / "app.py"
PORT = 8501
URL = f"http://localhost:{PORT}"


def find_python() -> str:
    """Find the real Python interpreter (not the PyInstaller exe).

    When running from a PyInstaller bundle, sys.executable points to the
    .exe itself, not to a Python interpreter. We need to find the actual
    Python installation.
    """
    # If running as normal .py script, sys.executable is fine
    if not getattr(sys, "frozen", False):
        return sys.executable

    # Running from PyInstaller exe — search for real Python
    # Check known install locations FIRST (before PATH, which may find WindowsApps stub)
    candidates = []

    appdata = os.environ.get("LOCALAPPDATA", "")
    if appdata:
        # pythoncore (python.org installer, newer style)
        for d in sorted(Path(appdata).glob("Python/pythoncore-*/python.exe"), reverse=True):
            candidates.append(str(d))

    for base in [
        Path(os.environ.get("LOCALAPPDATA", "")) / "Programs" / "Python",
        Path(os.environ.get("LOCALAPPDATA", "")) / "Python",
        Path("C:/Python"),
        Path("C:/Program Files/Python"),
    ]:
        if base.exists():
            for d in sorted(base.iterdir(), reverse=True):
                exe = d / "python.exe"
                if exe.exists():
                    candidates.append(str(exe))

    # PATH-based as fallback (may find WindowsApps stub)
    for which_name in ["python", "python3"]:
        found = shutil.which(which_name)
        if found:
            candidates.append(found)

    for c in candidates:
        if (c and Path(c).exists()
                and "scRNAseq_Hodge_Pipeline" not in str(c)
                and "WindowsApps" not in str(c)):
            # Verify it's a real Python
            try:
                result = subprocess.run(
                    [c, "--version"],
                    capture_output=True, text=True, timeout=5,
                )
                if result.returncode == 0 and "Python" in result.stdout:
                    return c
            except Exception:
                continue

    # No Python found — offer to install automatically
    print("=" * 50)
    print("  Python not found on this PC.")
    print("=" * 50)
    print()
    print("Options:")
    print("  1. Auto-install Python now (recommended)")
    print("  2. Exit and install manually from https://python.org")
    print()
    choice = input("Enter 1 or 2: ").strip()

    if choice == "1":
        installed_python = _auto_install_python()
        if installed_python:
            return installed_python

    print("\nPlease install Python 3.10+ from https://python.org")
    print("Make sure to check 'Add Python to PATH' during installation.")
    input("Press Enter to exit...")
    sys.exit(1)


def _auto_install_python() -> str:
    """Download and install Python automatically."""
    import urllib.request
    import tempfile

    PYTHON_VERSION = "3.12.8"
    INSTALLER_URL = f"https://www.python.org/ftp/python/{PYTHON_VERSION}/python-{PYTHON_VERSION}-amd64.exe"

    print(f"\nDownloading Python {PYTHON_VERSION}...")

    try:
        installer_path = Path(tempfile.gettempdir()) / f"python-{PYTHON_VERSION}-installer.exe"

        # Download with progress
        def _report(block, block_size, total):
            downloaded = block * block_size
            if total > 0:
                pct = min(downloaded / total * 100, 100)
                bar = "#" * int(pct // 5)
                print(f"\r  [{bar:<20}] {pct:.0f}%", end="", flush=True)

        urllib.request.urlretrieve(INSTALLER_URL, str(installer_path), _report)
        print("\n  Download complete.")

        # Run installer silently with PATH enabled
        print("  Installing Python (this may take a minute)...")
        result = subprocess.run(
            [
                str(installer_path),
                "/quiet",
                "InstallAllUsers=0",
                "PrependPath=1",
                "Include_pip=1",
                "Include_test=0",
            ],
            timeout=300,
        )

        # Clean up installer
        try:
            installer_path.unlink()
        except Exception:
            pass

        if result.returncode == 0:
            print("  Python installed successfully!")
            print()

            # Find the newly installed Python
            appdata = os.environ.get("LOCALAPPDATA", "")
            new_candidates = []
            if appdata:
                for d in sorted(Path(appdata).glob("Programs/Python/Python*/python.exe"), reverse=True):
                    new_candidates.append(str(d))
                for d in sorted(Path(appdata).glob("Python/pythoncore-*/python.exe"), reverse=True):
                    new_candidates.append(str(d))

            for c in new_candidates:
                try:
                    r = subprocess.run([c, "--version"], capture_output=True, text=True, timeout=5)
                    if r.returncode == 0 and "Python" in r.stdout:
                        return c
                except Exception:
                    continue

            # Try PATH refresh
            refreshed = shutil.which("python")
            if refreshed and "WindowsApps" not in refreshed:
                return refreshed

            print("  Python installed but could not locate it.")
            print("  Please restart this application.")
            input("Press Enter to exit...")
            sys.exit(0)
        else:
            print(f"  Installation failed (exit code {result.returncode}).")
            return None

    except Exception as e:
        print(f"\n  Download/install failed: {e}")
        return None


def check_and_install_deps(python: str):
    """Check and install required packages."""
    required = {
        "streamlit": "streamlit",
        "plotly": "plotly",
        "pandas": "pandas",
        "numpy": "numpy",
        "yaml": "pyyaml",
        "anndata": "anndata",
        "scipy": "scipy",
        "sklearn": "scikit-learn",
        "statsmodels": "statsmodels",
    }

    missing = []
    for import_name, pip_name in required.items():
        result = subprocess.run(
            [python, "-c", f"import {import_name}"],
            capture_output=True, timeout=10,
        )
        if result.returncode != 0:
            missing.append(pip_name)

    if missing:
        print(f"\nInstalling: {', '.join(missing)}")
        subprocess.check_call(
            [python, "-m", "pip", "install"] + missing,
        )
        print("Installation complete.\n")
    else:
        print("All dependencies OK.\n")


def kill_existing_streamlit():
    """Kill any existing Streamlit process on the port."""
    try:
        if sys.platform == "win32":
            result = subprocess.run(
                ["netstat", "-ano"],
                capture_output=True, text=True,
            )
            for line in result.stdout.split("\n"):
                if f":{PORT}" in line and "LISTENING" in line:
                    parts = line.split()
                    pid = parts[-1]
                    try:
                        subprocess.run(["taskkill", "/F", "/PID", pid],
                                     capture_output=True)
                        print(f"Killed existing process on port {PORT}")
                    except Exception:
                        pass
    except Exception:
        pass


def start_streamlit(python: str):
    """Start the Streamlit server."""
    print(f"Starting Streamlit on {URL} ...")
    print("(Close this window to stop the server)\n")

    proc = subprocess.Popen(
        [
            python, "-m", "streamlit", "run",
            str(APP_PY),
            "--server.port", str(PORT),
            "--server.headless", "true",
            "--browser.gatherUsageStats", "false",
        ],
        cwd=str(PIPELINE_ROOT),
    )

    time.sleep(3)
    webbrowser.open(URL)
    print(f"Browser opened: {URL}")
    print("Press Ctrl+C or close this window to stop.\n")

    try:
        proc.wait()
    except KeyboardInterrupt:
        print("\nShutting down...")
        proc.terminate()
        proc.wait(timeout=5)
        print("Done.")


def main():
    print("=" * 50)
    print("  scRNAseq Hodge Decomposition Pipeline")
    print("=" * 50)
    print()

    python = find_python()
    print(f"Python: {python}")

    check_and_install_deps(python)
    kill_existing_streamlit()
    start_streamlit(python)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        input("\nPress Enter to exit...")
        sys.exit(1)
