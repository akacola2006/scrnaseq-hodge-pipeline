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
    candidates = [
        # Common Windows Python locations
        shutil.which("python"),
        shutil.which("python3"),
    ]

    # Check known paths
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

    # Also check PATH-based pythoncore
    appdata = os.environ.get("LOCALAPPDATA", "")
    if appdata:
        for d in sorted(Path(appdata).glob("Python/pythoncore-*/python.exe"), reverse=True):
            candidates.append(str(d))

    for c in candidates:
        if c and Path(c).exists() and not str(c).endswith("scRNAseq_Hodge_Pipeline.exe"):
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

    print("ERROR: Could not find Python installation.")
    print("Please install Python 3.10+ from https://python.org")
    input("Press Enter to exit...")
    sys.exit(1)


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
    main()
