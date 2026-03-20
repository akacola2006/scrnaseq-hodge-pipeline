"""
scRNAseq Hodge Pipeline — Windows Launcher
============================================
Double-click this (or the .exe) to start the Streamlit app.
Automatically checks dependencies, installs missing ones, and opens the browser.
"""
import os
import subprocess
import sys
import time
import webbrowser
from pathlib import Path


PIPELINE_ROOT = Path(__file__).resolve().parent
APP_PY = PIPELINE_ROOT / "app.py"
PORT = 8501
URL = f"http://localhost:{PORT}"


def check_python():
    """Verify Python is available."""
    print(f"Python: {sys.executable} ({sys.version})")
    return True


def check_and_install_deps():
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
        try:
            __import__(import_name)
        except ImportError:
            missing.append(pip_name)

    if missing:
        print(f"\nInstalling missing packages: {', '.join(missing)}")
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install"] + missing,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
        )
        print("Installation complete.")
    else:
        print("All dependencies OK.")


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
                        print(f"Killed existing process on port {PORT} (PID {pid})")
                    except Exception:
                        pass
    except Exception:
        pass


def start_streamlit():
    """Start the Streamlit server."""
    print(f"\nStarting Streamlit on {URL} ...")
    print("(Close this window to stop the server)\n")

    proc = subprocess.Popen(
        [
            sys.executable, "-m", "streamlit", "run",
            str(APP_PY),
            "--server.port", str(PORT),
            "--server.headless", "true",
            "--browser.gatherUsageStats", "false",
        ],
        cwd=str(PIPELINE_ROOT),
    )

    # Wait for server to start, then open browser
    time.sleep(3)
    webbrowser.open(URL)

    print(f"Browser opened: {URL}")
    print("Press Ctrl+C to stop.\n")

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

    check_python()
    check_and_install_deps()
    kill_existing_streamlit()
    start_streamlit()


if __name__ == "__main__":
    main()
