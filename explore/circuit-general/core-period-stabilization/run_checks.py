"""Replay the independent exact theorem controls and their source digests."""

from pathlib import Path
import subprocess
import sys


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    subprocess.run([sys.executable, str(here / "verify.py")], check=True)
