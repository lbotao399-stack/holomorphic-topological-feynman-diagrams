"""Replay the stored complete identities and independent production controls."""
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
for path in (ROOT / 'code/tests/test_core_sector_periods.py', HERE / 'certify.py'):
    subprocess.run([sys.executable, str(path)], cwd=ROOT, check=True)
