"""Run focused checks and the complete source/row-witness certificate replay."""
import argparse
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--focused-only', action='store_true', help='omit the complete source-and-row replay')
options = parser.parse_args()
subprocess.run([sys.executable, str(HERE / 'check_sources.py')], check=True)
if not options.focused_only:
    subprocess.run([sys.executable, str(HERE / 'certify_sources.py')], check=True)
