"""Focused controls; --examples replays the full 7,368-coefficient target."""
import argparse
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--examples', action='store_true')
parser.add_argument('--dimension-shift', action='store_true',
                    help='also replay the larger unsuccessful q5-to-q7 closure')
options = parser.parse_args()
subprocess.run([sys.executable, str(HERE / 'check_target.py')], check=True)
if options.examples:
    subprocess.run([sys.executable, str(HERE / 'certify_target.py')], check=True)
if options.dimension_shift:
    subprocess.run([sys.executable, str(HERE / 'certify_target.py'), '--dimension-shift'], check=True)
