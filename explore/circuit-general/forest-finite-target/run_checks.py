"""Independent controls and the fast complete finite witness replay."""
import argparse
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--recompute', action='store_true',
                    help='also repeat bounded row generation and elimination')
options = parser.parse_args()
for script in ('check_witness.py', 'check_annihilators.py', 'certify_finite.py'):
    subprocess.run([sys.executable, str(HERE / script)], check=True)
if options.recompute:
    subprocess.run([sys.executable, str(HERE / 'certify_finite.py'), '--recompute'], check=True)
