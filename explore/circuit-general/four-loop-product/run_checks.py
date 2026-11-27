"""Replay the complete fourth-order product and its independent controls."""
from pathlib import Path
import argparse
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--recompute', action='store_true',
                        help='repeat all 350 exact sourced graph evaluations')
    options = parser.parse_args()
    commands = [('check_wick.py',), ('check_geometries.py',),
                ('check_dependencies.py',), ('render_tensor.py', '--check')]
    if options.recompute:
        commands.append(('quartic_b_squared.py', '--recompute'))
    for command in commands:
        print('\nRunning ' + ' '.join(command), flush=True)
        subprocess.run([sys.executable, str(HERE / command[0]), *command[1:]], check=True)
    print('\nComplete four-loop product checks passed.', flush=True)


if __name__ == '__main__':
    main()
