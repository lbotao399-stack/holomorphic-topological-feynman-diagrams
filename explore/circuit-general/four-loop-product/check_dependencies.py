"""Verify the complete recorded runtime and exact period-data fingerprints."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def main():
    data = json.loads((HERE / 'dependency_manifest.json').read_text())
    for relative, expected in data['sha256'].items():
        actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        if actual != expected:
            raise ArithmeticError('a recorded dependency changed: ' + relative)
    print('Verified', len(data['sha256']), 'runtime, exact-data, and complete-certificate fingerprints.')


if __name__ == '__main__':
    main()
