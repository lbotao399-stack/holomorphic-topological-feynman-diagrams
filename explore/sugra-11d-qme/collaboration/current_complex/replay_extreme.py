"""Reconstruct and verify the integral extreme-profile cycle in a temporary directory."""

import argparse
from pathlib import Path
import os
import shlex
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent


def run(args, **kwargs):
    return subprocess.run(args, check=True, **kwargs)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--rank", action="store_true", help="also verify the 174 independent boundaries"
    )
    parser.add_argument(
        "--save-cycle",
        type=Path,
        help="copy the two-column integral highest vector to this file",
    )
    args = parser.parse_args()
    with TemporaryDirectory(prefix="twist-extreme-") as directory:
        root = Path(directory)
        flags = [
            *shlex.split(os.environ.get("CXX", "c++")),
            "-O3",
            "-std=c++17",
            "-fsanitize=signed-integer-overflow",
            "-fno-sanitize-recover=all",
        ]
        fusion = root / "dual.dat"
        with fusion.open("w") as f:
            run([sys.executable, HERE / "extreme_fusion.py"], stdout=f)
        residues = []
        primes = (251, 257)
        for p in primes:
            exe = root / f"cycle-{p}"
            run([*flags, f"-DEXTREME_PRIME={p}", HERE / "extreme_cycle.cpp", "-o", exe])
            lines = run(
                [exe, "residue"], capture_output=True, text=True
            ).stdout.splitlines()
            residues.append(dict(tuple(map(int, row.split())) for row in lines))
        first, second = residues
        assert min(first) == min(second) == 4537660642477
        n = primes[0] * primes[1]
        integral = []
        for key in sorted(first.keys() | second.keys()):
            a, b = first.get(key, 0), second.get(key, 0)
            x = (
                (a + primes[0] * ((b - a) * pow(primes[0], -1, primes[1]) % primes[1]))
                * 1568
                % n
            )
            if 2 * x > n:
                x -= n
            if x:
                integral.append((key, x))
        assert len(integral) == 32676 and max(abs(x) for _, x in integral) == 27840
        text = "".join(f"{k} {x}\n" for k, x in integral)
        run([root / "cycle-251", "verify", fusion], input=text, text=True)
        if args.save_cycle:
            args.save_cycle.write_text(text)
        if args.rank:
            exe = root / "rank"
            run([*flags, HERE / "extreme_profile.cpp", "-o", exe])
            run([exe, fusion])


if __name__ == "__main__":
    main()
