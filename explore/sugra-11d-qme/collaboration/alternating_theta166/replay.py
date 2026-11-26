"""Rebuild the polarized theta166 graph from the graph and vertex rules."""
from fractions import Fraction
from pathlib import Path
import os
import shlex
import shutil
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
EXPECTED = (
    Fraction(65807437, 44157785817600),
    Fraction(-113467453, 44157785817600),
    Fraction(-3773629, 766628226000),
)


def run(args):
    result = subprocess.run(args, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode:
        raise RuntimeError(f"{args[0]} failed ({result.returncode})\n{result.stderr}\n{result.stdout}")
    return result.stdout


def gmp_flags():
    if shutil.which("pkg-config"):
        info = subprocess.run(["pkg-config", "--cflags", "--libs", "gmpxx"],
                              text=True, capture_output=True)
        if info.returncode == 0:
            return shlex.split(info.stdout)
    if shutil.which("brew"):
        result = subprocess.run(["brew", "--prefix", "gmp"], text=True, capture_output=True)
        if result.returncode == 0:
            prefix = Path(result.stdout.strip())
            return [f"-I{prefix / 'include'}", f"-L{prefix / 'lib'}", "-lgmpxx", "-lgmp"]
    return ["-lgmpxx", "-lgmp"]


def calculate(scratch):
    cxx = shlex.split(os.environ.get("CXX", "c++"))
    gmp = gmp_flags()
    for name in ("assemble", "moments", "independent_moments", "symmetry_numeric"):
        args = [*cxx, "-std=c++17", "-O3", HERE / f"{name}.cpp"]
        if "moments" in name:
            args += gmp
        run([*args, "-o", scratch / name])
    values = []
    for slot, expected in enumerate(EXPECTED, start=1):
        poly = scratch / f"B{slot}.poly"
        run([scratch / "assemble", str(slot-1), "0", "120", poly])
        first = run([scratch / "moments", poly]).splitlines()
        second = run([scratch / "independent_moments", poly]).splitlines()
        if len(first) != 121 or len(second) != 122:
            raise ValueError("incomplete integral output")
        summed = Fraction(0)
        for order in range(120):
            parts, independent = first[order].split(), second[order].split()
            if (int(parts[1]), int(parts[3])) != (order, int(independent[1])) or int(independent[0]) != order:
                raise ValueError(f"B{slot}, incompatible order or polynomial term count")
            value = Fraction(parts[-1])
            if value != Fraction(independent[-1]) / 158400:
                raise ValueError(f"B{slot}, order {order}, integrals disagree")
            summed += value
        if any(Fraction(text.split()[-1]) != expected for text in (first[-1], second[-1])) or summed != expected:
            raise ValueError(f"B{slot}, integral differs from recorded result")
        for seed in (7, 13):
            check = run([scratch / "symmetry_numeric", str(slot), poly, str(seed), "-1"]).strip()
            if ("120 blocks x 384 assignments x 4 Wick orientations" not in check
                    or "dropping exterior sign fails 120 blocks" not in check):
                raise ValueError("direct Wick evaluation or wrong-sign control failed")
            print(check, flush=True)
        print(f"A_B{slot}/pi = {summed}", flush=True)
        values.append(summed)
    orbit_sum = 4*values[0] + 4*values[1] + 2*values[2]
    if orbit_sum != Fraction(-12644, 892812375):
        raise ValueError("unexpected gamma orbit sum")
    print(f"All placements / pi = {orbit_sum}")
    print(f"All placements / (4*pi) = {orbit_sum/4}")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory(prefix="theta166-") as scratch:
        calculate(Path(scratch))
