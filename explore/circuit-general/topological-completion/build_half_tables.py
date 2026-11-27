"""Build complete half-integer scalar tables and retain their exact row systems."""
from fractions import Fraction
from hashlib import sha256
import argparse
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CORE = ROOT / "explore/circuit-general/cores"
sys.path.insert(0, str(CORE))
import cores
import msys


def build(name, pole):
    started = time.monotonic()
    system = msys.MultiSystem([cores.FAMILIES[name]], pole)
    system.build(verbose=True)
    values, pivots = system.solve(verbose=True)
    if len(values) != len(system.reps):
        raise ArithmeticError("the complete level is not determined")
    # Verify every original, uneliminated row in exact arithmetic.
    for row, expected in system.rows:
        actual = {}
        for column, coefficient in row.items():
            for label, value in values[column].items():
                actual[label] = actual.get(label, Fraction(0)) + coefficient * value
        if {k: v for k, v in actual.items() if v} != expected:
            raise ArithmeticError("an original scalar relation failed")
    encode = lambda value: {cores.enc(k): str(v) for k, v in sorted(value.items()) if v}
    data = {"schema": "complete-half-integer-core-periods-v1", "core": name,
            "q": str(pole), "edges": cores.FAMILIES[name],
            "orbit_count": len(system.reps), "rank": len(pivots),
            "row_counts": system.counts, "original_rows_verified": len(system.rows),
            "values": [[list(rep), encode(values[i])]
                       for i, (_, rep) in enumerate(system.reps)],
            "rows": [[[ [int(k), str(v)] for k, v in sorted(row.items())], encode(const)]
                     for row, const in system.rows],
            "seconds": round(time.monotonic() - started, 3),
            "dependencies": {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest()
                             for p in [Path(__file__)] + sorted(CORE.glob("*.py"))}}
    output = HERE / "tables" / (name + "_q" + cores.qname(pole) + ".json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(data, separators=(",", ":")) + "\n")
    print(name, str(pole), "complete", len(values), "orbits", len(system.rows),
          "exact rows", "seconds", data["seconds"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("core", choices=tuple(cores.FAMILIES))
    parser.add_argument("poles", nargs="+", type=Fraction)
    args = parser.parse_args()
    for pole in args.poles:
        build(args.core, pole)
