"""Exact full-source unisolvence, before integration or covariance reduction.

Every point and a nonsingular integer minor are retained.  The only input
is the already specified set of complete prism source evaluations.
"""

from hashlib import sha256
from itertools import permutations
from math import prod
from pathlib import Path
import argparse
import json

from flint import fmpq_mat, fmpz_mat

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / "explore/circuit-general/topological-completion/prism_source_covariance.json"
OUTPUT = HERE / "source_orbit_certificate.json"


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def homogeneous(variables, degree):
    if variables == 1:
        yield (degree,)
    else:
        for first in range(degree+1):
            for tail in homogeneous(variables-1, degree-first):
                yield (first,)+tail


def calculate():
    source = json.loads(SOURCE.read_text())
    edges = {tuple(sorted(e)) for e in source["edges"]}
    actions = []
    for permutation in permutations(range(6)):
        if {tuple(sorted((permutation[a], permutation[b]))) for a, b in edges} == edges:
            actions.append(permutation)
    if len(actions) != 12:
        raise ArithmeticError("the independently enumerated prism automorphism group changed")
    points = set()
    for sample in source["samples"]:
        full = tuple(sample["sources"])+(-sum(sample["sources"]),)
        for action in actions:
            inverse = tuple(action.index(i) for i in range(6))
            points.add(tuple(full[inverse[i]] for i in range(5)))
    points = sorted(points)
    if len(points) != 162:
        raise ArithmeticError("the complete prism source orbit has changed")
    records = []
    for degree, size in ((4, 70), (5, 126)):
        monomials = list(homogeneous(5, degree))
        entries = [[prod(x**p for x, p in zip(point, exponent)) for exponent in monomials]
                   for point in points]
        reduced, rank = fmpq_mat(entries).transpose().rref()
        if rank != size:
            raise ArithmeticError("the source orbit does not span the complete homogeneous polynomial space")
        indices = [next(j for j in range(len(points)) if reduced[i, j]) for i in range(rank)]
        determinant = fmpz_mat([entries[i] for i in indices]).det()
        if not determinant:
            raise ArithmeticError("the recorded full-source minor is singular")
        records.append({"degree": degree, "all_monomials": monomials, "rank": rank,
                        "independent_source_row_indices": indices,
                        "exact_integer_minor_determinant": str(determinant),
                        "complete_evaluation_matrix_sha256": sha256(canonical(entries)).hexdigest()})
    return {"schema": "full-source-orbit-unisolvence-v1", "core": "prism",
            "base_source_points": [s["sources"] for s in source["samples"]],
            "all_vertex_automorphisms": actions, "all_independent_source_rows": points,
            "root_source": "minus the sum of all five independent original sources",
            "complete_homogeneous_spaces": records,
            "sources": {str(SOURCE.relative_to(ROOT)): sha256(SOURCE.read_bytes()).hexdigest(),
                        str(Path(__file__).relative_to(ROOT)): sha256(Path(__file__).read_bytes()).hexdigest()}}


def verify(write=False):
    data = calculate()
    if write:
        OUTPUT.write_bytes(canonical(data)+b"\n")
    elif json.loads(OUTPUT.read_text()) != json.loads(canonical(data)):
        raise ArithmeticError("the full-source convergence rank certificate changed")
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = verify(args.write)
    print("Prism: 14 original source points; 12 automorphisms; 162 distinct balanced source rows.")
    print("Complete degree-four rank 70 and degree-five rank 126; both exact integer minors nonzero.")
    print("canonical_certificate_sha256:", sha256(canonical(data)).hexdigest())
