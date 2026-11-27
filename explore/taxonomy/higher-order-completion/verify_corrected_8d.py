"""Replay the corrected finite-field functional systems without omissions.

The color contraction is integral. The kernel values in this audit are
modulo the stated prime; this is not a universal rational nilpotency proof.
"""
import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PRIME = 2147483647


def residue(value):
    value = Fraction(value)
    return value.numerator * pow(value.denominator, -1, PRIME) % PRIME


def rank(rows):
    rows = [list(row) for row in rows]
    count = 0
    for column in range(len(rows[0])):
        pivot = next((i for i in range(count, len(rows)) if rows[i][column] % PRIME), None)
        if pivot is None:
            continue
        rows[count], rows[pivot] = rows[pivot], rows[count]
        inverse = pow(rows[count][column], -1, PRIME)
        rows[count] = [value*inverse % PRIME for value in rows[count]]
        for i in range(count+1, len(rows)):
            factor = rows[i][column]
            rows[i] = [(x-factor*y) % PRIME for x, y in zip(rows[i], rows[count])]
        count += 1
    return count


def read(path, functionals):
    data = json.loads(path.read_text())
    assert data["P"] == PRIME and data["F"] == functionals
    assert all(len(row) == functionals for row in data["E"].values())
    return data["E"]


def order_two():
    rows, rhs, paths = [], [], []
    required = {"G1.1|G1.1"}
    required |= {"Q0|G2.%d" % i for i in (1, 2, 3)}
    required |= {"G2.%d|Q0" % i for i in (1, 2, 3)}
    for tag in ("R3", "R5"):
        path = HERE/("corrected_8d_order2_"+tag+".json")
        values = read(path, 6)
        assert set(values) == required, "the entire order-two composition list is required"
        paths.append(path)
        for f in range(6):
            rows.append([(values["Q0|G2.%d" % i][f]+values["G2.%d|Q0" % i][f]) % PRIME
                         for i in (1, 2, 3)])
            rhs.append(-residue(Fraction(1, 4))*values["G1.1|G1.1"][f] % PRIME)
    weights = [Fraction(-1, 2), Fraction(-1), Fraction(-1, 6)]
    residuals = [(sum(a*residue(x) for a, x in zip(row, weights))-b) % PRIME
                 for row, b in zip(rows, rhs)]
    assert residuals == [0]*12
    assert rank(rows) == rank([row+[b] for row, b in zip(rows, rhs)]) == 3
    return {"all_compositions_per_family": 7, "complete_functionals": 12,
            "matrix": rows, "rhs": rhs, "rank": 3, "augmented_rank": 3,
            "weights": list(map(str, weights)), "residuals": residuals}, paths


def order_three():
    original = ROOT/"explore/taxonomy/b4-round2/o3g4/A_all.json"
    old = read(original, 3)
    # G3 consumes all four original inputs in Q0|G3. Thus this branch has
    # no untouched original color, and is unchanged by the repaired rule.
    values = {key: value for key, value in old.items() if key.startswith("Q0|G3.")}
    assert len(values) == 12
    paths = [original]
    corrected = {}
    for path in sorted(HERE.glob("corrected_8d_order3_*.json")):
        if path.name.endswith("_all.json"):
            continue
        rows = read(path, 3)
        for key, value in rows.items():
            if key in corrected:
                assert corrected[key] == value
            corrected[key] = value
        paths.append(path)
    required = {"G3.%d|Q0" % i for i in range(1, 13)}
    required |= {"G1.1|G2.%d" % i for i in (1, 2, 3)}
    required |= {"G2.%d|G1.1" % i for i in (1, 2, 3)}
    assert set(corrected) == required, ("incomplete corrected compositions", sorted(required-set(corrected)))
    values.update(corrected)
    names = ["G3.%d" % i for i in range(1, 13)]
    automorphisms = [2, 1, 2, 1, 1, 2, 2, 1, 2, 2, 8, 1]
    rows = [[(values["Q0|"+name][f]+values[name+"|Q0"][f]) % PRIME for name in names]
            for f in range(3)]
    lower = {"G1.1": Fraction(1, 2), "G2.1": Fraction(-1, 2),
             "G2.2": Fraction(-1), "G2.3": Fraction(-1, 6)}
    rhs = [sum(-residue(lower[key.split("|")[0]]*lower[key.split("|")[1]])*value[f]
               for key, value in values.items() if "Q0" not in key) % PRIME for f in range(3)]
    survivors = []
    for signs in product((1, -1), repeat=12):
        weights = [residue(Fraction(s, a)) for s, a in zip(signs, automorphisms)]
        residuals = [(sum(a*x for a, x in zip(row, weights))-b) % PRIME for row, b in zip(rows, rhs)]
        if not any(residuals):
            survivors.append(list(signs))
    assert rank(rows) == rank([row+[b] for row, b in zip(rows, rhs)])
    return {"all_compositions": len(values), "corrected_compositions": len(corrected),
            "unaffected_original_input_compositions": 12,
            "complete_functionals": 3, "matrix": rows, "rhs": rhs,
            "rank": rank(rows), "augmented_rank": rank([row+[b] for row, b in zip(rows, rhs)]),
            "complete_sign_patterns_checked": 4096, "surviving_sign_patterns": survivors,
            "changed_complete_compositions": sorted(key for key in corrected if corrected[key] != old[key]),
            "unrestricted_weights_remain_underdetermined": rank(rows) < 12}, paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order-two-only", action="store_true")
    args = parser.parse_args()
    second, paths = order_two()
    report = {"passed": True, "prime": PRIME, "kernel_arithmetic": "finite field",
              "color_arithmetic": "exact integers with prechecked intermediate bounds",
              "universal_rational_nilpotency_claimed": False, "order_two": second}
    if not args.order_two_only:
        third, extra = order_three()
        paths.extend(extra)
        report["order_three"] = third
    paths.extend(HERE/name for name in ("exact_colour.py", "safe_integer_contraction.py",
                                        "b4_modular_probe.py", "b4_bounded_integer_probe.py",
                                        "verify_corrected_8d.py"))
    report["input_sha256"] = {str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest()
                              for path in paths}
    target = "corrected_8d_order2_verification.json" if args.order_two_only else "corrected_8d_verification.json"
    (HERE/target).write_text(json.dumps(report, indent=2)+"\n")
    print("PASS:", target, "second-order rank", second["rank"], flush=True)
    if "order_three" in report:
        print("third-order rank", third["rank"], "sign patterns", third["surviving_sign_patterns"],
              "changed compositions", third["changed_complete_compositions"], flush=True)


if __name__ == "__main__":
    main()
