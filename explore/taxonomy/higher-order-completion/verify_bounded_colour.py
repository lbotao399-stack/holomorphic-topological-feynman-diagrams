"""Actual spectator tensors and overflow-fallback controls, all exact."""
from hashlib import sha256
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "qn-holomorphic"))
import numpy as np
import audit_spectators
import exact_colour
import safe_integer_contraction as safe


def main():
    samples = [
        ("ij,jk->ik", [np.array([[1, 2], [-3, 4]], dtype=object),
                       np.array([[2, -4], [3, 5]], dtype=object)]),
        ("i,i->", [np.array([2**40, 2**40], dtype=object),
                   np.array([2**30, 2**30], dtype=object)]),
        ("ii->", [np.array([[2**63, 1], [1, 2**63]], dtype=object)]),
        ("ij,jk,ki->", [np.array([[1, 2], [-3, 4]], dtype=object)]*3),
    ]
    for specification, arrays in samples:
        expected = np.einsum(specification, *arrays, optimize=False)
        assert np.array_equal(safe.einsum(specification, *arrays), expected)
    bundle = audit_spectators.weights.inputs(1)
    expected = [audit_spectators.history(bundle, "G2.3u", True),
                audit_spectators.history(bundle, "G2.1", False)]
    safe.install(exact_colour)
    actual = [audit_spectators.history(bundle, "G2.3u", True),
              audit_spectators.history(bundle, "G2.1", False)]
    assert actual == expected
    assert safe.STATISTICS["bounded_integer_steps"] and safe.STATISTICS["python_integer_steps"]
    report = {"passed": True, "synthetic_exact_controls": len(samples),
              "complete_original_spectator_tensor_matches": True,
              "complete_tensor_without_spectator_matches": True,
              "statistics": safe.STATISTICS,
              "sha256": {name: sha256((HERE/name).read_bytes()).hexdigest() for name in
                         ("safe_integer_contraction.py", "exact_colour.py", "verify_bounded_colour.py")}}
    (HERE/"bounded_colour_verification.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
