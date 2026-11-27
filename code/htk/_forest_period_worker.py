"""Private read-only finite-period worker; historical research imports stay here."""
from fractions import Fraction
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "explore/circuit-general/fiveloop"))
import per5
import gutil

_original_table = per5.table


def _read_only_table(name, q, build=True, verbose=False):
    return _original_table(name, q, build=False, verbose=False)


per5.table = _read_only_table


def evaluate(request):
    edges = tuple(tuple(map(int, edge)) for edge in request["edges"])
    q = int(request["q"])
    result = {}
    periods = 0
    for powers, coefficient in request["coefficients"]:
        powers = tuple(map(int, powers))
        coefficient = Fraction(coefficient)
        if not coefficient:
            continue
        if not gutil.convergent(edges, powers, q):
            raise ArithmeticError("nonconvergent weighted core monomial: %s, q=%d" % (powers, q))
        try:
            value = per5.period(edges, powers, q)
        except per5.UnknownCore as error:
            # Convergence was checked independently above. No missing entry
            # is interpreted as a zero, and table generation is disabled.
            return {"missing": str(error)}
        except ValueError as error:
            if str(error).startswith("divergent or unavailable:"):
                return {"missing": str(error)}
            raise
        for label, scalar in value.items():
            result[label] = result.get(label, Fraction(0)) + coefficient * scalar
        periods += 1
    return {"coefficients": {label: str(value) for label, value in result.items() if value},
            "stats": {"q": q, "weighted_periods": periods, "read_only_tables": True}}


for line in sys.stdin:
    try:
        answer = evaluate(json.loads(line))
    except Exception as error:
        answer = {"error": type(error).__name__ + ": " + str(error)}
    print(json.dumps(answer, separators=(",", ":")), flush=True)
