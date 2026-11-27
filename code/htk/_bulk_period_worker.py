"""Isolated, read-only scalar lookup for the general bulk dispatcher."""
from fractions import Fraction
from functools import lru_cache
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "explore/circuit-general/fiveloop"),
               str(ROOT / "explore/circuit-general/cores")]
import per5
import cores
import gutil
import prop
import laurent


class MissingPeriod(Exception):
    pass


old_five_table = per5.table
per5.table = lambda name, q, build=True, verbose=False: old_five_table(name, q, build=False, verbose=False)


@lru_cache(maxsize=None)
def read_core_table(name, q, verbose=False):
    filename = name + "_q" + cores.qname(q) + ".json"
    choices = (ROOT / "explore/circuit-general/topological-completion/tables" / filename,
               ROOT / "explore/circuit-general/cores/tables" / filename)
    path = next((path for path in choices if path.is_file()), None)
    if path is None:
        raise MissingPeriod(filename)
    data = json.loads(path.read_text())
    return {tuple(a): {cores.dec(k, Fraction(q).denominator != 1): Fraction(v)
                      for k, v in value.items()} for a, value in data["values"]}


cores.table = read_core_table
original_reduce, original_bubble = prop.reduce_graph.__wrapped__, prop.G_value
parameter_caches = {}
last_parameter = None


def set_q(q):
    global last_parameter
    parameter = Fraction(q), laurent.ORDER
    if parameter != last_parameter:
        if last_parameter is not None and parameter[1] != last_parameter[1]:
            laurent.gamma_ratio.cache_clear()
        if parameter not in parameter_caches:
            parameter_caches[parameter] = (lru_cache(maxsize=20000)(original_reduce),
                                           lru_cache(maxsize=20000)(original_bubble))
        prop.Q = parameter[0]
        prop.reduce_graph, prop.G_value = parameter_caches[parameter]
        last_parameter = parameter


prop.set_q = set_q


def evaluate(request):
    edges = tuple(tuple(map(int, edge)) for edge in request["edges"])
    q = Fraction(request["q"])
    integer = q.denominator == 1
    # The inherited helper is explicitly for simple graphs: on parallel
    # edges it repeats the last edge index and is not a permutation.
    # Keeping every labelled scalar term is exact for multigraphs.
    simple = len({tuple(sorted(edge)) for edge in edges}) == len(edges)
    autos = gutil.edge_automorphisms(edges) if simple else (tuple(range(len(edges))),)
    coefficients = {}
    for powers, coefficient in request["coefficients"]:
        representative = gutil.orbit_rep(tuple(map(int, powers)), autos)
        coefficients[representative] = coefficients.get(representative, Fraction(0)) + Fraction(coefficient)
    value = {}
    periods = 0
    for alpha, coefficient in coefficients.items():
        if not coefficient:
            continue
        if not gutil.convergent(edges, list(alpha), q):
            raise ArithmeticError("a complete scalar orbit is not convergent")
        try:
            answer = per5.period(edges, alpha, int(q)) if integer else cores.period(edges, alpha, q)
        except (MissingPeriod, per5.UnknownCore, prop.NotReducible) as error:
            return {"missing": str(error)}
        except ValueError as error:
            if str(error).startswith(("divergent or unavailable:", "divergent or unrecognized core")):
                return {"missing": str(error)}
            raise
        for label, scalar in answer.items():
            label = label if isinstance(label, str) else cores.enc(label)
            value[label] = value.get(label, Fraction(0)) + coefficient*scalar
        periods += 1
    return {"coefficients": {k: str(v) for k, v in value.items() if v},
            "stats": {"weighted_scalar_orbits": periods, "read_only_tables": True, "pole": str(q)}}


for line in sys.stdin:
    try:
        answer = evaluate(json.loads(line))
    except Exception as error:
        answer = {"error": type(error).__name__ + ": " + str(error)}
    print(json.dumps(answer, separators=(",", ":")), flush=True)
