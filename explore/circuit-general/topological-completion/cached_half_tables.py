"""Parameter-separated caches for the unchanged exact scalar recurrences.

Each (dimension, Laurent precision) owns different bounded caches.  The
construction records cold-route controls and restores the original module
bindings before recomputing those controls without the cache adapter.
"""
from contextlib import contextmanager
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import argparse
import json
from pathlib import Path

from build_half_tables import build, HERE
import cores
import laurent
import prop


@contextmanager
def parameter_caches(maxsize=20000):
    original = prop.set_q, prop.reduce_graph, prop.G_value, prop.cut_period, prop.Q
    caches = {}
    controls = {}
    previous = None
    stats = {"maxsize_per_parameter": maxsize, "parameter_changes": []}
    def set_q(q):
        nonlocal previous
        key = Fraction(q), laurent.ORDER
        if key != previous:
            if previous is not None and key[1] != previous[1]:
                laurent.gamma_ratio.cache_clear()
            if key not in caches:
                caches[key] = (lru_cache(maxsize=maxsize)(original[1].__wrapped__),
                               lru_cache(maxsize=maxsize)(original[2]))
            prop.Q = key[0]
            prop.reduce_graph, prop.G_value = caches[key]
            previous = key
            stats["parameter_changes"].append([str(key[0]), key[1]])
    def cut_period(edges, powers, q, cut=None, order=None):
        value = original[3](edges, powers, q, cut=cut, order=order)
        key = str(q), laurent.ORDER
        if key not in controls:
            controls[key] = (tuple(map(tuple, edges)), tuple(map(Fraction, powers)), Fraction(q), cut, order, dict(value))
        return value
    prop.set_q, prop.cut_period = set_q, cut_period
    try:
        yield stats, controls
    finally:
        stats["caches"] = {str(q)+":"+str(order): {"reduction": a.cache_info()._asdict(),
                                                     "bubble": b.cache_info()._asdict()}
                           for (q, order), (a, b) in caches.items()}
        for a, b in caches.values():
            a.cache_clear()
            b.cache_clear()
        prop.set_q, prop.reduce_graph, prop.G_value, prop.cut_period, prop.Q = original
        original[1].cache_clear()


def run(name, pole):
    with parameter_caches() as (stats, controls):
        build(name, pole)
    cold = []
    for edges, powers, q, cut, order, expected in controls.values():
        actual = prop.cut_period(edges, powers, q, cut=cut, order=order)
        if actual != expected:
            raise ArithmeticError("dimension-separated cache changed a scalar period")
        cold.append({"edges": edges, "powers": list(map(str, powers)), "q": str(q),
                     "cut": cut, "order": order,
                     "value": {cores.enc(k): str(v) for k, v in actual.items()}})
    path = HERE / "tables" / (name + "_q" + cores.qname(pole) + ".json")
    data = json.loads(path.read_text())
    data["parameter_cache"] = stats
    data["original_cache_clearing_controls"] = cold
    data["cache_adapter_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    path.write_text(json.dumps(data, separators=(",", ":")) + "\n")
    print(name, pole, "original cold controls", len(cold), "all equal", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("core", choices=("prism", "K33"))
    parser.add_argument("pole", type=Fraction)
    args = parser.parse_args()
    run(args.core, args.pole)
