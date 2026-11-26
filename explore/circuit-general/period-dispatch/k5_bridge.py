"""Read-only access to the existing exact K5 weighted-period certificates.

The K5 table format is [powers, kind, value], unlike per5's [powers, value].
This module never constructs a table, runs a solver, or substitutes zero for
an unavailable or divergent period.  Its default files are the tracked
``k5/tables/K5_q{3,...,8}.json`` certificates.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import hashlib
from itertools import combinations, groupby, permutations, product
import json
from numbers import Integral
from pathlib import Path
from types import MappingProxyType

E5 = tuple(combinations(range(5), 2))
_INDEX = {frozenset(edge): i for i, edge in enumerate(E5)}
_PERMUTATIONS = tuple(tuple(_INDEX[frozenset((p[u], p[v]))] for u, v in E5)
                      for p in permutations(range(5)))
_PERMUTATION_BY_ORDER = dict(zip(permutations(range(5)), _PERMUTATIONS))
_INCIDENT = tuple(tuple(e for e, edge in enumerate(E5) if vertex in edge) for vertex in range(5))
_TRIANGLES = tuple(tuple(_INDEX[frozenset(edge)] for edge in combinations(vertices, 2))
                   for vertices in combinations(range(5), 3))
_FOURS = tuple(tuple(_INDEX[frozenset(edge)] for edge in combinations(vertices, 2))
               for vertices in combinations(range(5), 4))
TABLE_DIRECTORY = Path(__file__).resolve().parents[1] / "k5" / "tables"
_TABLE_DIRECTORY_KEY = str(TABLE_DIRECTORY)


class K5Unavailable(LookupError):
    """A convergent period has no available certificate in the selected table."""


class K5NotConvergent(ValueError):
    """The positive Schwinger integral fails a strict convergence inequality."""


class K5TableError(ValueError):
    """An existing table does not satisfy the declared certificate format."""


def _integer(value, description):
    if isinstance(value, bool):
        raise ValueError(description + " must be an integer")
    if isinstance(value, Integral):
        return int(value)
    if isinstance(value, Fraction) and value.denominator == 1:
        return value.numerator
    raise ValueError(description + " must be an integer")


def is_k5_graph(edges):
    """Recognize a complete simple five-vertex graph, ignoring orientations."""
    try:
        edges = tuple(tuple(edge) for edge in edges)
        if len(edges) != 10 or any(len(edge) != 2 or edge[0] == edge[1] for edge in edges):
            return False
        return (len({vertex for edge in edges for vertex in edge}) == 5
                and len({frozenset(edge) for edge in edges}) == 10)
    except TypeError:
        return False


@lru_cache(maxsize=128)
def _edge_order(edges):
    if not is_k5_graph(edges):
        raise ValueError("the period bridge requires a simple complete five-vertex graph")
    try:
        vertices = sorted({vertex for edge in edges for vertex in edge})
    except TypeError as error:
        raise ValueError("vertex labels must admit a common ordering") from error
    labels = {vertex: i for i, vertex in enumerate(vertices)}
    inverse = [0] * 10
    for e, (u, v) in enumerate(edges):
        inverse[_INDEX[frozenset((labels[u], labels[v]))]] = e
    return tuple(inverse)


def _normalized_powers(edges, powers):
    edges = tuple(tuple(edge) for edge in edges)
    order = _edge_order(edges)
    powers = tuple(powers)
    if len(powers) != 10:
        raise ValueError("exactly ten edge powers are required")
    if not all(type(value) is int for value in powers):
        powers = tuple(_integer(value, "edge power") for value in powers)
    return tuple(powers[e] for e in order)


@lru_cache(maxsize=65536)
def canonical_powers(powers):
    """Original S5 minimum, retained as an independent canonicalization oracle."""
    return min(tuple(powers[i] for i in permutation) for permutation in _PERMUTATIONS)


@lru_cache(maxsize=65536)
def signature_key(powers):
    """Canonical key using exact weighted vertex signatures.

    A signature is the sorted list of the four incident weights.  Signatures
    order the vertex blocks invariantly; only permutations within equal
    blocks remain.  Minimizing over all of those permutations is a complete
    canonical form, not an isomorphism hash.  Tables are reindexed by this
    key, while their original S5 representatives are retained as provenance.
    """
    if len(set(powers)) == 1:
        return powers
    signatures = tuple(tuple(sorted(powers[e] for e in edges)) for edges in _INCIDENT)
    order = sorted(range(5), key=signatures.__getitem__)
    blocks = [tuple(vertices) for _, vertices in groupby(order, key=signatures.__getitem__)]
    if len(blocks) == 5:
        permutation = _PERMUTATION_BY_ORDER[tuple(order)]
        return tuple(powers[e] for e in permutation)
    best = None
    for orders in product(*(permutations(block) for block in blocks)):
        permutation = _PERMUTATION_BY_ORDER[tuple(vertex for block in orders for vertex in block)]
        candidate = tuple(powers[e] for e in permutation)
        if best is None or candidate < best:
            best = candidate
    return best


def convergence_failure(powers, q):
    """Return a violated K5 inequality, or None for strict convergence.

    The total weight must be 6q.  For positive integer weights, convergence
    is equivalent to a_e < q, every triangle having weight > q, and every
    four-vertex clique having weight > 3q.  Completing a connected subgraph
    on at most four vertices only increases q*h - sum(a), because q-a_e>0.
    A proper connected subgraph on all five vertices obeys the inequality
    from a_e<q and total weight 6q.  Disconnected subsets add these bounds.
    """
    if len(powers) != 10 or sum(powers) != 6 * q:
        return "the projective homogeneity condition is sum(a)=6q"
    for edge, power in enumerate(powers):
        if not 0 < power < q:
            return "edge %d must have 0 < a_e < q" % edge
    for edges in _TRIANGLES:
        if sum(powers[e] for e in edges) <= q:
            return "triangle %s must have total weight greater than q" % (edges,)
    for edges in _FOURS:
        if sum(powers[e] for e in edges) <= 3 * q:
            return "four-vertex clique %s must have total weight greater than 3q" % (edges,)
    return None


@lru_cache(maxsize=6)
def _table(q, directory):
    path = Path(directory) / ("K5_q%d.json" % q)
    if not path.is_file():
        raise K5Unavailable("no K5 period table at q=%d: %s" % (q, path))
    data = path.read_bytes()
    raw = json.loads(data)
    if (raw.get("core") != "K5" or int(raw.get("q", -1)) != q
            or tuple(tuple(edge) for edge in raw.get("edges", ())) != E5):
        raise K5TableError("K5 table metadata mismatch: %s" % path)
    values = {}
    kinds = Counter()
    for row in raw.get("values", ()):
        if not isinstance(row, list) or len(row) != 3:
            raise K5TableError("expected [powers, kind, value] in %s" % path)
        powers, kind, value = row
        powers = tuple(_integer(x, "table edge power") for x in powers)
        if len(powers) != 10 or kind not in {"T1", "T2", "G"}:
            raise K5TableError("invalid or duplicate K5 row in %s" % path)
        if convergence_failure(powers, q) is not None:
            raise K5TableError("nonconvergent K5 table row in %s" % path)
        if not isinstance(value, dict) or not value:
            raise K5TableError("a convergent positive K5 period cannot be an empty value")
        if any(not isinstance(x, (str, int)) or isinstance(x, bool) for x in value.values()):
            raise K5TableError("table coefficients must be exact integer or rational strings")
        exact = {label: Fraction(coefficient) for label, coefficient in value.items() if Fraction(coefficient)}
        if not exact:
            raise K5TableError("a convergent positive K5 period cannot be zero")
        key = signature_key(powers)
        if key in values:
            raise K5TableError("duplicate isomorphism orbit in %s" % path)
        values[key] = (kind, exact, powers)
        kinds[kind] += 1
    provenance = dict(q=q, file=str(path), sha256=hashlib.sha256(data).hexdigest(),
                      orbits=len(values), kinds=dict(sorted(kinds.items())),
                      index="exact weighted vertex signatures")
    return values, provenance


def _lookup(edges, powers, q, table_directory):
    q = _integer(q, "q")
    if q < 1:
        raise ValueError("q must be positive")
    powers = _normalized_powers(edges, powers)
    failure = convergence_failure(powers, q)
    if failure is not None:
        raise K5NotConvergent(failure)
    directory = _TABLE_DIRECTORY_KEY if table_directory is None else _directory_key(str(table_directory))
    values, provenance = _table(q, directory)
    row = values.get(signature_key(powers))
    if row is None:
        raise K5Unavailable("convergent K5 orbit %s is absent at q=%d" % (canonical_powers(powers), q))
    return row, provenance


@lru_cache(maxsize=128)
def _directory_key(directory):
    return str(Path(directory).resolve())


def k5_period_with_provenance(edges, powers, q, *, table_directory=None):
    """Return exact coefficients, orbit representative, and the table fingerprint."""
    (kind, value, representative), provenance = _lookup(edges, powers, q, table_directory)
    return dict(value=dict(value), representative=representative, kind=kind,
                table=dict(provenance))


def k5_period(edges, powers, q, *, table_directory=None):
    """I_q(K5; powers), with exact labels and Fraction coefficients.

    Missing data raises K5Unavailable; divergence raises K5NotConvergent.
    Neither case is represented by an empty coefficient dictionary.
    """
    (_, value, _), _ = _lookup(edges, powers, q, table_directory)
    return dict(value)


class K5Lookup:
    """Fixed-graph, fixed-q orbit index for bulk numerator contraction.

    ``key(powers)`` validates strict convergence and returns a collision-free
    orbit key.  ``values[key]`` is an immutable label-to-Fraction mapping.
    Integer or rational numerator coefficients may first be accumulated by
    this key, before any multiplication by the period coefficients.
    Construction loads one existing table and never attempts a solver.
    """
    def __init__(self, edges, q, *, table_directory=None):
        self.q = _integer(q, "q")
        if self.q < 1:
            raise ValueError("q must be positive")
        self.edges = tuple(tuple(edge) for edge in edges)
        self._order = _edge_order(self.edges)
        directory = _TABLE_DIRECTORY_KEY if table_directory is None else _directory_key(str(table_directory))
        rows, metadata = _table(self.q, directory)
        self.values = MappingProxyType({key: MappingProxyType(value)
                                        for key, (_, value, _) in rows.items()})
        self.metadata = MappingProxyType(dict(metadata))

    def key(self, powers):
        powers = tuple(powers)
        if len(powers) != 10:
            raise ValueError("exactly ten edge powers are required")
        if not all(type(value) is int for value in powers):
            powers = tuple(_integer(value, "edge power") for value in powers)
        powers = tuple(powers[e] for e in self._order)
        failure = convergence_failure(powers, self.q)
        if failure is not None:
            raise K5NotConvergent(failure)
        key = signature_key(powers)
        if key not in self.values:
            raise K5Unavailable("convergent K5 orbit %s is absent at q=%d" %
                                (canonical_powers(powers), self.q))
        return key

    def period(self, powers):
        return dict(self.values[self.key(powers)])


def inventory(*, table_directory=None):
    """Metadata for available tables, without any attempt to generate missing levels."""
    directory = TABLE_DIRECTORY if table_directory is None else Path(table_directory)
    out = []
    for path in sorted(directory.glob("K5_q*.json"), key=lambda p: int(p.stem.split("q")[-1])):
        q = int(path.stem.split("q")[-1])
        _, provenance = _table(q, str(directory.resolve()))
        out.append(dict(provenance))
    return out


if __name__ == "__main__":
    print(json.dumps(inventory(), indent=2))
