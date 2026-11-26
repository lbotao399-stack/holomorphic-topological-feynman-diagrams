"""Ordered BF kernels, with the graph orientation of the local reduced kernel engine."""

from collections import Counter
from fractions import Fraction
from itertools import permutations, product
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "code"))
from htk import Graph, theta, sp_reduce
import networkx as nx
import numpy as np

EPS = np.zeros((3, 3, 3), dtype=int)
for p in permutations(range(3)):
    EPS[p] = (-1) ** sum(a > b for i, a in enumerate(p) for b in p[i + 1 :])


def parity(p):
    return (-1) ** sum(a > b for i, a in enumerate(p) for b in p[i + 1 :])


def reroot(graph, origin):
    order = [v for v in range(graph.nv) if v != origin] + [origin]
    return Graph(
        len(order),
        [(order.index(u), order.index(v)) for u, v in graph.edges],
        graph.name,
    ), order


class Rooted:
    def __init__(self, graph, origin):
        self.G, self.order = reroot(graph, origin)
        self.edges = self.G.edges
        self.root = self.G.root
        self.inc = [
            [e for e, edge in enumerate(self.edges) if v in edge]
            for v in range(self.G.nv)
        ]
        self.outputs = [v for v in range(self.root) if len(self.inc[v]) == 2]
        assert all(len(self.inc[v]) in (2, 3) for v in range(self.root))
        self.geom = sp_reduce(self.G)
        net = nx.Graph()
        net.add_edges_from(self.edges)
        self.aut = sum(
            p[self.root] == self.root
            for p in nx.algorithms.isomorphism.GraphMatcher(
                net, net
            ).isomorphisms_iter()
        )
        self.basis = np.eye(3, dtype=int)

    def wick(self, root_legs):
        # Canonical vertex order: incident edge numbers, then a surviving field.
        # The root letters occur in their prescribed operator order.
        flags = [
            (v, e)
            for v in range(self.root)
            for e in self.inc[v] + ([-1] if v in self.outputs else [])
        ]
        flags.extend((self.root, e) for e in root_legs)
        pulled = [
            flags.index((v, e)) for e, ends in enumerate(self.edges) for v in ends
        ]
        pulled.extend(flags.index((v, -1)) for v in self.outputs)
        return parity(pulled)

    def orientations(self, root_legs, root_roles, ext_roles):
        # An edge is directed from its multiplier to its ordinary field.
        for tails in product(*self.edges):
            if any(
                (tails[e] == self.root) != role
                for e, role in zip(root_legs, root_roles)
            ):
                continue
            if all(
                sum(tails[e] == v for e in self.inc[v]) + int(ext_roles.get(v, False))
                == 1
                for v in range(self.root)
            ):
                yield tails

    def color(self, root_legs, root_colors, colors):
        operands = []
        for v in range(self.root):
            labels = self.inc[v][:]
            if v in self.outputs:
                labels.append(len(self.edges) + v)
            operands.extend((EPS, labels))
        for e, c in zip(root_legs, root_colors):
            operands.extend((self.basis[c], [e]))
        for v, c in colors.items():
            operands.extend((self.basis[c], [len(self.edges) + v]))
        return int(np.einsum(*operands, [], optimize=True))

    def sources(self, waves):
        p = dict(zip(self.outputs, waves))
        p[self.root] = tuple(-sum(x[k] for x in waves) for k in range(3))
        return p

    def gauge(self, root_colors, derivatives, external, root_roles=None):
        if root_roles is None:
            root_roles = [True] * len(root_colors)
        answer = Fraction()
        for perm in permutations(range(len(external))):
            # The three flags at each cubic vertex and the BF operator letters
            # are the odd holomorphic superfield components, before descendants.
            waves = [external[k][0] for k in perm]
            colors = dict(zip(self.outputs, (external[k][1] for k in perm)))
            roles = dict(zip(self.outputs, (external[k][2] for k in perm)))
            p = self.sources(waves)
            density = self.geom.density(p, check=1)
            for seq in permutations(self.inc[self.root]):
                count = sum(1 for _ in self.orientations(seq, root_roles, roles))
                coef = (
                    count
                    * self.wick(seq)
                    * parity(perm)
                    * self.color(seq, root_colors, colors)
                )
                if not coef:
                    continue
                powers = Counter()
                sg = 1
                for e, axes in zip(seq, derivatives):
                    for k in axes:
                        powers[e, k] += 1
                        if self.edges[e][0] == self.root:
                            sg = -sg
                answer += (
                    coef * sg * self.geom.integrate(p, powers=powers, dens=density)
                )
        return answer / self.aut

    def hamiltonian(
        self,
        waves,
        root_legs,
        root_roles,
        ext_roles=None,
        vertex_set=None,
        derivatives=None,
    ):
        # Ordered polynomial times reduced density. A mixed h,a,b vertex has
        # the same momentum determinant, and one Kronecker color contraction.
        p = self.sources(waves)
        density = self.geom.density(p, check=1)
        vertices = range(self.root) if vertex_set is None else vertex_set

        def numerator(y):
            val = 1
            for v in vertices:
                local = [
                    tuple((1 if self.edges[e][1] == v else -1) * t for t in y[e])
                    for e in self.inc[v]
                ]
                first, second = local[:2]
                val *= first[0] * second[1] - first[1] * second[0]
            return val

        powers = Counter()
        sign = 1
        for e, axes in zip(root_legs, derivatives or [()] * len(root_legs)):
            for k in axes:
                powers[e, k] += 1
                if self.edges[e][0] == self.root:
                    sign = -sign
        count = sum(
            1 for _ in self.orientations(root_legs, root_roles, ext_roles or {})
        )
        return (
            count
            * sign
            * self.wick(root_legs)
            * self.geom.integrate(p, powers=powers, insertion=numerator, dens=density)
        )


def rooted_two():
    return [
        Rooted(theta(*s)[0], v)
        for s, v in [((1, 3, 3), 1), ((2, 2, 3), 1), ((2, 2, 3), 3)]
    ]


def typed_color(graph, root_legs, root_letters, external, roles, pure_poisson=False):
    """All gauge/Poisson assignments with a given multiplier orientation."""
    edges = graph.edges
    rule = {
        ("a", "a", "b"): "gauge",
        ("g", "h", "h"): "gravity",
        ("a", "b", "h"): "mixed",
    }
    root_at = dict(zip(root_legs, root_letters))
    internal = [e for e in range(len(edges)) if e not in root_at]
    configurations = Counter()
    for choice in product(
        (True,) if pure_poisson else (False, True), repeat=len(internal)
    ):
        is_grav = {e: (name in ("g", "h")) for e, (name, _) in root_at.items()}
        is_grav.update(zip(internal, choice))
        insertions = []
        operands = []
        valid = True
        for v in range(graph.root):
            labels = graph.inc[v][:]
            kinds = [
                ("g" if roles[e] == v else "h")
                if is_grav[e]
                else ("b" if roles[e] == v else "a")
                for e in labels
            ]
            if v in external:
                labels.append(len(edges) + v)
                kinds.append(external[v][0])
            category = rule.get(tuple(sorted(kinds)))
            if category is None:
                valid = False
                break
            if category != "gauge":
                insertions.append(v)
            if category == "gauge":
                operands.extend((EPS, labels))
            if category == "mixed":
                operands.extend(
                    (
                        np.eye(3, dtype=int),
                        [e for e, k in zip(labels, kinds) if k in ("a", "b")],
                    )
                )
        if not valid:
            continue
        for e, (kind, col) in root_at.items():
            if kind in ("a", "b"):
                operands.extend((graph.basis[col], [e]))
        for v, (kind, col) in external.items():
            if kind in ("a", "b"):
                operands.extend((graph.basis[col], [len(edges) + v]))
        coefficient = int(np.einsum(*operands, [], optimize=True)) if operands else 1
        configurations[tuple(insertions)] += coefficient
    return configurations


def typed_symbol(graph, root_letters, derivatives, external, *, pure_poisson=False):
    if len(root_letters) != len(graph.inc[graph.root]) or len(derivatives) != len(
        root_letters
    ):
        raise ValueError("supply one field and derivative multi-index per root edge")
    if len(external) != len(graph.outputs):
        raise ValueError("supply one external field per uncontracted interaction slot")
    if pure_poisson and any(
        kind not in ("g", "h")
        for kind in [k for k, c in root_letters] + [k for p, k, c in external]
    ):
        raise ValueError("pure Poisson BF has only g and h fields")
    total = Fraction()
    for perm in permutations(range(len(external))):
        waves = [external[i][0] for i in perm]
        ext = {v: external[i][1:] for v, i in zip(graph.outputs, perm)}
        root_roles = [kind in ("b", "g") for kind, _ in root_letters]
        ext_roles = {v: kind in ("b", "g") for v, (kind, _) in ext.items()}
        sources = graph.sources(waves)
        density = graph.geom.density(sources, check=1)
        for seq in permutations(graph.inc[graph.root]):
            configurations = Counter()
            for tails in graph.orientations(seq, root_roles, ext_roles):
                configurations.update(
                    typed_color(graph, seq, root_letters, ext, tails, pure_poisson)
                )
            powers = Counter()
            sg = parity(perm) * graph.wick(seq)
            for e, axes in zip(seq, derivatives):
                for axis in axes:
                    powers[e, axis] += 1
                    if graph.edges[e][0] == graph.root:
                        sg = -sg
            for vertices, c in configurations.items():
                if not c:
                    continue

                def numerator(y):
                    val = 1
                    for v in vertices:
                        i, j = graph.inc[v][:2]
                        sign = (
                            1
                            if (graph.edges[i][0] == v) == (graph.edges[j][0] == v)
                            else -1
                        )
                        val *= sign * (y[i][0] * y[j][1] - y[i][1] * y[j][0])
                    return val

                total += (
                    sg
                    * c
                    * graph.geom.integrate(
                        sources, powers=powers, insertion=numerator, dens=density
                    )
                )
    return total / graph.aut
