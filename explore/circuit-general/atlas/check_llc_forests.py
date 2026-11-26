"""Cubic product forests of A*B on C^2 (LLC Thm thm-sum-all-kernels (iii), Thm b-00405):
admissible D = 2 graphs with marked edge AB, Gamma - {A, B} a forest of trees meeting A and B, every
other vertex of degree <= 3, one degree-2 vertex (the source) per tree.
Count (Gamma, A, B) up to automorphisms (ordered and unordered {A, B}); reduce with the tight polynomial
moves (definition A: series-parallel, source-free single and paired stars), sources at the tree roots
and at A, B."""
import itertools
import time

import networkx as nx

from graphs_enum import admissible
from reduce_search import tight_search, Network, search


def forests(G):
    out = []
    for A, B in itertools.permutations(G.nodes(), 2):
        if not G.has_edge(A, B):
            continue
        H = G.copy()
        H.remove_nodes_from([A, B])
        if not nx.is_forest(H):
            continue
        if any(G.degree(v) > 3 for v in H):
            continue
        ok = True
        roots = []
        for comp in nx.connected_components(H):
            r = [v for v in comp if G.degree(v) == 2]
            if len(r) != 1:
                ok = False
                break
            if not any(G.has_edge(v, A) for v in comp) or not any(G.has_edge(v, B) for v in comp):
                ok = False
                break
            roots += r
        if ok:
            out.append((A, B, roots))
    return out


def orbit_reps(G, marks):
    """marks: list of (A, B, roots); representatives up to Aut(G) acting on ordered (A, B)."""
    GM = nx.algorithms.isomorphism.GraphMatcher(G, G)
    autos = list(GM.isomorphisms_iter())
    seen, reps = set(), []
    for A, B, roots in marks:
        if (A, B) in seen:
            continue
        orb = {(g[A], g[B]) for g in autos}
        seen |= orb
        reps.append((A, B, roots, len(autos) // len(orb)))
    return reps, autos


if __name__ == "__main__":
    for n in range(1, 7):
        t0 = time.time()
        ordered = unordered = red_ok = 0
        fails = []
        for i, G in enumerate(admissible(2, n)):
            marks = forests(G)
            if not marks:
                continue
            reps, autos = orbit_reps(G, marks)
            ordered += len(reps)
            seenu = set()
            for A, B, roots, aut in reps:
                key = frozenset((A, B))
                orbu = {frozenset((g[A], g[B])) for g in autos}
                if not (orbu & seenu):
                    unordered += 1
                seenu |= orbu
                ts = tight_search(G, roots + [A, B], 2)
                if ts["reducible"]:
                    red_ok += 1
                else:
                    fails.append((i, A, B))
        print("n=%d: product forests (ordered A,B) %d, (unordered) %d; reducible by series-parallel + source-free single/paired stars: %d/%d %s | %.1fs" % (
            n, ordered, unordered, red_ok, ordered, fails[:5], time.time() - t0), flush=True)
