"""Task 3: graph sets selected by theory rules, classified by circuit reduction.

Sets (vertices of the interaction are cubic unless stated; sources = vertices with outputs or legs):
  bulk3(D, n)   cubic bulk graphs: max degree <= 3, legs at the bivalent vertices,
                L = (D-3) n + 4 legs; nonzero only if L - 1 >= H (sources span C^H).
  cone(D, n)    Q_n with one insertion O: Gamma - O a forest of p trees, deg O = n + p, others deg <= 3;
                sources: bivalent vertices and O.  (6d hCS on C^{3|2}, 8d SYM, MSYM (H,1).)
  twistor(n)    D = 3, Q_n of twistor BF: Gamma - O connected, cycle rank 0 or 1, others deg <= 3.
  fan(n)        D = 2, Q_n of 4d N=1 with cubic W: Gamma - O a path, deg O = n + 1.
Counts are orbits of (Gamma, O) under Aut(Gamma).  Classes: SP / B (polynomial: S, P, source-free
stars) / R_k (k sourced stars, rational transfer) / I (irreducible).  Output atlas.log, atlas.json.
"""
import json
import sys
import time

import networkx as nx

from graphs_enum import admissible
from reduce_search import Network, search, tight_search
from htk import catalog


def classify(G, S0, D):
    net = Network(G, S0)
    rp = search(net, ("S", "Y0"))
    if rp["reducible"]:
        k = 0
    else:
        rr = search(net, ("S", "Y0", "Y1"))
        k = rr["min_Y1"]
    sp = bool(catalog.is_sp(G))
    ta = tight_search(G, S0, D)["reducible"]
    return {"sp": sp, "k": k, "tight": ta}


def rooted_orbits(G, roots):
    GM = nx.algorithms.isomorphism.GraphMatcher(G, G)
    autos = list(GM.isomorphisms_iter())
    seen, reps = set(), []
    for O in roots:
        if O in seen:
            continue
        orb = {g[O] for g in autos}
        seen |= orb
        reps.append(O)
    return reps


def set_bulk3(D, n):
    out = []
    for G in admissible(D, n):
        if max(d for _, d in G.degree()) <= 3:
            biv = [v for v in G if G.degree(v) == 2]
            out.append((G, None, biv))
    return out


def set_cone(D, n):
    out = []
    for G in admissible(D, n):
        cand = []
        for O in G:
            H = G.copy()
            H.remove_node(O)
            if not nx.is_forest(H):
                continue
            p = nx.number_connected_components(H)
            if G.degree(O) != n + p or any(G.degree(v) > 3 for v in H):
                continue
            cand.append(O)
        for O in rooted_orbits(G, cand):
            S0 = sorted(set(v for v in G if G.degree(v) == 2) | {O})
            out.append((G, O, S0))
    return out


def set_twistor(n):
    out = {"tree": [], "unicyclic": []}
    for G in admissible(3, n):
        cand = {"tree": [], "unicyclic": []}
        for O in G:
            H = G.copy()
            H.remove_node(O)
            if not nx.is_connected(H) or any(G.degree(v) > 3 for v in H):
                continue
            r = H.number_of_edges() - H.number_of_nodes() + 1
            if r == 0:
                cand["tree"].append(O)
            elif r == 1:
                cand["unicyclic"].append(O)
        for key in cand:
            for O in rooted_orbits(G, cand[key]):
                S0 = sorted(set(v for v in G if G.degree(v) == 2) | {O})
                out[key].append((G, O, S0))
    return out


def set_fan(n):
    out = []
    for G in admissible(2, n):
        cand = []
        for O in G:
            H = G.copy()
            H.remove_node(O)
            if nx.is_connected(H) and H.number_of_edges() == H.number_of_nodes() - 1 and \
                    max((d for _, d in H.degree()), default=0) <= 2 and G.degree(O) == n + 1 and \
                    all(G.degree(v) <= 3 for v in H):
                cand.append(O)
        for O in rooted_orbits(G, cand):
            S0 = sorted(set(v for v in G if G.degree(v) == 2) | {O})
            out.append((G, O, S0))
    return out


def summarize(items, D):
    res = {"count": len(items), "SP": 0, "B_nonSP": 0, "tight": 0, "k": {}, "nonSP_names": []}
    for G, O, S0 in items:
        c = classify(G, S0, D)
        res["SP"] += c["sp"]
        res["tight"] += c["tight"]
        if c["k"] == 0 and not c["sp"]:
            res["B_nonSP"] += 1
        key = str(c["k"])
        res["k"][key] = res["k"].get(key, 0) + 1
        if not c["sp"]:
            res["nonSP_names"].append("%s k=%s" % (catalog.describe(G).split(", degrees")[0], c["k"]))
    return res


if __name__ == "__main__":
    log = open("atlas.log", "w")
    out = {}

    def P(*a):
        s = " ".join(str(x) for x in a)
        print(s)
        log.write(s + "\n")
        log.flush()

    t0 = time.time()
    P("== bulk3: cubic bulk graphs (legs at bivalent vertices), pattern = legs ==")
    for D, ns in [(2, range(1, 7)), (3, range(1, 5)), (4, range(1, 4)), (5, range(1, 4)), (6, range(1, 4))]:
        for n in ns:
            items = set_bulk3(D, n)
            s = summarize(items, D)
            L = (D - 3) * n + 4
            s["legs"] = L
            out["bulk3 D=%d n=%d" % (D, n)] = s
            P("bulk3 D=%d n=%d legs %d: graphs %d | SP %d | tight-poly %d | B nonSP %d | k_Y1 %s" % (
                D, n, L, s["count"], s["SP"], s["tight"], s["B_nonSP"], s["k"]))
            if s["nonSP_names"] and n <= 4:
                from collections import Counter
                for nm, c in sorted(Counter(s["nonSP_names"]).items()):
                    P("      %dx %s" % (c, nm))
    P("== cone: Q_n with Gamma - O a forest, deg O = n + p (pattern bivalent + O) ==")
    for D, ns in [(3, range(1, 5)), (4, range(1, 4)), (5, range(1, 4)), (6, range(1, 4)), (2, range(1, 7))]:
        for n in ns:
            items = set_cone(D, n)
            s = summarize(items, D)
            out["cone D=%d n=%d" % (D, n)] = s
            P("cone D=%d n=%d: (Gamma,O) %d | SP %d | tight-poly %d | k_Y1 %s" % (D, n, s["count"], s["SP"], s["tight"], s["k"]))
    P("== fan: D=2, Gamma - O a path ==")
    for n in range(1, 7):
        items = set_fan(n)
        s = summarize(items, 2)
        out["fan n=%d" % n] = s
        P("fan n=%d: %d | SP %d | k_Y1 %s" % (n, s["count"], s["SP"], s["k"]))
    P("== twistor BF Q_n: D=3, Gamma - O connected, cycle rank 0 / 1 ==")
    for n in range(1, 5):
        sets = set_twistor(n)
        for key in ("tree", "unicyclic"):
            s = summarize(sets[key], 3)
            out["twistor %s n=%d" % (key, n)] = s
            P("twistor %s n=%d: (Gamma,O) %d | SP %d | tight-poly %d | B nonSP %d | k_Y1 %s" % (
                key, n, s["count"], s["SP"], s["tight"], s["B_nonSP"], s["k"]))
    json.dump(out, open("atlas.json", "w"), indent=1)
    P("time %.1fs" % (time.time() - t0))
