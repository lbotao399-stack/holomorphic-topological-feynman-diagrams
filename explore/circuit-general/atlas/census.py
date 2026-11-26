"""Task 2: census of admissible graphs by circuit reducibility.

Patterns: gen = bivalent vertices sourced, all others source-free (sum of sources zero);
          all = every vertex sourced.
Classes (pattern gen):
  SP     series-parallel (no K4 minor): S moves only
  A      tight polynomial reduction to an edge (tight chain moves; D = 2: LLC (b), (c))
  B      polynomial reduction to an edge with S, P, Y0 (source-free stars, any number of merges)
  R_k    reduction needs k >= 1 sourced stars Y1 (rational transfer), k minimal
  I      no reduction even with Y1 (stuck at minimum degree >= 4)
Cores: B-core = minimal stuck state under S, Y0 (graph with source flags);
       A-core = minimal stuck tight state under tight polynomial moves.
Output: census.log, census.json.
"""
import json
import sys
import time

import networkx as nx

from graphs_enum import admissible
from reduce_search import Network, search, core_graph, tight_search, tight_core_graph, Classes, describe_core
from htk import catalog

RANGES = [(2, n) for n in range(1, 7)] + [(3, n) for n in range(1, 5)] + [(4, n) for n in range(1, 4)] + \
         [(5, n) for n in range(1, 4)] + [(6, n) for n in range(1, 4)]


def classify(G, D):
    rec = {}
    nv = G.number_of_nodes()
    biv = sorted(v for v in G if G.degree(v) == 2)
    rec["nv"], rec["ne"] = nv, G.number_of_edges()
    rec["sp"] = bool(catalog.is_sp(G))
    rec["n_biv"] = len(biv)
    rec["desc"] = catalog.describe(G)
    rec["edges"] = sorted(tuple(sorted(e)) for e in G.edges())
    out_cores = {}
    # generic pattern
    net = Network(G, biv)
    rp = search(net, ("S", "Y0"))
    rec["B"] = rp["reducible"]
    if rp["reducible"]:
        rec["k_gen"] = 0
        rec["seq_gen"] = rp["sequence"]
    else:
        out_cores["B"] = core_graph(net, rp["stuck"][0])
        rr = search(net, ("S", "Y0", "Y1"))
        rec["k_gen"] = rr["min_Y1"]
        rec["seq_gen"] = rr["sequence"]
        if not rr["reducible"]:
            out_cores["I"] = core_graph(net, rr["stuck"][0])
    # tight polynomial (definition A)
    ta = tight_search(G, biv, D)
    rec["A"] = ta["reducible"]
    if not ta["reducible"]:
        out_cores["A"] = tight_core_graph(ta["stuck"][0])
    if D == 2:
        ta2 = tight_search(G, biv, D, sourced_star=True)
        rec["A_rat"] = ta2["n_sourced"]
    # gen + one insertion vertex v* of degree >= 3
    ins = []
    for v in sorted(G):
        if G.degree(v) < 3:
            continue
        netv = Network(G, biv + [v])
        r1 = search(netv, ("S", "Y0"))
        if r1["reducible"]:
            ins.append(0)
        else:
            r2 = search(netv, ("S", "Y0", "Y1"))
            ins.append(r2["min_Y1"])
    rec["k_ins"] = ins
    # all vertices sourced
    neta = Network(G, list(G.nodes()))
    ra = search(neta, ("S", "Y0", "Y1"))
    rec["k_all"] = ra["min_Y1"]
    if not ra["reducible"]:
        out_cores["I_all"] = core_graph(neta, ra["stuck"][0])
    return rec, out_cores


def main(ranges, out_json="census.json", out_log="census.log"):
    log = open(out_log, "w")

    def P(*a):
        s = " ".join(str(x) for x in a)
        print(s)
        log.write(s + "\n")
        log.flush()

    P("census: patterns gen (bivalent sourced) and all; classes SP, A, B, R_k, I (see census.py)")
    data = {}
    table = []
    for D, n in ranges:
        t0 = time.time()
        gs = admissible(D, n)
        cls = {"B": Classes(), "A": Classes(), "I": Classes(), "I_all": Classes()}
        recs = []
        for i, G in enumerate(gs):
            rec, cores = classify(G, D)
            rec["index"] = i
            for key, H in cores.items():
                rec["core_" + key] = cls[key].add(H, i)
            recs.append(rec)
        n_sp = sum(r["sp"] for r in recs)
        n_A = sum(r["A"] for r in recs)
        n_B = sum(r["B"] for r in recs)
        kg = {}
        for r in recs:
            kg[r["k_gen"]] = kg.get(r["k_gen"], 0) + 1
        ka = {}
        for r in recs:
            ka[r["k_all"]] = ka.get(r["k_all"], 0) + 1
        few = sum(1 for r in recs if r["n_biv"] < 2)
        ins_all = sum(1 for r in recs if r["k_ins"] and all(k == 0 for k in r["k_ins"]))
        ins_some = sum(1 for r in recs if r["k_ins"] and any(k == 0 for k in r["k_ins"]))
        ins_none = sum(1 for r in recs if r["k_ins"] and not any(k == 0 for k in r["k_ins"]))
        row = {"D": D, "n": n, "graphs": len(recs), "SP": n_SP if False else n_sp, "A": n_A, "B": n_B,
               "B_nonSP": n_B - n_sp, "gen_k": {str(k): v for k, v in sorted(kg.items(), key=lambda x: (x[0] is None, x[0] or 0))},
               "all_k": {str(k): v for k, v in sorted(ka.items(), key=lambda x: (x[0] is None, x[0] or 0))},
               "fewer_than_2_sources_gen": few,
               "ins_poly_all": ins_all, "ins_poly_some": ins_some, "ins_poly_none": ins_none}
        if D == 2:
            ar = {}
            for r in recs:
                ar[r["A_rat"]] = ar.get(r["A_rat"], 0) + 1
            row["A_rat"] = {str(k): v for k, v in ar.items()}
        table.append(row)
        P("D=%d n=%d: graphs %d | SP %d | A %d | B %d (nonSP %d) | gen: k_Y1 %s | all: k_Y1 %s | gen pattern with <2 sources: %d | gen+1 insertion (deg>=3): poly for all/some/no choice %d/%d/%d%s | %.1fs" % (
            D, n, len(recs), n_sp, n_A, n_B, n_B - n_sp, row["gen_k"], row["all_k"], few, ins_all, ins_some, ins_none,
            (" | D=2 tight with sourced stars (LLC App E): n_sourced %s" % row["A_rat"]) if D == 2 else "",
            time.time() - t0))
        for key, title in (("A", "A-cores (tight polynomial reduction stuck), gen"),
                           ("B", "B-cores (S, Y0 stuck), gen"),
                           ("I", "irreducible with Y1, gen"),
                           ("I_all", "irreducible with Y1, all sourced")):
            reps = cls[key].reps
            if not reps:
                continue
            P("   %s: %d classes, %d graphs" % (title, len(reps), sum(r[2] for r in reps)))
            for k, H, cnt, members in sorted(reps, key=lambda r: (-r[2], r[1].number_of_nodes())):
                P("      x%-4d %s  members %s%s" % (cnt, describe_core(H), members[:12], " ..." if len(members) > 12 else ""))
        data["%d,%d" % (D, n)] = {"summary": row, "graphs": recs,
                                  "classes": {key: [{"desc": describe_core(H), "count": cnt, "members": members,
                                                     "edges": sorted(tuple(sorted(e)) for e in H.edges()),
                                                     "src": sorted(v for v in H if H.nodes[v].get("src"))}
                                                    for k, H, cnt, members in cls[key].reps]
                                              for key in cls}}
    json.dump({"table": table, "data": data}, open(out_json, "w"), default=str)
    P("done")


if __name__ == "__main__":
    rng = RANGES
    if len(sys.argv) > 1:
        rng = [tuple(map(int, a.split(","))) for a in sys.argv[1:]]
    main(rng)
