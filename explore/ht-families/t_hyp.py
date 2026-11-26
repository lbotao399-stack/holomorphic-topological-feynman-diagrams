"""Test: is each insertion I_m(J, Theta) a polynomial in hat J = J + (dt_m / 2 t_m) Theta ?"""
import sys, random
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import htk
from htk import spred
import enumerate_graphs as eg
import hf_recursion as hr
from hf_brute import wedge, fadd

def split(I):
    c, d = {}, {}
    for (odd, (p,)), v in I.items():
        if odd and odd[-1] == hr.THETA:
            d.setdefault(p, {})[odd[:-1]] = v
        else:
            c.setdefault(p, {})[odd] = v
    return c, d

gs = [Gx for Gx in eg.enumerate_tight(2, 4) if eg.is_sp(Gx)]
rng = random.Random(5)
for k, Gx in enumerate(gs):
    G, _ = htk.from_networkx(Gx, None)
    for prefer in ("large", "small"):
        red = spred.sp_reduce(G, prefer=prefer)
        plan = [(mv.p, mv.q, mv.chains, mv.signs) for mv in red.moves]
        rec = hr.Recursion(G.nv, G.edges, plan)
        t = [Fr(rng.randint(1, 30), rng.randint(1, 12)) for _ in G.edges]
        rho, U = rec.run(t)
        out = []
        for m in range(len(plan)):
            e = G.ne + m
            c, d = split(rec.ins[e])
            eps = {key: v / (2 * rec.tt[e]) for key, v in rec.dt[e].items()}
            good = True
            for p in set(c) | set(d):
                lhs = d.get(p, {})
                rhs = {kk: (p + 1) * v for kk, v in wedge(c.get(p + 1, {}), eps).items()}
                if fadd(lhs, rhs, -1):
                    good = False
            out.append((m + 1, sorted(c), good))
        print("graph", k, prefer, "moves (index, J-degrees, I = f(hatJ)):", out)
