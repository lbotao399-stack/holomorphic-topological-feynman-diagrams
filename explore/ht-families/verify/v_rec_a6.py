"""a6 exceptions: nested (2,1) graphs, small plan: merges whose chains are a plain chain and a single
reduced edge; is the consumed single edge's own insertion polynomial in its own Jhat?"""
import sys, random
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/ht-families/verify")
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import htk
from htk import spred
from v_rec import run_recursion, is_poly_in_jhat
from v_b2 import nested
rng = random.Random(9)
for inner in [(1, 3, 3), (2, 2, 3)]:
    nv, edges = nested(inner)
    red = spred.sp_reduce(htk.Graph(nv, edges), D=3, prefer="small")
    plan = [(mv.p, mv.q, mv.chains, mv.signs) for mv in red.moves]
    t = [Fr(rng.randint(1, 30), rng.randint(1, 30)) for _ in edges]
    rho, U, created = run_recursion(len(edges), plan, t)
    print("nested%s small plan:" % (inner,))
    for c in created:
        poly = is_poly_in_jhat(c["ins"], c["tm"], c["dtm"], len(edges))
        degs = sorted(set(p for (m, (p,)) in c["ins"]))
        print("   move %d (new edge %d): chains %s ; J-degrees %s ; polynomial in own Jhat: %s" % (
            c["move"], len(edges) + c["move"], c["chains"], degs, poly))
