import sys
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import htk, sympy as sp
from htk import spred
import enumerate_graphs as eg
gs = [Gx for Gx in eg.enumerate_tight(2, 4) if eg.is_sp(Gx)]
for k in (3, 1, 4):
    G, _ = htk.from_networkx(gs[k], None)
    red = spred.sp_reduce(G, prefer="small")
    print(red.describe())
    rho = red._rho(red.ring, red._bpoly)
    for S, lf in rho.items():
        mono = sp.Mul(*[sp.Symbol(b) ** sp.Rational(e.numerator, e.denominator) for b, e in zip(red.bary, lf.exp) if e != 0])
        print("  d", [red.coords[j] for j in S], ":", mono, "*", sp.factor(lf.num.as_expr()))
