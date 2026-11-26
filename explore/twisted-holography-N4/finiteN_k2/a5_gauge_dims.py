"""Dimensions of the twisted gauge-theory cochain spaces Sym(g (x) letters)^G per sector and letter number Y,
g = u(2) or su(2) (Molien-Weyl with the dimension plethystic exponential; fermionic letters exterior).
usage: python3 a5_gauge_dims.py 't1,t2,t3,n+,n-'"""
import sys
from math import comb
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw
lam = tuple(int(x) for x in sys.argv[1].split(','))
L = 2 * sum(lam[:3]) + 3 * sum(lam[3:])
alg, _ = bw.n4_letters(L)
letters = [(tuple(alg.w[i]), alg.par[i]) for i in range(alg.n) if all(a <= b for a, b in zip(alg.w[i], lam))]
def le(a, b): return all(x <= y for x, y in zip(a, b))
def run(adj):
    # series: key (weight tuple, Y, zpow) -> count
    S = {((0,) * 5, 0, 0): 1}
    for (w, par) in letters:
        for zp, mult in adj.items():
            for _ in range(mult):
                new = dict(S)
                if par:
                    for (om, Y, z), c in S.items():
                        no = tuple(a + b for a, b in zip(om, w))
                        if le(no, lam):
                            key = (no, Y + 1, z + zp); new[key] = new.get(key, 0) + c
                else:
                    for (om, Y, z), c in S.items():
                        n = 1
                        while True:
                            no = tuple(a + n * b for a, b in zip(om, w))
                            if not le(no, lam): break
                            key = (no, Y + n, z + n * zp); new[key] = new.get(key, 0) + c
                            n += 1
                S = new
    return S
out = {}
for name, adj, haar, norm in [("U(2)", {0: 2, 1: 1, -1: 1}, {0: 2, 1: -1, -1: -1}, 2),
                              ("SU(2)", {0: 1, 1: 1, -1: 1}, {0: 2, 1: -1, -1: -1}, 2)]:
    S = run(adj)
    dims = {}
    for (om, Y, z), c in S.items():
        if om != lam: continue
        for hz, hc in haar.items():
            if z + hz == 0:
                dims[Y] = dims.get(Y, 0) + c * hc
    dims = {Y: v // norm for Y, v in sorted(dims.items()) if v}
    print(f"{lam} {name}: dim cochains per Y = {dims}", flush=True)
