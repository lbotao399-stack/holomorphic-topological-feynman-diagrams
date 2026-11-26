"""Planar duality of core periods: I_q(G*; a) = I_q(G; q - a o pi), from U_{G*}(t) = (prod t) U_G(1/t).
Exact check with explore/circuit-general/cores/cores.period on self-dual K4 and W4.  usage: python3 d0_planar_duality.py"""
import itertools, os, random, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../cores"))
import cores, gutil

def check(name, E, dual, nloop, qs, per_q, seed):
    idx = {e: k for k, e in enumerate(E)}
    rng = random.Random(seed); tested = bad = 0
    for q in qs:
        cands = []
        for a in itertools.product(range(1, q + 1), repeat=len(E)):
            if sum(a) != nloop * q or not gutil.convergent(E, list(a), q): continue
            b = [0] * len(E)
            for e, k in idx.items():
                b[idx[dual(e)]] = q - a[k]
            if min(b) < 1 or not gutil.convergent(E, b, q): continue
            cands.append((a, b))
        rng.shuffle(cands)
        for a, b in cands[:per_q]:
            va, vb = cores.period(E, list(a), q), cores.period(E, b, q)
            tested += 1; bad += (va != vb)
            if va != vb: print("MISMATCH", name, q, a, va, b, vb)
        print("%s q=%d: %d dual-convergent pairs, %d tested" % (name, q, len(cands), min(per_q, len(cands))), flush=True)
    print("%s: %d pairs, %d mismatches" % (name, tested, bad), flush=True)

K4 = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
check("K4", K4, lambda e: tuple(sorted(set(range(4)) - set(e))), 3, (3, 4, 5), 12, 1)
n = 4
W4 = [(i, 4) for i in range(n)] + [tuple(sorted((i, (i + 1) % n))) for i in range(n)]
def dW(e):
    if e[1] == 4:
        i = e[0]; return tuple(sorted(((i - 1) % n, i)))
    i, j = e; i0 = i if (i + 1) % n == j else j
    return (i0, 4)
check("W4", W4, dW, 4, (3, 4), 10, 2)
