"""Direct numerics for the C^2 prism / K33 coefficient (no current formulas):
t(f, R): undo the two source-free star-triangle moves for resistances only, top level = K4 resistances R
(log-lattice, R_last = 1), fractions f in (0,1) (Gauss-Legendre, exact for polynomials in f; checked by
comparing m and m+1).  Density D(t) from htk.oracle (brute force) times |det[t; dt/d(f, log R)]|."""
import sys, os, time, itertools
import numpy as np
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code'))
from htk import oracle
from c2case import prism_case, k33_case, fe

def tmap(c):
    K4v = [c.Dc] + c.tri
    K4e = [fe(a, b) for a, b in itertools.combinations(sorted(K4v), 2)]
    Rs = {e: sp.Symbol('R%d%d' % tuple(sorted(e))) for e in K4e}
    net = {e: Rs[e] for e in K4e}
    for li in reversed(range(len(c.levels))):
        v, us, merged = c.levels[li]
        tau = {}
        new = dict(net)
        for p, m in merged.items():
            a, b = sorted(p)
            T = net[p]; del new[p]
            if m:
                f = c.fr[(li, a, b)]
                new[p] = T / f; tau[p] = T / (1 - f)
            else:
                tau[p] = T
        S = sum(tau.values())
        for u in us:
            o = [w for w in us if w != u]
            new[fe(v, u)] = tau[fe(u, o[0])] * tau[fe(u, o[1])] / S
        net = new
    t = [net[fe(a, b)] for a, b in c.edges]
    fs = [c.fr[k] for k in sorted(c.fr)]
    Rl = [Rs[e] for e in K4e]
    L = sp.symbols('L0:6')
    sub = {Rl[i]: sp.exp(L[i]) for i in range(6)}
    tL = [ti.subs(sub) for ti in t]
    coords = fs + list(L[:5])
    Jm = [tL] + [[sp.diff(ti, x) for ti in tL] for x in coords]
    fJ = sp.lambdify(fs + list(L), Jm, 'numpy')
    return fJ, len(fs)

def run(c, SA, SB, h, tol, m=4, shift=0.0, chunk=3000):
    fJ, nf = tmap(c)
    st = oracle._Setup(c.G)
    lam_all = np.array([[float(x) for x in v] for v in c.lam_num(SA, SB)])
    lam = lam_all[:5]
    s, w = np.polynomial.legendre.leggauss(m); s = (s + 1) / 2; w = w / 2
    fn = np.array(list(itertools.product(s, repeat=nf))); fw = np.prod(np.array(list(itertools.product(w, repeat=nf))), axis=1)
    def F(keys):
        Lv = np.array(keys, dtype=float) * h + shift
        B = len(keys); Q = len(fw)
        LL = np.repeat(Lv, Q, axis=0); FF = np.tile(fn, (B, 1))
        args = [FF[:, i] for i in range(nf)] + [LL[:, i] for i in range(5)] + [np.zeros(B * Q)]
        rows = fJ(*args)
        M = np.zeros((B * Q, 9, 9))
        for i in range(9):
            for j in range(9):
                M[:, i, j] = np.broadcast_to(np.asarray(rows[i][j], dtype=float), (B * Q,))
        t = M[:, 0, :]
        Sd = oracle.scaled_density_batch(c.G, t, lam, 0, _setup=st)        # D * prod t
        val = Sd / np.prod(t, axis=1) * np.abs(np.linalg.det(M))
        val = np.where(np.isfinite(val), val, 0.0)
        return (val.reshape(B, Q) @ fw)
    seen = {}; front = [tuple([0] * 5)]; seen[front[0]] = None
    total = 0.0; vmax = 0.0; below = 0.0; n = 0
    while front:
        nxt = []
        for i in range(0, len(front), chunk):
            blk = front[i:i + chunk]
            vals = F(blk); n += len(blk)
            for k, v in zip(blk, vals):
                total += v; av = abs(v); vmax = max(vmax, av)
                if av > tol * vmax:
                    for d in range(5):
                        for sg in (-1, 1):
                            kk = list(k); kk[d] += sg; kk = tuple(kk)
                            if kk not in seen:
                                seen[kk] = None; nxt.append(kk)
                else:
                    below += av
        front = nxt
    return total * h ** 5, below * h ** 5, n

if __name__ == "__main__":
    from fractions import Fraction as Fr
    which = sys.argv[1]; h = float(sys.argv[2]); tol = float(sys.argv[3])
    SA, SB = Fr(sys.argv[4]), Fr(sys.argv[5])
    m = int(sys.argv[6]) if len(sys.argv) > 6 else 4
    from c2case import C2Case
    from graphs import K33
    if which == 'prism':
        c = prism_case(0)
    elif which == 'K33':
        c = k33_case(0)
    else:   # 'K33x': sources at 0,1,2,3 (sourceless 4 in {0,2,4}, 5 in {1,3,5}); resistance map only:
            # source-free move at 4 (nbrs 1,3,5; no merge), then move at 0 (nbrs 1,3,5; three merges), K4 {1,2,3,5}
        c = C2Case.__new__(C2Case)
        import sympy as sp
        from htk.graph import Graph
        c.name, c.edges, c.nv, c.G = 'K33x', K33, 6, Graph(6, K33)
        c.SA, c.SB = sp.symbols('SA SB')
        lam = {0: sp.Matrix([1, 0]), 1: sp.Matrix([0, 1]), 2: sp.Matrix([c.SA, c.SB])}
        lam[3] = -(lam[0] + lam[1] + lam[2])
        c.lamsym = {v: lam.get(v, sp.zeros(2, 1)) for v in range(6)}
        c.Dc, c.tri = 1, [2, 3, 5]
        c.levels = [(4, (1, 3, 5), {fe(1, 3): False, fe(1, 5): False, fe(3, 5): False}),
                    (0, (1, 3, 5), {fe(1, 3): True, fe(1, 5): True, fe(3, 5): True})]
        c.fr = {(1, 1, 3): sp.Symbol('f1_13'), (1, 1, 5): sp.Symbol('f1_15'), (1, 3, 5): sp.Symbol('f1_35')}
    t0 = time.time()
    v0, b0, n0 = run(c, SA, SB, h, tol, m=m)
    v1, b1, n1 = run(c, SA, SB, h, tol, m=m, shift=h / 2)
    print("%s SA=%s SB=%s h=%g tol=%g m=%d: K(shift 0) = %.10e  K(shift h/2) = %.10e  points %d %d  sub-threshold %.2e  time %.0fs"
          % (which, SA, SB, h, tol, m, v0, v1, n0, n1, b0, time.time() - t0))
