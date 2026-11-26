"""Chain elimination + sourced star-triangle (Y-Delta) currents for the K4 subdivisions (my derivation):
  K4 sources  mu_p += sum_k s_{E,k} lam_{E,k},  mu_q += sum_k (1 - s_{E,k}) lam_{E,k}   (chain E = p..q)
  chain currents (direction p->q)  y_{E,j} = Y_E - mhat_E + Q_{E,j},  mhat_E = sum_k s_{E,k} lam_{E,k},  Q_{E,j} = sum_{k<j} lam_{E,k}
  Y-Delta at Dc:  p_v = (1-r_v) g_v,  Q = pA+pB+pC,  mu'_v = mu_v + g_v mu_Dc,
     merged triangle  Y_{XY} = (p_Y mu'_X - p_X mu'_Y)/Q,
     original triangle edge X->Y: r_Z Y_{XY};  star edge Dc->v: g_v mu_Dc - sum_{w != v} (1-r_{Z(vw)}) Y_{vw}.
Checked here against the brute-force Kirchhoff currents (htk.brute.currents) at random rational xi."""
import sys, random
from fractions import Fraction as Fr
import sympy as sp
from k4exact import K4Case

def formula_currents(c):
    A, B, C = c.tri
    g = {A: c.gA, B: c.gB, C: 1 - c.gA - c.gB}
    r = {A: c.rA, B: c.rB, C: c.rC}
    lamv = {v: sp.Matrix(c.lam.get(v, [0, 0, 0])) for v in range(c.nv)}
    mu = {v: sp.zeros(3, 1) for v in c.cubic}
    mhat, Qs = {}, {}
    for ci, (path, es) in enumerate(c.chains):
        p, q = path[0], path[-1]
        a = len(es)
        mh = sp.zeros(3, 1)
        for k in range(1, a):
            s = c.sig[(ci, k)]
            mu[p] += s * lamv[path[k]]
            mu[q] += (1 - s) * lamv[path[k]]
            mh += s * lamv[path[k]]
        mhat[ci] = mh
        Qs[ci] = [sum((lamv[path[k]] for k in range(1, j)), sp.zeros(3, 1)) for j in range(1, a + 1)]
    pv = {v: (1 - r[v]) * g[v] for v in c.tri}
    Q = pv[A] + pv[B] + pv[C]
    mup = {v: mu[v] + g[v] * mu[c.Dc] for v in c.tri}
    def Ytri(X, Y):
        return (pv[Y] * mup[X] - pv[X] * mup[Y]) / Q
    def third(X, Y):
        return [u for u in c.tri if u not in (X, Y)][0]
    def YK4(p, q):
        if p == c.Dc:
            v = q
            return g[v] * mu[c.Dc] - sum(((1 - r[third(v, w)]) * Ytri(v, w) for w in c.tri if w != v), sp.zeros(3, 1))
        if q == c.Dc:
            return -YK4(q, p)
        return r[third(p, q)] * Ytri(p, q)
    y = [None] * 10
    for ci, (path, es) in enumerate(c.chains):
        Y = YK4(path[0], path[-1])
        for j, e in enumerate(es, start=1):
            val = Y - mhat[ci] + Qs[ci][j - 1]
            sgn = 1 if c.edges[e] == (path[j - 1], path[j]) else -1
            y[e] = sgn * val
    return y, Q

def check_currents(c, npts=3, seed=1):
    y, Q = formula_currents(c)
    f = sp.lambdify(c.xi, [list(v) for v in y], 'sympy')
    rng = random.Random(seed)
    ok = True
    for _ in range(npts):
        pt = c.random_point(rng)
        P = [sp.Rational(p.numerator, p.denominator) for p in pt]
        t = [Fr(int(sp.numer(x)), int(sp.denom(x))) for x in c.t_f(*P)]
        yb = c.currents_brute(t)
        yf = f(*P)
        for e in range(10):
            for i in range(3):
                v = sp.Rational(yf[e][i])
                if Fr(int(sp.numer(v)), int(sp.denom(v))) != yb[e][i]:
                    ok = False
    return ok

if __name__ == "__main__":
    for w in range(6):
        for dci in range(4):
            c = K4Case(w, Dc_index=dci)
            print(c.name, "Dc", c.Dc, "currents formula == brute Kirchhoff:", check_currents(c))
