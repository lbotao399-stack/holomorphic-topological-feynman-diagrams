"""g^1, nu^0 sector of {I^12}_0 (11d, S3): flows, vertex tensors, Grassmann signs, integrand.

Conventions (S3):
  W = beta d_i mu^i + f(nu) mu^i mu^j d_i gamma_j + (g/6) eps^{ijklm} gamma_i d_j gamma_k d_l gamma_m   (R-action)
  Delta_vw = dL/dgamma_{i,v} dL/dmu^i_w + dL/dmu^i_v dL/dgamma_{i,w} + (beta,nu)                       (eq-Delta)
  W_Gamma = Delta_{e_1} ... Delta_{e_13} applied to W_0 W_1 ... W_11 (vertex order), edge order of htk theta
  (= S3 order A_1..A_p, B_1..B_q, C_1..C_r, kernel orientation N -> N').
  d_a on the contracted letter at e(0) inserts -y_e[a], at e(1) inserts +y_e[a]   (S3 Section 4.3).
  External letter at v: plane wave, mu -> u_v[a], gamma -> vg[a], d_a on it -> lam_v[a].
Sector: one (c) vertex, eleven (b_0) vertices, external mu^9 gamma (structure.py).
"""
import sys
from fractions import Fraction as Fr
from itertools import product, permutations
import numpy as np
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
from htk.graph import theta
from htk import brute
from htk.linalg import perm_sign

PQR = [(1, 6, 6), (2, 5, 6), (3, 4, 6), (3, 5, 5), (4, 4, 5)]
AUT = {(1, 6, 6): 4, (2, 5, 6): 2, (3, 4, 6): 2, (3, 5, 5): 4, (4, 4, 5): 4}
H = 5

# letters: (field, index variable, derivative variables)
TERMS = {
    'b': (('a', 'b'), [('mu', 'a', ()), ('mu', 'b', ()), ('ga', 'b', ('a',))]),
    'c': (('i', 'j', 'k', 'l', 'm'), [('ga', 'i', ()), ('ga', 'k', ('j',)), ('ga', 'm', ('l',))]),
}

def eps5(vals):
    if len(set(vals)) < 5:
        return 0
    return perm_sign(list(vals))

def coef(typ, env):
    if typ == 'b':
        return Fr(1)
    return Fr(eps5([env[s] for s in 'ijklm']), 6)

def m1_flows(G):
    """Flow orientations with exactly one (c) vertex (K = 0).  fl[e] = (tail, head) of the flow."""
    out = []
    for o in product((0, 1), repeat=G.ne):
        fl = [(a, b) if o[e] == 0 else (b, a) for e, (a, b) in enumerate(G.edges)]
        ind = [0] * G.nv; outd = [0] * G.nv
        for (t, h) in fl:
            outd[t] += 1; ind[h] += 1
        typ = {}; ok = True
        for v in range(G.nv):
            io = (ind[v], outd[v]); d = G.degree(v)
            if d == 2:
                typ[v] = {(1, 1): 'b', (0, 2): 'b', (2, 0): 'c'}[io]
            elif io == (1, 2):
                typ[v] = 'b'
            elif io == (3, 0):
                typ[v] = 'c'
            else:
                ok = False
        if ok and sum(1 for v in typ if typ[v] == 'c') == 1:
            out.append((fl, typ))
    return out

def patterns(G, fl, typ, order=None, vorder=None):
    """All letter-to-edge maps.  Returns list of (sign, maps) with maps[v] = dict letter position ->
    edge or 'ext'.  The sign is the Grassmann sign of Delta_{e_1}..Delta_{e_13} on W_0..W_11 relative
    to the external odd letters left in vertex order."""
    inc = {v: [e for e, (a, b) in enumerate(G.edges) if v in (a, b)] for v in range(G.nv)}
    per_v = []
    for v in range(G.nv):
        L = TERMS[typ[v]][1]
        opts = []
        for perm in permutations(range(len(L)), len(inc[v])):
            good = True
            for e, pos in zip(inc[v], perm):
                want = 'mu' if fl[e][0] == v else 'ga'
                if L[pos][0] != want:
                    good = False
            if good:
                mp = {pos: e for e, pos in zip(inc[v], perm)}
                for pos in range(len(L)):
                    mp.setdefault(pos, 'ext')
                opts.append(mp)
        per_v.append(opts)
    order = list(range(G.ne)) if order is None else list(order)
    vorder = list(range(G.nv)) if vorder is None else list(vorder)
    res = []
    for cidx in product(*[range(len(o)) for o in per_v]):
        choice = [per_v[v][c] for v, c in enumerate(cidx)]
        letters = [(v, pos, TERMS[typ[v]][1][pos][0]) for v in vorder for pos in range(len(TERMS[typ[v]][1]))]
        sign = 1
        for e in reversed(order):
            t, h = fl[e]
            # dL/dmu at the tail (odd), then dL/dgamma at the head (even)
            k = next(i for i, (v, pos, f) in enumerate(letters) if v == t and choice[v][pos] == e)
            assert letters[k][2] == 'mu'
            nodd = sum(1 for (_, _, f) in letters[:k] if f == 'mu')
            sign *= (-1) ** nodd
            letters.pop(k)
            k = next(i for i, (v, pos, f) in enumerate(letters) if v == h and choice[v][pos] == e)
            assert letters[k][2] == 'ga'
            letters.pop(k)
        # remaining odd letters: bring them to increasing vertex order
        rem = [v for (v, pos, f) in letters if f == 'mu']
        sign *= perm_sign(sorted(range(len(rem)), key=lambda i: rem[i])) if rem else 1
        res.append((sign, cidx))
    return res, per_v

def vertex_tensor(G, v, typ, mp, y, ext):
    """Tensor over the incident edges of v (in increasing edge order), entries Fractions."""
    inc = [e for e, (a, b) in enumerate(G.edges) if v in (a, b)]
    var, L = TERMS[typ]
    T = np.empty((H,) * len(inc), dtype=object); T.fill(Fr(0))
    for vals in product(range(H), repeat=len(var)):
        env = dict(zip(var, vals))
        c = coef(typ, env)
        if c == 0:
            continue
        idx = {}
        for pos, (f, iv, ders) in enumerate(L):
            e = mp[pos]
            if e == 'ext':
                c *= (ext['u'] if f == 'mu' else ext['vg'])[env[iv]]
                for dv in ders:
                    c *= ext['lam'][env[dv]]
            else:
                idx[e] = env[iv]
                for dv in ders:
                    c *= (-y[e][env[dv]] if G.edges[e][0] == v else y[e][env[dv]])
            if c == 0:
                break
        if c == 0:
            continue
        T[tuple(idx[e] for e in inc)] += c
    return T

def network(G, interior, tensors):
    """Contract the theta network: N = 0, N' = root, paths along `interior`."""
    paths, off = [], 0
    for ints in interior:
        n_ = len(ints) + 1
        M = None
        for v in ints:
            M = tensors[v] if M is None else M.dot(tensors[v])
        if M is None:
            M = np.empty((H, H), dtype=object); M.fill(Fr(0))
            for i in range(H):
                M[i, i] = Fr(1)
        paths.append(M)
    TN, TNp = tensors[0], tensors[G.root]
    tot = Fr(0)
    MA, MB, MC = paths
    X = np.tensordot(TN, MA, axes=([0], [0]))          # (iB, iC, jA)
    X = np.tensordot(X, MB, axes=([0], [0]))           # (iC, jA, jB)
    X = np.tensordot(X, MC, axes=([0], [0]))           # (jA, jB, jC)
    return sum((X * TNp).flatten(), Fr(0))

def flow_integrand(G, interior, fl, typ, pats, y, ext_at):
    """Sum over patterns of sign * contracted network; pats = patterns(...); ext_at[v] = dict(u, vg, lam)."""
    res, per_v = pats
    cache = [[vertex_tensor(G, v, typ[v], mp, y, ext_at[v]) for mp in per_v[v]] for v in range(G.nv)]
    tot = Fr(0)
    for sign, cidx in res:
        tot += sign * network(G, interior, [cache[v][c] for v, c in enumerate(cidx)])
    return tot

# ---------------------------------------------------------------- polynomial currents (S3 eq-series, eq-parallel)
def poly_currents(pqr, interior, x, w, lam):
    """x[P] = (x^P_1..x^P_{n_P}) summing to 1, w = (w_A,w_B,w_C); lam[v] for v interior, lam_N."""
    mh, Qs = [], []
    for P, ints in enumerate(interior):
        n_ = pqr[P]
        Q = [[Fr(0)] * H]
        for k in range(1, n_):
            Q.append([Q[-1][i] + lam[ints[k - 1]][i] for i in range(H)])
        Qs.append(Q)
        mh.append([sum(x[P][j] * Q[j][i] for j in range(n_)) for i in range(H)])
    lamN = lam[0]
    S = [lamN[i] + sum(m[i] for m in mh) for i in range(H)]
    y = []
    for P in range(3):
        l = [w[P] * S[i] - mh[P][i] for i in range(H)]
        for j in range(pqr[P]):
            y.append([l[i] + Qs[P][j][i] for i in range(H)])
    return y

def tmap(pqr, x, w):
    return [x[P][j] / w[P] for P in range(3) for j in range(pqr[P])]
