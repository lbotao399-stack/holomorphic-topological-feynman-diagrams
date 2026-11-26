"""Colour weights W_Gamma of cubic bulk graphs for hCS (vertex f(a,b,c) = kappa(a,[b,c]), edge kappa^{-1}),
with the orientation rule of explore/anomaly-structure (lie_eval.orientation_rule, c3):
  L  = (external letters in ext order) + ((e,0),(e,1) for e in edge order),
  L' = concatenation over vertices v = 0..nv-1 of [external of v if any] + half-edges at v sorted by (e,s),
  W  = sgn(L -> L') * contraction(f over the vertex words, kappa^{-1} over the edge pairs), output in ext order.
Two evaluators:
  ribbon(...)  gl(t) in the formal trace basis ('t Hooft double lines, kappa = tr, kappa^{-1} = sum E_ij (x) E_ji,
               f(a,b,c) = tr(abc) - tr(acb)); result: dict {(t-power, canonical multiset of cyclic words): int}.
               gl(N|M) is the same with str and t = N - M (Deligne category Rep(GL_t)).
  einsum_W(...) any matrix Lie algebra (basis, trace form), exact integers via scaled kappa^{-1}."""
import itertools
from fractions import Fraction as Fr
import numpy as np


def rule_words(nv, edges, ext):
    L = [('x', v) for v in ext]
    for e in range(len(edges)):
        L += [('e', e, 0), ('e', e, 1)]
    words = []
    for v in range(nv):
        hs = [('x', v)] if v in ext else []
        hs += sorted([('e', e, s) for e, ab in enumerate(edges) for s in (0, 1) if ab[s] == v], key=lambda h: (h[1], h[2]))
        assert len(hs) == 3, (v, hs)
        words.append(hs)
    Lp = [h for hs in words for h in hs]
    pos = {h: i for i, h in enumerate(L)}
    perm = [pos[h] for h in Lp]
    inv = sum(1 for i in range(len(perm)) for j in range(i + 1, len(perm)) if perm[i] > perm[j])
    return words, (-1 if inv % 2 else 1)


def canon_cyc(word):
    if not word:
        return ()
    rots = [tuple(word[i:] + word[:i]) for i in range(len(word))]
    return min(rots)


def ribbon(nv, edges, ext):
    """W_Gamma for gl(t): {(tpow, faces): coeff}; faces = sorted tuple of canonical cyclic words of ext slot indices."""
    words, sgn = rule_words(nv, edges, ext)
    slot = {('x', v): i for i, v in enumerate(ext)}
    partner = {}
    for e in range(len(edges)):
        partner[('e', e, 0)] = ('e', e, 1)
        partner[('e', e, 1)] = ('e', e, 0)
    out = {}
    for orient in itertools.product((0, 1), repeat=nv):
        nxt = {}
        sign = sgn
        for v in range(nv):
            w0, w1, w2 = words[v]
            cyc = (w0, w1, w2) if orient[v] == 0 else (w0, w2, w1)
            if orient[v] == 1:
                sign = -sign
            for k in range(3):
                nxt[cyc[k]] = cyc[(k + 1) % 3]
        # nodes: (h,'L'), (h,'R'); traversal: at R(h) -> L(nxt[h]); at L(h): external -> R(h) (record), internal -> R(partner)
        seen = set()
        faces = []
        empty = 0
        allnodes = [(h, side) for hs in words for h in hs for side in ('L', 'R')]
        for start in allnodes:
            if start in seen:
                continue
            word = []
            node = start
            while node not in seen:
                seen.add(node)
                h, side = node
                if side == 'R':
                    node = (nxt[h], 'L')
                else:
                    if h[0] == 'x':
                        word.append(slot[h])
                        node = (h, 'R')
                    else:
                        node = (partner[h], 'R')
            if word:
                faces.append(canon_cyc(word))
            else:
                empty += 1
        key = (empty, tuple(sorted(faces)))
        out[key] = out.get(key, 0) + sign
    return {k: v for k, v in out.items() if v != 0}


def ribbon_single(rib):
    """evaluate at x_1 = .. = x_4 = x: {(tpow, partition): coeff}, partition = sorted face sizes (p_k = tr x^k)"""
    out = {}
    for (tp, faces), c in rib.items():
        key = (tp, tuple(sorted((len(f) for f in faces), reverse=True)))
        out[key] = out.get(key, 0) + c
    return {k: v for k, v in out.items() if v != 0}


def ribbon_eval(rib, mats, t):
    """numerical value of a ribbon dict at matrices mats (list in ext order), t = N"""
    tot = Fr(0)
    for (tp, faces), c in rib.items():
        val = Fr(t) ** tp
        for f in faces:
            m = None
            for i in f:
                m = mats[i] if m is None else m @ mats[i]
            val *= int(np.trace(m))
        tot += c * val
    return tot


# ------------------------------------------------------------------ explicit matrix Lie algebras
class MatAlg:
    def __init__(self, name, basis):
        self.name = name
        self.B = [np.array(b, dtype=object) for b in basis]
        d = len(basis)
        G = [[Fr(int(np.trace(self.B[a].dot(self.B[b])))) for b in range(d)] for a in range(d)]
        import sympy as sp
        Gi = sp.Matrix(G).inv()
        den = 1
        for x in Gi:
            den = den * sp.Rational(x).q // np.gcd(den, sp.Rational(x).q)
        self.kinv_scale = int(den)
        self.kinv = np.array([[int(Gi[a, b] * den) for b in range(d)] for a in range(d)], dtype=np.int64)
        f = np.zeros((d, d, d), dtype=np.int64)
        for a in range(d):
            for b in range(d):
                c = self.B[a].dot(self.B[b]) - self.B[b].dot(self.B[a])
                for e in range(d):
                    f[e, a, b] = int(np.trace(self.B[e].dot(c)))   # f(e,a,b) = tr(e [a,b])
        self.f = f
        self.G = np.array([[int(x) for x in row] for row in G], dtype=np.int64)
        self.dim = d

    def coords(self, m):
        """coordinates of a matrix in the basis (exact, via Gram)"""
        import sympy as sp
        v = sp.Matrix([int(np.trace(np.array(m, dtype=object).dot(b))) for b in self.B])
        return [Fr(int(x.p), int(x.q)) for x in (sp.Matrix(self.G.tolist()).inv() * v)]


def gl(N):
    return MatAlg('gl%d' % N, [[[int(i == a and j == b) for j in range(N)] for i in range(N)] for a in range(N) for b in range(N)])


def sl(N):
    B = [[[int(i == a and j == b) for j in range(N)] for i in range(N)] for a in range(N) for b in range(N) if a != b]
    for k in range(N - 1):
        B.append([[int(i == j == k) - int(i == j == N - 1) for j in range(N)] for i in range(N)])
    return MatAlg('sl%d' % N, B)


def so(N):
    B = [[[int(i == a and j == b) - int(i == b and j == a) for j in range(N)] for i in range(N)] for a in range(N) for b in range(a + 1, N)]
    return MatAlg('so%d' % N, B)


def sp_(N2):
    n = N2 // 2
    J = np.zeros((N2, N2), dtype=object)
    for i in range(n):
        J[i, n + i] = 1
        J[n + i, i] = -1
    B = []
    # X = J S with S symmetric
    for a in range(N2):
        for b in range(a, N2):
            S = np.zeros((N2, N2), dtype=object)
            S[a, b] += 1
            if a != b:
                S[b, a] += 1
            B.append((J.dot(S)).tolist())
    return MatAlg('sp%d' % N2, B)


def einsum_W(alg, nv, edges, ext, xs):
    """exact W_Gamma(x_1..x_4) * kinv_scale^{|E|} as Python int; xs: coordinate vectors (ints) in ext order"""
    words, sgn = rule_words(nv, edges, ext)
    letters = {}
    alph = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    k = 0
    for hs in words:
        for h in hs:
            if h[0] == 'e':
                letters[h] = alph[k]; k += 1
    ops, terms = [], []
    for v, hs in enumerate(words):
        T = alg.f
        idx = ''
        # contract the external vector into the first slot if present
        if hs[0][0] == 'x':
            vec = np.array(xs[ext.index(v)], dtype=np.int64)
            T = np.tensordot(vec, T, axes=([0], [0]))
            idx = ''.join(letters[h] for h in hs[1:])
        else:
            idx = ''.join(letters[h] for h in hs)
        ops.append(T); terms.append(idx)
    for e in range(len(edges)):
        ops.append(alg.kinv); terms.append(letters[('e', e, 0)] + letters[('e', e, 1)])
    val = np.einsum(','.join(terms) + '->', *ops, optimize='greedy')
    return sgn * int(val)


def einsum_W_tensor(alg, nv, edges, ext, X4):
    """exact sum_{abcd} W_Gamma^{abcd} X4[a,b,c,d] * kinv_scale^{|E|}; X4 an int64 dim^4 array, slots in ext order"""
    words, sgn = rule_words(nv, edges, ext)
    alph = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    letters, k = {}, 0
    for hs in words:
        for h in hs:
            letters[h] = alph[k]; k += 1
    ops, terms = [], []
    for v, hs in enumerate(words):
        ops.append(alg.f); terms.append(''.join(letters[h] for h in hs))
    for e in range(len(edges)):
        ops.append(alg.kinv); terms.append(letters[('e', e, 0)] + letters[('e', e, 1)])
    ops.append(X4); terms.append(''.join(letters[('x', v)] for v in ext))
    return sgn * int(np.einsum(','.join(terms) + '->', *ops, optimize='greedy'))


def ribbon_alt(rib):
    """alternation over the 4 external slots of a multilinear ribbon dict: {(tpow, faces): coeff} (formal basis)"""
    out = {}
    for s in itertools.permutations(range(4)):
        sg = 1
        for i in range(4):
            for j in range(i + 1, 4):
                if s[i] > s[j]:
                    sg = -sg
        for (tp, faces), c in rib.items():
            nf = tuple(sorted(canon_cyc([s[i] for i in f]) for f in faces))
            out[(tp, nf)] = out.get((tp, nf), 0) + sg * c
    return {k: v for k, v in out.items() if v != 0}


def ribbon_sym(rib):
    out = {}
    for s in itertools.permutations(range(4)):
        for (tp, faces), c in rib.items():
            nf = tuple(sorted(canon_cyc([s[i] for i in f]) for f in faces))
            out[(tp, nf)] = out.get((tp, nf), 0) + c
    return {k: v for k, v in out.items() if v != 0}


def einsum_W_fast(alg, nv, edges, ext, xs_list):
    """W_Gamma at several external assignments: kappa^{-1} absorbed into the head-end vertex tensor, externals
    contracted first, one optimal path reused. Returns list of exact ints (scaled by kinv_scale^{|E|})."""
    words, sgn = rule_words(nv, edges, ext)
    alph = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    letters = {}
    for e in range(len(edges)):
        letters[('e', e, 0)] = alph[e]; letters[('e', e, 1)] = alph[e]
    # vertex tensors with raised head-end slots: T_v = f contracted with kinv on slots that are heads ('e',e,1)
    base = []
    for hs in words:
        T = alg.f.copy()
        for pos, h in enumerate(hs):
            if h[0] == 'e' and h[2] == 1:
                T = np.moveaxis(np.tensordot(alg.kinv, T, axes=([1], [pos])), 0, pos)
        base.append(T)
    out = []
    path = None
    for xs in xs_list:
        ops, terms = [], []
        for v, hs in enumerate(words):
            T = base[v]
            if hs[0][0] == 'x':
                T = np.tensordot(np.array(xs[ext.index(v)], dtype=np.int64), T, axes=([0], [0]))
                idx = ''.join(letters[h] for h in hs[1:])
            else:
                idx = ''.join(letters[h] for h in hs)
            ops.append(T); terms.append(idx)
        expr = ','.join(terms) + '->'
        if path is None:
            path = np.einsum_path(expr, *ops, optimize='optimal')[0]
        out.append(sgn * int(np.einsum(expr, *ops, optimize=path)))
    return out


def einsum_W_float(alg, nv, edges, ext, xs):
    """same as einsum_W in float64 (BLAS); exact when every partial sum stays below 2^53 -- callers compare with
    the int64 route on samples and require integral results"""
    words, sgn = rule_words(nv, edges, ext)
    alph = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    letters, k = {}, 0
    for hs in words:
        for h in hs:
            if h[0] == 'e':
                letters[h] = alph[k]; k += 1
    f = alg.f.astype(np.float64); ki = alg.kinv.astype(np.float64)
    ops, terms = [], []
    for v, hs in enumerate(words):
        if hs[0][0] == 'x':
            ops.append(np.tensordot(np.array(xs[ext.index(v)], dtype=np.float64), f, axes=([0], [0])))
            terms.append(''.join(letters[h] for h in hs[1:]))
        else:
            ops.append(f); terms.append(''.join(letters[h] for h in hs))
    for e in range(len(edges)):
        ops.append(ki); terms.append(letters[('e', e, 0)] + letters[('e', e, 1)])
    val = np.einsum(','.join(terms) + '->', *ops, optimize='greedy')
    return sgn * float(val)
