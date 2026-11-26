"""Vertex orders, ribbon orientations and trace contractions for the theta graphs."""
from fractions import Fraction as Fr

TENS = {
    'f': [(1, [[0, 1, 2]]), (-1, [[0, 2, 1]])],
    'd': [(1, [[0, 1, 2]]), (1, [[0, 2, 1]])],
    's1': [(1, [[0], [1, 2]]), (1, [[1], [0, 2]]), (1, [[2], [0, 1]])],
    's0': [(1, [[0], [1], [2]])],
}


def perm_sign(seq):
    seq = list(seq)
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def canon_cycle(c):
    c = tuple(c)
    i = c.index(min(c))
    return c[i:] + c[:i]


def canon_tau(cycles):
    return tuple(sorted(canon_cycle(c) for c in cycles))


def local_orders(nv, edges, ext):
    local = []
    for v in range(nv):
        hs = [('x', v)] if v in ext else []
        hs += sorted([('e', e, s) for e, ab in enumerate(edges) for s in (0, 1) if ab[s] == v],
                     key=lambda h: (h[1], h[2]))
        local.append(hs)
    return local


def orientation_sign(nv, edges, ext):
    """sgn(L -> L'), L = (x_v1..x_v4, (e0,0),(e0,1),...),  L' = concatenated local orders."""
    L = [('x', v) for v in ext]
    for e in range(len(edges)):
        L += [('e', e, 0), ('e', e, 1)]
    Lp = [h for hs in local_orders(nv, edges, ext) for h in hs]
    pos = {h: i for i, h in enumerate(L)}
    return perm_sign([pos[h] for h in Lp])


def colour_weight(nv, edges, ext, vtens=None):
    """dict tau (labels = positions in ext) -> {Npow: int}; vtens: {v: tensor name or list}."""
    vtens = vtens or {}
    local = local_orders(nv, edges, ext)
    sgn = orientation_sign(nv, edges, ext)
    terms = [(sgn, 0, [])]
    for v in range(nv):
        T = vtens.get(v, 'f')
        T = TENS[T] if isinstance(T, str) else T
        hs = local[v]
        new = []
        for (c, npw, ws) in terms:
            for coef, words in T:
                new.append((c * coef, npw, ws + [[hs[i] for i in w] for w in words]))
        terms = new
    for e in range(len(edges)):
        h, hp = ('e', e, 0), ('e', e, 1)
        new = []
        for (c, npw, ws) in terms:
            i = next(k for k, w in enumerate(ws) if h in w)
            j = next(k for k, w in enumerate(ws) if hp in w)
            ws = [list(w) for w in ws]
            if i != j:
                w1, w2 = ws[i], ws[j]
                a = w1.index(h)
                w1 = w1[a + 1:] + w1[:a]
                b = w2.index(hp)
                w2 = w2[b + 1:] + w2[:b]
                merged = w1 + w2
                rest = [w for k, w in enumerate(ws) if k not in (i, j)]
                if merged:
                    rest.append(merged)
                    new.append((c, npw, rest))
                else:
                    new.append((c, npw + 1, rest))
            else:
                w = ws[i]
                a = w.index(h)
                w = w[a:] + w[:a]
                b = w.index(hp)
                B = w[1:b]
                A = w[b + 1:]
                rest = [x for k, x in enumerate(ws) if k != i]
                add = 0
                for part in (A, B):
                    if part:
                        rest.append(part)
                    else:
                        add += 1
                new.append((c, npw + add, rest))
        terms = new
    out = {}
    for (c, npw, ws) in terms:
        cyc = [tuple(ext.index(h[1]) for h in w) for w in ws]
        tau = canon_tau(cyc)
        dd = out.setdefault(tau, {})
        dd[npw] = dd.get(npw, 0) + c
    return {t: {p: Fr(v) for p, v in dd.items() if v} for t, dd in out.items() if any(dd.values())}
