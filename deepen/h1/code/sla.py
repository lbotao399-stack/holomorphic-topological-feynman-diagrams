"""Sparse exact linear algebra on vectors given as dicts key -> int.

rank_modp(vecs, p): rank modulo p (sparse elimination).
rank_exact(vecs):   rank over Q (fraction-free sparse elimination, rows normalised by gcd).
kernel_exact(images, sub, nx): kernel of x_i -> images[i] in Y / span(sub), returned as integer vectors over range(nx).
"""
from math import gcd

P1, P2 = 2147483629, 2147483587


def _keyorder(vecs):
    keys = set()
    for v in vecs:
        keys.update(v.keys())
    order = {k: i for i, k in enumerate(sorted(keys, key=repr))}
    return order


def rank_modp(vecs, p, order=None):
    if order is None:
        order = _keyorder(vecs)
    piv = {}
    r = 0
    for v in vecs:
        w = {order[k]: c % p for k, c in v.items() if c % p}
        while w:
            c0 = min(w)
            if c0 in piv:
                row = piv[c0]
                f = w[c0]
                for k, c in row.items():
                    nv = (w.get(k, 0) - f * c) % p
                    if nv:
                        w[k] = nv
                    else:
                        w.pop(k, None)
            else:
                inv = pow(w[c0], p - 2, p)
                piv[c0] = {k: (c * inv) % p for k, c in w.items()}
                r += 1
                break
    return r


def _normalize(w):
    g = 0
    for c in w.values():
        g = gcd(g, c)
        if g == 1:
            break
    lead = w[min(w)]
    if lead < 0:
        g = -g
    if g not in (0, 1):
        for k in w:
            w[k] //= g
    return w


def rank_exact(vecs, order=None):
    if order is None:
        order = _keyorder(vecs)
    piv = {}
    r = 0
    for v in vecs:
        w = {order[k]: c for k, c in v.items() if c}
        while w:
            c0 = min(w)
            if c0 in piv:
                row = piv[c0]
                a = row[c0]
                b = w[c0]
                g = gcd(a, b)
                fa, fb = a // g, b // g
                nw = {}
                for k, c in w.items():
                    nw[k] = fa * c
                for k, c in row.items():
                    nv = nw.get(k, 0) - fb * c
                    if nv:
                        nw[k] = nv
                    else:
                        nw.pop(k, None)
                w = _normalize(nw) if nw else nw
            else:
                piv[c0] = _normalize(w)
                r += 1
                break
    return r


def ranks(vecs):
    """(rank mod P1, rank mod P2, exact rank)."""
    order = _keyorder(vecs)
    return rank_modp(vecs, P1, order), rank_modp(vecs, P2, order), rank_exact(vecs, order)


def kernel_exact(images, sub, nx):
    """kernel of x_i -> images[i] (dict over Y-keys) modulo span(sub) (dicts over Y-keys).
    returns list of integer vectors (dict i -> c) spanning the kernel."""
    vecs = []
    for i, im in enumerate(images):
        v = {("y", k): c for k, c in im.items()}
        v[("x", i)] = 1
        vecs.append(v)
    for s in sub:
        vecs.append({("y", k): c for k, c in s.items()})
    ykeys = sorted({k for v in vecs for k in v if k[0] == "y"}, key=repr)
    order = {k: i for i, k in enumerate(ykeys)}
    for i in range(nx):
        order[("x", i)] = len(ykeys) + i
    piv = {}
    for v in vecs:
        w = {order[k]: c for k, c in v.items() if c}
        while w:
            c0 = min(w)
            if c0 in piv:
                row = piv[c0]
                a = row[c0]; b = w[c0]
                g = gcd(a, b); fa, fb = a // g, b // g
                nw = {k: fa * c for k, c in w.items()}
                for k, c in row.items():
                    nv = nw.get(k, 0) - fb * c
                    if nv:
                        nw[k] = nv
                    else:
                        nw.pop(k, None)
                w = _normalize(nw) if nw else nw
            else:
                piv[c0] = _normalize(w)
                break
    ny = len(ykeys)
    out = []
    for c0, row in piv.items():
        if c0 >= ny:
            assert all(k >= ny for k in row)
            out.append({k - ny: c for k, c in row.items()})
    return out
