"""Colour weights W_133, W_223 from the index forms of run-1 claim c3, by my own Casimir contraction.

f(U,V,W) = kappa([U,V],W) = tr(UVW) - tr(VUW), kappa = tr.
W_133 = - f_xyz f_{a1 y'u} f_{a2 u'w} f_{a3 z's} f_{a4 s'r} f_{x'w'r'}
W_223 = - f_xyz f_{a1 x'u} f_{a2 y'w} f_{a3 z's} f_{a4 s'r} f_{u'w'r'}
Primed pairs contracted with kappa^{-1} (gl_N or sl_N Casimir).
Result: dict (sigma on 4 slots, npow) -> Fraction, slots 0..3 = a1..a4.
"""
from fractions import Fraction as Fr
from vce import perm_from_cycles, has_fixed


def f_terms(U, V, W):
    return [([(U, V, W)], 1), ([(V, U, W)], -1)]


def contract(terms, pairs, sl):
    """terms: list of (cycles(list of tuples of labels), npow, coef)."""
    for (h, hp) in pairs:
        new = []
        for cyc, npow, cf in terms:
            ci = [i for i, c in enumerate(cyc) if h in c]
            cj = [i for i, c in enumerate(cyc) if hp in c]
            assert len(ci) == 1 and len(cj) == 1
            ci, cj = ci[0], cj[0]
            rest = [c for i, c in enumerate(cyc) if i not in (ci, cj)]
            if ci != cj:
                c1 = cyc[ci]; a = c1.index(h); P = c1[a + 1:] + c1[:a]
                c2 = cyc[cj]; b = c2.index(hp); Q = c2[b + 1:] + c2[:b]
                outs = [([P + Q], 1, 0)]
                if sl:
                    outs.append(([P, Q], -1, -1))
                for cl, s, dn in outs:
                    cc = list(rest)
                    nn = npow + dn
                    for w in cl:
                        if len(w) == 0:
                            nn += 1
                        else:
                            cc.append(tuple(w))
                    new.append((cc, nn, cf * s))
            else:
                c = cyc[ci]
                a = c.index(h)
                c = c[a:] + c[:a]  # starts with h
                b = c.index(hp)
                A = c[1:b]
                B = c[b + 1:]
                outs = [([A, B], 1, 0)]
                if sl:
                    outs.append(([A + B], -1, -1))
                for cl, s, dn in outs:
                    cc = list(rest)
                    nn = npow + dn
                    for w in cl:
                        if len(w) == 0:
                            nn += 1
                        else:
                            cc.append(tuple(w))
                    new.append((cc, nn, cf * s))
        terms = new
    return terms


def weight(vertices, pairs, sl, ext=('a1', 'a2', 'a3', 'a4'), overall=-1):
    terms = [([], 0, Fr(overall))]
    for (U, V, W) in vertices:
        nt = []
        for cyc, npow, cf in terms:
            for cl, s in f_terms(U, V, W):
                nt.append((cyc + [tuple(x) for x in cl], npow, cf * s))
        terms = nt
    terms = contract(terms, pairs, sl)
    lab = {e: i for i, e in enumerate(ext)}
    out = {}
    for cyc, npow, cf in terms:
        cy = [tuple(lab[x] for x in c) for c in cyc]
        sig = perm_from_cycles(cy, len(ext))
        if sl and has_fixed(sig):
            continue
        key = (sig, npow)
        out[key] = out.get(key, 0) + cf
    return {k: v for k, v in out.items() if v != 0}


PAIRS = [(x, x + "'") for x in 'xyzuwsr']


def W133(sl):
    V = [('x', 'y', 'z'), ('a1', "y'", 'u'), ('a2', "u'", 'w'), ('a3', "z'", 's'), ('a4', "s'", 'r'),
         ("x'", "w'", "r'")]
    return weight(V, PAIRS, sl)


def W223(sl):
    V = [('x', 'y', 'z'), ('a1', "x'", 'u'), ('a2', "y'", 'w'), ('a3', "z'", 's'), ('a4', "s'", 'r'),
         ("u'", "w'", "r'")]
    return weight(V, PAIRS, sl)
