"""Minimal exterior / Grassmann algebra with sympy coefficients.

Element: dict {tuple of sorted generator indices: sympy expression}; all generators odd.
"""
import sympy as sp


def _merge(a, b):
    """sorted merge of two sorted tuples of odd generators: (tuple, sign) or None."""
    out, i, j, inv = [], 0, 0, 0
    la, lb = len(a), len(b)
    while i < la and j < lb:
        if a[i] < b[j]:
            out.append(a[i]); i += 1
        elif a[i] > b[j]:
            out.append(b[j]); j += 1
            inv += la - i
        else:
            return None
    out.extend(a[i:]); out.extend(b[j:])
    return tuple(out), (-1 if inv & 1 else 1)


def clean(p):
    out = {}
    for k, v in p.items():
        v = sp.expand(v)
        if v != 0:
            out[k] = v
    return out


def add(*ps):
    r = {}
    for p in ps:
        for k, v in p.items():
            r[k] = r.get(k, 0) + v
    return clean(r)


def scale(p, c):
    return clean({k: c * v for k, v in p.items()})


def wedge(p, q):
    r = {}
    for a, ca in p.items():
        for b, cb in q.items():
            m = _merge(a, b)
            if m is None:
                continue
            key, s = m
            r[key] = r.get(key, 0) + s * ca * cb
    return clean(r)


def wedge_all(ps):
    out = {(): sp.Integer(1)}
    for p in ps:
        out = wedge(out, p)
    return out


def gen(i, c=1):
    return {(i,): sp.sympify(c)}


def one(c=1):
    return {(): sp.sympify(c)}


def interior(p, vec):
    """iota_V p for V = sum_i vec[i] d/d(gen i) (odd derivation acting from the left)."""
    r = {}
    for key, c in p.items():
        for pos, g in enumerate(key):
            if g in vec:
                k2 = key[:pos] + key[pos + 1:]
                r[k2] = r.get(k2, 0) + (-1) ** pos * vec[g] * c
    return clean(r)


def coeff_apply(p, f):
    """apply the function f to every coefficient (even operator)."""
    return clean({k: f(v) for k, v in p.items()})


def d_op(p, pairs):
    """sum_g gen_g wedge (D_g p) for pairs = [(g, D_g)], D_g acting on coefficients."""
    out = {}
    for g, Dg in pairs:
        out = add(out, wedge(gen(g), coeff_apply(p, Dg)))
    return out


def iota_op(p, pairs):
    """sum_g iota_{gen g} (D_g p) for pairs = [(g, D_g)]."""
    out = {}
    for g, Dg in pairs:
        out = add(out, interior(coeff_apply(p, Dg), {g: 1}))
    return out


def is_zero(p):
    return all(sp.simplify(v) == 0 for v in p.values())
