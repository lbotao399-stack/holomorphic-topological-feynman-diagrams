"""Minimal exterior algebra with sympy (or Fraction) coefficients.

A form is a dict {tuple of generator indices (strictly increasing): coefficient}.
Generators are integers; their order defines the sign convention of the canonical monomial.
"""
import sympy as sp


def _sort_sign(seq):
    """sign of the permutation sorting seq (distinct entries), and the sorted tuple; 0 if repeated."""
    seq = list(seq)
    if len(set(seq)) != len(seq):
        return 0, None
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s, tuple(sorted(seq))


def wedge(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            s, k = _sort_sign(ka + kb)
            if s == 0:
                continue
            out[k] = out.get(k, 0) + s * va * vb
    return {k: v for k, v in out.items() if v != 0}


def add(a, b, c=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v != 0}


def scale(a, c):
    return {k: c * v for k, v in a.items() if c * v != 0}


def one_form(coeffs):
    """coeffs: dict generator -> coefficient."""
    return {(g,): c for g, c in coeffs.items() if c != 0}


def wedge_all(forms):
    out = {(): 1}
    for f in forms:
        out = wedge(out, f)
    return out


def interior(vec, a):
    """contraction iota_vec a, vec: dict generator -> component."""
    out = {}
    for k, v in a.items():
        for pos, g in enumerate(k):
            if g in vec and vec[g] != 0:
                kk = k[:pos] + k[pos + 1:]
                out[kk] = out.get(kk, 0) + (-1) ** pos * vec[g] * v
    return {k: v for k, v in out.items() if v != 0}


def d_coeff(a, var, gen):
    """the operator gen-differential times d/dvar acting on coefficients: dx^gen wedge d(coeff)/dvar."""
    out = {}
    for k, v in a.items():
        s, kk = _sort_sign((gen,) + k)
        if s == 0:
            continue
        dv = sp.diff(v, var)
        if dv != 0:
            out[kk] = out.get(kk, 0) + s * dv
    return {k: v for k, v in out.items() if v != 0}


def simplify_form(a):
    out = {}
    for k, v in a.items():
        w = sp.simplify(v)
        if w != 0:
            out[k] = w
    return out


def equal(a, b):
    diff = add(a, b, -1)
    return all(sp.simplify(v) == 0 for v in diff.values())
