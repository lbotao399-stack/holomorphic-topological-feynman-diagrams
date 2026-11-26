"""Minimal exterior algebra (verifier's own implementation).
A form is a dict {tuple of generator indices (strictly increasing): coefficient}.
Coefficients: sympy expressions or Fractions (anything supporting + * and == 0)."""


def perm_sign(seq):
    seq = list(seq)
    s = 1
    n = len(seq)
    for i in range(n):
        for j in range(i + 1, n):
            if seq[i] > seq[j]:
                s = -s
            elif seq[i] == seq[j]:
                return 0
    return s


def sort_word(word):
    """(sign, sorted tuple) of a word of odd generators; sign 0 if repeated."""
    s = perm_sign(word)
    return s, tuple(sorted(word))


def add(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + v
    return clean(out)


def scale(a, c):
    return clean({k: c * v for k, v in a.items()})


def clean(a, simp=None):
    out = {}
    for k, v in a.items():
        if simp is not None:
            v = simp(v)
        if v != 0:
            out[k] = v
    return out


def wedge(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            s, k = sort_word(ka + kb)
            if s == 0:
                continue
            out[k] = out.get(k, 0) + s * va * vb
    return clean(out)


def wedge_all(forms):
    out = {(): 1}
    for f in forms:
        out = wedge(out, f)
    return out


def gen(i, c=1):
    return {(i,): c}


def interior(vec_index, a):
    """iota_{d/dx^vec_index} with generator vec_index dual: removes the generator with sign."""
    out = {}
    for k, v in a.items():
        if vec_index in k:
            pos = k.index(vec_index)
            nk = k[:pos] + k[pos + 1:]
            out[nk] = out.get(nk, 0) + (-1) ** pos * v
    return clean(out)


def apply_coeff(a, f):
    return clean({k: f(v) for k, v in a.items()})
