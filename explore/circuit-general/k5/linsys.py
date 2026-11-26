"""Sparse linear algebra for relation systems: modular rank / nullspace analysis and exact solving."""
from fractions import Fraction as F

P = 2147483647


def inv(a, p=P):
    return pow(a % p, p - 2, p)


def to_mod(x, p=P):
    x = F(x)
    return x.numerator % p * inv(x.denominator, p) % p


def rref_mod(rows, ncols, p=P, col_order=None):
    """rows: list of dict col -> int (mod p).  Gaussian elimination with columns eliminated in col_order
    (pivot preference).  Returns (pivot_rows dict col -> row dict normalized with pivot 1)."""
    if col_order is None:
        col_order = list(range(ncols))
    rank_of = {c: i for i, c in enumerate(col_order)}
    piv = {}
    for r in rows:
        r = {c: v % p for c, v in r.items() if v % p}
        while r:
            c = min(r, key=lambda c: rank_of[c])
            if c in piv:
                f = r[c]
                pr = piv[c]
                for cc, vv in pr.items():
                    nv = (r.get(cc, 0) - f * vv) % p
                    if nv:
                        r[cc] = nv
                    else:
                        r.pop(cc, None)
            else:
                f = inv(r[c], p)
                r = {cc: vv * f % p for cc, vv in r.items()}
                # back-substitute into existing pivots to keep reduced form lazily (not needed)
                piv[c] = r
                break
    return piv


def full_reduce(piv, col_order, p=P):
    """make the echelon form reduced (each pivot column appears only in its pivot row)"""
    rank_of = {c: i for i, c in enumerate(col_order)}
    cols = sorted(piv, key=lambda c: rank_of[c], reverse=True)
    for c in cols:
        r = piv[c]
        # eliminate later pivot columns from r
        changed = True
        while changed:
            changed = False
            for cc in list(r):
                if cc != c and cc in piv and rank_of[cc] > rank_of[c]:
                    f = r[cc]
                    for c2, v2 in piv[cc].items():
                        nv = (r.get(c2, 0) - f * v2) % p
                        if nv:
                            r[c2] = nv
                        else:
                            r.pop(c2, None)
                    changed = True
        piv[c] = r
    return piv


def determined(piv, unknown_cols, col_order, p=P):
    """unknown columns fixed by the system once the remaining (known) columns are given:
    after full reduction with unknown columns first, unknown c is determined iff c is a pivot
    and its row contains no other unknown column."""
    U = set(unknown_cols)
    piv = full_reduce(piv, col_order, p)
    det = set()
    for c in U:
        if c in piv and all((cc == c or cc not in U) for cc in piv[c]):
            det.add(c)
    return det, piv


def solve_exact(rows, unknown, known_vals):
    """exact Fraction elimination.  rows: list of dict key -> Fraction coefficient (sum = 0).
    unknown: ordered list of keys to solve for; known_vals: key -> value (value = dict label -> Fraction).
    Returns dict key -> value for determined unknowns, plus list of inconsistencies (rows reducing to
    nonzero constant)."""
    U = {k: i for i, k in enumerate(unknown)}
    piv = {}
    bad = []
    for r in rows:
        # split into unknown part and constant (value) part
        ru = {}
        const = {}
        skip = False
        for k, c in r.items():
            if c == 0:
                continue
            if k in U:
                ru[U[k]] = ru.get(U[k], 0) + F(c)
            elif k in known_vals:
                for lab, x in known_vals[k].items():
                    const[lab] = const.get(lab, 0) + F(c) * x
            else:
                skip = True
                break
        if skip:
            continue
        ru = {k: v for k, v in ru.items() if v}
        const = {k: v for k, v in const.items() if v}
        while ru:
            c = min(ru)
            if c in piv:
                pr, pc = piv[c]
                f = ru[c]
                for cc, vv in pr.items():
                    nv = ru.get(cc, 0) - f * vv
                    if nv:
                        ru[cc] = nv
                    else:
                        ru.pop(cc, None)
                for lab, x in pc.items():
                    nv = const.get(lab, 0) - f * x
                    if nv:
                        const[lab] = nv
                    else:
                        const.pop(lab, None)
            else:
                f = 1 / ru[c]
                ru = {cc: vv * f for cc, vv in ru.items()}
                const = {lab: x * f for lab, x in const.items()}
                piv[c] = (ru, const)
                break
        else:
            if const:
                bad.append(const)
    # full back substitution: each pivot expressed through free columns + constants
    expr = {}
    for c in sorted(piv, reverse=True):
        ru, const = piv[c]
        free = {}
        val = {lab: -x for lab, x in const.items()}
        for cc, vv in ru.items():
            if cc == c:
                continue
            if cc in expr:
                fr, vl = expr[cc]
                for f2, w in fr.items():
                    nv = free.get(f2, 0) - vv * w
                    if nv:
                        free[f2] = nv
                    else:
                        free.pop(f2, None)
                for lab, x in vl.items():
                    nv = val.get(lab, 0) - vv * x
                    if nv:
                        val[lab] = nv
                    else:
                        val.pop(lab, None)
            else:
                nv = free.get(cc, 0) - vv
                if nv:
                    free[cc] = nv
                else:
                    free.pop(cc, None)
        expr[c] = (free, val)
    sol = {unknown[c]: {lab: x for lab, x in val.items() if x} for c, (free, val) in expr.items() if not free}
    return sol, bad, len(piv)
