"""Exact linear algebra for cochains (sympy sparse DomainMatrix over QQ and GF(p))."""
from fractions import Fraction as Fr
from sympy.polys.matrices.sdm import SDM
from sympy.polys.matrices import DomainMatrix
from sympy import QQ, GF
from vlib import LP

P1, P2 = 2147483629, 2147483587


def coeff_at(v, N):
    if isinstance(v, LP):
        return v.at(N)
    return Fr(v)


def to_rows(cochains, N, proj=None):
    """cochains -> list of sparse rows {colkey: Fr}; proj(cochain, N) -> {colkey: Fr} optional."""
    rows = []
    for w in cochains:
        if proj is not None:
            r = proj(w, N)
        else:
            r = {}
            for key, v in w.items():
                x = coeff_at(v, N)
                if x:
                    r[key] = x
        rows.append(r)
    return rows


def _matrix(rows, dom):
    cols = {}
    for r in rows:
        for k in r:
            if k not in cols:
                cols[k] = len(cols)
    data = {}
    for i, r in enumerate(rows):
        rr = {}
        for k, x in r.items():
            if x:
                if dom == QQ:
                    rr[cols[k]] = QQ(x.numerator, x.denominator)
                else:
                    rr[cols[k]] = dom(x.numerator) / dom(x.denominator)
        if rr:
            data[i] = rr
    return SDM(data, (len(rows), max(1, len(cols))), dom), cols


def rank(rows, exact=True):
    """exact rank over Q (sparse rref); also mod two primes (returned for the log)."""
    out = {}
    if exact:
        M, _ = _matrix(rows, QQ)
        out['Q'] = len(M.rref()[1])
    for p in (P1, P2):
        M, _ = _matrix(rows, GF(p))
        out[p] = len(M.rref()[1])
    return out


def rankQ(rows):
    M, _ = _matrix(rows, QQ)
    return len(M.rref()[1])


def rankp(rows, p=P1):
    M, _ = _matrix(rows, GF(p))
    return len(M.rref()[1])


def solve_in_span(rows, target):
    """exact rational x with sum x_i rows[i] = target, or None."""
    allr = rows + [target]
    cols = {}
    for r in allr:
        for k in r:
            if k not in cols:
                cols[k] = len(cols)
    import sympy as sp
    n = len(rows)
    A = sp.zeros(len(cols), n)
    b = sp.zeros(len(cols), 1)
    for i, r in enumerate(rows):
        for k, x in r.items():
            A[cols[k], i] = sp.Rational(x.numerator, x.denominator)
    for k, x in target.items():
        b[cols[k], 0] = sp.Rational(x.numerator, x.denominator)
    dA = DomainMatrix.from_Matrix(A).convert_to(QQ)
    db = DomainMatrix.from_Matrix(b).convert_to(QQ)
    aug = dA.hstack(db)
    R, piv = aug.rref()
    if n in piv:
        return None
    x = [Fr(0)] * n
    Rl = R.to_Matrix()
    for r_i, c in enumerate(piv):
        val = Rl[r_i, n]
        x[c] = Fr(int(val.p), int(val.q))
    return x


def indep_rows(rows):
    """indices of a maximal independent subset of rows (exact, over Q), via rref of the transpose."""
    cols = {}
    for r in rows:
        for k in r:
            if k not in cols:
                cols[k] = len(cols)
    data = {}
    for i, r in enumerate(rows):
        for k, x in r.items():
            if x:
                data.setdefault(cols[k], {})[i] = QQ(x.numerator, x.denominator)
    M = SDM(data, (max(1, len(cols)), len(rows)), QQ)
    _, piv = M.rref()
    return list(piv)
