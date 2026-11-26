"""Q_0-cohomology in graded sectors and the induced d_2 = [Q_2] (exact Fractions).
Sector: fixed (J0, f, af) with J0 = N_B + N_psi + D, f = N_phi - N_psi (per multiplet or total),
af = N_B + N_psi; complex graded by gh0 = N_c + N_psi.  Q_0 preserves (J0, f, af) when W = kappa = 0;
Q_2 preserves (J0, f), lowers af by 2, raises gh0 by 1."""
from fractions import Fraction as Fr
from itertools import combinations_with_replacement
from ht3d import mono, padd


def letters_upto(th, J0max):
    """all letters with J0-weight <= J0max: (letter, weight)"""
    out = []
    for d in range(J0max + 1):
        for a in range(th.G):
            out.append((('c', a, d), d))
        for a in range(th.G):
            if 1 + d <= J0max:
                out.append((('B', a, d), 1 + d))
        for k, m in enumerate(th.multiplets):
            for i in range(m[0]):
                out.append((('p%d' % k, i, d), d))
                if 1 + d <= J0max:
                    out.append((('s%d' % k, i, d), 1 + d))
    return out


def grading(m):
    J0 = f = af = gh = 0
    for (n, i, d) in m:
        if n == 'c':
            J0 += d; gh += 1
        elif n == 'B':
            J0 += 1 + d; af += 1
        elif n[0] == 'p':
            J0 += d; f += 1
        else:
            J0 += 1 + d; f -= 1; af += 1; gh += 1
    return J0, f, af, gh


def sector_basis(th, J0, f, af):
    """monomials with given (J0, f, af), grouped by gh0. N_phi bounded by f + N_psi <= f + J0."""
    A = th.A
    L = letters_upto(th, J0)
    res = {}
    # recursive enumeration with weight budget
    Ls = sorted(L)

    def rec(start, cur, J, nphi, npsi, nB):
        if J == J0 and nphi - npsi == f and nB + npsi == af:
            s, m = A.canon([x[0] for x in cur])
            if s:
                res.setdefault(grading(m)[3], set()).add(m)
        for idx in range(start, len(Ls)):
            (Lt, w) = Ls[idx]
            n = Lt[0]
            if J + w > J0:
                continue
            dphi = 1 if n[0] == 'p' else 0
            dpsi = 1 if n[0] == 's' else 0
            dB = 1 if n == 'B' else 0
            if nB + npsi + dB + dpsi > af:
                continue
            if nphi + dphi > f + af:  # N_phi = f + N_psi <= f + af
                continue
            if A.parity(Lt) and cur and cur[-1][0] == Lt:
                continue
            rec(idx, cur + [(Lt, w)], J + w, nphi + dphi, npsi + dpsi, nB + dB)

    rec(0, [], 0, 0, 0, 0)
    return {g: sorted(v) for g, v in res.items()}


# ---------------- sparse exact linear algebra
class Echelon:
    """incremental row echelon form over Q; rows as dicts {col: Fraction}"""

    def __init__(self):
        self.rows = {}  # pivot col -> row (normalized pivot 1)

    def reduce(self, v):
        v = dict(v)
        while v:
            p = min(v)
            if p in self.rows:
                c = v[p]
                for k, x in self.rows[p].items():
                    nv = v.get(k, 0) - c * x
                    if nv == 0:
                        v.pop(k, None)
                    else:
                        v[k] = nv
            else:
                return v
        return v

    def add(self, v):
        v = self.reduce(v)
        if not v:
            return False
        p = min(v)
        c = v[p]
        row = {k: x / c for k, x in v.items()}
        # keep fully reduced w.r.t. new pivot not required (reduce loops on min pivot)
        self.rows[p] = row
        return True

    def rank(self):
        return len(self.rows)


def kernel_basis(cols, image_vecs, ncols_target):
    """kernel of the linear map sending basis vector j -> image_vecs[j] (dict over target index).
    Returns list of dicts over source index."""
    # Gaussian elimination on augmented [image | identity]
    piv = {}  # pivot target col -> (img, src combination)
    ker = []
    for j, img in enumerate(image_vecs):
        v = dict(img)
        comb = {j: Fr(1)}
        while v:
            p = min(v)
            if p in piv:
                pi, pc = piv[p]
                c = v[p]
                for k, x in pi.items():
                    nv = v.get(k, 0) - c * x
                    if nv == 0:
                        v.pop(k, None)
                    else:
                        v[k] = nv
                for k, x in pc.items():
                    nv = comb.get(k, 0) - c * x
                    if nv == 0:
                        comb.pop(k, None)
                    else:
                        comb[k] = nv
            else:
                break
        if not v:
            ker.append(comb)
        else:
            p = min(v)
            c = v[p]
            piv[p] = ({k: x / c for k, x in v.items()}, {k: x / c for k, x in comb.items()})
    return ker


def poly_to_vec(P, index):
    v = {}
    for m, c in P.items():
        if m not in index:
            raise KeyError("monomial outside sector: %s" % (m,))
        v[index[m]] = c
    return v


def vec_to_poly(v, basis):
    return {basis[i]: c for i, c in v.items()}
