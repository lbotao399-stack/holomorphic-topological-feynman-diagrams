"""Finite-N tree-level cohomology H(Q_0) of the relative complex for U(N) or SU(N), per sector.

Invariants are spanned by multi-traces; at finite N they are evaluated on N x N matrices of letters
(graded-commutative component variables).  For a sector lam and letter number Y:
    dim Inv^Y      = rank ev(multi-traces with Y letters)
    rank Q_0^Y     = rank ev(Q_r(multi-traces with Y letters))
    dim H^Y        = dim Inv^Y - rank Q_0^Y - rank Q_0^{Y-1}.
SU(N): the last diagonal entry is replaced by minus the sum of the others.
"""
from itertools import product as iproduct
from . import rank_mod


def merge_sign(A, B, par):
    """Graded product of two sorted symbol tuples.  par(sym) -> 0/1.  Returns (sign, tuple) or (0, None)."""
    out = []
    i = j = 0
    sign = 1
    nA = len(A)
    # number of odd symbols remaining in A from position i
    oddA_suffix = [0] * (nA + 1)
    for k in range(nA - 1, -1, -1):
        oddA_suffix[k] = oddA_suffix[k + 1] + par(A[k])
    while i < nA and j < len(B):
        a, b = A[i], B[j]
        if a < b:
            out.append(a)
            i += 1
        elif b < a:
            if par(b) and (oddA_suffix[i] & 1):
                sign = -sign
            out.append(b)
            j += 1
        else:
            if par(a):
                return 0, None
            out.append(a)
            i += 1
    out.extend(A[i:])
    out.extend(B[j:])
    return sign, tuple(out)


def sort_sign(seq, par):
    """Sort a sequence of symbols with Koszul sign; repeated odd symbol -> 0."""
    seq = list(seq)
    sign = 1
    # insertion sort counting odd transpositions
    for k in range(1, len(seq)):
        x = seq[k]
        px = par(x)
        m = k - 1
        while m >= 0 and seq[m] > x:
            if px and par(seq[m]):
                sign = -sign
            seq[m + 1] = seq[m]
            m -= 1
        seq[m + 1] = x
    for k in range(1, len(seq)):
        if seq[k] == seq[k - 1] and par(seq[k]):
            return 0, None
    return sign, tuple(seq)


class MultiTraceSector:
    """Multi-trace monomials (sorted tuples of canonical cyclic words) of total weight lam."""

    def __init__(self, alg, lam):
        self.alg = alg
        self.lam = tuple(lam)
        cyc = {}
        for w in alg.words_upto(lam):
            s, cw = alg.canon_cyclic(w)
            if s != 0:
                cyc[cw] = True
        self.traces = sorted(cyc, key=lambda w: (len(w), w))
        self.monos = self._enumerate()

    def _enumerate(self):
        alg = self.alg
        gens = self.traces
        W = [alg.wword(w) for w in gens]
        P = [alg.pword(w) for w in gens]
        out = []

        def rec(start, cur, rem):
            if all(r == 0 for r in rem):
                out.append(tuple(cur))
                return
            for idx in range(start, len(gens)):
                wi = W[idx]
                if all(a <= b for a, b in zip(wi, rem)):
                    nrem = tuple(b - a for a, b in zip(wi, rem))
                    cur.append(gens[idx])
                    rec(idx + 1 if P[idx] else idx, cur, nrem)
                    cur.pop()
        rec(0, [], self.lam)
        return out

    def Y(self, mono):
        return sum(len(w) for w in mono)

    def by_Y(self):
        d = {}
        for m in self.monos:
            d.setdefault(self.Y(m), []).append(m)
        return d

    def Q_mono(self, mono):
        """Q_r on a multi-trace: dict mono -> coeff."""
        alg = self.alg
        out = {}
        pre = 0
        for pos, w in enumerate(mono):
            s0 = -1 if pre else 1
            for nw, c in alg.Qword(w).items():
                s, cw = alg.canon_cyclic(nw)
                if s == 0:
                    continue
                rest = list(mono[:pos]) + [cw] + list(mono[pos + 1:])
                sg, key = sort_sign_gens(rest, alg)
                if sg == 0:
                    continue
                out[key] = out.get(key, 0) + s0 * s * sg * c
            pre ^= alg.pword(w)
        return {k: v for k, v in out.items() if v}


def sort_sign_gens(gens, alg):
    """Sort generators (words) with Koszul sign by word parity; repeated odd -> 0."""
    return sort_sign(gens, lambda w: alg.pword(w))


class Evaluator:
    """Evaluate traces of words on N x N matrices (U(N) or SU(N)) into component polynomials."""

    def __init__(self, alg, N, special=False):
        self.alg = alg
        self.N = N
        self.special = special
        self.cache = {}

    def par(self, sym):
        return self.alg.par[sym[0]]

    def entry(self, x, i, j):
        """List of (coeff, symbol) for the (i,j) entry of letter x."""
        N = self.N
        if self.special and i == j == N - 1:
            return [(-1, (x, k, k)) for k in range(N - 1)]
        return [(1, (x, i, j))]

    def trace(self, w):
        if w in self.cache:
            return self.cache[w]
        N = self.N
        n = len(w)
        poly = {}
        for idx in iproduct(range(N), repeat=n):
            # product of entries x_k^{idx[k], idx[k+1]}
            terms = [(1, ())]
            for k in range(n):
                i, j = idx[k], idx[(k + 1) % n]
                new = []
                for c, seq in terms:
                    for c2, sym in self.entry(w[k], i, j):
                        new.append((c * c2, seq + (sym,)))
                terms = new
            for c, seq in terms:
                s, key = sort_sign(seq, self.par)
                if s:
                    poly[key] = poly.get(key, 0) + c * s
        poly = {k: v for k, v in poly.items() if v}
        self.cache[w] = poly
        return poly

    def mono(self, mono):
        acc = {(): 1}
        for w in mono:
            tw = self.trace(w)
            new = {}
            for k1, c1 in acc.items():
                for k2, c2 in tw.items():
                    s, k = merge_sign(k1, k2, self.par)
                    if s:
                        new[k] = new.get(k, 0) + s * c1 * c2
            acc = {k: v for k, v in new.items() if v}
            if not acc:
                break
        return acc


def cohomology_sector(alg, lam, N, special=False, field=rank_mod.P1):
    """Returns {Y: (dim Inv^Y, rank Q^Y, dim H^Y)} for U(N) (special=False) or SU(N)."""
    S = MultiTraceSector(alg, lam)
    byY = S.by_Y()
    ev = Evaluator(alg, N, special)
    dimInv = {}
    rankQ = {}
    for Y, monos in byY.items():
        sp = rank_mod.make_span(field)
        for m in monos:
            sp.add(ev.mono(m))
        dimInv[Y] = sp.dim()
    for Y, monos in byY.items():
        sp = rank_mod.make_span(field)
        for m in monos:
            q = S.Q_mono(m)
            v = {}
            for mm, c in q.items():
                for k, cc in ev.mono(mm).items():
                    v[k] = v.get(k, 0) + c * cc
            sp.add({k: x for k, x in v.items() if x})
        rankQ[Y] = sp.dim()
    out = {}
    for Y in sorted(byY):
        h = dimInv[Y] - rankQ.get(Y, 0) - rankQ.get(Y - 1, 0)
        out[Y] = (dimInv[Y], rankQ.get(Y, 0), h)
    return out


def cohomology_sector_stable(alg, lam, field=rank_mod.P1):
    """Large-N (no trace relations) cohomology of (C[T], Q_r) in sector lam: {Y: dim H^Y}."""
    S = MultiTraceSector(alg, lam)
    byY = S.by_Y()
    rankQ = {}
    for Y, monos in byY.items():
        sp = rank_mod.make_span(field)
        for m in monos:
            q = S.Q_mono(m)
            if q:
                sp.add(q)
        rankQ[Y] = sp.dim()
    return {Y: len(byY[Y]) - rankQ.get(Y, 0) - rankQ.get(Y - 1, 0) for Y in sorted(byY)}
