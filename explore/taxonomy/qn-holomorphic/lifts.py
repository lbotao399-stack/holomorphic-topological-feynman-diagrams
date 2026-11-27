"""Q0-cohomology of gauge-invariant local operators and the one-loop differential d_1 (modulo P).

Physical complex [src LLC sec-summary]: g-invariant polynomials in jet letters, independent of the undifferentiated gauge
ghost c_0 (letters of the gauge block 'g' with jet 0).  On such operators Q0 = c_0^a L_a + Q_rel, so Q_rel O = Q0 O with
c_0 terms dropped (they cancel on invariants).
Sectors: (q, d, k): q = block charges (integer gradings of g_tot constant on blocks, solved from the structure constants),
d = derivative vector (N^H), k = number of letters.  Q0: (q, d, k) -> (q, d, k+1);  Q_n: (q, d, k) -> (q + n dq, d + n, k + (H-3)n + 1).
Invariants: Cartan weight 0 and killed by the positive root vectors (weight bases: sl(N), sp(2N) matrix units).
d_1 rank = rank[Q1 Z_S ; B_T] - rank B_T  (Z_S = ker Q0 on Inv_S, B_T = Q0 Inv_{T-})."""
import itertools
import sys
import time
from fractions import Fraction as Fr

import numpy as np
import sympy as sp

import polyops as po


# =============================================================================== linear algebra mod P
def rref(rows, ncols):
    """rows: list of dict col -> val.  Returns list of pivot rows (dicts, normalized) and pivot cols."""
    P = po.P
    piv = {}
    for r in rows:
        r = {c: v % P for c, v in r.items() if v % P}
        while r:
            c = min(r)
            if c in piv:
                pr = piv[c]
                f = r[c]
                for cc, vv in pr.items():
                    nv = (r.get(cc, 0) - f * vv) % P
                    if nv:
                        r[cc] = nv
                    else:
                        r.pop(cc, None)
            else:
                inv = pow(r[c], P - 2, P)
                r = {cc: vv * inv % P for cc, vv in r.items()}
                piv[c] = r
                break
    return piv


def rank(rows):
    return len(rref(rows, None))


def kernel(rows_by_col, nvars):
    """kernel of the linear map given as columns: rows_by_col[j] = dict(target index -> value) image of basis j.
    Returns list of kernel vectors (dict j -> value)."""
    P = po.P
    # Gaussian elimination on the transpose: find relations sum_j x_j col_j = 0
    piv = {}          # target col -> (row dict over targets, combination dict over j)
    ker = []
    for j in range(nvars):
        r = dict(rows_by_col[j])
        comb = {j: 1}
        while r:
            c = min(r)
            if c in piv:
                pr, pc = piv[c]
                f = r[c]
                for cc, vv in pr.items():
                    nv = (r.get(cc, 0) - f * vv) % P
                    if nv:
                        r[cc] = nv
                    else:
                        r.pop(cc, None)
                for cc, vv in pc.items():
                    nv = (comb.get(cc, 0) - f * vv) % P
                    if nv:
                        comb[cc] = nv
                    else:
                        comb.pop(cc, None)
            else:
                inv = pow(r[c], P - 2, P)
                piv[c] = ({cc: vv * inv % P for cc, vv in r.items()}, {cc: vv * inv % P for cc, vv in comb.items()})
                break
        if not r:
            ker.append(comb)
    return ker


# =============================================================================== theory data
class LiftTheory:
    def __init__(self, alg, gname_block='g', cartan=None, posroots=None, charges=None):
        """alg: wick.Algebra; cartan: list of basis indices of g (Cartan elements); posroots: basis indices of the
        positive root vectors of g (both in the gauge block of the L-part)."""
        self.alg = alg
        sl = alg.sl
        self.H = alg.H
        self.N = alg.N
        self.pi = alg.pi
        self.gblock = sl.blocks[gname_block]
        self.blocks = sl.blocks
        # block partition (finest constant-parity)
        from blocksutil import blocks_for
        self.part = blocks_for(sl)
        self.blk_of = {}
        for b, ids in self.part.items():
            for i in ids:
                self.blk_of[i] = b
        self.cartan = cartan
        self.posroots = posroots
        self._weights()
        self.charges = charges if charges is not None else self._charges()
        self.wt = {A: tuple(int(c[self.blk_of[A]]) for c in self.charges) for A in range(self.N)}

    def _weights(self):
        """Cartan weights of the basis: [H_i, T_A] = w_i(A) T_A."""
        f = self.alg.sl.f
        self.hw = []
        for A in range(self.N):
            w = []
            for h in self.cartan:
                row = f[h, A]
                nz = np.nonzero(row)[0]
                if len(nz) == 0:
                    w.append(0)
                    continue
                assert list(nz) == [A], 'basis not a weight basis at %d' % A
                w.append(int(row[A]))
            self.hw.append(tuple(w))

    def _charges(self):
        """integer gradings constant on blocks with w(b1)+w(b2)=w(b3) for all nonzero brackets."""
        f = self.alg.sl.f
        bn = sorted(self.part)
        bi = {b: i for i, b in enumerate(bn)}
        eqs = set()
        for A, B, C in zip(*np.nonzero(f)):
            e = [0] * len(bn)
            e[bi[self.blk_of[A]]] += 1
            e[bi[self.blk_of[B]]] += 1
            e[bi[self.blk_of[C]]] -= 1
            eqs.add(tuple(e))
        M = sp.Matrix(sorted(eqs)) if eqs else sp.zeros(1, len(bn))
        ns = M.nullspace()
        out = []
        for v in ns:
            den = sp.ilcm(*[sp.fraction(x)[1] for x in v])
            vv = [int(x * den) for x in v]
            out.append({b: vv[bi[b]] for b in bn})
        return out

    # ------------------------------------------------------------------ letters
    def is_c0(self, letter):
        A, al = letter
        return A in self.gblock and not any(al)

    def letter_charge(self, letter):
        return self.wt[letter[0]]


def letters_upto(th, dmax):
    """all letters (A, alpha) with |alpha| <= dmax componentwise bounded by dmax vector, excluding c_0."""
    jets = [al for al in itertools.product(*[range(x + 1) for x in dmax])]
    out = []
    for A in range(th.N):
        for al in jets:
            if th.is_c0((A, al)):
                continue
            out.append((A, al))
    return out


def sector_monomials(th, q, d, k, wt0=True):
    """all monomials with k letters, total charge q, total jet d, total Cartan weight 0 (if wt0)."""
    H = th.H
    lets = letters_upto(th, d)
    lets.sort()
    nq = None if q is None else len(q)
    res = []
    hw = th.hw
    pi = th.pi

    def rec(start, kk, qrem, drem, wrem, cur):
        if kk == 0:
            if (qrem is None or all(x == 0 for x in qrem)) and all(x == 0 for x in drem) and \
                    (not wt0 or all(x == 0 for x in wrem)):
                res.append(tuple(cur))
            return
        for i in range(start, len(lets)):
            A, al = lets[i]
            if any(a > r for a, r in zip(al, drem)):
                continue
            if cur and lets[i] == cur[-1] and pi[A]:
                continue
            w = th.wt[A]
            nq_ = None if qrem is None else tuple(a - b for a, b in zip(qrem, w))
            nd = tuple(a - b for a, b in zip(drem, al))
            nw = tuple(a - b for a, b in zip(wrem, hw[A]))
            cur.append(lets[i])
            rec(i if not pi[A] else i + 1, kk - 1, nq_, nd, nw, cur)
            cur.pop()
    rec(0, k, None if q is None else tuple(q), tuple(d), tuple([0] * len(th.cartan)), [])
    return res


def g_action(th, a, mono):
    """L_a (even derivation, coadjoint on letters): L_a phi^A = - sum_C f[a,C,A] phi^C."""
    f = th.alg.sl.f
    out = {}
    for k, (A, al) in enumerate(mono):
        col = f[a, :, A]
        for C in np.nonzero(col)[0]:
            word = mono[:k] + ((int(C), al),) + mono[k + 1:]
            s, m = po.canon(word, th.pi)
            if s:
                po.padd(out, m, -s * int(col[C]))
    return out


def invariants(th, monos):
    """basis of g-invariants (vectors over monos) = weight-0 space killed by positive root vectors."""
    idx = {m: i for i, m in enumerate(monos)}
    cols = []
    for m in monos:
        img = {}
        for r, a in enumerate(th.posroots):
            for mm, v in g_action(th, a, m).items():
                img[(r, mm)] = (img.get((r, mm), 0) + v) % po.P
        cols.append({key: v for key, v in img.items() if v})
    # map keys to integers
    keymap = {}
    cols2 = []
    for c in cols:
        cols2.append({keymap.setdefault(key, len(keymap)): v for key, v in c.items()})
    return kernel(cols2, len(monos))


def vec_to_poly(vec, monos):
    return {monos[j]: v for j, v in vec.items()}


def drop_c0(th, poly):
    return {m: v for m, v in poly.items() if not any(th.is_c0(l) for l in m)}


def poly_to_row(poly, index):
    """poly over target monomials -> dict index -> value; raises if a monomial is outside the index."""
    out = {}
    for m, v in poly.items():
        if m not in index:
            raise KeyError('monomial outside target sector: %s' % (m,))
        out[index[m]] = v
    return out
