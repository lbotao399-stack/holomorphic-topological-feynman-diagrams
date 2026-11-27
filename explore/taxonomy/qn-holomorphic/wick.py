"""Algebraic (colour + Koszul) layer of the rooted graph operators, for a quadratic Lie superalgebra.

Labelled rooted graph: bulk vertices 0..V-1, root = V; edges (tail, head) in edge order, spokes (v, V).
Wick rule:  hatQ_Gamma O = K_Gamma * [prod_{e in edge order} G^{A_e B_e} d^{(tail)}_{A_e} d^{(head)}_{B_e}] (I_0 ... I_{V-1} O)
with I_v = (1/6) V_{ABC} phi^A phi^B phi^C; the 3! slot assignments at a vertex with distinct legs cancel 1/6.
Leg order at v: in-edges (edge order), outputs, out-edges (edge order).  Word = leg letters of I_0, ..., I_{V-1}, then O.
Removals: edges last to first; for each edge the head letter, then the tail letter; each removal of a letter of
parity q contributes (-1)^{q * (parity of the letters still to its left)}.  Remaining word = outputs (vertex order),
then the unused letters of O.
"""
import itertools
from functools import lru_cache

import numpy as np

import polyops as po


class Algebra:
    """g_tot with vertex V, propagator G = Gint/Gden; H fixes the theory dimension."""

    def __init__(self, sl, H, variant=0, check=True):
        self.sl = sl
        self.H = H
        if check:
            sl.check()
        assert sl.pomega % 2 == (H + 1) % 2, 'pairing parity must be H + 1'
        self.V = sl.vertex()
        self.G, self.Gden = sl.propagator(variant)
        self.pi = [int(x) for x in sl.pi]
        self.N = sl.N
        self.names = sl.names
        self.Vf = self.V.astype(np.float64)
        self.Gf = self.G.astype(np.float64)
        # sparse Q0 table: Q0 phi^D = 1/2 (-1)^{pi_D H} sum G[A,D] V[A,B,C] phi^B phi^C / Gden
        self._q0 = {}
        GV = np.einsum('ad,abc->dbc', self.G, self.V)
        for D in range(self.N):
            terms = []
            for B, C in zip(*np.nonzero(GV[D])):
                terms.append((int(B), int(C), int(GV[D, B, C])))
            self._q0[D] = terms
        self.inv2den = pow(2 * self.Gden, po.P - 2, po.P)

    def q0_letter(self, D, al):
        """hatQ0 on d^al phi^D (Poly)."""
        key = (D, al, po.P)
        if key in self._q0cache:
            return self._q0cache[key]
        s = -1 if (self.pi[D] * self.H) % 2 else 1
        out = {}
        for B, C, c in self._q0[D]:
            for be, ga, bc in po.binom_split(al):
                sg, m = po.canon(((B, be), (C, ga)), self.pi)
                if sg:
                    po.padd(out, m, s * sg * c * bc * self.inv2den)
        self._q0cache[key] = out
        return out

    _q0cache = {}

    def q0(self, poly):
        out = {}
        pi = self.pi
        for mono, coef in poly.items():
            par = 0
            for k, (D, al) in enumerate(mono):
                img = self.q0_letter(D, al)
                s = -1 if par % 2 else 1
                pre, post = mono[:k], mono[k + 1:]
                for m2, c2 in img.items():
                    sg, m = po.canon(pre + m2 + post, pi)
                    if sg:
                        po.padd(out, m, s * sg * coef * c2)
                par += pi[D]
        return out


def removal_sign(par, removals):
    """par: parities of the word letters; removals: positions in removal order."""
    present = [True] * len(par)
    s = 0
    for pos in removals:
        q = par[pos]
        if q:
            s += sum(par[i] for i in range(pos) if present[i])
        present[pos] = False
    return -1 if s % 2 else 1


class GraphAlg:
    """Leg structure and Koszul signs of a labelled rooted tree (or graph) for the Wick rule."""

    def __init__(self, nb, edges):
        self.nb = nb
        self.root = nb
        self.edges = [tuple(e) for e in edges]
        self.spokes = [i for i, e in enumerate(self.edges) if self.root in e]
        assert all(self.edges[i][1] == self.root for i in self.spokes)
        legs = {}
        nout = {}
        for v in range(nb):
            ins = [i for i, e in enumerate(self.edges) if e[1] == v]
            outs = [i for i, e in enumerate(self.edges) if e[0] == v]
            k = 3 - len(ins) - len(outs)
            assert k >= 0
            nout[v] = k
            legs[v] = [('e', i) for i in ins] + [('out', v, j) for j in range(k)] + [('e', i) for i in outs]
        self.legs, self.nout = legs, nout
        self.outputs = [(v, j) for v in range(nb) for j in range(nout[v])]
        # word positions of the bulk legs
        pos = {}
        p = 0
        for v in range(nb):
            for leg in legs[v]:
                if leg[0] == 'e':
                    pos[(leg[1], v)] = p
                else:
                    pos[('out', leg[1], leg[2])] = p
                p += 1
        self.pos, self.nbulk = pos, p

    def sign(self, leg_par, out_par, O_par, spoke_pos):
        """Koszul sign.  leg_par[(e, v)] parity of the letter at the end of edge e at bulk vertex v;
        out_par[(v, j)]; O_par: parities of the letters of O; spoke_pos[j] = position in O of the letter
        contracted by spoke j (spokes in edge order)."""
        par = [0] * (self.nbulk + len(O_par))
        for (e, v), q in leg_par.items():
            par[self.pos[(e, v)]] = q
        for (v, j), q in out_par.items():
            par[self.pos[('out', v, j)]] = q
        for k, q in enumerate(O_par):
            par[self.nbulk + k] = q
        sj = {e: j for j, e in enumerate(self.spokes)}
        rem = []
        for e in reversed(range(len(self.edges))):
            t, h = self.edges[e]
            if h == self.root:
                rem.append(self.nbulk + spoke_pos[sj[e]])
            else:
                rem.append(self.pos[(e, h)])
            rem.append(self.pos[(e, t)])
        return removal_sign(par, rem)
