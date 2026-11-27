"""Explicit graph operators hatQ_Gamma on polynomials in jet letters (exact modulo P).

hatQ_Gamma O = sum_{injective spoke maps} sum_{indices} [prod_v V_legs(v) prod_e G^{A_e B_e} / Gden^{|I|}] * koszul
               * [lambda^beta] I_Gamma(lambda; input jets) * (outputs d^{beta_v} phi^{B_v}, vertex order)(unused letters of O)
I_Gamma = (-1)^{n(n+1)/2} K_Gamma (htk K, S2 normalization), kernel polynomial from kern_h.GraphKernel.kernel_symbolic.
Q_n = c^n sum_Gamma kappa_Gamma hatQ_Gamma with kappa fixed by nilpotency (run scripts).
"""
import itertools
import time

import numpy as np

import polyops as po
import kern_h
from wick import GraphAlg


def _solve_gf2(rows, rhs, nun):
    """solve A x = b over GF(2); rows: list of sets of unknown indices.  Returns x or None."""
    A = [(set(r), b % 2) for r, b in zip(rows, rhs)]
    x = [None] * nun
    piv = []
    for col in range(nun):
        k = next((i for i, (r, _) in enumerate(A) if col in r and i not in [p for p, _ in piv]), None)
        if k is None:
            continue
        rk, bk = A[k]
        for i, (r, b) in enumerate(A):
            if i != k and col in r:
                A[i] = (r ^ rk, b ^ bk)
        piv.append((k, col))
    for i, (r, b) in enumerate(A):
        if not r and b:
            return None
    sol = [0] * nun
    for k, col in piv:
        r, b = A[k]
        sol[col] = b          # free unknowns set to 0 (trees: none)
    return sol


def find_loop_edge(nb, edges, internal):
    """an internal edge on the bulk cycle (unicyclic graphs), else None."""
    import networkx as nx
    Gb = nx.MultiGraph()
    Gb.add_nodes_from(range(nb))
    for e in internal:
        Gb.add_edge(*edges[e], key=e)
    try:
        cyc = nx.find_cycle(Gb)
    except nx.NetworkXNoCycle:
        return None
    return cyc[0][2]


def solve_legs(op, pin, out_par, q=None):
    """parities of all bulk legs from the spoke-letter parities pin, output parities out_par[(v, j)] and (unicyclic)
    the letter parity q at the tail of op.loop_edge.  Returns leg_par or None (inconsistent)."""
    pw = op.alg.sl.pomega % 2
    leg_par = {}
    for j, e in enumerate(op.spokes):
        t, _ = op.edges[e]
        leg_par[(e, t)] = (pin[j] + pw) % 2
    ui = {e: k for k, e in enumerate(op.internal)}
    rows, rhs = [], []
    for v in range(op.nb):
        r, b = set(), (pw + 1) % 2
        for leg in op.ga.legs[v]:
            if leg[0] == 'out':
                b += out_par[(v, leg[2])]
            else:
                e = leg[1]
                if e in ui:
                    t, h = op.edges[e]
                    r ^= {ui[e]}
                    if v == h:
                        b += pw
                else:
                    b += leg_par[(e, v)]
        rows.append(r)
        rhs.append(b)
    if op.loop_edge is not None:
        rows.append({ui[op.loop_edge]})
        rhs.append(q)
    sol = _solve_gf2(rows, rhs, len(op.internal))
    if sol is None:
        return None
    for e, k in ui.items():
        t, h = op.edges[e]
        leg_par[(e, t)] = sol[k]
        leg_par[(e, h)] = (sol[k] + pw) % 2
    return leg_par


class GraphOp:
    def __init__(self, alg, g, name=None):
        """alg: wick.Algebra; g: dict(nb, edges, aut, kind) from rgraphs."""
        self.alg = alg
        self.H = alg.H
        self.nb, self.edges = g['nb'], [tuple(e) for e in g['edges']]
        self.aut = g['aut']
        self.kind = g['kind']
        self.name = name or g.get('name', 'G')
        self.ga = GraphAlg(self.nb, self.edges)
        self.n = len(self.edges) - self.nb
        self.gk = kern_h.GraphKernel(self.H, self.nb, self.edges, name=self.name, aut=self.aut)
        self.spokes = self.ga.spokes
        self.internal = [e for e in range(len(self.edges)) if e not in self.spokes]
        self.outv = [(v, j) for (v, j) in self.ga.outputs]
        assert all(j == 0 for v, j in self.outv), 'one output per vertex (n >= 1)'
        self.loop_edge = find_loop_edge(self.nb, self.edges, self.internal) if self.kind == 'unicyclic' else None
        self.qvals = (0, 1) if self.loop_edge is not None else (None,)
        self.loop_proj = {q: np.array([1.0 if alg.pi[A] == q else 0.0 for A in range(alg.N)]) for q in (0, 1)}
        self.nG = len(self.edges)
        self.inv_den = pow(alg.Gden ** self.nG, po.P - 2, po.P)
        self._sign_cache = {}
        self._kcache = {}
        self._einsum_spec()

    # ---------------------------------------------------------------- colour contraction
    def _einsum_spec(self):
        sym = {}

        def s(key):
            if key not in sym:
                sym[key] = len(sym)
            return sym[key]
        ops = []           # ('V', idx list) / ('G', [t, h]) / ('g', j, [idx])
        for v in range(self.nb):
            idx = []
            for leg in self.ga.legs[v]:
                if leg[0] == 'out':
                    idx.append(s(('out', v)))
                else:
                    idx.append(s(('leg', leg[1], v)))
            ops.append(('V', idx))
        for e in self.internal:
            t, h = self.edges[e]
            ops.append(('G', [s(('leg', e, t)), s(('leg', e, h))]))
        for j, e in enumerate(self.spokes):
            t, _ = self.edges[e]
            ops.append(('g', j, [s(('leg', e, t))]))
        if self.loop_edge is not None:
            t, _ = self.edges[self.loop_edge]
            ops.append(('q', [s(('leg', self.loop_edge, t))]))
        letters = [chr(65 + i) if i < 26 else chr(97 + i - 26) for i in range(len(sym))]
        self._ops = ops
        self._outsyms = [sym[('out', v)] for v, _ in self.outv]
        self._spec = ','.join(''.join(letters[i] for i in o[-1]) for o in ops) + '->' + \
            ''.join(letters[i] for i in self._outsyms)
        self._path = None

    def colour(self, Bin, q=None):
        """unsigned tensor over output indices (vertex order) for spoke-letter indices Bin (spoke order);
        q: letter parity at the tail of the loop edge (unicyclic)."""
        alg = self.alg
        arrs = []
        for o in self._ops:
            if o[0] == 'V':
                arrs.append(alg.Vf)
            elif o[0] == 'G':
                arrs.append(alg.Gf)
            elif o[0] == 'q':
                arrs.append(self.loop_proj[q])
            else:
                arrs.append(alg.Gf[:, Bin[o[1]]])
        if self._path is None:
            self._path = np.einsum_path(self._spec, *arrs, optimize='greedy')[0]
        T = np.einsum(self._spec, *arrs, optimize=self._path)
        if not np.all(np.isfinite(T)) or np.abs(T).max(initial=0) > 2.0 ** 52:
            raise OverflowError('colour beyond float64 exactness')
        return np.rint(T).astype(np.int64)

    # ---------------------------------------------------------------- signs
    def sign(self, pin, pout, O_par, spoke_pos, q=None):
        """pin: parities of the spoke letters (spoke order); pout: output parities (vertex order)."""
        key = (pin, pout, O_par, spoke_pos, q)
        if key in self._sign_cache:
            return self._sign_cache[key]
        out_par = {(v, j): pout[k] for k, (v, j) in enumerate(self.outv)}
        leg_par = solve_legs(self, pin, out_par, q)
        if leg_par is None:
            self._sign_cache[key] = 0
            return 0
        sg = self.ga.sign(leg_par, out_par, list(O_par), list(spoke_pos))
        self._sign_cache[key] = sg
        return sg

    # ---------------------------------------------------------------- kernel polynomial
    def kernel_poly(self, jets):
        """jets: tuple over spokes of multi-indices.  Returns list of (beta per output (vertex order), coef)."""
        if jets in self._kcache:
            return self._kcache[jets]
        comps = tuple(tuple(a for a in range(self.H) for _ in range(al[a])) for al in jets)
        lexp, coef = self.gk.kernel_symbolic(comps)
        H = self.H
        out = []
        idx = {nm: i for i, nm in enumerate(self.gk.lam_names)}
        for t in range(lexp.shape[0]):
            if coef[t] % po.P == 0:
                continue
            betas = []
            for v, _ in self.outv:
                betas.append(tuple(int(lexp[t, idx['l%d_%d' % (v, a)]]) for a in range(H)))
            out.append((tuple(betas), int(coef[t]) % po.P))
        self._kcache[jets] = out
        return out

    # ---------------------------------------------------------------- action
    def apply(self, poly, target=None):
        """hatQ_Gamma on a Poly.  target: optional predicate on result monomials."""
        alg = self.alg
        pi = alg.pi
        out = {}
        ns = len(self.spokes)
        if not self.outv:
            return out          # no output sources: all currents vanish, K = 0 [drv]
        for mono, coef in poly.items():
            m = len(mono)
            O_par = tuple(pi[A] for A, _ in mono)
            for sel in itertools.permutations(range(m), ns):
                Bin = [mono[k][0] for k in sel]
                jets = tuple(mono[k][1] for k in sel)
                for q in self.qvals:
                  T = self.colour(Bin, q)
                  nz = np.nonzero(T)
                  if len(nz[0]) == 0:
                    continue
                  kp = self.kernel_poly(jets)
                  if not kp:
                    continue
                  rest = tuple(mono[k] for k in range(m) if k not in sel)
                  pin = tuple(pi[B] for B in Bin)
                  for outidx in zip(*nz):
                    val = int(T[outidx])
                    pout = tuple(pi[B] for B in outidx)
                    sg = self.sign(pin, pout, O_par, sel, q)
                    if sg == 0:
                        continue
                    base = coef * val * sg % po.P * self.inv_den % po.P
                    for betas, kc in kp:
                        word = tuple((int(B), be) for B, be in zip(outidx, betas)) + rest
                        s2, mm = po.canon(word, pi)
                        if s2 == 0:
                            continue
                        if target is not None and not target(mm):
                            continue
                        po.padd(out, mm, s2 * base * kc)
        return out


def density_vanishes(gk):
    """True if the reduced density is the zero polynomial (then K = 0 for all input jets) [drv]."""
    d = gk.density_factors()
    if isinstance(d, list):
        return any(len(f.keys) == 0 for f in d)
    return False


def graph_ops(alg, n, unicyclic=False, drop_zero=True):
    import rgraphs
    gs = rgraphs.rooted_trees(alg.H, n)
    if unicyclic:
        gs = gs + rgraphs.rooted_unicyclic(alg.H, n)
    ops = []
    for k, g in enumerate(gs, 1):
        g = dict(g)
        g['name'] = 'G%d.%d%s' % (n, k, 'u' if g['kind'] == 'unicyclic' else '')
        op = GraphOp(alg, g)
        if drop_zero and density_vanishes(op.gk):
            print('  [graph %s %s: density identically zero, dropped]' % (g['name'], g['edges']), flush=True)
            continue
        ops.append(op)
    return ops
