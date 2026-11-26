"""Compositions hatQ_a hatQ_b O paired with plane-wave test functionals, modulo a prime.

Grassmann layer and colour layer: conventions of explore/six-d-Q3/verify/vops.py (S2):
  C = c + th1 g1 + th2 g2 + th1 th2 b;  delta/delta C = th1 th2 d_c + th2 d_g1 - th1 d_g2 + d_b (left derivatives,
  the last printed derivative acts first);  int d^2 th th1 th2 = 1;  hatQ_0 C = (1/2) f C C, Q_0 th = - th Q_0;
  hatQ_Gamma O = int d^2th I_Gamma(d; z) prod_v f_{legs(v)} prod_{outputs, vertex order} C(x + X_v, th)
                 delta_1 ... delta_m O,  legs(v) = (in-edges in edge order, output, out-edges in edge order).
Pairing: letter (type t, colour D, momentum p) -> sum_{i: t_i = t} eps_i A_i^D (jet -> u_i), coefficient of
eps_1 ... eps_N.  Kernels: kern.py.
"""
import itertools
import math
import time
from fractions import Fraction as Fr

import numpy as np

import kern
from kern import P, XP, xp_const, xp_add, xp_mul, xp_value, eval_symbolic

PAR = {'c': 1, 'b': 1, 'g1': 0, 'g2': 0}
THETA_OUT = {'c': (), 'g1': (1,), 'g2': (2,), 'b': (1, 2)}
THETA_DEL = {'c': ((1, 2), 1), 'g1': ((2,), 1), 'g2': ((1,), -1), 'b': ((), 1)}
Q0_RULE = {'c': [(Fr(1, 2), 'c', 'c')], 'g1': [(Fr(1), 'c', 'g1')], 'g2': [(Fr(1), 'c', 'g2')],
           'b': [(Fr(1), 'c', 'b'), (Fr(1), 'g1', 'g2')]}


# =============================================================================== Lie algebras
def so_n(N):
    idx = {}
    for i in range(N):
        for j in range(i + 1, N):
            idx[(i, j)] = len(idx)
    dim = len(idx)
    gens = []
    for (i, j) in idx:
        M = np.zeros((N, N), dtype=int)
        M[i, j], M[j, i] = 1, -1
        gens.append(M)
    f = np.zeros((dim, dim, dim), dtype=np.int64)
    for B in range(dim):
        for C in range(dim):
            Z = gens[B] @ gens[C] - gens[C] @ gens[B]
            for (i, j), A in idx.items():
                f[A, B, C] = Z[i, j]
    return dim, f


def su2():
    f = np.zeros((3, 3, 3), dtype=np.int64)
    for a, b, c in itertools.permutations(range(3)):
        f[a, b, c] = 1 if (a, b, c) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)) else -1
    return 3, f


def check_lie(f):
    assert np.all(f == -f.transpose(1, 0, 2)) and np.all(f == -f.transpose(0, 2, 1))
    jac = (np.einsum('abe,ecd->abcd', f, f) + np.einsum('ace,edb->abcd', f, f) + np.einsum('ade,ebc->abcd', f, f))
    assert np.all(jac == 0)


# =============================================================================== operators
class Op:
    def __init__(self, g, prefix):
        self.g, self.prefix = g, prefix            # g: GraphKernel or None (hatQ_0)
        self.kind = 'q0' if g is None else g.name

    @property
    def nin(self):
        return 1 if self.g is None else len(self.g.spokes)

    def outtags(self):
        if self.g is None:
            return [(self.prefix, 1), (self.prefix, 2)]
        return [(self.prefix, v) for v in self.g.outputs]


def out_assignments(nout, ths):
    """output type tuples completing the theta word ths to {1, 2}."""
    s = sorted(ths)
    res = []
    if s == [1, 2]:
        res.append(('c',) * nout)
    elif s == [1]:
        for k in range(nout):
            res.append(tuple('g2' if i == k else 'c' for i in range(nout)))
    elif s == [2]:
        for k in range(nout):
            res.append(tuple('g1' if i == k else 'c' for i in range(nout)))
    elif s == []:
        for k in range(nout):
            res.append(tuple('b' if i == k else 'c' for i in range(nout)))
        for k1 in range(nout):
            for k2 in range(nout):
                if k1 != k2:
                    res.append(tuple('g1' if i == k1 else 'g2' if i == k2 else 'c' for i in range(nout)))
    return res


def apply_graph(op, mono, binding):
    letters = list(mono)
    ths = []
    sign = 1
    for tag in reversed(binding):
        pos = next(i for i, (t, _) in enumerate(letters) if t == tag)
        ty = letters[pos][1]
        if PAR[ty]:
            par = len(ths) + sum(PAR[t] for _, t in letters[:pos])
            if par & 1:
                sign = -sign
        del letters[pos]
        t, s = THETA_DEL[ty]
        sign *= s
        ths = list(t) + ths
    if len(ths) > 2 or len(set(ths)) < len(ths):
        return []
    outs = op.outtags()
    nout = len(outs)
    res = []
    for types in out_assignments(nout, ths):
        s = sign
        word = []
        nodd = 0
        for t in types:
            th = THETA_OUT[t]
            if len(th) and (nodd & 1) and (len(th) & 1):
                s = -s
            word += th
            nodd += PAR[t]
        if (len(ths) & 1) and (nodd & 1):
            s = -s
        word += ths
        if word == [2, 1]:
            s = -s
        elif word != [1, 2]:
            raise AssertionError(word)
        res.append((s, tuple(zip(outs, types)) + tuple(letters)))
    return res


def apply_q0(op, mono, binding):
    tag = binding[0]
    pos = next(i for i, (t, _) in enumerate(mono) if t == tag)
    s = -1 if sum(PAR[ty] for _, ty in mono[:pos]) & 1 else 1
    ty = mono[pos][1]
    out = []
    for c, t1, t2 in Q0_RULE[ty]:
        new = mono[:pos] + (((op.prefix, 1), t1), ((op.prefix, 2), t2)) + mono[pos + 1:]
        out.append((s * c, new))
    return out


def apply_op(op, mono, binding):
    if op.g is None:
        return apply_q0(op, mono, binding)
    return [(Fr(s), m) for s, m in apply_graph(op, mono, binding)]


# =============================================================================== functionals and rows
class Functionals:
    """F sets of N plane-wave functionals with common types."""

    def __init__(self, types, dim, F, rng, amp_col=2, amp_mom=4, nmom=0):
        """nmom > 0: functional i of each set carries momentum class i mod nmom (nmom random momenta)."""
        self.types = list(types)
        self.N = N = len(types)
        self.F = F
        self.A = np.array([[[rng.randint(-amp_col, amp_col) for _ in range(dim)] for _ in range(N)]
                           for _ in range(F)], dtype=np.int64)
        self.Af = self.A.astype(np.float64)
        if nmom:
            base = [[[rng.randint(-amp_mom, amp_mom) for _ in range(3)] for _ in range(nmom)] for _ in range(F)]
            self.u = np.array([[base[f][i % nmom] for i in range(N)] for f in range(F)], dtype=np.int64)
        else:
            self.u = np.array([[[rng.randint(-amp_mom, amp_mom) for _ in range(3)] for _ in range(N)]
                               for _ in range(F)], dtype=np.int64)
        # momentum class of each functional (equal momenta share a class)
        self.mcls = np.zeros((F, N), dtype=np.int64)
        for f in range(F):
            seen = {}
            for i in range(N):
                self.mcls[f, i] = seen.setdefault(tuple(self.u[f, i]), len(seen))
        perms = np.array(list(itertools.permutations(range(N))), dtype=np.int64)
        self.perms = perms
        ftype = np.array([self.types.index(t) for t in self.types])  # placeholder
        tname = self.types
        self.pat_of = {}
        pats = []
        for r, s in enumerate(perms):
            pats.append(tuple(tname[i] for i in s))
        self.pats = pats
        for r, pt in enumerate(pats):
            self.pat_of.setdefault(pt, []).append(r)
        self.pat_of = {k: np.array(v, dtype=np.int64) for k, v in self.pat_of.items()}
        odd = np.array([PAR[t] for t in tname])
        sg = []
        for s in perms:
            o = [i for i in s if odd[i]]
            inv = sum(1 for a in range(len(o)) for b in range(a + 1, len(o)) if o[a] > o[b])
            sg.append(-1 if inv & 1 else 1)
        self.sg = np.array(sg, dtype=np.int64)
        self.nrows = len(perms)


# =============================================================================== colour
_paths = {}


def colour_values(ops, bindings, O_cols, slots, slot_funcs, func, f):
    """Functional-space colour tensor: returns array (F, n_0, ..., n_{N-1}) where slot s runs over slot_funcs[s]."""
    names = {}

    def nm(x):
        if x not in names:
            names[x] = len(names)
        return names[x]
    operands, subs = [], []
    for k, (op, bnd) in enumerate(zip(ops, bindings)):
        if op.g is None:
            operands.append(f)
            subs.append([nm(('col', bnd[0])), nm(('col', (op.prefix, 1))), nm(('col', (op.prefix, 2)))])
            continue
        g = op.g
        spoke_of = {e: j for j, e in enumerate(g.spokes)}
        for v in g.verts:
            idx = []
            for leg in g.legs[v]:
                if leg[0] == 'out':
                    idx.append(nm(('col', (op.prefix, leg[1]))))
                elif leg[1] in spoke_of:
                    idx.append(nm(('col', bnd[spoke_of[leg[1]]])))
                else:
                    idx.append(nm(('int', k, leg[1])))
            operands.append(f)
            subs.append(idx)
    for tag, vec in O_cols.items():
        if ('col', tag) in names:
            operands.append(vec)
            subs.append([nm(('col', tag))])
    outsub = []
    slot_ops = []
    for s, tag in enumerate(slots):
        fi = nm(('fun', s))
        slot_ops.append((s, [fi, nm(('col', tag))]))
        outsub.append(fi)
    if len(names) > 52:
        raise ValueError('too many indices')
    letters = [chr(65 + i) if i < 26 else chr(97 + i - 26) for i in range(len(names))]
    allsubs = subs + [sb for _, sb in slot_ops]
    spec = ','.join(''.join(letters[i] for i in s) for s in allsubs) + '->' + ''.join(letters[i] for i in outsub)
    fl = [o.astype(np.float64) for o in operands]
    outs = []
    for b in range(func.F):
        sl = [func.Af[b][slot_funcs[s], :] for s, _ in slot_ops]
        ops_b = fl + sl
        key = (spec, tuple(o.shape for o in ops_b))
        if key not in _paths:
            _paths[key] = np.einsum_path(spec, *ops_b, optimize='greedy')[0]
        outs.append(np.einsum(spec, *ops_b, optimize=_paths[key]))
    out = np.array(outs)
    if np.abs(out).max() > 2.0 ** 52:
        raise OverflowError('colour value beyond float64 exactness')
    return np.rint(out).astype(np.int64)


# =============================================================================== kernels per key
def _sumxp(parts, L):
    xs = [kern._as_xp(p, L) for p in parts]
    return xp_add(*xs)


def key_kernel(key, info, mom, L, mode='symbolic'):
    """Kernel of a composition hatop2 hatop1 at the slot momenta mom (list over slots of (3, L) residues).
    info: ops, smap (op2 spoke -> op1 output or None), jets (O-tag -> components), jetloc ((tag, loc), ...)."""
    op1, op2 = info['ops']
    smap, jets, jetloc = info['smap'], info['jets'], info['jetloc']
    g1, g2 = op1.g, op2.g
    const = np.ones(L, dtype=np.int64)
    extra = []
    if g1 is None:                                          # hatQ_Gamma hatQ_0
        lam2 = {v: mom[s] for s, v in enumerate(g2.outputs)}
        p2 = g2.moment_xp(lam2, L) if jetloc else None
        q_slot = {}
        rest = iter(range(len(g2.outputs), len(g2.outputs) + 2))
        for k in (1, 2):
            if k not in smap:
                q_slot[k] = next(rest)
        for tag, loc in jetloc:
            for a in jets[tag]:
                if loc[0] == 'o':
                    extra.append(p2[loc[1]][a])
                elif loc[0] == 'i':
                    parts = [p2[smap.index(k)][a] if k in smap else mom[q_slot[k]][a] for k in (1, 2)]
                    extra.append(_sumxp(parts, L))
                else:
                    const = const * mom[loc[1]][a] % P
        return g2.kernel_numeric(lam2, L, extra=extra) * const % P
    if g2 is None:                                          # hatQ_0 hatQ_Gamma
        consumed = smap[0]
        lam1 = {}
        rest = iter(range(2, 2 + len(g1.outputs)))
        for v in g1.outputs:
            lam1[v] = (mom[0] + mom[1]) % P if v == consumed else mom[next(rest)]
        p1 = g1.moment_xp(lam1, L) if jetloc else None
        for tag, loc in jetloc:
            for a in jets[tag]:
                if loc[0] == 'i':
                    extra.append(p1[loc[1]][a])
                elif loc[0] == 'o':
                    const = const * ((mom[0][a] + mom[1][a]) % P) % P
                else:
                    const = const * mom[loc[1]][a] % P
        return g1.kernel_numeric(lam1, L, extra=extra) * const % P
    # two graph operators: g1 inner, g2 outer
    lam2 = {v: mom[s] for s, v in enumerate(g2.outputs)}
    p2 = g2.moment_xp(lam2, L)
    consumed = {v: j for j, v in enumerate(smap) if v is not None}
    rest = iter(range(len(g2.outputs), len(g2.outputs) + len(g1.outputs)))
    lam1 = {}
    for v in g1.outputs:
        lam1[v] = p2[consumed[v]] if v in consumed else mom[next(rest)]
    inner_jets = any(loc[0] == 'i' for _, loc in jetloc)
    if mode == 'symbolic' and not inner_jets:
        args = []
        for nm in g1.lam_names:
            v, a = map(int, nm[1:].split('_'))
            args.append(lam1[v][a])
        lexp, coef = g1.kernel_symbolic(((),) * len(g1.spokes))
        extra = [eval_symbolic(lexp, coef, args, L)]
        for tag, loc in jetloc:
            for a in jets[tag]:
                if loc[0] == 'o':
                    extra.append(p2[loc[1]][a])
                else:
                    const = const * mom[loc[1]][a] % P
        return g2.kernel_numeric(lam2, L, extra=extra) * const % P
    shift = g2.nc
    facs = g2.density_xp(lam2, L) + kern.density_mixed(g1, lam1, L, shift)
    p1 = kern.moments_mixed(g1, lam1, L, shift) if inner_jets else None
    for tag, loc in jetloc:
        for a in jets[tag]:
            if loc[0] == 'i':
                facs.append(p1[loc[1]][a])
            elif loc[0] == 'o':
                facs.append(p2[loc[1]][a])
            else:
                const = const * mom[loc[1]][a] % P
    blocks = list(g2.blocks) + [[c + shift for c in b] for b in g1.blocks]
    val = kern.integrate_factors(facs, blocks, L)
    return val * ((g1.sign_n * g2.sign_n) % P) % P * const % P


# =============================================================================== evaluation
def evaluate(op1, op2, O, O_cols, func, f, jets=None, mode='symbolic', chunk=1024, verbose=False, stats=None):
    """<Phi_f, hatop2 hatop1 O> for every functional set f.  jets: O-tag -> tuple of derivative components
    (0, 1, 2 for d_1, d_2, d_3).  Returns residues (F,)."""
    jets = {t: tuple(v) for t, v in (jets or {}).items() if v}
    t0 = time.time()
    O = tuple(O)
    Otags = [t for t, _ in O]
    Gsum = {}
    infos = {}
    nhist = 0
    F, N = func.F, func.N
    for bind1 in itertools.permutations(Otags, op1.nin):
        terms1 = {}
        for c, m in apply_op(op1, O, bind1):
            terms1[m] = terms1.get(m, 0) + c
        terms1 = [(c, m) for m, c in terms1.items() if c != 0]
        if not terms1:
            continue
        tags1 = [t for t, _ in terms1[0][1]]
        out1 = op1.outtags()
        for bind2 in itertools.permutations(tags1, op2.nin):
            terms2 = {}
            for c1, m1 in terms1:
                for c, m in apply_op(op2, m1, bind2):
                    terms2[m] = terms2.get(m, 0) + c1 * c
            # slots
            slots = op2.outtags() + [t for t in out1 if t not in bind2] + \
                [t for t in Otags if t not in bind1 and t not in bind2]
            assert len(slots) == N
            spos = {t: i for i, t in enumerate(slots)}
            byp = {}
            for m, c in terms2.items():
                if c == 0:
                    continue
                tps = dict(m)
                pat = tuple(tps[t] for t in slots)
                if pat not in func.pat_of:
                    continue
                order = [spos[t] for t, ty in m if PAR[ty]]
                inv = sum(1 for a in range(len(order)) for b in range(a + 1, len(order)) if order[a] > order[b])
                cc = -c if inv & 1 else c
                byp[pat] = byp.get(pat, 0) + cc
            byp = {k: v for k, v in byp.items() if v != 0}
            if not byp:
                continue
            nhist += 1
            # kernel key
            jetloc = []
            for t in Otags:
                if t in jets:
                    if t in bind1:
                        jetloc.append((t, ('i', bind1.index(t))))
                    elif t in bind2:
                        jetloc.append((t, ('o', bind2.index(t))))
                    else:
                        jetloc.append((t, ('s', slots.index(t))))
            jetloc = tuple(jetloc)
            smap = tuple(b[1] if b[0] == op1.prefix else None for b in bind2)
            if op1.g is None:
                key = ('Q0in', op2.kind, smap if jetloc else None, jetloc)
            elif op2.g is None:
                key = ('Q0out', op1.kind, smap, jetloc)
            else:
                key = ('GG', op1.kind, op2.kind, smap, jetloc)
            infos[key] = dict(ops=(op1, op2), smap=smap, jets=jets, jetloc=jetloc)
            # colour on the allowed functionals of each slot
            slot_types = [set(pt[s] for pt in byp) for s in range(N)]
            slot_funcs = [np.array([i for i in range(N) if func.types[i] in slot_types[s]], dtype=np.int64)
                          for s in range(N)]
            G = colour_values((op1, op2), (bind1, bind2), O_cols, slots, slot_funcs, func, f)
            acc = Gsum.setdefault(key, np.zeros((F, func.nrows), dtype=np.int64))
            for pat, cc in byp.items():
                rows = func.pat_of[pat]
                S = func.perms[rows]
                # position of S[:, s] inside slot_funcs[s]
                idx = [np.searchsorted(slot_funcs[s], S[:, s]) for s in range(N)]
                vals = G[(slice(None),) + tuple(idx)] % P
                acc[:, rows] = (acc[:, rows] + vals * kern.modq(cc)) % P
    tcol = time.time() - t0
    total = np.zeros(F, dtype=np.int64)
    nkeys = 0
    for key, acc in Gsum.items():
        acc = acc * (func.sg[None, :] % P) % P
        fi, ri = np.nonzero(acc)
        if len(fi) == 0:
            continue
        nkeys += 1
        # rows with equal slot momenta share the kernel
        S_all = func.perms[ri]
        code = fi.copy()
        for s in range(N):
            code = code * N + func.mcls[fi, S_all[:, s]]
        ucode, first, inv = np.unique(code, return_index=True, return_inverse=True)
        gs = np.zeros(len(ucode), dtype=np.int64)
        np.add.at(gs, inv, acc[fi, ri])
        gs %= P
        uf, ur = fi[first], ri[first]
        for c0 in range(0, len(ucode), chunk):
            fc, rc = uf[c0:c0 + chunk], ur[c0:c0 + chunk]
            L = len(fc)
            S = func.perms[rc]
            mom = [np.array([func.u[fc, S[:, s], a] % P for a in range(3)], dtype=np.int64) for s in range(N)]
            K = key_kernel(key, infos[key], mom, L, mode=mode)
            contrib = K * gs[c0:c0 + chunk] % P
            np.add.at(total, fc, contrib)
            total %= P
    if stats is not None:
        stats.update(hist=nhist, keys=nkeys, tcol=tcol, ttot=time.time() - t0)
    if verbose:
        print('    %s <- %s: histories %d, kernel keys %d, colour %.1fs, total %.1fs'
              % (op2.kind, op1.kind, nhist, nkeys, tcol, time.time() - t0), flush=True)
    return total % P


# =============================================================================== single operator
def evaluate1(op, O, O_cols, func, f, jets=None, chunk=1024):
    """<Phi_f, hatop O> for a single graph operator (jets on O letters allowed).  Returns residues (F,)."""
    jets = {t: tuple(v) for t, v in (jets or {}).items() if v}
    g = op.g
    O = tuple(O)
    Otags = [t for t, _ in O]
    F, N = func.F, func.N
    Gsum = {}
    for bind in itertools.permutations(Otags, op.nin):
        terms = {}
        for c, m in apply_op(op, O, bind):
            terms[m] = terms.get(m, 0) + c
        slots = op.outtags() + [t for t in Otags if t not in bind]
        assert len(slots) == N
        spos = {t: i for i, t in enumerate(slots)}
        byp = {}
        for m, c in terms.items():
            if c == 0:
                continue
            tps = dict(m)
            pat = tuple(tps[t] for t in slots)
            if pat not in func.pat_of:
                continue
            order = [spos[t] for t, ty in m if PAR[ty]]
            inv = sum(1 for a in range(len(order)) for b in range(a + 1, len(order)) if order[a] > order[b])
            byp[pat] = byp.get(pat, 0) + (-c if inv & 1 else c)
        byp = {k: v for k, v in byp.items() if v != 0}
        if not byp:
            continue
        jetloc = tuple((t, ('i', bind.index(t)) if t in bind else ('s', slots.index(t))) for t in Otags if t in jets)
        slot_types = [set(pt[s] for pt in byp) for s in range(N)]
        slot_funcs = [np.array([i for i in range(N) if func.types[i] in slot_types[s]], dtype=np.int64)
                      for s in range(N)]
        G = colour_values((op,), (bind,), O_cols, slots, slot_funcs, func, f)
        acc = Gsum.setdefault(jetloc, np.zeros((F, func.nrows), dtype=np.int64))
        for pat, cc in byp.items():
            rows = func.pat_of[pat]
            S = func.perms[rows]
            idx = [np.searchsorted(slot_funcs[s], S[:, s]) for s in range(N)]
            acc[:, rows] = (acc[:, rows] + G[(slice(None),) + tuple(idx)] % P * kern.modq(cc)) % P
    total = np.zeros(F, dtype=np.int64)
    for jetloc, acc in Gsum.items():
        acc = acc * (func.sg[None, :] % P) % P
        fi, ri = np.nonzero(acc)
        for c0 in range(0, len(fi), chunk):
            fc, rc = fi[c0:c0 + chunk], ri[c0:c0 + chunk]
            L = len(fc)
            S = func.perms[rc]
            mom = [np.array([func.u[fc, S[:, s], a] % P for a in range(3)], dtype=np.int64) for s in range(N)]
            lam = {v: mom[s] for s, v in enumerate(g.outputs)}
            extra, const = [], np.ones(L, dtype=np.int64)
            p = g.moment_xp(lam, L) if jetloc else None
            for t, loc in jetloc:
                for a in jets[t]:
                    if loc[0] == 'i':
                        extra.append(p[loc[1]][a])
                    else:
                        const = const * mom[loc[1]][a] % P
            K = g.kernel_numeric(lam, L, extra=extra) * const % P
            np.add.at(total, fc, K * acc[fc, rc] % P)
            total %= P
    return total % P
