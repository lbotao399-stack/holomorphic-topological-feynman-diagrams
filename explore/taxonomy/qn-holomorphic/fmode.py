"""Plane-wave functional evaluation <Phi, hatop2 hatop1 O> modulo P (generalizes 6d-nilpotency/comp.py to any H and
any quadratic Lie superalgebra).

Types = blocks of g_tot (basis subsets of fixed parity).  O = product of smeared letters (tag, block, jet) with colour
vectors; intermediate letters carry explicit block types (enumerated, filtered by a random structural probe); colour =
float64 einsum of V (vertices), G (edges), block projectors, O colours and functional colours (exactness checked);
Koszul signs of each Wick step from wick.GraphAlg; pairing: letter at slot s -> sum_i eps_i A_i u_i^alpha, coefficient of
eps_1..eps_N (odd eps anticommute).  Kernels: kern_h.GraphKernel (I = (-1)^{n(n+1)/2} K).
"""
import itertools
import time

import numpy as np

import kern_h as kern
from kern_h import XP, xp_const, xp_add, eval_symbolic
from explicit import GraphOp, _solve_gf2, solve_legs
from wick import GraphAlg

P = kern.P


# =============================================================================== operators
class FOp:
    def __init__(self, alg, gop, prefix):
        """gop: explicit.GraphOp or None (hatQ0)."""
        self.alg, self.gop, self.prefix = alg, gop, prefix
        if gop is None:
            self.ga = GraphAlg(1, [(0, 1)])
            self.kind = 'Q0'
            self.nG = 1
            self.outv = [(0, 0), (0, 1)]
            self.spokes = [0]
            self.edges = [(0, 1)]
            self.internal = []
            self.nb = 1
            self.loop_edge = None
            self.qvals = (None,)
        else:
            self.ga = gop.ga
            self.kind = gop.name
            self.nG = gop.nG
            self.outv = gop.outv
            self.spokes = gop.spokes
            self.edges = gop.edges
            self.internal = gop.internal
            self.nb = gop.nb
            self.loop_edge = gop.loop_edge
            self.qvals = gop.qvals
            self.loop_proj = gop.loop_proj
        self._sign_cache = {}
        self._allowed = {}

    @property
    def g(self):
        return None if self.gop is None else self.gop.gk

    @property
    def nin(self):
        return len(self.spokes)

    def outtags(self):
        return [(self.prefix, k) for k in range(len(self.outv))]

    # -------------------------------------------------------------- Koszul sign (as GraphOp.sign, two outputs ok)
    def sign(self, pin, pout, word_par, spoke_pos, q=None):
        key = (pin, pout, word_par, spoke_pos, q)
        if key in self._sign_cache:
            return self._sign_cache[key]
        out_par = {ov: pout[k] for k, ov in enumerate(self.outv)}
        leg_par = solve_legs(self, pin, out_par, q)
        if leg_par is None:
            self._sign_cache[key] = 0
            return 0
        sg = self.ga.sign(leg_par, out_par, list(word_par), list(spoke_pos))
        self._sign_cache[key] = sg
        return sg


# =============================================================================== context
class Context:
    def __init__(self, alg, blocks):
        """blocks: dict name -> list of basis indices (a partition, constant parity)."""
        self.alg = alg
        self.blocks = blocks
        self.bnames = list(blocks)
        self.bpar = {}
        N = alg.N
        cover = np.zeros(N, dtype=int)
        for b, ids in blocks.items():
            ps = {alg.pi[i] for i in ids}
            assert len(ps) == 1, 'block %s mixes parities' % b
            self.bpar[b] = ps.pop()
            cover[ids] += 1
        assert (cover == 1).all(), 'blocks must partition the basis'
        self.proj = {b: np.zeros(N) for b in blocks}
        for b, ids in blocks.items():
            self.proj[b][ids] = 1.0
        self.rng = np.random.RandomState(12345)
        self.Vf, self.Gf = alg.Vf, alg.Gf
        self.bd = BlockData(self)

    # -------------------------------------------------------------- allowed output blocks of one graph op
    def allowed_outputs(self, op, in_blocks):
        """output block tuples allowed by the block sparsity of V and G (superset of the nonzero ones)."""
        key = (op.kind, in_blocks)
        if key in op._allowed:
            return op._allowed[key]
        bd = self.bd
        ga = op.ga
        nout = len(op.outv)
        res = set()
        # enumerate blocks of spoke tails, internal legs and outputs vertex by vertex
        legkeys = []
        dom = {}
        for j, e in enumerate(op.spokes):
            t, _ = op.edges[e]
            dom[('leg', e, t)] = bd.Gdual_l[in_blocks[j]]
        for e in op.internal:
            t, h = op.edges[e]
            dom[('leg', e, t)] = list(self.bnames)
            dom[('leg', e, h)] = None          # determined by the tail through G
        for i in range(nout):
            dom[('out', i)] = list(self.bnames)

        def vlegs(v):
            out = []
            for leg in ga.legs[v]:
                out.append(('out', op.outv.index((v, leg[2]))) if leg[0] == 'out' else ('leg', leg[1], v))
            return out
        order = []
        for v in range(op.nb):
            for lk in vlegs(v):
                if lk not in order:
                    order.append(lk)
        heads = {('leg', e, op.edges[e][1]): ('leg', e, op.edges[e][0]) for e in op.internal}

        def check(a):
            for v in range(op.nb):
                bl = [a.get(lk) for lk in vlegs(v)]
                if None not in bl and tuple(bl) not in bd.Vnz:
                    return False
            for hk, tk in heads.items():
                if hk in a and tk in a and (a[tk], a[hk]) not in bd.Gnz:
                    return False
            return True

        def rec(i, a):
            if i == len(order):
                res.add(tuple(a[('out', k)] for k in range(nout)))
                return
            lk = order[i]
            cands = dom[lk]
            if cands is None:
                tk = heads[lk]
                cands = bd.Gdual[a[tk]] if tk in a else list(self.bnames)
            for b in cands:
                a[lk] = b
                if check(a):
                    rec(i + 1, a)
                del a[lk]
        rec(0, {})
        out = sorted(res)
        op._allowed[key] = out
        return out


def single_tensor(op, ctx, invecs, q=None):
    """colour tensor of one operator with input colour vectors invecs (spoke order) contracted through G."""
    sym = {}

    def s(k):
        if k not in sym:
            sym[k] = len(sym)
        return sym[k]
    operands, subs = [], []
    ga = op.ga
    for v in range(op.nb):
        idx = []
        for leg in ga.legs[v]:
            idx.append(s(('out', v, leg[2])) if leg[0] == 'out' else s(('leg', leg[1], v)))
        operands.append(ctx.Vf)
        subs.append(idx)
    for e in op.internal:
        t, h = op.edges[e]
        operands.append(ctx.Gf)
        subs.append([s(('leg', e, t)), s(('leg', e, h))])
    for j, e in enumerate(op.spokes):
        t, _ = op.edges[e]
        operands.append(np.einsum('ab,b->a', ctx.Gf, invecs[j]))
        subs.append([s(('leg', e, t))])
    if op.loop_edge is not None:
        t, _ = op.edges[op.loop_edge]
        operands.append(op.loop_proj[q])
        subs.append([s(('leg', op.loop_edge, t))])
    outs = [s(('out', v, k)) for v, k in op.outv]
    letters = [chr(65 + i) if i < 26 else chr(97 + i - 26) for i in range(len(sym))]
    spec = ','.join(''.join(letters[i] for i in sb) for sb in subs) + '->' + ''.join(letters[i] for i in outs)
    return np.einsum(spec, *operands, optimize='greedy')


# =============================================================================== type-level application
def apply_types(op, ctx, word, binding):
    """word: tuple of (tag, block); binding: tags bound to the spokes (spoke order).
    Returns list of (sign, new_word, out_blocks) with new_word = outputs (op order) + unused letters."""
    alg = ctx.alg
    tags = [t for t, _ in word]
    pos = [tags.index(t) for t in binding]
    in_blocks = tuple(word[p][1] for p in pos)
    rest = tuple(word[k] for k in range(len(word)) if k not in pos)
    word_par = tuple(ctx.bpar[b] for _, b in word)
    pin = tuple(ctx.bpar[b] for b in in_blocks)
    res = []
    for combo in ctx.allowed_outputs(op, in_blocks):
        pout = tuple(ctx.bpar[b] for b in combo)
        for q in op.qvals:
            sg = op.sign(pin, pout, word_par, tuple(pos), q)
            if sg == 0:
                continue
            new = tuple(zip(op.outtags(), combo)) + rest
            res.append((sg, new, combo, q))
    return res


# =============================================================================== functionals
class Functionals:
    def __init__(self, ctx, types, F, rng, amp_col=2, amp_mom=4):
        self.types = list(types)
        self.N = N = len(types)
        self.F = F
        H = ctx.alg.H
        dim = ctx.alg.N
        self.A = np.zeros((F, N, dim), dtype=np.int64)
        for f in range(F):
            for i, t in enumerate(types):
                for a in ctx.blocks[t]:
                    self.A[f, i, a] = rng.randint(-amp_col, amp_col)
        self.Af = self.A.astype(np.float64)
        self.u = np.array([[[rng.randint(-amp_mom, amp_mom) for _ in range(H)] for _ in range(N)]
                           for _ in range(F)], dtype=np.int64)
        self.mcls = np.zeros((F, N), dtype=np.int64)
        for f in range(F):
            seen = {}
            for i in range(N):
                self.mcls[f, i] = seen.setdefault(tuple(self.u[f, i]), len(seen))
        perms = np.array(list(itertools.permutations(range(N))), dtype=np.int64)
        self.perms = perms
        self.pat_of = {}
        for r, s in enumerate(perms):
            self.pat_of.setdefault(tuple(types[i] for i in s), []).append(r)
        self.pat_of = {k: np.array(v, dtype=np.int64) for k, v in self.pat_of.items()}
        odd = np.array([ctx.bpar[t] for t in types])
        sg = []
        for s in perms:
            o = [i for i in s if odd[i]]
            inv = sum(1 for a in range(len(o)) for b in range(a + 1, len(o)) if o[a] > o[b])
            sg.append(-1 if inv & 1 else 1)
        self.sg = np.array(sg, dtype=np.int64)
        self.nrows = len(perms)


# =============================================================================== colour of a two-step history
_paths = {}


class BlockData:
    """block slices of V and G and their nonzero patterns."""

    def __init__(self, ctx):
        self.ctx = ctx
        alg = ctx.alg
        bn = ctx.bnames
        self.ids = {b: np.array(ctx.blocks[b]) for b in bn}
        self.Vnz = set()
        self.Gnz = set()
        self.Vs, self.Gs = {}, {}
        for b1 in bn:
            for b2 in bn:
                sub = alg.Gf[np.ix_(self.ids[b1], self.ids[b2])]
                if np.any(sub):
                    self.Gnz.add((b1, b2))
                    self.Gs[(b1, b2)] = sub
                for b3 in bn:
                    sub = alg.Vf[np.ix_(self.ids[b1], self.ids[b2], self.ids[b3])]
                    if np.any(sub):
                        self.Vnz.add((b1, b2, b3))
                        self.Vs[(b1, b2, b3)] = sub
        self.Gdual = {b: [b2 for b2 in bn if (b, b2) in self.Gnz] for b in bn}
        self.Gdual_l = {b: [b1 for b1 in bn if (b1, b) in self.Gnz] for b in bn}
        self._assign = {}


def leg_assignments(bd, op, k, ext, q):
    """blocks of all bulk legs of op (index k in the history) given ext: ('col', tag) -> block of the letter at
    each output and at each bound spoke letter.  Returns list of dicts leg-key -> block."""
    key = (op.kind, k, tuple(sorted(ext.items())), q)
    if key in bd._assign:
        return bd._assign[key]
    ctx = bd.ctx
    ga = op.ga
    legs = []                       # leg keys in a fixed order: spoke tails, then internal (tail, head)
    choices = {}
    for j, e in enumerate(op.spokes):
        t, _ = op.edges[e]
        lk = ('leg', k, e, t)
        legs.append(lk)
        choices[lk] = bd.Gdual_l[ext[('spoke', j)]]
    for e in op.internal:
        t, h = op.edges[e]
        legs.append(('leg', k, e, t))
        legs.append(('leg', k, e, h))
    out = []

    def vblock(v, a):
        bl = []
        for leg in ga.legs[v]:
            if leg[0] == 'out':
                bl.append(ext[('out', op.outv.index((v, leg[2])))])
            else:
                bl.append(a.get(('leg', k, leg[1], v)))
        return bl

    def ok_partial(a):
        for v in range(op.nb):
            bl = vblock(v, a)
            if None in bl:
                continue
            if tuple(bl) not in bd.Vnz:
                return False
        return True

    internal = list(op.internal)

    def rec(i, a):
        if i < len(op.spokes):
            lk = legs[i]
            for b in choices[lk]:
                a[lk] = b
                if ok_partial(a):
                    rec(i + 1, a)
                del a[lk]
            return
        j = i - len(op.spokes)
        if j == len(internal):
            out.append(dict(a))
            return
        e = internal[j]
        t, h = op.edges[e]
        lt, lh = ('leg', k, e, t), ('leg', k, e, h)
        for bt in ctx.bnames:
            if op.loop_edge == e and ctx.bpar[bt] != q:
                continue
            for bh in bd.Gdual[bt]:
                a[lt], a[lh] = bt, bh
                if ok_partial(a):
                    rec(i + 1, a)
                del a[lt], a[lh]
    rec(0, {})
    bd._assign[key] = out
    return out


def colour_values(ctx, ops, bindings, wtypes, O_cols, slots, slot_funcs, func, qs=(None, None)):
    """ops = (op1 inner, op2 outer); bindings = (bind1, bind2); wtypes: tag -> block of every letter that occurs.
    Block-sparse contraction, batched over functional sets.  Returns array (F, n_0, ..., n_{N-1})."""
    bd = ctx.bd
    per_op = []
    for k, (op, bnd) in enumerate(zip(ops, bindings)):
        ext = {('out', i): wtypes[(op.prefix, i)] for i in range(len(op.outv))}
        for j in range(op.nin):
            ext[('spoke', j)] = wtypes[bnd[j]]
        per_op.append(leg_assignments(bd, op, k, ext, qs[k]))
    names = {}

    def nm(x):
        if x not in names:
            names[x] = len(names)
        return names[x]
    # index structure (fixed), operands depend on the assignment
    struct = []                     # (kind, data, subs)
    for k, (op, bnd) in enumerate(zip(ops, bindings)):
        ga = op.ga
        sp = {e: j for j, e in enumerate(op.spokes)}
        for v in range(op.nb):
            idx, lk = [], []
            for leg in ga.legs[v]:
                if leg[0] == 'out':
                    tg = ('col', (op.prefix, op.outv.index((v, leg[2]))))
                    idx.append(nm(tg))
                    lk.append(('fixed', wtypes[tg[1]]))
                else:
                    idx.append(nm(('leg', k, leg[1], v)))
                    lk.append(('var', ('leg', k, leg[1], v)))
            struct.append(('V', k, lk, idx))
        for e in op.internal:
            t, h = op.edges[e]
            struct.append(('G', k, [('var', ('leg', k, e, t)), ('var', ('leg', k, e, h))],
                           [nm(('leg', k, e, t)), nm(('leg', k, e, h))]))
        for e in op.spokes:
            t, _ = op.edges[e]
            tg = bnd[sp[e]]
            struct.append(('G', k, [('var', ('leg', k, e, t)), ('fixed', wtypes[tg])],
                           [nm(('leg', k, e, t)), nm(('col', tg))]))
    ocols = []
    for tag, vec in O_cols.items():
        if ('col', tag) in names:
            b = wtypes[tag]
            ocols.append((vec[bd.ids[b]], [nm(('col', tag))]))
    slot_ops = []
    outsub = []
    for s, tag in enumerate(slots):
        fi = nm(('fun', s))
        b = wtypes[tag]
        slot_ops.append((s, b, [fi, nm(('col', tag))]))
        outsub.append(fi)
    if len(names) > 52:
        raise ValueError('too many indices')
    letters = [chr(65 + i) if i < 26 else chr(97 + i - 26) for i in range(len(names))]
    allsubs = [st[3] for st in struct] + [sb for _, sb in ocols] + [sb for _, _, sb in slot_ops]
    spec = ','.join(''.join(letters[i] for i in sb) for sb in allsubs) + '->' + ''.join(letters[i] for i in outsub)
    total = None
    for combo in itertools.product(*per_op):
        a = {}
        for d in combo:
            a.update(d)
        arrs = []
        for kind, k, lk, idx in struct:
            bl = tuple(x[1] if x[0] == 'fixed' else a[x[1]] for x in lk)
            arrs.append(bd.Vs[bl] if kind == 'V' else bd.Gs[bl])
        arrs += [v for v, _ in ocols]
        res = []
        for fb in range(func.F):
            sl = [func.Af[fb][slot_funcs[s], :][:, bd.ids[b]] for s, b, _ in slot_ops]
            ab = arrs + sl
            key = (spec, tuple(o.shape for o in ab))
            if key not in _paths:
                _paths[key] = np.einsum_path(spec, *ab, optimize='greedy')[0]
            res.append(np.einsum(spec, *ab, optimize=_paths[key]))
        r = np.array(res)
        total = r if total is None else total + r
    if total is None:
        shape = [func.F] + [len(slot_funcs[s]) for s in range(len(slots))]
        return np.zeros(shape, dtype=np.int64)
    if not np.all(np.isfinite(total)) or np.abs(total).max(initial=0) > 2.0 ** 52:
        raise OverflowError('colour value beyond float64 exactness')
    return np.rint(total).astype(np.int64)


# =============================================================================== kernels per key
def _sumxp(parts, L):
    return xp_add(*[kern._as_xp(p, L) for p in parts])


def key_kernel(info, mom, L, H, mode='symbolic'):
    """as comp.key_kernel, H components; mom: list over slots of (H, L) residues."""
    op1, op2 = info['ops']
    smap, jets, jetloc = info['smap'], info['jets'], info['jetloc']
    g1, g2 = op1.g, op2.g
    const = np.ones(L, dtype=np.int64)
    extra = []
    if op1.gop is None and op2.gop is None:                 # Q0 Q0
        return const
    if op1.gop is None:                                      # hatQ_Gamma hatQ_0
        lam2 = {v: mom[s] for s, (v, _) in enumerate(op2.outv)}
        p2 = g2.moment_xp(lam2, L) if jetloc else None
        q_slot = {}
        rest = iter(range(len(op2.outv), len(op2.outv) + 2))
        for k in (0, 1):
            if k not in smap:
                q_slot[k] = next(rest)
        for tag, loc in jetloc:
            for a in jets[tag]:
                if loc[0] == 'o':
                    extra.append(p2[loc[1]][a])
                elif loc[0] == 'i':
                    parts = [p2[smap.index(k)][a] if k in smap else mom[q_slot[k]][a] for k in (0, 1)]
                    extra.append(_sumxp(parts, L))
                else:
                    const = const * mom[loc[1]][a] % P
        return g2.kernel_numeric(lam2, L, extra=extra) * const % P
    if op2.gop is None:                                      # hatQ_0 hatQ_Gamma
        consumed = smap[0]
        lam1 = {}
        rest = iter(range(2, 2 + len(op1.outv)))
        for k, (v, _) in enumerate(op1.outv):
            lam1[v] = (mom[0] + mom[1]) % P if k == consumed else mom[next(rest)]
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
    lam2 = {v: mom[s] for s, (v, _) in enumerate(op2.outv)}
    p2 = g2.moment_xp(lam2, L)
    consumed = {k: j for j, k in enumerate(smap) if k is not None}
    rest = iter(range(len(op2.outv), len(op2.outv) + len(op1.outv)))
    lam1 = {}
    for k, (v, _) in enumerate(op1.outv):
        lam1[v] = p2[consumed[k]] if k in consumed else mom[next(rest)]
    inner_jets = any(loc[0] == 'i' for _, loc in jetloc)
    if mode == 'symbolic' and not inner_jets and g1.red.tight_plan:
        args = []
        for nmv in g1.lam_names:
            v, a = map(int, nmv[1:].split('_'))
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
def evaluate(ctx, op1, op2, O, O_cols, func, jets=None, mode='symbolic', chunk=1024, stats=None):
    """<Phi_f, hatop2 hatop1 O> for every functional set.  O: tuple of (tag, block); jets: tag -> components.
    Includes 1/Gden^{|I|} and the 1/2 of hatQ0.  Returns residues (F,)."""
    alg = ctx.alg
    H = alg.H
    jets = {t: tuple(v) for t, v in (jets or {}).items() if v}
    t0 = time.time()
    O = tuple(O)
    Otags = [t for t, _ in O]
    Gsum, infos = {}, {}
    nhist = 0
    F, N = func.F, func.N
    norm = pow(alg.Gden, op1.nG + op2.nG, P)
    norm = pow(norm, P - 2, P)
    if op1.gop is None:
        norm = norm * pow(2, P - 2, P) % P
    if op2.gop is None:
        norm = norm * pow(2, P - 2, P) % P
    for bind1 in itertools.permutations(Otags, op1.nin):
        terms1 = apply_types(op1, ctx, O, bind1)
        if not terms1:
            continue
        for s1, w1, combo1, q1 in terms1:
            tags1 = [t for t, _ in w1]
            for bind2 in itertools.permutations(tags1, op2.nin):
                terms2 = apply_types(op2, ctx, w1, bind2)
                if not terms2:
                    continue
                for s2, w2, combo2, q2 in terms2:
                    slots = [t for t, _ in w2]
                    pat = tuple(b for _, b in w2)
                    if pat not in func.pat_of:
                        continue
                    nhist += 1
                    wtypes = dict(O)
                    wtypes.update(dict(zip(op1.outtags(), combo1)))
                    wtypes.update(dict(zip(op2.outtags(), combo2)))
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
                    if op1.gop is None:
                        key = ('Q0in', op2.kind, smap if jetloc else None, jetloc)
                    elif op2.gop is None:
                        key = ('Q0out', op1.kind, smap, jetloc)
                    else:
                        key = ('GG', op1.kind, op2.kind, smap, jetloc)
                    infos[key] = dict(ops=(op1, op2), smap=smap, jets=jets, jetloc=jetloc)
                    slot_funcs = [np.array([i for i in range(N) if func.types[i] == pat[s]], dtype=np.int64)
                                  for s in range(N)]
                    Gv = colour_values(ctx, (op1, op2), (bind1, bind2), wtypes, O_cols, slots, slot_funcs, func,
                                       qs=(q1, q2))
                    acc = Gsum.setdefault(key, np.zeros((F, func.nrows), dtype=np.int64))
                    rows = func.pat_of[pat]
                    S = func.perms[rows]
                    idx = [np.searchsorted(slot_funcs[s], S[:, s]) for s in range(N)]
                    vals = Gv[(slice(None),) + tuple(idx)] % P
                    acc[:, rows] = (acc[:, rows] + vals * (s1 * s2 % P)) % P
    tcol = time.time() - t0
    total = np.zeros(F, dtype=np.int64)
    nkeys = 0
    for key, acc in Gsum.items():
        acc = acc * (func.sg[None, :] % P) % P
        fi, ri = np.nonzero(acc)
        if len(fi) == 0:
            continue
        nkeys += 1
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
            mom = [np.array([func.u[fc, S[:, s], a] % P for a in range(H)], dtype=np.int64) for s in range(N)]
            K = key_kernel(infos[key], mom, L, H, mode=mode)
            contrib = K * gs[c0:c0 + chunk] % P
            np.add.at(total, fc, contrib)
            total %= P
    total = total * norm % P
    if stats is not None:
        stats.update(hist=nhist, keys=nkeys, tcol=tcol, ttot=time.time() - t0)
    return total % P
