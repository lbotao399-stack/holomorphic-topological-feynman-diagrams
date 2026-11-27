"""[fork A copy of run_lifts.py: flavour U(1)^Nf charges, U(Nf)-singlet option, printed representatives]
Scan of Q0-cohomology sectors of gauge-invariant, c_0-independent local operators and of d_1 = [Q_1] (mod P).
usage: python3 run_lifts.py THEORY GAUGE MATTER KMAX DMAX [--charges i,j] [--kappa1 k] [--d2]
Sector (q, d, k): q block charges, d in N^H with |d| <= DMAX, k <= KMAX letters.  Writes one line per sector with
dim H(Q0) and, for sources, rank d_1 into the target (q + dq, d + 1, k + H - 2).
"""
import argparse
import itertools
import json
import sys
import time
from fractions import Fraction as Fr

import numpy as np

import polyops as po
import theories
from wick import Algebra
from explicit import graph_ops
import lifts as lf

ap = argparse.ArgumentParser()
ap.add_argument('theory')
ap.add_argument('gauge')
ap.add_argument('matter')
ap.add_argument('kmax', type=int)
ap.add_argument('dmax', type=int)
ap.add_argument('--kmin', type=int, default=1)
ap.add_argument('--maxdim', type=int, default=20000, help='skip sectors with more monomials')
ap.add_argument('--only_d', default='', help='restrict to derivative vectors, e.g. 1.1;2.2')
ap.add_argument('--d2', action='store_true')
ap.add_argument('--only_q', default='')
ap.add_argument('--flavour', type=int, default=0, help='N_f: add U(1)^N_f flavour charges (U = T*(fund*N_f))')
ap.add_argument('--flavsinglet', action='store_true', help='only U(N_f)-singlet sectors/invariants')
ap.add_argument('--show', type=int, default=0, help='print up to this many lifted source representatives per sector')
args = ap.parse_args()

matter = None if args.matter == '-' else args.matter
H, sl = theories.build(args.theory, args.gauge, matter)
alg = Algebra(sl, H)
g = theories.gauge(args.gauge)
# Cartan and positive roots of g (weight basis)
mats = g.mats
cartan = [a for a, m in enumerate(mats) if not (m - np.diag(np.diag(m))).any()]
gidx = sl.blocks['g']
wts = []
for a in range(g.dim):
    w = []
    for h in cartan:
        row = g.f[h, a]
        nz = np.nonzero(row)[0]
        w.append(int(row[a]) if len(nz) else 0)
    wts.append(w)
rho = [1000 ** (len(cartan) - i) for i in range(len(cartan))]
pos = [a for a in range(g.dim) if a not in cartan and sum(x * y for x, y in zip(wts[a], rho)) > 0]
th = lf.LiftTheory(alg, cartan=[gidx[a] for a in cartan], posroots=[gidx[a] for a in pos])
NF = args.flavour
nflav_base = len(th.charges)
flav_raise = []          # list of dicts letter index -> list of (letter index, coeff)
if NF:
    uidx = sl.blocks['u']
    usidx = sl.blocks['u*'] if 'u*' in sl.blocks else sl.blocks['u^v']
    dimR = len(uidx) // 2
    comp = dimR // NF
    def upart(p):
        # (part, copy, component): part 0 = R, 1 = R^*
        if p < dimR:
            return 0, p // comp, p % comp
        return 1, (p - dimR) // comp, (p - dimR) % comp
    def upos(part, j, s_):
        return (0 if part == 0 else dimR) + j * comp + s_
    flw = {}
    for p, A in enumerate(uidx):
        part, j, s_ = upart(p)
        w = [0] * NF
        w[j] = -1 if part == 0 else 1          # letter phi^{u}: dual (U^*) weights
        flw[A] = w
    for p, A in enumerate(usidx):
        part, j, s_ = upart(p)
        w = [0] * NF
        w[j] = 1 if part == 0 else -1          # letter phi^{u*}: U weights
        flw[A] = w
    for A in range(alg.N):
        th.wt[A] = tuple(th.wt[A]) + tuple(flw.get(A, [0] * NF))
    th.charges = list(th.charges) + [None] * NF
    # U(N_f) raising operators E_ij (i < j) on letters: u* letters as U vectors, u letters as -rho^T (U^*)
    for i in range(NF):
        for j in range(i + 1, NF):
            act = {}
            for p, A in enumerate(usidx):
                part, jj, s_ = upart(p)
                if part == 0 and jj == j:
                    act.setdefault(A, []).append((usidx[upos(0, i, s_)], 1))
                if part == 1 and jj == i:
                    act.setdefault(A, []).append((usidx[upos(1, j, s_)], -1))
            for p, A in enumerate(uidx):
                part, jj, s_ = upart(p)
                # -rho(E_ij)^T = -rho(E_ji): R_i -> -R_j ; R*_j -> +R*_i
                if part == 0 and jj == i:
                    act.setdefault(A, []).append((uidx[upos(0, j, s_)], -1))
                if part == 1 and jj == j:
                    act.setdefault(A, []).append((uidx[upos(1, i, s_)], 1))
            flav_raise.append(act)


def flav_action(act, mono):
    out = {}
    for k, (A, al) in enumerate(mono):
        for B, c in act.get(A, []):
            word = mono[:k] + ((B, al),) + mono[k + 1:]
            sg, m = po.canon(word, th.pi)
            if sg:
                po.padd(out, m, sg * c)
    return out


def my_invariants(monos):
    cols = []
    keymap = {}
    for m in monos:
        img = {}
        for r, a in enumerate(th.posroots):
            for mm, v in lf.g_action(th, a, m).items():
                img[('g', r, mm)] = (img.get(('g', r, mm), 0) + v) % po.P
        if args.flavsinglet:
            for r, act in enumerate(flav_raise):
                for mm, v in flav_action(act, m).items():
                    img[('f', r, mm)] = (img.get(('f', r, mm), 0) + v) % po.P
        cols.append({keymap.setdefault(key, len(keymap)): v for key, v in img.items() if v})
    return lf.kernel(cols, len(monos))


print('theory %s gauge %s matter %s H=%d N=%d; blocks %s' % (args.theory, args.gauge, matter, H, alg.N,
      {b: len(v) for b, v in th.part.items()}), flush=True)
print('charges (per block):', th.charges, flush=True)
Q0 = alg.q0
G1 = graph_ops(alg, 1)
assert len(G1) == 1
Q1op = G1[0]
G2 = graph_ops(alg, 2) if args.d2 else []
kap2 = {2: [po.modq(Fr(1, 2))] * len(G2)} if H == 2 else None   # H=2 fans: kappa_n = 1/|Aut| [chk nilp_explicit]
kap1 = po.modq(Fr(1, Q1op.aut))


def Q1(poly):
    return po.pscale(Q1op.apply(poly), kap1)


def Q2(poly):
    out = {}
    for g, k in zip(G2, kap2[2]):
        for m, v in g.apply(poly).items():
            po.padd(out, m, k * v)
    return out


# charge shift of Q_1 from a random nonzero image
rng = __import__('random').Random(1)
dq = None
for _ in range(2000):
    lets = [(rng.randrange(alg.N), tuple(rng.randint(0, 1) for _ in range(H))) for _ in range(rng.randint(2, 3))]
    s, m = po.canon(lets, alg.pi)
    if not s:
        continue
    img = Q1({m: 1})
    if img:
        m2 = next(iter(img))
        q1 = tuple(sum(th.wt[A][i] for A, _ in m) for i in range(len(th.charges)))
        q2 = tuple(sum(th.wt[A][i] for A, _ in m2) for i in range(len(th.charges)))
        dq = tuple(b - a for a, b in zip(q1, q2))
        break
print('Q1 charge shift dq =', dq, ' letter shift', H - 2, flush=True)

cache = {}


def fmt_poly(p, maxterms=16):
    return po.show(p, sl.names, H, maxterms=maxterms)



def sector(q, d, k):
    key = (q, d, k)
    if key in cache:
        return cache[key]
    if k <= 0:
        cache[key] = None
        return None
    monos = lf.sector_monomials(th, q, d, k)
    if len(monos) > args.maxdim:
        cache[key] = 'big'
        return 'big'
    inv = my_invariants(monos) if monos else []
    polys = [lf.vec_to_poly(v, monos) for v in inv]
    imgs = [lf.drop_c0(th, Q0(p)) for p in polys]
    idx = {}
    cols = []
    for im in imgs:
        cols.append({idx.setdefault(m, len(idx)): v for m, v in im.items()})
    ker = lf.kernel(cols, len(polys))
    Z = []
    for kv in ker:
        z = {}
        for j, c in kv.items():
            for m, v in polys[j].items():
                po.padd(z, m, c * v)
        Z.append(z)
    res = dict(monos=monos, inv=polys, Q0img=imgs, Z=Z)
    cache[key] = res
    return res


def boundaries(q, d, k):
    """B = Q0(Inv_{k-1}) in sector (q, d, k)."""
    prev = sector(q, d, k - 1)
    if prev is None:
        return []
    if prev == 'big':
        return 'big'
    return [im for im in prev['Q0img'] if im]


def rank_polys(polys):
    idx = {}
    rows = []
    for p in polys:
        rows.append({idx.setdefault(m, len(idx)): v for m, v in p.items()})
    return lf.rank(rows)


def all_charges(d, k):
    monos = lf.sector_monomials_anyq(th, d, k) if hasattr(lf, 'sector_monomials_anyq') else None
    return monos


dvecs = []
for tot in range(args.dmax + 1):
    for c in itertools.combinations_with_replacement(range(H), tot):
        dv = [0] * H
        for a in c:
            dv[a] += 1
        dvecs.append(tuple(dv))
if args.only_d:
    dvecs = [tuple(int(x) for x in s.split('.')) for s in args.only_d.split(';')]
T0 = time.time()
found = []
for k in range(args.kmin, args.kmax + 1):
    for d in dvecs:
        # enumerate all charges occurring
        qs = set()
        for m in lf.sector_monomials(th, None, d, k):
            qs.add(tuple(sum(th.wt[A][i] for A, _ in m) for i in range(len(th.charges))))
        if args.only_q:
            qs = {tuple(int(x) for x in s.split('.')) for s in args.only_q.split(';')} & qs
        for q in sorted(qs):
            if args.flavsinglet and any(q[nflav_base:]):
                continue
            S = sector(q, d, k)
            if S is None or S == 'big':
                print('sector q=%s d=%s k=%d: %s' % (q, d, k, S), flush=True)
                continue
            B = boundaries(q, d, k)
            if B == 'big':
                print('sector q=%s d=%s k=%d: boundary sector too big' % (q, d, k), flush=True)
                continue
            dimZ, rkB = len(S['Z']), rank_polys(B)
            dimH = dimZ - rkB
            if dimH == 0:
                continue
            # d_1
            qT = tuple(a + b for a, b in zip(q, dq))
            dT = tuple(x + 1 for x in d)
            kT = k + H - 2
            imgs = [Q1(z) for z in S['Z']]
            c0 = sum(1 for im in imgs for mm in im if any(th.is_c0(l) for l in mm))
            imgs = [lf.drop_c0(th, im) for im in imgs]
            BT = boundaries(qT, dT, kT)
            if BT == 'big':
                print('sector q=%s d=%s k=%d: dim H = %d; target boundaries too big' % (q, d, k, dimH), flush=True)
                continue
            rB = rank_polys(BT)
            r1 = rank_polys(BT + imgs) - rB
            line = 'sector q=%s d=%s k=%d: monomials %d, invariants %d, dim Z %d, rank B %d, dim H %d; d_1 -> (q=%s, d=%s, k=%d): rank %d%s' % (
                q, d, k, len(S['monos']), len(S['inv']), dimZ, rkB, dimH, qT, dT, kT, r1,
                ('  [c0 terms in Q1 Z: %d]' % c0) if c0 else '')
            print(line, flush=True)
            if r1:
                found.append((q, d, k, dimH, r1))
                if args.show:
                    # representatives: Z elements whose images are independent modulo B_T
                    basis = list(BT)
                    shown = 0
                    for z, im in zip(S['Z'], imgs):
                        if rank_polys(basis + [im]) > rank_polys(basis):
                            basis.append(im)
                            print('    source cocycle:', ' + '.join(fmt_poly(z).split('\n')), flush=True)
                            print('    Q_1 image     :', ' + '.join(fmt_poly(im, 12).split('\n')), flush=True)
                            shown += 1
                            if shown >= args.show:
                                break
            if args.d2 and dimH - r1 > 0:
                # kernel of d_1: combinations with Q1 z = - Q0 o1
                prevT = sector(qT, dT, kT - 1)
                binv = [] if prevT is None else prevT['inv']
                bimg = [] if prevT is None else prevT['Q0img']
                cols_polys = imgs + bimg
                idx2 = {}
                cols = [{idx2.setdefault(m, len(idx2)): v for m, v in p.items()} for p in cols_polys]
                kv = lf.kernel(cols, len(cols_polys))
                nz = len(imgs)
                R2s = []
                for vec in kv:
                    z, o1 = {}, {}
                    for j, c in vec.items():
                        if j < nz:
                            for m, v in S['Z'][j].items():
                                po.padd(z, m, c * v)
                        else:
                            for m, v in binv[j - nz].items():
                                po.padd(o1, m, c * v)
                    if not z:
                        continue
                    R2 = po.psum(lf.drop_c0(th, Q2(z)), lf.drop_c0(th, Q1(o1)))
                    R2s.append(R2)
                q2 = tuple(a + 2 * b for a, b in zip(q, dq))
                d2 = tuple(x + 2 for x in d)
                k2 = k + 2 * (H - 3) + 1
                B2 = boundaries(q2, d2, k2)
                # Q1 Z of the preimage sector (q + dq, d + 1, k2 - (H - 3) - 1)
                pre = sector(tuple(a + b for a, b in zip(q, dq)), dT, k2 - (H - 3) - 1)
                QZ = [] if pre in (None, 'big') else [lf.drop_c0(th, Q1(zz)) for zz in pre['Z']]
                if B2 == 'big' or pre == 'big':
                    print('    d_2: target too big', flush=True)
                    continue
                den = rank_polys(B2 + QZ)
                r2 = rank_polys(B2 + QZ + R2s) - den
                print('    d_2 on ker d_1 (%d chains) -> (q=%s, d=%s, k=%d): rank %d' % (len(R2s), q2, d2, k2, r2),
                      flush=True)
                if r2:
                    found.append(('d2', q, d, k, r2))
print('nonzero d_1:', found, flush=True)
print('total time %.0fs' % (time.time() - T0), flush=True)
