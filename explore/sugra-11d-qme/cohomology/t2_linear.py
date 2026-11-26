"""T2: cohomology of the arity-preserving part Q_lin of Q on single letters of C_full (n = 2, 3).
Q_lin is read off from the arity-1 terms of Q(letter) (computed, not typed in).
Slots (G, t, N) of single letters; N=0: nu, gamma; N=1: X, beta.
Writes t2.log."""
import sys, itertools
from collections import defaultdict
from jets import Alg, multi_indices, Echelon, TYPES

out = open('t2.log', 'w')


def log(s):
    print(s)
    out.write(s + '\n')
    out.flush()


def run(n, M):
    A = Alg(n, quot=False, sX=-1)
    letters = []
    for m in range(M + 2):
        for a in multi_indices(n, m):
            for typ in TYPES:
                for i in (range(n) if typ in ('X', 'g') else [0]):
                    letters.append(A.L(typ, i, a))
    slots = defaultdict(lambda: {0: [], 1: []})
    for x in letters:
        typ, i, a = A.info[x]
        slots[(A.Gw[x], A.tw[x])][A.par[x]].append(x)
    # Q_lin
    qlin = {}
    lin_types = set()
    for x in letters:
        if A.par[x] == 0 and sum(A.info[x][2]) <= M:
            q = {m[0]: c for m, c in A.Qletter(x).items() if len(m) == 1}
            qlin[x] = q
            for y in q:
                lin_types.add((A.info[x][0], A.info[y][0]))
    for x in letters:
        if A.par[x] == 1 and sum(A.info[x][2]) <= M + 1:
            q = {m[0]: c for m, c in A.Qletter(x).items() if len(m) == 1}
            assert not q
    log('n=%d: arity-1 terms of Q occur only for (source, target) types %s; odd letters have none' % (n, sorted(lin_types)))
    # explicit form check
    ok_form = True
    for x, q in qlin.items():
        typ, i, a = A.info[x]
        if typ == 'n':
            exp = {}
            for j in range(n):
                b = list(a); b[j] += 1
                exp[A.L('X', j, tuple(b))] = 1
        else:
            b = list(a); b[i] += 1
            exp = {A.L('b', 0, tuple(b)): 1}
        if exp != q:
            ok_form = False
    log('  Q_lin nu_a = sum_i X^i_{a+e_i}, Q_lin gamma_{i,a} = beta_{a+e_i}: %s' % ok_form)
    # cohomology per slot, G <= M-1 so that all letters of the slot are present
    table = []
    allok = True
    for (G, t), d in sorted(slots.items()):
        if G > M - 1:
            continue
        ev, od = d[0], d[1]
        oidx = {y: j for j, y in enumerate(od)}
        E = Echelon()
        for x in ev:
            E.insert({oidx[y]: c for y, c in qlin[x].items()})
        r = E.rank
        H0 = len(ev) - r
        H1 = len(od) - r
        # predictions
        nX = sum(1 for y in od if A.info[y][0] == 'X')
        nB = sum(1 for y in od if A.info[y][0] == 'b')
        nN = sum(1 for x in ev if A.info[x][0] == 'n')
        nG = sum(1 for x in ev if A.info[x][0] == 'g')
        # rank of nu -> X (should be nN), gamma -> beta
        En = Echelon()
        for x in ev:
            if A.info[x][0] == 'n':
                En.insert({oidx[y]: c for y, c in qlin[x].items()})
        Eg = Echelon()
        for x in ev:
            if A.info[x][0] == 'g':
                Eg.insert({oidx[y]: c for y, c in qlin[x].items()})
        # p ^ gamma combos: gamma_{j,a+e_i} - gamma_{i,a+e_j}
        gidx = {x: j for j, x in enumerate(ev)}
        Ef = Echelon()
        fin_ker = True
        mG = None
        for x in ev:
            if A.info[x][0] == 'g':
                mG = sum(A.info[x][2])
        if mG is not None and mG >= 1:
            for a in multi_indices(n, mG - 1):
                for i in range(n):
                    for j in range(i + 1, n):
                        ai = list(a); ai[i] += 1
                        aj = list(a); aj[j] += 1
                        lj = A.L('g', j, tuple(ai))
                        li = A.L('g', i, tuple(aj))
                        if lj in gidx and li in gidx:
                            v = {gidx[lj]: 1, gidx[li]: -1}
                            img = {}
                            for x, c in ((lj, 1), (li, -1)):
                                for y, cc in qlin[x].items():
                                    img[y] = img.get(y, 0) + c * cc
                            if any(img.values()):
                                fin_ker = False
                            Ef.insert(v)
        kerg = nG - Eg.rank
        pred_H1 = (nX - nN) + (1 if (G, t) == (-1, (0,) * n) else 0)
        pred_H0 = kerg
        ok = (En.rank == nN) and (H1 == pred_H1) and (Ef.rank == kerg) and fin_ker and (nB - Eg.rank == (1 if (G, t) == (-1, (0,) * n) else 0))
        allok &= ok
        table.append((G, t, H0, H1, nX - nN, kerg, ok))
    log('  slots (G, t) with G <= %d: %d ; all predictions hold: %s' % (M - 1, len(table), allok))
    log('  [H^0 = ker(gamma->beta) = span{gamma_{j,a+e_i}-gamma_{i,a+e_j}}, H^1 = coker(nu->X) (+ C beta_0 at (G,t)=(-1,0)), nu->X injective]')
    # per (|a|, t) tables
    log('  X mod p.X : dim per (|a|, t)   [X^i_a has G=|a|-1, t=e_i-a]')
    byord = defaultdict(list)
    for (G, t, H0, H1, xm, kg, ok) in table:
        if xm:
            byord[G + 1].append((t, xm))
    for m in sorted(byord):
        tot = sum(v for _, v in byord[m])
        exp = n * comb(m + n - 1, n - 1) - (comb(m + n - 2, n - 1) if m >= 1 else 0)
        log('    |a|=%d total %d (formula n*C(m+n-1,n-1)-C(m+n-2,n-1) = %d): %s' % (m, tot, exp, ' '.join('%s:%d' % (t, v) for t, v in byord[m])))
    log('  Z(gamma) = functions of p^gamma : dim per (|a|, t)   [gamma_{i,a} has G=|a|, t=-e_i-a]')
    byord = defaultdict(list)
    for (G, t, H0, H1, xm, kg, ok) in table:
        if kg:
            byord[G].append((t, kg))
    for m in sorted(byord):
        tot = sum(v for _, v in byord[m])
        exp = n * comb(m + n - 1, n - 1) - comb(m + n, n - 1)
        log('    |a|=%d total %d (formula n*C(m+n-1,n-1)-C(m+n,n-1) = %d): %s' % (m, tot, exp, ' '.join('%s:%d' % (t, v) for t, v in byord[m])))
    nb = [(G, t, H1 - xm) for (G, t, H0, H1, xm, kg, ok) in table if H1 - xm]
    log('  beta part of H^1: %s' % nb)


from math import comb
if __name__ == '__main__':
    log('T2 (t2_linear.py)')
    run(2, 5)
    run(3, 4)
