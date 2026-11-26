"""T3: nu-elimination.  n = 2, w in {0,1}.
C_full^{<=A} := C_full / (arity > A)  (a quotient complex since Q raises arity by >= 0);
C_quot^{<=A} likewise.  Compare dim H(C_full^{<=A}), dim H(C_quot^{<=A}) and dim H(C_quot) per slot (w, t, N).
Also checks that pi: C_full -> C_quot (nu -> 0, X -> [X]) is a dg map on generators.
Writes t3.log."""
import sys, time, itertools
from collections import defaultdict
from jets import Alg, multi_indices, Echelon, vec, TYPES

out = open('t3_w1.log' if (len(sys.argv) > 1 and sys.argv[1] == 'w1') else 't3.log', 'w')


def log(s):
    print(s, flush=True)
    out.write(s + '\n')
    out.flush()


def letters_upto(A, maxord, quot):
    n = A.n
    L = []
    for m in range(maxord + 1):
        for a in multi_indices(n, m):
            for typ in TYPES:
                if quot and typ == 'n':
                    continue
                for i in (range(n) if typ in ('X', 'g') else [0]):
                    lid = A.L(typ, i, a)
                    if quot and not A.is_basis(lid):
                        continue
                    L.append(lid)
    return L


def enumerate_slot(A, w, t, Amax, quot):
    """all monomials of C_full (quot=False) or C_quot (quot=True) with weight (w,t) and arity <= Amax.
    returns dict N -> list."""
    n = A.n
    G0 = -w - sum(t)
    maxord = G0 + n + 2
    L = letters_upto(A, maxord, quot)
    L.sort(key=lambda x: -A.Gw[x])
    neg_after = [0] * (len(L) + 1)
    for j in range(len(L) - 1, -1, -1):
        neg_after[j] = neg_after[j + 1] + (1 if A.Gw[L[j]] < 0 else 0)
    res = defaultdict(list)
    isw = [A.info[x][0] in ('g', 'b') for x in L]

    def rec(j, G, ar, wc, tt, chosen):
        if j == len(L):
            if G == G0 and wc == w and tuple(tt) == t:
                m = tuple(sorted(chosen))
                res[A.Ndeg(m)].append(m)
            return
        # pruning on G
        if G - neg_after[j] > G0:
            return
        gj = A.Gw[L[j]]
        if gj >= 0 and G + (Amax - ar) * gj < G0 and gj > 0:
            # remaining letters have G <= gj
            return
        if G + (Amax - ar) * max(gj, 0) < G0:
            return
        x = L[j]
        rec(j + 1, G, ar, wc, tt, chosen)
        maxmult = 1 if A.par[x] else Amax - ar
        if isw[j]:
            maxmult = min(maxmult, w - wc)
        tx = A.tw[x]
        for k in range(1, maxmult + 1):
            if ar + k > Amax:
                break
            nt = [tt[i] + k * tx[i] for i in range(n)]
            rec(j + 1, G + k * gj, ar + k, wc + (k if isw[j] else 0), nt, chosen + [x] * k)

    rec(0, 0, 0, 0, [0] * n, [])
    return {N: sorted(v) for N, v in res.items()}


def H_trunc(A, sl, Amax):
    if not sl:
        return {}
    ranks = {}
    for N in sl:
        tgt = sl.get(N + 1, [])
        if not tgt:
            ranks[N] = 0
            continue
        idx = {m: j for j, m in enumerate(tgt)}
        E = Echelon()
        for m in sl[N]:
            q = {mm: c for mm, c in A.Qmono(m).items() if len(mm) <= Amax}
            E.insert(vec(q, idx))
        ranks[N] = E.rank
    return {N: len(sl[N]) - ranks[N] - ranks.get(N - 1, 0) for N in range(min(sl), max(sl) + 1) if N in sl}


def check_pi():
    """pi(Q phi) == Q_quot(pi phi) on generators, n = 2, 3, |a| <= 3."""
    res = []
    for n in (2, 3):
        Af = Alg(n, quot=False)
        Aq = Alg(n, quot=True)
        ok = True
        for m in range(4):
            for a in multi_indices(n, m):
                for typ in TYPES:
                    for i in (range(n) if typ in ('X', 'g') else [0]):
                        lf = Af.L(typ, i, a)
                        q = Af.Qletter(lf)
                        # transport to Aq and reduce
                        qq = {}
                        for mm, c in q.items():
                            p = {(): c}
                            for x in mm:
                                p = Aq.pmul(p, Aq.red_letter(Aq.L(*Af.info[x])))
                            qq = Aq.padd(qq, p)
                        lq = Aq.L(typ, i, a)
                        rhs = {}
                        for (x,), c in Aq.red_letter(lq).items():
                            rhs = Aq.padd(rhs, Aq.Qletter(x), c)
                        if Aq.padd(qq, rhs, -1):
                            ok = False
        res.append((n, ok))
    return res


if __name__ == '__main__':
    log('T3 (t3_nu.py)')
    log('pi: C_full -> C_quot (nu_a -> 0, X -> [X]) commutes with Q on generators |a|<=3: %s' % check_pi())
    n = 2
    Af = Alg(n, quot=False)
    Aq = Alg(n, quot=True)
    if len(sys.argv) > 1 and sys.argv[1] == 'w1':
        tests = [(1, (1, 1)), (1, (0, 0)), (1, (0, -1)), (1, (-1, -1))]
    else:
        tests = [(0, (1, 1)), (0, (0, 0)), (0, (0, -1)), (0, (-1, -1)), (0, (-2, -2)),
                 (1, (1, 1)), (1, (0, 0)), (1, (0, -1)), (1, (-1, -1))]
    CAP = 120000
    for w, t in tests:
        G = -w - sum(t)
        t0 = time.time()
        full = enumerate_slot(Aq, w, t, 99, True)
        amax = max((len(m) for v in full.values() for m in v), default=0)
        Hq = H_trunc(Aq, full, 99)
        log('slot w=%d G=%d t=%s : C_quot dims %s, max arity in C_quot slot a_max=%d' % (w, G, t, {N: len(v) for N, v in sorted(full.items())}, amax))
        log('   H(C_quot) = %s' % {N: h for N, h in sorted(Hq.items()) if h})
        for Amax in range(1, amax + 3):
            slf = enumerate_slot(Af, w, t, Amax, False)
            if sum(len(v) for v in slf.values()) > CAP:
                log('   A=%2d  dim C_full^{<=A}=%6d > cap %d: skipped' % (Amax, sum(len(v) for v in slf.values()), CAP))
                break
            slq = enumerate_slot(Aq, w, t, Amax, True)
            Hf = H_trunc(Af, slf, Amax)
            Hqt = H_trunc(Aq, slq, Amax)
            nzf = {N: h for N, h in sorted(Hf.items()) if h}
            nzq = {N: h for N, h in sorted(Hqt.items()) if h}
            log('   A=%2d  dim C_full^{<=A}=%6d  H(C_full^{<=A})=%s  H(C_quot^{<=A})=%s  equal:%s  equals H(C_quot):%s' % (
                Amax, sum(len(v) for v in slf.values()), nzf, nzq, nzf == nzq, nzf == {N: h for N, h in Hq.items() if h}))
        log('   (%.1fs)' % (time.time() - t0))
