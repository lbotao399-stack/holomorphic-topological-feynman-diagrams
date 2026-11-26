"""T5: classes of H(C_quot) and H(V), w = 1, G in [Gmin, Gmax], all torus weights (S_n-orbit representatives).
For each H(C) class: an exact cocycle representative z (reduced modulo Q C^{N-1}), and its beta_0-coefficient
c(z) := coefficient of beta_0 (beta_0 moved to the left).  Q z = 0 forces Q c = 0 (w = 0 complex, slot (t, N-1));
the class [c] in H^{N-1}(C_quot^{w=0})_t is the image under C^{w=1} -> C^{w=1}/C_F ~ beta_0 C^{w=0}.
For each H(V) class: the explicit lemma map  z = eta_K c,  eta_k = iota_k + h_k (h_k beta_a = -gamma_{k,a}),
c running over H(C) representatives at (t+e_K, N+|K|); checked: Q z in sum_k D_k C (V-cocycle) and the z's are
a basis of H(V)_{t,N} (rank of [Q C^{N-1} + I^N + span z] - rank[Q C^{N-1} + I^N] = dim H(V)).
usage: python3 t5_central.py n Gmin Gmax  -> t5_n{n}.log"""
import sys, itertools, time
from cohom import Space, Cohom, all_weights
from jets import Echelon, vec
from reps import reps_C, full_reduce

n, Gmin, Gmax = (int(x) for x in sys.argv[1:4])
supp_only = len(sys.argv) > 4 and sys.argv[4] == 'supp'
out = open('t5_n%d.log' % n, 'w')


def in_supp(t):
    for r in range(n + 1):
        for K in itertools.combinations(range(n), r):
            u = [t[j] + (j in K) for j in range(n)]
            if len(set(u)) == 1:
                return True
    return False


def log(s):
    print(s, flush=True)
    out.write(s + '\n')
    out.flush()


S1 = Space(n, 1, Gmax + n + 2)
S0 = Space(n, 0, Gmax + n + 3)
C1 = Cohom(S1)
C0 = Cohom(S0)
A = S1.A
A0 = S0.A


def show(Aa, p, maxt=24):
    return Aa.pstr(p, maxt) + '   [%d terms]' % len(p)


def beta0_coeff(z):
    b0 = A.L('b', 0, (0,) * n)
    c = {}
    for m, v in z.items():
        if b0 in m:
            j = m.index(b0)
            sgn = -1 if sum(1 for x in m[:j] if A.par[x]) % 2 else 1
            rest = m[:j] + m[j + 1:]
            c[rest] = c.get(rest, 0) + sgn * v
    return {m: v for m, v in c.items() if v}


def to_A0(p):
    r = {}
    for m, v in p.items():
        mm = tuple(sorted(A0.L(*A.info[x]) for x in m))
        # sign of re-sorting odd letters
        ids = [A0.L(*A.info[x]) for x in m]
        odd = [x for x in ids if A0.par[x]]
        inv = sum(1 for i in range(len(odd)) for j in range(i + 1, len(odd)) if odd[i] > odd[j])
        r[mm] = r.get(mm, 0) + (-1 if inv % 2 else 1) * v
    return r


def class_in_w0(c, t, N):
    """is c exact in the w=0 complex at (t, N)? returns (closed, exact, coordinates w.r.t. H reps)."""
    if not c:
        return True, True
    c0 = to_A0(c)
    q = A0.Q(c0)
    closed = not q
    src = S0.slot(t).get(N - 1, [])
    mid = S0.slot(t).get(N, [])
    idx = {m: j for j, m in enumerate(mid)}
    E = Echelon()
    for m in src:
        E.insert(vec(A0.Qmono(m), idx))
    r0 = E.rank
    E.insert(vec(c0, idx))
    return closed, E.rank == r0


def eta(k):
    def f(lid):
        typ, i, a = A.info[lid]
        if typ == 'X' and i == k and sum(a) == 0:
            return {(): 1}
        if typ == 'b':
            return {(A.L('g', k, a),): -1}
        return {}
    return f


def canon(t):
    return tuple(sorted(t, reverse=True))


log('T5 (t5_central.py) n=%d w=1 G=%d..%d%s' % (n, Gmin, Gmax, ' (H(C) on t=(k,..,k) only, H(V) on t+e_K=(k,..,k) only)' if supp_only else ''))
HCreps = {}
for G in range(Gmin - n, Gmax + 1):
    for t in all_weights(S1, G):
        if t != canon(t) and G >= Gmin:
            pass
    reps_t = sorted({canon(t) for t in all_weights(S1, G)}, reverse=True)
    if supp_only:
        reps_t = [t for t in reps_t if len(set(t)) == 1]
    for t in reps_t:
        h = C1.HC(t)
        for N in sorted(h):
            d, hh = h[N]
            if hh:
                zs = reps_C(C1, t, N)
                assert len(zs) == hh
                HCreps[(t, N)] = zs
                if G < Gmin:
                    continue
                log('H(C_quot) G=%d t=%s N=%d: dim %d (slot dim %d)' % (G, t, N, hh, d))
                for z in zs:
                    log('   z = ' + show(A, z))
                    c = beta0_coeff(z)
                    if c:
                        closed, exact = class_in_w0(c, t, N - 1)
                        log('   beta_0-coefficient c = ' + show(A, c, 12) + ' ; Qc=0: %s ; c exact in w=0: %s' % (closed, exact))
                    else:
                        log('   beta_0-coefficient c = 0')
    log('# G=%d H(C) done' % G)


def hc_reps(t, N):
    tc = canon(t)
    if (tc, N) not in HCreps:
        return []
    if tc == t:
        return HCreps[(tc, N)]
    # permute indices: t = sigma(tc); transport the representatives by the index permutation
    perm = None
    for p in itertools.permutations(range(n)):
        if tuple(tc[p[j]] for j in range(n)) == t:
            perm = p
            break
    # letter map: index i -> position j with p[j] = i
    inv = {p_: j for j, p_ in enumerate(perm)}
    res = []
    for z in HCreps[(tc, N)]:
        p = {}
        for m, v in z.items():
            acc = {(): v}
            for x in m:
                typ, i, a = A.info[x]
                ni = inv[i] if typ in ('X', 'g') else i
                na = tuple(a[perm[j]] for j in range(n))
                acc = A.pmul(acc, A.red_letter(A.L(typ, ni, na)) if True else None)
            p = A.padd(p, acc)
        res.append(p)
    return res


log('H(V) classes via the explicit lemma map z = eta_K c:')
for G in range(Gmin, Gmax + 1):
    reps_t = sorted({canon(t) for t in all_weights(S1, G)}, reverse=True)
    if supp_only:
        reps_t = [t for t in reps_t if in_supp(t)]
    for t in reps_t:
        hv = C1.HV(t)
        for N in sorted(hv):
            d, dV, hh = hv[N]
            zs = []
            labels = []
            for r in range(n + 1):
                for K in itertools.combinations(range(n), r):
                    tt = tuple(t[j] + (j in K) for j in range(n))
                    for c in hc_reps(tt, N + r):
                        z = c
                        for k in reversed(K):
                            z = A.apply_der(z, eta(k), True)
                        zs.append(z)
                        labels.append((K, tt, N + r))
            if not hh and not zs:
                continue
            mid = C1.basis(t, N)
            tgt = C1.basis(t, N + 1)
            midx = {m: j for j, m in enumerate(mid)}
            tidx = {m: j for j, m in enumerate(tgt)}
            EI = Echelon()
            for k in range(n):
                for m in C1.basis(C1.shiftw(t, k), N + 1):
                    EI.insert(vec(A.Dq({m: 1}, k), tidx))
            cocyc = all(not EI.reduce(vec(A.Q(z), tidx)) for z in zs) if tgt else True
            E = Echelon()
            for k in range(n):
                for m in C1.basis(C1.shiftw(t, k), N):
                    E.insert(vec(A.Dq({m: 1}, k), midx))
            for m in C1.basis(t, N - 1):
                E.insert(vec(A.Qmono(m), midx))
            r0 = E.rank
            Ecopy = Echelon()
            Ecopy.piv = {k: dict(v) for k, v in E.piv.items()}
            for z in zs:
                E.insert(vec(z, midx))
            log('H(V) G=%d t=%s N=%d: dim %d ; lemma images %d, V-cocycles: %s, independent rank %d -> basis: %s' % (
                G, t, N, hh, len(zs), cocyc, E.rank - r0, (E.rank - r0 == hh == len(zs)) and cocyc))
            for z, lab in zip(zs, labels):
                zr = full_reduce(Ecopy, vec(z, midx))
                zr = {mid[j]: v for j, v in zr.items()}
                log('   K=%s from H(C)%s N=%d :  z = %s' % (tuple(k + 1 for k in lab[0]), lab[1], lab[2], show(A, zr, 16)))
    log('# G=%d H(V) done' % G)
