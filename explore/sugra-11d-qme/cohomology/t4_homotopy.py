"""T4 part 1: exact operator identities on C_full and C_quot.
  [Q, D_k] = 0;  D_k preserves the div ideal;
  iota_k : X^i_a -> delta^i_k delta_{a,0} (odd derivation),  h_k^{(s)} : beta_a -> s gamma_{k,a} (odd derivation);
  [Q, iota_k + h_k^{(s)}] (graded commutator = anticommutator) compared with +-D_k for s = +1, -1;
  [eta_k, eta_l] = 0, [eta_k, D_l] = 0, eta_k preserves the ideal;
  Cartan homotopies for gl_n: xi = z^b d_a  (a != b) and xi = sum_j c_j z^j d_j:
     iota_xi(X^i_a) = d^a xi^i (0),  h_xi(beta_a) = -d^a(xi^k gamma_k)(0);
     [Q, iota_xi + h_xi] = +-L_xi, L_xi the gl_n action on jets (X upper, gamma lower, derivatives lower).
Checked on all generators with |a| <= 3 (n = 2, 3) and on random monomials (arity <= 5) of C_quot.
Writes t4_homotopy.log."""
import random, itertools
from jets import Alg, multi_indices, unit, TYPES

out = None


def log(s):
    print(s, flush=True)
    if out is not None:
        out.write(s + '\n')
        out.flush()


def gens(A, maxord, quot):
    n = A.n
    for m in range(maxord + 1):
        for a in multi_indices(n, m):
            for typ in TYPES:
                if quot and typ == 'n':
                    continue
                for i in (range(n) if typ in ('X', 'g') else [0]):
                    lid = A.L(typ, i, a)
                    if quot and not A.is_basis(lid):
                        continue
                    yield lid


def make_ops(A):
    n = A.n

    def iota(k, sgn=1):
        def f(lid):
            typ, i, a = A.info[lid]
            if typ == 'X' and i == k and sum(a) == 0:
                return {(): sgn}
            return {}
        return f

    def hk(k, s):
        def f(lid):
            typ, i, a = A.info[lid]
            if typ == 'b':
                return {(A.L('g', k, a),): s}
            return {}
        return f

    def plus(f, g):
        def h(lid):
            return A.padd(f(lid), g(lid))
        return h
    return iota, hk, plus


def anti(A, f, p, odd_f=True):
    """[Q, F] p for an odd derivation F given on letters: Q F p + F Q p."""
    Fp = A.apply_der(p, f, True)
    return A.padd(A.Q(Fp), A.apply_der(A.Q(p), f, True))


def check_translation(n, quot, maxord=3, nrand=150, seed=1):
    A = Alg(n, quot=quot, sX=-1)
    iota, hk, plus = make_ops(A)
    G = list(gens(A, maxord, quot))
    rnd = random.Random(seed)
    # random monomials from basis letters with |a| <= 2
    low = [x for x in G if sum(A.info[x][2]) <= 2]
    samples = []
    while len(samples) < nrand:
        k = rnd.randint(1, 5)
        m = rnd.sample(low, k)
        p = {(): 1}
        for x in m:
            p = A.pmul(p, {(x,): 1})
        if p:
            samples.append(p)
    tests = [{(x,): 1} for x in G] + samples
    res = {}
    for k in range(n):
        Dk = lambda p, k=k: A.reduce(A.D(p, k))
        # [Q, D_k] = 0
        res.setdefault('[Q,D_k]=0', True)
        for p in tests:
            if A.padd(A.Q(Dk(p)), Dk(A.Q(p)), -1):
                res['[Q,D_k]=0'] = False
        for s in (1, -1):
            eta = plus(iota(k), hk(k, s))
            okp = okm = True
            for p in tests:
                c = anti(A, eta, p)
                d = Dk(p)
                if A.padd(c, d, -1):
                    okp = False
                if A.padd(c, d, 1):
                    okm = False
            res.setdefault('[Q,iota_k+h_k]=+D_k, h_k beta=%+d gamma_k' % s, True)
            res.setdefault('[Q,iota_k+h_k]=-D_k, h_k beta=%+d gamma_k' % s, True)
            res['[Q,iota_k+h_k]=+D_k, h_k beta=%+d gamma_k' % s] &= okp
            res['[Q,iota_k+h_k]=-D_k, h_k beta=%+d gamma_k' % s] &= okm
    # commutators of eta's (with s = -1) and with D
    ok1 = ok2 = True
    for k in range(n):
        for l in range(n):
            ek = plus(iota(k), hk(k, -1))
            el = plus(iota(l), hk(l, -1))
            for p in tests:
                a = A.padd(A.apply_der(A.apply_der(p, el, True), ek, True), A.apply_der(A.apply_der(p, ek, True), el, True))
                if a:
                    ok1 = False
                b = A.padd(A.apply_der(A.reduce(A.D(p, l)), ek, True), A.reduce(A.D(A.apply_der(p, ek, True), l)), -1)
                if b:
                    ok2 = False
    res['[eta_k,eta_l]=0'] = ok1
    res['[eta_k,D_l]=0'] = ok2
    if quot:
        # eta_k and D_k preserve the ideal: evaluate on div_a in the full algebra then reduce
        okI = okD = True
        for m in range(maxord + 1):
            for a in multi_indices(n, m):
                div = {}
                for i in range(n):
                    b = list(a); b[i] += 1
                    div[(A.L('X', i, tuple(b)),)] = 1
                for k in range(n):
                    # iota_k(div_a) (full letters) : only X^k_0 contributes, never in div
                    v = sum(1 for (x,) in div if A.info[x][2] == (0,) * n)
                    if v:
                        okI = False
                    if A.reduce(A.D(div, k)):
                        okD = False
        res['eta_k(ideal) in ideal'] = okI
        res['D_k(ideal) in ideal'] = okD
    return res, len(G), len(samples)


def gl_action(A, a_, b_):
    """natural action of E_{ab} on letters, as an even derivation (letter -> poly), full letters."""
    n = A.n

    def f(lid):
        typ, i, al = A.info[lid]
        r = {}
        if typ == 'X' and i == b_:
            r = A.padd(r, {(A.L('X', a_, al),): 1})
        if typ == 'g' and i == a_:
            r = A.padd(r, {(A.L('g', b_, al),): -1})
        if al[a_] >= 1:
            be = list(al); be[a_] -= 1; be[b_] += 1
            r = A.padd(r, {(A.L(typ, i, tuple(be)),): -al[a_]})
        return A.reduce(r)
    return f


def check_cartan(n, maxord=3, nrand=150, seed=2):
    A = Alg(n, quot=True, sX=-1)
    G = list(gens(A, maxord, True))
    rnd = random.Random(seed)
    low = [x for x in G if sum(A.info[x][2]) <= 2]
    samples = []
    while len(samples) < nrand:
        k = rnd.randint(1, 5)
        p = {(): 1}
        for x in rnd.sample(low, k):
            p = A.pmul(p, {(x,): 1})
        if p:
            samples.append(p)
    tests = [{(x,): 1} for x in G] + samples
    res = {}
    # off-diagonal xi = z^b d_a
    ok = okQ = okp = True
    for a_ in range(n):
        for b_ in range(n):
            if a_ == b_:
                continue

            def eta(lid, a_=a_, b_=b_):
                typ, i, al = A.info[lid]
                if typ == 'X':
                    return {(): 1} if (i == a_ and al == unit(n, b_)) else {}
                if typ == 'b':
                    if al[b_] >= 1:
                        be = list(al); be[b_] -= 1
                        return {(A.L('g', a_, tuple(be)),): -al[b_]}
                    return {}
                return {}
            L = gl_action(A, b_, a_)  # E_{ba}: weight shift e_b - e_a, same as xi = z^b d_a
            for p in tests:
                c = anti(A, eta, p)
                Lp = A.apply_der(p, L, False)
                if A.padd(c, Lp, 1):
                    ok = False
                if A.padd(c, Lp, -1):
                    okp = False
                if A.padd(A.Q(Lp), A.apply_der(A.Q(p), L, False), -1):
                    okQ = False
    res['[Q, iota_xi + h_xi] = -E_ba, xi = z^b d_a (a!=b)'] = ok
    res['[Q, iota_xi + h_xi] = +E_ba, xi = z^b d_a (a!=b)'] = okp
    res['[Q, E_ab] = 0 (a!=b)'] = okQ
    # torus xi = sum c_j z^j d_j with sum c = 0  (c = e_a - e_b)
    ok = okp = True
    for a_ in range(n):
        for b_ in range(n):
            if a_ == b_:
                continue
            c = [0] * n
            c[a_] += 1
            c[b_] -= 1

            def eta(lid, c=c):
                typ, i, al = A.info[lid]
                if typ == 'X':
                    return {(): c[i]} if (al == unit(n, i) and c[i]) else {}
                if typ == 'b':
                    r = {}
                    for k in range(n):
                        if al[k] >= 1 and c[k]:
                            be = list(al); be[k] -= 1
                            r = A.padd(r, {(A.L('g', k, tuple(be)),): -c[k] * al[k]})
                    return r
                return {}
            for p in tests:
                cc = anti(A, eta, p)
                exp = {}
                for m, v in p.items():
                    G0, t = A.weight(m)
                    ct = sum(c[j] * t[j] for j in range(n))
                    if ct:
                        exp[m] = -ct * v
                if A.padd(cc, exp, -1):
                    ok = False
                if A.padd(cc, exp, 1):
                    okp = False
    res['[Q, iota_c + h_c] = -(c.t) Id, c = e_a - e_b'] = ok
    res['[Q, iota_c + h_c] = +(c.t) Id, c = e_a - e_b'] = okp
    # the trace part c=(1,..,1) does not preserve the ideal: iota_c(div_0) = n
    return res, len(G), len(samples)


if __name__ == '__main__':
    out = open('t4_homotopy.log', 'w')
    log('T4 operator identities (t4_homotopy.py); Q with s_X = -1')
    for n in (2, 3):
        for quot in (False, True):
            res, ng, ns = check_translation(n, quot)
            log('n=%d %s: %d generators |a|<=3 + %d random monomials' % (n, 'C_quot' if quot else 'C_full', ng, ns))
            for k, v in res.items():
                log('   %-45s %s' % (k, v))
    for n in (2, 3):
        res, ng, ns = check_cartan(n)
        log('n=%d C_quot gl_n Cartan: %d generators + %d random monomials' % (n, ng, ns))
        for k, v in res.items():
            log('   %-55s %s' % (k, v))
