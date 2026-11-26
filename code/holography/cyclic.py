"""Cyclic homology of A = C[z1,z2] (x) Lambda[theta1,theta2,theta3] (z even, theta odd), two ways.

(1) Connes complex.  Monomial key (m1, m2, mask), mask bit I <-> theta_{I+1}; parity |mask| mod 2;
    weight (m1, m2, bit0, bit1, bit2).  Chains x = (a0,...,an) in A^{(x)(n+1)} (reduced: all a_i != 1).
      b(x)  = sum_{i<n} (-1)^i (..., a_i a_{i+1}, ...) + s(x) (a_n a_0, a_1, ..., a_{n-1}),
      t(x)  = s(x) (a_n, a_0, ..., a_{n-1}),     s(x) = (-1)^{n + |a_n|(|a_0|+...+|a_{n-1}|)},
      C^lambda_n = A^{(x)(n+1)}/(1-t),  HC_n = H_n(C^lambda, b).
(2) Cyclic-word model of planar single-trace cochains of gl_N (x) A with the CE differential.
    Letter C_{m,J} (component of the superfield C = sum_{m,J} z^m/m! theta^J C_{m,J}), parity 1+|J|,
    weight (m; J).  Q C = C.C gives
      Q C_{M,K} = (-1)^{|K|} sum_{m+m'=M, J u J'=K} binom(M,m) eps(J,J') (-1)^{|J'|(1+|J|)} C_{m,J} C_{m',J'},
    theta^J theta^J' = eps(J,J') theta^K.  Q acts on tr(X_1...X_Y) as an odd derivation from the left;
    tr(X_1...X_Y) = (-1)^{|X_Y|(|X_1|+...+|X_{Y-1}|)} tr(X_Y X_1 ... X_{Y-1}).
    Relative complex (constant gl_N): words without c = C_{0,{}}; the c-terms of Q form ad_c and cancel
    on traces.  LLC sector (theta1,theta2,theta3; n+, n-)_Y = weight (n+, n-; theta) with Y letters.
"""
from fractions import Fraction
from itertools import combinations, product
from math import comb

from .linalg import rank, vadd

ONE = (0, 0, 0)


def popc(x):
    return bin(x).count("1")


def theta_sign(M1, M2):
    """theta^{M1} theta^{M2} = sign theta^{M1|M2} (increasing order); 0 if overlap."""
    if M1 & M2:
        return 0
    inv = 0
    for i in range(3):
        if M1 >> i & 1:
            inv += sum(1 for j in range(i) if M2 >> j & 1)
    return (-1) ** inv


def amul(a, b):
    s = theta_sign(a[2], b[2])
    if s == 0:
        return 0, None
    return s, (a[0] + b[0], a[1] + b[1], a[2] | b[2])


def apar(a):
    return popc(a[2]) % 2


def aweight(a):
    return (a[0], a[1], a[2] & 1, a[2] >> 1 & 1, a[2] >> 2 & 1)


# ----------------------------------------------------------------- chains
def compositions(total, parts):
    if parts == 0:
        if total == 0:
            yield ()
        return
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in compositions(total - first, parts - 1):
            yield (first,) + rest


def chains(beta, n, reduced=False):
    """All (n+1)-tuples of monomials with total weight beta."""
    L = n + 1
    out = []
    th_choices = [list(combinations(range(L), beta[2 + I])) for I in range(3)]
    for c1 in compositions(beta[0], L):
        for c2 in compositions(beta[1], L):
            for t0, t1, t2 in product(*th_choices):
                masks = [0] * L
                for I, slots in enumerate((t0, t1, t2)):
                    for s in slots:
                        masks[s] |= 1 << I
                x = tuple((c1[i], c2[i], masks[i]) for i in range(L))
                if reduced and any(a == ONE for a in x):
                    continue
                out.append(x)
    return out


def s_sign(x):
    n = len(x) - 1
    pre = sum(apar(a) for a in x[:-1])
    return (-1) ** (n + apar(x[-1]) * pre)


def canon(x, sign_fn=s_sign):
    """[x] = c [rep] in the quotient by (1 - t);  returns (rep, c) or None if [x] = 0."""
    cur, c = x, 1
    best = {}
    for _ in range(len(x)):
        if cur in best:
            if best[cur] != c:
                return None
        else:
            best[cur] = c
        c = c * sign_fn(cur)
        cur = (cur[-1],) + cur[:-1]
    # after a full turn we must return to x with the accumulated sign
    if c != 1:
        return None
    rep = min(best)
    # [rot^j x] = c_j [x]  =>  [x] = c_j [rep]
    return rep, best[rep]


def hochschild_b(x):
    n = len(x) - 1
    if n == 0:
        return {}
    out = {}
    for i in range(n):
        s, ab = amul(x[i], x[i + 1])
        if s:
            y = x[:i] + (ab,) + x[i + 2:]
            out[y] = out.get(y, 0) + (-1) ** i * s
    s, ab = amul(x[n], x[0])
    if s:
        y = (ab,) + x[1:n]
        out[y] = out.get(y, 0) + s_sign(x) * s
    return {k: v for k, v in out.items() if v}


def t_op(x):
    return s_sign(x), (x[-1],) + x[:-1]


def lam_basis(beta, n, reduced=False):
    reps = set()
    for x in chains(beta, n, reduced):
        r = canon(x)
        if r is not None:
            reps.add(r[0])
    return sorted(reps)


def b_on_lambda(rep, reduced=False):
    out = {}
    for y, c in hochschild_b(rep).items():
        if reduced and any(a == ONE for a in y):
            # in the complex of the augmentation ideal products never produce 1
            raise AssertionError("unit produced in reduced complex")
        r = canon(y)
        if r is not None:
            vadd(out, {r[0]: Fraction(r[1])}, c)
    return out


def HC(beta, n, reduced=False):
    """dim HC_n(A) (or of the augmentation ideal) at weight beta, via the Connes complex."""
    Bn = lam_basis(beta, n, reduced)
    Bn1 = lam_basis(beta, n + 1, reduced)
    rk_out = rank([b_on_lambda(r, reduced) for r in Bn]) if n >= 1 else 0
    rk_in = rank([b_on_lambda(r, reduced) for r in Bn1])
    return len(Bn) - rk_out - rk_in, len(Bn)


# ----------------------------------------------------------------- cyclic-word CE model
def letters_of_weight(w):
    """All letters (m1, m2, mask) of weight w (w = 5-tuple, theta entries 0/1)."""
    if min(w) < 0 or max(w[2:]) > 1:
        return []
    return [(w[0], w[1], w[2] | w[3] << 1 | w[4] << 2)]


def lpar(L):
    return (1 + popc(L[2])) % 2


def word_sign(x):
    pre = sum(lpar(a) for a in x[:-1])
    return (-1) ** (lpar(x[-1]) * pre)


def Q_letter(L):
    """Q C_{M,K} as dict {(L1, L2): coeff} (ordered matrix product L1 L2)."""
    M1, M2, K = L
    out = {}
    sub = [J for J in range(8) if J & ~K == 0]
    for m1 in range(M1 + 1):
        for m2 in range(M2 + 1):
            bc = comb(M1, m1) * comb(M2, m2)
            for J in sub:
                Jp = K & ~J
                eps = theta_sign(J, Jp)
                s = (-1) ** popc(K) * bc * eps * (-1) ** (popc(Jp) * (1 + popc(J)))
                key = ((m1, m2, J), (M1 - m1, M2 - m2, Jp))
                out[key] = out.get(key, 0) + s
    return {k: v for k, v in out.items() if v}


def Q_word(x, relative=False):
    out = {}
    pre = 0
    for i, L in enumerate(x):
        for (L1, L2), c in Q_letter(L).items():
            if relative and (L1 == ONE or L2 == ONE):
                continue
            y = x[:i] + (L1, L2) + x[i + 1:]
            r = canon(y, word_sign)
            if r is not None:
                vadd(out, {r[0]: Fraction(r[1])}, c * (-1) ** pre)
        pre += lpar(L)
    return out


def words(beta, Y, relative=False):
    reps = set()
    for x in chains(beta, Y - 1, reduced=relative):
        r = canon(x, word_sign)
        if r is not None:
            reps.add(r[0])
    return sorted(reps)


def CE_single_trace(beta, Y, relative=False):
    """dim of planar single-trace CE cohomology with Y letters at weight beta."""
    WY = words(beta, Y, relative)
    Wm = words(beta, Y - 1, relative) if Y >= 2 else []
    rk_out = rank([Q_word(x, relative) for x in WY])
    rk_in = rank([Q_word(x, relative) for x in Wm])
    return len(WY) - rk_out - rk_in, len(WY)
