"""Fixed small N: the 9d complex for sl_N (kappa = tr) with explicit evaluation coordinates.
Formal trace-basis cochains (o2_9d_slN conventions) are restricted to sl_N (a chain map onto the invariant
cochains); coordinates = values at random explicit points x_v = X_v e^{p_v z}, X_v in sl_N, sum p = 0.
B = d(C^7_inv) from a maximal independent subset of the formal generators; d evaluated through
(d eta)(x) = sum_{i<j} (-1)^{i+j} eta([x_i,x_j], ...) on explicit commutators.
Exactness claims are verified by the residual at fresh samples."""
import sys
import itertools
import numpy as np
from fractions import Fraction as Fr
from itertools import combinations, combinations_with_replacement
from common import Log, PRIMES, np_rank_mod, np_solve_mod, det_int, ratrec, crt_pair
import colour as col
import trace_complex as tc
import cartan_lib as cl

log = Log(__file__)
p = int(__import__('os').environ.get('SMALLN_P', '1048573'))   # < 2^20: batched int64 matmul mod p is exact
H = 4


def rand_sl(N, rng):
    X = rng.integers(-3, 4, size=(N, N)).astype(np.int64)
    X[N - 1, N - 1] -= np.trace(X)
    return X


class Samples:
    def __init__(self, k, N, S, rng, base=None):
        self.k, self.N, self.S = k, N, S
        if base is None:
            self.X = np.array([[rand_sl(N, rng) for _ in range(k)] for _ in range(S)], dtype=np.int64) % p
            mom = rng.integers(-3, 4, size=(S, k - 1, H))
            self.P = np.concatenate([mom, -mom.sum(axis=1, keepdims=True)], axis=1).astype(np.int64)
        else:
            self.X, self.P = base
        self.minors = {(): np.ones(S, dtype=np.int64)}
        for I in combinations(range(k), H):
            self.minors[I] = tc._det_int(self.P[:, list(I), :]) % p
        self._trcache = {}

    def trace_word(self, w):
        if w not in self._trcache:
            M = self.X[:, w[0]]
            for a in w[1:]:
                M = np.matmul(M, self.X[:, a]) % p
            self._trcache[w] = np.trace(M, axis1=1, axis2=2) % p
        return self._trcache[w]

    def T(self, rho):
        v = np.ones(self.S, dtype=np.int64)
        for c in tc.cycles(rho):
            m = min(range(len(c)), key=lambda i: c[i:] + c[:i])
            w = tuple(c[m:] + c[:m])
            v = v * self.trace_word(w) % p
        return v

    def derived(self):
        """for d: the S*C(k,2) explicit (k-1)-point samples ([x_i,x_j], rest), with signs."""
        Xs, Ps, sg = [], [], []
        for i in range(self.k):
            for j in range(i + 1, self.k):
                rest = [a for a in range(self.k) if a not in (i, j)]
                A, B = self.X[:, i], self.X[:, j]
                C = (np.matmul(A, B) - np.matmul(B, A)) % p
                Xs.append(np.concatenate([C[:, None], self.X[:, rest]], axis=1))
                Ps.append(np.concatenate([(self.P[:, i] + self.P[:, j])[:, None], self.P[:, rest]], axis=1))
                sg.append(np.full(self.S, (-1) ** (i + j)))
        # order: sample-major
        X = np.stack(Xs, axis=1).reshape(self.S * len(Xs), self.k - 1, self.N, self.N)
        P = np.stack(Ps, axis=1).reshape(self.S * len(Ps), self.k - 1, H)
        sg = np.stack(sg, axis=1).reshape(-1)
        return Samples(self.k - 1, self.N, X.shape[0], None, base=(X, P)), sg


FPF = {}


def fpf_perms(k):
    if k not in FPF:
        FPF[k] = [r for r in itertools.permutations(range(k)) if all(r[i] != i for i in range(k))]
    return FPF[k]


def evaluate(om, smp):
    """omega(x) at all samples (mod p)."""
    out = np.zeros(smp.S, dtype=np.int64)
    for rho in fpf_perms(om.k):
        lam, sig0 = tc.conjugator(rho)
        P = om.data.get(lam)
        if not P:
            continue
        Q = tc.poly_pullback(P, sig0, H)
        if not Q:
            continue
        c = np.zeros(smp.S, dtype=np.int64)
        for (I, J), v in Q.items():
            c = (c + cl.fmod(v, p) * (smp.minors[I] * smp.minors[J] % p)) % p
        out = (out + tc.sgn(sig0) * c % p * smp.T(rho)) % p
    return out % p


def evaluate_d(om, smp, der=None):
    if der is None:
        der = smp.derived()
    dsmp, sg = der
    vals = evaluate(om, dsmp)
    npair = smp.k * (smp.k - 1) // 2
    return (vals * sg % p).reshape(smp.S, npair).sum(axis=1) % p


def theta_cochain(N):
    o = None
    for shape in cl.SHAPES:
        om = tc.Cochain(8, H)
        for perm, w in col.colour_perm_dict(shape, N).items():
            om.add_term(perm, {(I, J): e for I, J, e in cl.DTH[shape]}, cl.SIGMA * w)
        om = om.scaled(cl.KC[shape] / cl.AUT[shape])
        o = om if o is None else o.plus(om)
    return o


def run(N, S7=40, S8=80, seed=0):
    rng = np.random.default_rng(seed)
    s7 = Samples(7, N, S7, rng)
    s8 = Samples(8, N, S8, rng)
    s8b = Samples(8, N, 30, rng)                  # fresh samples for verification
    subs = list(combinations(range(6), 4))
    gens = []
    for lam in tc.fpf_types(7):
        for I, J in combinations_with_replacement(subs, 2):
            g = tc.Cochain(7, H)
            g.add_term(tc.canonical(lam), {(I, J): 1}, 1)
            if g.data:
                gens.append(g)
    # maximal independent subset on C^7 coordinates
    rows, basis = [], []
    r = 0
    for g in gens:
        v = evaluate(g, s7)
        if not np.any(v):
            continue
        trial = rows + [v]
        rr = np_rank_mod(np.array(trial), p)
        if rr > r:
            rows, r = trial, rr
            basis.append(g)
    log('sl_%d: dim C^7_inv >= %d (rank of %d formal generators at %d samples)' % (N, r, len(gens), S7))
    der8, der8b = s8.derived(), s8b.derived()
    Bm = np.array([evaluate_d(b, s8, der8) for b in basis]) if basis else np.zeros((0, S8), dtype=np.int64)
    Bmb = np.array([evaluate_d(b, s8b, der8b) for b in basis]) if basis else np.zeros((0, 30), dtype=np.int64)
    o = theta_cochain(N)
    vo, vob = evaluate(o, s8), evaluate(o, s8b)
    T23 = tc.sym_trace_dict((2, 3))
    T5 = tc.sym_trace_dict((5,))
    Qs = {'Q2323': tc.half_bracket_sym(T23, H, 120, N)}
    if N >= 5:
        Qs['Q55'] = tc.half_bracket_sym(T5, H, 120, N)
        mix = {k: T5.get(k, 0) + T23.get(k, 0) for k in set(T5) | set(T23)}
        Qs['Qmix'] = tc.half_bracket_sym(mix, H, 120, N)
    vQ = {k: evaluate(v, s8) for k, v in Qs.items()}
    vQb = {k: evaluate(v, s8b) for k, v in Qs.items()}
    rB = np_rank_mod(Bm, p) if len(Bm) else 0
    rBo = np_rank_mod(np.vstack([Bm, vo[None]]), p)
    rQ = np_rank_mod(np.vstack([Bm] + [vQ[k][None] for k in vQ]), p)
    rQo = np_rank_mod(np.vstack([Bm] + [vQ[k][None] for k in vQ] + [vo[None]]), p)
    log('sl_%d: rank B = %d, rank[B;o\'] = %d, rank[B;Q] = %d, rank[B;Q;o\'] = %d  (S8 = %d), o\' nonzero: %s'
        % (N, rB, rBo, rQ, rQo, S8, bool(np.any(vo))))
    # solve o' = sum x_Q Q + B y, verify on fresh samples
    A = np.vstack([vQ[k][None] for k in vQ] + [Bm]).T
    ok, x, _ = np_solve_mod(A, vo, p)
    if ok:
        Ab = np.vstack([vQb[k][None] for k in vQb] + [Bmb]).T
        resid = (Ab @ x - vob) % p
        xs = [ratrec(int(v), p) for v in x[:len(vQ)]]
        log('sl_%d: o\' = sum x_Q Q + B y solvable; x_Q =' % N, dict(zip(vQ, xs)), ' residual at 30 fresh samples zero:',
            not np.any(resid))
    else:
        log('sl_%d: o\' not in span(Q, B) at these coordinates' % N)
    okB, xB, _ = np_solve_mod(Bm.T, vo, p) if len(Bm) else (not np.any(vo), None, 0)
    if okB and xB is not None:
        resid = (Bmb.T @ xB - vob) % p
        log('sl_%d: o\' in B; residual at fresh samples zero:' % N, not np.any(resid))
    return rB, rBo, rQ, rQo


Ns = [int(a) for a in sys.argv[1:]] or [2, 3]
S7 = int(__import__('os').environ.get('S7', '40'))
S8 = int(__import__('os').environ.get('S8', '80'))
log('prime p =', p, ' S7 =', S7, ' S8 =', S8)
for N in Ns:
    run(N, S7=S7, S8=S8)
log.close()
