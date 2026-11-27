"""Order-hbar counterterm J_1 = a j_T, j_T = alt([p_1..p_4] T), T in (Sym^5 g^*)^g, on Cartan cycles:
does o'|_C lie in span{1/2{j,j}|_C} (necessary for [o'] in the image of the quadratic map)?
sl_N: T_5 = Sym tr X^5, T_23 = Sym tr X^2 tr X^3 (defining rep, kappa = tr); E_6: T_5 = Sym tr_27 X^5."""
import sys
import numpy as np
from fractions import Fraction as Fr
from common import Log, PRIMES, np_rank_mod
import lie
import cartan_lib as cl

log = Log(__file__)
p = PRIMES[0]
NS = 6


def orbit(lam, A):
    seen = {tuple(lam)}
    todo = [tuple(lam)]
    while todo:
        x = todo.pop()
        for i in range(len(x)):
            if x[i] != 0:
                y = tuple(x[j] - x[i] * A[j][i] for j in range(len(x)))
                if y not in seen:
                    seen.add(y); todo.append(y)
    return sorted(seen)


def simple_root_matrix(L):
    """A[i][j] = alpha_j(h_i): weights (under the Cartan basis) of the simple root vectors."""
    r = len(L.cartan)
    wts = [tuple(int(x) for x in L.wt[a]) for a in range(L.dim)]
    nz = [w for w in wts if any(w)]
    gen = [10 ** (3 * i) + i for i in range(r)]
    pos = [w for w in nz if sum(a * g for a, g in zip(w, gen)) > 0]
    ps = set(pos)
    simple = [w for w in pos if not any(tuple(a - b for a, b in zip(w, u)) in ps for u in pos if u != w)]
    order = []
    for i in range(r):
        order.append(next(w for w in simple if w[i] == 2 and w not in order))
    return [[order[j][i] for j in range(r)] for i in range(r)]


def run(key, Tpolys):
    L = lie.build(key)
    r = len(L.cartan)
    detc = cl.Detector(L, p)
    Ts = {name: cl.sym_tensor_from_poly(P, r, p) for name, P in Tpolys.items()}
    names = list(Ts)
    rng = np.random.default_rng(31)
    rows = []
    for t in range(NS):
        hc, mom = detc.sample(rng)
        o = detc.o_prime(hc, mom)
        q = {}
        for a in names:
            q[(a, a)] = cl.half_bracket_tensor(L, Ts[a], hc, mom, p)
        if len(names) == 2:
            a, b = names
            Tab = {k: (Ts[a][k] + Ts[b][k]) % p for k in Ts[a]}
            q[(a, b)] = (cl.half_bracket_tensor(L, Tab, hc, mom, p) - q[(a, a)] - q[(b, b)]) % p
        rows.append([o] + [q[k] for k in sorted(q)])
    M = np.array(rows, dtype=np.int64)
    rq = np.linalg.matrix_rank(M[:, 1:].astype(float)) if False else np_rank_mod(M[:, 1:].T, p)
    ro = np_rank_mod(M.T, p)
    log('%s (%s): samples %d; rank{Q|_C} = %d, rank{Q|_C, o\'|_C} = %d  -> o\'|_C in span: %s; Q keys %s'
        % (key, L.name, NS, rq, ro, rq == ro, sorted(q)))
    return rq, ro


def sl_weights(N):
    # mu_k(h_i) = delta_{k,i} - delta_{k,i+1} in the Kac basis (h_i = alpha_i), verified via the trace form below
    return [[(1 if k == i else 0) - (1 if k == i + 1 else 0) for i in range(N - 1)] for k in range(N)]


for N in ([int(x) for x in sys.argv[1:]] or range(3, 10)):
    W = sl_weights(N)
    L = lie.build('A%d' % (N - 1))
    # check: kappa(h_i,h_j) = sum_k mu_k(h_i) mu_k(h_j) (trace form) for the Kac basis
    K = L.kap.toarray()
    ok = all(int(K[L.cartan[i], L.cartan[j]]) == sum(w[i] * w[j] for w in W) for i in range(N - 1) for j in range(N - 1))
    def P5(x, W=W): return sum(sum(w[i] * x[i] for i in range(len(x))) ** 5 for w in W)
    def P23(x, W=W): return sum(sum(w[i] * x[i] for i in range(len(x))) ** 2 for w in W) * sum(sum(w[i] * x[i] for i in range(len(x))) ** 3 for w in W)
    log('sl_%d: Kac Cartan form = trace form of C^N: %s' % (N, ok))
    polys = {'T5': P5, 'T23': P23} if N >= 5 else {'T23': P23}
    run('A%d' % (N - 1), polys)

L6 = lie.build('E6')
A = [list(map(int, row)) for row in np.array(lie.cartan_E(6))]   # (alpha_i|alpha_j), Kac node order = Cartan basis order
for i in range(6):
    lam = [0] * 6; lam[i] = 1
    orb = orbit(lam, A)
    if len(orb) == 27:
        W27 = orb
        log('E6: 27 = W-orbit of fundamental weight', i, ' size', len(orb))
        break
def P5e(x): return sum(sum(w[i] * x[i] for i in range(6)) ** 5 for w in W27)
def P2e(x): return sum(sum(w[i] * x[i] for i in range(6)) ** 2 for w in W27)
log('E6: tr_27 X^5 on h nonzero:', any(P5e(list(np.random.default_rng(s).integers(-3, 4, 6))) for s in range(3)))
run('E6', {'T5': P5e})
# so(10): (Sym^5 g*)^g = C Pf (degree-5 Pfaffian); Pf|_h = prod of one weight from each pair +-mu of the vector rep
L10 = lie.build('D5')
A10 = [list(map(int, row)) for row in np.array(lie.cartan_D(5))]
for i in range(5):
    lam = [0] * 5; lam[i] = 1
    orb = orbit(lam, A10)
    if len(orb) == 10:
        V10 = orb
        log('so10: vector rep = W-orbit of fundamental weight', i, ' size', len(orb))
        break
half = []
for w in V10:
    if tuple(-x for x in w) not in [tuple(h) for h in half]:
        half.append(w)
assert len(half) == 5
def Pf(x): return int(np.prod([sum(w[i] * x[i] for i in range(5)) for w in half]))
log('so10: Pf|_h nonzero:', any(Pf(list(np.random.default_rng(s).integers(-3, 4, 5))) for s in range(3)))
run('D5', {'Pf': Pf})
log.close()
