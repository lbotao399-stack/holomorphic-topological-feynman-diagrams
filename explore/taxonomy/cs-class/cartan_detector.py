"""9d two-loop obstruction on Cartan cycles for all simple g (lie.build constructions, exact integer structure
constants), mod two primes.  o'|_C != 0  =>  [O_2] != 0 in H^8 of sl... g[[z_1..z_4]] (coinvariants).
Rank 2: o'|_C = c_g det(h_v (x) p_v) (Lambda^8 of h (x) C^4 is one-dimensional): ratio printed.
sl_N, E_6: J_1 brackets on the same cycles (T = Sym tr_V X^5, V = C^N, 27)."""
import sys
import numpy as np
from fractions import Fraction as Fr
from common import Log, PRIMES, ratrec, crt_pair, det_int
import lie
import cartan_lib as cl

log = Log(__file__)
keys = sys.argv[1:] or ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'B2', 'B3', 'B4', 'C3', 'C4', 'D4', 'D5', 'G2', 'F4', 'E6', 'E7', 'E8']
rng = np.random.default_rng(11)


def kron_det(hc, mom, r):
    """det of the 8 vectors h_v (x) p_v in C^r (x) C^4 (only r = 2)."""
    rows = []
    for v in range(8):
        rows.append([int(hc[v, i]) * int(mom[v, k]) for i in range(r) for k in range(4)])
    return det_int(rows)


def minuscule_weights(L, n_target):
    """weights (Dynkin labels w.r.t. the Cartan basis of L) of the W-orbit of a fundamental weight of size n_target."""
    r = len(L.cartan)
    A = np.array([[int(L.wt[L.cartan[j], i]) for i in range(r)] for j in range(r)])  # not used
    # Cartan matrix from simple root weights: simple root vectors are the basis elements whose weight is a row
    # of the Cartan matrix; obtain A_ij = alpha_j(h_i) from L.wt of the simple root vectors
    return None


def orbit(lam, Acart):
    seen = {tuple(lam)}
    todo = [tuple(lam)]
    while todo:
        x = todo.pop()
        for i in range(len(x)):
            if x[i] != 0:
                y = tuple(x[j] - x[i] * Acart[i][j] for j in range(len(x)))
                if y not in seen:
                    seen.add(y); todo.append(y)
    return sorted(seen)


def cartan_matrix(L):
    """A[i][j] = alpha_j(h_i) for simple roots alpha_j (positive roots with minimal height: weights that are
    not sums of two positive weights)."""
    r = len(L.cartan)
    wts = [tuple(int(x) for x in L.wt[a]) for a in range(L.dim)]
    nz = [w for w in wts if any(w)]
    # choose positivity by a generic functional
    gen = [1000 ** i + 7 * i for i in range(r)]
    pos = [w for w in nz if sum(a * g for a, g in zip(w, gen)) > 0]
    pos_set = set(pos)
    simple = [w for w in pos if not any(tuple(a - b for a, b in zip(w, u)) in pos_set for u in pos if u != w)]
    assert len(simple) == r, (len(simple), r)
    # A[i][j] = alpha_j(h_i) where alpha_j weight vector gives alpha_j(h_i) = w[i]
    # reorder simple roots so that alpha_j(h_j) = 2 matches h_j
    order = []
    for i in range(r):
        cand = [w for w in simple if w[i] == 2 and w not in order]
        order.append(cand[0])
    A = [[order[j][i] for j in range(r)] for i in range(r)]
    return A, order


for key in keys:
    L = lie.build(key)
    r = len(L.cartan)
    out = {}
    for p in PRIMES[:2]:
        det_ = cl.Detector(L, p)
        rng_p = np.random.default_rng(100)
        vals = []
        for t in range(3):
            hc, mom = det_.sample(rng_p)
            v = det_.o_prime(hc, mom)
            ratio = None
            if r == 2:
                kd = kron_det(hc, mom, r) % p
                ratio = v * pow(kd, p - 2, p) % p if kd else None
            vals.append((v, ratio, hc, mom))
        out[p] = vals
    v0 = [x[0] for x in out[PRIMES[0]]]
    nonzero = any(v0) or any(x[0] for x in out[PRIMES[1]])
    msg = '%s %s dim %d rank %d: o\'|_C mod p (3 samples): %s  nonzero: %s' % (key, L.name, L.dim, r, v0, nonzero)
    if r == 2:
        rr = []
        for (a0, b0, _, _), (a1, b1, _, _) in zip(out[PRIMES[0]], out[PRIMES[1]]):
            if b0 is not None and b1 is not None:
                vv, m = crt_pair(b0, PRIMES[0], b1, PRIMES[1])
                rr.append(ratrec(vv, m))
        msg += '  o\'|_C / det(h (x) p) = %s (kappa of the construction)' % rr
    log(msg)
log.close()
