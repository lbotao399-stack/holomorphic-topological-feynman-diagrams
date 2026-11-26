"""Q_1^2 + {Q_0, Q_2} = 0 for holomorphic CS on C^{4|1} (adapted from S2 code/test_q2.py).

Solve  hat Q_1^2 + sum_Gamma k_Gamma {hat Q_0, hat Q_Gamma} = 0  for (k_144, k_234, k_333) on test operators;
prediction from the assembly (S2 Section 5.5 with 3n+1 vertices, raw factor (16 pi^4)^{-n}):
|k_Gamma| = 4/|Aut Gamma| = 2, 4, 2/3, with J_Gamma = the density g of K_Gamma (K_Gamma = int g).
Usage: python3 test_q2_c41.py ALGEBRA TEST [TEST ...]
"""
import sys, time, pickle
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/code')
import sympy as sp
from htk import sp_reduce
from graphs_c4 import theta_O, THETAS, AUT
from ops_c41 import structure_constants, add, mul, Theory, TwoLoop, FastTwoLoop, L

H = 4
import glob
res = {}
for PKL in sorted(glob.glob((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/code/kernels_c4_*.pkl')):
    res.update(pickle.load(open(PKL, 'rb')))


class KernelData:
    """Weights {m: Fraction} of the theta kernels (Cauchy-Binet, weights_c4.py, stored by kernels_c4.py):
    ks = 0 and first order in one z_j."""
    def __init__(self, pqr):
        r = res[pqr]
        self.pqr = pqr
        self.outputs = r['outputs']
        self.paths = r['paths']
        self.w = r['weights']

    def __call__(self, ks):
        return self.w[ks]


def main():
    alg = sys.argv[1] if len(sys.argv) > 1 else 'su2'
    which = sys.argv[2:] or ['g0 g1 g2']
    n, f = structure_constants(alg)
    T = Theory(n, f)
    kd = {}
    for pqr in THETAS:
        K = KernelData(pqr)
        arms = [list(K.paths[j][1]) for j in range(3)]
        kd[pqr] = dict(outputs=K.outputs, arms=arms, K=K)
    T2 = FastTwoLoop(T, kd)
    e = lambda a: tuple(1 if i == a else 0 for i in range(H))
    tests = {
        'g0 g1 g2': mul(mul(L('g', 0), L('g', 1)), L('g', 2)),
        'g0 g1 c2': mul(mul(L('g', 0), L('g', 1)), L('c', 2)),
        'g0 g0 g1': mul(mul(L('g', 0), L('g', 0)), L('g', 1)),
        'g0 g1 g3': mul(mul(L('g', 0), L('g', 1)), L('g', 3)),
        'g0 g1 c3': mul(mul(L('g', 0), L('g', 1)), L('c', 3)),
        'd1g0 g1 g2': mul(mul(L('g', 0, e(0)), L('g', 1)), L('g', 2)),
        'g0 g1 d2g2': mul(mul(L('g', 0), L('g', 1)), L('g', 2, e(1))),
        'd1g0 g1 c2': mul(mul(L('g', 0, e(0)), L('g', 1)), L('c', 2)),
        'g0 g1 d3c2': mul(mul(L('g', 0), L('g', 1)), L('c', 2, e(2))),
        'g0 g1 g2 c0': mul(mul(mul(L('g', 0), L('g', 1)), L('g', 2)), L('c', 0)),
    }
    print('algebra', alg, 'dim', n, flush=True)
    rows = []
    for nm in which:
        p = tests[nm]
        t0 = time.time()
        lhs = T.Q1(T.Q1(p))
        t1 = time.time()
        print('%s: |hatQ1^2 p| = %d (%.1fs)' % (nm, len(lhs), t1 - t0), flush=True)
        anti = {}
        for pqr in THETAS:
            ta = time.time()
            q2p = T2.Q2(pqr, p)
            a = add(T.Q0(q2p), T2.Q2(pqr, T.Q0(p)))
            anti[pqr] = a
            print('   theta%s: |hatQ_Gamma p| = %d, |{hatQ0, hatQ_Gamma} p| = %d (%.1fs)' % (pqr, len(q2p), len(a), time.time() - ta), flush=True)
        rows.append((nm, lhs, anti))
    eqs = []
    for nm, lhs, anti in rows:
        mons = set(lhs)
        for a_ in anti.values():
            mons |= set(a_)
        for m in mons:
            eqs.append([Fr(anti[pqr].get(m, 0)) for pqr in THETAS] + [-Fr(lhs.get(m, 0))])
    rank_x, sol, consistent = solve_small(eqs, 3)
    # proportionality constants (rank-deficient case): X_Gamma = a_Gamma X_ref, hatQ1^2 = b X_ref on each operator
    for nm, lhs, anti in rows:
        ref = anti[THETAS[0]]
        if not ref:
            continue
        m0 = next(iter(ref))
        a_ = {pqr: Fr(anti[pqr].get(m0, 0)) / Fr(ref[m0]) for pqr in THETAS}
        b_ = Fr(lhs.get(m0, 0)) / Fr(ref[m0])
        prop = all(all(Fr(anti[pqr].get(m, 0)) == a_[pqr] * Fr(ref.get(m, 0)) for m in set(ref) | set(anti[pqr])) for pqr in THETAS)
        propl = all(Fr(lhs.get(m, 0)) == b_ * Fr(ref.get(m, 0)) for m in set(ref) | set(lhs))
        print('   %s: {hatQ0,hatQ_Gamma} = a_Gamma X_ref with a = %s (exact on all monomials: %s); hatQ1^2 = b X_ref, b = %s (%s)'
              % (nm, [str(a_[q]) for q in THETAS], prop, b_, propl))
        import itertools as it
        for sg in it.product((1, -1), repeat=3):
            if sum(sg[i] * Fr(4, AUT[q]) * a_[q] for i, q in enumerate(THETAS)) + b_ == 0:
                print('      sign pattern (eps_144, eps_234, eps_333) = %s with |k| = 4/|Aut| solves this operator' % (sg,))
    import pickle as pk_
    pk_.dump({nm: (lhs, anti) for nm, lhs, anti in rows}, open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/code/q2_rows_%s_%s.pkl' % (alg, len(rows)), 'wb'))
    print('equations:', len(eqs), ' rank [X_144 X_234 X_333]:', rank_x, ' consistent:', consistent, flush=True)
    print('solution (k144, k234, k333):', sol)
    print('prediction |k| = 4/|Aut| :', [Fr(4, AUT[pqr]) for pqr in THETAS])
    for i, pqr in enumerate(THETAS):
        print('   column theta%s nonzero entries: %d' % (pqr, sum(1 for r in eqs if r[i] != 0)))
    nz_lhs = sum(1 for r in eqs if r[3] != 0)
    print('   nonzero entries of hatQ1^2: %d' % nz_lhs)
    if sol is not None and all(x is not None for x in sol):
        res_ = [sum(r[i] * sol[i] for i in range(3)) - r[3] for r in eqs]
        print('   residual of all equations at the solution is zero:', all(x == 0 for x in res_))


def solve_small(eqs, n):
    """Exact Gaussian elimination for rows [a_1..a_n | b]; returns (rank of A, solution or None, consistent)."""
    piv = []          # list of (col, row) reduced rows
    basis = []
    for r in eqs:
        r = list(r)
        for (c, br) in basis:
            if r[c] != 0:
                f = r[c] / br[c]
                r = [x - f * y for x, y in zip(r, br)]
        c = next((i for i in range(n) if r[i] != 0), None)
        if c is None:
            if r[n] != 0:
                return len(basis), None, False
            continue
        basis.append((c, r))
        if len(basis) == n:
            pass
    rank = len(basis)
    # back substitution
    sol = [None] * n
    if rank == n:
        basis.sort(key=lambda cb: cb[0])
        for c, br in reversed(basis):
            val = br[n] - sum(br[j] * sol[j] for j in range(c + 1, n) if sol[j] is not None)
            sol[c] = val / br[c]
        # consistency of all rows
        cons = all(sum(r[i] * sol[i] for i in range(n)) == r[n] for r in eqs)
        return rank, sol, cons
    return rank, None, True


if __name__ == '__main__':
    main()
