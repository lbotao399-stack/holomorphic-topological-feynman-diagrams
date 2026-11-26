"""(f) sl_N at fixed N = 3, 4, 5 (kappa = tr): faithful evaluation of trace-basis cochains
(formal coordinates for N >= k; explicit.Evaluator tensor mode for N < k except (N,k) = (4,5): sample mode with
the antisymmetrizer identity certifying the rank, as in deform/code/algebras.py).
  H^3 at (2,2) and (3,3), H^4 at (3,3);  exactness of theta', {D, iota}, o_3' = {D,iota} + theta' at N;
  values of the formal primitives: o_3' == d xi at N.
Usage: python3 o9_fixedN.py N"""
import sys, time
from common import *
import explicit

N = int(sys.argv[1])
sl = True
t0 = time.time()


class Coord:
    def __init__(self, N, sl, k):
        self.N, self.sl, self.k = N, sl, k
        self.ev = None
        if N < k:
            if k == 5 and N == 4:
                ev = explicit.Evaluator(N, sl, k, mode='samples', nsamp=(160 if not sl else 70), seed=11)
                rel = explicit.antisym_relation(k, sl)
                okrel = explicit.check_relation_on_samples(rel, N, sl, k)
                ev.exact = okrel and ev.rank == len(ev.taus) - 1
                print("    antisymmetrizer identity (k=5) at N=4 on 6 fresh samples: %s; rank %d = %d - 1: exact %s"
                      % (okrel, ev.rank, len(ev.taus), ev.exact), flush=True)
                assert ev.exact
            else:
                ev = explicit.Evaluator(N, sl, k, mode='tensor')
            self.ev = ev

    def vec(self, c):
        if self.ev is not None:
            return self.ev.vec(c)
        out = {}
        for (t, m), v in c.items():
            if self.sl and ce.has_singleton(t):
                continue
            x = ce.np_eval(v, self.N) if isinstance(v, dict) else Fr(v)
            if x != 0:
                out[(t, m)] = x
        return out


co = {k: Coord(N, sl, k) for k in (3, 4, 5)}


def npl(c):
    return {k: (d if isinstance(d, dict) else {0: Fr(d)}) for k, d in c.items()}


def basis_from(cochains, coord):
    vecs = [coord.vec(c) for c in cochains]
    vecs_nz = [(c, v) for c, v in zip(cochains, vecs) if v]
    if not vecs_nz:
        return [], 0
    keys, M = la.matrix_from([v for _, v in vecs_nz])
    cert = la.certified_rank(M)
    return [vecs_nz[j][0] for j in cert['pivcols']], cert['rank']


def rank(cochains, coord):
    vecs = [coord.vec(c) for c in cochains]
    vecs = [v for v in vecs if v]
    if not vecs:
        return 0
    keys, M = la.matrix_from(vecs)
    return la.certified_rank(M)['rank']


def cohom(deg, want_Z4=False):
    C2, _ = basis_from([npl(c) for c in spaces.spanning(2, deg, sl)], co[3]) if False else ([], 0)
    C3, _ = basis_from([npl(c) for c in spaces.spanning(3, deg, sl)], co[3])
    C4, _ = basis_from([npl(c) for c in spaces.spanning(4, deg, sl)], co[4])
    B4 = [ce.d_ce(c, 3, sl) for c in C3]
    r3 = rank(B4, co[4])
    r4 = rank([ce.d_ce(c, 4, sl) for c in C4], co[5])
    return dict(C3=C3, C4=C4, B4=B4, r3=r3, r4=r4)


print("sl%d" % N)
c2 = spaces.spanning(2, (2, 2), sl)
print("  C^2 spanning set at (2,2), (3,3): %d, %d elements (alt of symmetric 2-tensors x even monomials)"
      % (len(c2), len(spaces.spanning(2, (3, 3), sl))))
for deg in [(2, 2), (3, 3)]:
    R = cohom(deg)
    print("  weight %s: dim C^3 = %d, dim C^4 = %d, rank d3 = %d, rank d4 = %d: H^3 = %d, H^4 = %d   (t=%.0f)"
          % (deg, len(R['C3']), len(R['C4']), R['r3'], R['r4'], len(R['C3']) - R['r3'], len(R['C4']) - R['r4'] - R['r3'],
             time.time() - t0), flush=True)
    if deg == (3, 3):
        R33 = R
H = load('h4w3_sl.pkl')
P_ = load('primitives_sl.pkl')
B4 = R33['B4']
rB = R33['r3']
for nm, c in (("theta'", H['thp']), ("{D, iota}", H['br']), ("o_3' = {D,iota} + theta'", H['o3p'])):
    vc = co[4].vec(c)
    closed = not co[5].vec(ce.d_ce(c, 4, sl))
    r = rank(B4 + [c], co[4])
    print("  %s at N=%d: nonzero %s, closed %s, rank[B^4 | c] = %d = rank B^4 (%d): exact %s"
          % (nm, N, bool(vc), closed, r, rB, r == rB), flush=True)
diff = ce.lin((1, ce.d_ce(P_['xi'], 3, sl)), (-1, H['o3p']))
print("  formal primitive xi at N=%d: d xi - o_3' vanishes at N: %s" % (N, not co[4].vec(diff)))
print("time %.0fs" % (time.time() - t0))
