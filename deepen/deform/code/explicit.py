"""Faithful evaluation of trace-basis cochains at a given N (gl_N or sl_N).

For arity k the trace monomials (fixed-point-free ones for sl_N) are evaluated either on all
entries of the explicit tensors on the basis of lie.py (mode 'tensor') or at random integer
matrix tuples (mode 'samples').  A set P of evaluation points with rank E_P = rank E is chosen;
then c -> (sum_tau c_{tau,mono} E[tau, P])_{mono, P} is injective on functions.
'tensor': rank E exact (certified).  'samples': rank E_P is a lower bound; it is exact when it equals
#taus - (number of exactly verified independent trace identities) (antisymmetrizer identities)."""
import sys, random
import numpy as np
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/deform/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
import ce, la, spaces
import lie, formal_gl as fg


class Evaluator:
    def __init__(self, N, sl, k, mode='tensor', nsamp=None, seed=1, verbose=True):
        self.N, self.sl, self.k = N, sl, k
        self.taus = spaces.taus(k, sl)
        self.tidx = {t: i for i, t in enumerate(self.taus)}
        if mode == 'tensor':
            g = lie.sl(N) if sl else lie.gl(N)
            rows = [fg.tau_tensor(t, g, k).reshape(-1) for t in self.taus]
            E = np.array(rows, dtype=np.int64)          # taus x entries
            keep = np.nonzero(np.any(E != 0, axis=0))[0]
            E = E[:, keep]
            Ml = [[int(x) for x in col] for col in E.T]   # entries x taus  (rows x cols)
            cert = la.certified_rank(Ml)
            self.rank = cert['rank']
            self.relations = cert['kernel']                # exact left-kernel vectors (tau coefficients)
            P = cert['pivrows']
            self.EP = [[int(E[i, p]) for p in P] for i in range(len(self.taus))]
            self.exact = True
            self.note = "tensor entries %d, rank %d (certified), relations %d" % (E.shape[1], self.rank, len(self.relations))
        else:
            rng = random.Random(seed)
            pts = []
            for s in range(nsamp):
                mats = []
                for _ in range(k):
                    A = [[rng.randint(-4, 4) for _ in range(N)] for _ in range(N)]
                    if sl:
                        A[N - 1][N - 1] -= sum(A[i][i] for i in range(N))
                    mats.append(np.array(A, dtype=object))
                pts.append(mats)
            E = [[self._tau_val(t, m) for m in pts] for t in self.taus]
            Ml = [[E[i][s] for i in range(len(self.taus))] for s in range(nsamp)]
            cert = la.certified_rank(Ml)
            self.rank = cert['rank']
            self.relations = cert['kernel']
            P = cert['pivrows']
            self.EP = [[E[i][p] for p in P] for i in range(len(self.taus))]
            self.exact = None   # decided by the caller (trace identities)
            self.note = "samples %d, rank %d (lower bound), sample-kernel %d" % (nsamp, self.rank, len(self.relations))
        if verbose:
            print("  [eval N=%d %s k=%d] %s" % (N, "sl" if sl else "gl", k, self.note), flush=True)

    @staticmethod
    def _tau_val(tau, mats):
        N = mats[0].shape[0]
        v = 1
        for c in tau:
            M = np.eye(N, dtype=object)
            for i in c:
                M = M.dot(mats[i])
            v *= int(sum(M[i, i] for i in range(N)))
        return v

    def vec(self, c, N=None):
        """compressed vector of a cochain (Npoly or Fraction coefficients): dict (mono, j) -> Fraction."""
        N = self.N if N is None else N
        out = {}
        for (t, m), v in c.items():
            if self.sl and ce.has_singleton(t):
                continue
            x = ce.np_eval(v, N) if isinstance(v, dict) else Fr(v)
            if x == 0:
                continue
            row = self.EP[self.tidx[t]]
            for j, e in enumerate(row):
                if e:
                    key = (m, j)
                    out[key] = out.get(key, 0) + x * e
        return {k_: x for k_, x in out.items() if x != 0}


def antisym_relation(k, sl):
    """sum_sigma sgn(sigma) tau_sigma (trace identity for N < k), restricted for sl."""
    rel = {}
    for s, sg in ce.perms(k):
        seen, cyc = set(), []
        for i in range(k):
            if i in seen:
                continue
            c, j = [], i
            while j not in seen:
                seen.add(j)
                c.append(j)
                j = s[j]
            cyc.append(c)
        t = ce.canon_tau(cyc)
        if sl and ce.has_singleton(t):
            continue
        rel[t] = rel.get(t, 0) + sg
    return rel


def check_relation_on_samples(rel, N, sl, k, nsamp=6, seed=5):
    rng = random.Random(seed)
    ok = True
    for _ in range(nsamp):
        mats = []
        for _ in range(k):
            A = [[rng.randint(-6, 6) for _ in range(N)] for _ in range(N)]
            if sl:
                A[N - 1][N - 1] -= sum(A[i][i] for i in range(N))
            mats.append(np.array(A, dtype=object))
        ok &= sum(c * Evaluator._tau_val(t, mats) for t, c in rel.items()) == 0
    return ok
