"""Order-hbar deformations J_1(x_0,x_1,x_2) = [q_0 q_1] T(X_0,X_1,X_2), T in (S^3 g^*)^g (setting (i)).
 - (S^3 g^*)^g: exact nullspace (sympy) of the invariance system on h-weight-zero multisets; full invariance of
   the resulting tensor re-checked on all index values.
 - J_1 alternating: [q_s0 q_s1] = sgn(s) [q_0 q_1] mod q_0+q_1+q_2 = 0 (checked for all s in S_3).
 - d J_1 (x_0..x_3) = sum_{i<j} (-1)^{i+j} [p_i+p_j, p_k] T([X_i,X_j], X_k, X_l), coefficients over (u, v, w),
   evaluated at X_0 = b_theta and all basis triples (injective on invariant tensors); exact integers.
 - control: the same formula with a random symmetric (non-invariant) T."""
import sys, itertools, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-g/code')
import numpy as np, sympy as sp
import algebras as al, core
from algebras import perm_sign


def sym3_invariants(L):
    d, wt = L.dim, L.wt
    zero = tuple(0 for _ in L.cartan)
    add = lambda *ws: tuple(sum(x) for x in zip(*ws))
    idxs = [t for t in itertools.combinations_with_replacement(range(d), 3) if add(*[wt[i] for i in t]) == zero]
    unk = {t: n for n, t in enumerate(idxs)}
    rows = set()
    for y in range(d):
        for t in itertools.combinations_with_replacement(range(d), 3):
            if add(wt[y], *[wt[i] for i in t]) != zero:
                continue
            row = {}
            for i in range(3):
                for c in np.nonzero(L.fu[y, t[i]])[0]:
                    tt = list(t); tt[i] = int(c)
                    n = unk[tuple(sorted(tt))]
                    row[n] = row.get(n, 0) + int(L.fu[y, t[i], c])
            row = tuple(sorted((k, v) for k, v in row.items() if v))
            if row:
                rows.add(row)
    M = sp.zeros(len(rows), len(idxs))
    for r, row in enumerate(sorted(rows)):
        for k, v in row:
            M[r, k] = v
    ns = M.nullspace()
    out = []
    for vec in ns:
        den = sp.ilcm(*[sp.Rational(x).q for x in vec])
        T = np.zeros((d, d, d), dtype=np.int64)
        for t, n in unk.items():
            val = int(vec[n] * den)
            for p in set(itertools.permutations(t)):
                T[p] = val
        out.append(T)
    return out, len(idxs), len(rows)


def check_invariant3(L, T):
    fu = L.fu
    r = (np.einsum('yae,ebc->yabc', fu, T) + np.einsum('ybe,aec->yabc', fu, T) + np.einsum('yce,abe->yabc', fu, T))
    return not r.any()


def alt_check():
    # [q_a q_b] over the single basis [q0 q1] with q2 = -(q0+q1)
    base = {(0, 1): 1, (0, 2): -1, (1, 2): 1}
    b = lambda a, c: 0 if a == c else (base[(a, c)] if a < c else -base[(c, a)])
    return all(b(s[0], s[1]) == perm_sign(s) for s in itertools.permutations(range(3)))


def dJ1(L, T, theta):
    d = L.dim
    DD = core.ct('abe,ecd->abcd', L.fu, T)        # DD[a,b,c,d] = T([b_a,b_b], b_c, b_d)
    out = np.zeros((3, d, d, d), dtype=np.int64)
    for i, j in itertools.combinations(range(4), 2):
        k, l = [r for r in range(4) if r not in (i, j)]
        lin = core.br(i, k) + core.br(j, k)
        Tvar = np.transpose(DD, list(np.argsort([i, j, k, l])))[theta]
        for c in range(3):
            if lin[c]:
                out[c] += (-1) ** (i + j) * int(lin[c]) * Tvar
    return out


if __name__ == "__main__":
    print("J_1 alternating ([q_s0 q_s1] = sgn(s)[q0 q1] for all s in S_3):", alt_check())
    rng = np.random.default_rng(7)
    for key in sys.argv[1:]:
        t0 = time.time()
        L = al.build(key)
        Ts, nu, ne = sym3_invariants(L)
        line = "%-14s (S^3 g*)^g exact nullspace: unknowns %d, distinct equations %d, dim %d" % (L.name, nu, ne, len(Ts))
        for T in Ts:
            inv_ok = check_invariant3(L, T)
            sym_ok = all(np.array_equal(T, np.transpose(T, p)) for p in itertools.permutations(range(3)))
            dj = dJ1(L, T, L.theta)
            line += " ; T: symmetric %s, invariant (all entries) %s, nonzero %s ; d J_1 = 0: %s" % (sym_ok, inv_ok, T.any(), not dj.any())
        # control
        R = rng.integers(-2, 3, size=(L.dim,) * 3)
        R = sum(np.transpose(R, p) for p in itertools.permutations(range(3))).astype(np.int64)
        line += " ; control (random symmetric non-invariant T): d J = 0: %s" % (not dJ1(L, R, L.theta).any())
        print(line + "  (%.1fs)" % (time.time() - t0), flush=True)
