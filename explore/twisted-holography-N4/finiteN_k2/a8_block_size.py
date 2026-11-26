"""Size of one skeleton block of the k=2, N=2 complex in the O_0 sector: words gamma^1 x4, gamma^2 x4, gamma^3 x4
(each a separate even word of gh 0) plus n0 empty-word mesons M_0, at charge T = (2,2;2,2) and antifield number A
(g = -A).  Configurations: Sym^4(F)^{(x)3} (x) Sym^{n0}(F_0), F = {M,P,MM,R}^{ab}, F_0 = {M^{ab}}."""
from itertools import product
gens = []   # (parity, A, dLc(2), dRc(2))
for t in range(4):
    for a in range(2):
        for b in range(2):
            Lc = [0, 0]; Rc = [0, 0]
            if t == 0: Lc[a] += 1; Rc[b] += 1
            elif t == 1: Lc[a] += 1; Lc[b] -= 1
            elif t == 2: Rc[a] -= 1; Rc[b] += 1
            else: Rc[a] -= 1; Lc[b] -= 1
            gens.append(((0, 1, 1, 0)[t], (0, 1, 1, 2)[t], tuple(Lc), tuple(Rc)))
def sym(m):
    """generating dict (A, Lc, Rc) -> count for Sym^m(F) (odd generators exterior)."""
    S = {(0, (0, 0), (0, 0), 0): 1}     # key (A, Lc, Rc, degree)
    for (par, A, L, R) in gens:
        new = dict(S)
        for (a, lc, rc, deg), c in S.items():
            j = 1
            while deg + j <= m and (j <= 1 or not par):
                key = (a + j * A, (lc[0] + j * L[0], lc[1] + j * L[1]), (rc[0] + j * R[0], rc[1] + j * R[1]), deg + j)
                new[key] = new.get(key, 0) + c
                j += 1
        S = new
    return {(a, lc, rc): c for (a, lc, rc, deg), c in S.items() if deg == m}
def mul(X, Y):
    Z = {}
    for (a, l, r), c in X.items():
        for (a2, l2, r2), c2 in Y.items():
            key = (a + a2, (l[0] + l2[0], l[1] + l2[1]), (r[0] + r2[0], r[1] + r2[1]))
            Z[key] = Z.get(key, 0) + c * c2
    return Z
S4 = sym(4)
B = mul(mul(S4, S4), S4)
for A in range(0, 9):
    tot = 0
    for (a, l, r), c in B.items():
        if a != A: continue
        # M_0 completion: n0 = 4 - (Lc total); contingency count with row sums 2-l, col sums 2-r
        rr = (2 - l[0], 2 - l[1]); cc = (2 - r[0], 2 - r[1])
        if min(rr + cc) < 0 or sum(rr) != sum(cc): continue
        n = 0
        for m11 in range(0, min(rr[0], cc[0]) + 1):
            m12 = rr[0] - m11; m21 = cc[0] - m11; m22 = rr[1] - m21
            if m12 >= 0 and m21 >= 0 and m22 >= 0 and m12 + m22 == cc[1]: n += 1
        tot += c * n
    print(f"A={A} (g={-A}, Y={12-A}): block size {tot}")
